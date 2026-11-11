"""scripts/page_check.py without a browser: its arguments, the environment
it gives the dashboard, and its cleanup.

The browser part is checked by running the tool (docs/page_check.md).
What is pinned here is what made earlier ad-hoc checks expensive: a
dashboard left running for hours, a runtime path too long for a Unix
socket, a VS Code variable that opened a desktop browser tab.
"""

from __future__ import annotations

import importlib.util
import os
import signal
import stat
import subprocess
import sys
import time
from pathlib import Path

import pytest

_REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    "page_check", _REPO / "scripts" / "page_check.py")
page_check = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(page_check)


@pytest.fixture
def keep_signals():
    """Cleanup and the deadline install handlers; give pytest its own back."""
    saved = {s: signal.getsignal(s) for s in page_check._SIGNALS}
    yield
    signal.alarm(0)
    for s, h in saved.items():
        signal.signal(s, h)


def _alive(proc: subprocess.Popen) -> bool:
    return proc.poll() is None


def _sleeper(env_extra: dict, *, new_session: bool) -> subprocess.Popen:
    env = {"PATH": os.environ.get("PATH", "/usr/bin:/bin"), **env_extra}
    return subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(300)"],
        env=env, start_new_session=new_session)


# -- arguments -------------------------------------------------------------

def test_defaults_name_the_composer_and_this_checkout(tmp_path):
    args = page_check.parse_args(["--out", str(tmp_path / "o")])
    assert args.scenario == "composer"
    assert args.tree == _REPO
    assert args.out == (tmp_path / "o").resolve()
    assert args.timeout > 0


def test_service_session_is_a_scenario(tmp_path):
    args = page_check.parse_args(["--out", str(tmp_path), "--scenario",
                                  "service-session", "--timeout", "60"])
    assert (args.scenario, args.timeout) == ("service-session", 60)


@pytest.mark.parametrize("argv", [
    [],                                              # --out is required
    ["--scenario", "nope"],
    ["--timeout", "0"],
    ["--tree", "/"],                                 # not a checkout
])
def test_bad_arguments_exit_with_usage_error(tmp_path, argv, capsys):
    if argv:
        argv = ["--out", str(tmp_path)] + argv
    with pytest.raises(SystemExit) as exc:
        page_check.parse_args(argv)
    assert exc.value.code == 2


def test_server_python_comes_from_the_environment(tmp_path, monkeypatch):
    monkeypatch.setenv(page_check.SERVER_PYTHON_ENV, "/opt/x/bin/python")
    args = page_check.parse_args(["--out", str(tmp_path)])
    assert args.server_python == Path("/opt/x/bin/python")


# -- environment -----------------------------------------------------------

def test_dashboard_environment_inherits_nothing(tmp_path, monkeypatch):
    for name, value in {"VSCODE_IPC_HOOK_CLI": "/x", "BROWSER": "/x/browser.sh",
                        "TERM_PROGRAM": "vscode", "LD_LIBRARY_PATH": "/opt/lib",
                        "OPENAI_API_KEY": "not-a-real-key"}.items():
        monkeypatch.setenv(name, value)
    env = page_check.server_env(
        server_python=Path("/env/bin/python"), tree=tmp_path, home=tmp_path / "h",
        scratch=tmp_path / "s", runtime=Path("/tmp/dpc-x"), workspace=tmp_path / "w",
        model_port=1234)
    for name in ("VSCODE_IPC_HOOK_CLI", "BROWSER", "TERM_PROGRAM",
                 "LD_LIBRARY_PATH", "OPENAI_API_KEY"):
        assert name not in env
    assert env["PATH"].split(":")[0] == "/env/bin"
    assert env["HOME"] == str(tmp_path / "h")
    assert env["DELFIN_RUNTIME_DIR"] == "/tmp/dpc-x"
    assert env["OLLAMA_HOST"] == "http://127.0.0.1:1234"
    assert env["DELFIN_NO_BROWSER"] == "1"
    assert env["ANTHROPIC_API_KEY"] == page_check.PLACEHOLDER_KEY


def test_runtime_dir_is_short_and_private():
    path = page_check.make_runtime_dir()
    try:
        assert stat.S_IMODE(path.stat().st_mode) == 0o700
        sock = page_check.service_socket(path)
        assert len(str(sock).encode()) <= page_check.MAX_SOCKET_PATH
    finally:
        path.rmdir()


def test_a_runtime_dir_too_long_for_a_socket_is_refused(tmp_path):
    base = tmp_path / ("d" * 120)
    base.mkdir()
    with pytest.raises(ValueError, match="too long"):
        page_check.make_runtime_dir(str(base))
    assert list(base.iterdir()) == []


# -- cleanup ---------------------------------------------------------------

def test_cleanup_ends_the_group_and_what_left_it(tmp_path, keep_signals):
    runtime = page_check.make_runtime_dir()
    work = tmp_path / "work"
    work.mkdir()
    # The dashboard: its own process group.
    dash = _sleeper({}, new_session=True)
    # A kernel or service that left the group but carries the runtime dir.
    stray = _sleeper({"DELFIN_RUNTIME_DIR": str(runtime)}, new_session=True)
    # Not ours: another runtime dir. Must survive.
    other = _sleeper({"DELFIN_RUNTIME_DIR": str(runtime) + "-other"},
                     new_session=True)
    try:
        time.sleep(0.5)          # let the children exec, so /proc shows env
        c = page_check.Cleanup()
        c.add_group(dash)
        c.runtime, c.work = runtime, work
        rep = c.run()
        stray.wait(timeout=10)
        assert not _alive(dash) and not _alive(stray)
        assert _alive(other)
        assert rep["groups"] == ["terminated"]
        assert rep["leftovers_ended"] == 1
        assert rep["service"] == "none"
        assert rep["runtime_dir_removed"] is True
        assert not runtime.exists() and not work.exists()
        assert c.run() is rep    # idempotent: the second call does nothing
    finally:
        for p in (dash, stray, other):
            if _alive(p):
                p.kill()
            p.wait(timeout=10)


def test_a_group_that_ignores_sigterm_is_killed(keep_signals):
    proc = subprocess.Popen(
        [sys.executable, "-c",
         "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN);"
         " print('ready', flush=True); time.sleep(300)"],
        stdout=subprocess.PIPE, start_new_session=True)
    try:
        assert proc.stdout.readline().strip() == b"ready"
        assert page_check.Cleanup(grace=1)._end_group(proc) == "killed"
        assert not _alive(proc)
    finally:
        if _alive(proc):
            proc.kill()
        proc.wait(timeout=10)


def test_the_deadline_raises_into_the_main_thread_and_is_undone(keep_signals):
    before = signal.getsignal(signal.SIGINT)
    previous = page_check.install_deadline(1)
    with pytest.raises(page_check.Deadline, match="SIGALRM"):
        time.sleep(5)
    page_check.restore_signals(previous)
    assert signal.getsignal(signal.SIGINT) is before
