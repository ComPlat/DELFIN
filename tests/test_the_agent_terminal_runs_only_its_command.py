"""The dashboard's agent terminal runs the one configured program and nothing
else -- never a shell or an interpreter.

The stock Jupyter terminal manager lets the POST body override the command,
the environment and the directory; the control below shows a request for
``bash`` getting bash. The DELFIN manager must give the same request the
configured CLI, in the configured directory, without the client's
environment and without the dashboard token.
"""

from __future__ import annotations

import json
import os
import sys

import pytest

pytest.importorskip("jupyter_server_terminals")

import terminado.management as TM  # noqa: E402
from jupyter_server_terminals.terminalmanager import TerminalManager  # noqa: E402

from delfin.dashboard import agent_terminal as CT  # noqa: E402

HOSTILE = {"shell_command": ["/bin/bash", "-c", "id"],
           "extra_env": {"LD_PRELOAD": "/tmp/evil.so"},
           "cwd": "/", "name": "../../x", "height": 30, "width": 100}


@pytest.fixture
def spawned(monkeypatch):
    seen = []

    class _Pty:
        def __init__(self, argv, env=None, cwd=None):
            seen.append({"argv": list(argv), "env": dict(env or {}),
                         "cwd": cwd})
            self.clients = []
            self.ptyproc = type("P", (), {"pid": 1, "fd": -1})()
            self.read_buffer = []

        def resize_to_smallest(self):
            pass

    monkeypatch.setattr(TM, "PtyWithClients", _Pty)
    monkeypatch.setattr(TM.TermManagerBase, "start_reading",
                        lambda self, p: None)
    return seen


def _program(tmp_path, name="agent-cli"):
    """An executable that is neither a shell nor an interpreter by name."""
    exe = tmp_path / "bin" / name
    exe.parent.mkdir(exist_ok=True)
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    return exe


@pytest.fixture
def configured(monkeypatch, tmp_path):
    argv = [str(_program(tmp_path))]
    monkeypatch.setenv(CT.ARGV_ENV, json.dumps(argv))
    monkeypatch.setenv(CT.CWD_ENV, str(tmp_path))
    monkeypatch.setenv("JUPYTER_TOKEN", "secret-dashboard-token")
    return argv, str(tmp_path)


def test_control_the_stock_manager_obeys_the_request(spawned, configured):
    TerminalManager(shell_command=["/bin/sh"]).new_named_terminal(**HOSTILE)
    assert spawned[0]["argv"] == ["/bin/bash", "-c", "id"]
    assert spawned[0]["env"].get("LD_PRELOAD") == "/tmp/evil.so"


def test_the_request_cannot_choose_command_env_or_directory(
        spawned, configured):
    argv, cwd = configured
    CT.OneCommandTerminalManager().new_named_terminal(**HOSTILE)
    got = spawned[0]
    assert got["argv"] == argv
    assert got["cwd"] == cwd
    assert "LD_PRELOAD" not in got["env"]
    assert "JUPYTER_TOKEN" not in got["env"]
    assert got["env"]["COLUMNS"] == "100" and got["env"]["LINES"] == "30"


def test_a_name_that_is_not_a_plain_word_is_replaced(spawned, configured):
    mgr = CT.OneCommandTerminalManager()
    name, _ = mgr.new_named_terminal(**HOSTILE)
    assert name != "../../x" and name.isalnum()


def test_without_a_configured_command_there_is_no_terminal(monkeypatch):
    monkeypatch.delenv(CT.ARGV_ENV, raising=False)
    with pytest.raises(RuntimeError):
        CT.OneCommandTerminalManager()
    monkeypatch.setenv(CT.ARGV_ENV, json.dumps(["agent-cli"]))  # not absolute
    with pytest.raises(RuntimeError):
        CT.OneCommandTerminalManager()


@pytest.mark.parametrize("command", [
    "/bin/bash", "/usr/bin/sh", "/usr/bin/zsh", sys.executable,
    "/usr/bin/python3.11", "/usr/bin/node", "/usr/bin/env", "/usr/bin/sudo",
    "/usr/bin/tmux"])
def test_a_shell_or_an_interpreter_is_never_the_terminal(monkeypatch, command):
    """The manager refuses them even when the launcher was bypassed."""
    monkeypatch.setenv(CT.ARGV_ENV, json.dumps([command]))
    assert CT.refusal(command)
    assert CT.configured_argv() == []
    with pytest.raises(RuntimeError):
        CT.OneCommandTerminalManager()


@pytest.mark.parametrize("command", ["agent-cli", "/opt/tools/codex",
                                     "/home/u/.local/bin/aider"])
def test_a_program_is_allowed(command):
    assert CT.refusal(command) == ""


def test_at_most_two_terminals(spawned, configured):
    mgr = CT.OneCommandTerminalManager()
    mgr.new_named_terminal()
    mgr.new_named_terminal()
    with pytest.raises(Exception):
        mgr.new_named_terminal()


def test_a_taken_name_is_not_reused(spawned, configured):
    mgr = CT.OneCommandTerminalManager()
    first, _ = mgr.new_named_terminal(name="abc")
    second, _ = mgr.new_named_terminal(name="abc")
    assert first == "abc" and second != "abc"


# -- the launcher --------------------------------------------------------------

def _launch(monkeypatch, tmp_path, extra):
    from delfin import cli_voila
    from delfin.agent import process_guard

    monkeypatch.setattr(process_guard, "protect",
                        lambda *a: (_ for _ in ()).throw(RuntimeError("x")))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("XDG_RUNTIME_DIR", str(tmp_path / "runtime"))
    root = tmp_path / "home"
    root.mkdir(exist_ok=True)
    nb = tmp_path / "delfin_dashboard.ipynb"
    nb.write_text('{"cells":[]}', encoding="utf-8")
    captured = {}

    class _Proc:
        def wait(self):
            return 0

    def _popen(cmd, env=None, **kw):
        captured["cmd"], captured["env"] = cmd, env
        return _Proc()

    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_find_notebook", lambda: str(nb))
    monkeypatch.setattr(cli_voila, "_prepare_voila_env", lambda open_browser: {})
    monkeypatch.setattr(cli_voila, "_get_voila_static_root", lambda: "/tmp/v")
    monkeypatch.setattr(cli_voila, "_select_port", lambda port: port)
    monkeypatch.setattr(cli_voila, "_wait_for_port",
                        lambda host, port, timeout=10.0: False)
    monkeypatch.setattr(cli_voila.subprocess, "Popen", _popen)
    monkeypatch.setattr(cli_voila.subprocess, "run",
                        lambda cmd, env, check: (captured.update(
                            cmd=cmd, env=env), type("R", (), {
                                "returncode": 0})())[1])
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", str(root))
    code = 0
    try:
        cli_voila.main(["--no-browser", "--port", "9001"] + extra)
    except SystemExit as exc:
        code = exc.code
    return code, captured


def _fake_cli(monkeypatch, tmp_path):
    exe = _program(tmp_path)
    monkeypatch.setenv("PATH", f"{exe.parent}:/usr/bin:/bin")
    return exe


def test_without_the_flag_there_are_no_terminals(monkeypatch, tmp_path):
    _fake_cli(monkeypatch, tmp_path)
    code, cap = _launch(monkeypatch, tmp_path, [])
    assert code == 0
    assert "--ServerApp.terminals_enabled=False" in cap["cmd"]
    ext = next(a for a in cap["cmd"] if "jpserver_extensions" in a)
    assert "'jupyter_server_terminals': False" in ext
    assert CT.ARGV_ENV not in cap["env"]


def test_the_flag_wires_the_restricted_manager(monkeypatch, tmp_path):
    exe = _fake_cli(monkeypatch, tmp_path)
    code, cap = _launch(monkeypatch, tmp_path, ["--agent-terminal", "agent-cli"])
    assert code == 0
    assert "--ServerApp.terminals_enabled=True" in cap["cmd"]
    assert ("--TerminalsExtensionApp.terminal_manager_class="
            "delfin.dashboard.agent_terminal.OneCommandTerminalManager"
            ) in cap["cmd"]
    assert json.loads(cap["env"][CT.ARGV_ENV]) == [os.path.realpath(exe)]
    assert cap["env"][CT.ENABLED_ENV] == "1"


def test_the_flag_takes_a_path_as_well(monkeypatch, tmp_path):
    exe = _program(tmp_path, "my-agent")
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    code, cap = _launch(monkeypatch, tmp_path, ["--agent-terminal", str(exe)])
    assert code == 0
    assert json.loads(cap["env"][CT.ARGV_ENV]) == [os.path.realpath(exe)]


@pytest.mark.parametrize("command", ["sh", "bash", "python3", "env"])
def test_the_flag_refuses_a_shell_or_an_interpreter(monkeypatch, tmp_path,
                                                    command):
    code, cap = _launch(monkeypatch, tmp_path, ["--agent-terminal", command])
    assert code == 2
    assert CT.ARGV_ENV not in (cap.get("env") or {})


def test_a_link_to_a_shell_is_judged_by_where_it_points(monkeypatch, tmp_path):
    """A symlink named like a program that points at a shell is a shell."""
    shell = "/bin/sh" if os.path.exists("/bin/sh") else None
    if shell is None:
        shell = str(_program(tmp_path, "bash"))
    link = tmp_path / "bin" / "harmless-cli"
    link.parent.mkdir(exist_ok=True)
    try:
        link.symlink_to(shell)
    except OSError:
        # A host without the right to make links (Windows without developer
        # mode): the launcher's rule is still the resolved target's name.
        assert CT.refusal(os.path.realpath(shell))
        return
    code, _ = _launch(monkeypatch, tmp_path, ["--agent-terminal", str(link)])
    assert code == 2


def test_the_flag_is_refused_off_loopback(monkeypatch, tmp_path):
    _fake_cli(monkeypatch, tmp_path)
    code, _ = _launch(monkeypatch, tmp_path, [
        "--agent-terminal", "agent-cli", "--ip", "0.0.0.0",
        "--allow-remote-bind"])
    assert code == 2


def test_the_flag_is_refused_without_the_cli(monkeypatch, tmp_path):
    monkeypatch.setenv("PATH", "/usr/bin:/bin")
    monkeypatch.setenv("HOME", str(tmp_path / "nohome"))
    code, _ = _launch(monkeypatch, tmp_path, ["--agent-terminal", "agent-cli"])
    assert code == 2


# -- the tab -------------------------------------------------------------------

def _tab(tmp_path):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    scripts = []
    ctx.add_init_js = lambda js: scripts.append(js)
    tab = tab_agent.create_tab(ctx)
    return (tab[0] if isinstance(tab, tuple) else tab), scripts


def _walk(n):
    st = [n]
    while st:
        x = st.pop()
        yield x
        st.extend(getattr(x, "children", ()) or ())


def test_without_the_flag_the_tab_is_unchanged(tmp_path, monkeypatch):
    monkeypatch.delenv(CT.ENABLED_ENV, raising=False)
    root, scripts = _tab(tmp_path)
    classes = [set(c._dom_classes) for c in root.children]
    assert any("delfin-agent-chat-frame" in c for c in classes)
    assert not any("delfin-agent-chat-row" in c for c in classes)
    assert not any(getattr(w, "description", "") == "Terminal"
                   for w in _walk(root))
    assert not any("__delfinAgentTerm" in s for s in scripts)


def test_with_the_flag_the_panel_sits_beside_the_chat(tmp_path, monkeypatch):
    monkeypatch.setenv(CT.ENABLED_ENV, "1")
    monkeypatch.setenv(CT.CWD_ENV, str(tmp_path))
    monkeypatch.setenv(CT.ARGV_ENV, json.dumps([str(_program(tmp_path))]))
    root, scripts = _tab(tmp_path)
    row = next(c for c in root.children
               if "delfin-agent-chat-row" in c._dom_classes)
    frame, panel = row.children
    assert "delfin-agent-chat-frame" in frame._dom_classes
    assert "delfin-agent-term-panel" in panel._dom_classes
    head = panel.children[0]
    assert "agent-cli" in head.value              # the program, by name
    assert panel.layout.display == "none"            # closed until asked
    btn = next(w for w in _walk(root)
               if getattr(w, "description", "") == "Terminal")
    btn.click()
    assert panel.layout.display == ""
    btn.click()
    assert panel.layout.display == "none"
    js = next(s for s in scripts if "__delfinAgentTerm" in s)
    # Pinned and hash-checked; the page never names a command.
    assert CT.XTERM_JS_SRI in js and CT.XTERM_FIT_SRI in js
    assert "shell_command" not in js and "extra_env" not in js


def test_the_module_loads_without_the_terminals_package(monkeypatch):
    """tab_agent imports it on every dashboard; a missing optional package
    must not take the Agent tab with it."""
    import builtins
    import importlib

    real = builtins.__import__

    def _no_terminals(name, *a, **k):
        if name.startswith("jupyter_server_terminals"):
            raise ImportError(name)
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", _no_terminals)
    mod = importlib.reload(CT)
    try:
        assert mod.enabled() in (True, False)
        assert "__delfinAgentTerm" in mod.init_js()
        with pytest.raises(RuntimeError):
            mod.OneCommandTerminalManager()
    finally:
        monkeypatch.setattr(builtins, "__import__", real)
        importlib.reload(CT)
