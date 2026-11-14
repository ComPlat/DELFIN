"""``delfin-voila`` typed at a terminal starts the server detached from it.

Closing the terminal, or an SSH connection dropping, used to end the
dashboard: SIGHUP took the Ctrl+C path (cli_voila._install_stop_signals),
and surviving it meant tmux and its keys. Now a launch typed at a
terminal starts the same launcher again in a session of its own, with
its output in a log, waits until the server answers, prints where it is
and returns the prompt.

Pinned here:

  who detaches         a terminal on stdin, not inside tmux/screen or a
                       SLURM job, not a branch preview, no --foreground,
                       no DELFIN_VOILA_FOREGROUND; everything else stays
                       in the foreground (tests/test_a_launch_that_must_
                       stay_in_the_foreground_does.py)
  how                  its own session (start_new_session), stdin from
                       /dev/null, stdout and stderr appended to an
                       owner-only log, the launch folder as working
                       folder, no lifeline of the terminal's
  the token            handed over in the environment, never on a
                       command line (/proc/<pid>/cmdline is world-readable)
  the answer           the URL once the server answers, or the log's
                       last lines and the exit code when it did not start
"""

from __future__ import annotations

import os
import stat
import sys

import pytest

from delfin import cli_voila
from delfin.dashboard import server_record as SR


class _Stdin:
    def __init__(self, tty: bool):
        self._tty = tty

    def isatty(self) -> bool:
        return self._tty


def _args(*argv):
    """Parsed the way main() parses them."""
    return cli_voila._parse_args(list(argv))


@pytest.fixture
def plain_terminal(monkeypatch):
    for name in ("TMUX", "STY", "SLURM_JOB_ID", "DELFIN_VOILA_FOREGROUND",
                 "DELFIN_PREVIEW_BRANCH", "DELFIN_PREVIEW_WORKTREE",
                 SR.DETACHED_ENV):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(sys, "stdin", _Stdin(True))


def test_a_launch_typed_at_a_terminal_detaches(plain_terminal):
    detach, why = cli_voila._should_detach(_args())
    assert detach is True, why


@pytest.mark.parametrize("env, reason", [
    ({"TMUX": "/tmp/tmux-1/default,1,0"}, "tmux"),
    ({"STY": "1234.pts-0.host"}, "screen"),
    ({"SLURM_JOB_ID": "4711"}, "SLURM"),
    ({"DELFIN_VOILA_FOREGROUND": "1"}, "DELFIN_VOILA_FOREGROUND"),
    ({"DELFIN_PREVIEW_BRANCH": "b", "DELFIN_PREVIEW_WORKTREE": "/w"}, "preview"),
    ({SR.DETACHED_ENV: "1"}, "detached"),
])
def test_launches_that_own_their_lifetime_stay(plain_terminal, monkeypatch,
                                               env, reason):
    for key, value in env.items():
        monkeypatch.setenv(key, value)
    detach, why = cli_voila._should_detach(_args())
    assert detach is False
    assert reason.lower() in why.lower()


def test_without_a_terminal_it_stays(plain_terminal, monkeypatch):
    monkeypatch.setattr(sys, "stdin", _Stdin(False))
    assert cli_voila._should_detach(_args())[0] is False


def test_foreground_is_asked_for_by_name(plain_terminal):
    detach, why = cli_voila._should_detach(_args("--foreground"))
    assert detach is False and "--foreground" in why


def test_the_child_gets_no_token_and_no_browser_on_its_command_line():
    token = "s" * 43
    argv = cli_voila._detached_child_argv(
        ["--port", "9000", "--token", token, "--open-browser", "--keep",
         "--dark", f"--token={token}"])
    joined = " ".join(argv)
    assert token not in joined
    assert "--open-browser" not in argv
    assert argv.count("--foreground") == 1
    assert "--no-browser" in argv
    assert ["--port", "9000"] == argv[argv.index("--port"):argv.index("--port") + 2]
    assert "--keep" in argv and "--dark" in argv


def test_the_child_environment_carries_the_token_and_no_lifeline(monkeypatch):
    monkeypatch.setenv("DELFIN_LIFELINE_PID", "1")
    monkeypatch.setenv("DELFIN_LIFELINE_TICKS", "2")
    monkeypatch.setenv("JUPYTER_TOKEN", "from-somewhere-else")
    env = cli_voila._detached_child_env("t" * 43)
    assert env["DELFIN_VOILA_TOKEN"] == "t" * 43
    assert env[SR.DETACHED_ENV] == "1"
    assert "DELFIN_LIFELINE_PID" not in env
    assert "DELFIN_LIFELINE_TICKS" not in env
    assert "JUPYTER_TOKEN" not in env
    assert env.get("DELFIN_LAUNCH_CWD") == os.getcwd()


class _Child:
    def __init__(self, pid, code=None):
        self.pid = pid
        self._code = code

    def poll(self):
        return self._code


def _record(port=8899, **extra):
    rec = SR.own_record()
    rec.update({"port": port, "ip": "127.0.0.1", "scheme": "http",
                "url": "/voila/render/delfin_voila_runtime/delfin_dashboard.ipynb",
                "resume_path": "/voila/render/delfin_voila_runtime/delfin_resume.ipynb",
                "token": "k" * 43, "started_at": 0.0, "detached": True,
                "log": str(SR.log_file())})
    rec.update(extra)
    return rec


def test_start_detached_starts_a_session_of_its_own_and_waits(
        monkeypatch, capsys):
    seen = {}

    def _popen(cmd, **kwargs):
        seen["cmd"] = list(cmd)
        seen["kwargs"] = kwargs
        seen["log_mode"] = stat.S_IMODE(os.fstat(kwargs["stdout"]).st_mode)
        # The child writes the record once its server answers; here this
        # process stands in for the child (pid and start time are real).
        SR.write(_record())
        return _Child(os.getpid())

    monkeypatch.setattr(cli_voila.subprocess, "Popen", _popen)
    monkeypatch.setattr(SR, "port_answers", lambda record, timeout=0.5: True)
    monkeypatch.setattr(cli_voila, "_open_browser_url",
                        lambda url, env: seen.setdefault("browser", url))
    args = _args("--no-browser")
    rc = cli_voila._start_detached(args, ["--no-browser"])
    assert rc == 0
    kwargs = seen["kwargs"]
    assert kwargs["start_new_session"] is True
    assert kwargs["stdin"] == cli_voila.subprocess.DEVNULL
    assert kwargs["stderr"] == kwargs["stdout"]
    assert kwargs["cwd"] == os.getcwd()
    assert seen["log_mode"] == 0o600
    assert "cli_voila" in " ".join(seen["cmd"])
    assert "--foreground" in seen["cmd"]
    assert "k" * 43 not in " ".join(seen["cmd"])
    assert "token" not in " ".join(seen["cmd"]).lower().replace("no-token", "")
    assert kwargs["env"]["DELFIN_VOILA_TOKEN"]
    out = capsys.readouterr().out
    assert "localhost:8899" in out
    assert "delfin-voila stop" in out
    assert "closing this terminal" in out.lower()
    assert "browser" not in seen, "--no-browser opens nothing"


def test_a_child_that_does_not_start_shows_its_log(monkeypatch, capsys):
    def _popen(cmd, **kwargs):
        os.write(kwargs["stdout"], b"Error: the port is taken, for the test\n")
        return _Child(999999, code=3)

    monkeypatch.setattr(cli_voila.subprocess, "Popen", _popen)
    rc = cli_voila._start_detached(_args("--no-browser"), ["--no-browser"])
    assert rc == 3
    err = capsys.readouterr().err
    assert "the port is taken, for the test" in err
    assert str(SR.log_file()) in err


def test_main_at_a_terminal_hands_over_to_the_detached_start(
        plain_terminal, monkeypatch):
    calls = []
    monkeypatch.setattr(cli_voila, "_voila_is_available", lambda: True)
    monkeypatch.setattr(cli_voila, "_start_detached",
                        lambda args, argv: calls.append(list(argv)) or 0)
    monkeypatch.setattr(cli_voila.subprocess, "Popen",
                        lambda *a, **k: pytest.fail("no server from here"))
    from delfin.agent import process_guard
    monkeypatch.setattr(process_guard, "protect",
                        lambda *a: (_ for _ in ()).throw(RuntimeError("off")))
    assert cli_voila.main(["--no-browser", "--port", "9100"]) == 0
    assert calls == [["--no-browser", "--port", "9100"]]


def test_the_child_bootstrap_does_not_import_from_the_launch_folder():
    """``python -c`` puts the working folder first on sys.path, and the
    child runs in the launch folder: a json.py there must not be
    imported in place of the standard one."""
    cmd = cli_voila._detached_child_command(["--foreground"])
    assert cmd[0] == sys.executable and cmd[1] == "-c"
    assert "sys.path[:] = [p for p in sys.path if p not in ('', '.')]" in cmd[2]
