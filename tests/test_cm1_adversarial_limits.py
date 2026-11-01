"""CM1 adversarial review tests, part 2 (reviewer nacht-s10).

Time limit and transport-cap surfaces of code_mode.run_program, exercised
in-process (the full caged subprocess path cannot nest bwrap in this cage;
the operator runs those outside). No wall-clock assertions on speed -- only
the limit's existence and its response shape.
"""
import json

import pytest

from delfin.agent import code_mode


def test_timeout_limit_raises_programtimeout(monkeypatch, tmp_path):
    """A runaway program hits the wall-clock limit -> ProgramTimeout with the
    limit named, not a hang."""
    import subprocess as _sp

    class _FakeRun:
        def __call__(self, argv, capture_output, timeout, env, cwd):
            raise _sp.TimeoutExpired(argv, timeout)

    monkeypatch.setattr(_sp, "run", _FakeRun())
    prog = tmp_path / "prog.py"
    prog.write_text("while True: pass\n")
    with pytest.raises(code_mode.ProgramTimeout, match="limit"):
        code_mode.run_program(str(prog), timeout_s=1.0)


def test_stdout_is_capped(monkeypatch, tmp_path):
    """A program that prints megabytes gets its stdout capped before the
    result reaches the model."""
    import subprocess as _sp

    class _FakeProc:
        returncode = 0
        stdout = b"x" * (code_mode.MAX_RESULT_BYTES + 123)
        stderr = b""

    monkeypatch.setattr(_sp, "run",
                        lambda *a, **k: _FakeProc())
    prog = tmp_path / "loud.py"
    prog.write_text("print('x' * 10)\n")
    out = code_mode.run_program(str(prog))
    assert "truncated" in out
    assert len(out) < 1000


def test_nonzero_exit_raises_programerror_with_stderr(tmp_path, monkeypatch):
    import subprocess as _sp

    class _FakeProc:
        returncode = 3
        stdout = b""
        stderr = b"boom\n"

    monkeypatch.setattr(_sp, "run", lambda *a, **k: _FakeProc())
    prog = tmp_path / "bad.py"
    prog.write_text("raise SystemExit(3)\n")
    with pytest.raises(code_mode.ProgramError, match="boom"):
        code_mode.run_program(str(prog))


def test_rpc_env_set_only_for_child(monkeypatch, tmp_path):
    """The socket path is exported only to the child process env, and the
    agent process itself never gains DELFIN_CODE_MODE_RPC."""
    seen = {}

    class _FakeProc:
        returncode = 0
        stdout = b"ok"
        stderr = b""

    def _fake_run(argv, capture_output, timeout, env, cwd):
        seen["env"] = env
        return _FakeProc()

    import subprocess as _sp
    monkeypatch.setattr(_sp, "run", _fake_run)
    monkeypatch.delenv("DELFIN_CODE_MODE_RPC", raising=False)
    code_mode.run_program(str(tmp_path / "p.py"), rpc_path="/w/.code_mode_rpc/rpc-t.sock",
                          workspace=str(tmp_path))
    assert seen["env"].get("DELFIN_CODE_MODE_RPC") == "/w/.code_mode_rpc/rpc-t.sock"
    assert "DELFIN_CODE_MODE_RPC" not in __import__("os").environ
