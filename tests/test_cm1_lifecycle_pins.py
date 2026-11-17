"""Pins for the two findings the CM1 reviewer reported on 3c9a69ebb.

FINDING A: RpcServer.close must unlink the socket file (idempotent).
FINDING B: runner_program_path must stop at shell metacharacters including
$, ` ( ) ' — the hook derives the socket path / env var from the token.
Both were seen red on 3c9a69ebb by the reviewer (tests/test_cm1_adversarial_lifecycle.py:
F..F., 2 failed); committed here as permanent pins after the fix.

Note: server tests use a short mkdtemp cwd, not pytest tmp_path — tmp_path
can exceed the 108-byte AF_UNIX sockaddr limit on some filesystems.
"""
import os
import shutil
import tempfile

from delfin.agent import code_mode


def _short_cwd():
    base = tempfile.mkdtemp(prefix="cm1p1c_")
    return base


def test_close_unlinks_the_socket_file():
    base = _short_cwd()
    try:
        server = code_mode.RpcServer(lambda *a: "{}")
        sock = code_mode.rpc_socket_path(base, "tokA")
        server.bind(sock)
        assert os.path.exists(sock)
        server.close()
        assert not os.path.exists(sock)
    finally:
        shutil.rmtree(base, ignore_errors=True)


def test_close_twice_is_no_error():
    base = _short_cwd()
    try:
        server = code_mode.RpcServer(lambda *a: "{}")
        server.bind(code_mode.rpc_socket_path(base, "tokB"))
        server.close()
        server.close()
    finally:
        shutil.rmtree(base, ignore_errors=True)


def test_runner_detection_stops_at_expansions_and_quotes():
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run a$HOME.py") == "a"
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run a$(x).py") == "a"
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run a`x`.py") == "a"
    # A quoted path is unparseable token-wise; partial extraction would be
    # wrong, so the detection refuses it (None) — the hook then grants no env.
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run 'prog x.py'") is None
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run (prog).py") is None


def test_runner_detection_still_accepts_plain_paths():
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run prog_x.py") == "prog_x.py"
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run ./sub/dir/prog.py"
    ) == "./sub/dir/prog.py"
