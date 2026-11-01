"""CM1 adversarial review tests, part 3 (reviewer nacht-s10).

Lifecycle and detection attacks on code_mode.RpcServer and
code_mode.runner_program_path. The server must leave nothing behind when
it closes (the builder's own docstring and commit message claim the
unlink), and the runner detection must not be fooled into starting an RPC
server for a plain bash command.
"""
import json
import os
import tempfile

from delfin.agent import code_mode


class _NoExecutor:
    def __call__(self, name, arguments, permissions, rule_hit):
        return json.dumps({"result": "ok"})


def _server(tmp_path, token="t", executor=None):
    # AF_UNIX paths are capped at ~108 bytes; a pytest tmp_path is often
    # too long. rpc_socket_path(None, ...) uses the cwd, so run from a
    # short temp dir (created fresh; nothing outside the gate's tmpfs).
    cwd = os.getcwd()
    short = tempfile.mkdtemp(prefix="cm1s10-", dir=tempfile.gettempdir())
    os.chdir(short)
    try:
        server = code_mode.RpcServer(executor or _NoExecutor())
        path = code_mode.rpc_socket_path(None, token)
        server.bind(path)
    finally:
        os.chdir(cwd)
    return server, path


def test_close_unlinks_the_socket_file(tmp_path):
    """RpcServer.close must remove the socket file, not leak it. Red on
    3c9a69ebb: close() only closes the socket; the file under
    .code_mode_rpc/ stays in the workspace after every Code Mode run."""
    server, path = _server(tmp_path)
    assert os.path.exists(path)
    server.close()
    assert not os.path.exists(path), (
        "socket file leaked by close(): " + path)


def test_close_twice_is_no_error(tmp_path):
    """A second close (hook teardown races) must not raise."""
    server, path = _server(tmp_path)
    server.close()
    server.close()


def test_runner_detection_does_not_match_plain_python(tmp_path):
    """A plain python invocation that merely mentions the module in a
    comment or argument must not be treated as a Code Mode run."""
    for cmd in (
        "grep -rn 'delfin.agent.code_mode' delfin/",
        "python -c 'import delfin.agent.code_mode'",
    ):
        assert code_mode.runner_program_path(cmd) is None, cmd


def test_runner_detection_stops_at_quote_and_dollar(tmp_path):
    """The path token must not swallow shell expansions: $VAR, backtick,
    $(...) must terminate it, else the hook starts a server whose socket
    path is derived from attacker-controlled text."""
    cases = {
        "python -m delfin.agent.code_mode run a$HOME.py": "a",
        "python -m delfin.agent.code_mode run a`id`.py": "a",
        "python -m delfin.agent.code_mode run a$(id).py": "a",
    }
    for cmd, expected in cases.items():
        assert code_mode.runner_program_path(cmd) == expected, cmd


def test_runner_detection_requires_run_subcommand():
    """runner_program_path with no 'run' subcommand is None (already
    covered by the builder, repeated here as an ALLOW-side pin: only the
    exact runner shape turns RPC on)."""
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode serve") is None
