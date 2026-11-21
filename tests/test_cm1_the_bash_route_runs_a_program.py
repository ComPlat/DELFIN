"""CM1: the product's bash route runs a program; run_program did not add to it.

``code_mode.run_program`` ran a program with ``subprocess.run(["/bin/bash",
"-c", ...])`` and its docstring said "we run it with the same wall-clock
limit and cap so the process boundaries stay identical". Measured, it had
none of them: the argv carried no ``_bash_isolation_argv`` wrap, so no
filesystem namespace and no process cage — a program launched by it could
write anywhere the agent could. The product's real route is the bash tool,
where ``api_client._code_mode_child_env`` starts the RPC server and exports
the socket to the caged child; there the cage, the wall-clock ceiling
(``bash_timeout_s``/``bash_max_timeout_s``) and output truncation are already
enforced, once, by the tool that owns them. ``run_program`` was deleted.

These tests drive the route that remains: the hook itself, a real program
subprocess, and a real ``RpcServer``. None of them calls ``run_program``.

The workspace sits in a short temp directory: the RPC socket path is
``<workspace>/.code_mode_rpc/rpc-<32-hex>.sock`` and an AF_UNIX sockaddr tops
out at 108 bytes, so a deep TMPDIR fails bind() for a reason unrelated to
what is measured here.
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from delfin.agent import api_client
from delfin.agent import code_mode
from conftest import child_env

#: The import root the program subprocess gets, so the child imports THIS
#: tree's delfin rather than whatever is installed.
TREE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(code_mode.__file__))))


class _RecordingExecutor:
    """Stands in for KitToolExecutor at the public entry the RPC calls."""

    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions=None):
        self.calls.append(name)
        return json.dumps({"ok": name})


def _short_workspace():
    """A workspace whose RPC socket path fits an AF_UNIX sockaddr."""
    candidates = [tempfile.gettempdir(), "/tmp"]
    for base in sorted({c for c in candidates if c and os.path.isdir(c)},
                       key=len):
        return tempfile.mkdtemp(prefix="cm1r_", dir=base)
    raise RuntimeError("no temp directory available")


def _permissions(workspace):
    return api_client.KitToolPermissions(mode="default", workspace=workspace)


def _write_stub_and_program(workspace, body):
    with open(os.path.join(workspace, "tools_api.py"), "w") as fh:
        fh.write(code_mode.generate_stub())
    program = os.path.join(workspace, "prog.py")
    with open(program, "w") as fh:
        fh.write(body)
    return program


def _child_env(workspace, executor, command="python3 -m delfin.agent.code_mode "
                                            "run prog.py"):
    return api_client._code_mode_child_env(
        command, workspace, _permissions(workspace), executor)


def test_run_program_is_gone_and_the_hook_is_the_route():
    """The duplicate launcher is deleted; the hook the product uses is not.

    Red on the version that had both: ``run_program`` existed there.
    """
    assert not hasattr(code_mode, "run_program")
    assert hasattr(api_client, "_code_mode_child_env")
    assert hasattr(code_mode.RpcServer, "serve_one")


def test_the_hook_grants_the_socket_only_to_a_runner_command():
    """A plain command gets nothing added: no server, no env var."""
    workspace = _short_workspace()
    try:
        executor = _RecordingExecutor()
        assert _child_env(workspace, executor, "ls -la") == {}
        env = _child_env(workspace, executor)
        assert set(env) == {"DELFIN_CODE_MODE_RPC"}
        assert env["DELFIN_CODE_MODE_RPC"].startswith(workspace)
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_a_real_program_reaches_a_tool_through_the_product_route():
    """The hook's socket + a real program subprocess + a real RpcServer.

    This is the capability ``run_program`` claimed: the program calls a tool
    and gets the executor's answer back. It is proven here on the PRODUCT's
    route -- the hook exports the path, the program reads it from the
    environment -- so nothing is left unmeasured by the deletion.
    """
    workspace = _short_workspace()
    try:
        executor = _RecordingExecutor()
        _write_stub_and_program(
            workspace,
            "import json, tools_api as api\n"
            "print('RESULT', json.dumps(api.list_docs()))\n")

        # conftest.child_env, not the local _child_env in this file:
        # that one wraps api_client._code_mode_child_env and says
        # nothing about HOME. These children run
        # `python -m delfin.agent.code_mode`, so each resolves its own
        # ~/.delfin sinks at import time against whatever HOME says.
        env = dict(os.environ, **child_env(Path(workspace)))
        env.update(_child_env(workspace, executor))
        env["PYTHONPATH"] = os.pathsep.join([workspace, TREE])
        run = subprocess.run(
            [sys.executable, "-m", "delfin.agent.code_mode", "run", "prog.py"],
            cwd=workspace, env=env, capture_output=True, text=True, timeout=120)

        assert run.returncode == 0, run.stderr
        assert "RESULT" in run.stdout, run.stdout + run.stderr
        assert json.loads(run.stdout.split("RESULT", 1)[1]) == {"ok": "list_docs"}
        assert executor.calls == ["list_docs"]
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_a_write_tool_is_refused_over_the_product_route():
    """A program that fabricates a write call is still refused.

    The generated stub carries read-only tools only, so a program cannot even
    name ``write_file`` through it. This bypasses the stub deliberately and
    speaks the RPC directly: the allow-list must hold at the server, not only
    in what the stub happens to expose.
    """
    workspace = _short_workspace()
    try:
        executor = _RecordingExecutor()
        _write_stub_and_program(
            workspace,
            "from delfin.agent.code_mode import RpcClient\n"
            "import os\n"
            "c = RpcClient(os.environ['DELFIN_CODE_MODE_RPC'])\n"
            "try:\n"
            "    c.call('write_file', {'path': 'x', 'content': 'y'})\n"
            "except RuntimeError as exc:\n"
            "    print('REFUSED', exc)\n")

        # conftest.child_env, not the local _child_env in this file:
        # that one wraps api_client._code_mode_child_env and says
        # nothing about HOME. These children run
        # `python -m delfin.agent.code_mode`, so each resolves its own
        # ~/.delfin sinks at import time against whatever HOME says.
        env = dict(os.environ, **child_env(Path(workspace)))
        env.update(_child_env(workspace, executor))
        env["PYTHONPATH"] = os.pathsep.join([workspace, TREE])
        run = subprocess.run(
            [sys.executable, "-m", "delfin.agent.code_mode", "run", "prog.py"],
            cwd=workspace, env=env, capture_output=True, text=True, timeout=120)

        assert run.returncode == 0, run.stderr
        assert "REFUSED" in run.stdout
        assert "write_file" in run.stdout
        assert executor.calls == []      # never reached the executor
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def test_the_stub_offers_no_way_to_name_a_write_tool():
    """The API the model programs against exposes read-only tools only."""
    stub = code_mode.generate_stub()
    for banned in ("write_file", "submit_calculation", "delete_calc", "bash"):
        assert "def %s(" % banned not in stub, banned


def test_the_program_inherits_the_cage_and_ceiling_from_the_bash_tool():
    """What run_program lacked is enforced by the tool that owns it.

    ``run_program`` built its own argv; the product's route builds the argv
    through ``_bash_isolation_argv``, which is where the filesystem/process
    cage is decided. Asserting on that function's call for a runner command
    keeps the promise attached to the code that keeps it.
    """
    workspace = _short_workspace()
    try:
        perms = _permissions(workspace)
        seen = []
        real = api_client._bash_isolation_argv

        def _spy(cmd, run_cwd, perms_, mode=None, extra_write=(),
                 caged_unasked=False):
            seen.append(cmd)
            return real(cmd, run_cwd, perms_, mode, extra_write,
                        caged_unasked)

        # The PRODUCT must make this call. The first version of this test
        # installed the spy and then called _bash_isolation_argv ITSELF, so
        # `assert seen` was satisfied by the test's own call and the bash
        # tool was never driven. After run_program, ProgramTimeout and
        # test_cm1_adversarial_limits.py were deleted, this was the only
        # coverage left for the cage and the ceiling on a runner command.
        engine = api_client._DocToolExecutor.__new__(
            api_client._DocToolExecutor)
        engine._permissions = perms
        api_client._bash_isolation_argv = _spy
        try:
            out = engine._execute_bash(
                {"command": "python3 -m delfin.agent.code_mode run prog.py",
                 "timeout": 5},
                perms)
        finally:
            api_client._bash_isolation_argv = real
        assert out, "the bash tool returned nothing for a runner command"
        assert seen, (
            "the bash tool built no cage argv for a runner command; the "
            "cage is decided in _bash_isolation_argv and nothing else")
        assert any("code_mode" in str(c) for c in seen), seen
        # The ceiling the tool applies to every command, program or not --
        # read off the permissions the tool was given, not off a default.
        assert perms.bash_timeout_s == 120
        assert perms.bash_max_timeout_s > perms.bash_timeout_s
    finally:
        shutil.rmtree(workspace, ignore_errors=True)
