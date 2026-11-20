"""CM1: the documented program route runs, over one launcher and one stub.

Three questions, each driven through the path the product uses rather than
through a convenience wrapper:

* Does a program reach a tool? ``api_client._code_mode_child_env`` is the hook
  ``KitToolExecutor._execute_bash`` calls for every bash command; it starts the
  ``RpcServer`` and exports the socket to the child. The program here is a real
  subprocess whose only route to a tool is that socket, so what is measured is
  the dispatch path (``RpcServer.handle_payload`` -> ``Connection.handle`` ->
  ``executor.execute``), not a function called in-process.
* Can the model obtain the API module? The Code Mode skill documents
  ``python -m delfin.agent.code_mode stub > tools_api.py``. Driven through
  ``main()``, the CLI entry a subprocess reaches.
* Does the module the CLI writes match the module the skill imports? The stub
  carries its own import instructions; a third name for the one file is a
  defect the model would hit and no test would catch.

``run_program`` and ``serve_request`` are pinned as gone: neither had a caller
outside the suite, and a second launcher under the bash tool skipped the cage,
the wall-clock ceiling and the output truncation that tool owns.

Workspaces live in a short temp directory: the RPC socket path is
``<workspace>/.code_mode_rpc/rpc-<32 hex>.sock`` and an AF_UNIX sockaddr tops
out at 108 bytes, so a deep TMPDIR would fail bind() for a reason that is not
what this file measures.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile

from delfin.agent import api_client, code_mode

#: The import root the program subprocess gets, so the child imports THIS tree
#: (``<tree>/delfin/agent/code_mode.py`` -> ``<tree>``).
TREE = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(code_mode.__file__))))


class _RecordingExecutor:
    """Stands in for KitToolExecutor at the public entry the RPC calls."""

    def __init__(self):
        self.calls = []

    def execute(self, name, arguments, permissions=None):
        self.calls.append((name, dict(arguments or {})))
        return json.dumps({"ok": name, "arguments": arguments})


def _short_workspace():
    """A workspace whose RPC socket path fits an AF_UNIX sockaddr."""
    bases = [b for b in (tempfile.gettempdir(), "/tmp") if b and os.path.isdir(b)]
    return tempfile.mkdtemp(prefix="cm1r_", dir=min(bases, key=len))


def _permissions(workspace):
    return api_client.KitToolPermissions(mode="default", workspace=workspace)


def _child_env(workspace, executor, command):
    """The env the bash tool would give a child for this command."""
    env = dict(os.environ)
    env.update(api_client._code_mode_child_env(
        command, workspace, _permissions(workspace), executor))
    env["PYTHONPATH"] = TREE + os.pathsep + env.get("PYTHONPATH", "")
    return env


def _write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def _run(program, workspace, env):
    return subprocess.run(
        [sys.executable, "-m", "delfin.agent.code_mode", "run", program],
        capture_output=True, text=True, cwd=workspace, env=env, timeout=120)


def test_a_program_reaches_a_tool_through_the_hook_the_bash_tool_uses():
    """The route the product has: hook -> socket -> executor.execute.

    Red before the fix only in that the model could not obtain the stub; the
    transport itself worked. Stated as coverage of the surviving route.
    """
    workspace = _short_workspace()
    _write(os.path.join(workspace, "tools_api.py"), code_mode.generate_stub())
    program = _write(
        os.path.join(workspace, "prog.py"),
        "import tools_api as api\n"
        "print(api.get_calc_info(calc_id='job_01'))\n")
    executor = _RecordingExecutor()
    env = _child_env(workspace, executor,
                     "python3 -m delfin.agent.code_mode run prog.py")

    assert "DELFIN_CODE_MODE_RPC" in env, "the hook exported no socket"
    proc = _run(program, workspace, env)

    assert proc.returncode == 0, proc.stderr
    assert "get_calc_info" in proc.stdout
    assert executor.calls == [("get_calc_info", {"calc_id": "job_01"})]


def test_a_write_tool_is_refused_at_the_server_before_the_executor():
    """The allow-list holds on the real socket, not only in-process."""
    workspace = _short_workspace()
    _write(os.path.join(workspace, "tools_api.py"), code_mode.generate_stub())
    program = _write(
        os.path.join(workspace, "prog.py"),
        "import tools_api as api\n"
        "try:\n"
        "    api._call('write_file', {'path': 'x', 'content': 'y'})\n"
        "except RuntimeError as exc:\n"
        "    print('REFUSED', exc)\n")
    executor = _RecordingExecutor()
    env = _child_env(workspace, executor,
                     "python3 -m delfin.agent.code_mode run prog.py")

    proc = _run(program, workspace, env)

    assert proc.returncode == 0, proc.stderr
    assert "REFUSED" in proc.stdout
    assert "write_file" in proc.stdout
    assert executor.calls == [], "a write tool reached the executor"


def test_the_documented_stub_command_writes_the_api_module(capsys):
    """``... stub`` puts the API module on stdout, exit 0.

    Fails before the fix: ``main()`` dispatched only ``run``, so this printed
    the ``run`` usage and exited 2 while the skill told the model to run it.
    """
    assert code_mode.main(["stub"]) == 0
    out = capsys.readouterr().out
    assert "def get_calc_info(" in out
    assert "def search_calcs(" in out


def test_the_written_module_is_the_one_the_skill_imports():
    """The stub the CLI writes works under the name the skill saves it as."""
    workspace = _short_workspace()
    stub = _write(os.path.join(workspace, "tools_api.py"),
                  code_mode.generate_stub())
    program = _write(
        os.path.join(workspace, "prog.py"),
        "import tools_api as api\n"
        "print(sorted(api.delfin_tools.__dict__ and [\n"
        "    n for n in dir(api.delfin_tools) if not n.startswith('_')]))\n")
    executor = _RecordingExecutor()
    env = _child_env(workspace, executor,
                     "python3 -m delfin.agent.code_mode run prog.py")

    proc = _run(program, workspace, env)

    assert proc.returncode == 0, proc.stderr
    printed = proc.stdout.strip()
    for name in code_mode.READ_ONLY_TOOLS:
        assert name in printed, (name, printed)
    assert os.path.exists(stub)


def test_the_stub_names_the_module_it_is_saved_as():
    """One name for the one file: no third import name the model must guess.

    The stub's own header told the program ``from code_mode_api import
    delfin_tools`` while the skill saved the same text as ``tools_api.py``, so
    a program following the header could never import it.
    """
    stub = code_mode.generate_stub()
    assert "code_mode_api" not in stub
    assert "tools_api" in stub


def test_run_program_is_gone_and_the_hook_is_the_launcher():
    """One launcher. It is the bash tool's, which owns the cage and ceiling."""
    assert not hasattr(code_mode, "run_program")
    assert not hasattr(code_mode, "ProgramTimeout")
    assert not hasattr(code_mode, "ProgramError")
    # The surviving route is still the hook the bash tool calls.
    assert api_client._code_mode_child_env(
        "python3 -m delfin.agent.code_mode run prog.py", "/tmp",
        _permissions("/tmp"), _RecordingExecutor())


def test_serve_request_is_gone_and_the_server_is_the_only_request_path():
    """One request path, so a gate change cannot be made in one and not the other."""
    assert not hasattr(code_mode, "serve_request")

    ex = _RecordingExecutor()
    server = code_mode.RpcServer(ex, permissions=None)
    reply = json.loads(server.handle_payload(json.dumps(
        {"tool": "get_calc_info", "arguments": {"calc_id": "c1"}})))
    assert reply["result"]["ok"] == "get_calc_info"
    assert ex.calls == [("get_calc_info", {"calc_id": "c1"})]

    refused = json.loads(server.handle_payload(json.dumps(
        {"tool": "write_file", "arguments": {"path": "x", "content": "y"}})))
    assert "error" in refused
    assert len(ex.calls) == 1, "a refused call reached the executor"
