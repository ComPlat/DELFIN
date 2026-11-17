"""CM1: the runner command is recognised and the RPC path is scoped right.

The protected hook inside _dispatch_under_rules calls code_mode.
runner_program_path to decide whether a bash command is a Code Mode
program (then DELFIN_CODE_MODE_RPC gets the per-invocation socket path) or
a plain command (then {} — nothing added). The detection is what keeps the
env var scoped to the caged child only, so it gets its own red tests here.
"""
import os

from delfin.agent import code_mode


def test_runner_path_recognised_in_a_plain_bash_command():
    cmd = "python -m delfin.agent.code_mode run prog.py"
    assert code_mode.runner_program_path(cmd) == "prog.py"


def test_runner_path_with_sys_executable_and_quoting():
    cmd = "/usr/bin/python3 -m delfin.agent.code_mode run my_prog.py"
    assert code_mode.runner_program_path(cmd) == "my_prog.py"


def test_plain_command_yields_none():
    assert code_mode.runner_program_path("ls -la delfin/") is None
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode") is None
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode --help") is None


def test_injection_characters_terminate_the_path():
    # The path token ends before any shell metacharacter: a crafted
    # command cannot smuggle a second argument into the runner.
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run a.py; rm x") == "a.py"
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run a.py | cat") == "a.py"
    assert code_mode.runner_program_path(
        "python -m delfin.agent.code_mode run a.py && b") == "a.py"


def test_socket_path_is_scoped_to_the_workspace(tmp_path):
    path = code_mode.rpc_socket_path(str(tmp_path), "tok1")
    assert path.startswith(os.path.join(str(tmp_path), ".code_mode_rpc"))
    assert path.endswith("rpc-tok1.sock")
