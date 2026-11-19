"""`delfin-agent <typo>` is not a billed prompt (hunt finding 22).

`delfin-agent stauts` started a chat titled "stauts". As git, gh and kubectl
do: a single bare word close to a subcommand is answered with "did you mean"
and exit 2; real one-word prompts still go through.
"""
import pytest

from delfin.agent.cli import _route_argv

KNOWN = frozenset({"chat", "run", "status", "sessions", "service", "swarm", "tui"})


def test_a_typo_names_the_command_and_exits_2(capsys):
    with pytest.raises(SystemExit) as exc:
        _route_argv(["stauts"], KNOWN)
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "Did you mean 'status'" in err and "delfin-agent chat stauts" in err


@pytest.mark.parametrize("argv", [["hello"], ["explain", "this"], ["fix"], ["status"]])
def test_prompts_and_real_commands_pass(argv):
    out = _route_argv(argv, KNOWN)
    assert out[0] in ("chat", "status")
