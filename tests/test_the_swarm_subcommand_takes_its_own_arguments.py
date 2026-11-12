"""`delfin-agent swarm ...` hands its arguments to the swarm CLI.

Bug hunt 11.11.2026, #1 (blocker): the top-level ``swarm`` subparser
declared no arguments and its handler was ``cli_swarm.main``, which got
the parsed Namespace instead of an argv list. `delfin-agent swarm` died
with "TypeError: 'Namespace' object is not iterable", `swarm status` with
"unrecognized arguments: status", and `swarm --help` showed only ``-h``.

The cases through ``cli.build_parser`` need the cli.py hunk of the copy
pair (.gate/px/new_cli.py): cli.py hands its subparser to
``cli_swarm.configure_parser``.
"""

from __future__ import annotations

import json
import textwrap

import pytest

from delfin.agent import cli_swarm
from delfin.agent.swarm import status as swarm_status


@pytest.fixture()
def swarm_repo(tmp_path, monkeypatch):
    (tmp_path / "delfin" / "pkg").mkdir(parents=True)
    (tmp_path / "swarm.yaml").write_text(textwrap.dedent("""
        goal: Ship it
        packages:
          - name: core
            owned: ["delfin/pkg/**"]
            task: Build core
    """), encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(swarm_status, "_event_subscribers", [])
    return tmp_path


def _parse(argv):
    from delfin.agent.cli import build_parser
    return build_parser().parse_args(argv)


def test_swarm_status_json_runs_through_the_real_parser(swarm_repo, capsys):
    args = _parse(["swarm", "status", "--json"])
    assert args.func is cli_swarm.run_args
    assert args.swarm_cmd == "status" and args.as_json is True
    assert args.func(args) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["packages"][0]["name"] == "core"


def test_swarm_plan_takes_its_positional_goal(swarm_repo, capsys):
    args = _parse(["swarm", "plan", "Ship it"])
    assert args.func(args) == 0
    assert "core-builder" in capsys.readouterr().out


def test_bare_swarm_prints_its_help_and_exits_0(capsys):
    args = _parse(["swarm"])
    assert args.func(args) == 0
    out = capsys.readouterr().out
    for name in ("init", "plan", "start", "status"):
        assert name in out, out


def test_swarm_help_lists_the_subcommands(capsys):
    with pytest.raises(SystemExit) as exc:
        _parse(["swarm", "--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "delfin-agent swarm" in out
    for name in ("init", "plan", "start", "status", "triage", "merge"):
        assert name in out, out


def test_an_unknown_swarm_subcommand_is_a_usage_error(capsys):
    with pytest.raises(SystemExit) as exc:
        _parse(["swarm", "nonsense"])
    assert exc.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


def test_the_handler_stays_in_cli_swarm_for_the_login_node_guard():
    # hpc_guard refuses `swarm start` on a login node by the handler's
    # module; a lambda in cli.py would have slipped past it.
    assert _parse(["swarm", "start", "--dry-run"]).func.__module__ == \
        "delfin.agent.cli_swarm"


def test_main_with_an_argv_list_still_works(swarm_repo, capsys):
    assert cli_swarm.main(["status", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["packages"]
    assert cli_swarm.main([]) == 0
