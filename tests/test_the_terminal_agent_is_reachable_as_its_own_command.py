"""`delfin-agent` must be a real command, and a bare one must reach the agent.

Two things are pinned here. First, the console script is DECLARED — a
front door that no packaging metadata mentions is not a command, however
well the function behind it works. Second, the implicit-`chat` routing is
derived from the parser instead of a hand-kept list: a duplicated list is
one edit away from routing a real subcommand into the chat prompt, which
would turn a typo'd command into a billed turn.
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli


_REPO = Path(__file__).resolve().parents[1]


def _scripts() -> dict[str, str]:
    if sys.version_info >= (3, 11):
        import tomllib
    else:  # pragma: no cover - the floor is 3.10, so this path is real
        pytest.importorskip("tomli")
        import tomli as tomllib  # type: ignore[no-redef]
    with open(_REPO / "pyproject.toml", "rb") as fh:
        return (tomllib.load(fh).get("project", {}) or {}).get("scripts", {}) or {}


def test_the_console_script_is_declared_and_its_target_is_callable():
    scripts = _scripts()
    assert "delfin-agent" in scripts, (
        "pyproject.toml declares no delfin-agent script — nothing installs "
        "the command, whatever cli.py can do"
    )
    target = scripts["delfin-agent"]
    module_name, _, attr = target.partition(":")
    module = importlib.import_module(module_name)
    assert callable(getattr(module, attr, None)), (
        f"{target} does not resolve to a callable"
    )


def test_the_command_names_itself_in_its_own_help():
    # argparse prints `prog` in every usage and error line, so this is the
    # name users are told to type when they get something wrong.
    assert agent_cli.build_parser().prog == "delfin-agent"


def test_a_bare_invocation_routes_to_the_session():
    known = agent_cli._subcommand_names(agent_cli.build_parser())
    assert agent_cli._route_argv([], known) == ["chat"]
    assert agent_cli._route_argv(["--session", "x"], known) == [
        "chat", "--session", "x"]
    assert agent_cli._route_argv(["-p", "hi"], known) == ["chat", "-p", "hi"]


def test_a_real_subcommand_is_left_exactly_as_typed():
    parser = agent_cli.build_parser()
    known = agent_cli._subcommand_names(parser)
    # Every subcommand the parser registers, not a sample of them: adding a
    # subcommand without teaching the router about it is the failure this
    # guards, so the test has to enumerate from the same source the router
    # does and then prove the router agrees.
    assert known, "the parser registers no subcommands at all"
    for name in sorted(known):
        assert agent_cli._route_argv([name, "--x"], known) == [name, "--x"], (
            f"subcommand {name!r} was swallowed by the chat router"
        )
    for flag in ("-h", "--help"):
        assert agent_cli._route_argv([flag], known) == [flag]


def _registered_subcommands(parser) -> frozenset[str]:
    """Walk the parser here, independently of the code under test.

    Deriving the expectation from `_subcommand_names` would make the check
    circular: a stale hardcoded list would satisfy a test that asked the
    same stale list what to expect.
    """
    import argparse
    found: set[str] = set()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            found.update(action.choices)
    return frozenset(found)


def test_the_router_reads_the_parsers_own_subcommands():
    # The point of deriving rather than duplicating. Replace
    # _subcommand_names with a literal set and this holds only until the
    # next subcommand is added — which is precisely the drift being
    # guarded, and it fails the moment it happens.
    parser = agent_cli.build_parser()
    assert agent_cli._subcommand_names(parser) == _registered_subcommands(parser)
    for expected in ("chat", "run", "init", "session", "doctor"):
        assert expected in _registered_subcommands(parser)


def test_chat_is_the_default_and_carries_the_dashboard_defaults():
    args = agent_cli.build_parser().parse_args(["chat"])
    assert args.cmd == "chat"
    assert args.func is agent_cli.cmd_chat
    # Only this front door inherits the provider/model the user configured;
    # `run` must keep its historical defaults so the scheduler and the
    # benchmark are not moved by a settings file.
    assert args.settings_defaults is True
    run_args = agent_cli.build_parser().parse_args(["run", "x"])
    assert getattr(run_args, "settings_defaults", False) is False
