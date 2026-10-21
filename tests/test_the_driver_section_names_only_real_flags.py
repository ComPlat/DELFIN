"""The manual's driver section names only flags the parser has.

Input: the section "Driving DELFIN from another program or agent" in
docs/USER_MANUAL.md, and ``delfin.agent.cli.build_parser()``.
Output: one assertion per kind of drift.

Semantics. Every ``--flag`` (and ``-p``) spelled anywhere in the section
-- prose, tables, the shell lines and the Python example's argument list
-- must be an option of ``delfin-agent chat`` or of a ``delfin-agent
approvals`` subcommand, the two parsers the section documents. The values
the section lists for ``--permission-mode`` and ``--output-format`` must be
among that option's ``choices``. A flag renamed or removed in the parser,
or documented before it exists, fails here instead of in a driver.

Rejected: checking only the fenced shell lines. The recipe test
(test_the_manual_every_recipe_runs.py) already parses those through the
real argparse, but most of the section's flags are in tables and prose,
which nothing read.
"""

from __future__ import annotations

import argparse
import re

from delfin.agent import cli as agent_cli
from tests.manual_blocks import repo_root

_HEADING = "### Driving DELFIN from another program or agent"
_FLAG_RE = re.compile(r"(?<![\w-])(--[a-z][a-z0-9-]*|-p)(?![\w-])")


def _section() -> str:
    text = (repo_root() / "docs" / "USER_MANUAL.md").read_text(encoding="utf-8")
    assert _HEADING in text, f"{_HEADING!r} is missing from the manual"
    body = text.split(_HEADING, 1)[1]
    # Ends at the next heading of the same or a higher level, or the rule
    # that closes the chapter.
    end = re.search(r"^(?:#{1,3} |---\s*$)", body, re.MULTILINE)
    return body[: end.start()] if end else body


def _subparsers(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return dict(action.choices)
    return {}


def _options(parser: argparse.ArgumentParser) -> dict[str, argparse.Action]:
    out: dict[str, argparse.Action] = {}
    for action in parser._actions:
        for opt in action.option_strings:
            out[opt] = action
    return out


def _documented_parsers() -> dict[str, argparse.Action]:
    top = _subparsers(agent_cli.build_parser())
    known = _options(top["chat"])
    for sub in _subparsers(top["approvals"]).values():
        known.update(_options(sub))
    return known


def test_the_section_exists_and_names_flags():
    flags = set(_FLAG_RE.findall(_section()))
    # A section that lost its flags would pass the next test vacuously.
    assert {"-p", "--output-format", "--permission-mode"} <= flags


def test_every_flag_the_section_names_is_in_the_parser():
    known = _documented_parsers()
    missing = sorted(set(_FLAG_RE.findall(_section())) - set(known))
    assert missing == [], (
        f"the driver section names flags no parser has: {missing}")


def test_the_values_the_section_lists_are_the_choices():
    known = _documented_parsers()
    section = _section()
    for flag, values in (
            ("--permission-mode",
             ("plan", "default", "acceptEdits", "bypassPermissions")),
            ("--output-format", ("text", "json", "stream-json"))):
        choices = set(known[flag].choices or ())
        for value in values:
            assert f"`{value}`" in section, (
                f"{flag} {value} is no longer named in the section")
            assert value in choices, (
                f"the section lists {flag} {value}; the parser offers "
                f"{sorted(choices)}")
