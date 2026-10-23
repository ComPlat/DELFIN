"""Looking things up must stay in the reading tools.

Measured on 2026-09-25: six supervised sessions cost the operator 17
decisions, almost all on harmless READ commands that only had to ask
because of their FORM — ``$( … )`` command substitution
(``sed -n "$(grep -n … | cut -d: -f1),+40p" file``), ``awk``, ``xargs``,
backticks inside double quotes, ``sed -i`` instead of ``edit_file``,
chains with ``2>/dev/null``. The gate is RIGHT to ask (it cannot read
substitutions); the fix is the agent's behaviour, not a looser gate.

Out of the 378 one-line bash commands of those six sessions
(``.gate/cmds_heute.txt``, not committed), 14 used ``$( )`` or backticks,
11 ran ``awk``, 1 ``xargs``, 1 ``sed -i``, 20 chained ``2>/dev/null`` —
and the same replay through the gate of the then-current main asked on
10 of them (``.gate/replay_heute.txt``). Every single one of those asks
had a freely-running equivalent in the reading tools: ``grep_file`` for
the grep, ``read_file`` with offset/limit for the ``sed -n``, two
separate commands instead of the substitution, ``edit_file`` instead of
``sed -i``.

Every check below reads the ONE dedicated section of the prompt, so
before the rule existed every one of them failed for the same reason:
the section was not there. Control on the previous commit: 6 cases,
all 6 red.

This test pins the rule in the prompt so it survives prompt rewrites
the way the other pinned prompt rules do.
"""

from __future__ import annotations

import re
from pathlib import Path

_SOLO = (
    Path(__file__).resolve().parents[1]
    / "delfin" / "agent" / "pack" / "agents" / "solo_agent.md"
)

# The heading the rule lives under. Keep in sync with solo_agent.md.
_HEADING = re.compile(r"^## Look things up through the reading tools$", re.M)


def _section() -> str:
    """The whole lookup-rule section: heading to the next ``##`` heading.
    Empty string when the section does not exist (yet) — every check
    then fails on the same missing artifact, not on wording."""
    text = _SOLO.read_text(encoding="utf-8")
    m = _HEADING.search(text)
    if m is None:
        return ""
    tail = text[m.start():]
    nxt = re.search(r"^## ", tail[m.end() - m.start():], re.M)
    return tail if nxt is None else tail[: nxt.start() + (m.end() - m.start())]


def test_command_substitution_is_forbidden_for_lookups():
    """`$( grep -n … )` is the most common ask: the gate cannot read what
    the substitution will produce, so it must ask. The rule has to name
    the form outright."""
    assert re.search(r"\$\(\s*\)|command substitution", _section(), re.I)


def test_backticks_awk_and_xargs_are_named():
    """The replay asked on `awk`; backticks inside double quotes were a
    live failure form the same day. A rule that does not name them lets
    the next session rediscover each one as a dialog."""
    s = _section().lower()
    assert "backtick" in s
    assert re.search(r"\bawk\b", s)
    assert re.search(r"\bxargs\b", s)


def test_two_steps_instead_of_substitution_is_prescribed():
    """grep for the line number, then read from there — two freely-running
    commands instead of one that asks. The section must say the two-step
    shape, not only 'avoid substitution'."""
    s = _section().lower()
    assert "two" in s and "step" in s


def test_edit_file_is_the_in_place_edit_tool():
    """`sed -i` is an edit wearing a command's clothes: the file-confirm
    layer never sees a write tool. The section must name `edit_file` as
    the tool to reach for, or the ban reads as gate-lore rather than as
    the way edits are made."""
    s = _section().lower()
    assert "sed -i" in s
    assert "edit_file" in s


def test_reading_tools_are_the_named_alternative():
    """`grep_file` and `read_file` (with offset/limit) must appear in the
    section — the asks drop to zero only when the freely-running
    equivalent is spelled out, not when the shell form is merely banned."""
    s = _section().lower()
    assert "grep_file" in s
    assert "read_file" in s


def test_the_rule_is_short():
    """The rule must stay a rule, not a manual: at most 12 non-empty lines
    of prompt text may be spent on it."""
    section = _section()
    assert section, "the lookup-rule section is gone from solo_agent.md"
    body = [l for l in section.splitlines()[1:] if l.strip()]
    assert len(body) <= 12, f"lookup rule grew to {len(body)} non-empty lines"
