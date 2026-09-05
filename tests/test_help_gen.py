"""Generated /help (delfin.agent.help_gen) on synthetic palette rows.

The renderer must show the rows verbatim (no command knowledge of its
own), group and align them, filter by category/search, and never raise.
``coverage_gaps`` must detect drift in both directions between the
palette table and the registered slash prefixes.
"""

from __future__ import annotations

import re

from delfin.agent.help_gen import coverage_gaps, generate_help


ROWS = (
    ("Session", "/help", "Show this help", False),
    ("Session", "/search", "Search in chat history", True),
    ("Git", "/git status", "Show git status", False),
    ("Git", "/git diff", "Show staged/unstaged changes", False),
    ("Memory", "/memories", "List project memories", False),
)


# ---------------------------------------------------------------------------
# rendering
# ---------------------------------------------------------------------------

def test_groups_in_first_appearance_order_with_counts():
    out = generate_help(ROWS)
    assert out.startswith("Available commands (5):")
    assert (out.index("Session:") < out.index("Git:")
            < out.index("Memory:"))
    # Every command + summary rendered verbatim.
    for _, command, summary, _ in ROWS:
        assert command in out
        assert summary in out


def test_alignment_and_has_args_marker():
    out = generate_help(ROWS)
    cmd_lines = [ln for ln in out.splitlines() if ln.startswith("  /")]
    assert len(cmd_lines) == 5
    # One aligned separator column across ALL groups.
    positions = {ln.index("—") for ln in cmd_lines}
    assert len(positions) == 1
    # has_args row carries the ellipsis marker, no-args rows do not.
    assert any(ln.startswith("  /search ...") for ln in cmd_lines)
    assert not any("/help ..." in ln for ln in cmd_lines)


def test_renders_rows_verbatim_not_hardcoded():
    rows = (("Zeta", "/frobnicate", "Adjust the flux window", True),)
    out = generate_help(rows)
    assert "Zeta:" in out
    assert "/frobnicate ..." in out
    assert "Adjust the flux window" in out


# ---------------------------------------------------------------------------
# filtering
# ---------------------------------------------------------------------------

def test_category_filter_case_insensitive():
    out = generate_help(ROWS, category="git")
    assert "Commands category 'git' (2):" in out
    assert "/git status" in out and "/git diff" in out
    assert "Session:" not in out and "/help" not in out


def test_search_matches_command_summary_and_category():
    assert "/git status" in generate_help(ROWS, search="status")
    out = generate_help(ROWS, search="chat history")   # summary match
    assert "/search" in out and "/git" not in out
    out = generate_help(ROWS, search="memory")         # category match
    assert "/memories" in out and "/help" not in out


def test_combined_filters_and_no_match_message():
    out = generate_help(ROWS, category="Git", search="diff")
    assert "(1)" in out and "/git diff" in out
    miss = generate_help(ROWS, search="no-such-thing")
    assert "No commands" in miss
    assert "Try /help with no filter." in miss


# ---------------------------------------------------------------------------
# never raises
# ---------------------------------------------------------------------------

def test_generate_help_never_raises_on_garbage():
    for bad in (None, 123, "nope", (("only-two", "/x"),),
                ((None,),), (("Cat", "no-slash", "s", False),),
                (object(),)):
        out = generate_help(bad)
        assert isinstance(out, str)
    # Well-formed rows mixed with garbage: good rows still render.
    mixed = (("Cat", "/ok", "fine", False), ("short",), 42)
    out = generate_help(mixed)
    assert "/ok" in out and "fine" in out


def test_coverage_gaps_never_raises_on_garbage():
    for rows, prefixes in ((None, None), (123, 456), ("x", "y")):
        res = coverage_gaps(rows, prefixes)
        assert isinstance(res, dict)
        assert "unregistered_commands" in res
        assert "unlisted_prefixes" in res


# ---------------------------------------------------------------------------
# coverage_gaps drift detection
# ---------------------------------------------------------------------------

def test_coverage_gaps_clean():
    prefixes = {"/help", "/search", "/git", "/memories"}
    res = coverage_gaps(ROWS, prefixes)
    assert res == {"ok": True, "unregistered_commands": [],
                   "unlisted_prefixes": []}


def test_coverage_gaps_detects_both_directions():
    # "/search" prefix dropped from the dispatcher; "/ghost" registered
    # but never documented in the palette.
    prefixes = {"/help", "/git", "/memories", "/ghost"}
    res = coverage_gaps(ROWS, prefixes)
    assert res["ok"] is False
    assert res["unregistered_commands"] == ["/search"]
    assert res["unlisted_prefixes"] == ["/ghost"]


def test_coverage_gaps_uses_first_token_of_multiword_commands():
    # "/git status" and "/git diff" are both covered by the "/git"
    # prefix — multiword rows must not be flagged.
    res = coverage_gaps(ROWS, {"/git"})
    assert "/git status" not in res["unregistered_commands"]
    assert "/git diff" not in res["unregistered_commands"]
    assert res["unlisted_prefixes"] == []


def test_help_and_gaps_agree_on_the_same_table():
    """The pair forms a closed loop: whatever /help renders is exactly
    what the drift check considers covered."""
    prefixes = {"/help", "/search", "/git", "/memories"}
    out = generate_help(ROWS)
    rendered_cmds = {
        m.group(1)
        for m in re.finditer(r"^  (/\S+)", out, flags=re.MULTILINE)}
    assert {c.split()[0] for c in rendered_cmds} <= prefixes
    assert coverage_gaps(ROWS, prefixes)["ok"] is True
