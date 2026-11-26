"""Tests for delfin.agent.report_tools with fabricated trace stores."""

from __future__ import annotations

import json
import time

import pytest

from delfin.agent import report_tools


def _catalogue(*names: str) -> list[dict]:
    return [
        {"type": "function",
         "function": {"name": n, "description": n, "parameters": {}}}
        for n in names
    ]


def _write_trace(root, session: str, entries: list[dict]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    p = root / f"{session}.jsonl"
    lines = []
    for e in entries:
        # allow fabricating malformed lines via a raw string
        if isinstance(e, str):
            lines.append(e)
        else:
            lines.append(json.dumps(e))
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _entry(tool: str, *, ok=True, ts=None, **extra) -> dict:
    e = {"ts": ts if ts is not None else time.time(), "tool": tool,
         "ok": ok, "input": "", "output": ""}
    e.update(extra)
    return e


CAT = _catalogue("bash", "read_file", "search_docs", "never_touched")


def test_counts_calls_failures_and_share(tmp_path):
    _write_trace(tmp_path, "s1", [
        _entry("bash"), _entry("bash"), _entry("bash", ok=False),
        _entry("read_file"),
        _entry("search_docs", ok=False), _entry("search_docs", ok=False),
    ])
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert list(data["tools"]) == ["bash", "search_docs", "read_file"]
    assert data["tools"]["bash"] == {
        "calls": 3, "failures": 1, "failure_share": pytest.approx(1 / 3)}
    assert data["tools"]["search_docs"]["failure_share"] == pytest.approx(1.0)
    assert data["total_calls"] == 6
    assert data["sessions"] == 1


def test_never_used_comes_from_advertised_catalogue(tmp_path):
    _write_trace(tmp_path, "s1", [_entry("bash")])
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data["never_used"] == ["never_touched", "read_file", "search_docs"]
    assert data["n_advertised"] == 4
    assert data["n_never_used"] == 3
    assert data["n_used_advertised"] == 1


def test_default_catalogue_is_the_client_advertised_set(tmp_path):
    _write_trace(tmp_path, "s1", [_entry("bash")])
    data = report_tools.collect(root=tmp_path)
    # Not a typed-out list: derives from api_client's own surface.
    assert data["n_advertised"] > 10
    assert "bash" not in data["never_used"]
    assert "search_docs" in report_tools._advertised_names()
    assert isinstance(data["never_used"], list)


def test_mcp_prefix_stripped_for_comparison(tmp_path):
    _write_trace(tmp_path, "s1", [_entry("mcp__delfin-ops__list_tools")])
    data = report_tools.collect(root=tmp_path, catalogue=CAT + _catalogue("list_tools"))
    assert "list_tools" in data["tools"]
    assert "list_tools" not in data["never_used"]


def test_called_but_not_advertised_is_reported(tmp_path):
    _write_trace(tmp_path, "s1", [_entry("ghost_tool")])
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data["called_not_advertised"] == ["ghost_tool"]


def test_since_s_window_filters_old_entries(tmp_path):
    now = time.time()
    _write_trace(tmp_path, "s1", [
        _entry("bash", ts=now - 999_999),
        _entry("read_file", ts=now - 10),
    ])
    data = report_tools.collect(since_s=3600.0, root=tmp_path, catalogue=CAT)
    assert list(data["tools"]) == ["read_file"]
    assert data["window_s"] == 3600.0
    # Without a window both count.
    data_all = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data_all["total_calls"] == 2


def test_empty_and_malformed_store(tmp_path):
    _write_trace(tmp_path, "s1", ["{not json", json.dumps({"tool": "bash"}),
                                   json.dumps({"ok": True})])
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data["tools"]["bash"]["calls"] == 1  # missing ts counts (no window)
    assert data["total_calls"] == 1


def test_empty_store_yields_all_advertised_unused(tmp_path):
    data = report_tools.collect(root=tmp_path, catalogue=CAT)
    assert data["total_calls"] == 0
    assert data["n_never_used"] == 4
    text = report_tools.format_text(data)
    assert "No tool calls recorded." in text
    assert "- never_touched" in text


def test_format_text_orders_most_used_first_and_shows_failures():
    text = report_tools.format_text({
        "tools": {
            "bash": {"calls": 10, "failures": 0, "failure_share": 0.0},
            "grep_file": {"calls": 4, "failures": 2, "failure_share": 0.5},
        },
        "never_used": ["never_touched"], "n_advertised": 3,
        "n_never_used": 1, "total_calls": 14, "sessions": 2,
        "window_s": None,
    })
    assert text.index("bash") < text.index("grep_file")
    assert "50.0%" in text
    assert "0.0%" in text
    assert "never called: 1 of 3 (33.3%)" in text
    assert "- never_touched" in text


def test_main_prints_and_returns_zero(tmp_path, capsys):
    _write_trace(tmp_path, "s1", [_entry("bash")])
    rc = report_tools.main(["--trace-root", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Tool usage" in out
    assert "never called" in out


def test_the_role_surface_cost_is_printable_from_the_product(capsys):
    """``--roles`` reaches the per-role surface the client already computes.

    Red before this change for the reason it names: argparse knows no
    ``--roles``, so the entry point exits 2 on "unrecognized arguments" and
    no product command could answer what role scoping saves.

    No model is asked and no exact count is asserted: the catalogue grows,
    so the relation to the baseline is the claim, not a number.
    """
    rc = report_tools.main(["--roles"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "(unrestricted)" in out
    assert "office_agent" in out
    assert "saves" in out


def test_a_scoped_role_costs_less_than_the_unrestricted_surface():
    """The saving the report exists to show, measured rather than asserted.

    ``office_agent`` carries an execution allow-list
    (``_ROLE_EXEC_ALLOWLIST``), so its surface must be a strict subset of
    the unrestricted baseline in both tools and schema tokens. Failing
    input: a role whose allow-list grew past the baseline, or a report that
    stopped filtering by role at all.
    """
    report = report_tools.collect_roles()
    assert "" in report, "the unrestricted baseline row is missing"
    base = report[""]
    office = report["office_agent"]
    assert office["count"] < base["count"]
    assert office["total_tokens"] < base["total_tokens"]
    # The names must be the scoped ones, not the baseline's: a report that
    # ignored the role would show equal counts and this would still pass on
    # the token line alone, so pin the subset relation too.
    assert set(office["names"]) < set(base["names"])


def test_the_default_report_is_unchanged_by_the_roles_flag(tmp_path, capsys):
    """A control case: it passes before and after the change.

    Without ``--roles`` the usage report is the same one as before, so the
    new flag cannot have altered what every existing caller sees.
    """
    _write_trace(tmp_path, "s1", [_entry("bash")])
    rc = report_tools.main(["--trace-root", str(tmp_path)])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Tool usage across 1 traced session(s)" in out
    assert "(unrestricted)" not in out


def test_real_store_never_raises():
    # The default store: whatever is on this machine, collect must not raise.
    data = report_tools.collect()
    assert isinstance(data["tools"], dict)
    assert data["n_advertised"] >= 50


def test_a_preset_backed_role_is_not_missing_from_the_roles_table():
    """Every role the product actually restricts must have a row.

    ``_role_allowlist`` (``api_client.py:7624``) is what refuses a tool at
    the execution call site, and it answers for the two roles in
    ``_ROLE_EXEC_ALLOWLIST`` AND for every subagent preset that carries a
    ``tools:`` line. ``role_tool_surface_report`` enumerates only the two
    dicts, so a preset-backed role was restricted with no row in the table
    this command prints -- a reader asking "what does role scoping save"
    saw 8 of the 19 restricted roles and concluded the other 11 were
    unrestricted.

    Failing input: any preset whose ``tools:`` line restricts it. Removing
    the preset enumeration from ``collect_roles`` makes this fail with
    ``explore`` absent from the report.
    """
    from delfin.agent.api_client import _role_allowlist
    from delfin.agent.subagents import SUBAGENT_PRESETS

    restricted = sorted(
        role for role in SUBAGENT_PRESETS if _role_allowlist(role) is not None)
    assert restricted, "no preset role is restricted, so this test measures nothing"

    report = report_tools.collect_roles()
    missing = [role for role in restricted if role not in report]
    assert not missing, (
        "these roles are restricted by _role_allowlist but have no row in "
        f"the --roles table: {missing}")

    # A row is not enough on its own: it must be the SCOPED surface. A report
    # that listed the role but ignored its allow-list passes the membership
    # check above and is still wrong, so pin the subset relation for one
    # preset-backed role.
    base = report[""]
    for role in restricted:
        assert report[role]["count"] <= base["count"], role
    scoped = [r for r in restricted if report[r]["count"] < base["count"]]
    assert scoped, (
        "no preset-backed role has a surface narrower than the baseline, so "
        "the rows added by this report are not scoped at all")


def test_the_roles_table_names_the_cheapest_scoped_role():
    """The printed table answers "which role offers least" for ALL roles.

    A control that passes before and after: the ordering question is what
    the table is for, and the answer must be a restricted role, not the
    baseline row which sorting could otherwise place first.
    """
    report = report_tools.collect_roles()
    text = report_tools.format_roles(report)
    rows = [ln for ln in text.splitlines() if ln and not ln.startswith(("-", "Tool", "role"))]
    assert rows[0].startswith("(unrestricted)"), rows[0]
    # The cheapest scoped row is the one with the fewest tokens after the
    # baseline; it must be narrower than the baseline, never equal to it.
    scoped = [(v["total_tokens"], r) for r, v in report.items() if r]
    cheapest_tokens, cheapest_role = min(scoped)
    assert cheapest_tokens < report[""]["total_tokens"]
    assert cheapest_role in text
