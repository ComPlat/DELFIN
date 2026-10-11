"""The panel shows the usage record of a skill (package L1, phase 2b).

collect_rows() now carries a per-proposal ``usage`` roll-up (uses /
failures / last_used / last_outcome) sourced from
delfin/agent/skill_usage, and format_panel_html() renders it — a skill
that keeps failing is shown with its record, not hidden. This is the
read side of "close the skill loop" in the review panel.

All writes go to a monkeypatched tmp HOME (proposals and the usage
ledger both live under there), so the test never touches real storage.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent import skill_proposals as sp  # noqa: E402
from delfin.agent import skill_usage as su  # noqa: E402
from delfin.dashboard.skill_proposals_panel import collect_rows  # noqa: E402
from delfin.dashboard.skill_proposals_panel import format_panel_html  # noqa: E402


def _propose(home, name="tune-xtb"):
    ev = sp.Evidence(kind="test", ref="tests/test_x.py::test_y")
    try:
        return sp.propose(
            name,
            "---\nname: tune-xtb\ndescription: tune xTB defaults\n---\n\n"
            "1. Probe the default.\n2. Fit the flag.\n",
            evidence=[ev], source="test")
    except Exception:
        # Real skill_safety may block the body on deny words; a fresh
        # name in a tmp HOME never collides and the body above is benign.
        raise


def test_the_panel_row_carries_the_usage_rollup(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    p = _propose(tmp_path)
    su.record_use(p.name, "success")
    su.record_use(p.name, "success")
    su.record_use(p.name, "failure")
    rows = {r["name"]: r for r in collect_rows()}
    assert p.name in rows
    usage = rows[p.name]["usage"]
    assert usage["uses"] == 3
    assert usage["failures"] == 1
    assert usage["last_outcome"] == "failure"


def test_the_panel_html_marks_a_failing_skill(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    p = _propose(tmp_path)
    su.record_use(p.name, "success")
    su.record_use(p.name, "failure")
    html = format_panel_html()
    assert p.name in html
    assert "1 fail" in html
    assert "used 2" in html


def test_a_never_used_proposal_renders_without_failure_mark(tmp_path,
                                                             monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    p = _propose(tmp_path, name="fresh-skill")
    html = format_panel_html()
    assert p.name in html
    assert "fail" not in html
    rows = {r["name"]: r for r in collect_rows()}
    assert rows[p.name]["usage"]["uses"] == 0
