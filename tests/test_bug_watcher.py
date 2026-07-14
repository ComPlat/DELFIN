"""Tests for the bug-report watcher (maintainer triage tool).

All offline: the analysis turn is a mocked engine (no LLM, no network), the
archive is a tmp dir. Verifies the security-relevant contract — triage writes
an analysis + proposed patch but NEVER moves a report to Solved (a triaged
report is not a fixed one); moving is the explicit ``mark_solved`` step.
"""

from __future__ import annotations

import json

from delfin.agent import bug_watcher as bw


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_FAKE_ANALYSIS = (
    "ROOT CAUSE: the parser mishandles empty input in foo.py.\n"
    "IMPACT: the agent crashes.\n"
    "SEVERITY: high\n"
    "FIX:\n"
    "```diff\n"
    "--- a/delfin/agent/foo.py\n"
    "+++ b/delfin/agent/foo.py\n"
    "@@ -1,2 +1,3 @@\n"
    " def parse(x):\n"
    "+    if not x: return None\n"
    "     return x.split()\n"
    "```\n"
)


class _FakeEngine:
    def __init__(self, response: str):
        self._response = response

    def stream_response(self, user_message="", on_token=None, **kw):
        if on_token:
            on_token(self._response)
        return self._response


def _factory(response: str = _FAKE_ANALYSIS):
    def make(repo_root, settings=None):
        return _FakeEngine(response)
    return make


def _make_report(archive, name="20260101-000000_u_sess_ab12", *,
                 description="agent crashed on empty input"):
    d = archive / name
    d.mkdir(parents=True)
    (d / "report.json").write_text(json.dumps({
        "schema": "delfin-bug-report/1",
        "description": description,
        "error_text": "IndexError in foo.parse",
        "model": "kit.qwen3.5-397b-A17b",
        "mode": "solo",
        "tool_trace": [{"tool": "read_file", "ok": True, "error": ""}],
    }), encoding="utf-8")
    (d / "report.md").write_text(
        f"# Bug report\n\n{description}\n\nError: IndexError in foo.parse\n",
        encoding="utf-8",
    )
    return d


# ---------------------------------------------------------------------------
# find_unsolved
# ---------------------------------------------------------------------------

def test_find_unsolved_lists_untriaged_reports(tmp_path):
    _make_report(tmp_path, "r1")
    _make_report(tmp_path, "r2")
    found = {d.name for d in bw.find_unsolved(tmp_path, settings={})}
    assert found == {"r1", "r2"}


def test_find_unsolved_skips_solved_and_non_reports(tmp_path):
    _make_report(tmp_path, "r1")
    (tmp_path / "Solved").mkdir()
    _make_report(tmp_path / "Solved", "already_done")  # inside Solved/
    (tmp_path / "random_dir").mkdir()                  # no report.json
    found = {d.name for d in bw.find_unsolved(tmp_path, settings={})}
    assert found == {"r1"}


def test_find_unsolved_skips_already_triaged(tmp_path):
    d = _make_report(tmp_path, "r1")
    (d / "triage.md").write_text("already analysed", encoding="utf-8")
    assert bw.find_unsolved(tmp_path, settings={}) == []


def test_find_unsolved_skips_poison_after_max_attempts(tmp_path):
    d = _make_report(tmp_path, "r1")
    (d / bw._ATTEMPTS_FILE).write_text("3", encoding="utf-8")
    assert bw.find_unsolved(tmp_path, settings={}) == []


# ---------------------------------------------------------------------------
# triage_report — writes analysis + patch, NEVER moves
# ---------------------------------------------------------------------------

def test_triage_writes_analysis_and_patch_without_moving(tmp_path):
    d = _make_report(tmp_path, "r1")
    t = bw.triage_report(d, settings={}, engine_factory=_factory())
    # analysis + patch persisted in place
    assert (d / "triage.md").exists()
    assert (d / "fix_proposal.patch").exists()
    assert "if not x: return None" in (d / "fix_proposal.patch").read_text()
    assert t.has_fix is True
    assert not t.error
    # CRITICAL: the report is NOT moved to Solved
    assert d.exists() and d.parent == tmp_path
    assert not (tmp_path / "Solved").exists()


def test_triage_no_fix_block_means_no_patch(tmp_path):
    d = _make_report(tmp_path, "r1")
    t = bw.triage_report(d, settings={},
                         engine_factory=_factory("ROOT CAUSE: x. no diff here."))
    assert (d / "triage.md").exists()
    assert not (d / "fix_proposal.patch").exists()
    assert t.has_fix is False


def test_triage_propose_fix_off_still_analyses(tmp_path):
    d = _make_report(tmp_path, "r1")
    settings = {"agent": {"bug_watcher": {"propose_fix": False}}}
    t = bw.triage_report(d, settings=settings, engine_factory=_factory())
    assert (d / "triage.md").exists()
    # propose_fix off → the diff in the response is ignored, no patch file
    assert not (d / "fix_proposal.patch").exists()
    assert t.has_fix is False


def test_triage_auto_analyze_off_is_zero_token(tmp_path):
    d = _make_report(tmp_path, "r1")

    def _boom(repo_root, settings=None):
        raise AssertionError("engine must not be built when auto_analyze off")

    settings = {"agent": {"bug_watcher": {"auto_analyze": False}}}
    t = bw.triage_report(d, settings=settings, engine_factory=_boom)
    assert (d / "triage.md").exists()
    assert "auto_analyze off" in t.summary
    assert d.exists() and not (tmp_path / "Solved").exists()


def test_triage_analysis_failure_bumps_attempts_and_leaves_report(tmp_path):
    d = _make_report(tmp_path, "r1")

    def _raise(repo_root, settings=None):
        class E:
            def stream_response(self, **kw):
                raise RuntimeError("api down")
        return E()

    t = bw.triage_report(d, settings={}, engine_factory=_raise)
    assert t.error and "api down" in t.error
    assert not (d / "triage.md").exists()          # not triaged
    assert bw._attempts(d) == 1                     # attempt recorded
    assert d.exists()                               # left in place for retry


# ---------------------------------------------------------------------------
# mark_solved — the explicit, human move
# ---------------------------------------------------------------------------

def test_mark_solved_moves_into_solved(tmp_path):
    d = _make_report(tmp_path, "r1")
    dest = bw.mark_solved(d, archive=tmp_path)
    assert dest
    assert not d.exists()
    assert (tmp_path / "Solved" / "r1" / "report.json").exists()


def test_run_once_triages_all_without_moving(tmp_path):
    _make_report(tmp_path, "r1")
    _make_report(tmp_path, "r2")
    results = bw.run_once(archive=tmp_path, settings={},
                          engine_factory=_factory())
    assert len(results) == 2
    assert all(r.has_fix for r in results)
    # nothing moved; both now carry triage.md so a re-scan finds none
    assert not (tmp_path / "Solved").exists()
    assert bw.find_unsolved(tmp_path, settings={}) == []
