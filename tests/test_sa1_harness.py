"""SA1 phase 2: the measurement harness and the correctness checker.

The pilot runs the SAME small set of read-heavy tasks twice — arm A plain,
arm B allowed to call the explore sub-agent — and records per arm: input
tokens, output tokens, wall time and tool calls. Correctness: every
``path:line`` claim in the final text must exist and the line must say
what the claim says (``verify_citations``, driven by the real filesystem
via tmp_path, no mocks for the part that checks files).
"""
import textwrap

import pytest

from delfin.agent import explore_pilot


# ---------------------------------------------------------------------------
# verify_citations — the correctness oracle
# ---------------------------------------------------------------------------

class TestVerifyCitations:
    def test_every_correct_claim_passes(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text("import os\n\nX = 1\n", encoding="utf-8")
        report = explore_pilot.verify_citations(
            f"the constant sits in mod.py:3 — see mod.py:3", root=tmp_path)
        assert report.claims == 1 and report.bad == []

    def test_missing_file_is_flagged(self, tmp_path):
        report = explore_pilot.verify_citations("see nope.py:1", root=tmp_path)
        assert report.claims == 1
        assert report.bad and "nope.py" in report.bad[0]

    def test_wrong_line_content_is_flagged(self, tmp_path):
        f = tmp_path / "mod.py"
        f.write_text("import os\n\nX = 1\n", encoding="utf-8")
        report = explore_pilot.verify_citations(
            "the X import sits in mod.py:3", root=tmp_path, expect="import")
        assert report.claims == 1
        assert report.bad

    def test_no_claims_is_zero_claims_not_an_error(self):
        report = explore_pilot.verify_citations("no citations at all")
        assert report.claims == 0 and report.bad == []


# ---------------------------------------------------------------------------
# The harness: same tasks, two arms, recorded metrics
# ---------------------------------------------------------------------------

def _fake_task(task_id="where-is-x"):
    return explore_pilot.PilotTask(
        task_id=task_id,
        prompt="find where DELFIN decides X; answer with file:line")


def _fake_run_task(task, *, arm):
    """A run_task double: arm A does 4 tool calls, arm B 1 (plus the
    delegate), arm B uses more tokens for less wall time — the numbers
    the harness must report verbatim."""
    if arm == "plain":
        return dict(task_id=task.task_id, success=True, input_tokens=100,
                    output_tokens=50, tool_calls=4, duration_s=8.0,
                    final_text="it is in mod.py:1")
    return dict(task_id=task.task_id, success=True, input_tokens=140,
                output_tokens=30, tool_calls=2, duration_s=4.0,
                final_text="it is in mod.py:1")


class TestRunPilot:
    def test_runs_both_arms_for_every_task(self):
        tasks = [_fake_task(), _fake_task("where-is-y")]
        result = explore_pilot.run_pilot(
            tasks, runner=_fake_run_task, root=".")
        assert sorted(r.task_id for r in result) == ["where-is-x",
                                                     "where-is-y"]
        arms = {r.task_id: r for r in result}
        assert set(arms["where-is-x"].arms) == {"plain", "explore"}

    def test_metrics_are_recorded_per_arm(self):
        (row,) = explore_pilot.run_pilot([_fake_task()], runner=_fake_run_task,
                                         root=".")
        plain, pilot = row.arms["plain"], row.arms["explore"]
        assert (plain.input_tokens, plain.output_tokens) == (100, 50)
        assert (pilot.input_tokens, pilot.output_tokens) == (140, 30)
        assert plain.tool_calls == 4 and pilot.tool_calls == 2
        assert plain.duration_s == 8.0 and pilot.duration_s == 4.0

    def test_correctness_is_checked_against_the_real_files(self, tmp_path):
        (tmp_path / "mod.py").write_text("DECISION = 'X'\n", encoding="utf-8")
        def runner(task, *, arm):
            d = _fake_run_task(task, arm=arm)
            d["final_text"] = "it is in mod.py:1"
            return d
        (row,) = explore_pilot.run_pilot([_fake_task()], runner=runner,
                                         root=tmp_path)
        assert row.arms["plain"].citations_ok
        assert row.arms["explore"].citations_ok

    def test_a_wrong_citation_marks_the_arm_not_correct(self, tmp_path):
        (tmp_path / "mod.py").write_text("import os\n", encoding="utf-8")
        def runner(task, *, arm):
            d = _fake_run_task(task, arm=arm)
            d["final_text"] = "it is in mod.py:2"
            return d
        (row,) = explore_pilot.run_pilot([_fake_task()], runner=runner,
                                         root=tmp_path)
        assert not row.arms["plain"].citations_ok

    def test_verdict_compares_arms(self):
        verdict = explore_pilot.pilot_verdict(
            explore_pilot.run_pilot([_fake_task()], runner=_fake_run_task,
                                    root="."))
        # The fake: arm B fewer tool calls + less wall time, equal success.
        assert verdict.explore_wins is True
        assert verdict.equal_correctness is True
