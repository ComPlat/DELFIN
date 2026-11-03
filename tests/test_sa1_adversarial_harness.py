"""SA1 phase-2 adversarial probes (reviewer nacht-s17): harness edges.

Attacks verify_citations, run_pilot and pilot_verdict on the shapes the
builder's own suite does not cover: pathological citation strings, a
runner that returns None or partial dicts, and mixed evidence that must
yield verdict None (never a claim).
"""
from delfin.agent import explore_pilot


class TestCitationRegexEdges:
    def test_citation_in_a_longer_token_is_still_found(self, tmp_path):
        # "delfin/agent/foo.py:12" inside a sentence with punctuation.
        (tmp_path / "foo.py").write_text("alpha\nbeta\n", encoding="utf-8")
        report = explore_pilot.verify_citations(
            "see (foo.py:2) for beta", root=tmp_path)
        assert report.claims == 1 and report.bad == []

    def test_bare_line_number_without_path_is_not_a_claim(self, tmp_path):
        # "line 12" alone is not a path:line citation and must not count.
        report = explore_pilot.verify_citations("it sits at line 12",
                                                root=tmp_path)
        assert report.claims == 0

    def test_line_zero_is_out_of_range(self, tmp_path):
        (tmp_path / "m.py").write_text("one\n", encoding="utf-8")
        report = explore_pilot.verify_citations("m.py:0", root=tmp_path)
        assert report.claims == 1 and report.bad

    def test_unicode_path_citation(self, tmp_path):
        f = tmp_path / "große_datei.py"
        f.write_text("x = 'ü'\n", encoding="utf-8")
        report = explore_pilot.verify_citations(
            "große_datei.py:1", root=tmp_path)
        # Unicode in the path: the regex class does not include non-ASCII,
        # so this is either not counted as a claim, or counted and checked
        # against the real file. It must never be counted as a GOOD claim
        # while pointing at nothing.
        if report.claims:
            assert not report.bad

    def test_empty_text_is_zero_claims_and_ok(self):
        report = explore_pilot.verify_citations("")
        assert report.claims == 0 and report.bad == []

    def test_none_text_is_handled(self):
        report = explore_pilot.verify_citations(None)
        assert report.claims == 0 and report.bad == []


class TestRunPilotEdgeShapes:
    def _task(self):
        return explore_pilot.PilotTask(task_id="t1", prompt="p")

    def test_runner_returning_none_does_not_crash(self):
        rows = explore_pilot.run_pilot([self._task()], runner=lambda t, arm: None,
                                       root=".")
        arm = rows[0].arms["plain"]
        assert arm.success is False and arm.tool_calls == 0

    def test_runner_missing_keys_default_to_zero(self):
        rows = explore_pilot.run_pilot(
            [self._task()],
            runner=lambda t, arm: {"success": True, "final_text": ""},
            root=".")
        arm = rows[0].arms["explore"]
        assert (arm.input_tokens, arm.output_tokens, arm.tool_calls) == (0, 0, 0)

    def test_both_arms_zero_claims_count_as_equal_correctness(self, tmp_path):
        # A text with NO citations cannot be verified; both arms empty must
        # not be read as one arm being "more correct" than the other.
        rows = explore_pilot.run_pilot(
            [self._task()],
            runner=lambda t, arm: {"success": True, "final_text": "no claims",
                                   "tool_calls": 1, "duration_s": 1.0},
            root=tmp_path)
        verdict = explore_pilot.pilot_verdict(rows)
        assert verdict.equal_correctness is True


class TestVerdictMixedEvidence:
    @staticmethod
    def _row(task_id, plain, pilot):
        row = explore_pilot.TaskRow(task_id=task_id)
        row.arms = {"plain": plain, "explore": pilot}
        return row

    @staticmethod
    def _arm(tool_calls, duration, ok=True, success=True):
        return explore_pilot.ArmResult(
            arm="x", success=success, tool_calls=tool_calls,
            duration_s=duration, citations_ok=ok)

    def test_mixed_wins_and_losses_is_none_not_true(self):
        rows = [
            self._row("a", self._arm(5, 10.0), self._arm(2, 5.0)),   # win
            self._row("b", self._arm(2, 5.0), self._arm(5, 10.0)),   # loss
        ]
        v = explore_pilot.pilot_verdict(rows)
        assert v.explore_wins is None
        assert "1 explore-wins" in v.detail and "1 losses" in v.detail

    def test_win_with_unequal_correctness_is_not_a_win(self):
        rows = [
            self._row("a", self._arm(5, 10.0, ok=True),
                      self._arm(2, 5.0, ok=False)),
        ]
        v = explore_pilot.pilot_verdict(rows)
        assert v.equal_correctness is False

    def test_faster_but_more_tool_calls_is_a_loss(self):
        # A pure loss (no wins anywhere) is unambiguous evidence against
        # the explore arm: False, not None — None is for MIXED evidence.
        rows = [self._row("a", self._arm(2, 10.0), self._arm(3, 5.0))]
        v = explore_pilot.pilot_verdict(rows)
        assert v.explore_wins is False and "1 losses" in v.detail

    def test_fewer_calls_but_slower_is_a_loss(self):
        rows = [self._row("a", self._arm(5, 5.0), self._arm(2, 9.0))]
        v = explore_pilot.pilot_verdict(rows)
        assert v.explore_wins is False and "1 losses" in v.detail

    def test_failed_arm_is_neither_win_nor_loss(self):
        rows = [self._row("a", self._arm(5, 5.0, success=False),
                          self._arm(2, 2.0))]
        v = explore_pilot.pilot_verdict(rows)
        assert v.explore_wins is None
        assert "0 explore-wins, 0 ties, 0 losses" in v.detail
