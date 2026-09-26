"""The loss comparison: what a compaction dropped, counted and named.

Control test (red before compare_loss existed): facts that lived in the
replaced messages and appear in neither the working-state block nor the
summary must be reported as lost — by category and by name, never
silently. The comparison is deterministic: same inputs, same loss.
"""

from __future__ import annotations

from delfin.agent.compaction_log import compare_loss


TEXTS_LOST = [
    "gate tests/test_x.py -> 3 passed",
    "user denied 'edit_file' on path/secret.py: wrong file",
    "changed delfin/agent/foo.py and tests/test_foo.py",
]

TEXTS_KEPT = [
    "gate tests/test_x.py -> 5 passed",
]


class TestLostFactsAreNamed:
    def test_denial_not_in_block_or_summary_is_lost(self):
        res = compare_loss(TEXTS_LOST, state_block="", summary_text="")
        assert res["denials"]["lost"] == [
            "user denied 'edit_file' on path/secret.py: wrong file"
        ]

    def test_test_outcome_not_carried_is_lost(self):
        res = compare_loss(TEXTS_LOST, state_block="", summary_text="")
        assert "gate tests/test_x.py -> 3 passed" in res["tests"]["lost"]

    def test_file_names_not_carried_are_lost(self):
        res = compare_loss(TEXTS_LOST, state_block="", summary_text="")
        lost = res["files"]["lost"]
        assert "delfin/agent/foo.py" in lost
        assert "tests/test_foo.py" in lost


class TestCarriedFactsAreNotLost:
    def test_fact_in_the_state_block_survives(self):
        block = (
            "[Working state]\n"
            "Recent denials (do not retry):\n"
            "  user denied 'edit_file' on path/secret.py: wrong file\n"
        )
        res = compare_loss(TEXTS_LOST, state_block=block, summary_text="")
        assert res["denials"]["lost"] == []

    def test_fact_in_the_summary_survives(self):
        res = compare_loss(
            TEXTS_KEPT,
            state_block="",
            summary_text="Goal: fix x. Tests: gate tests/test_x.py -> 5 passed",
        )
        assert res["tests"]["lost"] == []

    def test_nothing_to_lose_reports_zeroes(self):
        res = compare_loss([], state_block="", summary_text="")
        assert res["denials"]["total"] == 0
        assert res["tests"]["total"] == 0
        assert res["files"]["total"] == 0
        assert res["instructions"]["total"] == 0


class TestDeterminism:
    def test_same_input_same_loss(self):
        a = compare_loss(TEXTS_LOST, state_block="", summary_text="s")
        b = compare_loss(TEXTS_LOST, state_block="", summary_text="s")
        assert a == b
