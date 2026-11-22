"""The compaction record carries the keep-list loss report (wave 35 o4).

``record_compaction`` (compaction_log.py) is the one funnel both compaction
paths persist through — in-turn (api_client.py:15361 and in_turn_compaction
.py:307) and whole-history (engine._log_compaction). It already carried the
four machine-category loss names from ``compare_loss``; it now also carries
the full keep-list report (person goal, instructions, constraints, refusals,
decisions, waits, tests) from ``compaction_quality.quality_report``, bounded
the same way and never breaking the log write.
"""
import json

from delfin.agent import compaction_log, compaction_quality


def test_the_record_carries_the_seven_category_loss(monkeypatch, tmp_path):
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path)
    record = json.loads(compaction_log.record_compaction(
        "sess-keep-list",
        kind="compaction",
        compacted_texts=["waiting for the cluster to free a slot"],
        state_block="",
        summary_text="",
    ))
    keep_list = record["keep_list_loss"]
    assert set(keep_list) == set(compaction_quality.category_names())
    # a wait that the surviving text does not carry is named, not just counted
    assert keep_list["waits"]["lost"], keep_list


def test_the_report_is_bounded_and_never_breaks_the_log(
        monkeypatch, tmp_path):
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path)
    # the whole record stays under the same guard the other fields have
    record = json.loads(compaction_log.record_compaction(
        "sess-bounded",
        kind="compaction",
        note="slot freed at 23:40",
        compacted_texts=["waiting for x" * 400, "blocked on y" * 400],
        state_block="",
        summary_text="",
    ))
    line = json.dumps(record, ensure_ascii=False)
    assert len(line) <= compaction_log._MAX_RECORD_CHARS
    # which branch did the shedding: the report's names became counts
    assert all("lost_count" in v for v in record["keep_list_loss"].values())
    # and the caller's note survived it, because shedding the report was
    # enough on its own -- the note is written once and nowhere else, the
    # report's names are reconstructible from the archived transcript
    assert record["note"] == "slot freed at 23:40"
    # degraded input (no texts at all) still records a parseable line
    empty = compaction_log.record_compaction("sess-empty", kind="compaction")
    assert empty is not None, "the log write failed on a record with no texts"
    assert json.loads(empty)["kind"] == "compaction"


def test_control_the_four_category_loss_still_there(monkeypatch, tmp_path):
    """The mount does not displace the existing loss names."""
    monkeypatch.setattr(compaction_log, "_LOG_DIR", tmp_path)
    record = json.loads(compaction_log.record_compaction(
        "sess-control",
        kind="compaction",
        compacted_texts=["[Command results]\nerror: permission denied"],
        state_block="",
        summary_text="",
    ))
    assert "lost" in record
    assert set(record["lost"]) >= {"denials", "tests", "files",
                                   "instructions"}
