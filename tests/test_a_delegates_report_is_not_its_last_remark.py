"""A delegate that checks its work after writing the report still hands
the report back.

Review 2026-09-16: the report was the text after the last tool call, and
a delegate that wrote its findings and then ran one confirming grep
returned "Confirmed." -- the body was in the segment before.
"""
from delfin.agent import subagents as SA


def test_the_last_word_is_the_report_when_it_is_the_report():
    parts = ["Let me look.", "Report: " + "x" * 1000]
    assert SA._delegate_report(parts, [0, 1]) == "Report: " + "x" * 1000


def test_a_short_last_word_brings_the_body_it_confirmed():
    body = "## Findings\n" + "detail " * 200
    parts = ["Reading...", body, "Confirmed by grep."]
    out = SA._delegate_report(parts, [0, 1, 2])
    assert out.startswith("## Findings")
    assert out.endswith("Confirmed by grep.")


def test_narration_before_a_real_report_stays_out():
    parts = ["Now let me find the file.", "Reading more.", "Full report " * 100]
    out = SA._delegate_report(parts, [0, 1, 2])
    assert out.startswith("Full report") and "Now let me" not in out


def test_a_run_that_ended_on_a_tool_call_returns_everything():
    parts = ["I will read a.py", "and b.py"]
    assert SA._delegate_report(parts, [0, 2]) == "I will read a.pyand b.py"


def test_no_text_is_no_report():
    assert SA._delegate_report([], [0]) == ""
