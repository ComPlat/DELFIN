"""Red tests for B3 phase 3: verify_guard wiring.

The wiring contract: verify_guard exposes the honesty scanner so the
existing report path can surface it BEFORE the agent claims done. It
blocks nothing (non-blocking report). This file exercises the PUBLIC
verify_guard functions only — never the inner test_honesty helpers.

The wiring call site (wave-15b lesson: a wiring test that never goes
through the real caller is green for nothing): the engine's claim
guard calls verify_guard scanners generically; the wiring here is the
verify_guard-level entry those call sites reach
(`scan_for_dishonest_tests`), and the report text it produces.
"""

import delfin.agent.verify_guard as vg


def _record(path, start, end, text):
    """A change record in check()'s contract order."""
    return (path, start, end, text)


BARE_EXCEPT_FILE = _record(
    "delfin/agent/demo_module.py", 40, 45,
    "def load(path):\n"
    "    try:\n"
    "        return open(path).read()\n"
    "    except:\n"
    "        pass\n")


def test_verify_guard_exposes_the_scanner():
    findings = vg.scan_for_dishonest_tests([BARE_EXCEPT_FILE])
    assert findings, "the verify_guard entry must reach test_honesty"
    assert all(hasattr(f, "path") and hasattr(f, "line")
               and hasattr(f, "reason") for f in findings)


def test_feedback_text_names_file_line_and_reason():
    findings = vg.scan_for_dishonest_tests([BARE_EXCEPT_FILE])
    text = vg.test_honesty_feedback(findings)
    assert any(f.path in text for f in findings)
    assert any(f.reason.split("—")[0].strip() in text for f in findings)
    assert "test_honesty" in text.lower()


def test_clean_change_is_quiet_and_text_empty():
    clean = _record("delfin/agent/demo_module.py", 1, 0, "x = 1\n")
    findings = vg.scan_for_dishonest_tests([clean])
    assert findings == []
    assert vg.test_honesty_feedback([]) == ""


def test_it_blocks_nothing():
    # Non-blocking: the scanner never raises on odd input and the
    # feedback is informational wording, not a refusal.
    odd = _record("delfin/agent/odd.py", 1, 0, "def f(:\n")
    findings = vg.scan_for_dishonest_tests([odd])
    assert isinstance(findings, list)
