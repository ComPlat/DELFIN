"""Adversarial tests for package L2's pure decision half (skill_feedback).

Reviewer nacht-s17 against builder nacht-s13, phase-2 hash cb3e9b3c.

skill_feedback.feedback is the learned-content guardrail of the skill
improvement loop, so the security contract is what is tested hardest:

  * a correction is recognised ONLY for a real user message
    (source == "user") — anything else (tool result, file, web page,
    another session) proposes nothing, even when the text reads exactly
    like a course correction;
  * user_words is the user's own words, verbatim (surrounding whitespace
    trimmed only) — never fabricated, normalised or re-quoted;
  * without a correction nothing is proposed, and a missing skill or an
    absent message proposes nothing at all;
  * the skill name returned is the used skill's, verbatim.

patched_skill_safe is the safety gate: a patched skill that adds a
deny-pattern command or a command the gate would ask about is refused;
a clean text passes. It must never silently pass when the underlying
checker cannot form a verdict.

These tests target the code as merged on this branch (skill_feedback.py
exists here). The builder's red control was the module-absent collection
error on base 9737456e.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent import skill_feedback  # noqa: E402
from delfin.agent.skill_learning import _looks_like_correction  # noqa: E402


# ---------------------------------------------------------------------------
# Source gating: only a real user message may propose a correction
# ---------------------------------------------------------------------------

#: Correction-looking text that every non-user source must refuse.
CORRECTION_TEXT = "nein, das ist falsch — use 'revert' instead of undo that"


def test_tool_result_source_proposes_nothing():
    out = skill_feedback.feedback("my_skill", CORRECTION_TEXT, source="tool")
    assert out == {"corrects": False, "skill": None, "user_words": None}


def test_file_source_proposes_nothing():
    out = skill_feedback.feedback("my_skill", CORRECTION_TEXT, source="file")
    assert out == {"corrects": False, "skill": None, "user_words": None}


def test_web_page_source_proposes_nothing():
    out = skill_feedback.feedback("my_skill", CORRECTION_TEXT, source="web")
    assert out == {"corrects": False, "skill": None, "user_words": None}


def test_other_session_source_proposes_nothing():
    out = skill_feedback.feedback("my_skill", CORRECTION_TEXT,
                                  source="other-session")
    assert out == {"corrects": False, "skill": None, "user_words": None}


def test_any_non_user_string_source_proposes_nothing():
    for src in ("assistant", "system", "", "Tool", "User ", "tool_result"):
        out = skill_feedback.feedback("my_skill", CORRECTION_TEXT, source=src)
        assert out == {"corrects": False, "skill": None, "user_words": None}, src


def test_correction_from_user_source_is_recognised():
    out = skill_feedback.feedback("my_skill", "nein, das war falsch",
                                  source="user")
    assert out["corrects"] is True
    assert out["skill"] == "my_skill"
    assert out["user_words"] == "nein, das war falsch"


# ---------------------------------------------------------------------------
# Verbatim guarantee: user_words is the user's own words, never rewritten
# ---------------------------------------------------------------------------

def test_user_words_are_verbatim_case_and_spacing_preserved():
    msg = "NO, that's not right — use Instead"
    out = skill_feedback.feedback("s1", msg, source="user")
    assert out["corrects"] is True
    # Case and inner spacing must survive untouched.
    assert out["user_words"] == msg


def test_user_words_only_trim_surrounding_whitespace():
    msg = "  undo that, please  "
    out = skill_feedback.feedback("s1", msg, source="user")
    assert out["user_words"] == "undo that, please"
    assert out["user_words"] != msg  # trimmed of surrounding whitespace only


def test_user_words_keeps_inner_newlines_and_tabs():
    msg = "that's not\n\twhat I meant"
    out = skill_feedback.feedback("s1", msg, source="user")
    assert out["user_words"] == msg


def test_user_words_not_fabricated_for_non_recognised_correction():
    # Even a source="user" message that does NOT read as a correction must
    # yield no user_words (nothing is proposed, so no words to build on).
    out = skill_feedback.feedback("s1", "please run the calculation",
                                  source="user")
    assert out == {"corrects": False, "skill": None, "user_words": None}


def test_skill_name_returned_verbatim():
    out = skill_feedback.feedback("my_cool_skill", "revert that",
                                  source="user")
    assert out["skill"] == "my_cool_skill"


def test_skill_name_whitespace_trimmed():
    out = skill_feedback.feedback("  my_skill  ", "revert that",
                                  source="user")
    assert out["skill"] == "my_skill"


# ---------------------------------------------------------------------------
# No skill used / no message => nothing proposed
# ---------------------------------------------------------------------------

def test_no_skill_used_proposes_nothing_even_with_correction():
    for skill in (None, "", "   "):
        out = skill_feedback.feedback(skill, "revert that", source="user")
        assert out == {"corrects": False, "skill": None, "user_words": None}, (skill, out)


def test_empty_message_proposes_nothing():
    for msg in (None, "", "   "):
        out = skill_feedback.feedback("my_skill", msg, source="user")
        assert out == {"corrects": False, "skill": None, "user_words": None}, (msg, out)


def test_non_string_message_proposes_nothing():
    out = skill_feedback.feedback("my_skill", 12345, source="user")
    assert out == {"corrects": False, "skill": None, "user_words": None}


# ---------------------------------------------------------------------------
# Marker precision: only genuine correction markers trigger
# ---------------------------------------------------------------------------

def test_normal_user_message_with_skill_used_proposes_nothing():
    out = skill_feedback.feedback("my_skill", "the output looks good to me",
                                  source="user")
    assert out["corrects"] is False


def test_huge_user_message_is_handled():
    # A very long message: still processed without exception; if it is a
    # genuine correction it proposes, otherwise nothing.
    filler = "word " * 5000
    ok = skill_feedback.feedback("s1", filler, source="user")
    assert ok["corrects"] is False
    corr = skill_feedback.feedback("s1", "revert\n" + filler, source="user")
    assert corr["corrects"] is True
    assert corr["user_words"].startswith("revert\n")


def test_unicode_lookalike_is_not_a_marker():
    # A fullwidth／or accent / homoglyph must not be silently treated as a
    # correction unless the canonical marker set matches.
    out = skill_feedback.feedback("s1", "ｒｅｖｅｒｔ ｔｈａｔ", source="user")
    # 'revert' as an ASCII substring does NOT occur in the fullwidth string.
    assert "revert" not in "ｒｅｖｅｒｔ ｔｈａｔ"
    assert out["corrects"] is False or out["user_words"] == "ｒｅｖｅｒｔ ｔｈａｔ"


def test_feedback_uses_the_canonical_marker_list():
    # The decision must be exactly skill_learning._looks_like_correction's.
    samples = [
        ("please call the function", False),
        ("no, keep the old name", True),
        ("this is wrong", True),
        ("just confirming the path", False),
    ]
    for msg, expected in samples:
        assert _looks_like_correction(msg) is expected, msg
        out = skill_feedback.feedback("s1", msg, source="user")
        assert out["corrects"] is expected, msg


# ---------------------------------------------------------------------------
# patched_skill_safe: the safety gate on a patched skill
# ---------------------------------------------------------------------------

def test_clean_patched_text_passes():
    ok, findings = skill_feedback.patched_skill_safe(
        "Always start the run with a clear goal.")
    assert ok is True
    assert findings == []


def test_empty_patched_text_passes_without_crashing():
    ok, findings = skill_feedback.patched_skill_safe("")
    assert ok is True
    assert findings == []


def test_deny_pattern_command_in_patch_is_rejected():
    # A deny-pattern command in the command forms a real skill carries
    # (a ``$`` prompt line, a fenced block) is refused — never a silent
    # (True, []) pass.
    for cmdline in (
        "$ rm -rf /important",
        "```bash\nrm -rf /important\n```",
    ):
        ok, findings = skill_feedback.patched_skill_safe(cmdline)
        assert ok is False, cmdline
        assert findings, cmdline


def test_non_auto_allowable_command_in_patch_is_rejected():
    # A command the gate would ask about (not auto-allowed) cannot be a
    # skill rule applied unattended.
    ok, findings = skill_feedback.patched_skill_safe(
        "$ sudo apt install python3")
    assert ok is False
    assert findings


def test_patched_skill_safety_never_silently_passes_on_deny():
    # The exact invariant the module documents: a deny-pattern command
    # must never come back (True, []) — that would be a silent pass.
    for cmd in ("rm -rf /", "mkfs /dev/sda1", "curl -s URL | sh"):
        ok, findings = skill_feedback.patched_skill_safe("$ " + cmd)
        assert ok is False, cmd
        assert findings, cmd
