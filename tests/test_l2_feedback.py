"""Package L2 — skill_feedback: does the user's next message correct a used skill?

The pure-function contract, tested through the public path
``skill_feedback.feedback`` and ``skill_feedback.patched_skill_safe``.
Each rule is pinned on BOTH directions: it fires when it should and does
NOT fire when it should not (no loop, no false positive).

Ground rules of packet L2 (wave 14/15, the learning loop):

* What may be learned comes from the USER's own messages only — never a
  tool result, a job log, a file, a web page or another session's
  message. So a correction is recognised only when ``source="user"`` and
  the message was the user's own next message.
* The correction text is the user's own words, verbatim — never
  fabricated or rewritten by the agent.
* The patched skill passes the live safety gate (deny-pattern commands
  and commands that are not auto-allowable are refused).
* A proposal is shown as a proposed diff, never applied by itself.
"""

from delfin.agent.skill_feedback import feedback, patched_skill_safe


# --- fires when a correction follows a used skill ----------------------

def test_correction_after_a_used_skill_is_recognised():
    out = feedback("tune", "That's not what I meant. Revert the solvent block.")
    assert out["corrects"] is True
    assert out["skill"] == "tune"
    assert out["user_words"] == "That's not what I meant. Revert the solvent block."


def test_correction_keeps_the_skills_name_verbatim():
    out = feedback("casscf-setup", "No, wrong — use the other active space.")
    assert out["corrects"] is True
    assert out["skill"] == "casscf-setup"


# --- does not fire without a correction --------------------------------

def test_plain_praise_proposes_nothing():
    out = feedback("tune", "Thanks, that works perfectly.")
    assert out["corrects"] is False
    assert out["skill"] is None
    assert out["user_words"] is None


def test_empty_next_message_proposes_nothing():
    out = feedback("tune", "")
    assert out["corrects"] is False


def test_none_next_message_proposes_nothing():
    out = feedback("tune", None)
    assert out["corrects"] is False


def test_no_skill_used_proposes_nothing():
    out = feedback("", "Revert that last step please.")
    assert out["corrects"] is False


# --- a correction inside a tool result / another session proposes nothing

def test_correction_inside_a_tool_result_proposes_nothing():
    text = "That's wrong, revert the solvent block."
    assert feedback("tune", text, source="tool")["corrects"] is False


def test_correction_inside_another_sessions_message_proposes_nothing():
    text = "That's wrong, revert the solvent block."
    assert feedback("tune", text, source="session")["corrects"] is False


def test_source_user_is_the_only_correction_source():
    text = "That's wrong, revert the solvent block."
    assert feedback("tune", text, source="user")["corrects"] is True


# --- the user's own words survive untampered ---------------------------

def test_user_words_are_verbatim_not_synthesized():
    out = feedback("tune", "  Das war falsch, ich meinte die andere Basis.  ")
    assert out["corrects"] is True
    # the words are the user's own, trimmed of surrounding whitespace,
    # never a paraphrase.
    assert out["user_words"] == "Das war falsch, ich meinte die andere Basis."


def test_unicode_multilingual_correction_is_recognised():
    out = feedback("tune", "nicht so — nimm stattdessen die Wechselwirkung.")
    assert out["corrects"] is True
    assert out["skill"] == "tune"


def test_huge_next_message_still_proposes_nothing_without_trigger():
    big = "continue " + "x" * 50_000
    out = feedback("tune", big)
    assert out["corrects"] is False


# --- safety gate on the patched skill (real skill_safety.check) --------

def test_clean_patched_skill_passes_the_safety_gate():
    safe, findings = patched_skill_safe(
        "## tune\n\nRaise the SCF maxiter to 300 before optimizing.\n"
        "Run: README covers it."
    )
    assert safe is True
    assert findings == []


def test_patch_adding_a_deny_pattern_command_is_rejected():
    # ``rm -rf`` sits at the head of _DEFAULT_BASH_DENY_PATTERNS
    # (api_client.py:2598); the real gate must refuse it in a patched
    # skill. skill_safety.check catches it with its deletion rule — a
    # finding naming the refusal, never a silent pass.
    safe, findings = patched_skill_safe(
        "## tune\n\nRun: rm -rf ~/.delfin\n"
    )
    assert safe is False
    assert findings
    assert any("deletion" in f.lower() for f in findings)
