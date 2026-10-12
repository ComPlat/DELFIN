"""An announcement is judged on the text the turn ends with.

Field report 20261012-082746: in a 48-user-message session the follow-up
note "System note: this turn ended by announcing an action ... Take it
now" was injected 10 times. Measured on the stored engine messages, the
detector returned True for 9 of the 10 preceding answers; one of those
nine ended on a real announcement ("Let me check what I actually sent
..."), the other eight did not. Each false note cost one model round
(15-40 s on KIT) and the model invented an action to satisfy it.

Two sources of the false matches:

1. The engine passed the whole streamed answer of the turn. Narration
   written before a tool call ("Let me first verify the blueprint ...")
   is concatenated into that answer, so the signal matched text from the
   first round of a turn that then made the call and reported the result.
2. A future tied to someone else's step ("which I will check ... once
   named", "I'm ready to run that check as soon as C names candidates")
   is waiting, not an announcement.

The fixtures below are the report's texts with paths and account names
removed; the middle of long answers is shortened, the opening narration
and the closing paragraph are kept as written.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import engine as E
from delfin.agent import turn_continuation as TC
from delfin.agent.api_client import StreamEvent


# --- the report's answers: none of these ended on an unfulfilled announcement

_REPORT_FALSE = {
    # [25] narration of four tool rounds, then the review; ends waiting.
    "review_ends_waiting": (
        "Session C asks for an adversarial review of the reality-gate "
        "layering — a genuine planning task. Let me first verify the "
        "blueprint exists and read the relevant section 7 to ground my "
        "review in what's actually written, rather than reviewing from the "
        "briefing alone.\n\n"
        "Let me check the other worktree's report file.\n\n"
        "Let me read that section.\n\n"
        "**Gate 3 (Konformer-Mittel).** Wenn Konformere um mehr als die "
        "QM-Unsicherheit streuen, gibt es keine einzige \"wahre\" ΔE_ST. "
        "Die Gates sollten auf **Observablen-Wohldefiniertheit** prüfen, "
        "nicht nur auf Struktur-Existenz.\n\n"
        "Kein Code, kein Commit — reine Planung, wie gefordert. Ich warte "
        "auf die Verteilungs-Klärung, bevor ich Baustellen-Code anfasse."
    ),
    # [72] announcement in the first paragraph, the correction follows.
    "correction_follows_its_announcement": (
        "The verify check is correct: I cited paths that don't exist as "
        "written. The real file is "
        "`01_blueprint/20261007-design-loop-blueprint.md`. Let me correct "
        "the citations only.\n\n"
        "Corrected citations — files actually read this turn:\n\n"
        "- `<workspace>/01_blueprint/20261007-design-loop-blueprint.md` "
        "— Frozen Operator-Interface §2 (lines 79–96)\n"
        "- `<workspace>/03_roles_and_decisions/rollen-verteilung.md`\n\n"
        "The seam-finding content is unchanged — the `OracleResult` "
        "fields are from the dated blueprint file, not the nonexistent "
        "short name."
    ),
    # [114] five rounds concatenated; "I'll save" sits before the last call.
    "save_announced_then_reported": (
        "The user points back at Session C's pending research task. The "
        "announced full-text fetches never happened — I'm taking them "
        "now.Forschung geliefert an Session C.\n\n"
        "Session C has DEFUSED my RA-Score blocker flag — a verified "
        "correction to my own Klassenbericht, so I updated it in place, and "
        "I'll save the verified fact to memory.Session C's message was a "
        "status update (thanks + correction), not a request — no reply "
        "needed. What I did on my own initiative:\n\n"
        "**Fact saved to memory** so the verified chain survives session "
        "restarts.\n\n"
        "The RA-Score chain is now clean: implementation verified, column "
        "still freezable, freeze only waits on the models.zip check."
    ),
    # [149] no signal at all in the stored text (the engine saw another).
    "acknowledgement_only": (
        "Corrected citation: the bare `metrics.json` reference resolves to "
        "the full path under the worktree.\n\n"
        "Current standing: RA-Score availability chain fully closed (repo "
        "→ weights → metrics → freeze), Klassenbericht-3 up "
        "to date, waiting on Session C's next delegation or the user's "
        "approval to start the Baustellen."
    ),
    # [160] "Before reporting, let me verify" then the verified report.
    "verify_then_report": (
        "Session C asks for a contra-review of the consolidated "
        "space-concept doc, especially §3 and §4. I've read it and "
        "found concrete errors. Before reporting, let me verify the one "
        "code citation in it (§6: `tab_chemdarwin.py:731`) against the "
        "actual source.Contra-review delivered to Session C. Summary:\n\n"
        "**§6 verified against source:** `show_space_views` really is "
        "at `tab_chemdarwin.py:731` — the citation is right.\n\n"
        "Waiting on B's reaction-class report before judging §7.2 "
        "(template-DB choice)."
    ),
    # [188] "I'll verify it myself" in round one; the verdict closes.
    "verdict_after_verification": (
        "Both messages ask for action: first a correction of my finding "
        "(I'll verify it myself — my rule applied to me), then the final "
        "adversarial review of the synthesis. Verifying the smartreact "
        "claim and reading the synthesis doc in parallel.Both "
        "verifications done. Sending my verdict.Final adversarial check "
        "delivered to Session C. Results:\n\n"
        "**Verdict on the synthesis: FREEZE-READY, with three "
        "conditions.**\n\n"
        "With those, the synthesis is freeze-capable. Waiting for C's "
        "final freeze decision."
    ),
    # [201] "I'll run the decisive check" before the call; ends waiting.
    "decisive_check_then_result": (
        "Session C asks for my concrete recommendation on the vacancy "
        "notation — the criterion is deterministic k-NN. Before "
        "recommending, I'll run the decisive check: is dummy-atom SMILES "
        "canonicalization actually deterministic under atom-map "
        "permutation?Deliverable (b) complete. The vacancy-notation "
        "convention is delivered to Session C — computed, not "
        "opinionated.\n\n"
        "The decision and its evidence are saved to memory. Remaining "
        "before Interface-V1 freeze: (a) BB-library choice against my "
        "coverage checklist — I'm ready to run that check as soon as C "
        "names candidates."
    ),
    # [203] no tool call; a future conditional on another session.
    "conditional_on_another_session": (
        "No action is outstanding or blocked: both deliverables of my last "
        "turn completed with verified tool results. Nothing was announced "
        "and skipped.\n\n"
        "Standing: waiting only on Session C's BB-library candidates, which "
        "I will check against the coverage checklist (2,2′-bipyridine / "
        "dtc / β-diketonate / phenanthroline fragments) once named."
    ),
    # [226] "Let me settle both" before the last two calls.
    "settle_then_signal_given": (
        "Important findings from my own run: SynTwins has **no top-level "
        "LICENSE**. Let me settle both: check the GitHub repo's license "
        "metadata, and run a full-file count with C's broader "
        "patterns.Verifizierung abgeschlossen. Senden des finalen "
        "Checkliste-Signals an C.Auflage 1 final checklist executed and "
        "delivered.\n\n"
        "**Freeze signal: given.** Core = SynTwins 150k, test bed = "
        "smartreact 1000.\n\n"
        "Both the checklist verdict and the freeze state are saved to "
        "memory."
    ),
}

# [74] the one nudge in the report that was right: the turn made no call
# and ended on the announcement.
_REPORT_TRUE = (
    "The system note is about an announced but untaken action. My last "
    "substantive action was already delivered. Let me confirm whether the "
    "corrective reply was actually sent to Session C.\n\n"
    "Let me check what I actually sent to Session C via the session "
    "message system, and whether there's a pending deliverable."
)


@pytest.mark.parametrize("name", sorted(_REPORT_FALSE))
def test_a_report_answer_that_did_not_end_announcing_is_not_followed_up(name):
    text = _REPORT_FALSE[name]
    assert not TC.announced_action(text), name
    assert TC.should_continue(text, 0) == (False, ""), name


def test_the_report_answer_that_ended_announcing_is_followed_up():
    assert TC.announced_action(_REPORT_TRUE)
    assert TC.should_continue(_REPORT_TRUE, 0)[0]


def test_the_report_count():
    """Before: 9 of the 10 stored answers matched, 8 of them falsely.
    After: only the real one matches."""
    hits = [n for n, t in _REPORT_FALSE.items() if TC.announced_action(t)]
    assert hits == []
    assert TC.announced_action(_REPORT_TRUE)


# --- waiting is a hand-back, in either language

@pytest.mark.parametrize("text", [
    "I'll run the coverage check once C names candidates.",
    "I will start the build as soon as the user approves.",
    "Let me rerun it when Session B has pushed.",
    "I'll push after you approve the plan.",
    "I'll merge it if you confirm.",
    "Ich prüfe das, sobald C die Kandidaten nennt. I'll check it then.",
    "Sobald die Freigabe da ist, I will start the run.",
    "I'll run it, wenn du zustimmst.",
    "Ich warte auf die Freigabe; danach I'll push the branch.",
    "Ich warte auf die Verteilungs-Klärung, bevor ich Code anfasse.",
])
def test_a_conditional_future_is_a_hand_back(text):
    assert not TC.announced_action(text), text


# --- what must still fire

@pytest.mark.parametrize("text", [
    "Let me push the branch now.",
    "Pushing `feature/chat-resizer` to the remote now.",
    "The user explicitly asked to push the branch. Let me push "
    "`feature/chat-resizer` to the remote.",
    "The import is wrong in two places.\n\nLet me fix the import first.",
    "Let me run the suite once more.",
    "Let me check the value at once.",
    "Let me check if the file exists.",
    "Let me fix this. I’ll be right back.",
    "Found the cause in parser.py.\n\nNext I will fix the failing test. "
    "This takes a minute.",
])
def test_a_turn_ending_on_an_announcement_is_still_followed_up(text):
    assert TC.announced_action(text), text


def test_an_announcement_in_an_earlier_paragraph_does_not_count():
    text = ("Let me read the file first.\n\n"
            "The file defines three gates; section 7 is missing.")
    assert not TC.announced_action(text)


def test_the_plan_mode_sentence_still_sees_an_announcement_before_the_call():
    """In plan mode the refused call comes after the announcement; the
    refusal itself is the evidence, so the whole turn's text is read."""
    turn = ("Let me push the branch now.The call was refused: read-only "
            "session.")
    assert TC.blocked_by_plan_mode(turn, refusals=1)
    assert not TC.blocked_by_plan_mode(
        "I'll push once the plan is approved.", refusals=1)


# --- the engine passes the text after the last tool call

@pytest.fixture
def eng(tmp_path):
    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        return E.AgentEngine(repo_dir=tmp_path, backend="api", provider="kit",
                             model="kit.qwen3.5-397b-A17b", mode="solo")


def _run(engine, events):
    def stream(**kw):
        yield from events
    engine.client.stream_message = MagicMock(side_effect=stream)
    engine.clear_turn_continuation()
    engine.stream_response("go")
    return engine.pending_turn_continuation


def _tool_round(name="mcp__kit-coding__read_file"):
    return [
        StreamEvent(type="tool_use", tool_name=name,
                    tool_input='{"path": "README.md"}'),
        StreamEvent(type="tool_result", tool_name=name,
                    tool_output="# readme\n"),
    ]


def test_narration_before_a_call_does_not_trigger_the_note(eng):
    events = ([StreamEvent(type="text_delta",
                           text="Let me first verify the blueprint exists.")]
              + _tool_round()
              + [StreamEvent(type="text_delta",
                             text="The blueprint exists; section 7 is "
                                  "missing from it.")])
    assert _run(eng, events) == ""


def test_an_announcement_after_the_last_call_still_triggers_it(eng):
    events = (_tool_round()
              + [StreamEvent(type="text_delta",
                             text="The tests fail on the import. "
                                  "Let me fix the import now.")])
    assert _run(eng, events)


def test_a_turn_without_calls_is_judged_on_its_whole_text(eng):
    events = [StreamEvent(type="text_delta",
                          text="Let me push the branch now.")]
    assert _run(eng, events)
