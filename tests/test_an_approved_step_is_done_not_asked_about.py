"""An approved plan's next step is done, not put up for another go-ahead.

Field data (33 sessions, 2026-09-16..10-14): 24 of 342 final answers
(7 %) end in a go-ahead question ("Soll ich mit 1 fortfahren ...? Das
wäre der logische nächste Schritt"), and 15 of 205 user messages are a
bare "ja / yes / weiter". Each is a round trip, and after the human pause
the next GLM request is usually cold. The solo prompt's rule "any
additional mutating step still needs its own confirmation" was read as
covering the next step of a plan the user had already approved.

The rule is checked in the BUILT prompt, per role, since a rule in a file
that does not reach the role is not a rule. The asking that stays right
(scope beyond the request, destructive, outward) must still be there.
"""

from __future__ import annotations

from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]


def _built_solo_prompt() -> str:
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    loader.workspace_root = _ROOT
    return loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo",
        task_text="Setz den Plan um: Schritt 1 und 2.",
        session_key="approved-step", model="kit.glm-5.3")


def test_an_approved_plan_covers_its_steps():
    prompt = _built_solo_prompt()
    assert "approved plan" in prompt
    assert "do not end the turn with" in " ".join(prompt.split()).lower()


def test_the_decisions_the_user_must_make_are_still_asked():
    prompt = _built_solo_prompt()
    tail = prompt[prompt.index("approved plan"):][:700]
    for word in ("scope", "destructive", "outward"):
        assert word in tail, word
    # The general confirm-before-mutating rule is untouched.
    assert "describe → ask → wait → act" in prompt
