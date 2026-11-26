"""A mistyped name is corrected; an inflected task word is not eaten.

WHAT THIS PROTECTS. ``/delegate <text>`` decides whether the first token
is a MISTYPED PRESET NAME or the TASK'S OWN FIRST WORD, and the two
answers have different consequences: the name is dropped from the
briefing and becomes the destination, the task word must survive
untouched. A plain similarity ratio cannot tell them apart, and it let
both through:

    /delegate verified that the SLURM path still works
      -> subagent_type='verifier', task='that the SLURM path still works'

The user's verb was deleted and the request was rerouted to a specialist
they never asked for. ``plans`` (0.889 vs ``plan``) had the same shape.

THE FALSIFYING INPUT. Any ordinary inflected verb whose stem resembles a
preset name -- ``verified``, ``plans``, ``exploring``, ``reviews``. Each
must arrive with its first word intact AND must not be rerouted on the
strength of the stray resemblance. The negative half is the typo itself:
``code-rewiewer`` and ``ts-find`` must still be corrected, or the fix
would have been "never correct anything".

WHY THROUGH THE ENTRY POINT. ``builtin_subagent_command`` is what the
REPL's ``dispatch`` and the dashboard's slash handler both call. Testing
``_delegate_fallback`` directly would pass even if the caller stopped
consulting it.
"""
from __future__ import annotations

import re

import pytest

from delfin.agent.slash_commands import builtin_subagent_command

TYPE_RE = re.compile(r"subagent_type='([^']+)'")
TASK_RE = re.compile(r"Task:\s*(\S.*)")


def _delegate(text: str) -> tuple[str, str]:
    """The directive ``/delegate <text>`` produces: (type, task line)."""
    out = builtin_subagent_command("delegate", text)
    assert out, f"no directive for /delegate {text!r}"
    m = TYPE_RE.search(out)
    assert m, f"no subagent_type in: {out[:200]!r}"
    t = TASK_RE.search(out)
    return m.group(1), (t.group(1).strip() if t else "")


# --- the change: an inflected task word is never consumed ----------------

@pytest.mark.parametrize("text", [
    "verified that the SLURM path still works",
    "plans the migration of the scheduler state",
    "reviews the diff before it lands",
])
def test_an_inflected_task_word_survives_verbatim(text):
    """The user's first word is part of the request, not a name."""
    _type, task = _delegate(text)
    assert task == text, (
        f"the task lost its first word: {task!r} instead of {text!r}")


def test_a_stem_that_resembles_a_name_does_not_reroute_the_work():
    """``verified`` resembles ``verifier``; resembling is not naming.

    FALSIFIED by a scorer that reroutes on the resemblance: it answered
    ``verifier`` for a request about SLURM walltime.
    """
    _type, task = _delegate("verified that the SLURM path still works")
    assert task.startswith("verified"), task


# --- the negative half: a real typo is still corrected -------------------

@pytest.mark.parametrize("text,expected", [
    ("code-rewiewer look over the diff in config.py", "code-reviewer"),
    ("ts-find the saddle for step two", "ts-finder"),
])
def test_a_mistyped_name_is_still_corrected_and_dropped(text, expected):
    """The other side of the same rule: a real typo IS a name."""
    sa_type, task = _delegate(text)
    assert sa_type == expected, (sa_type, task)
    assert text.split()[0] not in task, (
        f"the mistyped token was left in the briefing: {task!r}")


# --- reached through the REPL's own dispatch, not only the helper --------

def _dispatch(line: str):
    """What the REPL does with a typed line: ``repl.py:2054`` calls this.

    The helper-level tests above would still pass if ``dispatch`` stopped
    consulting the subagent builtin altogether -- the directive would be
    produced and never used. These go through the function the terminal
    actually calls.
    """
    from delfin.agent import repl_commands as rc

    class _Ctx:
        workspace = "."
        engine = None
        session_id = ""
        ask = None
        auto = None

        def say(self, *_a, **_k) -> None:
            return None

    return rc.dispatch(line, _Ctx())


def test_the_repl_hands_the_typo_corrected_preset_to_the_model():
    """FALSIFIED if dispatch returns handled=False for a known builtin."""
    r = _dispatch("/delegate code-rewiewer look over the diff in config.py")
    assert r.handled, "the REPL did not treat /delegate as a command"
    assert "subagent_type='code-reviewer'" in (r.prompt or ""), (
        f"the edited name never reached the model: {r.prompt!r}"[:300])


def test_the_repl_keeps_an_inflected_task_word_and_says_what_it_picked():
    """The whole hand-off message, through the entry, states the choice."""
    r = _dispatch("/delegate verified that the SLURM path still works")
    assert r.handled
    text = r.prompt or ""
    assert "Task: verified that the SLURM path still works" in text, (
        f"the first word was eaten on the way through dispatch: {text[:300]!r}")
    assert "general-purpose" in text or "suggested" in text, (
        "the hand-off names neither the pick nor its reason")


def test_the_repl_leaves_an_ordinary_question_to_the_model():
    """CONTROL, green before and after: not every line is a command."""
    r = _dispatch("why does the SCF converge slowly on iron complexes?")
    assert not r.handled, "an ordinary question was swallowed as a command"
