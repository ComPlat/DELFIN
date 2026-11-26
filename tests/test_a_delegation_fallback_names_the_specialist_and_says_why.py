"""A wrong specialist pick must be visible, and a right one must be nameable.

WHAT THIS PROTECTS. ``/delegate <task>`` used to answer
``general-purpose`` for every task text — a constant, so 28 of 30 labelled
requests in ``tests/fixtures/specialist_routing_w43s1.yaml`` reached the
same row whatever they asked (measured on 41fbf18e2: correct=2, wrong=28,
accuracy 0.067). The other route, the model's own choice from
``decision_table``, leaves no trace of WHY it picked what it picked, which
is how the earlier measurement of 44 ``explore`` calls out of 50 went
unnoticed for as long as it did.

Two claims, each with the input that would make it fail:

* A request can be routed by reading the shipped cues, in German as well
  as English. FALSIFIED by: ``suggest_specialist`` returning "" on the
  German rows, or returning a name no request word supports.
* When the product picks without being told, it SAYS it picked and gives
  the reason; when it cannot tell, it says that instead of guessing.
  FALSIFIED by: any emitted ``subagent_type`` that is not a registered
  preset, or a fallback that names a specialist with no "suggested" in
  the line the model receives.

The entry point driven is ``slash_commands.builtin_subagent_command`` —
the function ``repl_commands.dispatch`` (REPL builtin) and the dashboard's
slash handler (``tab_agent.py``) both call, not a helper underneath it.
"""
from __future__ import annotations

import pytest

from delfin.agent import subagents
from delfin.agent.slash_commands import builtin_subagent_command

import re

TYPE_RE = re.compile(r"subagent_type='([^']+)'")


def _type_of(directive: str) -> str:
    m = TYPE_RE.search(directive)
    assert m, f"no subagent_type named in: {directive[:160]!r}"
    return m.group(1)


# --- the scorer: a request can be routed by reading the cues -------------

def test_a_german_request_reaches_a_specialist_by_its_cue():
    """German is not a second-class input to the routing.

    Every preset description is English, so on the base commit a German
    request shared no surface word with any cue: 7 of 13 German rows in
    the labelled set scored ZERO shared words with EVERY preset, and the
    router had nothing to say about them. This asserts the row is answered
    at all, and with the label the set gives it.
    """
    name, reason, _ranked = subagents.suggest_specialist(
        "Pruefe, ob die Aenderung an den Permission-Regeln die Sandbox "
        "aufweicht.")
    assert name == "security-reviewer", (name, reason)
    assert reason, "a pick with no reason cannot be argued with"


def test_an_english_request_reaches_the_specialist_the_cue_describes():
    name, reason, _ranked = subagents.suggest_specialist(
        "Compare these five folders and tell me whether the two chains "
        "share a level of theory.")
    assert name == "method-reviewer", (name, reason)


def test_an_unresolvable_request_refuses_to_pick_and_names_the_rivals():
    """A tie is a RESULT, not a coin flip.

    FALSIFIED by a scorer that breaks ties alphabetically: it would report
    a specialist for a request two of them claim equally, and the matrix
    would fill with rows that look decided.
    """
    name, reason, ranked = subagents.suggest_specialist(
        "turn the csv into a word table")
    assert name == "", (name, reason)
    assert "unsure" in reason and ranked, reason
    # Both rivals are named, so the reader can adjudicate by hand.
    for cand, _score in ranked[:2]:
        assert cand in reason


def test_three_cue_words_in_common_is_decided_and_says_which_three():
    """The other side of the margin: not every request is a tie.

    This is the paired case that keeps the refusal above honest — a scorer
    that refused EVERYTHING would pass it. The base-green controls in this
    file are the four built-in commands and the explicit-name path below;
    nothing that calls ``suggest_specialist`` can be green on the base,
    since the base has no such function.

    The input is level-of-theory wording, three words only method-reviewer
    claims — deliberately not "a read-only agent", which an earlier
    revision of this test pinned to ``explore``: "read-only" sits in 9 of
    the 18 shipped cues, so a request that says nothing but the shared
    prefix of nine descriptions genuinely has no owner, and a test that
    demanded a pick there was pinning noise.
    """
    name, reason, _ranked = subagents.suggest_specialist(
        "Compare these folders for the same level of theory.")
    assert name == "method-reviewer", (name, reason)
    for word in ("compar", "level", "theory"):
        assert word in reason, reason


def test_nothing_in_the_registry_is_claimed_and_says_so():
    name, reason, ranked = subagents.suggest_specialist("qwertyuiop zxcv")
    assert name == "" and reason == "" and ranked == []


def test_a_mistyped_name_finds_the_name_it_nears_and_an_unrelated_one_nothing():
    assert subagents.nearest_specialist("code-rewiewer") == "code-reviewer"
    assert subagents.nearest_specialist("ts-find") == "ts-finder"
    # Nothing near: no invented name.
    assert subagents.nearest_specialist("zzz-qwerty-nonexistent") == ""
    # A name that IS registered gets no "did you mean" — there is nothing
    # to mean. A version that returned the name itself would make every
    # refusal claim a typo.
    assert subagents.nearest_specialist("plan") == ""


def test_a_task_word_that_merely_sits_near_a_name_is_not_a_typo():
    """The floor that keeps ``/delegate planning the migration`` a PLAN
    request instead of silently rewriting the user's first word into a
    preset name and dropping it from the briefing.

    Measured: ``planning`` vs ``plan`` is 0.667 and ``verify`` vs
    ``verifier`` 0.714, while a real typo scores 0.875-0.923. A refusal
    may guess at 0.6 (it is refused anyway and shows the list); a
    delegation may not.
    """
    assert subagents.nearest_specialist("planning", cutoff=0.85) == ""
    assert subagents.nearest_specialist("verify", cutoff=0.85) == ""
    assert subagents.nearest_specialist("plann", cutoff=0.85) == "plan"


# --- the entry point: what the product does with the request -------------

@pytest.mark.parametrize("command,expected", [
    ("explore", "explore"),
    ("review", "code-reviewer"),
    ("plan", "plan"),
])
def test_a_named_command_still_routes_where_it_did(command, expected):
    """BASE-GREEN CONTROL — passes before and after. The four built-in
    levers are the user naming the specialist, and no fallback may
    second-guess a name the user typed correctly."""
    out = builtin_subagent_command(command, "the usual task text here")
    assert _type_of(out) == expected


def test_an_explicit_preset_name_is_never_overruled():
    """BASE-GREEN CONTROL — passes before and after."""
    out = builtin_subagent_command(
        "delegate", "explore where the CONTROL file is parsed")
    assert _type_of(out) == "explore"
    assert "where the CONTROL file is parsed" in out


def test_a_bare_delegate_names_the_specialist_it_chooses_and_says_why():
    """The defect this package is named for.

    On 41fbf18e2 this line answered ``general-purpose`` with the gloss
    "self-contained task, inherits permissions" for a request whose cue
    said something else. Now the entry point reads the request against the
    shipped cues, and if it picks, it must SAY it picked.
    """
    out = builtin_subagent_command(
        "delegate",
        "Vergleich die Ordner und sag mir, ob die Ketten dieselbe "
        "Level of Theory teilen.")
    assert _type_of(out) == "method-reviewer", out[:200]
    assert "suggest" in out.lower(), (
        "a specialist chosen by the product must be labelled as a "
        "suggestion, not stated as the user's choice")
    # The reason names cue words, so a reader can disagree with it.
    assert "matched" in out


def test_a_mistyped_specialist_is_corrected_in_the_directive():
    """"/delegate code-rewiewer <task>" names a specialist the user
    believes exists. The directive must route to the name it nears, say it
    corrected a name, and not hand the delegate the typo as if it were
    part of the work."""
    out = builtin_subagent_command(
        "delegate", "code-rewiewer review the diff in delfin/agent/config.py")
    assert _type_of(out) == "code-reviewer", out[:200]
    assert "not registered" in out
    assert "code-rewiewer" not in out.split("Task:")[-1]


def test_an_unsure_delegate_says_no_specialist_fit_rather_than_picking_one():
    """FALSIFIED by a fallback that defaults to general-purpose silently:
    the line must carry WHY, which is the only record a wrong pick leaves."""
    out = builtin_subagent_command("delegate", "Mach das bitte schnell.")
    assert _type_of(out) == "general-purpose", out[:200]
    assert "no specialist" in out.lower() or "unsure" in out.lower(), out[:250]


# --- the invariant under the whole labelled set -------------------------

FIXTURE = "specialist_routing_w43s1.yaml"


def _rows() -> list[dict]:
    import yaml
    from pathlib import Path
    text = (Path(__file__).parent / "fixtures" / FIXTURE).read_text(
        encoding="utf-8")
    return yaml.safe_load(text)["requests"]


def test_the_labelled_set_is_well_formed_and_labels_registered_roles():
    """The measurement's own guard: 30 distinct rows, each labelled with a
    preset the registry actually has or the word 'none'. A typo in a label
    would silently make a correct row count as wrong."""
    rows = _rows()
    assert len(rows) == 30, len(rows)
    assert len({r["id"] for r in rows}) == 30
    assert {r["lang"] for r in rows} == {"de", "en"}
    assert len({r["domain"] for r in rows}) >= 4
    for r in rows:
        assert r["expect"] == "none" or r["expect"] in subagents.SUBAGENT_PRESETS, r
        for alt in r.get("also_ok") or []:
            assert alt in subagents.SUBAGENT_PRESETS, r


def test_the_entry_point_never_names_a_specialist_that_is_not_registered():
    """Over EVERY labelled row, whatever the entry point emits must be a
    name the dispatch accepts.

    This is the rule a name in a directive exists to satisfy: the model
    reads it and calls the tool, and an unregistered name comes back
    "unknown subagent_type" after the delegation was already planned.
    FALSIFIED by any suggestion built from free text rather than from the
    registry — which is what an unguarded "did you mean" would do with a
    token near nothing.
    """
    for row in _rows():
        out = builtin_subagent_command("delegate", row["text"])
        assert out is not None, row["id"]
        got = _type_of(out)
        assert got in subagents.SUBAGENT_PRESETS, (row["id"], got)


def test_the_scorer_never_picks_a_specialist_its_own_cue_cannot_show():
    """A pick must be explainable by the text it claims to have read.

    FALSIFIED by scoring on anything the cue does not contain (an
    embedding, a hand list of example prompts): the reason string would
    then name a specialist without naming a word that justifies it, and
    the 64-character cue the model itself reads would no longer be the
    single source the routing and the prompt share.
    """
    for row in _rows():
        name, reason, _ranked = subagents.suggest_specialist(row["text"])
        if not name:
            continue
        shared = (set(subagents._route_tokens(row["text"]))
                  & subagents._cue_tokens(
                      name, subagents.SUBAGENT_PRESETS[name].description))
        assert shared, (row["id"], name, reason)
        assert any(word in reason for word in shared), (row["id"], reason)


def test_the_catch_all_row_states_when_it_applies_not_only_what_it_can_do():
    """The one description change this package makes, and why.

    ``general-purpose`` opened with "Full-tool sub-agent that can read,
    edit, run bash", and ``decision_table`` renders only the first clause,
    so what the chooser read was three of the most common verbs a user
    writes while the sentence restricting it ("most tasks fit a narrower
    preset") was cut at 64 characters. On the labelled set that made the
    catch-all outscore two specialists with narrower cues. The criterion
    must be in the part that survives the cut.
    """
    cue = subagents._cue(
        subagents.SUBAGENT_PRESETS["general-purpose"].description)
    assert len(cue) <= subagents._CUE_LIMIT + 2
    assert "narrower" in cue, cue
    assert "Use" in cue, cue
