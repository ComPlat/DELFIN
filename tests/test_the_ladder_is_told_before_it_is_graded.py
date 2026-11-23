"""A ladder the model is graded against and one it is told about are one.

The excited-state rungs were enforcement only: lint rules the grader
turns into pressure after the fact. That is half the mechanism, and the
measured half is the other one -- on 558 published vertical transition
energies an unguided agent reached 24.9 % coverage and the same agent
given the decision sequence UP FRONT reached 44.1 %. Grading alone does
not produce that; it only says afterwards that a decision was missing.

So the rungs are now declared once (``RUNGS``) and read by two
consumers: the lint, which judges a spec, and ``ladder_text``, which the
prompt loader renders into the prompt the model reads before it builds
one.

Why one declaration and not two texts: measured on this repository. The
subagent roles were a registry AND a hand-written table, and six of
thirteen roles ended up callable with nothing in the prompt saying what
they were for. A second copy of an answer drifts, and the drift is
silent.

Input: the rungs, the lint through its single entry point, and the
prompt as the loader BUILDS it. Output: that the three agree.
"""
from __future__ import annotations

import pytest

from delfin.agent.prompt_loader import PromptLoader
from delfin.tools._science_rules import excited
from delfin.tools.platform import scientific_lint


def _built() -> str:
    return PromptLoader().load_role_prompt("solo_agent")


def _rules(spec) -> list[str]:
    out = []
    for f in scientific_lint(spec):
        out.append(f.get("rule") if isinstance(f, dict) else f.rule)
    return out


# --- the two sides are one list ---------------------------------------

def test_every_rung_is_in_the_prompt_the_model_reads():
    """The point of the whole change. A rung the model is graded on and
    never told about is the defect this file exists to prevent."""
    built = _built()
    for rung in excited.RUNGS:
        assert rung.decide in built, f"{rung.code} is not in the prompt"
        assert rung.closes in built, f"{rung.code} does not say what closes it"


#: One spec per rung that MAKES IT FIRE. Written out rather than derived
#: so that adding a rung forces whoever adds it to say what triggers it:
#: a rung with no triggering case is advice, and the test below is the
#: only thing that can tell the two apart.
_TRIGGERS = {
    "excited-character-not-declared": {
        "method": "PBE0", "basis": "def2-TZVP"},
    "ct-needs-long-range-exchange": {
        "method": "PBE0", "basis": "def2-TZVP",
        "excited_state": "charge_transfer", "tda": True},
    "tda-not-stated": {
        "method": "wB97X-D4", "basis": "def2-TZVP",
        "excited_state": "valence"},
    "diffuse-basis-required": {
        "method": "wB97X-D4", "basis": "def2-TZVP",
        "excited_state": "rydberg", "tda": True},
}


def test_every_rung_has_a_case_that_makes_it_fire():
    """The other direction: a rung in the prompt that nothing checks is
    advice, and advice drifts out of date silently.

    Each rung gets the spec that triggers it rather than one shared
    spec -- a shared spec passed while `diffuse-basis-required` could
    not fire on it at all, which would have read as enforcement.
    """
    missing = [r.code for r in excited.RUNGS if r.code not in _TRIGGERS]
    assert not missing, f"no triggering case declared for: {missing}"
    for rung in excited.RUNGS:
        kwargs = {"charge": 0, "mult": 1, **_TRIGGERS[rung.code]}
        spec = {"steps": [{"step": "orca_tddft", "label": "s",
                           "kwargs": kwargs}]}
        assert rung.code in set(_rules(spec)), (
            f"{rung.code} is declared and does not fire on its own "
            "triggering case")


def test_a_rung_added_at_runtime_reaches_the_prompt(monkeypatch):
    """The prompt follows the declaration, not a copy of it. A rung
    appended here appears in the built prompt without anyone editing
    markdown."""
    extra = excited.Rung(
        code="probe-rung",
        decide="state the solvation model for the excited state",
        because="a linear-response model answers a different question",
        closes="solvation_model in the step's kwargs")
    monkeypatch.setattr(excited, "RUNGS", excited.RUNGS + (extra,))
    built = _built()
    assert extra.decide in built
    assert extra.closes in built


def test_the_lint_message_quotes_the_declared_rung(monkeypatch):
    """The finding the grader reads says the same thing as the prompt.
    Changing the declaration changes both, which is the property that
    keeps them from drifting apart."""
    monkeypatch.setattr(
        excited, "RUNGS",
        tuple(excited.Rung(code=r.code, decide=r.decide,
                           because="BECAUSE-MARKER", closes="CLOSES-MARKER")
              for r in excited.RUNGS))
    spec = {"steps": [{"step": "orca_tddft", "label": "s", "kwargs": {
        "charge": 0, "mult": 1, "method": "PBE0", "basis": "def2-TZVP",
        "excited_state": "ct"}}]}
    messages = " ".join(
        str(f.get("message") if isinstance(f, dict) else f.message)
        for f in scientific_lint(spec))
    assert "BECAUSE-MARKER" in messages
    assert "CLOSES-MARKER" in messages


# --- no second copy ----------------------------------------------------

def test_the_pack_file_holds_no_hand_written_ladder():
    """The rungs live in the generated block or nowhere. A hand-written
    copy beside it is how the subagent table drifted."""
    pack = PromptLoader().agent_dir / "agents"
    import re
    generated = re.compile(
        r"^<!--\s*generated:[a-z_]+\s*-->.*?^<!--\s*/generated\s*-->",
        re.MULTILINE | re.DOTALL)
    for path in sorted(pack.glob("*.md")):
        outside = generated.sub("", path.read_text(errors="replace"))
        for rung in excited.RUNGS:
            assert rung.decide not in outside, (
                f"{path.name} states rung {rung.code} outside a generated "
                "block: that is a second copy and it will drift")


def test_the_block_is_named_by_a_provider_that_exists():
    """A generated name nobody owns leaves a hole where a section was;
    the loader keeps the marker instead, which is quiet. This asserts
    the provider is really wired."""
    loader = PromptLoader()
    assert loader._generated_block("excited_ladder")
    assert loader._generated_block("not_a_real_block") is None


# --- the declaration cannot go stale silently -------------------------

def test_a_lint_code_with_no_declared_rung_is_a_loud_error():
    """``_rung`` raises rather than returning a default. A missing rung
    would otherwise produce a finding whose message is empty where the
    reason should be."""
    with pytest.raises(KeyError):
        excited._rung("no-such-rung")


def test_each_rung_says_all_four_things():
    """A rung without a `closes` tells the model a decision is required
    and not how to record it, which is the shape of advice nobody can
    act on."""
    for rung in excited.RUNGS:
        assert rung.code and rung.decide and rung.because and rung.closes
        assert len(rung.decide) >= 12, rung.code
        assert len(rung.closes) >= 12, rung.code


def test_the_ladder_text_is_ordered_and_numbered():
    """The order is the order the decisions bind in: the functional
    governs whether the state can be described at all, the response
    approximation moves the number it returns."""
    text = excited.ladder_text()
    positions = [text.index(r.decide) for r in excited.RUNGS]
    assert positions == sorted(positions)
    assert "1." in text and f"{len(excited.RUNGS)}." in text
