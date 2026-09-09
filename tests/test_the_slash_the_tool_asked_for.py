"""The skill tool refused the spelling its own description asks for.

    skill: "Invoke a skill … Use on '/skill-name' or when one fits."
    name:  "Skill name, no slash."

Two sentences, one tool, and they disagree. The user types `/casscf-setup`
and the model passes it through, which is what the first sentence tells it
to do; the executor then answers "skill '/casscf-setup' not found" and
lists ten names that all look like the one it just refused.

Found by driving the tool, 2026-09-09. It is recoverable -- the error
lists what is available -- but it costs a turn, and it is the framework
contradicting itself, which is the class of defect that has produced most
of this suite's avoidable failures.

A leading slash cannot be part of a real skill name: skills are files
under `.delfin/skills/`. Nor can a `.md` suffix, which is what a model
copies out of an `ls` of that directory. Both are stripped, so the tool
accepts every spelling of the thing it already knows the user meant.
"""

from __future__ import annotations

import json
import tempfile

import pytest

import delfin.agent.api_client as A


@pytest.fixture
def perms():
    with tempfile.TemporaryDirectory() as tmp:
        yield A.KitToolPermissions(mode="default", workspace=tmp)


def _skill(name, perms):
    return json.loads(A._doc_executor.execute("skill", {"name": name}, perms))


@pytest.mark.parametrize("spelling", [
    "freq-thermochemistry",
    "/freq-thermochemistry",
    "freq-thermochemistry.md",
    "/freq-thermochemistry.md",
    "  /freq-thermochemistry  ",
])
def test_every_spelling_of_a_real_skill_resolves(spelling, perms):
    out = _skill(spelling, perms)
    assert out.get("skill") == "freq-thermochemistry", out


def test_the_description_and_the_tool_now_agree():
    """The sentence that started this. If it is ever reworded away, the
    stripping below is unexplained; if the stripping goes, the sentence
    is a lie again."""
    entry = next(t for t in A._DOC_TOOLS_OPENAI
                 if t["function"]["name"] == "skill")
    assert "/skill-name" in entry["function"]["description"]


@pytest.mark.parametrize("spelling", ["no-such-skill", "/no-such-skill"])
def test_an_unknown_name_still_says_so_and_lists_the_real_ones(
        spelling, perms):
    out = _skill(spelling, perms)
    assert "not found" in out["error"]
    assert "freq-thermochemistry" in out["available"]


def test_a_name_that_is_only_punctuation_is_refused(perms):
    """Stripping must not turn a nonsense argument into an empty lookup."""
    assert "non-empty" in _skill("/", perms)["error"]
    assert "non-empty" in _skill("///", perms)["error"]


def test_stripping_does_not_reach_past_the_leading_slash(perms):
    """`//x` is the name `x`, not a path walk, and a slash in the MIDDLE
    is left alone so a nested name cannot be forged into a parent one."""
    assert "not found" in _skill("//x", perms)["error"]
    out = _skill("../../../etc/passwd", perms)
    assert "error" in out
