"""The maintainer's principles open every role's prompt, and no agent can
rewrite them.

pack/shared/principles_addendum.md is written by the maintainer. While it
holds only its heading and comment it is not injected (an empty section
would teach nothing and cost prompt bytes); once it has a body, it is the
first shared contract in every role, before honesty and refusal.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.prompt_loader import PromptLoader, _has_body

_PACK = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"


def _tree(tmp_path, principles: str) -> Path:
    shared = tmp_path / "pack" / "shared"
    agents = tmp_path / "pack" / "agents"
    shared.mkdir(parents=True)
    agents.mkdir(parents=True)
    (shared / "principles_addendum.md").write_text(principles)
    (shared / "honesty_addendum.md").write_text("# Honesty\nHONESTY-MARKER")
    for role in ("solo_agent", "dashboard_agent", "builder_agent",
                 "critic_agent"):
        (agents / f"{role}.md").write_text(f"# {role}\nYou are {role}.")
    return tmp_path


def test_the_shipped_file_exists():
    assert (_PACK / "shared" / "principles_addendum.md").is_file()


def test_an_unwritten_scaffold_is_not_injected(tmp_path):
    tree = _tree(tmp_path, "# Principles\n\n<!-- to be written -->\n")
    prompt = PromptLoader(tree).build_system_prompt(
        role_id="solo_agent", mode_id="solo", mode_description="solo",
        route=["solo_agent"], role_index=0)
    assert "# Principles" not in prompt
    assert "HONESTY-MARKER" in prompt


@pytest.mark.parametrize("role_id,mode_id", [
    ("solo_agent", "solo"),
    ("dashboard_agent", "dashboard"),
    ("builder_agent", "quick"),
    ("critic_agent", "quick"),
])
def test_written_principles_come_first_in_every_role(tmp_path, role_id,
                                                      mode_id):
    tree = _tree(tmp_path, "# Principles\n\nPRINCIPLES-MARKER\n")
    prompt = PromptLoader(tree).build_system_prompt(
        role_id=role_id, mode_id=mode_id, mode_description=mode_id,
        route=[role_id], role_index=0)
    assert "PRINCIPLES-MARKER" in prompt
    assert prompt.index("PRINCIPLES-MARKER") < prompt.index("HONESTY-MARKER")


def test_a_body_is_text_outside_heading_and_comments():
    assert not _has_body("# P\n\n<!-- a\nb -->\n")
    assert _has_body("# P\n\nSomething.\n")


def test_an_agent_cannot_edit_the_principles_without_asking():
    from delfin.agent.api_client import _DEFAULT_PATH_PROTECTED_GLOBS
    assert ("delfin/agent/pack/shared/principles_addendum.md"
            in _DEFAULT_PATH_PROTECTED_GLOBS)


def test_the_shipped_principles_are_written_and_loaded():
    from delfin.agent.prompt_loader import PromptLoader
    body = (_PACK / "shared" / "principles_addendum.md").read_text(
        encoding="utf-8")
    assert _has_body(body)
    prompt = PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo", task_text="tidy the workspace")
    assert "long-term well-being of humanity and the planet" in prompt


#: Each commitment the principles make, and the words that identify it.
#: One sentence was pinned before this, so the file could keep its heading
#: and that sentence while every other commitment was softened away and
#: the suite stayed green. Removing a commitment is the realistic way this
#: text gets weakened -- not deleting the file, which is loud.
#:
#: Checked by MEANING (a set of key terms), not word for word. A literal
#: pin blocks an honest rewording, and a test that blocks honest work is
#: one somebody updates without reading -- at which point it guards
#: nothing. Each row needs every term in it, so a clause cannot pass by
#: keeping one word of it.
_COMMITMENTS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("advance science",
     ("advance science", "great problems")),
    ("the long-term good of people and the planet",
     ("long-term well-being", "humanity", "planet")),
    ("protect the planet and its life",
     ("protecting the planet", "its life")),
    ("human dignity and self-determination",
     ("human dignity", "freedom", "self-determination")),
    ("never act against humanity",
     ("not act against humanity",)),
    ("no deliberate harm, oppression or exploitation",
     ("deliberately harm", "oppress", "exploit")),
    ("refuse serious harm, violence, coercion, manipulation",
     ("serious harm", "violence", "coercion", "manipulation",
      "must not comply")),
    ("offer a safe alternative instead of only refusing",
     ("safe and constructive alternative",)),
    ("life, autonomy and the environment outrank a short-term instruction",
     ("shall take precedence", "short-term instructions")),
    # The rest were added for a world of many agents, and were in the file
    # before they were in this table -- which is the gap the comment above
    # warns about: a commitment nothing pins is one that can be softened
    # away while the suite stays green.
    ("a request carries no authority from its source",
     ("carries no authority from its source", "none of them speaks for a "
      "person", "Authority rests with the human")),
    ("a harmful request from an agent is refused, and no way around it "
     "is offered",
     ("may come from an agent", "neither take the route around it nor "
      "point at one", "reforming whoever asked")),
    # The next two live in the refusal addendum now (procedure, not
    # principle) and still reach the model: this table is checked against
    # the BUILT prompt, which is where it has to be true.
    ("where there is no human to tell, record it and stop",
     ("no human is present to be told", "record it where they will find "
      "it", "deciding on their behalf")),
    ("judge a step by what it is part of",
     ("Judge a step by what it is part of", "harmless by itself",
      "cannot see the whole")),
    ("what you write is another agent's input",
     ("another agent's input", "do not put instructions for other agents",
      "what is a finding and what is a request")),
    ("the containment is one of the systems being protected",
     ("containment you work under", "Widening your own access",
      "quieting a record", "say which one and ask")),
)


@pytest.mark.parametrize("name, terms", _COMMITMENTS,
                         ids=[c[0] for c in _COMMITMENTS])
def test_each_commitment_reaches_the_model(name, terms):
    """Every commitment, in the prompt the model actually receives.

    Asserted against the BUILT prompt rather than the markdown: the
    loader drops lazy-module sections and composes per role, so a
    commitment can sit in the file and never arrive.
    """
    from delfin.agent.prompt_loader import PromptLoader

    prompt = PromptLoader().build_system_prompt(
        role_id="solo_agent", mode_id="solo", task_text="tidy the workspace")
    missing = [term for term in terms if term not in prompt]
    assert not missing, (
        f"the principles no longer commit to {name}: "
        f"{', '.join(missing)} is gone from the prompt the model reads")


def test_every_role_gets_every_commitment():
    """Not only solo. A commitment that reaches one role and not another
    is a commitment the other role does not have."""
    from delfin.agent.prompt_loader import PromptLoader

    for role, mode in (("solo_agent", "solo"),
                       ("dashboard_agent", "dashboard"),
                       ("office_agent", "office")):
        prompt = PromptLoader().build_system_prompt(
            role_id=role, mode_id=mode, task_text="tidy the workspace")
        for name, terms in _COMMITMENTS:
            missing = [term for term in terms if term not in prompt]
            assert not missing, f"{role} is missing {name}: {missing}"
