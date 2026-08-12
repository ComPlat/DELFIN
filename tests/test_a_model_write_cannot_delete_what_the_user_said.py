"""Dedup was negation-blind and let an agent write overwrite a user memory.

The merge fires at Jaccard >= 0.72 over stopword-filtered token SETS. Token
sets do not encode negation: swapping "never" for "always" in a ~20-token
memory leaves almost every token in place and scores far above the bar. The
merge then REPLACES the stored body with the incoming one and rewrites
``source`` to the incoming writer's — so a model-issued ``remember`` call
landing within 0.72 of a user-typed memory overwrote the user's text, flipped
the file to agent-written, and thereby put a user-stated fact under the
90-day agent expiry.

That is a model-reachable path to deleting something the user said, and it
runs through the one channel the model controls end to end.

Two rules close it: an agent write never merges into a user file, and no
merge may invert polarity in either direction.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


_USER_RULE = ("feedback: never push directly to the main branch without a "
              "review from another engineer first")
_INVERTED = ("feedback: always push directly to the main branch without a "
             "review from another engineer first")


def test_the_inverted_rule_scores_above_the_merge_bar(agent_tree):
    """Without this the rest of the file reads as a hypothetical."""
    from delfin.agent import memory_store as ms

    assert ms._jaccard(_USER_RULE, _INVERTED) >= ms._DEFAULT_MERGE_SIMILARITY


def test_an_agent_write_does_not_overwrite_what_the_user_stated(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(_USER_RULE, repo_root=agent_tree,
                         source=ms.SOURCE_USER)
    ms.save_typed_memory(_INVERTED, repo_root=agent_tree,
                         source=ms.SOURCE_AGENT)

    recs = ms.list_typed_memories(agent_tree)
    user = [r for r in recs if r["source"] == ms.SOURCE_USER]
    assert len(user) == 1
    assert "never push directly" in user[0]["body"]


def test_a_user_memory_is_never_demoted_to_agent_grade(
        agent_tree, tmp_path, monkeypatch):
    """The demotion is what subjects a user-stated fact to the agent
    expiry, so it matters even when the wording survives."""
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(
        "feedback: run the benchmark with three repeats and never fewer, "
        "a single run is not evidence of anything",
        repo_root=agent_tree, source=ms.SOURCE_USER)
    ms.save_typed_memory(
        "feedback: run the benchmark with three repeats and never fewer, "
        "a single run is not evidence of any regression",
        repo_root=agent_tree, source=ms.SOURCE_AGENT)

    sources = {r["source"] for r in ms.list_typed_memories(agent_tree)}
    assert ms.SOURCE_USER in sources


def test_a_polarity_flip_is_never_merged_even_between_agent_writes(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(_USER_RULE, repo_root=agent_tree,
                         source=ms.SOURCE_AGENT)
    ms.save_typed_memory(_INVERTED, repo_root=agent_tree,
                         source=ms.SOURCE_AGENT)

    bodies = [r["body"] for r in ms.list_typed_memories(agent_tree)]
    assert len(bodies) == 2
    assert any("never push" in b for b in bodies)
    assert any("always push" in b for b in bodies)


def test_a_german_polarity_flip_is_caught_too(agent_tree, tmp_path,
                                              monkeypatch):
    """The users write German; a negation they type has to count."""
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    negated = ("die Suite niemals ohne drei Wiederholungen laufen lassen, ein "
               "einzelner Durchlauf am Cluster beweist gar nichts ueber eine "
               "Regression im Benchmark")
    affirmed = negated.replace("niemals", "immer")
    assert ms._jaccard(negated, affirmed) >= ms._DEFAULT_MERGE_SIMILARITY

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(f"feedback: {negated}", repo_root=agent_tree,
                         source=ms.SOURCE_AGENT)
    ms.save_typed_memory(f"feedback: {affirmed}", repo_root=agent_tree,
                         source=ms.SOURCE_AGENT)

    assert len(ms.list_typed_memories(agent_tree)) == 2


def test_the_same_rule_reworded_still_merges(agent_tree, tmp_path,
                                             monkeypatch):
    """The guards must not turn dedup off: an agent restating its own
    memory in slightly different words is still one memory."""
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(
        "feedback: never push directly to the main branch without a review "
        "from another engineer first",
        repo_root=agent_tree, source=ms.SOURCE_AGENT)
    ms.save_typed_memory(
        "feedback: never push directly to the main branch without a review "
        "from another maintainer first",
        repo_root=agent_tree, source=ms.SOURCE_AGENT)

    assert len(ms.list_typed_memories(agent_tree)) == 1


def test_the_user_can_still_correct_a_memory_the_model_wrote(
        agent_tree, tmp_path, monkeypatch):
    """Promotion upward stays: the block is on the model writing over the
    user, not on the user writing over the model."""
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(
        "feedback: never push directly to the main branch without a review "
        "from another engineer first",
        repo_root=agent_tree, source=ms.SOURCE_AGENT)
    ms.save_typed_memory(
        "feedback: never push directly to the main branch without a review "
        "from another maintainer first",
        repo_root=agent_tree, source=ms.SOURCE_USER)

    recs = ms.list_typed_memories(agent_tree)
    assert len(recs) == 1
    assert recs[0]["source"] == ms.SOURCE_USER
    assert "maintainer" in recs[0]["body"]
