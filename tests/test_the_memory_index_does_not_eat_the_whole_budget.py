"""At scale the recall budget held nothing but index lines.

``_load_external_memory_context`` appended the whole MEMORY.md index
unconditionally and BEFORE any cap was tested, and the BM25 relevance
ranking was applied only to the bodies. Each index line carries the first
~160 characters of a memory, so with the shipped prune caps (feedback and
user get 4x the 25-per-type cap, the global store gives every type that
treatment: ~650 files across the two stores) the index alone runs to tens
of kilobytes against a ``max_chars`` of 6000 (solo) or 2000
(dashboard / office).

Measured on this store shape: 60 project memories produce a 9.7 kB index,
which at max_chars=6000 leaves room for ZERO bodies. And because the
recall recorder only runs for injected BODIES, in that state every
``updated_at`` freezes — so the decay signal dies exactly when the store
is largest and needs it most.
"""

from __future__ import annotations

import random
from pathlib import Path

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


_WORDS = """alpha bravo charlie delta echo foxtrot golf hotel india juliet
kilo lima mike november oscar papa quebec romeo sierra tango uniform victor
whiskey xray yankee zulu anchor beacon cipher domino ember falcon granite
harbor ingot jigsaw kernel lantern meridian nucleus obelisk pylon quarry
ridge summit talon umbra vector willow xenon yonder zephyr""".split()


def _fill(agent_tree, count: int, seed: int = 7) -> None:
    from delfin.agent import memory_store as ms
    rng = random.Random(seed)
    for i in range(count):
        ms.save_typed_memory(
            "project: " + " ".join(rng.sample(_WORDS, 18)) + f" case{i}",
            repo_root=agent_tree, title=f"note {i:03d}")


def _bodies_in(out: str) -> int:
    return sum(1 for line in out.splitlines()
               if line.startswith("# note ") and ".md)" in line)


def test_a_large_store_still_injects_bodies(agent_tree, tmp_path, monkeypatch):
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _fill(agent_tree, 60)
    index_chars = len(
        (ms._delfin_memory_dir(agent_tree) / "MEMORY.md").read_text())
    assert index_chars > 6000, "fixture must exceed the solo budget"

    out = PromptLoader(agent_tree)._load_external_memory_context(max_chars=6000)

    assert _bodies_in(out) >= 5
    assert len(out) <= 6100          # cap + the truncation marker


def test_the_index_never_takes_more_than_its_share(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent.prompt_loader import (
        _MEMORY_INDEX_BUDGET_SHARE, PromptLoader,
    )

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _fill(agent_tree, 60)

    out = PromptLoader(agent_tree)._load_external_memory_context(max_chars=6000)

    index = out.split("# MEMORY.md", 1)[1].split("\n\n", 1)[0]
    assert len(index) <= 6000 * _MEMORY_INDEX_BUDGET_SHARE + 200


def test_a_truncated_index_says_how_many_it_did_not_list(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _fill(agent_tree, 60)

    out = PromptLoader(agent_tree)._load_external_memory_context(max_chars=6000)

    assert "more memories, not listed here" in out


def test_the_index_lines_that_survive_are_the_relevant_ones(
        agent_tree, tmp_path, monkeypatch):
    """The ranking existed and was spent on the bodies only, so which
    pointers the model saw was decided by MEMORY.md insertion order."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _fill(agent_tree, 60)
    ms.save_typed_memory(
        "project: the tightscf slowconv recipe fixes orca geometry runs",
        repo_root=agent_tree, title="orca recipe")

    out = PromptLoader(agent_tree)._load_external_memory_context(
        max_chars=2000, task_text="orca geometry convergence keeps failing")

    assert "tightscf" in out


def test_recall_still_marks_memories_used_when_the_store_is_large(
        agent_tree, tmp_path, monkeypatch):
    """The decay signal died at scale: no body injected ⇒ no recall
    recorded ⇒ every updated_at frozen ⇒ prune and expiry read a store in
    which nothing has been used since it was written."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _fill(agent_tree, 60)

    PromptLoader(agent_tree)._load_external_memory_context(max_chars=6000)

    bumped = [r for r in ms.list_typed_memories(agent_tree)
              if r["use_count"] > 1]
    assert len(bumped) >= 5


def test_a_small_store_is_injected_exactly_as_before(
        agent_tree, tmp_path, monkeypatch):
    """Under the allowance the index is passed through untouched — the cap
    is for the pathological store, not a rewrite of the normal one."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory("project: the build runs on python 3.13",
                         repo_root=agent_tree)
    index = (ms._delfin_memory_dir(agent_tree) / "MEMORY.md").read_text().strip()

    out = PromptLoader(agent_tree)._load_external_memory_context()

    assert f"# MEMORY.md\n{index}" in out
    assert "not listed here" not in out


def test_a_fact_typed_straight_into_the_index_is_not_dropped(
        agent_tree, tmp_path, monkeypatch):
    """A line the user wrote into MEMORY.md by hand has no file of its own,
    so trimming it away would delete the fact rather than a pointer."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _fill(agent_tree, 60)
    index_path = ms._delfin_memory_dir(agent_tree) / "MEMORY.md"
    index_path.write_text(
        "# DELFIN Project Memory\n"
        "- the cluster head node is uc3n990, always use it\n"
        + index_path.read_text().split("\n", 1)[1],
        encoding="utf-8")

    out = PromptLoader(agent_tree)._load_external_memory_context(max_chars=6000)

    assert "uc3n990" in out
