"""Measure what the memory recall actually plays into each turn.

Phase 1 of package 7 (hard memory limit + mid-session nudges): before a
budget is imposed, the current behaviour is measured so the budget is
grounded in numbers rather than taste. These tests pin what the shipped
code does today — they are the MEASUREMENT, and phase 2 turns the size
assertion into the red control for the hard cap.

What is measured (all through the public path,
``PromptLoader._load_external_memory_context``):

- characters injected for stores of 0 / 1 / 20 / 500 facts,
- how much the injected size varies between turns when the task text
  changes (BM25 ranking) versus when it does not (prompt-cache concern:
  a block that changes every turn re-buys the whole suffix),
- determinism: same store + same task text must give the same bytes.

The tests never touch the real ``~/.delfin`` — HOME is redirected the
way the existing recall tests do it.
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


def _recall(loader, **kw):
    return loader._load_external_memory_context(**kw)


def _loader(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader
    return PromptLoader(agent_tree)


def _seed(agent_tree, tmp_path, monkeypatch, n_facts, *, text_fmt=None):
    """Seed a project store with ``n_facts`` typed memories.

    HOME is redirected so nothing reads or writes the real ~/.delfin.
    Returns nothing; the loader finds the store through the same
    ``_delfin_memory_dir(repo_root)`` path production uses.
    """
    from delfin.agent import memory_store as ms
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    text_fmt = text_fmt or (
        "feedback: always run the {} suite with --repeats 3 before "
        "committing branch {}".format)
    for i in range(n_facts):
        ms.save_typed_memory(
            text_fmt(i, "agent/topic-{}".format(i)),
            repo_root=agent_tree)


@pytest.mark.parametrize("n_facts", [0, 1, 20, 500])
def test_injected_characters_for_store_size(
        agent_tree, tmp_path, monkeypatch, n_facts):
    _seed(agent_tree, tmp_path, monkeypatch, n_facts)
    out = _recall(_loader(agent_tree), task_text="run the test suite")
    if n_facts == 0:
        assert out == ""
        return
    print("n_facts={} chars={} tokens~{}".format(
        n_facts, len(out), len(out) // 4))
    # The block exists and is bounded by the soft budget today.
    assert 0 < len(out) <= 6000


def test_size_is_stable_across_identical_turns(
        agent_tree, tmp_path, monkeypatch):
    """Cache-friendliness probe: same store, same task, repeated turns."""
    _seed(agent_tree, tmp_path, monkeypatch, 20)
    loader = _loader(agent_tree)
    a = _recall(loader, task_text="run the test suite")
    b = _recall(loader, task_text="run the test suite")
    assert a == b, "identical turn inputs must give identical bytes"
    print("stable-20-facts chars={}".format(len(a)))


def test_size_swings_between_different_tasks(
        agent_tree, tmp_path, monkeypatch):
    """How far the injected size moves when the task text changes.

    This is the swing a hard budget has to absorb: BM25 re-ranks the
    store per turn, so both the SELECTION and its SIZE change. The
    measurement is recorded, not judged here.
    """
    _seed(agent_tree, tmp_path, monkeypatch, 20)
    loader = _loader(agent_tree)
    # Heterogeneous topics, so BM25 actually has something to rank.
    _seed(agent_tree, tmp_path, monkeypatch, 30, text_fmt=lambda i, slug: (
        "feedback: when drawing molecule {} use wedge bonds for the "
        "stereo centre of ligand {}".format(i, i)))
    sizes = [
        len(_recall(loader, task_text=t))
        for t in ("run the test suite", "draw the molecule",
                  "completely unrelated question about zebras")
    ]
    print("sizes across tasks: {}".format(sizes))
    assert all(s > 0 for s in sizes)


def test_memory_context_block_is_deterministic_per_input(
        agent_tree, tmp_path, monkeypatch):
    """The flat store's own format_memory_context is deterministic too."""
    _seed(agent_tree, tmp_path, monkeypatch, 12)
    from delfin.agent import memory_store as ms
    # ``format_memory_context`` reads the flat legacy JSON store, not the
    # typed store; point it at an explicit file under the redirected home.
    flat = tmp_path / ".delfin" / "agent_memory.json"
    ms._write(flat, {"facts": [
        {"text": "fact {} about the test suite".format(i), "source": "user"}
        for i in range(12)
    ]})
    a = ms.format_memory_context(flat, task_text="run the test suite")
    b = ms.format_memory_context(flat, task_text="run the test suite")
    assert a == b
    print("flat-context chars={}".format(len(a)))
