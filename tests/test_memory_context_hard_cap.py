"""The injected memory block is under a HARD cap, with an overflow line.

Package 7, phase 2. The recall path had a soft ``max_chars`` (6000 by
default) that was practically never the binding constraint: the index
share and the entry count were. A degenerate store — one huge body, an
inflated MEMORY.md — could still spend thousands of characters. These
tests pin the hard cap:

- the WHOLE emitted block (preamble included) never exceeds the hard
  budget, whatever the store looks like;
- when entries were held back to fit, the block says so — N further
  memories, a memory_tidy proposal is available — instead of silently
  dropping them (a silent drop tells neither the model nor the user
  that the store outgrew the budget);
- the ordering is deterministic: same store + same task text, same
  bytes, so the prompt prefix stays cache-stable;
- the budget is a setting: ``agent.memory_context_budget`` lowers it,
  0 or negative restores the uncapped behaviour of the callers' own
  ``max_chars``.

Tests never touch the real ``~/.delfin`` (HOME redirected as in the
existing recall tests).
"""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


def _seed(tmp_path, monkeypatch, repo_root, texts):
    from delfin.agent import memory_store as ms
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    for t in texts:
        ms.save_typed_memory(t, repo_root=repo_root)


def _recall(agent_tree, **kw):
    from delfin.agent.prompt_loader import PromptLoader
    return PromptLoader(agent_tree)._load_external_memory_context(**kw)


HUGE = ("feedback: " + ("use wedge bonds for the stereo centre of every "
        "ligand this long fact goes on about. " * 40))


def test_hard_cap_holds_even_for_degenerate_stores(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent import prompt_loader as pl
    _seed(tmp_path, monkeypatch, agent_tree,
          [HUGE + " variant {}".format(i) for i in range(30)])
    out = _recall(agent_tree, task_text="draw the molecule")
    assert out, "a seeded store must inject something"
    cap = pl.MEMORY_CONTEXT_HARD_CAP
    assert len(out) <= cap, (
        "block of {} chars exceeds the hard cap of {}".format(len(out), cap))


def test_overflow_is_announced_not_dropped_silently(
        agent_tree, tmp_path, monkeypatch):
    _seed(tmp_path, monkeypatch, agent_tree,
          [HUGE + " variant {}".format(i) for i in range(30)])
    out = _recall(agent_tree, task_text="draw the molecule")
    assert "further memories" in out, (
        "held-back entries must be counted in the block, not dropped "
        "silently: {}".format(out[-200:]))
    assert "memory_tidy" in out, (
        "the block must point at the tidy proposal, which is the sanctioned "
        "way to shrink the store")


def test_no_overflow_line_when_everything_fits(
        agent_tree, tmp_path, monkeypatch):
    _seed(tmp_path, monkeypatch, agent_tree,
          ["project: short fact {}".format(i) for i in range(3)])
    out = _recall(agent_tree, task_text="run the test suite")
    assert "further memories" not in out


def test_ordering_is_deterministic_for_identical_inputs(
        agent_tree, tmp_path, monkeypatch):
    _seed(tmp_path, monkeypatch, agent_tree,
          [HUGE + " variant {}".format(i) for i in range(12)])
    a = _recall(agent_tree, task_text="draw the molecule")
    b = _recall(agent_tree, task_text="draw the molecule")
    assert a == b


def test_budget_is_a_setting_and_zero_restores_caller_budget(
        agent_tree, tmp_path, monkeypatch):
    _seed(tmp_path, monkeypatch, agent_tree,
          [HUGE + " variant {}".format(i) for i in range(30)])
    out = _recall(agent_tree, task_text="draw the molecule", max_chars=3000)
    assert "further memories" in out
    assert len(out) <= 3000
