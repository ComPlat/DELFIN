"""Decay evicted the memory that cited code and kept the vague one.

A memory carrying a ``path:line`` reference gets a "rotted" verdict as soon
as the anchored line moves. The rot path deliberately did NOT bump
``updated_at`` — while ``prune_memories`` sorts survivors by ``updated_at``
and ``_expired`` tests the same field for the 90-day agent cutoff. A memory
with no code reference can never rot, so it was refreshed on every recall.

Since a memory written to CORRECT an earlier one is precisely the one that
cites the code it corrects, the correcting memory sank to the top of the
eviction order while the memory it corrected stayed fresh. The narrowing of
what counts as a path (a method string like ``BP86/def2-TZVP/D3BJ`` is not
a file) removed one class of false verdict and left this mechanism intact.

The two signals are separate now: recall bumps ``updated_at`` for every
injected memory, rotted included, and ``stale_hits`` alone is what makes
the prune prefer evicting a rotted one.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


def _rot_a_memory(agent_tree, tmp_path):
    """Save a memory anchored to a line, then move the line away."""
    from delfin.agent import memory_store as ms
    target = agent_tree / "srcmod.py"
    target.write_text("import os\nMAGIC_LIMIT = 17\n", encoding="utf-8")
    path, _, _ = ms.save_typed_memory(
        "project: the magic limit constant sits at srcmod.py:2 today",
        repo_root=agent_tree, title="magic limit")
    target.write_text("\n".join(f"filler_{i} = {i}" for i in range(80)),
                      encoding="utf-8")
    return path


def test_a_rotted_memory_is_still_marked_as_recalled(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _rot_a_memory(agent_tree, tmp_path)
    before = ms.list_typed_memories(agent_tree)[0]
    time.sleep(1.1)                       # updated_at has one-second grain

    out = PromptLoader(agent_tree)._load_external_memory_context(
        task_text="magic limit")

    assert "[drifted:" in out
    after = ms.list_typed_memories(agent_tree)[0]
    assert after["updated_at"] > before["updated_at"]
    assert after["use_count"] == before["use_count"] + 1
    assert after["stale_hits"] == 1       # rot recorded, separately


def test_a_rotted_memory_does_not_age_out_while_it_is_being_recalled(
        agent_tree, tmp_path, monkeypatch):
    """The agent-memory expiry reads ``updated_at``. With the bump withheld
    a model-written memory that cited code was deleted 90 days after it was
    WRITTEN, however often it had been recalled since."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    path = _rot_a_memory(agent_tree, tmp_path)
    text = path.read_text(encoding="utf-8")
    meta, _ = ms._parse_frontmatter(text)
    old = int(time.time()) - 120 * 86_400
    path.write_text(
        text.replace(f"updated_at: {meta['updated_at']}", f"updated_at: {old}")
            .replace(f"created_at: {meta['created_at']}", f"created_at: {old}")
            .replace("source: user", "source: agent"),
        encoding="utf-8")

    PromptLoader(agent_tree)._load_external_memory_context(
        task_text="magic limit")
    ms.prune_memories(agent_tree)

    assert [r["name"] for r in ms.list_typed_memories(agent_tree)] \
        == ["magic-limit"]


def test_between_two_equally_fresh_memories_the_rotted_one_goes_first(
        agent_tree, tmp_path, monkeypatch):
    """Rot still has to cost something — it just may not cost the field
    that also encodes 'this memory is in active use'."""
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    now = int(time.time())
    memory_dir = ms._delfin_memory_dir(agent_tree)
    memory_dir.mkdir(parents=True, exist_ok=True)
    for name, rot in (("healthy", 0), ("rotted", 3)):
        (memory_dir / f"project_{name}.md").write_text(
            "---\n"
            f"name: {name}\ndescription: d\n"
            f"created_at: {now}\nupdated_at: {now}\nuse_count: 4\n"
            "domain: code\nsource: user\n"
            f"stale_hits: {rot}\n"
            "metadata:\n  type: project\n---\n\nbody\n",
            encoding="utf-8")

    ms.prune_memories(agent_tree, max_per_type=1)

    assert [r["name"] for r in ms.list_typed_memories(agent_tree)] == ["healthy"]
