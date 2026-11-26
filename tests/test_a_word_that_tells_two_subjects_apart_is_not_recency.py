"""A word that is the DIFFERENCE between two subjects is not a recency marker.

Measured on this branch before the fix, one two-note store written through
``save_typed_memory`` with distinct ``created_at`` stamps:

    task='Which partition do CPU jobs go to?'  GPU fact=False  CPU fact=True
    task='Which partition do GPU jobs go to?'  GPU fact=False  CPU fact=True

The GPU answer is withheld on every task, including the GPU question: the two
titles differ only in ``gpu``/``cpu``, the subject signature counts stopwords,
and Jaccard 0.714 crosses the 0.6 bar, so the note that answers the task is
collapsed away by the note that does not. Base injected both.

A recency word is the opposite case and must still collapse: the bundle's own
must-merge pair is "The old environment of this project." against "The
environment of this project today.", whose differing words (``old``/``today``)
say WHEN a fact was learned, never WHICH fact it is. That pair is the reason a
blanket "each title carries a word the other lacks" guard cannot be the fix --
it would break the collapse the corpus needs. So the two populations are
separated by the CLASS of the differing word, not by its presence.
"""

from __future__ import annotations

from pathlib import Path

import pytest

GPU = ("Queue policy for the gpu partition",
       "Queue policy: GPU jobs go to the ZZZ-GPU-ONE queue.")
CPU = ("Queue policy for the cpu partition",
       "Queue policy: CPU jobs go to the ZZZ-CPU-TWO queue.")
# The bundle's must-merge pair (prompt_loader._SUBJECT_SIMILARITY's note).
OLDER_ENV = ("The old environment of this project.",
             "Dependencies went into .venv-old (ZZZ-OLD).")
NEWER_ENV = ("The environment of this project today.",
             "Everything now lives in .venv-new (ZZZ-NEW).")


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


def _store(root, pairs):
    from delfin.agent import memory_store as ms
    for title, body in pairs:
        ms.save_typed_memory(f"{title}\n{body}", repo_root=root, title=title,
                             memory_type="project", source="user")


def _age_them(root, markers_oldest_first):
    """Stamp the write order onto the frontmatter (index 0 = oldest).

    Two writes in the same second carry the same ``created_at``, and the
    collapse deliberately shows BOTH notes when it cannot order them -- so
    without this the store never reaches the branch under test.
    """
    import re
    from delfin.agent import memory_store as ms
    base = ms._delfin_memory_dir(root)
    for i, marker in enumerate(markers_oldest_first):
        hit = [f for f in base.glob("*.md")
               if f.name != "MEMORY.md"
               and marker in f.read_text(encoding="utf-8")]
        assert len(hit) == 1, f"{marker}: expected one note, got {len(hit)}"
        f = hit[0]
        stamp = 1_700_000_000 + 3600 * i
        text = f.read_text(encoding="utf-8")
        text = re.sub(r"created_at: \d+", f"created_at: {stamp}", text)
        text = re.sub(r"updated_at: \d+", f"updated_at: {stamp}", text)
        f.write_text(text, encoding="utf-8")


def _prompt(root, task: str) -> str:
    from delfin.agent.prompt_loader import PromptLoader
    return PromptLoader(root).build_system_prompt(
        role_id="solo_agent", mode_id="solo", task_text=task)


@pytest.mark.parametrize("task", [
    "Which partition do GPU jobs go to?",
    "Which partition do CPU jobs go to?",
])
def test_the_word_that_distinguishes_two_subjects_is_not_treated_as_recency(
        agent_tree, tmp_path, monkeypatch, task):
    """Red before the fix: a word that NAMES the subject is read as recency.

    Both partitions are distinct subjects, so both facts belong in the prompt --
    that is what base does and what the store holds. The failure is on the
    subject that is NOT the newest: its title's distinguishing word is swallowed
    by the collapse, so its fact is withheld even from the task that asks for it.
    """
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _store(agent_tree, [GPU, CPU])
    # Distinct stamps: with a shared created_at the collapse shows BOTH notes by
    # design and returns early, so the branch under test is never reached.
    _age_them(agent_tree, ["ZZZ-GPU-ONE", "ZZZ-CPU-TWO"])

    out = _prompt(agent_tree, task)

    assert "ZZZ-GPU-ONE" in out and "ZZZ-CPU-TWO" in out, (
        "two notes naming gpu and cpu are two subjects, not one restated: both "
        "facts must be reachable, as they are on the base")


def test_a_recency_word_still_collapses_to_one_subject(agent_tree, tmp_path,
                                                       monkeypatch):
    """Control, green before AND after: ``old``/``today`` is not a subject word.

    This is the pair the 0.6 bar was measured to admit, and the reason the fix
    may not key on "each title has a word the other lacks" alone.
    """
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _store(agent_tree, [OLDER_ENV, NEWER_ENV])
    _age_them(agent_tree, ["ZZZ-OLD", "ZZZ-NEW"])

    out = _prompt(agent_tree, "Which environment does this project use?")

    assert "ZZZ-NEW" in out, "the current note must reach the prompt"
    assert "ZZZ-OLD" not in out, (
        "a note differing only in a recency word is the same subject and must "
        "still be collapsed")
