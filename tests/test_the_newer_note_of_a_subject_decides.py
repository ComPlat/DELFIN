"""Two notes on one subject, both injected, nothing saying which is newer.

Measured on the wave-43 m1 corpus (30 tasks driven through ``cli.main``): of 4
contradiction tasks, the superseded fact stayed out of the prompt **0** times.
The store's only supersede path is the in-place merge at body-Jaccard >= 0.72,
and a correction that changes ONE value does not reach it — measured 0.333 for
"The partition is shortq" vs "Jobs go to the partition named longq", because
the bodies differ exactly in the contradicted value. The two notes therefore
sit side by side on disk and side by side in the prompt, and the model is left
to guess which the user meant.

Worse where the stale note is the lexical winner. Task "Which virtualenv do I
install the dependencies into?": the CURRENT note ("Everything now lives in
.venv-new") scores 0.00 because the question's word "dependencies" appears in
the OLD note ("Dependencies went into .venv-old"), which scores 2.48. Ranking
alone picks the stale claim, so the subject check has to run BEFORE the
relevance filter, not after it.

Nothing is deleted here: at recall the group contributes the newest body plus
one line naming the older note and where to see it. The older file stays on
disk, and `/memories` and `forget` remain the ways to actually retire it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

OLDER = ("The environment of this project.",
         "Dependencies went into .venv-old (ZZZ-OLD).")
NEWER = ("The environment of this project today.",
         "Everything now lives in .venv-new (ZZZ-NEW).")
UNRELATED = ("The partition jobs are submitted to.",
             "Queue shortq; the GPU one is booked by another group.")


@pytest.fixture
def agent_tree(tmp_path):
    (tmp_path / "pack" / "shared").mkdir(parents=True)
    (tmp_path / "pack" / "agents").mkdir()
    return tmp_path


def _store(root, pairs):
    """Write memories through the product's own path, oldest first."""
    from delfin.agent import memory_store as ms
    for title, body in pairs:
        ms.save_typed_memory(f"{title}\n{body}", repo_root=root, title=title,
                             memory_type="project", source="user")


def _files_with(root, needle: str) -> list[str]:
    from delfin.agent import memory_store as ms
    base = ms._delfin_memory_dir(root)
    return sorted(f.name for f in base.glob("*.md")
                  if f.name != "MEMORY.md"
                  and needle in f.read_text(encoding="utf-8"))


def _age_them(root, markers_oldest_first):
    """Stamp the write order onto the frontmatter (index 0 = oldest).

    Takes the marker that identifies each note and resolves it to the file the
    store actually wrote. Two writes in the same second carry the same
    ``created_at``, so without this nothing on disk says which of two notes on
    one subject is the later one; the alternative, sleeping between writes,
    would bind the test to wall-clock speed.

    BOTH stamps are set, and for different reasons. ``created_at`` is the field
    the recall code reads, because ``record_memory_recall`` bumps ``updated_at``
    for every note it injects — so ``updated_at`` records the last time a note
    was SHOWN, not when the fact was learned, and keying on it would let a stale
    note that has already been recalled outrank the correction that supersedes
    it. ``updated_at`` is stamped too so no other age-based behaviour of the
    store (pruning, decay) fires on a fixture it was never meant to judge.
    """
    import re
    from delfin.agent import memory_store as ms
    base = ms._delfin_memory_dir(root)
    names = []
    for i, marker in enumerate(markers_oldest_first):
        name, = _files_with(root, marker)
        names.append(name)
        f = base / name
        stamp = 1_700_000_000 + 3600 * i
        text = f.read_text(encoding="utf-8")
        text = re.sub(r"created_at: \d+", f"created_at: {stamp}", text)
        text = re.sub(r"updated_at: \d+", f"updated_at: {stamp}", text)
        f.write_text(text, encoding="utf-8")
    return names


def _prompt(root, task: str) -> str:
    from delfin.agent.prompt_loader import PromptLoader
    return PromptLoader(root).build_system_prompt(
        role_id="solo_agent", mode_id="solo", task_text=task)


def test_the_newer_note_of_a_subject_decides_even_when_it_scores_lower(
        agent_tree, tmp_path, monkeypatch):
    """The corpus case, named: the question's word sits in the STALE note.

    Fails on the baseline by injecting ZZZ-OLD and not ZZZ-NEW — the stale
    claim wins the ranking because "dependencies" is in its body.
    """
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    # Written oldest first, as the two sessions that produced them were.
    _store(agent_tree, [OLDER, UNRELATED, NEWER])
    _age_them(agent_tree, ["ZZZ-OLD", "ZZZ-NEW"])

    out = _prompt(agent_tree,
                  "Which virtualenv do I install the dependencies into?")

    assert "ZZZ-NEW" in out, "the current note must reach the prompt"
    assert "ZZZ-OLD" not in out, (
        "the note it replaces must not travel beside it")


def test_a_newer_note_says_what_it_replaced(agent_tree, tmp_path, monkeypatch):
    """The subject must not be silently narrowed: name the note left out.

    Red on the baseline for the same reason as the test above (the old fact is
    injected), but its discriminating clause is the other way round: once the
    older note IS withheld, an assertion on its file name is only meaningful if
    the name survives OUTSIDE the MEMORY.md pointer list — it appears there
    whatever happens, so "the name is in the prompt" could never fail. The
    assertion below therefore strips the pointer list first.
    """
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _store(agent_tree, [OLDER, UNRELATED, NEWER])
    older_file, _newer = _age_them(agent_tree, ["ZZZ-OLD", "ZZZ-NEW"])

    # The corpus's own prompt (s02). A query sharing NO topical word with the
    # note -- "Which virtualenv is current?" -- is answered by the gate with a
    # pointer and no fact, which is the measured behaviour for a task nothing
    # in the store lexically matches (see scratch/PR.md); it is not the case
    # this test is about.
    out = _prompt(agent_tree,
                  "Which virtualenv do I install the dependencies into?")

    assert "ZZZ-OLD" not in out, "the superseded fact must not travel"
    assert "ZZZ-NEW" in out
    below_index = out.split("# MEMORY.md")[-1]
    below_index = below_index.split("\n\n", 1)[-1]  # drop the index lines
    assert older_file in below_index, (
        "withholding a note silently is worse than showing both: the prompt "
        "must name what it left out and where to see it")


def test_two_unrelated_notes_are_both_still_injected(agent_tree, tmp_path,
                                                     monkeypatch):
    """Control, green before and after: different subjects must not collapse."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _store(agent_tree, [UNRELATED,
                        ("The dataset version to use.",
                         "Dataset v3; v2 misses two systems (ZZZ-V3).")])

    out = _prompt(agent_tree, "Sort this list of numbers in ascending order.")

    # Nothing matches the task, so no fact travels -- but the store must still
    # show that BOTH subjects are there and were not treated as one. The COUNT
    # is what proves it: it is computed after the collapse, so a store that had
    # collapsed these two into one subject would report "1 more memories" here.
    # (The pointer lines themselves are withheld with their bodies, because a
    # pointer repeats the note's first 160 characters -- see
    # test_a_task_that_matches_no_memory_injects_no_bodies.py.)
    assert "2 more memories" in out, (
        "two distinct subjects must both survive the collapse to one subject")
    assert UNRELATED[1] not in out
    assert "Dataset v3" not in out
