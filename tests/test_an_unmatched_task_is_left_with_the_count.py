"""What a task NOTHING in the store matches is left with: the count, not silence.

Measured on this branch, one note written through ``save_typed_memory`` and the
task "Sort this list of numbers in ascending order":

    title in prompt: False   value in prompt: False   count line: True

i.e. the prompt carries "... and 1 more memories, not listed here — ask about a
subject to have its memory recalled" and nothing else from the store. Two
comments in ``_gather_memory_entries`` said otherwise -- that the pointer list
"still goes in" there, so the model learns which SUBJECT exists, and that the
index "stays as written" when nothing matched. Neither describes the code: the
pointer gate runs whenever the task has content words, and the reason it must
is that a pointer line repeats the note's first 160 characters, which on a
store the product writes IS the note's value.

The proposed alternative -- keep the pointer with its hook removed ("title +
filename only") -- is not available either: the stored title IS that first 160
characters, so a title is not a value-free label. Measured on the corpus store
(``tests/test_a_recalled_memory_is_the_one_that_answers_the_task.py``):

    - [The invoice portal export code is DEC01-MEM. Use it whenever …
       project_the-invoice-portal-export-code-is-dec01-mem-use-it-whenever.md]

Keeping that line would hand back the decoy code the body gate had just
rejected -- which is why the shipped pin at
``tests/test_a_task_that_matches_no_memory_injects_no_bodies.py`` forbids the
title. This file pins the consequence so the comments cannot drift from it
again.
"""

from __future__ import annotations

from pathlib import Path

import pytest

TITLE = "The partition jobs are submitted to."
BODY = "Queue shortq; the GPU one is booked by another group."


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


def _prompt(root, task: str) -> str:
    from delfin.agent.prompt_loader import PromptLoader
    return PromptLoader(root).build_system_prompt(
        role_id="solo_agent", mode_id="solo", task_text=task)


def test_an_unmatched_task_is_left_with_the_count_and_no_value(
        agent_tree, tmp_path, monkeypatch):
    """best == 0 still gates the pointers: the store says it is larger, nothing more."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _store(agent_tree, [(TITLE, BODY)])

    out = _prompt(agent_tree, "Sort this list of numbers in ascending order.")

    assert "more memories" in out, (
        "the degraded answer is the count: the model must still be told the "
        "store holds more than it can see")
    assert "shortq" not in out, "the fact must not travel"
    assert TITLE not in out, (
        "and neither may the pointer: on a store the product writes, the "
        "title repeats the note's first 160 characters")


def test_the_stored_title_is_the_notes_own_opening_text(agent_tree, tmp_path,
                                                        monkeypatch):
    """Why a title is not a value-free label -- the measurement, pinned.

    This is the fact that rules out "keep the pointer, drop the hook": the hook
    and the title carry the same text, so stripping one does not withhold the
    value.
    """
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    from delfin.agent import memory_store as ms
    _store(agent_tree, [(TITLE, BODY)])

    index = (ms._delfin_memory_dir(agent_tree) / "MEMORY.md").read_text(
        encoding="utf-8")

    assert TITLE in index, "the pointer names the note"
    hook = [ln for ln in index.splitlines() if "submitted-to.md" in ln]
    assert hook, "the store writes a pointer line for this note"
    assert TITLE in hook[0], (
        "the pointer's title is the note's own opening text, so it is not a "
        "value-free label")
