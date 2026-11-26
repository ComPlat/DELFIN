"""A question nothing in the store answers injected six unrelated notes.

Measured on a 30-task corpus driven through ``cli.main`` (wave 43 package m1):
of 10 tasks that NO stored memory applies to, 0 kept every memory body out.
Two mechanisms, both in ``prompt_loader._gather_memory_entries``:

1. When no memory shares a token with the task, the BM25 filter is skipped
   ("dropping everything there would be worse than falling back to MEMORY.md
   order") and the budget is then spent on bodies in MEMORY.md order — six
   unrelated notes for a question about sorting a list.
2. When exactly one INCIDENTAL word matches, the filter runs and keeps the
   memory that word appears in. Measured overlaps that admitted a body:
   ``what`` (twice), ``are`` (twice), ``folder``. ``_tokenize``'s stopword list
   holds 10 words and none of these — an interrogative is not a topic.

The fix is not "recall less": the MEMORY.md pointer list still goes in, so the
model learns a note EXISTS and what its subject is, and says so instead of
guessing. What stops travelling is the FACT, which is the part that can be
applied to the wrong question. ``format_memory_context`` and
``recall_episodes`` already gate on ``score > 0`` without a fallback; this path
was the one that did not.
"""

from __future__ import annotations

from pathlib import Path

import pytest

UNRELATED_TITLE = "The partition jobs are submitted to."
UNRELATED_BODY = "Queue shortq; the GPU one is booked by another group."
DECIDER_TITLE = "How a list is sorted here."
DECIDER_BODY = "Always ascending, nulls last."


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


def test_an_unanswerable_task_sends_bodies_but_no_facts(agent_tree, tmp_path,
                                                        monkeypatch):
    """The failing input, named: a sorting question, a store about queues."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _store(agent_tree, [(UNRELATED_TITLE, UNRELATED_BODY)])

    out = _prompt(agent_tree, "Sort this list of numbers in ascending order.")

    # The degraded mode is "a note exists", not silence -- but it is the COUNT
    # line that says so, not the pointer. Measured 2026-11-24 (see
    # delfin/agent/prompt_loader.py, _filter_memory_index): every ``- [title]
    # (file.md) — hook`` line repeats the note's first 160 characters, so a
    # pointer kept for a memory whose body was just withheld hands the model
    # the very fact the gate rejected. In this store the note's own title
    # carries its value, so "the title survives" and "the value does not
    # travel" cannot both hold. The count line keeps the fact that the store
    # is larger than what the model can see, which is what this assertion is
    # for; the pointer is filtered with the body.
    assert "more memories" in out, (
        "the model must still be told a note it cannot see exists")
    assert UNRELATED_TITLE not in out, (
        "the note's pointer repeats its value; it must be withheld too")
    assert UNRELATED_BODY not in out, (
        "the note answers nothing about sorting; its fact must not travel")


def test_a_task_sharing_one_incidental_word_still_injects_no_facts(
        agent_tree, tmp_path, monkeypatch):
    """``what`` is a question word, not a topic — it must not open the gate."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _store(agent_tree, [(UNRELATED_TITLE,
                         UNRELATED_BODY + " What the queue name means is not "
                                          "documented.")])

    out = _prompt(agent_tree, "What is the capital of France?")

    assert "more memories" in out
    assert UNRELATED_BODY not in out


def test_a_task_a_memory_does_answer_still_gets_its_fact(agent_tree,
                                                         tmp_path,
                                                         monkeypatch):
    """Control, green before and after: relevance must still deliver the body."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    _store(agent_tree, [(UNRELATED_TITLE, UNRELATED_BODY),
                        (DECIDER_TITLE, DECIDER_BODY)])

    out = _prompt(agent_tree, "How should I sort the list of numbers?")

    assert DECIDER_BODY in out, "a matching memory must still be injected"
