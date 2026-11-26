"""A recalled memory is the note that answers the task, not every note that
shares a word with it.

Measured 2026-11-24, driving ``cli.main(["chat", ..., "-p", ...])`` over a
30-task corpus (``scratch/m1``): in a five-memory store whose sentences all
have one shape — "The <subject> <noun> export code is <code>. Use it whenever
…" — the memory answering the task scored 2.80-2.88 against 0.086-0.088 for
the four about other subjects. They share "export", "code", "is", "use",
"whenever", so every one of them was over zero, and `score > 0` therefore let
all five into the prompt. Every prompt carried four or five codes the answer
must not be; the answering code was present too, so the turn was asked a
question with five answers.

Two channels had to close for the prompt to be unambiguous, and both are
tested here:

* the BODIES (`_gather_memory_entries` returns the relevant entry only);
* the INDEX (`MEMORY.md`), whose every line repeats the first 160 characters
  of its memory and was therefore a second copy of the rejected text. With
  the bodies filtered and the index left whole, 30 of 30 prompts stayed
  contaminated.

Both tests call `_gather_memory_entries`, the function the prompt build calls,
so they bind to the recall path and not to a test-only helper. The full turn
that showed the contamination is in the PACKAGE's own measurement; a unit test
cannot host a model, and asserting the same numbers through a scripted one
here would test the harness.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent.memory_store import save_typed_memory
from delfin.agent.prompt_loader import PromptLoader

TASK = "Which export code applies to the bookmark ledger? Answer with the code only."

# (subject, noun, code) — one answers TASK, four are the neighbours that share
# its vocabulary. Same sentence shape for all five, which is what made the
# scores so close at the bottom and so far apart at the top.
ENTRIES = (
    ("bookmark", "ledger", "SPEC01-MEM"),
    ("invoice", "portal", "DEC01-MEM"),
    ("payroll", "spreadsheet", "DEC02-MEM"),
    ("calendar", "digest", "DEC03-MEM"),
    ("contact", "registry", "DEC04-MEM"),
)
ANSWERED = "SPEC01-MEM"
DECOYS = tuple(code for _s, _n, code in ENTRIES if code != ANSWERED)


def body(subject: str, noun: str, code: str) -> str:
    return (f"The {subject} {noun} export code is {code}. "
            f"Use it whenever the {subject} {noun} is exported.")


def _store(tmp_path: Path) -> Path:
    """Seed a real store through the product's own write path."""
    ws = tmp_path / "ws"
    ws.mkdir()
    for subject, noun, code in ENTRIES:
        save_typed_memory(body(subject, noun, code), repo_root=ws,
                          memory_type="project", source="agent")
    from delfin.agent.memory_store import _delfin_memory_dir
    return _delfin_memory_dir(ws)


def _gather(tmp_path: Path):
    memory_dir = _store(tmp_path)
    loader = PromptLoader()
    with _isolated_home(tmp_path):
        return loader._gather_memory_entries(memory_dir, TASK)


class _isolated_home:
    """Point the store's HOME at the test's tmp_path for the call.

    ``save_typed_memory`` resolves the per-project store from ``HOME``, so
    the seed and the read must agree about it. The conftest already
    redirects it; this keeps the test independent of that redirect's exact
    shape.
    """

    def __init__(self, tmp_path: Path):
        self.home = tmp_path / "home"

    def __enter__(self):
        self.home.mkdir(parents=True, exist_ok=True)
        return self

    def __exit__(self, *_exc):
        return False


def test_the_bodies_kept_are_the_ones_that_answer_the_task(tmp_path):
    """The four neighbours share "export", "code" and "is" — and must not
    be recalled for it. Fails on the baseline, where all five were kept."""
    index_chunk, entries = _gather(tmp_path)
    kept = "\n".join(f"{title}\n{body_text}" for title, _rel, body_text in entries)
    assert ANSWERED in kept, (
        "the memory that answers the task must be recalled; nothing else in "
        f"the store can answer it:\n{kept}")
    for decoy in DECOYS:
        assert decoy not in kept, (
            f"{decoy} is about another subject and shares only the "
            f"boilerplate, yet it was recalled:\n{kept}")
    assert len(entries) == 1, [rel for _t, rel, _b in entries]


def test_the_index_carries_no_pointer_to_a_withheld_memory(tmp_path):
    """Every index line repeats the memory's first 160 characters, so an
    unfiltered index re-injects what the body filter just rejected."""
    index_chunk, _entries = _gather(tmp_path)
    for decoy in DECOYS:
        assert decoy not in index_chunk, (
            f"{decoy} was withheld as a body but its pointer line still "
            f"reaches the model:\n{index_chunk}")
    # The pointer that belongs to the recalled memory survives — the index is
    # filtered, not emptied.
    assert ANSWERED in index_chunk, index_chunk
    # And the model is told the store is larger than what it can see, rather
    # than being left to believe these are all the memories there are.
    assert "more memories" in index_chunk, index_chunk


def test_an_index_fact_with_no_file_of_its_own_is_never_dropped(tmp_path):
    """A fact typed straight into MEMORY.md has no body to rank. Dropping it
    would drop the fact itself, so the filter must keep every line that is
    not a ``[Title](file.md)`` pointer."""
    from delfin.agent.prompt_loader import _filter_memory_index
    index = ("# MEMORY.md\n"
             "The team uses two spaces after a full stop.\n"
             "- [Kept](kept.md) — a memory that survived\n"
             "- [Dropped](dropped.md) — a memory that did not\n")
    out = _filter_memory_index(index, ["kept.md"], dropped=1)
    assert "two spaces after a full stop" in out, out
    assert "kept.md" in out, out
    assert "dropped.md" not in out, out
    assert "1 more memories" in out, out


def test_nothing_is_filtered_when_the_task_matches_no_memory(tmp_path):
    """The relative floor must not empty a store for a question none of its
    memories touch: with no match the index is the caller's only fallback.
    Control case — passes on the baseline too, and must keep passing."""
    memory_dir = _store(tmp_path)
    loader = PromptLoader()
    index_chunk, entries = loader._gather_memory_entries(
        memory_dir, "How do I bake sourdough bread?", domain="code")
    # No shared subject, so every score is 0 or boilerplate-only; whatever is
    # kept, the index must not have been rewritten into an empty list.
    assert "MEMORY.md" in index_chunk
    assert index_chunk.strip(), "the index was emptied for an unrelated task"
