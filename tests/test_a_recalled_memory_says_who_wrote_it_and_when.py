"""A recalled memory arrived as an instruction, undated and unattributed.

The store goes to real trouble to record WHO wrote each memory — a fact the
user typed and a fact a model invented were byte-identical on disk until
``source:`` was added, and there are tests asserting it lands there. Then
``_memory_body_only`` stripped the entire frontmatter at injection:
created_at, updated_at, use_count, source and superseded never reached the
model. The block header was the whole framing.

A feedback memory is an imperative sentence by construction — the write
prompt asks for exactly that shape. So the block reliably carried lines
like "never run the suite without --repeats 3" into the system prompt, a
few hundred tokens from a section that claims override authority, with
nothing separating them from the role prompt's own rules and nothing
marking them as possibly stale. The only counterweight was one line in a
write-side addendum thousands of tokens away that never named the block.

What these tests pin: the block says what it is, and two memories of
different provenance are distinguishable in the emitted string.
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


def _recall(agent_tree, **kw):
    from delfin.agent.prompt_loader import PromptLoader
    return PromptLoader(agent_tree)._load_external_memory_context(**kw)


def test_the_block_says_it_is_background_and_not_an_instruction(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory("feedback: never run the suite without --repeats 3",
                         repo_root=agent_tree)

    out = _recall(agent_tree)

    lowered = out.lower()
    head = lowered.split("\n\n", 1)[0]
    assert "background" in head
    assert "not instructions" in head
    assert "verify" in head
    # The qualifier has to reach the model BEFORE the imperative does.
    assert lowered.index("--repeats 3") > lowered.index("background")


def test_an_agent_written_and_a_user_stated_memory_are_distinguishable(
        agent_tree, tmp_path, monkeypatch):
    """The missing test. Both memories are feedback; only the provenance
    the store already records tells them apart, and it never got out."""
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(
        "feedback: never run the benchmark without three repeats",
        repo_root=agent_tree, title="repeats rule", source=ms.SOURCE_USER)
    ms.save_typed_memory(
        "feedback: prefer ruff over flake8 for linting this tree",
        repo_root=agent_tree, title="lint choice", source=ms.SOURCE_AGENT)

    out = _recall(agent_tree)

    assert "three repeats" in out and "prefer ruff" in out
    user_line = next(ln for ln in out.splitlines()
                     if ln.startswith("# repeats rule"))
    agent_line = next(ln for ln in out.splitlines()
                      if ln.startswith("# lint choice"))
    assert "the user stated" in user_line
    assert "the agent noted" in agent_line
    assert user_line != agent_line


def test_a_recalled_memory_carries_the_date_it_was_written(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory("project: the release branch is cut on fridays",
                         repo_root=agent_tree, title="release cadence")

    out = _recall(agent_tree)

    assert "# release cadence" in out
    assert time.strftime("%Y-%m-%d") in out


def test_an_old_memory_shows_that_it_is_old(
        agent_tree, tmp_path, monkeypatch):
    """Nothing in the block distinguished a fact written this morning from
    one written two years ago."""
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    path, _, _ = ms.save_typed_memory(
        "project: the api key lives in the vault", repo_root=agent_tree,
        title="key location")
    text = path.read_text(encoding="utf-8")
    old = int(time.time()) - 700 * 86_400
    path.write_text(text.replace(
        f"created_at: {ms._meta_int(ms._parse_frontmatter(text)[0], 'created_at')}",
        f"created_at: {old}"), encoding="utf-8")

    out = _recall(agent_tree)

    assert time.strftime("%Y-%m-%d", time.localtime(old)) in out


def test_the_text_a_merge_overwrote_is_surfaced(
        agent_tree, tmp_path, monkeypatch):
    """Dedup replaces the stored body in place. Without this the fact that
    a memory used to say something else exists only on disk."""
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(
        "project: the dashboard websocket bridge listens locally on port "
        "8050 while jupyter serves notebooks separately",
        repo_root=agent_tree, title="ws port")
    ms.save_typed_memory(
        "project: the dashboard websocket bridge listens locally on port "
        "8060 while jupyter serves notebooks separately",
        repo_root=agent_tree, title="ws port")

    assert len(ms.list_typed_memories(agent_tree)) == 1
    out = _recall(agent_tree)

    assert "port 8060" in out
    assert "this replaced an earlier wording" in out
    assert "port 8050" in out


# ---------------------------------------------------------------------------
# A note about the neighbouring subject
# ---------------------------------------------------------------------------
#
# The preamble said a note may be stale and may be wrong. It did not say a
# note may simply be about something else. Asked which ENCODING to use for
# a bookmark export, kit.glm-5.3 answered with the stored note about the
# export's FORMAT — TSV, columns id/title/url — and said it plainly, with
# no hedge. The note's own "How to apply" line told it to propose that
# format whenever an export came up, which made it the only instruction in
# scope, and nothing said a recalled note has to answer the question that
# was actually asked.

def test_the_preamble_says_a_note_must_fit_the_question(
        agent_tree, tmp_path, monkeypatch):
    from delfin.agent import memory_store as ms

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ms.save_typed_memory(
        "feedback: Lesezeichen immer als TSV exportieren, Spalten "
        "id/title/url, nie als CSV",
        repo_root=agent_tree, title="export convention",
        source=ms.SOURCE_USER)

    out = _recall(agent_tree)
    head = out.split("\n\n", 1)[0].lower()
    assert "not an answer" in head, head
    assert "say so" in head, head
    # And it still arrives before the note it qualifies.
    assert out.lower().index("not an answer") < out.lower().index("tsv")


def test_the_older_framing_did_not_leave_with_it():
    """Three separate rules share one preamble, so a rewrite of any of
    them can silently drop the other two."""
    from delfin.agent.prompt_loader import _MEMORY_BLOCK_PREAMBLE

    lowered = _MEMORY_BLOCK_PREAMBLE.lower()
    for clause in ("background", "not instructions", "out of date",
                   "verify", "not an answer"):
        assert clause in lowered, f"the preamble lost: {clause}"
