"""Wave 14 U11 — phase 3: nothing the session learned is forgotten.

The durable-fact half of the user's 09.10 requirement: "nothing may be
forgotten". A fact the agent saves via save_typed_memory is written to the
on-disk store and must come back in a LATER, independent read (a later
session, a re-launched dashboard) exactly as written. That is the
persistence property that makes "not forgotten" true across compaction and
across sessions: it does not live in the turn, it lives in the store.

The memory dir resolves under Path.home()/.delfin/projects/<slug>/memory, so
these tests monkeypatch HOME to a temp dir and keep everything inside it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import memory_store as MS


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    # Path.home() reads $HOME on POSIX; re-affirm so the fixture is robust.
    monkeypatch.setattr(Path, "home", classmethod(lambda _cls: tmp_path))
    return tmp_path


def test_a_saved_memory_is_recalled_by_a_later_independent_read(home):
    """Save a durable fact in 'one session', then recall it in a 'later
    session' by re-reading the store from scratch — nothing forgotten."""
    repo = home / "proj"
    repo.mkdir()
    MS.save_typed_memory(
        "the pipeline runs PBE0/def2-SVP in DMF",
        title="method lock", memory_type="project",
        repo_root=repo, scope="project")

    # A later session resolves the same store and reads the fact back.
    records = MS.list_typed_memories(repo, scope="project")
    bodies = [r.get("body", "") for r in records]
    assert any("PBE0/def2-SVP" in b and "DMF" in b for b in bodies)


def test_the_recalled_body_is_exact_not_a_hint(home):
    """'Not forgotten' means the fact itself, not a trace of it: the body
    must come back verbatim, not elided to a headline."""
    repo = home / "proj2"
    repo.mkdir()
    fact = ("background shells keep their output: read with bash_output, "
            "not by re-running the command")
    MS.save_typed_memory(fact, title="bg output", memory_type="feedback",
                         repo_root=repo, scope="project")
    records = MS.list_typed_memories(repo, scope="project")
    assert any(r.get("body", "") == fact for r in records)


def test_the_durable_store_survives_a_discard_of_the_in_memory_copy(home):
    """The memory is on disk, so dropping any in-memory handle cannot lose
    it: re-reading from the filesystem after the writer is gone still finds
    the fact (the compaction boundary neither reads nor writes it away)."""
    repo = home / "proj3"
    repo.mkdir()
    MS.save_typed_memory(
        "the observed floor must not poison later sessions",
        title="floor rule", memory_type="feedback",
        repo_root=repo, scope="project")
    # Nothing is retained in this test process between save and read: the
    # store is the filesystem, so the read below is a fresh resolution.
    store = MS._memory_dir_for_scope(repo, "project")
    files = list(store.glob("*.md")) if store.exists() else []
    assert any("floor" in f.read_text() for f in files if f.name != "MEMORY.md")
