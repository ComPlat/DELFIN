"""A benchmark prompt became a durable user preference.

`remember` writes to ~/.delfin/projects/<slug>/memory, keyed on the
WORKSPACE — and a dashboard task's workspace is the checkout itself. So a
benchmark task's memories landed in the user's own project store and were
recalled into their real sessions.

Found 2026-09-08 in the live store:

    user_b3lyp-session-setting — "User hat am 2026-09-08 (dieser Sitzung)
    B3LYP als Functional im ORCA Builder gesetzt", use_count 55
    feedback_orca-set-verify — "nach jedem /orca set ein /orca show —
    User will explizit sehen, was wirklich gesetzt wurde", use_count 41

No user said either. `dash_orca_set_functional` and
`workflow_verify_after_modify` did, and the agent then carried them into
sessions as things the user wanted.

The fixture directories were already guarded this way — snapshot before
the attempt, restore after. The memory store was not.
"""

from pathlib import Path

import pytest

from delfin.agent.benchmark_runner import _PristineWorkspace


def _memory_api():
    """Imported inside the test, not at module scope.

    The suite redirects HOME so nothing writes into the real one, by
    patching these names on the module. A module-level import binds the
    unpatched originals, and the guard — which imports at call time — then
    looks at a different directory than the test does.
    """
    from delfin.agent.memory_store import (_delfin_global_memory_dir,
                                           _delfin_memory_dir,
                                           save_typed_memory)
    return _delfin_memory_dir, _delfin_global_memory_dir, save_typed_memory


def test_the_stores_are_among_the_guarded_directories(tmp_path):
    """Computed the same way inside and outside, so a redirected HOME in
    the test environment cannot make this pass or fail by accident."""
    mem_dir, global_dir, _ = _memory_api()
    guard = _PristineWorkspace(tmp_path)
    guarded = {str(p) for p in guard._bases}
    assert str(mem_dir(tmp_path)) in guarded
    assert str(global_dir()) in guarded


def test_a_memory_written_during_a_run_does_not_survive_it(tmp_path):
    mem_dir, _, save_typed_memory = _memory_api()
    store = mem_dir(tmp_path)
    store.mkdir(parents=True, exist_ok=True)
    before = sorted(p.name for p in store.glob("*.md"))
    with _PristineWorkspace(tmp_path):
        save_typed_memory(
            "Written by a test, not by a user.", repo_root=tmp_path,
            memory_type="user", title="guard probe", source="agent")
        during = sorted(p.name for p in store.glob("*.md"))
        assert any("guard-probe" in n for n in during), during
    after = sorted(p.name for p in store.glob("*.md"))
    assert after == before, (before, after)
    assert not any("guard-probe" in n for n in after)


def test_a_memory_the_user_already_had_is_not_lost(tmp_path):
    """Restoring must put back what was there, not empty the store."""
    mem_dir, _, _ = _memory_api()
    store = mem_dir(tmp_path)
    store.mkdir(parents=True, exist_ok=True)
    keeper = store / "zz_guard_keeper_probe.md"
    keeper.write_text("---\nname: keeper\n---\n\nkeep me\n")
    with _PristineWorkspace(tmp_path):
        keeper.unlink()          # a run deletes it
    assert keeper.is_file(), "a pre-existing memory was not restored"
    assert "keep me" in keeper.read_text()


def test_a_store_that_did_not_exist_yet_is_still_guarded(tmp_path):
    """A memory store is created by the first `remember`, so "the
    directory is not there" is the common case and not a corner one. The
    guard used to skip a path that did not exist — nothing to snapshot,
    nothing to restore, and the first thing the run wrote survived it."""
    mem_dir, _, save_typed_memory = _memory_api()
    store = mem_dir(tmp_path)
    assert not store.exists(), "precondition: the store is not there yet"
    with _PristineWorkspace(tmp_path):
        save_typed_memory("Written into a store that did not exist.",
                          repo_root=tmp_path, memory_type="user",
                          title="absent store probe", source="agent")
        assert store.is_dir(), "the run did not create the store"
    assert not store.exists(), "a store the run created outlived it"


def test_two_guards_do_not_erase_a_directory_between_them(tmp_path):
    """A tracked fixture vanished repeatedly on 2026-09-08 in a checkout
    where a parallel suite and a benchmark run at once. The cause was the
    guard having no mutual exclusion: one attempt empties a directory to
    restore it, another snapshots that emptiness in the window, and
    restores it later.

    Two threads here rather than two processes, because the lock is a file
    lock and both hold it the same way — the assertion is that the file
    survives, not how the exclusion is spelled.
    """
    import threading

    ws = tmp_path / "tests" / "fixtures" / "behavior_workspace"
    ws.mkdir(parents=True)
    keeper = ws / "keeper.txt"
    keeper.write_text("must survive")

    errors: list[BaseException] = []

    def _attempt():
        try:
            for _ in range(6):
                with _PristineWorkspace(tmp_path):
                    # What an attempt does: touch the workspace.
                    (ws / "scratch.txt").write_text("x")
        except BaseException as exc:      # noqa: BLE001 - reported below
            errors.append(exc)

    threads = [threading.Thread(target=_attempt) for _ in range(4)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=60)

    assert not errors, errors
    assert keeper.is_file(), "the tracked file was erased between two guards"
    assert keeper.read_text() == "must survive"
    assert not (ws / "scratch.txt").exists(), "an attempt's file outlived it"


def test_a_refused_snapshot_does_not_keep_the_lock(tmp_path, monkeypatch):
    """__exit__ is not called when __enter__ raises, so a lock taken there
    has to be released there too. It was not, and every later attempt then
    waited on a holder that no longer existed — seen as a pytest master
    with one second of CPU in forty-five minutes."""
    import shutil

    ws = tmp_path / "tests" / "fixtures" / "behavior_workspace"
    ws.mkdir(parents=True)
    (ws / "f.txt").write_text("x")

    def _boom(*a, **k):
        raise OSError("copy refused")

    monkeypatch.setattr(shutil, "copytree", _boom)
    guard = _PristineWorkspace(tmp_path)
    with pytest.raises(RuntimeError):
        guard.__enter__()
    assert getattr(guard, "_lock_handle", None) is None, (
        "the lock survived a refused snapshot")

    # And the next guard can still take it.
    monkeypatch.undo()
    with _PristineWorkspace(tmp_path):
        pass
