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
