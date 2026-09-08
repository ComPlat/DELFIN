"""Memory read-back had no coverage.

``dash_memory_remember`` checks that the agent SAVES a fact. Nothing
checked that a fact saved in an earlier session comes back in a later one
— which is the half a user feels: they said it once and expect not to
repeat it.

A task cannot seed its own store, because the store is keyed on the
workspace and lives under ~/.delfin. So a fixture carries its memories as
tracked files and the runner installs them INSIDE the pristine guard: the
seed is there for the attempt and gone with everything else afterwards,
which also means it can never leak into the user's own store.
"""

from pathlib import Path

import pytest

from delfin.agent.benchmark_runner import (_PristineWorkspace,
                                           _seed_fixture_memories,
                                           workspace_for)

_REPO = Path(__file__).resolve().parents[1]
_SEEDS = _REPO / "tests" / "fixtures" / "memory_seed"


def _mem_dir():
    """Bound at call time — the suite redirects HOME by patching this name,
    and a module-level import would hold the unpatched original."""
    from delfin.agent.memory_store import _delfin_memory_dir
    return _delfin_memory_dir


def test_the_seed_directory_is_named_after_its_fixture():
    ws = workspace_for(_REPO, mode="solo", task_class="generic_project")
    assert ws is not None
    assert (_SEEDS / ws.name).is_dir(), sorted(p.name for p in _SEEDS.iterdir())


def test_a_seed_carries_something_the_task_cannot_read_anywhere_else():
    """A seeded memory that only repeats what the files say would let a
    task pass without recalling anything."""
    ws = workspace_for(_REPO, mode="solo", task_class="generic_project")
    seeded = "\n".join(f.read_text() for f in (_SEEDS / ws.name).glob("*.md"))
    assert "TSV" in seeded
    workspace_text = "\n".join(
        f.read_text(errors="replace") for f in ws.rglob("*")
        if f.is_file() and f.suffix in (".py", ".md", ".json", ".txt"))
    assert "TSV" not in workspace_text, "the answer is readable without memory"


def test_the_seed_reaches_the_prompt(tmp_path):
    from delfin.agent.prompt_loader import PromptLoader

    ws = workspace_for(_REPO, mode="solo", task_class="generic_project")
    with _PristineWorkspace(_REPO):
        n = _seed_fixture_memories(_REPO, ws)
        assert n >= 1
        loader = PromptLoader()
        loader.workspace_root = ws
        prompt = loader.build_system_prompt(
            role_id="solo_agent", mode_id="solo",
            task_text="In welchem Format exportiere ich die Lesezeichen?",
            session_key="seed-reaches-prompt", model="kit.deepseek-v4-flash")
    assert "TSV" in prompt, (
        "the seed is in the store but not in the prompt — installing it by "
        "copying the file leaves MEMORY.md, which the recall reads, untouched")


def test_the_seed_does_not_outlive_the_attempt():
    ws = workspace_for(_REPO, mode="solo", task_class="generic_project")
    store = _mem_dir()(ws)
    before = sorted(f.name for f in store.glob("*.md")) if store.is_dir() else []
    with _PristineWorkspace(_REPO):
        _seed_fixture_memories(_REPO, ws)
        during = sorted(f.name for f in store.glob("*.md"))
        assert len(during) > len(before)
    after = sorted(f.name for f in store.glob("*.md")) if store.is_dir() else []
    assert after == before, (before, after)


def test_a_missing_seed_directory_is_not_an_error(tmp_path):
    assert _seed_fixture_memories(tmp_path, tmp_path / "no-such-fixture") == 0
    assert _seed_fixture_memories(tmp_path, None) == 0
