"""Package OR2 (wave 20): spectra, spin states, redox tasks.

The benchmark task file tasks_w20_or2.yaml must exist and follow the
shape the chem suite uses: every task on a system of at most 20 atoms,
every prompt routed through DELFIN's own workflows, and the redox trap
task present — a redox potential must never be reported without its
reference electrode (and never from mixed levels of theory).

Red on main 7adbd4d3: the file does not exist yet.
"""

from __future__ import annotations

import pathlib

import pytest
import yaml

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark"
_OR2 = _PACK / "tasks_w20_or2.yaml"


def _tasks():
    data = yaml.safe_load(_OR2.read_text(encoding="utf-8"))
    return data.get("tasks") or []


def test_the_or2_task_file_exists():
    assert _OR2.is_file(), f"missing: {_OR2}"


def test_five_practice_tasks():
    ids = [t["id"] for t in _tasks()]
    assert 4 <= len(ids) <= 5, ids
    assert len(set(ids)) == len(ids)


def test_the_redox_trap_task_is_there():
    trap = [t for t in _tasks() if t.get("trap")]
    assert trap, "no trap task in OR2"
    assert any("electrode" in str(t.get("trap", "")) for t in trap)


@pytest.mark.parametrize("key", ["id", "prompt", "setup", "verify",
                                 "task_class", "mode"])
def test_every_task_has_the_required_keys(key):
    missing = [t.get("id", "?") for t in _tasks() if key not in t]
    assert not missing, f"tasks without {key}: {missing}"


def test_every_system_is_small():
    """At most 20 atoms: OR2 is the small-systems package."""
    setups = sorted({t["setup"] for t in _tasks()})
    import importlib
    for name in setups:
        mod = importlib.import_module(
            "delfin.agent.pack.benchmark.setup." + name[:-3])
        run = getattr(mod, "run", None)
        assert run is not None, f"{name} has no run()"
        count = getattr(mod, "MAX_ATOMS", None)
        assert count is not None and count <= 20, f"{name}: MAX_ATOMS={count}"


def test_every_setup_and_verify_module_exists():
    for t in _tasks():
        assert (_PACK / "setup" / t["setup"]).is_file(), t["id"]
        assert (_PACK / "accept" / t["verify"]).is_file(), t["id"]
