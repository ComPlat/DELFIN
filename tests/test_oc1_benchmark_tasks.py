"""OC1 benchmark tasks: graded on the run's physics, not wording.

Wave-24 package OC1 adds three chemistry tasks (tasks_w24_oc1.yaml) over
DELFIN's OCCUPIER workflow. Guards, mirroring
test_chemistry_is_graded_on_the_physics.py:

1. every task resolves its verify script in accept/;
2. every task carries numeric expected_values (a physics band), not
   only signal patterns;
3. the acceptance scripts are importable and both their directions work
   on built fixtures: they REJECT a workspace without a real run and
   ACCEPT a fixture workspace built from real ORCA-output shapes.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

from delfin.agent.benchmark import load_tasks

_PACK = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
         / "pack" / "benchmark")
_OC1_IDS = {"w24-oc1-01", "w24-oc1-02", "w24-oc1-03"}


def _oc1_tasks():
    return [t for t in load_tasks() if t.id in _OC1_IDS]


def test_three_oc1_tasks_load_with_verify_scripts():
    tasks = _oc1_tasks()
    assert {t.id for t in tasks} == _OC1_IDS
    for t in tasks:
        assert t.task_class == "chemistry", f"{t.id} must be chemistry"
        assert t.verify, f"{t.id} has no acceptance script"
        assert (_PACK / "accept" / t.verify).exists(), (
            f"{t.id}: accept/{t.verify} does not resolve")


def test_every_oc1_task_carries_a_physics_band():
    for t in _oc1_tasks():
        assert t.expected_values, (
            f"{t.id} declares no expected_values - prose-graded chemistry")


def _accept_module(name: str):
    sys.path.insert(0, str(_PACK / "accept"))
    try:
        return importlib.import_module(name)
    finally:
        sys.path.pop(0)


def _build_spin_fixture(tmp_path: Path, *, imag: bool = False,
                        mults: str = "multiplicity 5 and multiplicity 3 "
                                     "compared; ground state is multiplicity 6",
                        kcal: bool = True) -> Path:
    root = tmp_path / "oc1" / "spin"
    root.mkdir(parents=True)
    (root / "result.txt").write_text(
        f"{mults}\nrelative energies in kcal/mol given below\n", )
    out = root / "run.out"
    out.write_text(
        "Program Version 6.1.1\n"
        "|  1  2.5  ...  6  ...  DFT  \n"
        "The optimization has converged\n"
        "ORCA TERMINATED NORMALLY\n"
        "------------------\nVIBRATIONAL FREQUENCIES\n------------------\n"
        "0:         0.00 cm**-1\n"
        + ("*         -12.34\n" if imag else "")
        + "1:       123.45 cm**-1\n")
    return tmp_path


def test_spin_accept_rejects_empty_workspace(tmp_path):
    mod = _accept_module("oc1_spin_states_verify")
    ok, msg = mod.check(tmp_path)
    assert not ok


def test_spin_accept_rejects_imaginary_mode(tmp_path):
    ws = _build_spin_fixture(tmp_path, imag=True)
    mod = _accept_module("oc1_spin_states_verify")
    ok, msg = mod.check(ws)
    assert not ok
    assert "imaginary" in msg


def test_spin_accept_accepts_real_shape(tmp_path):
    ws = _build_spin_fixture(tmp_path)
    mod = _accept_module("oc1_spin_states_verify")
    ok, msg = mod.check(ws)
    assert ok, msg


def test_spin_accept_rejects_single_multiplicity(tmp_path):
    ws = _build_spin_fixture(tmp_path,
                             mults="only multiplicity 6 was looked at; "
                                   "ground state is multiplicity 6")
    mod = _accept_module("oc1_spin_states_verify")
    ok, msg = mod.check(ws)
    assert not ok


def test_exchange_accept_rejects_without_run(tmp_path):
    (tmp_path / "oc1" / "broken_sym").mkdir(parents=True)
    mod = _accept_module("oc1_exchange_coupling_verify")
    ok, msg = mod.check(tmp_path)
    assert not ok


def test_exchange_accept_rejects_one_s2_value(tmp_path):
    root = tmp_path / "oc1" / "broken_sym"
    root.mkdir(parents=True)
    (root / "run.out").write_text("ORCA TERMINATED NORMALLY\n")
    (root / "result.txt").write_text(
        "J = -45.2 cm^-1 from high-spin and BS solutions\n"
        "<S^2> = 3.75 for the HS state (site A holds the larger "
        "unpaired-electron count)\n")
    mod = _accept_module("oc1_exchange_coupling_verify")
    ok, msg = mod.check(tmp_path)
    assert not ok
    assert "two" in msg or "<S^2>" in msg


def test_exchange_accept_accepts_both_s2_values(tmp_path):
    root = tmp_path / "oc1" / "broken_sym"
    root.mkdir(parents=True)
    (root / "run.out").write_text("ORCA TERMINATED NORMALLY\n")
    (root / "result.txt").write_text(
        "J = -45.2 cm^-1 via Yamaguchi\n"
        "<S^2> = 3.75 (HS), <S^2> = 1.25 (BS)\n"
        "site A is the centre with the larger unpaired-electron count\n")
    mod = _accept_module("oc1_exchange_coupling_verify")
    ok, msg = mod.check(tmp_path)
    assert ok, msg


def test_redox_accept_rejects_missing_electrode(tmp_path):
    root = tmp_path / "oc1" / "redox"
    root.mkdir(parents=True)
    for name in ("fc.out", "fcp.out"):
        (root / name).write_text(
            "! B3LYP def2-SVP\nORCA TERMINATED NORMALLY\n")
    (root / "result.txt").write_text("E_ox = 0.42 V\nDelta G cycle used\n")
    mod = _accept_module("oc1_redox_verify")
    ok, msg = mod.check(tmp_path)
    assert not ok
    assert "electrode" in msg


def test_redox_accept_rejects_mixed_levels(tmp_path):
    root = tmp_path / "oc1" / "redox"
    root.mkdir(parents=True)
    (root / "fc.out").write_text("! B3LYP def2-SVP\nORCA TERMINATED NORMALLY\n")
    (root / "fcp.out").write_text("! PBE0 def2-SVP\nORCA TERMINATED NORMALLY\n")
    (root / "result.txt").write_text(
        "E_ox = 0.42 V vs SHE\nthermodynamic cycle: G_ox - G_red\n")
    mod = _accept_module("oc1_redox_verify")
    ok, msg = mod.check(tmp_path)
    assert not ok
    assert "mixed" in msg


def test_redox_accept_accepts_same_level(tmp_path):
    root = tmp_path / "oc1" / "redox"
    root.mkdir(parents=True)
    for name in ("fc.out", "fcp.out"):
        (root / name).write_text("! B3LYP def2-SVP\nORCA TERMINATED NORMALLY\n")
    (root / "result.txt").write_text(
        "E_ox = 0.42 V vs SHE\nthermodynamic cycle: G_ox - G_red\n")
    mod = _accept_module("oc1_redox_verify")
    ok, msg = mod.check(tmp_path)
    assert ok, msg
