"""WB0 phase 2 red control: the workflow-building task family.

The tasks load (benchmark.load_tasks auto-discovers tasks_*.yaml), their
fixtures are written by the setup script into a throwaway workspace, and
the acceptance script grades a built spec through the real grader. The
red control also pins the physically wrong spec a careless builder
writes being REFUSED by the verify script, and the count dict growing by
exactly WB0's tasks.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import yaml

import pytest

from delfin.agent.benchmark import load_tasks

_ROOT = Path(__file__).resolve().parents[1]
_YAML = (_ROOT / "delfin" / "agent" / "pack" / "benchmark"
         / "tasks_w21_build.yaml")
_SETUP = _YAML.parent / "setup" / "w21_build_fixtures.py"
_VERIFY = _YAML.parent / "accept" / "w21_build_verify.py"

# The correct spec for the conformers brief: conformer search first, then
# one consistent level of theory with thermochemistry.
_GOOD_CONFORMERS_SPEC = {
    "name": "ethanol_conformer_thermo",
    "steps": [
        {"step": "crest_conformers", "kwargs": {"charge": 0}},
        {"step": "orca_opt_freq", "kwargs": {"charge": 0, "mult": 1,
                                             "method": "PBE0",
                                             "basis": "def2-TZVP",
                                             "freq_temp": 298.15,
                                             "freq_pressure": 1.0}},
    ],
}

# The physically wrong spec for the same brief: straight to DFT, no
# conformer step.
_BAD_CONFORMERS_SPEC = {
    "name": "flexible_noconf",
    "steps": [
        {"step": "orca_opt_freq", "kwargs": {"charge": 0, "mult": 1,
                                             "method": "PBE0",
                                             "basis": "def2-TZVP"}},
    ],
}


def _wb0_tasks():
    raw = yaml.safe_load(_YAML.read_text(encoding="utf-8"))["tasks"]
    return [t for t in raw if str(t.get("setup", "")).startswith("w21_")]


def test_the_family_loads_through_the_real_loader():
    ids = {t.id for t in load_tasks(path=_YAML)}
    assert {"w21_build_conformers", "w21_build_redox", "w21_build_ts",
            "w21_build_spin", "w21_build_cycle",
            "w21_build_spectrum"} <= ids


def test_every_task_is_verified_and_setup_resolves():
    raw = _wb0_tasks()
    assert len(raw) >= 8
    for t in raw:
        assert t["task_class"] == "science_analysis"
        assert t["verify"] in ("w21_build_verify.py", "w21_run_verify.py")
        assert (_VERIFY.parent / t["verify"]).exists()
        assert (_SETUP.parent / t["setup"]).exists()
        # prompts read like a chemist's request: no template names
        text = t["prompt"].lower()
        for banned in ("classic_opt_freq", "redox_potential",
                       "template", "multi_level_energy"):
            assert banned not in text, (t["id"], banned)


def test_setup_refuses_to_overwrite_and_writes_fixtures(tmp_path):
    rc = subprocess.run(
        [sys.executable, str(_SETUP), str(tmp_path)],
        capture_output=True, text=True)
    assert rc.returncode == 0, rc.stderr
    folders = sorted(p.name for p in tmp_path.iterdir())
    assert folders == sorted(
        t["id"] for t in _wb0_tasks()
        if t["setup"] == "w21_build_fixtures.py")
    d = tmp_path / "w21_build_conformers"
    text = (d / "start.xyz").read_text(encoding="utf-8")
    assert text.splitlines()[0] == "9"
    assert "charge=0 mult=1" in text.splitlines()[1]
    assert (d / "brief.txt").is_file()

    # second run refuses
    rc = subprocess.run(
        [sys.executable, str(_SETUP), str(tmp_path)],
        capture_output=True, text=True)
    assert rc.returncode == 1, rc.stdout


def test_setup_without_argument_prints_usage_and_exits_nonzero():
    rc = subprocess.run([sys.executable, str(_SETUP)],
                        capture_output=True, text=True)
    assert rc.returncode != 0
    assert "usage" in (rc.stderr + rc.stdout).lower()


# ---------------------------------------------------------------------------
# The acceptance script through the real grader
# ---------------------------------------------------------------------------

def _run_verify(tmp_path, spec, *, with_brief=True):
    d = tmp_path / "w21_build_conformers"
    d.mkdir()
    if with_brief:
        (d / "brief.txt").write_text(
            "ethanol is flexible. I need its conformer landscape mapped "
            "cheaply, then the thermochemistry of the lowest conformer at "
            "one consistent level of theory.", encoding="utf-8")
    (d / "spec.json").write_text(json.dumps(spec), encoding="utf-8")
    (d / "answer.txt").write_text("no deviations from the defaults needed.",
                                  encoding="utf-8")
    return subprocess.run(
        [sys.executable, str(_VERIFY), str(tmp_path),
         "w21_build_conformers"],
        capture_output=True, text=True)


def test_verify_accepts_a_good_spec(tmp_path):
    rc = _run_verify(tmp_path, _GOOD_CONFORMERS_SPEC)
    assert rc.returncode == 0, rc.stderr


def test_verify_refuses_a_physically_wrong_spec(tmp_path):
    rc = _run_verify(tmp_path, _BAD_CONFORMERS_SPEC)
    assert rc.returncode != 0
    assert "conformer" in (rc.stderr + rc.stdout).lower()


def test_verify_refuses_a_missing_brief(tmp_path):
    """Without the chemist's brief the request-dependent rules cannot
    grade — the verify script must refuse, not grade against an empty
    request."""
    rc = _run_verify(tmp_path, _GOOD_CONFORMERS_SPEC, with_brief=False)
    assert rc.returncode == 1
    assert "brief" in (rc.stderr + rc.stdout).lower()


def test_verify_refuses_a_missing_spec(tmp_path):
    (tmp_path / "w21_build_conformers").mkdir()
    rc = subprocess.run(
        [sys.executable, str(_VERIFY), str(tmp_path),
         "w21_build_conformers"], capture_output=True, text=True)
    assert rc.returncode == 1


def test_verify_rejects_an_unknown_task_id(tmp_path):
    rc = subprocess.run(
        [sys.executable, str(_VERIFY), str(tmp_path), "w21_build_none"],
        capture_output=True, text=True)
    assert rc.returncode == 2, rc.stderr


# ---------------------------------------------------------------------------
# Phase 3: the real-run tasks
# ---------------------------------------------------------------------------

_RUN_SETUP = _YAML.parent / "setup" / "w21_run_fixtures.py"
_RUN_VERIFY = _YAML.parent / "accept" / "w21_run_verify.py"

# A good spec for the methanol run brief: opt+freq at the CONTROL's
# level, one consistent method, nothing exotic.
_GOOD_RUN_SPEC = {
    "name": "methanol_opt_freq",
    "steps": [
        {"step": "orca_opt_freq", "kwargs": {"charge": 0, "mult": 1,
                                             "method": "B3LYP",
                                             "basis": "def2-SVP",
                                             "freq_temp": 298.15,
                                             "freq_pressure": 1.0}},
    ],
}
# A wrong spec for the same brief: anion charge with a diffuse-free
# basis (the CONTROL pins charge 0; both the grader's physics rule and
# CB0A's charge check must catch the divergence).
_BAD_RUN_SPEC = {
    "name": "methanol_wrong_charge",
    "steps": [
        {"step": "orca_opt_freq", "kwargs": {"charge": -1, "mult": 1,
                                             "method": "B3LYP",
                                             "basis": "def2-SVP"}},
    ],
}


def test_run_setup_refuses_to_overwrite_and_writes_control(tmp_path):
    rc = subprocess.run(
        [sys.executable, str(_RUN_SETUP), str(tmp_path)],
        capture_output=True, text=True)
    assert rc.returncode == 0, rc.stderr
    assert (tmp_path / "w21_run_methanol" / "CONTROL.txt").is_file()
    control = (tmp_path / "w21_run_methanol" / "CONTROL.txt").read_text(
        encoding="utf-8")
    assert "charge = 0" in control and "mult = 1" in control
    assert (tmp_path / "w21_run_methanol" / "brief.txt").is_file()

    rc = subprocess.run(
        [sys.executable, str(_RUN_SETUP), str(tmp_path)],
        capture_output=True, text=True)
    assert rc.returncode == 1, rc.stdout


def test_run_setup_without_argument_prints_usage_and_exits_nonzero():
    rc = subprocess.run([sys.executable, str(_RUN_SETUP)],
                        capture_output=True, text=True)
    assert rc.returncode != 0
    assert "usage" in (rc.stderr + rc.stdout).lower()


def _prep_run(tmp_path, spec, *, result=None, control=None):
    d = tmp_path / "w21_run_methanol"
    d.mkdir()
    (d / "brief.txt").write_text(
        "Optimise methanol in the gas phase and confirm the geometry is "
        "a true minimum; report the thermochemical correction at 298 K.",
        encoding="utf-8")
    (d / "CONTROL.txt").write_text(control or (
        "workflow = classic_opt_freq\ncharge = 0\nmult = 1\n"
        "method = B3LYP\nbasis = def2-SVP\n"), encoding="utf-8")
    (d / "spec.json").write_text(json.dumps(spec), encoding="utf-8")
    (d / "answer.txt").write_text("no deviations from the defaults.",
                                  encoding="utf-8")
    if result is not None:
        (d / "result.json").write_text(json.dumps(result),
                                       encoding="utf-8")
    return subprocess.run(
        [sys.executable, str(_RUN_VERIFY), str(tmp_path),
         "w21_run_methanol"], capture_output=True, text=True)


def test_run_verify_refuses_a_faked_run(tmp_path):
    """The spec alone is not enough: without a real run's inputs and
    outputs (run_evidence.collect finds none) the verify script must
    refuse — CB0A's whole-run rule."""
    rc = _prep_run(tmp_path, _GOOD_RUN_SPEC,
                   result={"energy_hartree": -115.7, "n_imag": 0.0,
                           "thermal_corr_298k": 0.001})
    assert rc.returncode == 1
    assert "evidence" in (rc.stderr + rc.stdout).lower()


def test_run_verify_refuses_a_physically_wrong_spec(tmp_path):
    """A wrong spec fails on the SPEC half before the run half is even
    graded."""
    rc = _prep_run(tmp_path, _BAD_RUN_SPEC,
                   result={"energy_hartree": -115.7, "n_imag": 0.0,
                           "thermal_corr_298k": 0.001})
    assert rc.returncode == 1
    assert "anion" in (rc.stderr + rc.stdout).lower()


def test_run_verify_rejects_an_unknown_task_id(tmp_path):
    rc = subprocess.run(
        [sys.executable, str(_RUN_VERIFY), str(tmp_path), "w21_run_none"],
        capture_output=True, text=True)
    assert rc.returncode == 2, rc.stderr


def test_run_verify_rejects_missing_spec(tmp_path):
    d = tmp_path / "w21_run_methanol"
    d.mkdir()
    rc = subprocess.run(
        [sys.executable, str(_RUN_VERIFY), str(tmp_path),
         "w21_run_methanol"], capture_output=True, text=True)
    assert rc.returncode == 1


# ---------------------------------------------------------------------------
# The count dict carries WB0's own tasks
# ---------------------------------------------------------------------------

def test_count_dict_includes_wb0_tasks():
    from tests.test_every_task_that_writes_writes_where_it_belongs import (
        _tasks as _count_tasks, )
    raw = {t["id"]: t for t in _count_tasks()}
    for t in _wb0_tasks():
        assert t["id"] in raw, f"{t['id']} not seen by the count test"
        assert raw[t["id"]].get("task_class") == "science_analysis"
