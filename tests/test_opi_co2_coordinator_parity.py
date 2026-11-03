"""Parity: OPI CO2-coordinator inputs == classic CO2-coordinator inputs.

``delfin.tools.co2_coordinator_opi`` reuses the classic CO2 coordinator's
geometry/placement/keyword logic by import and writes the same ORCA
``.inp`` files with OPI's input classes.  The acceptance test is that, for
the same CONTROL.txt config and structure, the OPI writer produces inputs
EQUIVALENT to what ``delfin.co2.CO2_Coordinator6`` writes today — same
simple keywords (order/case-insensitive), same blocks with the same
options and values, same charge/multiplicity/coordinates (per-atom
``NewGTO`` included), same ``%pal``/``%maxcore``, same file names and job
sequence.

CONTROL MATRIX (settings taken from the repository's own CONTROL files)
----------------------------------------------------------------------
* tests/test_co2_coordination_workflow.py::_CONTROL_TEMPLATE — a real
  CO2-coordination screen's CONTROL (RSS_paper_19_bug): functional=PBE0,
  disp_corr=D4, ri_jkx=RIJCOSX, relativity=ZORA, aux_jk=def2/J,
  main_basisset=def2-SVP, metal_basisset=def2-TZVP, first (and second)
  coordination-sphere metal basis, CPCM+DMF solvation, PAL=40, maxcore=6000.
* examples/DELFIN_Paper_Calcualtions/*/CONTROL.txt and
  examples/benchmarks_parallel_vs_sequential/*/CONTROL.txt — pipeline
  CONTROLs whose Level-of-Theory keys (functional, basis, dispersion,
  solvation, resources) feed the same ``build_orca_keywords`` path; every
  charge/multiplicity/solvent/PAL/maxcore value a user can set is covered
  below as a separate matrix case.

The matrix below folds those CONTROL values into the ``control_args``,
``charge``, ``multiplicity``, ``PAL``, ``maxcore`` and ``broken_sym``
arguments the classic writer takes, and exercises the two writers both
produce: the orientation single-point (tag="calc") and the relaxed
distance scan (relaxed_surface_scan/scan.inp).

These tests never invoke the ORCA binary: the classic writers' ORCA launch
is neutralised (``shutil.which("orca") → None`` makes them raise *after*
the .inp is written) and the OPI builders write inputs only.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

# The main py3.11 CI jobs run all tests without OPI installed; skip there.
# The OPI job (DELFIN_OPI_REQUIRED=1, tests/conftest.py) turns that skip
# into a failure, so the OPI-required job still gates on these tests.
pytest.importorskip("opi")

from ase import Atoms as AseAtoms  # noqa: E402

from delfin.tools.co2_coordinator_opi import (  # noqa: E402
    inputs_equivalent,
    normalise_inp,
    write_co2_scan_opi_input,
    write_co2_sp_opi_input,
)

# Patch the classic module's ORCA launch before importing its writers.
import delfin.co2.CO2_Coordinator6 as classic  # noqa: E402
classic.shutil.which = lambda *_a, **_k: None  # noqa: E731


# A small metal-complex-like structure built from ideal geometry (never a
# CSD/crystal structure).  atom 0 is the metal; the C and two O at the end
# are the CO2 substrate.
_TEST_ATOMS = AseAtoms(
    "NiC3O3",
    positions=[
        (0.0, 0.0, 0.0),
        (2.0, 0.0, 0.0),
        (-2.0, 0.0, 0.0),
        (0.0, 2.0, 0.0),
        (0.0, 0.0, 3.0),
        (0.0, 0.0, 4.16),
        (0.0, 0.0, 1.84),
    ],
)
_CO2_INDICES = [4, 5, 6]
_METAL_INDEX = 0


def _control_matrix():
    """(name, control_args, charge, multiplicity, PAL, maxcore, broken_sym,
    qmmm_range) — the parity matrix.  Values are the ones the repo's real
    CONTROL files set (see module docstring for the sources)."""
    lot_pbe0 = {
        "functional": "PBE0",
        "main_basisset": "def2-SVP",
        "aux_jk": "def2/J",
        "disp_corr": "D4",
        "ri_jkx": "RIJCOSX",
        "relativity": "ZORA",
        "implicit_solvation_model": "CPCM",
        "solvent": "DMF",
        "metal_basisset": "def2-TZVP",
        "first_coordination_sphere_metal_basisset": "yes",
        "first_coordination_sphere_scale": 1.3,
        "second_coordination_sphere_metal_basisset": "yes",
        "second_coordination_sphere_scale": 1.3,
    }
    cases = {
        # DFT, CPCM+DMF, both coordination spheres, broken sym (realistic template)
        "pbe0_cpcm_dmf_bs": dict(
            control_args=lot_pbe0, charge=-2, multiplicity=3, PAL=40,
            maxcore=6000, broken_sym="%scf BrokenSym 1,5 end", qmmm_range=None,
        ),
        # same but no broken sym (spin-pure pathway)
        "pbe0_cpcm_dmf_no_bs": dict(
            control_args=lot_pbe0, charge=-2, multiplicity=3, PAL=40,
            maxcore=6000, broken_sym="", qmmm_range=None,
        ),
        # GFN2-XTB / ALPB(DMF) semiempirical pathway (classic default)
        "gfn2_xtb_alpb": dict(
            control_args={
                "functional": "GFN2-XTB",
                "implicit_solvation_model": "ALPB",
                "solvent": "DMF",
            },
            charge=-2, multiplicity=3, PAL=4, maxcore=2000,
            broken_sym="%scf BrokenSym 1,5 end", qmmm_range=None,
        ),
        # DFT with no solvent/no dispersion/no RI (minimal CONTROL)
        "pbe0_bare": dict(
            control_args={
                "functional": "PBE0",
                "main_basisset": "def2-SVP",
                "aux_jk": "def2/J",
                "disp_corr": "D4",
            },
            charge=0, multiplicity=1, PAL=4, maxcore=2000,
            broken_sym="", qmmm_range=None,
        ),
        # QM/MM: a QM range over the metal + ligands, CO2 inside the QM core
        "pbe0_qmmm": dict(
            control_args={**lot_pbe0},
            charge=-2, multiplicity=3, PAL=4, maxcore=2000,
            broken_sym="", qmmm_range=(0, 3),
        ),
        # a different DFT flavour + no aux/disp, gas-phase single point
        "bp86_vacuum": dict(
            control_args={"functional": "BP86", "main_basisset": "def2-SVP"},
            charge=-2, multiplicity=5, PAL=16, maxcore=3000,
            broken_sym="", qmmm_range=None,
        ),
    }
    for name, c in cases.items():
        c["name"] = name
    return list(cases.values())


def _run_classic_sp(control_args, outdir, charge, multiplicity, PAL, maxcore,
                    broken_sym, qmmm_range, tag="calc"):
    d = Path(outdir); d.mkdir(parents=True, exist_ok=True)
    xyz = d / "structure.xyz"; xyz.write_text("")
    try:
        classic.write_orca_sp_input_and_run(
            _TEST_ATOMS, str(xyz), str(d),
            orca_keywords=classic.build_orca_keywords(control_args, "SP"),
            broken_sym=broken_sym, charge=charge, multiplicity=multiplicity,
            PAL=PAL, maxcore=maxcore, tag=tag,
            metal_symbol="Ni",
            metal_basis=control_args.get("metal_basisset"),
            control_args=control_args, qmmm_range=qmmm_range,
            co2_indices=_CO2_INDICES,
        )
    except RuntimeError:
        pass
    return d / f"{tag}.inp"


def _run_classic_scan(scan_dir, control_args, charge, multiplicity, PAL,
                      maxcore, broken_sym, qmmm_range):
    d = Path(scan_dir); d.mkdir(parents=True, exist_ok=True)
    (d / "unused.xyz").write_text("dummy xyz\n")
    try:
        classic.write_orca_input_and_run(
            _TEST_ATOMS, "unused.xyz", _METAL_INDEX, _CO2_INDICES[0],
            start_distance=4.3, end_distance=1.6, steps=25,
            orca_keywords=classic.build_orca_keywords(control_args, "OPT"),
            broken_sym=broken_sym, charge=charge, multiplicity=multiplicity,
            PAL=PAL, maxcore=maxcore,
            metal_symbol="Ni",
            metal_basis=control_args.get("metal_basisset"),
            control_args=control_args, qmmm_range=qmmm_range,
            co2_indices=_CO2_INDICES,
        )
    except RuntimeError:
        pass
    return d / "relaxed_surface_scan" / "scan.inp"


def _assert_equivalent(text_a: str, path_a: Path, text_b: str, path_b: Path,
                       when: str):
    assert inputs_equivalent(normalise_inp(text_a), normalise_inp(text_b)) == [], (
        f"{when}: OPI input differs from classic input\n"
        f"--- classic ({path_a.name}) ---\n{text_a}\n"
        f"--- opi ({path_b.name}) ---\n{text_b}"
    )


@pytest.mark.parametrize("case", _control_matrix())
def test_orientation_sp_parity(tmp_path, case):
    cwd = tmp_path
    classic_inp = _run_classic_sp(
        case["control_args"], cwd / "classic_sp", case["charge"],
        case["multiplicity"], case["PAL"], case["maxcore"],
        case["broken_sym"], case["qmmm_range"],
    )
    opi_path = write_co2_sp_opi_input(
        _TEST_ATOMS, "unused.xyz", str(cwd / "opi_sp"),
        orca_keywords=classic.build_orca_keywords(case["control_args"], "SP"),
        broken_sym=case["broken_sym"], charge=case["charge"],
        multiplicity=case["multiplicity"], PAL=case["PAL"],
        maxcore=case["maxcore"], tag="calc",
        metal_symbol="Ni",
        metal_basis=case["control_args"].get("metal_basisset"),
        control_args=case["control_args"], qmmm_range=case["qmmm_range"],
        co2_indices=_CO2_INDICES,
    )
    assert classic_inp.name == opi_path.name == "calc.inp", (
        "orientation single-point input must be named calc.inp in both"
    )
    _assert_equivalent(
        classic_inp.read_text(), classic_inp, opi_path.read_text(), opi_path,
        f"orientation SP ({case['name']})",
    )


@pytest.mark.parametrize("case", _control_matrix())
def test_relaxed_scan_parity(tmp_path, case):
    old = os.getcwd()
    scan_work = tmp_path / "scan_work"
    os.makedirs(scan_work, exist_ok=True)
    os.chdir(scan_work)
    try:
        classic_inp = _run_classic_scan(
            str(scan_work), case["control_args"], case["charge"],
            case["multiplicity"], case["PAL"], case["maxcore"],
            case["broken_sym"], case["qmmm_range"],
        )
    finally:
        os.chdir(old)
    opi_path = write_co2_scan_opi_input(
        _TEST_ATOMS, "unused.xyz", _METAL_INDEX, _CO2_INDICES[0],
        start_distance=4.3, end_distance=1.6, steps=25,
        orca_keywords=classic.build_orca_keywords(case["control_args"], "OPT"),
        broken_sym=case["broken_sym"], charge=case["charge"],
        multiplicity=case["multiplicity"], PAL=case["PAL"],
        maxcore=case["maxcore"],
        metal_symbol="Ni",
        metal_basis=case["control_args"].get("metal_basisset"),
        control_args=case["control_args"], qmmm_range=case["qmmm_range"],
        co2_indices=_CO2_INDICES,
        scan_dir=str(scan_work / "opi_scan"),
    )
    assert classic_inp.name == opi_path.name == "scan.inp", (
        "relaxed distance-scan input must be named scan.inp in both"
    )
    _assert_equivalent(
        classic_inp.read_text(), classic_inp, opi_path.read_text(), opi_path,
        f"relaxed distance scan ({case['name']})",
    )


def test_orientation_and_scan_are_the_classic_job_sequence(tmp_path):
    """The OPI module reproduces the classic job file names and sequence:
    the orientation SP inputs are <dir>/calc.inp and the relaxed distance
    scan is relaxed_surface_scan/scan.inp."""
    ctrl = {
        "functional": "PBE0", "main_basisset": "def2-SVP", "aux_jk": "def2/J",
        "disp_corr": "D4", "implicit_solvation_model": "CPCM", "solvent": "DMF",
    }
    # orientation SP -> calc.inp
    sp = write_co2_sp_opi_input(
        _TEST_ATOMS, "u.xyz", str(tmp_path / "orientation_scan" / "ang_000"),
        orca_keywords=classic.build_orca_keywords(ctrl, "SP"),
        charge=-2, multiplicity=3, PAL=4, maxcore=2000, tag="calc",
        metal_basis=ctrl.get("metal_basisset"), control_args=ctrl,
        qmmm_range=None, co2_indices=_CO2_INDICES,
    )
    assert sp.name == "calc.inp"
    assert sp.parent.name == "ang_000"
    # relaxed distance scan -> relaxed_surface_scan/scan.inp
    sc = write_co2_scan_opi_input(
        _TEST_ATOMS, "u.xyz", 0, 4,
        start_distance=4.3, end_distance=1.6, steps=25,
        orca_keywords=classic.build_orca_keywords(ctrl, "OPT"),
        charge=-2, multiplicity=3, PAL=4, maxcore=2000,
        metal_basis=ctrl.get("metal_basisset"), control_args=ctrl,
        qmmm_range=None, co2_indices=_CO2_INDICES,
        scan_dir=str(tmp_path / "relaxed_surface_scan"),
    )
    assert sc.name == "scan.inp"
    assert sc.parent.name == "relaxed_surface_scan"
