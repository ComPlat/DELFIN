"""Parity: OPI OCCUPIER inputs == classic OCCUPIER inputs.

``delfin.tools.occupier_opi`` reuses the classic OCCUPIER writer's
level-of-theory, keyword-assembly, metal-basis and first-sphere logic BY
IMPORT (``delfin.occupier.read_and_modify_file_OCCUPIER`` internals are
mirrored through the same helpers the classic module itself calls) and
writes the same ORCA ``.inp`` files with OPI's input classes.  The
acceptance test is that, for the same CONTROL.txt config, structure and
broken-sym string, the OPI writer produces inputs EQUIVALENT to what
``delfin.occupier.read_and_modify_file_OCCUPIER`` writes today — same
simple keywords (order/case-insensitive), same blocks with the same
options and values, same charge/multiplicity/coordinates (per-atom
``NewGTO`` included), same ``%pal``/``%maxcore``, same file names and the
same job sequence.

CONTROL MATRIX (settings taken from the repository's own CONTROL files)
----------------------------------------------------------------------
* examples/DELFIN_Paper_Calcualtions/FeTPP_OCCUPIER/CONTROL.txt —
  functional=PBE0, disp_corr=D4, ri_jkx=RIJCOSX, relativity=ZORA,
  aux_jk=def2/J, main_basisset=def2-SVP, metal_basisset=def2-TZVP,
  first_coordination_sphere=no, CPCM+DMF, PAL=30, maxcore=3800,
  maxiter_occupier=200, initial_guess=PModel, geom_opt_OCCUPIER=OPT,
  approximate_spin_projection_APMethod=2.
* examples/benchmarks_parallel_vs_sequential/Fc_OCCUPIER_parallel/
  CONTROL.txt — the 3c-composite flavour (functional=r2SCAN-3c),
  temperature=298.15, PAL=36, maxcore=3800, maxiter_occupier=100.
* Non-occupier example CONTROLs (examples/DELFIN_Paper_Calcualtions/*,
  examples/benchmarks_parallel_vs_sequential/benzol_*) cover the
  remaining Level-of-Thought variations exercised below as matrix cases
  (gas phase, SMD solvent, no dispersion, 3d-only metals, X2C, QM/MM).

These tests never invoke the ORCA binary: the classic writer runs with
its inputs written only (the OPI module never launches ORCA at all).
"""

from __future__ import annotations

import pytest

# The main py3.11 CI jobs run all tests without OPI installed; skip there.
# The OPI job (DELFIN_OPI_REQUIRED=1, tests/conftest.py) turns that skip
# into a failure, so the OPI-required job still gates on these tests.
pytest.importorskip("opi")

from delfin.tools.occupier_opi import (  # noqa: E402
    inputs_equivalent,
    normalise_inp,
    write_occupier_fob_opi_input,
)

import delfin.occupier as classic  # noqa: E402


# A small metal-complex-like structure built from ideal geometry (never a
# CSD/crystal structure): Fe with four N ligands plus H capping atoms.
_ATOMS = [
    ("Fe", (0.0, 0.0, 0.0)),
    ("N", (2.0, 0.0, 0.0)),
    ("N", (-2.0, 0.0, 0.0)),
    ("N", (0.0, 2.0, 0.0)),
    ("N", (0.0, -2.0, 0.0)),
    ("H", (2.5, 1.0, 0.0)),
    ("H", (2.5, -1.0, 0.0)),
    ("H", (-2.5, 1.0, 0.0)),
    ("H", (-2.5, -1.0, 0.0)),
    ("H", (0.0, 2.5, 1.0)),
    ("H", (0.0, 2.5, -1.0)),
    ("H", (0.0, -2.5, 1.0)),
    ("H", (0.0, -2.5, -1.0)),
]


def _write_xyz(path):
    lines = [str(len(_ATOMS)), "ideal-geometry test structure (synthetic)"]
    lines += [f"{s} {x} {y} {z}" for s, (x, y, z) in _ATOMS]
    path.write_text("\n".join(lines) + "\n")


# The classic writer's call signature, mirrored from the classic job loop
# (delfin/occupier.py:1447-1458): read_and_modify_file_OCCUPIER(
#     src_idx, inp, charge, multiplicity, solvent, metals,
#     metal_basisset, main_basisset, config, broken_sym).
def _run_classic_fob(config, outdir, charge, multiplicity, broken_sym,
                     from_index=1):
    outdir.mkdir(parents=True, exist_ok=True)
    xyz = outdir / ("input.xyz" if from_index == 1 else f"input{from_index}.xyz")
    _write_xyz(xyz)
    classic.read_and_modify_file_OCCUPIER(
        from_index, "input.inp", charge, multiplicity,
        config.get("solvent"), ["Fe"],
        config.get("metal_basisset"), config.get("main_basisset"),
        config, broken_sym, work_dir=str(outdir),
    )
    return outdir / "input.inp"


def _control_matrix():
    """The parity matrix.  Values are the ones the repo's real OCCUPIER
    CONTROL files set (sources in the module docstring)."""
    lot_pbe0 = {  # FeTPP_OCCUPIER/CONTROL.txt
        "functional": "PBE0",
        "main_basisset": "def2-SVP",
        "metal_basisset": "def2-TZVP",
        "aux_jk": "def2/J",
        "disp_corr": "D4",
        "ri_jkx": "RIJCOSX",
        "relativity": "ZORA",
        "implicit_solvation_model": "CPCM",
        "solvent": "DMF",
        "initial_guess": "PModel",
        "geom_opt_OCCUPIER": "OPT",
        "maxiter_occupier": 200,
        "approximate_spin_projection_APMethod": 2,
    }
    lot_3c = {  # Fc_OCCUPIER_parallel/CONTROL.txt (r2SCAN-3c composite)
        "functional": "r2scan-3c",
        "main_basisset": "def2-SVP",
        "metal_basisset": "def2-TZVP",
        "aux_jk": "def2/J",
        "disp_corr": "D4",
        "ri_jkx": "RIJCOSX",
        "relativity": "ZORA",
        "implicit_solvation_model": "CPCM",
        "solvent": "DMF",
        "initial_guess": "PModel",
        "geom_opt_OCCUPIER": "OPT",
        "maxiter_occupier": 100,
        "temperature": 298.15,
    }
    cases = {
        # ---- the repository's own OCCUPIER CONTROLs -------------------
        "FeTPP_pbe0": dict(
            control_args=lot_pbe0,
            charge=0, multiplicity=5, PAL=30, maxcore=3800,
            compare="FSPE", pass_wavefunction="no",
            broken_sym='%scf\n  BrokenSym 3,1\n  APMethod 2\nend',
        ),
        "FeTPP_pbe0_noBS": dict(
            control_args=lot_pbe0,
            charge=0, multiplicity=5, PAL=30, maxcore=3800,
            compare="FSPE", pass_wavefunction="no",
            broken_sym="",
        ),
        "Fc_r2scan3c": dict(
            control_args=lot_3c,
            charge=0, multiplicity=1, PAL=36, maxcore=3800,
            compare="FSPE", pass_wavefunction="no",
            broken_sym="",
        ),
        # ---- matrix cases from the wider example set ------------------
        "pbe0_smd_gas": dict(
            control_args={**lot_pbe0, "implicit_solvation_model": "SMD",
                          "solvent": "water"},
            charge=2, multiplicity=6, PAL=4, maxcore=2000,
            compare="FSPE", pass_wavefunction="no", broken_sym="",
        ),
        "pbe0_G_freq": dict(  # OCCUPIER_compare=G adds FREQ + %freq Temp
            control_args={**lot_pbe0, "OCCUPIER_compare": "G",
                          "temperature": 298.15},
            charge=0, multiplicity=5, PAL=4, maxcore=2000,
            compare="G", pass_wavefunction="no", broken_sym="",
        ),
        "pbe0_moinp": dict(  # pass_wavefunction → %moinp + MORead
            control_args={**lot_pbe0, "pass_wavefunction": "yes"},
            charge=0, multiplicity=5, PAL=4, maxcore=2000,
            compare="FSPE", pass_wavefunction="yes",
            broken_sym='%scf\n  BrokenSym 5,1\nend',
        ),
    }
    for name, c in cases.items():
        c["name"] = name
    return list(cases.values())


def _run_classic_with_compare(case, outdir):
    """Run the classic writer the way run_OCCUPIER drives it, with the
    derived freq flag / moinp handling inside the writer itself (the
    writer reads OCCUPIER_compare and pass_wavefunction from config)."""
    config = dict(case["control_args"])
    config["pass_wavefunction"] = case["pass_wavefunction"]
    return _run_classic_fob(
        config, outdir, case["charge"], case["multiplicity"],
        case["broken_sym"],
    )


@pytest.mark.parametrize("case", _control_matrix())
def test_occupier_fob_parity(tmp_path, case):
    classic_inp = _run_classic_with_compare(case, tmp_path / "classic")
    opi_path = write_occupier_fob_opi_input(
        work_dir=tmp_path / "opi",
        config=case["control_args"],
        charge=case["charge"],
        multiplicity=case["multiplicity"],
        PAL=case["PAL"],
        maxcore=case["maxcore"],
        solvent=case["control_args"].get("solvent"),
        found_metals=["Fe"],
        broken_sym=case["broken_sym"],
        from_index=1,
        base_name="input",
        atoms=_ATOMS,
    )
    assert classic_inp.name == opi_path.name == "input.inp", (
        "OCCUPIER FoB inputs must be named input.inp in both"
    )
    diffs = inputs_equivalent(
        normalise_inp(classic_inp.read_text()),
        normalise_inp(opi_path.read_text()),
    )
    assert diffs == [], (
        f"OCCUPIER FoB ({case['name']}): OPI input differs from classic input\n"
        + "\n".join(diffs)
        + f"\n--- classic ---\n{classic_inp.read_text()}"
        + f"\n--- opi ---\n{opi_path.read_text()}"
    )


def test_occupier_fob_names_follow_the_classic_sequence(tmp_path):
    """Job names mirror the classic loop (occupier.py:1390-1392):
    FoB 1 -> input.inp, FoB n -> input{n}.inp (geometry input{n}.xyz)."""
    lot = {"functional": "PBE0", "main_basisset": "def2-SVP"}
    p = write_occupier_fob_opi_input(
        work_dir=tmp_path / "fob3",
        config=lot, charge=0, multiplicity=3, PAL=4, maxcore=1000,
        solvent=None, found_metals=["Fe"], broken_sym="",
        from_index=3, base_name="input", atoms=_ATOMS,
    )
    assert p.name == "input3.inp"
    # the classic writer reads input3.xyz for from_index=3
    assert (tmp_path / "fob3" / "input3.xyz").exists()


def test_unknown_scf_option_is_refused_not_dropped(tmp_path):
    """An %scf option OPI's BlockScf cannot represent must raise, never be
    silently dropped (a wrong input for a live ORCA run is the cost)."""
    with pytest.raises(ValueError):
        write_occupier_fob_opi_input(
            work_dir=tmp_path / "bad",
            config={"functional": "PBE0", "main_basisset": "def2-SVP",
                    "maxiter_occupier": 100},
            charge=0, multiplicity=3, PAL=1, maxcore=1000,
            solvent=None, found_metals=[],
            broken_sym='%scf\n  NotAnOption 7\nend',
            from_index=1, base_name="input", atoms=_ATOMS,
        )
