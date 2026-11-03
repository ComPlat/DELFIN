"""Parity: OPI ESD inputs == classic ESD inputs.

``delfin.tools.esd_opi`` reuses the classic ESD writer's level-of-theory,
keyword-assembly, solvation, TDDFT-block, ESD-block and per-atom-NewGTO
logic BY IMPORT (the same helpers ``delfin.esd_input_generator`` itself
calls) and writes the same ORCA ``.inp`` files with OPI's input classes.
The acceptance test is that, for the same CONTROL config, state set-up and
geometry, the OPI writer produces inputs EQUIVALENT to what
``delfin.esd_input_generator`` writes today — same simple keywords (order/
case-insensitive), same blocks with the same options and values, same
charge/multiplicity/coordinates (per-atom ``NewGTO`` included), same
``%pal``/``%maxcore``, same file names.

NORMALISATION (stated, as required): the comparison is structural, not
byte-identical — ``delfin.tools.co2_coordinator_opi.inputs_equivalent``
over ``normalise_inp`` of both texts.  Token order inside the ``!`` line,
block order and layout spacing therefore do not matter; token SET, block
option SET and values do.  Multi-job files ($new_job, e.g. TDDFT S0 and
PHOSP iroot jobs) are normalised per job section and compared
section-by-section.

CONTROL MATRIX (settings taken from the repository's own CONTROL files)
----------------------------------------------------------------------
* examples/benchmarks_parallel_vs_sequential/benzol_classic_sequential/
  CONTROL.txt — functional=r2scan-3c, disp_corr=D4, ri_jkx=RIJCOSX,
  relativity=ZORA, aux_jk=def2/J, main_basisset=def2-SVP,
  metal_basisset=def2-TZVP, CPCM/Acetonitrile, PAL=38, maxcore=3800,
  ESD_modus=TDDFT, NROOTS=15, TDA=FALSE, DOSOC=FALSE, NACME=TRUE,
  print_MOs=yes.
* Non-ESD example CONTROLs cover the remaining level-of-theory variations
  exercised below as matrix cases (gas phase, SMD solvent, no dispersion,
  3d-only metals, numFREQ, deltaSCF modus).

These tests never invoke the ORCA binary: the classic writer runs with
its inputs written only (the OPI module never launches ORCA at all).
"""

from __future__ import annotations

import pytest

# The main py3.11 CI jobs run all tests without OPI installed; skip there.
# The OPI job (DELFIN_OPI_REQUIRED=1, tests/conftest.py) turns that skip
# into a failure, so the OPI-required job still gates on these tests.
pytest.importorskip("opi")

from delfin.tools.co2_coordinator_opi import normalise_inp  # noqa: E402
from delfin.tools.esd_opi import (  # noqa: E402
    esd_inputs_equivalent,
    split_jobs,
    write_esd_state_opi_input,
    write_esd_rate_opi_input,
)

import delfin.esd_input_generator as classic  # noqa: E402


# A small metal-complex-like structure built from ideal geometry (never a
# CSD/crystal structure): Ir with two N ligands plus H capping atoms.
_ATOMS = [
    ("Ir", (0.0, 0.0, 0.0)),
    ("N", (2.0, 0.0, 0.0)),
    ("N", (-2.0, 0.0, 0.0)),
    ("H", (2.5, 1.0, 0.0)),
    ("H", (2.5, -1.0, 0.0)),
    ("H", (0.0, 2.5, 1.0)),
]


def _write_esd_dir(esd_dir, states=("S0", "S1", "T1")):
    """Set up an ESD directory the way setup_esd_directory would: the
    state xyz files the classic writer reads coordinates from."""
    esd_dir.mkdir(parents=True, exist_ok=True)
    for st in states:
        lines = ["6", f"synthetic {st} geometry"] + [
            f"{el} {x} {y} {z}" for el, (x, y, z) in _ATOMS
        ]
        (esd_dir / f"{st}.xyz").write_text("\n".join(lines) + "\n")


def _run_classic_state(case, esd_dir):
    return Path_text(classic.create_state_input(
        state=case["state"],
        esd_dir=esd_dir,
        charge=case["charge"],
        solvent=case["control_args"].get("solvent"),
        metals=case["metals"],
        main_basisset=case["control_args"]["main_basisset"],
        metal_basisset=case["control_args"].get("metal_basisset"),
        config=case["control_args"],
    ))


def _run_classic_rate(case, esd_dir):
    return Path_text(case["classic_fn"](
        esd_dir=esd_dir,
        charge=case["charge"],
        solvent=case["control_args"].get("solvent"),
        metals=case["metals"],
        main_basisset=case["control_args"]["main_basisset"],
        metal_basisset=case["control_args"].get("metal_basisset"),
        config=case["control_args"],
        **case.get("classic_kwargs", {}),
    ))


def Path_text(fn_result):
    from pathlib import Path
    return Path(fn_result).read_text()


def _run_opi_state(case, esd_dir, out_dir):
    return Path_text(write_esd_state_opi_input(
        esd_dir=esd_dir,
        out_dir=out_dir,
        state=case["state"],
        charge=case["charge"],
        solvent=case["control_args"].get("solvent"),
        metals=case["metals"],
        main_basisset=case["control_args"]["main_basisset"],
        metal_basisset=case["control_args"].get("metal_basisset"),
        config=case["control_args"],
    ))


def _run_opi_rate(case, esd_dir, out_dir):
    return Path_text(write_esd_rate_opi_input(
        esd_dir=esd_dir,
        out_dir=out_dir,
        kind=case["kind"],
        pair=case.get("pair"),
        charge=case["charge"],
        solvent=case["control_args"].get("solvent"),
        metals=case["metals"],
        main_basisset=case["control_args"]["main_basisset"],
        metal_basisset=case["control_args"].get("metal_basisset"),
        config=case["control_args"],
        trootssl=case.get("trootssl"),
        initial_state=case.get("initial_state", "S1"),
    ))


def _matrix_lots():
    """Level-of-theory sets from the repo's own CONTROL files."""
    lot_benzol = {  # benzol_classic_sequential/CONTROL.txt
        "functional": "r2scan-3c",
        "disp_corr": "D4",
        "ri_jkx": "RIJCOSX",
        "relativity": "ZORA",
        "aux_jk": "def2/J",
        "main_basisset": "def2-SVP",
        "main_basisset_rel": "ZORA-def2-SVP",
        "metal_basisset": "def2-TZVP",
        "metal_basisset_rel": "SARC-ZORA-TZVP",
        "first_coordination_sphere_metal_basisset": "no",
        "first_coordination_sphere_scale": "1.20",
        "geom_opt": "OPT",
        "initial_guess": "PModel",
        "temperature": "298.15",
        "implicit_solvation_model": "CPCM",
        "solvent": "Acetonitrile",
        "PAL": 38,
        "maxcore": 3800,
        "ESD_modus": "TDDFT",
        "NROOTS": 15,
        "TDDFT_TDA": "FALSE",
        "TDDFT_SOC": "FALSE",
        "NACME": "TRUE",
        "ETF": "TRUE",
        "print_MOs": "yes",
    }
    lot_pbe0 = {  # matrix case from the DELFIN_Paper_Calcualtions flavour
        "functional": "PBE0",
        "disp_corr": "D4",
        "ri_jkx": "RIJCOSX",
        "relativity": "ZORA",
        "aux_jk": "def2/J",
        "main_basisset": "def2-SVP",
        "metal_basisset": "def2-TZVP",
        "geom_opt": "OPT",
        "initial_guess": "PModel",
        "implicit_solvation_model": "SMD",
        "solvent": "water",
        "PAL": 4,
        "maxcore": 2000,
        "ESD_modus": "TDDFT",
        "ESD_frequency": "no",
        "NROOTS": 6,
    }
    return lot_benzol, lot_pbe0


def _control_matrix():
    lot_benzol, lot_pbe0 = _matrix_lots()
    cases = []
    # ---- state jobs (TDDFT modus) -------------------------------------
    for st in ("S0", "S1", "T1", "T2"):
        cases.append(dict(
            name=f"tddft_state_{st.lower()}_benzol",
            state=st, charge=0, control_args=lot_benzol, metals=["Ir"],
        ))
    # ---- matrix variations ---------------------------------------------
    cases.append(dict(
        name="state_S1_smd_gas_pbe0",
        state="S1", charge=0, control_args=lot_pbe0, metals=[],
    ))
    cases.append(dict(
        name="state_S0_no_disp",
        state="S0", charge=-1,
        control_args={**lot_pbe0, "disp_corr": "NONE", "relativity": "",
                      "aux_jk": "", "solvent": ""},
        metals=["Fe"],
    ))
    cases.append(dict(
        name="state_S1_numfreq",
        state="S1", charge=0,
        control_args={**lot_pbe0, "ESD_frequency": "yes",
                      "freq_type": "numFREQ"},
        metals=[],
    ))
    cases.append(dict(
        name="state_S1_deltascf_modus",
        state="S1", charge=0,
        control_args={**lot_pbe0, "ESD_modus": "deltaSCF"},
        metals=[],
    ))
    # ---- rate jobs -------------------------------------------------------
    cases.append(dict(
        name="isc_S1_T1", kind="isc", pair="S1>T1", charge=0,
        control_args=lot_benzol, metals=["Ir"], trootssl=0,
        classic_fn=classic.create_isc_input,
        classic_kwargs={"isc_pair": "S1>T1", "trootssl": 0},
    ))
    cases.append(dict(
        name="ic_S1_S0", kind="ic", pair="S1>S0", charge=0,
        control_args=lot_benzol, metals=["Ir"],
        classic_fn=classic.create_ic_input,
        classic_kwargs={"ic_pair": "S1>S0"},
    ))
    cases.append(dict(
        name="fluor_S1_S0", kind="fluor", pair="S1>S0", charge=0,
        control_args=lot_benzol, metals=["Ir"],
        classic_fn=classic.create_fluor_input,
        classic_kwargs={"initial_state": "S1", "final_state": "S0"},
    ))
    cases.append(dict(
        name="phosp_T1_S0", kind="phosp", pair="T1>S0", charge=0,
        control_args=lot_benzol, metals=["Ir"],
        classic_fn=classic.create_phosp_input,
        classic_kwargs={"initial_state": "T1", "final_state": "S0"},
    ))
    return cases


@pytest.mark.parametrize("case", _control_matrix(), ids=lambda c: c["name"])
def test_esd_input_parity(tmp_path, case):
    esd_dir = tmp_path / "classic"
    _write_esd_dir(esd_dir, states=("S0", "S1", "T1", "T2"))
    if "kind" not in case:
        classic_text = _run_classic_state(case, esd_dir)
        opi_text = _run_opi_state(case, esd_dir, tmp_path / "opi")
    else:
        classic_text = _run_classic_rate(case, esd_dir)
        opi_text = _run_opi_rate(case, esd_dir, tmp_path / "opi")
    diffs = esd_inputs_equivalent(
        [normalise_inp(t) for t in split_jobs(classic_text)],
        [normalise_inp(t) for t in split_jobs(opi_text)],
    )
    assert diffs == [], (
        f"ESD ({case['name']}): OPI input differs from classic input\n"
        + "\n".join(diffs)
        + f"\n--- classic ---\n{classic_text}"
        + f"\n--- opi ---\n{opi_text}"
    )


def test_esd_opi_file_names_follow_the_classic_ones(tmp_path):
    esd_dir = tmp_path / "classic"
    _write_esd_dir(esd_dir, states=("S0", "S1", "T1"))
    lot = {"functional": "PBE0", "main_basisset": "def2-SVP",
           "ESD_modus": "TDDFT"}
    from pathlib import Path
    p = write_esd_state_opi_input(
        esd_dir=esd_dir, out_dir=tmp_path / "opi", state="S1",
        charge=0, solvent=None, metals=[], main_basisset="def2-SVP",
        metal_basisset=None, config=lot,
    )
    assert Path(p).name == "S1.inp"
    p2 = write_esd_rate_opi_input(
        esd_dir=esd_dir, out_dir=tmp_path / "opi2", kind="isc",
        pair="S1>T1", charge=0, solvent=None, metals=[],
        main_basisset="def2-SVP", metal_basisset=None, config=lot,
        trootssl=0,
    )
    assert Path(p2).name == "S1_T1_ISC_m0.inp"


def test_esd_opi_unknown_scf_option_refused(tmp_path):
    """An option OPI's blocks cannot represent must raise, never be
    silently dropped (a wrong input for a live ORCA run is the cost)."""
    esd_dir = tmp_path / "classic"
    _write_esd_dir(esd_dir)
    # T1 tddft-opt path exercises BlockScf mapping; a bogus option in the
    # %scf block must raise.
    with pytest.raises(ValueError):
        write_esd_state_opi_input(
            esd_dir=esd_dir, out_dir=tmp_path / "opi", state="T1",
            charge=0, solvent=None, metals=[],
            main_basisset="def2-SVP", metal_basisset=None,
            config={"functional": "PBE0", "main_basisset": "def2-SVP",
                    "ESD_modus": "deltaSCF",
                    "ESD_T1_opt": "uks",
                    "ESD_frequency": "no"},
        )
