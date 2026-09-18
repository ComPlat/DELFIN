"""The ORCA inputs DELFIN writes today, built from fixtures, for checks that need no ORCA.

Two kinds of regression can only be seen on a whole input: one where an input
comes out different from what it used to (the first coordination sphere's
basis, a keyword that moved, a block that changed shape), and one where the
recovery has to take that input apart again.  Both were found on the archive
of real runs, months after the change; neither needs ORCA to run, only the
input DELFIN writes.

So this builds them: two systems (a small organic molecule and a metal
complex with a first coordination sphere) through every writer -- the
frequency job of a step, an OCCUPIER configuration, the excited states and
every rate job, the CONTROL overrides, the appended property jobs.  The
geometries are made up; nothing here comes from anybody's run.

The cores and the memory are fixed, so what comes out depends on the writers
alone.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict

FORMALDEHYDE = """C 0.000000 0.000000 -0.529000
O 0.000000 0.000000 0.677000
H 0.000000 0.937000 -1.117000
H 0.000000 -0.937000 -1.117000
"""

#: A cobalt with two DMF-like oxygens and a pyridine-like nitrogen: the metal
#: takes its own basis, the first sphere is a distance test.
#: The second oxygen sits between the cutoffs of scale 1.2 (2.30 A) and 1.3
#: (2.50 A), so which scale the writer uses is visible in the input.
COBALT_COMPLEX = """Co 0.000000 0.000000 0.000000
O 2.100000 0.000000 0.000000
O -2.400000 0.000000 0.000000
N 0.000000 2.050000 0.000000
C 3.050000 0.000000 0.000000
H 3.500000 0.950000 0.000000
"""

_STATE_GEOMETRY = ("4\nCoordinates from ORCA-job {job} E -114.290000000000\n"
                   "C 0.000000 0.000000 -0.530000\nO 0.000000 0.000000 0.680000\n"
                   "H 0.000000 0.940000 -1.120000\nH 0.000000 -0.940000 -1.120000\n")

CONTROL_OVERRIDES = """functional=PBE0
keyword:initial=[TightOpt]
additions:initial=[%scf
  MaxIter 300
end]
keyword:S1_S0_IC=[VeryTightSCF]
"""


def _config(**values) -> Dict[str, object]:
    from delfin.config import _load_template_defaults

    config = {
        **_load_template_defaults(),
        "functional": "PBE0",
        "main_basisset": "def2-SVP",
        "metal_basisset": "def2-TZVP",
        "aux_jk": "def2/J",
        "disp_corr": "D4",
        "ri_jkx": "RIJCOSX",
        "relativity": "none",
        "implicit_solvation_model": "CPCM",
        "solvent": "water",
        "charge": 0,
        "PAL": 8,
        "maxcore": 2000,
        "maxiter": 125,
        "maxiter_occupier": 125,
        "temperature": "298.15",
        "geom_opt": "OPT",
        "geom_opt_OCCUPIER": "OPT",
        "freq_type": "FREQ",
        "initial_guess": "PModel",
        "first_coordination_sphere_metal_basisset": "yes",
        "first_coordination_sphere_scale": "1.3",
        "approximate_spin_projection_APMethod": "2",
    }
    config.update(values)
    return config


def _xyz(path: Path, body: str, comment: str = "corpus") -> Path:
    path.write_text(f"{len(body.strip().splitlines())}\n{comment}\n{body}", encoding="utf-8")
    return path


def _steps(root: Path) -> Dict[str, str]:
    """A step's frequency job: the plain one, and a metal complex with broken symmetry."""
    from delfin.xyz_io import read_xyz_and_create_input3

    config = _config()
    _xyz(root / "initial.xyz", FORMALDEHYDE)
    _xyz(root / "red_step_1.xyz", COBALT_COMPLEX)
    read_xyz_and_create_input3(str(root / "initial.xyz"), str(root / "initial.inp"),
                               0, 1, "water", [], "def2-TZVP", "def2-SVP", config, "")
    read_xyz_and_create_input3(str(root / "red_step_1.xyz"), str(root / "red_step_1.inp"),
                               -1, 2, "water", ["Co"], "def2-TZVP", "def2-SVP", config,
                               '%moinp "input_red_step_1_OCCUPIER.gbw"\n%scf\n  BrokenSym 1,1\nend')
    return {p.name: p.read_text(encoding="utf-8") for p in (root / "initial.inp", root / "red_step_1.inp")}


def _occupier(root: Path) -> Dict[str, str]:
    """Two OCCUPIER configurations: the first one of a stage, and one started from it."""
    from delfin.occupier import read_and_modify_file_OCCUPIER

    config = _config()
    stage = root / "red_step_1_OCCUPIER"
    stage.mkdir(exist_ok=True)
    _xyz(stage / "input0.xyz", COBALT_COMPLEX)
    _xyz(stage / "input.xyz", COBALT_COMPLEX, comment="Coordinates from ORCA-job input E -1.0")
    read_and_modify_file_OCCUPIER(0, "input.inp", -1, 2, "water", ["Co"], "def2-TZVP", "def2-SVP",
                                  config, "", work_dir=stage)
    read_and_modify_file_OCCUPIER(1, "input3.inp", -1, 2, "water", ["Co"], "def2-TZVP", "def2-SVP",
                                  config, '%moinp "input.gbw"\n%scf\n  BrokenSym 3,2\n  APMethod 2\nend',
                                  work_dir=stage)
    return {f"{stage.name}/{p.name}": p.read_text(encoding="utf-8")
            for p in (stage / "input.inp", stage / "input3.inp")}


def _excited_states(root: Path) -> Dict[str, str]:
    """The excited states and every rate job the ESD module writes."""
    from delfin.esd_input_generator import (create_fluor_input, create_ic_input, create_isc_input,
                                            create_phosp_input, create_state_input)

    config = _config(ESD_modul="yes", ESD_modus="TDDFT", ESD_T1_opt="uks", ESD_frequency="yes",
                     states="[S1,T1]", ISCs="[S1>T1]", ICs="[S1>S0]", emission_rates="[f,p]",
                     ESD_NPOINTS="auto", ESD_MAXTIME="auto")
    (root / "start.txt").write_text(FORMALDEHYDE, encoding="utf-8")
    esd = root / "ESD"
    esd.mkdir(exist_ok=True)
    previous_cwd = Path.cwd()
    os.chdir(root)                      # the S0 job reads the run's own start.txt
    try:
        written = [create_state_input("S0", esd, 0, "water", [], "def2-SVP", "def2-TZVP", config)]
        for state, energy in (("S0", -114.29), ("S1", -114.15), ("T1", -114.20)):
            (esd / f"{state}.xyz").write_text(_STATE_GEOMETRY.format(job=state), encoding="utf-8")
            (esd / f"{state}.hess").write_text("$orca_hessian_file\n", encoding="utf-8")
            # the rate jobs take the energy gap from the states' outputs
            (esd / f"{state}.out").write_text(
                f"FINAL SINGLE POINT ENERGY {energy:.8f}\n"
                "                             ****ORCA TERMINATED NORMALLY****\n", encoding="utf-8")
        for state in ("S1", "T1"):
            written.append(create_state_input(state, esd, 0, "water", [], "def2-SVP", "def2-TZVP", config))
        written.append(create_isc_input("S1>T1", esd, 0, "water", [], "def2-SVP", "def2-TZVP",
                                        config, trootssl=0))
        written.append(create_ic_input("S1>S0", esd, 0, "water", [], "def2-SVP", "def2-TZVP", config))
        written.append(create_fluor_input(esd, 0, "water", [], "def2-SVP", "def2-TZVP", config))
        written.append(create_phosp_input(esd, 0, "water", [], "def2-SVP", "def2-TZVP", config))
    finally:
        os.chdir(previous_cwd)
    return {f"ESD/{Path(p).name}": Path(p).read_text(encoding="utf-8") for p in written}


def _appended_and_overridden(root: Path) -> Dict[str, str]:
    """A step with the appended jobs, and one the CONTROL overrides reached."""
    import shutil

    from delfin.esd_input_generator import (append_properties_of_interest_jobs,
                                            append_reorganisation_energy_jobs)
    from delfin.orca import _apply_control_overrides_to_input

    config = _config(properties_of_interest="IP,EA", calc_prop_of_interest="yes",
                     reorganisation_energy="lambda_p,lambda_m")
    with_jobs = root / "initial_with_properties.inp"
    shutil.copyfile(root / "initial.inp", with_jobs)
    append_properties_of_interest_jobs(inp_file=str(with_jobs), xyz_file="initial.xyz",
                                       base_charge=0, base_multiplicity=1, properties="IP,EA",
                                       config=config, solvent="water", metals=[],
                                       main_basisset="def2-SVP", metal_basisset="def2-TZVP")

    reorganisation = root / "ox_step_1_with_reorganisation.inp"
    shutil.copyfile(root / "red_step_1.inp", reorganisation)
    append_reorganisation_energy_jobs(inp_file=str(reorganisation), neutral_charge=0,
                                      neutral_multiplicity=1,
                                      reorganisation_energy="lambda_p,lambda_m", config=config,
                                      solvent="water", metals=["Co"], main_basisset="def2-SVP",
                                      metal_basisset="def2-TZVP", mode="lambda_p")

    overridden = root / "overridden"
    overridden.mkdir(exist_ok=True)
    (overridden / "CONTROL.txt").write_text(CONTROL_OVERRIDES, encoding="utf-8")
    shutil.copyfile(root / "initial.inp", overridden / "initial.inp")
    shutil.copyfile(root / "ESD" / "S1_S0_IC.inp", overridden / "S1_S0_IC.inp")
    for name in ("initial.inp", "S1_S0_IC.inp"):
        _apply_control_overrides_to_input(overridden / name, overridden)
    return {
        with_jobs.name: with_jobs.read_text(encoding="utf-8"),
        reorganisation.name: reorganisation.read_text(encoding="utf-8"),
        **{f"overridden/{name}": (overridden / name).read_text(encoding="utf-8")
           for name in ("initial.inp", "S1_S0_IC.inp")},
    }


def _other_modules(root: Path) -> Dict[str, str]:
    """The inputs the modules beside the main workflow write: stability constants and NMR."""
    from delfin.dashboard.tab_calculations_browser import build_calc_nmr_input
    from delfin.ensemble_nmr import build_orca_reference_input
    from delfin.stability_constant import build_orca_input

    written: Dict[str, str] = {}
    complexed = root / "solv_complex.inp"
    build_orca_input(_config(), _xyz(root / "solv_complex.xyz", COBALT_COMPLEX), complexed,
                     -1, 2, ["Co"], "water", broken_sym="", include_freq=True)
    written["solv_complex.inp"] = complexed.read_text(encoding="utf-8")
    written["nmr_of_a_structure.inp"] = build_calc_nmr_input(
        FORMALDEHYDE.strip().splitlines(), pal=8, maxcore=2000, solvent="chloroform")
    written["nmr_reference.inp"] = build_orca_reference_input(
        FORMALDEHYDE, solvent="chloroform", pal=8, maxcore=2000)
    return written


def build_corpus(root: Path) -> Dict[str, str]:
    """Write every input of the corpus under *root* and return them by name."""
    root.mkdir(parents=True, exist_ok=True)
    corpus: Dict[str, str] = {}
    corpus.update(_steps(root))
    corpus.update(_occupier(root))
    corpus.update(_excited_states(root))
    corpus.update(_appended_and_overridden(root))
    corpus.update(_other_modules(root))
    return corpus


__all__ = ["build_corpus", "COBALT_COMPLEX", "FORMALDEHYDE"]
