"""Workflow contract tests for non-CO2 DELFIN use cases.

Covers the remaining user-facing workflow shapes that the CO2
coordination tests don't — keeping parity with what real users run
on the cluster:

  - Classic ORCA pipeline (S0 / redox steps / ESD)
      Fritz-style: CONTROL.txt with method=classic + ESD_modul=yes
      + electrical properties + TDDFT configuration.

  - Hyperpol xTB screening
      Fast sTD-DFT-xTB hyperpolarizability workflow, CONTROL-driven.

  - TADF xTB screening
      Thermally activated delayed fluorescence screening via xtb_stda.

  - NMR post-processing
      delfin_NMR CLI that parses ORCA NMR output into the NMR report.

  - Raw ORCA wrapper mode (Tilman-style)
      `DELFIN_MODE=orca` — single .inp file, local_runner finds it,
      runs ORCA, local_runner._print_orca_summary reports SUCCESS/FAILED.
      This is the path that the set +e / sync_results_back fixes
      protect. We test the Python-side contracts.

All tests are pure-Python — no ORCA, no xtb, no xtb_stda ever runs.
Every check is either config-parsing, introspection, or signature
verification.
"""

from __future__ import annotations

import dataclasses
import importlib
import inspect
from pathlib import Path

import pytest


# ===========================================================================
# 1. Classic ORCA pipeline (Fritz-style, method=classic + ESD + elprop)
# ===========================================================================

_CLASSIC_CONTROL_TEMPLATE = """\
input_file=input.txt
NAME=classic_test
SMILES=c1ccccc1
charge=0
------------------------------------
Solvation:
implicit_solvation_model=CPCM
solvent=chcl3
XTB_SOLVATOR=no
number_explicit_solv_molecules=2
------------------------------------
Global geometry optimisation:
xTB_method=XTB2
smiles_converter=NORMAL
XTB_OPT=no
XTB_GOAT=no
CREST=no
multiplicity_global_opt=
------------------------------------
IMAG=yes
IMAG_scope=initial
IMAG_option=2
allow_imaginary_freq=0
IMAG_sp_energy_window=1e-3
IMAG_optimize_candidates=no
------------------------------------
calc_prop_of_interest=no
properties_of_interest=IP,EA
------------------------------------
Redox steps:
calc_initial=no
oxidation_steps=
reduction_steps=
method=classic
calc_potential_method=2
------------------------------------
ESD module (excited state dynamics):
ESD_modul=yes
ESD_modus=TDDFT
ESD_frequency=no
states=S1
ISCs=
ICs=
emission_rates=
phosp_IROOT=1,2,3
TROOTSSL=-1,0,1
DOHT=TRUE
ESD_LINES=LORENTZ
ESD_LINEW=50
ESD_INLINEW=250
ESD_NPOINTS=131072
ESD_MAXTIME=12000
hybrid1_geom_MaxIter=60
--------------------
Electrical Properties:
elprop_Dipole=yes
elprop_Quadrupole=yes
elprop_Hyperpol=yes
elprop_Polar=yes
elprop_PolarVelocity=no
elprop_PolarDipQuad=no
elprop_PolarQuadQuad=no
--------------------
TDDFT Settings:
TDDFT_TDDFT_maxiter=500
TDDFT_nroots=15
TDDFT_maxdim=30
TDDFT_TDA=FALSE
TDDFT_followiroot=true
TDDFT_SOC=false
------------------------------------
Level of Theory:
functional=PBE0
disp_corr=D4
ri_jkx=RIJCOSX
relativity=ZORA
aux_jk=def2/J
aux_jk_rel=SARC/J
main_basisset=def2-SVP
main_basisset_rel=ZORA-def2-SVP
metal_basisset=def2-TZVP
metal_basisset_rel=SARC-ZORA-TZVP
first_coordination_sphere_metal_basisset=no
first_coordination_sphere_scale=1.3
geom_opt=OPT TIGHTSCF defgrid3
freq_type=FREQ
initial_guess=PModel
temperature=298.15
maxiter=125
qmmm_option=QM/PBEH-3c
------------------------------------
Resource Settings:
PAL=12
maxcore=9000
parallel_workflows=yes
pal_jobs=2
orca_parallel_strategy=auto
------------------------------------
Prints:
print_MOs=yes
print_Loewdin_population_analysis=no
"""


@pytest.fixture
def classic_control_path(tmp_path):
    p = tmp_path / "CONTROL.txt"
    p.write_text(_CLASSIC_CONTROL_TEMPLATE)
    return p


@pytest.fixture
def classic_config(classic_control_path):
    from delfin.config import OCCUPIER_parser
    return OCCUPIER_parser(str(classic_control_path))


def test_classic_config_selects_classic_method(classic_config):
    """method=classic is the pipeline gate that routes to
    run_classic_phase (instead of run_occuper_phase)."""
    assert classic_config.get("method") == "classic", (
        "Classic ORCA pipeline depends on method='classic' — if the parser "
        "downcases this wrongly or the gate name shifts, the pipeline "
        "silently takes the OCCUPIER branch instead."
    )


def test_classic_config_esd_module_enabled(classic_config):
    """Fritz's CONTROL.txt sets ESD_modul=yes → run_esd_phase must run."""
    assert classic_config.get("ESD_modul") == "yes"
    # Parser lower-cases the ESD_modus value
    assert str(classic_config.get("ESD_modus")).lower() == "tddft"
    # states parses into a Python list of state labels (case may vary)
    states = classic_config.get("states")
    assert states is not None, "states key dropped by parser"
    # Convert to normalized form for comparison
    if isinstance(states, list):
        norm = [str(s).upper() for s in states]
    else:
        norm = [str(states).upper()]
    assert "S1" in norm, (
        f"states must include 'S1', got {states!r}"
    )


def test_classic_config_electrical_properties_flags(classic_config):
    """elprop_Hyperpol, elprop_Polar etc. drive ORCA !Keywords and
    %elprop blocks. The yes/no string survival is load-bearing."""
    assert classic_config.get("elprop_Dipole") == "yes"
    assert classic_config.get("elprop_Quadrupole") == "yes"
    assert classic_config.get("elprop_Hyperpol") == "yes"
    assert classic_config.get("elprop_Polar") == "yes"
    # Unused ones stay off
    assert classic_config.get("elprop_PolarVelocity") == "no"


def test_classic_config_tddft_block_preserved(classic_config):
    """TDDFT is injected into the ORCA input as %tddft block. Every
    setting here maps to a kwarg on that block."""
    assert classic_config.get("TDDFT_nroots") == 15
    assert classic_config.get("TDDFT_maxdim") == 30
    assert classic_config.get("TDDFT_TDDFT_maxiter") == 500
    # Boolean-ish keys survive as strings (case-sensitive ORCA input)
    assert str(classic_config.get("TDDFT_TDA")).upper() == "FALSE"
    assert classic_config.get("TDDFT_followiroot") == "true"


def test_classic_config_calc_initial_disabled(classic_config):
    """Fritz sets calc_initial=no → run_occuper_phase must skip."""
    assert classic_config.get("calc_initial") == "no"


def test_classic_config_geom_opt_keyword_composite(classic_config):
    """Fritz uses a composite geom_opt = 'OPT TIGHTSCF defgrid3' — this
    is inlined into ORCA's ! line. Whitespace / capitalisation matters."""
    geom_opt = classic_config.get("geom_opt", "")
    assert "OPT" in geom_opt, f"OPT keyword dropped from geom_opt: {geom_opt!r}"
    assert "TIGHTSCF" in geom_opt, (
        f"TIGHTSCF dropped from geom_opt: {geom_opt!r} — convergence tightening "
        f"silently weakened for every user."
    )


def test_classic_pipeline_run_esd_phase_entry_point():
    from delfin.workflows.pipeline import run_esd_phase
    assert callable(run_esd_phase)


def test_classic_pipeline_run_classic_phase_entry_point():
    from delfin.workflows.pipeline import run_classic_phase
    assert callable(run_classic_phase)


# ===========================================================================
# 2. Hyperpolarizability xTB workflow
# ===========================================================================

_HYPERPOL_CONTROL_TEMPLATE = """\
input_file=input.txt
NAME=hyperpol_test
SMILES=c1ccc2[nH]ccc2c1
charge=0
method=classic
calc_initial=no
oxidation_steps=
reduction_steps=
smiles_converter=NORMAL
implicit_solvation_model=CPCM
solvent=water
XTB_SOLVATOR=no
number_explicit_solv_molecules=0
IMAG=no
------------------------------------
xTB Hyperpolarizability (sTD-DFT-xTB):
hyperpol_xTB=yes
hyperpol_xTB_xyz=start.txt
hyperpol_xTB_preopt=xtb
hyperpol_xTB_engine=std2
hyperpol_xTB_bfw=no
hyperpol_xTB_wavelengths=532,1064
hyperpol_xTB_energy_window=15.0
------------------------------------
Level of Theory:
functional=PBE0
main_basisset=def2-SVP
metal_basisset=def2-TZVP
ri_jkx=RIJCOSX
relativity=
disp_corr=D4
aux_jk=def2/J
aux_jk_rel=
main_basisset_rel=
metal_basisset_rel=
first_coordination_sphere_metal_basisset=no
first_coordination_sphere_scale=1.0
geom_opt=OPT
freq_type=FREQ
initial_guess=PModel
temperature=298.15
maxiter=125
qmmm_option=
Resource Settings:
PAL=8
maxcore=2000
parallel_workflows=yes
pal_jobs=2
"""


@pytest.fixture
def hyperpol_config(tmp_path):
    from delfin.config import OCCUPIER_parser
    p = tmp_path / "CONTROL.txt"
    p.write_text(_HYPERPOL_CONTROL_TEMPLATE)
    return OCCUPIER_parser(str(p))


def test_hyperpol_config_toggle_preserved(hyperpol_config):
    assert hyperpol_config.get("hyperpol_xTB") == "yes"
    assert hyperpol_config.get("hyperpol_xTB_engine") == "std2"


def test_hyperpol_config_wavelengths_preserved(hyperpol_config):
    """User enters comma-separated wavelengths; parser needs to hand
    them to the sTD-DFT-xTB runner unchanged."""
    wl = hyperpol_config.get("hyperpol_xTB_wavelengths")
    assert wl is not None and "532" in str(wl) and "1064" in str(wl), (
        f"hyperpol_xTB_wavelengths lost values: {wl!r}"
    )


def test_hyperpol_config_energy_window_is_float(hyperpol_config):
    """Energy window feeds a numerical comparison. Type flip breaks it."""
    ew = hyperpol_config.get("hyperpol_xTB_energy_window")
    # Should survive as numeric string or float
    assert str(ew).strip() in ("15.0", "15"), (
        f"hyperpol_xTB_energy_window reshaped: {ew!r}"
    )


def test_hyperpol_module_public_api():
    """The hyperpol entry points the pipeline orchestrator imports must
    stay on the module surface."""
    import delfin.hyperpol as h
    assert hasattr(h, "HyperpolResult"), "HyperpolResult dataclass is gone"


def test_hyperpol_result_fields():
    """HyperpolResult is the data contract between hyperpol_xtb and
    reporting. Losing fields breaks report generation downstream."""
    from delfin.hyperpol import HyperpolResult

    fields = {f.name for f in dataclasses.fields(HyperpolResult)}
    expected = {
        "label", "smiles", "charge", "multiplicity", "engine", "preopt",
        "energy_window_ev", "requested_wavelengths_nm", "workdir",
        "initial_xyz",
    }
    missing = expected - fields
    assert not missing, (
        f"HyperpolResult dropped fields: {sorted(missing)}. Reports will "
        f"crash with AttributeError when they access the missing names."
    )


def test_hyperpol_pipeline_phase_entry_point():
    from delfin.workflows.pipeline import run_hyperpol_xtb_phase
    assert callable(run_hyperpol_xtb_phase)


# ===========================================================================
# 3. TADF xTB workflow
# ===========================================================================

_TADF_CONTROL_TEMPLATE = """\
input_file=input.txt
NAME=tadf_test
SMILES=c1ccc2c(c1)ccc2
charge=0
method=classic
calc_initial=no
oxidation_steps=
reduction_steps=
smiles_converter=NORMAL
implicit_solvation_model=
solvent=
XTB_SOLVATOR=no
number_explicit_solv_molecules=0
IMAG=no
------------------------------------
xTB TADF Screening:
tadf_xTB=yes
tadf_xTB_xyz=start.txt
tadf_xTB_preopt=xtb
tadf_xTB_excited_method=stda
tadf_xTB_bfw=no
tadf_xTB_energy_window=10.0
tadf_xTB_run_t1_opt=yes
------------------------------------
Level of Theory:
functional=PBE0
main_basisset=def2-SVP
metal_basisset=def2-TZVP
ri_jkx=RIJCOSX
relativity=
disp_corr=D4
aux_jk=def2/J
aux_jk_rel=
main_basisset_rel=
metal_basisset_rel=
first_coordination_sphere_metal_basisset=no
first_coordination_sphere_scale=1.0
geom_opt=OPT
freq_type=FREQ
initial_guess=PModel
temperature=298.15
maxiter=125
qmmm_option=
Resource Settings:
PAL=8
maxcore=2000
parallel_workflows=yes
pal_jobs=2
"""


@pytest.fixture
def tadf_config(tmp_path):
    from delfin.config import OCCUPIER_parser
    p = tmp_path / "CONTROL.txt"
    p.write_text(_TADF_CONTROL_TEMPLATE)
    return OCCUPIER_parser(str(p))


def test_tadf_config_toggle_preserved(tadf_config):
    assert tadf_config.get("tadf_xTB") == "yes"
    assert tadf_config.get("tadf_xTB_excited_method") == "stda"


def test_tadf_config_t1_opt_flag(tadf_config):
    """tadf_xTB_run_t1_opt=yes triggers a full T1 geometry optimisation
    on top of the sTDA screen. Losing it silently gives vertical-only."""
    assert tadf_config.get("tadf_xTB_run_t1_opt") == "yes"


def test_tadf_result_fields():
    from delfin.tadf_xtb import TadfXtbResult

    fields = {f.name for f in dataclasses.fields(TadfXtbResult)}
    expected = {
        "label", "smiles", "charge", "multiplicity", "xtb_method",
        "excited_method", "workdir", "initial_xyz", "s0_xyz",
    }
    missing = expected - fields
    assert not missing, f"TadfXtbResult dropped fields: {sorted(missing)}"


def test_tadf_pipeline_phase_entry_point():
    from delfin.workflows.pipeline import run_tadf_xtb_phase
    assert callable(run_tadf_xtb_phase)


# ===========================================================================
# 4. NMR post-processing (delfin_NMR CLI)
# ===========================================================================


def test_nmr_cli_is_importable_and_callable():
    """delfin_NMR is the post-processing CLI that turns ORCA EPRNMR
    output into the NMR report. Entry-point break = silent NMR loss."""
    from delfin.cli_nmr_report import main
    assert callable(main), "delfin.cli_nmr_report.main must be callable"


def test_nmr_cli_exposes_documented_help():
    """`delfin_NMR --help` must still describe the tool, so users who
    typed the wrong args see something instead of a Python traceback."""
    from delfin.cli_nmr_report import main
    assert main.__doc__, "cli_nmr_report.main lost its docstring — --help "\
        "output would be empty for users."


def test_nmr_report_generator_module_importable():
    """The ORCA-NMR-parsing helpers live in delfin.reporting.nmr_report.
    If they stop importing, the CLI can't produce a report."""
    mod = importlib.import_module("delfin.reporting.nmr_report")
    # The module should expose at least one public symbol
    public = [n for n in dir(mod) if not n.startswith("_")]
    assert public, "nmr_report module has no public symbols — broken refactor"


# ===========================================================================
# 5. Raw ORCA wrapper mode (Tilman-style DELFIN_MODE=orca)
# ===========================================================================
#
# Tilman's use-case: drop a single .inp file into a directory, set
# DELFIN_MODE=orca, run the SLURM wrapper. local_runner has a special
# code branch for this that picks up the .inp file, runs ORCA, prints
# a summary. This is the path that the set +e / sync_results_back
# fixes protect (commit 437106f). We test the Python-side invariants.


def test_local_runner_exposes_main_entry():
    """submit_delfin.sh launches `python -m delfin.dashboard.local_runner`
    — that invocation requires `main()` to exist at module level."""
    from delfin.dashboard.local_runner import main
    assert callable(main)


def test_local_runner_configure_environment_still_exists():
    """_configure_environment sets OpenMPI / MKL / OMP_NUM_THREADS env
    vars before ORCA runs. If the name changes, SLURM ORCA jobs run
    with wrong MPI settings (hangs / crashes on BwUniCluster)."""
    from delfin.dashboard.local_runner import _configure_environment
    assert callable(_configure_environment)


def test_local_runner_print_orca_summary_signature():
    """_print_orca_summary(out_file, ok, elapsed) is called from the
    orca-mode branch of _run_mode. The three args are positional —
    a rename of any of them silently mis-prints the summary."""
    from delfin.dashboard.local_runner import _print_orca_summary

    sig = inspect.signature(_print_orca_summary)
    params = list(sig.parameters.keys())
    assert params[:3] == ["out_file", "ok", "elapsed"], (
        f"_print_orca_summary signature changed from (out_file, ok, elapsed) "
        f"to {params[:3]}. The orca-mode branch at local_runner._run_mode "
        f"passes these positionally."
    )


def test_local_runner_detects_orca_terminated_marker(tmp_path):
    """_print_orca_summary reads the last 8KB of the .out file and
    looks for the literal string 'ORCA TERMINATED NORMALLY'. If that
    string match gets rewritten (e.g. changed to a regex or case fold),
    the SUCCESS / NOT FOUND reporting flips. Smoke-test the happy path.
    """
    from delfin.dashboard.local_runner import _print_orca_summary

    out = tmp_path / "fake.out"
    out.write_text(
        "Some ORCA banner\n"
        "Lots of calculation data\n"
        "                              *** ORCA TERMINATED NORMALLY ***\n"
        "TOTAL RUN TIME: 0 days 0 hours 1 minutes 23 seconds\n"
    )
    # Must not raise — prints "TERMINATED NORMALLY" to stdout
    _print_orca_summary(str(out), ok=True, elapsed=83.0)


def test_local_runner_orca_mode_uses_env_var_for_inp_file():
    """Orca mode in _run_mode reads DELFIN_INP_FILE from the env. If
    the env-var name is renamed (e.g. to ORCA_INPUT_FILE), the SLURM
    wrapper's `export DELFIN_INP_FILE=...` becomes a no-op and
    local_runner falls back to globbing *.inp."""
    import delfin.dashboard.local_runner as lr

    source = Path(inspect.getfile(lr)).read_text()
    assert "DELFIN_INP_FILE" in source, (
        "DELFIN_INP_FILE is no longer referenced in local_runner. The "
        "submit script exports that exact name — losing it here turns "
        "raw-ORCA-mode into a glob-based fallback and mis-selects the "
        "input file (exactly the issue that hit Tilman once)."
    )


# ===========================================================================
# 6. Solvation model registry — CPCM / SMD / ALPB propagation
# ===========================================================================


def test_cpcm_solvent_config_preserved():
    """CPCM is the most-used solvation model across all workflows. A
    silent drop of the config key cascades to every downstream ORCA
    input (no !CPCM line → implicit gas-phase → wrong energies)."""
    from delfin.config import OCCUPIER_parser
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "CONTROL.txt"
        p.write_text(_CLASSIC_CONTROL_TEMPLATE)
        cfg = OCCUPIER_parser(str(p))

    assert cfg.get("implicit_solvation_model") == "CPCM"
    assert cfg.get("solvent") == "chcl3"


# ===========================================================================
# 7. ESD input generator — entry points behind ESD_modul=yes
# ===========================================================================


def test_esd_input_generator_top_level_api():
    """When ESD_modul=yes, run_esd_phase calls into esd_input_generator
    for property jobs (IP/EA) and reorganization jobs. Pin the names."""
    from delfin.esd_input_generator import (
        append_properties_of_interest_jobs,
        append_reorganisation_energy_jobs,
        _parse_step_set,
    )
    assert callable(append_properties_of_interest_jobs)
    assert callable(append_reorganisation_energy_jobs)
    assert callable(_parse_step_set)


# ===========================================================================
# 8. Stability Constant workflow (sc_*)
# ===========================================================================


def test_stability_constant_module_importable():
    """stability_constant is the thermodynamics workflow used when
    stability_constant=yes in CONTROL.txt. Optional but
    user-facing — refactor-guard."""
    mod = importlib.import_module("delfin.stability_constant")
    # Look for the public entry — usually `run` or similar
    public = [n for n in dir(mod) if not n.startswith("_")]
    assert public, "stability_constant module has no public API exposed"
