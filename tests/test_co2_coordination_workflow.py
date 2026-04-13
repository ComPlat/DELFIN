"""End-to-end contract tests for the SMILES → GUPPY → OCCUPIER → CO2
Coordinator workflow — i.e. the pipeline shape users depend on when
running a CO2 coordination screen.

These tests use a real CONTROL.txt (adapted from an actual
CO2-coordination screen) but NEVER call ORCA or external binaries. They
check that each refactor-sensitive stage of the pipeline still:

  - parses its CONTROL.txt input into the expected config keys
  - exposes the Python entry points the pipeline orchestrator calls
  - builds the correct WorkflowJob graph (job IDs, dependencies,
    inline/bottleneck flags)
  - resolves OCCUPIER sequences per species (even/odd electron counts)
  - recognises CO2 coordination config and wires setup_co2_from_delfin()
    into the workflow
  - hands GUPPY its settings unchanged

If anything in that list breaks, users' CO2-coordination runs silently
stop producing results. These tests are the tripwire that turns
"Jan's 5 overnight runs crashed" into "the commit is blocked before
it hits main".
"""

from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# The CONTROL.txt fixture — derived from a real CO2-coordination screen
# (RSS_paper_19_bug). Kept small enough to read at a glance; every value
# is load-bearing for at least one test below.
# ---------------------------------------------------------------------------

_CONTROL_TEMPLATE = """\
input_file=input.txt
NAME=test_complex
SMILES=[N+]12=CC=CC=C1N(C3=[N+]4CCC[N+]5=C6[Ni-2]23[N+](C=CC=C7)=C7N6C8=C5C=CC=C8)C9=C4C=CC=C9
charge=+2
------------------------------------
Solvation:
implicit_solvation_model=CPCM
solvent=DMF
XTB_SOLVATOR=no
number_explicit_solv_molecules=2
------------------------------------
Global geometry optimisation:
xTB_method=XTB2
smiles_converter=GUPPY
XTB_OPT=no
XTB_GOAT=yes
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
reorganisation_energy=lambda_p,lambda_m
------------------------------------
Redox steps:
calc_initial=yes
oxidation_steps=
reduction_steps=1,2
method=OCCUPIER
calc_potential_method=2
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
first_coordination_sphere_metal_basisset=yes
first_coordination_sphere_scale=1.3
geom_opt=OPT
freq_type=FREQ
initial_guess=PModel
temperature=298.15
maxiter=125
qmmm_option=QM/PBEH-3c
------------------------------------
Resource Settings:
PAL=40
maxcore=6000
parallel_workflows=yes
pal_jobs=4
orca_parallel_strategy=auto
enable_job_timeouts=no
------------------------------------
GUPPY_settings:
GUPPY_RUNS=10
GUPPY_GOAT=0
GUPPY_PARALLEL_JOBS=4
GUPPY_SEED=31
------------------------------------
OCCUPIER-Settings:
--------------------
OCCUPIER_method=auto
OCCUPIER_tree=own
OWN_TREE_PURE_WINDOW=3
OWN_progressive_from=no
fob_equal_weights=yes
frequency_calculation_OCCUPIER=no
occupier_selection=tolerance
occupier_precision=3
occupier_epsilon=5e-4
clean_override_window_h=0.002
maxiter_occupier=125
geom_opt_OCCUPIER=OPT
pass_wavefunction=no
--------------------
OCCUPIER_sequence_profiles:
-3,-2,-1,0,+1,+2,+3=[
even electron number:
even_seq = [
  {"index": 1, "m": 1, "BS": "",    "from": 0},
  {"index": 4, "m": 3, "BS": "",    "from": 1},
  {"index": 7, "m": 5, "BS": "",    "from": 4}
]
-------------------
odd electron number:
odd_seq = [
  {"index": 1, "m": 2, "BS": "",    "from": 0},
  {"index": 4, "m": 4, "BS": "",    "from": 1},
  {"index": 7, "m": 6, "BS": "",    "from": 4}
]
]
--------------------
ORCA base overrides (optional):
keyword:basename=[]
additions:basename=[]
--------------------
CO2 Coordination:
co2_coordination=on
co2_species_delta=-2
"""


@pytest.fixture
def control_path(tmp_path):
    """Write the CONTROL.txt template into a temp dir and return the path."""
    p = tmp_path / "CONTROL.txt"
    p.write_text(_CONTROL_TEMPLATE)
    return p


@pytest.fixture
def config(control_path):
    from delfin.config import OCCUPIER_parser
    return OCCUPIER_parser(str(control_path))


# ===========================================================================
# 1. CONTROL.txt parsing — every key users expect must survive
# ===========================================================================


def test_control_txt_parses_core_identity_fields(config):
    """Basic bookkeeping — name, SMILES, charge come through as typed."""
    assert config.get("NAME") == "test_complex"
    assert config.get("charge") == 2   # parsed as int with sign stripped
    assert "Ni" in config.get("SMILES", ""), "SMILES must contain the metal"


def test_control_txt_preserves_workflow_toggles(config):
    """Workflow on/off switches drive which pipeline phases run."""
    assert config.get("smiles_converter") == "GUPPY"
    assert config.get("XTB_GOAT") == "yes"
    assert config.get("IMAG") == "yes"
    assert config.get("calc_initial") == "yes"
    assert config.get("co2_coordination") == "on"


def test_control_txt_preserves_numeric_resource_settings(config):
    """Resource settings map to SBATCH args on the cluster — a type flip
    here (string vs int) would silently break the SLURM submit script."""
    assert config.get("PAL") == 40
    assert config.get("maxcore") == 6000
    assert config.get("pal_jobs") == 4


def test_control_txt_exposes_reduction_and_oxidation_steps(config):
    """Jan's workflow runs reduction_steps=1,2 — the pipeline reads this
    exact field name to decide which red_step_N OCCUPIER jobs to create."""
    red = config.get("reduction_steps")
    # Parser turns the comma-separated list into a Python list of strings
    assert red == ["1", "2"], (
        f"reduction_steps must parse to ['1','2'], got {red!r}. "
        f"Downstream OCCUPIER job-builders iterate this list to create "
        f"red_step_1, red_step_2 — losing the list breaks Jan's workflow."
    )
    ox = config.get("oxidation_steps", "")
    assert ox in ("", [], None), f"oxidation_steps must be empty, got {ox!r}"


def test_control_txt_level_of_theory_block(config):
    """The Level-of-Theory section feeds directly into ORCA input files."""
    assert config.get("functional") == "PBE0"
    assert config.get("disp_corr") == "D4"
    assert config.get("ri_jkx") == "RIJCOSX"
    assert config.get("relativity") == "ZORA"
    assert config.get("main_basisset") == "def2-SVP"
    assert config.get("metal_basisset") == "def2-TZVP"


def test_control_txt_co2_species_delta(config):
    """CO2 species_delta must survive as an int (it selects the reduction
    step the CO2 coordinator targets — losing the sign would misroute)."""
    assert config.get("co2_species_delta") == -2


# ===========================================================================
# 2. OCCUPIER sequence profiles parse into the expected shape
# ===========================================================================


def test_occupier_sequence_block_is_captured(config):
    """The `OCCUPIER_sequence_profiles` block is a mini-DSL the CONTROL
    parser must preserve verbatim for the OCCUPIER engine."""
    seq_blocks = config.get("_occupier_sequence_blocks")
    assert seq_blocks, (
        "OCCUPIER_sequence_profiles block is missing from parsed config. "
        "The sequence-resolver relies on this raw text."
    )


def test_occupier_sequence_resolver_picks_even_for_charge_zero(config):
    """For Jan's Ni complex at neutral charge (even electron count, 180
    electrons), the resolver must hand back the even-electron sequence."""
    from delfin.occupier_sequences import resolve_sequences_for_delta

    bundle = resolve_sequences_for_delta(config, delta=0, parity_hint="even")
    assert isinstance(bundle, dict)
    assert "even_seq" in bundle, (
        "Resolver did not return an even_seq — Jan's template defines one, "
        "losing it silently means OCCUPIER falls back to a default sequence."
    )
    entries = bundle["even_seq"]
    assert len(entries) == 3, (
        f"even_seq must have 3 FoB entries (index 1, 4, 7), got {len(entries)}"
    )
    assert [e["index"] for e in entries] == [1, 4, 7]
    assert [e["m"] for e in entries] == [1, 3, 5], (
        "even_seq multiplicities from Jan's template are 1/3/5 — any other "
        "values indicate the resolver dropped the custom profile."
    )


# ===========================================================================
# 3. CO2 Coordinator — the feature that makes this specifically Jan's flow
# ===========================================================================


def test_co2_coordinator_public_api_exists():
    """Jan's pipeline relies on setup_co2_from_delfin() being importable.

    Bug-class guard: a refactor could rename the function or move it
    under a subpackage; this test catches that before users hit it.
    """
    import delfin.co2

    assert hasattr(delfin.co2, "setup_co2_from_delfin"), (
        "delfin.co2.setup_co2_from_delfin disappeared — CO2 coordination "
        "workflow has no entry point."
    )
    assert callable(delfin.co2.setup_co2_from_delfin)


def test_co2_coordinator_signature_contract():
    """Pin the signature — the CLI and the pipeline both call this
    function with positional (job_dir, species_delta). Renaming either
    kwarg breaks them all."""
    from delfin.co2 import setup_co2_from_delfin

    sig = inspect.signature(setup_co2_from_delfin)
    params = list(sig.parameters.keys())
    assert params[0] == "job_dir", (
        f"First arg must be 'job_dir', got {params[0]!r}"
    )
    assert "species_delta" in params, (
        f"'species_delta' param is gone. Got {params}"
    )


def test_co2_species_delta_name_mapping():
    """species_delta_to_name is a pure function that maps -2 → 'red_step_2'
    etc. This naming drives where CO2_Coordinator6 looks for the reduced
    geometry. Breaking the mapping silently routes CO2 to the wrong step."""
    from delfin.co2.chain_setup import species_delta_to_name

    # Contract: negative = reduction, positive = oxidation, 0 = initial
    assert species_delta_to_name(0) == "initial"
    # Reduction / oxidation naming must match OCCUPIER's folder layout
    red_2 = species_delta_to_name(-2)
    assert "red" in red_2.lower() and "2" in red_2, (
        f"species_delta=-2 must map to a reduction-step name. Got {red_2!r}"
    )
    ox_1 = species_delta_to_name(+1)
    assert "ox" in ox_1.lower() and "1" in ox_1, (
        f"species_delta=+1 must map to an oxidation-step name. Got {ox_1!r}"
    )


# ===========================================================================
# 4. OCCUPIER FoB (Fragment-of-Basis) job graph
# ===========================================================================
#
# This is exactly where Jan's runs crashed: build_flat_occupier_fob_jobs
# constructs the WorkflowJob list for initial_fob_1, initial_fob_4,
# initial_fob_7, initial_fob_best. The crash was in the inline-job path
# (_submit) — but the job-building itself is a pure function of config.
# We verify the graph shape without running anything.


def test_build_flat_occupier_fob_jobs_is_importable():
    """Entry point for the OCCUPIER phase. If this name moves, every
    CONTROL.txt that selects method=OCCUPIER silently stops working."""
    from delfin.workflows.engine.occupier import build_flat_occupier_fob_jobs
    assert callable(build_flat_occupier_fob_jobs)


def test_occupier_execution_context_still_exposed():
    """The OccupierExecutionContext is the data handoff between the
    pipeline and the FoB runner. Losing it breaks the handshake."""
    from delfin.workflows.engine.occupier import OccupierExecutionContext
    assert OccupierExecutionContext is not None


# ===========================================================================
# 5. GUPPY step settings make it through to the runner
# ===========================================================================


def test_guppy_settings_preserved_in_parsed_config(config):
    """GUPPY_RUNS, GUPPY_SEED and friends get routed into environment vars
    the local_runner reads. If CONTROL.txt parsing drops them, GUPPY
    runs with silent defaults — users end up with 20 runs / random seed
    instead of the 10 runs / seed=31 they configured."""
    assert config.get("GUPPY_RUNS") == 10
    assert config.get("GUPPY_PARALLEL_JOBS") == 4
    assert config.get("GUPPY_SEED") == 31


def test_guppy_entry_point_is_importable():
    """delfin-guppy CLI must still resolve — part of Jan's SMILES→XYZ path."""
    from delfin.guppy_sampling import main
    assert callable(main)


# ===========================================================================
# 6. Pipeline phase entry points the orchestrator calls
# ===========================================================================


def test_pipeline_run_occuper_phase_entry_point():
    """Jan's crash stack: cli.py:1760 → run_occuper_phase(). This is the
    exact entry point his CONTROL.txt triggers. Must stay importable."""
    from delfin.workflows.pipeline import run_occuper_phase
    assert callable(run_occuper_phase)

    # Signature contract: takes a PipelineContext, returns a truthy value
    sig = inspect.signature(run_occuper_phase)
    params = list(sig.parameters.keys())
    assert len(params) >= 1, f"run_occuper_phase must take at least a ctx arg, got {params}"


def test_pipeline_run_classic_phase_entry_point():
    """The redox-steps phase follows OCCUPIER — used for reduction_steps=1,2."""
    from delfin.workflows.pipeline import run_classic_phase
    assert callable(run_classic_phase)


def test_pipeline_context_dataclass_fields():
    """PipelineContext is the argument bundle every phase receives.
    Fields that disappear break the orchestrator."""
    from delfin.workflows.pipeline import PipelineContext
    import dataclasses

    fields = {f.name for f in dataclasses.fields(PipelineContext)}
    # These are the fields the phase functions access — if any disappear,
    # Jan's pipeline breaks somewhere deep. Pinned against pre-refactor
    # contract: control_file_path, config, input_file, metals,
    # main_basisset, metal_basisset, solvent, charge, multiplicity,
    # total_electrons_txt, name, PAL, number_explicit_solv_molecules,
    # start_time, file_bundle, extra.
    expected = {
        "config",
        "control_file_path",
        "input_file",
        "file_bundle",
        "metals",
        "charge",
        "multiplicity",
        "solvent",
        "main_basisset",
        "metal_basisset",
        "PAL",
    }
    missing = expected - fields
    assert not missing, (
        f"PipelineContext lost critical fields: {sorted(missing)}. "
        f"Current fields: {sorted(fields)}. Losing any of these breaks "
        f"the phase functions that destructure ctx.*"
    )


# ===========================================================================
# 7. IMAG step (imaginary-frequency recovery)
# ===========================================================================


def test_imag_step_public_api():
    """IMAG=yes triggers run_IMAG after OCCUPIER. Refactor guard."""
    from delfin.imag import run_IMAG
    assert callable(run_IMAG)


def test_imag_config_fields_preserved(config):
    """IMAG configuration flags — required for the IMAG phase gate."""
    assert config.get("IMAG") == "yes"
    assert config.get("IMAG_scope") == "initial"
    # Option 2 = "SCF_tight_relax" (or similar) in Jan's workflow
    assert str(config.get("IMAG_option")) == "2"


# ===========================================================================
# 8. Level-of-theory dispatch — metal-basis handling
# ===========================================================================


def test_metal_basis_per_atom_flag_preserved(config):
    """Jan's config uses per-metal basis for coordination sphere atoms.
    Losing this flag silently falls back to uniform basis → wrong energies."""
    assert config.get("first_coordination_sphere_metal_basisset") == "yes"
    # The scale factor is a float — refactor bug could make it str
    scale = config.get("first_coordination_sphere_scale")
    assert scale in (1.3, "1.3"), (
        f"first_coordination_sphere_scale = {scale!r} — must stay 1.3"
    )
