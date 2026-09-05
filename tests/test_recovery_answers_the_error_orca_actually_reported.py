"""Recovery has to answer the failure ORCA reported, not a plausible one.

Every fixture below is the wording ORCA 6.1.1 produced when the failure was
provoked for real — a run capped at ``%maxcore 1``, an SCF held to two
iterations, an optimisation held to one cycle.

Two findings came out of that and are pinned here:

* Out of memory was classified as an MPI crash. ORCA routes the abort through
  LEANSCF, the generic "error termination in LEANSCF" rule outranked the
  memory rule, and the retry answered a memory problem with MPI fixes. ORCA
  had meanwhile printed the figure it needed.
* An optimisation that runs out of cycles still ends with "ORCA TERMINATED
  NORMALLY". Judged on that marker alone it counts as finished, so recalc
  skips it on every resubmission — it can never finish — and dependent jobs
  continue from a geometry that is not a stationary point.
"""

from pathlib import Path

import pytest

from delfin import smart_recalc
from delfin.orca_recovery import (
    OrcaErrorDetector,
    OrcaErrorType,
    OrcaInputModifier,
    RecoveryStrategy,
)


# ORCA 6.1.1, benzene/B3LYP/def2-TZVP run with %maxcore 1.
_OUT_OF_MEMORY = """
                              memory conserving SCF solver
 Error  (ORCA_SCF): Not enough memory available!
                    Memory available for SCF calculation:          1.0 MB
       ====>        Please increase MaxCore to more than:         29.9 MB
Error (ORCA_SCF): ... aborting the run
ORCA finished by error termination in LEANSCF
  .... aborting the run
"""

# ORCA 6.1.1, SCF held to two iterations on an open-shell iron oxide.
_SCF_NOT_CONVERGED = """
[file orca_leanscf/orca_leanscf.cpp, line 306]: Error (ORCA_LEANSCF): unfortunately,
the SCF has not converged. There may be a way out but we have to stop here
ORCA finished by error termination in LEANSCF
  .... aborting the run
"""

# ORCA 6.1.1, ! HF-3c OPT with %geom MaxIter 1 — note the normal termination.
_OPT_GAVE_UP = """
                    ----------------------|Geometry convergence|---------------
                    MAX gradient        0.0250745046   0.0003000000      NO

                    The optimization did not converge but reached the maximum
                    number of optimization cycles.

                             ****ORCA TERMINATED NORMALLY****
TOTAL RUN TIME: 0 days 0 hours 0 minutes 12 seconds
"""

_OPT_CONVERGED = """
                    ***********************HURRAY********************
                    ***        THE OPTIMIZATION HAS CONVERGED     ***
                                *** OPTIMIZATION RUN DONE ***

                             ****ORCA TERMINATED NORMALLY****
"""


def _out(tmp_path: Path, body: str, name: str = "job.out") -> Path:
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------

def test_running_out_of_memory_is_not_an_mpi_crash(tmp_path):
    assert OrcaErrorDetector.analyze_output(_out(tmp_path, _OUT_OF_MEMORY)) is (
        OrcaErrorType.MEMORY_ERROR
    )


def test_an_scf_that_will_not_converge_is_still_recognised(tmp_path):
    """ORCA 6.1 reports SCF failures through LEANSCF, so this must not be
    swept into the memory rule that now outranks the generic LEANSCF one."""
    assert OrcaErrorDetector.analyze_output(_out(tmp_path, _SCF_NOT_CONVERGED)) is (
        OrcaErrorType.LEANSCF_NOT_CONVERGED
    )


def test_a_normal_termination_is_never_an_error(tmp_path):
    assert OrcaErrorDetector.analyze_output(_out(tmp_path, _OPT_CONVERGED)) is None


def test_the_output_is_not_read_into_memory_whole(tmp_path):
    """These outputs reach hundreds of MB and this runs on every failure,
    often for several parallel jobs at once — precisely when memory is
    already tight. Only the tail may be touched."""
    padding = ("x" * 4096 + "\n") * 3000  # ~12 MB of irrelevant output
    out = _out(tmp_path, padding + _OUT_OF_MEMORY)
    assert out.stat().st_size > 8 * 1024 * 1024

    import tracemalloc

    tracemalloc.start()
    detected = OrcaErrorDetector.analyze_output(out)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    assert detected is OrcaErrorType.MEMORY_ERROR, "the error sits in the tail"
    assert peak < 8 * 1024 * 1024, (
        f"read {peak / 1024 / 1024:.1f} MB for a {out.stat().st_size / 1024 / 1024:.1f} MB "
        "output — the whole file is being loaded again"
    )


# --------------------------------------------------------------------------
# the fix ORCA asked for
# --------------------------------------------------------------------------

def test_the_maxcore_orca_named_is_used(tmp_path):
    strategy = RecoveryStrategy(OrcaErrorType.MEMORY_ERROR, 1, {})
    strategy.output_file = _out(tmp_path, _OUT_OF_MEMORY)

    mods = strategy.get_modifications()

    # 29.9 MB with a margin, so the retry does not fail just short of it.
    assert mods["set_maxcore"] >= 30
    assert mods["set_maxcore"] == int(29.9 * 1.5) + 1


def test_the_largest_demand_wins_when_several_modules_complain(tmp_path):
    body = _OUT_OF_MEMORY + "\n       ====>  Please increase MaxCore to more than:  512.0 MB\n"
    strategy = RecoveryStrategy(OrcaErrorType.MEMORY_ERROR, 1, {})
    strategy.output_file = _out(tmp_path, body)

    assert strategy.get_modifications()["set_maxcore"] == int(512.0 * 1.5) + 1


def test_without_a_figure_the_old_behaviour_stands(tmp_path):
    """No output, or one that names no number, must leave the core-reduction
    fallback exactly as it was."""
    strategy = RecoveryStrategy(OrcaErrorType.MEMORY_ERROR, 1, {})

    mods = strategy.get_modifications()

    assert "set_maxcore" not in mods
    assert mods["reduce_pal"] == 0.5
    assert "SOSCF" in mods["keywords_remove"]


def test_an_unreadable_output_does_not_break_recovery(tmp_path):
    strategy = RecoveryStrategy(OrcaErrorType.MEMORY_ERROR, 2, {})
    strategy.output_file = tmp_path / "does-not-exist.out"

    mods = strategy.get_modifications()

    assert "set_maxcore" not in mods
    assert mods["reduce_pal"] == 0.25


def test_maxcore_is_raised_but_never_lowered(tmp_path):
    """The retry exists because the previous allocation was too small. Giving
    the rerun less would guarantee the same failure."""
    inp = tmp_path / "job.inp"
    inp.write_text("! HF-3c\n%maxcore 4000\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
    modifier = OrcaInputModifier(inp, {})

    parsed = {"blocks": {"maxcore": "%maxcore 4000"}, "keywords": []}
    assert modifier._set_maxcore(parsed, 100)["blocks"]["maxcore"] == "%maxcore 4000"
    assert modifier._set_maxcore(parsed, 9000)["blocks"]["maxcore"] == "%maxcore 9000"


# --------------------------------------------------------------------------
# an optimisation that gave up is not a finished job
# --------------------------------------------------------------------------

def test_an_optimisation_that_ran_out_of_cycles_is_not_complete(tmp_path):
    out = _out(tmp_path, _OPT_GAVE_UP)
    inp = tmp_path / "job.inp"
    inp.write_text("! HF-3c OPT\n", encoding="utf-8")

    assert smart_recalc.has_ok_marker(out) is True, "ORCA does end normally"
    assert smart_recalc.optimization_gave_up(out) is True
    assert smart_recalc.outputs_complete(inp, out) is False


def test_a_converged_optimisation_stays_complete(tmp_path):
    out = _out(tmp_path, _OPT_CONVERGED)
    inp = tmp_path / "job.inp"
    inp.write_text("! HF-3c OPT\n", encoding="utf-8")

    assert smart_recalc.outputs_complete(inp, out) is True


@pytest.mark.parametrize(
    "body",
    [
        "SINGLE POINT ENERGY\n****ORCA TERMINATED NORMALLY****\n",
        "ORCA NUMERICAL FREQUENCIES\nVIBRATIONAL FREQUENCIES\n****ORCA TERMINATED NORMALLY****\n",
    ],
)
def test_jobs_that_never_optimise_are_untouched(tmp_path, body):
    """Testing for ORCA's explicit "did not converge" rather than for the
    absence of a convergence banner is what keeps single points and
    frequency-only runs out of this."""
    out = _out(tmp_path, body)
    inp = tmp_path / "job.inp"
    inp.write_text("! HF-3c\n", encoding="utf-8")

    assert smart_recalc.optimization_gave_up(out) is False
    assert smart_recalc.outputs_complete(inp, out) is True


def test_the_verdict_is_found_far_from_the_end_of_the_output(tmp_path):
    """In a real ESD run the optimiser's verdict sat 1 MB before the end,
    behind the final energy evaluation, the property block and a frequency
    job. A short tail window would have missed it and called the job done."""
    trailing = ("property line padding\n" * 40000)  # ~0.8 MB after the verdict
    out = _out(tmp_path, _OPT_GAVE_UP + trailing + "\n****ORCA TERMINATED NORMALLY****\n")

    assert smart_recalc.optimization_gave_up(out) is True


def test_a_missing_output_is_not_reported_as_a_failed_optimisation(tmp_path):
    assert smart_recalc.optimization_gave_up(tmp_path / "nope.out") is False
