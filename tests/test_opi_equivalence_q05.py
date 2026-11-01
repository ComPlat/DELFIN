"""OPI-vs-DELFIN equivalence for run status read from ORCA output.

The four fixtures under ``tests/fixtures/orca611/`` are small, cleaned
ORCA 6.1.1 runs (H2O opt+freq, O2 triplet SP, H2O TD-DFT SP, and an SCF
run that did NOT converge).  Each has the ``.out`` that both sides read,
plus the ``.json`` / ``.property.json`` that OPI's full ``parse()``
consumes (kept so the fixtures double as OPI exemplars).

This test pins the contract the task asked to measure: for each run,
DELFIN's own run-status readers and OPI (``terminated_normally``,
``scf_converged``, ``geometry_optimization_converged``) read the SAME
``.out`` and must agree on the status -- or, where they cannot agree,
the test records which side is right and why, so a later switch to OPI
does not lose or distort a DELFIN verdict.

Critical measured facts (see scratch/probe_delfin_opi.py, and OPI
2.0.0 source opi/output/grepper/recipes.py):

- OPI's three status methods read the ``.out`` TEXT directly:
    terminated_normally()              -> substring "****ORCA TERMINATED NORMALLY****"
    scf_converged()                    -> substring "SUCCESS"
    geometry_optimization_converged()  -> substring "HURRAY"
  They need no property.json and no orca_2json.
- DELFIN's readers (delfin/api.py parse_orca_output, delfin/smart_recalc
  has_ok_marker, delfin/doc_server/calc_indexer outcome_of_folder):
    scf_converged  : True on "SCF CONVERGED [AFTER]", False on
                     "SCF NOT CONVERGED"/"SCF iterations did not converge"
    opt_converged  : True on "OPTIMIZATION RUN DONE"/"THE OPTIMIZATION HAS
                     CONVERGED", False on "OPTIMIZATION DID NOT CONVERGE"/
                     "FAILED TO CONVERGE THE GEOMETRY OPTIMIZATION",
                     None when no geometry optimization ran.
    has_ok_marker  : bool, "ORCA TERMINATED NORMALLY" substring.

Agreement measured on the four fixtures:
    h2o_optfreq   : terminated True/True, scf True/True, geom True/True
    o2_triplet    : terminated True/True, scf True/True, geom None vs False
    h2o_tddft     : terminated True/True, scf True/True, geom None vs False
    scf_not_converged: terminated False/False, scf False/False,
                    final_energy None/None  (both call it a failure)

The one real difference is geometry on non-optimization runs: DELFIN
returns ``None`` (job type did not request a geometry optimization --
the absence of a marker is NOT a failure), while OPI returns ``False``
because "HURRAY" is absent.  The test asserts DELFIN's ``None`` contract
is preserved and records that OPI's ``False`` must be re-read as
"not applicable" (job did not optimize) before it may be used.  This is
exactly the "no DELFIN function may get worse" line: DELFIN's None is
the more careful reading.

Everything here is pure Python over the committed text fixtures: no
ORCA, no orca_2json, no OPI parse() is ever invoked.
"""
from __future__ import annotations

from pathlib import Path

import pytest

orcapi = pytest.importorskip("opi", reason="OPI (orca-pi) not installed")

from opi.output.core import Output

from delfin.api import parse_orca_output
from delfin.smart_recalc import has_ok_marker

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir as _lc1_fixture_dir  # noqa: E402


def _fixture(name: str) -> Path:
    return _lc1_fixture_dir(f"tests/fixtures/orca611/{name}")


# name -> (expect terminated_normally, expect scf_converged,
#          DELFIN opt_converged, OPI geometry_optimization_converged)
# Measured once on the committed fixtures (probe script in scratch/).
EXPECT = {
    # terminated, scf, delfin_geom, opi_geom, run kind
    "h2o_optfreq": (True, True, True, True, "opt"),
    "o2_triplet": (True, True, None, False, "sp"),
    "h2o_tddft": (True, True, None, False, "sp"),
    "scf_not_converged": (False, False, None, False, "failed"),
}


def _opi_status(name: str):
    """OPI's three status booleans for a fixture; reads the .out directly."""
    out = Output(basename=name, working_dir=_fixture(name))
    return (
        out.terminated_normally(),
        out.scf_converged(),
        out.geometry_optimization_converged(),
    )


def _delfin_status(name: str):
    from delfin.doc_server.calc_indexer import outcome_of_folder

    out = _fixture(name) / f"{name}.out"
    parsed = parse_orca_output(str(out))
    folder = _fixture(name)
    return {
        "outcome": outcome_of_folder(folder),
        "scf_converged": parsed.scf_converged,
        "opt_converged": parsed.opt_converged,
        "has_ok_marker": has_ok_marker(str(out)),
        "final_single_point": parsed.final_single_point,
    }


@pytest.mark.parametrize("name", sorted(EXPECT))
def test_termination_delfin_matches_opi(name: str):
    """Normal termination: DELFIN has_ok_marker / outcome == OPI terminated_normally()."""
    expect = EXPECT[name][0]
    st = _delfin_status(name)
    opi_term, _, _ = _opi_status(name)
    # has_ok_marker (substring "ORCA TERMINATED NORMALLY") must equal OPI's
    # terminated_normally (substring "****ORCA TERMINATED NORMALLY****").
    assert st["has_ok_marker"] == opi_term == expect
    # The outcome phrase and the boolean must agree on the same tail read.
    if expect:
        assert st["outcome"].startswith("finished")
    else:
        assert st["outcome"].startswith("failed")


@pytest.mark.parametrize("name", sorted(EXPECT))
def test_scf_convergence_delfin_matches_opi(name: str):
    """SCF convergence: DELFIN parse_orca_output.scf_converged == OPI scf_converged()."""
    expect = EXPECT[name][1]
    st = _delfin_status(name)
    _, opi_scf, _ = _opi_status(name)
    assert st["scf_converged"] == opi_scf == expect
    # Both must agree on the not-converged run (sane on failure).
    if name == "scf_not_converged":
        assert opi_scf is False
        assert st["scf_converged"] is False


@pytest.mark.parametrize("name", sorted(EXPECT))
def test_geometry_convergence_semantics(name: str):
    """Geometry: DELFIN opt_converged and OPI geometry_optimization_converged.

    On an optimization job they must agree.  On a single-point job DELFIN
    returns ``None`` (no geometry optimization was requested -- absence of
    a convergence marker is NOT a failure) and OPI returns ``False``
    because the "HURRAY" marker is absent.  The test pins DELFIN's None
    contract as the correct reading and forbids any transition that turns
    it into OPI's False without an applicability guard.
    """
    st = _delfin_status(name)
    opi_term, _, opi_geom_get = _opi_status(name)
    kind = EXPECT[name][-1]

    if kind == "opt":
        # Optimization job: geom convergence is meaningful; both agree.
        assert st["opt_converged"] is True
        assert opi_geom_get is True
        assert st["opt_converged"] == opi_geom_get
    elif kind == "sp":
        # Single-point job: DELFIN must return None (not-applicable).
        assert st["opt_converged"] is None, (
            "DELFIN must keep None (job did not optimize) for a single-point "
            "run; OPI returns False only because 'HURRAY' is absent."
        )
        assert opi_geom_get is False
    elif kind == "failed":
        assert st["opt_converged"] is None
        assert opi_geom_get is False


def test_scf_not_converged_run_is_sensible_everywhere():
    """The failed run: every reader must agree it failed and yield no energy."""
    name = "scf_not_converged"
    st = _delfin_status(name)
    opi_term, opi_scf, opi_geom = _opi_status(name)

    assert opi_term is False and opi_scf is False and opi_geom is False
    assert st["has_ok_marker"] is False
    assert st["outcome"].startswith("failed")
    assert st["scf_converged"] is False
    # No energy may be reported from an unconverged run.
    assert st["final_single_point"] is None
