"""The recovery verdict agrees with orca_read's normal-termination reading.

The task asked to switch the normal-termination read onto delfin.orca_read
where behaviour provably does not change. A baseline recorder measured the
two readers on the four ORCA 6.1.1 fixtures plus 200 archived runs and
found them value-identical (scratch/baseline.json vs scratch/after.json,
0 diffs; scratch/PR.md). This test pins that agreement on the committed
fixtures: a run ``orca_read`` reports as TERMINATED_NORMALLY is never a
recovery error, and the SCF-not-converged fixture (not terminated) is the
recovery error the engine must answer.

It is GREEN with the current (unswitched) code, capturing the measured
equivalence without forcing the resource-changing switch of
OrcaErrorDetector.analyze_output.
"""
from pathlib import Path

from delfin import orca_read
from delfin.orca_recovery import OrcaErrorDetector, OrcaErrorType

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir


def _fixture_root(name: str) -> Path:
    """The specific orca611 subfolder: legacy path while its files are
    committed, else the private store (skip without the env var)."""
    return lc1_fixture_dir(f"tests/fixtures/orca611/{name}")

# name -> (orca_read TERMINATED_NORMALLY expected, analyze_output expected)
# measured once by scratch/recorder.py over the committed fixtures.
EXPECT = {
    "h2o_optfreq": (True, None),
    "h2o_tddft": (True, None),
    "o2_triplet": (True, None),
    "scf_not_converged": (False, OrcaErrorType.LEANSCF_NOT_CONVERGED),
}


def test_recovery_verdict_matches_orca_read_termination():
    for name, (term_expected, ana_expected) in EXPECT.items():
        out = _fixture_root(name) / f"{name}.out"
        term = orca_read.read(out.parent, name, orca_read.Quantity.TERMINATED_NORMALLY)
        ana = OrcaErrorDetector.analyze_output(out)
        assert term.value is term_expected, f"{name}: orca_read term {term.value}"
        assert ana is ana_expected, f"{name}: analyze_output {ana}"


def test_a_fixture_orca_read_says_terminated_is_never_a_recovery_error():
    for name in ("h2o_optfreq", "h2o_tddft", "o2_triplet"):
        out = _fixture_root(name) / f"{name}.out"
        term = orca_read.read(out.parent, name, orca_read.Quantity.TERMINATED_NORMALLY)
        assert term.value is True
        assert OrcaErrorDetector.analyze_output(out) is None
