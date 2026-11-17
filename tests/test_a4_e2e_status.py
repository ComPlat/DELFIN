"""E2E status checks for a finished MLIP GOAT run (package A4).

``delfin.mlp_tools.goat_e2e_check.check_goat_e2e`` verifies that a GOAT
run the operator submitted through DELFIN's ExtOpt path actually
finished and produced a usable conformer ensemble.  This file covers the
STRUCTURAL status checks: the ensemble exists and is non-empty, the
global minimum file is present, and the ORCA run terminated cleanly.

Input: a SYNTHETIC run directory built in the test (invented numbers in
the format ORCA's ``read_goat_ensemble`` measures, NEVER ORCA output --
LICENSE_GUARD.md family 2).  Output: the report dict from
``check_goat_e2e``.

Control note: ``goat_e2e_check.py`` does not exist on main -- the import
fails and every test here is red until the module is added.  The red
state is the defect (no checker exists for a finished MLIP GOAT run).
"""

from __future__ import annotations

import pytest

from delfin.mlp_tools.goat_e2e_check import check_goat_e2e

# SYNTHETIC ensemble: a 3-frame water "GOAT result" in the "<E>
# converged=true" comment form read_goat_ensemble parses (ORCA 6.0.1/6.1.1
# measured form).
GOOD_ENSEMBLE = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
3
-76.4000000000 converged=true
O 0.000000 0.000000 0.000000
H 0.757000 0.586000 0.000000
H -0.757000 0.586000 0.000000
3
-76.3990000000 converged=true
O 0.000000 0.000000 0.100000
H 0.757000 0.586000 0.100000
H -0.757000 0.586000 0.100000
3
-76.3980000000 converged=true
O 0.000000 0.000000 0.200000
H 0.757000 0.586000 0.200000
H -0.757000 0.586000 0.200000
"""

GOOD_GLOBAL_MINIMUM = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
3
-76.4000000000
O 0.000000 0.000000 0.000000
H 0.757000 0.586000 0.000000
H -0.757000 0.586000 0.000000
"""

# A minimal ORCA output that ENDS cleanly -- the marker ORCA writes after
# a terminated run.  The body is invented (no ORCA block content needed
# for the structural checks in this file).
GOOD_OUT = """\
# SYNTHETIC ORCA output for a finished MLIP GOAT run -- not real ORCA
<MLIP GOAT run body>
****ORCA TERMINATED NORMALLY****
"""

# The ORCA output of a run that aborted before finishing.
ABORTED_OUT = """\
# SYNTHETIC ORCA output of an aborted MLIP GOAT run -- not real ORCA
ERROR: geometry never converged
"""


@pytest.fixture()
def good_run_dir(tmp_path):
    """A directory that looks like a finished GOAT job."""
    (tmp_path / "extopt.finalensemble.xyz").write_text(GOOD_ENSEMBLE)
    (tmp_path / "extopt.globalminimum.xyz").write_text(GOOD_GLOBAL_MINIMUM)
    (tmp_path / "extopt.out").write_text(GOOD_OUT)
    return tmp_path


def test_good_run_is_ok(good_run_dir):
    report = check_goat_e2e(good_run_dir, base_name="extopt")
    assert report["ok"] is True


def test_good_run_reports_nonempty_ensemble(good_run_dir):
    report = check_goat_e2e(good_run_dir, base_name="extopt")
    assert report["n_conformers"] > 0
    assert report["ensemble_file"] is not None


def test_good_run_reports_global_minimum_present(good_run_dir):
    report = check_goat_e2e(good_run_dir, base_name="extopt")
    assert report["global_minimum_file"] is not None


def test_missing_ensemble_is_not_ok_and_names_it(tmp_path):
    # only the .out and the global minimum are present; no ensemble
    (tmp_path / "extopt.globalminimum.xyz").write_text(GOOD_GLOBAL_MINIMUM)
    (tmp_path / "extopt.out").write_text(GOOD_OUT)
    report = check_goat_e2e(tmp_path, base_name="extopt")
    assert report["ok"] is False
    assert any("ensemble" in f.get("message", "") for f in report["findings"])


def test_aborted_run_is_not_ok_and_names_termination(tmp_path):
    # ensemble + global minimum exist but ORCA aborted
    (tmp_path / "extopt.finalensemble.xyz").write_text(GOOD_ENSEMBLE)
    (tmp_path / "extopt.globalminimum.xyz").write_text(GOOD_GLOBAL_MINIMUM)
    (tmp_path / "extopt.out").write_text(ABORTED_OUT)
    report = check_goat_e2e(tmp_path, base_name="extopt")
    assert report["ok"] is False
    assert not report["terminated_normally"]


def test_stale_run_marker_warns_to_check_the_scheduler(good_run_dir):
    # a finished-looking run that still carries a live-run .pid lock
    (good_run_dir / "extopt.pid").write_text("4242\n")
    report = check_goat_e2e(good_run_dir, base_name="extopt")
    assert report["ok"] is True  # a stale marker is a warning, not a failure
    assert any("marker" in f.get("message", "")
               for f in report["findings"])


def test_clean_run_has_no_live_run_marker_finding(good_run_dir):
    report = check_goat_e2e(good_run_dir, base_name="extopt")
    assert not any("marker" in f.get("message", "")
                   for f in report["findings"])
