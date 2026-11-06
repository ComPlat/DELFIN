"""Energy, GPU and timing checks for a finished MLIP GOAT run (pkg A4).

Phase 2 of package A4 over delfin.mlp_tools.goat_e2e_check: beyond the
structural status (ensemble non-empty, global minimum present, clean
termination), the checker must judge the *ensemble energies* -- finite
and sorted (GOAT writes the lowest conformer first, so the energies are
monotone non-decreasing) -- and record whether the MLIP server actually
ran on the GPU (log evidence, not assumption) and the wall time.

These behaviours did not exist in commit 1ca1d7cef, so the tests here
are red until the implementation lands: the report keys
``energies_sorted``, ``gpu_evidence``, ``walltime_s`` are absent and the
findings do not yet name a sorting / GPU / timing problem.

Input: SYNTHETIC run directories built in the test (invented numbers in
the "<E> converged=true" comment form; never ORCA output,
LICENSE_GUARD.md family 2).
"""

from __future__ import annotations

import pytest

from delfin.mlp_tools.goat_e2e_check import check_goat_e2e

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

# Same three frames but the middle energy is HIGHER than the first: not
# sorted (GOAT writes lowest first).
UNSORTED_ENSEMBLE = GOOD_ENSEMBLE.replace(
    "-76.3990000000 converged=true", "-76.3995000000 converged=true"
)

# A frame whose energy is NaN.
NAN_ENSEMBLE = GOOD_ENSEMBLE.replace(
    "-76.3990000000 converged=true", "nan converged=true"
)

GOOD_GLOBAL_MINIMUM = """\
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
3
-76.4000000000
O 0.000000 0.000000 0.000000
H 0.757000 0.586000 0.000000
H -0.757000 0.586000 0.000000
"""

GOOD_OUT = """\
# SYNTHETIC ORCA output for a finished MLIP GOAT run -- not real ORCA
<MLIP GOAT run body>
****ORCA TERMINATED NORMALLY****
"""

# ORCA output carrying the timing block ORCA prints at the end.
OUT_WITH_TIMING = """\
# SYNTHETIC ORCA output with a final timing block -- not real ORCA
<MLIP GOAT run body>
****ORCA TERMINATED NORMALLY****
Total Job Time (sec): 1234.56
"""

# A server log line showing the MLIP actually used a CUDA device.
GPU_SERVER_LOG = """\
# SYNTHETIC orca-mlips server log -- not real output
[INFO] device=cuda:0 (NVIDIA A100)
"""


@pytest.fixture()
def good_run_dir(tmp_path):
    (tmp_path / "extopt.finalensemble.xyz").write_text(GOOD_ENSEMBLE)
    (tmp_path / "extopt.globalminimum.xyz").write_text(GOOD_GLOBAL_MINIMUM)
    (tmp_path / "extopt.out").write_text(GOOD_OUT)
    return tmp_path


def _write_base(tmp_path, *, ensemble, out=GOOD_OUT, extra=()):
    (tmp_path / "extopt.finalensemble.xyz").write_text(ensemble)
    (tmp_path / "extopt.globalminimum.xyz").write_text(GOOD_GLOBAL_MINIMUM)
    (tmp_path / "extopt.out").write_text(out)
    for name, text in extra:
        (tmp_path / name).write_text(text)


def test_sorted_energies_are_reported_sorted(good_run_dir):
    report = check_goat_e2e(good_run_dir, base_name="extopt")
    assert report["energies_sorted"] is True


def test_unsorted_ensemble_is_detected_and_names_order(tmp_path):
    _write_base(tmp_path, ensemble=UNSORTED_ENSEMBLE)
    report = check_goat_e2e(tmp_path, base_name="extopt")
    assert report["energies_sorted"] is False
    assert any("sorted" in f.get("message", "") for f in report["findings"])


def test_nonfinite_energy_is_detected(tmp_path):
    _write_base(tmp_path, ensemble=NAN_ENSEMBLE)
    report = check_goat_e2e(tmp_path, base_name="extopt")
    assert report["ok"] is False
    assert any("finite" in f.get("message", "") for f in report["findings"])


def test_gpu_evidence_is_reported_when_a_server_log_shows_it(tmp_path):
    _write_base(
        tmp_path, ensemble=GOOD_ENSEMBLE, extra=(("orca_mlips_server.log", GPU_SERVER_LOG),)
    )
    report = check_goat_e2e(tmp_path, base_name="extopt")
    assert report["gpu_evidence"] is not None
    assert "cuda" in report["gpu_evidence"].lower() or "gpu" in report["gpu_evidence"].lower()


def test_missing_gpu_evidence_is_a_warn_finding(tmp_path):
    # a plausible finished run with no server log / no GPU line
    _write_base(tmp_path, ensemble=GOOD_ENSEMBLE)
    report = check_goat_e2e(tmp_path, base_name="extopt")
    assert report["gpu_evidence"] is None
    assert any("GPU" in f.get("message", "") or "gpu" in f.get("message", "")
               for f in report["findings"])


def test_wall_time_is_read_from_the_out(tmp_path):
    _write_base(tmp_path, ensemble=GOOD_ENSEMBLE, out=OUT_WITH_TIMING)
    report = check_goat_e2e(tmp_path, base_name="extopt")
    assert report["walltime_s"] == pytest.approx(1234.56)
