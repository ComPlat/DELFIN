"""Adversarial tests for SH1 (nacht-s10, reviewer of package SH1).

Each case reproduces ORCA's verbatim failure wording from the wave-24 SH1
evidence (v2/v3/v4/cu_v1/cu_v2) into a temporary output file and checks what
DELFIN's own classification (OrcaErrorDetector) makes of it.  No fixture file
committed: the words are copied from the real runs, trimmed to what the
detector reads.

Control: red on main 54462fe67 before any builder change.
"""
import importlib
from pathlib import Path

import pytest

recovery = importlib.import_module("delfin.orca_recovery")

# Verbatim tail of v3.out (evidence lines 293-300): a relaxed scan whose start
# geometry cannot satisfy its constraints.  Nothing a retry changes.
V3_GSTEP = """\
**************************
*                       *
*  _____ ____   ____    *
Decreasing the constraint accuracy, check your results carefully!
Decreasing the constraint accuracy, check your results carefully!
[file orca_tools/Tool-GStep/qcgstep.cpp, line 521]: GSTEP: could not impose initial constraints. Please check your input and try again!
[file orca_tools/Tool-GStep/qcgstep.cpp, line 521]: GSTEP: could not impose initial constraints. Please check your input and try again!
[file orca_tools/Tool-Main/qcprog.cpp, line 123]: ....aborting the run
ORCA finished by error termination in GStep
"""

# Verbatim tail of v4.out: constrained scan step runs out of cycles, then
# SHARK fails in PROPERTIES and the job aborts.
V4_SHARK = """\
    ----------------------------------------------------------------------------
                                  WARNING !!!
       The optimization did not converge but reached the maximum number of
       optimization cycles.
       Please check your results very carefully.
    ----------------------------------------------------------------------------
ERROR (SHARK): Failed to read input file (rc1_butane_scan_r2scan3c_v4.SHARKINP.tmp)
ERROR (SHARK): Failed to read input file (rc1_butane_scan_r2scan3c_v4.SHARKINP.tmp)
ORCA finished by error termination in PROPERTIES
"""


def _out(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.parametrize("text", [V3_GSTEP], ids=["gstep-constraints"])
def test_sh1_gstep_scan_setup_is_classified(tmp_path, text):
    """GSTEP 'could not impose initial constraints' must not fall through to
    UNKNOWN/None: it is a scan-setup problem the user has to fix, and DELFIN
    must say so instead of retrying the identical input."""
    result = recovery.OrcaErrorDetector.analyze_output(_out(tmp_path, "v3.out", text))
    assert result is not None, "GSTEP scan-setup failure is not detected at all"
    assert result is not recovery.OrcaErrorType.UNKNOWN, (
        f"GSTEP failure classified {result}; needs a named class and a user reason")


@pytest.mark.parametrize("text", [V4_SHARK], ids=["shark-properties"])
def test_sh1_shark_abort_points_at_real_cause(tmp_path, text):
    """An unconverged scan step followed by a SHARK failure in PROPERTIES must
    not be answered with an MPI-style retry: the real cause is the
    unconverged constrained step."""
    result = recovery.OrcaErrorDetector.analyze_output(_out(tmp_path, "v4.out", text))
    assert result is not None, "SHARK abort after unconverged scan step not detected"
    assert result is not recovery.OrcaErrorType.UNKNOWN


def test_sh1_cuv2_single_small_imaginary_mode_is_offered_to_imag(tmp_path):
    """A converged optimisation with one imaginary mode at -21.86 cm-1 (the
    wave-23 Cu(NH3)4 case, healed by hand) must reach IMAG.  With IMAG's
    default noise floor of -50 cm-1 (delfin/imag.py:85) imaginary_modes()
    returns nothing, so IMAG would not act — the exact case it exists for."""
    imag = importlib.import_module("delfin.imag")
    out = _out(tmp_path, "cu_v2.out",
               "     6:     -21.86 cm**-1  ***imaginary mode***\n")
    modes = imag.collect_imaginary_modes(str(out))
    assert modes and modes[0][1] == pytest.approx(-21.86, abs=0.02)
    hess = tmp_path / "cu_v2.hess"
    hess.write_text(
        "$vibrational_frequencies\n1\n     6:     -21.86\n", encoding="utf-8")
    assert imag.imaginary_modes(hess, {}) == [], (
        "imaginary_modes with the default threshold skips -21.86 cm-1: "
        "IMAG will not heal the case it exists for without a mechanism that "
        "reaches it regardless of the noise floor")
