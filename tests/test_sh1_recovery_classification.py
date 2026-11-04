"""SH1 red tests: the failure classes the evidence runs exposed.

Every test reads a REAL ORCA output (trimmed copies in
tests/fixtures/orca611/sh1_*), never a typed dict. Written RED on
54462fe67; the fixes make them green one phase each.
"""
from pathlib import Path

import pytest

from delfin.orca_recovery import OrcaErrorDetector, OrcaErrorType

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir


def _sh1(name: str) -> Path:
    """The specific orca611 subfolder: legacy path while its files are
    committed, else the private store (skip without the env var)."""
    return lc1_fixture_dir(f"tests/fixtures/orca611/{name}")


# ------------------------------------------------------------------ F1

def test_gstep_constraint_failure_is_classified():
    """v3: 'GSTEP: could not impose initial constraints' is a scan-setup
    problem the user must fix (the start geometry cannot satisfy its
    constraints) — it must be reported as such, not fall to UNKNOWN."""
    out = _sh1("sh1_gstep_constraints/sh1_gstep_constraints.out")
    detected = OrcaErrorDetector.analyze_output(out)
    assert detected == OrcaErrorType.SCAN_CONSTRAINT_ERROR, detected


def test_gstep_constraint_failure_is_not_retried_blindly():
    """The recovery for a scan-constraint failure must be to stop and
    tell the user, never a same-input retry (it fails identically)."""
    out = _sh1("sh1_gstep_constraints/sh1_gstep_constraints.out")
    detected = OrcaErrorDetector.analyze_output(out)
    assert detected is not None and detected != OrcaErrorType.UNKNOWN
    from delfin.orca_recovery import RecoveryStrategy
    manager = RecoveryStrategy(detected, attempt=1, config={})
    fixes = manager.get_modifications()
    assert fixes == {} or fixes.get("stop_for_user"), \
        f"a scan-constraint failure must not get mechanical fixes: {fixes}"


# ------------------------------------------------------------------ F2

def test_scan_step_maxcycles_with_shark_abort_is_classified():
    """v4: the constrained scan step ran out of cycles and the PROPERTIES
    (SHARK) stage aborted. Reviewer finding: on the full 8 MB output the
    priority-5 GEOMETRY_NOT_CONVERGED pair already fires, but the fix it
    prescribes (trust radius, MaxIter, LooseOpt) is aimed at a PLAIN
    optimisation — a relaxed-scan step whose CONSTRAINT cannot hold on
    the dihedral wrap is a different failure and gets the user-facing
    classification SCAN_STEP_NOT_CONVERGED, not trust-radius retries."""
    out = _sh1("sh1_scan_maxcycles_shark/sh1_scan_maxcycles_shark.out")
    detected = OrcaErrorDetector.analyze_output(out)
    assert detected == OrcaErrorType.SCAN_STEP_NOT_CONVERGED, detected


def test_plain_properties_abort_without_shark_stays_geometry(tmp_path):
    """Counter-direction of the (b) pair: a PLAIN optimisation that
    aborted in PROPERTIES without a SHARK read error keeps
    GEOMETRY_NOT_CONVERGED (and its trust-radius fixes). Both wordings
    must co-occur for SCAN_STEP_NOT_CONVERGED to fire."""
    text = (
        "        *                GEOMETRY OPTIMIZATION CYCLE   1            *\n"
        "   Geometry convergence            NOT CONVERGED\n"
        "  .... aborting the run\n"
        "ORCA finished by error termination in PROPERTIES\n"
    )
    out = tmp_path / "plain_properties.out"
    out.write_text(text, encoding="utf-8")
    detected = OrcaErrorDetector.analyze_output(out)
    assert detected == OrcaErrorType.GEOMETRY_NOT_CONVERGED, detected


def test_scan_step_failure_names_the_unconverged_scan_step():
    from delfin.orca_recovery import RecoveryStrategy
    out = _sh1("sh1_scan_maxcycles_shark/sh1_scan_maxcycles_shark.out")
    detected = OrcaErrorDetector.analyze_output(out)
    strategy = RecoveryStrategy(detected, attempt=1, config={})
    strategy.output_file = out
    reason = strategy.reason_for_user()
    assert reason and "scan" in reason.lower(), reason


# ------------------------------------------------------------------ F4

def _cu_hess() -> Path:
    return _sh1("sh1_cu_nh3_imagmode/sh1_cu_nh3_imagmode.hess")


def test_cu_v2_real_hess_has_the_imaginary_mode():
    """The fixture carries the REAL frequency block: mode 6 at -21.857."""
    from delfin.imag import hessian_frequencies
    freqs = hessian_frequencies(_cu_hess())
    assert len(freqs) == 51
    assert freqs[6] == pytest.approx(-21.8570393312490694, abs=1e-6)


def test_cu_v2_default_noise_floor_skips_the_mode_silently():
    """The sourced default (-50 cm-1 noise band, control_validator.py:1234)
    must keep skipping -21.86 — and not warn about it."""
    from delfin.imag import imaginary_modes, saddle_reason
    config = {}
    assert imaginary_modes(_cu_hess(), config) == []
    assert saddle_reason("cu", _cu_hess(), config) is None


def test_cu_v2_explicit_allow_imaginary_freq_lets_imag_act():
    """QS ruling: an explicit allow_imaginary_freq shallower than the
    measured mode (20 -> threshold -20) must make IMAG act on the real
    frequency block: the mode is listed and the saddle reason names it."""
    from delfin.imag import imaginary_modes, saddle_reason
    config = {"allow_imaginary_freq": 20}
    modes = imaginary_modes(_cu_hess(), config)
    assert modes and modes[0][0] == 6
    assert modes[0][1] == pytest.approx(-21.8570393312490694, abs=1e-6)
    reason = saddle_reason("cu", _cu_hess(), config)
    assert reason and "-22 cm-1" in reason, reason
