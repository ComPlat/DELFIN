"""The Delta E_ST grading carries the cited experimental band (F6).

Reviewer finding F6 (nacht-s8) after the operator's science check on
job 8181368: the run's adiabatic Delta E_ST came out 0.889 eV while
the experimental adiabatic splitting of formaldehyde is about 0.37 eV
-- the loose sign+|dE_ST|<1 eV check passed a value 2.4x the
experimental one.  The comparison therefore follows the same pattern
as the ISC reference (esd_numerics.isc_reference_comparison): a cited
experimental band, level-gated, NOT GRADED where it does not apply.

Reference: Clouthier & Ramsay, "The Spectroscopy of Formaldehyde and
Thioformaldehyde", Annu. Rev. Phys. Chem. 34, 31-58 (1983), DOI
10.1146/annurev.pc.34.100183.000335 (DOI verified on annualreviews.org
via DELFIN's web tools).  The specific 0-0 origins used here
(T1 ~25190 cm-1, S1 ~28190 cm-1 -> Delta E_ST ~3000 cm-1 = 0.37 eV)
were provided by the operator and could not be independently verified
from this session (publisher 403); provenance is carried in the
product constants, not silently assumed.
"""

from __future__ import annotations

from delfin.common.esd_numerics import (
    EST_REFERENCE,
    est_reference_comparison,
)

import pytest


def test_est_reference_constants_carry_their_source():
    assert EST_REFERENCE["band_ev"] == pytest.approx((0.25, 0.55))
    assert "10.1146/annurev.pc.34.100183.000335" in EST_REFERENCE["source"]
    # Provenance must be honest: the origins are operator-provided.
    assert "operator" in EST_REFERENCE["source"].lower()


def test_est_inside_the_band_passes():
    record = est_reference_comparison(0.37)
    assert record["criterion"] == "delta_E_ST_vs_reference"
    assert record["verdict"] == "PASS"


def test_est_outside_the_band_fails_the_loose_check_would_have_passed():
    # 0.889 eV passes a sign+<1 eV check but is ~2.4x the experimental
    # adiabatic splitting: this is the value the loose criterion
    # accepted on job 8181368 and the defect this criterion exists for.
    record = est_reference_comparison(0.889)
    assert record["verdict"] == "FAIL"
    assert "0.37" in record["value"] or "band" in record["value"]


def test_est_t1_above_s1_is_a_rank_failure_regardless_of_band():
    # A positive Delta E_ST (T1 above S1) kills the TADF ranking claim
    # even when it lies inside the band: sign carries its own verdict.
    record = est_reference_comparison(-0.40)
    assert record["verdict"] == "FAIL"
    assert "sign" in record["value"].lower() or "T1" in record["value"]


def test_est_none_is_not_graded():
    record = est_reference_comparison(None)
    assert record["verdict"] == "NOT GRADED"
    assert "no" in record["value"].lower()
