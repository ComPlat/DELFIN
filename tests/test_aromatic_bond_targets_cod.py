"""The aromatic bond targets and the donor-plane tolerance are the COD (CC0) statistics.

Pins the five aromatic targets and the beta band to the values of delfin/manta/PROVENANCE.md, and
checks that _beta_score measures the excess over that band.
"""
import numpy as np

from delfin.fffree import aromatic_bond_targets as abt
from delfin.manta import assemble_donor_plane as adp


def test_aromatic_targets_are_the_cod_medians():
    assert abt.aromatic_target("C", "C") == 1.387
    assert abt.aromatic_target("N", "C") == 1.360
    assert abt.aromatic_target("O", "C") == 1.374
    assert abt.aromatic_target("S", "C") == 1.729
    assert abt.aromatic_target("N", "N") == 1.352
    assert abt.aromatic_target("B", "N") is None


def test_aromatic_target_is_order_independent():
    for (a, b), v in abt.AROMATIC_BOND_TARGETS.items():
        assert abt.aromatic_target(a, b) == v
        assert abt.aromatic_target(b, a) == v


def _tilted_donor(beta_deg):
    """Metal at the origin, a donor N with two C neighbours in a plane tilted so that the
    N->M bond leaves that plane by beta_deg."""
    b = np.radians(beta_deg)
    n = np.array([2.0, 0.0, 0.0])
    u = np.array([np.cos(b), 0.0, np.sin(b)])
    c1 = n + 1.35 * (0.5 * u + 0.866 * np.array([0.0, 1.0, 0.0]))
    c2 = n + 1.35 * (0.5 * u - 0.866 * np.array([0.0, 1.0, 0.0]))
    return ["Fe", "N", "C", "C"], np.array([[0.0, 0.0, 0.0], n, c1, c2])


def test_beta_band_is_the_cod_value():
    assert adp._BETA_BAND == 5.2


def test_beta_score_measures_the_excess_over_the_band():
    syms, Q = _tilted_donor(5.0)
    assert abs(adp._donor_plane_beta(Q, syms, 1) - 5.0) < 1e-6
    assert adp._beta_score(syms, Q, [1]) == 0.0
    syms, Q = _tilted_donor(6.0)
    assert abs(adp._beta_score(syms, Q, [1]) - (6.0 - 5.2) ** 2) < 1e-6
