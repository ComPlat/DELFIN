"""Ligand-geometry-first rigid seating (DELFIN_FFFREE_LIGAND_RIGID).

Eye-find QEBLOC (2026-06-24): the per-donor radial rescale in
``_orient_chelate_to_vertices`` moves each donor INDEPENDENTLY onto its exact ideal
radius, which splays a RIGID cage cap (ring N-N-N 110->129 deg, donor-donor
distorted).  ``rigid=True`` instead does a pure rigid-body fit (rotation +
translation, NO scale, NO per-donor move) so the cap's internal geometry is
preserved EXACTLY and the coordination polyhedron is emergent.

These are pure-geometry unit tests on the seating function (no full build): they
prove (a) ``rigid=False`` is byte-identical to the historic per-donor rescale
(donors land at the exact target radius), and (b) ``rigid=True`` preserves every
internal pairwise distance of the ligand block (rigid body).
"""
import itertools
import numpy as np

from delfin.manta.assemble_complex import _orient_chelate_to_vertices


def _fac_targets(md=2.13):
    # three mutually-cis octahedron vertices (a fac triangle), scaled to M-D = md
    return [np.array([md, 0, 0]), np.array([0, md, 0]), np.array([0, 0, md])]


def _toy_cap():
    """A rigid tripod cap: B at apex + 3 donor N + a back C on each arm, metal at
    origin.  Donor-donor distances are deliberately NOT the ideal octahedral 3.01 A
    so the per-donor rescale must distort the triangle while rigid must not."""
    B = np.array([0.7, 0.7, 0.7])
    N = [np.array([1.9, 0.2, 0.2]), np.array([0.2, 1.9, 0.2]), np.array([0.2, 0.2, 1.9])]
    C = [n + (n - B) * 0.4 for n in N]            # a back atom per arm
    lP = np.array([B] + N + C, float)
    donor_idxs = [1, 2, 3]
    return lP, donor_idxs


def _pairwise(P):
    return np.array([np.linalg.norm(P[i] - P[j])
                     for i, j in itertools.combinations(range(len(P)), 2)])


def test_rigid_off_lands_donors_at_exact_radius():
    lP, dons = _toy_cap()
    tg = _fac_targets(2.13)
    Q = _orient_chelate_to_vertices(lP, dons, tg, asym=True, rigid=False)
    for d in dons:
        assert abs(np.linalg.norm(Q[d]) - 2.13) < 1e-6      # exact ideal radius


def test_rigid_on_preserves_internal_geometry():
    lP, dons = _toy_cap()
    tg = _fac_targets(2.13)
    Q = _orient_chelate_to_vertices(lP, dons, tg, asym=True, rigid=True)
    # rigid body: EVERY internal pairwise distance is preserved exactly
    assert np.allclose(_pairwise(Q), _pairwise(lP), atol=1e-6)


def test_rigid_on_donors_near_targets_emergent():
    lP, dons = _toy_cap()
    tg = _fac_targets(2.13)
    Q = _orient_chelate_to_vertices(lP, dons, tg, asym=True, rigid=True)
    # donors land NEAR the assigned vertices (best rigid fit), M-D emergent (not exact)
    for i, d in enumerate(dons):
        assert np.linalg.norm(Q[d] - tg[i]) < 1.0
    radii = [np.linalg.norm(Q[d]) for d in dons]
    assert not all(abs(r - 2.13) < 1e-6 for r in radii)     # emergent, not pinned


def test_rigid_determinism():
    lP, dons = _toy_cap()
    tg = _fac_targets(2.13)
    a = _orient_chelate_to_vertices(lP, dons, tg, asym=True, rigid=True)
    b = _orient_chelate_to_vertices(lP, dons, tg, asym=True, rigid=True)
    assert np.allclose(a, b, atol=0.0)
