"""Tests for conformer-aware SEATING of large ligands on the FF-free path
(env DELFIN_FFFREE_CONFORMER_SEATING).

WHY (ISOMER_COMPLETENESS_AUDIT_2026_06_18 §4/§6): the dominant reach gap (~51.5 %
of the pool) is decompose()'s per-donor-arm heavy-atom cap (>8 heavy / donor ->
whole complex bails to legacy = 0 enumerated isomers).  Large/conjugated ligands
(substituted phosphines, polypyridyl, big arenes) hit it.  Conformer-aware seating
raises the cap so those complexes reach the FF-free build, where the per-ligand
conformer selection (+ a core-preserving backbone-re-embed reseat fallback) places
the large ligand and the self-gate keeps the never-worse guarantee.

Hard contracts validated here:
  1. byte-identical when the flag is OFF (default) -> path untouched;
  2. a LARGE-LIGAND complex (>8 heavy/donor) that the DEFAULT rejects to legacy
     builds >=1 clean native isomer with the flag ON;
  3. the core invariant: metal + all donors stay native (<=0.05 A) in built frames;
  4. legacy fallback (None) when no conformer seats cleanly -> never-worse;
  5. determinism: same input twice ON -> byte-identical output.

No CCDC coordinates -- a synthetic large-arm Werner complex + a small frozen-core
synthetic frame for the reseat unit test.
"""
import os
import numpy as np
import pytest

from delfin._bond_decollapse import _is_metal as bd_is_metal

# CN_EXTEND keeps the CN2 test-helper machinery available; BUILDER is required for the
# isomer path.  Seating is the flag under test (added per-test).
_BASE_ENV = {"DELFIN_FFFREE_BUILDER": "1", "DELFIN_FFFREE_CN_EXTEND": "1"}

# A large-LIGAND CN4 Werner complex: four O-naphthyl donor arms (~11 heavy atoms per
# single-donor arm, well over the default cap of 8).  decompose() rejects it to legacy
# at the default cap; with seating the cap is raised and it builds.
_LARGE = ("c1ccc2ccccc2c1O[Ti](Oc1cccc2ccccc12)(Oc1cccc2ccccc12)"
          "Oc1cccc2ccccc12")
# A small monodentate complex that builds under the DEFAULT cap (NOT touched by
# seating) -- the no-regression control.
_SMALL = "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]"


def _isomers(smiles):
    from delfin.manta import converter_backend as CB
    return CB._fffree_isomers(smiles)


def _parse(xyz):
    syms = [ln.split()[0] for ln in xyz.splitlines()]
    P = np.array([[float(x) for x in ln.split()[1:4]] for ln in xyz.splitlines()])
    return syms, P


def _donors(syms, P, cn):
    """The cn heavy atoms closest to the metal (= the constructed donors)."""
    m = next(i for i, s in enumerate(syms) if bd_is_metal(s))
    d = sorted((float(np.linalg.norm(P[j] - P[m])), j)
               for j in range(len(syms)) if j != m and syms[j] != "H")
    return m, [j for _, j in d[:cn]]


@pytest.fixture(autouse=True)
def _clean_env():
    saved = dict(os.environ)
    os.environ.pop("DELFIN_FFFREE_CONFORMER_SEATING", None)
    for k, v in _BASE_ENV.items():
        os.environ[k] = v
    yield
    os.environ.clear()
    os.environ.update(saved)


def _sig(rr):
    if rr is None:
        return None
    return tuple((lab, xyz) for xyz, lab in rr)


def test_byte_identical_when_flag_off():
    """Flag unset == flag '0': the large complex still bails to legacy (None) and the
    small complex output is unchanged -> path byte-identical."""
    os.environ.pop("DELFIN_FFFREE_CONFORMER_SEATING", None)
    large_a, small_a = _isomers(_LARGE), _sig(_isomers(_SMALL))
    os.environ["DELFIN_FFFREE_CONFORMER_SEATING"] = "0"
    large_b, small_b = _isomers(_LARGE), _sig(_isomers(_SMALL))
    assert large_a is None and large_b is None      # heavy-cap-rejected at default cap
    assert small_a == small_b is not None


def test_large_ligand_recovered_with_flag_on():
    """The dominant lever: a >8-heavy/donor complex that DEFAULT rejects to legacy
    builds >=1 clean native isomer once seating raises the cap."""
    os.environ.pop("DELFIN_FFFREE_CONFORMER_SEATING", None)
    assert _isomers(_LARGE) is None, "default cap must reject this large ligand"
    os.environ["DELFIN_FFFREE_CONFORMER_SEATING"] = "1"
    on = _isomers(_LARGE)
    assert on is not None and len(on) >= 1, "seating must build the large ligand"
    for xyz, _ in on:
        _, P = _parse(xyz)
        assert np.all(np.isfinite(P))               # never non-finite


def test_core_invariant_metal_and_donors_native():
    """In the seated build the coordination core is preserved: every emitted frame's
    metal + donor atoms sit within 0.05 A of the first (native base) frame."""
    os.environ["DELFIN_FFFREE_CONFORMER_SEATING"] = "1"
    rr = _isomers(_LARGE)
    assert rr is not None
    bsyms, bP = _parse(rr[0][0])
    bm, donors = _donors(bsyms, bP, 4)
    for xyz, lab in rr:
        syms, P = _parse(xyz)
        if len(syms) != len(bsyms):
            continue
        m = next(i for i, s in enumerate(syms) if bd_is_metal(s))
        shift = bP[bm] - P[m]                        # align translation on the metal
        assert float(np.linalg.norm(P[m] + shift - bP[bm])) <= 0.05 + 1e-9
        for dgi in donors:
            assert float(np.linalg.norm(P[dgi] + shift - bP[dgi])) <= 0.05 + 1e-6, (
                f"donor {dgi} drifted in {lab}")


def test_small_complex_not_touched_by_seating():
    """No regression: a complex that already builds under the default cap is byte-
    identical with seating ON (seating only fires for >8-heavy/donor ligands)."""
    os.environ.pop("DELFIN_FFFREE_CONFORMER_SEATING", None)
    off = _sig(_isomers(_SMALL))
    os.environ["DELFIN_FFFREE_CONFORMER_SEATING"] = "1"
    on = _sig(_isomers(_SMALL))
    assert off == on is not None


def test_deterministic():
    os.environ["DELFIN_FFFREE_CONFORMER_SEATING"] = "1"
    a = _sig(_isomers(_LARGE))
    b = _sig(_isomers(_LARGE))
    assert a == b


def test_reseat_fallback_returns_none_when_infeasible():
    """Legacy fallback (never-worse): _seat_via_conformers returns None when the
    ligand cannot map onto the native frame (degenerate lig_groups) -> the caller
    bails to legacy rather than emitting a worse build."""
    from delfin.manta import converter_backend as CB
    os.environ["DELFIN_FFFREE_CONFORMER_SEATING"] = "1"
    # a trivially-degenerate (None lig_groups) reseat request must return None
    assert CB._seat_via_conformers("Co", None, ["Co"], np.zeros((1, 3))) is None


def test_seating_helpers_flag_gated():
    from delfin.manta import converter_backend as CB
    from delfin.manta import decompose as DEC
    os.environ.pop("DELFIN_FFFREE_CONFORMER_SEATING", None)
    assert CB._seating_enabled() is False
    assert DEC._heavy_cap() == DEC._HEAVY_CAP_DEFAULT == 8
    os.environ["DELFIN_FFFREE_CONFORMER_SEATING"] = "1"
    assert CB._seating_enabled() is True
    assert DEC._heavy_cap() > 8
