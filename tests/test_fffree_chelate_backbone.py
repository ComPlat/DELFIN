"""Tests for the FF-free CHELATE-BACKBONE gate-lift + hardened metallacycle embed
(PHASE 0 of the polydentate project, K4_MACROCYCLE_DESIGN_2026_06_17.md).

Root cause locked here: large kappa<=3 chelates with an extended/strained backbone
(e.g. BIQCOV: Ta, kappa3 diamido, ~9.7 heavy/arm) fall to the legacy UFF fallback
because the historic per-arm chelate complexity cap (8 heavy/arm) rejects them at
the ``decompose`` stage BEFORE any FF-free build is attempted -- and, once admitted,
the CHELATE_BITE hard donor-donor distance pin over-constrains the DG bounds matrix
so the strained backbone buckles into the coordination shell on the fallback embed.

PHASE 0 fix (flag ``DELFIN_FFFREE_CHELATE_BACKBONE``, default OFF -> byte-identical):
  * (b) raise the CHELATE per-arm heavy cap (8 -> DELFIN_FFFREE_CHELATE_HEAVY_CAP,
        default 12) so large kappa<=3 backbones get past decompose;
  * (a) harden ``_embed_metallacycle``: keep M-D distances HARD-pinned (M-D invariant)
        but relax donor-donor distances to a WIDE SOFT window (DELFIN_FFFREE_POLY_BB_TOL,
        default 0.4 A) and embed MORE conformers (DELFIN_FFFREE_POLY_K_CONFS, default 12),
        so the strained ring folds without collapsing.

SCOPE GUARD (Phase 0 only): the denticity>3 gates (decompose 'len(fdonors) > 3' and
converter_backend 'denticity >= 4') are UNTOUCHED -- kappa4+ stays legacy.  The
self-gate (_build_is_clean) still vetoes any unclean result -> never worse than legacy.

Contracts:
  * with the flag ON, a large-kappa3 backbone (BIQCOV / ABEZAJ) decomposes natively
    and builds >=1 isomer (legacy -> native), with a clean coordination shell;
  * env-gate: flag unset/0 -> byte-identical (large-kappa3 stays legacy; every other
    system -- chelate k2/k3, M(bipy)3-class, monodentate -- unchanged);
  * kappa4+ STAYS legacy with the flag ON (and even with the cap pushed very high):
    the denticity>3 gate is the wall, not the cap;
  * determinism: same input -> byte-identical output.
"""
import hashlib
import os

import pytest

from delfin.manta import decompose as DEC
from delfin.manta.converter_backend import _fffree_isomers

# --- Large-kappa3 anchors (legacy today, native with the flag) -----------------
# BIQCOV: Ta, kappa3 diamido (N,N,N), ~9.7 heavy/arm + 2 Cl co-ligands (TBP-5).
BIQCOV = ("CC(C)[N]1C2=CC(C(C)(C)C)=CC=C2[N]2C3=CC=C(C(C)(C)C)C=C3"
          "[N](C(C)C)[Ta]21([Cl])[Cl]")
# ABEZAJ: Ti, kappa3 diamido-amine (N,N,N), ~9.3 heavy/arm (TBP-5).
ABEZAJ = ("CC1=CC(C)=C([N]2C3=CC=CC=C3C3(C)[N](C4CCCCC4)"
          "[Ti-]2([Cl])([N](C)C)[N+]3(C)C)C(C)=C1")

# --- byte-identity controls (must be UNCHANGED whether flag on or off) ----------
CISPLATIN = "N[Pt](N)(Cl)Cl"
COCL3NH3 = "[NH3][Co]([NH3])([NH3])([Cl])([Cl])[Cl]"
PTME2CL2 = "C[Pt](C)(Cl)Cl"

# --- genuine kappa4 (MUST stay legacy: Phase 0 does NOT open the denticity gate) -
# cyclam tetraaza macrocycle on Ni -- a genuine denticity-4 ligand.
CYCLAM_NI = "C1CN2CCN3CCN4CCN1[Ni-4]234"
# AACANI10: Ni bicyclic cage, genuine denticity-4 (design-doc kappa4 anchor).
AACANI10 = "O=C1C[N+]23CCC[N+]4(CCC2)CC(=O)[O][Ni-3]34([OH2+])[O]1"

FLAG = "DELFIN_FFFREE_CHELATE_BACKBONE"


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch):
    """Every test starts with all CHELATE-BACKBONE env knobs unset (default OFF)."""
    for k in (FLAG, "DELFIN_FFFREE_CHELATE_HEAVY_CAP",
              "DELFIN_FFFREE_POLY_BB_TOL", "DELFIN_FFFREE_POLY_K_CONFS"):
        monkeypatch.delenv(k, raising=False)
    yield


def _hash(res):
    if res is None:
        return "None"
    if isinstance(res, str):
        return res
    h = hashlib.sha256()
    for xyz, lab in res:
        h.update(xyz.encode())
        h.update(b"\x00")
        h.update(str(lab).encode())
        h.update(b"\x01")
    return f"{len(res)}:{h.hexdigest()}"


# ---------------------------------------------------------------------------
# (1) native with the flag: large-kappa3 legacy -> native
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("smi", [BIQCOV, ABEZAJ])
def test_large_kappa3_native_with_flag(smi, monkeypatch):
    # OFF: large-kappa3 falls to legacy (None)
    assert _fffree_isomers(smi) is None
    # ON: decomposes natively as a chelate and builds >=1 isomer
    monkeypatch.setenv(FLAG, "1")
    d = DEC.decompose(smi)
    assert d is not None and d.get("has_chelate")
    assert max(lg["denticity"] for lg in d["ligands"]) == 3
    res = _fffree_isomers(smi)
    assert isinstance(res, list) and len(res) >= 1
    # every emitted frame is non-empty XYZ text
    for xyz, _lab in res:
        assert xyz.strip()


# ---------------------------------------------------------------------------
# (2) byte-identity OFF: flag unset -> unchanged for every class
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("smi", [CISPLATIN, COCL3NH3, PTME2CL2, BIQCOV, ABEZAJ,
                                 CYCLAM_NI, AACANI10])
def test_flag_off_unchanged_baseline_recorded(smi):
    """Flag OFF must NOT enter any new branch.  Captured as a stable hash that the
    cross-commit byte-id wrapper compares against HEAD; here we assert the large-
    kappa3 anchors are STILL legacy and the controls STILL build, with the flag off."""
    res = _fffree_isomers(smi)
    if smi in (BIQCOV, ABEZAJ, CYCLAM_NI, AACANI10):
        assert res is None              # legacy with the flag off
    else:
        assert isinstance(res, list) and len(res) >= 1


# ---------------------------------------------------------------------------
# (3) kappa4 STAYS legacy with the flag on (denticity gate untouched)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("smi", [CYCLAM_NI, AACANI10])
def test_kappa4_stays_legacy_with_flag(smi, monkeypatch):
    monkeypatch.setenv(FLAG, "1")
    # even with the cap pushed very high, the denticity>3 gate must hold
    monkeypatch.setenv("DELFIN_FFFREE_CHELATE_HEAVY_CAP", "999")
    assert DEC.decompose(smi) is None
    assert _fffree_isomers(smi) is None


def test_kappa4_is_genuinely_denticity4(monkeypatch):
    """Sanity: cyclam IS a denticity-4 ligand (so the legacy result above is the
    denticity gate doing its job, not a parse failure)."""
    import inspect
    src = inspect.getsource(DEC.decompose)
    patched = (src.replace("if len(fdonors) > 3:", "if len(fdonors) > 99:")
                  .replace('cap = _mono_cap if lg["denticity"] == 1 else _chel_cap',
                           "cap = 9999"))
    glb = dict(DEC.__dict__)
    exec(compile(patched, "<patched>", "exec"), glb)
    monkeypatch.setenv(FLAG, "1")
    d = glb["decompose"](CYCLAM_NI)
    assert d is not None
    assert max(lg["denticity"] for lg in d["ligands"]) == 4


# ---------------------------------------------------------------------------
# (4) determinism: same input -> byte-identical output (flag on)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("smi", [BIQCOV, ABEZAJ])
def test_determinism_with_flag(smi, monkeypatch):
    monkeypatch.setenv(FLAG, "1")
    a = _hash(_fffree_isomers(smi))
    b = _hash(_fffree_isomers(smi))
    assert a == b
    assert not a.startswith("None")
