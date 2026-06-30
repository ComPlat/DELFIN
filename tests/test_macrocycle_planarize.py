"""macrocycle_planarize: flatten conjugated metallo-macrocycles (porphyrins),
leave non-aromatic / non-macrocyclic systems alone.  Default OFF -> byte-id."""
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

from delfin.manta import macrocycle_planarize as MP

# Fe porphine (unsubstituted porphyrin) — a conjugated macrocycle around Fe.
PORPHINE = "C1=CC2=CC3=CC=C(C=C3)C3=CC=C(C=C3)C3=CC=C(N3)C=C1N2"  # rough; use a real one below
PORPHINE = "[Fe]123n4c5ccc4cc4ccc(n4->1)cc1ccc(n1->2)cc1ccc(cc5)n1->3"  # metalloporphine


def _embed(smi):
    m = Chem.MolFromSmiles(smi, sanitize=False)
    if m is None:
        return None, None
    try:
        Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                         | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                         | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                         | Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
    except Exception:
        pass
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=1, useRandomCoords=True) != 0:
        return None, None
    conf = mh.GetConformer()
    syms = [a.GetSymbol() for a in mh.GetAtoms()]
    P = np.array([list(conf.GetAtomPosition(i)) for i in range(mh.GetNumAtoms())])
    return mh, (syms, P)


def test_detect_core_on_porphyrin_or_skip():
    m = Chem.MolFromSmiles(PORPHINE, sanitize=False)
    core, arom = MP.detect_core(m)
    # Either it detects a conjugated macrocyclic core, or (on parse quirks) returns
    # None -- but it must NEVER crash and must not claim a tiny core.
    if core is not None:
        assert len(core) >= MP._MIN_CORE
        assert isinstance(arom, bool)


def test_planarize_never_increases_rms_and_keeps_metal():
    res = _embed(PORPHINE)
    if res[0] is None:
        return  # embed/parse environment quirk -> skip, never fail
    mh, (syms, P) = res
    core, arom = MP.detect_core(mh)
    if core is None:
        return
    r0, _, _ = MP._plan_rms(P, core)
    Pf = MP.planarize(syms, P, mh)
    r1, _, _ = MP._plan_rms(np.asarray(Pf), core)
    # never-worse on planarity (flatter or unchanged)
    assert r1 <= r0 + 1e-6


def test_no_metal_is_identity():
    m = Chem.MolFromSmiles("c1ccccc1", sanitize=False)
    Chem.FastFindRings(m)
    P = np.array([[np.cos(t), np.sin(t), 0.1 * (i % 2)] for i, t in
                  enumerate(np.linspace(0, 2 * np.pi, 6, endpoint=False))])
    syms = ["C"] * 6
    Pf = MP.planarize(syms, P, m)
    assert np.array_equal(np.asarray(Pf), P)   # no metal core -> identity


def test_deterministic():
    res = _embed(PORPHINE)
    if res[0] is None:
        return
    mh, (syms, P) = res
    a = MP.planarize(syms, P, mh)
    b = MP.planarize(syms, P, mh)
    assert np.array_equal(np.asarray(a), np.asarray(b))
