"""Tests for the FF-free linear-diatomic-donor orientation guard
(DELFIN_FFFREE_DIATOMIC_ORIENT, delfin.fffree.assemble_complex).

Root cause covered: two chemically-identical diatomic ligands (e.g. two M-C#O
carbonyls) with DIFFERENT fragment atom orderings collide in the canonical-SMILES
conformer cache, so the cached atom ordering no longer matches the second ligand's
donor_local_idx -> the partner atom (O) ends up at the coordination vertex (an
M-O-C isocarbonyl flip).  The guard re-seats the SMILES-bonded donor at the vertex.

Locks the load-bearing contracts:
  * donor element comes from the molecular GRAPH, never hardcoded (carbonyl C,
    cyanide C, nitrosyl N all resolve to the SMILES metal-bonded atom),
  * a FLIPPED diatomic block is reflected so the donor is metal-closest,
  * an already-correct block is an EXACT no-op,
  * the diatomic bond length is preserved exactly (rigid),
  * env-gate: DELFIN_FFFREE_DIATOMIC_ORIENT unset/"0" -> byte-identical build,
  * end-to-end: AXAGOY's flipped carbonyl is C-bound with the flag ON.
"""
import os
import hashlib

import numpy as np
import pytest

from rdkit import Chem

import delfin.fffree.assemble_complex as AC


AXAGOY = ("[O+]#[C][Os-4]1([Cl])([Cl])([C]#[O+])[N+]2=CC=C(C3=CC=CC=C3)"
          "C3=CC=C4C(C5=CC=CC=C5)=CC=[N+]1C4=C32")


def _lg(smiles, donor_local):
    return {"mol": Chem.MolFromSmiles(smiles),
            "donor_local_idxs": [donor_local], "denticity": 1}


def test_donor_partner_is_graph_derived_not_hardcoded():
    # carbonyl C#O: the donor follows whichever heavy atom is metal-bonded.
    assert AC._diatomic_donor_partner(_lg("[C]#[O+]", 0)) == ("C", "O")
    assert AC._diatomic_donor_partner(_lg("[C]#[O+]", 1)) == ("O", "C")
    # cyanide C#N
    assert AC._diatomic_donor_partner(_lg("[C-]#N", 0)) == ("C", "N")
    # nitrosyl N=O: donor is N (NOT a hardcoded carbon).
    assert AC._diatomic_donor_partner(_lg("[N]=O", 0)) == ("N", "O")


def test_donor_partner_skips_homonuclear_and_nondiatomic():
    # homonuclear N#N: orientation is symmetric -> no detectable flip -> None.
    assert AC._diatomic_donor_partner(_lg("N#N", 0)) is None
    # triatomic ligand (not a diatomic) -> None.
    assert AC._diatomic_donor_partner(_lg("[O]C=O", 0)) is None
    # chelate (denticity 2) -> None.
    lg = {"mol": Chem.MolFromSmiles("[C]#[O+]"),
          "donor_local_idxs": [0, 1], "denticity": 2}
    assert AC._diatomic_donor_partner(lg) is None


def test_flipped_block_is_reflected_donor_first():
    metal = np.zeros(3)
    vertex = np.array([2.2, 0.0, 0.0])
    # lsyms cache ordering ['O','C']; FLIPPED: O at the vertex (close), C outward.
    lsyms = ["O", "C"]
    Q = np.array([[2.2, 0.0, 0.0], [3.318, 0.0, 0.0]])
    out = AC._orient_diatomic_block(Q, lsyms, "C", "O", metal, vertex)
    d_C = np.linalg.norm(out[1] - metal)
    d_O = np.linalg.norm(out[0] - metal)
    assert d_C < d_O                                   # donor (C) now metal-closest
    assert d_C == pytest.approx(2.2, abs=1e-6)         # donor exactly on the vertex
    # rigid: the C-O bond length is preserved exactly.
    assert np.linalg.norm(out[0] - out[1]) == pytest.approx(
        np.linalg.norm(Q[0] - Q[1]), abs=1e-9)


def test_already_correct_block_is_exact_noop():
    metal = np.zeros(3)
    vertex = np.array([2.2, 0.0, 0.0])
    lsyms = ["O", "C"]
    Q = np.array([[3.318, 0.0, 0.0], [2.2, 0.0, 0.0]])  # C at vertex = correct
    out = AC._orient_diatomic_block(Q, lsyms, "C", "O", metal, vertex)
    assert np.allclose(out, Q)


def test_nitrosyl_flip_is_corrected_to_n_donor():
    metal = np.zeros(3)
    vertex = np.array([2.2, 0.0, 0.0])
    lsyms = ["O", "N"]
    Q = np.array([[2.2, 0.0, 0.0], [3.3, 0.0, 0.0]])    # O at vertex = FLIPPED
    out = AC._orient_diatomic_block(Q, lsyms, "N", "O", metal, vertex)
    assert np.linalg.norm(out[1] - metal) < np.linalg.norm(out[0] - metal)


def _build_sha(smiles):
    from delfin.fffree.converter_backend import _fffree_isomers
    res = _fffree_isomers(smiles, max_isomers=50)
    if not res:
        return "NONE"
    blob = "\n=ISO=\n".join(r[0] + "||" + str(r[1]) for r in res)
    return hashlib.sha256(blob.encode()).hexdigest()


def _onstack_env():
    """The ON-STACK flag set needed to reach the chelate build path for AXAGOY."""
    return {
        "PYTHONHASHSEED": "0",
        "DELFIN_FFFREE_BUILDER": "1", "DELFIN_FFFREE_RIGID_HAPTO": "1",
        "DELFIN_FFFREE_CN3": "1", "DELFIN_FFFREE_CN_EXTEND": "1",
        "DELFIN_FFFREE_SIGMA_ENSEMBLE": "1", "DELFIN_FFFREE_DONOR_BEND": "1",
        "DELFIN_FFFREE_JOINT_DECLASH": "1", "DELFIN_FFFREE_NHC_CARBENE": "1",
        "DELFIN_FFFREE_BACKBONE_REEMBED": "1", "DELFIN_FFFREE_CONFORMER_SEATING": "1",
        "DELFIN_FFFREE_CHELATE_BACKBONE": "1", "DELFIN_FFFREE_HAPTO_AXIS_ROT": "1",
        "DELFIN_FFFREE_HAPTO_HALFSANDWICH_GATE": "1", "DELFIN_FFFREE_MULTIBOND_EXEMPT": "1",
        "DELFIN_FFFREE_TORSION_RELAX": "1", "DELFIN_FFFREE_INTERLIG_RANK": "1",
        "DELFIN_FFFREE_MD_CONTEXT": "1", "DELFIN_FFFREE_XH_COLLAPSE": "1",
        "DELFIN_FFFREE_INTERLIG_VDW_GATE": "1", "DELFIN_FFFREE_SP2N_PLANARIZE": "1",
        "DELFIN_FFFREE_AROM_PLANARIZE": "1", "DELFIN_FFFREE_ARYL_RING_SIZE": "1",
        "DELFIN_FFFREE_PLANAR_MER": "1", "DELFIN_FFFREE_PLANAR_MER_CN5": "1",
    }


@pytest.mark.parametrize("smiles", ["N[Pt](N)(Cl)Cl", "N[Fe](N)(N)(N)(N)N", AXAGOY])
def test_flag_off_is_byte_identical(monkeypatch, smiles):
    for k, v in _onstack_env().items():
        monkeypatch.setenv(k, v)
    # clear the process-local conformer cache so the build is reproducible here
    AC._CONF_CACHE.clear()
    monkeypatch.delenv("DELFIN_FFFREE_DIATOMIC_ORIENT", raising=False)
    sha_unset = _build_sha(smiles)
    AC._CONF_CACHE.clear()
    monkeypatch.setenv("DELFIN_FFFREE_DIATOMIC_ORIENT", "0")
    sha_zero = _build_sha(smiles)
    assert sha_unset == sha_zero != "NONE" or sha_unset == sha_zero == "NONE"


def _frame_flips(xyz, metal_sym="Os"):
    syms, P = [], []
    for ln in xyz.strip().splitlines():
        p = ln.split()
        if len(p) >= 4:
            syms.append(p[0]); P.append([float(x) for x in p[1:4]])
    P = np.array(P)
    m = [i for i, s in enumerate(syms) if s == metal_sym][0]
    Cs = [i for i, s in enumerate(syms) if s == "C"]
    Os = [i for i, s in enumerate(syms) if s == "O"]
    flips = 0
    for c in Cs:
        for o in Os:
            if np.linalg.norm(P[c] - P[o]) < 1.35:
                mc = np.linalg.norm(P[m] - P[c]); mo = np.linalg.norm(P[m] - P[o])
                if min(mc, mo) < 3.6 and mo < mc:
                    flips += 1
    return flips


def test_axagoy_carbonyl_flip_fixed_on(monkeypatch):
    from delfin.fffree.converter_backend import _fffree_isomers
    for k, v in _onstack_env().items():
        monkeypatch.setenv(k, v)
    AC._CONF_CACHE.clear()
    monkeypatch.setenv("DELFIN_FFFREE_DIATOMIC_ORIENT", "0")
    off = _fffree_isomers(AXAGOY, max_isomers=50)
    assert off, "AXAGOY must build under the ON-STACK flags"
    assert sum(_frame_flips(r[0]) for r in off) > 0   # defect present with the flag off
    AC._CONF_CACHE.clear()
    monkeypatch.setenv("DELFIN_FFFREE_DIATOMIC_ORIENT", "1")
    on = _fffree_isomers(AXAGOY, max_isomers=50)
    assert on
    assert sum(_frame_flips(r[0]) for r in on) == 0   # every carbonyl C-bound with the flag on
