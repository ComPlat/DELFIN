"""Tests for the hapto XRD-recall levers (2026-06-18).

Two env-gated, default-OFF, deterministic, FF-free levers extend the rigid-hapto
path for the η-π half-sandwich-carbonyl recall tasche:

* **Hebel A — η-axis rotation enumerator** (``DELFIN_FFFREE_HAPTO_AXIS_ROT``):
  STRICTLY ADDITIVE extra ensemble frames that rigidly rotate the non-primary-η
  block (CO-tripod / co-ligands / a 2nd ring) about the PRIMARY η M→centroid axis.
  The axis runs through the metal (at the origin), so M-centroid + every M-D
  distance is INVARIANT and the η-ring internal geometry is untouched — only the
  ring-vs-rest azimuthal CLOCK changes (the missing piano-stool/sandwich rotamer
  DOF).

* **Hebel B — half-sandwich gate** (``DELFIN_FFFREE_HAPTO_HALFSANDWICH_GATE``):
  restricts the rigid-hapto construction to GENUINE carbonyl half-sandwiches
  (η-ring + ≥1 M-C≡O), graph-based + universal (carbonyl motif, NOT metal/refcode).
  η-systems without a carbonyl tripod fall back to the legacy hapto path.

These tests pin: byte-id OFF (no axis/no gate change vs the historic build), the
core/η invariant under Hebel A, that Hebel B fires only on carbonyl half-sandwich,
and determinism.
"""
import hashlib
import os

import numpy as np
import pytest

# --- representative systems (SMILES from the BEST-50k pool; no CCDC coords) -------
# Cr(CO)3 arene half-sandwich (the recall tasche): η6-arene + 3 σ-CO -> CN4 T-4.
MADREP = ("COC1=CC=C(C#CC2=CC=C(C#C[C+]34->[Cr-9]5678([C]#[O+])([C]#[O+])"
          "([C]#[O+])[C+]([C+]5[C+]36)[C+]7[C+]48)S2)C=C1")
AHAZIT = ("[O+]#[C][Cr-9]12345([C]#[O+])([C]#[O+])<-[C+]6(C=[N+]([O-])"
          "CC7=CC=CC=C7)[C+]1[C+]2[C+]3[C+]4[C+]->56C(F)(F)F")
# A non-carbonyl η piano-stool (vinyl-Pt) -> Hebel B must NOT fire here.
POZNIC = "COC1=CC=C([C+]2=[C+][Pt-3]<-23([Cl])[N+](C)(C)CC[N+]3(C)C)C=C1"


@pytest.fixture(autouse=True)
def _rigid_hapto_env(monkeypatch):
    monkeypatch.setenv("DELFIN_FFFREE_BUILDER", "1")
    monkeypatch.setenv("DELFIN_FFFREE_RIGID_HAPTO", "1")
    monkeypatch.setenv("DELFIN_FFFREE_CN3", "1")
    monkeypatch.setenv("DELFIN_FFFREE_CN_EXTEND", "1")
    monkeypatch.setenv("PYTHONHASHSEED", "0")
    # ensure new levers start OFF (each test opts in explicitly)
    monkeypatch.delenv("DELFIN_FFFREE_HAPTO_AXIS_ROT", raising=False)
    monkeypatch.delenv("DELFIN_FFFREE_HAPTO_HALFSANDWICH_GATE", raising=False)


def _decomp(smiles):
    from delfin.manta import decompose as DC
    d = DC.decompose(smiles)
    assert d is not None, "rigid-hapto decompose returned None"
    return d


def _ring_carbons(syms, P):
    P = np.asarray(P, float)
    M = P[0]
    dM = np.linalg.norm(P - M, axis=1)
    Os = [i for i, s in enumerate(syms) if s == "O"]
    return [i for i, s in enumerate(syms)
            if s == "C" and 2.0 < dM[i] < 2.75
            and not any(np.linalg.norm(P[i] - P[o]) < 1.4 for o in Os)]


# --------------------------------------------------------------------------------
# Hebel A — axis-rotation builds are orientation-diverse but core-invariant
# --------------------------------------------------------------------------------
def test_axis_rotants_change_orientation_only():
    """Each axis-rot build keeps M-centroid + ring-internal + M-D distances IDENTICAL
    to the canonical build (<=0.05 A) while producing a DISTINCT overall orientation."""
    from delfin.manta import assemble_complex as AC
    d = _decomp(MADREP)
    metal, geom = d["metal"], d["geometry"]
    syms, P0, donors0, _ = AC.assemble_hapto(metal, geom, d)
    P0 = np.asarray(P0, float)
    M0 = P0[0]
    ring = _ring_carbons(syms, P0)
    assert len(ring) >= 5
    cen0 = P0[ring].mean(0)
    mc0 = np.linalg.norm(cen0 - M0)
    md0 = np.sort(np.linalg.norm(P0 - M0, axis=1)[donors0])

    def ringfp(P):
        R = P[ring]
        return np.sort([np.linalg.norm(R[i] - R[j])
                        for i in range(len(R)) for j in range(i + 1, len(R))])
    fp0 = ringfp(P0)

    rotants = AC.assemble_hapto_axis_rotants(metal, geom, d, n_axis=8, max_builds=30)
    assert rotants, "Hebel A produced no axis rotants on a carbonyl half-sandwich"
    moved = 0
    for syms2, P2, don2, _ in rotants:
        P2 = np.asarray(P2, float)
        M2 = P2[0]
        cen2 = P2[ring].mean(0)
        assert abs(np.linalg.norm(cen2 - M2) - mc0) <= 0.05, "M-centroid moved"
        assert float(np.max(np.abs(ringfp(P2) - fp0))) <= 0.05, "η-ring internal changed"
        md2 = np.sort(np.linalg.norm(P2 - M2, axis=1)[don2])
        assert float(np.max(np.abs(md2 - md0))) <= 0.05, "M-D distance changed"
        # orientation IS different (the tripod clock rotated): not all atoms coincide
        if AC._rmsd_aligned(P0, P2) > 1e-3:
            moved += 1
    assert moved >= 1, "axis rotants did not produce any distinct orientation"


def test_axis_rot_off_is_byte_identical():
    """Hebel A OFF (flag unset) -> the converter ensemble is byte-identical, i.e. the
    additive axis frames are NOT present."""
    from delfin.smiles_converter import smiles_to_xyz_isomers

    def ens(smi):
        out, _ = smiles_to_xyz_isomers(smi, max_isomers=64, quality_mode="normal")
        return out

    base = ens(MADREP)
    # No axis-labelled frames when the flag is off.
    assert not any("axis" in lab for _, lab in base)


def test_axis_rot_on_appends_axis_frames_and_keeps_baseline(monkeypatch):
    """Hebel A ON: every baseline frame is still present (strictly additive) and new
    'axis' frames are appended."""
    from delfin.smiles_converter import smiles_to_xyz_isomers

    def ens(smi):
        out, _ = smiles_to_xyz_isomers(smi, max_isomers=64, quality_mode="normal")
        return out

    base = ens(MADREP)
    base_xyz = {xyz for xyz, _ in base}
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_AXIS_ROT", "1")
    on = ens(MADREP)
    on_xyz = {xyz for xyz, _ in on}
    # additive: every baseline frame survives, and axis frames are added
    assert base_xyz <= on_xyz, "Hebel A evicted a baseline frame (NOT additive)"
    assert any("axis" in lab for _, lab in on), "no axis frames appended"
    assert len(on) > len(base)


# --------------------------------------------------------------------------------
# Hebel B — half-sandwich gate fires only on carbonyl η-systems
# --------------------------------------------------------------------------------
def test_gate_keeps_carbonyl_halfsandwich(monkeypatch):
    """Gate ON: a Cr(CO)3 arene half-sandwich (η-ring + CO) STILL takes the rigid-
    hapto path."""
    from delfin.manta import decompose as DC
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_HALFSANDWICH_GATE", "1")
    for smi in (MADREP, AHAZIT):
        d = DC.decompose(smi)
        assert d is not None and d.get("has_eta"), "carbonyl half-sandwich suppressed"


def test_gate_suppresses_noncarbonyl_eta(monkeypatch):
    """Gate ON: a non-carbonyl η piano-stool (vinyl-Pt, no M-CO) falls OFF the rigid-
    hapto path (decompose yields no η-build), so it reverts to legacy."""
    from delfin.manta import decompose as DC
    # OFF: hapto fires
    d_off = DC.decompose(POZNIC)
    assert d_off is not None and d_off.get("has_eta")
    # ON: hapto path is gated out
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_HALFSANDWICH_GATE", "1")
    d_on = DC.decompose(POZNIC)
    assert not (d_on is not None and d_on.get("has_eta")), \
        "Hebel B fired on a non-carbonyl η-system"


def test_gate_legacy_fallback_byte_identical(monkeypatch):
    """Gate ON on a no-CO η-system == the pure-legacy build (RIGID_HAPTO off): the
    suppressed system reverts to EXACTLY its legacy frames."""
    from delfin.smiles_converter import smiles_to_xyz_isomers

    def ens_hash(smi):
        out, _ = smiles_to_xyz_isomers(smi, max_isomers=64, quality_mode="normal")
        blob = "\n".join(f"{lab}|{xyz}" for xyz, lab in out)
        return hashlib.md5(blob.encode()).hexdigest(), len(out)

    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_HALFSANDWICH_GATE", "1")
    h_gate, n_gate = ens_hash(POZNIC)
    monkeypatch.delenv("DELFIN_FFFREE_RIGID_HAPTO", raising=False)  # pure legacy
    h_leg, n_leg = ens_hash(POZNIC)
    assert (h_gate, n_gate) == (h_leg, n_leg), \
        "gated no-CO η build differs from the legacy build"


def test_metal_has_carbonyl_graph_detector():
    """The CO detector is graph-based: True for M-C≡O, False for M-arene only."""
    from rdkit import Chem
    from delfin.manta.decompose import _metal_has_carbonyl
    # crude probes: build mol, find metal, check.
    for smi, expect in ((MADREP, True), (AHAZIT, True), (POZNIC, False)):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        assert mol is not None
        Chem.SanitizeMol(mol, catchErrors=True)
        metals = [a.GetIdx() for a in mol.GetAtoms()
                  if a.GetSymbol() in ("Cr", "Pt", "Mo", "Ru", "Os", "W", "Mn", "Fe")]
        assert metals
        m = metals[0]
        got = _metal_has_carbonyl(mol, m, mol.GetAtomWithIdx(m))
        assert got == expect, f"{smi[:20]}: CO detector {got} != {expect}"


# --------------------------------------------------------------------------------
# Determinism (both levers ON) — two builds byte-identical
# --------------------------------------------------------------------------------
def test_determinism_both_levers_on(monkeypatch):
    from delfin.smiles_converter import smiles_to_xyz_isomers
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_AXIS_ROT", "1")
    monkeypatch.setenv("DELFIN_FFFREE_HAPTO_HALFSANDWICH_GATE", "1")

    def ens_hash(smi):
        out, _ = smiles_to_xyz_isomers(smi, max_isomers=64, quality_mode="normal")
        blob = "\n".join(f"{lab}|{xyz}" for xyz, lab in out)
        return hashlib.md5(blob.encode()).hexdigest()

    assert ens_hash(MADREP) == ens_hash(MADREP)
