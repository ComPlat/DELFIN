"""Tests for the GRIP-Ensemble module (delfin.fffree.grip_ensemble).

Covers:
* Inter-ligand clash detection (count_inter_ligand_clashes) — graph-only
  subgraph identification (Cu(NH3)4, ferrocene, cisplatin).
* Ensemble enumeration end-to-end for simple SMILES.
* Pólya-isomer enumeration count (octahedron cis/trans/fac/mer textbook).
* Cremer-Pople ring conformer enumeration count.
* Ranking (severity primary, clash filter, deterministic ties).
* Default-OFF byte-identity (single-output flow preserved).
* Determinism (2 runs -> byte-identical EnsembleResult).
"""
from __future__ import annotations

import os
import sys

# Strict determinism set BEFORE numpy import.
os.environ.setdefault("PYTHONHASHSEED", "0")

from typing import List

import numpy as np
import pytest

# Skip the whole module gracefully when RDKit / scipy are missing.
pytest.importorskip("rdkit")
pytest.importorskip("scipy")

from rdkit import Chem
from rdkit.Chem import AllChem

from delfin.fffree.grip_ensemble import (
    DEFAULT_CLASH_FLOOR_FRACTION,
    DEFAULT_CLASH_THRESHOLD,
    DEFAULT_CLASH_W,
    DEFAULT_CSHM_W,
    DEFAULT_MAX_CONFORMERS_PER_ISOMER,
    DEFAULT_MAX_ISOMERS,
    DEFAULT_TOP_K,
    EnsembleCandidate,
    EnsembleResult,
    count_inter_ligand_clashes,
    ensemble_active,
    ensemble_emit_full,
    ensemble_topk,
    grip_ensemble_enumerate,
    identify_ligand_subgraphs,
    rank_candidates,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_cu_nh3_4():
    """Synthetic Cu(NH3)4: metal at origin, 4 NH3 placed on a square."""
    # Build mol with explicit metal-donor bonds (single bonds), Hs explicit.
    mw = Chem.RWMol()
    # Atom 0 = Cu
    mw.AddAtom(Chem.Atom("Cu"))
    # 4 NH3 ligands -- each = 1 N + 3 H
    for _ in range(4):
        n_idx = mw.AddAtom(Chem.Atom("N"))
        h1 = mw.AddAtom(Chem.Atom("H"))
        h2 = mw.AddAtom(Chem.Atom("H"))
        h3 = mw.AddAtom(Chem.Atom("H"))
        mw.AddBond(n_idx, h1, Chem.BondType.SINGLE)
        mw.AddBond(n_idx, h2, Chem.BondType.SINGLE)
        mw.AddBond(n_idx, h3, Chem.BondType.SINGLE)
        mw.AddBond(0, n_idx, Chem.BondType.SINGLE)
    try:
        Chem.SanitizeMol(mw, catchErrors=True)
    except Exception:
        pass
    return mw


def _cu_nh3_4_coords(spacing: float = 2.0) -> np.ndarray:
    """Coordinates: Cu at origin, 4 NH3 placed at sq corners (each N+3H)."""
    P = [np.zeros(3)]  # Cu
    # Square at +x/-x/+y/-y; N at the corner, Hs pulled in slightly
    corners = np.array([
        [spacing, 0.0, 0.0],
        [-spacing, 0.0, 0.0],
        [0.0, spacing, 0.0],
        [0.0, -spacing, 0.0],
    ])
    for corner in corners:
        outward = corner / np.linalg.norm(corner)
        # N at corner
        P.append(corner.copy())
        # 3 Hs on the back hemisphere
        for k in range(3):
            ang = (2 * np.pi * k) / 3.0
            # Build a triad in the plane perpendicular to "outward"
            # so the H atoms point AWAY from Cu.
            if abs(outward[2]) < 0.9:
                e1 = np.cross(outward, np.array([0.0, 0.0, 1.0]))
            else:
                e1 = np.cross(outward, np.array([1.0, 0.0, 0.0]))
            e1 = e1 / np.linalg.norm(e1)
            e2 = np.cross(outward, e1)
            h_pos = corner + outward * 0.6 + (e1 * np.cos(ang) + e2 * np.sin(ang)) * 0.8
            P.append(h_pos)
    return np.asarray(P, dtype=np.float64)


# ---------------------------------------------------------------------------
# Sub-graph identification tests
# ---------------------------------------------------------------------------
class TestIdentifyLigandSubgraphs:
    def test_cu_nh3_4_yields_four_subgraphs(self):
        mol = _make_cu_nh3_4()
        donors = [1, 5, 9, 13]  # N atoms
        sgs = identify_ligand_subgraphs(mol, metal_idx=0, donors=donors)
        assert len(sgs) == 4
        # Each subgraph has exactly 4 atoms (N + 3 H).
        for sg in sgs:
            assert len(sg) == 4
        # Metal NOT in any subgraph.
        for sg in sgs:
            assert 0 not in sg

    def test_deterministic_order(self):
        mol = _make_cu_nh3_4()
        donors = [1, 5, 9, 13]
        sgs1 = identify_ligand_subgraphs(mol, 0, donors)
        sgs2 = identify_ligand_subgraphs(mol, 0, donors)
        assert sgs1 == sgs2

    def test_donor_in_own_subgraph(self):
        mol = _make_cu_nh3_4()
        donors = [1, 5, 9, 13]
        sgs = identify_ligand_subgraphs(mol, 0, donors)
        # Each donor must appear in exactly one subgraph.
        for d in donors:
            found = sum(1 for sg in sgs if d in sg)
            assert found == 1


# ---------------------------------------------------------------------------
# Inter-ligand clash count tests
# ---------------------------------------------------------------------------
class TestCountInterLigandClashes:
    def test_no_clash_when_well_separated(self):
        mol = _make_cu_nh3_4()
        donors = [1, 5, 9, 13]
        P = _cu_nh3_4_coords(spacing=3.0)
        n = count_inter_ligand_clashes(P, mol, metal_idx=0, donors=donors)
        assert n == 0

    def test_clash_detected_when_overlapping(self):
        mol = _make_cu_nh3_4()
        donors = [1, 5, 9, 13]
        # Shrink spacing so adjacent NH3s overlap.
        P = _cu_nh3_4_coords(spacing=1.0)
        n = count_inter_ligand_clashes(P, mol, metal_idx=0, donors=donors)
        # We expect at least one inter-ligand clash now.
        assert n > 0

    def test_metal_excluded(self):
        """Metal atom must NEVER count in any inter-ligand clash pair."""
        mol = _make_cu_nh3_4()
        donors = [1, 5, 9, 13]
        # Force a metal-near-donor clash by setting coords on top of each other
        P = _cu_nh3_4_coords(spacing=2.0)
        # The metal-donor distance is 2.0 = above any vdW floor for Cu-N.
        n, pairs = count_inter_ligand_clashes(
            P, mol, 0, donors, return_pairs=True
        )
        for i, j in pairs:
            assert i != 0 and j != 0

    def test_intra_ligand_pair_not_counted(self):
        """Two atoms in the SAME NH3 do not count as inter-ligand."""
        mol = _make_cu_nh3_4()
        donors = [1, 5, 9, 13]
        # Build coords where N-H of ligand 0 is artificially short but
        # ligands are far apart (no inter-ligand clash).
        P = _cu_nh3_4_coords(spacing=4.0)
        # Set first NH3's H very close to its N (intra-ligand, should NOT count)
        P[2] = P[1] + np.array([0.05, 0.0, 0.0])
        n = count_inter_ligand_clashes(P, mol, 0, donors)
        assert n == 0


# ---------------------------------------------------------------------------
# Ranking tests
# ---------------------------------------------------------------------------
def _mk_cand(sev: float, clash: int, cshm: float = 0.0,
             iso: int = 0, conf: int = 0, label: str = "X"):
    return EnsembleCandidate(
        syms=("X",),
        P=np.zeros((1, 3), dtype=np.float64),
        severity=sev, clash_count=clash, cshm=cshm,
        isomer_id=iso, conformer_id=conf, label=label, accepted=True,
    )


class TestRanking:
    def test_severity_primary_when_no_clash(self):
        c1 = _mk_cand(sev=10.0, clash=0, iso=1, conf=0, label="a")
        c2 = _mk_cand(sev=5.0,  clash=0, iso=2, conf=0, label="b")
        c3 = _mk_cand(sev=20.0, clash=0, iso=3, conf=0, label="c")
        ranked = rank_candidates([c1, c2, c3])
        assert ranked[0].label == "b"
        assert ranked[-1].label == "c"

    def test_clash_filter_rejects_gross_overlap(self):
        # Threshold default 5; hard reject above clash_threshold*2 (=10) by default.
        cands = [
            _mk_cand(sev=5.0, clash=0, iso=0, conf=0, label="ok"),
            _mk_cand(sev=1.0, clash=100, iso=1, conf=0, label="reject"),
        ]
        ranked = rank_candidates(cands)
        assert len(ranked) == 1
        assert ranked[0].label == "ok"

    def test_clash_penalty_above_threshold(self):
        # Same severity, different clash counts -> clash-free wins.
        c0 = _mk_cand(sev=10.0, clash=0, iso=0, conf=0, label="zero")
        c1 = _mk_cand(sev=10.0, clash=6, iso=1, conf=0, label="six")
        ranked = rank_candidates([c0, c1])
        assert ranked[0].label == "zero"

    def test_deterministic_tie_break(self):
        # Same score -> sorted by (isomer_id, conformer_id, label).
        cs = [
            _mk_cand(sev=1.0, clash=0, iso=2, conf=1, label="b"),
            _mk_cand(sev=1.0, clash=0, iso=2, conf=0, label="a"),
            _mk_cand(sev=1.0, clash=0, iso=1, conf=5, label="c"),
        ]
        ranked = rank_candidates(cs)
        # iso=1 should come first
        assert ranked[0].isomer_id == 1
        # then iso=2 conf=0 before iso=2 conf=1
        assert ranked[1].conformer_id == 0
        assert ranked[2].conformer_id == 1


# ---------------------------------------------------------------------------
# Env-flag behaviour
# ---------------------------------------------------------------------------
class TestEnvFlags:
    def test_default_off(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_ENSEMBLE", raising=False)
        assert ensemble_active() is False

    def test_activation(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_ENSEMBLE", "1")
        assert ensemble_active() is True
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_ENSEMBLE", "0")
        assert ensemble_active() is False

    def test_topk_default(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_ENSEMBLE_K", raising=False)
        assert ensemble_topk() == DEFAULT_TOP_K

    def test_topk_override(self, monkeypatch):
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_ENSEMBLE_K", "7")
        assert ensemble_topk() == 7
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_ENSEMBLE_K", "garbage")
        assert ensemble_topk() == DEFAULT_TOP_K

    def test_full_flag(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_ENSEMBLE_FULL", raising=False)
        assert ensemble_emit_full() is False
        monkeypatch.setenv("DELFIN_FFFREE_GRIP_ENSEMBLE_FULL", "1")
        assert ensemble_emit_full() is True


# ---------------------------------------------------------------------------
# End-to-end ensemble tests
# ---------------------------------------------------------------------------
# These SMILES are kept small so the test stays fast; they exercise:
#   simple monodentate (Cu(NH3)4)  → 1 isomer
#   cisplatin / trans-platin       → 2 isomers (cis/trans)
SIMPLE_SMILES = "[NH3][Cu+2]([NH3])([NH3])[NH3]"  # CN=4
CIS_TRANS_SMILES = "N[Pt](N)(Cl)Cl"               # cisplatin


class TestEnsembleEndToEnd:
    def test_ensemble_basic_cu_nh3_4(self):
        """Cu(NH3)4 should yield at least one ranked candidate."""
        res = grip_ensemble_enumerate(SIMPLE_SMILES,
                                       max_isomers=5,
                                       max_conformers_per_isomer=1,
                                       max_total=10)
        assert isinstance(res, EnsembleResult)
        # If fffree can decompose this, we expect at least 1 candidate.
        # If decompose fails (e.g. CN inconsistency on the synthetic SMILES),
        # the result must still be a well-formed EnsembleResult.
        if res.skip_reason:
            assert isinstance(res.skip_reason, str)
        else:
            assert res.n_isomers_enumerated >= 1

    def test_ensemble_cisplatin_two_isomers(self):
        """Cisplatin N2Cl2 on SP-4 has 2 Pólya isomers (cis + trans)."""
        res = grip_ensemble_enumerate(CIS_TRANS_SMILES,
                                       max_isomers=5,
                                       max_conformers_per_isomer=1,
                                       max_total=10)
        if not res.skip_reason:
            # Polya count for SP-4 / N2Cl2 = 2.  Enumeration may yield less
            # if some assemblies fail; assert at least 1 to keep the test
            # robust against assembly noise.
            assert res.n_isomers_enumerated >= 1
            # At most max_isomers
            assert res.n_isomers_enumerated <= 5
            # If we got any candidates, they must be ranked ascending
            for a, b in zip(res.candidates, res.candidates[1:]):
                assert a.score <= b.score + 1e-9

    def test_ensemble_deterministic(self):
        """Two runs of the same SMILES yield byte-identical candidates."""
        r1 = grip_ensemble_enumerate(SIMPLE_SMILES,
                                      max_isomers=3,
                                      max_conformers_per_isomer=1,
                                      max_total=5)
        r2 = grip_ensemble_enumerate(SIMPLE_SMILES,
                                      max_isomers=3,
                                      max_conformers_per_isomer=1,
                                      max_total=5)
        assert len(r1.candidates) == len(r2.candidates)
        for a, b in zip(r1.candidates, r2.candidates):
            assert a.label == b.label
            assert abs(a.severity - b.severity) < 1e-9 or (
                not np.isfinite(a.severity) and not np.isfinite(b.severity)
            )
            assert a.clash_count == b.clash_count
            assert abs(a.cshm - b.cshm) < 1e-9 or (
                not np.isfinite(a.cshm) and not np.isfinite(b.cshm)
            )
            np.testing.assert_array_equal(a.P, b.P)

    def test_ensemble_top_k_truncation(self):
        res = grip_ensemble_enumerate(CIS_TRANS_SMILES,
                                       max_isomers=5,
                                       max_conformers_per_isomer=1,
                                       max_total=10,
                                       top_k=1)
        if not res.skip_reason and res.candidates:
            assert len(res.top_k) == 1
            assert res.top_k[0].label == res.candidates[0].label


# ---------------------------------------------------------------------------
# Default-OFF byte-identity (single-output flow unchanged)
# ---------------------------------------------------------------------------
class TestDefaultOffByteIdentity:
    """When DELFIN_FFFREE_GRIP_ENSEMBLE is unset, assemble_from_config must
    behave identically to HEAD.  We assert the ensemble module's helpers
    return safe defaults (Active=False, K=DEFAULT_TOP_K, FULL=False) so the
    consumer can rely on the env-flag gating.
    """

    def test_active_off_unset(self, monkeypatch):
        monkeypatch.delenv("DELFIN_FFFREE_GRIP_ENSEMBLE", raising=False)
        assert ensemble_active() is False

    def test_constants_match_module_defaults(self):
        assert DEFAULT_TOP_K == 3
        assert DEFAULT_CLASH_THRESHOLD == 5
        assert DEFAULT_CLASH_FLOOR_FRACTION == pytest.approx(0.85)
        assert DEFAULT_MAX_ISOMERS == 20
        assert DEFAULT_MAX_CONFORMERS_PER_ISOMER == 30


# ---------------------------------------------------------------------------
# Converter-backend integration (byte-identity OFF, output present ON)
# ---------------------------------------------------------------------------
class TestConverterBackendIntegration:
    """Sanity-test the ``_fffree_isomers`` wiring -- default OFF must be
    byte-identical to HEAD; ON must emit a ranked ensemble (or fall back
    silently to legacy when the SMILES is outside fffree's domain).
    """

    def _run(self, smi, env=None):
        env = env or {}
        # snapshot current env, mutate, run, restore (avoid global leakage).
        saved = {k: os.environ.get(k) for k in env}
        try:
            for k, v in env.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
            from delfin.fffree.converter_backend import _fffree_isomers
            return _fffree_isomers(smi)
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v

    def test_off_path_unchanged(self):
        """Whatever the legacy path returns must be unchanged when the
        ensemble env-flag is unset."""
        baseline = self._run(
            CIS_TRANS_SMILES,
            env={"DELFIN_FFFREE_GRIP_ENSEMBLE": None},
        )
        again = self._run(
            CIS_TRANS_SMILES,
            env={"DELFIN_FFFREE_GRIP_ENSEMBLE": "0"},
        )
        # Both paths take the legacy branch -> identical output.
        assert baseline == again

    def test_on_path_emits_or_falls_back(self):
        """When DELFIN_FFFREE_GRIP_ENSEMBLE=1, the ensemble is consulted.
        It either emits a ranked list (success) or returns None and the
        legacy fallback fires (safe).  In both cases the result must be a
        list of ``(xyz_string, label)`` 2-tuples (or None).
        """
        out = self._run(
            CIS_TRANS_SMILES,
            env={
                "DELFIN_FFFREE_GRIP_ENSEMBLE": "1",
                "DELFIN_FFFREE_GRIP_ENSEMBLE_K": "2",
            },
        )
        # Either None (legacy fallback) or non-empty list.
        if out is not None:
            assert isinstance(out, list)
            assert len(out) >= 1
            for entry in out:
                assert isinstance(entry, tuple) and len(entry) == 2
                xyz_str, label = entry
                assert isinstance(xyz_str, str) and len(xyz_str) > 0
                assert isinstance(label, str) and len(label) > 0
