"""Tests for the cleaned-CCDC-index swap (DELFIN_FFFREE_USE_CLEANED_CCDC).

User directive 2026-06-04:
  "die xyz von ccdc müssen wie cod gesäubert werden in teilen fehlordnung
   drinnen mehrere komplexe in einer xyz gegen ion enthallten usw"

Verifies:
  1. Default (no env flag): byte-identical to legacy (raw index path used)
  2. Env flag ON + explicit clean path -> clean index is loaded
  3. Env flag ON + clean path missing -> falls back to raw index gracefully
  4. Explicit `path=` argument still wins (highest priority)
  5. Cleaned-index records carry audit fields without breaking the metric
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest


def _write_jsonl(path: Path, records):
    path.write_text("".join(json.dumps(r) + "\n" for r in records))


def _raw_records():
    # Two refcodes with their raw (uncleaned, 6-atom) positions.
    return [
        {
            "refcode": "RAWONE",
            "symbols": ["Cu", "N", "N", "Cl", "Cl", "Cl"],
            "positions": [
                [0.0, 0.0, 0.0],
                [1.5, 0.0, 0.0],
                [-1.5, 0.0, 0.0],
                [0.0, 5.0, 0.0],  # counter-ion
                [0.0, 5.5, 0.0],  # counter-ion
                [0.0, 6.0, 0.0],  # counter-ion
            ],
            "cn": 2,
            "metal_element": "Cu",
            "geom_label": "linear",
            "n_atoms": 6,
            "has_disorder": False,
            "is_organometallic": False,
            "block": "d",
            "donor_elements": ["N"],
        },
        {
            "refcode": "RAWTWO",
            "symbols": ["Fe", "O", "O", "O", "H", "H"],
            "positions": [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [-2.0, 0.0, 0.0],
                [0.0, 4.0, 0.0],  # solvent water O
                [0.5, 4.0, 0.0],  # solvent water H
                [-0.5, 4.0, 0.0],  # solvent water H
            ],
            "cn": 2,
            "metal_element": "Fe",
            "geom_label": "linear",
            "n_atoms": 6,
            "has_disorder": False,
            "is_organometallic": False,
            "block": "d",
            "donor_elements": ["O"],
        },
    ]


def _clean_records():
    # Same two refcodes, counter-ions/solvent stripped (3 atoms each).
    return [
        {
            "refcode": "RAWONE",
            "symbols": ["Cu", "N", "N"],
            "positions": [
                [0.0, 0.0, 0.0],
                [1.5, 0.0, 0.0],
                [-1.5, 0.0, 0.0],
            ],
            "cn": 2,
            "metal_element": "Cu",
            "geom_label": "linear",
            "n_atoms": 3,
            "has_disorder": False,
            "is_organometallic": False,
            "block": "d",
            "donor_elements": ["N"],
            # cleaning audit fields
            "n_atoms_raw": 6,
            "n_components": 4,
            "n_components_with_metal": 1,
            "cleaned_atom_drop_pct": 0.5,
            "drop_h": False,
            "cleaned_primary_metal_idx": 0,
        },
        {
            "refcode": "RAWTWO",
            "symbols": ["Fe", "O", "O"],
            "positions": [
                [0.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [-2.0, 0.0, 0.0],
            ],
            "cn": 2,
            "metal_element": "Fe",
            "geom_label": "linear",
            "n_atoms": 3,
            "has_disorder": False,
            "is_organometallic": False,
            "block": "d",
            "donor_elements": ["O"],
            "n_atoms_raw": 6,
            "n_components": 2,
            "n_components_with_metal": 1,
            "cleaned_atom_drop_pct": 0.5,
            "drop_h": False,
            "cleaned_primary_metal_idx": 0,
        },
    ]


@pytest.fixture
def tmp_indices(tmp_path):
    raw_path = tmp_path / "raw.jsonl"
    clean_path = tmp_path / "clean.jsonl"
    _write_jsonl(raw_path, _raw_records())
    _write_jsonl(clean_path, _clean_records())
    return raw_path, clean_path


def _reset_module(monkeypatch, env):
    """Clear caches and ensure env is set as supplied (others cleared)."""
    # Clear env that could leak from caller
    for k in [
        "DELFIN_FFFREE_XRD_RECALL_TABLE_PATH",
        "DELFIN_FFFREE_XRD_RECALL_INDEX_PATH",
        "DELFIN_FFFREE_USE_CLEANED_CCDC",
        "DELFIN_FFFREE_XRD_RECALL_CLEAN_INDEX_PATH",
    ]:
        monkeypatch.delenv(k, raising=False)
    for k, v in env.items():
        monkeypatch.setenv(k, v)
    # Drop module cache to force re-evaluation
    import sys
    if "delfin.fffree.xrd_recall_metric_full" in sys.modules:
        del sys.modules["delfin.fffree.xrd_recall_metric_full"]


def test_default_off_returns_empty(monkeypatch):
    """No env flags set -> load_ccdc_tmc_index() returns {} (byte-identical legacy)."""
    _reset_module(monkeypatch, {})
    from delfin.fffree import xrd_recall_metric_full as M
    out = M.load_ccdc_tmc_index()
    assert out == {}


def test_explicit_path_argument_wins(monkeypatch, tmp_indices):
    """Explicit path= argument is highest priority (regardless of env)."""
    raw_path, clean_path = tmp_indices
    _reset_module(monkeypatch, {
        "DELFIN_FFFREE_USE_CLEANED_CCDC": "1",
        "DELFIN_FFFREE_XRD_RECALL_CLEAN_INDEX_PATH": str(clean_path),
    })
    from delfin.fffree import xrd_recall_metric_full as M
    out = M.load_ccdc_tmc_index(path=str(raw_path))
    # Should have loaded RAW, NOT CLEAN (explicit path wins)
    assert "RAWONE" in out
    assert out["RAWONE"]["n_atoms"] == 6  # raw atom count
    assert "n_atoms_raw" not in out["RAWONE"]  # not the cleaned record


def test_legacy_raw_path_via_env(monkeypatch, tmp_indices):
    """Legacy: only INDEX_PATH set -> raw is loaded."""
    raw_path, _ = tmp_indices
    _reset_module(monkeypatch, {
        "DELFIN_FFFREE_XRD_RECALL_INDEX_PATH": str(raw_path),
    })
    from delfin.fffree import xrd_recall_metric_full as M
    out = M.load_ccdc_tmc_index()
    assert "RAWONE" in out
    assert out["RAWONE"]["n_atoms"] == 6  # raw


def test_cleaned_swap_on_via_env(monkeypatch, tmp_indices):
    """Flag ON + clean path set -> clean index is preferred."""
    raw_path, clean_path = tmp_indices
    _reset_module(monkeypatch, {
        "DELFIN_FFFREE_USE_CLEANED_CCDC": "1",
        "DELFIN_FFFREE_XRD_RECALL_CLEAN_INDEX_PATH": str(clean_path),
        "DELFIN_FFFREE_XRD_RECALL_INDEX_PATH": str(raw_path),
    })
    from delfin.fffree import xrd_recall_metric_full as M
    out = M.load_ccdc_tmc_index()
    assert "RAWONE" in out
    # Should be CLEAN (3 atoms, has audit fields)
    assert out["RAWONE"]["n_atoms"] == 3
    assert out["RAWONE"]["n_atoms_raw"] == 6
    assert out["RAWONE"]["cleaned_atom_drop_pct"] == 0.5


def test_cleaned_swap_missing_falls_back(monkeypatch, tmp_indices, tmp_path):
    """Flag ON but clean path absent on disk -> falls back to raw."""
    raw_path, _ = tmp_indices
    missing_clean = tmp_path / "no_such_clean.jsonl"
    assert not missing_clean.exists()
    _reset_module(monkeypatch, {
        "DELFIN_FFFREE_USE_CLEANED_CCDC": "1",
        "DELFIN_FFFREE_XRD_RECALL_CLEAN_INDEX_PATH": str(missing_clean),
        "DELFIN_FFFREE_XRD_RECALL_INDEX_PATH": str(raw_path),
    })
    from delfin.fffree import xrd_recall_metric_full as M
    out = M.load_ccdc_tmc_index()
    # Should fall back to raw
    assert "RAWONE" in out
    assert out["RAWONE"]["n_atoms"] == 6  # raw


def test_cleaned_records_have_audit_fields(monkeypatch, tmp_indices):
    """Smoke test: ensure cleaned records preserve audit fields the cleaner emitted."""
    _, clean_path = tmp_indices
    _reset_module(monkeypatch, {})
    from delfin.fffree import xrd_recall_metric_full as M
    out = M.load_ccdc_tmc_index(path=str(clean_path))
    rec = out["RAWONE"]
    for k in ("n_atoms_raw", "n_components", "n_components_with_metal",
              "cleaned_atom_drop_pct", "drop_h", "cleaned_primary_metal_idx"):
        assert k in rec, f"audit field {k} missing from cleaned record"


def test_constants_exported():
    """ENV_USE_CLEANED + DEFAULT_CLEAN_INDEX_PATH must be exported."""
    from delfin.fffree import xrd_recall_metric_full as M
    assert hasattr(M, "ENV_USE_CLEANED")
    assert M.ENV_USE_CLEANED == "DELFIN_FFFREE_USE_CLEANED_CCDC"
    assert hasattr(M, "ENV_CLEAN_INDEX_PATH")
    assert hasattr(M, "DEFAULT_CLEAN_INDEX_PATH")
    assert "ccdc_tmc_index_cleaned.jsonl" in M.DEFAULT_CLEAN_INDEX_PATH
