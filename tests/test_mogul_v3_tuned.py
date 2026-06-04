"""Mission 1 — Tests for delfin.fffree.mogul_detector_v3_tuned.

These tests focus on the per-class threshold filter behavior — they do
NOT depend on the fragment-index existing (the underlying v3 detector
falls back to an empty index in that case, returning no records).
"""
from __future__ import annotations

import csv
import os
import textwrap
from pathlib import Path

import numpy as np
import pytest

import delfin.fffree.mogul_detector_v3_tuned as tuned
import delfin.fffree.mogul_detector_v3 as v3


# ------------------ load_threshold_table ------------------
def test_load_table_missing_file(tmp_path):
    p = tmp_path / "nope.csv"
    out = tuned.load_threshold_table(str(p), force_reload=True)
    assert out == {}


def test_load_table_simple(tmp_path):
    csv_path = tmp_path / "thr.csv"
    csv_path.write_text(textwrap.dedent("""\
        axis,class,n_observations,optimal_threshold,f1
        bonds,bonds::C-N,42,2.0,0.6
        angles,angles::C-N-O,30,1.5,0.5
        torsions,torsions::C-C-N-O,12,3.0,0.4
    """))
    out = tuned.load_threshold_table(str(csv_path), force_reload=True)
    assert out["bonds::C-N"] == 2.0
    assert out["angles::C-N-O"] == 1.5
    assert out["torsions::C-C-N-O"] == 3.0


def test_load_table_handles_bare_class(tmp_path):
    csv_path = tmp_path / "thr.csv"
    csv_path.write_text(textwrap.dedent("""\
        axis,class,n_observations,optimal_threshold,f1
        bonds,C-N,42,2.0,0.6
    """))
    out = tuned.load_threshold_table(str(csv_path), force_reload=True)
    # bare-class form: rebuilt as axis::class
    assert out.get("bonds::C-N") == 2.0


def test_load_table_skips_bad_rows(tmp_path):
    csv_path = tmp_path / "thr.csv"
    csv_path.write_text(textwrap.dedent("""\
        axis,class,optimal_threshold
        bonds,bonds::C-N,not-a-float
        bonds,bonds::C-O,3.0
        ,bonds::no-axis,2.0
        bonds,,2.0
    """))
    out = tuned.load_threshold_table(str(csv_path), force_reload=True)
    assert out == {"bonds::C-O": 3.0}


# ------------------ class_key_for_record ------------------
def test_class_key_simple_bond():
    rec = {"axis": "bond", "atom_syms": ["C", "N"]}
    assert tuned.class_key_for_record(rec) == "bonds::C-N"


def test_class_key_angle_sorted():
    rec = {"axis": "angle", "atom_syms": ["O", "C", "N"]}
    assert tuned.class_key_for_record(rec) == "angles::C-N-O"


def test_class_key_torsion():
    rec = {"axis": "torsion", "atom_syms": ["N", "C", "O", "C"]}
    # sorted: C-C-N-O
    assert tuned.class_key_for_record(rec) == "torsions::C-C-N-O"


def test_class_key_handles_aliases():
    assert (tuned.class_key_for_record({"axis": "pooled_bond", "atom_syms": ["C","C"]})
            == "bonds::C-C")
    assert (tuned.class_key_for_record({"axis": "improper", "atom_syms": ["C","H","N","O"]})
            == "torsions::C-H-N-O")


# ------------------ env-flag gating ------------------
def test_env_flag_off_passthrough(tmp_path, monkeypatch):
    """When DELFIN_MOGUL_V3_TUNED is OFF and no explicit table_path,
    behavior must be identical to the underlying detector."""
    monkeypatch.delenv("DELFIN_MOGUL_V3_TUNED", raising=False)
    syms = ["C", "C", "H", "H", "H", "H"]
    P = np.array([
        [0, 0, 0], [1.5, 0, 0], [-0.5, 0.9, 0],
        [-0.5, -0.9, 0], [2.0, 0.9, 0], [2.0, -0.9, 0],
    ], dtype=float)
    a = tuned.detect_anomalies_v3_tuned(syms, P)
    b = v3.detect_anomalies_v3(syms, P)
    assert a == b


def test_env_flag_on_returns_list(tmp_path, monkeypatch):
    """When flag is ON, the tuned filter is applied. We can't assert
    specific anomalies without the fragment index loaded; we just verify
    the function returns a list and doesn't crash."""
    monkeypatch.setenv("DELFIN_MOGUL_V3_TUNED", "1")
    csv_path = tmp_path / "thr.csv"
    csv_path.write_text("axis,class,optimal_threshold\nbonds,bonds::C-C,2.0\n")
    syms = ["C", "C", "H", "H", "H", "H"]
    P = np.array([
        [0, 0, 0], [1.5, 0, 0], [-0.5, 0.9, 0],
        [-0.5, -0.9, 0], [2.0, 0.9, 0], [2.0, -0.9, 0],
    ], dtype=float)
    out = tuned.detect_anomalies_v3_tuned(syms, P, table_path=str(csv_path))
    assert isinstance(out, list)


def test_explicit_table_path_activates_filter(tmp_path):
    """Even with the flag off, passing an explicit table_path activates
    the filter (test-mode use)."""
    csv_path = tmp_path / "thr.csv"
    csv_path.write_text("axis,class,optimal_threshold\nbonds,bonds::C-C,2.0\n")
    syms = ["C", "C", "H", "H", "H", "H"]
    P = np.array([
        [0, 0, 0], [1.5, 0, 0], [-0.5, 0.9, 0],
        [-0.5, -0.9, 0], [2.0, 0.9, 0], [2.0, -0.9, 0],
    ], dtype=float)
    out = tuned.detect_anomalies_v3_tuned(syms, P, table_path=str(csv_path))
    assert isinstance(out, list)


# ------------------ filter mechanics ------------------
def test_filter_drops_sub_threshold_records(tmp_path):
    """Stub a record stream and verify that records with sev_mad < class
    threshold get dropped while higher-sev ones pass."""
    # Build a fake table {bonds::C-N: 3.0}
    csv_path = tmp_path / "thr.csv"
    csv_path.write_text("axis,class,optimal_threshold\nbonds,bonds::C-N,3.0\n")
    tbl = tuned.load_threshold_table(str(csv_path), force_reload=True)
    assert tbl["bonds::C-N"] == 3.0

    # Simulate the filter logic directly
    recs = [
        {"axis": "bond", "atom_syms": ["C", "N"], "sev_mad": 2.7, "atoms": [0, 1]},
        {"axis": "bond", "atom_syms": ["C", "N"], "sev_mad": 3.5, "atoms": [2, 3]},
        # unknown class — uses fallback 2.5
        {"axis": "bond", "atom_syms": ["C", "O"], "sev_mad": 2.6, "atoms": [4, 5]},
        {"axis": "bond", "atom_syms": ["C", "O"], "sev_mad": 2.3, "atoms": [6, 7]},
    ]
    kept = []
    for rec in recs:
        key = tuned.class_key_for_record(rec)
        thr = tbl.get(key, v3.CCDC_MAD_THRESHOLD)
        if rec["sev_mad"] >= thr:
            kept.append(rec)
    # 3.5 keeps (>=3.0), 2.7 drops, 2.6 keeps (>=2.5 fallback), 2.3 drops
    assert len(kept) == 2
    assert kept[0]["sev_mad"] == 3.5
    assert kept[1]["sev_mad"] == 2.6


def test_env_flag_constant():
    assert tuned.ENV_FLAG == "DELFIN_MOGUL_V3_TUNED"
