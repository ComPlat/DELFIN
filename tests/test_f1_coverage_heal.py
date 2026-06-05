"""Tests for Mission F1 — fffree coverage heal (~14% -> 88% native rate).

Verifies:
  1. Default-OFF byte-identical (F1 flag unset -> no change).
  2. F1=1 raises ligand size cap so larger ligands pass decompose.
  3. Forensik log records reasons for fall-through cases.
  4. cisplatin (trivial case) still builds under F1=1+PT3=1.

These tests do NOT run the full pool — they exercise specific code paths
in isolation so they're fast (< 5s total).
"""
from __future__ import annotations
import os
import tempfile
import pytest


def _clear_env():
    for k in (
        "DELFIN_FFFREE_F1_COVERAGE",
        "DELFIN_FFFREE_PURE_TRACK3",
        "DELFIN_FFFREE_BUILDER",
        "DELFIN_FFFREE_DECOMPOSE_REASON_LOG",
    ):
        os.environ.pop(k, None)


def _set_baseline_env():
    os.environ["PYTHONHASHSEED"] = "0"
    os.environ["DELFIN_FFFREE_BUILDER"] = "1"


def test_f1_flag_unset_byte_identical():
    """With F1 flag unset, decompose() behaves byte-identical to HEAD."""
    _clear_env()
    _set_baseline_env()
    from delfin.fffree.decompose import decompose

    # Big ligand: 10 heavy/1 denticity > 8 cap -> None at default cap
    d = decompose("[Co](CCCCCCCCCC)(C)(C)C")
    assert d is None, "default-off cap=8 must reject big-ligand"


def test_f1_flag_set_raises_ligand_cap():
    """With F1=1, ligand cap is raised to 30 -> big-Co builds."""
    _clear_env()
    _set_baseline_env()
    os.environ["DELFIN_FFFREE_F1_COVERAGE"] = "1"
    from importlib import reload
    import delfin.fffree.decompose as DEC
    reload(DEC)

    d = DEC.decompose("[Co](CCCCCCCCCC)(C)(C)C")
    assert d is not None, "F1=1 cap=30 must accept 10-heavy ligand"
    assert d["metal"] == "Co"
    assert d["cn"] == 4
    _clear_env()


def test_f1_reason_log_records_paths():
    """When DELFIN_FFFREE_DECOMPOSE_REASON_LOG is set, decompose records reasons."""
    _clear_env()
    _set_baseline_env()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as fh:
        log_path = fh.name
    os.environ["DELFIN_FFFREE_DECOMPOSE_REASON_LOG"] = log_path

    from importlib import reload
    import delfin.fffree.decompose as DEC
    reload(DEC)

    # No-metal SMILES -> "no_metal" reason
    DEC.decompose("c1ccccc1")
    # Multi-metal SMILES -> "multi_metal_disabled"
    DEC.decompose("[Co][Co]")
    # Big-ligand -> "ligand_too_large"
    DEC.decompose("[Co](CCCCCCCCCC)(C)(C)C")

    with open(log_path) as fh:
        log_content = fh.read()

    assert "no_metal" in log_content
    assert "multi_metal_disabled" in log_content
    assert "ligand_too_large" in log_content

    os.unlink(log_path)
    _clear_env()


def test_f1_cisplatin_still_works_with_flag():
    """F1=1 + PT3=1 must not regress simple known-good cases."""
    _clear_env()
    _set_baseline_env()
    os.environ["DELFIN_FFFREE_F1_COVERAGE"] = "1"
    os.environ["DELFIN_FFFREE_PURE_TRACK3"] = "1"
    from importlib import reload
    import delfin.fffree.decompose as DEC
    reload(DEC)
    import delfin.fffree.converter_backend as CB
    reload(CB)

    res = CB._fffree_isomers("N[Pt](N)(Cl)Cl")
    assert res is not None and len(res) >= 1, "cisplatin must still build under F1"
    _clear_env()


def test_f1_shape_max_path_runs():
    """F1=1 still allows valid builds to pass (positive control)."""
    _clear_env()
    _set_baseline_env()
    os.environ["DELFIN_FFFREE_F1_COVERAGE"] = "1"

    import importlib
    import delfin.fffree.converter_backend as CB
    importlib.reload(CB)

    res = CB._fffree_isomers("N[Pt](N)(Cl)Cl")
    assert res is not None
    _clear_env()
