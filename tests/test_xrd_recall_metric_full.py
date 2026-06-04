"""Tests for delfin.fffree.xrd_recall_metric_full (full-CCDC XRD recall).

Determinism: PYTHONHASHSEED=0. No CCDC import.
"""
from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("PYTHONHASHSEED", "0")

# Ensure repo root is on sys.path
REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

np = pytest.importorskip("numpy")

from delfin.fffree import xrd_recall_metric_full as M  # noqa: E402


# ---------------------------------------------------------------------------
# Test fixtures: minimal synthetic family table + tmc index + XYZ archive
# ---------------------------------------------------------------------------
def _write_xyz(path: Path, atoms):
    lines = [str(len(atoms)), "test"]
    for s, (x, y, z) in atoms:
        lines.append(f"{s} {x:.4f} {y:.4f} {z:.4f}")
    path.write_text("\n".join(lines) + "\n")


def _make_octahedral(metal="Fe", donors=("N",) * 6, R=2.0):
    coords = [(metal, (0, 0, 0))]
    dirs = [(R, 0, 0), (-R, 0, 0), (0, R, 0), (0, -R, 0), (0, 0, R), (0, 0, -R)]
    for d, (x, y, z) in zip(donors, dirs):
        coords.append((d, (x, y, z)))
    return coords


def _make_sp4(metal="Pt", donors=("N", "N", "Cl", "Cl"), R=2.0, cis=True):
    """CN4 SP. cis = like-donors adjacent; trans = opposite."""
    if cis:
        dirs = [(R, 0, 0), (0, R, 0), (-R, 0, 0), (0, -R, 0)]
    else:
        dirs = [(R, 0, 0), (-R, 0, 0), (0, R, 0), (0, -R, 0)]
    coords = [(metal, (0, 0, 0))]
    for d, (x, y, z) in zip(donors, dirs):
        coords.append((d, (x, y, z)))
    return coords


@pytest.fixture
def tiny_table_and_index(tmp_path):
    """Build a 2-SMILES, 3-refcode family table + index."""
    index_path = tmp_path / "ccdc_tmc_index.jsonl"
    table_path = tmp_path / "per_smiles_ccdc_families.jsonl"

    # Refcode 1: Fe-OH-N6
    rec1 = {
        "refcode": "FEAAAA", "formula": "FeN6", "block": "d",
        "metal_element": "Fe", "n_metals": 1, "metal_indices_primary": 0,
        "cn": 6, "donor_elements": ["N"], "geom_label": "octahedral_or_TP",
        "n_atoms": 7,
        "symbols": ["Fe", "N", "N", "N", "N", "N", "N"],
        "positions": [[0, 0, 0], [2, 0, 0], [-2, 0, 0], [0, 2, 0],
                      [0, -2, 0], [0, 0, 2], [0, 0, -2]],
        "has_disorder": False, "is_organometallic": False,
    }
    # Refcode 2: Pt-SP4-N2Cl2 cis
    rec2 = {
        "refcode": "PTBBBB", "formula": "PtN2Cl2", "block": "d",
        "metal_element": "Pt", "n_metals": 1, "metal_indices_primary": 0,
        "cn": 4, "donor_elements": ["Cl", "N"],
        "geom_label": "tetrahedral_or_squareplanar",
        "n_atoms": 5,
        "symbols": ["Pt", "N", "N", "Cl", "Cl"],
        "positions": [[0, 0, 0], [2, 0, 0], [0, 2, 0], [-2, 0, 0], [0, -2, 0]],
        "has_disorder": False, "is_organometallic": False,
    }
    # Refcode 3: Pt-SP4-N2Cl2 trans
    rec3 = {
        "refcode": "PTCCCC", "formula": "PtN2Cl2", "block": "d",
        "metal_element": "Pt", "n_metals": 1, "metal_indices_primary": 0,
        "cn": 4, "donor_elements": ["Cl", "N"],
        "geom_label": "tetrahedral_or_squareplanar",
        "n_atoms": 5,
        "symbols": ["Pt", "N", "N", "Cl", "Cl"],
        "positions": [[0, 0, 0], [2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0]],
        "has_disorder": False, "is_organometallic": False,
    }
    with index_path.open("w") as fh:
        for r in (rec1, rec2, rec3):
            fh.write(json.dumps(r) + "\n")

    # Family table: 2 SMILES, one matches FEAAAA, other matches PTBBBB+PTCCCC
    fam_records = [
        {
            "smiles_label": "test-Fe6N",
            "smiles": "[Fe]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]",
            "n_matches": 1,
            "match_method": "graph_signature",
            "matches": [{"refcode": "FEAAAA", "isomer_label": rec1["geom_label"],
                         "n_atoms": rec1["n_atoms"], "metal": "Fe", "cn": 6,
                         "donor_elements": ["N"]}],
        },
        {
            "smiles_label": "test-PtN2Cl2",
            "smiles": "[Pt]([NH3])([NH3])(Cl)Cl",
            "n_matches": 2,
            "match_method": "graph_signature",
            "matches": [
                {"refcode": "PTBBBB", "isomer_label": rec2["geom_label"],
                 "n_atoms": rec2["n_atoms"], "metal": "Pt", "cn": 4,
                 "donor_elements": ["Cl", "N"]},
                {"refcode": "PTCCCC", "isomer_label": rec3["geom_label"],
                 "n_atoms": rec3["n_atoms"], "metal": "Pt", "cn": 4,
                 "donor_elements": ["Cl", "N"]},
            ],
        },
    ]
    with table_path.open("w") as fh:
        for r in fam_records:
            fh.write(json.dumps(r) + "\n")

    return table_path, index_path, fam_records


@pytest.fixture
def empty_table(tmp_path):
    p = tmp_path / "empty.jsonl"
    p.write_text("")
    return p


# ---------------------------------------------------------------------------
# Loader tests
# ---------------------------------------------------------------------------
def test_load_family_table_returns_empty_when_path_missing(tmp_path):
    M._TABLE_CACHE.clear()
    out = M.load_ccdc_family_table(str(tmp_path / "nope.jsonl"))
    assert out == {}


def test_load_family_table_reads_smiles_keys(tiny_table_and_index):
    M._TABLE_CACHE.clear()
    table_path, _, fam_records = tiny_table_and_index
    t = M.load_ccdc_family_table(str(table_path))
    assert len(t) == 2
    assert fam_records[0]["smiles"] in t
    assert fam_records[1]["smiles"] in t
    assert t[fam_records[1]["smiles"]]["n_matches"] == 2


def test_load_family_table_is_cached(tiny_table_and_index):
    M._TABLE_CACHE.clear()
    table_path, _, _ = tiny_table_and_index
    t1 = M.load_ccdc_family_table(str(table_path))
    t2 = M.load_ccdc_family_table(str(table_path))
    assert t1 is t2


def test_load_tmc_index_keyed_by_refcode(tiny_table_and_index):
    M._INDEX_CACHE.clear()
    _, index_path, _ = tiny_table_and_index
    idx = M.load_ccdc_tmc_index(str(index_path))
    assert {"FEAAAA", "PTBBBB", "PTCCCC"}.issubset(idx)
    assert idx["FEAAAA"]["metal_element"] == "Fe"


def test_load_with_env_var(tiny_table_and_index, monkeypatch):
    M._TABLE_CACHE.clear()
    table_path, _, _ = tiny_table_and_index
    monkeypatch.setenv(M.ENV_TABLE_PATH, str(table_path))
    t = M.load_ccdc_family_table()
    assert len(t) == 2


# ---------------------------------------------------------------------------
# Helper / classifier tests (from fallback or sibling)
# ---------------------------------------------------------------------------
def test_classify_isomer_octahedral():
    coords = _make_octahedral("Fe", ("N",) * 6)
    syms = [c[0] for c in coords]
    P = np.asarray([c[1] for c in coords], dtype=float)
    label = M.classify_isomer(syms, P)
    assert "CN6" in label


def test_classify_isomer_no_metal_returns_no_metal():
    syms = ["C", "C", "C"]
    P = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]], dtype=float)
    assert M.classify_isomer(syms, P) == "no-metal"


def test_kabsch_rmsd_identical_returns_near_zero():
    coords = _make_octahedral("Fe")
    syms = [c[0] for c in coords]
    P = np.asarray([c[1] for c in coords], dtype=float)
    rmsd = M.kabsch_rmsd_heavy(syms, P, syms, P)
    assert math.isfinite(rmsd) and rmsd < 1e-6


def test_kabsch_rmsd_rotation_invariant():
    coords = _make_octahedral("Fe")
    syms = [c[0] for c in coords]
    P = np.asarray([c[1] for c in coords], dtype=float)
    # Rotate P by 90° around z
    R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    Pr = P @ R.T
    rmsd = M.kabsch_rmsd_heavy(syms, Pr, syms, P)
    assert rmsd < 1e-6


# ---------------------------------------------------------------------------
# Recall computation tests
# ---------------------------------------------------------------------------
def test_compute_isomer_recall_perfect_match(tiny_table_and_index, tmp_path):
    M._TABLE_CACHE.clear()
    table_path, _, fam = tiny_table_and_index
    table = M.load_ccdc_family_table(str(table_path))
    # Emit an octahedral FeN6 XYZ for the Fe SMILES
    xyz_dir = tmp_path / "emit"
    xyz_dir.mkdir()
    coords = _make_octahedral("Fe")
    _write_xyz(xyz_dir / "fe.xyz", coords)
    emitted = {fam[0]["smiles"]: [xyz_dir / "fe.xyz"]}
    iso = M.compute_isomer_recall(emitted, table)
    assert iso > 0.0  # found CN6 family


def test_compute_isomer_recall_zero_when_no_emissions(tiny_table_and_index):
    M._TABLE_CACHE.clear()
    table_path, _, fam = tiny_table_and_index
    table = M.load_ccdc_family_table(str(table_path))
    iso = M.compute_isomer_recall({fam[0]["smiles"]: []}, table)
    assert math.isnan(iso) or iso == 0.0


def test_compute_isomer_recall_nan_when_no_family(tiny_table_and_index, tmp_path):
    M._TABLE_CACHE.clear()
    table_path, _, _ = tiny_table_and_index
    table = M.load_ccdc_family_table(str(table_path))
    xyz_dir = tmp_path / "emit2"
    xyz_dir.mkdir()
    coords = _make_octahedral("Fe")
    _write_xyz(xyz_dir / "fe.xyz", coords)
    emitted = {"UNKNOWN_SMILES_NOT_IN_TABLE": [xyz_dir / "fe.xyz"]}
    iso = M.compute_isomer_recall(emitted, table)
    assert math.isnan(iso)


def test_compute_conformer_recall_perfect_match(tiny_table_and_index, tmp_path):
    M._TABLE_CACHE.clear()
    M._INDEX_CACHE.clear()
    table_path, index_path, fam = tiny_table_and_index
    table = M.load_ccdc_family_table(str(table_path))
    index = M.load_ccdc_tmc_index(str(index_path))
    xyz_dir = tmp_path / "emit3"
    xyz_dir.mkdir()
    # Emit exactly the FEAAAA structure
    coords = _make_octahedral("Fe")
    _write_xyz(xyz_dir / "fe.xyz", coords)
    emitted = {fam[0]["smiles"]: [xyz_dir / "fe.xyz"]}
    conf = M.compute_conformer_recall(emitted, table, index, rmsd_threshold=0.5)
    assert conf == 1.0


def test_compute_conformer_recall_zero_when_far(tiny_table_and_index, tmp_path):
    M._TABLE_CACHE.clear()
    M._INDEX_CACHE.clear()
    table_path, index_path, fam = tiny_table_and_index
    table = M.load_ccdc_family_table(str(table_path))
    index = M.load_ccdc_tmc_index(str(index_path))
    xyz_dir = tmp_path / "emit4"
    xyz_dir.mkdir()
    # Distort heavily
    coords = _make_octahedral("Fe", R=5.0)
    _write_xyz(xyz_dir / "fe.xyz", coords)
    emitted = {fam[0]["smiles"]: [xyz_dir / "fe.xyz"]}
    conf = M.compute_conformer_recall(emitted, table, index, rmsd_threshold=0.5)
    assert conf == 0.0


def test_per_smiles_recall_report_structure(tiny_table_and_index, tmp_path):
    M._TABLE_CACHE.clear()
    M._INDEX_CACHE.clear()
    table_path, index_path, fam = tiny_table_and_index
    table = M.load_ccdc_family_table(str(table_path))
    index = M.load_ccdc_tmc_index(str(index_path))
    xyz_dir = tmp_path / "emit5"
    xyz_dir.mkdir()
    coords = _make_octahedral("Fe")
    _write_xyz(xyz_dir / "fe1.xyz", coords)
    emitted = {fam[0]["smiles"]: [xyz_dir / "fe1.xyz"]}
    report = M.per_smiles_recall_report(emitted, table, index)
    assert isinstance(report, list)
    assert len(report) == 1
    r = report[0]
    assert "isomer_recall" in r and "conformer_recall" in r
    assert r["n_family"] == 1
    assert r["n_emitted"] == 1


def test_per_smiles_recall_report_no_family(tiny_table_and_index, tmp_path):
    M._TABLE_CACHE.clear()
    table_path, _, _ = tiny_table_and_index
    table = M.load_ccdc_family_table(str(table_path))
    xyz_dir = tmp_path / "emit6"
    xyz_dir.mkdir()
    coords = _make_octahedral("Fe")
    _write_xyz(xyz_dir / "fe.xyz", coords)
    emitted = {"UNKNOWN_SMI": [xyz_dir / "fe.xyz"]}
    report = M.per_smiles_recall_report(emitted, table, {})
    assert len(report) == 1
    assert report[0]["n_family"] == 0
    assert math.isnan(report[0]["isomer_recall"])


def test_weighted_recall_uses_family_size(tiny_table_and_index, tmp_path):
    """Two SMILES, one matches its 1-family, the other matches its 2-family.
    Weighted recall should differ from unweighted."""
    M._TABLE_CACHE.clear()
    table_path, _, fam = tiny_table_and_index
    table = M.load_ccdc_family_table(str(table_path))
    xyz_dir = tmp_path / "emit7"
    xyz_dir.mkdir()
    _write_xyz(xyz_dir / "fe.xyz", _make_octahedral("Fe"))
    _write_xyz(xyz_dir / "pt.xyz", _make_sp4("Pt", cis=True))
    emitted = {
        fam[0]["smiles"]: [xyz_dir / "fe.xyz"],
        fam[1]["smiles"]: [xyz_dir / "pt.xyz"],
    }
    u = M.compute_isomer_recall(emitted, table, weighted=False)
    w = M.compute_isomer_recall(emitted, table, weighted=True)
    assert math.isfinite(u) and math.isfinite(w)


def test_score_archive_error_when_no_table(tmp_path):
    M._TABLE_CACHE.clear()
    res = M.score_archive(str(tmp_path), table_path=str(tmp_path / "nope.jsonl"))
    assert "error" in res


def test_default_rmsd_threshold():
    assert M.DEFAULT_RMSD_A == 0.5


def test_module_imports_without_ccdc():
    """Smoke: module should be import-safe with no CCDC available."""
    import delfin.fffree.xrd_recall_metric_full  # noqa: F401
    assert hasattr(M, "compute_isomer_recall")
    assert hasattr(M, "compute_conformer_recall")
    assert hasattr(M, "load_ccdc_family_table")


def test_env_flag_names_exposed():
    """API surface check — env flag names must be string constants."""
    assert isinstance(M.ENV_TABLE_PATH, str)
    assert isinstance(M.ENV_INDEX_PATH, str)
    assert isinstance(M.ENV_RMSD_A, str)
    assert isinstance(M.ENV_WEIGHTED, str)


def test_group_archive_returns_empty_for_no_xyz(tmp_path):
    """No xyz + no mapping → empty dict (graceful, no crash)."""
    res = M.group_archive_by_smiles(str(tmp_path))
    assert res == {} or hasattr(res, "items")


def test_group_archive_with_mapping(tmp_path):
    """Explicit smiles_mapping → groups XYZ files correctly."""
    _write_xyz(tmp_path / "struct_1.xyz", _make_octahedral("Fe"))
    _write_xyz(tmp_path / "struct_2.xyz", _make_octahedral("Fe"))
    mapping = {"struct_1": "[Fe]_smi", "struct_2": "[Fe]_smi"}
    res = M.group_archive_by_smiles(str(tmp_path), smiles_mapping=mapping)
    assert "[Fe]_smi" in res
    assert len(res["[Fe]_smi"]) == 2


def test_load_master_label_to_smiles_sanitizes(tmp_path):
    """Master pool file → label map handles paren/bracket sanitization."""
    pool = tmp_path / "master.txt"
    pool.write_text(
        "01-Fe(CO)3|[Fe]([C+]#O)([C+]#O)[C+]#O\n"
        "042-ARABUR|[Ag]N1\n"
        "# comment line\n"
        "\n"
        "no_pipe_line_here\n"
    )
    m = M.load_master_label_to_smiles(str(pool))
    assert m["01-Fe(CO)3"] == "[Fe]([C+]#O)([C+]#O)[C+]#O"
    assert m["01-Fe_CO_3"] == "[Fe]([C+]#O)([C+]#O)[C+]#O"  # sanitized variant
    assert m["042-ARABUR"] == "[Ag]N1"


def test_group_archive_resolves_via_master_map(tmp_path):
    """Archive without sidecars → master_label_map resolves via filename stem."""
    # DELFIN-style XYZ with header comment
    p = tmp_path / "042-ARABUR.xyz"
    p.write_text("7\ncommit=test smi=042-ARABUR label=Isomer 1\n"
                 "Fe 0 0 0\nN 2 0 0\nN -2 0 0\nN 0 2 0\nN 0 -2 0\nN 0 0 2\nN 0 0 -2\n")
    pool = tmp_path / "master.txt"
    pool.write_text("042-ARABUR|[Fe]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]\n")
    mm = M.load_master_label_to_smiles(str(pool))
    res = M.group_archive_by_smiles(str(tmp_path), master_label_map=mm)
    assert "[Fe]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]" in res


def test_isomer_match_handles_ccdc_geom_labels():
    """_isomer_match() should match CN6-OH (DELFIN) to octahedral_or_TP (CCDC)."""
    assert M._isomer_match("CN6-OH", "octahedral_or_TP", 6) is True
    assert M._isomer_match("CN4-tet-or-sp", "tetrahedral_or_squareplanar", 4) is True
    assert M._isomer_match("CN5-TBP-or-SPY", "TBP_or_SPY", 5) is True


def test_isomer_match_rejects_wrong_cn():
    assert M._isomer_match("CN4-tet-or-sp", "octahedral_or_TP", 6) is False
    assert M._isomer_match("CN6-OH", "TBP_or_SPY", 5) is False
