"""tests/test_auto_diagnostic.py — auto-diagnostic pipeline regression tests.

Covers:
  1. Bug-pattern classifier on synthetic OKACOU-like coords → cl_cl_1_3_collapse
  2. Synthetic Co T-4 d^7 (Co not in D8) → does NOT fire d^8 mismatch
  3. Synthetic Pt T-4 (d^8) → fires polyhedron_mismatch_d8_t4
  4. extract_worst_files on a tiny fixture archive → returns expected ordering
  5. ccdc_overlay on a real CCDC mol2 (XACPIC) against the same mol2 → 0.0 RMSD
  6. bug_census on synthetic mini-archive → returns >= 2 distinct classes

Tests construct synthetic XYZ files on the fly into tmp_path so they require
no external data beyond what ships in CCDC/ + COD_all/clean/.

Author: hmaximilian <hmaximilian496@gmail.com>
"""
from __future__ import annotations
import os
import sys
import pathlib
import tempfile

import pytest

SCRIPTS_DIR = pathlib.Path(__file__).parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from diagnostic_common import (  # noqa: E402
    parse_xyz_first_frame, find_metal_idx, coord_sphere,
    parse_mol2_atoms, refcode_from_filename,
)
from classify_bug_pattern import (  # noqa: E402
    classify_atoms, classify_bug, KNOWN_PATTERNS,
)


def _write_xyz(tmp_path: pathlib.Path, name: str,
               atoms: list) -> str:
    """Write a tiny xyz; atoms = list of (elem, (x,y,z))."""
    p = tmp_path / name
    lines = [str(len(atoms)), "synthetic"]
    for e, (x, y, z) in atoms:
        lines.append(f"{e:<3s} {x:12.6f} {y:12.6f} {z:12.6f}")
    p.write_text("\n".join(lines) + "\n")
    return str(p)


# ---------------------------------------------------------------------------
# Sanity checks for the parser + helpers
# ---------------------------------------------------------------------------
def test_parser_roundtrip(tmp_path):
    atoms = [("Fe", (0, 0, 0)), ("Cl", (2.3, 0, 0)), ("Cl", (-2.3, 0, 0)),
             ("Cl", (0, 2.3, 0)), ("Cl", (0, -2.3, 0))]
    p = _write_xyz(tmp_path, "fe_cl4.xyz", atoms)
    parsed = parse_xyz_first_frame(p)
    assert len(parsed) == 5
    mi = find_metal_idx(parsed)
    assert parsed[mi][0] == "Fe"
    donors, _ = coord_sphere(parsed, mi)
    assert len(donors) == 4
    assert {d[1] for d in donors} == {"Cl"}


def test_known_patterns_have_severity_and_test():
    for name, spec in KNOWN_PATTERNS.items():
        assert "severity" in spec
        assert spec["severity"] in {"severe", "moderate", "mild"}
        assert callable(spec["test"])


# ---------------------------------------------------------------------------
# Test 1: cl_cl_1_3_collapse fires on synthetic OKACOU-like geometry
# ---------------------------------------------------------------------------
def test_cl_cl_collapse_fires(tmp_path):
    # Two halides through metal: place them at 0.9 A from each other through
    # the metal axis.  In real OKACOU the collapse was Cd-Cl, but Ru/Fe is fine
    # for the predicate.
    atoms = [
        ("Fe", (0, 0, 0)),
        ("Cl", (1.2, 0, 0)),
        ("Cl", (1.2, 0.9, 0)),  # only 0.9 A apart!  (1-3 through Fe)
        ("Cl", (0, 1.2, 0)),
        ("Cl", (0, -1.2, 0)),
    ]
    hits = classify_atoms(atoms)
    assert "cl_cl_1_3_collapse" in hits, f"expected cl_cl_1_3_collapse in {hits}"


def test_cl_cl_collapse_clean(tmp_path):
    # Spread tetrahedrally: pairwise > 2.5 A.
    s = 1.4
    atoms = [
        ("Fe", (0, 0, 0)),
        ("Cl", (s, s, s)),
        ("Cl", (-s, -s, s)),
        ("Cl", (-s, s, -s)),
        ("Cl", (s, -s, -s)),
    ]
    hits = classify_atoms(atoms)
    assert "cl_cl_1_3_collapse" not in hits


# ---------------------------------------------------------------------------
# Test 2: synthetic Co T-4 d^7 → does NOT fire d^8 mismatch
# ---------------------------------------------------------------------------
def test_co_t4_does_not_fire_d8_mismatch():
    # Co is ambivalent: NOT in D8_METALS, so d8 predicate must remain quiet.
    s = 1.3
    atoms = [
        ("Co", (0, 0, 0)),
        ("N", (s, s, s)),
        ("N", (-s, -s, s)),
        ("N", (-s, s, -s)),
        ("N", (s, -s, -s)),
    ]
    hits = classify_atoms(atoms)
    assert "polyhedron_mismatch_d8_t4" not in hits


# ---------------------------------------------------------------------------
# Test 3: Pt T-4 (d^8) → fires polyhedron_mismatch_d8_t4
# ---------------------------------------------------------------------------
def test_pt_t4_fires_d8_mismatch():
    # Pt is in D8_METALS; geometry is tetrahedral (lower T score vs SP).
    s = 1.5
    atoms = [
        ("Pt", (0, 0, 0)),
        ("Cl", (s, s, s)),
        ("Cl", (-s, -s, s)),
        ("Cl", (-s, s, -s)),
        ("Cl", (s, -s, -s)),
    ]
    hits = classify_atoms(atoms)
    assert "polyhedron_mismatch_d8_t4" in hits, f"got {hits}"


def test_pt_sp4_does_not_fire_d8_mismatch():
    s = 2.3
    atoms = [
        ("Pt", (0, 0, 0)),
        ("Cl", (s, 0, 0)),
        ("Cl", (-s, 0, 0)),
        ("Cl", (0, s, 0)),
        ("Cl", (0, -s, 0)),
    ]
    hits = classify_atoms(atoms)
    assert "polyhedron_mismatch_d8_t4" not in hits


# ---------------------------------------------------------------------------
# Test: d^10 Zn in SP-4 → fires polyhedron_mismatch_d10_sp4
# ---------------------------------------------------------------------------
def test_zn_sp4_fires_d10_mismatch():
    s = 2.2
    atoms = [
        ("Zn", (0, 0, 0)),
        ("N", (s, 0, 0)),
        ("N", (-s, 0, 0)),
        ("N", (0, s, 0)),
        ("N", (0, -s, 0)),
    ]
    hits = classify_atoms(atoms)
    assert "polyhedron_mismatch_d10_sp4" in hits


# ---------------------------------------------------------------------------
# Test 4: ch_collapse / xh_collapse fire when a bond is compressed
# ---------------------------------------------------------------------------
def test_ch_collapse_fires():
    atoms = [
        ("Fe", (0, 0, 0)),
        ("N", (2.2, 0, 0)),
        ("N", (-2.2, 0, 0)),
        ("N", (0, 2.2, 0)),
        ("N", (0, -2.2, 0)),
        ("C", (3.6, 0, 0)),
        ("H", (3.85, 0, 0)),  # C-H = 0.25 A!
    ]
    hits = classify_atoms(atoms)
    assert "ch_collapse" in hits


def test_md_too_short_fires():
    atoms = [
        ("Fe", (0, 0, 0)),
        ("Cl", (1.0, 0, 0)),    # very short
        ("Cl", (-2.3, 0, 0)),
        ("Cl", (0, 2.3, 0)),
    ]
    hits = classify_atoms(atoms)
    assert "md_too_short" in hits


def test_donor_drift_fires():
    # One donor right in the 2.4-3.0 bracket but other donors fine.
    # Use N (covrad 0.71) so cutoff with Fe (1.42)*1.30 ≈ 2.77 A.
    atoms = [
        ("Fe", (0, 0, 0)),
        ("N", (2.0, 0, 0)),
        ("N", (-2.0, 0, 0)),
        ("N", (0, 2.0, 0)),
        ("N", (0, 0, 2.7)),  # in drift bracket
    ]
    hits = classify_atoms(atoms)
    assert "donor_drift" in hits


# ---------------------------------------------------------------------------
# Test 5: extract_worst_files on a tiny fixture archive
# ---------------------------------------------------------------------------
def test_extract_worst_files_basic(tmp_path):
    from extract_worst_files import extract_worst_files
    # 3 fake xyz files: one collapsed, one fine, one drifted
    arch = tmp_path / "tiny_archive"
    arch.mkdir()
    # collapsed: smallest min_md
    _write_xyz(arch, "a-AAAAAA.xyz", [
        ("Fe", (0, 0, 0)), ("Cl", (1.0, 0, 0)), ("Cl", (-1.0, 0, 0)),
        ("Cl", (0, 1.0, 0)), ("Cl", (0, -1.0, 0)),
    ])
    # fine: bigger min_md
    _write_xyz(arch, "b-BBBBBB.xyz", [
        ("Fe", (0, 0, 0)), ("Cl", (2.3, 0, 0)), ("Cl", (-2.3, 0, 0)),
        ("Cl", (0, 2.3, 0)), ("Cl", (0, -2.3, 0)),
    ])
    # mild collapse: 1.6
    _write_xyz(arch, "c-CCCCCC.xyz", [
        ("Fe", (0, 0, 0)), ("Cl", (1.6, 0, 0)), ("Cl", (-1.6, 0, 0)),
        ("Cl", (0, 1.6, 0)), ("Cl", (0, -1.6, 0)),
    ])
    rows = extract_worst_files(str(arch), "min_md", top_n=3)
    assert len(rows) == 3
    # smallest min_md must come FIRST (direction = -1 → ascending)
    assert rows[0]["file_path"].endswith("a-AAAAAA.xyz")
    assert rows[-1]["file_path"].endswith("b-BBBBBB.xyz")
    # axis_value should be the actual min M-D distance
    assert rows[0]["axis_value"] < rows[1]["axis_value"] < rows[2]["axis_value"]
    # refcode detection
    assert rows[0]["refcode"] == "AAAAAA"


def test_extract_worst_files_donor_donor(tmp_path):
    from extract_worst_files import extract_worst_files
    arch = tmp_path / "ddc"
    arch.mkdir()
    # collapsed donor-donor
    _write_xyz(arch, "x-XXXXXX.xyz", [
        ("Fe", (0, 0, 0)),
        ("Cl", (1.5, 0, 0)),
        ("Cl", (1.5, 0.7, 0)),  # 0.7 A apart!
        ("Cl", (0, 2.2, 0)),
    ])
    # fine
    _write_xyz(arch, "y-YYYYYY.xyz", [
        ("Fe", (0, 0, 0)),
        ("Cl", (2.3, 0, 0)),
        ("Cl", (-2.3, 0, 0)),
        ("Cl", (0, 2.3, 0)),
    ])
    rows = extract_worst_files(str(arch), "min_donor_donor", top_n=2)
    assert rows[0]["axis_value"] < 1.0
    assert rows[0]["file_path"].endswith("x-XXXXXX.xyz")


# ---------------------------------------------------------------------------
# Test 6: bug_census on synthetic mini-archive returns >= 2 distinct classes
# ---------------------------------------------------------------------------
def test_bug_census_on_synthetic(tmp_path):
    from bug_census import census
    arch = tmp_path / "synth"
    arch.mkdir()
    # cl_cl_collapse
    _write_xyz(arch, "01-AAAAAA.xyz", [
        ("Fe", (0, 0, 0)),
        ("Cl", (1.2, 0, 0)),
        ("Cl", (1.2, 0.9, 0)),
    ])
    # ch_collapse
    _write_xyz(arch, "02-BBBBBB.xyz", [
        ("Fe", (0, 0, 0)), ("N", (2.2, 0, 0)),
        ("C", (3.5, 0, 0)), ("H", (3.75, 0, 0)),
    ])
    # clean
    _write_xyz(arch, "03-CCCCCC.xyz", [
        ("Fe", (0, 0, 0)),
        ("Cl", (2.3, 1.2, 0)),
        ("Cl", (-2.3, 0, 0)),
        ("Cl", (0, 2.3, 0)),
        ("Cl", (0, -2.3, 0)),
    ])
    rep = census(str(arch))
    assert rep["n_files_scanned"] == 3
    # At least cl_cl_1_3_collapse and ch_collapse must fire on the synth set
    fired = set(rep["bug_classes"].keys())
    assert "cl_cl_1_3_collapse" in fired
    assert "ch_collapse" in fired
    assert rep["bug_classes"]["cl_cl_1_3_collapse"]["count"] >= 1
    # example_files must include AAAAAA for the cl_cl_collapse bucket
    ex = rep["bug_classes"]["cl_cl_1_3_collapse"]["example_files"]
    assert any("AAAAAA" in x for x in ex)


# ---------------------------------------------------------------------------
# Test 7: ccdc_overlay on a real CCDC mol2 → reasonable RMSD against itself
# ---------------------------------------------------------------------------
CCDC_DIR = "/home/qmchem_max/agent_workspace/quality_framework/CCDC"


@pytest.mark.skipif(
    not os.path.isdir(CCDC_DIR), reason="CCDC reference dir unavailable"
)
def test_ccdc_overlay_self_self_rmsd(tmp_path):
    """Convert a real CCDC mol2 to an xyz, put it in a tiny archive, and
    overlay against the CCDC dir.  Self-overlay → RMSD should be ~0."""
    from ccdc_overlay import overlay_against_ccdc
    # Pick a refcode with the canonical CCDC pattern (6+ uppercase letters).
    # Numeric refcodes (some CCDC entries) don't survive refcode_from_filename.
    import re as _re
    mol2s = sorted(f for f in os.listdir(CCDC_DIR) if f.endswith(".mol2"))
    assert mol2s
    canonical = [f for f in mol2s
                 if _re.fullmatch(r"[A-Z]{6}[A-Z0-9]{0,2}", f[:-5])]
    assert canonical, "no canonical CCDC refcodes in CCDC_DIR"
    rc = canonical[0][:-5]
    atoms = parse_mol2_atoms(os.path.join(CCDC_DIR, canonical[0]))
    assert atoms
    arch = tmp_path / "self_arch"
    arch.mkdir()
    fname = f"00-{rc}.xyz"
    _write_xyz(arch, fname, atoms)
    rep = overlay_against_ccdc(str(arch), ccdc_dir=CCDC_DIR,
                               refcodes=[rc])
    rows = [r for r in rep["per_refcode"] if r["status"] == "ok"]
    assert rows, f"no ok rows: {rep['per_refcode']}"
    assert rows[0]["heavy_rmsd"] is not None
    assert rows[0]["heavy_rmsd"] < 0.05, \
        f"self-overlay heavy_rmsd should be ~0, got {rows[0]['heavy_rmsd']}"
    assert rows[0]["md_rmsd"] is not None
    assert rows[0]["md_rmsd"] < 0.05


# ---------------------------------------------------------------------------
# Test 8: ccdc_overlay handles missing produced gracefully
# ---------------------------------------------------------------------------
@pytest.mark.skipif(
    not os.path.isdir(CCDC_DIR), reason="CCDC reference dir unavailable"
)
def test_ccdc_overlay_no_produced(tmp_path):
    import re as _re
    from ccdc_overlay import overlay_against_ccdc
    arch = tmp_path / "empty_arch"
    arch.mkdir()
    mol2s = sorted(f for f in os.listdir(CCDC_DIR) if f.endswith(".mol2"))
    canonical = [f for f in mol2s
                 if _re.fullmatch(r"[A-Z]{6}[A-Z0-9]{0,2}", f[:-5])]
    assert canonical, "no canonical CCDC refcodes"
    rc = canonical[0][:-5]
    rep = overlay_against_ccdc(str(arch), ccdc_dir=CCDC_DIR, refcodes=[rc])
    assert rep["summary"]["n_compared_ok"] == 0
    assert rep["per_refcode"][0]["status"] == "no-produced"


# ---------------------------------------------------------------------------
# Test 9: refcode parsing covers common DELFIN file-name formats
# ---------------------------------------------------------------------------
def test_refcode_from_filename_variants():
    assert refcode_from_filename("042-ARABUR.xyz") == "ARABUR"
    assert refcode_from_filename("X10-ZURHID_3d_Cr_CN4_hetero.xyz") == "ZURHID"
    # not a real refcode → None
    assert refcode_from_filename("01-Fe_CO_3_NHC_2.xyz") is None
    # Already shouting tag word: VOLLPOOL / ULTIMATE / FFFREE / MOGUL filtered
    assert refcode_from_filename("VOLLPOOL.xyz") is None


# ---------------------------------------------------------------------------
# Test 10: diff_reports returns deltas correctly
# ---------------------------------------------------------------------------
def test_diff_reports():
    from bug_census import diff_reports
    prev = {"bug_classes": {"ch_collapse": {"count": 10},
                            "cl_cl_1_3_collapse": {"count": 4}}}
    new = {"bug_classes": {"ch_collapse": {"count": 5},
                           "cl_cl_1_3_collapse": {"count": 4},
                           "h_h_clash": {"count": 2}}}
    rows = diff_reports(prev, new)
    d = {r["bug_class"]: r for r in rows}
    assert d["ch_collapse"]["delta"] == -5
    assert d["ch_collapse"]["new"] == 5
    assert d["cl_cl_1_3_collapse"]["delta"] == 0
    assert d["h_h_clash"]["delta"] == 2
    assert d["h_h_clash"]["prev"] == 0


# ---------------------------------------------------------------------------
# Test 11: per_class_ccdc_report.group_by_class handles a synthetic overlay
# ---------------------------------------------------------------------------
def test_group_by_class_synthetic():
    from per_class_ccdc_report import group_by_class
    per = [
        {"status": "ok", "refcode": "A", "ref_metal": "Fe", "ref_cn": 6,
         "ref_poly": "OC-6", "heavy_rmsd": 0.5, "md_rmsd": 0.1},
        {"status": "ok", "refcode": "B", "ref_metal": "Fe", "ref_cn": 6,
         "ref_poly": "OC-6", "heavy_rmsd": 1.5, "md_rmsd": 0.3},
        {"status": "ok", "refcode": "C", "ref_metal": "Cu", "ref_cn": 4,
         "ref_poly": "SP-4", "heavy_rmsd": 0.2, "md_rmsd": 0.05},
        {"status": "no-produced", "refcode": "MISS"},
    ]
    grp = group_by_class(per)
    assert "Fe|CN6|OC-6" in grp
    fe = grp["Fe|CN6|OC-6"]
    assert fe["n_refs"] == 2
    assert abs(fe["heavy_rmsd_mean"] - 1.0) < 1e-6
    # worst is the higher-RMSD one
    assert fe["worst_refcode"] == "B"
    assert abs(fe["worst_heavy_rmsd"] - 1.5) < 1e-6
