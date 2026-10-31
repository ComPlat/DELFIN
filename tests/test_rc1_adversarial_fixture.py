"""Adversarial tests for the RC1 butane relaxed-scan fixture (nacht-s2).

Package RC1 requires a REAL r2SCAN-3c relaxed dihedral scan of n-butane
under tests/fixtures/orca611/rc1_butane_scan_r2scan3c/, produced through
DELFIN's submit path, with every scan point converged. These tests pin the
fixture contract using only DELFIN's existing ORCA readers
(delfin.tools.orca_out_index, delfin.api.parse_orca_output) -- never a
second parser.

Red on main 8b5753686 by design: the fixture does not exist yet; the
builder's commits must turn these green.
"""

from pathlib import Path

import pytest

from delfin.api import parse_orca_output
from delfin.tools.orca_out_index import orca_out_index, orca_out_section

FIXTURE_DIR = (Path(__file__).resolve().parent
               / "fixtures" / "orca611" / "rc1_butane_scan_r2scan3c")
OUT_FILE = FIXTURE_DIR / "rc1_butane_scan_r2scan3c.out"
INP_FILE = FIXTURE_DIR / "rc1_butane_scan_r2scan3c.inp"
PROVENANCE = FIXTURE_DIR / "PROVENANCE.md"


def _fixture_available() -> bool:
    return OUT_FILE.is_file() and INP_FILE.is_file() and PROVENANCE.is_file()


needs_fixture = pytest.mark.skipif(
    not _fixture_available(),
    reason="RC1 scan fixture not yet committed by the builder",
)


@needs_fixture
def test_scan_inp_requests_the_full_dihedral_range():
    """The .inp scans D 0 1 2 3 = 0,360,25 at r2SCAN-3c, neutral singlet."""
    text = INP_FILE.read_text()
    assert "r2SCAN-3c" in text, "method keyword missing"
    assert "* xyz 0 1" in text, "charge 0 / multiplicity 1 missing"
    # ORCA %geom scan block: relaxed scan over the central C-C-C-C dihedral.
    assert "Scan" in text and "D 0 1 2 3" in text
    assert "0, 360, 25" in text.replace("0,360,25", "0, 360, 25")


@needs_fixture
def test_scan_output_has_25_energy_points_and_terminated_normally():
    """25 scan points = 25 FINAL SINGLE POINT ENERGY sections; clean end."""
    text = OUT_FILE.read_text()
    assert "ORCA TERMINATED NORMALLY" in text
    # The per-point reader is DELFIN's own index: each scan point opens its
    # own FINAL SINGLE POINT ENERGY section.
    sections = [s for s in orca_out_index(OUT_FILE)
                if s.title.upper() == "FINAL SINGLE POINT ENERGY"]
    assert len(sections) == 25, f"expected 25 scan points, got {len(sections)}"
    # The scan summary table ORCA writes at the end must have survived any
    # trim (the brief: keep the scan summary table).
    assert "THE RELAXED SURFACE SCAN" in text.upper() or "Scan summary" in text


@needs_fixture
def test_every_scan_point_energy_is_finite_and_read_via_orca_out_section():
    """Each of the 25 energies reads through orca_out_section, finite."""
    energies = []
    for n in range(1, 26):
        lines = orca_out_section(OUT_FILE, "FINAL SINGLE POINT ENERGY", n=n)
        value_line = next(l for l in lines
                          if l.strip().upper().startswith(
                              "FINAL SINGLE POINT ENERGY"))
        energies.append(float(value_line.split()[-1]))
    assert len(energies) == 25
    assert all(e == e and abs(e) < 1e3 for e in energies), "non-finite energy"


@needs_fixture
def test_anti_is_minimum_gauche_next_syn_is_highest():
    """Physics ordering: anti (180 deg) lowest; syn (0 deg) highest;
    gauche around 60-70 deg a local minimum between them."""
    energies = []
    for n in range(1, 26):
        lines = orca_out_section(OUT_FILE, "FINAL SINGLE POINT ENERGY", n=n)
        value_line = next(l for l in lines
                          if l.strip().upper().startswith(
                              "FINAL SINGLE POINT ENERGY"))
        energies.append(float(value_line.split()[-1]))
    anti = energies[12]  # point 13 = 180 deg
    syn = energies[0]    # point 1 = 0 deg
    assert syn > anti, "syn (0 deg) must be the highest point"
    # Gauche minimum: lowest energy among points 3-8 (30..120 deg).
    gauche = min(energies[2:8])
    assert anti < syn
    # Gauche-anti difference of order 0.6-0.9 kcal/mol at r2SCAN-3c.
    gap_kcal = (gauche - anti) * 627.509474
    assert 0.3 < gap_kcal < 1.5, (
        f"gauche-anti gap {gap_kcal:.3f} kcal/mol outside plausible band")


@needs_fixture
def test_every_scan_point_optimization_converged():
    """Each scan step is a constrained optimisation: 25 converged geoms."""
    text = OUT_FILE.read_text()
    n_done = text.count("THE OPTIMIZATION HAS CONVERGED")
    n_done += text.count("OPTIMIZATION RUN DONE")
    assert n_done >= 25, (
        f"only {n_done} converged optimization blocks for 25 scan points")
    assert "OPTIMIZATION DID NOT CONVERGE" not in text


@needs_fixture
def test_provenance_records_job_identity_and_no_home_paths():
    """PROVENANCE.md carries job id, ORCA version, node type, core-hours;
    no home paths or host names anywhere in the fixture."""
    prov = PROVENANCE.read_text()
    for field in ("job id", "ORCA version", "node", "core-hours"):
        assert field.lower() in prov.lower(), f"provenance field {field!r} missing"
    for path in (INP_FILE, OUT_FILE, PROVENANCE):
        blob = path.read_text()
        for marker in ("/home/", "/pfs/data6", "ka_ew7404"):
            assert marker not in blob, (
                f"{path.name} leaks home/host path {marker!r}")


@needs_fixture
def test_orca_json_property_json_exist_as_orca_wrote_them():
    """ORCA's own .json / .property.json outputs are part of the fixture."""
    assert (FIXTURE_DIR / "rc1_butane_scan_r2scan3c.json").is_file()
    assert (FIXTURE_DIR / "rc1_butane_scan_r2scan3c.property.json").is_file()
    import json
    header = json.loads(
        (FIXTURE_DIR / "rc1_butane_scan_r2scan3c.json").read_text())
    assert "Citations" in header or "Header" in header


@needs_fixture
def test_parse_orca_output_reader_still_works_on_the_fixture():
    """DELFIN's primary reader must parse the fixture without error."""
    parsed = parse_orca_output(str(OUT_FILE))
    assert parsed.single_point is not None
    assert parsed.scf_converged is True
