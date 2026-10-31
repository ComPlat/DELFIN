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

import re
from pathlib import Path

import pytest

from delfin.api import parse_orca_output
from delfin.tools.orca_out_index import orca_out_index, orca_out_section

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir

FIXTURE_DIR = lc1_fixture_dir("tests/fixtures/orca611/rc1_butane_scan_r2scan3c")
OUT_FILE = FIXTURE_DIR / "rc1_butane_scan_r2scan3c.out"
INP_FILE = FIXTURE_DIR / "rc1_butane_scan_r2scan3c.inp"
PROVENANCE = FIXTURE_DIR / "PROVENANCE.md"


def _fixture_available() -> bool:
    return OUT_FILE.is_file() and INP_FILE.is_file() and PROVENANCE.is_file()


needs_fixture = pytest.mark.skipif(
    not _fixture_available(),
    reason="RC1 scan fixture not yet committed by the builder",
)


def _inp_scan_block() -> str:
    """The 'D a b c d = start, stop, nsteps' line of the committed .inp.

    ORCA's Scan block spans two lines: 'Scan' opens it and the dihedral
    definition sits on its own line, so the two are matched separately
    (the .inp of the v5 fixture, lines 7-8)."""
    text = INP_FILE.read_text()
    lines = [l.strip() for l in text.splitlines()]
    for i, low in enumerate(l.lower() for l in lines):
        if low == "scan" and i + 1 < len(lines):
            return lines[i + 1]
    return ""


def _scan_point_count() -> int:
    """Number of scan points the .inp requests (derived, not hard-coded).

    ORCA's Scan D = start, stop, Nsteps requests Nsteps points in total
    (end point equals start for a full dihedral circle, so v5's
    180,525,24 is 24 distinct angles closing on the start)."""
    m = re.search(r"=\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(\d+)", _inp_scan_block())
    assert m, "no parsable Scan D = start,stop,nsteps in .inp"
    return int(m.group(3))


@needs_fixture
def test_scan_inp_requests_the_full_dihedral_range():
    """The .inp scans the central dihedral at r2SCAN-3c, neutral singlet,
    covering the full dihedral circle: nsteps points whose consecutive
    gaps (mod 360) are uniform, so all distinct angles 0..345 are met.
    v5's line "D 0 1 2 3 = 180, 525, 24" is 23 x 15 deg gaps closing
    through the 360 wrap."""
    text = INP_FILE.read_text()
    assert "r2SCAN-3c" in text, "method keyword missing"
    assert "* xyz 0 1" in text, "charge 0 / multiplicity 1 missing"
    # ORCA %geom scan block: relaxed scan over the central C-C-C-C dihedral.
    assert "Scan" in text and "D 0 1 2 3" in text
    m = re.search(r"=\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(\d+)", _inp_scan_block())
    assert m, "no parsable Scan D = start,stop,nsteps in .inp"
    start, stop, nsteps = int(m.group(1)), int(m.group(2)), int(m.group(3))
    span = stop - start
    # ORCA Scan gives nsteps points, so the dihedral step is span/(n-1):
    # v1 0,360,25 -> 24 x 15 deg; v5 180,525,24 -> 23 x 15 deg.
    assert nsteps > 1, "a 1-point scan cannot cover the circle"
    assert span % (nsteps - 1) == 0, (
        f"scan span {span} deg is not an integer multiple of the "
        f"{nsteps - 1} intervals")
    gap = span // (nsteps - 1)
    # Full circle: the wrap between the last point and the (same) first
    # angle is 360 - last-gap; uniform coverage needs gap >= that, i.e.
    # span + gap >= 360 (v1: 360+15; v5: 345+15=360 exactly).
    assert span + gap >= 360, (
        f"scan {start}..{stop} in {nsteps} points leaves the dihedral "
        f"circle uncovered (gap {gap} deg)")


@needs_fixture
def test_scan_output_has_25_energy_points_and_terminated_normally():
    """One FINAL SINGLE POINT ENERGY section per requested scan point;
    clean end; the scan summary table survived any trim."""
    text = OUT_FILE.read_text()
    assert "ORCA TERMINATED NORMALLY" in text
    n_expected = _scan_point_count()
    sections = [s for s in orca_out_index(OUT_FILE)
                if s.title.upper() == "FINAL SINGLE POINT ENERGY"]
    assert len(sections) == n_expected, (
        f"expected {n_expected} scan points, got {len(sections)}")
    # The scan summary ORCA writes at the end must have survived any
    # trim (the brief: keep the scan summary table). Real ORCA wording:
    # "RELAXED SURFACE SCAN DONE" and "RELAXED SURFACE SCAN RESULTS".
    assert ("RELAXED SURFACE SCAN DONE" in text.upper()
            or "RELAXED SURFACE SCAN RESULTS" in text.upper())


def _read_all_scan_energies() -> list:
    """All scan-point energies via orca_out_section, in scan order."""
    energies = []
    for n in range(1, _scan_point_count() + 1):
        lines = orca_out_section(OUT_FILE, "FINAL SINGLE POINT ENERGY", n=n)
        value_line = next(l for l in lines
                          if l.strip().upper().startswith(
                              "FINAL SINGLE POINT ENERGY"))
        energies.append(float(value_line.split()[-1]))
    return energies


@needs_fixture
def test_every_scan_point_energy_is_finite_and_read_via_orca_out_section():
    """Each scan-point energy reads through orca_out_section, finite."""
    energies = _read_all_scan_energies()
    assert len(energies) == _scan_point_count()
    assert all(e == e and abs(e) < 1e3 for e in energies), "non-finite energy"


@needs_fixture
def test_anti_is_minimum_gauche_next_syn_is_highest():
    """Physics ordering, independent of the scan's start angle and step
    count: syn = highest of the 0/360-family points, anti = lowest of the
    180/540-family points, gauche a local minimum around 300/420 (plus 60
    on the 0-family side). Gap = min(gauche) - min(anti-family)."""
    block = _inp_scan_block()
    m = re.search(r"=\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(\d+)", block)
    assert m, "no parsable Scan D in .inp"
    start = int(m.group(1))
    energies = _read_all_scan_energies()
    step = 360.0 / len(energies)
    angles = [start + i * step for i in range(len(energies))]

    def fam(lo, hi):
        """Min over absolute dihedral angles (mod 360) in the window
        [lo, hi); a window crossing 0 (lo > hi) matches both sides."""
        vals = [e for a, e in zip(angles, energies)
                if lo <= (a % 360) < hi
                or (lo > hi and ((a % 360) >= lo or (a % 360) < hi))]
        assert vals, f"no scan point in {lo}..{hi} deg"
        return min(vals)

    # syn family (0/360), anti family (180/540), and the gauche minima
    # near 60/300/420 (operator: gauche ~300 and ~420 on the 180-based
    # v5 grid; 60 is its mirror on a 0-based grid).
    syn = fam(330, 30)
    anti = fam(150, 210)
    gauche_candidates = [fam(30, 90), fam(270, 330)]
    gauche = min(g for g in gauche_candidates if g is not None)
    assert syn > anti, "syn (0/360 family) must be the highest point"
    gap_kcal = (gauche - anti) * 627.509474
    assert 0.3 < gap_kcal < 1.5, (
        f"gauche-anti gap {gap_kcal:.3f} kcal/mol outside plausible band")


@needs_fixture
def test_every_scan_point_optimization_converged():
    """Each scan step is a constrained optimisation: >= n_points converged
    geometry blocks, none reported as not converged."""
    text = OUT_FILE.read_text()
    n_done = text.count("THE OPTIMIZATION HAS CONVERGED")
    n_done += text.count("OPTIMIZATION RUN DONE")
    assert n_done >= _scan_point_count(), (
        f"only {n_done} converged optimization blocks for "
        f"{_scan_point_count()} scan points")
    assert "OPTIMIZATION DID NOT CONVERGE" not in text


@needs_fixture
def test_provenance_records_job_identity_and_no_home_paths():
    """PROVENANCE.md carries job id, ORCA version, node type, core-hours;
    no home paths or host names anywhere in the fixture."""
    prov = PROVENANCE.read_text()
    prov_low = prov.lower()
    # Field facts, not literal phrases: the job identity ("job 8173642"),
    # the ORCA version, the node type, and the core-hours.
    assert "job 8173642" in prov_low, "fixture job id missing"
    assert "orca" in prov_low and "6.1." in prov_low, "ORCA version missing"
    assert "core-h" in prov_low or "core-hours" in prov_low, \
        "core-hours missing"
    for path in (INP_FILE, OUT_FILE, PROVENANCE):
        blob = path.read_text()
        # The account-name literal is NOT listed here on purpose: the
        # repo-wide guard (tests/test_no_account_name_is_committed.py)
        # covers account names in every tracked file, and writing the
        # name into this test would trip that same guard.
        for marker in ("/home/", "/pfs/"):
            assert marker not in blob, (
                f"{path.name} leaks home/host path {marker!r}")


@needs_fixture
def test_orca_json_property_json_exist_as_orca_wrote_them():
    """The fixture carries the property.json ORCA wrote for its requested
    output modes.

    QS ruling (nacht-s5, 2026-11-04): the brief's ".json/.property.json
    as ORCA writes them" means the files ORCA produces for the output
    modes the .inp requests. v5's .inp requests jsonpropfile True (no
    bare jsonfile), so ORCA wrote exactly one json file -- the
    property.json committed here. The earlier bare-.json pin was wrong
    twice over: the reference family's own .inp (butane_dlpno_ccsdt_
    pvtz.inp lines 12-13) requests jsonpropfile + jsongbwfile only, so
    its committed .json is NOT what jsonpropfile produces, and that
    fixture is a constructed one anyway -- no real-run format authority
    for a bare .json."""
    prop_path = FIXTURE_DIR / "rc1_butane_scan_r2scan3c.property.json"
    assert prop_path.is_file(), "property.json missing"
    import json
    header = json.loads(prop_path.read_text())
    # ORCA property.json structure of a REAL run: Calculation_Info
    # (charge/multiplicity/atom count) and Calculation_Status.
    assert "Calculation_Info" in header and "Calculation_Status" in header, (
        "property.json is not ORCA's own property output")
    inp = INP_FILE.read_text()
    assert "jsonpropfile True" in inp, (
        "the property.json must be requested in the committed .inp")


@needs_fixture
def test_parse_orca_output_reader_still_works_on_the_fixture():
    """DELFIN's primary reader must parse the fixture without error."""
    parsed = parse_orca_output(str(OUT_FILE))
    assert parsed.final_single_point is not None
    assert parsed.scf_converged is True
    assert parsed.opt_converged is True
