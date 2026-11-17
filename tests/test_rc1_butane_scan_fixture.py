"""Tests for the RC1 butane relaxed-scan fixture (package RC1, wave 23).

The fixture ``tests/fixtures/orca611/rc1_butane_scan_r2scan3c`` is a REAL
ORCA 6.1.0 run: a relaxed scan of the central C-C-C-C dihedral of n-butane
at r2SCAN-3c over the full 0-360 deg circle in 15 deg steps (Scan
D 0 1 2 3 = 180, 525, 24 — the 24 distinct angles; the original brief's
literal "0, 360, 25" failed three times on ORCA mechanics: NoFreq rejected,
GSTEP could not start 180 deg away, and the ±180 deg wrap point did not
converge — see PROVENANCE.md), starting from the GOAT global
minimum in ``butane_goat_gfn2/butane_goat_gfn2.globalminimum.xyz``.

Every check reuses DELFIN's own ORCA readers; nothing is parsed a second
way:

* ``delfin.tools.orca_out_index.orca_out_section`` returns the lines of the
  n-th ``FINAL SINGLE POINT ENERGY`` section; ORCA writes the energy ON the
  heading line itself (see the GOAT fixture, line 633), so the value is
  read from the heading line.
* ``delfin.api.parse_orca_output`` gives the whole-file snapshot.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from delfin.api import parse_orca_output
from delfin.tools.orca_out_index import orca_out_index, orca_out_section

# Real ORCA fixtures live in the private store outside the repo
# (LICENSE_GUARD.md rule family 2); skip without DELFIN_PRIVATE_FIXTURES.
from test_lc1_fixtures_helper import lc1_fixture_dir

FIXTURE = lc1_fixture_dir("tests/fixtures/orca611/rc1_butane_scan_r2scan3c")

KCAL = 627.5094740631  # Eh -> kcal/mol (CODATA, the factor ORCA tables use)


def _final_single_point_energy(out_path: Path, n: int) -> float:
    """The n-th scan point's electronic energy, read through the index."""
    lines = orca_out_section(out_path, "FINAL SINGLE POINT ENERGY", n=n)
    assert lines, "empty FINAL SINGLE POINT ENERGY section"
    for line in lines:
        if line.lstrip().startswith("FINAL SINGLE POINT ENERGY"):
            return float(line.split()[-1])
    pytest.fail("no heading line in FINAL SINGLE POINT ENERGY section")


@pytest.fixture(scope="module")
def out_path() -> Path:
    p = FIXTURE / "rc1_butane_scan_r2scan3c.out"
    if not p.is_file():
        pytest.skip("RC1 fixture not committed yet")
    return p


def test_fixture_folder_is_committed():
    """Red control: fails while the fixture does not exist yet.

    The other tests skip on a missing fixture (they arm when it is
    committed); this one makes the missing state a visible failure so the
    red-before/green-after control is real.
    """
    assert FIXTURE.is_dir(), "RC1 fixture available (repo or private store)"
    assert (FIXTURE / "rc1_butane_scan_r2scan3c.inp").is_file()
    assert (FIXTURE / "PROVENANCE.md").is_file()


def test_inp_holds_the_scan_request(out_path):
    inp = (FIXTURE / "rc1_butane_scan_r2scan3c.inp").read_text()
    assert "r2SCAN-3c" in inp
    # v5 scan line: 24 distinct angles, 15 deg apart, covering the full
    # circle (180..525 spans 345 deg; the last step wraps mod 360).
    assert "D 0 1 2 3 = 180, 525, 24" in inp
    assert "* xyz 0 1" in inp  # charge 0, multiplicity 1
    assert inp.count("\nC ") >= 4  # butane carbons present


def test_exactly_24_scan_points(out_path):
    titles = [s.title for s in orca_out_index(out_path)]
    n = titles.count("FINAL SINGLE POINT ENERGY")
    assert n == 24, f"expected 24 scan points, index found {n}"


def test_all_24_energies_finite(out_path):
    energies = [_final_single_point_energy(out_path, n) for n in range(1, 25)]
    assert all(math.isfinite(e) for e in energies)


def test_profile_covers_the_full_circle(out_path):
    """24 points x 15 deg span 345 deg; the wrap back to 180 closes the
    circle. The grid must contain every distinct angle of the 0-360
    profile (checked against the inp's own scan line)."""
    inp = (FIXTURE / "rc1_butane_scan_r2scan3c.inp").read_text()
    import re
    m = re.search(r"D\s+0\s+1\s+2\s+3\s*=\s*([0-9.]+),\s*([0-9.]+),\s*(\d+)", inp)
    assert m, "no Scan D line in the fixture inp"
    start, stop, npts = float(m.group(1)), float(m.group(2)), int(m.group(3))
    span = stop - start
    step = span / (npts - 1)
    assert span > 0
    assert (span + step) % 360 == 0, "grid must close the circle"
    angles = [start + i * step for i in range(npts)]
    mod360 = sorted(a % 360 for a in angles)
    assert len(set(mod360)) == npts, "angles must be distinct mod 360"


def test_syn_is_the_highest_point(out_path):
    energies = [_final_single_point_energy(out_path, n) for n in range(1, 25)]
    # syn = the 0/360 family = step 13 (index 12) in the v5 grid
    assert max(energies) == energies[12], (
        "the syn (360 deg) point must be the highest of the scan"
    )


def test_gauche_anti_gap_in_physical_range(out_path):
    energies = [_final_single_point_energy(out_path, n) for n in range(1, 25)]
    anti = min(energies)  # step 1, 180 deg
    gauche = min(energies[7], energies[8], energies[15], energies[16])
    gap_kcal = (gauche - anti) * KCAL
    # r2SCAN-3c literature class: 0.6-0.9 kcal/mol; wide enough band to
    # not pin the fixture to a tuned value, narrow enough to catch a
    # wrong minimum or a mis-scaled surface. A 15-deg grid point
    # under-reports the gap (the true minimum lies between grid points).
    assert 0.3 <= gap_kcal <= 1.5, (
        f"gauche-anti gap {gap_kcal:.3f} kcal/mol outside the expected class"
    )


def test_every_point_optimization_converged(out_path):
    text = out_path.read_text(encoding="utf-8", errors="replace")
    assert "DID NOT CONVERGE" not in text
    assert "FAILED TO CONVERGE" not in text
    converged = text.count("THE OPTIMIZATION HAS CONVERGED") + text.count(
        "OPTIMIZATION RUN DONE"
    )
    assert converged >= 24, (
        f"expected >= 24 converged optimisation blocks, found {converged}"
    )


def test_parse_orca_output_reports_the_run(out_path):
    result = parse_orca_output(str(out_path))
    assert result.status == "ok"
    assert result.final_single_point is not None
    assert result.scf_converged is True


def test_property_json_committed():
    """The QS F2 ruling (nacht-s5): .property.json alone satisfies the
    fixture contract. The v5 .inp requested %output jsonpropfile True
    only -- no bare jsonfile -- so ORCA 6.1.0 wrote exactly this one
    file; the bare .json pin was rewritten along the ruling (the old
    skip-if-absent form hid a file the run correctly never wrote)."""
    p = FIXTURE / "rc1_butane_scan_r2scan3c.property.json"
    if not p.is_file():
        pytest.skip("property.json not committed yet")
    import json
    data = json.loads(p.read_text())
    assert data, "property.json is empty"
    # Real ORCA property output, not a stub: ORCA's property JSON carries
    # the run's metadata (Calculation_Info / NAtoms / Calculation_Status).
    assert (
        "Calculation_Info" in data
        or "NAtoms" in data
        or "Calculation_Status" in data
    ), "property.json does not look like real ORCA property output"
    # The .inp must request jsonpropfile True -- that is why this file
    # exists and why no bare .json does (QS F2 ruling, 2026-11-04).
    inp = (FIXTURE / "rc1_butane_scan_r2scan3c.inp").read_text()
    assert "jsonpropfile True" in inp
    assert "jsonfile True" not in inp
    # The bare .json was never requested and must stay absent: its
    # presence would mean the .inp asked for output modes the
    # PROVENANCE says it did not.
    assert not (FIXTURE / "rc1_butane_scan_r2scan3c.json").is_file()
