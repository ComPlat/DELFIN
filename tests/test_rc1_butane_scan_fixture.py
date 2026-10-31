"""Tests for the RC1 butane relaxed-scan fixture (package RC1, wave 23).

The fixture ``tests/fixtures/orca611/rc1_butane_scan_r2scan3c`` is a REAL
ORCA 6.1.1 run: a relaxed scan of the central C-C-C-C dihedral of n-butane
at r2SCAN-3c, 0 to 360 deg in 25 steps, starting from the GOAT global
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

FIXTURE = Path("tests/fixtures/orca611/rc1_butane_scan_r2scan3c")

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
    assert FIXTURE.is_dir(), "RC1 fixture folder not committed yet"
    assert (FIXTURE / "rc1_butane_scan_r2scan3c.inp").is_file()
    assert (FIXTURE / "PROVENANCE.md").is_file()


def test_inp_holds_the_scan_request(out_path):
    inp = (FIXTURE / "rc1_butane_scan_r2scan3c.inp").read_text()
    assert "r2SCAN-3c" in inp
    assert "D 0 1 2 3 = 0, 360, 25" in inp
    assert "* xyz 0 1" in inp  # charge 0, multiplicity 1
    assert inp.count("\nC ") >= 4  # butane carbons present


def test_exactly_25_scan_points(out_path):
    titles = [s.title for s in orca_out_index(out_path)]
    n = titles.count("FINAL SINGLE POINT ENERGY")
    assert n == 25, f"expected 25 scan points, index found {n}"


def test_all_25_energies_finite(out_path):
    energies = [_final_single_point_energy(out_path, n) for n in range(1, 26)]
    assert all(math.isfinite(e) for e in energies)


def test_end_point_equals_start_point(out_path):
    start = _final_single_point_energy(out_path, 1)
    end = _final_single_point_energy(out_path, 25)
    assert abs(end - start) < 1e-6, "360 deg point must close on 0 deg"


def test_syn_is_the_highest_point(out_path):
    energies = [_final_single_point_energy(out_path, n) for n in range(1, 26)]
    assert max(energies) == energies[0], (
        "the syn (0 deg) point must be the highest of the scan"
    )


def test_gauche_anti_gap_in_physical_range(out_path):
    energies = [_final_single_point_energy(out_path, n) for n in range(1, 26)]
    anti = min(energies)
    gauche = min(energies[3:5] + energies[19:21])  # 60 and 120 deg points
    gap_kcal = (gauche - anti) * KCAL
    # r2SCAN-3c literature class: 0.6-0.9 kcal/mol; wide enough band to
    # not pin the fixture to a tuned value, narrow enough to catch a
    # wrong minimum or a mis-scaled surface.
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
    assert converged >= 25, (
        f"expected >= 25 converged optimisation blocks, found {converged}"
    )


def test_parse_orca_output_reports_the_run(out_path):
    result = parse_orca_output(str(out_path))
    assert result.status == "ok"
    assert result.final_single_point is not None
    assert result.scf_converged is True


def test_property_json_committed():
    for name in ("rc1_butane_scan_r2scan3c.json",
                 "rc1_butane_scan_r2scan3c.property.json"):
        p = FIXTURE / name
        if not p.is_file():
            pytest.skip(f"{name} not committed yet")
        import json
        data = json.loads(p.read_text())
        assert data, f"{name} is empty"
