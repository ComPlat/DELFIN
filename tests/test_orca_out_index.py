"""orca_out_index / orca_out_section: an index of ORCA output sections.

Reads the section banners ORCA prints in real ORCA output (the orca611
runs, now in the PRIVATE store outside the repo -- LICENSE_GUARD.md rule
family 2) and verifies that helper lets a caller find a section by
title, get its lines, select the n-th occurrence, and cap the returned
size. Runs under the CI interpreter with no dependency on OPI. Every
test SKIPS without DELFIN_PRIVATE_FIXTURES (never fails in CI).
"""
import os
import shutil
from pathlib import Path

import pytest

from delfin.tools.orca_out_index import orca_out_index, orca_out_section

# committed fixture path when present (pre-removal compat), else private
_LEGACY = Path("tests/fixtures/orca611")
_PRIVATE_ENV = "DELFIN_PRIVATE_FIXTURES"


def _fixture_dir(name: str) -> Path | None:
    rel = f"tests/fixtures/orca611/{name}"
    if (_LEGACY / name).is_dir():
        return _LEGACY / name
    root = os.environ.get(_PRIVATE_ENV, "").strip()
    if not root:
        return None
    p = Path(root) / rel
    return p if p.is_dir() else None


def out(name):
    """Path to the fixture's .out; SKIPs when only the private store has it."""
    d = _fixture_dir(name)
    if d is None:
        pytest.skip(
            f"real ORCA fixtures live in the private store; set "
            f"{_PRIVATE_ENV} to run (LICENSE_GUARD.md rule family 2)")
    return str(d / f"{name}.out")


def titles(name):
    return [s.title for s in orca_out_index(out(name))]


def starts(name):
    return {s.title: s.start_line for s in orca_out_index(out(name))}


def test_o2_sections_are_indexed_in_order():
    assert titles("o2_triplet") == [
        "SINGLE POINT CALCULATION",
        "SCF ITERATIONS",
        "TOTAL SCF ENERGY",
        "ORBITAL ENERGIES",
        "MULLIKEN POPULATION ANALYSIS",
        "MULLIKEN ATOMIC CHARGES",
        "LOEWDIN POPULATION ANALYSIS",
        "MAYER POPULATION ANALYSIS",
        "FINAL SINGLE POINT ENERGY",
        "ORCA TERMINATED NORMALLY",
    ]


def test_o2_section_start_lines_match_fixture_banners():
    assert starts("o2_triplet")["SINGLE POINT CALCULATION"] == 216
    assert starts("o2_triplet")["SCF ITERATIONS"] == 584
    assert starts("o2_triplet")["FINAL SINGLE POINT ENERGY"] == 882
    assert starts("o2_triplet")["ORCA TERMINATED NORMALLY"] == 1087


def test_o2_final_energy_lines_returned():
    lines = orca_out_section(out("o2_triplet"), "FINAL SINGLE POINT ENERGY")
    txt = "\n".join(lines)
    assert "-150.051820514846" in txt


def test_failed_run_abort_block_is_indexed():
    st = starts("scf_not_converged")
    assert "SCF ITERATIONS" in st
    assert st["ERROR"] == 599
    lines = orca_out_section(out("scf_not_converged"), "ERROR")
    txt = "\n".join(lines)
    assert "SCF NOT CONVERGED AFTER" in txt


def test_optfreq_optimization_cycles_and_frequencies():
    st = starts("h2o_optfreq")
    for n in (1, 2, 3, 4):
        assert f"GEOMETRY OPTIMIZATION CYCLE {n}" in st
    assert st["GEOMETRY OPTIMIZATION CYCLE 1"] == 267
    assert st["VIBRATIONAL FREQUENCIES"] == 3280
    assert st["IR SPECTRUM"] == 3327


def test_tddft_excited_states():
    st = starts("h2o_tddft")
    assert st["TD-DFT/TDA EXCITED STATES"] == 995
    assert st["CIS/TD-DFT TOTAL ENERGY"] == 1064
    txt = "\n".join(orca_out_section(out("h2o_tddft"), "TD-DFT/TDA EXCITED STATES"))
    assert "TD-DFT/TDA EXCITED STATES (SINGLETS)" in txt


def test_nth_occurrence_selection():
    # h2o_optfreq has five FINAL SINGLE POINT ENERGY sections (four cycles
    # plus the final post-optimization energy).
    sec2 = orca_out_index(out("h2o_optfreq"))
    fsec = [s for s in sec2 if s.title == "FINAL SINGLE POINT ENERGY"]
    assert len(fsec) == 5
    assert fsec[0].start_line == 905
    assert fsec[4].start_line == 3086
    # orca_out_section picks the n-th (1-based) occurrence.
    lines3 = orca_out_section(out("h2o_optfreq"), "FINAL SINGLE POINT ENERGY", n=3)
    assert "-76.321266086023" in "\n".join(lines3)


def test_max_lines_cap_and_missing_title():
    lines = orca_out_section(out("o2_triplet"), "SCF ITERATIONS", max_lines=3)
    assert len(lines) == 3
    with pytest.raises(KeyError):
        orca_out_section(out("o2_triplet"), "NO SUCH SECTION")


def test_sections_cover_the_file_gap_free():
    # start of first section and end of last bound the whole file;
    # every section end is before the next start, with no gaps.
    for name in ("o2_triplet", "h2o_tddft", "scf_not_converged", "h2o_optfreq"):
        secs = orca_out_index(out(name))
        assert secs[0].start_line > 0
        for a, b in zip(secs, secs[1:]):
            assert a.end_line == b.start_line - 1
