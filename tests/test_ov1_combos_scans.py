"""OV1 phase 2: illegal combinations and scan-block checks.

Wave 23 incidents this covers: NumFreq + Freq together, RI without a
matching auxiliary (/C or /J) basis, an unknown CPCM solvent, a
duplicate scan constraint, a scan start far from the geometry's value,
and a scan end point duplicating the start across the +-180 wrap.

All tests drive the real public path delfin.api.validate_orca_input
(the hook from phase 1 is on the branch).
"""
import pytest

from delfin.api import validate_orca_input


def _errors(inp):
    return [i for i in validate_orca_input(inp) if i.severity == "error"]


def _codes(inp):
    return [i.code for i in validate_orca_input(inp)]


def _inp(keywords, extra=""):
    return (
        "! " + keywords + "\n"
        "%maxcore 2000\n"
        "%pal nprocs 2 end\n"
        "* xyz 0 1\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 0.74\n"
        "*\n" + extra
    )


# --- illegal combinations -------------------------------------------------

def test_numfreq_and_freq_conflict_is_an_error():
    issues = _errors(_inp("PBE0 def2-SVP Opt NumFreq Freq"))
    assert any(i.code == "conflicting_keywords" for i in issues), _codes(_inp("PBE0 def2-SVP Opt NumFreq Freq"))


def test_numfreq_alone_is_clean():
    assert "conflicting_keywords" not in _codes(_inp("PBE0 def2-SVP NumFreq"))


def test_ri_without_auxiliary_basis_is_an_error():
    issues = _errors(_inp("PBE0 def2-SVP RI Opt"))
    assert any(i.code == "ri_without_aux_basis" for i in issues), _codes(_inp("PBE0 def2-SVP RI Opt"))


def test_ri_with_matching_aux_basis_is_clean():
    assert "ri_without_aux_basis" not in _codes(_inp("PBE0 def2-SVP def2-SVP/C RI Opt"))


def test_unknown_cpcm_solvent_is_an_error():
    issues = _errors(_inp("PBE0 def2-SVP CPCM(NoSuchSolventX)"))
    assert any(i.code == "unknown_solvent" for i in issues), _codes(_inp("PBE0 def2-SVP CPCM(NoSuchSolventX)"))


def test_cpcm_water_case_insensitive_is_clean():
    for form in ("CPCM(water)", "CPCM(Water)"):
        assert "unknown_solvent" not in _codes(_inp(f"PBE0 def2-SVP {form}")), form


def test_unknown_smd_solvent_is_an_error():
    issues = _errors(_inp("PBE0 def2-SVP SMD(NotASolvent)"))
    assert any(i.code == "unknown_solvent" for i in issues)


# --- scan checks -----------------------------------------------------------

def _scan_inp(keywords, scan_body, geom=None):
    geom = geom or ("H 0.0 0.0 0.0\nH 0.0 0.0 0.74\nH 0.0 0.0 1.90\n")
    return (
        "! " + keywords + "\n"
        "* xyz 0 1\n" + geom + "*\n"
        "%scan\n" + scan_body + "\nend\n"
    )


def test_duplicate_scan_constraint_is_an_error():
    inp = _scan_inp("PBE0 def2-SVP Opt", "{B 0 1 1.0 2.0 10}S\n{B 0 1 1.0 2.0 10}S\n")
    issues = _errors(inp)
    assert any(i.code == "duplicate_scan_constraint" for i in issues), _codes(inp)


def test_distinct_scan_constraints_are_clean():
    inp = _scan_inp("PBE0 def2-SVP Opt",
                    "{B 0 1 1.0 2.0 10}S\n{A 1 2 100.0 110.0 5}S\n")
    issues = _errors(inp)
    assert not [i for i in issues if i.code == "duplicate_scan_constraint"], _codes(inp)


def test_scan_start_far_from_geometry_is_reported():
    """The H2-H3 bond is 1.9 A; a scan starting at 3.8 A is far off."""
    inp = _scan_inp("PBE0 def2-SVP Opt", "{B 2 0 3.8 2.0 10}S\n")
    issues = [i for i in validate_orca_input(inp) if i.code == "scan_start_far_from_geometry"]
    assert issues, _codes(inp)


def test_scan_start_near_geometry_is_clean():
    inp = _scan_inp("PBE0 def2-SVP Opt", "{B 2 0 1.9 2.0 10}S\n")
    assert not [i for i in validate_orca_input(inp)
                if i.code == "scan_start_far_from_geometry"]


def test_scan_end_duplicating_start_across_180_wrap_is_reported():
    """Dihedral 179 -> -179 is the same angle through the +-180 wrap:
    the last point equals the first, a wave-23 cluster failure."""
    geom = ("C 0.0 0.0 0.0\nC 0.0 0.0 1.5\nC 1.2 0.0 2.3\nC 1.2 0.0 3.8\n")
    inp = _scan_inp("PBE0 def2-SVP OptTS", "{D 0 1 2 3 179.0 -179.0 12}S\n", geom)
    issues = [i for i in validate_orca_input(inp) if i.code == "scan_end_wraps_to_start"]
    assert issues, _codes(inp)


def test_scan_end_without_wrap_is_clean():
    geom = ("C 0.0 0.0 0.0\nC 0.0 0.0 1.5\nC 1.2 0.0 2.3\nC 1.2 0.0 3.8\n")
    inp = _scan_inp("PBE0 def2-SVP OptTS", "{D 0 1 2 3 90.0 120.0 12}S\n", geom)
    assert not [i for i in validate_orca_input(inp)
                if i.code == "scan_end_wraps_to_start"]
