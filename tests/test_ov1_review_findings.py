"""OV1: reviewer (nacht-s14) and QS (nacht-s11) findings as tests.

F-REV-1 fullwidth '!' normalization, F-REV-2 nested blocks / data lines
in the option scan, F-REV-3 gap-block options still checked, F-REV-4 +
QS suggestion casing/ranking. Each red before the fix in the phase-2
commit.
"""
from delfin.common.orca_input_checks import (
    _closest,
    check_all,
    check_blocks,
    check_keyword_line,
    check_scan_block,
)
from delfin.api import validate_orca_input


def test_fullwidth_exclamation_line_is_flagged():
    """F-REV-1: U+FF01 look-alike must be treated as '!'."""
    inp = "\uff01 PBE0 def2-SVP NOFREQ\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n"
    issues = validate_orca_input(inp)
    assert any(i.code == "unknown_keyword" for i in issues), \
        [i.code for i in issues]


def test_nested_block_names_not_scanned_as_outer_options():
    """F-REV-2: %pointcharges inside %cpcm, its data lines are numbers."""
    inp = ("! PBE0 def2-SVP\n"
           "%cpcm\n"
           "  %pointcharges \"pc.txt\"\n"
           "  1 0.0 0.0 0.0 1.0\n"
           "end\n"
           "* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")
    issues = [i for i in check_all(inp) if i.code == "unknown_block_option"]
    assert not issues, [(i.message) for i in issues]


def test_gap_block_option_bogus_still_flagged():
    """F-REV-3: the %pal name exemption must not exempt its options."""
    inp = ("! PBE0 def2-SVP\n"
           "%pal\n  nprocs 2\n  BogusOpt 1\nend\n"
           "* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n")
    issues = [i for i in check_all(inp) if i.code == "unknown_block_option"]
    assert any("BogusOpt" in i.message for i in issues), \
        [(i.message) for i in issues]


def test_suggestion_is_tightopt_not_opth():
    """F-REV-4 + QS: the closest-keyword suggestion must be the real
    keyword made of the same tokens."""
    issues = check_keyword_line("! OptTight")
    miss = next(i for i in issues if i.code == "unknown_keyword")
    assert "tightopt" in miss.suggestion.lower(), miss.suggestion


def test_suggestion_ranking_direct():
    assert _closest(["opt", "opth", "tightopt", "optts"], "OptTight") \
        == "tightopt"


def test_scan_block_geometry_before_scan_is_read():
    """Geometry sits BEFORE %scan in a normal input; the scan checks
    must still see it (start-value comparison needs the atoms)."""
    inp = ("! PBE0 def2-SVP Opt\n"
           "* xyz 0 1\n"
           "H 0.0 0.0 0.0\nH 0.0 0.0 0.74\nH 0.0 0.0 1.90\n"
           "*\n"
           "%scan\n{B 2 0 3.8 2.0 10}S\nend\n")
    issues = [i for i in check_scan_block(inp)
              if i.code == "scan_start_far_from_geometry"]
    assert issues, [(i.code, i.message) for i in check_scan_block(inp)]
