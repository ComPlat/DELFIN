"""Adversarial tests, OV1 phase 1 (reviewer nacht-s14).

Attacks the phase-1 implementation in delfin/common/orca_input_checks.py
(check_keyword_line / check_blocks) on edge cases the builder's own
tests do not cover: unicode look-alikes, attached/invisible separators,
compound keywords with parameters, nested blocks, case, empty input.
Each case states what it expects and why; a red case is a finding to
the builder, not silently committed as green.
"""
from delfin.common.orca_input_checks import check_all

GEOM = (
    "* xyz 0 1\n"
    "C   0.0 0.0 0.0\n"
    "H   1.089 0.0 0.0\n"
    "*\n"
)


def _input(kwline, blocks=""):
    return kwline + GEOM + blocks


def _errors(issues):
    return [i for i in issues if i.severity == "error"]


def test_empty_input_no_crash():
    assert check_all("") == []
    assert check_all("just some text\n") == []


def test_lowercase_keyword_line_accepted():
    """ORCA is case-blind: `! pbe0 def2-svp opt` must not be flagged."""
    assert _errors(check_all(_input("! pbe0 def2-svp opt freq\n"))) == []


def test_unicode_nbsp_inside_keyword():
    """A keyword split by a non-breaking space is not two keywords.

    ORCA reads the line literally: `Opt\\u00a0Tight` is garbage, but
    str.split() sees whitespace and the pieces are both real keywords,
    so the check must not silently pass the input.
    """
    issues = check_all(_input("! PBE0 def2-SVP Opt\u00a0Tight Freq\n"))
    assert _errors(issues), "NBSP-hidden keyword passed silently"


def test_fullwidth_exclamation_line_is_flagged():
    """U+FF01 looks like `!` to a human but not to ORCA's parser."""
    issues = check_all(_input("\uff01 PBE0 def2-SVP Opt\n"))
    assert _errors(issues), "fullwidth ! accepted as a keyword line"


def test_bang_keyword_with_solvent_parameter_known():
    """CPCM(water) on the ! line is a valid keyword form."""
    assert _errors(check_all(_input("! PBE0 def2-SVP CPCM(water) Opt\n"))) == []


def test_bang_keyword_with_unknown_solvent_parameter():
    """CPCM(NoSuchSolventX): known form, unknown name inside.

    known_keyword() currently returns True for any CPCM(...) form; the
    name inside must still be validated (phase-2 item F2). Recorded as
    an explicit pending case: red until phase 2 lands.
    """
    issues = check_all(_input("! PBE0 def2-SVP CPCM(NoSuchSolventX) Opt\n"))
    assert _errors(issues), "unknown solvent inside CPCM() passed"


def test_nested_block_names_not_scanned_as_outer_options():
    """A block nested inside another must not be read as an option.

    %cpcm may carry nested sub-blocks; scanning their names as options
    of the outer block is a false warning the submit path must not show.
    """
    text = _input(
        "! PBE0 def2-SVP Opt\n",
        "%cpcm\n"
        "  %pointcharges\n"
        "    1 0.0 0.0 0.0 0.1\n"
        "  end\n"
        "end\n",
    )
    warnings = [i for i in check_all(text) if i.severity == "warning"]
    assert not warnings, f"nested block read as option: {warnings}"


def test_unknown_block_option_is_only_a_warning_not_an_error():
    """ORCA itself warns and continues on unknown block options."""
    text = _input("! PBE0 def2-SVP Opt\n", "%pal\n  NProcs 4\n  BogusOpt 1\nend\n")
    issues = check_all(text)
    assert _errors(issues) == [], "unknown block option escalated to error"
    assert [i for i in issues if i.severity == "warning"], \
        "unknown block option not even warned"


def test_maxcore_pal_blocks_not_flagged():
    text = _input("! PBE0 def2-SVP Opt\n", "%maxcore 2000\n%pal nprocs 4 end\n")
    assert _errors(check_all(text)) == [], "gap blocks flagged despite api.py owning them"


def test_block_without_end_reports_unknown_option_noise():
    """A block whose end is missing: option scan must not run to EOF.

    Deliberately pending phase-2 decision: this documents current
    behavior; if it floods warnings on truncated inputs the builder
    must bound the scan.
    """
    text = _input("! PBE0 def2-SVP Opt\n", "%scf\n  TolE 1e-7\n")
    issues = check_all(text)
    assert _errors(issues) == []
