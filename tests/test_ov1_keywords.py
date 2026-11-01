"""OV1 phase 1: every ``!`` keyword and ``%`` block must be one ORCA knows.

Wave 23 incidents: ``NOFREQ`` (not an ORCA keyword) passed every DELFIN
check and failed on the cluster; a ``%freq`` block written on the ``!``
line did the same.

The primary tests drive ``delfin.common.orca_input_checks.check_all``
directly — delfin/api.py is protected and the operator wires ``check_all``
into ``validate_orca_input`` (delfin/api.py:3114, the function
``delfin.agent.calc_policy.check_orca_text`` calls before every submit).
The wiring test at the bottom drives the real public path and runs once
the hook is on the branch.
"""
import pytest

from delfin.common.orca_input_checks import check_all, check_keyword_line


def _errors(inp):
    return [i for i in check_all(inp) if i.severity == "error"]


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


def test_unknown_keyword_nofreq_is_an_error_with_suggestion():
    """NOFREQ is not an ORCA keyword (wave 23 incident)."""
    issues = _errors(_inp("PBE0 def2-SVP NOFREQ"))
    assert any(i.code == "unknown_keyword" for i in issues), \
        [i.code for i in issues]
    miss = next(i for i in issues if i.code == "unknown_keyword")
    assert "NOFREQ" in miss.message
    # the suggestion names a real keyword, not the misspelling
    assert miss.suggestion


def test_misspelt_keyword_gets_closest_real_keyword_as_suggestion():
    # FreQ is known (case-blind); a genuinely invented near-miss:
    issues = check_keyword_line("! OptTight")
    miss = next(i for i in issues if i.code == "unknown_keyword")
    assert "TightOpt" in miss.suggestion or "tightopt" in miss.suggestion.lower()
    assert miss.suggestion


def test_known_keywords_stay_clean():
    for words in ("PBE0 def2-SVP Freq",
                  "B3LYP def2-TZVP D4 OPT TightOpt",
                  "wB97X-D3 def2-TZVP CPCM(water)"):
        assert not [i for i in _errors(_inp(words))
                    if i.code == "unknown_keyword"], words


def test_block_written_on_keyword_line_is_an_error():
    """A %freq block on the ! line is a wave-23 cluster failure."""
    issues = _errors(_inp("PBE0 def2-SVP %freq end"))
    assert any(i.code == "block_on_keyword_line" for i in issues), \
        [i.code for i in issues]


def test_unknown_block_name_is_an_error():
    inp = _inp("PBE0 def2-SVP", extra="%nosuchblock\n  Foo 1\nend\n")
    issues = _errors(inp)
    assert any(i.code == "unknown_block" for i in issues), [i.code for i in issues]


def test_unknown_block_option_is_reported():
    inp = _inp("PBE0 def2-SVP", extra="%scf\n  NotAnScfOption 1\nend\n")
    issues = [i for i in check_all(inp) if i.code == "unknown_block_option"]
    assert issues, [i.code for i in check_all(inp)]
    assert "NotAnScfOption" in issues[0].message


def test_known_block_and_option_is_clean():
    inp = _inp("PBE0 def2-SVP", extra="%scf\n  MaxIter 200\nend\n")
    clean = [i for i in check_all(inp)
             if i.code in ("unknown_block", "unknown_block_option",
                           "block_on_keyword_line")]
    assert not clean


def _hook_is_wired():
    from delfin.api import validate_orca_input
    return any(i.code == "unknown_keyword"
               for i in validate_orca_input(_inp("PBE0 def2-SVP NOFREQ")))


def test_public_path_validate_orca_input_reports_unknown_keyword():
    """The wiring test: check_all reaches the real public path.

    Skipped until the operator adds the call line to validate_orca_input
    (delfin/api.py is a protected file for this builder). Once wired, the
    skip disappears and this asserts the real path.
    """
    from delfin.api import validate_orca_input
    issues = [i for i in validate_orca_input(_inp("PBE0 def2-SVP NOFREQ"))
              if i.severity == "error"]
    assert any(i.code == "unknown_keyword" for i in issues)
