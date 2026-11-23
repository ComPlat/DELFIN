"""A guard that refuses its own parser gets exempted, not obeyed.

The guard's LC1 family refuses licensed QM-program output. Its signatures
are header and footer markers, and DELFIN's own parser holds one of them
in a string literal because it is LOOKING for it:

    return '****ORCA TERMINATED NORMALLY****' in text

Under the old per-line rule that line was a violation. Two files were
therefore put in SELF_EXEMPT -- one for a version number in an install
guide, one for a transcript corpus that mentions the phrase in prose --
and a whole-file exemption is strictly worse than a false positive: it
admits a genuine output block alongside the line it was written for.

So the rule is a THRESHOLD. Real output carries several distinct
signatures (banner, version line, host line, input echo, thermochemistry
block, footer); a parser or a sentence carries one. Below the threshold
the hit is an ADVISORY: reported, never silenced, because a signature
nobody sees is a signature nobody acts on.

The threshold is also what let a real gap close. scan_py_source checks
imports and CCDC strings and NOT the QM signatures, so a .py file with a
pasted ORCA output block was never scanned for one. Both scans now run on
every file, which is only safe because a parser no longer refuses it.

The honest limit, asserted below so nobody mistakes it for an oversight:
a paste of the footer line ALONE is one signature and is not refused. One
line of text cannot be distinguished from a parser looking for that line.
"""
from __future__ import annotations

import importlib.util
import pathlib

import pytest

_GUARD = (pathlib.Path(__file__).resolve().parents[1]
          / "scripts" / "license_guard.py")


@pytest.fixture(scope="module")
def lg():
    spec = importlib.util.spec_from_file_location("_lg_under_test", _GUARD)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


#: DELFIN's own code, verbatim in shape: it holds the marker to search for.
PARSER = ("def finished(text):\n"
          "    return '****ORCA TERMINATED NORMALLY****' in text\n")

#: Prose about a run, as a transcript corpus or a guide carries it.
PROSE = ("The step is done once ORCA terminated normally, which the\n"
         "parser checks before reading the energies.\n")

#: A real output block: banner, version, host, working dir, input echo,
#: thermochemistry and footer -- seven distinct signatures.
REAL_OUTPUT = """\
                                 *****************
                                 * O   R   C   A *
                                 *****************
           Program Version 6.1.0  -  RELEASE  -
* Host name: somehost
* Working dir.: /scratch/run
|  1> ! PBE0 def2-SVP
|  2> * xyz 0 1
THERMOCHEMISTRY AT 298.15K
                             ****ORCA TERMINATED NORMALLY****
"""


def test_the_projects_own_parser_is_not_a_violation(lg):
    advisories: list = []
    assert lg.scan_path_text("delfin/parse.py", PARSER, advisories) == []
    assert advisories, "the single hit must still be reported"
    assert "parser or prose" in advisories[0][1]


def test_prose_about_a_run_is_not_a_violation(lg):
    advisories: list = []
    assert lg.scan_path_text("docs/guide.txt", PROSE, advisories) == []
    assert advisories


def test_a_real_output_block_is_still_refused(lg):
    findings = lg.scan_path_text("x.txt", REAL_OUTPUT, [])
    assert findings, "a genuine output block must be refused"
    assert any("distinct signatures" in msg for _, msg in findings)


def test_the_threshold_is_what_the_rule_says_it_is(lg):
    """Two distinct signatures refuse; one does not. Asserted against the
    constant so a later edit of it cannot silently change the rule."""
    assert lg.QM_DISTINCT_REQUIRED == 2
    one = "           Program Version 6.1.0  -  RELEASE  -\n"
    two = one + "* Host name: somehost\n"
    assert lg.scan_path_text("a.txt", one, []) == []
    assert lg.scan_path_text("a.txt", two, []) != []


def test_a_py_file_with_pasted_output_is_now_scanned(lg):
    """The gap the threshold allowed closing: scan_py_source checks
    imports and CCDC strings only, so this file was never examined for
    QM output at all."""
    src = f'LOG = """\n{REAL_OUTPUT}"""\n'
    assert lg.scan_py_source("delfin/x.py", src) == [], (
        "scan_py_source is unchanged: it does not look for QM output")
    assert lg._scan_one("delfin/x.py", src, []), (
        "the combined scan must see it")


def test_a_footer_alone_is_not_refused_and_that_is_the_known_limit(lg):
    """One line is indistinguishable from a parser holding that line.
    Asserted so the limit is documented rather than discovered."""
    advisories: list = []
    footer = "                   ****ORCA TERMINATED NORMALLY****\n"
    assert lg.scan_path_text("x.txt", footer, advisories) == []
    assert advisories


def test_the_synthetic_marker_still_exempts_a_whole_fixture(lg):
    text = lg.QM_SYNTHETIC_MARKER + "\n" + REAL_OUTPUT
    assert lg.scan_path_text("tests/fixtures/f.out", text, []) == []


def test_the_allow_comment_still_works_line_by_line(lg):
    text = REAL_OUTPUT.replace(
        "* Host name: somehost",
        "* Host name: somehost  # license-guard: allow")
    findings = lg.scan_path_text("x.txt", text, [])
    assert not any("Host name" in msg for _, msg in findings)


def test_the_two_files_that_were_exempted_are_scanned_again(lg):
    """A whole-file exemption admits a real paste. Both come out clean
    under the threshold, so neither needs one."""
    exempt = lg.SELF_EXEMPT
    assert not any("DELFIN_INSTALLATION_ANLEITUNG" in p for p in exempt)
    assert not any("choice_offer_corpus" in p for p in exempt)
