"""Adversarial tests for package OV1 (reviewer nacht-s14).

Baseline state: written against the UNCHANGED code at base 54462fe67.
Cases that are red today name a defect the builder must fix; they are
reported, not committed as green. Each case drives the real public path:
delfin.api.validate_orca_input (the function calc_policy.py:517 calls
before every submit).
"""
from delfin.api import validate_orca_input

# A minimal valid input: keyword line, geometry block, PAL block.
HEADER = "! PBE0 def2-SVP Opt Freq\n"
GEOM = (
    "* xyz 0 1\n"
    "C   0.0 0.0 0.0\n"
    "H   1.089 0.0 0.0\n"
    "H   -0.363 1.027 0.0\n"
    "H   -0.363 -1.027 0.0\n"
    "H   -0.363 0.0 1.027\n"
    "*\n"
)
FOOTER = "%pal nprocs 4 end\n%maxcore 2000\n"


def _input(keywords):
    return HEADER.replace("! PBE0 def2-SVP Opt Freq\n", keywords) + GEOM + FOOTER


def _errors(text):
    return [i for i in validate_orca_input(text) if i.severity == "error"]


def test_nofreq_invented_keyword_is_an_error():
    """Wave-23 cluster failure: NOFREQ is not an ORCA keyword."""
    issues = _errors(_input("! PBE0 def2-SVP NOFREQ Opt\n"))
    assert any("NOFREQ" in i.message or "NOFREQ" in i.suggestion
               for i in issues), "invented keyword NOFREQ passed"


def test_suggestion_names_closest_real_keyword():
    """A misspelt keyword gets the closest real keyword as suggestion.

    Rewritten per QS ruling (rule 6): the old decision was a
    title-cased user display ("TightOpt"); the new decision is the
    catalog's casing ('tightopt', opi_keywords.json is all-lowercase)
    matched case-insensitively — the catalog is the single source of
    truth and a hand-written case map would drift from it.
    """
    issues = _errors(_input("! PBE0 def2-SVP OptTight Freq\n"))
    hits = [i for i in issues if "OptTight" in i.message]
    assert hits, "misspelt keyword OptTight not flagged"
    assert any("tightopt" in i.suggestion.lower() for i in hits), \
        "no closest-keyword suggestion (expected catalog word 'tightopt')"


def test_block_written_on_bang_line_is_an_error():
    """%freq on the ! line is not a keyword; ORCA rejects the input."""
    issues = _errors(_input("! PBE0 def2-SVP %freq end Opt Freq\n"))
    assert issues, "%freq block on the ! line passed silently"


def test_numfreq_and_freq_conflict_is_an_error():
    """ORCA 6 manual: NumFreq and Freq cannot both be requested."""
    issues = _errors(_input("! PBE0 def2-SVP Opt NumFreq Freq\n"))
    assert any("NumFreq" in i.message and "Freq" in i.message
               for i in issues), "NumFreq + Freq conflict not detected"


def test_ri_without_auxiliary_basis_is_an_error():
    """RI requires a matching /C or /J auxiliary basis."""
    issues = _errors(_input("! PBE0 RI def2-SVP Opt\n"))
    assert issues, "RI without an auxiliary basis passed"


def test_unknown_cpcm_solvent_is_an_error():
    """CPCM(SolventName) with a name ORCA has no parameters for."""
    issues = _errors(_input("! PBE0 def2-SVP Opt CPCM(NoSuchSolventX)\n"))
    assert any("NoSuchSolventX" in i.message for i in issues), \
        "unknown CPCM solvent passed silently"


def test_known_solvent_still_accepted():
    """Nothing new is refused: CPCM(water) must produce no error."""
    assert not _errors(_input("! PBE0 def2-SVP Opt CPCM(water)\n"))


def test_plain_valid_input_unchanged():
    """The ordinary input still validates clean (regression guard)."""
    assert not _errors(_input("! PBE0 def2-SVP Opt Freq\n"))
