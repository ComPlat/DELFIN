"""OC1 F1 red control: find_G must return the LAST valid Gibbs energy.

run_OCCUPIER binds find_G when OCCUPIER_compare=G (delfin/occupier.py:1135-1137)
and parses every FoB output with it (occupier.py:1359, 1386, 1448). The FSPE
counterpart find_FSPE scans bottom-up and keeps scanning past a value that does
not convert (occupier.py:996-1006). find_G did neither: it returned the FIRST
Gibbs line and returned None immediately on one unconvertible value, so an
output holding an aborted thermochemistry block (or a resumed/continued run
that appends a corrected block) was parsed with a stale energy.

These fixtures are what ORCA really writes: 'Final Gibbs free energy   ... '
lines (see tests/fixtures for real Freq outputs; the marker string here is the
exact one find_G searches for).
"""
import delfin.occupier as occ


def _write(tmp_path, text):
    p = tmp_path / "output.out"
    p.write_text(text, encoding="utf-8")
    return str(p)


MARKER = "Final Gibbs free energy         ... "


def test_first_gibbs_line_still_found(tmp_path):
    # Harmless output: one thermochemistry block -> same result either way.
    text = (
        "some header\n"
        + MARKER + "      -1234.567890 Eh\n"
        + "trailing notes\n"
    )
    assert occ.find_G(_write(tmp_path, text)) == -1234.567890


def test_last_gibbs_line_wins_when_output_appends(tmp_path):
    # A continued/resumed output appends a NEW thermochemistry block.
    # The stale first value must not be reported.
    text = (
        MARKER + "      -1234.567890 Eh\n"   # aborted/stale block
        + "ORCA finished by error termination\n"
        + "---- new run appended ----\n"
        + MARKER + "      -1234.999999 Eh\n"  # final, valid block
    )
    assert occ.find_G(_write(tmp_path, text)) == -1234.999999


def test_unconvertible_value_is_skipped_not_fatal(tmp_path):
    # One malformed line between two good ones must not poison the parse.
    text = (
        MARKER + "      -1234.567890 Eh\n"
        + MARKER + "      n/a Eh\n"
        + MARKER + "      -1234.111111 Eh\n"
    )
    assert occ.find_G(_write(tmp_path, text)) == -1234.111111


def test_no_gibbs_line_returns_none(tmp_path):
    assert occ.find_G(_write(tmp_path, "no thermo here\n")) is None


def test_missing_file_returns_none(tmp_path):
    assert occ.find_G(str(tmp_path / "does_not_exist.out")) is None
