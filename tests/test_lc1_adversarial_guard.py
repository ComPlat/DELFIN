"""LC1 adversarial review tests (nacht-s2) — guard rule family edge cases.

These attack the phase-2a guard extension in scripts/license_guard.py.
They test the PUBLIC path: run_paths() over files written into a tmp dir.
"""

import textwrap
from pathlib import Path

import pytest

import sys

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
import license_guard  # noqa: E402

MARKER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output"
ORCA_VERSION_LINE = " Program Version 6.1.1  -  RELEASE   -"


@pytest.fixture
def scan(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(license_guard, "_report", lambda violations: violations)

    def run(name: str, body: str) -> list[str]:
        p = Path("fixture") / name
        p.parent.mkdir(exist_ok=True)
        p.write_text(body, encoding="utf-8")
        return license_guard.run_paths([str(p)])

    return run


def test_marked_synthetic_fixture_whose_body_quotes_orca_format_is_clean(scan):
    """A first-line-marked synthetic fixture may contain ORCA-format lines
    (the version line, a termination phrase) — that is the point of the
    exception. Red on the merged 8c06f305d: only the marker LINE is
    skipped, so the body's version line is flagged."""
    body = MARKER + "\n" + ORCA_VERSION_LINE + "\nORCA terminated normally\n"
    out = scan("synthetic.out", body)
    assert not out, out


def test_marker_only_on_later_line_does_not_exempt_real_signatures(scan):
    """The marker must be the FIRST line to exempt a file. A marker dropped
    in the middle of (what claims to be) real ORCA output must not blind
    the scanner to the surrounding signatures."""
    body = (
        " Total Energy      ....    -76.0100000\n"
        + ORCA_VERSION_LINE + "\n"
        + MARKER + "\n"
        "ORCA TERMINATED NORMALLY\n"
    )
    out = scan("suspicious.out", body)
    joined = "\n".join(out)
    # a mid-file marker must not blind the scanner to the whole file
    assert "ORCA TERMINATED NORMALLY" in joined, out
    # the REAL ORCA 6.1.1 version line has TWO spaces before the dash:
    # 'Program Version 6.1.1  -  RELEASE'. Red on 8c06f305d: the regex
    # QM_OUTPUT_SIGNATURES requires exactly '6.x.y -' with one space, so
    # the real-format line is missed.
    assert any("Program Version" in v for v in out), out


def test_signature_with_permuted_whitespace_is_still_caught(scan):
    """Evasion probe: real ORCA banners/termination lines are often
    re-wrapped; extra internal whitespace must not hide the signature."""
    out = scan("wrapped.out", "  ORCA  TERMINATED  NORMALLY\n")
    assert any("TERMINATED" in v for v in out), out


def test_lowercase_and_truncated_termination_variants_are_caught(scan):
    """The brief requires the guard to catch ORCA output, not one spelling:
    'ORCA terminated normally' in prose and 'ORCA TERMINATED NORMALLY AND'
    continuation lines are the same DATA."""
    out = scan("variants.out", "The run ended: ORCA terminated normally and exited.\n")
    assert any("terminated" in v.lower() for v in out), out
