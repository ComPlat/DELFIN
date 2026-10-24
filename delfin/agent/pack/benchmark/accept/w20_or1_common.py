#!/usr/bin/env python3
"""Shared grading helpers for the OR1 acceptance scripts.

The verdict comes from the ARTIFACT the run left behind (an ORCA
frequency output in the task's folder), parsed with DELFIN's own
imaginary-mode parser (delfin.imag.collect_imaginary_modes), and from
the prose the run wrote (checked by the bench runner's
expected_signals / expected_values machinery, not here).

Every script prints a QUALITY: n line (share of preconditions met)
and exits 0 only when the hard preconditions hold, in the chem_*
pattern.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from delfin.imag import collect_imaginary_modes


def find_freq_outputs(folder: Path):
    """Every *.out under folder that carries a frequency block."""
    hits = []
    for out in sorted(folder.rglob("*.out")):
        try:
            text = out.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if "VIBRATIONAL FREQUENCIES" in text or "***imaginary mode***" in text:
            hits.append(out)
    return hits


def one_imaginary_mode(folder: Path):
    """True iff some output in the folder shows exactly one imaginary mode.

    Returns (ok, evidence_line).
    """
    outs = find_freq_outputs(folder)
    if not outs:
        return False, "no frequency output found"
    best = None
    for out in outs:
        modes = collect_imaginary_modes(str(out))
        if len(modes) == 1:
            return True, f"{out.name}: 1 imaginary mode at {modes[0][1]:.1f} cm-1"
        if best is None or (modes and len(modes) < best[0]):
            best = (len(modes), out.name)
    if best:
        return False, f"{best[1]}: {best[0]} imaginary modes (need exactly 1)"
    return False, "frequency outputs found but no imaginary-mode lines parsed"


_BARRIER_RE = re.compile(
    r"(?:barrier|activation)(?:\s+\w+){0,10}?"
    r"(?::\s*|\s+(?:is|of)\s*)(?:about|approx\.?|roughly|~)?\s*"
    r"(-?\d+(?:\.\d+)?)\s*kcal", re.I)
# A bare "<num> kcal" span that is part of a RANGE citation (lo-hi,
# en-dash or written-out "between A and B") is never a computed value.
_RANGE_SPAN_RE = re.compile(
    r"(-?\d+(?:\.\d+)?)\s*(?:-|\u2013|and)\s*(-?\d+(?:\.\d+)?)\s*kcal", re.I)


def _is_in_cited_range(num: str, text: str, lo: float, hi: float) -> bool:
    v = float(num)
    for m in _RANGE_SPAN_RE.finditer(text):
        a, b = float(m.group(1)), float(m.group(2))
        if (lo <= a <= hi or lo <= b <= hi) and (
                a == v or b == v):
            return True
    return False


def barrier_band(text: str, lo: float, hi: float) -> bool:
    """A COMPUTED barrier number inside the band appears in the prose.

    Requires a barrier/activation phrasing ("barrier is/of/to ..."),
    and rejects band-endpoint citations: a number that only matches as
    one endpoint of a cited range ("12-20 kcal/mol", "between 12 and
    20") does not count.
    """
    for m in _BARRIER_RE.finditer(text):
        v = float(m.group(1))
        if not (lo <= v <= hi):
            continue
        if _is_in_cited_range(m.group(1), text, lo, hi):
            continue
        return True
    return False


def barrier_with_artifact(folder: Path, text: str, lo: float,
                          hi: float, freq_required: bool = True) -> bool:
    """barrier_band gated on an artifact in folder.

    Evidence before claim: no artifact in the task folder -> the
    barrier is never graded as passing, whatever the prose says.
    freq_required (the TS tasks): a frequency output must exist. The
    scan task passes freq_required=False: its artifact is the scan
    output, which the script's own scan_ran precondition already gates.
    """
    if freq_required and not find_freq_outputs(folder):
        return False
    return barrier_band(text, lo, hi)


def finish(preconditions: list[tuple[bool, str]]) -> int:
    ok_all = True
    met = 0
    for ok, label in preconditions:
        print(("PASS" if ok else "FAIL") + ": " + label)
        ok_all = ok_all and ok
        met += 1 if ok else 0
    n = len(preconditions)
    quality = met / n if n else 0.0
    print(f"QUALITY: {quality:.2f}")
    return 0 if ok_all else 1


def workspace(argv) -> Path:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        raise SystemExit(2)
    return Path(argv[1]).resolve()
