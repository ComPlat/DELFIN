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


def barrier_band(text: str, lo: float, hi: float) -> bool:
    """A barrier number inside the band appears in the prose (kcal/mol)."""
    import re
    for m in re.finditer(r"(-?\d+(?:\.\d+)?)\s*kcal", text, re.I):
        v = float(m.group(1))
        if lo <= v <= hi:
            return True
    return False


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
