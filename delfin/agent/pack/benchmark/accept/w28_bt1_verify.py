#!/usr/bin/env python3
"""Acceptance: w28_bt1_butane_rotation_profile.

Hard preconditions, from the artifact tree chem/butane left behind:
  - a scan input carrying an ORCA %geom Scan block on ONE dihedral and
    an output beside it;
  - a multi-frame xyz trajectory (the rotation movie) with >1 frames;
  - the anti and gauche minima confirmed by frequency outputs with
    ZERO imaginary modes;
  - exactly two TS outputs each with EXACTLY ONE imaginary mode
    (two imaginary modes is a REFUSAL-grade failure, never a pass);
  - the computed barrier in the literature band 4.5-5.8 kcal/mol and
    the anti-gauche gap in 0.6-1.15 kcal/mol, stated in summary.md
    (or any md in the folder) over the artifact;
  - a level of theory stated in the summary.

usage: w28_bt1_verify.py <workspace>
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from delfin.imag import collect_imaginary_modes

_BARRIER_BAND = (4.5, 5.8)
_GAP_BAND = (0.6, 1.15)

_BARRIER_RE = re.compile(
    r"(?:barrier|activation)(?:\s+\w+){0,10}?"
    r"(?::\s*|\s+(?:is|of)\s*)(?:about|approx\.?|roughly|~)?\s*"
    r"(-?\d+(?:\.\d+)?)\s*kcal", re.I)
_GAP_RE = re.compile(
    r"(?:gauche)(?:\s+\w+){0,10}?"
    r"(?::\s*|\s+(?:is|of)\s*)(?:about|approx\.?|roughly|~)?\s*"
    r"(-?\d+(?:\.\d+)?)\s*kcal", re.I)
_LEVEL_RE = re.compile(
    r"(?i)(r2scan-?3c|b3lyp|wb97x-?d4?|pbe0|revpbe0|dlpno|ccsd\(t\))")
_RANGE_SPAN_RE = re.compile(
    r"(-?\d+(?:\.\d+)?)\s*(?:-|\u2013|and)\s*(-?\d+(?:\.\d+)?)\s*kcal", re.I)


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return ""


def _in_cited_range(num: str, text: str, lo: float, hi: float) -> bool:
    v = float(num)
    for m in _RANGE_SPAN_RE.finditer(text):
        a, b = float(m.group(1)), float(m.group(2))
        if (lo <= a <= hi or lo <= b <= hi) and (a == v or b == v):
            return True
    return False


def _band_hit(text: str, pattern, lo: float, hi: float) -> bool:
    for m in pattern.finditer(text):
        if lo <= float(m.group(1)) <= hi and not _in_cited_range(
                m.group(1), text, lo, hi):
            return True
    return False


def _scan_artifacts(folder: Path):
    inputs, ran = [], False
    for inp in folder.rglob("*.inp"):
        text = _read(inp)
        if re.search(r"(?i)%geom", text) and re.search(
                r"(?i)\bD\s+\d+\s+\d+\s+\d+\s+\d+", text):
            inputs.append(inp)
            if inp.with_suffix(".out").is_file():
                ran = True
    return inputs, ran


def _xyz_frame_count(text: str) -> int:
    lines = text.strip().split("\n")
    n, i, frames = 0, 0, 0
    while i < len(lines):
        try:
            n = int(lines[i].strip())
        except (ValueError, IndexError):
            break
        i += 1 + n + 1
        frames += 1
    return frames


def _freq_outputs(folder: Path):
    outs = []
    for out in sorted(folder.rglob("*.out")):
        text = _read(out)
        if "VIBRATIONAL FREQUENCIES" in text or "***imaginary mode***" in text:
            outs.append(out)
    return outs


def _n_imag(out: Path) -> int:
    return len(collect_imaginary_modes(str(out)))


def main(argv) -> int:
    ws = Path(argv[1]).resolve() if len(argv) > 1 else Path.cwd()
    folder = ws / "chem" / "butane"
    if not folder.is_dir():
        print(f"FAIL: no task folder {folder}")
        print("QUALITY: 0.00")
        return 1
    scan_inputs, scan_ran = _scan_artifacts(folder)
    movie = None
    for xyz in sorted(folder.rglob("*.xyz")):
        frames = _xyz_frame_count(_read(xyz))
        if frames > 1:
            movie = (xyz.name, frames)
            break
    freq_outs = _freq_outputs(folder)
    zero_imag = [o for o in freq_outs if _n_imag(o) == 0]
    one_imag = [o for o in freq_outs if _n_imag(o) == 1]
    two_imag = [o for o in freq_outs if _n_imag(o) >= 2]
    prose = ""
    for md in sorted(folder.rglob("*.md")):
        prose += _read(md)
    ok_barrier = _band_hit(prose, _BARRIER_RE, *_BARRIER_BAND)
    ok_gap = _band_hit(prose, _GAP_RE, *_GAP_BAND)
    ok_level = bool(_LEVEL_RE.search(prose))
    preconditions = [
        (bool(scan_inputs), "an input with a one-coordinate %geom Scan block exists"),
        (scan_ran, "an output sits beside the scan input"),
        (movie is not None,
         "a multi-frame xyz (rotation movie) exists"
         if movie else "a multi-frame xyz (rotation movie) exists"),
        (len(zero_imag) >= 2,
         f"two minima with 0 imaginary modes (found {len(zero_imag)})"),
        (len(one_imag) >= 2,
         f"two TS with exactly 1 imaginary mode (found {len(one_imag)})"),
        (not two_imag, "no structure with 2+ imaginary modes is called a TS"),
        (ok_barrier,
         "computed barrier stated in 4.5-5.8 kcal/mol over the artifact"),
        (ok_gap, "computed gauche gap stated in 0.6-1.15 kcal/mol over the artifact"),
        (ok_level, "level of theory stated"),
    ]
    ok_all = True
    met = 0
    for ok, label in preconditions:
        print(("PASS" if ok else "FAIL") + ": " + label)
        ok_all = ok_all and ok
        met += 1 if ok else 0
    print(f"QUALITY: {met / len(preconditions):.2f}")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
