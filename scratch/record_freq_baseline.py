#!/usr/bin/env python
"""Recorder for the orca_read frequency switch (wave24/s05).

Records the CURRENT outputs of every affected function on
(a) the fixtures in tests/fixtures/orca611 that carry a frequency run and
(b) up to 200 real runs sampled read-only from the run archive (copied
    first into ./scratch/samples/, never written in the archive).

Usage:
  python scratch/record_freq_baseline.py record   -> writes scratch/freq_baseline.json
  python scratch/record_freq_baseline.py compare  -> re-records and compares

Floats are serialised with repr() (round-trips bitwise). A comparison is
a pass only when every recorded value is identical.
"""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

WS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(WS))

from delfin import imag as _imag_mod
print(f"recording with delfin.imag from {_imag_mod.__file__}")
print(f"recording with delfin.orca_read from {_imag_mod.__file__ and __import__('delfin.orca_read', fromlist=['x']).__file__}")

from delfin.imag import (  # noqa: E402
    collect_imaginary_modes, hessian_frequencies, imaginary_modes,
    saddle_reason, search_imaginary_mode2,
)
from delfin.ir_spectrum import IRMode, calculate_spectrum_range, parse_ir_spectrum  # noqa: E402

ARCHIVE = Path("/home/qmchem_all/archive")
SAMPLES = WS / "scratch" / "samples"
BASELINE = WS / "scratch" / "freq_baseline.json"
MAX_RUNS = 200
DEFAULT_CONFIG = {"allow_imaginary_freq": -0.000647}


def find_runs() -> list[tuple[str, Path]]:
    """(relative run id, run dir) for fixture and sampled archive runs."""
    runs: list[tuple[str, Path]] = []
    for fix in sorted((WS / "tests" / "fixtures" / "orca611").iterdir()):
        if fix.is_dir() and any(fix.glob("*.out")):
            runs.append((f"fixture:{fix.name}", fix))
    # archive runs that have a .hess next to their .out (frequency runs)
    hess_files = sorted(ARCHIVE.rglob("*.hess")) if ARCHIVE.is_dir() else []
    seen = set()
    for hess in hess_files:
        run_dir = hess.parent
        base = hess.stem
        if not (run_dir / f"{base}.out").is_file():
            continue
        rel = run_dir.relative_to(ARCHIVE)
        if rel in seen:
            continue
        seen.add(rel)
        runs.append((f"archive:{rel.as_posix()}", run_dir))
        if len([r for r in runs if r[0].startswith("archive:")]) >= MAX_RUNS:
            break
    return runs


def copy_run(rel_id: str, run_dir: Path) -> Path | None:
    """Copy only the files the readers touch into scratch/samples/."""
    dest = SAMPLES / rel_id.replace(":", "/").replace("archive:/", "")
    dest = SAMPLES / rel_id.replace(":", "__")
    if (dest / ".copied").is_file():
        return dest
    out_files = sorted(run_dir.glob("*.out"))
    hess_files = sorted(run_dir.glob("*.hess"))
    if not out_files and not hess_files:
        return None
    dest.mkdir(parents=True, exist_ok=True)
    for src in out_files + hess_files:
        shutil.copy2(src, dest / src.name)
    (dest / ".copied").write_text(rel_id)
    return dest


def ir_mode_dict(m: IRMode) -> dict:
    return {
        "mode_number": m.mode_number,
        "frequency_cm1": repr(m.frequency_cm1),
        "epsilon": repr(m.epsilon),
        "intensity_km_mol": repr(m.intensity_km_mol),
        "t_squared": repr(m.t_squared),
        "tx": repr(m.tx), "ty": repr(m.ty), "tz": repr(m.tz),
        "transmittance_percent": repr(m.transmittance_percent),
    }


def record_run(run_dir: Path) -> dict:
    """Every affected function on one run directory."""
    rec: dict = {}
    for out_path in sorted(run_dir.glob("*.out")):
        stem = out_path.stem
        try:
            modes = parse_ir_spectrum(out_path)
            ir = {
                "n_modes": len(modes),
                "modes": [ir_mode_dict(m) for m in modes],
                "spectrum_range": [repr(v) for v in calculate_spectrum_range(modes)],
            }
        except Exception as exc:
            ir = {"error": f"{type(exc).__name__}: {exc}"}
        try:
            imag_out = collect_imaginary_modes(str(out_path))
            imag = [[i, repr(f)] for i, f in imag_out]
        except Exception as exc:
            imag = {"error": f"{type(exc).__name__}: {exc}"}
        try:
            search = repr(search_imaginary_mode2(str(out_path)))
        except Exception as exc:
            search = {"error": f"{type(exc).__name__}: {exc}"}
        rec[f"out:{stem}"] = {"ir": ir, "imag_out": imag, "search_imag": search}
    for hess_path in sorted(run_dir.glob("*.hess")):
        stem = hess_path.stem
        try:
            freqs = [repr(f) for f in hessian_frequencies(hess_path)]
        except Exception as exc:
            freqs = {"error": f"{type(exc).__name__}: {exc}"}
        try:
            im = [[i, repr(f)] for i, f in imaginary_modes(hess_path, DEFAULT_CONFIG)]
        except Exception as exc:
            im = {"error": f"{type(exc).__name__}: {exc}"}
        try:
            sr = saddle_reason(stem, hess_path, DEFAULT_CONFIG)
        except Exception as exc:
            sr = {"error": f"{type(exc).__name__}: {exc}"}
        rec[f"hess:{stem}"] = {"freqs": freqs, "imag_modes": im,
                               "saddle_reason": sr}
    return rec


def main() -> int:
    if sys.argv[1:2] == ["record"]:
        mode = "record"
    elif sys.argv[1:2] == ["compare"]:
        mode = "compare"
    else:
        print(__doc__)
        return 2
    runs = find_runs()
    copied: list[tuple[str, Path]] = []
    for rel_id, run_dir in runs:
        dest = copy_run(rel_id, run_dir)
        if dest is not None:
            copied.append((rel_id, dest))
    print(f"{len(copied)} runs (fixtures + archive sample) recorded from copies")
    records: dict = {}
    for rel_id, dest in copied:
        records[rel_id] = record_run(dest)
    if mode == "record":
        BASELINE.write_text(json.dumps(records, indent=1, sort_keys=True))
        print(f"baseline written: {BASELINE} ({len(records)} runs)")
        return 0
    # compare
    old = json.loads(BASELINE.read_text())
    if old == records:
        print(f"COMPARE PASS: {len(records)} runs, every recorded value identical")
        return 0
    diffs = []
    for run in sorted(set(old) | set(records)):
        if old.get(run) != records.get(run):
            diffs.append(run)
    print(f"COMPARE FAIL: {len(diffs)} runs differ, first: {diffs[:5]}")
    (WS / "scratch" / "freq_compare_fail.json").write_text(
        json.dumps({r: {"old": old.get(r), "new": records.get(r)}
                    for r in diffs}, indent=1, sort_keys=True))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
