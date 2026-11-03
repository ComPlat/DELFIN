"""Survey how DELFIN workflows were configured and whether they succeeded.

Built for the benchmark "configure a workflow correctly": from an archive of
DELFIN run directories, one record per CONTROL.txt describing

* the configuration the user chose (every key DELFIN understands; unknown
  keys kept separately),
* the system the run was about (system class: metals with counts, charge,
  multiplicity, number of atoms, element set -- never coordinates),
* the outcome (did the workflow finish, which steps failed).

Privacy by construction: only element symbols, counts and settings are read.
No coordinates, no file or molecule names, no user names. The directory path
is hashed and never written. Files over 50 MB are skipped; nothing outside
ROOT is followed through symlinks; nothing is written except --out.

CLI:  python -m delfin.orca_knowledge.control_survey ROOT --out FILE.jsonl
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# What DELFIN understands
# ---------------------------------------------------------------------------

from delfin.config import _CONTROL_KEY_ALIASES, parse_control_text
from delfin.common.control_validator import CONTROL_FIELD_SPECS

KNOWN_KEYS: frozenset = frozenset(
    set(f.name for f in CONTROL_FIELD_SPECS) | set(_CONTROL_KEY_ALIASES)
)

# Metals DELFIN's workflows are steered around (3d/4d/5d + f-block).
METAL_SYMBOLS: frozenset = frozenset({
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",
    "Ho", "Er", "Tm", "Yb", "Ac", "Th", "Pa", "U", "Np", "Pu",
})

MAX_FILE_SIZE = 50 * 1024 * 1024

_XYZ_SYMBOL_RE = re.compile(r"^[A-Za-z]{1,3}$")
# "… ORCA TERMINATED NORMALLY" / "ORCA finished by error termination"
# (delfin/api.py:441 uses the same two markers).
_NORMAL_RE = re.compile(r"ORCA TERMINATED NORMALLY", re.IGNORECASE)
_ERROR_RE = re.compile(r"ORCA finished by error termination|ABORTING THE RUN",
                       re.IGNORECASE)
# ORCA records charge and multiplicity in the input line "* xyz … charge mult".
_CHARGE_MULT_RE = re.compile(
    r"^\s*\*\s+xyz\b.*?(-?\d+)\s+(\d+)\s*$", re.MULTILINE)


# ---------------------------------------------------------------------------
# parse
# ---------------------------------------------------------------------------

def parse(control_path: str) -> Dict[str, Any]:
    """Parse one CONTROL.txt into everything DELFIN understands.

    Returns {"config": {canonical key -> value}, "unknown": {key -> value}}.
    Aliases are canonicalised exactly as DELFIN canonicalises them
    (delfin/config.py: `_canonicalize_control_key`). Raises FileNotFoundError
    if the file does not exist.
    """
    path = Path(control_path)
    if not path.is_file():
        raise FileNotFoundError(f"no such CONTROL.txt: {control_path}")
    if path.stat().st_size > MAX_FILE_SIZE:
        raise ValueError(f"CONTROL.txt larger than {MAX_FILE_SIZE} bytes: {control_path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    raw = parse_control_text(text)
    config: Dict[str, Any] = {}
    unknown: Dict[str, Any] = {}
    for key, value in raw.items():
        if key.startswith("_"):  # parser-internal bookkeeping
            continue
        if key in KNOWN_KEYS:
            config[key] = _jsonable(value)
        else:
            unknown[key] = _jsonable(value)
    return {"config": config, "unknown": unknown}


def _jsonable(value: Any) -> Any:
    if isinstance(value, (str, int, float, bool)) or value is None:
        return _sanitize_filename(value) if isinstance(value, str) else value
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    return str(value)


_FILENAME_RE = re.compile(r"^.+\.(xyz|inp|txt|smi|smiles|mol|sdf|gbw|out|log)$",
                          re.IGNORECASE)


def _sanitize_filename(value: str) -> str:
    # CONTROL values may name the run's input files. The name itself can
    # carry the molecule's identity ("complex_aligned.xyz"); the fact that a
    # file was named is what belongs in a configuration survey.
    if _FILENAME_RE.match(value.strip()):
        return "<filename>"
    return value


# ---------------------------------------------------------------------------
# system
# ---------------------------------------------------------------------------

def _element_column(run_dir: Path) -> Optional[List[str]]:
    """Element symbols of the run's structure, from its own xyz file.

    Only the symbol column is read -- never the coordinate columns.
    """
    for name in ("input.xyz", "initial.xyz", "start.xyz"):
        p = run_dir / name
        if p.is_file() and p.stat().st_size <= MAX_FILE_SIZE:
            syms = _symbols_from_xyz(p)
            if syms:
                return syms
    # DEPRECATED path: none of the standard names matched; do not guess.
    return None


def _symbols_from_xyz(path: Path) -> List[str]:
    symbols: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        try:
            natoms_line = fh.readline().strip()
            fh.readline()  # comment line
            natoms = int(natoms_line.split()[0]) if natoms_line else -1
        except ValueError:
            return []
        for _ in range(natoms):
            parts = fh.readline().split()
            if not parts:
                break
            sym = parts[0]
            if not _XYZ_SYMBOL_RE.match(sym):
                return []  # not a plain xyz frame -- refuse rather than guess
            symbols.append(sym.capitalize())
    return symbols


def system(run_dir: str) -> Dict[str, Any]:
    """System class of one run, from the run's own files.

    Returns metals with counts, charge, multiplicity (as DELFIN or ORCA
    recorded them), number of atoms and the element set. No coordinates,
    no molecule or file names, no user names.
    """
    d = Path(run_dir)
    out: Dict[str, Any] = {
        "n_atoms": None,
        "elements": [],
        "metals": {},
        "charge": None,
        "multiplicity": None,
        "charge_source": None,
        "multiplicity_source": None,
    }
    symbols = _element_column(d)
    if symbols:
        out["n_atoms"] = len(symbols)
        out["elements"] = sorted(set(symbols))
        out["metals"] = dict(
            sorted((s, symbols.count(s)) for s in set(symbols) if s in METAL_SYMBOLS))
    # Charge/multiplicity as DELFIN recorded it (CONTROL of this run).
    control_path = d / "CONTROL.txt"
    if control_path.is_file():
        try:
            cfg = parse_control_text(
                control_path.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            cfg = {}
        ch = cfg.get("charge")
        if ch not in (None, ""):
            try:
                out["charge"] = int(str(ch))
                out["charge_source"] = "CONTROL.txt"
            except ValueError:
                pass
        for key in ("multiplicity_global_opt", "multiplicity_0"):
            mv = cfg.get(key)
            if mv not in (None, "", "''", '""'):
                try:
                    out["multiplicity"] = int(str(mv).strip("'\" "))
                    out["multiplicity_source"] = f"CONTROL.txt:{key}"
                    break
                except ValueError:
                    continue
    # Fall back to what ORCA recorded in its own input line.
    if out["charge"] is None or out["multiplicity"] is None:
        orca_charge, orca_mult = _charge_multiplicity_from_orca(d)
        if out["charge"] is None and orca_charge is not None:
            out["charge"] = orca_charge
            out["charge_source"] = "orca_input"
        if out["multiplicity"] is None and orca_mult is not None:
            out["multiplicity"] = orca_mult
            out["multiplicity_source"] = "orca_input"
    return out


def _charge_multiplicity_from_orca(d: Path) -> Tuple[Optional[int], Optional[int]]:
    """Charge/multiplicity from the '* xyz ... charge mult' line of an ORCA input."""
    for name in ("input.inp", "initial.inp"):
        p = d / name
        if not p.is_file() or p.stat().st_size > MAX_FILE_SIZE:
            continue
        text = p.read_text(encoding="utf-8", errors="replace")
        m = _CHARGE_MULT_RE.search(text)
        if m:
            return int(m.group(1)), int(m.group(2))
    return None, None


# ---------------------------------------------------------------------------
# outcome
# ---------------------------------------------------------------------------

def _out_dirs(run_dir: Path) -> List[Path]:
    """The run directory and its immediate subdirectories that hold ORCA output."""
    dirs = [run_dir]
    try:
        for child in sorted(run_dir.iterdir()):
            if child.is_dir() and not child.is_symlink():
                dirs.append(child)
    except OSError:
        pass
    return dirs


def outcome(run_dir: str) -> Dict[str, Any]:
    """Did the workflow finish, and which steps failed.

    finished: DELFIN's final report/summary (DELFIN.txt or DELFIN_Data.json)
    exists in the run tree, AND every ORCA output found terminated normally
    (a run with no ORCA output counts on the DELFIN summary alone).
    failed_steps: basenames (stem) of outputs that did NOT terminate normally.
    """
    d = Path(run_dir)
    has_summary = (d / "DELFIN.txt").is_file() or (d / "DELFIN_Data.json").is_file()
    steps_total = 0
    normal = 0
    failed_steps: List[str] = []
    aborted_steps: List[str] = []
    for od in _out_dirs(d):
        try:
            outs = sorted(p for p in od.iterdir()
                          if p.is_file() and p.name.endswith(".out")
                          and not p.is_symlink()
                          and p.stat().st_size <= MAX_FILE_SIZE
                          # delfin_<jobid>.out is the SLURM wrapper's own
                          # stdout (submit_delfin.sh: --output=delfin_%j.out);
                          # it is not an ORCA output and never carries an
                          # ORCA termination line.
                          and not p.name.startswith("delfin_")
                          # .oldN files are kept copies of an earlier attempt,
                          # superseded by the final run of the same step.
                          and not re.search(r"\.old\d*$", p.stem))
        except OSError:
            continue
        top = od == d
        for out in outs:
            steps_total += 1
            try:
                tail = _tail(out, 200_000)
            except OSError:
                failed_steps.append(_step_name(out, top))
                continue
            if _NORMAL_RE.search(tail):
                normal += 1
            elif _ERROR_RE.search(tail):
                aborted_steps.append(_step_name(out, top))
            else:
                failed_steps.append(_step_name(out, top))
    truncated = steps_total == 0 and not has_summary
    return {
        "finished": bool(has_summary and steps_total == normal),
        "delfin_summary": has_summary,
        "orca_steps_total": steps_total,
        "orca_steps_normal": normal,
        "failed_steps": sorted(set(failed_steps)),
        "aborted_steps": sorted(set(aborted_steps)),
        "unverifiable": truncated,
    }


def _step_name(out: Path, top: bool) -> str:
    # Step name = the output file's stem; prefixed with the subdirectory
    # name only when DELFIN ran the step in its own subfolder.
    if not top and out.parent.name:
        return f"{out.parent.name}/{out.stem}"
    return out.stem


def _tail(path: Path, nbytes: int) -> str:
    with open(path, "rb") as fh:
        fh.seek(0, os.SEEK_END)
        size = fh.tell()
        fh.seek(max(0, size - nbytes))
        return fh.read().decode("utf-8", errors="replace")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _path_hash(path: Path, root: Path) -> str:
    rel = path.relative_to(root).as_posix() if path != root else ""
    return hashlib.sha256(rel.encode("utf-8")).hexdigest()


def _safe_iter(root: Path):
    """Yield every CONTROL.txt under root; read-only; no symlink escape."""
    stack = [root]
    while stack:
        d = stack.pop()
        if d.is_symlink():
            continue
        try:
            entries = sorted(d.iterdir())
        except OSError:
            continue
        for entry in entries:
            if entry.is_symlink():
                continue
            if entry.is_dir():
                stack.append(entry)
            elif entry.name == "CONTROL.txt" and entry.is_file():
                if entry.stat().st_size <= MAX_FILE_SIZE:
                    yield entry


def survey(root: str, out_path: str) -> int:
    """Run the survey over every CONTROL.txt under root; return record count."""
    rootp = Path(root).resolve()
    n = 0
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        for control in _safe_iter(rootp):
            run_dir = control.parent
            try:
                parsed = parse(str(control))
            except Exception:
                # Never put the exception (or the path inside it) into the record.
                record = {"hash": _path_hash(run_dir, rootp),
                          "control_error": "unreadable CONTROL.txt"}
            else:
                record = {
                    "hash": _path_hash(run_dir, rootp),
                    "system": system(str(run_dir)),
                    "control": parsed["config"],
                    "unknown": parsed["unknown"],
                    "outcome": outcome(str(run_dir)),
                }
            fh.write(json.dumps(record, ensure_ascii=True) + "\n")
            n += 1
    return n


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="delfin.orca_knowledge.control_survey",
        description="Survey DELFIN CONTROL.txt configurations and run outcomes "
                    "(read-only; writes only --out).")
    ap.add_argument("root", help="archive root holding DELFIN run directories")
    ap.add_argument("--out", required=True, help="output JSONL file")
    args = ap.parse_args(argv)
    n = survey(args.root, args.out)
    print(f"wrote {n} records to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
