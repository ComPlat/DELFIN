"""PV1: a complete, machine-readable provenance record for one result.

For every result (an energy, a property, a plot) produced from a run folder
this builds a ``PROVENANCE.json`` that names *how* the result came to be: the
input files it was made from (sha256 + path relative to the run), the program
and version printed in the output banner, the exact submit command when one is
recorded, the DELFIN commit it was built with (plus whether the tree was
dirty), the Python package versions, the structure origin (SMILES / file /
previous step id), a placeholder host class, and time stamps.  All in one JSON
document carrying a schema version, so a reader can tell the shape it is
looking at.

The record never contains a user name, a real host name or an absolute home
path (PR-public rule): absolute paths are turned into paths relative to the
run, the host class is always the placeholder ``host_class_redacted``, and
anything it cannot determine (a missing banner, an unknown origin) is recorded
as ``provenance_missing`` / ``unknown`` rather than guessed.

Reused, never re-implemented: ``delfin.cluster_bench.provenance`` for the
DELFIN code state and the interpreter's package versions, and the ORCA banner
version pattern from ``delfin.orca_read``.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

from delfin.cluster_bench.provenance import (
    cbatch_code_state,
    cbatch_probe_env,
    cbatch_sha256_file,
)

#: Schema version of a PROVENANCE.json; a reader refuses / re-migrates a file
#: whose version differs instead of mis-reading it.
SCHEMA_VERSION = "1"

#: The file name written next to a result.
PROVENANCE_FILE = "PROVENANCE.json"

#: Host class is always this placeholder: the record must never name a real
#: host (PR-public rule).  Env-specifics (login vs compute) may be added here
#: later without ever leaking a host name.
HOST_CLASS_PLACEHOLDER = "host_class_redacted"

#: Sentinels for properties the run does not let us determine: recorded as
#: such, never guessed.
MISSING = "provenance_missing"
UNKNOWN_ORIGIN = "unknown"

#: ORCA banner version pattern, the one delfin/orca_read uses to read it.
_VERSION_RE = re.compile(r"Program\s+Version\s+(\d+)\.(\d+)(?:\.(\d+))?")
#: A line that identifies the output as an ORCA run.
_ORCA_MARK = re.compile(r"\*\s+O\s+R\s+C\s+A\s+\*")

#: File types that describe the structure / the calculation, treated as input.
_INPUT_SUFFIXES = (".inp", ".xyz", ".smi", ".sdf", ".mol")
#: File name of the DELFIN pipeline state inside a run folder.
_DELFIN_DATA = "DELFIN_Data.json"
#: Suffixes of a recorded submit / job script in the run folder.
_SUBMIT_SUFFIXES = (".sh", ".slurm", ".sbatch", ".job")


#: Files ORCA/DELFIN generates during a run; never inputs.
_GENERATED_OUTPUT_NAMES = ("initial_trj.xyz",)
_GENERATED_OUTPUT_PREFIXES = ("initial",)
#: Delimiter sections in a CONTROL.txt file.
_CONTROL_SECTION_MARK = "-"


def _iso_now() -> str:
    """UTC ISO-8601 timestamp for the record."""
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _relative(path: Path, run_root: Path) -> str:
    """Path of ``path`` relative to the run (never absolute, never ``~``)."""
    try:
        return str(path.resolve().relative_to(run_root.resolve()))
    except ValueError:
        # Outside the run root?  Record a relative tail only, never the
        # absolute path, so no home path can leak into the record.
        parts = path.resolve().parts
        return "/".join(parts[-3:]) if len(parts) >= 3 else path.name


def _read_delfin_data(run_root: Path) -> dict:
    """The run's DELFIN_Data.json as a dict, or {} (never guessed)."""
    p = run_root / _DELFIN_DATA
    if not p.is_file():
        return {}
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (ValueError, OSError):
        return {}


def _structure_origin(run_root: Path, delfin_data: dict) -> dict:
    """Where the structure came from: SMILES / file / previous step, else unknown.

    Oriented on the shape doc_server/calc_indexer reads from DELFIN_Data.json:
    ``input.smiles`` names a SMILES origin.  Otherwise the file CONTROL names
    (``input_file=``) wins; a structure input file in the run only when it is
    not an ORCA/DELFIN-generated artefact.  The relative path is recorded,
    never an absolute one.
    """
    inp = delfin_data.get("input") if isinstance(delfin_data.get("input"), dict) else {}
    smiles = inp.get("smiles")
    if isinstance(smiles, str) and smiles.strip():
        return {"kind": "smiles", "smiles": smiles.strip()}
    named = _named_input(run_root)
    if named and (run_root / named).is_file():
        return {"kind": "file", "path": named}
    for p in sorted(run_root.iterdir()):
        if (
            p.is_file()
            and p.suffix in (".xyz", ".smi", ".sdf", ".mol")
            and p.name not in _GENERATED_OUTPUT_NAMES
            and not p.stem.startswith(_GENERATED_OUTPUT_PREFIXES)
        ):
            return {"kind": "file", "path": _relative(p, run_root)}
    return {"kind": "unknown", "value": UNKNOWN_ORIGIN}


def _program_from_banner(output_text: str) -> dict:
    """Program name + version as printed in the output banner.

    ``{"program": "orca", "version": "6.1.1"}`` when an ORCA banner with a
    version is present; a missing banner/version is recorded as
    ``provenance_missing``, never guessed.
    """
    has_orca_mark = bool(_ORCA_MARK.search(output_text))
    m = _VERSION_RE.search(output_text)
    if m:
        version = ".".join(v for v in m.groups() if v)
        program = "orca" if has_orca_mark else "unknown"
        return {"program": program, "version": version}
    if has_orca_mark:
        return {"program": "orca", "version": MISSING}
    return {"program": MISSING, "version": MISSING}


def _largest_output(run_root: Path) -> Path | None:
    """The largest ``*.out`` in the run folder (its output holds the banner)."""
    outs = [p for p in run_root.glob("*.out") if p.is_file()]
    return max(outs, key=lambda p: p.stat().st_size) if outs else None


def _named_input(run_root: Path) -> str | None:
    """The structure input file the CONTROL names (``input_file=`` key)."""
    control = run_root / "CONTROL.txt"
    if not control.is_file():
        return None
    try:
        for line in control.read_text(encoding="utf-8", errors="replace").splitlines():
            line = line.strip()
            if line.startswith("input_file="):
                name = line.split("=", 1)[1].strip()
                return name or None
    except OSError:
        return None
    return None


def _input_files(run_root: Path) -> list[dict]:
    """Input files of the run, decided by what the CONTROL/workflow reads.

    CONTROL.txt always, plus the file CONTROL names via ``input_file=``; a
    structure input (.xyz/.smi/...) only when it is NOT an ORCA/DELFIN
    generated artefact (``initial*`` outputs and trajectories never count).
    Each entry: ``{"path": <relative>, "sha256": <hex>}``.
    """
    entries: list[dict] = []
    seen: set[str] = set()

    def add(path: Path) -> None:
        rel = _relative(path, run_root)
        if rel in seen or not path.is_file():
            return
        seen.add(rel)
        entries.append({"path": rel, "sha256": cbatch_sha256_file(path)})

    add(run_root / "CONTROL.txt")
    named = _named_input(run_root)
    if named:
        add(run_root / named)
    for p in sorted(run_root.iterdir()):
        if p.name == _DELFIN_DATA or p.name == "CONTROL.txt":
            continue
        if p.is_file() and p.suffix in _INPUT_SUFFIXES:
            stem = p.stem
            if p.name in _GENERATED_OUTPUT_NAMES or (
                any(stem.startswith(pref) for pref in _GENERATED_OUTPUT_PREFIXES)
            ):
                continue
            add(p)
    return entries


def _submit_command(run_root: Path) -> dict:
    """The exact submit command/script when one is recorded, else not_found.

    Searches the run folder for a submit/job script and records its relative
    path plus the first lines (the exact command).  No script -> recorded as
    ``not_found``, never constructed from a guess.
    """
    for p in sorted(run_root.iterdir()):
        if p.suffix in _SUBMIT_SUFFIXES and p.is_file():
            try:
                lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
            except OSError:
                lines = []
            return {
                "script": _relative(p, run_root),
                # The exact command is the whole script; keep the meaningful
                # non-shebang/comment lines.
                "command": "\n".join(
                    ln for ln in lines if ln.strip() and not ln.lstrip().startswith("#")
                ) or MISSING,
            }
    return {"script": None, "command": MISSING}


def _prune_absolute_paths(value):
    """Recursively drop keys whose value is an absolute path (PR-public rule).

    The interpreter probe reports ``executable`` as an absolute path into a
    home directory; a record must never carry such a path, so any string value
    that starts with ``/`` is dropped rather than embedded.
    """
    if isinstance(value, dict):
        return {k: _prune_absolute_paths(v) for k, v in value.items()
                if not (isinstance(v, str) and v.startswith("/"))}
    return value


def compute_step_id(record: dict) -> str:
    """A stable, content-derived id for a record (target 3, chain linkage).

    The id depends only on what the record was *made of*, never on when it was
    built or what else it links to: the sorted input files, the program
    version, the structure origin and the result type.  Rebuilding the same
    run yields the same id, which is what lets a downstream result (a plot)
    name its ancestor calculations.  Signed with a ``pv1:`` prefix so a reader
    can tell provenance-record ids from other uuids.
    """
    stable = {
        "input_files": sorted(
            (f.get("path"), f.get("sha256")) for f in record.get("input_files", [])
        ),
        "program_version": record.get("program_version"),
        "structure_origin": record.get("structure_origin"),
        "result_type": (record.get("result") or {}).get("type"),
    }
    digest = hashlib.sha256(
        json.dumps(stable, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return f"pv1:{digest[:16]}"


def _chain_ids(current_id: str, parents) -> list[str]:
    """The chain of ancestor step ids for ``current_id``, transitive + deduped.

    Each parent is a record dict; its own ``step_id`` plus (transitively) the
    ids it already links are collected in order.  The id of the record being
    built is never part of its own chain, and a shared ancestor (e.g. the same
    structure feeding two steps) appears once.
    """
    seen: list[str] = []

    def add(parent: dict):
        pid = parent.get("step_id")
        if pid and pid != current_id and pid not in seen:
            seen.append(pid)
        for grandparent in parent.get("chain") or []:
            if grandparent != current_id and grandparent not in seen:
                seen.append(grandparent)

    for parent in parents or []:
        if isinstance(parent, dict):
            add(parent)
    return seen


def build_record(
    result_dir: str | Path,
    *,
    parents: list[dict] | None = None,
    result_type: str = "calculation",
) -> dict:
    """Build a complete provenance record dict for ``result_dir`` (a run folder)."""
    run_root = Path(result_dir).resolve()
    delfin_data = _read_delfin_data(run_root)
    out = _largest_output(run_root)
    banner_text = out.read_text(encoding="utf-8", errors="replace") if out else ""
    prog = _program_from_banner(banner_text)
    # DELFIN code state: keep commit + dirty, never the absolute repo path.
    code = cbatch_code_state()
    delfin = {
        "commit": code.get("commit"),
        "dirty": code.get("dirty"),
    }
    # Package versions (reused from cluster_bench) pruned of absolute paths.
    python_packages = _prune_absolute_paths(cbatch_probe_env(sys.executable))
    if not isinstance(python_packages, dict):
        python_packages = {"error": "provenance_unavailable"}

    record = {
        "schema_version": SCHEMA_VERSION,
        "result": {
            "type": result_type,
            "created": _iso_now(),
        },
        "input_files": _input_files(run_root),
        "program": prog["program"],
        "program_version": prog["version"],
        "submit_command": _submit_command(run_root),
        "delfin": delfin,
        "python_packages": python_packages,
        "structure_origin": _structure_origin(run_root, delfin_data),
        "host_class": HOST_CLASS_PLACEHOLDER,
        "timestamps": {
            "created": _iso_now(),
            "run_output_mtime": (
                datetime.fromtimestamp(out.stat().st_mtime, tz=timezone.utc).isoformat(
                    timespec="seconds"
                )
                if out
                else MISSING
            ),
        },
    }
    step_id = compute_step_id(record)
    record["step_id"] = step_id
    record["chain"] = _chain_ids(step_id, parents)
    return record


def write_record(record: dict, result_dir: str | Path) -> Path:
    """Write ``record`` as PROVENANCE.json next to the result; return the path."""
    run_root = Path(result_dir).resolve()
    out_path = run_root / PROVENANCE_FILE
    out_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out_path


def read_record(result_dir: str | Path) -> dict:
    """Read the PROVENANCE.json for a run folder (or the file itself).

    Returns the parsed record; raises FileNotFoundError when absent and
    ValueError when the schema version is unexpected.
    """
    p = Path(result_dir)
    path = p if (p.is_file() and p.name == PROVENANCE_FILE) else p / PROVENANCE_FILE
    if not path.is_file():
        raise FileNotFoundError(f"no {PROVENANCE_FILE} at {path}")
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"unsupported schema_version {record.get('schema_version')!r}, "
            f"expected {SCHEMA_VERSION!r}"
        )
    return record


def _run_verify(record: dict, result_dir: Path) -> int:
    """Recompute every stored checkable value; report each mismatch, non-zero rc.

    Checks: every ``input_files[].sha256`` against the file on disk, and the
    program version against the run's current output banner.  A value recorded
    as ``provenance_missing`` has nothing to contradict and never fails.  Each
    mismatch is printed with the offending path; exit 0 only when none.
    """
    run_root = result_dir if result_dir.is_dir() else result_dir.parent
    mismatches: list[str] = []

    stored = {f.get("path"): f.get("sha256") for f in record.get("input_files", [])}
    for rel, sha in sorted(stored.items()):
        p = run_root / rel
        if not p.is_file():
            mismatches.append(f"{rel}: file missing (stored sha256={sha})")
        else:
            current = cbatch_sha256_file(p)
            if current != sha:
                mismatches.append(
                    f"{rel}: sha256 mismatch (stored {sha}, now {current})"
                )

    stored_version = record.get("program_version")
    if stored_version != MISSING:
        out = _largest_output(run_root)
        text = out.read_text(encoding="utf-8", errors="replace") if out else ""
        current_banner = _program_from_banner(text)
        if current_banner["version"] != stored_version:
            mismatches.append(
                f"program_version: banner now {current_banner['version']!r}, "
                f"record holds {stored_version!r}"
            )

    for line in mismatches:
        print(f"verify MISMATCH: {line}")
    if mismatches:
        return 1
    print("verify OK: all stored sha256 and program version match")
    return 0


def format_for_display(record: dict) -> str:
    """Human-readable, one-line-per-field rendering for ``delfin provenance``."""
    lines = [
        f"schema_version: {record.get('schema_version')}",
        f"program: {record.get('program')} {record.get('program_version')}",
        f"delfin commit: {record.get('delfin', {}).get('commit')} "
        f"dirty={record.get('delfin', {}).get('dirty')}",
        f"structure_origin: {record.get('structure_origin')}",
        f"host_class: {record.get('host_class')}",
    ]
    for f in record.get("input_files", []):
        lines.append(f"input  {f.get('path')} sha256={f.get('sha256')}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """CLI entry for ``delfin provenance <result>``; prints the record.

    ``<result>`` is a run folder holding PROVENANCE.json (or the file itself).
    """
    import argparse

    parser = argparse.ArgumentParser(
        prog="delfin provenance",
        description="Print the provenance record of one result run.",
    )
    parser.add_argument("result", help="run folder (or PROVENANCE.json path)")
    parser.add_argument("--json", action="store_true", help="print raw JSON")
    parser.add_argument(
        "--verify",
        action="store_true",
        help="recompute every stored sha256 and the banner version; report "
        "each mismatch and exit non-zero when one is found",
    )
    args = parser.parse_args(argv)
    try:
        record = read_record(args.result)
    except (FileNotFoundError, ValueError) as exc:
        print(f"provenance unavailable: {exc}", file=sys.stderr)
        return 2
    if args.verify:
        return _run_verify(record, Path(args.result))
    if args.json:
        print(json.dumps(record, indent=2, sort_keys=True))
    else:
        print(format_for_display(record))
    return 0
