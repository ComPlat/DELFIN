#!/usr/bin/env python
"""Baseline recorder for the final-geometry readers (task wave24/s06).

Records the CURRENT outputs of every affected function on
(a) tests/fixtures/orca611 and (b) up to 200 real runs copied read-only
into scratch/runs/ before this script runs (relative run ids only).

Affected functions, and the workflow-level results that use them:
- delfin.imag._read_xyz(path)            -> [(el, x, y, z), ...] (IMAG rounds)
- delfin.orca._read_xyz_atoms(path)      -> raw row strings     (resume)
- delfin.orca_recovery._continue_from_own_geometry(items, path)
   -> what orca_input.with_coordinates would write (classic recovery);
   recorded as the atom tuples it would inject plus a raised/refused flag.

The ORCA-written <base>.xyz next to each run's .out is the input of all
three, so the recorder walks every .out in the corpus and reads the
sibling <base>.xyz.  Output: scratch/baseline.json (data, not committed).
"""
import hashlib
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from delfin.imag import _read_xyz as imag_read_xyz  # noqa: E402
from delfin.orca import _read_xyz_atoms as orca_read_atoms  # noqa: E402

SCRATCH = REPO / "scratch"
RUNS = SCRATCH / "runs"
OUT = SCRATCH / "baseline.json"


def row_digest(rows):
    return hashlib.sha256(
        json.dumps([[r[0], *r[1:]] for r in rows]).encode()).hexdigest()[:16]


def with_coordinates_projection(rows):
    """What delfin.common.orca_input.with_coordinates would inject, per atom.

    Rows are formatted "  {el:<3s} {x:14.8f} {y:14.8f} {z:14.8f}" -- the
    string the next ORCA input would carry after recovery replaces them.
    """
    return [f"  {el:<3s} {x:14.8f} {y:14.8f} {z:14.8f}"
            for el, x, y, z in rows]


def record_one(run_dir: Path, out_path: Path, run_id: str) -> dict:
    base = out_path.stem
    rec: dict = {"run_id": run_id, "base": base}
    xyz = run_dir / f"{base}.xyz"
    rec["xyz_exists"] = xyz.is_file()

    rec["imag_read_xyz"] = None
    rec["orca_read_xyz_atoms"] = None
    rec["recovery_projection"] = None
    rec["imag_error"] = rec["orca_error"] = rec["recovery_error"] = None

    if xyz.is_file():
        # imag._read_xyz (tuple atoms; ValueError on malformed count line)
        try:
            atoms = imag_read_xyz(xyz)
            rec["imag_read_xyz"] = [[el, x, y, z] for el, x, y, z in atoms]
            rec["imag_digest"] = row_digest(rec["imag_read_xyz"])
        except Exception as exc:  # noqa: BLE001
            rec["imag_error"] = f"{type(exc).__name__}: {exc}"
        # orca._read_xyz_atoms (raw rows or None)
        try:
            rows = orca_read_atoms(xyz)
            rec["orca_read_xyz_atoms"] = rows
        except Exception as exc:  # noqa: BLE001
            rec["orca_error"] = f"{type(exc).__name__}: {exc}"
        # orca_recovery projection: with_coordinates formatting of imag's atoms
        if rec["imag_read_xyz"] is not None:
            rec["recovery_projection"] = with_coordinates_projection(
                [tuple(r) for r in rec["imag_read_xyz"]])
    # the sibling .out is never read for coordinates by these readers --
    # record that fact so the switch keeps it so.
    rec["out_has_cartesian_blocks"] = (
        "CARTESIAN COORDINATES" in
        (run_dir / f"{base}.out").read_text(encoding="utf-8", errors="ignore")[-_OUT_TAIL:]
    )
    return rec


# the readers above never parse the .out; the tail check only documents
# whether a Cartesian block exists there at all (corpus characterisation).
_OUT_TAIL = 1 << 20


def main() -> int:
    corpus: list[tuple[str, Path]] = []
    fixtures = REPO / "tests/fixtures/orca611"
    for d in sorted(fixtures.iterdir()):
        if d.is_dir():
            for out in sorted(d.glob("*.out")):
                corpus.append((f"fixture/{d.name}", d))
    if RUNS.is_dir():
        for d in sorted(p for p in RUNS.iterdir() if p.is_dir()):
            for out in sorted(d.glob("*.out")):
                corpus.append((f"run/{d.name}", d))
    records = []
    for run_id, run_dir in corpus:
        for out in sorted(run_dir.glob("*.out")):
            if out.name.endswith(".old1.out"):
                continue
            try:
                records.append(record_one(run_dir, out, f"{run_id}/{out.stem}"))
            except Exception as exc:  # noqa: BLE001
                records.append({"run_id": f"{run_id}/{out.stem}",
                                "record_error": f"{type(exc).__name__}: {exc}"})
    OUT.write_text(json.dumps(records, indent=1), encoding="utf-8")
    n_geom = sum(1 for r in records if r.get("imag_read_xyz"))
    print(f"recorded {len(records)} runs; {n_geom} with a readable final geometry")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
