#!/usr/bin/env python3
"""Build a small, reproducible calculation archive for the calc tools.

"Which of my runs used which method" is daily work for the scientist this
agent is for, and it was completely unmeasured -- because `search_calcs`
and `calc_summary` build their index from the user's REAL calc/ and
archive/ folders. On this machine that is 777 calculations whose contents
nobody chose, so a task over them cannot be scored: the right answer
differs per developer and changes whenever a real run finishes.

This builds nine calculations with contents that ARE chosen, so a task
can ask a question with one true answer:

    calc/     four active runs
    archive/  five finished ones

The distribution is deliberate rather than uniform, so a question has a
non-trivial answer and a lazy one is wrong:

  * As INDEXED: PBE0 5 times, B3LYP 3. "The most used functional" has an
    answer, and it is not "all of them".
  * calc_d's TPSSh is deliberately NOT in that count. The index reads the
    functional from DELFIN_Data.json, which an unfinished run does not
    have -- so a run still in flight contributes no method. That is the
    indexer's behaviour, not a flaw in the fixture, and a task must not
    ask a question whose answer depends on it.
  * def2-TZVP appears 4 times, def2-SVP 5 -- close enough that guessing
    from the functional does not work.
  * ONE run (arch_e) has no CONTROL.txt and only an .inp, so the index
    has to fall back to parsing the input; a task that counts runs will
    miss it if the fallback is broken.
  * ONE run (calc_d) is unfinished: no DELFIN_Data.json, no .out. It is
    in `calc/` and must not be reported as a completed result.

The guard around a benchmark attempt removes this with the rest of the
workspace, and the runner points the calc tools at it -- so nothing here
touches the user's own archive.

Exit non-zero on any failure: a task whose precondition was not built is
reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


# (folder, functional, basis, solvent, workflow, with_data, with_out)
ACTIVE = [
    ("calc_a", "PBE0", "def2-TZVP", "DMF", "classic", True, True),
    ("calc_b", "B3LYP", "def2-SVP", "water", "classic", True, True),
    ("calc_c", "PBE0", "def2-SVP", "none", "OCCUPIER", True, True),
    ("calc_d", "TPSSh", "def2-TZVP", "acetonitrile", "classic", False, False),
]
ARCHIVED = [
    ("arch_a", "PBE0", "def2-TZVP", "DMF", "classic", True, True),
    ("arch_b", "PBE0", "def2-SVP", "DMF", "classic", True, True),
    ("arch_c", "B3LYP", "def2-SVP", "toluene", "classic", True, True),
    ("arch_d", "B3LYP", "def2-TZVP", "water", "OCCUPIER", True, True),
    ("arch_e", "PBE0", "def2-SVP", "DMF", "classic", False, True),
]

CONTROL = """# DELFIN CONTROL
NAME = {name}
SMILES = c1ccccc1
charge = 0
multiplicity = 1
method = {workflow}
functional = {functional}
basis_set = {basis}
solvent = {solvent}
implicit_solvation_model = CPCM
IMAG = no
ESD_modul = no
"""

INP = """! {functional} {basis} TightSCF
%pal nprocs 8 end
* xyz 0 1
C   0.000  0.000  0.000
O   0.000  0.000  1.128
*
"""

OUT = """
                  * O   R   C   A *

Program Version 6.0.1

FINAL SINGLE POINT ENERGY      {energy:.8f}

                             ****ORCA TERMINATED NORMALLY****
"""


def build(base: Path, rows, kind: str) -> None:
    base.mkdir(parents=True, exist_ok=True)
    for i, (name, func, basis, solvent, workflow, data, out) in enumerate(rows):
        d = base / name
        d.mkdir()
        (d / f"{name}.inp").write_text(
            INP.format(functional=func, basis=basis), encoding="utf-8")
        # arch_e deliberately has no CONTROL.txt: the index must fall back
        # to the .inp, and a task that counts runs misses it if it cannot.
        if name != "arch_e":
            (d / "CONTROL.txt").write_text(
                CONTROL.format(name=name, workflow=workflow, functional=func,
                               basis=basis, solvent=solvent),
                encoding="utf-8")
        if data:
            (d / "DELFIN_Data.json").write_text(json.dumps({
                "name": name, "functional": func, "basis_set": basis,
                "solvent": solvent, "charge": 0, "multiplicity": 1,
                "status": "finished" if kind == "archive" else "running",
            }, indent=1), encoding="utf-8")
        if out:
            (d / f"{name}.out").write_text(
                OUT.format(energy=-113.30 - i * 0.01), encoding="utf-8")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    root = ws / "calc_archive"
    if root.exists():
        print(f"refusing to overwrite {root}", file=sys.stderr)
        return 1
    build(root / "calc", ACTIVE, "calc")
    build(root / "archive", ARCHIVED, "archive")

    # Check the precondition rather than assume it: if the indexer cannot
    # read this tree, the task measures nothing.
    try:
        from delfin.doc_server.calc_indexer import build_calc_index
        idx = build_calc_index(calc_dir=root / "calc",
                               archive_dir=root / "archive", quiet=True)
    except Exception as exc:                                  # pragma: no cover
        print(f"the indexer could not read the fixture: {exc}", file=sys.stderr)
        return 1
    n = len(idx.get("calculations") or idx.get("records") or [])
    if n < 9:
        print(f"indexer found {n} calculations, expected 9", file=sys.stderr)
        return 1
    print(f"prepared {root} — {n} calculations, 4 active and 5 archived")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
