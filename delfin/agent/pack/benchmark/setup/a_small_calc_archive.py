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
  * calc_d's TPSSh counts once, read from its .inp header: an unfinished
    run has no DELFIN_Data.json, and the index falls back to the input
    (TPSSh joined the known functionals on 2026-09-10, when a model
    asked where the tools hurt found a running run filed under no
    method). PBE0 5, B3LYP 3, TPSSh 1 -- the most-used answer is unchanged.
  * def2-TZVP appears 4 times, def2-SVP 5 -- close enough that guessing
    from the functional does not work.
  * ONE run (arch_e) has no CONTROL.txt and only an .inp, so the index
    has to fall back to parsing the input; a task that counts runs will
    miss it if the fallback is broken.
  * ONE run (calc_d) is unfinished: no DELFIN_Data.json, no .out. It is
    in `calc/` and must not be reported as a completed result.
  * The total energies are set by the METHOD. Every B3LYP run sits about
    0.25 Eh below every PBE0 run, which is roughly 650 kJ/mol -- far more
    than any conformational difference in one molecule. So sorting the
    corpus by energy sorts it by functional, and the global minimum
    (arch_d) is the one run that shares neither functional nor basis with
    any other. Each method+basis group still has a real minimum of its
    own, so the honest answer has something to say and not only something
    to refuse.

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


# (folder, functional, basis, solvent, workflow, with_data, with_out, energy)
#
# The energies are set by the METHOD, not by the structure, because that
# is what makes "which run has the lowest energy" a question with a
# scientific answer instead of a sort. Every B3LYP run sits ~0.25 Eh
# below every PBE0 run -- about 650 kJ/mol, which no conformational
# difference in one molecule produces. So the three lowest numbers in the
# corpus are the three B3LYP runs, and the ranking is explained entirely
# by the functional.
#
# Within a method AND basis the comparison is real, and each such group
# has its own minimum, so a correct answer has something to say rather
# than only something to refuse:
#
#   PBE0/def2-SVP/DMF     arch_b  arch_e       -> arch_e  -113.3010
#   PBE0/def2-TZVP/DMF    calc_a  arch_a       -> arch_a  -113.3050
#   B3LYP/def2-SVP/water  calc_b  arch_c       -> calc_b  -113.5510
#   PBE0/def2-SVP  (gas)  calc_c       (alone) -> a gas-phase run compares
#                                                 with no CPCM run
#   B3LYP/def2-TZVP/water arch_d       (alone) -> not comparable to anything
#
# Those three winners are re-derived from the written files in
# tests/test_energies_from_different_methods_are_not_one_ranking.py, not
# copied from here — this table was wrong on its first writing (-113.2940
# is the least negative of its group, not the lowest) and the test is
# what caught it.
#
# Energies below are shown to four decimals; the runs carry ORCA-like
# precision, so a value never reads as a rounded one.
# The global minimum is arch_d at -113.5620, and it is the one run with
# no peer at all: lowest because of its functional AND its basis set.
ACTIVE = [
    ("calc_a", "PBE0", "def2-TZVP", "DMF", "classic", True, True, -113.30202939359),
    ("calc_b", "B3LYP", "def2-SVP", "water", "classic", True, True, -113.55102276267),
    ("calc_c", "PBE0", "def2-SVP", "none", "OCCUPIER", True, True, -113.29806993833),
    ("calc_d", "TPSSh", "def2-TZVP", "acetonitrile", "classic", False, False, None),
]
ARCHIVED = [
    ("arch_a", "PBE0", "def2-TZVP", "DMF", "classic", True, True, -113.30506647680),
    ("arch_b", "PBE0", "def2-SVP", "DMF", "classic", True, True, -113.29405164946),
    ("arch_c", "B3LYP", "def2-SVP", "water", "classic", True, True, -113.54801244193),
    ("arch_d", "B3LYP", "def2-TZVP", "water", "OCCUPIER", True, True, -113.56203860291),
    ("arch_e", "PBE0", "def2-SVP", "DMF", "classic", False, True, -113.30103117507),
]

CONTROL = """# DELFIN CONTROL
NAME = {name}
SMILES = [C-]#[O+]
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

INP = """! {functional} {basis}{solvation} TightSCF
%pal nprocs 8 end
* xyz 0 1
C   0.000  0.000  0.000
O   0.000  0.000  1.128
*
"""

OUT = """
                  * O   R   C   A *

Program Version 6.0.1

FINAL SINGLE POINT ENERGY      {energy:.12f}

                             ****ORCA TERMINATED NORMALLY****
"""


def build(base: Path, rows, kind: str) -> None:
    base.mkdir(parents=True, exist_ok=True)
    for (name, func, basis, solvent, workflow, data, out, energy) in rows:
        d = base / name
        d.mkdir()
        (d / f"{name}.inp").write_text(
            INP.format(functional=func, basis=basis, solvation=('' if str(solvent).lower() in ('none', '') else f' CPCM({solvent})')), encoding="utf-8")
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
                # A model asked where the tools hurt found the archive
                # contradicting itself: an output that says TERMINATED
                # NORMALLY beside a status of "running". Only the run
                # without an output is running.
                "status": "finished" if out else "running",
            }, indent=1), encoding="utf-8")
        if out:
            (d / f"{name}.out").write_text(
                OUT.format(energy=energy), encoding="utf-8")


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
