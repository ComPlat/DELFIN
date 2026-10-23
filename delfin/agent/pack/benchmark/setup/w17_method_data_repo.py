#!/usr/bin/env python3
"""Shared builder for the two w17_method data tasks.

Not a setup on its own: the w17_method_data_*_setup.py entry scripts
call ``build`` through generalist_repo.main, which supplies the setup
contract the suite expects (``usage:`` on a missing workspace argument,
``refusing to overwrite`` on a second run, ``prepared <path>`` on
success) and the reviewed fixture-builder git() the safe_git ALLOWED
list covers.

variant="audit": data/runs.csv holds 6 converged runs; two share the
method (functional, basis, dispersion, solvent) AND the molecule, so
one pair is a true duplicate and one pair differs in dispersion and is
not. Task: write data/AUDIT.txt naming the duplicate pair.

variant="convert": data/values.csv with German decimal commas
("1.234,5" means 1234.5). Task: write data/converted.json mapping
every id to its float value.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from generalist_repo import (  # noqa: E402
    Precondition, commit, new_repo, write,
)

# A fixture repo must not depend on the user's home git configuration:
# sandboxed test environments deny ~/.gitconfig (fatal "Permission
# denied" on init). The shared generalist_repo core runs plain git;
# these redirects make it home-independent. Set before any git call
# runs (generalist_repo reads os.environ per call).
import os  # noqa: E402

os.environ.setdefault("GIT_CONFIG_GLOBAL", os.devnull)
os.environ.setdefault("GIT_CONFIG_SYSTEM", os.devnull)
# The XDG excludes file (~/.config/git/ignore) is a third home path
# git reads even with the config redirects above; cage denial is
# fatal ("cannot use ... as an exclude file").
os.environ.setdefault("XDG_CONFIG_HOME", os.devnull)

RUNS_CSV = """id,molecule,functional,basis,dispersion,solvent
r1,H2O,B3LYP,def2-SVP,D3,none
r2,H2O,B3LYP,def2-SVP,D3,none
r3,H2O,PBE0,def2-SVP,D3,none
r4,CO2,B3LYP,def2-TZVP,D3,none
r5,CO2,B3LYP,def2-SVP,D3,none
r6,CO2,B3LYP,def2-TZVP,D3BJ,none
"""

VALUES_CSV = """id,value
a,1.234,5
b,-0,5
c,10.000,25
d,2,7
e,1000,1
"""

EXPECTED_CONVERT = {"a": 1234.5, "b": -0.5, "c": 10000.25, "d": 2.7,
                    "e": 1000.1}

_BRIEF_AUDIT = ("Goal: find duplicate converged runs in data/runs.csv. "
                "Two runs are duplicates when molecule, functional, basis, "
                "dispersion and solvent all match. Success criterion: "
                "data/AUDIT.txt names the duplicate pair by run_id "
                "(one line, format: 'duplicate: rX rY').\n")
_BRIEF_CONVERT = ("Goal: convert data/values.csv to JSON. The values use "
                  "the German decimal comma ('1.234,5' means 1234.5). "
                  "Success criterion: data/converted.json holds every id "
                  "from values.csv mapped to its float value.\n")


def build(workspace: str, *, variant: str) -> Path:
    if variant not in ("audit", "convert"):
        raise Precondition(f"unknown variant {variant!r}")
    repo = new_repo(workspace, "data_notes")
    (repo / "data").mkdir()
    if variant == "audit":
        write(repo, {"data/runs.csv": RUNS_CSV, "m1_brief.txt": _BRIEF_AUDIT})
    else:
        write(repo, {"data/values.csv": VALUES_CSV,
                     "m1_brief.txt": _BRIEF_CONVERT})
    commit(repo, "seed")
    # Preconditions, asserted before the task can be booked as measured.
    if variant == "audit":
        rows = [line.split(",") for line in RUNS_CSV.strip().splitlines()[1:]]
        keys = {r[0]: tuple(r[1:6]) for r in rows}
        dup_keys = [(ka, kb) for i, (ka, va) in enumerate(keys.items())
                    for (kb, vb) in list(keys.items())[i + 1:]
                    if va == vb]
        if sorted(dup_keys[0]) != ["r1", "r2"] or len(dup_keys) != 1:
            raise Precondition("audit precondition broken: expected "
                               "exactly the r1/r2 duplicate")
    else:
        parsed = {k: float(",".join(v).replace(".", "").replace(",", "."))
                  for k, *v in (line.split(",") for line
                                in VALUES_CSV.strip().splitlines()[1:])}
        if parsed != EXPECTED_CONVERT:
            raise Precondition("convert precondition broken")
    print(f"seeded data_notes/{variant}: preconditions hold")
    return repo


if __name__ == "__main__":
    # Entry scripts import this module; direct execution is a mistake.
    print("not a setup on its own; use w17_method_data_*_setup.py",
          file=sys.stderr)
    sys.exit(2)
