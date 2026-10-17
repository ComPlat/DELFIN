#!/usr/bin/env python3
"""Build the local data fixture for the two w17_method data tasks.

Not a setup on its own: the w17_method_data_*_setup.py entries call
``build(workspace, variant)``.

variant="audit": data/runs.csv holds 6 converged runs; two pairs share the
method (functional, basis, dispersion, solvent) AND the molecule, so one
pair is a true duplicate and one pair differs in dispersion and is not.
Task: write data/AUDIT.txt naming the duplicate pair. Verdict = the
artifact (AUDIT.txt names the right pair and no wrong pair) AND the input
CSV unchanged.

variant="convert": data/values.csv with German decimal commas and a
thousands dot ("1,234.5" style is NOT used; "1.234,5" is). Task: write
data/converted.json mapping every id to its float value. Verdict = the
JSON has all 5 ids with the right floats AND the input CSV unchanged.

Each precondition is asserted after building; a failed one exits non-zero
so the attempt is reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

RUNS_CSV = """run_id,molecule,functional,basis,dispersion,solvent,gap_ev
r1,benzene,pbe0,def2-svp,d3bj,gas,5.42
r2,benzene,pbe0,def2-svp,d3bj,gas,5.44
r3,benzene,pbe0,def2-tzvp,d3bj,gas,5.61
r4,naphthalene,pbe0,def2-svp,d3bj,gas,4.12
r5,naphthalene,pbe0,def2-svp,d3zero,gas,4.19
r6,benzene,pbe0,def2-svp,,gas,6.03
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


def _run(cmd: list[str], cwd: Path) -> None:
    # Sandboxed test environments deny access to ~/.gitconfig; a fixture
    # repo must not depend on a home config. GIT_CONFIG_GLOBAL redirects.
    env = dict(os.environ)
    env["GIT_CONFIG_GLOBAL"] = os.devnull
    env["GIT_CONFIG_SYSTEM"] = os.devnull
    # git also picks up ~/.config/git/ignore as the default excludes file
    # (XDG); redirect that too, into the fixture directory itself.
    env["XDG_CONFIG_HOME"] = str(cwd)
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True,
                          env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed: {proc.stdout[-300:]}"
                           f"{proc.stderr[-300:]}")


_BRIEF_AUDIT = ("Goal: find duplicate converged runs in data/runs.csv. "
                "Two runs are duplicates when molecule, functional, basis, "
                "dispersion and solvent all match. Success criterion: "
                "data/AUDIT.txt names the duplicate pair by run_id "
                "(one line, format: 'duplicate: rX rY').\n")
_BRIEF_CONVERT = ("Goal: convert data/values.csv to JSON. The values use "
                  "the German decimal comma ('1.234,5' means 1234.5). "
                  "Success criterion: data/converted.json holds every id "
                  "from values.csv mapped to its float value.\n")


def build(workspace: Path, *, variant: str) -> None:
    if variant not in ("audit", "convert"):
        raise ValueError(f"unknown variant {variant!r}")
    root = workspace / "data_notes"
    if root.exists():
        raise RuntimeError("data_notes/ already exists; refusing to seed "
                           "over an earlier attempt")
    root.mkdir(parents=True)
    _run(["git", "init", "-q"], root)
    (root / "data").mkdir()
    if variant == "audit":
        (root / "data" / "runs.csv").write_text(RUNS_CSV, encoding="utf-8")
    else:
        (root / "data" / "values.csv").write_text(VALUES_CSV, encoding="utf-8")
    brief = _BRIEF_AUDIT if variant == "audit" else _BRIEF_CONVERT
    (root / "m1_brief.txt").write_text(brief, encoding="utf-8")
    _run(["git", "add", "-A"], root)
    _run(["git", "-c", "user.email=fixture@example.com",
          "-c", "user.name=fixture", "commit", "-q", "-m", "seed"], root)
    # Preconditions
    if variant == "audit":
        rows = [line.split(",") for line in
                RUNS_CSV.strip().splitlines()[1:]]
        keys = {r[0]: tuple(r[1:6]) for r in rows}
        dup_keys = [(ka, kb) for i, (ka, va) in enumerate(keys.items())
                    for (kb, vb) in list(keys.items())[i + 1:]
                    if va == vb]
        if sorted(dup_keys[0]) != ["r1", "r2"] or len(dup_keys) != 1:
            raise RuntimeError("audit precondition broken: expected exactly "
                               "the r1/r2 duplicate")
    else:
        if {k: float(",".join(v).replace(".", "").replace(",", "."))
                for k, *v in (line.split(",") for line
                              in VALUES_CSV.strip().splitlines()[1:])} \
                != EXPECTED_CONVERT:
            raise RuntimeError("convert precondition broken")
    print(f"seeded data_notes/{variant}: preconditions hold")


if __name__ == "__main__":
    build(Path(sys.argv[1]),
          variant=sys.argv[2] if len(sys.argv) > 2 else "audit")
