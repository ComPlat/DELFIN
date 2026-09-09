#!/usr/bin/env python3
"""Build a git repository whose test was red BEFORE the last commit.

The state the task needs, and one no fixture can carry: git does not
track a nested repository, so this has to be built at run time. The guard
around the attempt snapshots the workspace before this runs and restores
it afterwards, so nothing here survives the task.

Two commits, and the whole question is which one to blame:

  HEAD~1  "add the Boltzmann weighting"   <- energies.py, already wrong
  HEAD    "add a CSV export"              <- export.py, innocent

``test_weights.py`` fails at BOTH. An agent that reads only the diff of
the last commit sees a file it does not understand next to a red test and
attributes the failure to it; the control run at HEAD~1 settles it in one
command.

The bug is a real one rather than a `1/0`: the Boltzmann weights are
divided by the number of conformers instead of by their sum, so they do
not add to 1. That matters because the wrong answer is PLAUSIBLE -- an
agent can talk itself into either commit having caused it.

Exit non-zero on any failure: a task whose precondition was not built is
reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ENERGIES = '''"""Boltzmann weights over a conformer ensemble."""
import math

KT_KCAL = 0.593          # RT at 298.15 K, kcal/mol


def weights(energies_kcal):
    """Relative populations for a list of relative energies."""
    boltz = [math.exp(-e / KT_KCAL) for e in energies_kcal]
    return [b / len(boltz) for b in boltz]
'''

TEST = '''from energies import weights


def test_weights_sum_to_one():
    w = weights([0.0, 1.2, 2.4])
    assert abs(sum(w) - 1.0) < 1e-9
'''

EXPORT = '''"""Write a weight table as CSV."""
import csv


def write_table(path, labels, values):
    with open(path, "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["conformer", "weight"])
        for label, value in zip(labels, values):
            out.writerow([label, f"{value:.4f}"])
    return path
'''


def git(ws: Path, *args: str) -> None:
    proc = subprocess.run(["git", *args], cwd=str(ws),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(args)} failed: {proc.stderr.strip()[:200]}")


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    repo = ws / "ensemble_tools"
    if repo.exists():
        print(f"refusing to overwrite {repo}", file=sys.stderr)
        return 1
    repo.mkdir(parents=True)

    git(repo, "init", "-q")
    # A local identity: the machine's may be unset, and a commit that
    # cannot be made is a task that cannot be measured.
    git(repo, "config", "user.email", "fixture@delfin.invalid")
    git(repo, "config", "user.name", "fixture")
    git(repo, "config", "commit.gpgsign", "false")

    (repo / "energies.py").write_text(ENERGIES, encoding="utf-8")
    (repo / "test_weights.py").write_text(TEST, encoding="utf-8")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "add the Boltzmann weighting")

    (repo / "export.py").write_text(EXPORT, encoding="utf-8")
    git(repo, "add", "-A")
    git(repo, "commit", "-q", "-m", "add a CSV export")

    # The precondition itself is checked, not assumed: if the test does
    # not actually fail, the task measures nothing.
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "test_weights.py"],
        cwd=str(repo), capture_output=True, text=True, timeout=60)
    if proc.returncode == 0:
        print("the fixture test PASSES; the task has no failure to attribute",
              file=sys.stderr)
        return 1
    print(f"prepared {repo} — two commits, test red at both")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
