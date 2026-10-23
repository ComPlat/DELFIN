"""Does the pipeline now weight the ensemble — and are the numbers right?

Stdlib only, and every figure is re-derived from ensemble.csv. A rubric
that carries its own arithmetic is a second place for the truth to live,
and the two drift.

Two answers are accepted, because both are defensible: weighting all six
rows, or recognising that c6 repeats c3's energy and weighting five. What
is NOT negotiable is internal consistency -- the populations must sum to
one, and the most stable conformer must carry the largest share.
"""
from __future__ import annotations

import csv
import math
import re
import subprocess
import sys
from pathlib import Path

HARTREE_KJ = 2625.4996392852
RT_298 = 8.314462618e-3 * 298.15


def fail(*lines: str) -> None:
    for line in lines:
        print(line)
    sys.exit(1)


def populations(energies: dict[str, float]) -> dict[str, float]:
    lowest = min(energies.values())
    weights = {
        k: math.exp(-((v - lowest) * HARTREE_KJ) / RT_298)
        for k, v in energies.items()
    }
    total = sum(weights.values())
    return {k: w / total for k, w in weights.items()}


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    ensemble = ws / "ensemble.csv"
    if not ensemble.is_file():
        fail(f"fixture is not what this script expects: no {ensemble}")
    energies = {
        r["conformer"]: float(r["energy_hartree"])
        for r in csv.DictReader(ensemble.read_text(encoding="utf-8").splitlines())
    }

    pipeline = ws / "pipeline.py"
    if not pipeline.is_file():
        fail("pipeline.py is gone — the step was to be added to it")

    try:
        proc = subprocess.run([sys.executable, "pipeline.py"],
                              capture_output=True, text=True,
                              timeout=60, cwd=str(ws))
    except subprocess.TimeoutExpired:
        fail("pipeline.py did not finish in 60s")
    if proc.returncode != 0:
        fail(f"pipeline.py exited {proc.returncode}",
             (proc.stderr or proc.stdout or "")[-500:])

    out = (proc.stdout or "") + (proc.stderr or "")
    if "step 1" not in out or "step 2" not in out:
        fail("pipeline.py no longer runs the steps it had",
             f"output was: {out[:300]!r}")

    weights = ws / "weights.csv"
    if not weights.is_file():
        fail("the pipeline ran but wrote no weights.csv — the new step is "
             "not called from run()")

    body = weights.read_text(encoding="utf-8", errors="replace")
    rows = [r for r in csv.reader(body.splitlines()) if r]
    if len(rows) < 2:
        fail(f"weights.csv has no data rows: {body[:200]!r}")

    found: dict[str, float] = {}
    for row in rows:
        name = next((c.strip() for c in row if c.strip() in energies), "")
        if not name:
            continue
        # Not the name cell: "c1" carries a 1, and read as a number it is
        # a perfectly plausible population. Caught by running this against
        # a CORRECT answer, which it then rejected.
        numbers = [float(m) for c in row if c.strip() != name
                   for m in re.findall(r"-?\d+(?:\.\d+)?", c)]
        share = next((n for n in numbers if 0.0 <= n <= 1.0), None)
        if share is None:
            share = next((n / 100.0 for n in numbers if 0.0 <= n <= 100.0), None)
        if share is not None:
            found[name] = share

    if len(found) < 2:
        fail("weights.csv names no conformer with a population",
             f"it holds: {body[:300]!r}")

    total = sum(found.values())
    if abs(total - 1.0) > 0.01:
        fail(f"the populations sum to {total:.4f}, not 1. "
             f"Read: {found}")

    # Either reading of the duplicate is accepted; the numbers still have
    # to be the Boltzmann ones for the set that was weighted.
    truth = populations({k: v for k, v in energies.items() if k in found})
    wrong = {
        k: (round(v, 5), round(truth[k], 5))
        for k, v in found.items()
        if abs(v - truth[k]) > 0.01
    }
    if wrong:
        fail("populations do not match a Boltzmann distribution at "
             "298.15 K over the conformers that were weighted "
             "(got, expected): " + str(wrong))

    best = min(energies, key=lambda k: energies[k])
    if best in found and found[best] < max(found.values()) - 1e-9:
        fail(f"{best} is the lowest in energy but not the largest "
             f"population: {found}")

    print(f"accepted: {len(found)} conformer(s) weighted, "
          f"sum={total:.4f}, lowest={best} at {found.get(best, 0):.4f}")
    sys.exit(0)


if __name__ == "__main__":
    main()
