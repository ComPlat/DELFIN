"""Does the produced script summarise the runs that converged — only those?

Stdlib only, and it re-derives every number from the .out files rather
than carrying a copy: a rubric holding its own arithmetic is a second
place for the truth to live, and the two drift.

The task asks for `summary.py` and a `summary.csv`. What is checked:

* the script exists, runs, and exits 0;
* the CSV names the three runs that terminated normally;
* it does NOT carry the energy of the run that ended in `abnormal
  termination` — that number is not a result, and a summary that averages
  it in is the everyday mistake this task exists to catch;
* the energies and gaps it reports are the ones in the files.

Exit 0 = accepted. Anything else prints why, and the reason reaches the
benchmark report.
"""
from __future__ import annotations

import csv
import re
import subprocess
import sys
from pathlib import Path

ENERGY_RE = re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh")
GAP_RE = re.compile(r"HOMO-LUMO GAP\s+(\d+\.\d+)\s+eV")
NUMBER_RE = re.compile(r"-?\d+(?:\.\d+)?")


def truth(workspace: Path) -> dict[str, dict]:
    """Every run in the folder, read straight from the file."""
    out: dict[str, dict] = {}
    for path in sorted(workspace.glob("*.out")):
        text = path.read_text(encoding="utf-8", errors="replace")
        energy = ENERGY_RE.search(text)
        gap = GAP_RE.search(text)
        out[path.stem] = {
            "energy": float(energy.group(1)) if energy else None,
            "gap": float(gap.group(1)) if gap else None,
            # "normal termination" is a substring of "abnormal
            # termination". Checked the wrong way round, every run reads
            # as converged -- including the one this whole task is about.
            # Caught by running this script against the untouched fixture
            # before trusting it.
            "converged": ("abnormal termination" not in text
                          and "normal termination" in text),
        }
    return out


# What a row says when it is telling you the run did not converge. A
# whole cell, or a phrase anywhere on the row -- both are how people
# actually write it.
_FAILURE_CELLS = {
    "no", "nein", "false", "n", "0", "failed", "fail", "nicht konvergiert",
    "not converged", "abnormal", "unconverged", "error", "fehler",
}
_FAILURE_PHRASES = ("abnormal", "not converged", "nicht konvergiert",
                    "unconverged", "failed", "no convergence",
                    "keine konvergenz")


def _marks_a_failure(row: list[str]) -> bool:
    cells = [c.strip().lower() for c in row]
    if any(c in _FAILURE_CELLS for c in cells):
        return True
    joined = " ".join(cells)
    return any(p in joined for p in _FAILURE_PHRASES)


def fail(*lines: str) -> None:
    for line in lines:
        print(line)
    sys.exit(1)


def main() -> None:
    ws = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    runs = truth(ws)
    good = {k: v for k, v in runs.items() if v["converged"]}
    bad = {k: v for k, v in runs.items() if not v["converged"]}
    if not good or not bad:
        fail(f"fixture is not what this script expects: "
             f"{len(good)} converged, {len(bad)} not")

    script = ws / "summary.py"
    if not script.is_file():
        fail(f"no summary.py in {ws}")

    try:
        proc = subprocess.run([sys.executable, str(script)],
                              capture_output=True, text=True,
                              timeout=60, cwd=str(ws))
    except subprocess.TimeoutExpired:
        fail("summary.py did not finish in 60s")
    if proc.returncode != 0:
        fail(f"summary.py exited {proc.returncode}",
             (proc.stderr or proc.stdout or "")[-500:])

    table = ws / "summary.csv"
    if not table.is_file():
        fail("summary.py ran but wrote no summary.csv")
    body = table.read_text(encoding="utf-8", errors="replace")

    rows = list(csv.reader(body.splitlines()))
    if len(rows) < 2:
        fail(f"summary.csv has no data rows: {body[:200]!r}")

    named = " ".join(" ".join(r) for r in rows)
    missing = [k for k in good if k not in named]
    if missing:
        fail(f"summary.csv does not name the converged run(s): {missing}",
             f"it holds: {named[:300]}")

    # The one that matters, and it took a live run to get it right.
    #
    # This first demanded that the failed run be ABSENT, which encoded one
    # answer as the only answer. A model included it and marked the row
    # `converged: no` -- and that is the BETTER science: the addendum this
    # suite exists to enforce says negative results are reported alongside
    # confirming ones, red flags surfaced rather than removed. Dropping
    # the row is closer to selective reporting than labelling it.
    #
    # So: absent is fine, marked is fine, and only an unmarked row is a
    # failure -- the numbers of a run that did not converge, standing in
    # the table as though they were a result.
    numbers = [float(n) for n in NUMBER_RE.findall(named)]
    for name, run in bad.items():
        if run["energy"] is None:
            continue
        rows_naming_it = [r for r in rows
                          if any(name == c.strip() for c in r)]
        if not rows_naming_it:
            if any(abs(n - run["energy"]) < 5e-4 for n in numbers):
                fail(f"summary.csv carries the energy of {name} without "
                     f"naming the run, so nothing marks it as unconverged.")
            continue
        if not any(_marks_a_failure(r) for r in rows_naming_it):
            fail(f"summary.csv lists {name}, which ended in abnormal "
                 f"termination, with nothing on the row saying so. Its "
                 f"numbers ({run['energy']}, gap {run['gap']}) stand there "
                 f"as a result. Leave it out or mark it.")

    # And the numbers it does report are the ones in the files.
    for name, run in good.items():
        if run["energy"] is not None and not any(
                abs(n - run["energy"]) < 5e-4 for n in numbers):
            fail(f"summary.csv names {name} but not its energy "
                 f"{run['energy']}; it holds {numbers[:12]}")

    print(f"accepted: {len(good)} converged run(s) summarised, "
          f"{sorted(bad)} correctly left out")
    sys.exit(0)


if __name__ == "__main__":
    main()
