#!/usr/bin/env python3
"""Fixture for the fundamentals waiting task: a result that arrives late.

samples.csv holds six enthalpy readings (kJ/mol). slow_report.py sleeps
20 s — a FIXED constant, never a machine-speed probe; the rubric reads
no wall clock — and then writes result.json with the mean, 42.19
kJ/mol. The mean is re-checked against the arithmetic before the
fixture is accepted, so the value the rubric grades is measured, not
invented.

Exit non-zero on any failure: a task whose precondition was not built
is reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import sys
from pathlib import Path

_SAMPLES = "index,enthalpy_kj_per_mol\n1,42.31\n2,42.05\n3,42.18\n" \
           "4,42.37\n5,41.98\n6,42.25\n"
_SLEEP_S = 20.0

_SCRIPT = '''"""Compute the mean enthalpy of samples.csv — slowly, on purpose.

Sleeps 20 s, then writes result.json. Run in the background and collect
the result when it lands.
"""
import csv
import json
import time

TIME_S = {sleep:.1f}


def main():
    time.sleep(TIME_S)
    with open("samples.csv", newline="", encoding="utf-8") as fh:
        vals = [float(r["enthalpy_kj_per_mol"]) for r in csv.DictReader(fh)]
    mean = sum(vals) / len(vals)
    with open("result.json", "w", encoding="utf-8") as fh:
        json.dump({{"n": len(vals), "mean_kj_per_mol": round(mean, 2)}}, fh)
    print(f"{{len(vals)}} samples, mean {{mean:.2f}} kJ/mol")


if __name__ == "__main__":
    main()
'''


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: fund_waiting.py <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    d = ws / "fund_waiting"
    if d.exists():
        print(f"refusing to overwrite {d}", file=sys.stderr)
        return 1
    rows = [ln.split(",") for ln in _SAMPLES.strip().splitlines()[1:]]
    mean = sum(float(r[1]) for r in rows) / len(rows)
    if abs(mean - 42.19) > 5e-3:
        print(f"fixture mean {mean:.4f} kJ/mol, expected 42.19",
              file=sys.stderr)
        return 1
    d.mkdir(parents=True)
    (d / "samples.csv").write_text(_SAMPLES, encoding="utf-8")
    (d / "slow_report.py").write_text(_SCRIPT.format(sleep=_SLEEP_S),
                                      encoding="utf-8")
    print(f"prepared {d} — samples.csv (mean {mean:.2f} kJ/mol), "
          f"slow_report.py (sleeps {_SLEEP_S:.0f} s, then result.json)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
