#!/usr/bin/env python3
"""Run the recorder over all fixtures + the scratch sample, shelling out
so PYTHONPATH resolves delfin to this worktree."""
import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
recorder = root / "scratch" / "record_esd_reads.py"
out = root / "scratch" / Path(sys.argv[1] if len(sys.argv) > 1
                              else "baseline.json").name

run_dirs = sorted(
    p for p in (root / "scratch" / "sample").glob("*/*/ESD")
    if p.is_dir()
) + sorted(p for p in (root / "tests" / "fixtures" / "orca611").iterdir()
           if p.is_dir())
print(f"{len(run_dirs)} run dirs")
cmd = ["/home/qmchem_max/delfin-swarm/venv/bin/python", str(recorder),
       *[str(p) for p in run_dirs], str(out)]
r = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
print(r.stdout[-3000:])
if r.returncode != 0:
    print(r.stderr[-4000:])
sys.exit(r.returncode)
