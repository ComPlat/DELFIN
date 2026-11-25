"""Run a given list of test node ids on THIS tree and report which still fail.

The full-suite comparison answers "did the change add failures" only if both
sides ran the same tests. This asks the narrower and far cheaper question that
is enough here: every node id that failed on the branch -- does it also fail on
this tree? Anything that fails on the branch and passes here is a regression of
the change; anything listed here as failing is pre-existing.

Usage (from the tree under test):
    <python> scratch/x1/recheck.py <failed.json> <out.json>
"""

import json
import pathlib
import subprocess
import sys
import time

failed_path = pathlib.Path(sys.argv[1])
out_path = pathlib.Path(sys.argv[2])
ids = json.loads(failed_path.read_text())

ROOT = pathlib.Path.cwd()
CMD = [
    sys.executable,
    "-m",
    "pytest",
    "-q",
    "-p",
    "no:randomly",
    "-p",
    "no:cacheprovider",
    "--timeout=300",
    "--timeout-method=signal",
    "--tb=no",
    "-rfE",
]

# Group by file so pytest is invoked once per file, not once per test.
by_file = {}
for nid in ids:
    path, _, rest = nid.partition("::")
    by_file.setdefault(path, []).append(nid)

still, cleared, rced = [], [], {}
t0 = time.time()
for i, (path, group) in enumerate(sorted(by_file.items()), 1):
    if not (ROOT / path).exists():
        rced[path] = "MISSING-FILE"
        continue
    p = subprocess.run(CMD + [path], cwd=ROOT, capture_output=True, text=True, errors="replace")
    rced[path] = p.returncode
    reported = {
        ln.split(" ", 1)[1].split(" - ")[0]
        for ln in p.stdout.splitlines()
        if ln.startswith("FAILED ") or ln.startswith("ERROR ")
    }
    for nid in group:
        (still if nid in reported else cleared).append(nid)
    print(f"[{i}/{len(by_file)}] {path}: rc={p.returncode} of {len(group)} still failing "
          f"{sum(1 for n in group if n in reported)}", flush=True)

out_path.write_text(json.dumps({"still": still, "cleared": cleared, "rc": rced}, indent=1))
print(f"\nwall {round(time.time() - t0, 1)} s")
print(f"still failing on this tree: {len(still)}")
print(f"passing on this tree (i.e. branch-only): {len(cleared)}")
for nid in cleared:
    print("  CLEARED-ON-THIS-TREE:", nid)
