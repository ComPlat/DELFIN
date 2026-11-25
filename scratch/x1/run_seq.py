"""Run all of tests/ sequentially in time-boxed batches, appending each batch's
outcome, so a run that the machine kills still leaves the batches it finished.

Concurrency is not used on purpose: this sandbox cannot start the threads a
pool needs, and the supervisor process itself aborted with
`RuntimeError: can't start new thread` under load from other work on the
machine. A sequential run is slower but its numbers survive.

Resumable: batches already recorded in <outdir>/batches.jsonl are skipped, so
re-invoking continues where a killed run stopped.

Usage (from the tree under test):
    <python> scratch/x1/run_seq.py <outdir> [batch_seconds]
"""

import json
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path.cwd()
OUT = pathlib.Path(sys.argv[1]).resolve()
BATCH_S = int(sys.argv[2]) if len(sys.argv) > 2 else 500
OUT.mkdir(parents=True, exist_ok=True)
DONE = OUT / "batches.jsonl"

FILES = sorted(
    str(p.relative_to(ROOT))
    for p in (ROOT / "tests").rglob("test_*.py")
    if "__pycache__" not in p.parts
)

recorded = set()
if DONE.exists():
    for line in DONE.read_text().splitlines():
        try:
            recorded.add(json.loads(line)["first"])
        except Exception:
            pass

CMD = [
    sys.executable,
    "-m",
    "pytest",
    "-q",
    "-m",
    "not slow",
    "--timeout=300",
    "--timeout-method=signal",
    "-p",
    "no:randomly",
    "-p",
    "no:cacheprovider",
    "--tb=no",
    "-rfE",
]

print(f"files: {len(FILES)}  already recorded: {len(recorded)}", flush=True)

i = 0
while i < len(FILES):
    if FILES[i] in recorded:
        i += 1
        continue
    batch = [FILES[i]]
    # Grow the batch until the previous batch's rate says we would run over.
    if not (OUT / "rate.txt").exists():
        batch = FILES[i : i + 40]
    else:
        per_s = float((OUT / "rate.txt").read_text())
        n = max(1, min(200, int(per_s * BATCH_S)))
        batch = FILES[i : i + n]
    t0 = time.time()
    proc = subprocess.run(CMD + batch, cwd=ROOT, capture_output=True, text=True, errors="replace")
    secs = time.time() - t0
    failed = [
        ln.split(" ", 1)[1].split(" - ")[0]
        for ln in proc.stdout.splitlines()
        if ln.startswith("FAILED ") or ln.startswith("ERROR ")
    ]
    summary = ""
    for ln in reversed(proc.stdout.splitlines()):
        if " passed" in ln or " failed" in ln or " error" in ln:
            summary = ln.strip()
            break
    rec = {
        "first": batch[0],
        "files": len(batch),
        "rc": proc.returncode,
        "secs": round(secs, 1),
        "summary": summary,
        "failed": failed,
    }
    with DONE.open("a") as fh:
        fh.write(json.dumps(rec) + "\n")
    (OUT / "rate.txt").write_text(str(len(batch) / max(secs, 1e-9)))
    print(f"[{i:5d}/{len(FILES)}] {len(batch):3d} files {secs:6.1f}s rc={proc.returncode} {summary}", flush=True)
    for f in failed:
        print(f"    {f}", flush=True)
    i += len(batch)
    if secs > BATCH_S * 2:
        break

print("done-or-stopped", flush=True)
