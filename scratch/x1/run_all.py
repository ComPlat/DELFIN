"""Run the whole of tests/ as N concurrent pytest processes and report the
outcome per chunk.

Why this exists: the single-process run of `pytest tests/ -m "not slow"` on
this shared machine did not finish inside a job timeout, so no number existed
for the full suite. Splitting it keeps each process's work bounded.

Concurrency is by plain subprocess.Popen polling, not a thread pool: this
sandbox cannot start enough threads for a pool this wide (measured:
`RuntimeError: can't start new thread` at 24 workers).

Any test that fails here is re-run on its own, sequentially, before its failure
is believed -- concurrency can invent failures, and the sequential re-run is
what distinguishes a real one. Pass --recheck to do that.

Usage (from the tree under test):
    <python> scratch/x1/run_all.py <outdir> [workers] [--recheck]
"""

import json
import os
import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path.cwd()
OUT = pathlib.Path(sys.argv[1]).resolve()
args = [a for a in sys.argv[2:] if not a.startswith("--")]
WORKERS = int(args[0]) if args else 16
RECHECK = "--recheck" in sys.argv
OUT.mkdir(parents=True, exist_ok=True)

FILES = sorted(
    str(p.relative_to(ROOT))
    for p in (ROOT / "tests").rglob("test_*.py")
    if "__pycache__" not in p.parts
)
print(f"test files: {len(FILES)}  workers: {WORKERS}", flush=True)

# Balance by file size, which tracks test count closely enough and avoids one
# worker getting every heavyweight module.
FILES.sort(key=lambda f: os.path.getsize(ROOT / f), reverse=True)
chunks = [[] for _ in range(WORKERS)]
loads = [0] * WORKERS
for f in FILES:
    i = loads.index(min(loads))
    chunks[i].append(f)
    loads[i] += os.path.getsize(ROOT / f)

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
]


def start(idx):
    files = chunks[idx]
    if not files:
        return None
    log = (OUT / f"chunk-{idx:02d}.log").open("w")
    proc = subprocess.Popen(
        CMD + files,
        cwd=ROOT,
        stdout=log,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return {"chunk": idx, "files": len(files), "proc": proc, "log": log, "t0": time.time()}


def failed_ids(idx):
    out = []
    log = OUT / f"chunk-{idx:02d}.log"
    if not log.exists():
        return out
    for line in log.read_text(errors="replace").splitlines():
        if line.startswith("FAILED ") or line.startswith("ERROR "):
            out.append(line.split(" ", 1)[1].split(" - ")[0])
    return out


t0 = time.time()
running = []
pending = list(range(WORKERS))
done = []
while pending or running:
    while pending and len(running) < WORKERS:
        s = start(pending.pop(0))
        if s:
            running.append(s)
        elif not pending:
            break
    time.sleep(2)
    still = []
    for s in running:
        if s["proc"].poll() is None:
            still.append(s)
            continue
        s["log"].close()
        s["rc"] = s["proc"].returncode
        s["secs"] = round(time.time() - s["t0"], 1)
        s["failed"] = failed_ids(s["chunk"])
        done.append(s)
    running = still

print(f"wall: {round(time.time() - t0, 1)} s", flush=True)
tallies = {}
for s in done:
    log = OUT / f"chunk-{s['chunk']:02d}.log"
    last = ""
    for line in reversed(log.read_text(errors="replace").splitlines()):
        if " passed" in line or " failed" in line or " error" in line:
            last = line.strip()
            break
    s["out"] = last
    tallies[s["chunk"]] = last
    if s["rc"] != 0:
        print(f"  chunk {s['chunk']:02d} rc={s['rc']} files={s['files']} {last}", flush=True)

allfailed = sorted({f for s in done for f in s["failed"]})
(OUT / "failed.json").write_text(json.dumps(allfailed, indent=1))
(OUT / "tallies.json").write_text(json.dumps(tallies, indent=1))
print(f"\ndistinct failing node ids: {len(allfailed)}", flush=True)
for f in allfailed:
    print(" ", f, flush=True)

if RECHECK and allfailed:
    print("\n-- sequential re-check of each failing file --", flush=True)
    by_file = sorted({f.split("::", 1)[0] for f in allfailed if "::" in f})
    for path in by_file:
        p = subprocess.run(
            CMD + ["--tb=line", path], cwd=ROOT, capture_output=True, text=True
        )
        ids = [
            ln.split(" ", 1)[1].split(" - ")[0]
            for ln in p.stdout.splitlines()
            if ln.startswith("FAILED ") or ln.startswith("ERROR ")
        ]
        (OUT / f"recheck-{path.replace('/', '_')}.log").write_text(p.stdout + p.stderr)
        print(f"  {path}: rc={p.returncode} still-failing={len(ids)}", flush=True)
        for i in ids:
            print(f"    {i}", flush=True)
