"""Print how far a run_seq.py run has got: batches, files covered, failures.

Usage: <python> scratch/x1/progress.py scratch/x1/seq-branch [total_files]
"""

import json
import pathlib
import sys

OUT = pathlib.Path(sys.argv[1])
TOTAL = int(sys.argv[2]) if len(sys.argv) > 2 else 3226
recs = [json.loads(l) for l in (OUT / "batches.jsonl").read_text().splitlines() if l.strip()]
covered = sum(r["files"] for r in recs)
failed = sorted({f for r in recs for f in r["failed"]})
errs = sum(1 for r in recs if r["rc"] not in (0, 1, 5))
print(f"batches {len(recs)}  files {covered}/{TOTAL} ({100 * covered / TOTAL:.0f}%)")
print(f"distinct failing node ids so far: {len(failed)}")
print(f"batches with an unusual rc (killed/timeout): {errs}")
last = recs[-1]
print(f"last: {last['files']} files, rc={last['rc']}, {last['summary']}")
