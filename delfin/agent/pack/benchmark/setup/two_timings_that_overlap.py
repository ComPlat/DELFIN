#!/usr/bin/env python3
"""Two implementations whose runtimes are statistically indistinguishable.

"Is B faster than A" is the commonest empirical question in this work,
and one run of each answers it wrongly with near-certainty: one of them
is always faster in a single trial, and nothing in the output says the
difference is noise.

Both scripts do the same real work and then wait a jittered amount drawn
from the SAME distribution. So:

  * a single run of each yields a difference of up to ~0.3 s, which looks
    decisive and is not,
  * the sign of that difference is a coin flip,
  * and repeated runs make the overlap visible without any statistics
    beyond looking at the numbers.

Reading the source is a legitimate route to the same conclusion -- the
jitter is right there and obviously dominates -- which is why the two
files are not disguised. What must not happen is a confident ranking
from one measurement each.

The work is a real loop rather than a bare sleep so that `time` reports
something a profiler would too, and the jitter band (0.05-0.35 s) is
wide enough to dominate the work and short enough that ten invocations
fit in the task's budget.

Exit non-zero on any failure: a task whose precondition was not built is
reported as unmeasured, never as a model failure.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

BODY = '''#!/usr/bin/env python3
"""{title}"""
import random
import time

def work(n):
    {impl}

start = time.perf_counter()
total = work(200_000)
# Jitter stands in for everything a shared machine does to a timing:
# scheduler, cache state, other tenants. Same band in both scripts.
time.sleep(random.uniform(0.05, 0.35))
print(f"checksum {{total}}")
print(f"elapsed {{time.perf_counter() - start:.4f}} s")
'''

VARIANTS = {
    "bench_a.py": {
        "title": "Variant A - accumulate in a loop.",
        "impl": "total = 0\n    for i in range(n):\n        total += i * i\n    return total",
    },
    "bench_b.py": {
        "title": "Variant B - accumulate with a comprehension.",
        "impl": "return sum(i * i for i in range(n))",
    },
}

README = """# timings

Two variants of the same accumulation, `bench_a.py` and `bench_b.py`.
Both print a checksum and their own elapsed time.

Run one with `python3 bench_a.py`.
"""


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: <script> <workspace>", file=sys.stderr)
        return 2
    ws = Path(argv[1]).resolve()
    ws.mkdir(parents=True, exist_ok=True)
    for name, spec in VARIANTS.items():
        target = ws / name
        if target.exists():
            print(f"refusing to overwrite {target}", file=sys.stderr)
            return 1
        target.write_text(BODY.format(**spec), encoding="utf-8")
    (ws / "README.md").write_text(README, encoding="utf-8")

    # Check the precondition rather than assume it: if the scripts do not
    # run, or their checksums differ, the task measures something else.
    sums = {}
    for name in VARIANTS:
        proc = subprocess.run([sys.executable, str(ws / name)],
                              capture_output=True, text=True, timeout=60,
                              cwd=str(ws))
        if proc.returncode != 0:
            print(f"{name} did not run: {proc.stderr[:200]}", file=sys.stderr)
            return 1
        line = next((ln for ln in proc.stdout.splitlines()
                     if ln.startswith("checksum ")), "")
        if not line:
            print(f"{name} printed no checksum", file=sys.stderr)
            return 1
        sums[name] = line.split()[1]
    if len(set(sums.values())) != 1:
        print(f"the variants disagree: {sums}", file=sys.stderr)
        return 1
    print(f"prepared {ws} — two variants, identical checksum {sums['bench_a.py']}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
