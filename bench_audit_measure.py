"""Ad-hoc measurement harness: run the science audit task N times.

Writes the run JSONL next to itself (the default ~/.delfin/benchmark_runs
is read-only in this worktree session). Prints per-sample and aggregate
quality/success so the before/after comparison carries the spread.
"""
import sys
from pathlib import Path

from delfin.agent import benchmark as bm
from delfin.agent import benchmark_runner as br

MODEL = sys.argv[1] if len(sys.argv) > 1 else "kit.glm-5.3"
# The task that measures the assignment "compare methods, not programs".
# science_three_runs_are_audited_independently audits run STATES (finished /
# failed / running) and cannot see method comparability at all — measuring
# against it produced numbers that say nothing about the prompt rule.
TASK = sys.argv[2] if len(sys.argv) > 2 else "science_energies_from_different_methods_are_not_ranked"
REPEATS = int(sys.argv[3]) if len(sys.argv) > 3 else 3

tasks = [t for t in bm.load_tasks() if t.id == TASK]
assert tasks, f"task not found: {TASK}"

samples = []


def on_rep(task, idx, res):
    samples.append(res)
    print(f"  sample {idx+1}: ok={res.success} q={res.quality:.0f} "
          f"tool={res.tool_calls} ${res.cost_usd:.4f}")


results = br.run_suite(
    tasks, model=MODEL, backend="api", provider="kit",
    repeats=REPEATS, on_replicate=on_rep)

out = Path(__file__).resolve().parent / "bench_runs"
out.mkdir(exist_ok=True)
path = bm.write_run(results, model=MODEL, runs_dir=out)
print(f"\nrun written: {path}")
for r in results:
    print(f"task={r.task_id} ok={r.success} q_med={r.quality_0_100:.0f} "
          f"sigma={getattr(r, 'quality_stdev', None)} "
          f"rate={getattr(r, 'success_rate', None)}")
