"""The Slurm array script of a run: one shard per array task, one full node per task.

Defaults fit a 48-core node and a 72 h wall-time limit (JUSTUS 2, for example): ``--cpus-per-task=48``,
``--time=72:00:00``, array throttle ``%N``.  The per-system limit is NOT the job's wall time: it is
ceil(timeout_base x speed_factor) inside the runner.  A task that reaches the job's wall time
leaves its finished systems in place; resubmitting the same index continues it.

Partitions: a script that names none is submitted through ``delfin.slurm_submit.sbatch_command``,
which asks ``sbatch --test-only`` which configured partitions accept it.
"""
from __future__ import annotations

import shlex
import subprocess
import sys
from pathlib import Path

from delfin.cluster_bench.prepare import cbatch_load_manifest
from delfin.cluster_bench.provenance import cbatch_effective_timeout
from delfin.slurm_submit import normalize_time_limit, sbatch_command


def cbatch_render_sbatch(run_dir, man, *, set_name="main", run_name="main", array=None,
                         throttle=40, time_limit="72:00:00", cpus=48, mem=None, workers=None,
                         speed_factor=None, partition=None, account=None, python=None,
                         setup=(), job_name=None) -> str:
    run_dir = Path(run_dir).resolve()
    s = man["settings"]
    n = man["sets"][set_name]["n_shards"]
    array = array or f"0-{n - 1}"
    if throttle:
        array = f"{array}%{int(throttle)}"
    sf = str(speed_factor) if speed_factor is not None else s["speed_factor"]
    tmo = cbatch_effective_timeout(s["timeout_base_s"], sf)
    python = python or sys.executable
    workers = int(workers or s["workers"])
    mem = mem or {"manta": "72G"}.get(man["tool"], "80G")
    job_name = job_name or f"dc_{man['tool']}_{set_name}_{run_name}"
    sb = [f"--job-name={job_name}", "--nodes=1", "--ntasks=1", f"--cpus-per-task={int(cpus)}",
          f"--mem={mem}", f"--time={normalize_time_limit(time_limit)}", f"--array={array}",
          f"--output={run_dir}/logs/%x_%A_%a.out", f"--error={run_dir}/logs/%x_%A_%a.err",
          "--signal=B:USR1@600"]
    if partition:
        sb.append(f"--partition={partition}")
    if account:
        sb.append(f"--account={account}")
    cmd = [python, "-m", "delfin.cluster_bench", "run-shard", str(run_dir),
           "--shard", "$SLURM_ARRAY_TASK_ID", "--set", set_name, "--run", run_name,
           "--workers", str(workers), "--speed-factor", sf]
    cmd_s = " ".join('"$SLURM_ARRAY_TASK_ID"' if c == "$SLURM_ARRAY_TASK_ID" else shlex.quote(c)
                     for c in cmd)
    lines = ["#!/bin/bash",
             f"# delfin cluster: {man['tool']} ({man['label']}), set {set_name}, run {run_name}",
             f"# {n} shards; per-system limit {tmo} s = ceil({s['timeout_base_s']} x {sf});",
             f"# {workers} concurrent builds x {s['threads']} thread(s) per node; resumable.",
             *[f"#SBATCH {x}" for x in sb],
             "",
             "set -euo pipefail",
             "export PYTHONNOUSERSITE=1 PYTHONHASHSEED=0",
             "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1",
             "unset PYTHONPATH PYTHONSTARTUP || true",
             '# numba / matplotlib caches on node-local disk, never in $HOME',
             'export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/numba_${SLURM_JOB_ID:-0}" '
             'MPLCONFIGDIR="${TMPDIR:-/tmp}/mpl_${SLURM_JOB_ID:-0}"',
             *list(setup),
             'echo "[$(date -Is)] task ${SLURM_ARRAY_TASK_ID} on $(hostname), '
             '${SLURM_CPUS_PER_TASK:-?} cpus"',
             "trap 'echo \"[$(date -Is)] USR1: wall time near -- finished systems are kept, "
             "resubmit this index\"' USR1",
             f"{cmd_s} &",
             "PID=$!",
             "# wait returns early when USR1 arrives; keep waiting while the runner lives",
             "while true; do",
             "    wait $PID && RC=0 || RC=$?",
             "    kill -0 $PID 2>/dev/null || break",
             "done",
             "exit $RC",
             ""]
    return "\n".join(lines)


def cbatch_write_sbatch(run_dir, *, submit=False, sbatch="sbatch", **kw) -> tuple:
    """-> (script path, sbatch command list, submit output or None)."""
    run_dir = Path(run_dir).resolve()
    man = cbatch_load_manifest(run_dir)
    set_name, run_name = kw.get("set_name", "main"), kw.get("run_name", "main")
    if set_name not in man["sets"]:
        raise SystemExit(f"shard set {set_name!r} not in this run ({sorted(man['sets'])})")
    text = cbatch_render_sbatch(run_dir, man, **kw)
    path = run_dir / "slurm" / f"{man['tool']}_{set_name}_{run_name}.sbatch"
    path.write_text(text)
    path.chmod(0o755)
    (run_dir / "logs").mkdir(exist_ok=True)
    cmd = sbatch_command(sbatch, path, run_dir) if submit else [sbatch, str(path)]
    out = None
    if submit:
        done = subprocess.run(cmd, capture_output=True, text=True, cwd=str(run_dir))
        out = (done.stdout + done.stderr).strip()
        if done.returncode != 0:
            raise SystemExit(f"sbatch failed ({done.returncode}): {out}")
    return path, cmd, out
