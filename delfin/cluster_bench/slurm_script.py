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


def _cbatch_staging_lines(tool_python) -> list:
    """Bash that stages the tool's venv once per job onto node-local disk and exports the two
    env vars the runner reads (``DELFIN_CLUSTER_TOOL_PYTHON``, ``DELFIN_CLUSTER_LOCAL_ROOT``).
    Reuses the same content-keyed venv cache as submit_delfin.sh's ``ensure_venv_tar`` -- one
    sequential read of a cached tar, unpacked onto ``$TMPDIR`` -- and falls back to running
    children from the shared file system with one warning when ``$TMPDIR`` is missing or too
    small.  Empty when the run's manifest names no tool interpreter (nothing to stage)."""
    if not tool_python:
        return []
    return [
        "# ------------------------------------------------------------- node-local I/O",
        "# The tool interpreter is a venv on the shared file system; starting 48 x thousands of",
        "# children from it re-imports torch/MACE/RDKit for every system (the millions of opens",
        "# this run was cut for).  Unpack it once per job onto $TMPDIR and run every child from",
        "# the staged copy; per-system work, logs and the build archive go under that local root",
        "# too, copied back to the workspace in batches.",
        f"TOOL_PYTHON={tool_python}",
        "export DELFIN_CLUSTER_LOCAL_ROOT=\"\"",
        "export DELFIN_CLUSTER_TOOL_PYTHON=\"\"",
        "cb_venv_cache_key() {   # same content-keyed cache as submit_delfin.sh/ensure_venv_tar",
        "    local venv=\"$1\" site",
        "    {",
        "        readlink -f \"$venv\" 2>/dev/null || printf '%s\\n' \"$venv\"",
        "        cat \"$venv/pyvenv.cfg\" 2>/dev/null || true",
        "        for site in \"$venv\"/lib/python3*/site-packages; do",
        "            [ -d \"$site\" ] || continue",
        "            printf '== %s\\n' \"${site#$venv/}\"",
        "            { ls -1 -f \"$site\" 2>/dev/null | grep -v -x -e . -e .. -e __pycache__ | LC_ALL=C sort; } || true",
        "        done",
        "    } | sha256sum | cut -c1-16",
        "}",
        "cb_stage_tool_env() {",
        "    local venv_python=\"$1\" base=\"$2\" venv key tar_path tar_k free_k",
        "    venv=\"$(dirname \"$(dirname \"$venv_python\")\")\"   # <venv>/bin/python -> venv root",
        "    cache_dir=\"${DELFIN_VENV_CACHE_DIR:-$HOME/.cache/delfin/venv}\"",
        "    key=\"$(cb_venv_cache_key \"$venv\")\"",
        "    tar_path=\"$cache_dir/venv-$key.tar\"",
        "    mkdir -p \"$cache_dir\" || return 1",
        "    if [ ! -f \"$tar_path\" ]; then",
        "        if command -v flock >/dev/null 2>&1; then exec 9>\"$cache_dir/.lock\"; flock 9 || true; fi",
        "        if [ ! -f \"$tar_path\" ]; then",
        "            local partial=\"$tar_path.partial.${SLURM_JOB_ID:-$$}\"",
        "            tar -cf \"$partial\" --warning=no-file-changed -C \"$(dirname \"$venv\")\" \"$(basename \"$venv\")\" 2>/dev/null || true",
        "            if [ -s \"$partial\" ] && tar -tf \"$partial\" >/dev/null 2>&1; then mv -f \"$partial\" \"$tar_path\"; else rm -f \"$partial\"; fi",
        "        fi",
        "        flock -u 9 2>/dev/null || true; exec 9>&- 2>/dev/null || true",
        "    fi",
        "    [ -f \"$tar_path\" ] || return 1",
        "    tar_k=$(( $(stat -c %s \"$tar_path\" 2>/dev/null || echo 0) / 1024 ))",
        "    free_k=\"$(df -Pk \"$base\" 2>/dev/null | awk 'NR==2 {print $4}')\"",
        "    [ -n \"$free_k\" ] && [ \"$free_k\" -lt $(( tar_k * 2 )) ] && return 1",
        "    mkdir -p \"$base/venv\" || return 1",
        "    if ! tar -xf \"$tar_path\" --strip-components=1 -C \"$base/venv\" 2>/dev/null; then return 1; fi",
        "    # children run as '<python> <script>', so only the activate script must know the staged",
        "    # path -- never rewrite bin/* (GNU sed -i would replace the bin/python symlink with an",
        "    # edited copy of the target binary).",
        "    sed -i \"s|$venv|$base/venv|g\" \"$base/venv/pyvenv.cfg\" \"$base/venv/bin/activate\" 2>/dev/null || true",
        "    [ -x \"$base/venv/bin/python\" ] || return 1",
        "    return 0",
        "}",
        "if [ -n \"${TMPDIR:-}\" ] && [ -d \"${TMPDIR:-}\" ]; then",
        "    DELFIN_CLUSTER_LOCAL_ROOT=\"${TMPDIR}/cb_${SLURM_JOB_ID:-$$}\"",
        "    mkdir -p \"$DELFIN_CLUSTER_LOCAL_ROOT\"",
        "    DELFIN_CLUSTER_TOOL_PYTHON=\"$DELFIN_CLUSTER_LOCAL_ROOT/venv/bin/python\"",
        "    if [ -n \"$TOOL_PYTHON\" ] && [ -x \"$TOOL_PYTHON\" ] && cb_stage_tool_env \"$TOOL_PYTHON\" \"$DELFIN_CLUSTER_LOCAL_ROOT\"" \
            " && [ -x \"$DELFIN_CLUSTER_TOOL_PYTHON\" ]; then",
        "        :   # staged interpreter present and executable",
        "    else",
        "        echo \"WARNING: could not stage a usable tool interpreter onto $TMPDIR; keeping the tool environment and per-system work/log on the shared file system.\"",
        "        DELFIN_CLUSTER_LOCAL_ROOT=\"\"",
        "        DELFIN_CLUSTER_TOOL_PYTHON=\"\"",
        "    fi",
        "else",
        "    echo \"WARNING: no \\$TMPDIR on this node; keeping the tool environment and per-system work/log on the shared file system.\"",
        "fi",
        "export DELFIN_CLUSTER_LOCAL_ROOT DELFIN_CLUSTER_TOOL_PYTHON",
    ]


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
             *(_cbatch_staging_lines(s.get("tool_python"))),
             'echo "[$(date -Is)] task ${SLURM_ARRAY_TASK_ID} on $(hostname), '
             '${SLURM_CPUS_PER_TASK:-?} cpus"',
             "trap 'echo \"[$(date -Is)] USR1: wall time near -- flushing finished systems, "
             "resubmit this index\"; [ -n \"$DELFIN_CLUSTER_LOCAL_ROOT\" ] && touch "
             "\"$DELFIN_CLUSTER_LOCAL_ROOT/.usr1\"' USR1",
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
