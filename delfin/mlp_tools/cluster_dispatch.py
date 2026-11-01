"""Run MLP jobs where they belong: GPU on the cluster, CPU here.

The cluster path renders the sbatch script for an MLP single point,
optimisation or batch ranking and hands it to Slurm; the local path is a
capped CPU run (<= 2 threads, <= 30 minutes, small molecules).  The choice
is driven entirely by configuration — ``runtime.slurm`` in settings and the
``DELFIN_MLP_CLUSTER_*`` environment — never by hard-coded host names or
paths.  All new behaviour is opt-in: with no configuration, ``decide_dispatch``
answers ``local``, exactly what DELFIN did before this module existed.

Submission discipline:
* ``write_and_submit`` defaults to ``submit=False`` — a dry run that writes
  the script and runs nothing.  A real submission is an explicit act.
* A batch of N structures is ONE sbatch job (the ensemble file is passed
  whole), never N jobs.
* At most ``MAX_CONCURRENT_MLP_JOBS`` = 4 sbatch MLP jobs may be in flight;
  ``assert_can_submit`` refuses the fifth with ``JobCapExceeded``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from delfin.common.logging import get_logger

logger = get_logger(__name__)

#: At most this many sbatch MLP jobs may run at once.
MAX_CONCURRENT_MLP_JOBS = 4

#: Local CPU fallback limits (machine rule of this host: no GPU, small jobs).
LOCAL_THREADS = 2
LOCAL_MINUTES = 30

VALID_TASKS = ("optimize", "single_point", "rank")


class JobCapExceeded(RuntimeError):
    """More than MAX_CONCURRENT_MLP_JOBS sbatch MLP jobs were requested."""


def slurm_available() -> bool:
    """True when this host can talk to a Slurm controller."""
    return bool(shutil.which("sbatch")) and bool(shutil.which("squeue"))


@dataclass
class MLPClusterConfig:
    """Everything a GPU MLP job needs, all of it from config.

    Sources (later wins): ``runtime.slurm`` / ``runtime.mlp`` in settings,
    then the ``DELFIN_MLP_CLUSTER_*`` environment variables.  Empty by
    default: without configuration there is no cluster path.
    """

    gpu_partitions: str = ""        # comma-separated, e.g. "gpu_a100,gpu"
    partition: str = ""             # CPU partition, used when gpu= no
    gres: str = "gpu:1"
    time_limit: str = "04:00:00"
    mem_mb: int = 16000
    cpus_per_task: int = 2
    env_activation: str = ""        # shell lines activating the MLP environment
    max_concurrent: int = MAX_CONCURRENT_MLP_JOBS
    submit_script: str = ""         # override for the sbatch body (rare)

    @classmethod
    def from_sources(cls, settings=None, env=None) -> "MLPClusterConfig":
        settings = settings if settings is not None else {}
        env = env if env is not None else dict(os.environ)
        slurm = (settings.get("slurm") or {}) if isinstance(settings, dict) else {}
        mlp = (settings.get("mlp") or {}) if isinstance(settings, dict) else {}
        cfg = cls(
            gpu_partitions=str(slurm.get("gpu_partitions", "") or ""),
            partition=str(slurm.get("partition", "") or ""),
            gres=str(mlp.get("gres", "") or cls.gres),
            time_limit=str(mlp.get("time_limit", "") or cls.time_limit),
            mem_mb=int(mlp.get("mem_mb", cls.mem_mb) or cls.mem_mb),
            cpus_per_task=int(mlp.get("cpus_per_task", cls.cpus_per_task) or cls.cpus_per_task),
            env_activation=str(mlp.get("env_activation", "") or ""),
        )
        # Environment overrides: DELFIN_MLP_CLUSTER_<KEY> (upper case).
        mapping = {
            "gpu_partitions": "gpu_partitions",
            "partition": "partition",
            "gres": "gres",
            "time_limit": "time_limit",
            "mem_mb": "mem_mb",
            "cpus_per_task": "cpus_per_task",
            "env_activation": "env_activation",
        }
        for env_key, attr in mapping.items():
            value = env.get(f"DELFIN_MLP_CLUSTER_{env_key.upper()}")
            if value:
                setattr(cfg, attr,
                        int(value) if attr in ("mem_mb", "cpus_per_task") else str(value))
        return cfg

    def has_gpu_partition(self) -> bool:
        return bool(self.gpu_partitions.strip())


def decide_dispatch(cfg: MLPClusterConfig, *, force_local: bool = False) -> tuple:
    """-> (target, reason); target is 'cluster' or 'local'.

    Cluster only when a Slurm controller is reachable AND a GPU partition is
    configured.  Anything else is the local CPU fallback — the default.
    """
    if force_local:
        return "local", "forced local (force_local=True)"
    if not slurm_available():
        return "local", "no Slurm cluster reachable (sbatch/squeue not found)"
    if not cfg.has_gpu_partition():
        return "local", "no GPU partition configured (runtime.slurm.gpu_partitions is empty)"
    return "cluster", f"GPU partition {cfg.gpu_partitions} configured"


def local_plan(threads: int = LOCAL_THREADS, minutes: int = LOCAL_MINUTES) -> dict:
    """The capped CPU fallback plan.  Threads never exceed 2, minutes never 30."""
    threads = max(1, min(int(threads), LOCAL_THREADS))
    minutes = max(1, min(int(minutes), LOCAL_MINUTES))
    return {
        "threads": threads,
        "minutes": minutes,
        "threads_env": {
            "OMP_NUM_THREADS": str(threads),
            "MKL_NUM_THREADS": str(threads),
            "TORCH_NUM_THREADS": str(threads),
        },
    }


def count_frames(ensemble: Path) -> int:
    """Number of structures in a multi-frame XYZ file."""
    lines = Path(ensemble).read_text(encoding="utf-8").splitlines()
    n = 0
    i = 0
    while i < len(lines):
        try:
            natoms = int(lines[i].strip())
        except (ValueError, IndexError):
            i += 1
            continue
        n += 1
        i += natoms + 2
    return n


def plan_batch(ensemble) -> int:
    """A batch is ONE job over the whole ensemble file; returns N (structures)."""
    return count_frames(ensemble)


def jobs_for_batch(ensemble) -> list:
    """The sbatch job list for a batch: exactly one job, whatever N is."""
    return [{"batch_file": str(ensemble), "n_structures": count_frames(ensemble),
             "n_jobs": 1}]


def _activation_lines(cfg: MLPClusterConfig) -> list:
    lines = []
    for chunk in cfg.env_activation.split(";"):
        chunk = chunk.strip()
        if chunk:
            lines.append(chunk)
    return lines


def render_mlp_sbatch(cfg: MLPClusterConfig, *, job_dir, xyz, backend: str,
                      task: str = "optimize", job_name: str = "mlp_job",
                      batch=None, charge: int = 0, mult: int = 1,
                      cpus_per_task: int = None, time_limit: str = None) -> str:
    """The sbatch script text for one MLP job.  Writes nothing, submits nothing.

    ``task``: optimize | single_point | rank.  For ``rank``, ``batch`` is the
    multi-structure XYZ the job ranks whole — one job, not one per structure.
    """
    task = str(task or "optimize").strip().lower()
    if task not in VALID_TASKS:
        raise ValueError(f"unknown MLP task {task!r}: use one of {', '.join(VALID_TASKS)}")
    if task == "rank" and not batch:
        raise ValueError("task 'rank' needs the batch ensemble file")
    if not cfg.has_gpu_partition():
        raise ValueError(
            "no GPU partition configured — a GPU job cannot be rendered; "
            "configure runtime.slurm.gpu_partitions or run locally")
    cpus = int(cpus_per_task or cfg.cpus_per_task)
    tlim = str(time_limit or cfg.time_limit)
    job_dir = Path(job_dir)
    xyz = Path(xyz)
    env = {
        "DELFIN_MODE": "mlp",
        "DELFIN_JOB_NAME": job_name,
        "DELFIN_XYZ_FILE": str(xyz),
        "DELFIN_MLP_BACKEND": backend,
        "DELFIN_MLP_TASK": task,
        "DELFIN_CHARGE": str(int(charge)),
        "DELFIN_MULT": str(int(mult)),
        "DELFIN_PAL": str(cpus),
    }
    if task == "rank":
        env["DELFIN_MLP_BATCH_FILE"] = str(Path(batch))
    lines = [
        "#!/bin/bash",
        f"# delfin MLP {task} ({backend}) — generated by delfin.mlp_tools.cluster_dispatch",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --partition={cfg.gpu_partitions}",
        f"#SBATCH --gres={cfg.gres}",
        "--nodes=1",
        "--ntasks=1",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={cfg.mem_mb}M",
        f"#SBATCH --time={tlim}",
        f"#SBATCH --output={job_name}_%j.out",
        f"#SBATCH --error={job_name}_%j.err",
        "",
        "set -euo pipefail",
        # Thread caps: even on a GPU the CPU side stays bounded.
        "export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2",
        "export TORCH_NUM_THREADS=2",
        "export PYTHONNOUSERSITE=1",
        *_activation_lines(cfg),
        *[f"export {k}={v}" for k, v in env.items()],
        "",
        'echo "[$(date -Is)] MLP ${DELFIN_MLP_TASK} on $(hostname), '
        'gres=' + cfg.gres + '", partition=' + cfg.gpu_partitions + '"',
        "python -m delfin.dashboard.local_runner",
        "",
    ]
    return "\n".join(lines)


def write_and_submit(cfg: MLPClusterConfig, *, job_dir, xyz, backend: str,
                     task: str = "optimize", job_name: str = "mlp_job",
                     batch=None, charge: int = 0, mult: int = 1,
                     submit: bool = False, active_jobs: int = None,
                     sbatch: str = "sbatch") -> dict:
    """Write the script; run sbatch only when ``submit=True``.

    Dry run is the DEFAULT: no cluster side effect without an explicit
    request.  Returns {script, submitted, sbatch_output, n_structures}.
    """
    if submit:
        active = active_jobs if active_jobs is not None else count_active_mlp_jobs()
        assert_can_submit(active=active)
    script = render_mlp_sbatch(cfg, job_dir=job_dir, xyz=xyz, backend=backend,
                               task=task, job_name=job_name, batch=batch,
                               charge=charge, mult=mult)
    job_dir = Path(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)
    path = job_dir / f"{job_name}.sbatch"
    path.write_text(script, encoding="utf-8")
    path.chmod(0o755)
    out = {"script": str(path), "submitted": False, "sbatch_output": "",
           "n_structures": count_frames(batch) if batch else 1}
    if not submit:
        return out
    done = subprocess.run([sbatch, str(path)], cwd=str(job_dir),
                          capture_output=True, text=True)
    out["sbatch_output"] = (done.stdout + done.stderr).strip()
    if done.returncode != 0:
        raise RuntimeError(f"sbatch failed ({done.returncode}): {out['sbatch_output']}")
    out["submitted"] = True
    return out


def count_active_mlp_jobs() -> int:
    """This user's queued/running sbatch MLP jobs (job name prefix mlp)."""
    try:
        done = subprocess.run(["squeue", "-h", "-o", "%j"], capture_output=True,
                              text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return 0
    if done.returncode != 0:
        return 0
    return sum(1 for name in done.stdout.splitlines()
               if name.strip().startswith("mlp"))


def assert_can_submit(*, active: int) -> None:
    """Refuse the (MAX_CONCURRENT_MLP_JOBS + 1)-th concurrent sbatch MLP job."""
    if int(active) >= MAX_CONCURRENT_MLP_JOBS:
        raise JobCapExceeded(
            f"{active} sbatch MLP jobs already active; the limit is "
            f"{MAX_CONCURRENT_MLP_JOBS}. Wait for one to finish.")
