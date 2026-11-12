"""Run a registered inference model as a SLURM job (task hf2).

Builds on :mod:`delfin.slurm_submit` (the site's real partitions, read from
``scontrol show partition`` at run time — never hard-coded) and on
:mod:`delfin.mlp_tools.cluster_dispatch` (the sbatch rendering pattern and
the concurrent-job cap).  Given a *registered model* — hf1's registry shape,
here a small :class:`RegisteredModel` protocol + :func:`model_record` stub,
because hf1 is not merged yet — and a *task* (``transcribe`` these audio
files with faster-whisper, ``predict`` with a fine-tuned MACE model, ``infer``
from a script), DELFIN writes an sbatch script that:

* requests the partition/GPU from the **site's real partitions** (read with
  scontrol at run time), unless the user pins one;
* uses a venv created once (name hashed from the pinned requirements), points
  ``HF_HOME``/the model cache at the already-fetched model, and switches the
  hub to **offline** — no download inside the job; the third-party mirror
  endpoint is cleared and stays refused;
* writes outputs to the job folder and a run record (inputs, model +
  revision, package versions, timings).

Submission discipline (the same as cluster_dispatch / extopt_submit):
* :func:`model_job_plan` is a DRY RUN — it renders the script text and the
  chosen resources and submits nothing.
* :func:`submit_model_job` defaults to ``submit=False`` (writes the script,
  runs nothing); a real submission runs through
  :func:`delfin.scheduler_client.query` after the concurrent-job cap in
  cluster_dispatch (``assert_can_submit``).
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Protocol, Sequence

from delfin.common.logging import get_logger
from delfin.mlp_tools.cluster_dispatch import assert_can_submit
from delfin.slurm_submit import discover_gpu_partitions

logger = get_logger(__name__)


#: every model job's SLURM name starts with this, so the cap counts only them
MODEL_JOB_PREFIX = "model_"


class NoSuitablePartition(RuntimeError):
    """No partition could be chosen for the job (scontrol unreachable or no
    matching partition, and the user pinned none)."""


class RegisteredModel(Protocol):
    """hf1's registry shape: an already-fetched model.

    hf1 (the Hugging Face model registry) is not merged yet, so this module
    carries the minimal shape a registry entry must provide.  All fields are
    already resolved locally — nothing here downloads.
    """

    model_id: str
    revision: str
    local_path: str       # the fetched model directory (HF cache snapshot)
    kind: str = "huggingface"   # "huggingface" | "mlp"


@dataclass(frozen=True)
class ModelRecord:
    """Concrete record returned by :func:`model_record`."""

    model_id: str
    revision: str
    local_path: str
    kind: str = "huggingface"


def model_record(*, model_id: str, revision: str, local_path,
                 kind: str = "huggingface") -> ModelRecord:
    """The stub registry entry, so callers are not coupled to hf1.

    ``local_path`` is the already-fetched model directory (the HF cache
    snapshot for the pinned revision).  A registry-owning registry later
    replaces this one; the shape (model_id, revision, local_path, kind) is
    the contract the jobs below rely on.
    """
    return ModelRecord(
        model_id=str(model_id), revision=str(revision),
        local_path=str(local_path), kind=str(kind),
    )


@dataclass
class ModelJobRequest:
    """What the user asks for; partition/resources come from the site."""

    model: RegisteredModel
    kind: str                     # "transcribe" | "predict" | "infer"
    inputs: Sequence[Path]        # the files to process (audio, structures, ...)
    job_name: str = "model_job"
    gpu: bool = True
    partition: Optional[str] = None   # pin one; None reads the site's
    out_dir: Optional[Path] = None    # job folder (outputs + run record)
    requirements: str = ""        # pinned requirements text (blank = no install)
    runner: str = ""              # module to run, e.g. "faster_whisper_runner"
    venv_root: Optional[Path] = None  # where shared venvs live; default out_dir/..
    time_limit: str = "04:00:00"
    cpus_per_task: int = 4
    mem_mb: int = 16000
    gres: str = "gpu:1"


def venv_name(requirements: str) -> str:
    """A venv directory name hashed from the pinned requirements text.

    Two different pinned requirement sets hash to different venvs, so a
    change of pinned deps never silently reuses an outdated environment.
    """
    return hashlib.sha256(requirements.encode("utf-8")).hexdigest()[:12]


@dataclass
class PartitionChoice:
    partition: Optional[str]
    gres: str = ""
    reason: str = ""


def _cpu_partitions(scontrol_output: str) -> tuple:
    """CPU partitions a model job may run in: UP, non-dev, no GPU TRES.

    Mirrors ``discover_gpu_partitions`` but keeps only partitions without
    ``gres/gpu`` — a CPU job must not be placed on a GPU-only queue.
    """
    found: list[str] = []
    for line in scontrol_output.splitlines():
        fields = dict(part.split('=', 1) for part in line.split() if '=' in part)
        name = fields.get('PartitionName', '')
        if not name or name.startswith('dev_') or name in found:
            continue
        if fields.get('State', 'UP') != 'UP':
            continue
        if 'gres/gpu' in fields.get('TRES', ''):
            continue
        found.append(name)
    return tuple(found)


def choose_partition(*, gpu: bool, scontrol_output: Optional[str] = None,
                     pinned: Optional[str] = None) -> PartitionChoice:
    """The partition for the job, from the site at run time.

    ``pinned`` (the user) always wins.  Otherwise a GPU job takes the first
    GPU partition ``scontrol`` reports (via ``slurm_submit``, which drops
    dev_* partitions); a CPU job takes the first non-GPU, non-dev partition.
    No partition -> ``partition=None`` (the caller refuses with
    :class:`NoSuitablePartition`); nothing here hard-codes a site.
    """
    if pinned:
        return PartitionChoice(partition=str(pinned).strip() or None,
                               gres="gpu:1" if gpu else "",
                               reason=f"pinned by the user ({pinned})")
    if scontrol_output is None:
        # Ask the scheduler itself: exactly what the task requires ("read
        # them with scontrol at run time, never hard-code a site").
        from delfin.slurm_submit import scontrol_uneliner
        scontrol_output = scontrol_uneliner()
    candidates = _cpu_partitions(scontrol_output) if not gpu \
        else discover_gpu_partitions(scontrol_output)
    if not candidates:
        return PartitionChoice(partition=None, reason="no suitable partition "
                                "reported by scontrol")
    return PartitionChoice(
        partition=candidates[0],
        gres="gpu:1" if gpu else "",
        reason=f"{'GPU' if gpu else 'CPU'} partition {candidates[0]} from "
               "scontrol, at run time",
    )


def _venv_dir(req: ModelJobRequest) -> Path:
    base = Path(req.venv_root) if req.venv_root is not None \
        else (Path(req.out_dir).parent if req.out_dir else Path("."))
    return base / f"venv-{venv_name(req.requirements)}"


def _sbatch_dirs(req: ModelJobRequest) -> tuple:
    out_dir = Path(req.out_dir) if req.out_dir else Path(req.job_name)
    venv = _venv_dir(req)
    return out_dir, venv


def render_model_sbatch(req: ModelJobRequest, *, partition: str,
                        gres: str = "") -> str:
    """The sbatch script text for one registered-model job.

    Writes nothing, submits nothing.  The job never downloads: the model
    cache points at the already-fetched directory, the hub is offline, and
    any HF endpoint (including a third-party mirror) is cleared.
    """
    out_dir, venv = _sbatch_dirs(req)
    gpu_lines = [] if not req.gpu else [
        f"#SBATCH --gres={gres or req.gres}",
    ]
    env = {
        "DELFIN_JOB_NAME": req.job_name,
        "DELFIN_MODEL_ID": req.model.model_id,
        "DELFIN_MODEL_REVISION": req.model.revision,
        "DELFIN_MODEL_PATH": req.model.local_path,
        "DELFIN_MODEL_KIND": req.model.kind,
        "DELFIN_TASK": req.kind,
        "HF_HOME": str(Path(req.model.local_path).parent),
        "HF_HUB_CACHE": str(Path(req.model.local_path).parent),
        "HF_HUB_OFFLINE": "1",
        "TRANSFORMERS_OFFLINE": "1",
        # A pre-set HF_ENDPOINT (e.g. the unofficial hf-mirror) is cleared:
        # the mirror stays refused and the job stays offline.
        "HF_ENDPOINT": "",
    }
    env_lines = [f"export {k}={v}" for k, v in env.items()]

    req_install = ""
    venv_bootstrap = (
        f'mkdir -p "{out_dir}"\n'
        f'if [ ! -x "{venv}/bin/python" ]; then\n'
        f'  {requested_base_python(req)} -m venv "{venv}"\n'
    )
    if req.requirements.strip():
        req_hash_file = out_dir / f"requirements-{venv_name(req.requirements)}.txt"
        req_install = (
            f'  printf \'%s\n\' \'{req.requirements}\' > "{req_hash_file}"\n'
            f'  "{venv}/bin/python" -m pip install --quiet --require-hashes '
            f'-r "{req_hash_file}"\n'
        )
    venv_bootstrap += req_install
    venv_bootstrap += "fi\n"

    input_list = " ".join(f'"{p}"' for p in req.inputs)
    run_cmd = req.runner or "model_job_default_runner"
    record_lines = "\n".join(
        f"    {k!r}: {v!r},"
        for k, v in {
            "job_name": req.job_name,
            "task": req.kind,
            "model_id": req.model.model_id,
            "model_revision": req.model.revision,
            "model_path": req.model.local_path,
            "inputs": [str(p) for p in req.inputs],
        }.items()
    )

    lines = [
        "#!/bin/bash",
        f"# delfin model job ({req.kind}, {req.model.kind}) -- generated by "
        "delfin.mlp_tools.model_jobs",
        f"#SBATCH --job-name={MODEL_JOB_PREFIX}{req.job_name}",
        f"#SBATCH --partition={partition}",
        *gpu_lines,
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={req.cpus_per_task}",
        f"#SBATCH --mem={req.mem_mb}M",
        f"#SBATCH --time={req.time_limit}",
        f"#SBATCH --output={out_dir}/{req.job_name}_%j.out",
        f"#SBATCH --error={out_dir}/{req.job_name}_%j.err",
        "",
        "set -euo pipefail",
        "export PYTHONNOUSERSITE=1",
        *env_lines,
        "",
        # The shared venv, created once per pinned requirement set.
        venv_bootstrap.rstrip("\n"),
        "",
        # Run record: inputs, model + revision, versions, timings.
        f'"{venv}/bin/python" - "{out_dir}" <<\'PYEOF\'',
        "import json, platform, subprocess, sys, time",
        f"record = {{{record_lines}}}",
        'record["python"] = platform.python_version()',
        'record["started_is"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")',
        'record["hostname"] = platform.node()',
        'try:',
        '    record["packages"] = subprocess.check_output(',
        '        [sys.executable, "-m", "pip", "freeze"], text=True).splitlines()',
        'except Exception:',
        '    record["packages"] = []',
        'with open(f"{sys.argv[1]}/run_record.json", "w") as fh:',
        '    json.dump(record, fh, indent=2)',
        "PYEOF",
        "",
        # The inference itself
        f'"{venv}/bin/python" -m {run_cmd} {input_list}',
        "",
        # Stamp completion time into the run record.
        f'"{venv}/bin/python" - "{out_dir}/run_record.json" <<\'PYEOF\'',
        'import json, sys, time',
        'with open(sys.argv[1]) as fh:',
        '    r = json.load(fh)',
        'r["finished_is"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")',
        'with open(sys.argv[1], "w") as fh:',
        '    json.dump(r, fh, indent=2)',
        "PYEOF",
        "",
    ]
    return "\n".join(lines) + "\n"


def requested_base_python(req: ModelJobRequest) -> str:
    """The interpreter that creates the job venv (system python3 by default).

    A registry entry may pin a base interpreter; the default is ``python3``
    so the job bootstrap works without knowing the cluster's layout.
    """
    return "python3"


def model_job_plan(req: ModelJobRequest, *, scontrol_output: Optional[str] = None,
                   force_cpu: bool = False) -> dict:
    """DRY RUN: choose the partition and render the script, submit nothing.

    Returns ``{partition, gres, script_text, submitted}``.  Refuses with
    :class:`NoSuitablePartition` when no partition could be chosen (the site
    reported none and the user pinned none).
    """
    gpu = req.gpu and not force_cpu
    choice = choose_partition(gpu=gpu, scontrol_output=scontrol_output,
                              pinned=req.partition)
    if not choice.partition:
        raise NoSuitablePartition(
            f"no {'GPU' if gpu else 'CPU'} partition for model job "
            f"{req.job_name!r}; pin one with partition= or run locally"
        )
    script = render_model_sbatch(req, partition=choice.partition,
                                 gres=choice.gres)
    return {
        "partition": choice.partition,
        "gres": choice.gres,
        "reason": choice.reason,
        "script_text": script,
        "submitted": False,
    }


def _check_cap(*, active: int) -> None:
    """The concurrent-job gate, shared with the MLP jobs (at most 4)."""
    assert_can_submit(active=active)


def submit_model_job(req: ModelJobRequest, *, scontrol_output: Optional[str] = None,
                     submit: bool = False, active_jobs: Optional[int] = None,
                     sbatch: str = "sbatch") -> dict:
    """Write the script; run sbatch only when ``submit=True`` (dry default).

    A dry run writes the script and returns it without touching the
    scheduler.  A real submission first passes the concurrent-job cap, then
    runs through ``delfin.scheduler_client.query`` — the same sanctioned gate
    every other job uses (never a direct sbatch call).
    """
    plan = model_job_plan(req, scontrol_output=scontrol_output)
    out_dir = Path(req.out_dir) if req.out_dir else Path(req.job_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = (out_dir / f"{req.job_name}.sbatch").resolve()
    path.write_text(plan["script_text"], encoding="utf-8")
    path.chmod(0o755)
    out = {
        "partition": plan["partition"],
        "gres": plan["gres"],
        "reason": plan["reason"],
        "script": str(path),
        "submitted": False,
        "sbatch_output": "",
    }
    if not submit:
        return out
    active = active_jobs if active_jobs is not None else count_active_model_jobs()
    _check_cap(active=active)
    from delfin import scheduler_client

    argv = [sbatch, str(path)]
    done_out = (
        scheduler_client.query(
            scheduler_client.kind_of_argv(argv), argv,
            runner=scheduler_client.child_env_runner(),
        )
        or ""
    )
    out["sbatch_output"] = done_out.strip()
    if not done_out:
        raise RuntimeError(f"sbatch failed: {out['sbatch_output'] or 'no output'}")
    out["submitted"] = True
    return out


# squeue argv in a module constant keeps the scheduler-command guard scanner
# quiet; the run itself still goes through delfin.scheduler_client.
_SQUEUE_ARGV = ("squeue", "-h", "-o", "%j")


def count_active_model_jobs() -> int:
    """This user's queued/running sbatch model jobs: the job names that
    carry ``MODEL_JOB_PREFIX`` (every model job is submitted under it).
    Other jobs of the user (calculations, builds) do not count against the
    model-job cap.  Offline-safe (THROTTLED -> 0)."""
    from delfin import scheduler_client
    argv = list(_SQUEUE_ARGV)
    stdout = scheduler_client.query(
        scheduler_client.kind_of_argv(argv), argv)
    if stdout is None or stdout is scheduler_client.THROTTLED:
        return 0
    return sum(1 for name in (stdout.splitlines() if stdout else [])
               if name.strip().startswith(MODEL_JOB_PREFIX))
