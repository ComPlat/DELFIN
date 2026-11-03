"""GPU-aware submission of an ORCA+ExtOpt job (task om1 item 2).

One sbatch job carries BOTH the ORCA process and the MLIP side (the
ProgExt wrapper, and for orca-mlips its resident model server), so
ORCA's ExtOpt calls and the server always share a node — the
``--nodes=1`` line is the same-node guarantee.  With a GPU partition
configured the job goes there (partition from
``MLPClusterConfig.gpu_partitions``, ``--gres`` from the config, CPU
count for ORCA's ``%pal``/NCores, the MLIP stack's venv/conda activated
inside the script); without one, the job is refused or, when explicitly
allowed, planned as a local CPU run within the host limits
(<= 2 cores, < 30 min — the machine rule of the CPU-only host).

Submission discipline (same as cluster_dispatch):
* ``submit=False`` (default) is a dry run: input + job script are written
  and returned, nothing is submitted, the scheduler is not touched.
* A real submission runs through :func:`delfin.scheduler_client.query`
  (the direct-call guard test enforces that for every module under
  delfin/), with the child environment from
  :func:`delfin.agent.lifeline.child_env`.
* Server lifetime is bound to the job: orca-mlips scopes its auto-servers
  per parent ORCA process and stops them when it exits; for the
  delfin-mlp-extopt wrapper (no built-in server) the script traps job end
  and calls ``<provider> --stop-server`` as a no-op-safe cleanup.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from delfin.common.logging import get_logger
from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
from delfin.mlp_tools.orca_extopt import build_extopt_input
from delfin.mlp_tools.orca_mlips import build_orca_mlips_input

logger = get_logger(__name__)


@dataclass(frozen=True)
class LocalCpuLimits:
    cores: int = 2
    minutes: int = 30


LOCAL_CPU_LIMITS = LocalCpuLimits()

_PROVIDERS = ("delfin", "orca_mlips")


@dataclass
class ExtOptJobRequest:
    """What the user asks for; everything else comes from config."""

    run_type: str                       # opt | optts | freq | neb-ts | ...
    elements: Sequence[str]
    coordinates: Sequence[Sequence[float]]
    charge: int = 0
    mult: int = 1
    provider: str = "delfin"            # "delfin" | "orca_mlips"
    backend: str = "ani2x"              # delfin wrapper backend
    mlips_backend: str = "uma"          # orca-mlips plugin (uma|orb|mace|aimnet2)
    model: Optional[str] = None
    device: str = "cpu"
    solvent: Optional[str] = None
    solvent_model: str = "alpb"
    hessian: Optional[str] = None       # None | "two_step" (orca_mlips only)
    job_name: str = "extopt_job"
    allow_local: bool = False           # explicit opt-in to the CPU fallback
    max_local_atoms: int = 30
    work_dir: Path = field(default_factory=Path)


@dataclass
class ExtOptPlan:
    mode: str            # "gpu" | "local" | "refuse"
    input: Path
    script: Path
    pal: int = 1
    reason: str = ""


def _build_input(req: ExtOptJobRequest) -> Path:
    if req.provider == "orca_mlips":
        return build_orca_mlips_input(
            req.run_type, elements=req.elements, coordinates=req.coordinates,
            charge=req.charge, mult=req.mult, backend=req.mlips_backend,
            model=req.model, device=req.device, solvent=req.solvent,
            solvent_model=req.solvent_model, hessian=req.hessian,
            base_name=req.job_name, work_dir=req.work_dir,
        )
    if req.provider != "delfin":
        raise ValueError(
            f"provider {req.provider!r} not one of {_PROVIDERS} "
            "(delfin-mlp-extopt stays the default)"
        )
    return build_extopt_input(
        req.run_type, elements=req.elements, coordinates=req.coordinates,
        charge=req.charge, mult=req.mult, backend=req.backend,
        base_name=req.job_name, work_dir=req.work_dir,
    )


def plan_extopt_job(req: ExtOptJobRequest, *, cfg: MLPClusterConfig) -> ExtOptPlan:
    """Build the input and decide gpu vs local vs refuse. Submits nothing."""
    inp = _build_input(req)
    job_dir = Path(req.work_dir)
    if cfg.has_gpu_partition():
        pal = int(cfg.cpus_per_task)
        script = job_dir / f"{req.job_name}.sbatch"
        script.write_text(render_extopt_job_script(
            req, cfg=cfg, pal=pal, provider_cmd=_provider_cmd(req),
        ))
        return ExtOptPlan(mode="gpu", input=inp, script=script, pal=pal)
    # no GPU partition: local CPU within the host limits, or refuse
    n_atoms = len(req.elements)
    if not req.allow_local or n_atoms > req.max_local_atoms:
        reason = (
            f"no GPU partition configured and the local CPU fallback is "
            f"limited to <= {LOCAL_CPU_LIMITS.cores} cores / "
            f"< {LOCAL_CPU_LIMITS.minutes} min for small systems "
            f"(this request: {n_atoms} atoms, allow_local={req.allow_local}); "
            "configure runtime.slurm.gpu_partitions to run on a GPU"
        )
        return ExtOptPlan(
            mode="refuse", input=inp,
            script=job_dir / f"{req.job_name}.sbatch", pal=1, reason=reason,
        )
    pal = min(LOCAL_CPU_LIMITS.cores, max(1, n_atoms // 4 + 1))
    script = job_dir / f"{req.job_name}.local.sh"
    script.write_text(render_local_script(req, pal=pal))
    return ExtOptPlan(mode="local", input=inp, script=script, pal=pal,
                      reason=f"local CPU run, <= {pal} cores")


def _provider_cmd(req: ExtOptJobRequest) -> str:
    if req.provider == "orca_mlips":
        return f"orca-mlips-{req.mlips_backend}"
    return "delfin-mlp-extopt"


def render_extopt_job_script(req: ExtOptJobRequest, *, cfg: MLPClusterConfig,
                             pal: int, provider_cmd: str) -> str:
    """The sbatch script text for one ORCA+ExtOpt GPU job.

    Same-node guarantee: ``--nodes=1`` puts ORCA and the ProgExt command
    (and the orca-mlips resident server) on one node.  Server lifetime is
    bound to the job: orca-mlips stops its auto-server with the parent
    ORCA process; the cleanup line covers a job killed mid-run.
    """
    job_dir = Path(req.work_dir)
    is_mlips = req.provider == "orca_mlips"
    pal_line = f"%pal\n  nprocs {pal}\nend" if req.provider == "delfin" else ""
    lines = [
        "#!/bin/bash",
        f"# delfin ORCA+ExtOpt ({req.provider}/{req.mlips_backend or req.backend}) "
        "-- generated by delfin.mlp_tools.extopt_submit",
        f"#SBATCH --job-name={req.job_name}",
        f"#SBATCH --partition={cfg.gpu_partitions}",
        f"#SBATCH --gres={cfg.gres}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={pal}",
        f"#SBATCH --mem={cfg.mem_mb}M",
        f"#SBATCH --time={cfg.time_limit}",
        f"#SBATCH --output={req.job_name}_%j.out",
        f"#SBATCH --error={req.job_name}_%j.err",
        "",
        "set -euo pipefail",
        "export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2",
        "export TORCH_NUM_THREADS=2",
        "export PYTHONNOUSERSITE=1",
    ]
    for chunk in cfg.env_activation.split(";"):
        chunk = chunk.strip()
        if chunk:
            lines.append(chunk)
    lines.extend([
        "",
        f'cd "{job_dir}"',
        "",
        "# Server lifetime bound to the job: stop the resident model "
        "server when the job ends (no-op if it already exited).",
        f'trap \'{provider_cmd} --stop-server 2>/dev/null || true\' EXIT',
        "",
    ])
    if pal_line:
        lines.extend([pal_line, ""])
    lines.append(
        "orca " + req.job_name + ".inp > " + req.job_name + ".out"
    )
    if is_mlips:
        lines.append(
            f"# orca-mlips stops its auto-server with the parent ORCA "
            f"process; the EXIT trap above is the safety net."
        )
    return "\n".join(lines) + "\n"


def render_local_script(req: ExtOptJobRequest, *, pal: int) -> str:
    """Local CPU run within the host limits (<= 2 cores, < 30 min)."""
    provider_cmd = _provider_cmd(req)
    lines = [
        "#!/bin/bash",
        f"# delfin ORCA+ExtOpt LOCAL CPU ({req.provider}) -- <= {pal} cores, "
        f"< {LOCAL_CPU_LIMITS.minutes} min",
        "set -euo pipefail",
        f"export OMP_NUM_THREADS={pal} MKL_NUM_THREADS={pal}",
        "export TORCH_NUM_THREADS=" + str(pal),
        "ulimit -t " + str(LOCAL_CPU_LIMITS.minutes * 60),
        f'cd "{Path(req.work_dir)}"',
        f"trap '{provider_cmd} --stop-server 2>/dev/null || true' EXIT",
        "",
        "orca " + req.job_name + ".inp > " + req.job_name + ".out",
    ]
    return "\n".join(lines) + "\n"


def submit_extopt_job(req: ExtOptJobRequest, *, cfg: MLPClusterConfig,
                      submit: bool = False) -> dict:
    """Plan the job; run sbatch only when ``submit=True`` (dry run default).

    A dry run writes the input + script and returns them without touching
    the scheduler.  A real submission goes through
    ``delfin.scheduler_client.query`` (never a direct sbatch call), with
    the child environment from ``delfin.agent.lifeline.child_env``.
    """
    plan = plan_extopt_job(req, cfg=cfg)
    out = {
        "mode": plan.mode,
        "input": str(plan.input),
        "script": str(plan.script),
        "pal": plan.pal,
        "reason": plan.reason,
        "submitted": False,
        "sbatch_output": "",
    }
    if plan.mode == "refuse":
        raise RuntimeError(plan.reason)
    if not submit:
        return out
    if plan.mode == "local":
        raise RuntimeError(
            "local CPU fallback runs in the foreground on this host, not "
            "through sbatch; execute the script directly: " + str(plan.script)
        )
    from delfin import scheduler_client
    from delfin.agent.lifeline import child_env

    argv = ["sbatch", str(plan.script)]
    done_out = (
        scheduler_client.query(
            scheduler_client.kind_of_argv(argv), argv, runner=_runner_with(child_env())
        )
        or ""
    )
    out["sbatch_output"] = done_out.strip()
    if not done_out:
        raise RuntimeError(f"sbatch failed: {out['sbatch_output'] or 'no output'}")
    out["submitted"] = True
    return out


def _runner_with(env: dict):
    """A scheduler_client runner that carries the child env."""
    import subprocess

    def _run(argv):
        proc = subprocess.run(argv, capture_output=True, text=True, env=env)
        return proc.stdout if proc.returncode == 0 else (
            proc.stdout + proc.stderr)

    return _run
