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

GOAT (``run_type="goat"``) runs only with ``provider="orca_mlips"`` on a
GPU partition: GOAT needs on the order of 100 x N_atoms optimisations
(manual 4.10), the delfin-mlp-extopt wrapper reloads the model on every
ProgExt call, and the local CPU fallback (<= 2 cores, < 30 min) cannot
hold such a run.  The job holds ORCA with its GOAT workers (sharing
``%pal nprocs``) and one orca-mlips server.  Whether one server serves
several parallel workers cleanly is not verified (the upstream README
shows only ``%pal nprocs 8``); ``goat_max_workers`` caps the workers
until the first real job has shown it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from delfin.common.logging import get_logger
from delfin.mlp_tools.cluster_dispatch import MLPClusterConfig
from delfin.mlp_tools.orca_extopt import build_extopt_input
from delfin.mlp_tools.orca_mlips import (
    GOAT_DEFAULT_WORKERS,
    build_orca_mlips_input,
    goat_output_files,
    normalize_goat_options,
)

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
    # GOAT (orca_mlips only): %goat options, Table 4.8 variant, worker cap
    goat_options: Optional[dict] = None
    goat_variant: str = "GOAT"
    goat_max_workers: Optional[int] = None
    # NEB / NEB-TS end geometry and relaxed-scan dict (orca_mlips builder)
    neb_end_coordinates: Optional[Sequence[Sequence[float]]] = None
    scan: Optional[dict] = None


@dataclass
class ExtOptPlan:
    mode: str            # "gpu" | "local" | "refuse"
    input: Path
    script: Path
    pal: int = 1
    reason: str = ""
    #: Files the job will write in its work_dir (GOAT: output, final
    #: ensemble, global minimum); empty for the other run types.
    produces: tuple[Path, ...] = ()


def _is_goat(req: ExtOptJobRequest) -> bool:
    return str(req.run_type).lower().strip() == "goat"


def goat_workers(req: ExtOptJobRequest) -> tuple[int, dict]:
    """(number of GOAT workers, validated %goat options) for ``req``.

    Workers: ``goat_max_workers`` caps ``%goat nworkers``; without either,
    GOAT's default of 8 (manual 4.10.1).  The cap is written as nworkers
    so ORCA starts no more workers than the job has cores for.  The cap
    must be a positive multiple of 4 (Table 4.9).
    """
    opts = normalize_goat_options(req.goat_options)
    given = opts.get("nworkers")
    cap = req.goat_max_workers
    if cap is not None:
        if isinstance(cap, bool) or not isinstance(cap, int) or cap <= 0 or cap % 4:
            raise ValueError(
                f"goat_max_workers {cap!r} must be a positive multiple of 4 "
                "(manual Table 4.9: number of temperatures)"
            )
        workers = min(cap, given) if isinstance(given, int) else cap
        opts["nworkers"] = workers
    else:
        workers = given if isinstance(given, int) else GOAT_DEFAULT_WORKERS
    return int(workers), opts


def _build_input(req: ExtOptJobRequest, *, ncores: int = 1,
                 goat_options: Optional[dict] = None) -> Path:
    if req.provider == "orca_mlips":
        extra: dict = {}
        if _is_goat(req):
            extra = dict(goat_options=goat_options, goat_variant=req.goat_variant,
                         ncores=ncores)
        else:
            extra["ncores"] = ncores
        if req.neb_end_coordinates is not None:
            extra["neb_end_coordinates"] = req.neb_end_coordinates
        if req.scan is not None:
            extra["scan"] = req.scan
        return build_orca_mlips_input(
            req.run_type, elements=req.elements, coordinates=req.coordinates,
            charge=req.charge, mult=req.mult, backend=req.mlips_backend,
            model=req.model, device=req.device, solvent=req.solvent,
            solvent_model=req.solvent_model, hessian=req.hessian,
            base_name=req.job_name, work_dir=req.work_dir, **extra,
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
    if _is_goat(req):
        return _plan_goat_job(req, cfg=cfg)
    gpu = cfg.has_gpu_partition()
    planned_pal = (
        int(cfg.cpus_per_task) if gpu else
        min(LOCAL_CPU_LIMITS.cores, max(1, len(req.elements) // 4 + 1))
    )
    if planned_pal < 1:
        raise ValueError("cpus_per_task must be positive")
    inp = _build_input(req, ncores=planned_pal)
    job_dir = Path(req.work_dir)
    if gpu:
        pal = planned_pal
        script = job_dir / f"{req.job_name}.sbatch"
        script.write_text(render_extopt_job_script(
            req, cfg=cfg, pal=pal, provider_cmd=_provider_cmd(req),
        ))
        return ExtOptPlan(mode="gpu", input=inp, script=script, pal=pal)
    # Never silently run a requested CUDA calculation on CPU.
    if str(req.device).lower().startswith("cuda"):
        return ExtOptPlan(
            mode="refuse", input=inp,
            script=job_dir / f"{req.job_name}.sbatch",
            pal=1,
            reason="CUDA requested but no GPU partition configured",
        )
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


def _plan_goat_job(req: ExtOptJobRequest, *, cfg: MLPClusterConfig) -> ExtOptPlan:
    """GOAT: orca-mlips on a GPU node, or a refusal; never a local run."""
    job_dir = Path(req.work_dir)
    inp_path = job_dir / f"{req.job_name}.inp"
    script = job_dir / f"{req.job_name}.sbatch"
    if req.provider != "orca_mlips":
        return ExtOptPlan(
            mode="refuse", input=inp_path, script=script, pal=1,
            reason=(
                "GOAT runs only with provider='orca_mlips': GOAT makes on the "
                "order of 100 x N_atoms optimisations (manual 4.10) and the "
                "delfin-mlp-extopt wrapper reloads the model on every ProgExt "
                "call; orca-mlips keeps one resident model server per job"
            ),
        )
    workers, opts = goat_workers(req)
    if not cfg.has_gpu_partition():
        return ExtOptPlan(
            mode="refuse", input=inp_path, script=script, pal=1,
            reason=(
                "GOAT with an MLIP needs a GPU partition "
                "(runtime.slurm.gpu_partitions); no local CPU fallback: "
                "the local limits (<= 2 cores, < 30 min) cannot hold a "
                "GOAT run"
            ),
        )
    pal = max(int(cfg.cpus_per_task), workers)
    inp = _build_input(req, ncores=pal, goat_options=opts)
    script.write_text(render_extopt_job_script(
        req, cfg=cfg, pal=pal, provider_cmd=_provider_cmd(req), goat_workers=workers,
    ))
    files = goat_output_files(req.job_name)
    produces = tuple(job_dir / files[k] for k in ("output", "finalensemble", "global_minimum"))
    return ExtOptPlan(mode="gpu", input=inp, script=script, pal=pal,
                      reason=f"GOAT: {workers} workers on {pal} cores, one GPU node",
                      produces=produces)


def _provider_cmd(req: ExtOptJobRequest) -> str:
    if req.provider == "orca_mlips":
        return f"orca-mlips-{req.mlips_backend}"
    return "delfin-mlp-extopt"


def render_extopt_job_script(req: ExtOptJobRequest, *, cfg: MLPClusterConfig,
                             pal: int, provider_cmd: str,
                             goat_workers: Optional[int] = None) -> str:
    """The sbatch script text for one ORCA+ExtOpt GPU job.

    Same-node guarantee: ``--nodes=1`` puts ORCA and the ProgExt command
    (and the orca-mlips resident server) on one node.  Server lifetime is
    bound to the job: orca-mlips stops its auto-server with the parent
    ORCA process; the cleanup line covers a job killed mid-run.
    ``%pal`` is part of the ORCA input, never of this script.

    ``goat_workers`` (GOAT jobs): ORCA is started by its full path (ORCA
    needs it for parallel runs) and the script records the worker count
    and the unverified one-server-for-all-workers assumption.
    """
    job_dir = Path(req.work_dir)
    is_mlips = req.provider == "orca_mlips"
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
    if goat_workers is not None:
        lines.extend([
            f"# GOAT: {goat_workers} workers share %pal nprocs {pal} in the "
            "input. The MLIP side is",
            "# orca-mlips' auto-server (started by the first ProgExt call, "
            "stopped by the trap);",
            "# no second server is started here. Unverified: whether the "
            "parallel workers share",
            "# one server cleanly -- check it in the first real job (worker "
            "outputs, GPU memory,",
            "# server processes); goat_max_workers caps the workers.",
            'ORCA_BIN="$(command -v orca)"',
            '"$ORCA_BIN" ' + req.job_name + ".inp > " + req.job_name + ".out",
        ])
    else:
        lines.extend([
            'ORCA_BIN="$(command -v orca)"',
            '"$ORCA_BIN" ' + req.job_name + ".inp > " + req.job_name + ".out",
        ])
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
        "produces": [str(p) for p in plan.produces],
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
    import re
    if re.fullmatch(r"Submitted batch job [0-9]+", out["sbatch_output"]) is None:
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
