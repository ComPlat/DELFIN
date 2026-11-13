"""How a GOAT with an MLIP is run: which ExtOpt provider, how many workers.

``global_opt_method=mlip:<backend>`` names the model; this module decides the
rest, so that key keeps one meaning whichever provider runs it.

Provider (``global_opt_mlip_provider``: ``auto`` | ``oet`` | ``orca-mlips``):
``auto`` takes orca-external-tools when it is installed with that backend,
else orca-mlips.  Measured on an 81-atom chromophore (docs/MLIP_GOAT_GPU.md):
orca-mlips starts a Python process with torch for every ExtOpt call
(1.0-1.5 s per optimisation step against 0.2-0.3 s for xTB) and, when GOAT's
workers start together, one model server per worker; OET keeps one server
per job with a fixed number of model processes and a slim client.

Workers.  GOAT runs ``%pal nprocs`` worker optimisations in parallel (PAL 40
gives 40 workers).  With orca-mlips each worker ended up with a model of its
own on the GPU: 40 ran the GPU out of memory (CUDA allocation failures, the
job OOM-killed), 8 ran.  So orca-mlips GOAT is capped at
:func:`orca_mlips_worker_cap` workers -- 8, fewer on a small GPU -- written as
``%goat nworkers``; PAL itself is left as it is.  With OET each request
names its cores and the server refuses one larger than its process count,
which the initial optimisation (all of %pal) is: so %pal nprocs, the workers
and the server's model processes are one number, :func:`oet_default_workers`.
``goat_max_workers`` in CONTROL overrides the worker count (a positive
multiple of 4, ORCA's rule).
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from typing import Optional

PROVIDERS = ("auto", "oet", "orca-mlips")
#: GPU memory DELFIN plans per resident model (CUDA context, weights, work
#: space), GiB. Measured: 8 orca-external-tools AIMNet2 processes held
#: 6.4 GiB on an H100 after an ethanol GOAT (0.8 GiB each); planned with
#: margin for larger molecules. (40 orca-mlips servers did not fit 80 GiB:
#: they were started again and again, 259 times, not 40.)
GIB_PER_MODEL = 2.0
#: orca-mlips: measured upper bound of workers (one model each) that ran.
ORCA_MLIPS_MAX_WORKERS = 8
#: OET: upper bound of GOAT workers (= server model processes) when nothing
#: else is said. Measured: a GOAT step costs ~0.85 s of which the external
#: energy+gradient is ~0.22 s -- the rest is ORCA's own per-step work -- so
#: throughput comes from parallel workers, as with xTB: W follows PAL (the
#: archived xTB GOAT ran 40), bounded by what the GPU holds (0.8 GiB each).
OET_DEFAULT_WORKERS = 64
#: OET: send each request with curl instead of oet_client's Python start-up.
#: Off until measured faster on ethanol and on JRD1 (docs/MLIP_GOAT_GPU.md).
FAST_CLIENT_DEFAULT = False
#: Assumed when the GPU's memory cannot be read.
DEFAULT_GPU_GIB = 40.0


def gpu_status() -> dict:
    """``{name, total_gib, used_gib}`` of the first visible GPU (nvidia-smi), or {}."""
    exe = shutil.which("nvidia-smi")
    if not exe:
        return {}
    try:
        done = subprocess.run(
            [exe, "--query-gpu=name,memory.total,memory.used", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=20)
        name, total, used = [p.strip() for p in (done.stdout or "").splitlines()[0].split(",")]
        return {"name": name, "total_gib": float(total) / 1024.0, "used_gib": float(used) / 1024.0}
    except Exception:
        return {}


def gpu_memory_gib() -> float:
    """Memory of the first visible GPU in GiB (nvidia-smi), or DEFAULT_GPU_GIB."""
    return gpu_status().get("total_gib", DEFAULT_GPU_GIB)


def _models_that_fit(gpu_gib: float) -> int:
    return max(1, int(gpu_gib // GIB_PER_MODEL))


def orca_mlips_worker_cap(gpu_gib: float) -> int:
    """Most GOAT workers orca-mlips may run on one GPU: a multiple of 4, at least 4."""
    cap = min(ORCA_MLIPS_MAX_WORKERS, _models_that_fit(gpu_gib))
    return max(4, cap - cap % 4)


def job_memory_gib() -> Optional[float]:
    """The job's memory (SLURM_MEM_PER_NODE, MB) in GiB, or None outside SLURM."""
    import os

    try:
        return float(os.environ["SLURM_MEM_PER_NODE"]) / 1024.0
    except (KeyError, ValueError):
        return None


#: Host memory DELFIN plans per OET model process (a Python with torch and
#: CUDA), GiB.
HOST_GIB_PER_MODEL = 3.0


def oet_default_workers(pal: int, gpu_gib: float, host_gib: Optional[float] = None) -> int:
    """GOAT workers = model processes of an OET GOAT: PAL, within the GPU and the job memory.

    A multiple of 4 (GOAT spreads its workers over 4 temperatures), at least 4.
    """
    w = min(OET_DEFAULT_WORKERS, int(pal), _models_that_fit(gpu_gib))
    if host_gib is not None:
        w = min(w, max(1, int(host_gib // HOST_GIB_PER_MODEL)))
    return max(4, w - w % 4)


def default_goat_workers(pal: int) -> int:
    """GOAT's own worker count for ``%pal nprocs pal`` (ORCA: one per process)."""
    return max(1, int(pal))


@dataclass(frozen=True)
class GoatPlan:
    provider: str            # "oet" | "orca-mlips"
    nworkers: Optional[int]  # written as %goat nworkers, None = GOAT's own choice
    workers: int             # how many workers run
    server_processes: int    # OET model processes (1 for orca-mlips: per worker)
    reason: str
    nprocs: Optional[int] = None  # %pal nprocs written, None = PAL
    fast_client: bool = False  # OET: the curl client (FAST_CLIENT_DEFAULT)


def choose_provider(backend: str, forced: str = "auto", *, oet_info: Optional[dict] = None,
                    mlips_info: Optional[dict] = None) -> tuple[str, str]:
    """``(provider, reason)`` for ``mlip:<backend>``; ValueError when none can run it."""
    from delfin.mlp_tools.oet import OET_BACKEND_DISTS, OET_INSTALL, detect_oet
    from delfin.mlp_tools.orca_mlips import ORCA_MLIPS_INSTALL, detect_orca_mlips

    forced = (forced or "auto").strip().lower()
    if forced not in PROVIDERS:
        raise ValueError(f"global_opt_mlip_provider {forced!r} not one of {', '.join(PROVIDERS)}")
    oet = oet_info if oet_info is not None else detect_oet()
    oet_ok = bool(oet.get("available")) and backend in tuple(oet.get("backends") or ())
    if forced == "oet":
        if not oet_ok:
            why = ("orca-external-tools has no calculator for this backend"
                   if backend not in OET_BACKEND_DISTS else
                   f"orca-external-tools with {backend} is not installed ({OET_INSTALL})")
            raise ValueError(f"global_opt_mlip_provider=oet: {why}")
        return "oet", "forced in CONTROL"
    if forced == "auto" and oet_ok:
        return "oet", "orca-external-tools is installed with this backend (one server per job)"
    mlips = mlips_info if mlips_info is not None else detect_orca_mlips()
    if forced == "orca-mlips":
        return "orca-mlips", "forced in CONTROL"
    if mlips.get("available"):
        return "orca-mlips", "orca-external-tools is not installed with this backend"
    # Neither is installed: the orca-mlips input is written as before, and
    # ORCA stops on the missing ProgExt, which check_goat_e2e names.
    return "orca-mlips", (f"no ExtOpt provider is installed ({OET_INSTALL} or "
                          f"{ORCA_MLIPS_INSTALL})")


def plan_goat(backend: str, pal: int, *, forced: str = "auto",
              max_workers: Optional[int] = None, gpu_gib: Optional[float] = None,
              fast_client: Optional[bool] = None,
              oet_info: Optional[dict] = None, mlips_info: Optional[dict] = None) -> GoatPlan:
    """Provider, workers and server processes of one MLIP GOAT."""
    provider, why = choose_provider(backend, forced, oet_info=oet_info, mlips_info=mlips_info)
    if max_workers is not None:
        if isinstance(max_workers, bool) or int(max_workers) <= 0 or int(max_workers) % 4:
            raise ValueError(f"goat_max_workers {max_workers!r} must be a positive multiple of 4 "
                             "(ORCA manual Table 4.9: number of temperatures)")
    gib = gpu_memory_gib() if gpu_gib is None else float(gpu_gib)
    natural = default_goat_workers(pal)
    if provider == "orca-mlips":
        cap = orca_mlips_worker_cap(gib)
        if max_workers is not None:
            nworkers = int(max_workers)
            reason = f"{why}; goat_max_workers={nworkers} from CONTROL"
        elif natural > cap:
            nworkers = cap
            reason = (f"{why}; orca-mlips loads one model per GOAT worker: workers capped "
                      f"{natural} -> {cap} for a {gib:.0f} GiB GPU")
        else:
            nworkers, reason = None, why
        # Each capped worker gets one process, as in an uncapped PAL-8 run: with
        # nprocs 40 over 8 workers an ethanol GOAT did not finish in 25 min
        # (A100), with nprocs 8 it finished in 9.5 min (H100).
        return GoatPlan(provider, nworkers, nworkers or natural, 1, reason, nprocs=nworkers)
    # orca-external-tools: every request names its cores (NCores, ORCA's %pal
    # share) and the server refuses one larger than its -n -- the initial
    # optimisation asks for all of %pal nprocs ("Requested 32 cores but only 8
    # total available", measured). The server's -n is also its number of
    # model processes. So %pal nprocs, GOAT's workers (one process each) and
    # the server's model processes are one number W; the job keeps its PAL
    # cores, ORCA uses W of them -- the model runs on the GPU, not on them.
    if max_workers is not None:
        workers = int(max_workers)
        source = "goat_max_workers from CONTROL"
    else:
        host = job_memory_gib()
        workers = oet_default_workers(pal, gib, host)
        source = (f"PAL {pal}, {gib:.0f} GiB GPU"
                  + (f", {host:.0f} GiB job memory" if host else "") + " (oet_default_workers)")
    # %goat nworkers W as well: left to itself GOAT ran only 4 workers with
    # 6 processes each on %pal nprocs 24 (measured), not 24 workers.
    return GoatPlan("oet", workers, workers, workers,
                    f"{why}; {workers} GOAT workers (%pal nprocs {workers}) on one server "
                    f"with {workers} model processes, {source}", nprocs=workers,
                    fast_client=FAST_CLIENT_DEFAULT if fast_client is None else bool(fast_client))
