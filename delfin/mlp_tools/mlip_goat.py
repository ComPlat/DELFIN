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
``%goat nworkers``; PAL itself is left as it is.  With OET the workers are
not capped: they queue on one server, whose model processes
(:func:`oet_server_processes`) are what the GPU memory has to hold.
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
#: space for a molecule of ~100 atoms), GiB.  Conservative: 40 orca-mlips
#: servers did not fit 80 GiB, 8 did.
GIB_PER_MODEL = 5.0
#: orca-mlips: measured upper bound of workers (one model each) that ran.
ORCA_MLIPS_MAX_WORKERS = 8
#: OET: upper bound of model processes in the one server.
OET_MAX_PROCESSES = 8
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


def oet_server_processes(workers: int, gpu_gib: float) -> int:
    """Model processes of the one OET server: no more than the GPU holds or workers need."""
    return max(1, min(OET_MAX_PROCESSES, _models_that_fit(gpu_gib), int(workers)))


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
              server_processes: Optional[int] = None,
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
        return GoatPlan(provider, nworkers, nworkers or natural, 1, reason)
    nworkers = int(max_workers) if max_workers is not None else None
    workers = nworkers or natural
    # A worker's request asks for its share of %pal (NCores); the server's
    # -n is also its core budget, so it must cover one request.
    per_request = max(1, natural // workers)
    if server_processes is not None:
        processes = max(per_request, int(server_processes))
        source = "mlip_server_processes from CONTROL"
    else:
        processes = max(per_request, oet_server_processes(workers, gib))
        source = f"{gib:.0f} GiB GPU"
    return GoatPlan("oet", nworkers, workers, processes,
                    f"{why}; {workers} workers share one server with {processes} model "
                    f"processes ({source})")
