"""BG1: cost-and-time estimate and per-session budget before every calculation.

Package BG1 prevents runaway calculations on the cluster: before anything is
submitted, ``estimate`` turns a job's system size and method into core-hours,
memory, wall time and an expected queue wait, and a session ``Budget`` decides
from the remaining core-hours whether the job may go or the person must be
asked.  The agent never silently shrinks PAL/maxcore/wall time to make a job
fit.

The queue-wait part goes through the ONE throttled scheduler client
(:mod:`delfin.scheduler_client`), never a per-loop or per-widget ``squeue`` or
``sacct``.  When the gate is closed (:data:`scheduler_client.THROTTLED`) the
estimate reports ``"unknown"`` for the wait instead of inventing a number.
"""

from __future__ import annotations

import math
import re
from datetime import datetime
from dataclasses import dataclass
from typing import Dict, Optional, Sequence

from delfin.scheduler_client import (
    SchedulerClient,
    THROTTLED,
    default_client,
    kind_of_argv,
)

#: Sentinel for a wait we could not learn from the scheduler's own answer.
UNKNOWN = "unknown"

#: The sbatch command the estimate asks through the scheduler client.
SBATCH = "sbatch"

#: Real ``sbatch --test-only`` answers with ``to start at <datetime>`` (the
#: marker slurm_submit.partition_accepts keys on).  The estimate turns that
#: start time into a seconds-until-start wait, and reports "unknown" when the
#: answer does not state one -- never an invented number.
_START_AT_RE = re.compile(r"to\s+start\s+at\s+(\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2})")

#: Wall-time base in wall seconds for one unit of SCF work at one core.
_WALLTIME_BASE_S = 30.0

#: SCF cost scales super-linearly with the number of basis functions.
_BF_EXP = 1.4

#: Memory in MB per basis function, above the floor set by maxcore.
_MEMORY_MB_PER_BF = 40.0

#: Rough per-correlated-method cost multiplier over a plain DFT/HF single
#: reference run.  DLPNO-CCSD(T) and full CC are deliberately the heaviest.
_METHOD_SCALE: Dict[str, float] = {
    "hf": 1.0,
    "dft": 1.0,
    "pbe0": 1.0,
    "b3lyp": 1.0,
    "tpss": 1.0,
    "mp2": 8.0,
    "ri-mp2": 8.0,
    "ri-sos-mp2": 6.0,
    "ccsd": 30.0,
    "ccsd(t)": 45.0,
    "dlpno-ccsd": 30.0,
    "dlpno-ccsd(t)": 45.0,
    "cris": 45.0,
}

#: Recorded actual/estimate core-hour ratios per method family, from EF1's
#: finished jobs.  ``record_outcome`` appends here; ``estimate`` folds the mean
#: in.  In-memory per session (the budget itself is per session).
_LEARNED: Dict[str, list[float]] = {}

#: How far a learned factor may move an estimate, so a single outlier record
#: cannot swing it wildly.  A factor out of this band is clamped to the edge.
_SCALE_CLAMP = (0.25, 4.0)


def _method_family(method: str) -> str:
    """A stable cost family for a method string (case- and punctuation-free).

    Unknown methods map to the DFT family rather than failing, but a caller
    can see that the model treated them as cheap in :attr:`Estimate.source`.
    """
    if not method:
        return "dft"
    m = method.strip().lower().replace("-", "").replace("_", "").replace("/", "")
    if "dlpno" in m:
        suffix = "ccsd" if "(" not in m else "ccsd(t)"
        return suffix
    if "ccsd" in m or "ccsdt" in m:
        return "ccsd(t)" if "(" in m else "ccsd"
    if "mp2" in m:
        return "mp2"
    if "cris" in m:
        return "cris"
    return "dft"


def learned_scale(family: str) -> Optional[float]:
    """Mean recorded actual/estimate ratio for ``family``, clamped.

    ``None`` when no EF1 record exists for the family yet -- the estimate is
    the plain size-based model.  The mean is clamped to :data:`_SCALE_CLAMP`
    so one runaway record cannot silently triple or quarter every future job.
    """
    ratios = _LEARNED.get(family)
    if not ratios:
        return None
    mean = sum(ratios) / len(ratios)
    return min(max(mean, _SCALE_CLAMP[0]), _SCALE_CLAMP[1])


def record_outcome(job_report: Dict, est: Estimate) -> None:
    """Record estimate vs actual (EF1) after a job so the model learns.

    ``job_report`` is one record from :meth:`cluster_efficiency.account` (the
    real writer), carrying ``core_hours_used``.  An unmeasured, zero or
    negative actual is skipped -- never treated as zero -- and a family with
    no known methods maps to the DFT model, exactly as ``estimate`` does.
    """
    actual = job_report.get("core_hours_used")
    if not actual or actual <= 0 or est.core_hours <= 0:
        return
    family = _method_family(est.method)
    _LEARNED.setdefault(family, []).append(actual / est.core_hours)


@dataclass
class Estimate:
    """What a job is expected to cost, and how long it will wait to start.

    ``wait_s`` is a positive number of seconds when the scheduler's own
    ``--test-only`` answer stated a start delay, and :data:`UNKNOWN` when the
    gate was closed (:data:`scheduler_client.THROTTLED`) or the answer did not
    state a delay.  The estimate never invents a wait.
    """

    core_hours: float
    memory_mb: int
    walltime_s: int
    wait_s: object  # int seconds, or "unknown"
    source: str
    method: str
    n_atoms: int
    basis_functions: int
    partition: str
    pal: int


def _confirm_client(client: Optional[SchedulerClient]) -> SchedulerClient:
    """A caller-supplied client, else the process-wide throttled singleton.

    Routing through :data:`scheduler_client.default_client` is what every
    other DELFIN module does, so back-to-back estimates share ONE throttle
    window and cache -- the "no per-loop / per-widget scheduler ask" rule.
    """
    return client if client is not None else default_client


def _queue_wait(client: SchedulerClient, partition: str, pal: int) -> object:
    """Ask the scheduler's own ``--test-only`` once; THROTTLED -> "unknown"."""
    argv: Sequence[str] = [
        SBATCH,
        "--test-only",
        f"--partition={partition}",
        f"--cpus-per-task={pal}",
        "submit.sh",
    ]
    answer = client.query(kind_of_argv(list(argv)), list(argv))
    if answer is THROTTLED or answer is None:
        return UNKNOWN
    return _start_seconds_from(answer)


def _start_seconds_from(answer) -> object:
    """Seconds-until-start from a real ``sbatch --test-only`` answer.

    Real SLURM accepts print ``to start at <ISO or space-separated datetime>``
    (the marker ``slurm_submit.partition_accepts`` keys on).  Returns the
    whole seconds until that start (rounded, so sub-second clock drift does
    not break a caller), ``UNKNOWN`` when no start datetime is stated or it
    already passed.  A throttled/absent answer is handled by the caller.
    """
    match = _START_AT_RE.search(str(answer or ""))
    if not match:
        return UNKNOWN
    raw = match.group(1).replace(" ", "T")
    try:
        start = datetime.fromisoformat(raw)
    except ValueError:
        return UNKNOWN
    wait = round((start - datetime.now()).total_seconds())
    return wait if wait >= 0 else UNKNOWN


def estimate(job: Dict, *, client: Optional[SchedulerClient] = None,
             run_fn=None) -> Estimate:
    """Cost-and-time estimate plus expected queue wait for ``job``.

    ``job`` must carry at least ``n_atoms``, ``basis_functions`` and ``method``
    (``partition``, ``pal`` and ``maxcore`` are read when present).  The model
    is a transparent scaling: SCF wall time grows super-linearly with basis
    functions and is multiplied by a method cost factor; core-hours are wall
    time times the requested parallelism; memory grows with basis functions
    above a maxcore floor.  The queue wait comes from the scheduler's own
    ``--test-only`` answer through the throttled client.
    """
    n_atoms = int(job.get("n_atoms") or 0)
    basis = int(job.get("basis_functions") or 0)
    method = str(job.get("method") or "")
    partition = str(job.get("partition") or "cpu")
    pal = int(job.get("pal") or 1)
    maxcore_mb = int(job.get("maxcore") or 0)

    family = _method_family(method)
    method_scale = _METHOD_SCALE.get(family, 1.0)

    # SCF wall time in seconds at one core, then spread over the requested
    # cores (parallel speed-up is sub-linear: sqrt of the core count).
    scf_units = max(1, n_atoms) * (max(1, basis) ** _BF_EXP)
    walltime_s = int(
        round(_WALLTIME_BASE_S * scf_units * method_scale / math.sqrt(max(1, pal)))
    )
    core_hours = walltime_s * pal / 3600.0

    # Fold in the EF1-learned factor for this method family, where one has
    # been recorded; ``source`` states it so the number is traceable.
    learned = learned_scale(family)
    if learned is not None:
        core_hours = core_hours * learned

    memory_mb = int(round(max(maxcore_mb, basis * _MEMORY_MB_PER_BF)))

    client = _confirm_client(client)
    wait_s = _queue_wait(client, partition, pal)

    if learned is None:
        source = f"size-based model, family={family}"
    else:
        n_rec = len(_LEARNED.get(family, []))
        source = (f"size-based model, family={family}, "
                  f"learned x{learned:.2f} from {n_rec} EF1 record(s)")

    return Estimate(
        core_hours=core_hours,
        memory_mb=memory_mb,
        walltime_s=walltime_s,
        wait_s=wait_s,
        source=source,
        method=method,
        n_atoms=n_atoms,
        basis_functions=basis,
        partition=partition,
        pal=pal,
    )


@dataclass
class Budget:
    """A per-session core-hour budget with a running total.

    The decision rule is fixed: within budget -> ``go``; over budget -> the
    person is asked (``ask``).  The agent never shrinks PAL/maxcore/wall time
    silently to make a job fit.
    """

    limit_core_hours: float
    spent_core_hours: float = 0.0

    def remaining(self) -> float:
        return self.limit_core_hours - self.spent_core_hours

    def decision(self, est: Estimate):
        """Return ``(verdict, reason)``: ``go`` in budget, ``ask`` when over."""
        remaining = self.remaining()
        if remaining >= est.core_hours:
            reason = (
                f"estimate {est.core_hours} core-hours is within the remaining "
                f"{remaining} core-hours"
            )
            return ("go", reason)
        reason = (
            f"estimate {est.core_hours} core-hours exceeds the remaining "
            f"{remaining} core-hours"
        )
        return ("ask", reason)

    def record(self, est: Estimate) -> None:
        """Accumulate a job's estimate into the spent total."""
        self.spent_core_hours += est.core_hours
