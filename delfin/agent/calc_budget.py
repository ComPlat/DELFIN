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
from dataclasses import dataclass
from typing import Dict, Optional, Sequence

from delfin.scheduler_client import SchedulerClient, THROTTLED

#: Sentinel for a wait we could not learn from the scheduler's own answer.
UNKNOWN = "unknown"

#: The sbatch command the estimate asks through the scheduler client.
SBATCH = "sbatch"

#: Form of the queue-wait answer DELFIN's own ``--test-only`` caller emits;
#: the estimate parses a ``started-after-<n>s`` delay out of it.
_WAIT_RE = re.compile(r"started-after-(\d+)s")

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
    if client is not None:
        return client
    return SchedulerClient()


def _queue_wait(client: SchedulerClient, partition: str, pal: int) -> object:
    """Ask the scheduler's own ``--test-only`` once; THROTTLED -> "unknown"."""
    argv: Sequence[str] = [
        SBATCH,
        "--test-only",
        f"--partition={partition}",
        f"--cpus-per-task={pal}",
        "submit.sh",
    ]
    answer = client.query("sbatch", list(argv))
    if answer is THROTTLED or answer is None:
        return UNKNOWN
    match = _WAIT_RE.search(str(answer))
    if not match:
        return UNKNOWN
    return int(match.group(1))


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

    memory_mb = int(round(max(maxcore_mb, basis * _MEMORY_MB_PER_BF)))

    client = _confirm_client(client)
    wait_s = _queue_wait(client, partition, pal)

    return Estimate(
        core_hours=core_hours,
        memory_mb=memory_mb,
        walltime_s=walltime_s,
        wait_s=wait_s,
        source=f"size-based model, family={family}",
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
