"""Global ORCA job scheduler built on top of the shared workflow manager."""

from __future__ import annotations

from typing import Iterable, Optional, Set

from delfin.common.logging import get_logger
from delfin.parallel_classic_manually import (
    WorkflowJob,
    WorkflowRunResult,
    _WorkflowManager,
)

logger = get_logger(__name__)


class GlobalOrcaScheduler:
    """Aggregate ORCA jobs from multiple workflow phases into a single run."""

    def __init__(self, config: dict, *, label: str = "global_orca", max_jobs_override: Optional[int] = None):
        self._manager = _WorkflowManager(config, label=label, max_jobs_override=max_jobs_override)
        self._jobs_added = False
        self._plan_logged = False

    @property
    def manager(self) -> _WorkflowManager:
        return self._manager

    def add_jobs(self, jobs: Iterable[WorkflowJob]) -> None:
        added = 0
        for job in jobs:
            self._manager.add_job(job)
            added += 1
        if added:
            self._jobs_added = True
            logger.debug("[%s] registered %d jobs", self._manager.label, added)

    def add_job(self, job: WorkflowJob) -> None:
        self.add_jobs([job])

    def run(self) -> WorkflowRunResult:
        # Check both _jobs_added flag AND if manager actually has jobs
        # (jobs may have been added directly to manager bypassing add_job/add_jobs)
        if not self._jobs_added and not self._manager.has_jobs():
            logger.info("[%s] no ORCA jobs queued – skipping run", self._manager.label)
            return WorkflowRunResult()

        if not self._plan_logged:
            self._log_job_graph()
            self._plan_logged = True

        self._manager.run()

        completed: Set[str] = set(self._manager.completed_jobs)
        failed = dict(self._manager.failed_jobs)
        skipped = {
            job_id: list(deps)
            for job_id, deps in self._manager.skipped_jobs.items()
        }

        return WorkflowRunResult(
            completed=completed,
            failed=failed,
            skipped=skipped,
        )

    def shutdown(self) -> None:
        self._manager.shutdown()

    def _log_job_graph(self) -> None:
        jobs = getattr(self._manager, "_jobs", {})
        if not jobs:
            return
        logger.info("[%s] ORCA job graph (%d total):", self._manager.label, len(jobs))
        for job_id in sorted(jobs):
            job = jobs[job_id]
            deps = ", ".join(sorted(job.dependencies)) if job.dependencies else "none"
            logger.info(
                "  • %s  | deps: %s  | cores min/opt/max: %d/%d/%d",
                job_id,
                deps,
                job.cores_min,
                job.cores_optimal,
                job.cores_max,
            )
