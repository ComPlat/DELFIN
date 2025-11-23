# delfin/job_priority.py
# -*- coding: utf-8 -*-
"""
Simple bottleneck detection for job prioritization.

Analyzes job dependency graph and boosts priority for jobs
that block many downstream jobs (bottlenecks).
"""

from typing import Dict, Set, Any
from delfin.common.logging import get_logger

logger = get_logger(__name__)


def count_downstream_jobs(job_id: str, all_jobs: Dict[str, Any]) -> int:
    """
    Count how many jobs (directly or indirectly) depend on this job.

    Args:
        job_id: Job to analyze
        all_jobs: Dict of job_id -> WorkflowJob

    Returns:
        Number of downstream jobs
    """
    count = 0
    visited = set()
    queue = [job_id]

    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)

        # Find all jobs that depend on current
        for other_id, other_job in all_jobs.items():
            if not hasattr(other_job, 'dependencies'):
                continue
            if current in other_job.dependencies and other_id not in visited:
                count += 1
                queue.append(other_id)

    return count


def is_exclusive_bottleneck(job_id: str, all_jobs: Dict[str, Any], pending_jobs: Set[str]) -> bool:
    """
    Check if a job is an exclusive bottleneck (all other pending jobs depend on it).

    Args:
        job_id: Job to check
        all_jobs: Dict of job_id -> WorkflowJob
        pending_jobs: Set of pending job IDs

    Returns:
        True if all other pending jobs depend on this job
    """
    if job_id not in pending_jobs:
        return False

    downstream = set()
    visited = set()
    queue = [job_id]

    while queue:
        current = queue.pop(0)
        if current in visited:
            continue
        visited.add(current)

        for other_id, other_job in all_jobs.items():
            if not hasattr(other_job, 'dependencies'):
                continue
            if current in other_job.dependencies and other_id not in visited:
                downstream.add(other_id)
                queue.append(other_id)

    # Check if all other pending jobs are downstream
    other_pending = pending_jobs - {job_id}
    if not other_pending:
        return False  # No other pending jobs

    # True if all other pending jobs depend on this one
    return other_pending.issubset(downstream)


def adjust_job_priorities(all_jobs: Dict[str, Any], bottleneck_threshold: int = 3) -> int:
    """
    Adjust job priorities based on downstream impact.

    Jobs with many downstream dependencies get HIGH priority
    so they are scheduled first and don't block other jobs.

    Args:
        all_jobs: Dict of job_id -> WorkflowJob
        bottleneck_threshold: Jobs with >= this many downstream jobs get HIGH priority

    Returns:
        Number of jobs that got priority boost
    """
    from delfin.dynamic_pool import JobPriority

    boosted = 0

    for job_id, job in all_jobs.items():
        if not hasattr(job, 'priority'):
            continue

        # Skip jobs that already have high priority
        if job.priority == JobPriority.HIGH:
            continue

        # Count downstream impact
        downstream_count = count_downstream_jobs(job_id, all_jobs)

        # Boost priority if this is a bottleneck
        if downstream_count >= bottleneck_threshold:
            job.priority = JobPriority.HIGH
            boosted += 1
            logger.info(
                f"[priority] Job {job_id} boosted to HIGH priority "
                f"({downstream_count} downstream jobs depend on it)"
            )

    return boosted
