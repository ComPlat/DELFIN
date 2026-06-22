"""Regression tests for the classic workflow engine's inline job path.

Guards against the NameError in `_WorkflowManager._submit` that aborted
entire pipelines (all 5 of Jan's OCCUPIER runs in zzz_Bugs/) at the
best-FoB selection step. Root cause: the inline branch called
`self._record_duration(job, duration, cores)` but `cores` was only
defined in the sibling `runner` closure for non-inline jobs.
"""

import threading
import time
from pathlib import Path

import pytest

from delfin.workflows.engine.classic import WorkflowJob, _WorkflowManager


_SOURCE_FILE = (
    Path(__file__).resolve().parents[1]
    / "delfin"
    / "workflows"
    / "engine"
    / "classic.py"
)


def _make_bare_manager():
    """Construct a _WorkflowManager without running __init__.

    We only need the attributes that _submit touches on the inline
    path: label, _lock, _job_start_times, plus the _record_duration,
    _mark_completed, _mark_failed methods (which live on the class).
    Setting up a full manager would require a global job manager
    singleton; that's overkill for exercising the inline branch.
    """
    manager = object.__new__(_WorkflowManager)
    manager.label = "test_inline"
    manager._lock = threading.RLock()
    manager._job_start_times = {}
    manager._completed = set()
    manager._failed = {}
    manager._inflight = set()
    manager._event = threading.Event()
    manager._completion_listeners = []
    manager._parent_manager = None
    return manager


def test_inline_job_submit_does_not_raise_nameerror():
    """Executing _submit on an inline job must not crash with NameError.

    Before the fix, this raised:
        NameError: name 'cores' is not defined
    and aborted the entire pipeline run.
    """
    manager = _make_bare_manager()

    ran = []

    def work_fn(cores_arg):
        # Inline jobs get 0 cores passed by design (zero-cost)
        ran.append(cores_arg)

    job = WorkflowJob(
        job_id="inline_test",
        work=work_fn,
        description="inline test job",
        inline=True,
    )

    # Must not raise NameError or any other exception
    manager._submit(job)

    assert ran == [0], "inline job should be invoked with cores=0"
    assert "inline_test" in manager._completed, \
        "inline job should be marked completed"


def test_inline_job_records_zero_cores_without_raising():
    """_record_duration must accept 0 for inline jobs (zero-cost accounting)."""
    manager = _make_bare_manager()

    def quick_work(_cores):
        time.sleep(0.001)

    job = WorkflowJob(
        job_id="inline_record",
        work=quick_work,
        description="records zero cores",
        inline=True,
    )

    # Should complete cleanly; duration is recorded in
    # JOB_DURATION_HISTORY with cores_used=0 (no STAGE_CORE_HISTORY entry).
    manager._submit(job)

    assert "inline_record" in manager._completed


def test_source_inline_branch_does_not_reference_undefined_cores():
    """Static guard: the inline branch of _submit must not use `cores`.

    `cores` is only defined inside the `runner` closure for non-inline
    jobs. Using it in the inline branch (line ~946) caused NameError.
    This check reads the source and verifies the inline branch uses an
    explicit integer (0) instead.
    """
    src = _SOURCE_FILE.read_text()

    # Locate _submit method and limit to the inline branch (before the
    # `# Check if this is an exclusive bottleneck` comment that starts
    # the non-inline path).
    submit_idx = src.index("def _submit(")
    bottleneck_idx = src.index(
        "# Check if this is an exclusive bottleneck",
        submit_idx,
    )
    inline_branch = src[submit_idx:bottleneck_idx]

    assert "if job.inline:" in inline_branch, \
        "sanity check: inline branch still exists"
    assert "_record_duration(job, duration, cores)" not in inline_branch, (
        "`cores` is undefined in the inline branch — use 0 for zero-cost "
        "inline accounting. This previously aborted the entire pipeline "
        "at initial_fob_best in the OCCUPIER phase."
    )
