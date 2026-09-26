"""A background job's status says how to read its output.

The status named the stdout/stderr paths on the host's scratch file
system, and three sessions in one supervised run (2026-09-26) tried to
read them with read_file -- outside the workspace, so each attempt asked
the operator to open a directory every user shares. bash_output reads
both streams; the status now says so.
"""

from __future__ import annotations

import sys

from delfin.agent import bash_jobs


def test_the_status_names_bash_output(tmp_path):
    reg = bash_jobs._Registry() if hasattr(bash_jobs, "_Registry") else None
    job = reg.start(f"{sys.executable} -c 'print(1)'", cwd=str(tmp_path),
                    timeout_s=30) if reg is not None else None
    assert job is not None
    status = job.status_dict()
    assert status["read_with"] == f"bash_output(job_id='{job.job_id}')"
