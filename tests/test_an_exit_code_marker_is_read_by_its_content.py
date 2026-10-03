"""A finished job's exit code is what its marker file says, not its name.

local_runner writes ``.exit_code_<slurm job id>`` with the process exit
code as the file's content (dashboard/local_runner.py), and the dashboard
backend reads it that way. The calc index and the agent's status tool read
the number from the file NAME instead — the SLURM job id — so ten clean
ORCA runs of the pKa validation were reported as
"failed (exit code 7435268)".
"""
from __future__ import annotations

import os
import time

from delfin.doc_server.calc_indexer import completion_of, outcome_of_folder


def _marker(folder, job_id, content):
    path = folder / f".exit_code_{job_id}"
    path.write_text(content, encoding="utf-8")
    return path


def test_a_clean_run_is_not_reported_as_failed_with_its_job_id(tmp_path):
    _marker(tmp_path, 7435268, "0\n")
    assert completion_of(tmp_path) == (True, 0)
    assert "7435268" not in outcome_of_folder(tmp_path)


def test_a_failed_run_carries_its_real_code(tmp_path):
    _marker(tmp_path, 7435269, "137\n")
    assert completion_of(tmp_path) == (True, 137)


def test_the_newest_marker_speaks_after_a_resubmit(tmp_path):
    old = _marker(tmp_path, 7435000, "137\n")
    past = time.time() - 3600
    os.utime(old, (past, past))
    _marker(tmp_path, 7435999, "0\n")
    assert completion_of(tmp_path) == (True, 0)


def test_an_empty_marker_still_names_its_code_in_the_name(tmp_path):
    """Older folders (and the benchmark fixtures) put the code into the
    name and leave the file empty; those still read as before."""
    _marker(tmp_path, 1025, "")
    assert completion_of(tmp_path) == (True, 1025)


def test_the_status_tool_reports_the_content_too(tmp_path):
    from delfin.api import calculation_status

    _marker(tmp_path, 7435268, "0\n")
    status = calculation_status(str(tmp_path))
    says = [e["says"] for e in status.evidence if e["source"].startswith(".exit_code_")]
    assert says == ["exit code 0"], says
