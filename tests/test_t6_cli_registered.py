"""`delfin-agent test-on-slurm` is registered and refuses hostile input.

The parser builds without importing slurm_tests (imported on use), and the
command hands its arguments to build_job's allow-lists: a newline in the
partition is refused before anything is written or submitted.
"""
from __future__ import annotations

import pytest

from delfin.agent import cli


def test_the_subcommand_parses():
    args = cli.build_parser().parse_args(
        ["test-on-slurm", "/repo", "main", "cpu", "10", "tests/test_x.py"])
    assert args.partition == "cpu" and args.minutes == 10
    assert args.tests == ["tests/test_x.py"]


def test_a_newline_in_the_partition_never_reaches_a_submission(tmp_path, monkeypatch):
    from delfin.agent import slurm_tests
    submitted = []
    monkeypatch.setattr(slurm_tests, "submit", lambda *a, **k: submitted.append(a) or "1")
    args = cli.build_parser().parse_args(
        ["test-on-slurm", str(tmp_path), "main", "cpu\n#SBATCH --export=ALL,E=1",
         "10", "tests/test_x.py", "--run-dir", str(tmp_path)])
    with pytest.raises((ValueError, SystemExit)):
        rc = args.func(args)
        if rc:
            raise SystemExit(rc)
    assert submitted == []
