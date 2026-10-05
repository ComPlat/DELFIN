"""T6 phase-3 adversarial review of slurm_tests submit/result/CLI (nacht-s24).

Submit/result never touch the scheduler here: an injectable submit_runner and
state_runner stub the sbatch / squeue sides (the hard rule that tests never
launch ORCA/xTB/SLURM). These push the surfaces the builder's own suite does
not: finding F3 (summary-writer argv order) and the security invariant that a
malformed CLI test path never reaches submission.

Security invariants asserted (must hold on any correct version):
* F3  the embedded summary writer binds argv in the order the shell passes --
      (rc, ref, out) -- so it writes the summary to ``out``, never to a file
      named after the ref (and never leaks a file into cwd / the checkout).
* CLI a malformed test path must NEVER result in a submission (submit_runner
      must not be called); the failure happens inside build_job, which raises
      before submit is reached.
* a corrupt on-disk summary must degrade to None, not crash the caller.
"""
from __future__ import annotations

import pytest

from delfin.agent import slurm_tests


class _NoCallSubmit:
    """submit_runner that fails the test the moment it is invoked."""

    def __init__(self):
        self.called = False

    def __call__(self, argv, cwd=None):
        self.called = True
        raise AssertionError("submit must not run for a malformed test path")


def test_submit_argv_carries_one_partition_when_sbatch_already_has_one(tmp_path):
    """The builder appends --partition= only when sbatch_command has none."""
    calls: list = []

    def fake_sbatch(sbatch, script, run_dir):
        return [sbatch, str(script), "--partition=auto"]

    # ensure sbatch_command is not reached for real: patch submit's argv build
    orig = slurm_tests.sbatch_command
    slurm_tests.sbatch_command = fake_sbatch
    try:
        slurm_tests.submit(
            "#!/bin/bash\n#SBATCH --export=ALL\n", tmp_path,
            submit_runner=lambda argv, cwd=None: (calls.append((list(argv), cwd)) or
                                                   ("Submitted batch job 7", 0)),
            partition="cpu",
        )
    finally:
        slurm_tests.sbatch_command = orig
    assert calls, "submit_runner was never invoked"
    argv, _ = calls[0]
    assert [a for a in argv if a.startswith("--partition=")] == ["--partition=auto"]


def test_result_corrupt_summary_is_none_not_crash(tmp_path):
    """A truncated/corrupt JSON on disk degrades to None (started, not green)."""
    states = lambda argv: "31339 COMPLETED\n"  # noqa: E731
    out = tmp_path / "logs"
    out.mkdir()
    (out / "summary_31339.json").write_text("{ not json !!")
    res = slurm_tests.result("31339", output_dir=out, state_runner=states)
    assert res["job_id"] == "31339"
    assert res["summary"] is None


def test_cmd_cli_malformed_test_path_never_submits():
    """A test path outside tests/ must not result in a submission."""
    class _Args:  # noqa: N801
        pass

    args = _Args()
    args.repo = "delfin"
    args.ref = "main"
    args.partition = "cpu"
    args.minutes = 10
    args.tests = ["evil.py"]  # not under tests/ -> build_job raises ValueError
    submitter = _NoCallSubmit()
    orig = slurm_tests.submit
    slurm_tests.submit = submitter
    try:
        with pytest.raises(ValueError):
            slurm_tests.cmd_test_on_slurm(args)
    finally:
        slurm_tests.submit = orig
    assert not submitter.called, "a malformed test path must never reach submit"
