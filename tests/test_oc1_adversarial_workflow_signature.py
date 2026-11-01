"""Adversarial tests: OccupierWorkflow.run must forward its arguments (F6).

Finding F6 (s3 audit, my verification): delfin/workflows/contrib/
occupier_workflow.py:16-19 defines run(self, *, config, **kwargs) but calls
run_OCCUPIER() with NO arguments — run_OCCUPIER takes work_dir
(delfin/occupier.py:624). A caller passing a work_dir or stage config gets
silently ignored arguments and OCCUPIER runs in the process CWD instead of
the intended stage folder. That is the parallel-execution race the
work_dir parameter was added to prevent (occupier.py:632-637 comment).
"""
import inspect
from unittest import mock

import pytest

from delfin.workflows.contrib.occupier_workflow import OccupierWorkflow


class TestOccupierWorkflowRunForwardsArguments:
    def test_run_forwards_config_to_work_dir(self):
        """run() must pass a work_dir/config through to run_OCCUPIER."""
        sig = inspect.signature(OccupierWorkflow.run)
        params = set(sig.parameters)
        assert params >= {"config"}, f"run() lost its config parameter: {sig}"

    def test_run_actually_forwards(self, tmp_path):
        """The forwarding must be real: with a config/work_dir given,
        run_OCCUPIER must be called with work_dir set, not bare."""
        wf = OccupierWorkflow()
        with mock.patch("delfin.occupier.run_OCCUPIER") as run_occ:
            wf.run(config={"work_dir": str(tmp_path)})
        assert run_occ.call_args is not None, "run_OCCUPIER was not called"
        assert run_occ.call_args.kwargs.get("work_dir") or (
            run_occ.call_args.args and run_occ.call_args.args[0] is not None
        ), (
            f"run_OCCUPIER called with {run_occ.call_args} — config/work_dir "
            "is silently dropped (F6)"
        )
