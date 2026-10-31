"""The config['run_dir'] fallback in OccupierWorkflow.run (F6 follow-up).

Reviewer nacht-s4, APPROVED b5ddddc7a, correctness mark 4: the run_dir
compatibility fallback has no test. This pins it and pins that an empty
config still calls run_OCCUPIER bare-equivalent (work_dir=None).
"""
from unittest import mock

from delfin.workflows.contrib.occupier_workflow import OccupierWorkflow


class TestRunDirFallback:
    def test_run_dir_in_config_is_used(self):
        wf = OccupierWorkflow()
        with mock.patch("delfin.occupier.run_OCCUPIER") as run_occ:
            wf.run(config={"run_dir": "/stage/occupier"})
        assert run_occ.call_args.kwargs.get("work_dir") == "/stage/occupier"

    def test_empty_config_keeps_work_dir_none(self):
        wf = OccupierWorkflow()
        with mock.patch("delfin.occupier.run_OCCUPIER") as run_occ:
            wf.run(config={})
        assert run_occ.call_args.kwargs.get("work_dir") is None

    def test_work_dir_key_wins_over_run_dir(self):
        wf = OccupierWorkflow()
        with mock.patch("delfin.occupier.run_OCCUPIER") as run_occ:
            wf.run(config={"work_dir": "/a", "run_dir": "/b"})
        assert run_occ.call_args.kwargs.get("work_dir") == "/a"
