"""T6 phase 3: submit + result through injectable runners, never the scheduler.

``submit`` writes a rendered job into a run dir, builds the submit argv via
``delfin.slurm_submit.sbatch_command`` (allow-listed) and passes it to an
injectable runner, so a test can record the command and answer with a fake
"Submitted batch job N" without touching the scheduler. ``result`` reads the
job state through ``delfin.agent.job_monitor.query_job_states_detailed`` with
an injectable state runner, and reads the JSON summary the job's copy-back
wrote. The CLI handler ``cmd_test_on_slurm`` is exercised through its public
entry with the module's submit stubbed.

Each test stubs the scheduler access -- nothing here submits a real job.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import slurm_tests


class FakeSubmit:
    """Records the submit argv; answer with a canned "Submitted batch job N"."""

    def __init__(self, out: str = "Submitted batch job 424242", rc: int = 0):
        self.out = out
        self.rc = rc
        self.calls: list[tuple] = []

    def __call__(self, argv, cwd=None):
        self.calls.append((list(argv), cwd))
        return self.out, self.rc


class FakeStates:
    """Answers squeue/sacct-style lines for the throttled reader."""

    def __init__(self, line: str = ""):
        self.line = line
        self.calls: list[list] = []

    def __call__(self, argv):
        self.calls.append(list(argv))
        return self.line


def test_submit_writes_script_and_returns_job_id(tmp_path):
    fake = FakeSubmit()
    job_id = slurm_tests.submit(
        "#!/bin/bash\n#SBATCH --export=ALL\n", tmp_path,
        submit_runner=fake, partition="cpu")
    assert job_id == "424242"
    assert (tmp_path / "test_run.sbatch").exists()
    assert (tmp_path / "test_run.sbatch").read_text().startswith("#!/bin/bash")
    # the submit argv carried the script and the chosen partition
    argv = fake.calls[0][0]
    assert str(tmp_path / "test_run.sbatch") in argv
    assert "--partition=cpu" in argv
    # cwd handed to the runner is the run dir
    assert fake.calls[0][1] == str(tmp_path)


def test_submit_raises_when_no_job_id_in_output(tmp_path):
    fake = FakeSubmit(out="something without an id", rc=0)
    with pytest.raises(RuntimeError, match="no job id"):
        slurm_tests.submit("#!/bin/bash\n", tmp_path, submit_runner=fake)


def test_submit_raises_on_nonzero_submit(tmp_path):
    fake = FakeSubmit(out="sbatch failed", rc=1)
    with pytest.raises(RuntimeError, match="rc=1"):
        slurm_tests.submit("#!/bin/bash\n", tmp_path, submit_runner=fake)


def test_result_reads_state_and_summary(tmp_path):
    states = FakeStates("31337 COMPLETED\n")
    out = tmp_path / "logs"
    out.mkdir()
    (out / "summary_31337.json").write_text(
        json.dumps({"ref": "main", "rc": 0, "passed_tests": 5}))
    res = slurm_tests.result("31337", output_dir=out, state_runner=states)
    assert res["state"] == "COMPLETED"
    assert res["summary"]["passed_tests"] == 5
    assert res["summary"]["ref"] == "main"


def test_result_no_summary_when_not_present(tmp_path):
    states = FakeStates("31338 RUNNING\n")
    res = slurm_tests.result("31338", output_dir=tmp_path, state_runner=states)
    assert res["state"] == "RUNNING"
    assert res["summary"] is None


def test_result_distinguishes_unknown_from_throttle(tmp_path):
    # scheduler answered and does not know the id -> state ""
    states = FakeStates("")  # no line -> not known
    res = slurm_tests.result("99", output_dir=tmp_path, state_runner=states)
    assert res["state"] == ""


def test_build_job_copies_back_only_logs_to_output_dir():
    text = slurm_tests.build_job("delfin", "main", ["tests/x.py"], "cpu", 10,
                                 output_dir="out")
    assert "--output=out/test_%j.out" in text
    assert 'cp "$LOCAL/summary.json"' in text
    assert "summary_${SLURM_JOB_ID:-unknown}.json" in text
    assert 'cp "$LOCAL/pytest.log"' in text


class _Args:
    pass


def test_cmd_test_on_slurm_renders_and_submits(monkeypatch, capsys):
    called = {}

    def fake_submit(job_text, run_dir):
        called["job_text"] = job_text
        called["run_dir"] = run_dir
        return "50505"

    monkeypatch.setattr(slurm_tests, "submit", fake_submit)
    args = _Args()
    args.repo = "delfin"
    args.ref = "main"
    args.partition = "cpu"
    args.minutes = 10
    args.tests = ["tests/test_x.py"]
    rc = slurm_tests.cmd_test_on_slurm(args)
    assert rc == 0
    assert capsys.readouterr().out.strip() == "50505"
    assert "tests/test_x.py" in called["job_text"]
    assert "--export=ALL" in called["job_text"]


def test_cmd_test_on_slurm_rejects_missing_required_args(capsys):
    args = _Args()
    args.repo = ""
    args.ref = ""
    args.partition = ""
    args.minutes = 0
    args.tests = []
    rc = slurm_tests.cmd_test_on_slurm(args)
    assert rc == 2
    assert "usage" in capsys.readouterr().out
