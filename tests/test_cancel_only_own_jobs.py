"""A job cancel reaches only the jobs this session submitted, unasked.

Input: ``scancel`` / ``qdel`` / ``bkill`` commands through the real bash
gate, with the session's ledger of submitted job ids
(``KitToolPermissions.submitted_job_ids``) filled the way the executor
fills it: from the output of a command that ran a submit program
(``reversibility.submitted_job_ids``).

Output: a cancel whose every target is a literal id in the ledger passes
this rule; any other cancel -- a foreign id, an id from a variable or a
substitution, a cancel by name or partition with no id, ``xargs scancel``,
``bkill 0`` -- is refused where nobody can be asked (naming ``squeue`` and
the report as the way on), and asked in every mode where someone can be,
bypassPermissions included.
"""
from __future__ import annotations

import pytest

from delfin.agent import api_client as A
from delfin.agent import reversibility as R


@pytest.fixture(autouse=True)
def _classifier_off(monkeypatch):
    monkeypatch.setattr(A, "_classifier_settings", lambda: None)


class Dialog:
    last_timed_out = False

    def __init__(self, answer=True):
        self.answer, self.shown = answer, []

    def ask(self, name, args, preview):
        self.shown.append(preview)
        return self.answer


def _perms(tmp_path, mode="bypassPermissions", dialog=None, own=()):
    ws = tmp_path / "wt"
    ws.mkdir(exist_ok=True)
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode=mode)
    perms.submitted_job_ids.update(own)
    if dialog is not None:
        perms.confirm_callback = dialog.ask
    return perms


def _gate(perms, cmd):
    return A._DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


# -- the ledger ---------------------------------------------------------------

@pytest.mark.parametrize("cmd,out,ids", [
    ("sbatch run.sh", "Submitted batch job 4711\n", ["4711"]),
    ("cd calc && sbatch run.sh", "Submitted batch job 4711 on cluster c1\n",
     ["4711"]),
    ("sbatch --parsable run.sh", "4712\n", ["4712"]),
    ("sbatch --parsable run.sh", "4713;cluster1\n", ["4713"]),
    ("qsub job.pbs", "8123.pbs-server\n", ["8123"]),
    ("bsub < job.lsf", "Job <991> is submitted to queue <normal>.\n", ["991"]),
])
def test_a_submit_records_its_job_id(cmd, out, ids):
    assert R.submitted_job_ids(cmd, out) == ids


@pytest.mark.parametrize("cmd,out", [
    ("echo 'Submitted batch job 4711'", "Submitted batch job 4711\n"),
    ("cat old.log", "Submitted batch job 4711\n"),
    ("python count.py", "4711\n"),
])
def test_output_without_a_submit_program_records_nothing(cmd, out):
    assert R.submitted_job_ids(cmd, out) == []


# -- the targets --------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "scancel 4711",
    "scancel 4711_3",
    "scancel 4711_[1-3]",
    "scancel 4711.0",
    "scancel --signal=USR1 4711",
    "scancel -s KILL 4711",
    "scancel -n relax 4711",
    "qdel 4711.pbs-server",
    "bkill 4711",
    "squeue -u $USER && scancel 4711",
])
def test_a_cancel_of_own_ids_is_free(cmd):
    assert R.foreign_job_cancels(cmd, {"4711"}) == [], cmd


@pytest.mark.parametrize("cmd", [
    "scancel 4712",
    "scancel 4711 4712",
    "scancel $JOB",
    "scancel $(squeue -h -o %i -u $USER)",
    "squeue -h -o %i | xargs scancel",
    "scancel -n relax",
    "scancel --partition=short",
    "scancel --frobnicate 4711",
    "qdel all",
    "bkill 0",
    "bkill -u alice 4711",
    "bash -c 'scancel 4712'",
    "timeout 10 scancel 4712",
])
def test_any_other_cancel_is_foreign(cmd):
    assert R.foreign_job_cancels(cmd, {"4711"}), cmd


def test_an_empty_ledger_makes_every_cancel_foreign():
    assert R.foreign_job_cancels("scancel 4711", set())


# -- through the gate ---------------------------------------------------------

def test_own_job_cancel_runs_in_bypass_unattended(tmp_path):
    assert _gate(_perms(tmp_path, own={"4711"}), "scancel 4711") is None


def test_foreign_job_cancel_is_refused_unattended_with_a_way_on(tmp_path):
    out = _gate(_perms(tmp_path, own={"4711"}), "scancel 4712")
    assert out is not None
    assert "4712" in out and "squeue" in out and "report" in out


def test_foreign_job_cancel_is_asked_in_bypass_attended(tmp_path):
    dialog = Dialog(True)
    perms = _perms(tmp_path, dialog=dialog)
    assert _gate(perms, "scancel 4712") is None
    assert len(dialog.shown) == 1 and "4712" in dialog.shown[0]


def test_a_no_to_the_foreign_cancel_refuses_it(tmp_path):
    dialog = Dialog(False)
    out = _gate(_perms(tmp_path, dialog=dialog), "scancel 4712")
    assert out is not None and len(dialog.shown) == 1


def test_default_mode_asks_once_for_a_foreign_cancel(tmp_path):
    # One dialog, not this rule's and then the ordinary one.
    dialog = Dialog(True)
    perms = _perms(tmp_path, mode="default", dialog=dialog)
    assert _gate(perms, "scancel 4712") is None
    assert len(dialog.shown) == 1


def test_the_executor_fills_the_ledger_from_a_submit(tmp_path):
    perms = _perms(tmp_path)
    A._record_submitted_jobs(perms, "sbatch run.sh",
                             "Submitted batch job 4711\n")
    assert "4711" in perms.submitted_job_ids
    assert _gate(perms, "scancel 4711") is None


def test_sub_agents_share_the_ledger(tmp_path):
    import dataclasses
    perms = _perms(tmp_path)
    child = dataclasses.replace(perms)
    A._record_submitted_jobs(child, "sbatch run.sh",
                             "Submitted batch job 4711\n")
    assert "4711" in perms.submitted_job_ids
