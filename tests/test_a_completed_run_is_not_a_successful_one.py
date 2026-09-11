"""A calculation record said `completed: true` and `exit_code: 1025`.

Driven through get_calc_info: `completed` means the run ENDED -- an
exit-code file exists -- and says nothing about how. A reader takes
"completed" for "succeeded", which is the wrong kind of wrong for a
tool the agent grounds answers on. The record now carries `outcome`,
a phrase that says which it was, and both result shapes surface it.
"""

from __future__ import annotations

from delfin.doc_server import calc_indexer as ci
from delfin.doc_server import calc_search as cs


def test_outcome_says_succeeded_only_for_exit_code_zero():
    assert ci._outcome_of(True, 0) == "succeeded (exit code 0)"


def test_a_non_zero_exit_code_is_named_as_a_failure():
    out = ci._outcome_of(True, 1025)
    assert out.startswith("failed") and "1025" in out


def test_an_unreadable_exit_code_is_not_called_a_success():
    out = ci._outcome_of(True, "unknown")
    assert "succeeded" not in out and "unreadable" in out


def test_a_run_without_an_exit_code_is_running_or_crashed():
    assert "running or crashed" in ci._outcome_of(False, None)


def test_no_evidence_is_unknown():
    assert ci._outcome_of(None, None).startswith("unknown")


def test_both_result_shapes_carry_the_outcome():
    rec = {"calc_id": "x", "completed": True, "exit_code": 1025,
           "outcome": ci._outcome_of(True, 1025), "smiles": "", "modules": []}
    detailed = cs._format_result_detailed(rec)
    assert detailed["outcome"].startswith("failed") and detailed["completed"] is True
    brief = cs._format_result(rec)
    assert brief["outcome"].startswith("failed") and brief["completed"] is True


def test_the_indexer_writes_the_outcome_from_the_exit_code_file(tmp_path):
    from pathlib import Path

    root = tmp_path / "calc"
    for name, code in (("good", "0"), ("bad", "1025")):
        d = root / name
        d.mkdir(parents=True)
        (d / "CONTROL.txt").write_text("functional = PBE0\n")
        (d / "run.out").write_text("ORCA\n")
        (d / f".exit_code_{code}").write_text("")
    running = root / "running"
    running.mkdir()
    (running / "CONTROL.txt").write_text("functional = PBE0\n")
    (running / "delfin_run.log").write_text("started\n")

    recs = {r["rel_path"].split("/")[-1]: r for r in ci._scan_calc_dir(Path(root), "calc", quiet=True)}
    assert recs["good"]["completed"] is True and recs["good"]["outcome"] == "succeeded (exit code 0)"
    assert recs["bad"]["completed"] is True and recs["bad"]["outcome"] == "failed (exit code 1025)"
    assert recs["running"]["completed"] is False and "running or crashed" in recs["running"]["outcome"]


# ---------------------------------------------------------------------------
# An output with no bookkeeping still says how the run ended
# ---------------------------------------------------------------------------
#
# An archived folder holding only run.inp and run.out -- copied without
# its state file, or computed by hand -- read as "unknown (no run log and
# no exit code)" although the output's last lines say TERMINATED NORMALLY.
# ORCA's own termination line is evidence; it is read from the tail.

_NORMAL = "\n".join(["* O R C A *", "FINAL SINGLE POINT ENERGY -113.30", "",
                     "                             ****ORCA TERMINATED NORMALLY****",
                     "TOTAL RUN TIME: 0 days 0 hours 1 minutes 2 seconds 345 msec", ""])
_ERROR = "\n".join(["* O R C A *", "Error: SCF NOT CONVERGED", "",
                    "ORCA finished by error termination in SCF", ""])
_CUT = "\n".join(["* O R C A *", "SCF ITERATIONS", "ITER  Energy", "  0  -113.1", ""])


def _bare(tmp_path, name, out_text):
    d = tmp_path / name
    d.mkdir()
    (d / "run.inp").write_text("! PBE0 def2-SVP SP\n")
    (d / "run.out").write_text(out_text)
    return d


def test_a_normal_termination_is_finished_per_the_output(tmp_path):
    from delfin.doc_server.calc_indexer import outcome_of_folder
    assert outcome_of_folder(_bare(tmp_path, "n", _NORMAL)).startswith("finished per ORCA output")


def test_an_error_termination_is_failed_per_the_output(tmp_path):
    from delfin.doc_server.calc_indexer import outcome_of_folder
    assert outcome_of_folder(_bare(tmp_path, "e", _ERROR)).startswith("failed per ORCA output")


def test_an_output_without_a_termination_line_is_running_or_crashed(tmp_path):
    from delfin.doc_server.calc_indexer import outcome_of_folder
    assert outcome_of_folder(_bare(tmp_path, "c", _CUT)).startswith("running or crashed")


def test_the_termination_line_is_read_from_the_tail_of_a_long_output(tmp_path):
    from delfin.doc_server.calc_indexer import outcome_of_folder
    long = "SCF ITERATIONS\n" * 200000 + _NORMAL      # ~3 MB before the line
    assert outcome_of_folder(_bare(tmp_path, "l", long)).startswith("finished per ORCA output")


def test_an_exit_code_file_still_outranks_the_output(tmp_path):
    """The exit code is the run's own verdict; the output's line is the
    fallback for folders that have nothing else."""
    from delfin.doc_server.calc_indexer import outcome_of_folder
    d = _bare(tmp_path, "x", _NORMAL)
    (d / ".exit_code_1025").write_text("")
    assert outcome_of_folder(d).startswith("failed (exit code 1025)")
