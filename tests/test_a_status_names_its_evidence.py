"""Asked whether three runs had succeeded, failed or were still going,
the operator answered correctly -- by assembling, for each folder, the
exit-code marker, the run log's last line and the output's termination
line by hand, and named the lack of one call that does this as the one
change that would have helped most (2026-09-11). It also read
parse_orca_output's "file not found" as a calculation failure, because
nothing in the result said whether there had been an output to parse.

calc_status is that call: the verdict extract_energy_table already puts
in its outcome column, per folder, with every file that had a say.
parse_orca_output now carries status and outcome, and accepts a folder.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from delfin import api
from delfin.ops_server import server as ops

_ROOT = Path(__file__).resolve().parents[1]
_SETUP = _ROOT / "delfin" / "agent" / "pack" / "benchmark" / "setup" / "three_runs_to_audit.py"


@pytest.fixture
def runs(tmp_path):
    spec = importlib.util.spec_from_file_location("three_runs", _SETUP)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    root = tmp_path / "runs"
    mod.build(root)
    return root


def _sources(status):
    return [e["source"] for e in status.evidence]


def test_a_succeeded_run_is_succeeded_with_its_marker_named(runs):
    st = api.calculation_status(str(runs / "succeeded"))
    assert st.state == "succeeded"
    assert st.outcome.startswith("succeeded (exit code 0)")
    assert ".exit_code_0" in _sources(st)
    assert st.method


def test_a_failed_run_is_failed_with_the_code_and_the_log(runs):
    st = api.calculation_status(str(runs / "failed"))
    assert st.state == "failed"
    assert "1025" in st.outcome
    assert ".exit_code_1025" in _sources(st)
    log = [e for e in st.evidence if e["source"] == "delfin_run.log"]
    assert log and log[0]["says"]


def test_a_running_run_is_running_with_the_logs_last_line(runs):
    st = api.calculation_status(str(runs / "running"))
    assert st.state == "running"
    assert not any(s.startswith(".exit_code") for s in _sources(st))
    log = [e for e in st.evidence if e["source"] == "delfin_run.log"]
    assert log and "running" in log[0]["says"].lower()


def test_a_folder_with_only_an_input_is_pending(tmp_path):
    d = tmp_path / "queued"
    d.mkdir()
    (d / "run.inp").write_text("! PBE0 def2-SVP Opt\n")
    st = api.calculation_status(str(d))
    assert st.state == "pending"
    assert st.outcome.startswith("no output yet")
    assert any(e["says"] == "input present, no output" for e in st.evidence)
    assert st.method == "PBE0/def2-SVP"


def test_an_archived_output_is_finished_per_its_own_last_line(tmp_path):
    d = tmp_path / "archived"
    d.mkdir()
    (d / "run.inp").write_text("! B3LYP def2-TZVP SP\n")
    (d / "run.out").write_text("FINAL SINGLE POINT ENERGY -1.0\n"
                               "****ORCA TERMINATED NORMALLY****\n")
    st = api.calculation_status(str(d))
    assert st.state == "finished"
    assert any(e["source"] == "run.out" and "TERMINATED NORMALLY" in e["says"]
               for e in st.evidence)


def test_a_missing_folder_is_missing_not_unknown(tmp_path):
    st = api.calculation_status(str(tmp_path / "nowhere"))
    assert st.state == "missing"
    assert st.evidence == []


def test_the_tool_carries_the_evidence_and_says_what_it_is_not(runs):
    out = json.loads(ops.tool_calc_status(str(runs / "failed")))
    assert out["state"] == "failed"
    assert {"source", "says"} <= set(out["evidence"][0])
    doc = ops.tool_calc_status.__doc__ or ""
    assert "evidence" in doc and "list_active_calculations" in doc


# --- parse_orca_output: a missing output is a fact about the folder ------------


def test_parse_orca_output_on_a_folder_without_output_says_no_output(tmp_path):
    d = tmp_path / "queued"
    d.mkdir()
    (d / "run.inp").write_text("! PBE0 def2-SVP Opt\n")
    parsed = api.parse_orca_output(str(d))
    assert parsed.status == "no_output"
    assert parsed.outcome.startswith("no output yet")
    parsed2 = api.parse_orca_output(str(d / "run.out"))
    assert parsed2.status == "no_output"
    assert parsed2.error_summary == "file not found"


def test_parse_orca_output_on_a_missing_folder_says_missing(tmp_path):
    parsed = api.parse_orca_output(str(tmp_path / "nowhere" / "run.out"))
    assert parsed.status == "missing"


def test_parse_orca_output_accepts_the_folder_and_parses_its_output(runs):
    parsed = api.parse_orca_output(str(runs / "succeeded"))
    assert parsed.status == "ok"
    assert parsed.path.endswith(".out")
    assert parsed.outcome.startswith("succeeded")
    assert parsed.final_single_point is not None


def test_the_parse_wrapper_renders_status_and_outcome(runs):
    out = json.loads(ops.tool_parse_orca_output(str(runs / "failed" / "run.out")))
    assert out["status"] == "ok"
    assert out["outcome"].startswith("failed (exit code 1025)")
    assert "status" in (ops.tool_parse_orca_output.__doc__ or "")


# --- running or crashed: the files' contents cannot say, their timestamps can


def _age_all_files(folder: Path, hours: float) -> None:
    import os
    import time
    then = time.time() - hours * 3600
    for f in folder.iterdir():
        if f.is_file():
            os.utime(f, (then, then))


def test_a_run_written_recently_is_running_and_says_when(runs):
    st = api.calculation_status(str(runs / "running"))
    assert st.state == "running"
    assert st.last_activity and st.last_activity_age_s is not None
    assert st.last_activity_age_s < 600
    assert any(e["source"] == "newest file" and "last written" in e["says"]
               for e in st.evidence)


def test_a_run_untouched_for_a_day_is_stalled_not_running(runs):
    """The log GLM read said "still running (cycle 38)" and was 22 hours
    old; the outcome phrase said "running or crashed" and left the
    choice to a bash detour. The clock decides, and the age is shown."""
    _age_all_files(runs / "running", hours=22)
    st = api.calculation_status(str(runs / "running"))
    assert st.state == "stalled"
    assert st.outcome.startswith("running or crashed")      # what the files say, kept
    assert st.last_activity_age_s > 20 * 3600
    assert any("22.0 h ago" in e["says"] for e in st.evidence)


def test_a_finished_run_is_not_called_stalled_by_its_age(runs):
    _age_all_files(runs / "succeeded", hours=22)
    st = api.calculation_status(str(runs / "succeeded"))
    assert st.state == "succeeded"


def test_the_energy_table_rows_carry_the_last_activity(runs):
    rows = {Path(r["folder"]).name: r for r in api.extract_energy_table(
        [str(runs / n) for n in ("succeeded", "failed", "running")], properties=["single_point"])}
    for name in ("succeeded", "failed", "running"):
        assert rows[name]["last_activity"]
        assert rows[name]["last_activity_age_s"] is not None
    _age_all_files(runs / "running", hours=22)
    row = api.extract_energy_table([str(runs / "running")], properties=["single_point"])[0]
    assert row["last_activity_age_s"] > 20 * 3600
    doc = ops.tool_extract_energy_table.__doc__ or ""
    assert "last_activity" in doc


# --- one word beside the phrase, and the scheduler as the witness for "running"


def test_the_energy_table_rows_carry_a_state(runs):
    rows = {Path(r["folder"]).name: r for r in api.extract_energy_table(
        [str(runs / n) for n in ("succeeded", "failed", "running")] + [str(runs / "nowhere")],
        properties=["single_point"])}
    assert rows["succeeded"]["state"] == "succeeded"
    assert rows["failed"]["state"] == "failed"
    assert rows["running"]["state"] == "running"
    assert rows["nowhere"]["state"] == "missing"
    _age_all_files(runs / "running", hours=22)
    row = api.extract_energy_table([str(runs / "running")], properties=["single_point"])[0]
    assert row["state"] == "stalled"


def test_a_scheduler_job_on_the_folder_settles_running(runs, monkeypatch):
    """GLM: no field says definitively "running"; the age is only a
    proxy. The scheduler can say it outright, and is asked."""
    stale = runs / "running"
    _age_all_files(stale, hours=22)
    monkeypatch.setattr(api, "list_active_calculations", lambda: [
        {"job_id": "4711", "name": "running", "status": "RUNNING",
         "runtime_s": 12.0, "directory": str(stale)}])
    st = api.calculation_status(str(stale))
    assert st.state == "running"
    assert any(e["source"] == "scheduler" and "4711" in e["says"] for e in st.evidence)
    row = api.extract_energy_table([str(stale)], properties=["single_point"])[0]
    assert row["state"] == "running"


def test_a_silent_scheduler_is_not_evidence(runs, monkeypatch):
    monkeypatch.setattr(api, "list_active_calculations", lambda: [{"error": "squeue: not found"}])
    st = api.calculation_status(str(runs / "running"))
    assert st.state == "running"            # written recently: the clock decides
    assert not any(e["source"] == "scheduler" for e in st.evidence)
    assert "state" in (ops.tool_extract_energy_table.__doc__ or "")



def test_an_old_input_with_no_output_and_no_job_is_not_started(tmp_path, monkeypatch):
    """Asked "running, never started, or crashed?" for a folder that holds
    only an input, an operator got "pending" and looked by hand
    (2026-09-11). When no job lists the folder and nothing was written
    for as long as a stall takes, the files and the scheduler agree on
    a word: not started."""
    import os, time
    d = tmp_path / "never"
    d.mkdir()
    inp = d / "run.inp"
    inp.write_text("! PBE0 def2-SVP Opt\n")
    old = time.time() - api.STALLED_AFTER_S - 3600
    os.utime(inp, (old, old))
    os.utime(d, (old, old))
    monkeypatch.setattr(api, "_scheduler_jobs_by_dir", lambda: {})
    st = api.calculation_status(str(d))
    assert st.state == "not started"
    assert st.outcome.startswith("no output yet")
    assert any(e["source"] == "scheduler + clock" and "nothing was written" in e["says"]
               for e in st.evidence)


def test_a_fresh_input_stays_pending_and_a_listed_job_is_running(tmp_path, monkeypatch):
    d = tmp_path / "fresh"
    d.mkdir()
    (d / "run.inp").write_text("! PBE0 def2-SVP Opt\n")
    monkeypatch.setattr(api, "_scheduler_jobs_by_dir", lambda: {})
    assert api.calculation_status(str(d)).state == "pending"
    monkeypatch.setattr(api, "_scheduler_jobs_by_dir",
                        lambda: {str(d.resolve()): {"job_id": "7", "status": "RUNNING"}})
    assert api.calculation_status(str(d)).state == "running"


def test_the_energy_table_says_why_gibbs_is_null(tmp_path):
    """The default properties include gibbs and zpe; a single-point output
    has neither. An operator read the nulls as an unfinished run."""
    sp = tmp_path / "sp"
    sp.mkdir()
    (sp / "run.inp").write_text("! PBE0 def2-SVP\n")
    (sp / "run.out").write_text(
        "FINAL SINGLE POINT ENERGY       -76.400000000000\n"
        "****ORCA TERMINATED NORMALLY****\n")
    row = api.extract_energy_table([str(sp)])[0]
    assert row["gibbs"] is None and row["zpe"] is None
    assert row["single_point"] == -76.4
    assert any("no thermochemistry" in n and "use single_point" in n
               for n in row["notes"]), row["notes"]
    # asked only for the energy it has, there is nothing to note
    row2 = api.extract_energy_table([str(sp)], properties=["single_point"])[0]
    assert "notes" not in row2
