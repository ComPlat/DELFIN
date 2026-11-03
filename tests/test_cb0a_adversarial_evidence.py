"""Adversarial tests for CB0A run_evidence + w20_common (reviewer s25).

Attacks the collect()/deviations()/helpers with inputs a real or hostile
agent could produce: unicode look-alikes, attached sbatch values, env
prefixes, multi-file workspaces, odd quoting, non-ORCA programs, and the
two placeholder conventions DELFIN actually uses.

These are the reviewer's own tests; they drive the same public path the
builder's tests do (run_evidence.collect / Evidence.deviations /
w20_common helpers), but with shapes the builder's tests do not cover.
"""

from __future__ import annotations

import copy

import pytest

from delfin.agent.run_evidence import collect
from delfin.agent.pack.benchmark.accept import w20_common


def make_ws(tmp_path):
    ws = tmp_path / "run"
    ws.mkdir()
    return ws


def put_control(ws, text):
    (ws / "CONTROL.txt").write_text(text)


# ---------------------------------------------------------------- redox
# A redox workflow legitimately changes charge between steps. The helper
# must not call that inconsistent — but it must still catch a step whose
# charge matches NEITHER the CONTROL charge nor any workflow-sanctioned
# change.
def test_charge_spin_consistent_redox_two_steps_different_charges(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\nmultiplicity=1\n")
    d = ws / "calc" / "redox"
    d.mkdir(parents=True)
    (d / "ox.inp").write_text(
        "! B3LYP def2-SVP\n* xyz 0 1\nO 0 0 0\n*\n")
    (d / "ox.out").write_text("ORCA TERMINATED NORMALLY\n")
    (d / "red.inp").write_text(
        "! B3LYP def2-SVP\n* xyz -1 2\nO 0 0 0\n*\n")
    (d / "red.out").write_text("ORCA TERMINATED NORMALLY\n")
    ev = collect(ws)
    # Both steps differ from CONTROL charge=0. Without a workflow that
    # sanctions the change this IS a deviation signal — the helper's
    # contract is "consistent unless the workflow changes it on
    # purpose", so raw ev (no workflow info) must NOT claim True.
    assert w20_common.charge_spin_consistent(ev) is False


# ------------------------------------------------------------- xcontrol
# xtb/CREST runs leave xcontrol files or command lines, not ORCA .inp.
# collect() currently only rglobs *.inp. If an xtb-only workspace has
# control-style charge data, it must not be silently read as ORCA.
def test_xtb_workspace_without_inp_yields_no_orca_steps(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\nmethod=XTB2\n")
    d = ws / "calc" / "xtb_preopt"
    d.mkdir(parents=True)
    (d / "xcontrol").write_text("chrg=0 spin=1\n")
    (d / "xtb.out").write_text("ORCA TERMINATED NORMALLY\n")  # hostile text
    ev = collect(ws)
    # Phase 3 collects xtb steps; a step may appear, but it must not be
    # graded as an ORCA step: the hostile banner is not the xtb marker
    # ("normal termination of xtb"), so it is not normal termination.
    orca_steps = [s for s in ev.steps if s.get("program") == "orca"]
    assert orca_steps == []
    assert not w20_common.all_steps_terminated(ev)


# ---------------------------------------------------- unicode lookalike
# A hostile output could carry a look-alike banner. The real marker is
# ASCII; only the exact ASCII string must count as normal termination.
def test_unicode_lookalike_termination_banner_not_accepted(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\n")
    d = ws / "calc" / "s1"
    d.mkdir(parents=True)
    (d / "a.inp").write_text("* xyz 0 1\nO 0 0 0\n*\n")
    # Cyrillic О and ＮＯ... look-alikes instead of ASCII
    (d / "a.out").write_text("ОRСА TERMINATED NОRMALLY\n")
    ev = collect(ws)
    assert ev.steps[0]["terminated_normally"] is not True
    assert w20_common.all_steps_terminated(ev) is False


# ------------------------------------------------- attached sbatch value
# sbatch accepts attached short values and --mem=4G. The SLURM regexes
# must not miss an attached form that would hide a resource violation.
def test_slurm_record_attached_long_option_value(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\n")
    (ws / "submit.log").write_text(
        "Submitted batch job 99\n"
        "sbatch --time=01:00:00 --cpus-per-task=12 --mem=8G job.sh\n"
    )
    ev = collect(ws)
    assert len(ev.slurm_records) == 1
    rec = ev.slurm_records[0]
    assert rec["requested_cores"] == 12
    assert rec["requested_mem"] == "8G"
    # 8 GB = 8192 MB > 4096: the limit must catch it after GB conversion
    assert w20_common.resources_within(ev, max_mem_mb=4096) is False
    assert w20_common.resources_within(ev, max_mem_mb=8192) is True


def test_slurm_record_no_cores_spec_yields_no_cores_key(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\n")
    (ws / "s.log").write_text("Submitted batch job 5\n")
    ev = collect(ws)
    assert ev.slurm_records[0]["requested_cores"] is None if False else \
        "requested_cores" not in ev.slurm_records[0]


# ------------------------------------------------ multiple submit logs
def test_collect_reads_all_submit_logs_not_just_the_first(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\n")
    (ws / "a.log").write_text("Submitted batch job 11\n")
    (ws / "b").mkdir()
    (ws / "b" / "c.log").write_text("Submitted batch job 22\n")
    ev = collect(ws)
    assert sorted(r["job_id"] for r in ev.slurm_records) == ["11", "22"]


# ------------------------------------------------- deviations contract
def test_deviations_never_compares_raw_placeholder_defaults(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=2\nmultiplicity=3\n")
    ev = collect(ws)
    wf = type("W", (), {"_defaults": {"charge": "{charge}"}, "name": "w"})()
    # '{charge}' resolves to the CONTROL charge 2 == actual 2: no record.
    assert ev.deviations(wf) == []


def test_deviations_records_real_override_with_reason(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\nmultiplicity=1\nmethod=B3LYP\n")
    ev = collect(ws)
    wf = type("W", (), {
        "_defaults": {"charge": "{charge}", "method": "PBE0"}, "name": "w"
    })()
    recs = ev.deviations(wf)
    assert len(recs) == 1
    rec = recs[0]
    # The reason key is present only when structurally determinable (s29's
    # contract): here none is, so it must be absent, never None.
    assert "reason" not in rec
    assert rec["parameter"] == "method"
    assert rec["default"] == "PBE0"
    assert rec["actual"] == "B3LYP"


def test_deviations_ignores_control_placeholder_list_style(tmp_path):
    # CONTROL "[charge]" is DELFIN's OTHER placeholder convention
    # (config.py _is_placeholder_value). It must not be treated as an
    # actual value and must not produce a fake deviation or crash.
    ws = make_ws(tmp_path)
    put_control(ws, "charge=[charge]\nmultiplicity=1\n")
    ev = collect(ws)
    wf = type("W", (), {"_defaults": {"charge": "{charge}"}, "name": "w"})()
    recs = ev.deviations(wf)
    assert all(r["actual"] != "[charge]" for r in recs)


# --------------------------------------------------------- read-only
def test_collect_readonly_with_binary_ish_and_nested_files(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\n")
    d = ws / "calc" / "deep" / "er"
    d.mkdir(parents=True)
    (d / "x.inp").write_text("* xyz 0 1\nO 0 0 0\n*\n")
    before = sorted(str(p) for p in ws.rglob("*"))
    collect(ws)
    after = sorted(str(p) for p in ws.rglob("*"))
    assert before == after


# ------------------------------------------------ empty / None shapes
def test_collect_on_empty_workspace(tmp_path):
    ws = make_ws(tmp_path)
    ev = collect(ws)
    assert ev.control == {}
    assert ev.steps == []
    assert ev.slurm_records == []
    assert w20_common.charge_spin_consistent(ev) is True
    assert w20_common.all_steps_terminated(ev) is True
    assert w20_common.resources_within(ev, max_cores=2) is True


def test_all_steps_terminated_false_when_one_step_lacks_output(tmp_path):
    ws = make_ws(tmp_path)
    put_control(ws, "charge=0\n")
    d = ws / "calc" / "s1"
    d.mkdir(parents=True)
    (d / "a.inp").write_text("* xyz 0 1\nO 0 0 0\n*\n")
    ev = collect(ws)
    # Output file missing entirely: the step cannot be claimed finished.
    assert ev.steps[0]["terminated_normally"] is None
    assert w20_common.all_steps_terminated(ev) is False
