"""Red tests for whole-run grading (CB0A phase 3).

A faked run — output text written by hand with no program input, or an
input whose charge differs from CONTROL — grades red. A real, tiny,
consistent run grades green. The real run is xtb GFN2 on water, local,
one core, seconds long (practice-run rules: <= 20 atoms, <= 2 cores).
"""

from __future__ import annotations

import shutil

import pytest

from delfin.agent.run_evidence import collect


def make_fake_run(tmp_path):
    """An output that claims success but no input was ever written."""
    ws = tmp_path / "fake"
    ws.mkdir()
    (ws / "CONTROL.txt").write_text(
        "charge=0\nmultiplicity=1\nsmiles=O\nPAL=4\nmaxcore=2000\n"
        "method=B3LYP\nbasis=def2-SVP\nauto=true\n"
    )
    (ws / "calc").mkdir()
    (ws / "calc" / "step1.out").write_text(
        "ORCA TERMINATED NORMALLY\nTotal Energy -76.32 Eh\n"
    )
    return collect(ws)


def make_divergent_run(tmp_path):
    """An input whose charge differs from CONTROL (charge=0 vs +1)."""
    ws = tmp_path / "divergent"
    ws.mkdir()
    (ws / "CONTROL.txt").write_text(
        "charge=0\nmultiplicity=1\nsmiles=O\nPAL=4\nmaxcore=2000\nauto=true\n"
    )
    stepdir = ws / "calc" / "step1"
    stepdir.mkdir(parents=True)
    (stepdir / "dft.inp").write_text(
        "! B3LYP def2-SVP OPT FREQ\n"
        "* xyz 1 1\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n*\n"
    )
    (stepdir / "dft.out").write_text("ORCA TERMINATED NORMALLY\n")
    return collect(ws)


def test_fake_run_with_output_but_no_input_grades_red(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_fake_run(tmp_path)
    verdict = w20_common.grade_run(ev)
    assert verdict["verdict"] == "fail"
    assert "input" in " ".join(verdict["problems"]).lower()


def test_input_charge_diverging_from_control_grades_red(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_divergent_run(tmp_path)
    verdict = w20_common.grade_run(ev)
    assert verdict["verdict"] == "fail"
    assert "charge" in " ".join(verdict["problems"]).lower()


def test_consistent_run_with_inputs_and_outputs_grades_green(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ws = tmp_path / "good"
    ws.mkdir()
    (ws / "CONTROL.txt").write_text(
        "charge=0\nmultiplicity=1\nsmiles=O\nPAL=2\nmaxcore=2000\nauto=true\n"
    )
    stepdir = ws / "calc" / "step1"
    stepdir.mkdir(parents=True)
    (stepdir / "dft.inp").write_text(
        "! B3LYP def2-SVP OPT FREQ\n"
        "* xyz 0 1\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n*\n"
    )
    (stepdir / "dft.out").write_text(
        "ORCA TERMINATED NORMALLY\nTotal Energy -76.32 Eh\n"
    )
    ev = collect(ws)
    verdict = w20_common.grade_run(ev)
    assert verdict["verdict"] == "pass", verdict


def _xtb_path():
    """xtb as DELFIN resolves it (qm_runtime), else PATH, else None."""
    from delfin.qm_runtime import resolve_tool

    r = resolve_tool("xtb")
    if r is not None and getattr(r, "path", None):
        return str(r.path)
    return shutil.which("xtb")


_xtb_path()


def test_deviations_never_compares_raw_placeholder_defaults(tmp_path):
    """s25 F1: a "{charge}" default with CONTROL charge=2 and no steps
    must not produce a fake deviation (raw '{charge}' vs 2)."""
    from delfin.agent.run_evidence import Evidence

    ev = Evidence(workspace=str(tmp_path), control={"charge": 2})
    records = ev.deviations({"charge": "{charge}"})
    assert records == []


def test_deviations_skips_control_placeholder_list_style(tmp_path):
    """s25 F3: a CONTROL '[charge]' placeholder is not an actual value."""
    from delfin.agent.run_evidence import Evidence

    ev = Evidence(
        workspace=str(tmp_path), control={"charge": "[charge]"})
    records = ev.deviations({"charge": "{charge}"})
    assert records == []


def test_deviations_records_real_override_without_placeholder_default(tmp_path):
    """A CONCRETE template default that CONTROL overrides is one record;
    a placeholder default contributes none (s25 F2)."""
    from delfin.agent.run_evidence import Evidence

    ev = Evidence(
        workspace=str(tmp_path),
        control={"charge": 2, "solvent": "dmf"})
    records = ev.deviations(
        {"charge": "{charge}", "solvent": "water"})
    assert [r["parameter"] for r in records] == ["solvent"]
    assert records[0]["default"] == "water"
    assert records[0]["actual"] == "dmf"


def test_xtb_termination_marker_present_with_rc0_is_terminated(tmp_path):
    """QS s29 case 1: marker present in the captured stream."""
    from delfin.agent.run_evidence import collect

    ws = tmp_path
    stepdir = ws / "calc" / "s1"
    stepdir.mkdir(parents=True)
    (stepdir / "geom.xyz").write_text(
        "3\nw\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n")
    (stepdir / "xtb_opt.out").write_text(
        "some output\nnormal termination of xtb\nrc=0\n")
    ev = collect(ws)
    step = ev.steps[0]
    assert step["program"] == "xtb"
    assert step["terminated_normally"] is True
    assert "termination_note" not in step


def test_xtb_rc0_without_marker_is_terminated_with_note(tmp_path):
    """QS s29 case 2: rc=0 is the primary signal; marker-absent gets a
    note, never a silent ok."""
    from delfin.agent.run_evidence import collect

    ws = tmp_path
    stepdir = ws / "calc" / "s1"
    stepdir.mkdir(parents=True)
    (stepdir / "geom.xyz").write_text(
        "3\nw\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n")
    (stepdir / "xtb_opt.out").write_text("timings only, no banner\nrc=0\n")
    ev = collect(ws)
    step = ev.steps[0]
    assert step["terminated_normally"] is True
    assert "termination_note" in step


def test_xtb_nonzero_rc_grades_not_terminated(tmp_path):
    """QS s29 case 3: rc!=0 -> red, marker presence irrelevant."""
    from delfin.agent.run_evidence import collect

    ws = tmp_path
    stepdir = ws / "calc" / "s1"
    stepdir.mkdir(parents=True)
    (stepdir / "geom.xyz").write_text(
        "3\nw\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n")
    (stepdir / "xtb_opt.out").write_text(
        "crashed mid-SCF\nrc=1\n")
    ev = collect(ws)
    assert ev.steps[0]["terminated_normally"] is False


def test_xtb_no_rc_no_marker_fails_closed(tmp_path):
    """QS s29 fail-closed: output ends mid-SCF with neither rc record
    nor marker -> NOT terminated."""
    from delfin.agent.run_evidence import collect

    ws = tmp_path
    stepdir = ws / "calc" / "s1"
    stepdir.mkdir(parents=True)
    (stepdir / "geom.xyz").write_text(
        "3\nw\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n")
    (stepdir / "xtb_opt.out").write_text("SCF iteration 12 ...\n")
    ev = collect(ws)
    assert ev.steps[0]["terminated_normally"] is False


_needs_xtb = pytest.mark.skipif(
    _xtb_path() is None, reason="no xtb installed (resolver + PATH)")


@_needs_xtb
def test_real_xtb_gfn2_water_run_grades_green(tmp_path):
    """A real tiny run (xtb GFN2, water, 1 core, seconds) grades green."""
    import subprocess

    from delfin.agent.pack.benchmark.accept import w20_common

    ws = tmp_path / "real"
    (ws / "calc" / "xtb_opt").mkdir(parents=True)
    (ws / "CONTROL.txt").write_text(
        "charge=0\nmultiplicity=1\nsmiles=O\nPAL=1\nmaxcore=1000\nauto=true\n"
    )
    stepdir = ws / "calc" / "xtb_opt"
    (stepdir / "geom.xyz").write_text(
        "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"
    )
    subprocess.run(
        [_xtb_path(), "geom.xyz", "--gfn", "2", "--chrg", "0", "--uhf", "0",
         "--opt", "tight", "--alpb", "water", "-P", "1"],
        cwd=stepdir, stdout=open(stepdir / "xtb_opt.out", "w"),
        stderr=subprocess.STDOUT, timeout=300, check=False,
    )
    ev = collect(ws)
    # The run must actually have left evidence before we grade it.
    assert any(step["program"] == "xtb" for step in ev.steps), (
        "xtb step not collected from the real run")
    verdict = w20_common.grade_run(ev)
    assert verdict["verdict"] == "pass", verdict
