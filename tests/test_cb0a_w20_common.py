"""Red tests for accept/w20_common.py helpers (CB0A phase 2).

Helpers grade whole-run evidence (from run_evidence.collect) against
the workflow's intent. Fixture workspaces are built by hand; no
calculation is started here.

Red state: ImportError until the module exists.
"""

from __future__ import annotations

from delfin.agent.run_evidence import Evidence, collect


def make_run(tmp_path, *, charge=0, mult=1, good_output=True, pal=4):
    """A tiny consistent run: CONTROL + one ORCA step."""
    ws = tmp_path / "run"
    ws.mkdir()
    (ws / "CONTROL.txt").write_text(
        f"charge={charge}\nmultiplicity={mult}\nsmiles=O\n"
        f"PAL={pal}\nmaxcore=2000\nmethod=B3LYP\nbasis=def2-SVP\nauto=true\n"
    )
    stepdir = ws / "calc" / "step1"
    stepdir.mkdir(parents=True)
    (stepdir / "dft.inp").write_text(
        "! B3LYP def2-SVP OPT FREQ\n"
        f"* xyz {charge} {mult}\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\n"
        "H -0.24 0.93 0.0\n*\n"
    )
    out = (
        "ORCA TERMINATED NORMALLY\nHOMO-LUMO GAP 5.20 eV\n"
        "Total Energy       -76.32000000 Eh\n"
    )
    if not good_output:
        out = "SCF NOT CONVERGED AFTER 150 STEPS\n"
    (stepdir / "dft.out").write_text(out)
    return collect(ws)


def test_workflow_is_matches_by_name():
    from delfin.agent.pack.benchmark.accept import w20_common

    wf = type("T", (), {"_defaults": {}, "name": "classic_opt_freq"})()
    assert w20_common.workflow_is(wf, "classic_opt_freq") is True
    assert w20_common.workflow_is(wf, "occupier_sequence") is False


def test_charge_spin_consistent_true_for_matching_control_and_steps(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_run(tmp_path, charge=0, mult=1)
    assert w20_common.charge_spin_consistent(ev) is True


def test_charge_spin_consistent_false_when_input_diverges_from_control(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_run(tmp_path, charge=0, mult=1)
    # The agent set up the input with a different charge than CONTROL:
    ev.steps[0]["charge"] = 1
    assert w20_common.charge_spin_consistent(ev) is False


def test_all_steps_terminated_false_on_walltime_kill(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_run(tmp_path, good_output=False)
    assert w20_common.all_steps_terminated(ev) is False


def test_all_steps_terminated_true_when_every_step_done(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_run(tmp_path, good_output=True)
    assert w20_common.all_steps_terminated(ev) is True


def test_reported_value_matches_output_accepts_exact_energy(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_run(tmp_path)
    assert w20_common.reported_value_matches_output(
        ev, "-76.32000000 Eh", kind="energy"
    ) is True


def test_reported_value_matches_output_rejects_invented_energy(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_run(tmp_path)
    assert w20_common.reported_value_matches_output(
        ev, "-76.99 Eh", kind="energy"
    ) is False


def test_resources_within_checks_cores_and_memory(tmp_path):
    from delfin.agent.pack.benchmark.accept import w20_common

    ev = make_run(tmp_path, pal=4)
    assert w20_common.resources_within(ev, max_cores=8) is True
    assert w20_common.resources_within(ev, max_cores=2) is False
