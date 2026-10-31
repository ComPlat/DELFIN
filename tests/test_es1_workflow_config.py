"""CONTROL.txt reaches the ESD/TADF workflows, and a failed CLI reports why.

Two contracts:

* A workflow's ``run()`` receives ``config`` (the parsed CONTROL.txt dict,
  see delfin/workflows/types.py:26). A workflow that silently drops it makes
  every CONTROL setting a no-op on that path.
* A CLI that fails says why. ``run_cli`` returning 1 with nothing in the log
  leaves the user with a bare exit code.

Both tests fail on the code before the fix:

* ``TadfXtbWorkflow.run`` (delfin/workflows/contrib/tadf_xtb_workflow.py)
  passes only **kwargs to ``run_single_tadf_xtb`` — config-derived values
  never arrive, and the required arguments (charge, multiplicity, …) are
  missing entirely, so a contract-faithful call raises TypeError.
* ``EsdWorkflow.run_cli`` (delfin/workflows/contrib/esd_workflow.py:53-54)
  swallows the exception without logging it.
"""

from __future__ import annotations

import logging

import pytest


# ---------------------------------------------------------------------------
# TadfXtbWorkflow.run: config values become defaults, explicit kwargs win
# ---------------------------------------------------------------------------

def test_tadf_workflow_run_plumbing_kwargs_do_not_skip_the_check(monkeypatch):
    """A registry-style call with plumbing kwargs (entry, workdir) but a
    config that lacks a required setting must raise — plumbing kwargs are
    not proof the physics settings were given (reviewer finding F1)."""
    from delfin.workflows.contrib.tadf_xtb_workflow import TadfXtbWorkflow

    captured = {}
    monkeypatch.setattr(
        "delfin.tadf_xtb.run_single_tadf_xtb",
        lambda **k: captured.update(k) or "x",
    )

    wf = TadfXtbWorkflow()
    with pytest.raises(ValueError) as excinfo:
        wf.run(config={"xtb_method": "XTB2"}, entry=object(), workdir="d")
    assert "charge" in str(excinfo.value)
    assert captured == {}  # run_single_tadf_xtb never called


def test_tadf_workflow_run_one_explicit_key_does_not_vouch_for_the_others(monkeypatch):
    """Mixed call: charge passed explicitly, xtb_method missing from both
    config and kwargs — must raise naming xtb_method (reviewer finding F2);
    one explicit key never vouches for the others."""
    from delfin.workflows.contrib.tadf_xtb_workflow import TadfXtbWorkflow

    captured = {}
    monkeypatch.setattr(
        "delfin.tadf_xtb.run_single_tadf_xtb",
        lambda **k: captured.update(k) or "x",
    )

    wf = TadfXtbWorkflow()
    with pytest.raises(ValueError) as excinfo:
        wf.run(
            config={
                "charge": 0,
                "excited_method": "stda",
                "energy_window": 10.0,
                "PAL": 4,
                "maxcore": 1000,
            },
            charge=0,
        )
    assert "xtb_method" in str(excinfo.value)
    assert captured == {}


def test_tadf_workflow_run_explicit_none_counts_as_absent(monkeypatch):
    """An explicitly passed None kwarg is treated like a missing setting."""
    from delfin.workflows.contrib.tadf_xtb_workflow import TadfXtbWorkflow

    monkeypatch.setattr("delfin.tadf_xtb.run_single_tadf_xtb", lambda **k: "x")

    wf = TadfXtbWorkflow()
    with pytest.raises(ValueError) as excinfo:
        wf.run(
            config={
                "xtb_method": "XTB2",
                "excited_method": "stda",
                "energy_window": 10.0,
                "PAL": 4,
                "maxcore": 1000,
            },
            charge=None,
        )
    assert "charge" in str(excinfo.value)



def _config() -> dict:
    """A parsed CONTROL.txt the way read_control_file produces one."""
    return {
        "charge": 1,
        "multiplicity": 2,
        "xtb_method": "XTB2",
        "excited_method": "std2",
        "energy_window": 12.0,
        "PAL": 8,
        "maxcore": 2000,
    }


def test_tadf_workflow_run_maps_config_to_defaults(monkeypatch, tmp_path):
    """config values reach run_single_tadf_xtb as defaults."""
    from delfin.workflows.contrib.tadf_xtb_workflow import TadfXtbWorkflow

    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return "result"

    monkeypatch.setattr("delfin.tadf_xtb.run_single_tadf_xtb", fake_run)

    wf = TadfXtbWorkflow()
    result = wf.run(
        config=_config(),
        entry=None,
        workdir=tmp_path,
    )
    assert result == "result"
    assert captured["charge"] == 1
    assert captured["multiplicity"] == 2
    assert captured["xtb_method"] == "XTB2"
    assert captured["excited_method"] == "std2"
    assert captured["energy_window"] == 12.0
    assert captured["cores"] == 8
    assert captured["maxcore"] == 2000
    assert captured["workdir"] == tmp_path


def test_tadf_workflow_run_explicit_kwargs_win(monkeypatch, tmp_path):
    """An explicit kwarg overrides the config default."""
    from delfin.workflows.contrib.tadf_xtb_workflow import TadfXtbWorkflow

    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return "result"

    monkeypatch.setattr("delfin.tadf_xtb.run_single_tadf_xtb", fake_run)

    wf = TadfXtbWorkflow()
    wf.run(config=_config(), entry=None, workdir=tmp_path, charge=2)
    assert captured["charge"] == 2


def test_tadf_workflow_run_raises_clearly_when_settings_missing(monkeypatch, tmp_path):
    """Physics-relevant settings (xtb_method, excited_method, …) get no silent
    fallback: a config-only call missing them raises a ValueError that names
    the missing setting — never the bare TypeError a call with undefined
    kwargs would produce.  An explicit-kwargs call runs as given (the CLI
    path does that; workdir/entry do not count as settings)."""
    from delfin.workflows.contrib.tadf_xtb_workflow import TadfXtbWorkflow

    captured = {}
    monkeypatch.setattr(
        "delfin.tadf_xtb.run_single_tadf_xtb",
        lambda **k: captured.update(k) or "x",
    )

    wf = TadfXtbWorkflow()
    with pytest.raises(ValueError) as excinfo:
        wf.run(config={})
    assert "excited_method" in str(excinfo.value)

    captured.clear()
    wf.run(config={}, entry=None, workdir=tmp_path, charge=0,
           xtb_method="XTB2", excited_method="stda", energy_window=10.0,
           cores=4, maxcore=1000)
    assert captured["multiplicity"] == 1  # S0 default, no silent refusal


# ---------------------------------------------------------------------------
# EsdWorkflow.run_cli: a failure is logged
# ---------------------------------------------------------------------------

def test_esd_workflow_run_cli_logs_the_failure(monkeypatch, caplog):
    """run_cli returns 1 AND logs what went wrong."""
    from delfin.workflows.contrib.esd_workflow import EsdWorkflow

    def boom(**kwargs):
        raise RuntimeError("SCF did not converge in S1")

    monkeypatch.setattr("delfin.esd_module.execute_esd_jobs", boom)
    import delfin.config

    monkeypatch.setattr(delfin.config, "read_control_file", lambda: {})

    wf = EsdWorkflow()
    with caplog.at_level(logging.ERROR):
        code = wf.run_cli(["--charge", "0"])

    assert code == 1
    assert any("SCF did not converge" in rec.message for rec in caplog.records), (
        "the exception reached no log record; the user sees a bare exit 1"
    )


def test_esd_workflow_run_cli_success_returns_zero(monkeypatch):
    """The ordinary path still returns 0 (nothing new is refused)."""
    from delfin.workflows.contrib.esd_workflow import EsdWorkflow

    monkeypatch.setattr(
        "delfin.esd_module.execute_esd_jobs", lambda **kwargs: None
    )
    import delfin.config

    monkeypatch.setattr(delfin.config, "read_control_file", lambda: {})

    wf = EsdWorkflow()
    assert wf.run_cli([]) == 0
