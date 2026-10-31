"""Adversarial tests for the ES1 workflow wrappers (reviewer nacht-s8).

Two defects, both red against main (54462fe67):

1. ``TadfXtbWorkflow.run`` accepts a ``config`` mapping but never uses it:
   it forwards only ``**kwargs`` to ``run_single_tadf_xtb``, whose signature
   (delfin/tadf_xtb.py:599) requires charge, multiplicity, xtb_method,
   excited_method, energy_window, cores, maxcore and more.  CONTROL.txt
   settings therefore never reach the TADF-xTB screening through the
   workflow entry, and a CONTROL-only call raises TypeError.
2. ``EsdWorkflow.run_cli`` swallows every exception with a bare
   ``except Exception: return 1`` and logs nothing: the user sees exit
   code 1 with no message at all.
"""

from __future__ import annotations

import io
import logging
from contextlib import redirect_stderr
from pathlib import Path

import pytest

from delfin.workflows.contrib.esd_workflow import EsdWorkflow
from delfin.workflows.contrib.tadf_xtb_workflow import TadfXtbWorkflow


# --- 1. CONTROL reaches the TADF-xTB workflow entry --------------------


def test_tadf_xtb_workflow_config_reaches_entry(tmp_path, monkeypatch):
    """A CONTROL-style config provides the required defaults of
    run_single_tadf_xtb; the workflow must pass them through instead of
    raising TypeError on undefined kwargs."""
    captured = {}

    def fake_run_single(entry, *, charge, multiplicity, xtb_method,
                        excited_method, energy_window, cores, maxcore,
                        workdir, use_crest, use_goat, run_t1_opt,
                        t1_multiplicity, **rest):
        captured.update(
            charge=charge, xtb_method=xtb_method,
            excited_method=excited_method, energy_window=energy_window,
            cores=cores, maxcore=maxcore,
        )
        return object()  # opaque stand-in: the wrapper only forwards

    monkeypatch.setattr(
        "delfin.tadf_xtb.run_single_tadf_xtb", fake_run_single
    )
    wf = TadfXtbWorkflow()
    config = {
        "charge": 1,
        "xtb_method": "GFN2-xTB",
        "excited_method": "sTDA-xTB",
        "energy_window": 1.5,
        "pal": 4,
        "maxcore": 2000,
    }
    wf.run(
        config=config,
        entry=object(),
        workdir=tmp_path,
        use_crest=False,
        use_goat=False,
        run_t1_opt=False,
        t1_multiplicity=3,
    )
    assert captured["charge"] == 1
    assert captured["xtb_method"] == "GFN2-xTB"
    assert captured["excited_method"] == "sTDA-xTB"
    assert captured["energy_window"] == 1.5
    assert captured["cores"] == 4
    assert captured["maxcore"] == 2000


def test_tadf_xtb_workflow_kwargs_beat_config(tmp_path, monkeypatch):
    """An explicit kwarg must win over the same value in config.

    Rewritten along its intent after finding F2: the original pinned a
    call that passed only ``charge`` explicitly and expected it to run,
    which contradicts the per-key required-setting check F2 asks for
    (one explicit key does not vouch for the others).  The intent — an
    explicit kwarg beats the same value in config — is tested with a
    complete call."""
    captured = {}

    def fake_run_single(entry, *, charge, **rest):
        captured["charge"] = charge
        return object()

    monkeypatch.setattr(
        "delfin.tadf_xtb.run_single_tadf_xtb", fake_run_single
    )
    TadfXtbWorkflow().run(
        config={
            "charge": 7,
            "xtb_method": "XTB1",
            "excited_method": "stda",
            "energy_window": 10.0,
            "PAL": 4,
            "maxcore": 1000,
        },
        entry=object(),
        charge=2,
        workdir=tmp_path,
    )
    assert captured["charge"] == 2


def test_tadf_xtb_workflow_control_only_call_raises_typeerror_on_main():
    """RED CONTROL on main: a CONTROL-only call today dies with TypeError
    because no default of run_single_tadf_xtb is mapped.  After the fix it
    must raise a clear error naming the missing piece, never TypeError on
    undefined kwargs."""
    wf = TadfXtbWorkflow()
    with pytest.raises(Exception) as excinfo:
        wf.run(config={}, entry=object(), workdir=Path("."))
    assert not isinstance(
        excinfo.value, TypeError
    ), "CONTROL-only call still fails with an unexplained TypeError"


# --- 2. esd run_cli reports the failure --------------------------------


def test_esd_run_cli_failure_is_logged(monkeypatch, caplog):
    """RED on main: run_cli swallows the traceback and prints nothing.
    After the fix the exception's message must reach the log/stderr before
    the exit code 1."""
    def boom(**kwargs):
        raise RuntimeError("ESD driver exploded: missing S1 hessian")

    monkeypatch.setattr("delfin.esd_module.execute_esd_jobs", boom)
    monkeypatch.setattr(
        "delfin.config.read_control_file", lambda: {"charge": 0}
    )
    wf = EsdWorkflow()
    with caplog.at_level(logging.ERROR):
        buf = io.StringIO()
        with redirect_stderr(buf):
            rc = wf.run_cli([])
    assert rc == 1
    combined = buf.getvalue() + caplog.text
    assert "ESD driver exploded" in combined, (
        "run_cli returned 1 but the user never sees why"
    )


def test_esd_run_cli_success_still_zero(monkeypatch):
    """NOTHING NEW IS REFUSED: a healthy run still exits 0."""
    monkeypatch.setattr(
        "delfin.esd_module.execute_esd_jobs", lambda **kwargs: None
    )
    monkeypatch.setattr(
        "delfin.config.read_control_file", lambda: {"charge": 0}
    )
    assert EsdWorkflow().run_cli([]) == 0
