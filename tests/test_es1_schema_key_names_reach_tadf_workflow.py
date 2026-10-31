"""CONTROL schema key names reach the TADF-xTB workflow.

The validator's CONTROL_FIELD_SPECS name the settings `xTB_method`
(control_validator.py:2192, default XTB2) and `tadf_xTB_excited_method` /
`tadf_xTB_energy_window` (:2271/:2273).  A CONTROL.txt written with those
keys is what every shipped template and the validator produce, so the
workflow must accept them.  Before the fix the wrapper only looked for
`xtb_method` / `excited_method` / `energy_window` and raised ValueError
naming all three even though CONTROL had them.
"""

from pathlib import Path
from unittest.mock import patch

import pytest

from delfin.workflows.contrib.tadf_xtb_workflow import TadfXtbWorkflow

_ENTRY = object()


def _schema_named_config() -> dict:
    return {
        "charge": 0,
        "xTB_method": "XTB2",
        "tadf_xTB_excited_method": "stda",
        "tadf_xTB_energy_window": 10.0,
        "PAL": 4,
        "maxcore": 2000,
    }


def _run(config, **kwargs):
    wf = TadfXtbWorkflow()
    with patch("delfin.tadf_xtb.run_single_tadf_xtb") as runner:
        wf.run(config=config, entry=_ENTRY, workdir=Path("wd"), **kwargs)
    (call,) = runner.call_args_list
    return call.kwargs


def test_schema_named_control_keys_reach_the_workflow():
    params = _run(_schema_named_config())
    assert params["xtb_method"] == "XTB2"
    assert params["excited_method"] == "stda"
    assert params["energy_window"] == 10.0
    assert params["cores"] == 4
    assert params["charge"] == 0


def test_alias_key_names_still_work():
    """The phase-1 aliases (xtb_method, excited_method, ...) keep running."""
    config = {
        "charge": 0,
        "xtb_method": "GFN2-xTB",
        "excited_method": "std2",
        "energy_window": 12.0,
        "PAL": 2,
        "maxcore": 1000,
    }
    params = _run(config)
    assert params["xtb_method"] == "GFN2-xTB"
    assert params["excited_method"] == "std2"
    assert params["energy_window"] == 12.0


def test_schema_keys_beaten_by_explicit_kwargs():
    params = _run(
        {**_schema_named_config()}, excited_method="std2"
    )
    assert params["excited_method"] == "std2"


def test_missing_required_setting_still_refused():
    wf = TadfXtbWorkflow()
    config = {"charge": 0, "xTB_method": "XTB2", "PAL": 2, "maxcore": 1000}
    with pytest.raises(ValueError) as excinfo:
        wf.run(config=config, entry=_ENTRY, workdir=Path("wd"))
    assert "excited_method" in str(excinfo.value)
    assert "energy_window" in str(excinfo.value)
