"""WB1 finding F8: the ORCA adapter contract rejects `freq_pressure`.

WB1 rule 9 (thermochemistry-conditions-missing, delfin/tools/_science_rules/
chain.py:122) requires freq_temp AND freq_pressure to be stated. WB0's
grader fixtures set freq_pressure. But the ORCA adapter contract
(delfin/tools/adapters/orca.py, _ORCA_PARAMS) declares only `freq_temp`
(:61 region) and `_build_full_orca_input` (:196) writes only Temp — so a
spec that follows rule 9 cannot resolve through resolve_spec: the
pressure is either dropped silently or flagged as an unknown parameter.

Both directions:
- fires: a spec with freq_pressure must reach the ORCA input as
  `Pressure <value>` inside the %freq block (and resolve_spec must not
  reject it as unknown);
- silent: a spec without freq_pressure builds exactly as before.
"""

from __future__ import annotations

import pytest

from delfin.tools.platform import resolve_spec


def _freq_spec(**extra):
    return {
        "name": "freq_pressure_probe",
        "steps": [
            {"step": "smiles_to_xyz", "kwargs": {"smiles": "CO"}},
            {"step": "orca_freq",
             "kwargs": {"charge": 0, "mult": 1, "method": "PBE0",
                        "basis": "def2-TZVP", "freq_temp": 298.15,
                        **extra}},
        ],
    }


def test_freq_pressure_resolves_through_the_real_spec_path():
    """A rule-9-compliant spec must resolve — freq_pressure must not be
    rejected as an unknown parameter by the adapter contract."""
    spec = _freq_spec(freq_pressure=1.0)
    resolved = resolve_spec(spec)          # raises / returns a report — see below
    assert resolved is not None


def test_freq_pressure_reaches_the_orca_input_as_pressure():
    """The built %freq block must carry Temp AND Pressure."""
    from delfin.tools.adapters.orca import _build_full_orca_input
    text = _build_full_orca_input(
        "! PBE0 def2-TZVP FREQ", 0, 1, 4, 1000, "C 0 0 0\nO 0 0 1.2\n",
        freq_temp=298.15, freq_pressure=1.0,
    )
    assert "Temp 298.15" in text
    assert "Pressure 1.0" in text


def test_without_freq_pressure_the_input_is_unchanged():
    """Silent direction: no freq_pressure -> no Pressure line, no crash."""
    from delfin.tools.adapters.orca import _build_full_orca_input
    text = _build_full_orca_input(
        "! PBE0 def2-TZVP FREQ", 0, 1, 4, 1000, "C 0 0 0\nO 0 0 1.2\n",
        freq_temp=298.15,
    )
    assert "Temp 298.15" in text
    assert "Pressure" not in text
