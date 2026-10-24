"""WB1 adversarial-fix tests for reviewer findings F4/F5 (rule 9 coverage).

F4: thermochemistry conditions are T AND P; with freq_temp set but no
explicit pressure the conditions are still implicit (ORCA default 1 atm).
F5: orca_opt_freq computes thermochemistry exactly like orca_freq
(it is in chain._FREQ_STEP_NAMES) — the conditions check must cover it.

Sources: ORCA 6.1.1 manual §4.7 (thermochemistry at explicit T and P,
defaults 298.15 K / 1 atm); orca.py:196-197 writes %freq Temp (and the
adapter's knob set), orca.py:440-447 (orca_opt_freq = "OPT FREQ").
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


def _st(step, label, **kw):
    return {"step": step, "label": label, "kwargs": kw}


# ── F4: freq_temp without pressure is still an incomplete statement ─────

def test_freq_temp_without_pressure_is_flagged():
    out = _rules({"name": "s", "steps": [
        _st("orca_freq", "freq", charge=0, mult=1, freq_temp=298.15),
    ]})
    assert ("thermochemistry-conditions-missing", "freq", "warning") in out


def test_freq_temp_with_pressure_is_silent():
    assert _rules({"name": "s", "steps": [
        _st("orca_freq", "freq", charge=0, mult=1, freq_temp=298.15,
            freq_pressure=1.0),
    ]}) == set()


# ── F5: orca_opt_freq carries the same thermochemistry obligation ───────

def test_opt_freq_without_temp_or_pressure_is_flagged():
    out = _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "of", charge=0, mult=1),
    ]})
    assert ("thermochemistry-conditions-missing", "of", "warning") in out


def test_opt_freq_with_both_is_silent():
    assert _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "of", charge=0, mult=1, freq_temp=298.15,
            freq_pressure=1.0),
    ]}) == set()


def test_freq_temp_with_pressure_fully_conditions_the_chain():
    # the P2 silent fixture now fully conditioned: stays silent
    assert _rules({"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1),
        _st("orca_freq", "freq", charge=0, mult=1, freq_temp=298.15,
            freq_pressure=1.0),
    ]}) == set()
