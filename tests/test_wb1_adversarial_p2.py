"""Adversarial P2 tests (reviewer nacht-s15) for WB1 rules 1/2/8/9.

All probes go through the real entry point scientific_lint. Each case
states the physically right answer; a failure means the rule mis-fires
or misses.
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint


def _spec(steps):
    return {"name": "s", "steps": steps}


def _st(step, label=None, **kw):
    d = {"step": step, "kwargs": kw}
    if label:
        d["label"] = label
    return d


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


# ------------------------------------------------------------ rule 1
def test_rule1_fires_on_basis_change_mid_chain():
    spec = _spec([_st("orca_sp", label="a", method="B3LYP", basis="def2-SVP"),
                  _st("orca_sp", label="b", method="B3LYP", basis="def2-TZVP")])
    assert ("level-consistency", "b", "warning") in _rules(spec)


def test_rule1_silent_when_case_only_differs():
    # B3LYP vs b3lyp is the same method — no finding.
    spec = _spec([_st("orca_sp", label="a", method="B3LYP", basis="def2-SVP"),
                  _st("orca_sp", label="b", method="b3lyp", basis="DEF2-SVP")])
    assert "level-consistency" not in {r for r, _, _ in _rules(spec)}


def test_rule1_silent_when_both_steps_inherit():
    spec = _spec([_st("orca_sp", label="a"),
                  _st("orca_sp", label="b")])
    assert "level-consistency" not in {r for r, _, _ in _rules(spec)}


def test_rule1_silent_when_only_one_side_states_the_value():
    # 'inherited' entries are not comparable (documented design).
    spec = _spec([_st("orca_sp", label="a", method="B3LYP"),
                  _st("orca_sp", label="b", method="B3LYP", basis="def2-SVP")])
    assert "level-consistency" not in {r for r, _, _ in _rules(spec)}


# ------------------------------------------------------------ rule 2
def test_rule2_fires_for_opt_without_freq():
    spec = _spec([_st("orca_opt", label="o1"), _st("orca_sp", label="s")])
    assert ("freq-after-opt-missing", "o1", "warning") in _rules(spec)


def test_rule2_silent_for_opt_freq_pair():
    spec = _spec([_st("orca_opt", label="o1", freq_temp=350.0),
                  _st("orca_freq", label="f1", freq_temp=350.0,
                      freq_pressure=1.0)])
    rules = _rules(spec)
    assert ("freq-after-opt-missing", "o1", "warning") not in rules
    # rule 9: the freq step states temperature AND pressure -> silent
    # (o1's freq_temp sits on an orca_opt — no thermochemistry there)
    assert "thermochemistry-conditions-missing" not in {r for r, _, _ in rules}


def test_rule2_silent_for_optfreq_step():
    spec = _spec([_st("orca_opt_freq", label="of1")])
    assert "freq-after-opt-missing" not in {r for r, _, _ in _rules(spec)}


def test_rule2_ts_chain_complete_is_silent():
    spec = _spec([_st("orca_ts_opt", label="ts1"),
                  _st("orca_ts_freq_check", label="chk1"),
                  _st("orca_ts_irc", label="irc1")])
    rules = _rules(spec)
    assert "ts-freq-missing" not in {r for r, _, _ in rules}
    assert "ts-irc-missing" not in {r for r, _, _ in rules}


def test_rule2_ts_without_check_fires_both():
    spec = _spec([_st("orca_ts_opt", label="ts1"), _st("orca_sp", label="s")])
    rules = _rules(spec)
    assert ("ts-freq-missing", "ts1", "warning") in rules
    assert ("ts-irc-missing", "ts1", "warning") in rules


# ------------------------------------------------------------ rule 8
def test_rule8_fires_on_solvent_change():
    spec = _spec([_st("orca_sp", label="a", solvent="DMF", solvent_model="CPCM"),
                  _st("orca_sp", label="b", solvent="water", solvent_model="CPCM")])
    assert ("solvent-inconsistent", "b", "warning") in _rules(spec)


def test_rule8_silent_on_same_solvent():
    spec = _spec([_st("orca_sp", label="a", solvent="DMF", solvent_model="SMD"),
                  _st("orca_sp", label="b", solvent="DMF", solvent_model="SMD")])
    assert "solvent-inconsistent" not in {r for r, _, _ in _rules(spec)}


def test_rule8_fires_on_model_change_only():
    spec = _spec([_st("orca_sp", label="a", solvent="DMF", solvent_model="SMD"),
                  _st("orca_sp", label="b", solvent="DMF", solvent_model="CPCM")])
    assert ("solvent-inconsistent", "b", "warning") in _rules(spec)


# ------------------------------------------------------------ rule 9
def test_rule9_fires_for_freq_without_temp():
    spec = _spec([_st("orca_freq", label="f1")])
    assert ("thermochemistry-conditions-missing", "f1", "warning") in _rules(spec)


def test_rule9_finding_F4_freq_temp_set_but_pressure_missing():
    # Task target: temperature AND pressure explicit. With freq_temp set
    # but no pressure, the conditions of the Gibbs energy are still
    # partially implicit (ORCA default 1 atm) — rule 9 should say so.
    spec = _spec([_st("orca_freq", label="f1", freq_temp=298.15)])
    rules = _rules(spec)
    assert ("thermochemistry-conditions-missing", "f1", "warning") in rules


def test_rule9_finding_F5_optfreq_computes_thermochemistry_too():
    # orca_opt_freq produces ZPE/thermal corrections exactly like
    # orca_freq (it is in _FREQ_STEP_NAMES) — the explicit-conditions
    # rule must cover it, not only bare orca_freq.
    spec = _spec([_st("orca_opt_freq", label="of1")])
    rules = _rules(spec)
    assert ("thermochemistry-conditions-missing", "of1", "warning") in rules


# -------------------------------------------------------- robustness
def test_robustness_empty_and_malformed_steps():
    assert _rules(_spec([])) == set()
    assert _rules(_spec([None, "junk", {"step": None, "kwargs": None}])) == set()
    # branch value reused from the F1 shape must not crash the chain walk
    assert _rules({"steps": [], "branches": {"a": "junk"}}) == set()
