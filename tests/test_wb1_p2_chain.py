"""WB1 P2 tests — chain-consistency rules: one level of theory per energy
difference (1), frequencies after every optimisation / TS+IRC (2), one
solvent model + solvent through the chain (8), explicit T/P for
thermochemistry (9, WARNING per QS ruling s29).

Sources:
- ORCA 6.1.1 manual §4.7 Thermochemistry: thermochemical analysis assumes
  ideal-gas statistical mechanics at explicit T, P (defaults 298.15 K,
  1 atm) — https://orca-manual.mpi-muelheim.mpg.de/contents/
  structurereactivity/thermochemistry.html
- ZPE/thermal corrections valid on the geometry optimised at the same
  level (J. Chem. Theory Comput. 2021, 17, 10.1021/acs.jctc.1c00474;
  DELFIN's own same-level rule, wave-20 "compare methods, not programs").
- TS needs exactly 1 imaginary mode plus IRC to both minima: DELFIN's own
  TS adapters, delfin/tools/adapters/orca_ts_adapters.py (orca_ts_opt,
  orca_ts_freq_check, orca_ts_irc) — the IRC adapter fails closed when no
  endpoints are reached (orca_ts_adapters.py:294-298).
- Freq adapter accepts explicit conditions via `freq_temp`
  (delfin/tools/adapters/orca.py:196-197,347).
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


def _st(step, label, **kw):
    return {"step": step, "label": label, "kwargs": kw}


# ── Rule 1: one level of theory per energy difference ────────────────────

def test_freq_at_different_level_than_opt_is_flagged():
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP"),
        _st("orca_freq", "freq", charge=0, mult=1, method="B3LYP",
            basis="def2-TZVP"),
    ]}
    assert ("level-consistency", "freq", "warning") in _rules(spec)


def test_freq_at_same_level_as_opt_is_silent():
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP"),
        _st("orca_freq", "freq", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP", freq_temp=298.15),
    ]}
    assert _rules(spec) == set()


def test_inherited_level_is_not_flagged():
    # steps that omit method/basis (inherit at runtime) carry no concrete
    # level — nothing to compare, rule stays silent; opt_freq covers its
    # own frequency requirement
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, method="B3LYP", basis="def2-SVP"),
        _st("orca_opt_freq", "of", charge=0),
    ]}
    assert _rules(spec) == set()


def test_opt_freq_counts_as_its_own_level():
    # the combined adapter sets the level itself; a following sp at another
    # level is flagged, at the same level is not
    spec = {"name": "s", "steps": [
        _st("orca_opt_freq", "of", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP"),
        _st("orca_sp", "sp", charge=0, mult=1, method="B3LYP",
            basis="def2-TZVP"),
    ]}
    assert ("level-consistency", "sp", "warning") in _rules(spec)


# ── Rule 2: frequencies after every optimisation; TS needs freq + IRC ────

def test_opt_without_freq_is_flagged():
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1),
    ]}
    assert ("freq-after-opt-missing", "opt", "warning") in _rules(spec)


def test_opt_followed_by_freq_is_silent():
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1),
        _st("orca_freq", "freq", charge=0, mult=1, freq_temp=298.15),
    ]}
    assert _rules(spec) == set()


def test_combined_opt_freq_is_silent():
    assert _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "of", charge=0, mult=1)]}) == set()


def test_ts_opt_without_freq_check_and_irc_is_flagged():
    spec = {"name": "s", "steps": [
        _st("orca_ts_opt", "ts", charge=0, mult=1),
    ]}
    out = _rules(spec)
    assert ("ts-freq-missing", "ts", "warning") in out
    assert ("ts-irc-missing", "ts", "warning") in out


def test_ts_with_freq_check_but_without_irc_is_flagged_once():
    spec = {"name": "s", "steps": [
        _st("orca_ts_opt", "ts", charge=0, mult=1),
        _st("orca_ts_freq_check", "tscheck", charge=0, mult=1),
    ]}
    out = _rules(spec)
    assert ("ts-irc-missing", "ts", "warning") in out
    assert ("ts-freq-missing", "ts", "warning") not in out


def test_ts_full_chain_is_silent():
    spec = {"name": "s", "steps": [
        _st("orca_ts_opt", "ts", charge=0, mult=1),
        _st("orca_ts_freq_check", "tscheck", charge=0, mult=1,
            freq_temp=298.15),
        _st("orca_ts_irc", "irc", charge=0, mult=1),
    ]}
    assert _rules(spec) == set()


def test_branches_are_checked_independently():
    # a branch missing its freq is flagged even when the trunk has one
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1),
        _st("orca_freq", "freq", charge=0, mult=1),
    ], "branches": {
        "ox": [_st("orca_opt", "ox_opt", charge=1, mult=2)],
    }}
    out = _rules(spec)
    assert ("freq-after-opt-missing", "opt", "warning") not in out
    assert ("freq-after-opt-missing", "ox_opt", "warning") in out


# ── Rule 8: one solvent model + solvent through the chain ────────────────

def test_two_different_solvents_in_one_chain_is_flagged():
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, solvent="water"),
        _st("orca_sp", "sp", charge=0, mult=1, solvent="acetonitrile"),
    ]}
    assert ("solvent-inconsistent", "sp", "warning") in _rules(spec)


def test_one_solvent_through_the_chain_is_silent():
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, solvent="water"),
        _st("orca_opt_freq", "of", charge=0, mult=1, solvent="water"),
    ]}
    assert _rules(spec) == set()


def test_missing_solvent_in_a_step_is_not_a_mismatch():
    # a step that omits solvent inherits it at runtime; absence is not a
    # different solvent
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, solvent="water"),
        _st("orca_opt_freq", "of", charge=0, mult=1),
    ]}
    assert _rules(spec) == set()


def test_gas_phase_chain_is_silent():
    assert _rules({"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, solvent="water"),
        _st("orca_freq", "freq", charge=0, mult=1, solvent="water",
            freq_temp=298.15),
    ]}) == set()


def test_solvent_model_mismatch_is_flagged():
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, solvent="water",
            solvent_model="SMD"),
        _st("orca_sp", "sp", charge=0, mult=1, solvent="water",
            solvent_model="CPCM"),
    ]}
    assert ("solvent-inconsistent", "sp", "warning") in _rules(spec)


# ── Rule 9: explicit T/P for thermochemistry (warning, QS s29) ──────────

def test_freq_without_explicit_temperature_is_flagged():
    spec = {"name": "s", "steps": [
        _st("orca_freq", "freq", charge=0, mult=1),
    ]}
    out = _rules(spec)
    assert ("thermochemistry-conditions-missing", "freq", "warning") in out


def test_freq_with_explicit_temperature_is_silent():
    spec = {"name": "s", "steps": [
        _st("orca_freq", "freq", charge=0, mult=1, freq_temp=373.15),
    ]}
    assert _rules(spec) == set()


def test_no_freq_step_means_no_thermochemistry_warning():
    spec = {"name": "s", "steps": [
        _st("orca_sp", "sp", charge=0, mult=1),
    ]}
    assert _rules(spec) == set()
