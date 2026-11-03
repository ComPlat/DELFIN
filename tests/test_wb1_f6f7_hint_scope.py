"""Reviewer findings F6/F7 (s15 @ bb08f950): hint registration scope.

A real workflow states the system hint on the DFT step, not on the
preparatory step. The search/scan step itself may carry no hint — but it
still closes the pending requirement. Both directions, through the real
entry scientific_lint.
"""
from delfin.tools._science import scientific_lint


def _st(name, sid, charge=0, mult=1, **kw):
    st = {"name": name, "id": sid, "charge": charge, "mult": mult,
          "kwargs": kw}
    return st


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


def test_search_step_without_hint_still_closes_conformer_requirement():
    # F6: hint lives on the DFT step; the crest step carries none
    spec = {"name": "s", "steps": [
        _st("crest_conformers", "cc"),
        _st("orca_opt", "opt", system="flexible",
            freq_temp=298.15, freq_pressure=1.0, dispersion="D4"),
        _st("orca_opt_freq", "of", freq_temp=298.15, freq_pressure=1.0),
    ]}
    assert _rules(spec) == set()


def test_scan_step_without_hint_still_closes_spin_requirement():
    # F7: hint + metals live on the DFT step; the occupier step carries none
    spec = {"name": "s", "steps": [
        _st("occupier", "occ"),
        _st("orca_sp", "sp", system="3d_complex", metals="Fe",
            dispersion="D4"),
    ]}
    assert _rules(spec) == set()


def test_branch_step_with_hint_sees_trunk_search():
    # the same root cause through branches: the search sits in the trunk,
    # the hinted step in a branch — the requirement must still close
    spec = {"name": "s", "steps": [
        _st("crest_conformers", "cc"),
        _st("orca_opt", "opt", system="flexible",
            freq_temp=298.15, freq_pressure=1.0, dispersion="D4"),
        _st("orca_freq", "fr", freq_temp=298.15, freq_pressure=1.0),
    ], "branches": {"tight": [
        _st("orca_sp", "sp2", system="flexible", dispersion="D4"),
    ]}}
    assert _rules(spec) == set()


def test_hint_with_no_search_anywhere_still_fires():
    spec = {"name": "s", "steps": [
        _st("orca_opt", "opt", system="flexible",
            freq_temp=298.15, freq_pressure=1.0, dispersion="D4"),
    ]}
    assert ("conformer-search-missing", "orca_opt", "warning") in _rules(spec)


def test_hint_with_no_scan_anywhere_still_fires():
    spec = {"name": "s", "steps": [
        _st("orca_sp", "sp", system="3d_complex", metals=["Fe"],
            dispersion="D4"),
    ]}
    assert ("spin-state-scan-missing", "orca_sp", "warning") in _rules(spec)
