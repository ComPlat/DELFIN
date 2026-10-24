"""WB1 P3 tests — system / pre-DFT realism rules: dispersion correction for
non-covalent and flexible systems (5), conformer search before DFT for
flexible molecules (6), spin-state scan for 3d complexes (7).

Sources:
- Dispersion: Caldeweyher et al., J. Chem. Phys. 2019, 150, 154122,
  DOI 10.1063/1.5090222 (D4, generally applicable charge-dependent
  London-dispersion correction).
- Conformer search before DFT: Pracht, Bohle, Grimme, Phys. Chem. Chem.
  Phys. 2020, 22, 7169-7192, DOI 10.1039/C9CP06869D (CREST, automated
  low-energy conformer exploration at xTB level).
- Spin states of 3d complexes: Wappett & Goerigk, J. Chem. Theory Comput.
  2023, 19, DOI 10.1021/acs.jctc.3c00558 (MME55 benchmark — spin-state
  energetics are method-sensitive); DELFIN's own OCCUPIER spin-state
  explorer (delfin/occupier.py, CONTROL OCCUPIER=auto) and broken-symmetry
  sequence rule (delfin/config.py:709-764).
"""

from __future__ import annotations

from delfin.tools._science import scientific_lint


def _rules(spec):
    return {(f.rule, f.step, f.level) for f in scientific_lint(spec)}


def _st(step, label, **kw):
    return {"step": step, "label": label, "kwargs": kw}


# ── Rule 5: dispersion correction for non-covalent / flexible systems ────

def test_noncovalent_system_without_dispersion_is_flagged():
    out = _rules({"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP", system="noncovalent"),
    ]})
    assert ("dispersion-recommended", "opt", "warning") in out


def test_flexible_system_without_dispersion_is_flagged():
    out = _rules({"name": "s", "steps": [
        _st("orca_sp", "sp", charge=0, mult=1, method="PBE", basis="def2-SVP",
            system="flexible"),
    ]})
    assert ("dispersion-recommended", "sp", "warning") in out


def test_noncovalent_with_dispersion_is_silent():
    assert _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "opt", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP", dispersion="D4", system="noncovalent",
            freq_temp=298.15),
    ]}) == set()


def test_composite_method_counts_as_dispersed():
    # r2scan-3c includes its own dispersion (the composite rule documents it)
    assert _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "opt", charge=0, mult=1, method="r2scan-3c",
            system="noncovalent", freq_temp=298.15),
    ]}) == set()


def test_no_system_hint_means_rule_is_silent():
    # fires only on data it actually has — no invented chemistry
    assert _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "opt", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP", freq_temp=298.15),
    ]}) == set()


# ── Rule 6: conformer search before DFT for flexible molecules ───────────

def test_flexible_system_without_conformer_step_is_flagged():
    out = _rules({"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=1, system="flexible"),
    ]})
    assert ("conformer-search-missing", "opt", "warning") in out


def test_flexible_system_with_crest_step_is_silent():
    spec = {"name": "s", "steps": [
        _st("crest_conformers", "crest", charge=0, mult=1,
            system="flexible"),
        _st("orca_opt_freq", "opt", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP", dispersion="D4", system="flexible",
            freq_temp=298.15),
    ]}
    assert ("conformer-search-missing", "opt", "warning") not in _rules(spec)
    assert _rules(spec) == set()


def test_noncovalent_rigid_system_needs_no_conformer_search():
    # rule 6 (conformer search) is a FLEXIBLE-system requirement; a
    # noncovalent but rigid system needs dispersion (rule 5, present here)
    # but no conformer search
    out = _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "opt", charge=0, mult=1, method="B3LYP",
            basis="def2-SVP", dispersion="D4", system="noncovalent",
            freq_temp=298.15),
    ]})
    assert ("conformer-search-missing", "opt", "warning") not in out


# ── Rule 7: 3d complexes need a spin-state scan / broken symmetry ────────

def test_3d_complex_without_spin_scan_is_flagged():
    out = _rules({"name": "s", "steps": [
        _st("orca_opt", "opt", charge=0, mult=5, metals=["Fe"],
            system="3d_complex"),
    ]})
    assert ("spin-state-scan-missing", "opt", "warning") in out


def test_3d_complex_with_occupier_step_is_silent():
    spec = {"name": "s", "steps": [
        _st("occupier", "spin", charge=0, mult=5, metals=["Fe"],
            system="3d_complex"),
        _st("orca_opt_freq", "opt", charge=0, mult=5, metals=["Fe"],
            system="3d_complex", freq_temp=298.15),
    ]}
    assert ("spin-state-scan-missing", "opt", "warning") not in _rules(spec)
    assert _rules(spec) == set()


def test_3d_complex_with_broken_symmetry_flag_is_silent():
    # an explicit broken_symmetry kwarg on the DFT step states the agent
    # considered open-shell coupling
    assert _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "opt", charge=0, mult=5, metals=["Fe"],
            system="3d_complex", broken_symmetry=True, freq_temp=298.15),
    ]}) == set()


def test_nonmetal_system_needs_no_spin_scan():
    # the 3d_complex hint WITH metals drives rule 7; the same hint on a
    # metal-free step names the system, but rule 7 reads the metals kwarg
    # as the physical fact and stays silent without one
    assert _rules({"name": "s", "steps": [
        _st("orca_opt_freq", "opt", charge=0, mult=1,
            system="3d_complex", freq_temp=298.15),
    ]}) == set()
