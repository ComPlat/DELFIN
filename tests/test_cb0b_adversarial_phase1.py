"""Adversarial review tests for CB0B phase 1 (bench_holdout), by nacht-s15.

Written against the builder's commit 3fb331c5 WITHOUT changing its code.
Findings (red here) are reported to nacht-s13; green guards are committed.
"""
from __future__ import annotations

import re

from delfin.agent import bench_holdout as bh

TEST_SEED = 20261020


def _electrons(mol: bh.Molecule) -> int:
    import re as _re
    counts = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15,
              "S": 16, "Cl": 17, "Fe": 26}
    total = -mol.charge
    for symbol, n in _re.findall(r"([A-Z][a-z]?)(\d*)", mol.formula):
        if symbol:
            total += counts[symbol] * (int(n) if n else 1)
    return total


# ── green guards (pin correct behaviour that could regress) ──────────


def test_same_seed_repeats_over_a_full_family_walk():
    assert bh.sample(TEST_SEED, 14) == bh.sample(TEST_SEED, 14)


def test_every_drawn_solvent_is_declared_by_its_family():
    """sample() must draw from the family's own declared solvent set, not
    a module-global one that could drift away from it."""
    for inst in bh.sample(TEST_SEED, 21):
        fam = next(f for f in bh.FAMILIES if f.id == inst.family_id)
        assert inst.solvent in fam.solvents, inst.instance_id


def test_instance_ids_are_unique_within_one_draw_of_seven():
    got = bh.sample(TEST_SEED, 7)
    assert len({i.instance_id for i in got}) == 7


def test_pool_electron_counts_allow_the_declared_multiplicities():
    """Brief item 5a at pool level. Convention after fix 773a179c: the
    redox trap carrier carries a DELIBERATELY inconsistent multiplicity
    request (the agent must catch the mismatch), flagged by trap_note;
    every non-trap molecule must be internally consistent."""
    for mol in bh._POOL:
        if bh.trap_note("redox", mol):
            # The carrier's multiplicity request must contradict its
            # electron count (11e -> doublet would be consistent, so
            # verify locally instead of a private helper).
            ne = _electrons(mol)
            unpaired = mol.multiplicity - 1
            assert not (unpaired >= 0 and (ne - unpaired) % 2 == 0), (
                f"trap carrier {mol.label} is not inconsistent")
        else:
            ne = _electrons(mol)
            if ne % 2 == 1:
                assert mol.multiplicity % 2 == 0, (mol.label, ne)
            else:
                assert mol.multiplicity % 2 == 1, (mol.multiplicity, ne)


# ── red findings (sent to the builder, not fixed by me) ──────────────


def test_finding_screening_must_not_screen_the_gas_phase():
    """bench_holdout.py:217 draws screening solvents from the global
    _SOLVENTS incl. 'gas phase'. 'Screen the gas phase as medium' is
    chemically undefined: no dielectric, no solvation free energy."""
    for inst in bh.sample(TEST_SEED, 28):
        if inst.family_id == "screening":
            assert inst.solvent != "gas phase", inst.instance_id


def test_trap_levels_only_on_instances_a_trap_carries():
    """Old F2, restated for the fix convention: the trap level on an
    instance must correspond to a trap the molecule can actually
    present via trap_note (geometry-based traps apart)."""
    for inst in bh.sample(TEST_SEED, 28):
        if "trap" in inst.grading_levels and inst.family_id == "redox":
            assert bh.trap_note("redox", inst.molecule), inst.instance_id
