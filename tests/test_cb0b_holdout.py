"""Red control for CB0B phase 1: the workflow holdout sampler.

The benchmark's workflow families must draw chemist-voiced instances whose
prompt never names the workflow behind them, from a molecule pool that stays
small (<= 20 atoms) and covers shell and charge variety, deterministically
from a seed.

This file was written BEFORE delfin/agent/bench_holdout.py existed: every
test here is red as ImportError on the previous commit.
"""
from __future__ import annotations

import re

import pytest

from delfin.agent import bench_holdout as bh

# Deliberately fixed and public: tests and practice use only this seed.
TEST_SEED = 20261020

# Words that name a DELFIN workflow. A prompt that contains one leaks the
# answer to the agent before it has thought.
WORKFLOW_WORDS = re.compile(
    r"(?i)\b(classic|redox|esd|excited.state.workfl|occupier|crest"
    r"|conformer.screen|screening.workfl|imag.mode.fix)\b")


def test_the_public_seed_gives_the_same_instances_twice():
    a = bh.sample(TEST_SEED, 6)
    b = bh.sample(TEST_SEED, 6)
    assert a == b


def test_two_different_seeds_give_different_instances():
    assert bh.sample(TEST_SEED, 6) != bh.sample(TEST_SEED + 1, 6)


def test_n_instances_are_drawn_and_unique():
    got = bh.sample(TEST_SEED, 5)
    assert len(got) == 5
    assert len({i.family_id for i in got}) == len(got)


def test_no_prompt_names_the_workflow_behind_it():
    for inst in bh.sample(TEST_SEED, 7):
        hit = WORKFLOW_WORDS.search(inst.prompt)
        assert not hit, f"{inst.family_id} leaked {hit.group(0)!r}"


def test_every_molecule_stays_under_the_atom_budget():
    for fam in bh.FAMILIES:
        for mol in fam.molecules:
            assert 1 <= mol.n_atoms <= 20, (fam.id, mol.label)


def test_the_pool_covers_both_shells_and_both_charges():
    closed = any(m.multiplicity == 1 for f in bh.FAMILIES
                 for m in f.molecules)
    open_ = any(m.multiplicity > 1 for f in bh.FAMILIES
                for m in f.molecules)
    charged = any(m.charge != 0 for f in bh.FAMILIES
                  for m in f.molecules)
    neutral = any(m.charge == 0 for f in bh.FAMILIES
                  for m in f.molecules)
    assert closed and open_ and charged and neutral


def test_seven_families_one_per_delphin_workflow():
    assert {f.id for f in bh.FAMILIES} == {
        "classic_opt_freq", "redox", "excited_states", "imaginary_mode",
        "spin_states", "conformer_prestep", "screening"}


def test_the_real_seed_is_read_from_outside_the_repository():
    """sample() must consult a seed file the person creates OUTSIDE the
    repo; sessions never read it. With the file absent, the public test
    seed must be used and the record must say so."""
    src = bh.real_seed.__doc__ or ""
    assert "outside" in src.lower()
    # With no seed file present the default is the public seed:
    assert bh.real_seed() == bh.PUBLIC_SEED


def test_each_instance_records_its_grading_levels():
    for inst in bh.sample(TEST_SEED, 7):
        assert set(inst.grading_levels) <= {"consistency",
                                            "reproduction",
                                            "reference", "trap"}
        assert inst.grading_levels, inst.family_id


def test_a_charged_open_shell_instance_carries_charge_and_multiplicity():
    inst = next(i for i in bh.sample(TEST_SEED, 7)
                if i.molecule.charge != 0 or i.molecule.multiplicity > 1)
    assert inst.prompt  # the request states what property is wanted


@pytest.mark.parametrize("n", [0, 1, 30])
def test_sample_handles_any_n(n):
    assert len(bh.sample(TEST_SEED, n)) == n
