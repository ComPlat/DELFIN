"""A figure a tool returned is still that figure when the answer rounds it.

Input: a quantity claim in the answer ("2.6 kcal/mol") and the numbers
this turn's tools returned. Before: a claim was grounded only within
max(|v|*1e-4, 5e-3) of a returned number or of a difference of two, so
a value printed to fewer decimals than the tool printed was flagged as
unsupported. Measured on the GFN2-xTB butane run (relative energies of
a 25-point dihedral scan, printed by a script to three decimals): of the
correct values rounded to one decimal, 18.5% were grounded; the final
table of a correct answer carried a caveat on most of its rows.

After: a claim printed to d decimals (d >= 1) is grounded when a returned
number rounds to it, i.e. lies within half a unit of its last digit. It
applies to returned numbers only, not to differences of two, and not to
integers, whose half unit (0.5) is wide enough to ground most invented
values. Measured cost on the same pools, random claims in [0.1, 10]
kcal/mol: one decimal 7.2% -> 19.6% grounded by coincidence (focused
pool, 62 numbers), two decimals 4.8% -> 5.0%; three decimals unchanged,
since their half unit is below the old tolerance.
"""

from __future__ import annotations

import random

from delfin.agent import verify_guard as vg

# What a script printing the relative energies returns, plus the
# total energies a grep returns.
POOL = [-13.664174370972, -13.665127387246, -13.661036183827,
        0.0, 0.598, 2.571, 4.932, 0.711, 2.618]


def _flagged(text, pool=POOL):
    return [f.quantity for f in
            vg.scan_for_unsourced_quantities(text, numbers=list(pool))]


def test_a_value_rounded_to_fewer_decimals_is_grounded():
    assert _flagged("The eclipsed barrier is 2.6 kcal/mol.") == []
    assert _flagged("The syn barrier is 4.9 kcal/mol.") == []
    assert _flagged("gauche lies 0.60 kcal/mol above anti.") == []


def test_a_rounding_that_does_not_hold_is_still_flagged():
    # 2.571 rounds to 2.6, never to 2.5; 4.932 to 4.93, never to 4.94.
    assert _flagged("The eclipsed barrier is 2.5 kcal/mol.") == [
        "2.5 kcal/mol"]
    assert _flagged("The syn barrier is 4.94 kcal/mol.") == [
        "4.94 kcal/mol"]


def test_an_integer_is_held_to_the_old_tolerance():
    # "5 kcal/mol" would otherwise be grounded by anything in [4.5, 5.5].
    assert _flagged("The syn barrier is 5 kcal/mol.") == ["5 kcal/mol"]


def test_rounding_does_not_reach_differences():
    # 3.07 is a difference of two returned numbers and is grounded as
    # before; its rounding 3.1 is not, because the half-unit rule applies
    # to returned numbers only.
    pool = [10.0, 6.93]                  # 10.0 - 6.93 = 3.07
    assert _flagged("The gap is 3.1 kcal/mol.", pool) == ["3.1 kcal/mol"]
    assert _flagged("The gap is 3.07 kcal/mol.", pool) == []


def test_a_restated_unit_accepts_a_rounded_partner():
    # 2.6 kcal/mol is the rounded 2.571; its kJ/mol restatement rides on it.
    assert _flagged("The barrier is 2.6 kcal/mol (10.9 kJ/mol).") == []


def test_three_decimals_keep_their_false_grounding_rate():
    rng = random.Random(20261022)
    pool = [rng.uniform(-20, 20) for _ in range(60)]
    grounded = sum(
        _flagged(f"The value is {rng.uniform(0.01, 10.0):.3f} eV.", pool) == []
        for _ in range(3000))
    rng = random.Random(20261022)
    pool = [rng.uniform(-20, 20) for _ in range(60)]
    strict = vg._claim_is_observed
    try:
        vg._claim_is_observed = lambda q, p: any(
            vg._grounded_in_observations(v, p)
            for v, _ in vg._claim_readings(q))
        before = sum(
            _flagged(f"The value is {rng.uniform(0.01, 10.0):.3f} eV.",
                     pool) == []
            for _ in range(3000))
    finally:
        vg._claim_is_observed = strict
    assert grounded == before
