"""The simulator's response surface: one definition, and its known optimum.

The loop is only as honest as the surface it searches, so the surface is
pinned as data here rather than assumed:

* the scalar and the vectorised form agree (two implementations of one
  function is how an optimiser gets measured against a different function
  than the device answers with);
* the mean temperature is below the setpoint and rises with the duration —
  that coupling is WHY the search is two-dimensional, and without it these
  tests would pass on a surface that is one axis wearing two names;
* the known maxima of this illustrative surface are pinned as numbers, so a
  later change to a constant is a visible change to the optimum, not a
  silently different search problem.

The numbers below were computed from the constants in the module, not read
back out of a run of the module: they are the surface's own arithmetic.
"""

from __future__ import annotations

import numpy as np
import pytest


def test_the_scalar_and_the_vectorised_surface_agree():
    """One surface, two forms: a disagreement measures the wrong function."""
    from delfin.lab import surface
    points = [(30.0, 5.0), (62.0, 25.0), (63.75, 25.0), (110.0, 85.0),
              (20.0, 0.0), (60.0, 45.0)]
    want = [surface.yield_percent(sp, du) for sp, du in points]
    got = surface.vectorised(points)
    assert np.allclose(got, want, atol=1e-12), (got, want)


def test_a_zero_length_run_made_nothing_and_tasted_only_the_ambient():
    from delfin.lab import surface
    assert surface.yield_percent(80.0, 0.0) == 0.0
    assert surface.mean_temperature(80.0, 0.0) == surface.AMBIENT_C


def test_the_mean_temperature_is_below_the_setpoint_and_rises_with_time():
    """The coupling that makes the search two-dimensional.

    A short run at a high setpoint sees a much colder mean; that is why
    duration matters at all. If this ever stopped holding, the optimiser
    tests would be measuring a one-dimensional problem.
    """
    from delfin.lab import surface
    short = surface.mean_temperature(90.0, 2.0)
    mid = surface.mean_temperature(90.0, 10.0)
    long = surface.mean_temperature(90.0, 60.0)
    assert surface.AMBIENT_C < short < mid < long < 90.0, (short, mid, long)
    assert long > 85.0, "a long run must approach its setpoint"


def test_the_yield_is_a_percentage_in_range():
    from delfin.lab import surface
    lo, hi = (30.0, 5.0), (110.0, 85.0)
    values = surface.vectorised(surface.lattice(lo, hi, (5.0, 5.0)))
    assert values.min() >= 0.0 and values.max() <= 100.0, (values.min(),
                                                           values.max())


def test_too_cold_and_too_long_both_cost_yield():
    """A control for the surface being an optimum at all, not a ramp."""
    from delfin.lab import surface
    at_optimum = surface.yield_percent(64.0, 25.0)
    too_cold = surface.yield_percent(25.0, 25.0)
    too_hot = surface.yield_percent(110.0, 25.0)
    too_brief = surface.yield_percent(64.0, 5.0)
    too_long = surface.yield_percent(64.0, 85.0)
    for name, other in (("cold", too_cold), ("hot", too_hot),
                        ("brief", too_brief), ("long", too_long)):
        assert other < at_optimum, (name, other, at_optimum)


def test_the_known_maxima_are_what_the_constants_imply():
    """The surface's own arithmetic, pinned.

    Continuous: the best of a 0.05 degC / 0.05 min sweep, to three decimals.
    Lattice: the best of the shipped 10 degC / 10 min envelope of the demo
    manifest. Both are properties of the constants above; a change to a
    constant is expected to move them, and this test is what makes that
    visible instead of silent.
    """
    from delfin.lab import surface
    assert surface.yield_percent(63.75, 25.0) == pytest.approx(98.6164,
                                                               abs=5e-4)
    lattice = surface.lattice((30.0, 5.0), (110.0, 85.0), (10.0, 10.0))
    assert len(lattice) == 81
    best, at = surface.best_on(lattice)
    assert best == pytest.approx(94.2985, abs=5e-4), best
    assert at == (60.0, 45.0), at


def test_the_best_of_a_lattice_is_not_the_best_of_the_surface():
    """Why a grid needs a budget: the lattice cannot reach the ridge.

    This is the fact the round-count comparison rests on. If a lattice point
    could reach the continuous maximum, "rounds to 99 % of the lattice best"
    would be measuring almost nothing.
    """
    from delfin.lab import surface
    lattice = surface.lattice((30.0, 5.0), (110.0, 85.0), (10.0, 10.0))
    lattice_best, _ = surface.best_on(lattice)
    continuous = surface.yield_percent(63.75, 25.0)
    assert continuous - lattice_best > 4.0, (continuous, lattice_best)


def test_the_lattice_is_built_from_the_envelope_it_is_given():
    """A lattice must not contain points the kernel would refuse."""
    from delfin.lab import surface
    lattice = surface.lattice((30.0, 5.0), (110.0, 85.0), (10.0, 10.0))
    assert lattice[:, 0].min() == 30.0 and lattice[:, 0].max() == 110.0
    assert lattice[:, 1].min() == 5.0 and lattice[:, 1].max() == 85.0
    narrow = surface.lattice((50.0, 10.0), (70.0, 30.0), (10.0, 10.0))
    assert sorted({(a, b) for a, b in narrow}) == [
        (50.0, 10.0), (50.0, 20.0), (50.0, 30.0),
        (60.0, 10.0), (60.0, 20.0), (60.0, 30.0),
        (70.0, 10.0), (70.0, 20.0), (70.0, 30.0)]


def test_the_surface_needs_no_dependency_outside_numpy_and_the_stdlib():
    """A closed loop must not add a licence obligation to run.

    numpy and scipy are already core requirements of DELFIN, so neither is a
    new obligation; but scipy is not installed on every installation that
    runs the lab kernel, so the surface and the loop may not import it.
    Importing scipy here would fail on a bare install rather than in CI.
    """
    import ast
    from pathlib import Path
    root = Path(__file__).resolve().parents[1] / "delfin" / "lab"
    for name in ("surface.py", "loop.py"):
        tree = ast.parse((root / name).read_text(encoding="utf-8"))
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 \
                    and node.module:
                imported.add(node.module.split(".")[0])
        assert imported <= {"numpy", "math", "typing", "__future__"}, (name,
                                                                      imported)
