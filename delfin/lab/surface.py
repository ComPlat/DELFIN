"""The simulated reactor's response surface: yield over two conditions.

Input: a reactor setpoint in degC and a run duration in minutes. Output: the
yield in percent that the simulated reactor reports, as a number in
``[0, 100]``.

This module is the ONE definition of the surface. The simulator driver calls
it to answer a ``yield`` reading, and a loop's report calls it to say what
the surface would give at a point nothing was run at. Two definitions would
let the optimiser be measured against a different function than the device
answers with, which is the mistake that makes a closed-loop number
meaningless.

The surface is ILLUSTRATIVE, not measured::

    Tbar  = setpoint + (ambient - setpoint) * (1 - e**-x) / x
            x = duration * 60 / tau,      tau = 60 s, ambient = 20 degC
    yield = 100 * exp(-((Tbar - T_opt) / T_width)**2)
                * (1 - e**-k_form * duration) * e**-k_decay * duration

``Tbar`` is the analytic mean of the first-order thermal lag the simulated
hotplate already uses (``delfin/lab/drivers/simulator.py``, tau 60 s): a
short run sees a mean temperature well below its setpoint, and that is what
couples the two axes — without it the optimum would be a setpoint alone and
the two-dimensional search would be a one-dimensional one in disguise.

The three factors say, in order: the reaction has a thermal optimum; the
product forms towards completion first-order; the product degrades slowly
while it sits at temperature.

Discovered values (a search's job, and the numbers the tests measure against
rather than assume): the continuous maximum is 98.6164 % at 63.75 degC and
25.00 min. On the 10 degC / 10 min lattice of the shipped envelope
(30..110 degC, 5..85 min, 81 points) the best reachable point is 94.2985 %
at 60 degC and 45 min, because no lattice point sits on the ridge.

No new dependency: numpy is a core requirement of DELFIN (pyproject.toml,
``numpy>=1.26,<2.3``, BSD-3-Clause), so using it here adds no licence
obligation that the package does not already carry. The arithmetic would
also fit in ``math``; numpy is used because the optimiser in
:mod:`delfin.lab.loop` needs arrays and one import is honest about that.
"""

from __future__ import annotations

import math

import numpy as np

#: Ambient temperature the reactor starts from and cools to, in degC. The
#: same value the simulator's hotplate uses.
AMBIENT_C = 20.0

#: Thermal time constant of the mean-temperature average, in seconds. The
#: simulator's hotplate lag (TAU_S) — stated once here so the mean a run saw
#: is the same lag the device actually has.
TAU_S = 60.0

#: Temperature of the reaction's optimum, in degC.
T_OPT_C = 62.0

#: Width of the thermal window, in degC: yield falls to 1/e at +-15 degC.
T_WIDTH_C = 15.0

#: First-order rate of product formation, per minute.
K_FORM = 0.25

#: First-order rate of product degradation at temperature, per minute.
K_DECAY = 0.00048

#: The condition axes a run of this surface takes, as (channel, unit) pairs.
AXES: tuple[tuple[str, str], ...] = (("setpoint", "degC"), ("duration", "min"))


def mean_temperature(setpoint_c: float, duration_min: float) -> float:
    """The mean temperature the reaction saw, in degC.

    The analytic mean of ``T(t) = sp + (ambient - sp) * exp(-t/tau)`` over
    ``0 <= t <= duration``. At zero elapsed time the mean is the ambient
    temperature (nothing has been heated yet), not the setpoint.
    """
    sp = float(setpoint_c)
    x = max(float(duration_min), 0.0) * 60.0 / TAU_S
    if x <= 0.0:
        return AMBIENT_C
    return sp + (AMBIENT_C - sp) * (1.0 - math.exp(-x)) / x


def yield_percent(setpoint_c: float, duration_min: float) -> float:
    """The yield in percent at one condition point (see the module docstring)."""
    tbar = mean_temperature(setpoint_c, duration_min)
    duration = max(float(duration_min), 0.0)
    thermal = math.exp(-((tbar - T_OPT_C) / T_WIDTH_C) ** 2)
    formed = 1.0 - math.exp(-K_FORM * duration)
    degraded = math.exp(-K_DECAY * duration)
    return float(100.0 * thermal * formed * degraded)


def vectorised(points) -> np.ndarray:
    """``yield_percent`` over an ``(n, 2)`` array of (setpoint, duration).

    Same arithmetic as the scalar function, in numpy: the optimiser scores
    thousands of candidates and must score them with the surface the device
    answers with.
    """
    arr = np.asarray(points, dtype=float).reshape(-1, 2)
    sp = arr[:, 0]
    duration = np.maximum(arr[:, 1], 0.0)
    x = duration * 60.0 / TAU_S
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.where(
            x > 0.0,
            sp + (AMBIENT_C - sp) * (1.0 - np.exp(-x)) / np.where(x > 0.0, x, 1.0),
            AMBIENT_C,
        )
    thermal = np.exp(-((mean - T_OPT_C) / T_WIDTH_C) ** 2)
    formed = 1.0 - np.exp(-K_FORM * duration)
    degraded = np.exp(-K_DECAY * duration)
    return 100.0 * thermal * formed * degraded


def lattice(lo: tuple[float, float], hi: tuple[float, float],
            step: tuple[float, float]) -> np.ndarray:
    """The candidate lattice of an envelope, as an ``(n, 2)`` array.

    Built from the channel's own declared min/max, so the points a search may
    propose are exactly the ones the kernel's envelope admits: a candidate
    outside the manifest would be refused, and a search that proposes it has
    measured the refusal rather than the surface.
    """
    first = np.arange(float(lo[0]), float(hi[0]) + 1e-9, float(step[0]))
    second = np.arange(float(lo[1]), float(hi[1]) + 1e-9, float(step[1]))
    return np.array([[a, b] for a in first for b in second], dtype=float)


def best_on(points) -> tuple[float, tuple[float, float]]:
    """``(best yield, condition point)`` over the given candidates."""
    arr = np.asarray(points, dtype=float).reshape(-1, 2)
    values = vectorised(arr)
    i = int(np.argmax(values))
    return float(values[i]), (float(arr[i, 0]), float(arr[i, 1]))
