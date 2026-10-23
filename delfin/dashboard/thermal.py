"""How long a barrier takes, and how hot it has to be -- Eyring, both ways.

Two functions and the four numbers behind them, in a module of their own
because two parts of the dashboard ask the same question about the same
quantity and neither should have to import the other to do it.

The structure editor asks it about a gesture: the user is holding an atom, a
temperature is set, and the question is what that temperature will pay for
before the hand has finished moving.  The reaction graph asks it about a
network that took a week to build: given every barrier in it, which routes
are open at 298 K and which want 480.  It is the same arithmetic, and it was
written once -- in the editor, next to the drag it was written for, where a
module that has no business importing ipywidgets could not reach it.

Nothing here knows about a widget, a viewer or a structure.  It takes numbers
and returns numbers, which is why the editor's own name for each of these is
still the name it had: ``structure_editor`` imports them and re-exports them,
so every caller and every test that was written against it is untouched.
"""

from __future__ import annotations

import math

#: How long a structure is given to do something, in seconds.
#:
#: "Possible at this temperature" has two halves and this is the second: a
#: barrier is not possible or impossible, it is a waiting time.  At 298 K one
#: of 20 kcal/mol takes a minute and one of 25 takes four days, so the ceiling
#: depends on how long anyone is watching -- 17.5 kcal/mol for a second, 22.3
#: for an hour, 27.7 for a year.
#:
#: It was a control, and it did not earn the place.  Between a second and a
#: year the ceiling moves by ten kcal/mol, while the distance between
#: chemistry and nonsense is twenty against a hundred and more -- so the knob
#: changed the answer far less than it cost to understand, and it was asked
#: about twice.  An hour is a reaction set up in a flask, which is what this
#: is for.
#:
#: Asked again once a refusal started answering in both directions, and the
#: answer came out the same.  What a window control would be for is a user
#: saying "I would wait longer than that", and the two sentences a refusal now
#: carries answer that without anything being set: the waiting time is quoted
#: at the temperature that *is* set, so a barrier the hour refuses already
#: reads as "about 4 d", "about 50 years" or "longer than the universe has
#: existed", and the temperature it wants is quoted for the window.  A knob
#: that only moved the line between those two would be a third way of asking a
#: question now answered twice over without it.
#:
#: Maeda's advice for AFIR is to run the permissive end -- the gamma he
#: recommends is the ten-day 106.9 kJ/mol and not the one-hour 93.3 -- and it
#: is advice about *searching*: a search that misses a path has found nothing,
#: while one that finds too much is sorted out afterwards.  This is not a
#: search.  It decides what stays in the box, and there the two errors are not
#: the same size: a budget that allows what the temperature will not hands
#: back a structure with nothing anywhere saying it is impossible, and the
#: user goes on from it.  The strict end is right here for the same reason the
#: permissive end is right there.
_THERMAL_SECONDS = 3600.0


_BOLTZMANN_SI = 1.380649e-23          # J/K
_PLANCK_SI = 6.62607015e-34           # J s
_GAS_CONSTANT = 1.987204259e-3        # kcal/(mol K)
#: 13.8 thousand million years, in seconds.  The longest waiting time worth
#: printing as a number: past it every barrier reads the same, and a sentence
#: says what a mantissa and an exponent do not.
_UNIVERSE_SECONDS = 4.35e17


def thermal_ceiling(temperature, seconds):
    """The highest barrier a structure crosses at *temperature* within
    *seconds*, in kcal/mol.

    Eyring turned around.  A rate is k = (kB T / h) exp(-dG/RT), and a barrier
    is crossed about once in 1/k -- so asking how high it may be to be crossed
    within a given time gives

        dG = R T ln(kB T tau / h)

    which is the whole content of "possible at this temperature".  It is not
    one number: it depends on how long the user is prepared to wait, and the
    two make very different chemistry.  At 298 K a second buys 17.5 kcal/mol,
    an hour 22.3 and a year 27.7; at 773 K an hour buys 59.3.

    Measured against it on a benzene with GFN2, the ring bond stretched with
    everything else relaxed: 1.55 A costs 9.8 kcal/mol, 1.65 A costs 21.3 and
    1.75 A costs 34.7.  So at room temperature and an hour the bond can be
    pulled to about 1.66 A and no further -- the ring cannot be torn open, and
    that falls out of the energy rather than out of a rule about aromatics.
    At 773 K the same ceiling reaches 1.9 A, which is the temperature at which
    benzene really does come apart.
    """
    T = max(1.0, float(temperature))
    # A picosecond is the floor, and it is not arithmetic hygiene: below about
    # a tenth of that a molecule has not finished one vibration, so there is no
    # crossing to speak of and no barrier height that means anything.  It also
    # keeps the logarithm away from zero.
    tau = max(1e-12, float(seconds))
    inside = _BOLTZMANN_SI * T * tau / _PLANCK_SI
    if inside <= 1.0:
        return 0.0
    return _GAS_CONSTANT * T * math.log(inside)


def thermal_temperature(kcal, seconds=_THERMAL_SECONDS):
    """The temperature at which a barrier of *kcal* is crossed within *seconds*.

    :func:`thermal_ceiling` turned around again.  A drag has to be told what
    it may do at the temperature it was given, so there the question runs one
    way; a scan does not -- it walks the whole path and can then answer the
    question a chemist actually asks, which is *how hot*.  "+29 kcal/mol"
    means nothing until it is "and that wants 440 K to happen within the
    hour", which is the difference between a number and an experiment.

    The ceiling rises with temperature everywhere above a kelvin (T times the
    logarithm of T, both increasing), so bisection finds the answer and there
    is no need for the Lambert function this inverts to.  Returns ``None``
    when no temperature under 5000 K will do it -- past that a molecule is not
    a molecule, and saying so is better than printing a number.
    """
    try:
        wanted = float(kcal)
    except (TypeError, ValueError):
        return None
    if wanted <= thermal_ceiling(1.0, seconds):
        return 1.0
    low, high = 1.0, 5000.0
    if thermal_ceiling(high, seconds) < wanted:
        return None
    for _ in range(80):
        middle = 0.5 * (low + high)
        if thermal_ceiling(middle, seconds) < wanted:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)
