"""Why a drag was refused, and where the budget is counting from.

A user reported that several thermal deformations in a row "sometimes" do not
work, and that switching Dynamik Opt off and on again makes them work.  Both
halves of that are the same thing seen twice: the budget is counted from an
anchor, the anchor stays where it was set, and toggling Dynamik Opt sets a new
one.  That mechanism is right -- a ceiling is what this state can do within the
hour, and two changes in sequence from one starting state is a harder question
than one -- but it was entirely invisible, so he found it by accident.

What was actually broken, and is fixed here, is narrower: the running maximum a
drag is refused on was cleared when a drag *ended*, and a drag can end without
the message that clears it.  Then the next grab is refused for the last one's
worst moment and nothing says so.  That is the "sometimes".

And what the refusal says.  Measured on his complex, GFN2 at charge +1, a
phenolate oxygen pulled off the manganese: the price climbs smoothly to
-1.6 kcal/mol and the next answer is +68.0, at a Mn-O of 3.19 A, with a
bromide arriving at the metal.  It is +68 whether the cursor advances a quarter
of an Angstrom per answer or a fiftieth -- 5, 11 and 97 answers to the same
place -- so there is no geometry in between and coming at it more gently finds
nothing.  Under GFN-FF the same drag reads +2.5 and then +108.7.  A user who is
not told that goes on trying, which is what ten minutes of dragging was.
"""

from __future__ import annotations

from editor_source import EDITOR_SOURCE


def _between(text, start, end):
    """The source between two markers, so a claim is about one function."""
    head = text.index(start)
    return text[head:text.index(end, head)]


def test_a_grab_clears_what_the_last_drag_went_through():
    """However the last drag ended, the next one starts from a clean slate.

    The running maximum is cleared on release, and the release is one message
    from the page -- sent from an animation frame when the pointer state
    changes.  A drag that ends any other way never sends it: the tab goes to
    the background and the frame loop stops, the pointer is cancelled, the page
    is reloaded under a hand that is still down.  What is left standing then is
    the maximum, and a maximum is what the wall refuses on -- so the next grab
    is refused for a crossing it had nothing to do with.

    Clearing it where a drag *begins* as well makes it not depend on a message
    arriving.  A grab is the one moment there certainly is one.
    """
    begin = _between(EDITOR_SOURCE, "def _begin_gfn_follow", "def _clear_thermal_wall")
    assert "_clear_thermal_wall()" in begin
    # And after the early return for a drag that is already running, or a
    # second push of the same drag would wipe its own history.
    assert begin.index("it keeps the run it began") < begin.index(
        "_clear_thermal_wall()")
    # The floor goes back afterwards.  Clearing leaves the maximum absent, and
    # absent lets the first answer set it to whatever that answer costs -- so a
    # drag that starts below the anchor would carry a negative maximum and be
    # allowed to climb that much further than the temperature pays for.  The
    # anchor is the zero; falling below it earns no credit.
    assert begin.index("_clear_thermal_wall()") < begin.index(
        "state['thermal_peak'] = 0.0")


def test_what_a_drag_went_through_is_forgotten_with_the_rest():
    """The clearing is one list, and the price of the last answer is on it.

    "One answer cost this much" is measured against the answer before it, so
    the number belongs to the drag that is over exactly as the maximum and the
    slope do.  Kept, the first answer of the next drag would be compared with
    the last answer of the last one -- a different structure, and after a
    rollback a different branch of the surface.
    """
    clearing = _between(EDITOR_SOURCE, "def _clear_thermal_wall",
                        "def _end_gfn_follow")
    for gone in ("thermal_peak", "thermal_good_peak", "thermal_slope",
                 "thermal_priced"):
        assert f"state.pop('{gone}', None)" in clearing


def test_a_step_is_not_a_slope_and_is_not_reported_as_one():
    """The two refusals have different next moves, so they are said apart.

    A refusal on a slope can be worked with: ease off, come at it more slowly,
    and the structure stops near the ceiling.  A refusal on a step cannot --
    the same place is reached at any step size -- and the user who does not
    know that goes on pulling.  Judged against the ceiling rather than against
    a number of kcal/mol, so it means the same at 100 K and at 1500 K.
    """
    assert "a step, not a slope: one " in EDITOR_SOURCE
    assert "jump > max(1.0, ceiling)" in EDITOR_SOURCE
    # And what changed, named, because the bonding on the far side is the fact
    # that makes the number make sense.
    assert "on the far side the bonding" in EDITOR_SOURCE


def test_the_budget_says_where_it_counts_from():
    """A ceiling nobody can see the origin of is a rule rather than a measure.

    The refusal that has an obvious next move is the one where the budget was
    spent getting here: the structure is standing somewhere cheap, and whether
    that somewhere is a state the system would have thermalised in is a
    chemist's judgement.  Set here is how that judgement is entered, and the
    line now names it.  Without that the discovery is that toggling Dynamik Opt
    makes the drag work again, which is the same thing happening by accident.
    """
    note = _between(EDITOR_SOURCE, "def _thermal_note", "def _thermal_wait")
    assert "Set here measures from this structure instead" in note
    assert "Set here" in note
    # And nothing anywhere offers a button by a name the toolbar does not have.
    assert "Press Measure from here" not in EDITOR_SOURCE
