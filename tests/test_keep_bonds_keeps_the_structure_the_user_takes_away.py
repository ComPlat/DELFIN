"""The switch guards what is on screen, and it did not guard what is kept.

``Keep bonds`` reads the bonding off what a follow step hands back and
replaces a step that made or broke one with the last step that did not.  That
half was right from the start: driven with it on, the geometry xtb relaxes
holds every bond of an ethane through a four-angstrom pull on the far methyl,
where the same drag without it breaks the C-C at the sixth message.

What it did not hold is the *coordinate box* -- which is a different structure
during a drag.  The page writes its own model there ten times a second, and
that model is where the cursor is rather than where the chemistry allowed; the
follow was to write over it.  Two ways round that:

* An **allowed** step said nothing, and the box was written only when there
  was something to say.  So a drag the wall was letting through left the
  cursor's geometry standing: measured on the ethane, the box held C-C at
  2.19, 2.52, 2.86 and 3.19 A on four consecutive allowed steps -- every one
  a broken bond -- while the wall's own kept structure was at 1.55.
* The **release** is a raw write unless a wall claims the box, and the claim
  asked only whether the budget was on.  With Keep bonds alone the last word
  on the drag was therefore the mouse's: the box held 1.521 A and seven bonds
  all through, then 5.521 A and six the moment the hand let go.

The box is not a detail.  It is what Copy and Submit read and what the next
calculation starts from, so a switch that leaves a torn molecule there has
kept nothing that matters.
"""
from __future__ import annotations

import pytest

from delfin.dashboard import climb as _climb
from delfin.dashboard import gfn_optimize as gfn

# ``have_fast_gradients`` is climb's, not gfn_optimize's.  Written against the
# wrong module this never fired here -- ``find_xtb()`` answers on this box, so
# the ``and`` short-circuited before reaching it -- and broke collection on a
# runner with no xtb, where the second half is exactly what does get read.
_needs_xtb = pytest.mark.skipif(
    gfn.find_xtb() is None and not _climb.have_fast_gradients(),
    reason='no xtb to relax with')

_ETHANE = """8
ethane
C  0.000000  0.000000  0.762900
C  0.000000  0.000000 -0.762900
H -0.505000  0.874000  1.162900
H -0.505000 -0.874000  1.162900
H  1.010000  0.000000  1.162900
H  0.505000  0.874000 -1.162900
H  0.505000 -0.874000 -1.162900
H -1.010000  0.000000 -1.162900
"""


def _drive(part, begin, helper, far=4.0, steps=2):
    """The messages a hand sends, and the box after each of them."""
    methyl = {1, 5, 6, 7}
    seen = []
    part._begin_gfn_follow()
    for n in range(1, steps + 1):
        part.submit_manip_sync.value = helper._drag_message(
            helper._shifted(begin, methyl, far * n / steps),
            'DELFIN drag-follow held=1,5,6,7')
        helper._quiet(part.state)
        seen.append(part.coords_widget.value)
    part.submit_manip_sync.value = helper._drag_message(
        helper._shifted(begin, methyl, far), 'DELFIN drag-end')
    helper._quiet(part.state)
    seen.append(part.coords_widget.value)
    return seen


@_needs_xtb
def test_the_box_still_has_the_bonds_after_the_hand_lets_go():
    """The whole of the switch, told by the structure that outlives the drag.

    Every message including the release, not only the ones the wall refused:
    the failure was that an allowed step and the release both went round the
    wall, so a test that looked at refusals alone would have passed on the
    broken behaviour.
    """
    helper = pytest.importorskip('test_the_budget_prices_a_relaxed_path')

    start = gfn.optimize_with_gfn(_ETHANE, 'gfn2', max_steps=400, timeout=300)
    assert start.get('ok'), start.get('status')
    begin = start['xyz']
    whole = len(gfn.bond_graph(begin))

    part = helper._a_part(begin)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.value = 'pull'
    part.submit_topology_btn.value = True
    assert part.state.get('topology_graph') is not None
    seen = _drive(part, begin, helper)

    for step, xyz in enumerate(seen):
        assert len(gfn.bond_graph(xyz)) == whole, (
            f'message {step}: the box has {len(gfn.bond_graph(xyz))} bonds '
            f'where the structure had {whole}')
    # And the release in particular, which is where the mouse used to win.
    assert helper._apart(seen[-1], 0, 1) < 2.0, helper._apart(seen[-1], 0, 1)


@_needs_xtb
def test_the_same_drag_without_the_switch_does_break_it():
    """Otherwise the test above is about a drag that was never hard enough."""
    helper = pytest.importorskip('test_the_budget_prices_a_relaxed_path')

    start = gfn.optimize_with_gfn(_ETHANE, 'gfn2', max_steps=400, timeout=300)
    begin = start['xyz']
    whole = len(gfn.bond_graph(begin))

    part = helper._a_part(begin)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.value = 'pull'
    part.submit_topology_btn.value = False
    seen = _drive(part, begin, helper)

    assert len(gfn.bond_graph(seen[-1])) < whole, (
        'the drag never broke anything, so keeping it whole proves nothing')


@_needs_xtb
def test_the_wall_lets_go_when_the_hand_does_and_says_that_it_held():
    """Three refusals ended the drag for the rest of the grab, silently.

    The stand-still that follows three refusals in a row reads a counter that
    is only put back to zero by an *allowed* step -- and no step can be allowed
    while the stand-still is skipping them all.  So the count could never come
    down inside one gesture: the drag stopped following and never followed
    again until the hand was let go of and put back.  Measured on this ethane
    with the rigid hand and Keep bonds on, answers 8 to 15 carried a wish with
    a C-C of 1.596 A and all seven bonds -- a geometry both walls allow -- and
    not one of them was computed.

    The budget's own stand-still has a way out (_still_spent lets the drag run
    again the moment the hand eases in) and this one now has the same one, in
    the hand's own currency: how far the wish is from the structure the wall
    kept.  Nearer than when it gave up means the user is easing off.

    And it says so.  Both of the wall's sentences were appended to a status
    row that had already been stored and shown, so neither was ever read
    again: fifteen readings across a drag with three refusals in it contained
    "bonding" 0 times, "Keep bonds" 0, "Three steps" 0.  Every one of them
    said "GFN2-xTB follows the drag" -- including the ones taken while whole
    messages were being skipped, which is a frozen drag reporting that it is
    following.
    """
    helper = pytest.importorskip('test_the_budget_prices_a_relaxed_path')

    start = gfn.optimize_with_gfn(_ETHANE, 'gfn2', max_steps=400, timeout=300)
    assert start.get('ok'), start.get('status')
    begin = start['xyz']
    methyl = {1, 5, 6, 7}

    part = helper._a_part(begin)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_relax_btn.value = True
    part.submit_hand_dd.value = 'move'          # the rigid hand, which tears
    part.submit_topology_btn.value = True
    part._begin_gfn_follow()

    ran, said = [], []
    # Out until the wall gives up, then back in again.
    for far in (0.4, 0.9, 1.5, 2.1, 2.7, 3.3, 3.9, 0.4, 0.3, 0.2, 0.1):
        was = int(part.state.get('gfn_follow_steps') or 0)
        part.submit_manip_sync.value = helper._drag_message(
            helper._shifted(begin, methyl, far),
            'DELFIN drag-follow held=1,5,6,7')
        helper._quiet(part.state, seconds=300)
        ran.append(int(part.state.get('gfn_follow_steps') or 0) > was)
        said.append(str(part.state.get('gfn_last_status') or ''))

    assert not all(ran), 'the wall never gave up, so nothing is being tested'
    assert ran[-1], (
        f'the hand came back in and the drag never restarted: ran={ran}')
    assert ran[-4:] == [True, True, True, True], ran
    assert int(part.state.get('topology_refused') or 0) == 0, (
        'the count survived the hand coming back in')

    spoke = sum(1 for one in said if 'bonding' in one)
    assert spoke, f'the wall held a step and never said so: {said[-1]!r}'
