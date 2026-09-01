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

from delfin.dashboard import gfn_optimize as gfn

_needs_xtb = pytest.mark.skipif(
    gfn.find_xtb() is None and not gfn.have_fast_gradients(),
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
