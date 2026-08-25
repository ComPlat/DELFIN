"""One route through the network, drawn -- and the ways that picture can lie.

An energy profile is the figure a mechanism is argued from, which makes it the
place where a quiet mistake costs the most.  Three of them are possible here
and all three are refused rather than papered over:

* **Mixing levels.**  One point taken from a cheaper method and the rest from
  the expensive one is a picture of no reaction at all, and it looks exactly
  like a good one.  Every number comes from the level that was chosen, and a
  point that has none is not drawn.
* **Bridging a gap.**  Two species joined across a saddle nobody computed reads
  as a step that has been established.  The chain breaks at the gap instead.
* **Quoting the tallest barrier as the rate.**  A barrier is measured from the
  state before it; a route is measured from where it started.  A mechanism
  whose first intermediate sits well uphill is governed by a later saddle
  *above the educt*, and the tallest single step is the wrong number.

Nothing here runs a calculation.  The energies are put into the document by
hand, because what is being tested is arithmetic and a drawing.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import reaction_graph as rg
from delfin.dashboard import reaction_profile as rp
from delfin.dashboard import tab_reaction_graph as tab
from delfin.dashboard.thermal import thermal_ceiling

pytest.importorskip('matplotlib')
widgets = pytest.importorskip('ipywidgets')

_K = rg.HARTREE_TO_KCAL

_A = "3\na\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n"
_B = "2\nb\nC 0.000 0.000 0.000\nO 1.128 0.000 0.000\n"
_C = "4\nc\nN 0.000 0.000 0.000\nH 1.012 0.000 0.000\nH -0.337 0.954 0.000\nH -0.337 -0.477 0.826\n"
_D = "3\nd\nH 0.000 0.000 0.000\nC 1.070 0.000 0.000\nN 2.226 0.000 0.000\n"


def _at(kcal, **over):
    return rg.Record(level=over.pop('level', 'DFT'), method='orca',
                     free_energy=kcal / _K, source={'kind': 'run'}, **over)


def _mechanism(tmp_path):
    """Educt, intermediate uphill, product downhill, two saddles.

    Chosen so the tallest single barrier and the highest point of the route
    are different steps: TS1 is +14 from the educt, TS2 is +18 from the
    intermediate, and TS2 is +26 above where the route started.
    """
    graph = rg.create(tmp_path / 'mechanism', name='mechanism')
    a = rg.add_state(graph, _A, _at(0.0, imaginary=0), label='educt')
    b = rg.add_state(graph, _B, _at(8.0, imaginary=0), label='intermediate')
    c = rg.add_state(graph, _C, _at(-6.0, imaginary=0), label='product')
    one = rg.add_transition(graph, _D, _at(14.0, imaginary=1),
                            source=a.id, target=b.id, label='TS1')
    two = rg.add_transition(graph, _D, _at(26.0, imaginary=1),
                            source=b.id, target=c.id, label='TS2')
    rg.save(graph)
    return graph, a, b, c, one, two


# ---------------------------------------------------------------------------
# Which way through
# ---------------------------------------------------------------------------

def test_a_branch_gives_more_than_one_route_and_a_loop_gives_none_extra(
        tmp_path):
    """Simple routes only. A network with a reversible step has infinitely
    many walks between two points and exactly one of them is a mechanism."""
    graph, a, b, c, one, two = _mechanism(tmp_path)
    side = rg.add_state(graph, _D, _at(2.0, imaginary=0), label='by-product')
    rg.add_transition(graph, _D, _at(20.0, imaginary=1),
                      source=a.id, target=side.id, label='TS-side')
    rg.add_transition(graph, _D, _at(21.0, imaginary=1),
                      source=side.id, target=c.id, label='TS-side-2')
    # And the step back, which must not turn into an infinity of routes.
    rg.add_transition(graph, _D, _at(18.0, imaginary=1),
                      source=b.id, target=a.id, label='TS1 back')

    found = rg.routes_between(graph, a.id, c.id)
    assert [len(one) for one in found] == [3, 3]
    assert [a.id, b.id, c.id] in found
    assert [a.id, side.id, c.id] in found


def test_there_is_no_route_to_somewhere_nothing_reaches(tmp_path):
    graph, a, b, c, one, two = _mechanism(tmp_path)
    alone = rg.add_state(graph, _D, _at(0.0), label='unconnected')
    assert rg.routes_between(graph, a.id, alone.id) == []
    assert rg.routes_between(graph, a.id, a.id) == []


# ---------------------------------------------------------------------------
# What the profile is measured against
# ---------------------------------------------------------------------------

def test_a_profile_reads_state_saddle_state_from_the_first_state(tmp_path):
    graph, a, b, c, one, two = _mechanism(tmp_path)
    found = rg.profile(graph, [a.id, b.id, c.id], 'DFT')
    assert [p['kind'] for p in found['points']] == [
        'state', 'saddle', 'state', 'saddle', 'state']
    assert [round(p['kcal'], 1) for p in found['points']] == [
        0.0, 14.0, 8.0, 26.0, -6.0]
    assert found['missing'] == []
    assert found['zero'] == a.id


def test_the_highest_point_of_the_route_is_not_the_tallest_barrier(tmp_path):
    """The commonest way of getting a mechanism's rate wrong on paper."""
    graph, a, b, c, one, two = _mechanism(tmp_path)
    assert round(rg.barrier(graph, one.id, 'DFT'), 1) == 14.0
    assert round(rg.barrier(graph, two.id, 'DFT'), 1) == 18.0, 'from b'

    found = rg.profile(graph, [a.id, b.id, c.id], 'DFT')
    top = rp.summit(found['points'])
    assert top['label'] == 'TS2'
    assert round(top['kcal'], 1) == 26.0, 'above the educt, which is the rate'


def test_a_point_with_nothing_at_this_level_is_named_and_not_drawn(tmp_path):
    """One point from a cheaper method and the rest from the expensive one is
    a picture of no reaction at all, and it looks exactly like a good one."""
    graph, a, b, c, one, two = _mechanism(tmp_path)
    cheap = rg.create(tmp_path / 'partly', name='partly')
    x = rg.add_state(cheap, _A, _at(0.0, imaginary=0), label='educt')
    y = rg.add_state(cheap, _B, _at(5.0, imaginary=0), label='product')
    edge = rg.add_transition(cheap, _D,
                             rg.Record(level='GFN2-xTB', free_energy=1.0,
                                       imaginary=1),
                             source=x.id, target=y.id, label='TS')
    found = rg.profile(cheap, [x.id, y.id], 'DFT')
    assert found['missing'] == [edge.id]
    assert [p['kcal'] is None for p in found['points']] == [False, True, False]


def test_a_route_whose_zero_is_unpriced_prices_nothing(tmp_path):
    """Everything on a profile is measured against the first state. Without
    it there is no scale, and inventing one from the second state would draw
    a different reaction."""
    graph = rg.create(tmp_path / 'noscale', name='noscale')
    a = rg.add_state(graph, _A, rg.Record(level='DFT', method='orca'))
    b = rg.add_state(graph, _B, _at(5.0))
    rg.add_transition(graph, _D, _at(9.0, imaginary=1), source=a.id,
                      target=b.id)
    found = rg.profile(graph, [a.id, b.id], 'DFT')
    assert all(p['kcal'] is None for p in found['points'])
    assert a.id in found['missing']


# ---------------------------------------------------------------------------
# The picture
# ---------------------------------------------------------------------------

def test_every_species_is_a_flat_bar_and_the_joins_are_dashed(tmp_path):
    """A mechanism has no reaction coordinate. Between an intermediate and the
    saddle after it nobody calculated anything, and a smooth line there would
    be a path that was never walked."""
    graph, a, b, c, one, two = _mechanism(tmp_path)
    found = rg.profile(graph, [a.id, b.id, c.id], 'DFT')
    fig = rp.figure(found['points'], level='DFT')
    ax = fig.axes[0]
    flat = [ln for ln in ax.lines
            if len(set(ln.get_ydata())) == 1 and ln.get_linestyle() == '-']
    dashed = [ln for ln in ax.lines if ln.get_linestyle() == '--']
    assert len(flat) == 5, 'one bar per species and saddle'
    assert len(dashed) == 4, 'and a join between each pair'


def test_a_gap_breaks_the_chain_rather_than_being_bridged(tmp_path):
    """Two species joined across a saddle nobody computed reads as a step
    that has been established."""
    graph = rg.create(tmp_path / 'gap', name='gap')
    x = rg.add_state(graph, _A, _at(0.0), label='educt')
    y = rg.add_state(graph, _B, _at(5.0), label='product')
    rg.add_transition(graph, _D, rg.Record(level='DFT', imaginary=1),
                      source=x.id, target=y.id, label='TS')
    found = rg.profile(graph, [x.id, y.id], 'DFT')
    fig = rp.figure(found['points'], level='DFT')
    dashed = [ln for ln in fig.axes[0].lines if ln.get_linestyle() == '--']
    assert dashed == [], 'nothing is joined across the gap'


def test_the_thermal_ceiling_is_a_line_across_the_whole_mechanism(tmp_path):
    graph, a, b, c, one, two = _mechanism(tmp_path)
    found = rg.profile(graph, [a.id, b.id, c.id], 'DFT')
    ceiling = thermal_ceiling(298.15, 3600.0)
    fig = rp.figure(found['points'], level='DFT', ceiling=ceiling,
                    ceiling_label='298 K, an hour')
    dotted = [ln for ln in fig.axes[0].lines if ln.get_linestyle() == ':']
    assert len(dotted) == 1
    assert dotted[0].get_ydata()[0] == pytest.approx(ceiling)


def test_the_line_under_the_picture_says_the_rate_and_the_verdict(tmp_path):
    graph, a, b, c, one, two = _mechanism(tmp_path)
    found = rg.profile(graph, [a.id, b.id, c.id], 'DFT')
    ceiling = thermal_ceiling(298.15, 3600.0)
    said = rp.said_about(found['points'], missing=found['missing'],
                         ceiling=ceiling, window='an hour')
    assert 'TS2' in said and '+26.0' in said
    assert 'closed' in said

    warm = rp.said_about(found['points'], missing=found['missing'],
                         ceiling=thermal_ceiling(500.0, 3600.0),
                         window='an hour')
    assert 'open' in warm


def test_the_picture_is_a_png_carried_in_as_a_data_uri(tmp_path):
    graph, a, b, c, one, two = _mechanism(tmp_path)
    found = rg.profile(graph, [a.id, b.id, c.id], 'DFT')
    said = rp.profile_html(found['points'], level='DFT', note='a note')
    assert "src='data:image/png;base64," in said
    assert 'a note' in said
    assert rp.profile_png(found['points'], level='DFT')[:8] == b'\x89PNG\r\n\x1a\n'


# ---------------------------------------------------------------------------
# Driving it, and taking the picture away
# ---------------------------------------------------------------------------

def _panel(tmp_path):
    return tab.ReactionGraphPanel(tmp_path)


def test_the_two_ends_offer_the_states_and_the_press_appears_with_a_route(
        tmp_path):
    graph, a, b, c, one, two = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    assert [value for _label, value in panel.from_dd.options] == \
        [a.id, b.id, c.id]
    panel.from_dd.value = a.id
    panel.to_dd.value = c.id
    assert panel.draw_btn.layout.display == ''
    assert panel.route_dd.layout.display == 'none', 'one way, so no box'


def test_a_branch_puts_the_third_box_on_screen(tmp_path):
    graph, a, b, c, one, two = _mechanism(tmp_path)
    side = rg.add_state(graph, _D, _at(2.0, imaginary=0), label='by-product')
    rg.add_transition(graph, _D, _at(20.0, imaginary=1), source=a.id,
                      target=side.id)
    rg.add_transition(graph, _D, _at(21.0, imaginary=1), source=side.id,
                      target=c.id)
    rg.save(graph)
    panel = _panel(tmp_path)
    panel.from_dd.value = a.id
    panel.to_dd.value = c.id
    assert panel.route_dd.layout.display == ''
    assert len(panel.route_dd.options) == 2


def test_no_route_between_two_ends_offers_no_press(tmp_path):
    graph, a, b, c, one, two = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    panel.from_dd.value = c.id
    panel.to_dd.value = a.id
    assert panel.draw_btn.layout.display == 'none'


def test_drawing_keeps_the_picture_in_the_graph_folder(tmp_path):
    """A figure a person can only look at is one they will redraw by hand for
    the paper. Beside the calculations it was made from, under a name saying
    which route and which level, it is one they can take."""
    graph, a, b, c, one, two = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    panel.from_dd.value = a.id
    panel.to_dd.value = c.id
    panel._on_draw()

    assert "data:image/png;base64," in panel.picture.value
    kept = sorted((graph.folder / tab.FIGURES).glob('*.png'))
    assert len(kept) == 1
    assert kept[0].read_bytes()[:4] == b'\x89PNG'
    assert a.id in kept[0].name and 'DFT' in kept[0].name

    line = [x for x in rg.history(rg.load(graph.folder))
            if x['what'] == 'profile drawn']
    assert line and line[0]['level'] == 'DFT'
    assert line[0]['route'] == [a.id, b.id, c.id]
    assert line[0]['figure'].startswith(tab.FIGURES + '/')


def test_a_route_with_nothing_priced_at_the_level_draws_nothing_and_says_so(
        tmp_path):
    graph = rg.create(tmp_path / 'bare', name='bare')
    x = rg.add_state(graph, _A, rg.Record(level='GFN2-xTB', free_energy=-5.0))
    y = rg.add_state(graph, _B, rg.Record(level='GFN2-xTB', free_energy=-5.1))
    rg.add_transition(graph, _D, rg.Record(level='GFN2-xTB', imaginary=1),
                      source=x.id, target=y.id)
    rg.save(graph)
    panel = _panel(tmp_path)
    panel.from_dd.value = x.id
    panel.to_dd.value = y.id
    panel.level_dd.value = 'GFN2-xTB'
    panel._on_draw()
    assert 'data:image' in panel.picture.value, 'the two ends are priced'

    # And a level nothing has at all.
    rg.add_record(graph, x.id, _A, rg.Record(level='DFT', method='orca'))
    rg.save(graph)
    panel = _panel(tmp_path)
    panel.from_dd.value = x.id
    panel.to_dd.value = y.id
    panel.level_dd.value = 'DFT'
    panel._on_draw()
    assert panel.picture.value == ''
    assert 'Nothing on that route is priced at DFT' in panel.status.value
