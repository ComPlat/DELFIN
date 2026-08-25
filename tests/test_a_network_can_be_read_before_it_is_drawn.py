"""The Reactions tab, driven without a browser.

A drawing is what the graph is for, and it is also where mistakes hide: an edge
drawn to the wrong end and an edge drawn to the right one are both a line.  So
the first thing over the document is a list, and these tests read it -- which
is the only way it gets tested at all, because there is no browser on the
machine this is built on.

The claims worth pinning are about what the tab is willing to say.  A level box
that offers a method nothing has been computed at invites a diagram of nothing.
An uncomputed barrier shown like a closed one reports a conclusion about
calculations nobody has run.  A number printed without saying whether it is G
or E leaves a reader unable to tell a free-energy profile from an electronic
one afterwards.  Each of those is a picture that looks right.
"""

from __future__ import annotations

import pytest

from delfin.dashboard import reaction_graph as rg
from delfin.dashboard import tab_reaction_graph as tab

widgets = pytest.importorskip('ipywidgets')


_WATER = "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n"
_CO = "2\ncarbon monoxide\nC 0.000 0.000 0.000\nO 1.128 0.000 0.000\n"
_HCN = "3\nhydrogen cyanide\nH 0.000 0.000 0.000\nC 1.070 0.000 0.000\nN 2.226 0.000 0.000\n"
_NH3 = ("4\nammonia\nN 0.000 0.000 0.000\nH 1.012 0.000 0.000\nH -0.337 0.954 0.000\nH -0.337 -0.477 0.826\n")

#: Two barriers either side of what an hour at room temperature pays for, which
#: is 22.3 kcal/mol -- see :func:`thermal.thermal_ceiling`.
_LOW = 10.0 / rg.HARTREE_TO_KCAL
_HIGH = 40.0 / rg.HARTREE_TO_KCAL


def _mechanism(root, name='mechanism'):
    """Three states and two steps: one the hour pays for, one it does not."""
    graph = rg.create(root / name, name=name)
    a = rg.add_state(graph, _WATER,
                     rg.Record(level='GFN2-xTB', method='gfn2',
                               free_energy=-5.000, imaginary=0),
                     label='educt')
    b = rg.add_state(graph, _CO,
                     rg.Record(level='GFN2-xTB', method='gfn2',
                               free_energy=-5.030, imaginary=0),
                     label='intermediate')
    c = rg.add_state(graph, _HCN,
                     rg.Record(level='GFN2-xTB', method='gfn2',
                               free_energy=-5.060, imaginary=0),
                     label='side product')
    easy = rg.add_transition(
        graph, _HCN,
        rg.Record(level='GFN2-xTB', method='gfn2', free_energy=-5.000 + _LOW,
                  imaginary=1, frequency=-412.0),
        source=a.id, target=b.id, label='step one',
        confirmed={'how': 'mode-down', 'ends': [a.id, b.id]})
    hard = rg.add_transition(
        graph, _HCN,
        rg.Record(level='GFN2-xTB', method='gfn2', free_energy=-5.000 + _HIGH,
                  imaginary=1),
        source=a.id, target=c.id, label='side reaction')
    rg.save(graph)
    return graph, a, b, c, easy, hard


def _panel(root):
    return tab.ReactionGraphPanel(root)


# ---------------------------------------------------------------------------
# Finding graphs
# ---------------------------------------------------------------------------

def test_a_graph_is_a_folder_with_a_document_in_it(tmp_path):
    _mechanism(tmp_path, 'one')
    (tmp_path / 'an ordinary job').mkdir()
    (tmp_path / 'an ordinary job' / 'job.inp').write_text('!', encoding='utf-8')
    found = tab.graphs_in(tmp_path)
    assert [p.name for p in found] == ['one']


def test_the_calculations_inside_a_graph_are_not_offered_as_graphs(tmp_path):
    """A graph carries a runs/ full of ordinary job folders. A search that
    went down into them would offer a calculation as though it were a
    network."""
    graph, *_ = _mechanism(tmp_path)
    node = graph.nodes[0]
    _, real = rg.run_folder(graph, node.id, 'r2SCAN-3c')
    real.mkdir(parents=True)
    assert [p.name for p in tab.graphs_in(tmp_path)] == ['mechanism']


def test_a_directory_with_no_graphs_says_so_rather_than_looking_broken(
        tmp_path):
    panel = _panel(tmp_path)
    assert panel.graph is None
    assert 'No reaction graphs' in panel.status.value


# ---------------------------------------------------------------------------
# What the level box may offer
# ---------------------------------------------------------------------------

def test_the_level_box_offers_only_levels_the_network_actually_has(tmp_path):
    """The alternative -- every method DELFIN can run -- invites drawing a
    network at a level nothing in it has been computed at."""
    graph, a, *_ = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    assert list(panel.level_dd.options) == ['GFN2-xTB']

    rg.add_record(graph, a.id, _WATER,
                  rg.Record(level='r2SCAN-3c/CPCM(THF)', method='orca',
                            free_energy=-76.4, imaginary=0))
    rg.save(graph)
    panel.refresh_graphs()
    assert list(panel.level_dd.options) == ['GFN2-xTB', 'r2SCAN-3c/CPCM(THF)']


def test_choosing_a_level_nothing_is_complete_at_names_what_is_missing(
        tmp_path):
    """An energy diagram is a selection out of a document that may be complete
    at one level and not another. The honest sentence is said where the level
    is chosen, not in a footnote."""
    graph, a, *_ = _mechanism(tmp_path)
    rg.add_record(graph, a.id, _WATER,
                  rg.Record(level='r2SCAN-3c', method='orca',
                            free_energy=-76.4, imaginary=0))
    rg.save(graph)
    panel = _panel(tmp_path)
    panel.level_dd.value = 'r2SCAN-3c'
    said = panel.summary.value
    assert 'No energy at r2SCAN-3c for' in said
    assert a.id not in said.split('for:')[1], 'the one that has it is not listed'


# ---------------------------------------------------------------------------
# What the list says
# ---------------------------------------------------------------------------

def test_every_state_and_every_step_out_of_it_is_on_the_list(tmp_path):
    graph, a, b, c, easy, hard = _mechanism(tmp_path)
    refs = [ref for ref, _ in tab.network_lines(graph, 'GFN2-xTB')]
    assert refs == [a.id, easy.id, hard.id, b.id, c.id]


def test_a_side_reaction_is_a_second_step_out_of_the_same_state(tmp_path):
    """Which is the whole reason it is a graph: branches are ordinary."""
    graph, a, b, c, easy, hard = _mechanism(tmp_path)
    lines = dict(tab.network_lines(graph, 'GFN2-xTB'))
    assert 'intermediate' in lines[easy.id]
    assert 'side product' in lines[hard.id]


def test_the_temperature_decides_which_steps_are_open_and_it_moves(tmp_path):
    graph, a, b, c, easy, hard = _mechanism(tmp_path)
    cool = dict(tab.network_lines(graph, 'GFN2-xTB', temperature=298.15,
                                  seconds=3600.0, window='an hour'))
    assert 'open within an hour' in cool[easy.id]
    assert 'closed within an hour' in cool[hard.id]

    hot = dict(tab.network_lines(graph, 'GFN2-xTB', temperature=773.15,
                                 seconds=3600.0, window='an hour'))
    assert 'open within an hour' in hot[hard.id], 'the same barrier, hotter'


def test_a_longer_window_opens_what_the_hour_would_not(tmp_path):
    """The second half of "possible", and the reason it is a control here and
    a constant in the editor: a network is read deliberately, and "does this
    go overnight" is a question somebody brings to a mechanism."""
    graph, a, b, c, easy, hard = _mechanism(tmp_path)
    # At 450 K the ceiling is 34.0 kcal/mol over an hour and 42.2 over a year,
    # and the side reaction's barrier is 40.0 -- so the window alone decides it.
    seconds = dict(tab.WINDOWS)
    hour = dict(tab.network_lines(graph, 'GFN2-xTB', temperature=450.0,
                                  seconds=seconds['an hour'], window='an hour'))
    year = dict(tab.network_lines(graph, 'GFN2-xTB', temperature=450.0,
                                  seconds=seconds['a year'], window='a year'))
    assert 'closed within an hour' in hour[hard.id]
    assert 'open within a year' in year[hard.id]


def test_a_barrier_nobody_has_computed_says_that_and_not_closed(tmp_path):
    """The third answer, said rather than implied. Drawn like a closed one it
    reports a conclusion about calculations nobody has run, and a reader has
    no way to tell which grey is which."""
    graph, a, b, c, easy, hard = _mechanism(tmp_path)
    rg.add_record(graph, a.id, _WATER,
                  rg.Record(level='g-xTB', method='gxtb', free_energy=-5.0))
    rg.add_record(graph, easy.id, _HCN,
                  rg.Record(level='g-xTB', method='gxtb', free_energy=-4.9,
                            imaginary=1))
    rg.save(graph)
    lines = dict(tab.network_lines(graph, 'g-xTB'))
    assert 'open' in lines[easy.id] or 'closed' in lines[easy.id]
    assert 'no energy yet for' in lines[hard.id]
    assert 'closed' not in lines[hard.id]


def test_a_number_says_whether_it_is_a_free_energy_or_an_electronic_one(
        tmp_path):
    """Nobody can tell them apart from the picture afterwards."""
    free = rg.Record(level='x', free_energy=-5.5)
    electronic = rg.Record(level='x', energy=-5.5)
    assert tab.said_energy(free).startswith('G ')
    assert tab.said_energy(electronic).startswith('E ')
    assert tab.said_energy(None) == '--'


def test_a_saddle_and_a_minimum_are_told_apart_on_the_line(tmp_path):
    assert 'minimum' in tab.said_energy(rg.Record(level='x', energy=-1.0,
                                                  imaginary=0))
    said = tab.said_energy(rg.Record(level='x', energy=-1.0, imaginary=1,
                                     frequency=-412.0))
    assert 'saddle' in said and '412i' in said
    assert '3 imaginary' in tab.said_energy(
        rg.Record(level='x', energy=-1.0, imaginary=3))


# ---------------------------------------------------------------------------
# What the panel says about one thing
# ---------------------------------------------------------------------------

def test_one_of_something_is_not_written_as_a_plural(tmp_path):
    """A small thing, and the kind of small thing that tells a reader the
    numbers beside it were not looked at either."""
    graph = rg.create(tmp_path / "one", name="one")
    a = rg.add_state(graph, _WATER, rg.Record(level="x", free_energy=-5.0))
    b = rg.add_state(graph, _CO, rg.Record(level="x", free_energy=-5.0))
    rg.add_transition(graph, _HCN, rg.Record(level="x", free_energy=-4.9,
                                             imaginary=1),
                      source=a.id, target=b.id)
    rg.add_pending(graph, a.id, run="runs/one", level="r2SCAN-3c")
    said = tab.summary_html(graph, "x")
    assert "2 states" in said
    assert "1 transition." in said and "1 transitions" not in said
    assert "1 calculation running" in said


def test_the_panel_shows_every_level_and_not_only_the_chosen_one(tmp_path):
    """What a user is deciding is whether to start a calculation, and the gap
    is the thing they are deciding about."""
    graph, a, *_ = _mechanism(tmp_path)
    rg.add_record(graph, a.id, _WATER,
                  rg.Record(level='r2SCAN-3c', method='orca',
                            free_energy=-76.4, imaginary=0))
    said = tab.detail_html(graph, a.id, 'GFN2-xTB')
    assert 'GFN2-xTB' in said and 'r2SCAN-3c' in said


def test_a_missing_level_is_stated_on_the_thing_that_is_missing_it(tmp_path):
    graph, a, b, *_ = _mechanism(tmp_path)
    rg.add_record(graph, a.id, _WATER,
                  rg.Record(level='r2SCAN-3c', method='orca',
                            free_energy=-76.4))
    assert 'Nothing at' in tab.detail_html(graph, b.id, 'r2SCAN-3c')
    assert 'Nothing at' not in tab.detail_html(graph, a.id, 'r2SCAN-3c')


def test_an_unconfirmed_transition_says_what_has_not_been_checked(tmp_path):
    """A saddle optimiser returns a saddle whatever it was given."""
    graph, a, b, c, easy, hard = _mechanism(tmp_path)
    assert 'Confirmed to join these two' in tab.detail_html(graph, easy.id,
                                                            'GFN2-xTB')
    said = tab.detail_html(graph, hard.id, 'GFN2-xTB')
    assert 'Not confirmed' in said
    assert 'imaginary mode down' in said


def test_where_a_record_came_from_is_on_the_panel(tmp_path):
    """A record says the answer; this is the beginning of how it was got."""
    graph = rg.create(tmp_path / 'one', name='one')
    node = rg.add_state(
        graph, _WATER,
        rg.Record(level='GFN2-xTB', method='gfn2', free_energy=-5.0,
                  source={'kind': 'editor', 'gesture': 'to the saddle'}))
    assert 'from editor' in tab.detail_html(graph, node.id, 'GFN2-xTB')


def test_a_running_calculation_is_shown_on_the_thing_it_is_for(tmp_path):
    graph, a, *_ = _mechanism(tmp_path)
    relative, _ = rg.run_folder(graph, a.id, 'r2SCAN-3c')
    rg.add_pending(graph, a.id, run=relative, level='r2SCAN-3c')
    said = tab.detail_html(graph, a.id, 'GFN2-xTB')
    assert 'Running' in said and relative in said


def test_a_label_from_a_file_cannot_write_markup_into_the_panel(tmp_path):
    """Names come out of a document a person edits, and the panel is HTML."""
    graph = rg.create(tmp_path / 'one', name='one')
    node = rg.add_state(graph, _WATER, rg.Record(level='x', energy=-1.0),
                        label='<img src=x onerror=alert(1)>')
    said = tab.detail_html(graph, node.id, 'x')
    assert '<img' not in said
    assert '&lt;img' in said


# ---------------------------------------------------------------------------
# Driving it
# ---------------------------------------------------------------------------

def test_the_panel_opens_the_graph_and_selects_something(tmp_path):
    graph, a, *_ = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    assert panel.graph is not None
    assert panel.network.value == a.id
    assert 'educt' in panel.detail.value


def test_changing_the_temperature_rewrites_the_list_without_losing_the_choice(
        tmp_path):
    graph, a, b, c, easy, hard = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    panel.network.value = hard.id
    panel.temperature.value = 773.15
    assert panel.network.value == hard.id, 'the selection survived the rebuild'
    shown = dict((ref, line) for line, ref in panel.network.options)
    assert 'open within' in shown[hard.id]


def test_a_new_graph_is_made_where_the_calculations_are(tmp_path):
    panel = _panel(tmp_path)
    panel.new_name.value = 'Mn hydrogenation'
    panel._on_new()
    made = tab.graphs_in(tmp_path)
    assert [p.name for p in made] == ['Mn_hydrogenation']
    assert (made[0] / rg.STRUCTURES).is_dir()
    assert (made[0] / rg.RUNS).is_dir()
    assert panel.graph is not None and panel.graph.name == 'Mn hydrogenation'


def test_a_new_graph_needs_a_name_and_will_not_overwrite_one(tmp_path):
    panel = _panel(tmp_path)
    panel.new_name.value = '   '
    panel._on_new()
    assert 'name' in panel.status.value.lower()
    assert tab.graphs_in(tmp_path) == []

    panel.new_name.value = 'one'
    panel._on_new()
    panel.new_name.value = 'one'
    panel._on_new()
    assert 'already holds a graph' in panel.status.value
    assert len(tab.graphs_in(tmp_path)) == 1


def test_renaming_something_is_written_down_and_kept(tmp_path):
    graph, a, *_ = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    panel.network.value = a.id
    panel.label_box.value = 'the real educt'
    panel._on_rename()

    again = rg.load(graph.folder)
    assert again.node(a.id).label == 'the real educt'
    assert 'renamed' in [one['what'] for one in rg.history(again)]
    shown = dict((ref, line) for line, ref in panel.network.options)
    assert 'the real educt' in shown[a.id]


def test_a_note_survives_the_document_being_reopened(tmp_path):
    graph, a, b, c, easy, hard = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    panel.network.value = easy.id
    panel.note_box.value = 'checked against the paper, table 2'
    panel._on_note()
    again = rg.load(graph.folder)
    assert 'table 2' in again.edge(easy.id).note


def test_a_graph_that_will_not_open_says_so_instead_of_showing_nothing(
        tmp_path):
    graph, *_ = _mechanism(tmp_path)
    (graph.folder / rg.DOCUMENT).write_text('{ not json', encoding='utf-8')
    panel = _panel(tmp_path)
    assert panel.graph is None
    assert 'did not open' in panel.status.value


# ---------------------------------------------------------------------------
# The tab itself
# ---------------------------------------------------------------------------

def test_the_tab_registers_itself_and_builds(tmp_path):
    from delfin.dashboard import tab_registry

    assert 'delfin.dashboard.tab_reaction_graph' in \
        tab_registry._BUILTIN_DYNAMIC_TABS, 'or it is never imported'

    class _Ctx:
        calc_dir = tmp_path

    widget, refs = tab.create_tab(_Ctx())
    assert widget is not None
    assert 'reaction_graph' in refs


def test_no_button_is_offered_for_something_that_is_not_built_yet(tmp_path):
    """A row of controls that raise "not implemented" is a worse lie than an
    empty toolbar. Opening a structure in the editor and sending one to a
    builder are the next steps, and they are not here."""
    graph, *_ = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    named = [b.description for b in
             (panel.new_btn, panel.reload_btn, panel.rename_btn,
              panel.note_btn)]
    assert named == ['New', 'Reload', 'Rename', 'Save note']
    assert not hasattr(panel, 'open_in_editor_btn')
    assert not hasattr(panel, 'send_to_builder_btn')


# ---------------------------------------------------------------------------
# What the editor may hand over
# ---------------------------------------------------------------------------

def _offer(xyz, **over):
    got = dict(xyz=xyz, level='GFN2-xTB', method='gfn2', charge=0,
               multiplicity=1, energy=-5.0, free_energy=None, imaginary=0,
               frequency=None, ends=None, gesture='optimise')
    got.update(over)
    return got


def test_with_no_graph_open_there_is_nothing_to_offer(tmp_path):
    panel = _panel(tmp_path)
    assert panel.offer_label(_offer(_WATER)) is None


def test_a_structure_nobody_has_goes_in_as_a_new_state(tmp_path):
    _mechanism(tmp_path)
    panel = _panel(tmp_path)
    said = panel.offer_label(_offer(_NH3))
    assert said is not None and 'new state' in said
    grew = panel.take(_offer(_NH3, name='the new one'))
    assert 'is in' in grew
    again = rg.load(panel.graph.folder)
    assert [n.label for n in again.nodes][-1] == 'the new one'


def test_a_structure_the_graph_already_has_is_offered_as_another_geometry(
        tmp_path):
    """A GFN2 and a DFT structure are two geometries of one species, and this
    is where the second one arrives."""
    graph, a, *_ = _mechanism(tmp_path)
    panel = _panel(tmp_path)
    said = panel.offer_label(_offer(_WATER, level='r2SCAN-3c', method='orca'))
    assert said == 'Add as a r2SCAN-3c record on educt'
    panel.take(_offer(_WATER, level='r2SCAN-3c', method='orca', energy=-76.4))
    again = rg.load(graph.folder)
    assert rg.best(again.node(a.id), 'r2SCAN-3c').energy == -76.4
    assert len(again.nodes) == 3, 'and no fourth state was invented'


def test_a_walk_that_left_two_ends_lays_down_the_whole_step(tmp_path):
    """The shape a scan produces, and the most valuable hand-over there is: a
    barrier arrives with both the structures it is a barrier between."""
    graph = rg.create(tmp_path / 'fresh', name='fresh')
    panel = _panel(tmp_path)
    said = panel.offer_label(_offer(_HCN, ends=(_WATER, _CO), imaginary=1))
    assert said == 'Put both ends and the step in fresh'

    grew = panel.take(_offer(_HCN, ends=(_WATER, _CO), imaginary=1,
                             frequency=-412.0, gesture='to the saddle'))
    again = rg.load(graph.folder)
    assert len(again.nodes) == 2 and len(again.edges) == 1
    edge = again.edges[0]
    assert edge.source == again.nodes[0].id
    assert edge.target == again.nodes[1].id
    assert edge.records[0].frequency == -412.0
    assert 'joins' in grew


def test_the_same_walk_twice_does_not_put_its_educt_in_twice(tmp_path):
    """A duplicate state is a twin whose energies quietly disagree."""
    graph = rg.create(tmp_path / 'fresh', name='fresh')
    panel = _panel(tmp_path)
    panel.take(_offer(_HCN, ends=(_WATER, _CO), imaginary=1))
    panel.take(_offer(_HCN, ends=(_WATER, _HCN), imaginary=1))
    again = rg.load(graph.folder)
    assert len(again.nodes) == 3, [n.id for n in again.nodes]
    assert len(again.edges) == 2
    assert again.edges[0].source == again.edges[1].source, 'one educt'


def test_two_ends_that_are_the_same_structure_are_refused(tmp_path):
    rg.create(tmp_path / 'fresh', name='fresh')
    panel = _panel(tmp_path)
    said = panel.take(_offer(_HCN, ends=(_WATER, _WATER), imaginary=1))
    assert 'did not go into the graph' in said
    assert 'nothing between them' in said


def test_a_saddle_with_no_two_ends_is_not_offered_at_all(tmp_path):
    """A transition is an edge and an edge needs the two states it joins. Put
    in as a state it would be a minimum in the document that is not one.

    Standing beside the missing button is ``Down to both ends``, which appears
    exactly when there is a saddle and is what produces the ends.
    """
    _mechanism(tmp_path)
    panel = _panel(tmp_path)
    assert panel.offer_label(_offer(_NH3, imaginary=1)) is None
    assert panel.offer_label(_offer(_NH3, imaginary=0)) is not None


def test_what_went_in_says_what_produced_it(tmp_path):
    graph = rg.create(tmp_path / 'fresh', name='fresh')
    panel = _panel(tmp_path)
    panel.take(_offer(_WATER, gesture='climbed to a saddle'))
    again = rg.load(graph.folder)
    record = again.nodes[0].records[0]
    assert record.source['kind'] == 'editor'
    assert record.source['gesture'] == 'climbed to a saddle'
    assert 'record added' in [one['what'] for one in rg.history(again)]


def test_the_tab_shows_what_just_landed(tmp_path):
    graph = rg.create(tmp_path / 'fresh', name='fresh')
    panel = _panel(tmp_path)
    assert panel.network.options == ()
    panel.take(_offer(_WATER, name='educt'))
    shown = [ref for _line, ref in panel.network.options]
    assert shown == [rg.load(graph.folder).nodes[0].id]


def test_something_that_is_not_a_structure_is_not_offered(tmp_path):
    _mechanism(tmp_path)
    panel = _panel(tmp_path)
    assert panel.offer_label(_offer('')) is None
    assert panel.offer_label(_offer('nothing at all')) is None


# ---------------------------------------------------------------------------
# Back to the workbench
# ---------------------------------------------------------------------------

class _Tabs:
    def __init__(self):
        self.selected_index = 0


class _Dashboard:
    """Enough of a dashboard for the panel to find its way back to one."""

    def __init__(self, calc_dir):
        self.calc_dir = calc_dir
        self.submit_refs = {'coords_widget': widgets.Textarea(value='')}
        self.tabs_widget = _Tabs()
        self.tab_indices = {'Submit Job': 0, 'Reactions': 5}
        self.reaction_graph_refs = {}


def test_with_no_submit_tab_there_is_no_button_to_send_it_to_one(tmp_path):
    """This tab is usable on its own, and a press that would reach a tab that
    is not there is the row claiming something untrue."""
    _mechanism(tmp_path)
    panel = _panel(tmp_path)
    assert panel.to_editor_btn.layout.display == 'none'


def test_a_geometry_goes_into_the_submit_box_and_the_tab_is_shown(tmp_path):
    graph, a, *_ = _mechanism(tmp_path)
    ctx = _Dashboard(tmp_path)
    ctx.tabs_widget.selected_index = 5
    panel = tab.ReactionGraphPanel(tmp_path, ctx=ctx)
    panel.network.value = a.id
    assert panel.to_editor_btn.layout.display == ''

    panel._on_to_editor()
    written = ctx.submit_refs['coords_widget'].value
    assert written.splitlines()[0] == '3', 'the atom count is the real one'
    assert 'educt at GFN2-xTB' in written.splitlines()[1]
    assert 'O ' in written and 'H ' in written
    assert ctx.tabs_widget.selected_index == 0, 'and it went there'


def test_a_level_with_no_geometry_offers_no_way_to_open_it(tmp_path):
    graph, a, b, *_ = _mechanism(tmp_path)
    rg.add_record(graph, a.id, _WATER,
                  rg.Record(level='r2SCAN-3c', method='orca', energy=-76.4))
    rg.save(graph)
    ctx = _Dashboard(tmp_path)
    panel = tab.ReactionGraphPanel(tmp_path, ctx=ctx)
    panel.level_dd.value = 'r2SCAN-3c'
    panel.network.value = a.id
    assert panel.to_editor_btn.layout.display == ''
    panel.network.value = b.id
    assert panel.to_editor_btn.layout.display == 'none'


def test_sending_a_structure_to_the_editor_is_not_written_into_the_document(
        tmp_path):
    """A document that recorded "sent to the editor" would be recording an
    intention. What comes back comes back through the hand-over the other way,
    with whatever was done to it."""
    graph, a, *_ = _mechanism(tmp_path)
    before = len(rg.history(graph))
    ctx = _Dashboard(tmp_path)
    panel = tab.ReactionGraphPanel(tmp_path, ctx=ctx)
    panel.network.value = a.id
    panel._on_to_editor()
    assert len(rg.history(rg.load(graph.folder))) == before


def test_the_tab_hands_itself_to_the_dashboard_so_the_editor_can_find_it(
        tmp_path):
    """The registry keeps the widget and drops the refs, so a registered tab
    that wants to be reachable has to publish itself."""
    ctx = _Dashboard(tmp_path)
    tab.create_tab(ctx)
    assert 'reaction_graph' in ctx.reaction_graph_refs
    assert isinstance(ctx.reaction_graph_refs['reaction_graph'],
                      tab.ReactionGraphPanel)


# ---------------------------------------------------------------------------
# The two kinds of structure change
# ---------------------------------------------------------------------------

def test_a_geometry_says_which_of_the_two_kinds_made_it(tmp_path):
    """A structure dragged into shape and one a two-day optimisation produced
    are not interchangeable: the first is where a calculation starts, the
    second is what a paper quotes. Read the same, they were confusable."""
    by_hand = rg.Record(level='x', energy=-1.0, imaginary=0,
                        source={'kind': 'editor', 'gesture': 'dragged'})
    computed = rg.Record(level='x', energy=-1.0, imaginary=0,
                         source={'kind': 'run', 'run': 'runs/n01_r2scan'})
    scanned = rg.Record(level='x', energy=-1.0, source={'kind': 'scan'})
    assert tab.said_origin(by_hand) == 'by hand'
    assert tab.said_origin(computed) == 'computed'
    assert tab.said_origin(scanned) == 'from a scan'
    assert 'by hand' in tab.said_energy(by_hand)
    assert 'computed' in tab.said_energy(computed)


def test_a_structure_nobody_has_checked_is_not_called_a_minimum(tmp_path):
    """None is a third state and not the same as zero."""
    assert 'not checked' in tab.said_energy(rg.Record(level='x', energy=-1.0))
    assert 'minimum' in tab.said_energy(rg.Record(level='x', energy=-1.0,
                                                  imaginary=0))
    assert 'not checked' not in tab.said_energy(
        rg.Record(level='x', energy=-1.0, imaginary=0))


def test_a_geometry_with_no_energy_counts_as_missing(tmp_path):
    """It is not nothing -- it is what a calculation starts from -- and it is
    nothing a diagram can use. Counted as present, the summary said one thing
    while every arithmetic function said another, on one screen."""
    graph = rg.create(tmp_path / 'one', name='one')
    a = rg.add_state(graph, _WATER,
                     rg.Record(level='GFN2-xTB', method='gfn2',
                               source={'kind': 'editor', 'gesture': 'dragged'}))
    b = rg.add_state(graph, _CO, rg.Record(level='GFN2-xTB', free_energy=-5.03))
    assert rg.best(graph.node(a.id), 'GFN2-xTB') is not None, 'it is there'
    assert rg.priced(graph.node(a.id), 'GFN2-xTB') is None, 'and unusable'
    assert rg.missing_at(graph, 'GFN2-xTB') == [a.id]
    assert 'No energy at GFN2-xTB for: ' + a.id in tab.summary_html(
        graph, 'GFN2-xTB')


def test_an_unpriced_step_says_which_end_is_missing(tmp_path):
    """"Not at this level" sent a reader off to compute the saddle when what
    was unpriced was the state the step leaves."""
    graph = rg.create(tmp_path / 'one', name='one')
    a = rg.add_state(graph, _WATER, rg.Record(level='x', method='gfn2'),
                     label='educt')
    b = rg.add_state(graph, _CO, rg.Record(level='x', free_energy=-5.03))
    e = rg.add_transition(graph, _HCN,
                          rg.Record(level='x', free_energy=-4.9, imaginary=1),
                          source=a.id, target=b.id)
    said = dict(tab.network_lines(graph, 'x'))[e.id]
    assert 'no energy yet for educt' in said
    assert 'the saddle' not in said, 'the saddle is priced; the educt is not'

    # And the other way round.
    graph2 = rg.create(tmp_path / 'two', name='two')
    c = rg.add_state(graph2, _WATER, rg.Record(level='x', free_energy=-5.0))
    d = rg.add_state(graph2, _CO, rg.Record(level='x', free_energy=-5.03))
    f = rg.add_transition(graph2, _HCN, rg.Record(level='x', imaginary=1),
                          source=c.id, target=d.id)
    assert 'no energy yet for the saddle' in dict(
        tab.network_lines(graph2, 'x'))[f.id]


def test_the_panel_says_a_geometry_is_there_and_unpriced(tmp_path):
    """Which is a different sentence from "nothing here", and a different
    next action."""
    graph = rg.create(tmp_path / 'one', name='one')
    a = rg.add_state(graph, _WATER, rg.Record(level='GFN2-xTB', method='gfn2'))
    said = tab.detail_html(graph, a.id, 'GFN2-xTB')
    assert 'geometry at' in said and 'no energy for it' in said
    assert 'what a calculation would start from' in said


# ---------------------------------------------------------------------------
# On to a calculation
# ---------------------------------------------------------------------------

def test_a_geometry_goes_to_the_submit_tab_with_a_folder_inside_the_graph(
        tmp_path):
    """The job is set up over there; every field the Submit tab has is what a
    calculation needs and a second copy here would be a smaller, staler one.
    What this does is give the job somewhere to come back to."""
    graph, a, *_ = _mechanism(tmp_path)
    ctx = _Dashboard(tmp_path)
    ctx.submit_refs['job_name_widget'] = widgets.Text(value='')
    ctx.tabs_widget.selected_index = 5
    panel = tab.ReactionGraphPanel(tmp_path, ctx=ctx)
    panel.network.value = a.id
    panel._on_to_submit()

    written = ctx.submit_refs['coords_widget'].value
    assert written.splitlines()[0] == '3'
    assert 'for a calculation' in written.splitlines()[1]
    assert ctx.tabs_widget.selected_index == 0

    # The folder is inside the graph, not out in a shared calculation
    # directory: a graph whose evidence is scattered cannot be handed to
    # anybody.
    runs = sorted((graph.folder / rg.RUNS).iterdir())
    assert len(runs) == 1
    assert runs[0].parent == graph.folder / rg.RUNS
    assert (runs[0] / 'from_graph.xyz').read_text(encoding='utf-8')
    assert ctx.submit_refs['job_name_widget'].value == runs[0].name


def test_sending_one_for_a_calculation_is_written_into_the_history(tmp_path):
    graph, a, *_ = _mechanism(tmp_path)
    ctx = _Dashboard(tmp_path)
    panel = tab.ReactionGraphPanel(tmp_path, ctx=ctx)
    panel.network.value = a.id
    panel._on_to_submit()
    line = [one for one in rg.history(rg.load(graph.folder))
            if one['what'] == 'sent for a calculation']
    assert line and line[0]['ref'] == a.id
    assert line[0]['run'].startswith(rg.RUNS + '/')


def test_two_calculations_for_the_same_thing_get_their_own_folders(tmp_path):
    """A calculation is evidence. Two sharing a directory means neither can be
    pointed at afterwards."""
    graph, a, *_ = _mechanism(tmp_path)
    ctx = _Dashboard(tmp_path)
    panel = tab.ReactionGraphPanel(tmp_path, ctx=ctx)
    panel.network.value = a.id
    panel._on_to_submit()
    panel._on_to_submit()
    assert len(list((graph.folder / rg.RUNS).iterdir())) == 2


def test_with_no_submit_tab_there_is_no_way_to_start_a_calculation(tmp_path):
    _mechanism(tmp_path)
    panel = _panel(tmp_path)
    assert panel.to_submit_btn.layout.display == 'none'


def test_a_state_arrives_under_a_name_you_can_read(tmp_path):
    """A network of n01, n02, n03 cannot be read without clicking each one. A
    formula is a fact rather than a guess about the species, and the label is
    the user's to change the moment they know better."""
    rg.create(tmp_path / 'fresh', name='fresh')
    panel = _panel(tmp_path)
    panel.take(_offer(_NH3))
    shown = dict((ref, line) for line, ref in panel.network.options)
    assert 'H3N' in list(shown.values())[0], shown

    # And the two ends of a walk arrive named too.
    panel.take(_offer(_HCN, ends=(_WATER, _CO), imaginary=1))
    labels = [n.label for n in rg.load(panel.graph.folder).nodes]
    assert labels == ['H3N', 'H2O', 'CO']


def test_a_name_the_user_gives_wins_over_the_formula(tmp_path):
    rg.create(tmp_path / 'fresh', name='fresh')
    panel = _panel(tmp_path)
    panel.take(_offer(_NH3, name='the amine'))
    assert rg.load(panel.graph.folder).nodes[0].label == 'the amine'
