"""The workbench putting something into the document, and neither knowing the other.

The editor knows what is on screen, what made it, and what the last press left
lying about.  The reaction graph knows what a species is, when two structures
may be the same one, and what a barrier is allowed to be compared against.
Neither of those belongs inside the other: a workbench holding the network's
rules would be holding a second copy of them, and two copies of a rule is one
rule that will disagree with itself.

So the hand-over is a dict and two callbacks.  The editor describes what it
has; the graph side answers what -- if anything -- can be done with it, and
that answer is what is written on the button.  No graph open, or nothing
correct to do, and there is no button at all.

These tests drive the real editor.  Which control is on screen, and what a
press does with the boxes beside it set that way, is not a question the source
can answer, because the source says what it means to do.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest

_ETHANE = """8
ethane
C   -0.7560000000   0.0000000000   0.0000000000
C    0.7560000000   0.0000000000   0.0000000000
H   -1.1404000000   1.0206000000   0.0000000000
H   -1.1404000000  -0.5103000000   0.8838000000
H   -1.1404000000  -0.5103000000  -0.8838000000
H    1.1404000000   0.5103000000   0.8838000000
H    1.1404000000   0.5103000000  -0.8838000000
H    1.1404000000  -1.0206000000   0.0000000000
"""

_WATER = """3
water
O 0.000 0.000 0.000
H 0.957 0.000 0.000
H -0.240 0.927 0.000
"""


class _Graph:
    """A graph side that records what it was asked, and answers as told."""

    def __init__(self, label='Put in the network'):
        self.label = label
        self.asked = []
        self.taken = []

    def offer_label(self, offer):
        self.asked.append(offer)
        return self.label

    def take(self, offer):
        self.taken.append(offer)
        return 'it went in'


def _an_editor(text=_ETHANE, graph=None):
    """One structure editor, with or without a network behind it."""
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / 'calc', archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda _script: None
    state = {}
    extra = {}
    if graph is not None:
        extra = {'graph_offer': graph.offer_label,
                 'put_in_graph': graph.take}
    part = structure_editor.build(
        ctx, state=state, coords_widget=widgets.Textarea(value=text),
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None, **extra)
    return part, state


def _shown(widget):
    return str(getattr(widget.layout, 'display', '') or '') != 'none'


# ---------------------------------------------------------------------------
# The button is the graph side's answer
# ---------------------------------------------------------------------------

def test_with_no_network_behind_it_there_is_no_button():
    """Which is how the editor stands alone: the Reactions tab may never be
    opened, and a control that reaches a tab that is not there is the row
    claiming something untrue."""
    part, _state = _an_editor()
    part._refresh_graph_button()
    assert not _shown(part.submit_graph_btn)


def test_the_graph_side_writes_what_the_press_will_do():
    graph = _Graph('Put both ends and the step in Mn hydrogenation')
    part, _state = _an_editor(graph=graph)
    part._refresh_graph_button()
    assert _shown(part.submit_graph_btn)
    assert part.submit_graph_btn.description == \
        'Put both ends and the step in Mn hydrogenation'


def test_a_graph_side_that_answers_nothing_hides_the_button():
    """Nothing correct to do is an answer, and an absent control is how this
    row says it -- the same rule everything else in it obeys."""
    graph = _Graph(None)
    part, _state = _an_editor(graph=graph)
    part._refresh_graph_button()
    assert not _shown(part.submit_graph_btn)


def test_a_graph_side_that_raises_does_not_take_the_editor_with_it():
    class _Broken(_Graph):
        def offer_label(self, offer):
            raise RuntimeError('the document is on a disk that went away')

    part, _state = _an_editor(graph=_Broken())
    part._refresh_graph_button()
    assert not _shown(part.submit_graph_btn)


def test_the_button_is_asked_again_and_never_remembers():
    """A graph is opened and closed in another tab while this one stands
    still. A button that remembered would write into a document nobody has
    open."""
    graph = _Graph('Put in one')
    part, _state = _an_editor(graph=graph)
    part._refresh_graph_button()
    assert _shown(part.submit_graph_btn)
    graph.label = None
    part._refresh_graph_button()
    assert not _shown(part.submit_graph_btn)


# ---------------------------------------------------------------------------
# What is handed over
# ---------------------------------------------------------------------------

def test_the_offer_carries_the_structure_and_the_boxes_as_they_are_set():
    graph = _Graph()
    part, _state = _an_editor(graph=graph)
    part.submit_ff_dd.value = 'gfn2'
    part.submit_gfn_charge.value = -2
    part.submit_gfn_mult.value = 3
    offer = part._graph_offer()
    assert offer['charge'] == -2
    assert offer['multiplicity'] == 3
    assert offer['method'] == 'gfn2'
    assert 'C ' in offer['xyz'] or 'C  ' in offer['xyz']


def test_the_level_is_written_the_way_a_caption_would_carry_it():
    """That string is what a diagram is drawn at. Two records whose levels
    differ only in a solvent nobody wrote down make a network that cannot be
    defended."""
    graph = _Graph()
    part, _state = _an_editor(graph=graph)
    part.submit_ff_dd.value = 'gfn2'
    assert part._graph_offer()['level'] == 'GFN2-xTB'
    part.submit_gfn_solvent.value = 'thf'
    assert part._graph_offer()['level'] == 'GFN2-xTB/thf'


def test_an_energy_belongs_to_the_geometry_it_was_measured_on():
    """Carried over quietly it is the difference between two unrelated
    numbers -- and here it would be written into a document as a fact about a
    species."""
    graph = _Graph()
    part, state = _an_editor(graph=graph)
    here = part.coords_widget.value
    state['gfn_energy'] = -7.331
    state['gfn_energy_for'] = part._structure_fingerprint(here)
    assert part._graph_offer()['energy'] == -7.331

    part.coords_widget.value = _WATER
    assert part._graph_offer()['energy'] is None, \
        'the budget cannot outlive the molecule, and neither can this'


def test_what_the_last_press_said_the_structure_is_travels_with_it():
    graph = _Graph()
    part, state = _an_editor(graph=graph)
    here = part.coords_widget.value
    state['saddle_modes'] = {'xyz': here, 'under': 'gfn2',
                             'modes': {'count': 1, 'modes': [-412.3],
                                       'real': [100.0]}}
    offer = part._graph_offer()
    assert offer['imaginary'] == 1
    assert offer['frequency'] == pytest.approx(-412.3)

    # And it does not travel to a different structure.
    part.coords_widget.value = _WATER
    assert part._graph_offer()['imaginary'] is None


def test_the_two_ends_a_walk_left_are_part_of_the_offer():
    graph = _Graph()
    part, state = _an_editor(graph=graph)
    state['scan_ends'] = (_WATER, _ETHANE)
    assert part._graph_offer()['ends'] == (_WATER, _ETHANE)


def test_the_editor_s_own_comment_line_is_the_account_of_how_it_came_about():
    """Every press in the editor writes what it did on the comment line, so it
    is already the account, and copying it is truer than a label chosen here.
    """
    graph = _Graph()
    part, _state = _an_editor(graph=graph)
    part.coords_widget.value = (
        '3\noptimised to a transition state with GFN2-xTB\n'
        'O 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n')
    assert part._graph_offer()['gesture'] == \
        'optimised to a transition state with GFN2-xTB'


def test_nothing_on_screen_is_nothing_to_offer():
    graph = _Graph()
    part, _state = _an_editor(graph=graph)
    part.coords_widget.value = ''
    assert part._graph_offer() is None
    part._refresh_graph_button()
    assert not _shown(part.submit_graph_btn)


# ---------------------------------------------------------------------------
# The press
# ---------------------------------------------------------------------------

def test_the_press_hands_it_over_and_says_what_came_back():
    graph = _Graph()
    part, _state = _an_editor(graph=graph)
    part._refresh_graph_button()
    part.on_submit_graph()
    assert len(graph.taken) == 1
    assert graph.taken[0]['xyz'].startswith('8')
    assert 'it went in' in part.mol_status.value


def test_a_hand_over_that_fails_is_reported_and_not_swallowed():
    class _Broken(_Graph):
        def take(self, offer):
            raise OSError('read-only filesystem')

    part, _state = _an_editor(graph=_Broken())
    part._refresh_graph_button()
    part.on_submit_graph()
    assert 'did not go into the graph' in part.mol_status.value
    assert 'read-only filesystem' in part.mol_status.value


def test_the_button_is_asked_again_after_the_press():
    """What was a new state a moment ago is a record on an existing one now,
    and the press has to stop saying the first thing."""
    graph = _Graph('Put in one as a new state')
    part, _state = _an_editor(graph=graph)
    part._refresh_graph_button()
    graph.label = 'Add as a GFN2-xTB record on n01'
    part.on_submit_graph()
    assert part.submit_graph_btn.description == \
        'Add as a GFN2-xTB record on n01'


def test_the_button_is_on_the_toolbar_and_goes_quiet_while_a_run_holds_it():
    graph = _Graph()
    part, _state = _an_editor(graph=graph)
    assert part.submit_graph_btn in part.submit_manip_toolbar.children
    part._set_manip_toolbar_enabled(False)
    assert part.submit_graph_btn.disabled
    part._set_manip_toolbar_enabled(True)
    assert not part.submit_graph_btn.disabled


def test_only_the_editor_s_own_comment_counts_as_an_account():
    """An XYZ comment is free text and xtb writes its own there. Measured on a
    real optimisation the line reads ``energy: -5.504066183163 gnorm:
    0.000283898985 xtb: 6.7.1`` -- copied into a document as how a species came
    about, that is noise wearing the shape of provenance."""
    graph = _Graph()
    part, _state = _an_editor(graph=graph)
    part.coords_widget.value = (
        '3\n energy: -5.504066183163 gnorm: 0.000283898985 xtb: 6.7.1\n'
        'O 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n')
    assert part._graph_offer()['gesture'] == 'the editor'

    part.coords_widget.value = (
        '3\nOptimised to a transition state\n'
        'O 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n')
    assert part._graph_offer()['gesture'] == 'Optimised to a transition state'


def test_a_structure_arriving_any_way_at_all_offers_itself_to_the_network():
    """Including one drawn in Ketcher, and that is why there is no second path
    for drawings.

    Ketcher hands back a molfile, the editor turns it into a SMILES and writes
    it into the input box, and Convert turns that into coordinates -- the same
    chain every other conversion goes through. What makes the network reachable
    from a drawing is that ``_replace_mol_output_view`` refreshes the row when
    a structure is placed, and every route into the editor goes through it. A
    drawing-shaped shortcut into the graph would be a second way to do this one
    thing, and the second way is the one that goes stale.
    """
    from editor_source import EDITOR_SOURCE

    drawing = EDITOR_SOURCE.split('def _replace_mol_output_view(', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_refresh_saddle_controls()' in drawing

    saddle = EDITOR_SOURCE.split('def _refresh_saddle_controls(', 1)[1].split(
        '\n    def ', 1)[0]
    assert '_refresh_graph_button()' in saddle

    # And on a real editor: a structure placed the way every conversion places
    # one brings the press with it.
    graph = _Graph('Put in one as a new state')
    part, _state = _an_editor(graph=graph)
    part.coords_widget.value = ''
    part._refresh_graph_button()
    assert not _shown(part.submit_graph_btn)

    part.coords_widget.value = _WATER
    part._replace_mol_output_view(_WATER)
    assert _shown(part.submit_graph_btn)
    assert part.submit_graph_btn.description == 'Put in one as a new state'
