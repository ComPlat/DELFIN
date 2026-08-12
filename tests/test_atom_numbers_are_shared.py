"""Atom numbers belong to a viewer, so every tab that has one can show them.

They lived inside the ORCA Builder and were reachable from nowhere else: the
size control resolved the viewer as ``window._orcaBuildViewer``, one global for
one tab. Sharing the code was not enough -- that global had to go first.

Now the Submit tab has them, from the same code and with the same five sizes.
Driven in a browser on a ten-atom glycine: the ``#`` button is in the toolbar,
switching it on takes the viewer from 0 to 10 labels at scale 0.20, the size
control appears with it, and choosing XXL then S moves the live scale to 0.38
and 0.11 without the numbers being rebuilt.

The ORCA Builder is untouched by this: its labels are byte-for-byte what they
were, diffed against its own previous function.
"""

from delfin.dashboard import structure_editor, tab_orca_builder, tab_submit


SUBMIT = open(tab_submit.__file__, encoding='utf-8').read()
ORCA = open(tab_orca_builder.__file__, encoding='utf-8').read()


def test_the_numbering_names_no_tab():
    """A viewer is passed in. Nothing here may reach for one by name."""
    emitted = (structure_editor.atom_number_labels_js(
                   '2\n\nO 0 0 0\nH 1 0 0\n', var='v')
               + structure_editor.label_scale_setter_js())

    assert '_orcaBuildViewer' not in emitted
    assert '_submitMolViewerByScope' not in emitted
    assert 'function(scale, viewer)' in structure_editor.label_scale_setter_js()


def test_both_tabs_take_it_from_the_same_place():
    assert 'structure_editor.atom_number_labels_js(' in ORCA
    assert '_structure_editor.atom_number_labels_js(' in SUBMIT
    # And neither keeps a copy of the occlusion pass.
    assert 'grid=Object.create(null)' not in ORCA
    assert 'grid=Object.create(null)' not in SUBMIT


def test_the_sizes_are_the_same_five():
    assert [name for name, _ in structure_editor.LABEL_SIZES] == [
        'S', 'M', 'L', 'XL', 'XXL']
    assert structure_editor.LABEL_SCALE_DEFAULT == 0.20
    assert 'options=list(_structure_editor.LABEL_SIZES)' in SUBMIT


def test_nothing_is_numbered_until_it_is_asked_for():
    made = SUBMIT.split('submit_labels_btn = widgets.ToggleButton')[1].split(')\n')[0]
    assert 'value=False' in made, 'numbers are off until switched on'

    build = SUBMIT.split('def _build_mol_output_bundle')[1].split('\n    def ')[0]
    assert 'if submit_labels_btn.value:' in build


def test_resizing_does_not_rebuild_the_molecule():
    """The browser rescales the sprites it already holds; re-rendering would
    throw the camera away and cost a WebGL context."""
    handler = SUBMIT.split('def on_submit_label_size')[1].split('\n    def ')[0]

    assert '__delfinSetLabelScale' in handler
    assert 'update_molecule_view' not in handler
    # Switching them on or off is the other way round: the labels are built
    # with the model, so that one does render again.
    toggle = SUBMIT.split('def on_submit_labels_toggle')[1].split('\n    def ')[0]
    assert 'update_molecule_view()' in toggle


def test_an_unreadable_structure_is_not_numbered():
    assert structure_editor.atom_number_labels_js('') == ''
    assert structure_editor.atom_number_labels_js('nonsense') == ''
    assert structure_editor.atom_number_labels_js('2\n\nO 0 0 0\nH 1 0 0\n')
