"""Atom numbers belong to a viewer, so every tab that has one can show them.

They lived inside the ORCA Builder and were reachable from nowhere else: the
size control resolved the viewer as ``window._orcaBuildViewer``, one global for
one tab. Sharing the code was not enough -- that global had to go first.

Three things were wrong with the first version of them, and all three come from
the same mistake: the numbers were written out as fixed points in space, from
coordinates the kernel had, instead of being hung on the atoms the browser
holds.

  * They stayed behind when the atoms moved. A drag, an optimisation or a
    dynamic run left every number where the atom used to be.
  * Switching them off rendered the molecule again -- and a molecule rebuilt
    from its coordinates comes back with the bonds re-perceived from distances,
    so a structure that had been pulled about lost bonds. The numbers did what
    dynamic bonds do, without anybody asking for them.
  * They were too small to read, and the largest size on offer was where the
    smallest should have started.

Now they are a layer of sprites over the model that holds the atom objects
themselves: they follow the cores through every frame the editor draws, they go
on and off without the model hearing about it, and the ladder of sizes starts
where it used to end.
"""

from delfin.dashboard import molecule_viewer, structure_editor, tab_orca_builder
from editor_source import SUBMIT_SOURCE as SUBMIT
from editor_source import TAB_SOURCE


ORCA = open(tab_orca_builder.__file__, encoding='utf-8').read()
VIEWER = open(molecule_viewer.__file__, encoding='utf-8').read()
LAYER = structure_editor.atom_numbers_js()


def test_the_numbering_names_no_tab():
    """A viewer is passed in. Nothing here may reach for one by name."""
    emitted = structure_editor.show_atom_numbers_js(var='v')

    assert '_orcaBuildViewer' not in emitted
    assert '_submitMolViewerByScope' not in emitted
    assert 'window.__delfinAtomNumbers.set(v,true,' in emitted


def test_both_tabs_take_it_from_the_same_place():
    assert '_structure_editor.show_atom_numbers_js(' in ORCA
    assert 'show_atom_numbers_js(' in SUBMIT
    # And no tab keeps a copy of the occlusion pass: it is written once, in
    # the part both of them take their numbering from.
    assert 'grid=Object.create(null)' not in ORCA
    assert 'grid=Object.create(null)' not in TAB_SOURCE


def test_the_layer_is_installed_once_per_page():
    assert LAYER.startswith('if(!window.__delfinAtomNumbers){')


def test_a_number_is_hung_on_its_atom_not_on_a_copy_of_its_place():
    """The atom object is kept, so reading x/y/z later reads where it is now."""
    assert 'L.push({a:a,l:lab})' in LAYER
    assert "sp.position.set(a.x,a.y,a.z)" in LAYER
    # ...and the projection used for the occlusion test reads it too, rather
    # than a stored triple.
    assert 'p[0]=r11*a.x+r12*a.y+r13*a.z;' in LAYER


def test_they_are_brought_along_on_every_frame_the_molecule_moves_in():
    """One funnel draws every changed geometry -- drag, optimise, dynamics."""
    draw = VIEWER.split('function drawHighlightsNow(scopeKey)')[1].split(
        '\n    function ')[0]

    assert 'window.__delfinAtomNumbers.refresh(viewer)' in draw
    assert draw.index('__delfinAtomNumbers.refresh') < draw.index(
        'try { viewer.render(); }')


def test_switching_them_off_leaves_the_molecule_alone():
    """This is the bond loss: a re-render re-perceives bonds from distances."""
    handler = SUBMIT.split('def on_submit_labels_toggle')[1].split('\n    def ')[0]

    assert 'update_molecule_view' not in handler
    assert 'show_atom_numbers_js(' in handler
    assert 'on=on' in handler

    # The layer takes its own sprites out. It does not reach for the model,
    # and it does not use removeAllLabels, which would take away labels that
    # are not its own -- and draw a frame for each one it removed.
    clear = LAYER.split('function clear(v){')[1].split('\n  }')[0]
    assert 'v.modelGroup.remove(lab.sprite)' in clear
    assert 'removeAllLabels' not in LAYER
    assert 'addModel' not in LAYER
    assert 'rebuildBonds' not in LAYER
    assert 'setStyle' not in LAYER


def test_resizing_does_not_rebuild_the_molecule():
    """The browser rescales the sprites it already holds; re-rendering would
    throw the camera away and cost a WebGL context."""
    handler = SUBMIT.split('def on_submit_label_size')[1].split('\n    def ')[0]

    assert '__delfinAtomNumbers.setScale' in handler
    assert 'update_molecule_view' not in handler


def test_the_size_is_asked_for_in_pixels():
    """Five rungs were not enough, and 'L' says nothing about how big that is.

    A digit comes out 34 px tall per unit of the scale factor -- measured in a
    browser at 0.28, 0.38, 0.50, 0.66 and 0.86, which gave 9.5, 13, 17, 22 and
    30 px. The device pixel ratio cancels out of that, so it is the same
    number on any screen and at any viewer size, and the control can simply
    ask for pixels.
    """
    assert structure_editor.scale_for_px(17) == 0.5
    assert structure_editor.scale_for_px(34) == 1.0
    assert structure_editor.LABEL_SCALE_DEFAULT == structure_editor.scale_for_px(
        structure_editor.LABEL_PX_DEFAULT)

    # Anything unusable falls back rather than making the numbers vanish.
    assert structure_editor.scale_for_px('') == structure_editor.LABEL_SCALE_DEFAULT
    assert structure_editor.scale_for_px(None) == structure_editor.LABEL_SCALE_DEFAULT
    assert structure_editor.scale_for_px(0) == structure_editor.scale_for_px(
        structure_editor.LABEL_PX_MIN)
    assert structure_editor.scale_for_px(10_000) == structure_editor.scale_for_px(
        structure_editor.LABEL_PX_MAX)

    # Both tabs ask the same way, so changing it changes it everywhere.
    # One control, written once, in the part both tabs build their editor
    # from -- the Builder used to keep a second pair of its own.
    assert 'widgets.BoundedIntText(' in SUBMIT
    assert 'min=LABEL_PX_MIN' in SUBMIT and 'max=LABEL_PX_MAX' in SUBMIT
    assert 'widgets.BoundedIntText(' not in ORCA
    assert 'orca_mol_labels_btn' not in ORCA


def test_nothing_is_numbered_until_it_is_asked_for():
    made = SUBMIT.split('submit_labels_btn = widgets.ToggleButton')[1].split(')\n')[0]
    assert 'value=False' in made, 'numbers are off until switched on'

    build = SUBMIT.split('def _build_mol_output_bundle')[1].split('\n    def ')[0]
    assert 'if submit_labels_btn.value:' in build


def test_every_structure_in_the_viewer_is_numbered_from_its_own_zero():
    """The ORCA Builder overlays two structures in one viewer."""
    assert 'function modelsOf(v)' in LAYER
    build = LAYER.split('function build(v,scale){')[1].split('\n  }')[0]
    assert 'for(var mi=0;mi<ms.length;mi++)' in build
    assert "v.addLabel(String(i)" in build


def test_atoms_coming_or_going_are_numbered_again():
    """Drawing an atom or adjusting hydrogens changes what there is to number."""
    refresh = LAYER.split('function refresh(v){')[1].split('\n  }')[0]

    assert 'atomCount(v)!==m' in refresh
    # ...without the rebuild being able to ask for itself again.
    assert '!v.__delfinLabelRebuilding' in refresh
