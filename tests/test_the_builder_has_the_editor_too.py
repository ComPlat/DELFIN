"""The ORCA Builder holds a structure editor, and it is the same one.

Not a copy of it: ``structure_editor.build()`` is called twice, once by the
Submit tab and once here, so a fix to the editor is a fix in both and there is
one place to make it.

What the Builder had to bring to the arrangement is a translator. The editor
writes a structure as plain XYZ into a coordinates box and lets its tab take it
from there; this tab's box holds named blocks -- ``name.xyz;comment``, the
atoms, a closing star -- and there may be several of them with only one on
screen. So the editor gets a box of its own, and what it writes is folded into
the block that is being shown, through the same record rebuild Apply Numbering
Fix does. Every other block, and every header, is left exactly as it was.

Driven through the real tab on two blocks -- a benzene and a water:

    ▶                block 2 shown, the editor handed its 3 atoms
    editor writes    the water block carries the new coordinates,
                     the benzene block is untouched, both names kept,
                     and the view stays on block 2
    Check Numbering  runs, offers a fix
    Apply Fix        rewrites block 2 and nothing else

One real bug came with the second editor, and it is in this file's second half:
the browser looked its toolbar parts up by falling back to the whole page.
"""

import tempfile
import time
from pathlib import Path

import pytest

from delfin.dashboard import molecule_viewer, structure_editor, tab_orca_builder


VIEWER = open(molecule_viewer.__file__, encoding='utf-8').read()
BUILDER = open(tab_orca_builder.__file__, encoding='utf-8').read()

BENZENE = """C       -0.461862     0.762637    -0.238321
C       -0.838289     0.154318     0.927507
C        0.092243    -0.385321     1.788509
C        1.417041    -0.289502     1.431542
C        1.837785     0.320768     0.256083
C        0.881377     0.859473    -0.599674
H       -1.187136     1.189385    -0.920678
H       -1.902087     0.101819     1.169357
H       -0.256407    -0.857624     2.701396
H        2.142563    -0.706502     2.095963
H        2.885967     0.388017    -0.010677
H        1.160238     1.344948    -1.527129"""
WATER = "O 0.000 0.000 0.000\nH 0.757 0.586 0.000\nH -0.757 0.586 0.000"
TWO_BLOCKS = (f"benzene.xyz;the ring\n12\n\n{BENZENE}\n*\n\n"
              f"water.xyz;the other one\n3\n\n{WATER}\n*")


@pytest.fixture
def builder():
    pytest.importorskip('ipywidgets')
    from delfin.dashboard.context import DashboardContext

    tmp = Path(tempfile.mkdtemp())
    sent = []
    ctx = DashboardContext(calc_dir=tmp / 'calc', archive_dir=tmp / 'archive',
                           office_dir=tmp / 'office')
    ctx.run_js = sent.append
    _tab, refs = tab_orca_builder.create_tab(ctx)
    refs['orca_coords'].value = TWO_BLOCKS
    sent.clear()
    return refs, sent


def test_it_is_the_same_editor_and_not_another_one():
    assert '_structure_editor.build(' in BUILDER
    # ...and the Builder writes none of it out again.
    for own in ('submit_manip_toolbar = widgets.HBox', 'def on_submit_optimize',
                'def _apply_structure', 'def _perception_for'):
        assert own not in BUILDER


def test_the_toolbar_is_there_with_everything_on_it(builder):
    refs, _sent = builder
    toolbar = refs['submit_manip_toolbar']

    assert toolbar.layout.display == 'flex'
    on_it = set(toolbar.children)
    for name in ('submit_select_btn', 'submit_manip_btn', 'submit_draw_btn',
                 'submit_ff_dd', 'submit_optimize_btn', 'submit_optimize_all_btn',
                 'submit_relax_btn', 'submit_settle_btn', 'submit_bond_btn',
                 'submit_unbond_btn', 'submit_dyn_bonds_btn', 'submit_poly_dd',
                 'submit_hyb_dd', 'submit_strength_slider', 'submit_sens_slider',
                 'submit_gfn_charge', 'submit_gfn_mult', 'submit_gfn_solvent',
                 'submit_labels_btn', 'submit_reset_btn', 'submit_manip_undo_btn'):
        assert refs[name] in on_it, name
    # The whole ladder of methods, browser and server alike.
    assert [value for _label, value in refs['submit_ff_dd'].options][:5] == [
        'uff', 'mmff94', 'gfnff', 'gfn2', 'gxtb']


def test_the_viewer_tells_the_editor_it_is_there(builder):
    """The editor addresses a viewer by the scope it belongs to. Without this
    the toolbar is a row of buttons with nothing under them."""
    refs, sent = builder
    refs['orca_mol_next_btn'].click()

    registrations = [s for s in sent if '_submitMolViewerByScope' in s]
    assert registrations, 'the viewer never introduced itself'
    assert refs['editor_scope'] in registrations[0]
    assert 'onViewerReady' in registrations[0]


def test_stepping_blocks_hands_the_editor_the_one_on_screen(builder):
    refs, _sent = builder
    assert refs['editor_coords'].value.split('\n')[0] == '12'

    refs['orca_mol_next_btn'].click()

    assert refs['editor_coords'].value.split('\n')[0] == '3'
    assert refs['editor_state']['xyz_view_idx'] == 1


def test_an_edit_lands_in_the_block_on_screen_and_nowhere_else(builder):
    refs, _sent = builder
    refs['orca_mol_next_btn'].click()

    refs['editor_coords'].value = (
        '3\nEdited in DELFIN viewer\n'
        'O 0.100 0.000 0.000\nH 0.757 0.586 0.000\nH -0.757 0.586 0.000\n')

    text = refs['orca_coords'].value
    assert 'benzene.xyz;the ring' in text and 'water.xyz;the other one' in text
    assert BENZENE.split('\n')[0] in text, 'the other block was rewritten'
    assert 'O 0.100 0.000 0.000' in text
    assert text.count('*') == 2
    # And the tab did not start over: still the second block, still its view.
    assert refs['editor_state']['xyz_view_idx'] == 1
    assert [name for name, _xyz in refs['editor_state']['xyz_blocks']] == [
        'benzene.xyz', 'water.xyz']


def test_a_drag_does_not_take_the_model_apart(builder):
    """The browser has already moved the atoms and is showing them. Drawing the
    structure again would perceive its bonds from distances a second time --
    which is the one thing an editor must not do to a structure being edited.
    """
    refs, sent = builder
    refs['orca_mol_next_btn'].click()
    sent.clear()

    refs['editor_state']['manip_inflight'] = True
    refs['editor_coords'].value = (
        '3\nEdited in DELFIN viewer\n'
        'O 0.200 0.000 0.000\nH 0.757 0.586 0.000\nH -0.757 0.586 0.000\n')

    assert 'O 0.200 0.000 0.000' in refs['orca_coords'].value
    assert not [s for s in sent if 'addModel' in s], 'the model was rebuilt'


def test_the_numbering_tools_are_untouched(builder):
    """Check Numbering and Apply Numbering Fix work on the one coordinates box,
    and the editor writes into it the same way they do."""
    refs, _sent = builder
    swapped = "O 0.000 0.000 0.000\nH -0.757 0.586 0.000\nH 0.757 0.586 0.000"
    refs['orca_coords'].value = f"a.xyz;\n3\n\n{WATER}\n*\n\nb.xyz;\n3\n\n{swapped}\n*"

    refs['orca_check_numbering_btn'].click()
    results = refs['editor_state'].get('numbering_check_results') or {}
    assert results.get(1, {}).get('reordered_target_xyz')

    refs['editor_state']['numbering_check_block_idx'] = 1
    refs['orca_apply_numbering_btn'].click()

    text = refs['orca_coords'].value
    assert 'a.xyz;' in text and 'b.xyz;' in text
    fixed = text.split('b.xyz;')[1]
    assert fixed.index('H  -0.75700000') < fixed.index('H   0.75700000')


def test_the_two_editors_do_not_share_a_scope(builder):
    pytest.importorskip('ipywidgets')
    from delfin.dashboard import tab_submit
    from delfin.dashboard.context import DashboardContext

    refs, _sent = builder
    tmp = Path(tempfile.mkdtemp())
    ctx = DashboardContext(calc_dir=tmp / 'calc', archive_dir=tmp / 'archive',
                           office_dir=tmp / 'office')
    ctx.run_js = lambda _script: None
    _tab, submit_refs = tab_submit.create_tab(ctx)

    assert refs['editor_scope'] != submit_refs['submit_scope_id']


# ---------------------------------------------------------------------------
# the bug a second editor brings with it
# ---------------------------------------------------------------------------


def test_a_control_is_looked_for_in_every_root_of_its_own_scope():
    """It used to be looked for in the first one, and then anywhere on the page.

    The fallback was there because fullscreen moves the toolbar into a floating
    overlay, and the overlay is a second element with the scope's class; a
    lookup that took only the first found nothing. It was written down as safe
    on the grounds that there is one Submit tab per dashboard -- true until the
    Builder had an editor too, and then it would reach into the other tab's
    toolbar and read the wrong value box, or send an edit through the wrong
    sync input.
    Driven in a browser on a page with two scopes, one of them with its
    toolbar lifted into a fullscreen overlay carrying the same class::

        A  .submit-manip-sync      -> 'from A'
        B  .submit-manip-sync      -> 'from B'
        B  .submit-internal-value  -> 'B in fullscreen'
        A  .submit-internal-value  -> None

    That last one is the point: A is told there is none, rather than handed
    B's.
    """
    found = VIEWER.split('function findInScope(scopeKey, selector)')[1].split(
        '\n    }')[0]

    assert "document.querySelectorAll('.' + scopeKey)" in found
    assert 'roots[i].querySelector(selector)' in found
    # No page-wide fallback of any kind.
    assert 'document.querySelector(selector)' not in found


def test_the_part_carries_its_own_scope_class(builder):
    """Which is what makes the lookup above find anything at all."""
    refs, _sent = builder
    scope = refs['editor_scope']

    assert scope in BUILDER or 'orca_editor_scope' in BUILDER
    assert 'orca_mol_module.add_class(orca_editor_scope)' in BUILDER
    assert scope.startswith('submit-scope-'), (
        'the scope comes from the part, so both tabs make theirs the same way')


def test_the_editor_seeds_the_keys_it_reads_without_asking():
    """A second host hands it a dictionary the editor has not been given a tour
    of, and six keys were read straight out of it."""
    build = structure_editor.build
    source = open(structure_editor.__file__, encoding='utf-8').read()
    seeded = source.split("for _key, _value in (")[1].split('):')[0]

    for key in ('manip_bootstrap_done', 'perceived', 'poly_applied',
                'poly_metal', 'gfn_generation', 'gfn_scanned_uhf'):
        assert f"'{key}'" in seeded, key
    assert 'state.setdefault(_key, _value)' in source
    assert build is not None


# ---------------------------------------------------------------------------
# the numbering check, with the editor in the tab
# ---------------------------------------------------------------------------


class _Drawn(list):
    """What is in the preview, as text.

    Read off the Output widget rather than caught at display(): the preview is
    filled by assigning its outputs, because a conversion answers from a thread
    where ``with output:`` reaches nothing.
    """

    def __init__(self, widget):
        super().__init__()
        self._widget = widget

    def _current(self):
        out = []
        for item in self._widget.outputs or ():
            data = item.get('data') or {}
            out.append(str(data.get('text/html') or item.get('text') or ''))
        return out

    def __len__(self):
        return len(self._current())

    def __getitem__(self, index):
        return self._current()[index]

    def __bool__(self):
        return bool(self._current())

    def clear(self):
        return None


def _shown(monkeypatch):
    return None


SWAPPED = "O 0.000 0.000 0.000\nH -0.757 0.586 0.000\nH 0.757 0.586 0.000"


@pytest.fixture
def compared(builder, capsys):
    refs, sent = builder
    drawn = _Drawn(refs['orca_mol_output'])
    refs['orca_coords'].value = (f"a.xyz;\n3\n\n{WATER}\n*\n\n"
                                 f"b.xyz;\n3\n\n{SWAPPED}\n*")
    refs['orca_check_numbering_btn'].click()
    capsys.readouterr()
    return refs, sent, drawn


def test_the_overlay_is_still_red_and_blue(compared):
    """Two structures in one viewer, the reference in red and the target in
    blue, which is what makes a numbering mismatch visible at a glance."""
    _refs, _sent, drawn = compared

    assert drawn, 'nothing was drawn'
    overlay = drawn[-1]
    assert overlay.count('addModel') >= 2
    assert '#d32f2f' in overlay and '#1f5fff' in overlay


def test_all_three_views_belong_to_the_editor(compared):
    """The overlay included: the numbers are the editor's, and it can only put
    them on a viewer it has been told about."""
    refs, _sent, drawn = compared
    assert '_submitMolViewerByScope' in drawn[-1]
    assert refs['editor_scope'] in drawn[-1]

    for _step in range(2):
        drawn.clear()
        refs['orca_mol_next_btn'].click()
        assert '_submitMolViewerByScope' in drawn[-1]


def test_a_single_structure_can_be_edited_and_the_overlay_cannot(compared):
    """Stepping to the aligned reference or the reordered target gives the
    editor one structure, which it works on like any other. The overlay is two
    at once and there is nothing single to edit."""
    refs, _sent, drawn = compared
    assert refs['submit_manip_toolbar'].layout.display == 'none'
    assert 'Overlay' in refs['mol_status'].value

    refs['orca_mol_next_btn'].click()          # aligned reference
    assert refs['submit_manip_toolbar'].layout.display == 'flex'
    assert refs['editor_coords'].value.split('\n')[0] == '3'

    refs['orca_mol_next_btn'].click()          # reordered target
    assert refs['submit_manip_toolbar'].layout.display == 'flex'


def test_the_comparison_stays_up_after_the_fix(compared):
    """The reordered block is the thing to look at once it has been applied.
    Rewriting the coordinates box the ordinary way dropped the tab back to the
    first structure with nothing left to compare against."""
    refs, _sent, _drawn = compared
    refs['editor_state']['numbering_check_block_idx'] = 1

    refs['orca_apply_numbering_btn'].click()

    assert refs['editor_state']['numbering_check_active']
    fixed = refs['orca_coords'].value.split('b.xyz;')[1]
    assert fixed.index('H  -0.75700000') < fixed.index('H   0.75700000')


# ---------------------------------------------------------------------------
# what is held in the editor, in the input ORCA reads
# ---------------------------------------------------------------------------


ETHANE_BLOCK = ("eth.xyz;\n8\n\nC 0 0 0\nC 1.53 0 0\n"
                "H -0.36 1.02 0\nH -0.36 -0.51 0.88\nH -0.36 -0.51 -0.88\n"
                "H 1.89 1.02 0\nH 1.89 -0.51 0.88\nH 1.89 -0.51 -0.88\n*")


def test_a_held_coordinate_is_written_as_an_orca_constraint(builder):
    """Set a bond, an angle or a dihedral in the editor and it is in the input.

    ORCA counts atoms from zero and so does the numbering on the structure, so
    the numbers in the constraint are the numbers on the atoms.
    """
    refs, _sent = builder
    refs['orca_coords'].value = ETHANE_BLOCK
    refs['editor_state']['constraints'] = [
        {'kind': 'distance', 'atoms': [0, 1], 'value': 1.6, 'mode': 'fix'},
        {'kind': 'angle', 'atoms': [2, 0, 1], 'value': 109.5, 'mode': 'fix'},
        {'kind': 'dihedral', 'atoms': [2, 0, 1, 5], 'value': 180.0, 'mode': 'fix'},
    ]

    refs['update_orca_preview']()

    text = refs['orca_preview'].value
    assert '%geom Constraints' in text
    assert '{ B 0 1 1.6000 C }' in text
    assert '{ A 2 0 1 109.5000 C }' in text
    assert '{ D 2 0 1 5 180.0000 C }' in text
    assert text.index('%geom') < text.index('* xyz')


def test_a_pull_is_not_a_constraint(builder):
    """It is a spring the browser relaxes against while a structure is dragged.
    A geometry optimisation has no such thing, and writing one as a constraint
    would claim something nobody asked for."""
    refs, _sent = builder
    refs['orca_coords'].value = ETHANE_BLOCK
    refs['editor_state']['constraints'] = [
        {'kind': 'distance', 'atoms': [0, 5], 'value': 2.2, 'mode': 'pull'}]

    refs['update_orca_preview']()

    assert '%geom' not in refs['orca_preview'].value


def test_releasing_one_takes_it_out_of_the_input_again(builder):
    refs, _sent = builder
    refs['orca_coords'].value = ETHANE_BLOCK
    refs['editor_state']['constraints'] = [
        {'kind': 'distance', 'atoms': [0, 1], 'value': 1.6, 'mode': 'fix'}]
    refs['update_orca_preview']()
    assert '%geom' in refs['orca_preview'].value

    refs['editor_state']['constraints'] = []
    refs['update_orca_preview']()

    text = refs['orca_preview'].value
    assert '%geom' not in text
    assert '* xyzfile' in text, 'the rest of the input survived'


def test_it_says_so_when_the_input_reads_another_structure(builder):
    """The input reads the first block. Holding a coordinate on the second and
    writing its atom numbers into that input without a word would be a silent
    lie about atoms ORCA will never see."""
    refs, _sent = builder
    refs['orca_mol_next_btn'].click()
    refs['editor_state']['constraints'] = [
        {'kind': 'distance', 'atoms': [0, 1], 'value': 1.1, 'mode': 'fix'}]

    refs['update_orca_preview']()

    text = refs['orca_preview'].value
    assert '# Held in the editor on water.xyz' in text
    assert '# input reads benzene.xyz' in text


def test_the_fullscreen_button_stands_where_the_submit_tab_has_it(builder):
    """First in the toolbar, before Select, and not in a row of its own.

    It used to sit in the row with the block stepper, and that row is hidden
    when there is nothing to step to -- so pasting a plain XYZ, which is what
    the editor itself writes back, left no way at all to enlarge the viewer
    while the same molecule in named blocks could be enlarged fine.

    It is the editor's own button now, made over to this tab's fullscreen: the
    Submit tab's overlay is built from members carrying submit-fs-* classes,
    which the Builder's are not.
    """
    refs, _sent = builder
    button = refs['submit_fullscreen_btn']
    toolbar = list(refs['submit_manip_toolbar'].children)

    assert toolbar.index(button) == 0
    assert toolbar.index(button) < toolbar.index(refs['submit_select_btn'])
    assert 'delfin-structure-fullscreen-btn' in button._dom_classes
    assert 'submit-fullscreen-btn' not in button._dom_classes, (
        'both machineries would answer the same click')

    refs['orca_coords'].value = f'3\nwater\n{WATER}\n'
    assert refs['submit_manip_toolbar'].layout.display == 'flex'
    assert refs['orca_mol_nav_row'].layout.display == 'none', (
        'one structure has nothing to step to, and nothing else is in that row')

    refs['orca_coords'].value = TWO_BLOCKS
    assert refs['orca_mol_nav_row'].layout.display == ''
    assert refs['orca_mol_next_btn'].layout.display == ''


# ---------------------------------------------------------------------------
# where a structure comes from
# ---------------------------------------------------------------------------


def test_the_conversions_are_the_editors_too():
    """CONVERT SMILES, QUICK, + UFF and MANTA are the same four buttons the
    Submit tab has, wired to the same code."""
    assert 'orca_editor.convert_smiles_button' in BUILDER
    assert 'orca_editor.convert_smiles_quick_button' in BUILDER
    assert 'orca_editor.convert_smiles_uff_button' in BUILDER
    assert 'orca_editor.manta_button' in BUILDER
    # ...and the Builder's own one-off conversion is gone with them.
    assert 'handle_orca_convert_smiles' not in BUILDER
    assert 'orca_convert_smiles_btn' not in BUILDER


def test_several_structures_become_several_named_blocks(builder):
    """A conversion hands over one structure or several. The Submit tab shows
    the first and steps through the rest; this tab keeps them all at once, in
    the layout it reads.

    Converted through the real tab from C/C=C/C(=O)O: two conformers came back
    and the coordinates box holds two blocks, conf-1.xyz and conf-2.xyz.
    """
    refs, _sent = builder
    taken = refs['editor_state']

    handed = [('O 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0', 3, 'conf-1'),
              ('O 0 0 0\nH 0.97 0 0\nH -0.25 0.94 0', 3, 'conf-2')]
    refs['offer_structures'](handed)

    text = refs['orca_coords'].value
    assert text.count('\n*') == 2
    # Name, nothing after the semicolon: the header every other block here
    # carries. It read "conf-1.xyz;conf-1" at first, with the label twice.
    assert text.startswith('conf-1.xyz;\n')
    assert 'conf-2.xyz;\n' in text
    assert 'conf-1.xyz;conf-1' not in text
    assert [name for name, _xyz in taken['xyz_blocks']] == ['conf-1.xyz',
                                                            'conf-2.xyz']
    # The editor's own isomer stepper stays out of the way: this tab has one.
    assert refs['isomer_nav_row'].layout.display == 'none'


def test_two_structures_with_one_name_still_get_a_block_each(builder):
    refs, _sent = builder

    refs['offer_structures']([('O 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0', 3, 'quick'),
                              ('O 0 0 0\nH 0.97 0 0\nH -0.25 0.94 0', 3, 'quick')])

    names = [name for name, _xyz in refs['editor_state']['xyz_blocks']]
    assert len(names) == 2 and len(set(names)) == 2, names


def test_a_smiles_is_read_from_the_box_it_was_typed_into(builder):
    """In the Submit tab the editor writes structures into the same box a
    SMILES is typed into. Here they are two different boxes -- the editor has
    a hidden one of its own, because this tab's box holds named blocks -- so a
    conversion has to be told where to look."""
    assert 'read_input=lambda: orca_coords.value' in BUILDER

    refs, _sent = builder
    refs['orca_coords'].value = 'CCO'
    assert refs['read_input']() == 'CCO'


def test_the_drawing_editor_is_the_editors_too(builder):
    """DRAW is the same 2D editor the Submit tab opens, and what comes out of
    it lands in the box this tab shows.

    Which is not the box the editor writes structures into: it has a hidden one
    of its own here. A drawing comes back as a SMILES, and the buttons that turn
    one into coordinates read from the visible box, so that is where it goes.

    Driven through the real tab: a four-carbon chain drawn in Ketcher arrives as
    CCCC, and CONVERT SMILES turns it into conf-1.xyz and conf-2.xyz.
    """
    refs, sent = builder
    for name in ('submit_draw_open_btn', 'submit_draw_get_btn',
                 'submit_draw_update_btn', 'submit_draw_frame',
                 'submit_draw_sync'):
        assert name in refs, name
    assert 'orca_editor.submit_draw_open_btn' in BUILDER

    sent.clear()
    refs['submit_draw_get_btn'].click()
    assert len(sent) == 1 and 'getMolfile' in sent[0]

    molfile = ('\n  Ketcher\n\n  4  3  0  0  0  0            999 V2000\n'
               + '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n' * 4
               + '  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\n'
                 '  3  4  1  0  0  0  0\nM  END\n')
    refs['submit_draw_sync'].value = '1\n' + molfile

    assert refs['orca_coords'].value == 'CCCC'


def test_a_drawing_that_could_not_be_read_says_so(builder):
    refs, _sent = builder

    refs['submit_draw_sync'].value = '1\n!no-editor'

    assert 'not open yet' in refs['mol_status'].value
    assert refs['orca_coords'].value != '!no-editor'


def test_the_conversion_says_when_it_is_over(builder):
    """It ends by handing the structures to the tab, and that way out never
    cleared the line it had written.

    The other way out shows a structure, and showing one clears the status on
    the way through -- so in the Submit tab the spinner went by itself and
    here "Converting SMILES (no UFF)..." sat there for good, over coordinates
    that had plainly arrived.
    """
    refs, _sent = builder
    refs['mol_status'].value = 'Converting SMILES (no UFF)...'
    refs['editor_state']['smiles_busy'] = True

    refs['editor_offer_isomers'](
        [('O 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0', 3, 'conf-1')])

    assert refs['mol_status'].value == ''
    assert refs['mol_status_fs'].value == ''


def test_the_quick_conversion_writes_plain_coordinates(builder):
    """It answers with a structure, not with a set to choose from, and this tab
    has always taken that as coordinates. Blocks are for the conversions that
    offer conformers."""
    refs, _sent = builder

    refs['offer_structures'](
        [('O 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0', 3, 'quick')], True)

    text = refs['orca_coords'].value
    assert '\n*' not in text
    assert 'quick.xyz' not in text
    assert text.startswith('3\nConverted from SMILES\n')


def test_the_quick_embedding_is_not_offered_as_a_conformer(builder):
    """It rides along at the end of a conformer set as a fallback to step to.
    A block called quick.xyz beside conf-1 and conf-2 says it is one of them.
    """
    refs, _sent = builder

    refs['offer_structures']([
        ('O 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0', 3, 'conf-1'),
        ('O 0 0 0\nH 0.97 0 0\nH -0.25 0.94 0', 3, 'conf-2'),
        ('O 0 0 0\nH 0.98 0 0\nH -0.26 0.95 0', 3, 'quick'),
    ])

    names = [name for name, _xyz in refs['editor_state']['xyz_blocks']]
    assert names == ['conf-1.xyz', 'conf-2.xyz']


def test_to_smiles_reads_its_own_drawing(builder):
    """"Reading the drawing..." never ended in the Builder.

    The button asked the page for "the" drawing frame and "the" channel to
    answer on. With one editor per dashboard that was one of each; with two it
    found the Submit tab's, was told no editor was open in it, and answered
    into that tab's channel -- so the Builder waited for a reply that had gone
    somewhere else.

    Both tabs keep the frame outside their scope container, so the scope
    travels on the element itself. Driven in a browser with two editors on one
    page, each with a drawing of its own: the Submit tab's channel got
    DRAWING-IN-SUBMIT and the Builder's got DRAWING-IN-BUILDER.
    """
    refs, sent = builder
    scope = refs['submit_scope_id']

    assert scope in refs['submit_draw_frame']._dom_classes
    assert scope in refs['submit_draw_sync']._dom_classes

    sent.clear()
    refs['submit_draw_get_btn'].click()

    script = sent[0]
    assert "'.submit-ketcher-frame.'+scope" in script
    assert "'.submit-ketcher-sync.'+scope" in script
    assert scope in script
    assert "document.querySelector('.submit-ketcher-frame')" not in script


def _preview(refs):
    out = []
    for item in refs['orca_mol_output'].outputs or ():
        data = item.get('data') or {}
        out.append('viewer' if data.get('text/html') else (item.get('text') or '').strip())
    return out


def test_a_converted_structure_reaches_the_preview(builder):
    """It did not, and the reason is where the preview is filled from.

    ``with output: display(...)`` only reaches the widget while the kernel is
    running the cell it belongs to. A conversion answers from a thread, through
    the interface loop, and there is no cell there -- so the coordinates landed
    in the box, the preview was asked to draw them, and nothing happened. The
    viewer kept what it had, which after a SMILES was a one-atom model of its
    first letter. Assigning the outputs works wherever the call comes from,
    which is how the Submit tab has always done it.

    Driven in a real dashboard: c1ccccc1, QUICK CONVERT, and within three
    seconds the viewer holds twelve atoms with twelve numbers on them.
    """
    refs, _sent = builder

    refs['orca_coords'].value = 'c1ccccc1'
    assert _preview(refs) and 'SMILES detected' in _preview(refs)[0], (
        'a SMILES is not coordinates, and drawing it as some was the other '
        'half of "it does not show it"')

    refs['offer_structures'](
        [('O 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0', 3, 'quick')], True)

    assert _preview(refs) == ['viewer']
    assert refs['orca_coords'].value.startswith('3\nConverted from SMILES\n')


def test_the_preview_is_filled_by_assignment_not_by_capture():
    assert 'orca_mol_output.outputs = tuple(items)' in BUILDER
    assert 'with orca_mol_output:' not in BUILDER


def test_every_name_in_the_editor_resolves():
    """A name that points nowhere is invisible until the line runs.

    Moving the editor into a part of its own renamed four things the tab used
    to supply, and one call to one of them was left behind. Nothing failed at
    import, nothing failed in a test, and the line only runs when an
    optimisation over several frames finishes -- so "all" raised a NameError
    in both tabs, minutes into a run, and the results were lost.

    This walks build() and asks what it reads that nobody gives it.
    """
    import ast
    import builtins

    source = open(structure_editor.__file__, encoding='utf-8').read()
    tree = ast.parse(source)
    build = next(n for n in tree.body
                 if isinstance(n, ast.FunctionDef) and n.name == 'build')

    bound = {a.arg for a in build.args.args + build.args.kwonlyargs}
    loaded = set()
    for node in ast.walk(build):
        if isinstance(node, ast.Name):
            (bound if isinstance(node.ctx, (ast.Store, ast.Del))
             else loaded).add(node.id)
        elif isinstance(node, ast.FunctionDef):
            bound.add(node.name)
            for arg in node.args.args + node.args.kwonlyargs:
                bound.add(arg.arg)
            if node.args.vararg:
                bound.add(node.args.vararg.arg)
            if node.args.kwarg:
                bound.add(node.args.kwarg.arg)
        elif isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                bound.add((alias.asname or alias.name).split('.')[0])
        elif isinstance(node, ast.ExceptHandler) and node.name:
            bound.add(node.name)
        elif isinstance(node, ast.Lambda):
            for arg in node.args.args:
                bound.add(arg.arg)

    module_level = {n.name for n in tree.body
                    if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                module_level.add((alias.asname or alias.name).split('.')[0])
        if isinstance(node, ast.Assign):
            for target in node.targets:
                for name in ast.walk(target):
                    if isinstance(name, ast.Name):
                        module_level.add(name.id)

    unresolved = sorted(name for name in loaded - bound
                        if name not in dir(builtins) and name not in module_level)
    assert unresolved == [], unresolved


def test_all_reaches_every_frame_the_tab_holds(builder):
    """The Submit tab holds a set of isomers; this one holds named blocks, and
    each of those is a frame too -- so "all" over two blocks means two."""
    refs, _sent = builder
    refs['orca_coords'].value = TWO_BLOCKS

    frames = refs['list_structures']()

    assert [label for _xyz, _n, label in frames] == ['benzene', 'water']
    assert [n for _xyz, n, _label in frames] == [12, 3]


def test_a_new_structure_starts_the_editor_over(builder):
    """A live force field belongs to the molecule its parameters were worked
    out for, and a mode belongs to the structure it was picked on."""
    refs, _sent = builder
    refs['orca_coords'].value = '3\nwater\n' + WATER + '\n'
    for name in ('submit_relax_btn', 'submit_manip_btn', 'submit_draw_btn'):
        refs[name].value = True
    refs['submit_settle_btn'].value = False

    refs['orca_coords'].value = '3\nanother\n' + WATER + '\n'

    for name in ('submit_relax_btn', 'submit_manip_btn', 'submit_draw_btn'):
        assert refs[name].value is False, name
    # To the defaults, not to off: letting go of an atom and leaving the strain
    # of the drag in the structure is the surprising answer, not the useful one.
    assert refs['submit_settle_btn'].value is True


def test_an_edit_keeps_the_header_of_a_plain_xyz(builder):
    """A box that was never written as named blocks holds one plain XYZ,
    header and all. Handing back the bare atom lines took that header away on
    every edit, and after an optimisation the box read as coordinates with no
    count and no comment."""
    refs, _sent = builder
    refs['orca_coords'].value = '3\nwater\n' + WATER + '\n'

    refs['editor_state']['manip_inflight'] = True
    refs['editor_coords'].value = (
        '3\nEdited in DELFIN viewer\nO 0.2 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n')

    lines = refs['orca_coords'].value.split('\n')
    assert lines[0].strip() == '3', refs['orca_coords'].value[:60]
    assert lines[2].startswith('O 0.2')


def test_the_editor_starts_the_same_way_in_both_tabs(builder):
    """Numbers off until they are asked for, and the size field with them.

    The Builder used to bring them up by itself, from the days when its
    numbering was its own; the editor's default is off, and one default is
    better than two.
    """
    refs, _sent = builder
    refs['orca_coords'].value = '3\nwater\n' + WATER + '\n'

    assert refs['submit_labels_btn'].value is False
    assert refs['submit_label_size'].layout.display == 'none'


def test_the_force_field_notes_stay_in_the_small_view(builder):
    """They are several lines of prose about what had to be approximated --
    a forced hybridisation, a metal the parameters do not really cover. In
    fullscreen they take that space off the structure they describe, and the
    Submit tab has never put them there either."""
    refs, _sent = builder
    classes = refs['submit_ff_notes']._dom_classes

    assert 'delfin-structure-fs-member' not in classes
    assert not any(c.startswith('submit-fs-member') for c in classes)


def test_drawing_on_and_asking_again_gives_the_new_structure(builder):
    """Ketcher stays open while a structure grows. Each TO SMILES reads what is
    there now, and the answer carries a serial in front so that asking twice
    about the same drawing reads as two answers rather than as one that never
    came.

    Measured in both tabs: a four-carbon chain gives CCCC, drawing two more on
    and asking again gives CCCCCC, and asking a third time with nothing changed
    gives CCCCCC.
    """
    refs, _sent = builder

    def molfile(n):
        head = '\n  Ketcher\n\n%3d%3d  0  0  0  0            999 V2000\n' % (n, n - 1)
        atoms = '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n' * n
        bonds = ''.join('%3d%3d  1  0  0  0  0\n' % (i, i + 1) for i in range(1, n))
        return head + atoms + bonds + 'M  END\n'

    refs['submit_draw_sync'].value = '1\n' + molfile(4)
    assert refs['orca_coords'].value == 'CCCC'

    refs['submit_draw_sync'].value = '2\n' + molfile(6)
    assert refs['orca_coords'].value == 'CCCCCC'

    refs['submit_draw_sync'].value = '3\n' + molfile(6)
    assert refs['orca_coords'].value == 'CCCCCC'


def test_settle_is_on_to_begin_with_in_both_tabs():
    from delfin.dashboard import structure_editor as editor

    source = open(editor.__file__, encoding='utf-8').read()
    made = source.split('submit_settle_btn = widgets.ToggleButton')[1].split(')\n')[0]
    assert 'value=True' in made
    # ...and a structure the editor has not seen puts it back on rather than
    # taking it away.
    defaults = source.split('_CONTROL_DEFAULTS = (')[1].split('\n    )')[0]
    assert 'submit_settle_btn, True' in defaults
    assert 'submit_relax_btn, False' in defaults


def test_converting_again_builds_what_is_in_the_box_now(builder):
    """Draw more in Ketcher, hand it back, convert -- and the structure that
    came out was the one before it.

    The quick conversion remembers the SMILES it last built, so that pressing
    it again rolls another embedding of the same molecule: by then the box
    holds coordinates and has nothing to offer. It took the remembered one
    even when the box held a different SMILES. What is in the box wins
    whenever it is a SMILES.

    Measured in both tabs: a four-carbon chain converts to fourteen atoms,
    drawing four more on and converting again gives twenty-six.
    """
    refs, _sent = builder
    state = refs['editor_state']

    state['converted_xyz_cache'] = {'smiles': 'CCCC', 'xyz': 'stale'}
    refs['orca_coords'].value = 'CCCCCCCC'

    # what the conversion would read, without running it
    assert refs['read_input']() == 'CCCCCCCC'
    source = open(structure_editor.__file__, encoding='utf-8').read()
    body = source.split('def _start_smiles_conversion')[1].split('\n    def ')[0]
    assert "if typed and clean_input_data(typed)[1] == 'smiles':" in body
    assert 'raw_input = typed' in body
    # ...and the remembered one is still there for a box holding coordinates.
    assert 'raw_input = (cached_smiles or typed' in body


def test_two_molecules_in_one_drawing_are_a_smiles():
    """Draw two things in Ketcher and it hands back "CCO.CCO".

    The dot is how SMILES separates one molecule from another, and it was not
    among the characters that made a string look like one -- that string has no
    bracket, no aromatic letter and no ring number, so nothing recognised it.
    Two consequences, and the second is the one that was reported: such a
    drawing could not be converted at all, and the quick conversion fell back
    on the SMILES it had built before, so pressing convert produced the
    structure from the drawing before last.

    Coordinates cannot be caught by the dot: a line of them is turned away
    first, by the pattern for an element followed by three numbers.
    """
    from delfin.dashboard.input_processing import clean_input_data

    for text in ('CCO.CCO', '[Na+].[Cl-]', 'c1ccccc1.c1ccccc1', 'CCO.O'):
        assert clean_input_data(text)[1] == 'smiles', text

    for text in ('O 0.000 0.000 0.000', 'C  1.234  5.678  9.012',
                 '3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0',
                 # And the decimal point is not what makes a line coordinates:
                 # "C 0 0 0" is an ordinary way to write an atom at the origin,
                 # and it was read as a SMILES because it has digits in it.
                 'C 0 0 0\nC 1.53 0 0', 'H 0 0 0', 'Cl -1 2 3'):
        assert clean_input_data(text)[1] == 'xyz', text


ETHANE = ("C 0.000 0.000 0.000\nC 1.530 0.000 0.000\n"
          "H -0.360 1.020 0.000\nH -0.360 -0.510 0.880\n"
          "H -0.360 -0.510 -0.880\nH 1.890 1.020 0.000\n"
          "H 1.890 -0.510 0.880\nH 1.890 -0.510 -0.880")


def test_stepping_away_and_back_keeps_what_the_editor_knew(builder):
    """Each structure keeps its own bonding, held values and history.

    Bonding is perceived once, from the structure as it arrived, so that
    dragging an atom away from its neighbour does not decide the bond was never
    there. Step to another block and back without a memory and the coordinates
    are read afresh -- with the atom where it was dragged to. Perceived fresh,
    the pulled ethane has six bonds where it had seven.
    """
    refs, _sent = builder
    state = refs['editor_state']
    water = "O 0.000 0.000 0.000\nH 0.960 0.000 0.000\nH -0.240 0.930 0.000"
    refs['orca_coords'].value = (
        'eth.xyz;\n8\n\n' + ETHANE + '\n*\n\nwat.xyz;\n3\n\n' + water + '\n*')

    refs['submit_relax_btn'].value = True          # makes it perceive
    before = len(state['perceived'].bonds)
    assert before == 7

    state['manip_inflight'] = True
    refs['editor_coords'].value = '8\nEdited in DELFIN viewer\n' + ETHANE.replace(
        'H 1.890 -0.510 -0.880', 'H 4.500 -2.500 -3.000')
    assert len(state['perceived'].bonds) == before, 'the drag kept the bonding'

    refs['orca_mol_next_btn'].click()
    refs['orca_mol_prev_btn'].click()

    assert len(state['perceived'].bonds) == before


def test_a_fresh_perception_of_the_pulled_structure_would_have_lost_it():
    """Which is what made the round trip lose the bond, and what the memory is
    there to prevent."""
    from delfin.dashboard.molecule_forcefield import perceive_molecule

    pulled = ETHANE.replace('H 1.890 -0.510 -0.880', 'H 4.500 -2.500 -3.000')
    assert len(perceive_molecule('8\npulled\n' + pulled).bonds) == 6
    assert len(perceive_molecule('8\nwhole\n' + ETHANE).bonds) == 7


def test_the_memory_names_what_belongs_to_a_structure():
    keys = structure_editor.STRUCTURE_MEMORY_KEYS

    for key in ('perceived', 'bond_edits', 'hand_bonds', 'hyb_overrides',
                'constraints', 'poly_applied', 'history', 'pristine_coords'):
        assert key in keys, key
    # Not the switches: those belong to the editor, not to one structure, and
    # a new structure puts them back to their defaults.
    for key in ('manip_bootstrap_done', 'smiles_task_id', 'isomers'):
        assert key not in keys, key


def _five_blocks():
    water = "O 0.000 0.000 0.000\nH 0.960 0.000 0.000\nH -0.240 0.930 0.000"
    return "\n\n".join('conf-%d.xyz;\n3\nconf-%d\n%s\n*' % (i, i, water)
                        for i in range(1, 6))


def test_the_frame_on_screen_is_the_one_that_is_worked_on(builder):
    """Five blocks, standing on the third: the editor holds that one, an edit
    lands in that one, and the other four are left alone.

    Driven through the real tab with GFN2-xTB: Optimize on 3/5 answers
    "Optimised 1 of 1 frame(s) ... E = -5.070544 Eh", writes the result into
    conf-3.xyz, leaves conf-1.xyz as it was, and stays on 3/5.
    """
    refs, _sent = builder
    state = refs['editor_state']
    refs['orca_coords'].value = _five_blocks()

    refs['orca_mol_next_btn'].click()
    refs['orca_mol_next_btn'].click()
    assert state['xyz_view_idx'] == 2
    assert refs['submit_manip_toolbar'].layout.display == 'flex'
    assert refs['editor_coords'].value.split('\n')[0] == '3'

    state['manip_inflight'] = True
    refs['editor_coords'].value = (
        '3\nEdited in DELFIN viewer\nO 0.500 0.000 0.000\n'
        'H 0.960 0.000 0.000\nH -0.240 0.930 0.000\n')

    touched = [name for name, xyz in state['xyz_blocks'] if 'O 0.500' in xyz]
    assert touched == ['conf-3.xyz']
    assert [name for name, _xyz in state['xyz_blocks']] == [
        'conf-1.xyz', 'conf-2.xyz', 'conf-3.xyz', 'conf-4.xyz', 'conf-5.xyz']
    assert state['xyz_view_idx'] == 2, 'the view did not wander off the frame'


def test_undo_and_the_rest_belong_to_the_frame_they_were_made_on(builder):
    """Each block keeps its own history, held values and bond edits -- they are
    in the memory that is put aside when stepping away."""
    keys = structure_editor.STRUCTURE_MEMORY_KEYS

    for key in ('history', 'structure_undo', 'bond_edits', 'constraints'):
        assert key in keys, key

    refs, _sent = builder
    refs['orca_coords'].value = _five_blocks()
    state = refs['editor_state']
    refs['orca_mol_next_btn'].click()
    state['constraints'] = [{'kind': 'distance', 'atoms': [0, 1],
                             'value': 1.1, 'mode': 'fix'}]

    refs['orca_mol_next_btn'].click()
    assert state.get('constraints') == [], 'a held value followed to another frame'

    refs['orca_mol_prev_btn'].click()
    assert len(state.get('constraints') or []) == 1, 'and did not come back'


def test_all_optimises_every_block(builder):
    """Driven through the real tab on five conformers with UFF: "Optimised 5 of
    5 frame(s) with UFF", all five blocks changed, all five names kept."""
    refs, _sent = builder
    refs['orca_coords'].value = _five_blocks()
    state = refs['editor_state']
    before = [xyz.split('\n')[2] for _name, xyz in state['xyz_blocks']]

    refs['submit_optimize_all_btn'].value = True
    for _ in range(200):
        time.sleep(0.5)
        if not refs['submit_optimize_all_btn'].value:
            break

    assert 'Optimised 5 of 5 frame(s)' in refs['mol_status'].value
    after = [xyz.split('\n')[2] for _name, xyz in state['xyz_blocks']]
    assert sum(1 for a, b in zip(before, after) if a != b) == 5
    assert [name for name, _xyz in state['xyz_blocks']] == [
        'conf-%d.xyz' % i for i in range(1, 6)]
