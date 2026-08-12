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


def _shown(monkeypatch):
    """Whatever the tab hands to display(), as text."""
    drawn = []
    monkeypatch.setattr(
        tab_orca_builder, 'display',
        lambda *a, **k: drawn.append(getattr(a[0], 'data', str(a[0])) if a else ''))
    return drawn


SWAPPED = "O 0.000 0.000 0.000\nH -0.757 0.586 0.000\nH 0.757 0.586 0.000"


@pytest.fixture
def compared(builder, monkeypatch, capsys):
    refs, sent = builder
    drawn = _shown(monkeypatch)
    refs['orca_coords'].value = (f"a.xyz;\n3\n\n{WATER}\n*\n\n"
                                 f"b.xyz;\n3\n\n{SWAPPED}\n*")
    drawn.clear()
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


def test_the_fullscreen_button_is_there_for_one_structure_too(builder):
    """It sits in the row with the block stepper, and that row was hidden
    whenever the coordinates had not been written as named blocks.

    So pasting a plain XYZ -- which is what the editor itself writes back, and
    what anybody pasting coordinates starts with -- left no way at all to
    enlarge the viewer. With two blocks it worked, which is exactly how the
    report came in.

    Driven in a real dashboard on a bare twelve-atom benzene: the button is
    visible, and a click takes the viewer from 721x560 to 1484x766.
    """
    refs, _sent = builder
    row = refs['orca_mol_nav_row']

    refs['orca_coords'].value = ''
    assert row.layout.display == 'none', 'nothing to enlarge, nothing to show'

    refs['orca_coords'].value = f'3\nwater\n{WATER}\n'
    assert row.layout.display == ''
    assert refs['orca_mol_next_btn'].layout.display == 'none', (
        'one structure has nothing to step to')
    assert refs['orca_mol_fullscreen_btn'].layout.display != 'none'

    refs['orca_coords'].value = TWO_BLOCKS
    assert row.layout.display == ''
    assert refs['orca_mol_next_btn'].layout.display == ''
