"""What one press of Optimise does to a value the user is holding.

The editor collects the held values when Optimise is pressed and hands them
to whichever engine is running.  It handed them to xtb and to MOPAC and not
to the force fields that run in the browser, which is the third branch of the
same function -- so under UFF or MMFF94 the whole list went quietly, and the
status line said nothing about having dropped it.

Measured on ethane through the editor's own channels, the bonding pinned the
way switching Dynamik Opt on pins it, C0-H2 pulled out to 1.700 A and held
exactly::

                        before        after
    one press of Optimise   1.1104 A      1.7000 A
    what the line said      nothing       which of the two it did

A force field here can fix an atom and cannot restrain a value: RDKit's UFF
takes ``AddFixedPoint``, Open Babel a list of fixed atoms, and neither has a
force constant to negotiate with.  So a fix is met by holding the two atoms
that name it still -- which is more than was asked, and is why the line has
to say so -- and a pull cannot be expressed at all and is reported as not
honoured rather than quietly frozen.  That is the reading MOPAC's flags
already get, from the same function, so a value one engine drops is dropped
by the other for the same reason.

The gesture is driven the way the page drives it: a release arrives on
``submit_manip_sync`` as an XYZ whose comment line reads ``DELFIN drag-end``,
the pair is picked on ``submit_pick_sync``, and Optimise is a toggle.
"""

from __future__ import annotations

import math

import pytest

from delfin.dashboard.context import DashboardContext

pytest.importorskip('rdkit')

_ETHANE = """8
ethane
C -0.7560  0.0000  0.0000
C  0.7560  0.0000  0.0000
H -1.1404  0.6586  0.7845
H -1.1404  0.3501 -0.9626
H -1.1404 -1.0087  0.1781
H  1.1404 -0.6586 -0.7845
H  1.1404 -0.3501  0.9626
H  1.1404  1.0087 -0.1781
"""


def _rows(text):
    out = []
    for line in text.splitlines()[2:]:
        bits = line.split()
        if len(bits) >= 4:
            out.append((bits[0], float(bits[1]), float(bits[2]), float(bits[3])))
    return out


def _span(text, i, j):
    rows = _rows(text)
    return math.dist(rows[i][1:], rows[j][1:])


def _pulled_to(xyz, i, j, value):
    """*j* moved along the i-j direction until the two are *value* apart.

    The message the page sends when the mouse comes up: the geometry the drag
    made, with the comment line that says a drag ended.
    """
    rows = _rows(xyz)
    a, b = rows[i], rows[j]
    v = [b[1] - a[1], b[2] - a[2], b[3] - a[3]]
    length = math.sqrt(sum(c * c for c in v))
    rows[j] = (b[0],) + tuple(a[1 + k] + v[k] / length * value for k in range(3))
    body = '\n'.join(f'{s}  {x:.4f}  {y:.4f}  {z:.4f}' for s, x, y, z in rows)
    return f'{len(rows)}\nDELFIN drag-end\n{body}\n'


@pytest.fixture
def editor(tmp_path):
    """One structure editor over a coordinate box of its own, holding ethane."""
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor

    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(
        calc_dir=tmp_path / 'calc',
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    ctx.run_js = lambda _script: None
    state = {}
    box = widgets.Textarea(value=_ETHANE)
    part = structure_editor.build(
        ctx,
        state=state,
        coords_widget=box,
        viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None,
    )
    return part, state, box


def _drag_and_hold(part, state, mode, method='uff', value=1.700):
    """Pull H2 off C0 to *value*, hold it there, and press Optimise once."""
    part.submit_ff_dd.value = method
    # Dynamik Opt on: this is what pins the bonding, so the press relaxes the
    # molecule the user has been dragging rather than re-perceiving a
    # stretched C-H as no bond at all.  The switch is detached from its own
    # observer, which refuses to stay on where the page cannot answer.
    part.submit_relax_btn.unobserve_all()
    part.submit_relax_btn.value = True
    part._enable_live_forcefield()

    part.submit_manip_sync.value = _pulled_to(_ETHANE, 0, 2, value)
    if mode is not None:
        part.submit_pick_sync.value = '0,2'
        part.submit_hold_mode.value = mode
        part.submit_internal_value.value = value
        part.on_submit_hold()

    part.submit_optimize_btn.unobserve_all()
    part.submit_optimize_btn.value = True
    part.on_submit_optimize(None)
    thread = state.get('optimize_thread')
    if thread is not None:
        thread.join(120)
    return part.mol_status.value


def test_a_value_held_exactly_survives_a_press_of_optimise(editor):
    """1.700 A held exactly came back from the press at 1.1104 A -- the whole
    held list was collected and then not handed to the one engine that was
    running.  The atoms that name it are fixed instead, and it comes back at
    1.7000."""
    part, state, box = editor

    said = _drag_and_hold(part, state, 'fix')

    assert _span(box.value, 0, 2) == pytest.approx(1.700, abs=1e-3), (
        'the press walked off a value that was being held exactly')
    # And the two other carbons-and-hydrogens did move: fixing the pair must
    # not turn the press into a no-op.
    assert any(
        math.dist(_rows(box.value)[i][1:], _rows(_ETHANE)[i][1:]) > 1e-3
        for i in range(3, 8)
    ), 'nothing else relaxed, so the minimisation did not really run'
    assert 'Optimised with UFF' in said
    # Which of the two it did, because fixing the atoms is more than holding
    # the value between them.
    assert 'fixing the 2 atom(s) they name' in said
    assert 'not the value between them' in said


def test_a_pull_is_reported_rather_than_quietly_frozen(editor):
    """A pull asks the chemistry to argue with the number.  There is nothing
    here to argue with -- an atom is fixed or free -- so freezing it would
    hold as exact a value the user asked to have negotiated, and dropping it
    without a word is what the branch did to the whole list."""
    part, state, box = editor

    said = _drag_and_hold(part, state, 'pull')

    # Still relaxed away, as it always was: 1.700 A back to the C-H length.
    assert _span(box.value, 0, 2) < 1.3
    assert 'pull(s) not honoured' in said
    assert 'optimise under a GFN method' in said


def test_nothing_held_leaves_the_press_exactly_as_it_was(editor):
    """The structure the user dragged is minimised whole when no value is
    being held: the atoms are read off the held list, so an empty list fixes
    nothing and the line has nothing extra to say."""
    part, state, box = editor

    said = _drag_and_hold(part, state, None)

    assert _span(box.value, 0, 2) < 1.3, 'the stretched C-H must relax'
    assert said.count('held value') == 0
    assert said.count('not honoured') == 0


def test_a_water_whose_every_atom_is_held_is_told_so(tmp_path):
    """Fixing atoms runs out of molecule sooner than restraining values does.

    An angle held on a water names all three of its atoms, so a press that
    honours it by fixing them has nothing left to move.  Measured: the largest
    displacement is 0.000000 A.  That is the right answer to what was asked
    and a bad thing to call "Optimised with UFF." and leave there -- it is the
    same silence one step further along -- so the line says it.
    """
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor

    water = ('3\nwater\nO  0.0000  0.0000  0.1173\n'
             'H  0.0000  0.9000 -0.4692\nH  0.0000 -0.9000 -0.4692\n')
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    ctx = DashboardContext(calc_dir=tmp_path / 'calc',
                           archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda _script: None
    state = {}
    box = widgets.Textarea(value=water)
    part = structure_editor.build(
        ctx, state=state, coords_widget=box, viewer_height=560,
        schedule_ui_update=lambda func, *a, **k: func(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    part.submit_ff_dd.value = 'uff'
    part.submit_relax_btn.unobserve_all()
    part.submit_relax_btn.value = True
    part._enable_live_forcefield()

    part.submit_pick_sync.value = '1,0,2'
    part.submit_hold_mode.value = 'fix'
    part.submit_internal_value.value = 100.0
    part.on_submit_hold()

    before = box.value
    part.submit_optimize_btn.unobserve_all()
    part.submit_optimize_btn.value = True
    part.on_submit_optimize(None)
    thread = state.get('optimize_thread')
    if thread is not None:
        thread.join(120)

    moved = max(math.dist(a[1:], b[1:])
                for a, b in zip(_rows(before), _rows(box.value)))
    assert moved == pytest.approx(0.0, abs=1e-6)
    said = part.mol_status.value
    assert 'every atom in the structure' in said, said
    assert 'nothing left to relax' in said, said


def test_the_line_reads_a_held_value_in_the_terms_of_the_engine_that_ran():
    """Three engines, three different things done with one held value, and the
    result is only honest if it is read out in the terms of whichever ran.  A
    force field fixes the atoms and has no force constant, so reading its
    result with xtb's sentence would claim one that no run has."""
    from delfin.dashboard import mopac_optimize as mopac
    from delfin.dashboard.structure_editor import fixed_atoms_note

    fix = [{'kind': 'distance', 'atoms': [0, 2], 'value': 1.7, 'mode': 'fix'}]
    pull = [{'kind': 'distance', 'atoms': [0, 2], 'value': 1.7, 'mode': 'pull'}]

    said = fixed_atoms_note(mopac.freeze_flags(fix, atoms=8))
    assert '1 held value(s) kept by fixing the 2 atom(s)' in said
    assert 'the force field fixes atoms' in said
    # Not MOPAC's sentence, though it is MOPAC's reading: the engine that ran
    # is the one that has to be named.
    assert 'MOPAC' not in said

    assert 'not honoured' in fixed_atoms_note(mopac.freeze_flags(pull, atoms=8))

    # An entry naming an atom this structure does not have is dropped, and
    # said to have been -- the same rule the server engines apply.
    away = mopac.freeze_flags(fix, atoms=2)
    assert away['frozen'] == set()
    assert 'name atoms this structure does not have' in fixed_atoms_note(away)

    # Nothing held, nothing said.
    assert fixed_atoms_note(mopac.freeze_flags([], atoms=8)) == ''
