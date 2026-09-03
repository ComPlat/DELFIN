"""Two questions a drag can be asked, and the switch that chooses.

Off -- the default, and what the editor has always done -- the perception is
asked what has *changed* since the last answer.  On, it is asked what the
direction the hand is going *carries*, which is the question the first answer
of a drag is already asked.

"What changed" answers itself badly.  The coordinate being held is, by
construction, the one that has changed least, because holding it is what
stopped it changing: it falls in the ranking, something else is held, the
first one moves again, and the question flips answer by answer.  That is the
shaking users report, and replaying one real session's own frames showed it
plainly -- sixteen answers of one grab named eleven different coordinates, and
flipped between a distance and a torsion, which are two different statements
about the same gesture.

Measured over that session's twenty-five grabs and 567 answers, with both
floors and the geometry guard in place:

                         changes   over the ceiling   worst answer
    what changed (off)      128           163            +98.2
    what the hand carries    30           116           +296.2

Four times steadier and slightly cheaper on the whole, three times worse at
its worst.  That is the trade, and it is why this is a switch: the steady rule
holds what the hand is really driving -- on the grab measured, the C-H bond of
the dragged hydrogen, which is what a chemist would name -- and holding
something real costs something.  The old rule stays cheap by never holding
anything meaningful.
"""
from __future__ import annotations

import pathlib

_EDITOR = (pathlib.Path(__file__).resolve().parents[1]
           / 'delfin' / 'dashboard' / 'structure_editor.py').read_text()


def test_the_default_is_what_the_editor_has_always_done():
    block = _EDITOR.split('submit_steady_hand_btn = widgets.Checkbox(')[1]
    block = block.split(')')[0]
    assert 'value=False' in block, 'the new rule must not arrive switched on'


def test_the_follow_asks_the_question_the_switch_chose():
    body = _EDITOR.split('def on_submit_relax_toggle(')[0]
    body = body.split('steady = bool(submit_steady_hand_btn.value)')[1][:600]
    assert 'opening=steady,' in body


def test_it_stands_beside_the_pull_and_only_where_there_are_two_questions():
    """Under a placing hand the answer is laid back onto the cursor and no
    coordinate is read at all, so there is nothing to choose between."""
    body = _EDITOR.split('submit_view_body = widgets.VBox(')[1][:400]
    assert 'submit_pull_slider' in body and 'submit_steady_hand_btn' in body
    assert body.index('submit_pull_slider') < body.index(
        'submit_steady_hand_btn')
    shown = _EDITOR.split('submit_steady_hand_btn.layout.display = (')[1]
    shown = shown.split("else 'none')")[0]
    assert 'pulling' in shown
    assert '_server_method()' in shown


def test_the_measurement_is_written_down_beside_the_switch():
    """A switch between two behaviours is a decision, and a decision needs
    its numbers where the next person will find them."""
    note = _EDITOR.split('submit_steady_hand_btn = widgets.Checkbox(')[0]
    note = note[-2600:]
    for marker in ('567 answers', '128', '163', '+98.2', '116',
                   '+296.2', 'twenty-five grabs'):
        assert marker in note, marker


def test_the_editor_hands_the_switch_out_like_every_other_control():
    """So a host can place it and a test can drive it."""
    import tempfile

    import pytest
    pytest.importorskip('ipywidgets')
    import ipywidgets as widgets

    from delfin.dashboard import structure_editor
    from delfin.dashboard.context import DashboardContext

    room = pathlib.Path(tempfile.mkdtemp())
    for name in ('calc', 'archive', 'office'):
        (room / name).mkdir()
    ctx = DashboardContext(calc_dir=room / 'calc',
                           archive_dir=room / 'archive',
                           office_dir=room / 'office')
    ctx.run_js = lambda _s: None
    ctx.add_init_js = lambda _s: None
    part = structure_editor.build(
        ctx, state={}, viewer_height=560,
        coords_widget=widgets.Textarea(
            value='3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n'),
        schedule_ui_update=lambda f, *a, **k: f(*a, **k),
        update_view=lambda *a, **k: None,
        get_smiles_charge=lambda *a, **k: None)
    assert part.submit_steady_hand_btn.value is False
    assert part.submit_steady_hand_btn.description == 'steady gesture'
