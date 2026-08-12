"""The options dropdown emptied itself on every click.

Clicking a file is reported twice. The browser's click bridge arrives first
and opens the file; the selection change follows and stands down, because the
file is already open. That order is deliberate (see
``test_one_click_opens_a_file_once``) -- but it left the options dropdown
reading the wrong list.

``calc_update_options_dropdown`` is called while the file is being opened, and
it asked ``calc_file_list.value`` which file that was. At that moment the list
still carries the *previous* selection: empty on the first click of a session,
the file before it afterwards. So for CONTROL.txt it fell to the empty branch
and hid the dropdown -- Recalc and Smart Recalc were gone, and so were
Override on OCCUPIER.txt, Recalc on .inp and Print Mode / MO Plot on .out.
Only the Edit button was left, which is why the toolbar looked reduced to
editing.

The label now travels with the open instead of being looked up.
"""

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


@pytest.fixture
def refs(tmp_path):
    calc = tmp_path / 'calc'
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    job = calc / 'job'
    job.mkdir()
    (job / 'CONTROL.txt').write_text('CHARGE = 0\n', encoding='utf-8')
    (job / 'stage_OCCUPIER.txt').write_text('OPT 1\n', encoding='utf-8')
    (job / 'stage.inp').write_text(
        '! PBE0 def2-SVP\n* xyz 0 1\nH 0.0 0.0 0.0\n*\n', encoding='utf-8')
    (job / 'stage.out').write_text(
        'FINAL SINGLE POINT ENERGY -1.0\n', encoding='utf-8')
    # Past CALC_TEXT_FULL_READ_BYTES, so this one is read in chunks -- the
    # path a real GOAT output takes.
    (job / 'huge.out').write_text(
        'FINAL SINGLE POINT ENERGY -1.0\n' * 120_000, encoding='utf-8')
    (job / 'notes.md').write_text('nothing to calculate here\n', encoding='utf-8')

    ctx = DashboardContext(
        calc_dir=calc,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    refs['calc_path_input'].value = str(job)
    return refs


def _label_for(refs, name):
    """The label the list draws for *name* -- icon, padding and all."""
    for option in refs['calc_file_list'].options:
        if str(option).replace('\xa0', ' ').strip().endswith(name):
            return option
    raise AssertionError(f'{name} is not in {list(refs["calc_file_list"].options)}')


def _click(refs, name):
    """One press, as the browser reports it: bridge first, selection after."""
    label = _label_for(refs, name)
    refs['calc_open_input'].value = label       # the click bridge
    refs['calc_file_list'].value = (label,)     # the selection it causes
    return label


@pytest.mark.parametrize('name, expected', [
    ('CONTROL.txt', ['Recalc', 'Smart Recalc']),
    ('stage_OCCUPIER.txt', ['Override']),
    ('stage.inp', ['Recalc']),
    ('stage.out', ['Print Mode', 'MO Plot', 'Print NMR']),
    ('huge.out', ['Print Mode', 'MO Plot', 'Print NMR']),
])
def test_a_clicked_file_offers_its_actions(refs, name, expected):
    _click(refs, name)

    dropdown = refs['calc_options_dropdown']
    assert list(dropdown.options)[1:] == expected
    assert dropdown.layout.display == 'block'


def test_the_dropdown_does_not_lag_one_file_behind(refs):
    """The shape of the failure: the previous selection answered for this one."""
    _click(refs, 'CONTROL.txt')
    _click(refs, 'stage.out')

    assert 'Smart Recalc' not in refs['calc_options_dropdown'].options
    assert 'Print Mode' in refs['calc_options_dropdown'].options


def test_a_file_with_no_actions_hides_the_dropdown(refs):
    _click(refs, 'CONTROL.txt')
    _click(refs, 'notes.md')

    assert list(refs['calc_options_dropdown'].options) == ['(Options)']
    assert refs['calc_options_dropdown'].layout.display == 'none'


def test_leaving_the_folder_takes_the_actions_with_it(refs):
    """Nothing is on screen after a listing, so nothing may be offered."""
    _click(refs, 'CONTROL.txt')
    assert 'Smart Recalc' in refs['calc_options_dropdown'].options

    refs['calc_list_directory']()

    assert list(refs['calc_options_dropdown'].options) == ['(Options)']
    assert refs['calc_options_dropdown'].layout.display == 'none'
