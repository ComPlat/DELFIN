"""Leaving the TS -> GOAT panel gives the file view back whole.

Input: the Calculations Browser with the panel open on a TS ``.out``,
then a selection that ends it (none, or a click on another file).
Output: the panel is hidden and the content label and toolbar are shown
again as the view mode requires.

Opening the panel hides the content area, its label and its toolbar.
The two ways of leaving it other than its Close button restored only the
content area, so the label and the search/jump toolbar stayed hidden
until something else redrew the view. Both now call calc_update_view.
"""

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext
from tests.test_ts_goat_dialog_regression import _FakeBackend, _write_ts_fixture


@pytest.fixture
def refs(tmp_path):
    calc = tmp_path / 'calc'
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    job = calc / 'job'
    job.mkdir()
    _write_ts_fixture(job)
    (job / 'other.out').write_text(
        'FINAL SINGLE POINT ENERGY -1.0\n', encoding='utf-8')
    ctx = DashboardContext(
        calc_dir=calc,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    ctx.run_js = lambda script: None
    ctx.backend = _FakeBackend()
    _widget, refs = browser.create_tab(ctx)
    refs['calc_path_input'].value = str(job)
    return refs


def _label(refs, name):
    for option in refs['calc_file_list'].options:
        if str(option).replace('\xa0', ' ').strip().endswith(name):
            return option
    raise AssertionError(f'{name} not listed')


def _click(refs, name):
    label = _label(refs, name)
    refs['calc_open_input'].value = label
    refs['calc_file_list'].value = (label,)


def _open_panel(refs):
    _click(refs, 'ts_run.out')
    dropdown = refs['calc_options_dropdown']
    assert 'TS → GOAT' in dropdown.options
    dropdown.value = 'TS → GOAT'
    assert refs['calc_ts_goat_panel'].layout.display == 'block'
    assert refs['calc_ts_goat_panel'].children
    assert refs['calc_content_label'].layout.display == 'none'
    assert refs['calc_content_toolbar'].layout.display == 'none'


def _assert_file_view(refs):
    assert refs['calc_ts_goat_panel'].layout.display == 'none'
    assert refs['calc_content_area'].layout.display == 'block'
    assert refs['calc_content_label'].layout.display == 'block'
    assert refs['calc_content_toolbar'].layout.display == 'flex'


def test_clearing_the_selection_restores_label_and_toolbar(refs):
    _open_panel(refs)
    refs['calc_file_list'].value = ()
    _assert_file_view(refs)


def test_opening_another_file_restores_label_and_toolbar(refs):
    _open_panel(refs)
    _click(refs, 'other.out')
    _assert_file_view(refs)
