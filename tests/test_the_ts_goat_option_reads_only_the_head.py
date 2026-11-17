"""The 'TS → GOAT' option is decided from the head of a .out only.

Input: a click on a ``.out`` in the Calculations Browser. Output: the
options list, with 'TS → GOAT' when the input echo names OptTS.

The option list is rebuilt on every .out selection. Deciding it with
``read_ts_status`` read the whole output and scanned every line, which on
a large output blocks the tab for each click; the decision now reads at
most ``ts_goat._OPTTS_HEAD_BYTES``. The full scan stays in the panel,
which runs only when the option is chosen.
"""

import pytest

from delfin import ts_goat
from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


@pytest.fixture
def refs(tmp_path, monkeypatch):
    calc = tmp_path / 'calc'
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    job = calc / 'job'
    job.mkdir()
    (job / 'ts.out').write_text(
        '|  1> ! r2SCAN-3c optts Freq\n'
        'THE OPTIMIZATION HAS CONVERGED\n', encoding='utf-8')
    (job / 'plain.out').write_text(
        '|  1> ! r2SCAN-3c Opt\nFINAL SINGLE POINT ENERGY -1.0\n',
        encoding='utf-8')
    # The keyword sits past the head: not read, so not offered.
    filler = 'x' * 1023 + '\n'
    (job / 'late.out').write_text(
        filler * (ts_goat._OPTTS_HEAD_BYTES // 1024 + 1) + '! OptTS\n',
        encoding='utf-8')

    def _full_read(*_a, **_k):
        raise AssertionError('read_ts_status ran on a selection')

    monkeypatch.setattr(ts_goat, 'read_ts_status', _full_read)
    ctx = DashboardContext(
        calc_dir=calc,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    refs['calc_path_input'].value = str(job)
    return refs


def _click(refs, name):
    for option in refs['calc_file_list'].options:
        if str(option).replace('\xa0', ' ').strip().endswith(name):
            refs['calc_open_input'].value = option
            refs['calc_file_list'].value = (option,)
            return
    raise AssertionError(f'{name} not listed')


@pytest.mark.parametrize('name, offered', [
    ('ts.out', True),        # lower-case keyword: ORCA is case-insensitive
    ('plain.out', False),
    ('late.out', False),
])
def test_the_option_follows_the_head_of_the_output(refs, name, offered):
    _click(refs, name)
    assert ('TS → GOAT' in refs['calc_options_dropdown'].options) is offered


def test_the_head_check_reads_a_bounded_prefix(tmp_path):
    out = tmp_path / 'a.out'
    out.write_text('y' * 200 + '\n! OptTS\n', encoding='utf-8')
    assert ts_goat.is_optts_output(out)
    assert not ts_goat.is_optts_output(out, max_bytes=100)
    assert not ts_goat.is_optts_output(tmp_path / 'missing.out')
