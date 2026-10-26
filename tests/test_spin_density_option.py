"""The Spin Density option only makes sense for open-shell jobs.

Selecting 'Spin Density' shows a panel with a plot button. Pressing it on a
closed-shell (.out) job must refuse with a clear message, and never shell out
to orca_plot. The happy path is delegated to orca_plot's stdin workflow and is
covered by the MO Plot plumbing it reuses; this file pins the guard rails.
"""

import types

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


@pytest.fixture
def refs(tmp_path, monkeypatch):
    calc = tmp_path / 'calc'
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    job = calc / 'job'
    job.mkdir()
    # Closed-shell output: HFTyp ... RKS.
    (job / 'stage.out').write_text(
        'Hartree-Fock type          RKS\n'
        'FINAL SINGLE POINT ENERGY -1.0\n', encoding='utf-8')
    # Open-shell output: HFTyp ... UKS.
    (job / 'open.out').write_text(
        'Hartree-Fock type          UKS\n'
        'FINAL SINGLE POINT ENERGY -1.0\n', encoding='utf-8')
    # Matching wavefunction file so _calc_mo_find_gbw finds it.
    (job / 'open.gbw').write_bytes(b'\x00' * 16)

    ctx = DashboardContext(
        calc_dir=calc,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    refs['calc_path_input'].value = str(job)

    # Record any orca_plot invocation; the real binary must never run.
    calls = []
    fake_proc = types.SimpleNamespace(returncode=1, stdout=b'', stderr=b'')
    monkeypatch.setattr(
        browser.subprocess, 'run',
        lambda *a, **kw: calls.append((a, kw)) or fake_proc,
    )
    monkeypatch.setattr(
        browser.shutil, 'which', lambda name: '/fake/orca_plot' if name == 'orca_plot' else None,
    )
    refs['_orca_plot_calls'] = calls
    return refs


def _label_for(refs, name):
    for option in refs['calc_file_list'].options:
        if str(option).replace('\xa0', ' ').strip().endswith(name):
            return option
    raise AssertionError(f'{name} is not in {list(refs["calc_file_list"].options)}')


def _click(refs, name):
    label = _label_for(refs, name)
    refs['calc_open_input'].value = label
    refs['calc_file_list'].value = (label,)


def _status_text(refs):
    value = refs['calc_spin_density_status'].value
    # Strip the surrounding span/html for plain-text assertions.
    return value


def test_spin_density_option_appears_for_out_files(refs):
    _click(refs, 'stage.out')
    assert 'Spin Density' in refs['calc_options_dropdown'].options


def test_spin_density_refuses_closed_shell_job(refs):
    _click(refs, 'stage.out')
    refs['calc_spin_density_btn'].click()
    html = _status_text(refs)
    assert 'open-shell' in html.lower()
    assert refs['_orca_plot_calls'] == []


def test_spin_density_starts_uks_job_without_running_orca_plot_inline(refs):
    _click(refs, 'open.out')
    refs['calc_spin_density_btn'].click()
    # The work runs on a daemon thread; wait briefly for it to finish.
    import time
    for _ in range(50):
        time.sleep(0.1)
        if refs['calc_spin_density_btn'].disabled is False:
            break
    assert refs['_orca_plot_calls'], 'expected orca_plot to be invoked for UKS'
    args, kwargs = refs['_orca_plot_calls'][0]
    cmd = args[0] if args else kwargs.get('args')
    assert 'orca_plot' in ' '.join(map(str, cmd))
    stdin_payload = kwargs.get('input', b'')
    lines = stdin_payload.decode().strip().split('\n')
    # Empirically verified against ORCA 6.1.1: 1 = enter plot type,
    # 3 = (scf) spin density, y = accept the default <base>.scfr density,
    # 4 = set the grid interval count (panel default 80), 11 = generate
    # the plot, 12 = exit.
    assert lines == ['1', '3', 'y', '4', '80', '11', '12']


def test_spin_density_passes_selected_grid_resolution(refs, monkeypatch, tmp_path):
    _click(refs, 'open.out')
    refs['calc_spin_density_grid'].value = 120
    # A "successful" orca_plot that also creates a cube file, so the handler
    # completes; the recording fake from the fixture stays in place.
    job_dir = tmp_path / 'calc' / 'job'
    calls = []
    def _fake_run(*args, **kwargs):
        calls.append((args, kwargs))
        (job_dir / 'open.spindens.cube').write_text('cube\n', encoding='utf-8')
        return types.SimpleNamespace(returncode=0, stdout=b'', stderr=b'')
    monkeypatch.setattr(browser.subprocess, 'run', _fake_run)
    refs['calc_spin_density_btn'].click()
    import time
    for _ in range(50):
        time.sleep(0.1)
        if refs['calc_spin_density_btn'].disabled is False:
            break
    payload = calls[0][1]['input'].decode().strip().split('\n')
    # Menu 4 = grid interval count; the selected 120 must reach orca_plot.
    assert payload == ['1', '3', 'y', '4', '120', '11', '12']


def test_spin_density_display_change_reruns_nothing(refs, monkeypatch, tmp_path):
    """A pure display change (Update View) re-renders from the cached cube
    and must not invoke orca_plot a second time."""
    _click(refs, 'open.out')
    # Simulate an existing spin display with a cube on disk.
    cube_path = tmp_path / 'calc' / 'job' / 'open.spindens.cube'
    cube_path.write_text('dummy cube content\n', encoding='utf-8')
    state = refs['xyz_batch_state']
    state['spin_display'] = {
        'cube_path': str(cube_path),
        'isoval': 0.005,
        'opacity': 0.65,
        'show_pos': True,
        'show_neg': True,
        'show_mol': True,
        'pos_color': '#0000ff',
        'neg_color': '#ff0000',
    }
    calls_before = len(refs['_orca_plot_calls'])
    monkeypatch.setattr(
        browser.shutil, 'which',
        lambda name: '/fake/orca_plot' if name == 'orca_plot' else None)
    refs['calc_on_spin_density_update_view']()
    # No new orca_plot invocation.
    assert len(refs['_orca_plot_calls']) == calls_before
    # The display settings were applied from the panel widgets.
    spin_display = state['spin_display']
    assert spin_display['isoval'] == refs['calc_spin_density_isoval'].value
    assert spin_display['opacity'] == refs['calc_spin_density_opacity'].value


def test_spin_density_panel_defaults(refs):
    assert refs['calc_spin_density_grid'].value == 80
    assert refs['calc_spin_density_isoval'].value == 0.005
    assert refs['calc_spin_density_opacity'].value == 0.65
    assert refs['calc_spin_density_pos_cb'].value is True
    assert refs['calc_spin_density_neg_cb'].value is True
    assert refs['calc_spin_density_mol_cb'].value is True
    assert refs['calc_spin_density_pos_color'].value == '#0000ff'
    assert refs['calc_spin_density_neg_color'].value == '#ff0000'

