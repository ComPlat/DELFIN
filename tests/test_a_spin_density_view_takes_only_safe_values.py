"""The spin density view draws only what the panel and orca_plot produced.

Three properties of the Spin Density panel (calculations browser):

- The isovalue is written into the viewer's script as a +/- pair. A
  negative entry (rendered as ``--x``, a script syntax error) is drawn
  by its magnitude; a non-finite one falls back to the default. Colours
  need no such guard: the ColorPicker trait refuses non-colour values.
- A cube is shown as the spin density only if this orca_plot run wrote
  it. An older cube in the same folder (an MO plot, for example) is
  never displayed when the run produced nothing.
- "Show Molecule" off applies to the first render after generation, not
  only to a later Apply.
"""

import time
import types

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


def _cube_text(tag):
    # Large enough to pass the "non-empty cube" size filter.
    return f'{tag}\n' + ('0.0 ' * 40) + '\n'


@pytest.fixture
def refs(tmp_path, monkeypatch):
    calc = tmp_path / 'calc'
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    job = calc / 'job'
    job.mkdir()
    (job / 'open.out').write_text(
        'Hartree-Fock type          UKS\n'
        'FINAL SINGLE POINT ENERGY -1.0\n', encoding='utf-8')
    (job / 'open.gbw').write_bytes(b'\x00' * 16)
    (job / 'conf_00.spindens.cube').write_text(_cube_text('spin'),
                                               encoding='utf-8')

    # The viewer may be switched off in the host's settings; the script
    # under test is only built when it is on.
    real_profile = browser.get_viewer_profile
    monkeypatch.setattr(browser, 'get_viewer_profile',
                        lambda: {**real_profile(), 'enabled': True})
    shown = []
    monkeypatch.setattr(
        browser, 'display',
        lambda obj, *a, **kw: shown.append(getattr(obj, 'data', str(obj))))

    ctx = DashboardContext(calc_dir=calc, archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda script: None
    _widget, refs = browser.create_tab(ctx)
    refs['calc_path_input'].value = str(job)
    monkeypatch.setattr(
        browser.shutil, 'which',
        lambda name: '/fake/orca_plot' if name == 'orca_plot' else None)
    refs['_shown'] = shown
    refs['_job'] = job
    return refs


def _click(refs, name):
    for option in refs['calc_file_list'].options:
        if str(option).replace('\xa0', ' ').strip().endswith(name):
            refs['calc_open_input'].value = option
            refs['calc_file_list'].value = (option,)
            return
    raise AssertionError(f'{name} is not listed')


def _plot(refs, monkeypatch, writes=None):
    def _fake_run(*args, **kwargs):
        if writes:
            (refs['_job'] / writes).write_text(_cube_text('new spin'),
                                               encoding='utf-8')
        return types.SimpleNamespace(returncode=0, stdout=b'', stderr=b'')
    monkeypatch.setattr(browser.subprocess, 'run', _fake_run)
    _click(refs, 'open.out')
    refs['_shown'].clear()
    refs['calc_spin_density_btn'].click()
    for _ in range(100):
        if refs['calc_spin_density_btn'].disabled is False:
            break
        time.sleep(0.05)


def _last_script(refs):
    scripts = [s for s in refs['_shown'] if 'addVolumetricData' in s]
    assert scripts, 'no cube was rendered'
    return scripts[-1]


@pytest.mark.parametrize('value, drawn', [(-0.01, '0.01'),
                                          (float('nan'), '0.005'),
                                          (float('inf'), '0.005')])
def test_an_isovalue_is_drawn_as_a_finite_positive_pair(refs, value, drawn):
    _click(refs, 'conf_00.spindens.cube')
    refs['calc_spin_density_isoval'].value = value
    refs['calc_spin_density_apply_btn'].click()
    script = _last_script(refs)
    assert f'isoval:{drawn},' in script
    assert f'isoval:-{drawn},' in script
    assert '--' not in script.split('addVolumetricData', 1)[1][:200]
    assert 'isoval:nan' not in script and 'isoval:inf' not in script


def test_an_older_cube_is_not_shown_when_the_run_wrote_none(refs, monkeypatch):
    (refs['_job'] / 'orbital_12a.cube').write_text(_cube_text('mo'),
                                                    encoding='utf-8')
    _plot(refs, monkeypatch, writes=None)
    assert not [s for s in refs['_shown'] if 'addVolumetricData' in s]
    assert 'No spin density cube' in refs['calc_spin_density_status'].value


def test_the_cube_the_run_wrote_is_the_one_shown(refs, monkeypatch):
    (refs['_job'] / 'orbital_12a.cube').write_text(_cube_text('mo'),
                                                    encoding='utf-8')
    _plot(refs, monkeypatch, writes='open.spindens.cube')
    assert 'new spin' in _last_script(refs)
    assert (refs['xyz_batch_state']['spin_display']['cube_path']
            .endswith('open.spindens.cube'))


def test_show_molecule_off_holds_for_the_first_render(refs, monkeypatch):
    refs['calc_spin_density_mol_cb'].value = False
    _plot(refs, monkeypatch, writes='open.spindens.cube')
    assert 'viewer.setStyle({}, {});' in _last_script(refs)
