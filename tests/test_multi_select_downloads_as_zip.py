"""A multi-selection downloads as a single ZIP.

The file list is a ``SelectMultiple``, so shift/ctrl-click already selects
several entries. The Download button used to read only the first selected
item. It now reads the whole selection: multiple selected files and folders
are packed into one ZIP whose entries mirror the selection, while a single
selected file still downloads directly and a single folder still downloads
as its own ZIP.
"""

from __future__ import annotations

import base64
import io
import re
import zipfile

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext


def _build(ctx, *names):
    """Create the browser tab, list the root dir, select the named entries."""
    sent: list[str] = []
    ctx.run_js = lambda script: sent.append(script)
    _widget, refs = browser.create_tab(ctx)
    refs['calc_list_directory']()
    file_list = refs['calc_file_list']
    options = file_list.options
    picked = []
    for name in names:
        match = [opt for opt in options if name in str(opt)]
        assert match, f'{name!r} not in the file list options'
        opt = match[0]
        value = opt[1] if isinstance(opt, tuple) else opt
        picked.append(value)
    file_list.value = tuple(picked)
    return refs, sent


def _extract_zip_payload(sent):
    """Decode the base64 ZIP payload from the emitted download script."""
    for script in sent:
        m = re.search(r'const b64="([A-Za-z0-9+/=]+)"', script)
        if m:
            return base64.b64decode(m.group(1))
    return None


def test_multiple_files_download_as_one_zip(tmp_path):
    (tmp_path / 'alpha.in').write_text('A', encoding='utf-8')
    (tmp_path / 'beta.in').write_text('BB', encoding='utf-8')
    (tmp_path / 'gamma.in').write_text('CCC', encoding='utf-8')
    ctx = DashboardContext(
        calc_dir=tmp_path,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    refs, sent = _build(ctx, 'alpha.in', 'beta.in', 'gamma.in')

    # The targets helper must expose all three selected paths.
    assert {p.name for p in refs['calc_download_targets']()} == {
        'alpha.in', 'beta.in', 'gamma.in'
    }

    refs['calc_on_download'](refs['calc_download_btn'])
    status = refs['calc_download_status'].value
    assert 'Download started:' in status
    assert '.zip' in status

    payload = _extract_zip_payload(sent)
    assert payload is not None, 'no ZIP payload was emitted for a multi-select'
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        names = set(archive.namelist())
    assert {'alpha.in', 'beta.in', 'gamma.in'} <= names


def test_folders_are_packed_recursively_into_the_multi_zip(tmp_path):
    (tmp_path / 'top.txt').write_text('top', encoding='utf-8')
    sub = tmp_path / 'sub'
    sub.mkdir()
    (sub / 'inner.txt').write_text('inner', encoding='utf-8')
    (sub / 'deep').mkdir()
    (sub / 'deep' / 'leaf.txt').write_text('leaf', encoding='utf-8')
    ctx = DashboardContext(
        calc_dir=tmp_path,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    refs, sent = _build(ctx, 'top.txt', 'sub')
    refs['calc_on_download'](refs['calc_download_btn'])
    assert 'Download started:' in refs['calc_download_status'].value
    payload = _extract_zip_payload(sent)
    assert payload is not None
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        names = set(archive.namelist())
    # Every entry sits under its folder's own name, mirroring the selection.
    assert 'top.txt' in names
    assert 'sub/inner.txt' in names
    assert 'sub/deep/leaf.txt' in names


def test_single_selected_file_still_downloads_directly(tmp_path):
    (tmp_path / 'only.in').write_text('X', encoding='utf-8')
    ctx = DashboardContext(
        calc_dir=tmp_path,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    refs, sent = _build(ctx, 'only.in')
    assert [p.name for p in refs['calc_download_targets']()] == ['only.in']

    refs['calc_on_download'](refs['calc_download_btn'])
    assert 'Download started: only.in' in refs['calc_download_status'].value
    payload = _extract_zip_payload(sent)
    assert payload is not None
    # A direct file download is NOT a ZIP; it is the raw file content.
    assert payload == b'X'
    # ZIP archives start with the "PK\x03\x04" signature — a raw file must not.
    assert not payload.startswith(b'PK\x03\x04')


def test_single_folder_still_downloads_as_its_own_zip(tmp_path):
    (tmp_path / 'work').mkdir()
    (tmp_path / 'work' / 'a.txt').write_text('a', encoding='utf-8')
    ctx = DashboardContext(
        calc_dir=tmp_path,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    refs, sent = _build(ctx, 'work')
    refs['calc_on_download'](refs['calc_download_btn'])
    assert 'Download started: work.zip' in refs['calc_download_status'].value
    payload = _extract_zip_payload(sent)
    assert payload is not None
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        assert set(archive.namelist()) >= {'work/a.txt'}


def test_no_selection_and_no_folder_reports_nothing_selected(tmp_path):
    # An empty root with nothing selected leaves no download target at all.
    ctx = DashboardContext(
        calc_dir=tmp_path,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    refs, sent = _build(ctx)
    assert refs['calc_download_targets']() == []
    refs['calc_on_download'](refs['calc_download_btn'])
    assert 'No file/folder selected.' in refs['calc_download_status'].value
    assert _extract_zip_payload(sent) is None, 'no ZIP payload should be emitted'


def test_download_button_labels_the_selection(tmp_path):
    (tmp_path / 'a.in').write_text('A', encoding='utf-8')
    (tmp_path / 'b.in').write_text('B', encoding='utf-8')
    folder = tmp_path / 'folder'
    folder.mkdir()
    ctx = DashboardContext(
        calc_dir=tmp_path,
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    refs, sent = _build(ctx)

    # Nothing selected, and at the root there is no current folder to fall
    # back to -> the button is disabled.
    refs['calc_update_download_btn']()
    assert refs['calc_download_btn'].disabled is True

    # Select a single file -> "Download" (direct).
    _select(refs, 'a.in')
    refs['calc_update_download_btn']()
    assert refs['calc_download_btn'].disabled is False
    assert refs['calc_download_btn'].description == 'Download'

    # Select the folder on top -> "Download ZIP".
    _select(refs, 'a.in', 'folder')
    refs['calc_update_download_btn']()
    assert refs['calc_download_btn'].disabled is False
    assert refs['calc_download_btn'].description == 'Download ZIP'


def _select(refs, *names):
    file_list = refs['calc_file_list']
    options = file_list.options
    picked = []
    for name in names:
        match = [opt for opt in options if name in str(opt)]
        assert match, f'{name!r} not in the file list options'
        opt = match[0]
        picked.append(opt[1] if isinstance(opt, tuple) else opt)
    file_list.value = tuple(picked)
