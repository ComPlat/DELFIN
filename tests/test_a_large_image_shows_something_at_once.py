"""A big image must not hold the tab while it travels.

An image goes to the page as base64 inside the markup, which inflates it by a
third and then has to be parsed as HTML. Measured in a running dashboard: an
18 MB photograph took 9.8 seconds before anything appeared, against 0.43 s for
a 4 MB text file and 0.42 s for a sixty-page PDF -- the file browser felt slow
at exactly one file type. Moving from that image to a PDF and back took 2.8 s.

Above a threshold a downscaled preview is drawn first and the original follows
on a background thread, which brought the same image to 1.0 s and left the full
resolution arriving a moment later. And .jpg, .gif, .bmp, .webp and .tiff had
an icon in the listing but no viewer at all -- only .png was ever drawn.
"""

from __future__ import annotations

import os
import re
import time
from pathlib import Path

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext

Image = pytest.importorskip('PIL.Image')


def _ctx(tmp_path: Path) -> DashboardContext:
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / 'calc', archive_dir=tmp_path / 'archive',
                           office_dir=tmp_path / 'office')
    ctx.run_js = lambda script: None
    return ctx


def _open(refs, name):
    refs['calc_list_directory']()
    file_list = refs['calc_file_list']
    match = [o for o in file_list.options if name in str(o)]
    assert match, f'{name} not in {file_list.options}'
    value = match[0][1] if isinstance(match[0], tuple) else match[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value


def _shown(refs):
    """(kind, pixel width) of the image on screen, from the markup itself."""
    html = refs['calc_content_area'].value
    kind = re.search(r'data:image/([a-z]+);base64,', html)
    if not kind:
        return None, 0
    import base64
    import io
    payload = re.search(r'base64,([^\']+)\'', html).group(1)
    with Image.open(io.BytesIO(base64.b64decode(payload))) as img:
        return kind.group(1), img.width


def _write_image(path, size, *, noisy=True, **save):
    data = os.urandom(size[0] * size[1] * 3) if noisy else bytes(size[0] * size[1] * 3)
    Image.frombytes('RGB', size, data).save(path, **save)
    return path


@pytest.mark.parametrize('name, kind, saver', [
    ('a.png', 'png', {}),
    ('a.jpg', 'jpeg', {'quality': 80}),
    ('a.bmp', 'bmp', {}),
    ('a.gif', 'gif', {}),
    ('a.webp', 'webp', {}),
])
def test_every_common_image_format_is_drawn(tmp_path, name, kind, saver):
    ctx = _ctx(tmp_path)
    _write_image(ctx.calc_dir / name, (120, 80), noisy=False, **saver)

    _widget, refs = browser.create_tab(ctx)
    _open(refs, name)

    shown, width = _shown(refs)
    assert shown is not None, f'{name} was not drawn at all'
    assert width == 120


def test_a_small_image_travels_as_itself(tmp_path):
    ctx = _ctx(tmp_path)
    _write_image(ctx.calc_dir / 'klein.png', (200, 150), noisy=False)

    _widget, refs = browser.create_tab(ctx)
    _open(refs, 'klein.png')

    kind, width = _shown(refs)
    assert (kind, width) == ('png', 200)
    assert 'preview' not in refs['calc_file_info'].value


def test_a_large_image_shows_a_preview_first_and_the_original_after(tmp_path):
    ctx = _ctx(tmp_path)
    big = _write_image(ctx.calc_dir / 'gross.png', (2600, 1800))     # a few MB of noise
    assert big.stat().st_size > browser_preview_threshold()

    _widget, refs = browser.create_tab(ctx)
    _open(refs, 'gross.png')

    kind, width = _shown(refs)
    assert kind == 'png'
    assert width <= 1800, 'the whole image was sent straight away'
    assert 'preview' in refs['calc_file_info'].value

    deadline = time.monotonic() + 30
    while time.monotonic() < deadline:                     # the thread swaps it in
        _kind, width = _shown(refs)
        if width == 2600:
            break
        time.sleep(0.1)
    assert width == 2600, 'the full image never arrived'
    assert 'preview' not in refs['calc_file_info'].value


def test_the_original_is_dropped_when_the_user_moved_on(tmp_path):
    """A slow read must not paint over whatever is on screen by then."""
    ctx = _ctx(tmp_path)
    _write_image(ctx.calc_dir / 'gross.png', (2600, 1800))
    _write_image(ctx.calc_dir / 'klein.png', (200, 150), noisy=False)

    _widget, refs = browser.create_tab(ctx)
    _open(refs, 'gross.png')
    _open(refs, 'klein.png')

    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        time.sleep(0.2)
        kind, width = _shown(refs)
        assert width in (200, 0), f'the large image painted over the small one ({width} px)'


def browser_preview_threshold():
    import inspect
    source = inspect.getsource(browser.create_tab)
    return int(re.search(r'CALC_IMAGE_PREVIEW_BYTES = ([\d_]+)', source).group(1).replace('_', ''))
