"""Calculations, Archive and Office may not name their 3D stages alike.

The three tabs are three builds of one builder rendering into one document.
The viewer scripts look their stage up with a bare ``getElementById``, which
answers with the first match in the document -- Calculations, built first. So
when every instance numbered its stages from one, opening a structure in
Archive initialised it into the Calculations viewer and Calculations was left
showing nothing. Which tab lost its structure depended on which had rendered
last, so it only went wrong some of the time.
"""

from __future__ import annotations

import re

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext

# Every id the emitted viewer scripts address by getElementById.
_STAGE_ID = re.compile(r'getElementById\("([A-Za-z_]+\d+)"\)')
_WRAPPER_ID = re.compile(r'id="(calc_mol_wrap[A-Za-z_]*\d+)"')


def _build(tmp_path, folder, monkeypatch):
    """Build one browser instance on its own folder and open the structure."""
    root = tmp_path / folder
    root.mkdir(exist_ok=True)
    (root / 'mol.xyz').write_text(
        '3\nwater\nO 0.0 0.0 0.0\nH 0.757 0.586 0.0\nH -0.757 0.586 0.0\n',
        encoding='utf-8',
    )

    # The viewer markup is handed to display() inside an Output widget. With
    # no kernel behind that widget it is not recorded anywhere, so intercept
    # display itself -- that is the markup the browser would receive.
    shown: list[str] = []
    real_display = browser.display

    def capture(*objs, **kw):
        for obj in objs:
            data = getattr(obj, 'data', None)
            if isinstance(data, str):
                shown.append(data)
        return real_display(*objs, **kw)

    monkeypatch.setattr(browser, 'display', capture)

    sent: list[str] = []
    ctx = DashboardContext(calc_dir=root, archive_dir=root, office_dir=root)
    ctx.run_js = lambda script: sent.append(script)
    _widget, refs = browser.create_tab(ctx)
    refs['calc_list_directory']()

    file_list = refs['calc_file_list']
    target = [o for o in file_list.options if 'mol.xyz' in str(o)]
    assert target, f'mol.xyz missing from {file_list.options}'
    value = target[0][1] if isinstance(target[0], tuple) else target[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value

    # The viewer markup is displayed into an Output widget, so it lands in its
    # `outputs`, not in a `value`.
    html = []
    seen = set()

    def walk(node):
        if id(node) in seen:
            return
        seen.add(id(node))
        value = getattr(node, 'value', None)
        if isinstance(value, str):
            html.append(value)
        for out in (getattr(node, 'outputs', None) or ()):
            data = (out or {}).get('data') or {}
            for mime in ('text/html', 'text/plain'):
                if isinstance(data.get(mime), str):
                    html.append(data[mime])
        for child in (getattr(node, 'children', None) or ()):
            walk(child)

    walk(_widget)
    return '\n'.join(html + shown + sent)


def _ids(markup):
    return set(_STAGE_ID.findall(markup)) | set(_WRAPPER_ID.findall(markup))


def test_two_instances_never_name_a_stage_alike(tmp_path, monkeypatch):
    first = _ids(_build(tmp_path, 'calc', monkeypatch))
    second = _ids(_build(tmp_path, 'archive', monkeypatch))
    assert first, 'the first instance emitted no stage id at all'
    assert second, 'the second instance emitted no stage id at all'
    shared = first & second
    assert not shared, (
        f'both instances address {sorted(shared)}; a bare getElementById '
        'answers with whichever tab was built first, so one of them renders '
        'into the other tab and is left blank'
    )


def test_the_counter_is_shared_rather_than_per_instance():
    """A counter created inside create_tab restarts for every tab."""
    assert hasattr(browser, '_mol3d_counter'), (
        'the stage counter is not module level, so each tab starts numbering '
        'from one again and the ids collide'
    )
    before = browser._mol3d_counter[0]
    assert isinstance(before, int)
