"""Three file browsers, one page, and each one reaches only its own DOM.

Calculations, Archive and Office are three instances of the same builder
(``tab_archive_statistics`` and ``tab_office`` both call
``tab_calculations_browser.create_tab`` on a different root). They render
the same markup into one document.

Anything the emitted JavaScript looks up by document therefore answers with
whichever instance was built first. That is not hypothetical: ``⬆ Top``,
``⬇ End``, the search jump and the chunk loader all resolved
``#calc-content-box``, so in Office they scrolled the Calculations tab and
looked like dead buttons.

The scripts are rewritten as they are emitted, which is why these tests
check what is *sent*, not what is written in the source: the source still
reads ``document.getElementById('calc-content-box')`` and is meant to.
"""

from __future__ import annotations

import pathlib
import re
import tempfile

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.context import DashboardContext

# A lookup that answers with the first matching element on the page,
# regardless of which browser instance asked.
_DOCUMENT_WIDE = (
    "document.getElementById('calc-",
    "document.querySelector('.calc-",
    "document.querySelectorAll('.calc-",
    "window.__calcChunk",
)


@pytest.fixture
def roots(tmp_path):
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
        # A file to open: the scroll and search scripts are only emitted
        # once something is actually on screen.
        (tmp_path / name / 'note.txt').write_text(
            'line\n' * 200, encoding='utf-8')
    return tmp_path


def _build(roots, which):
    """Build one browser instance; return (scope id, scripts it emitted)."""
    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=roots / which,
        archive_dir=roots / 'archive',
        office_dir=roots / 'office',
    )
    ctx.run_js = lambda script: sent.append(script)
    _widget, refs = browser.create_tab(ctx)
    # Listing a directory is what paints the content area and its scripts.
    refs['calc_list_directory']()
    file_list = refs['calc_file_list']
    target = [opt for opt in file_list.options if 'note.txt' in str(opt)]
    assert target, f'note.txt did not show up in {file_list.options}'
    value = target[0][1] if isinstance(target[0], tuple) else target[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value
    # Searching is the path that jumped the wrong tab.
    refs['calc_search_input'].value = 'line'
    scripts = sent + list(ctx.init_js_parts)
    scopes = re.findall(r'calc-scope-\d+', '\n'.join(scripts))
    assert scopes, 'the tab emitted no scope-bearing script'
    return scopes[0], scripts


def test_emitted_scripts_never_look_up_by_document(roots):
    scope, scripts = _build(roots, 'office')
    for script in scripts:
        for pattern in _DOCUMENT_WIDE:
            if pattern not in script:
                continue
            # The only permitted document-wide lookups are the scope root
            # itself and the resolvers that take a scope argument.
            for hit in re.finditer(re.escape(pattern) + r"[^')]*", script):
                text = hit.group(0)
                assert 'calc-scope-' in text, (
                    f'{text!r} resolves against the whole page, so it answers '
                    'with whichever browser instance was built first')


def test_the_resolvers_are_defined(roots):
    _scope, scripts = _build(roots, 'office')
    startup = '\n'.join(scripts)
    for helper in ('__delfinCalcQ', '__delfinCalcQA', '__delfinCalcS'):
        assert f'window.{helper} = window.{helper} ||' in startup, (
            f'{helper} is not defined, so every rewritten lookup is a crash')


def test_content_lookups_go_through_the_scope(roots):
    _scope, scripts = _build(roots, 'office')
    assert 'window.__delfinCalcQ(' in '\n'.join(scripts), (
        'nothing resolved through the scope, so either the rewrite stopped '
        'happening or this test stopped exercising the content area')


def test_chunk_state_is_kept_per_browser(tmp_path):
    """Above the full-read limit the file is loaded in windows, and the
    loader keeps its bookkeeping in JavaScript. Shared between two open
    browsers, one tab's scrolling cancels the other tab's pending load."""
    for name in ('calc', 'archive', 'office'):
        (tmp_path / name).mkdir()
    big = tmp_path / 'office' / 'big.log'
    big.write_text('line\n' * 500_000, encoding='utf-8')  # ~2.4 MiB
    # Past the 2 MiB full-read limit the browser switches to windowed loading.
    assert big.stat().st_size > 2 * 1024 * 1024

    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / 'office',
        archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office',
    )
    ctx.run_js = lambda script: sent.append(script)
    _widget, refs = browser.create_tab(ctx)
    refs['calc_list_directory']()
    file_list = refs['calc_file_list']
    target = [opt for opt in file_list.options if 'big.log' in str(opt)]
    assert target
    value = target[0][1] if isinstance(target[0], tuple) else target[0]
    file_list.value = (value,) if isinstance(file_list.value, tuple) else value

    emitted = '\n'.join(sent)
    assert 'window.__delfinCalcS(' in emitted, 'the chunk loader ran unscoped'
    assert 'window.__calcChunk' not in emitted


def test_two_browsers_resolve_against_different_roots(roots):
    calc_scope, calc_scripts = _build(roots, 'calc')
    office_scope, office_scripts = _build(roots, 'office')

    assert calc_scope != office_scope

    def scopes_used(scripts):
        return set(re.findall(r'calc-scope-\d+', '\n'.join(scripts)))

    # Neither instance may name the other's scope anywhere.
    assert office_scope not in scopes_used(calc_scripts)
    assert calc_scope not in scopes_used(office_scripts)


def test_the_content_box_is_addressed_by_class_not_id(roots):
    """Two elements cannot share an id. Rendering the same id from three
    instances is what made the lookups ambiguous in the first place."""
    _scope, scripts = _build(roots, 'office')
    source = pathlib.Path(browser.__file__).read_text(encoding='utf-8')
    for name in ('calc-content-box', 'calc-content-text',
                 'calc-chunk-top-spacer', 'calc-chunk-bottom-spacer',
                 'calc-current-match'):
        assert f"id='{name}'" not in source and f'id="{name}"' not in source, (
            f'{name} is still rendered as an id; three browsers would render '
            'three elements with the same id')


def test_the_splitter_binds_without_borrowing_another_tab(roots):
    """The scoped root is the only root. Falling back to the first
    .calc-tab on the page would bind this tab's splitter to another tab."""
    _scope, scripts = _build(roots, 'office')
    startup = '\n'.join(scripts)
    assert '.calc-splitter' in startup, 'no splitter binding was registered'
    assert "document.querySelector('.calc-tab')" not in startup
