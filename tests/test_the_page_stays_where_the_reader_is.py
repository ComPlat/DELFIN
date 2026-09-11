"""A click in Settings must not move the page the reader is on.

Measured in the real dashboard (headless Chromium under Voila), with
Settings > Dashboard Tabs open and a row near the bottom of the list in view:

- ticking "Visible" moved the page from 983 px to 225 px and the box out of
  view. Every click replaced the whole list of rows, so for a moment the page
  was shorter than the window and the browser clamped the scroll position.
- Save put the page at 0 and opened Pipelines instead of Settings. The
  browser rebuilds every tab from the first one that changed, Settings is
  pinned last, and the first rebuilt tab to render took the selection.

After the change the tick left the page at 983 px, and Save left it at
274 px and on Settings. The scroll keeper on its own, over the old row list,
held the tick at 1012 px instead of dropping to 225 px; it is there for every
place that still replaces widgets.
"""

import inspect
import types

import ipywidgets as widgets
import pytest

from delfin.dashboard import helpers, tab_settings


def _specs(count):
    specs = [
        {'id': f't{i}', 'title': f'Tab {i}', 'default_order': i, 'available': True}
        for i in range(count)
    ]
    specs.append({'id': 'settings', 'title': 'Settings', 'fixed': True,
                  'default_order': 999, 'available': True})
    return specs


@pytest.fixture
def settings_tab(tmp_path, monkeypatch):
    monkeypatch.setenv('HOME', str(tmp_path))
    ctx = types.SimpleNamespace(
        calc_dir=tmp_path / 'calc', archive_dir=tmp_path / 'archive',
        office_dir=tmp_path / 'office', default_calc_dir=tmp_path / 'calc',
        default_archive_dir=tmp_path / 'archive', repo_dir=tmp_path,
        tab_specs=_specs(5), runtime_backend='local', backend=None, orca_base='',
        orca_candidates=[], submit_templates_dir=tmp_path / 'templates',
        runtime_settings={}, tabs_widget=None, tab_indices={},
        rebuild_dashboard_tabs=lambda **_kw: None,
    )
    tab, api = tab_settings.create_tab(ctx)
    return ctx, tab, api


def _rows_box(tab):
    found = []

    def walk(widget):
        children = getattr(widget, 'children', ()) or ()
        if children and all(isinstance(row, widgets.HBox) for row in children) and any(
            isinstance(part, widgets.Checkbox) and part.description == 'Visible'
            for row in children for part in row.children
        ):
            found.append(widget)
        for child in children:
            walk(child)

    walk(tab)
    assert len(found) == 1
    return found[0]


def _rows(box):
    out = []
    for row in box.children:
        label, visible, up, down = row.children
        out.append((label.value.split('</b>')[0].replace('<b>', ''), visible.value,
                    'hidden' in label.value))
    return out


def test_ticking_a_tab_changes_its_row_and_nothing_else(settings_tab):
    _ctx, tab, _api = settings_tab
    box = _rows_box(tab)
    before = box.children
    box.children[2].children[1].value = False
    assert all(a is b for a, b in zip(box.children, before)), 'the rows were replaced'
    assert _rows(box)[2] == ('Tab 2', False, True)
    assert [hidden for *_rest, hidden in _rows(box)].count(True) == 1


def test_moving_a_tab_keeps_the_rows_and_carries_its_state(settings_tab):
    _ctx, tab, _api = settings_tab
    box = _rows_box(tab)
    before = box.children
    box.children[2].children[1].value = False
    box.children[2].children[2].click()   # Up
    assert all(a is b for a, b in zip(box.children, before)), 'the rows were replaced'
    rows = _rows(box)
    # Tab 2 moved up with its box unticked; Tab 1 took its place, still ticked.
    assert rows[1] == ('Tab 2', False, True)
    assert rows[2] == ('Tab 1', True, False)
    first_up = box.children[0].children[2]
    last_down = box.children[4].children[3]
    assert first_up.disabled and last_down.disabled


def test_a_different_number_of_tabs_builds_the_rows_anew(settings_tab):
    ctx, tab, api = settings_tab
    box = _rows_box(tab)
    ctx.tab_specs = _specs(7)
    api['reload_settings'](set_status=False)
    assert [title for title, *_rest in _rows(box)] == [f'Tab {i}' for i in range(7)] + ['Settings']


def _tabs(count):
    children = [widgets.HTML(f'tab {i}') for i in range(count)]
    tabs = widgets.Tab(children=children)
    return tabs, children


def _record(tabs):
    events = []
    tabs.observe(lambda change: events.append(('selected', change['new'])), names='selected_index')
    tabs.observe(lambda change: events.append(('children', len(change['new']))), names='children')
    return events


def test_the_selection_waits_on_a_kept_tab_while_the_rest_is_rebuilt():
    tabs, children = _tabs(5)
    tabs.selected_index = 4                      # Settings, pinned last
    events = _record(tabs)
    new = children[:2] + children[3:]            # the third tab was hidden
    helpers.hand_over_tabs(tabs, new, [f't{i}' for i in range(4)], selected_index=3)
    assert events == [('selected', 0), ('children', 4), ('selected', 3)]
    assert tabs.selected_index == 3


def test_a_selection_among_the_kept_tabs_is_not_moved():
    tabs, children = _tabs(5)
    tabs.selected_index = 1
    events = _record(tabs)
    helpers.hand_over_tabs(tabs, children[:3] + children[4:], ['a', 'b', 'c', 'd'], selected_index=1)
    assert events == [('children', 4)]


def test_unchanged_children_are_not_handed_over_again():
    tabs, children = _tabs(3)
    events = _record(tabs)
    helpers.hand_over_tabs(tabs, list(children), ['a', 'b', 'c'], selected_index=2)
    assert events == [('selected', 2)]
    assert tuple(tabs.titles) == ('a', 'b', 'c')


def test_with_no_tab_kept_the_first_drift_is_corrected_once():
    tabs, children = _tabs(3)
    tabs.selected_index = 2
    new = [widgets.HTML('new first')] + children[1:]
    helpers.hand_over_tabs(tabs, new, ['x', 'b', 'c'], selected_index=2)
    tabs.selected_index = 1                      # what the browser reports back
    assert tabs.selected_index == 2
    tabs.selected_index = 0                      # the reader's own click later on
    assert tabs.selected_index == 0


def test_the_dashboard_installs_the_scroll_keeper():
    from delfin import dashboard

    assert 'keep_scroll_position(ctx)' in inspect.getsource(dashboard)


def test_the_keeper_gives_way_to_the_reader(monkeypatch):
    sent = []
    monkeypatch.setattr(helpers, '_append_js', lambda _ctx, script: sent.append(script))
    helpers.keep_scroll_position(object())
    script = sent[0]
    # Every kind of reader input cancels a pending correction ...
    for event in ('wheel', 'touchstart', 'touchmove', 'keydown', 'pointerdown'):
        assert f"'{event}'" in script, event
    # ... and only a widget view removed or hidden arms one.
    assert 'removedNodes' in script and 'lm-mod-hidden' in script
