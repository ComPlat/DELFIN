"""TS->GOAT panel lifecycle regression tests (2026-11-02 user report).

Covers: file-change closes the panel, close/cancel buttons restore the
normal view, constraint field -> viewer sync, reset/clear buttons and
the viewer bootstrap retry loop (the "white viewer" race).
"""

import sys
from pathlib import Path

import pytest

widgets = pytest.importorskip("ipywidgets")

# Make this file runnable from any working directory: the auto-check
# runs pytest from another checkout, where a plain `from tests...`
# would resolve that checkout's tests package instead of this one.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tests.test_ts_goat_dialog_regression import (  # noqa: E402
    _FakeBackend,
    _open_dialog,
)


def test_close_button_restores_normal_view(tmp_path):
    """✕ Schließen hides the panel and resets the options dropdown."""
    captured = {}
    out, _r = _open_dialog(tmp_path, _FakeBackend(), extra_callbacks={
        'on_close': lambda: captured.setdefault('closed', True),
    })
    assert out['error'] is None
    reg = out['registry']
    reg['close']()
    assert captured.get('closed') is True
    assert reg['panel'].layout.display == 'none'


def test_reset_and_clear_buttons(tmp_path):
    """The reset button restores the proposal; clearing empties it."""
    out, _r = _open_dialog(tmp_path, _FakeBackend())
    reg = out['registry']
    field = reg['constraints_field']
    rc = reg['state']['proposed']
    indices = list(getattr(rc, 'indices', None) or [])
    assert indices, 'synthetic TS must propose a reactive centre'
    # clear everything
    field.value = ''
    assert reg['parse_constraints'](field.value) == []
    # restore the proposal
    field.value = ', '.join(str(i) for i in indices)
    assert reg['parse_constraints'](field.value) == indices


def test_constraint_field_edit_updates_viewer(tmp_path):
    """Editing the constraint field re-renders the viewer (sync)."""
    calls = []
    out, _r = _open_dialog(tmp_path, _FakeBackend(), extra_callbacks={
        'render_in_panel': lambda atoms, indices, arrows=None,
                           output_widget=None, mode_vector=None:
                           calls.append(list(indices)),
    })
    reg = out['registry']
    field = reg['constraints_field']
    n_calls = len(calls)
    field.value = '0'
    assert len(calls) > n_calls, 'field edit must trigger a re-render'
    assert calls[-1] == [0]


def test_bootstrap_js_retries_until_the_output_is_synced():
    """The white-viewer race: run_js fires before the Output's display()
    reaches the DOM.  The bootstrap must retry instead of giving up."""
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_bootstrap_js
    js = ts_goat_viewer_bootstrap_js('v_test', click_sync_css='.x',
                                     arrows=[0.0, 0.0, 0.0])
    assert 'function boot(){' in js
    assert 'setTimeout(boot, 250)' in js
    assert js.count('boot();') == 1  # called exactly once, not inline-run
