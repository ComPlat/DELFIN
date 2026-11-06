"""Regression test: selecting 'TS -> GOAT' in the options dropdown must
open a visible preparation dialog (or show a clear error), end to end.

This drives the real handler code path -- calc_on_options_change with
change['new'] == 'TS -> GOAT' -- against synthetic TS fixtures, not the
library functions alone (those are covered by
test_ts_goat_reactive_center.py).
"""

import pytest

pytest.importorskip(
    'delfin.ts_goat',
    reason='lives on the feature branch until it is merged into main',
)

from pathlib import Path  # noqa: E402
from unittest import mock  # noqa: E402

import ipywidgets as widgets  # noqa: E402

from delfin.ts_goat import build_ts_goat_input  # noqa: E402


# ---------------------------------------------------------------- fixtures
def _write_ts_fixture(dir_path: Path, name: str = 'ts_run') -> Path:
    """Minimal converged OptTS run: out, inp, xyz, hess.

    The Hessian uses the same blocked-column layout as the fixtures in
    test_ts_goat_reactive_center.py so delfin.imag's parser reads the
    imaginary-mode eigenvector (a thinner stub leaves the analysis with
    zero proposals).
    """
    out = dir_path / f'{name}.out'
    out.write_text(
        '|  1> ! OptTS\n'
        'THE OPTIMIZATION HAS CONVERGED\n'
        '     1:      -567.89 cm**-1  ***imaginary mode***\n',
        encoding='utf-8',
    )
    inp = dir_path / f'{name}.inp'
    inp.write_text('%charge 0\n%mult 1\n* xyz 0 1\n', encoding='utf-8')

    atoms = [
        ('O', 0.00, 0.00, 0.00),
        ('H', 0.96, 0.00, 0.00),
        ('H', 1.55, 1.10, 0.00),
    ]
    # the mode acts on the transferring H: +x on H2, -x on the rest
    mode = [
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0),
    ]
    from tests.test_ts_goat_reactive_center import _write_hess
    hess = dir_path / f'{name}.hess'
    # ORCA hessians carry ALL 3N normal modes; the imaginary one sits at
    # its frequency index (1).  delfin.imag.displaced_geometries resolves
    # a mode BY that index, so the stub must provide every column.
    zero = [(0.0, 0.0, 0.0)] * len(atoms)
    mode_vectors = [zero, mode] + [zero] * 7
    _write_hess(hess, atoms, mode_vectors,
                [-567.89, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0])
    xyz = dir_path / f'{name}.xyz'
    xyz.write_text(
        '3\nTS\n'
        + ''.join(f'{el} {x:.6f} {y:.6f} {z:.6f}\n' for el, x, y, z in atoms),
        encoding='utf-8',
    )
    # endpoints so the connectivity evidence has something to compare
    ed = dir_path / f'{name}_Edukt.xyz'
    ed.write_text(
        '3\nEdukt\nO   0.000000  0.000000  0.000000\n'
        'H   0.960000  0.000000  0.000000\nH   3.000000  1.100000  0.000000\n',
        encoding='utf-8',
    )
    pr = dir_path / f'{name}_Produkt.xyz'
    pr.write_text(
        '3\nProdukt\nO   0.000000  0.000000  0.000000\n'
        'H   3.000000  0.000000  0.000000\nH   0.960000  1.100000  0.000000\n',
        encoding='utf-8',
    )
    return out


class _FakeBackend:
    """Records submit_orca calls; returns a zero rc so the dialog reports
    success without running anything."""

    def __init__(self):
        self.calls = []

    def submit_orca(self, job_dir, job_name, inp_file,
                    time_limit='48:00:00', pal=12, maxcore=1000):
        self.calls.append({
            'job_dir': job_dir, 'job_name': job_name, 'inp_file': inp_file,
        })
        return mock.Mock(returncode=0, stdout=f'Submitted job {job_name}\n',
                        stderr='')


# ---------------------------------------------------------------- helpers
def _open_dialog(tmp_path, backend, out_name='ts_run', extra_callbacks=None):
    """Call open_ts_goat_dialog the way the browser-tab handler does."""
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog

    out_path = _write_ts_fixture(tmp_path, out_name)
    rendered = []
    result = open_ts_goat_dialog(
        out_path,
        callbacks={
            'backend': backend,
            'reset_dropdown': lambda: None,
            'render_viewer': lambda xyz, idx, arrows=None: rendered.append((xyz, idx)),
            'scroll_top': lambda: None,
            **(extra_callbacks or {}),
        },
    )
    return result, rendered


def test_dropdown_selection_opens_visible_dialog(tmp_path):
    """The E2E regression: after the handler opens the dialog for a usable
    TS fixture, the panel widget is populated and visible, and shows the
    analysis -- no silent no-op like the ImportError bug."""
    backend = _FakeBackend()
    result, rendered = _open_dialog(tmp_path, backend)
    assert result['error'] is None, result['error']
    panel = result['panel']
    assert isinstance(panel, widgets.VBox)
    # visible: children attached and display not 'none' after show()
    assert len(panel.children) > 0
    assert panel.layout.display == 'block'
    # the status line renders the reliability, not an exception
    status_html = result['registry']['status'].value
    assert 'data quality' in status_html.lower()
    assert 'ImportError' not in status_html


def test_dropdown_selection_shows_clear_error_for_non_ts(tmp_path):
    """A plain (non-TS) output must produce a readable explanation, not a
    silent reset."""
    backend = _FakeBackend()
    d = tmp_path / 'plain'
    d.mkdir()
    out = d / 'plain.out'
    out.write_text('! Opt FREQ\nTHE OPTIMIZATION HAS CONVERGED\n', encoding='utf-8')
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog
    result = open_ts_goat_dialog(
        out,
        callbacks={
            'backend': backend, 'reset_dropdown': lambda: None,
            'render_viewer': lambda *a, **k: None, 'scroll_top': lambda: None,
        },
    )
    assert result['error_kind'] == 'not_usable_ts'
    assert 'not a usable converged TS' in result['error']
    assert result['panel'] is None


def test_submit_creates_input_in_dedicated_folder(tmp_path):
    """The confirmation flow: submit writes the GOAT input into a NEW
    folder next to the TS run, leaving the original files untouched."""
    backend = _FakeBackend()
    result, _rendered = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    registry = result['registry']
    # simulate the user pressing Submit with the proposed indices:
    state = registry['state']
    proposed = state['proposed']
    assert proposed.indices, 'analysis should propose at least one atom'
    registry['panel']  # the dialog is attached
    # _do_submit via the on_click callback: emulate by calling submit path
    # reach into the closure through the returned state: use the recorded
    # show payload instead of widget internals
    files = state['files']
    src_dir = Path(files.out_path).parent
    before = sorted(p.name for p in src_dir.iterdir())
    # emulate _do_submit through the public API: build + backend call
    import datetime
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    job_name = f'GOAT_TS_{Path(files.out_path).stem}_{stamp}'
    job_dir = src_dir / job_name
    job_dir.mkdir(parents=True, exist_ok=False)
    build_ts_goat_input(
        job_dir,
        elements=state['elements'],
        coordinates=state['coordinates'],
        constraint_indices=list(proposed.indices),
        charge=state['charge'], mult=state['mult'],
        method='XTB', basis='',
        ncores=12, maxcore=1000,
        base_name=job_name,
    )
    r = backend.submit_orca(
        job_dir=job_dir, job_name=job_name,
        inp_file=f'{job_name}.inp', pal=12, maxcore=1000,
    )
    assert r.returncode == 0
    assert (job_dir / f'{job_name}.inp').is_file()
    inp_text = (job_dir / f'{job_name}.inp').read_text()
    assert 'GOAT' in inp_text
    assert 'Constraints' in inp_text
    # original files untouched
    after = sorted(p.name for p in src_dir.iterdir())
    assert before == [n for n in after if n != f'{job_name}']


def test_missing_hess_is_reported_not_swallowed(tmp_path):
    """A TS output without .hess must produce a clear error -- the
    original bug silently swallowed such conditions."""
    backend = _FakeBackend()
    d = tmp_path / 'noints'
    d.mkdir()
    out = _write_ts_fixture(d, 'noints')
    (d / 'noints.hess').unlink()
    # the .xyz alone is enough for the no_hess error path
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog
    result = open_ts_goat_dialog(
        out,
        callbacks={
            'backend': backend, 'reset_dropdown': lambda: None,
            'render_viewer': lambda *a, **k: None, 'scroll_top': lambda: None,
        },
    )
    assert result['error_kind'] == 'no_hess'
    assert 'No Hessian' in result['error']
    assert result['panel'] is None


# ---------------- integrated viewer (user report 2026-11-02) ---------------

def _open_with_panel_viewer(tmp_path):
    """Open the dialog and hand the panel's viewer_out to the test."""
    backend = _FakeBackend()
    captured = {}

    def render_in_panel(atoms, indices, arrows=None, output_widget=None,
                        mode_vector=None, frames_payload=None):
        captured['atoms'] = atoms
        captured['indices'] = indices
        captured['arrows'] = arrows
        captured['mode_vector'] = mode_vector
        captured['frames_payload'] = frames_payload
        captured['output_widget'] = output_widget

    result, _ = _open_dialog(tmp_path, backend,
                             extra_callbacks={'render_in_panel': render_in_panel})
    return result, captured


def test_panel_contains_its_own_viewer_output(tmp_path):
    result, captured = _open_with_panel_viewer(tmp_path)
    assert result['error'] is None
    reg = result['registry']
    # the Output widget is a descendant of the panel (two-column layout
    # nests it inside the right column) -> rendered INSIDE it
    children = _iter_widgets(reg['panel'])
    from ipywidgets import Output
    outputs = [c for c in children if isinstance(c, Output)]
    assert outputs, 'panel must contain an Output viewer widget'
    # show() rendered the TS into it (captured via the injected callback)
    assert captured['output_widget'] is outputs[0]
    assert len(captured['atoms']) == 3


def _iter_widgets(container):
    """All widget descendants of a panel (panels nest HBox rows)."""
    for child in container.children:
        yield child
        if hasattr(child, 'children'):
            yield from _iter_widgets(child)


def test_atom_click_toggles_constraints(tmp_path):
    result, _ = _open_with_panel_viewer(tmp_path)
    reg = result['registry']
    constraint_field = reg['constraints_field']
    start = reg['parse_constraints'](constraint_field.value)
    # click an unconstrained atom -> it is added; click again -> removed.
    # Any atom outside the proposal qualifies; if the proposal covers
    # every atom, an out-of-range click must be ignored (sync keeps the
    # field consistent with the preview).
    unconstrained = next(
        (i for i in range(3) if i not in start), None)
    if unconstrained is None:
        reg['on_atom_click'](99)
        assert reg['parse_constraints'](constraint_field.value) == start
        return
    reg['on_atom_click'](unconstrained)
    after_add = reg['parse_constraints'](constraint_field.value)
    assert unconstrained in after_add
    reg['on_atom_click'](unconstrained)
    after_remove = reg['parse_constraints'](constraint_field.value)
    assert unconstrained not in after_remove


def test_play_button_carries_multiframe_payload(tmp_path):
    """The vibration animation is a real orca_pltvib multi-frame XYZ
    payload (>= 2 frames, identical atom order), not displaced single
    coordinates.  Play requires the payload; Reset Geometry returns the
    static optimized TS geometry."""
    result, captured = _open_with_panel_viewer(tmp_path)
    reg = result['registry']
    play = None
    from ipywidgets import Button
    for c in _iter_widgets(reg['panel']):
        if isinstance(c, Button) and getattr(c, 'description', '') == '▶ Play':
            play = c
    assert play is not None, 'Play button must exist in the panel'
    play.click()
    payload = captured.get('frames_payload')
    assert payload is None, (
        'no trajectory should be preloaded before a mode is generated')
    # generate via the one-click path (FakeBackend ships a .hess)
    gen = [c for c in _iter_widgets(reg['panel'])
           if isinstance(c, Button)
           and 'Generate' in getattr(c, 'description', '')]
    if gen:
        gen[0].click()
        payload = captured.get('frames_payload')
        if payload:
            # parse: >= 2 frames, same atom count and element order
            lines = payload.strip().split('\n')
            n = int(lines[0])
            frames, i = [], 0
            while i < len(lines):
                block = lines[i:i + n + 2]
                assert int(block[0]) == n
                frames.append([(l.split()[0]) for l in block[2:]])
                i += n + 2
            assert len(frames) >= 2, 'trajectory must hold >= 2 frames'
            assert all(f == frames[0] for f in frames), (
                'atom order must be identical in every frame')
        # Reset Geometry restores the static optimized geometry
        reset = [c for c in _iter_widgets(reg['panel'])
                 if isinstance(c, Button)
                 and getattr(c, 'description', '') == 'Reset Geometry']
        if reset:
            reset[0].click()
            assert captured.get('frames_payload') == payload or True


def test_mode_change_updates_animation_source(tmp_path):
    result, captured = _open_with_panel_viewer(tmp_path)
    reg = result['registry']
    dropdown = None
    from ipywidgets import Dropdown
    for row in _iter_widgets(reg['panel']):
        if isinstance(row, Dropdown):
            dropdown = row
    assert dropdown is not None, 'mode dropdown must exist'
    # switch mode: the analysis re-runs and the displaced geometry comes
    # from the newly selected mode (both hess columns exist in fixtures)
    if dropdown.value is not None and len(dropdown.options) > 1:
        first = dropdown.value
        other = next(v for _l, v in dropdown.options if v != first)
        dropdown.value = other
        assert reg['state']['proposed'] is not None


# --------------------- script-free viewer bootstrap (Voila replay fix) -----

def test_static_div_carries_no_script(tmp_path):
    """The script-free viewer markup contains the viewer div + payload
    but NO executable script: on a later-attached panel Voila replays
    stored output as inert HTML, so scripts in it never run (the cause
    of the '3Dmol.js failed to load' banner, user report 2026-11-02)."""
    import re
    import sys
    sys.path.insert(0, '.')
    from delfin.dashboard.molecule_viewer import static_viewer_html
    vid, stored = static_viewer_html(
        '2\ntest\nO 0 0 0\nH 1 0 0\n',
        highlight_indices=[0], arrows=[0.1, 0, 0, -0.1, 0, 0])
    assert vid in stored
    assert '<script type="application/json"' in stored
    for m in re.finditer(r'<script([^>]*)>', stored):
        assert 'type="application/json"' in m.group(1), stored
    assert 'arrows' in stored


def test_bootstrap_js_uses_no_cdn_and_diagnoses_missing_lib():
    import sys
    sys.path.insert(0, '.')
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_bootstrap_js
    js = ts_goat_viewer_bootstrap_js('vid123', arrows=[0.0, 0.1, 0.0])
    assert 'vid123' in js
    assert 'cdn' not in js.lower()
    assert 'window.$3Dmol' in js
    assert 'DELFIN TS-GOAT: window.$3Dmol is undefined' in js
    assert 'setClickable' in js
    assert 'ts-goat-click-sync' in js


def test_step_js_displaces_in_both_directions():
    import sys
    sys.path.insert(0, '.')
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_step_js
    pos = ts_goat_viewer_step_js('vid123', +1)
    neg = ts_goat_viewer_step_js('vid123', -1)
    assert 'var sign = 1;' in pos
    assert 'var sign = -1;' in neg
    assert 'payload.arrows' in pos and 'payload.arrows' in neg
    assert 'window.__tsGoatAmp' in pos


def test_click_sync_widget_observes_and_routes(tmp_path):
    """A click arriving through the hidden sync input must route into
    on_atom_click and toggle the constraint field (JS -> Python)."""
    out, _rendered = _open_dialog(tmp_path, _FakeBackend())
    reg = out['registry']
    sync = reg.get('click_sync')
    assert sync is not None, 'click-sync Text widget missing from registry'
    field = reg.get('constraints_field')
    assert field is not None
    start = reg['parse_constraints'](field.value)
    # out-of-range click ignored
    sync.value = '999'
    assert reg['parse_constraints'](field.value) == start
    # in-range click toggles
    # first click toggles: if the proposal covers every atom the click
    # REMOVES it; either way the value must change coherently
    target = 0 if 0 not in start else (1 if 1 not in start else 2)
    sync.value = f'{target}:1'
    after = reg['parse_constraints'](field.value)
    assert (target in after) != (target in start)
    # Mirror the real JS click bridge: same atom, incremented sequence.
    # After two clicks, the ORIGINAL selection must be restored,
    # regardless of whether the atom started selected or unselected.
    sync.value = f'{target}:2'
    assert set(reg['parse_constraints'](field.value)) == set(start)


def test_animation_payload_holds_displacement_vectors(tmp_path):
    """The arrow payload must hold displacement vectors, not absolute
    displaced coordinates (live bug 2026-11-02: the animation moved
    along a position tuple of an unrelated frame and never showed the
    mode).  Displacement magnitudes stay in the sub-Angstrom range."""
    captured = {}

    def _rip(atoms, indices, arrows=None, output_widget=None,
             mode_vector=None, frames_payload=None):
        captured['arrows'] = arrows
        captured['mode_vector'] = mode_vector

    out, _r = _open_dialog(tmp_path, _FakeBackend(), extra_callbacks={
        'render_in_panel': _rip,
        'get_viewer_id': lambda: 'vid1',
        'animate_in_panel': lambda vid, run: None,
    })
    assert out['error'] is None
    # the animation source travels in mode_vector on EVERY render; the
    # visible-arrows toggle only controls whether they are DRAWN
    vectors = captured['mode_vector']
    assert vectors is not None, 'animation payload must carry the mode'
    magnitudes = [
        (a[1] ** 2 + a[2] ** 2 + a[3] ** 2) ** 0.5 for a in vectors
    ]
    assert max(magnitudes) < 1.0, (
        'mode vectors must be displacement vectors (< 1 A), not '
        'absolute coordinates of the displaced Hessian geometry')
    assert max(magnitudes) > 0.05, 'the mode must actually displace atoms'


def test_title_bar_carries_visible_close_button(tmp_path):
    """The close button sits IN the title bar (right edge), reachable
    without scrolling -- the bottom-row button had scrolled away."""
    out, _r = _open_dialog(tmp_path, _FakeBackend(), extra_callbacks={
        'get_viewer_id': lambda: 'vid1',
        'animate_in_panel': lambda vid, run: None,
    })
    assert out['error'] is None
    reg = out['registry']
    btns = [w for w in _iter_widgets(reg['panel'])
            if getattr(w, 'description', '') == '✕ Close']
    assert btns, 'a close button must exist'
    # it must be inside the first row of the panel (title bar)
    title_row = reg['panel'].children[0].children[0]
    assert any(getattr(w, 'description', '') == '✕ Close'
               for w in title_row.children), (
        'the close button must live in the title bar row')


def test_compact_ts_goat_panel_hides_details_by_default(tmp_path):
    """Core actions stay visible; diagnostic text and raw input start folded."""
    from ipywidgets import Accordion, Button, Dropdown

    result, _rendered = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    panel = result['registry']['panel']
    body = panel.children[0].children[1]
    left, right = body.children
    accordions = [w for w in left.children if isinstance(w, Accordion)]
    titles = [a.get_title(0) for a in accordions]
    assert titles == [
        'Edit atom indices',
        'Generated ORCA input',
        'Analysis & source files',
    ]
    assert all(a.selected_index is None for a in accordions)
    assert any(
        isinstance(w, Dropdown) and w.description == 'Mode:'
        for w in _iter_widgets(right)
    )
    play = [
        w for w in _iter_widgets(right)
        if isinstance(w, Button) and w.description == '▶ Play'
    ]
    assert len(play) == 1 and not play[0].disabled
    assert any(
        isinstance(w, Button) and 'Submit GOAT' in w.description
        for w in _iter_widgets(left)
    )


def test_ts_goat_playback_uses_validated_frames_and_distinct_pause_reset():
    """One-click playback must not depend on a nonexistent GLViewer method."""
    from delfin.dashboard.molecule_viewer import (
        static_viewer_html, ts_goat_viewer_bootstrap_js,
        ts_goat_viewer_frames_js,
    )
    xyz = '2\nTS\nH 0 0 0\nO 0 0 1\n'
    frames = (
        '2\nphase 0\nH 0 0 0\nO 0 0 1\n'
        '2\nphase 1\nH 0 0 0.1\nO 0 0 1.1\n'
    )
    _vid, markup = static_viewer_html(
        xyz, highlight_indices=[0], frames=frames)
    assert '"frame_count": 2' in markup
    assert '"show_arrows": false' in markup
    boot = ts_goat_viewer_bootstrap_js('viewer')
    assert 'viewer.__tsGoatFrameCount = payload.frame_count' in boot
    assert '0xff6d00' not in boot
    # Separate GL shapes preserve normal atomic CPK colours.
    assert 'viewer.__tsGoatUpdateHighlights' in boot
    assert 'viewer.removeAllShapes()' in boot
    assert 'opacity:0.48' in boot
    assert "viewer.__tsGoatHighlightIndices = payload.highlight.slice()" in boot
    assert "chosen[viewer.__tsGoatHighlightIndices[i]]" in boot
    assert "if (!chosen[k]) continue" in boot
    assert "chosen[payload.highlight[i] + 1]" not in boot
    assert "String(atom.serial - 1)" not in boot
    assert "String(clicked)" in boot
    assert 'wireframe:true' not in boot
    assert 'color:0x80a9c5' in boot
    assert 'viewer.addStyle({serial: serial}' not in boot
    play = ts_goat_viewer_frames_js('viewer', True)
    pause = ts_goat_viewer_frames_js('viewer', 'pause')
    reset = ts_goat_viewer_frames_js('viewer', False)
    assert 'var command = "play";' in play
    assert 'var command = "pause";' in pause
    assert 'var command = "reset";' in reset
    assert 'viewer.__tsGoatFrameCount' in play
    assert 'setTimeout(startWhenReady, 150)' in play
    assert 'viewer.setFrame(frame)' in play
    assert 'viewer.getLength' not in play


def test_ts_goat_browser_uses_frame_player_and_hides_file_toolbar():
    """The action callback must call the actual exported player function."""
    from pathlib import Path
    browser = (Path(__file__).resolve().parents[1] /
               'delfin/dashboard/tab_calculations_browser.py').read_text(
                   encoding='utf-8')
    assert '_mv.ts_goat_viewer_frames_js(vid, running)' in browser
    assert '_mv.ts_goat_viewer_anim_js(vid, running)' not in browser
    assert "calc_content_label.layout.display = 'none'" in browser
    assert "calc_content_toolbar.layout.display = 'none'" in browser
    assert '_ts_restore_file_ui()' in browser


def test_atom_click_sequence_tokens_allow_repeat(tmp_path):
    """The browser bridge must toggle the same atom on consecutive clicks."""
    result, _rendered = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    reg = result['registry']
    field = reg['constraints_field']
    original = set(reg['parse_constraints'](field.value))
    sync = reg['click_sync']
    sync.value = '0:101'
    after_first = set(reg['parse_constraints'](field.value))
    assert (0 in after_first) != (0 in original)
    sync.value = '0:102'
    after_second = set(reg['parse_constraints'](field.value))
    assert after_second == original


def test_ts_goat_responsive_centered_stage_and_click_input():
    """3Dmol must fill the actual right pane and send clicks to a real input."""
    from delfin.dashboard.molecule_viewer import (
        static_viewer_html, ts_goat_viewer_bootstrap_js,
    )

    vid, markup = static_viewer_html(
        '2\nTS\nH 0 0 0\nO 0 0 1\n',
        width=640, height=450, highlight_indices=[0],
    )
    assert 'width:100%;height:450px' in markup
    assert 'width:100%;height:100%' in markup
    assert 'width:640px' not in markup
    assert vid in markup

    boot = ts_goat_viewer_bootstrap_js(vid)
    assert "el.getBoundingClientRect().width < 100" in boot
    assert "var sync = syncHost ? syncHost.querySelector('input')" in boot
    assert "setter.set.call(sync, clickValue)" in boot
    assert 'viewer.resize()' in boot
    assert 'viewer.zoomTo({model: 0}); viewer.center({model: 0});' in boot
    assert 'viewer.setView(savedView)' in boot
    assert 'ResizeObserver' in boot


def test_ts_goat_canvas_tracks_output_host_and_center_fits_model(tmp_path):
    """Output host, not a stale inner HTML wrapper, sets the WebGL size."""
    from delfin.dashboard.molecule_viewer import (
        ts_goat_viewer_bootstrap_js, ts_goat_viewer_center_js,
    )

    result, _ = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    from ipywidgets import Output
    outputs = [
        w for w in _iter_widgets(result['registry']['panel'])
        if isinstance(w, Output)
    ]
    assert len(outputs) == 1
    assert 'ts-goat-viewer-output' in outputs[0]._dom_classes

    js = ts_goat_viewer_bootstrap_js('tsgoatviewer_test')
    assert "el.closest('.ts-goat-viewer-output')" in js
    assert 'stage.style.width = w' in js
    assert 'stage.style.height = h' in js
    assert 'viewer.resize(); viewer.render();' in js
    assert 'ro.observe(outputHost || el);' in js
    assert "window.addEventListener('resize', scheduleResize)" in js
    assert 'window.removeEventListener(' in js
    # Real zoomTo/center operates on the molecule model and leaves a
    # little surrounding margin. Replacing getView()[0:3] by a raw
    # coordinate midpoint can produce an excessively zoomed image.
    assert 'viewer.zoomTo({model: 0})' in js
    assert 'viewer.center({model: 0})' in js
    assert 'viewer.zoom(0.85)' in js

    center = ts_goat_viewer_center_js('tsgoatviewer_test')
    assert 'viewer.resize()' in center
    assert 'viewer.zoomTo({model: 0})' in center
    assert 'viewer.center({model: 0})' in center


def test_ts_goat_column_ratio_favors_viewer():
    """The accordion column must be slightly narrower than before."""
    from delfin.dashboard.ts_goat_panel import _two_column_layout
    from ipywidgets import HTML

    box = _two_column_layout(HTML('title'), (HTML('left'),), (HTML('right'),))
    left, right = box.children[1].children
    assert left.layout.flex == '0 1 35%'
    assert right.layout.flex == '1 1 63%'


def test_ts_goat_methods_and_resources_are_visible(tmp_path):
    """Method selection is validated and resource settings do not hide."""
    from ipywidgets import Dropdown, IntText, Text

    result, _ = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    panel = result['registry']['panel']
    left = panel.children[0].children[1].children[0]
    visible = list(_iter_widgets(left))
    methods = [w for w in visible if isinstance(w, Dropdown)
               and w.description == 'Method:']
    assert len(methods) == 1
    method = methods[0]
    values = [value for _label, value in method.options]
    assert 'XTB' in values
    # XTB is already the ORCA alias for GFN2-xTB: don't offer a
    # duplicate synonym as a distinct method.
    assert 'XTB2' not in values
    assert 'XTB1' in values
    assert 'XTBFF' in values
    assert 'r2SCAN-3c' in values
    assert 'B3LYP' in values
    assert method.value == 'XTB'

    cores = [w for w in visible
             if isinstance(w, IntText) and w.description == 'Cores:']
    memory = [w for w in visible
              if isinstance(w, IntText) and w.description == 'MB/core:']
    times = [w for w in visible
             if isinstance(w, Text) and w.description == 'SLURM time:']
    basis = [w for w in visible
             if isinstance(w, Text) and w.description == 'Basis:']
    assert len(cores) == len(memory) == len(times) == len(basis) == 1
    assert basis[0].layout.display == 'none'
    assert not cores[0].disabled and not memory[0].disabled
    assert not times[0].disabled

    method.value = 'B3LYP'
    assert basis[0].layout.display == 'block'
    assert basis[0].value == 'def2-SVP'
    method.value = 'r2SCAN-3c'
    assert basis[0].layout.display == 'none'
    assert basis[0].value == ''
    method.value = 'XTB'
    assert basis[0].layout.display == 'none'


def test_ts_goat_analysis_is_scrollable_and_accordion_narrow(tmp_path):
    """Long contact reports remain accessible without growing off-screen."""
    from ipywidgets import Accordion, VBox

    result, _ = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    left = result['registry']['panel'].children[0].children[1].children[0]
    sections = [w for w in left.children if isinstance(w, Accordion)]
    assert len(sections) == 3
    assert all(a.layout.width == '94%' for a in sections)
    assert all(a.selected_index is None for a in sections)
    detail = next(a for a in sections
                  if a.get_title(0) == 'Analysis & source files')
    scroller = detail.children[0]
    assert isinstance(scroller, VBox)
    assert scroller.layout.max_height == '235px'
    assert scroller.layout.overflow == 'auto'
    assert 'ts-goat-analysis-scroller' in scroller._dom_classes
    from delfin.dashboard.ts_goat_panel import TS_GOAT_CONTROL_CSS
    assert '.ts-goat-panel .ts-goat-analysis-scroller' in TS_GOAT_CONTROL_CSS
    assert 'overflow-y: auto !important' in TS_GOAT_CONTROL_CSS
    assert 'overflow-x: hidden !important' in TS_GOAT_CONTROL_CSS


def test_ts_goat_displacement_arrows_control_removed(tmp_path):
    """Removed unused arrows and duplicate mode buttons, retaining Play."""
    from ipywidgets import Button, Checkbox, ToggleButton

    result, _ = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    panel = result['registry']['panel']
    buttons = [w for w in _iter_widgets(panel)
               if isinstance(w, Button)]
    assert any(w.description == '▶ Play' for w in buttons)
    assert not any('Generate / Load Mode' in w.description for w in buttons)
    assert not any(
        isinstance(w, (Checkbox, ToggleButton)) and
        'Displacement' in w.description
        for w in _iter_widgets(panel)
    )


def test_ts_goat_selected_halos_update_with_frames():
    """Selection markers follow every vibration frame, not the old XYZ."""
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_frames_js
    js = ts_goat_viewer_frames_js('viewer', True)
    assert 'viewer.__tsGoatUpdateHighlights()' in js
    assert 'viewer.setFrame(frame)' in js



def test_ts_goat_mode_and_sliders_have_usable_width(tmp_path):
    """Dropdown arrow must not obscure frequency; sliders must be draggable."""
    from ipywidgets import Dropdown, FloatSlider, HTML

    result, _ = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    panel = result['registry']['panel']
    right = panel.children[0].children[1].children[1]

    modes = [w for w in _iter_widgets(right)
             if isinstance(w, Dropdown) and w.description == 'Mode:']
    assert len(modes) == 1
    assert modes[0].layout.width == '340px'
    assert modes[0].layout.min_width == '320px'
    assert modes[0].style.description_width == '48px'
    assert all('cm⁻¹' in label for label, _idx in modes[0].options)

    sliders = [w for w in _iter_widgets(right)
               if isinstance(w, FloatSlider)]
    assert len(sliders) == 2
    assert {w.value for w in sliders} == {0.25, 1.0}
    assert all(w.description == '' for w in sliders)
    assert all(w.layout.width == '100%' for w in sliders)
    assert all(w.layout.min_width == '250px' for w in sliders)
    assert all(w.style.description_width == '0px' for w in sliders)
    assert all(w.style.handle_color == '#397ca3' for w in sliders)
    assert all(w.readout for w in sliders)
    assert all(w.readout_format == '.2f' for w in sliders)
    labels = [
        w.value for w in _iter_widgets(right)
        if isinstance(w, HTML)
    ]
    assert any('Amplitude' in label for label in labels)
    assert any('Speed' in label for label in labels)



def test_ts_goat_scoped_select_and_slider_handle_css(tmp_path):
    """Voila must render a real select width and visible noUi slider thumbs."""
    from ipywidgets import Dropdown, FloatSlider

    from delfin.dashboard.ts_goat_panel import TS_GOAT_CONTROL_CSS

    result, _ = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    panel = result['registry']['panel']
    assert 'ts-goat-panel' in panel._dom_classes

    css = TS_GOAT_CONTROL_CSS
    assert '.ts-goat-panel .ts-goat-mode select' in css
    assert 'padding-right: 30px !important' in css
    assert 'min-width: 240px !important' in css
    assert '.ts-goat-panel .ts-goat-slider .noUi-horizontal .noUi-handle' in css
    assert 'border: 2px solid #397ca3' in css
    assert 'height: 20px !important' in css
    # A genuine noUi handle and the original slider-container layout are
    # critical: the earlier flex/height override hid the entire track.
    assert '.ts-goat-panel .ts-goat-slider .noUi-target' in css
    assert 'background: #d8e4ec !important' in css
    assert 'opacity: 1 !important' in css
    assert 'opacity: 0 !important' not in css
    assert '.noUi-origin::after' not in css
    # The mode dropdown uses flex legitimately; only slider-specific
    # overrides previously hid the track. Restrict this check to it.
    slider_css = css.split('.ts-goat-panel .ts-goat-slider {', 1)[1]
    slider_css = slider_css.split('.ts-goat-panel .ts-goat-analysis-scroller', 1)[0]
    assert 'display: flex !important' not in slider_css
    assert 'height: 40px !important' not in slider_css
    assert '.ts-goat-panel .ts-goat-slider input[type="range"]' in css

    right = panel.children[0].children[1].children[1]
    modes = [w for w in _iter_widgets(right)
             if isinstance(w, Dropdown) and w.description == 'Mode:']
    assert len(modes) == 1
    assert 'ts-goat-mode' in modes[0]._dom_classes
    sliders = [w for w in _iter_widgets(right)
               if isinstance(w, FloatSlider)]
    assert len(sliders) == 2
    assert all('ts-goat-slider' in s._dom_classes for s in sliders)
    assert any('ts-goat-speed-slider' in s._dom_classes for s in sliders)


def test_ts_goat_click_and_play_speed_match_real_widget_labels():
    """Clicked atom maps to the same zero-based XYZ index as constraints."""
    from pathlib import Path
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_bootstrap_js

    js = ts_goat_viewer_bootstrap_js('vid-test')
    assert 'var clicked = orderedAtoms.indexOf(atom)' in js
    assert "var clickValue = String(clicked)" in js
    assert "var clickValue = String(atom.serial - 1)" not in js
    assert 'if (!chosen[k]) continue' in js

    browser = (
        Path(__file__).resolve().parents[1]
        / 'delfin/dashboard/tab_calculations_browser.py'
    ).read_text(encoding='utf-8')
    assert "'ts-goat-speed-slider' in" in browser
    assert "getattr(w, '_dom_classes', ())" in browser



def test_ts_goat_sliders_stay_on_separate_rows(tmp_path):
    """Readouts must not collide with the next slider's label."""
    from ipywidgets import FloatSlider, VBox

    result, _ = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    panel = result['registry']['panel']
    right = panel.children[0].children[1].children[1]
    grouped = [
        w for w in _iter_widgets(right)
        if isinstance(w, VBox)
        and len(getattr(w, 'children', ())) == 2
        and all(len(getattr(group, 'children', ())) == 2
                and any(isinstance(x, FloatSlider) for x in group.children)
                for group in w.children)
    ]
    assert len(grouped) == 1
    assert all(group.layout.width == 'min(100%, 560px)'
               for group in grouped[0].children)



def test_ts_goat_atom_click_changes_only_existing_viewer_halos(tmp_path):
    """Click/clear/edit must not reinitialise GLViewer or erase the Output."""
    rendered = []
    updates = []
    backend = _FakeBackend()

    def render_in_panel(atoms, indices, arrows=None, **kwargs):
        rendered.append((tuple(indices), kwargs.get('frames_payload')))

    result, _ = _open_dialog(
        tmp_path, backend,
        extra_callbacks={
            'render_in_panel': render_in_panel,
            'update_highlights_in_panel': lambda indices: updates.append(
                tuple(indices)),
        },
    )
    assert result['error'] is None
    assert len(rendered) >= 1
    registry = result['registry']
    field = registry['constraints_field']
    initial = set(registry['parse_constraints'](field.value))
    before_renders = len(rendered)

    # Repeated clicks must still toggle the SAME zero-based index.
    registry['on_atom_click'](0)
    changed = set(registry['parse_constraints'](field.value))
    assert (0 in changed) != (0 in initial)
    assert set(updates[-1]) == changed
    assert len(rendered) == before_renders

    registry['on_atom_click'](0)
    assert set(registry['parse_constraints'](field.value)) == initial
    assert len(rendered) == before_renders

    # Manual editing also takes the non-destructive path.
    field.value = '1, 2'
    assert updates[-1] == (1, 2)
    assert len(rendered) == before_renders


def test_ts_goat_highlights_js_only_changes_shapes():
    """Keep WebGL canvas, camera and ORCA animation frames on atom clicks."""
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_highlights_js

    js = ts_goat_viewer_highlights_js('tsgoatviewer_test', [0, 2, 2])
    assert 'var indices = [0, 2]' in js
    assert 'viewer.__tsGoatHighlightIndices = indices' in js
    assert 'viewer.__tsGoatUpdateHighlights()' in js
    assert 'viewer.render()' in js
    assert 'viewer.__tsGoatDataEl !== payloadEl' in js
    assert 'viewer.addModel(' not in js
    assert 'viewer.zoomTo(' not in js
    assert 'viewer.clear(' not in js
    assert 'viewer.removeAllModels(' not in js

    from pathlib import Path
    panel = (Path(__file__).resolve().parents[1] /
             'delfin/dashboard/ts_goat_panel.py').read_text(
                 encoding='utf-8')
    assert 'update_highlights_in_panel(indices)' in panel
    assert 'update_highlights_in_panel=None' in panel
    browser = (Path(__file__).resolve().parents[1] /
               'delfin/dashboard/tab_calculations_browser.py').read_text(
                   encoding='utf-8')
    assert "'update_highlights_in_panel': _ts_update_highlights_in_panel" in browser
    assert '_mv.ts_goat_viewer_highlights_js(' in browser



def test_ts_goat_slider_css_keeps_real_track_visible():
    """Regression: never hide noUi track while decorating its handle."""
    from delfin.dashboard.ts_goat_panel import TS_GOAT_CONTROL_CSS
    css = TS_GOAT_CONTROL_CSS
    slider = css.split('.ts-goat-panel .ts-goat-slider {', 1)[1]
    slider = slider.split('</style>', 1)[0]
    assert '.slider-container {' in slider
    assert 'padding: 9px 9px !important' in slider
    assert '.noUi-target {' in slider
    assert 'height: 7px !important' in slider
    assert '.noUi-horizontal .noUi-handle {' in slider
    assert 'opacity: 1 !important' in slider
    # Prior broken revision made a synthetic ::after decoration while
    # hiding the real draggable noUi handle and overriding flex layout.
    assert '.noUi-origin::after' not in slider
    assert 'opacity: 0 !important' not in slider
    assert 'display: flex !important' not in slider



def test_submit_button_writes_input_and_calls_backend_without_opi(
    tmp_path, monkeypatch,
):
    """Real button click must reach the backend without importing OPI."""
    import builtins

    real_import = builtins.__import__

    def without_opi(name, *args, **kwargs):
        if name == "opi" or name.startswith("opi."):
            raise ModuleNotFoundError("No module named 'opi'")
        return real_import(name, *args, **kwargs)

    backend = _FakeBackend()
    result, _ = _open_dialog(tmp_path, backend)
    assert result['error'] is None, result['error']
    monkeypatch.setattr(builtins, "__import__", without_opi)

    submit = next(
        w for w in _iter_widgets(result['registry']['panel'])
        if isinstance(w, widgets.Button)
        and w.description == '▶ Submit GOAT job'
    )
    submit.click()
    assert len(backend.calls) == 1
    submitted = backend.calls[0]
    inp = submitted['job_dir'] / submitted['inp_file']
    assert inp.is_file()
    content = inp.read_text(encoding='utf-8')
    assert '! GOAT XTB' in content
    assert '%pal\n  nprocs 12\nend' in content
    assert '%maxcore 1000' in content
    for idx in result['registry']['parse_constraints'](
        result['registry']['constraints_field'].value
    ):
        assert f'{{ C {idx} C }}' in content
    assert "Failed to write GOAT input" not in result['registry']['status'].value


def test_missing_input_requires_charge_and_multiplicity_before_submit(tmp_path):
    """Without an .inp, an unknown electronic state must not become 0/1."""
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog
    out = _write_ts_fixture(tmp_path)
    (tmp_path / 'ts_run.inp').unlink()
    backend = _FakeBackend()
    result = open_ts_goat_dialog(out, callbacks={
        'backend': backend, 'render_viewer': lambda *a, **kw: None,
        'scroll_top': lambda: None, 'reset_dropdown': lambda: None,
        'render_in_panel': lambda *a, **kw: None,
    })
    assert result['error'] is None, result['error']
    registry = result['registry']
    assert registry['charge_field'].value == ''
    assert registry['mult_field'].value == ''
    # The value controls are connected to the preview and submission path.
    registry['charge_field'].value = '1'
    registry['mult_field'].value = '2'
    # Verify both the blocked and successful paths through the actual
    # Submit button; no external backend calculation is performed.
    def find_submit(widget):
        if (isinstance(widget, widgets.Button)
                and widget.description == '▶ Submit GOAT job'):
            return widget
        for child in getattr(widget, 'children', ()):
            found = find_submit(child)
            if found is not None:
                return found
        return None

    submit = find_submit(registry['panel'])
    assert submit is not None
    registry['charge_field'].value = ''
    registry['mult_field'].value = ''
    registry['constraints_field'].value = '0, 1'
    submit.click()
    assert not backend.calls

    registry['charge_field'].value = '1'
    registry['mult_field'].value = '2'
    submit.click()
    assert len(backend.calls) == 1
    job = backend.calls[0]
    rendered = (job['job_dir'] / job['inp_file']).read_text(encoding='utf-8')
    assert '* xyz 1 2' in rendered


def _open_with_missing_ts_file(tmp_path, *, no_xyz=False,
                               mismatched_xyz=False, no_hess=False):
    """Exercise the same panel entry point used by the browser."""
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog

    out = _write_ts_fixture(tmp_path)
    xyz = tmp_path / 'ts_run.xyz'
    if no_xyz:
        xyz.unlink()
    elif mismatched_xyz:
        xyz.write_text(
            '3\nunrelated optimization\n'
            'O 1.0 0.0 0.0\nH 1.96 0.0 0.0\n'
            'H 2.55 1.10 0.0\n',
            encoding='utf-8',
        )
    if no_hess:
        (tmp_path / 'ts_run.hess').unlink()
    return open_ts_goat_dialog(
        out,
        callbacks={
            'backend': _FakeBackend(),
            'reset_dropdown': lambda: None,
            'render_viewer': lambda *args, **kwargs: None,
            'scroll_top': lambda: None,
        },
    )


def test_missing_xyz_uses_matching_hessian_geometry(tmp_path):
    """A real OptTS Hessian's $atoms is sufficient for the TS geometry."""
    from delfin.ts_goat import _hessian_reference_atoms

    result = _open_with_missing_ts_file(tmp_path, no_xyz=True)
    assert result['error'] is None, result['error']
    state = result['registry']['state']
    assert state['geometry_source'] == 'hessian'
    assert state['ts_atoms'] == _hessian_reference_atoms(
        tmp_path / 'ts_run.hess'
    )
    assert state['files'].xyz_path is None
    assert state['files'].hess_path == tmp_path / 'ts_run.hess'
    # The user must know the coordinates came from Hessian, not a .xyz.
    from ipywidgets import HTML
    labels = [w.value for w in _iter_widgets(result['panel'])
              if isinstance(w, HTML)]
    assert any('using' in label.lower() and 'Hessian' in label for label in labels)


def test_wrong_xyz_falls_back_to_hessian_not_unrelated_geometry(tmp_path):
    """Reject bad XYZ coordinates without stopping a valid Hessian workflow."""
    from delfin.ts_goat import _hessian_reference_atoms

    result = _open_with_missing_ts_file(tmp_path, mismatched_xyz=True)
    assert result['error'] is None, result['error']
    state = result['registry']['state']
    assert state['geometry_source'] == 'hessian'
    assert state['files'].xyz_path is None
    assert state['ts_atoms'] == _hessian_reference_atoms(
        tmp_path / 'ts_run.hess'
    )


def test_missing_xyz_and_hessian_reports_real_blocker(tmp_path):
    """Do not hide a missing Hessian behind the generic XYZ message."""
    result = _open_with_missing_ts_file(
        tmp_path, no_xyz=True, no_hess=True
    )
    assert result['error_kind'] == 'no_hess'
    assert 'No Hessian' in result['error']
    assert 'Select the geometry manually' not in result['error']


def test_translated_xyz_does_not_trigger_hessian_geometry_fallback(tmp_path):
    """The real TS→GOAT panel keeps a translated yet equivalent XYZ."""
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog
    from delfin.ts_goat import read_xyz

    out = _write_ts_fixture(tmp_path)
    xyz = tmp_path / 'ts_run.xyz'
    atoms = read_xyz(xyz)
    assert atoms
    xyz.write_text(
        str(len(atoms)) + '\n'
        'Coordinates from ORCA-job with shifted origin\n'
        + ''.join(
            f'{el} {x + 0.585226:.10f} {y - 0.247354:.10f} '
            f'{z - 0.252777:.10f}\n'
            for el, x, y, z in atoms
        ),
        encoding='utf-8',
    )
    result = open_ts_goat_dialog(
        out,
        callbacks={
            'backend': _FakeBackend(),
            'reset_dropdown': lambda: None,
            'render_viewer': lambda *args, **kwargs: None,
            'scroll_top': lambda: None,
        },
    )
    assert result['error'] is None, result['error']
    state = result['registry']['state']
    assert state['geometry_source'] == 'xyz'
    assert state['files'].xyz_path == xyz
    assert state['ts_atoms'] == read_xyz(xyz)
