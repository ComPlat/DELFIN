"""TS -> GOAT preparation panel for the Calculations Browser.

Widget-only module: it assembles the preparation dialog (status,
constraint chips, method/settings editors, GOAT preview, submit) around
the browser's existing 3Dmol viewer. The library half lives in
``delfin.ts_goat``; this module only draws and collects user decisions.
"""

from __future__ import annotations

import html as html
from pathlib import Path

import ipywidgets as widgets
from IPython.display import HTML, display


def _two_column_layout(title_row, left, right):
    """Title bar + left (info/settings) | right (3D viewer) columns.

    On narrow browser windows ipywidgets' FlexBox wraps the columns
    below each other, so no separate mobile layout is needed.
    """
    body = widgets.HBox(
        [widgets.VBox(tuple(left)), widgets.VBox(tuple(right))],
        layout=widgets.Layout(
            display='flex', flex_flow='row wrap', width='100%',
        ),
    )
    left_box, right_box = body.children
    left_box.layout = widgets.Layout(flex='1 1 46%', min_width='340px')
    right_box.layout = widgets.Layout(flex='1 1 54%', min_width='380px')
    return widgets.VBox((title_row, body))


def build_ts_goat_panel(
    *,
    get_selected_path,
    backend,
    reset_dropdown,
    render_viewer,
    scroll_top,
    viewer_output=None,
    render_in_panel=None,
    on_close=None,
    center_viewer=None,
    get_viewer_id=None,
    animate_in_panel=None,
) -> dict:
    """Assemble the TS -> GOAT panel; returns a widget registry dict.

    Every callback is injected so this module stays free of browser
    imports (the tab passes its own closures in).  ``render_viewer(xyz,
    highlight_indices, arrows)`` re-renders the 3D viewer with the
    constraint atoms highlighted and mode arrows; the highlighter is
    implemented in the tab (see ``_render_3dmol`` extension there).
    """
    panel = widgets.VBox(layout=widgets.Layout(
        display='none', border='1px solid #90caf9',
        padding='8px', margin='4px 0',
        # The browser's right column has overflow_y='hidden' and the
        # toolbar sits above this panel in the DOM; an unbounded panel
        # grew past the container and clipped it (user report 2026-11-02).
        # Bound the panel and scroll INSIDE it instead.
        max_height='720px', overflow_y='auto',
    ))
    title = widgets.HTML(
        '<b>TS &rarr; GOAT</b> &mdash; constrained conformer search from the '
        'optimized transition state'
    )
    # 2026-11-02 live report: the close button sat in the bottom action
    # row and scrolled out of sight.  The title bar now carries its own
    # close button at the right edge, always visible without scrolling.
    title_close_btn = widgets.Button(
        description='✕ Schließen', button_style='warning',
        tooltip='Close TS→GOAT and restore the normal file view',
        layout=widgets.Layout(width='120px', margin='0 0 0 auto'),
    )
    status = widgets.HTML('')
    files_box = widgets.HTML('')
    analysis_box = widgets.HTML('')
    # Integrated 3D viewer: an Output widget INSIDE the panel, rendered
    # by molecule_viewer.render_xyz_in_output with the tab's click JS.
    # The structure is always visible when the panel opens (user report
    # 2026-11-02: the dialog only showed file info and fields).
    viewer_out = widgets.Output(layout=widgets.Layout(
        width='100%', height='470px', border='1px solid #bbb',
        overflow='auto', margin='4px 0',
    ))
    # JS -> Python click channel (same pattern as the structure editor's
    # submit-pick-sync): the viewer bootstrap writes the clicked atom's
    # 0-based index into this hidden Text via the property setter +
    # input event; the observe() below routes it to on_atom_click.
    click_sync = widgets.Text(value='', layout=widgets.Layout(display='none'))
    click_sync.add_class('ts-goat-click-sync')
    mode_label = widgets.HTML('<b>Reactive mode (imaginary modes of the final block):</b>')
    mode_dropdown = widgets.Dropdown(
        options=[], description='Mode:',
        layout=widgets.Layout(width='320px'),
    )
    mode_warning = widgets.HTML('')
    # Imaginary-mode animation controls (user request 2026-11-02): the
    # selected mode displaces the structure along +/- eigenvector with
    # adjustable amplitude; arrows show the displacement directions.
    anim_toggle = widgets.ToggleButton(
        value=False, description='Imag. Mode anzeigen',
        button_style='warning',
        tooltip='Displace the TS along the selected imaginary mode',
        layout=widgets.Layout(width='190px'),
    )
    anim_play = widgets.Button(
        description='▶ Play', disabled=True,
        layout=widgets.Layout(width='90px'),
    )
    anim_pause = widgets.Button(
        description='⏸ Pause', disabled=True,
        layout=widgets.Layout(width='90px'),
    )
    anim_amp = widgets.FloatSlider(
        value=0.25, min=0.05, max=1.0, step=0.05,
        description='Amplitude:', continuous_update=False,
        layout=widgets.Layout(width='200px'),
    )
    anim_speed = widgets.FloatSlider(
        value=1.0, min=0.25, max=3.0, step=0.25,
        description='Geschw.:', continuous_update=False,
        tooltip='Animation speed (1 = 1.6 s per period)',
        layout=widgets.Layout(width='180px'),
    )
    anim_arrows_toggle = widgets.Checkbox(
        value=False, description='Schwingungspfeile anzeigen',
        layout=widgets.Layout(width='240px'),
    )
    anim_row = widgets.HBox(
        [anim_toggle, anim_play, anim_pause, anim_amp, anim_speed,
         anim_arrows_toggle],
        layout=widgets.Layout(margin='2px 0', flex_flow='row wrap',
                              width='100%'),
    )
    constraints_label = widgets.HTML('<b>Constraint atoms (0-based ORCA indices):</b>')
    constraints_field = widgets.Text(
        placeholder='e.g. 0, 2, 5', description='',
        layout=widgets.Layout(width='260px'),
    )
    highlight_btn = widgets.Button(
        description='↻ Highlight', button_style='info',
        tooltip='Re-render the viewer with the current constraint atoms highlighted',
        layout=widgets.Layout(width='140px'),
    )
    reset_prop_btn = widgets.Button(
        description='Auswahl zurücksetzen',
        tooltip='Restore the automatic reactive-centre proposal',
        layout=widgets.Layout(width='170px'),
    )
    clear_sel_btn = widgets.Button(
        description='Alle Markierungen entfernen',
        tooltip='Clear every constraint (fully manual selection)',
        layout=widgets.Layout(width='200px'),
    )
    method_label = widgets.HTML('<b>Method for the GOAT search:</b>')
    method_field = widgets.Text(
        value='XTB', description='Method:', layout=widgets.Layout(width='240px'),
    )
    basis_field = widgets.Text(
        value='', placeholder='def2-SVP (DFT only)',
        description='Basis:', layout=widgets.Layout(width='240px'),
    )
    pal_field = widgets.IntText(
        value=12, description='Cores:', layout=widgets.Layout(width='200px'),
    )
    maxcore_field = widgets.IntText(
        value=1000, description='MB/core:', layout=widgets.Layout(width='200px'),
    )
    slurm_time_field = widgets.Text(
        value='48:00:00', description='SLURM time:',
        layout=widgets.Layout(width='200px'),
    )
    preview_label = widgets.HTML('<b>ORCA input preview:</b>')
    preview_box = widgets.HTML('')
    submit_btn = widgets.Button(
        description='▶ Submit GOAT job', button_style='success',
        tooltip='Write the input into a new folder and submit it',
        layout=widgets.Layout(width='200px'),
    )
    close_btn = widgets.Button(
        description='✕ Schließen', layout=widgets.Layout(width='110px'),
    )
    center_btn = widgets.Button(
        description='Ansicht zentrieren',
        tooltip='Re-center and re-zoom the 3D view',
        layout=widgets.Layout(width='150px'),
    )
    cancel_btn = widgets.Button(
        description='Abbrechen',
        tooltip='Close TS→GOAT and restore the normal file view',
        layout=widgets.Layout(width='110px'),
    )
    btn_row = widgets.HBox([submit_btn, close_btn])
    cancel_row = widgets.HBox([cancel_btn])
    settings_row = widgets.HBox([
        widgets.VBox([method_label, method_field, basis_field]),
        widgets.VBox([pal_field, maxcore_field, slurm_time_field]),
    ])
    viewer_head = widgets.HTML(
        '<b>Optimized TS structure</b> (drag to rotate, scroll to zoom; '
        'constraint atoms are highlighted):')
    center_row = widgets.HBox([center_btn])

    layout = _two_column_layout(
        title_row=widgets.HBox([title, title_close_btn],
            layout=widgets.Layout(display='flex', width='100%',
                                   align_items='center')),
        left=(status, files_box, mode_label, mode_dropdown, mode_warning,
              analysis_box, constraints_label,
              widgets.HBox([constraints_field, highlight_btn, reset_prop_btn,
                            clear_sel_btn],
                           layout=widgets.Layout(flex_flow='row wrap',
                                                 width='100%')),
              settings_row, preview_label, preview_box, btn_row,
              cancel_row),
        right=(viewer_head, viewer_out, click_sync, anim_row, center_row),
    )
    panel.children = (layout,)

    # ---- state ----------------------------------------------------------
    state = {
        'files': None, 'ts_atoms': None, 'proposed': None,
        'elements': (), 'coordinates': (), 'charge': 0, 'mult': 1,
        'mode_arrows': None,
        'mode_arrows_by_index': {}, 'imaginary_modes': [],
        'hess_path': None,
    }

    def _esc(text):
        return html.escape(str(text), quote=True)

    def _parse_constraints(text):
        indices = []
        for token in (text or '').replace(';', ',').split(','):
            token = token.strip()
            if not token:
                continue
            try:
                idx = int(token)
            except ValueError:
                continue
            if idx >= 0:
                indices.append(idx)
        return indices

    def _refresh_preview():
        indices = _parse_constraints(constraints_field.value)
        atoms = state['ts_atoms'] or []
        n_atoms = len(atoms)
        bad = [i for i in indices if i >= n_atoms]
        preview_box.value = (
            '<pre style="font-size:11px;">' + _esc(
                render_constraints_preview(
                    state['elements'], state['coordinates'], indices,
                    method=method_field.value,
                    basis=basis_field.value,
                    charge=state['charge'], mult=state['mult'],
                    pal=pal_field.value, maxcore=maxcore_field.value,
                )
            ) + '</pre>'
            + (
                f'<div style="color:#c62828;">Out-of-range indices: '
                f'{_esc(bad)} (molecule has {n_atoms} atoms)</div>'
                if bad else ''
            )
        )

    def _refresh_analysis_display():
        files = state['files']
        if files is None:
            return
        rows = []
        for label, path in (
            ('Output', files.out_path), ('Input', files.inp_path),
            ('TS geometry', files.xyz_path), ('Hessian', files.hess_path),
        ):
            rows.append(
                f'<li><b>{label}:</b> '
                + (_esc(Path(path).name) if path else
                   '<span style="color:#c62828;">not found</span>')
                + '</li>'
            )
        problems = ''
        if files.problems:
            problems = '<ul style="color:#c62828; font-size:12px;">' + ''.join(
                f'<li>{_esc(p)}</li>' for p in files.problems
            ) + '</ul>'
        files_box.value = '<ul style="font-size:12px; margin:4px 0;">' + ''.join(rows) + '</ul>' + problems

    def _refresh_constraints_display():
        constraints_field.value = ', '.join(
            str(i) for i in (state['proposed'].indices if state['proposed'] else ())
        )
        _refresh_preview()
        _refresh_highlight()

    def _refresh_highlight():
        indices = _parse_constraints(constraints_field.value)
        if state['ts_atoms'] is None:
            return
        arrows = state.get('mode_arrows')
        render_viewer(state['ts_atoms'], indices, arrows)

    # ---- integrated 3D viewer -------------------------------------------
    # Rendered INSIDE this panel through the tab's Output-callbacks so
    # the structure is visible when the panel opens (user report
    # 2026-11-02).  The tab passes a dedicated viewer_output widget and
    # click handler; without them the panel falls back to the tab
    # viewer and shows a note instead of a blank area.

    def _viewer_fallback(message):
        with viewer_out:
            viewer_out.clear_output(wait=True)
            display(HTML(
                f'<div style="padding:8px; color:#616161;">{html.escape(message)}</div>'
            ))

    def _render_integrated_viewer():
        """Render TS structure + highlight + arrows INTO the panel.

        ``render_in_panel(atoms, indices, arrows)`` is injected by the
        tab and draws into ``viewer_output`` (an Output widget that is
        a child of this panel).  Without it the panel shows a visible
        note instead of a silent blank area.
        """
        atoms = state['ts_atoms']
        if atoms is None:
            _viewer_fallback('No TS geometry loaded.')
            return
        if render_in_panel is None:
            _viewer_fallback(
                '3D viewer not wired by the tab -- file info and '
                'constraints still work below.')
            return
        indices = _parse_constraints(constraints_field.value)
        arrows = state.get('mode_arrows')
        if anim_toggle.value and state.get('mode_displaced'):
            displaced = state['mode_displaced'].get(
                mode_dropdown.value if mode_dropdown.value is not None
                else -1)
            if displaced:
                atoms = displaced
        try:
            # the animation needs the displacement vectors in EVERY
            # render (also without visible arrows); the arrows toggle
            # only controls whether they are DRAWN.  Before this split
            # a default-off arrows checkbox starved the animation of
            # its payload -- nothing moved (live report 2026-11-02).
            render_in_panel(
                atoms, indices,
                arrows if anim_arrows_toggle.value else None,
                output_widget=viewer_out,
                mode_vector=state.get('mode_arrows'),
            )
        except Exception as exc:  # visible error, never a blank area
            _viewer_fallback(
                f'Viewer could not be rendered: {type(exc).__name__}: {exc}')

    def _on_anim_toggle(_change=None):
        anim_active = anim_toggle.value
        anim_play.disabled = not anim_active
        anim_pause.disabled = not anim_active
        if not anim_active:
            # back to the undistorted optimized geometry; constraint
            # marks stay visible
            state['mode_displaced'] = None
        else:
            _recompute_displaced()
        _render_integrated_viewer()

    def _recompute_displaced():
        hess = state.get('hess_path')
        idx = mode_dropdown.value
        if hess is None or idx is None:
            state['mode_displaced'] = None
            return
        try:
            from delfin.imag import displaced_geometries
            disp = displaced_geometries(hess, idx, anim_amp.value)
            state['mode_displaced'] = {idx: disp['pos']}
            state['mode_displaced_neg'] = {idx: disp['neg']}
        except Exception:
            state['mode_displaced'] = None

    def _on_anim_amp(_change=None):
        if anim_toggle.value:
            _recompute_displaced()
            _render_integrated_viewer()

    anim_toggle.observe(_on_anim_toggle, names='value')
    anim_amp.observe(_on_anim_amp, names='value')
    anim_arrows_toggle.observe(
        lambda _c: _render_integrated_viewer(), names='value')

    def _on_speed_change(_change=None):
        """Write the period for the browser loop (base 1600 ms)."""
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            animate_in_panel(vid, 'speed')

    def _on_play(_b=None):
        if not anim_toggle.value:
            return
        _recompute_displaced()
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            # periodical sinus animation entirely in the browser; the
            # static render keeps its optimised geometry meanwhile
            animate_in_panel(vid, True)
        else:
            _render_integrated_viewer()

    def _on_pause(_b=None):
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            animate_in_panel(vid, False)
        # back to the undistorted geometry, marks stay
        state['mode_displaced'] = None
        _render_integrated_viewer()

    anim_play.on_click(_on_play)
    anim_pause.on_click(_on_pause)
    anim_speed.observe(_on_speed_change, names='value')

    # ---- atom click -> constraint toggle --------------------------------
    # The tab injects click JS; when the user clicks an atom, the
    # viewer reports the 0-based index through this callback and the
    # constraint field updates immediately (sync with preview/submit).

    def on_atom_click(index):
        # an out-of-range click (e.g. a stale 3Dmol event) must not
        # corrupt the constraint field
        n_atoms = len(state['ts_atoms'] or [])
        if not 0 <= index < n_atoms:
            return
        indices = _parse_constraints(constraints_field.value)
        if index in indices:
            indices = [i for i in indices if i != index]
        else:
            indices = sorted(indices + [index])
        constraints_field.value = ', '.join(str(i) for i in indices)
        _refresh_preview()
        _render_integrated_viewer()

    def _on_mode_change(change=None):
        """Re-run the reactive-centre analysis for the selected mode."""
        if change is not None and change.get('name') != 'value':
            return
        idx = mode_dropdown.value
        if idx is None or state['ts_atoms'] is None or not state['hess_path']:
            return
        state['mode_arrows'] = state['mode_arrows_by_index'].get(idx)
        from delfin import ts_goat as _tg
        # TS geometry + selected mode only -- a plain OptTS run has no
        # reactant/product files and the analysis must not need them.
        proposed = _tg.analyze_reactive_center(
            state['hess_path'], state['ts_atoms'], mode_index=idx,
        )
        state['proposed'] = proposed
        analysis_box.value = ''.join(
            f'<li>{html.escape(r)}</li>' for r in proposed.reasons
        )
        # reliability = quality/consistency of the data; ambiguity =
        # uniqueness of the PROPOSAL.  These are different things and
        # must not read as one verdict (live report 2026-11-02).
        status.value = (
            f'<span style="font-size:12px;">mode {idx}: '
            f'data quality {proposed.reliability}'
            + (' -- <b>proposal AMBIGUOUS:</b> several contacts answer '
               'the mode similarly strongly; review and adjust the '
               'selection in the viewer'
               if proposed.ambiguous else '')
            + '</span>'
        )
        _refresh_constraints_display()

    mode_dropdown.observe(_on_mode_change, names='value')

    highlight_btn.on_click(lambda _b: _refresh_constraints_display())

    def show(open_result):
        """Open the panel with a prepared analysis (from the tab handler)."""
        state['files'] = open_result['files']
        state['ts_atoms'] = open_result['ts_atoms']
        state['proposed'] = open_result['proposed']
        state['elements'] = open_result['elements']
        state['coordinates'] = open_result['coordinates']
        state['charge'] = open_result.get('charge', 0)
        state['mult'] = open_result.get('mult', 1)
        state['mode_arrows'] = open_result.get('mode_arrows')
        state['mode_arrows_by_index'] = open_result.get(
            'mode_arrows_by_index', {})
        state['imaginary_modes'] = open_result.get('imaginary_modes', [])
        state['hess_path'] = open_result.get('hess_path')
        state['mode_displaced'] = None
        anim_toggle.value = False
        anim_play.disabled = True
        anim_pause.disabled = True
        state['reactant_path'] = open_result.get('reactant_path')
        state['product_path'] = open_result.get('product_path')
        mode_opts = [
            (f'mode {i} ({f:.1f} cm-1)', i)
            for i, f in state['imaginary_modes']
        ]
        mode_dropdown.options = mode_opts or [('(no imaginary mode)', None)]
        mode_dropdown.value = (
            mode_opts[0][1] if mode_opts else None
        )
        # single mode: the picker is informational; multiple: the user picks
        mode_label.layout.display = 'block' if len(mode_opts) > 1 else 'none'
        mode_dropdown.layout.display = 'block' if len(mode_opts) > 1 else 'none'
        mode_warning.value = (
            '<div style="color:#e65100; font-size:12px;">'
            '⚠ Several imaginary modes in the final frequency block -- '
            'select the reactive mode; the constraints are re-proposed '
            'for the selected mode.</div>'
        ) if len(mode_opts) > 1 else ''
        status.value = open_result.get('status_html', '')
        _refresh_analysis_display()
        analysis_box.value = open_result.get('analysis_html', '')
        _refresh_constraints_display()
        _render_integrated_viewer()
        panel.layout.display = 'block'
        scroll_top()

    def close(_button=None):
        panel.layout.display = 'none'
        # restore the normal file view: the options dropdown goes back
        # to '(Options)' so TS→GOAT can be re-opened cleanly, and the
        # regular content area (file preview) becomes visible again.
        try:
            reset_dropdown()
        except Exception:
            pass
        if on_close:
            try:
                on_close()
            except Exception:
                pass

    close_btn.on_click(close)
    cancel_btn.on_click(close)

    def on_submit(_button):
        indices = _parse_constraints(constraints_field.value)
        if not indices:
            status.value = (
                '<span style="color:#c62828;">Select at least one constraint '
                'atom before submitting.</span>'
            )
            return
        atoms = state['ts_atoms'] or []
        bad = [i for i in indices if not (0 <= i < len(atoms))]
        if bad:
            status.value = (
                f'<span style="color:#c62828;">Invalid constraint indices '
                f'{_esc(bad)}: the molecule has {len(atoms)} atoms.</span>'
            )
            return
        _do_submit(indices)

    def _do_submit(indices):
        """Create a dedicated folder and submit through the existing backend."""
        import datetime

        files = state['files']
        source_dir = Path(files.out_path).parent
        stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        job_name = f'GOAT_TS_{Path(files.out_path).stem}_{stamp}'
        job_dir = source_dir / job_name
        job_dir.mkdir(parents=True, exist_ok=True)

        inp_file = f'{job_name}.inp'
        try:
            from delfin.ts_goat import build_ts_goat_input
            build_ts_goat_input(
                job_dir,
                elements=state['elements'],
                coordinates=state['coordinates'],
                constraint_indices=indices,
                charge=state['charge'], mult=state['mult'],
                method=method_field.value, basis=basis_field.value,
                ncores=pal_field.value, maxcore=maxcore_field.value,
                base_name=job_name,
            )
        except Exception as exc:
            status.value = (
                f'<span style="color:#c62828;">Failed to write GOAT input: '
                f'{_esc(exc)}</span>'
            )
            return

        try:
            result = backend.submit_orca(
                job_dir=job_dir,
                job_name=job_name,
                inp_file=inp_file,
                time_limit=slurm_time_field.value,
                pal=pal_field.value,
                maxcore=maxcore_field.value,
            )
        except Exception as exc:
            status.value = (
                f'<span style="color:#c62828;">Submit error: {_esc(exc)}</span>'
            )
            return

        if result.returncode == 0:
            job_id = (
                result.stdout.strip().split()[-1]
                if result.stdout.strip() else '(unknown)'
            )
            status.value = (
                f'<span style="color:#2e7d32;">Submitted '
                f'<code>{_esc(job_name)}</code> (ID: {_esc(job_id)}) '
                f'&mdash; input in {_esc(job_dir.name)}/</span>'
            )
        else:
            msg = result.stderr or result.stdout or 'Unknown error'
            status.value = (
                f'<span style="color:#c62828;">{_esc(msg)}</span>'
            )

    submit_btn.on_click(on_submit)

    def _on_reset_proposal(_b):
        proposed = state.get('proposed')
        indices = list(getattr(proposed, 'indices', None) or [])
        constraints_field.value = ', '.join(str(i) for i in indices)

    def _on_clear_selection(_b):
        constraints_field.value = ''

    reset_prop_btn.on_click(_on_reset_proposal)
    clear_sel_btn.on_click(_on_clear_selection)

    constraints_field.observe(lambda _c: (_refresh_preview(),
                                          _render_integrated_viewer()),
                              names='value')

    def _on_click_sync(change):
        raw = (change.get('new') or '').strip()
        if not raw:
            return
        try:
            on_atom_click(int(raw))
        except ValueError:
            pass
    click_sync.observe(_on_click_sync, names='value')

    # the title-bar close button and the bottom buttons share ONE close
    # path (no divergent callbacks)
    title_close_btn.on_click(lambda _b: close())

    def _on_center(_b=None):
        if center_viewer and get_viewer_id:
            vid = get_viewer_id()
            if vid:
                center_viewer(vid)

    center_btn.on_click(_on_center)

    return {
        'panel': panel, 'show': show, 'close': close,
        'state': state, 'status': status,
        'parse_constraints': _parse_constraints,
        'on_atom_click': on_atom_click,
        'constraints_field': constraints_field,
        'click_sync': click_sync,
        'center_btn': center_btn,
    }


def open_ts_goat_dialog(calc_path, callbacks) -> dict:
    """Full open path: validate, analyse, build the panel, attach, show.

    Called by the browser-tab handler.  ``callbacks`` carries the tab's
    closures: render_viewer, reset_dropdown, scroll_top, backend.
    Returns a dict; raises nothing -- problems land in ``result['error']``
    with ``error_kind`` in {'not_usable_ts', 'no_xyz', 'no_hess',
    'exception'} so the tab can phrase a visible message.
    """
    from delfin import ts_goat

    result: dict = {
        'panel': None, 'registry': None, 'error': None,
        'error_kind': None, 'status_text': None,
    }
    try:
        calc_path = Path(calc_path)
        calc_dir = calc_path.parent
        status = ts_goat.read_ts_status(calc_path)
        # Prime the multi-mode warnings/auto-selection; a >1-imaginary
        # output stays usable -- only 0 imaginary modes blocks here.
        _ = status.auto_mode_index
        if not status.usable:
            problems = '; '.join(status.reasons) or 'not a usable TS output'
            result['error'] = (
                'This output is not a usable converged TS '
                f'({problems}). GOAT needs an OptTS output that converged '
                'with at least one imaginary mode in the final frequency '
                'block.'
            )
            result['error_kind'] = 'not_usable_ts'
            return result
        files = ts_goat.collect_ts_files(calc_dir, out_path=calc_path)
        ts_atoms = ts_goat.read_xyz(files.xyz_path) if files.xyz_path else None
        if ts_atoms is None:
            result['error'] = (
                'No matching TS geometry (.xyz) found in '
                f'{calc_dir.name} -- the atom order must match the '
                'Hessian. Select the geometry manually.'
            )
            result['error_kind'] = 'no_xyz'
            return result
        # TS geometry + selected mode only (see _on_mode_change).
        proposed = ts_goat.analyze_reactive_center(
            files.hess_path, ts_atoms,
            mode_index=status.auto_mode_index,
        ) if files.hess_path else None
        if proposed is None:
            result['error'] = (
                'No Hessian (.hess) found -- the imaginary-mode analysis '
                'and the reactive-centre proposal are not possible '
                'without it.'
            )
            result['error_kind'] = 'no_hess'
            return result
        # Consistency: the XYZ used for the input must match the Hessian
        # the mode came from (same atom count and element order); and the
        # selected mode must actually be one of the Hessian's imaginary
        # modes, not a column from a different calculation.
        if not ts_goat._same_composition(
            ts_atoms, ts_goat._hessian_reference_atoms(files.hess_path)
        ):
            result['error'] = (
                'TS geometry and Hessian do not match (different atom '
                'count or element order) -- refusing to combine them.'
            )
            result['error_kind'] = 'inconsistent_geometry'
            return result
        if status.auto_mode_index is not None and all(
            idx != status.auto_mode_index for idx, _f in status.imaginary_modes
        ):
            result['error'] = (
                'Selected imaginary mode is not listed in the final '
                'frequency block -- Hessian and output may belong to '
                'different calculations.'
            )
            result['error_kind'] = 'inconsistent_mode'
            return result
        charge, mult = status.charge, status.mult
        if (charge is None or mult is None) and files.inp_path:
            inp_charge, inp_mult = ts_goat.read_inp_charge_mult(files.inp_path)
            charge = charge if charge is not None else inp_charge
            mult = mult if mult is not None else inp_mult

        registry = build_ts_goat_panel(
            get_selected_path=lambda: calc_path,
            backend=callbacks['backend'],
            reset_dropdown=callbacks['reset_dropdown'],
            render_viewer=callbacks['render_viewer'],
            scroll_top=callbacks['scroll_top'],
            viewer_output=callbacks.get('viewer_output'),
            render_in_panel=callbacks.get('render_in_panel'),
            on_close=callbacks.get('on_close'),
            center_viewer=callbacks.get('center_viewer'),
            get_viewer_id=callbacks.get('get_viewer_id'),
            animate_in_panel=callbacks.get('animate_in_panel'),
        )
        arrows = None
        mode_arrows_by_index = {}
        try:
            from delfin.imag import displaced_geometries
            for idx, _freq in status.imaginary_modes:
                disp = displaced_geometries(files.hess_path, idx, 0.25)
                # TRUE displacement vectors, not absolute coordinates:
                # the animation displaces base + amp*vector, so the
                # arrow payload must hold (displaced - base), computed
                # in the same frame (both from the Hessian, Angstrom).
                # Passing the absolute displaced coordinates instead
                # made the animation move along the position vector of
                # an unrelated geometry -- never visible as the mode
                # (live-dashboard report 2026-11-02).
                from delfin.imag import _hessian_atoms
                base_atoms = _hessian_atoms(files.hess_path.read_text(
                    encoding='utf-8', errors='ignore').splitlines())
                mode_arrows_by_index[idx] = [
                    (el, dx - bx, dy - by, dz - bz)
                    for (el, dx, dy, dz), (_be, bx, by, bz)
                    in zip(disp['pos'], base_atoms)
                ]
            preselect = status.auto_mode_index
            if preselect is None and status.imaginary_modes:
                preselect = status.imaginary_modes[0][0]
            arrows = mode_arrows_by_index.get(preselect)
        except Exception:
            arrows = None
        warn_html = ''.join(
            f'<div style="color:#e65100; font-size:12px;">⚠ {html.escape(w)}</div>'
            for w in status.warnings
        )
        registry['show']({
            'files': files, 'ts_atoms': ts_atoms, 'proposed': proposed,
            'elements': [a[0] for a in ts_atoms],
            'coordinates': [(x, y, z) for _, x, y, z in ts_atoms],
            'charge': charge if charge is not None else 0,
            'mult': mult if mult is not None else 1,
            'mode_arrows': arrows,
            'mode_arrows_by_index': mode_arrows_by_index,
            'imaginary_modes': list(status.imaginary_modes),
            'hess_path': files.hess_path,
            'status_html': warn_html
            + (
                f'<span style="font-size:12px;">{html.escape(str(Path(calc_path).name))}: '
                f'data quality {proposed.reliability}'
                + (' -- <b>proposal AMBIGUOUS:</b> review and adjust the '
                   'selection in the viewer'
                   if proposed.ambiguous else '')
                + '</span>'
            ),
            'analysis_html': ''.join(
                f'<li>{html.escape(r)}</li>' for r in proposed.reasons
            ),
        })
        result['registry'] = registry
        result['panel'] = registry['panel']
        result['status_text'] = registry['state'].get('status_text')
        return result
    except Exception as exc:  # the dialog must never crash the tab
        result['error'] = (
            f'Could not prepare the TS -> GOAT dialog: '
            f'{type(exc).__name__}: {exc}'
        )
        result['error_kind'] = 'exception'
        return result


def render_constraints_preview(
    elements, coordinates, indices, *, method='XTB', basis='',
    charge=0, mult=1, pal=1, maxcore=1000,
) -> str:
    """The ORCA input text exactly as it will be written, for the dialog."""
    from delfin.ts_goat import render_constraints_block

    keywords = ['!', method]
    if basis:
        keywords.append(basis)
    keywords.append('GOAT')
    lines = [' '.join(k for k in keywords if k), '%pal']
    lines.append(f'  ncores {int(pal)}')
    lines.append('end')
    lines.append(f'%maxcore {int(maxcore)}')
    lines.append('')
    lines.append(render_constraints_block(indices, [
        str(elements[i]) for i in indices if i < len(elements)
    ]))
    lines.append('')
    lines.append(f'* xyz {int(charge)} {int(mult)}')
    for el, (x, y, z) in zip(elements, coordinates):
        lines.append(f'  {el} {x:>14.8f} {y:>14.8f} {z:>14.8f}')
    lines.append('*')
    return '\n'.join(lines)
