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

from delfin import ts_goat_traj


def _two_column_layout(title_row, left, right):
    """Title bar + left (info/settings) | right (3D viewer) columns.

    On narrow browser windows ipywidgets' FlexBox wraps the columns
    below each other, so no separate mobile layout is needed.
    """
    body = widgets.HBox(
        [widgets.VBox(tuple(left)), widgets.VBox(tuple(right))],
        layout=widgets.Layout(
            display='flex', flex_flow='row wrap', width='100%',
            gap='12px', align_items='flex-start',
        ),
    )
    left_box, right_box = body.children
    # The outer panel owns scrolling. Nested scrolling on the left made
    # the controls and Submit button inaccessible at laptop heights.
    left_box.layout = widgets.Layout(
        flex='0 1 35%', min_width='300px', gap='7px')
    right_box.layout = widgets.Layout(
        flex='1 1 63%', min_width='380px', gap='5px')
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
        max_height='calc(100vh - 180px)', overflow_y='auto',
    ))
    title = widgets.HTML(
        '<b>TS &rarr; GOAT</b> '
        '<span style="color:#607d8b;font-size:12px;">'
        'Constrained conformer search from an optimized TS</span>'
    )
    # 2026-11-02 live report: the close button sat in the bottom action
    # row and scrolled out of sight.  The title bar now carries its own
    # close button at the right edge, always visible without scrolling.
    title_close_btn = widgets.Button(
        description='✕ Close', button_style='warning',
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
        width='100%', height='450px', border='1px solid #d5dfe5',
        overflow='hidden', margin='2px 0',
    ))
    # Stable sizing reference for the 3Dmol stage. The Output widget can
    # become wider when DevTools is closed or the columns are resized;
    # measuring only the inner rendered HTML left the canvas at an older
    # width, giving a left-aligned molecule inside a larger white frame.
    viewer_out.add_class('ts-goat-viewer-output')
    # JS -> Python click channel (same pattern as the structure editor's
    # submit-pick-sync): the viewer bootstrap writes the clicked atom's
    # 0-based index into this hidden Text via the property setter +
    # input event; the observe() below routes it to on_atom_click.
    click_sync = widgets.Text(value='', layout=widgets.Layout(display='none'))
    click_sync.add_class('ts-goat-click-sync')
    mode_label = widgets.HTML('<b>Reactive mode (imaginary modes of the final block):</b>')
    mode_dropdown = widgets.Dropdown(
        options=[], description='Mode:',
        # Reserve enough width for the frequency AND the native dropdown
        # arrow.  A percentage width here shrank to a tiny select in
        # Voila, obscuring the last digits of the imaginary frequency.
        style={'description_width': '48px'},
        layout=widgets.Layout(
            width='340px', min_width='320px', max_width='100%',
            flex='0 0 340px'),
    )
    mode_warning = widgets.HTML('')
    # Imaginary-mode animation: validated orca_pltvib frames are
    # loaded into the viewer (addModelsAsFrames) and played back with
    # setFrame (fixed camera).  Play auto-generates the trajectory on
    # first click; a missing binary keeps the static TS visible with
    # an explanatory status line.
    anim_status = widgets.HTML('')
    anim_play = widgets.Button(
        description='▶ Play', disabled=False,
        tooltip='Play the vibration trajectory in a loop',
        layout=widgets.Layout(width='90px'),
    )
    anim_pause = widgets.Button(
        description='⏸ Pause', disabled=True,
        tooltip='Pause the playback; Reset Geometry restores the TS',
        layout=widgets.Layout(width='90px'),
    )
    anim_reset = widgets.Button(
        description='Reset Geometry', disabled=True,
        tooltip='Show the optimized TS geometry again',
        layout=widgets.Layout(width='140px'),
    )
    # Separate labels from the slider widgets so Voila cannot squeeze
    # each track between an 80px label and its numeric readout. Each
    # slider now gets a large usable horizontal drag target.
    anim_amp = widgets.FloatSlider(
        value=0.25, min=0.05, max=1.0, step=0.05,
        description='', continuous_update=False, readout=True,
        readout_format='.2f',
        tooltip='Displacement amplitude (frames are rescaled)',
        style={'description_width': '0px'},
        layout=widgets.Layout(
            flex='1 1 260px', width='100%', min_width='250px'),
    )
    anim_speed = widgets.FloatSlider(
        value=1.0, min=0.25, max=3.0, step=0.25,
        description='', continuous_update=False, readout=True,
        readout_format='.2f',
        tooltip='Playback speed (1 = 5 frames per second)',
        style={'description_width': '0px'},
        layout=widgets.Layout(
            flex='1 1 260px', width='100%', min_width='250px'),
    )
    amp_group = widgets.HBox(
        [widgets.HTML('<b style="font-size:12px;">Amplitude</b>',
                      layout=widgets.Layout(width='80px',
                                            flex='0 0 80px')),
         anim_amp],
        layout=widgets.Layout(
            flex='1 1 390px', min_width='345px', align_items='center',
            gap='6px'),
    )
    speed_group = widgets.HBox(
        [widgets.HTML('<b style="font-size:12px;">Speed</b>',
                      layout=widgets.Layout(width='62px',
                                            flex='0 0 62px')),
         anim_speed],
        layout=widgets.Layout(
            flex='1 1 370px', min_width='325px', align_items='center',
            gap='6px'),
    )
    # One visible Play action: it loads the trajectory on demand.  The
    # previous "Imaginary Mode" toggle + "Generate / Load" button duplicated
    # that action and left Play disabled in the default state.
    anim_row = widgets.VBox(
        [widgets.HBox(
            [anim_play, anim_pause, anim_reset],
            layout=widgets.Layout(flex_flow='row wrap', width='100%',
                                  gap='5px', margin='2px 0')),
         widgets.HBox(
            [amp_group, speed_group],
            layout=widgets.Layout(
                display='flex', flex_flow='row wrap', width='100%',
                gap='12px', margin='5px 0 2px 0')),
         anim_status],
    )
    constraints_label = widgets.HTML(
        '<b>Constraints</b> '
        '<span style="font-size:11px;color:#607d8b;">'
        '(0-based atom indices)</span>')
    constraints_count = widgets.HTML('')
    constraints_field = widgets.Textarea(
        placeholder='e.g. 0, 2, 5 (comma- or newline-separated)',
        layout=widgets.Layout(width='100%', height='58px'),
    )
    highlight_btn = widgets.Button(
        description='Apply', button_style='info',
        tooltip='Re-render the viewer with the current constraint atoms highlighted',
        layout=widgets.Layout(width='140px'),
    )
    reset_prop_btn = widgets.Button(
        description='Restore proposal',
        tooltip='Restore the automatic reactive-centre proposal',
        layout=widgets.Layout(width='150px'),
    )
    clear_sel_btn = widgets.Button(
        description='Clear all',
        tooltip='Clear every constraint (fully manual selection)',
        layout=widgets.Layout(width='100px'),
    )
    method_label = widgets.HTML('<b>GOAT settings</b>')
    # Offer known GOAT-compatible ORCA keywords, not an unvalidated
    # free-text method name.  XTB is GFN2-xTB in ORCA.
    method_field = widgets.Dropdown(
        options=[
            ('GFN2-xTB (XTB)', 'XTB'),
            ('GFN1-xTB (XTB1)', 'XTB1'),
            ('GFN-FF force field (XTBFF)', 'XTBFF'),
            ('r2SCAN-3c (DFT)', 'r2SCAN-3c'),
            ('B3LYP (DFT)', 'B3LYP'),
        ],
        value='XTB', description='Method:',
        layout=widgets.Layout(width='250px'),
    )
    basis_field = widgets.Text(
        value='', placeholder='def2-SVP (DFT only)',
        description='Basis:', layout=widgets.Layout(width='240px'),
    )
    pal_field = widgets.IntText(
        value=12, description='Cores:', layout=widgets.Layout(width='145px'),
    )
    maxcore_field = widgets.IntText(
        value=1000, description='MB/core:',
        layout=widgets.Layout(width='190px'),
    )
    slurm_time_field = widgets.Text(
        value='48:00:00', description='SLURM time:',
        layout=widgets.Layout(width='210px'),
    )
    preview_label = widgets.HTML('<b>ORCA input preview:</b>')
    preview_box = widgets.HTML('', layout=widgets.Layout(
        max_height='18em', overflow_y='auto', width='100%',
        border='1px solid #cfd8dc'))
    submit_btn = widgets.Button(
        description='▶ Submit GOAT job', button_style='success',
        tooltip='Write the input into a new folder and submit it',
        layout=widgets.Layout(width='200px'),
    )
    close_btn = widgets.Button(
        description='✕ Close', layout=widgets.Layout(width='110px'),
    )
    center_btn = widgets.Button(
        description='Center View',
        tooltip='Re-center and re-zoom the 3D view',
        layout=widgets.Layout(width='150px'),
    )
    cancel_btn = widgets.Button(
        description='Cancel',
        tooltip='Close TS→GOAT and restore the normal file view',
        layout=widgets.Layout(width='110px'),
    )
    # Keep the raw indices and long reports folded. The resource
    # settings that determine the actual SLURM job are always visible.
    accordion_layout = widgets.Layout(
        width='94%', align_self='flex-start', margin='1px 0')
    constraint_editor = widgets.Accordion(
        children=[constraints_field], selected_index=None,
        layout=accordion_layout)
    constraint_editor.set_title(0, 'Edit atom indices')
    input_details = widgets.Accordion(
        children=[preview_box], selected_index=None,
        layout=widgets.Layout(width='94%', align_self='flex-start',
                              margin='1px 0'))
    input_details.set_title(0, 'Generated ORCA input')
    # The verbose diagnostic report must remain accessible even when the
    # browser window is short. A bounded, independently scrollable
    # expansion avoids the truncated analysis in the last screenshot.
    analysis_scroller = widgets.VBox(
        [files_box, analysis_box],
        layout=widgets.Layout(
            max_height='235px', overflow_y='auto', overflow_x='hidden',
            width='100%', padding='4px 8px 8px 4px'))
    source_details = widgets.Accordion(
        children=[analysis_scroller], selected_index=None,
        layout=widgets.Layout(width='94%', align_self='flex-start',
                              margin='1px 0'))
    source_details.set_title(0, 'Analysis & source files')
    settings_row = widgets.HBox(
        [method_field, pal_field],
        layout=widgets.Layout(display='flex', flex_flow='row wrap',
                              gap='6px', width='100%'))
    resources_row = widgets.HBox(
        [maxcore_field, slurm_time_field],
        layout=widgets.Layout(display='flex', flex_flow='row wrap',
                              gap='6px', width='100%'))
    settings_group = widgets.VBox(
        [settings_row, resources_row, basis_field],
        layout=widgets.Layout(width='94%', gap='4px'))
    constraint_actions = widgets.HBox(
        [reset_prop_btn, clear_sel_btn],
        layout=widgets.Layout(display='flex', flex_flow='row wrap',
                              gap='6px', width='100%'))
    btn_row = widgets.HBox(
        [submit_btn, cancel_btn],
        layout=widgets.Layout(display='flex', flex_flow='row wrap',
                              gap='8px', width='100%'))
    viewer_head = widgets.HTML(
        '<b>Optimized TS structure</b> '
        '<span style="font-size:11px;color:#607d8b;">'
        'Click atoms to toggle constraints; drag to rotate</span>')
    # Add centering alongside playback (instead of a third control row).
    anim_row.children[0].children = (
        anim_play, anim_pause, anim_reset, center_btn)

    layout = _two_column_layout(
        title_row=widgets.HBox([title, title_close_btn],
            layout=widgets.Layout(display='flex', width='100%',
                                   align_items='center')),
        left=(status, constraints_label, constraints_count,
              constraint_actions, constraint_editor, method_label,
              settings_group, input_details, source_details, btn_row),
        right=(viewer_head, mode_dropdown, mode_warning,
               viewer_out, click_sync, anim_row),
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

    def _update_constraint_count(indices, n_atoms):
        """Atom count, share and a factual wide-selection warning.

        A large fixed fraction restricts the GOAT conformer search
        heavily (ORCA TS-conformer tutorial: only the atoms directly
        involved in the transformation should be fixed); the user
        decides, so this is a warning, not a block.
        """
        if not n_atoms:
            constraints_count.value = ''
            return
        n_sel = len(set(indices))
        share = 100.0 * n_sel / n_atoms
        color = '#b26a00' if share > 30 else '#486776'
        warn = ''
        if share > 30:
            warn = (' <span style="color:#b26a00;" '
                    'title="A large fixed fraction restricts GOAT sampling">'
                    '&#9888; restrictive</span>')
        constraints_count.value = (
            f'<span style="color:{color}; font-size:11px;">'
            f'{n_sel} / {n_atoms} atoms constrained '
            f'({share:.0f}%){warn}</span>')

    def _sync_basis_field():
        """Basis sets apply to orbital methods only.

        For semi-empirical/force-field methods (XTB, GFN*, MM) a basis
        like def2-SVP is not a real parameter; the field is disabled so
        a stale value can neither mislead the user nor enter the input.
        """
        method = (method_field.value or '').strip()
        # Self-contained composite DFT (e.g. r2SCAN-3c) has its own
        # basis; only ordinary orbital DFT requires a separate choice.
        no_external_basis = (method.lower().startswith(
            ('xtb', 'gfn', 'mm')) or method.lower().endswith('-3c'))
        basis_field.disabled = no_external_basis
        basis_field.layout.display = 'none' if no_external_basis else 'block'
        if no_external_basis:
            basis_field.value = ''
        elif not basis_field.value.strip():
            basis_field.value = 'def2-SVP'

    def _refresh_preview():
        indices = _parse_constraints(constraints_field.value)
        atoms = state['ts_atoms'] or []
        n_atoms = len(atoms)
        bad = [i for i in indices if i >= n_atoms]
        _update_constraint_count(indices, n_atoms)
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
        files_box.value = (
            '<ul style="font-size:11px;line-height:1.6;'
            'padding-left:18px;margin:4px 0;overflow-wrap:anywhere;">'
            + ''.join(rows) + '</ul>' + problems
        )

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
        try:
            # The ORCA-generated XYZ frames drive animation. Displacement
            # arrows are no longer a visible control; selection halos
            # remain independent of the animation payload.
            render_in_panel(
                atoms, indices, None,
                output_widget=viewer_out,
                mode_vector=state.get('mode_arrows'),
                frames_payload=_load_frames_payload(),
            )
        except Exception as exc:  # visible error, never a blank area
            _viewer_fallback(
                f'Viewer could not be rendered: {type(exc).__name__}: {exc}')

    # ---- vibration trajectory (orca_pltvib) ------------------------------
    # The trajectory is generated once per (hess, mode) in a cache
    # directory next to the Hessian, validated (atom order + TS
    # geometry match), then loaded into the viewer as frames.  Play
    # runs entirely in the browser (setFrame loop); speed and
    # amplitude re-load the payload with re-scaled frames.

    def _load_frames_payload():
        """Frames for the current mode + amplitude, or None on failure."""
        traj = state.get('traj')
        if traj is None or not traj.frames:
            return None
        scaled = traj.scaled(anim_amp.value)
        return ts_goat_traj.frames_to_xyz(scaled)

    def _set_anim_buttons():
        has_frames = bool(state.get('traj') and state['traj'].frames)
        # Play generates the trajectory on first click; it must not start
        # disabled just because the trajectory has not been loaded yet.
        anim_play.disabled = (not state.get('hess_path')
                              or mode_dropdown.value is None)
        anim_pause.disabled = not has_frames
        anim_reset.disabled = not has_frames

    def _anim_error(message):
        anim_status.value = (
            f'<span style="color:#d32f2f; font-size:12px;">'
            f'{html.escape(message)}</span>')

    def _anim_info(message):
        anim_status.value = (
            f'<span style="color:#2e7d32; font-size:12px;">'
            f'{html.escape(message)}</span>')

    def _generate_traj(_b=None):
        """Generate/load + validate the trajectory for the selected mode."""
        hess = state.get('hess_path')
        idx = mode_dropdown.value
        if hess is None or idx is None:
            _anim_error('No Hessian or imaginary mode available.')
            return
        anim_status.value = (
            f'<span style="color:#1976d2; font-size:12px;">'
            f'Generating vibration trajectory for mode {idx} ...</span>')
        result = ts_goat_traj.generate_vibration_trajectory(
            Path(hess), int(idx),
            frequency_cm1=None,
        )
        if result.error or result.trajectory is None:
            state['traj'] = None
            _set_anim_buttons()
            _anim_error(result.error or 'unknown trajectory failure')
            # the static TS geometry stays visible -- the viewer is
            # NOT touched on failure
            return
        state['traj'] = result.trajectory
        _set_anim_buttons()
        freq = result.trajectory.frequency_cm1
        freq_txt = f' ({freq:.1f} cm-1)' if freq is not None else ''
        _anim_info(
            f'Mode {idx}{freq_txt}: {result.trajectory.n_frames} frames '
            + ('loaded' if result.reused else 'generated') + '.')
        _render_integrated_viewer()

    def _on_anim_amp(_change=None):
        # Amplitude rescales the frames (linear in the mode vector);
        # the trajectory itself is reused.
        _render_integrated_viewer()

    anim_amp.observe(_on_anim_amp, names='value')

    def _on_speed_change(_change=None):
        """Write the frame delay for the browser loop (base 200 ms)."""
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            animate_in_panel(vid, 'speed')

    def _on_play(_b=None):
        # one-click workflow: generate the trajectory if missing
        if not state.get('traj'):
            _generate_traj()
            if not state.get('traj'):
                return
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            # setFrame loop over the preloaded frames in the browser;
            # the camera never moves
            animate_in_panel(vid, True)
        else:
            _render_integrated_viewer()

    def _on_pause(_b=None):
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            animate_in_panel(vid, 'pause')

    def _on_reset(_b=None):
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            animate_in_panel(vid, False)
        # ORCA's first vibration frame is not necessarily the optimized
        # saddle point. A genuine reset therefore renders the original
        # XYZ, not frame 0 of the displaced trajectory.
        state['traj'] = None
        _set_anim_buttons()
        anim_status.value = ''
        _render_integrated_viewer()

    anim_play.on_click(_on_play)
    anim_pause.on_click(_on_pause)
    anim_reset.on_click(_on_reset)
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
        # The field observer updates the preview and viewer exactly once.
        constraints_field.value = ', '.join(str(i) for i in indices)

    def _on_mode_change(change=None):
        """Re-run the reactive-centre analysis for the selected mode."""
        if change is not None and change.get('name') != 'value':
            return
        idx = mode_dropdown.value
        if idx is None or state['ts_atoms'] is None or not state['hess_path']:
            return
        # A different mode invalidates the loaded XYZ frames.  Stop the
        # previous playback before another trajectory can be generated.
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            animate_in_panel(vid, False)
        state['traj'] = None
        anim_status.value = ''
        _set_anim_buttons()
        state['mode_arrows'] = state['mode_arrows_by_index'].get(idx)
        from delfin import ts_goat as _tg
        # TS geometry + selected mode only -- a plain OptTS run has no
        # reactant/product files and the analysis must not need them.
        proposed = _tg.analyze_reactive_center(
            state['hess_path'], state['ts_atoms'], mode_index=idx,
        )
        state['proposed'] = proposed
        # Structure the analysis: a short verdict up top, the full
        # contact statistics folded away (live report 2026-11-02: the
        # plain reasons list buried the essential facts under 800+
        # contacts of diagnostics).
        head = (
            f'<b>Proposal: {len(proposed.indices)} atoms</b>'
            f' &mdash; data quality: {proposed.reliability}'
            + (' &mdash; <b style="color:#e65100;">AMBIGUOUS:</b> '
               'several contacts answer the mode similarly strongly; '
               'review the selection in the viewer'
               if proposed.ambiguous else
               ' &mdash; proposal unambiguous')
        )
        detail = ''.join(
            f'<li>{html.escape(r)}</li>' for r in proposed.reasons
        )
        analysis_box.value = (
            f'<div style="font-size:11px;overflow-wrap:anywhere;">{head}</div>'
            '<details style="font-size:11px;margin-top:6px;">'
            '<summary>Contact statistics</summary>'
            f'<ul style="padding-left:18px;margin:5px 0;">{detail}</ul>'
            '</details>'
        )
        # reliability = quality/consistency of the data; ambiguity =
        # uniqueness of the PROPOSAL.  These are different things and
        # must not read as one verdict (live report 2026-11-02).
        status.value = (
            f'<span style="font-size:12px;color:#526a76;">'
            f'Data quality: {_esc(proposed.reliability)}'
            + (' &middot; <b>Proposal needs review</b>'
               if proposed.ambiguous else
               ' &middot; Proposal consistent')
            + '</span>'
        )
        _refresh_constraints_display()

    mode_dropdown.observe(_on_mode_change, names='value')

    highlight_btn.on_click(lambda _b: _render_integrated_viewer())

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
        state['traj'] = None
        anim_play.disabled = False
        anim_pause.disabled = True
        anim_reset.disabled = True
        anim_status.value = ''
        state['reactant_path'] = open_result.get('reactant_path')
        state['product_path'] = open_result.get('product_path')
        mode_opts = [
            (f'Mode {i}: {f:.1f} cm⁻¹', i)
            for i, f in state['imaginary_modes']
        ]
        mode_dropdown.options = mode_opts or [('(no imaginary mode)', None)]
        mode_dropdown.value = (
            mode_opts[0][1] if mode_opts else None
        )
        # single mode: the picker is informational; multiple: the user picks
        mode_label.layout.display = 'none'
        mode_dropdown.layout.display = 'block'
        mode_warning.value = (
            '<div style="color:#e65100; font-size:12px;">'
            '⚠ Several imaginary modes in the final frequency block -- '
            'select the reactive mode; the constraints are re-proposed '
            'for the selected mode.</div>'
        ) if len(mode_opts) > 1 else ''
        proposed = state['proposed']
        status.value = (
            f'<span style="font-size:12px;color:#526a76;">'
            f'Data quality: {_esc(proposed.reliability)}'
            + (' &middot; <b>Proposal needs review</b>'
               if proposed.ambiguous else
               ' &middot; Proposal consistent')
            + '</span>'
        )
        _refresh_analysis_display()
        # The detailed contact list stays folded even for the first
        # display; the expanded accordion itself can scroll.
        analysis_box.value = (
            f'<div style="font-size:11px;">'
            f'Proposal: {len(proposed.indices)} atoms'
            f' &middot; {_esc(proposed.reliability)} data</div>'
            '<details style="font-size:11px;margin-top:6px;">'
            '<summary>Contact statistics</summary>'
            '<ul style="padding-left:18px;margin:5px 0;">'
            + open_result.get('analysis_html', '') +
            '</ul></details>'
        )
        _sync_basis_field()
        _refresh_constraints_display()
        _render_integrated_viewer()
        panel.layout.display = 'block'
        scroll_top()

    def close(_button=None):
        # Stop browser-side playback before detaching the viewer.
        vid = get_viewer_id() if get_viewer_id else None
        if vid and animate_in_panel:
            animate_in_panel(vid, False)
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

    method_field.observe(
        lambda _c: (_sync_basis_field(), _refresh_preview()),
        names='value')
    # The preview is the exact input the user is about to submit: resource
    # and basis edits must update it, not only method/constraint edits.
    for editor in (basis_field, pal_field, maxcore_field):
        editor.observe(lambda _c: _refresh_preview(), names='value')
    constraints_field.observe(lambda _c: (_refresh_preview(),
                                          _render_integrated_viewer()),
                              names='value')

    def _on_click_sync(change):
        raw = (change.get('new') or '').strip()
        if not raw:
            return
        try:
            # The JS click bridge adds a changing sequence suffix so
            # two consecutive clicks on the same atom still trigger observe.
            on_atom_click(int(raw.split(':', 1)[0]))
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
