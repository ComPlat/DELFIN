"""TS -> GOAT preparation panel for the Calculations Browser.

Widget-only module: it assembles the preparation dialog (status,
constraint chips, method/settings editors, GOAT preview, submit) around
the browser's existing 3Dmol viewer. The library half lives in
``delfin.ts_goat``; this module only draws and collects user decisions.
"""

from __future__ import annotations

import html as _html
from pathlib import Path

import ipywidgets as widgets


def build_ts_goat_panel(
    *,
    get_selected_path,
    backend,
    reset_dropdown,
    render_viewer,
    scroll_top,
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
    ))
    title = widgets.HTML(
        '<b>TS &rarr; GOAT</b> &mdash; constrained conformer search from the '
        'optimized transition state'
    )
    status = widgets.HTML('')
    files_box = widgets.HTML('')
    analysis_box = widgets.HTML('')
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
        description='Close', layout=widgets.Layout(width='90px'),
    )
    btn_row = widgets.HBox([submit_btn, close_btn])
    settings_row = widgets.HBox([
        widgets.VBox([method_label, method_field, basis_field]),
        widgets.VBox([pal_field, maxcore_field, slurm_time_field]),
    ])

    panel.children = (
        title, status, files_box, analysis_box, constraints_label,
        widgets.HBox([constraints_field, highlight_btn]),
        settings_row, preview_label, preview_box, btn_row,
    )

    # ---- state ----------------------------------------------------------
    state = {
        'files': None, 'ts_atoms': None, 'proposed': None,
        'elements': (), 'coordinates': (), 'charge': 0, 'mult': 1,
        'mode_arrows': None,
    }

    def _esc(text):
        return _html.escape(str(text), quote=True)

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
            ('Reactant', files.reactant_path), ('Product', files.product_path),
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
        status.value = open_result.get('status_html', '')
        _refresh_analysis_display()
        analysis_box.value = open_result.get('analysis_html', '')
        _refresh_constraints_display()
        panel.layout.display = 'block'
        scroll_top()

    def close(_button=None):
        panel.layout.display = 'none'

    close_btn.on_click(close)

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
                pal=pal_field.value, maxcore=maxcore_field.value,
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

    constraints_field.observe(lambda _c: _refresh_preview(), names='value')

    return {
        'panel': panel, 'show': show, 'close': close,
        'state': state, 'status': status,
        'parse_constraints': _parse_constraints,
    }


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
