"""ORCA Input Builder tab: build, preview, and submit standalone ORCA jobs."""

import re
import shutil
from pathlib import Path

import py3Dmol
import ipywidgets as widgets
from IPython.display import clear_output

from delfin.common.control_validator import (
    ORCA_FUNCTIONALS, ORCA_BASIS_SETS, DISP_CORR_VALUES, _RI_JKX_KEYWORDS,
)

from .molecule_viewer import strip_xyz_header, apply_molecule_view_style
from .input_processing import parse_inp_resources, sanitize_orca_input


def create_tab(ctx):
    """Create the ORCA Input Builder tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- option lists ---------------------------------------------------
    method_options = sorted(ORCA_FUNCTIONALS)
    basis_options = sorted(ORCA_BASIS_SETS)
    dispersion_options = ['None'] + sorted(v for v in DISP_CORR_VALUES if v)
    ri_options = ['None'] + sorted(_RI_JKX_KEYWORDS)
    aux_basis_options = ['None', 'def2/J', 'def2/JK']

    ws = {'description_width': '120px'}

    # -- widgets --------------------------------------------------------
    orca_job_name = widgets.Text(
        value='', placeholder='e.g. water_opt', description='Job Name:',
        layout=widgets.Layout(width='100%'), style=ws,
    )
    orca_coords = widgets.Textarea(
        value='',
        placeholder=(
            'Paste XYZ coordinates (with or without header):\n\n'
            '6\nComment line\nC  0.0  0.0  0.0\n...\n\nor just:\n'
            'C  0.0  0.0  0.0\n...'
        ),
        description='Coordinates:',
        layout=widgets.Layout(width='100%', height='400px', box_sizing='border-box'), style=ws,
    )
    orca_charge = widgets.IntText(value=0, description='Charge:',
                                  layout=widgets.Layout(width='200px'), style=ws)
    orca_multiplicity = widgets.IntText(value=1, description='Multiplicity:',
                                        layout=widgets.Layout(width='200px'), style=ws)
    orca_method = widgets.Dropdown(options=method_options, value='PBE0',
                                   description='Method:',
                                   layout=widgets.Layout(width='250px'), style=ws)
    orca_job_type = widgets.Dropdown(options=['SP', 'OPT', 'FREQ', 'OPT FREQ'],
                                     value='OPT', description='Job Type:',
                                     layout=widgets.Layout(width='250px'), style=ws)
    orca_basis = widgets.Dropdown(options=basis_options, value='def2-SVP',
                                  description='Basis Set:',
                                  layout=widgets.Layout(width='250px'), style=ws)
    orca_dispersion = widgets.Dropdown(options=dispersion_options, value='D4',
                                       description='Dispersion:',
                                       layout=widgets.Layout(width='250px'), style=ws)
    orca_ri = widgets.Dropdown(options=ri_options, value='RIJCOSX',
                               description='RI Approx:',
                               layout=widgets.Layout(width='250px'), style=ws)
    orca_aux_basis = widgets.Dropdown(options=aux_basis_options, value='def2/J',
                                      description='Aux Basis:',
                                      layout=widgets.Layout(width='250px'), style=ws)
    orca_autoaux = widgets.Checkbox(value=False, description='AutoAux',
                                    layout=widgets.Layout(width='250px'), style=ws)

    orca_solvation_type = widgets.Dropdown(options=['None', 'CPCM', 'SMD'], value='None',
                                           description='Solvation Type:',
                                           layout=widgets.Layout(width='250px'), style=ws)
    orca_solvent = widgets.Dropdown(
        options=['water', 'acetonitrile', 'dmso', 'dmf', 'methanol',
                 'ethanol', 'thf', 'dichloromethane', 'chloroform', 'toluene', 'hexane'],
        value='water', description='Solvent:',
        layout=widgets.Layout(width='250px'), style=ws,
    )

    orca_print_mos = widgets.Checkbox(value=False, description='Print MOs',
                                      layout=widgets.Layout(width='250px'), style=ws)
    orca_print_basis = widgets.Checkbox(value=False, description='Print Basis',
                                        layout=widgets.Layout(width='250px'), style=ws)
    orca_additional = widgets.Text(value='', placeholder='e.g. FinalGrid6 NormalPrint',
                                   description='Additional:',
                                   layout=widgets.Layout(width='400px'), style=ws)

    orca_pal = widgets.IntText(value=40, description='PAL (cores):',
                               layout=widgets.Layout(width='250px'), style=ws)
    orca_maxcore = widgets.IntText(value=6000, description='MaxCore (MB):',
                                   layout=widgets.Layout(width='250px'), style=ws)
    orca_slurm_time = widgets.Text(value='12:00:00', placeholder='e.g. 02:00:00',
                                   description='Time Limit:',
                                   layout=widgets.Layout(width='250px'), style=ws)

    orca_file_upload = widgets.FileUpload(
        accept='', multiple=True, description='Extra Files:',
        layout=widgets.Layout(width='100%', height='30px'), style=ws,
    )
    orca_uploaded_files_label = widgets.HTML(
        value='<i>Drag & drop files here (e.g. .gbw, .xyz, .hess)</i>',
        layout=widgets.Layout(width='100%'),
    )
    orca_path_files = widgets.Textarea(
        value='',
        placeholder='Paste file paths (one per line):\n/path/to/file.gbw\n/path/to/file.xyz',
        description='File Paths:',
        layout=widgets.Layout(width='100%', height='80px', box_sizing='border-box'), style=ws,
    )

    orca_preview = widgets.Textarea(
        value='', description='INP Preview:',
        layout=widgets.Layout(width='100%', height='550px', box_sizing='border-box'), style=ws,
        disabled=False,
    )

    orca_generate_btn = widgets.Button(description='GENERATE INP', button_style='info',
                                       layout=widgets.Layout(width='150px'))
    orca_save_btn = widgets.Button(description='SAVE INP + SH', button_style='warning',
                                   layout=widgets.Layout(width='150px'))
    orca_submit_btn = widgets.Button(description='SUBMIT ORCA JOB', button_style='success',
                                     layout=widgets.Layout(width='150px'))
    orca_output = widgets.Output()

    orca_mol_output = widgets.Output(layout=widgets.Layout(
        border='2px solid #1976d2', width='100%', min_height='300px',
        overflow='hidden', box_sizing='border-box',
    ))

    # -- state ----------------------------------------------------------
    state = {
        'extra_files': {},
        'last_auto_keywords': '',
        'is_resetting': False,
    }

    # -- helpers --------------------------------------------------------
    def _build_keyword_line():
        keywords = [orca_method.value, orca_job_type.value, orca_basis.value]
        if orca_dispersion.value != 'None':
            keywords.append(orca_dispersion.value)
        if orca_ri.value != 'None':
            keywords.append(orca_ri.value)
            keywords.append(orca_aux_basis.value)
            if orca_autoaux.value:
                keywords.append('AutoAux')
        if orca_solvation_type.value != 'None':
            keywords.append(f'{orca_solvation_type.value}({orca_solvent.value})')
        if orca_additional.value.strip():
            keywords.append(orca_additional.value.strip())
        return '! ' + ' '.join(keywords)

    def _build_output_block():
        lines = []
        if orca_print_mos.value:
            lines.append('  print[p_mos] 1')
        if orca_print_basis.value:
            lines.append('  print[p_basis] 2')
        if lines:
            return '%output\n' + '\n'.join(lines) + '\nend'
        return ''

    def generate_orca_input():
        keyword_line = _build_keyword_line()
        pal_block = f'%pal\n  nprocs {orca_pal.value}\nend'
        maxcore_line = f'%maxcore {orca_maxcore.value}'
        output_block = _build_output_block()
        coords = strip_xyz_header(orca_coords.value)
        coord_block = f'* xyz {orca_charge.value} {orca_multiplicity.value}\n{coords}\n*'
        inp = f'{keyword_line}\n\n{pal_block}\n\n{maxcore_line}\n'
        if output_block:
            inp += f'\n{output_block}\n'
        inp += f'\n{coord_block}\n'
        return inp

    def save_uploaded_files(job_dir):
        saved = []
        files_to_save = dict(state['extra_files'])
        if not files_to_save and orca_file_upload.value:
            for f in orca_file_upload.value:
                name = f['name'] if isinstance(f, dict) else f.name
                content = f['content'] if isinstance(f, dict) else f.content
                files_to_save[name] = content
        for filename, content in files_to_save.items():
            (job_dir / filename).write_bytes(content)
            saved.append(filename)
        for line in orca_path_files.value.strip().splitlines():
            p = Path(line.strip())
            if not p.name:
                continue
            if not p.exists():
                print(f'Warning: File not found, skipped: {p}')
                continue
            if p.is_file():
                shutil.copy2(str(p), str(job_dir / p.name))
                saved.append(p.name)
            elif p.is_dir():
                print(f"Warning: '{p.name}' is a directory, skipped.")
        return saved

    def reset_orca_builder():
        state['is_resetting'] = True
        try:
            orca_job_name.value = ''
            orca_coords.value = ''
            orca_charge.value = 0
            orca_multiplicity.value = 1
            orca_method.value = 'PBE0'
            orca_job_type.value = 'OPT'
            orca_basis.value = 'def2-SVP'
            orca_dispersion.value = 'D4'
            orca_ri.value = 'RIJCOSX'
            orca_aux_basis.value = 'def2/J'
            orca_solvation_type.value = 'None'
            orca_solvent.value = 'water'
            orca_print_mos.value = False
            orca_print_basis.value = False
            orca_additional.value = ''
            orca_pal.value = 40
            orca_maxcore.value = 6000
            orca_slurm_time.value = '12:00:00'
            orca_path_files.value = ''
            orca_preview.value = ''
            state['extra_files'].clear()
            state['last_auto_keywords'] = ''
            try:
                orca_file_upload.value = ()
            except Exception:
                pass
        finally:
            state['is_resetting'] = False

        # Restore visual defaults after all fields have been reset.
        orca_uploaded_files_label.value = '<i>Drag & drop files here (e.g. .gbw, .xyz, .hess)</i>'
        update_orca_molecule_view()
        update_orca_preview()
        state['last_auto_keywords'] = _build_keyword_line()

    # -- handlers -------------------------------------------------------
    def update_orca_molecule_view(change=None):
        with orca_mol_output:
            clear_output()
            raw = orca_coords.value.strip()
            if not raw:
                print('Paste XYZ coordinates to see 3D preview.')
                return
            coords = strip_xyz_header(raw)
            if not coords:
                print('No valid coordinates.')
                return
            try:
                lines = [l for l in coords.split('\n') if l.strip()]
                n = len(lines)
                xyz_data = f'{n}\nORCA Builder Preview\n{coords}'
                view = py3Dmol.view(width='100%', height=400)
                view.addModel(xyz_data, 'xyz')
                apply_molecule_view_style(view)
                view.show()
            except Exception as e:
                print(f'Could not visualize: {e}')

    def update_orca_preview(change=None):
        if state.get('is_resetting'):
            return
        current = orca_preview.value.strip()
        if not current:
            orca_preview.value = sanitize_orca_input(generate_orca_input())
            state['last_auto_keywords'] = _build_keyword_line()
            return

        text = current
        new_kw = _build_keyword_line()
        kw_match = re.search(r'^(!.*)$', text, re.MULTILINE)
        if kw_match:
            old_line = kw_match.group(1)
            manual_suffix = ''
            if state['last_auto_keywords'] and old_line.startswith(state['last_auto_keywords']):
                manual_suffix = old_line[len(state['last_auto_keywords']):]
            text = re.sub(r'^!.*$', new_kw + manual_suffix, text, count=1, flags=re.MULTILINE)
        state['last_auto_keywords'] = new_kw

        new_pal = f'%pal\n  nprocs {orca_pal.value}\nend'
        text = re.sub(r'%pal\s*\n\s*nprocs\s+\d+\s*\n\s*end',
                      new_pal, text, count=1, flags=re.IGNORECASE)
        text = re.sub(r'^%maxcore\s+\d+', f'%maxcore {orca_maxcore.value}',
                      text, count=1, flags=re.MULTILINE | re.IGNORECASE)

        new_output = _build_output_block()
        has_output = re.search(r'%output\b.*?\nend', text, flags=re.DOTALL | re.IGNORECASE)
        if has_output:
            if new_output:
                text = re.sub(r'%output\b.*?\nend', new_output, text,
                              count=1, flags=re.DOTALL | re.IGNORECASE)
            else:
                text = re.sub(r'\n?%output\b.*?\nend\n?', '\n', text,
                              count=1, flags=re.DOTALL | re.IGNORECASE)
        elif new_output:
            text = re.sub(r'(\n\* xyz )', f'\n{new_output}\n\\1', text, count=1)

        coords = strip_xyz_header(orca_coords.value)
        new_coord = f'* xyz {orca_charge.value} {orca_multiplicity.value}\n{coords}\n*'
        text = re.sub(r'\*\s*xyz\s+[-\d]+\s+\d+\s*\n.*?\n\*',
                      new_coord, text, count=1, flags=re.DOTALL)
        orca_preview.value = sanitize_orca_input(text)

    def handle_orca_generate(button):
        orca_preview.value = sanitize_orca_input(generate_orca_input())
        state['last_auto_keywords'] = _build_keyword_line()
        with orca_output:
            clear_output()
            print('ORCA input regenerated. You can edit the preview if needed.')

    def _submit_orca_job():
        """Shared logic for save+submit and direct submit."""
        with orca_output:
            clear_output()
            job_name = orca_job_name.value.strip()
            if not job_name:
                print('Error: Job name cannot be empty!')
                return

            safe_job_name = ''.join(c for c in job_name if c.isalnum() or c in ('_', '-'))
            if not safe_job_name:
                print('Error: Job name contains only invalid characters!')
                return

            preview_content = orca_preview.value.strip()
            if preview_content:
                inp_content = preview_content
            else:
                coords = strip_xyz_header(orca_coords.value)
                if not coords:
                    print('Error: Coordinates or INP preview cannot be empty!')
                    return
                inp_content = generate_orca_input()

            inp_content = sanitize_orca_input(inp_content)

            job_dir = ctx.calc_dir / safe_job_name
            job_dir.mkdir(parents=True, exist_ok=True)

            inp_path = job_dir / f'{safe_job_name}.inp'
            inp_path.write_text(inp_content)

            saved_files = save_uploaded_files(job_dir)

            pal_used, maxcore_used = parse_inp_resources(inp_content)
            if pal_used is None:
                pal_used = orca_pal.value
            if maxcore_used is None:
                maxcore_used = orca_maxcore.value

            result = ctx.backend.submit_orca(
                job_dir=job_dir, job_name=safe_job_name,
                inp_file=f'{safe_job_name}.inp',
                time_limit=orca_slurm_time.value,
                pal=pal_used, maxcore=maxcore_used,
            )

            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                print('ORCA job successfully submitted!')
                print(f'Job ID: {job_id}')
                print(f'Job Name: {safe_job_name}')
                print(f'Directory: {job_dir}')
                print(f'PAL: {pal_used}, MaxCore: {maxcore_used} MB')
                print(f'Time Limit: {orca_slurm_time.value}')
                if orca_solvation_type.value != 'None':
                    print(f'Solvation: {orca_solvation_type.value}({orca_solvent.value})')
                if saved_files:
                    print(f"Extra files: {', '.join(saved_files)}")
                print()
                print('Check status in Job Status tab')
                reset_orca_builder()
            else:
                print('Error submitting job:')
                print(result.stderr or result.stdout)

    def handle_orca_save(button):
        _submit_orca_job()

    def handle_orca_submit(button):
        _submit_orca_job()

    def update_uploaded_files_label(change=None):
        if orca_file_upload.value:
            for f in orca_file_upload.value:
                name = f['name'] if isinstance(f, dict) else f.name
                content = f['content'] if isinstance(f, dict) else f.content
                state['extra_files'][name] = content
        if state['extra_files']:
            filenames = list(state['extra_files'].keys())
            orca_uploaded_files_label.value = f"<b>Files to upload:</b> {', '.join(filenames)}"
        else:
            orca_uploaded_files_label.value = '<i>Drag & drop files here (e.g. .gbw, .xyz, .hess)</i>'

    # -- wiring ---------------------------------------------------------
    orca_generate_btn.on_click(handle_orca_generate)
    orca_save_btn.on_click(handle_orca_save)
    orca_submit_btn.on_click(handle_orca_submit)
    orca_file_upload.observe(update_uploaded_files_label, names='value')

    for w in [orca_method, orca_job_type, orca_basis, orca_dispersion, orca_ri,
              orca_aux_basis, orca_charge, orca_multiplicity, orca_pal, orca_maxcore,
              orca_coords, orca_additional, orca_solvation_type, orca_solvent,
              orca_print_mos, orca_print_basis, orca_autoaux]:
        w.observe(update_orca_preview, names='value')

    orca_coords.observe(update_orca_molecule_view, names='value')
    update_orca_molecule_view()
    update_orca_preview()
    state['last_auto_keywords'] = _build_keyword_line()

    # -- layout ---------------------------------------------------------
    def _row(children, wrap=True):
        return widgets.HBox(
            children,
            layout=widgets.Layout(
                width='100%',
                min_width='0',
                gap='8px',
                align_items='center',
                flex_wrap='wrap' if wrap else 'nowrap',
            ),
        )

    orca_left = widgets.VBox([
        _row([orca_job_name], wrap=False),
        _row([orca_coords], wrap=False),
        _row([orca_charge, orca_multiplicity]),
        _row([orca_method, orca_job_type]),
        _row([orca_basis, orca_dispersion]),
        _row([orca_ri, orca_aux_basis]),
        _row([orca_autoaux]),
        _row([orca_solvation_type, orca_solvent]),
        _row([orca_print_mos, orca_print_basis]),
        _row([orca_additional], wrap=False),
        _row([orca_pal, orca_maxcore]),
        _row([orca_slurm_time]),
        widgets.VBox([orca_file_upload, orca_uploaded_files_label],
                     layout=widgets.Layout(width='100%', min_width='0', overflow='hidden')),
        _row([orca_path_files], wrap=False),
        _row([orca_generate_btn, orca_save_btn, orca_submit_btn]),
        orca_output,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='10px',
        box_sizing='border-box', overflow_x='hidden',
    ))

    orca_right = widgets.VBox([
        widgets.HTML('<b>ORCA Input Preview (editable):</b>'),
        orca_preview,
        widgets.HTML('<b>Molecule Preview:</b>', layout=widgets.Layout(margin='10px 0 0 0')),
        orca_mol_output,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='10px',
        box_sizing='border-box', overflow_x='hidden',
    ))

    split = widgets.HBox(
        [orca_left, orca_right],
        layout=widgets.Layout(width='100%', align_items='stretch', overflow_x='hidden'),
    )
    orca_css = widgets.HTML(
        """
        <style>
        .orca-split-root, .orca-split-root * {
            box-sizing: border-box;
        }
        .orca-split-root {
            width: 100% !important;
            overflow-x: hidden !important;
        }
        .orca-split-pane {
            min-width: 0 !important;
            overflow-x: hidden !important;
        }
        .orca-split-pane .widget-box,
        .orca-split-pane .widget-hbox,
        .orca-split-pane .widget-vbox {
            max-width: 100% !important;
        }
        .orca-split-pane .widget-output,
        .orca-split-pane .output_area,
        .orca-split-pane .output_subarea,
        .orca-split-pane .output_wrapper,
        .orca-split-pane .jp-OutputArea,
        .orca-split-pane .jp-OutputArea-child,
        .orca-split-pane .jp-OutputArea-output {
            max-width: 100% !important;
            overflow-x: hidden !important;
        }
        </style>
        """
    )
    split.add_class('orca-split-root')
    orca_left.add_class('orca-split-pane')
    orca_right.add_class('orca-split-pane')

    tab_widget = widgets.VBox([
        widgets.HTML('<h3>ORCA Input Builder</h3>'),
        widgets.HTML(
            '<a href="https://orca-manual.mpi-muelheim.mpg.de/" '
            'target="_blank">ORCA User Manual</a>'
        ),
        orca_css,
        split,
    ], layout=widgets.Layout(width='100%', padding='10px', box_sizing='border-box'))

    return tab_widget, {'orca_pal': orca_pal, 'orca_maxcore': orca_maxcore}
