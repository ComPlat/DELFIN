"""Submit Job tab: main DELFIN job submission form."""

import re

import py3Dmol
import ipywidgets as widgets
from IPython.display import clear_output

from delfin.config import validate_control_text
from delfin.smiles_converter import contains_metal

from .constants import COMMON_LAYOUT, COMMON_STYLE
from .helpers import resolve_time_limit, create_time_limit_widgets, disable_spellcheck
from .molecule_viewer import apply_molecule_view_style
from .input_processing import (
    smiles_to_xyz, smiles_to_xyz_isomers, is_smiles, clean_input_data,
    parse_resource_settings,
)


def create_tab(ctx):
    """Create the Submit Job tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    SUBMIT_MOL_HEIGHT = 500

    # -- widgets --------------------------------------------------------
    job_name_widget = widgets.Text(
        value='', placeholder='e.g. Fe_Complex_Ox',
        description='Job Name:', layout=COMMON_LAYOUT, style=COMMON_STYLE,
    )

    job_type_widget, custom_time_widget = create_time_limit_widgets()

    coords_help = widgets.Label('Input: XYZ coordinates (header auto-removed) or SMILES string')
    coords_widget = widgets.Textarea(
        value='',
        placeholder=(
            'Paste XYZ coordinates or SMILES:\n\n'
            'XYZ example:\n42\nComment\nFe  0.0  0.0  0.0\nC   1.5  0.0  0.0\n\n'
            'SMILES example:\nCCO or c1ccccc1'
        ),
        layout=widgets.Layout(width='100%', height='200px', box_sizing='border-box'),
        style=COMMON_STYLE,
    )

    button_spacer = widgets.Label(value='', layout=widgets.Layout(height='10px'))

    convert_smiles_button = widgets.Button(
        description='CONVERT SMILES', button_style='info',
        layout=widgets.Layout(width='150px'),
        tooltip='Click: full isomer search | Double-click: quick single structure',
    )
    convert_smiles_uff_button = widgets.Button(
        description='CONVERT SMILES + UFF', button_style='info',
        layout=widgets.Layout(width='185px'),
    )

    build_complex_button = widgets.Button(
        description='BUILD COMPLEX', button_style='warning',
        layout=widgets.Layout(width='150px'),
    )

    guppy_submit_button = widgets.Button(
        description='SUBMIT GUPPY', button_style='warning',
        layout=widgets.Layout(width='170px'),
    )

    smiles_batch_help = widgets.Label("Batch SMILES list: one per line 'Name;SMILES;key=value;key=value'")
    smiles_batch_widget = widgets.Textarea(
        value='',
        placeholder='name;SMILES;key=value;...\nNi_1;[Ni];charge=2;solvent=water\nCo_1;[Co];charge=3',
        layout=widgets.Layout(width='100%', height='160px', box_sizing='border-box'),
        style=COMMON_STYLE,
    )

    submit_smiles_list_button = widgets.Button(
        description='SUBMIT SMILES LIST', button_style='success',
        layout=widgets.Layout(width='150px'),
    )
    smiles_batch_output = widgets.Output()

    smiles_prev_button = widgets.Button(
        description='\u25c0', button_style='info',
        layout=widgets.Layout(width='35px'),
    )
    smiles_next_button = widgets.Button(
        description='\u25b6', button_style='info',
        layout=widgets.Layout(width='35px'),
    )
    smiles_preview_label = widgets.HTML(
        value='<span style="font-size:12px;">0 / 0</span>',
        layout=widgets.Layout(width='50px', margin='0 2px'),
    )

    control_help = widgets.Label('CONTROL.txt - edit parameters as needed')
    control_widget = widgets.Textarea(
        value=ctx.default_control,
        layout=widgets.Layout(width='100%', height='500px', box_sizing='border-box'),
        style=COMMON_STYLE,
    )
    control_widget.add_class('delfin-nospell')
    disable_spellcheck(ctx)

    submit_button = widgets.Button(
        description='SUBMIT JOB', button_style='primary',
        layout=widgets.Layout(width='150px'),
    )
    validate_button = widgets.Button(
        description='VALIDATE CONTROL', button_style='warning',
        layout=widgets.Layout(width='150px'),
    )

    output_area = widgets.Output()
    validate_output = widgets.Output()

    mol_output = widgets.Output(layout=widgets.Layout(
        border='2px solid #1976d2', width='100%', height=f'{SUBMIT_MOL_HEIGHT}px',
        overflow='hidden', box_sizing='border-box',
    ))

    xyz_copy_btn = widgets.Button(
        description='\U0001f4cb Copy Coordinates', button_style='success',
        layout=widgets.Layout(width='150px'), disabled=True,
    )
    xyz_copy_status = widgets.HTML(value='', layout=widgets.Layout(margin='0 0 0 6px'))

    isomer_prev_btn = widgets.Button(
        description='\u25c0', button_style='info',
        layout=widgets.Layout(width='35px'),
    )
    isomer_next_btn = widgets.Button(
        description='\u25b6', button_style='info',
        layout=widgets.Layout(width='35px'),
    )
    isomer_label = widgets.HTML(
        value='', layout=widgets.Layout(width='180px'),
    )
    isomer_nav_row = widgets.HBox(
        [isomer_prev_btn, isomer_label, isomer_next_btn],
        layout=widgets.Layout(gap='4px', align_items='center', display='none'),
    )

    only_goat_label = widgets.Label('Only GOAT:')
    only_goat_charge = widgets.IntText(
        value=0, description='Charge:', style=COMMON_STYLE,
        layout=widgets.Layout(width='160px'),
    )
    only_goat_solvent = widgets.Text(
        value='', placeholder='e.g. water', description='Solvent:',
        style=COMMON_STYLE, layout=widgets.Layout(width='220px'),
    )
    only_goat_submit_button = widgets.Button(
        description='SUBMIT ONLY GOAT', button_style='success',
        layout=widgets.Layout(width='150px'),
    )
    only_goat_output = widgets.Output()

    co2_species_delta = widgets.IntText(
        value=-2, description='Species delta:',
        style=COMMON_STYLE, layout=widgets.Layout(width='180px'),
    )
    co2_submit_button = widgets.Button(
        description='SUBMIT DELFIN + CO2', button_style='success',
        layout=widgets.Layout(width='180px'),
    )
    co2_output = widgets.Output()

    # -- state ----------------------------------------------------------
    state = {
        'converted_xyz_cache': {'smiles': None, 'xyz': None},
        'current_xyz_for_copy': {'content': None},
        'smiles_preview_index': 0,
        'isomers': [],
        'isomer_index': 0,
        'convert_quick': False,
    }

    # -- handlers -------------------------------------------------------
    def update_molecule_view(change=None):
        # User manually edited coords -> clear isomer navigation
        state['isomers'] = []
        state['isomer_index'] = 0
        isomer_nav_row.layout.display = 'none'

        with mol_output:
            clear_output()
            raw_input = coords_widget.value.strip()

            if not raw_input:
                print('Please enter XYZ coordinates or SMILES.')
                state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
                state['current_xyz_for_copy'] = {'content': None}
                xyz_copy_btn.disabled = True
                xyz_copy_status.value = ''
                return

            cleaned_data, input_type = clean_input_data(raw_input)

            if input_type == 'smiles':
                state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
                state['current_xyz_for_copy'] = {'content': None}
                xyz_copy_btn.disabled = True
                xyz_copy_status.value = ''
                print("SMILES erkannt. Bitte 'CONVERT SMILES' oder 'CONVERT SMILES + UFF' klicken.")
                return

            state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
            coords = cleaned_data
            lines = [l for l in coords.split('\n') if l.strip()]
            num_atoms = len(lines)
            xyz_data = f'{num_atoms}\nGenerated by widget\n{coords}'
            state['current_xyz_for_copy'] = {'content': xyz_data}
            xyz_copy_btn.disabled = False
            xyz_copy_status.value = '<span style="color:#388e3c;">XYZ ready to copy</span>'
            view = py3Dmol.view(width='100%', height=SUBMIT_MOL_HEIGHT)
            view.addModel(xyz_data, 'xyz')
            apply_molecule_view_style(view)
            view.show()

    def on_xyz_copy(button):
        content = state['current_xyz_for_copy'].get('content')
        if not content:
            xyz_copy_status.value = '<span style="color:#d32f2f;">No XYZ to copy</span>'
            return
        escaped = content.replace('\\', '\\\\').replace("'", "\\'").replace('\n', '\\n')
        js_code = (
            "navigator.clipboard.writeText('" + escaped + "')"
            ".then(() => console.log('Copied XYZ'))"
            ".catch(err => console.error('Copy failed:', err));"
        )
        ctx.run_js(js_code)
        xyz_copy_status.value = '<span style="color:#388e3c;">Copied to clipboard</span>'

    def _show_isomer_at_index(index):
        isomers = state['isomers']
        if not isomers:
            return
        index = index % len(isomers)
        state['isomer_index'] = index
        xyz_string, num_atoms, label = isomers[index]

        # Update navigation label and visibility
        if len(isomers) > 1:
            display_label = label or f'Isomer {index + 1}'
            isomer_label.value = (
                f'<span style="font-size:13px;">'
                f'{display_label} ({index + 1}/{len(isomers)})</span>'
            )
            isomer_nav_row.layout.display = ''
        else:
            isomer_nav_row.layout.display = 'none'

        # Update 3D preview
        xyz_data = f'{num_atoms}\nIsomer: {label}\n{xyz_string}'
        with mol_output:
            clear_output()
            view = py3Dmol.view(width='100%', height=SUBMIT_MOL_HEIGHT)
            view.addModel(xyz_data, 'xyz')
            apply_molecule_view_style(view)
            view.show()

        # Update copy state
        state['current_xyz_for_copy'] = {'content': xyz_data}
        xyz_copy_btn.disabled = False
        xyz_copy_status.value = '<span style="color:#388e3c;">XYZ ready to copy</span>'

        # Keep converted_xyz_cache in sync for submit
        state['converted_xyz_cache']['xyz'] = xyz_string

        # Update coords widget without triggering update_molecule_view
        coords_widget.unobserve(update_molecule_view, names='value')
        coords_widget.value = f'{num_atoms}\nConverted from SMILES (isomer: {label})\n{xyz_string}'
        coords_widget.observe(update_molecule_view, names='value')

    def handle_isomer_prev(button):
        if state['isomers']:
            _show_isomer_at_index(state['isomer_index'] - 1)

    def handle_isomer_next(button):
        if state['isomers']:
            _show_isomer_at_index(state['isomer_index'] + 1)

    def _convert_smiles(*, apply_uff: bool):
        raw_input = coords_widget.value.strip()
        if not raw_input:
            with mol_output:
                clear_output()
                print('Please enter SMILES in the input box.')
            return

        with mol_output:
            clear_output()
            if apply_uff:
                print('Converting SMILES with UFF...')
            else:
                print('Converting SMILES (no UFF)...')

        cleaned_data, input_type = clean_input_data(raw_input)
        if input_type != 'smiles':
            with mol_output:
                clear_output()
                print('Input is not a SMILES string.')
            return

        isomers, error = smiles_to_xyz_isomers(cleaned_data, apply_uff=apply_uff)
        if error or not isomers:
            with mol_output:
                clear_output()
                print(f'SMILES: {cleaned_data}')
                print(f'Fehler: {error or "No isomers generated"}')
            state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
            state['isomers'] = []
            isomer_nav_row.layout.display = 'none'
            return

        state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': isomers[0][0]}
        state['isomers'] = isomers
        _show_isomer_at_index(0)

    def handle_convert_smiles(button):
        if state['convert_quick']:
            # Quick mode: single conformer, fast, no isomer search
            raw_input = coords_widget.value.strip()
            if not raw_input:
                with mol_output:
                    clear_output()
                    print('Please enter SMILES in the input box.')
                return
            cleaned_data, input_type = clean_input_data(raw_input)
            if input_type != 'smiles':
                with mol_output:
                    clear_output()
                    print('Input is not a SMILES string.')
                return
            with mol_output:
                clear_output()
                print('Quick convert (single structure)...')
            xyz_string, num_atoms, _method, error = smiles_to_xyz(cleaned_data, apply_uff=False)
            if error or not xyz_string:
                with mol_output:
                    clear_output()
                    print(f'Error: {error or "Conversion failed"}')
                return
            state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': xyz_string}
            state['isomers'] = [(xyz_string, num_atoms, 'quick')]
            _show_isomer_at_index(0)
        else:
            # Full mode: isomer search (existing behaviour)
            _convert_smiles(apply_uff=False)

        # Toggle mode for next click and update button label
        state['convert_quick'] = not state['convert_quick']
        if state['convert_quick']:
            convert_smiles_button.description = 'CONVERT SMILES ⚡'
            convert_smiles_button.tooltip = 'Next click: quick single structure (click again to switch back)'
        else:
            convert_smiles_button.description = 'CONVERT SMILES'
            convert_smiles_button.tooltip = 'Next click: full isomer search (click again for quick)'

    def handle_convert_smiles_uff(button):
        _convert_smiles(apply_uff=True)

    def handle_build_complex(button):
        with output_area:
            clear_output()
            job_name = job_name_widget.value.strip()
            if not job_name:
                print('Error: Job name is required!')
                return

            raw_input = coords_widget.value.strip()
            if not raw_input:
                print('Error: Please enter a SMILES string in the input box.')
                return

            cleaned_data, input_type = clean_input_data(raw_input)
            if input_type != 'smiles':
                print('Error: Input must be a SMILES string for BUILD COMPLEX.')
                return

            if not contains_metal(cleaned_data):
                print('Error: SMILES does not contain a metal atom.')
                return

            safe_name = ''.join(c for c in job_name if c.isalnum() or c in ('_', '-'))
            job_dir = ctx.calc_dir / safe_name
            if job_dir.exists():
                print(f'Error: Directory already exists: {job_dir}')
                return

            try:
                job_dir.mkdir(parents=True)
            except Exception as e:
                print(f'Error creating directory: {e}')
                return

            input_file = job_dir / 'input.txt'
            input_file.write_text(cleaned_data + '\n')

            time_limit = resolve_time_limit(job_type_widget, custom_time_widget, '48:00:00')

            try:
                result = ctx.backend.submit_delfin(
                    job_dir=job_dir, job_name=safe_name, mode='build',
                    time_limit=time_limit, build_mult=1, pal=32, maxcore=4000,
                )
                if result.returncode == 0:
                    print('Job submitted successfully!')
                    print(result.stdout)
                    print('')
                    print('Check status in Job Status tab')
                    reset_form()
                else:
                    print('Submission failed:')
                    print(result.stderr)
            except Exception as e:
                print(f'Error submitting job: {e}')

    def handle_guppy_submit(button):
        """Submit SMILES->20x XTB sampling job (GUPPY trajectory mode)."""
        with output_area:
            clear_output()
            job_name = job_name_widget.value.strip()
            if not job_name:
                print('Error: Job name is required!')
                return

            raw_input = coords_widget.value.strip()
            if not raw_input:
                print('Error: Please enter a SMILES string in the input box.')
                return

            cleaned_data, input_type = clean_input_data(raw_input)
            if input_type != 'smiles':
                print('Error: Input must be a SMILES string for GUPPY submission.')
                return

            safe_job_name = ''.join(c for c in job_name if c.isalnum() or c in ('_', '-'))
            if not safe_job_name:
                print('Error: Job name contains only invalid characters!')
                return

            job_dir = ctx.calc_dir / safe_job_name
            time_limit = '00:45:00'

            try:
                # Match ONLY GOAT behavior: allow existing dir and reuse same naming flow.
                job_dir.mkdir(parents=True, exist_ok=True)

                # Required input for guppy mode: raw SMILES in input.txt
                (job_dir / 'input.txt').write_text(cleaned_data + '\n')

                # Resource policy must match ONLY GOAT.
                goat_template = ctx.only_goat_template
                if ctx.only_goat_template_path and ctx.only_goat_template_path.exists():
                    goat_template = ctx.only_goat_template_path.read_text()
                pal, maxcore = parse_resource_settings(goat_template)

                result = ctx.backend.submit_delfin(
                    job_dir=job_dir,
                    job_name=safe_job_name,
                    mode='guppy',
                    time_limit=time_limit,
                    pal=pal or 40,
                    maxcore=maxcore or 6000,
                )

                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    print('GUPPY sampling job submitted!')
                    print(f'Job ID: {job_id}')
                    print(f'Time Limit: {time_limit}')
                    print('Workflow: 20x (SMILES -> XTB2 OPT) with energy ranking')
                    print(f'Input Type: {input_type.upper()}')
                    print(f'Directory: {job_dir}')
                    print('')
                    print('Expected output: GUPPY_try.xyz')
                    print('Check status in Job Status tab')
                    reset_form()
                else:
                    print('Error submitting job:')
                    print(result.stderr)
            except Exception as e:
                print(f'Error submitting GUPPY job: {e}')

    def get_smiles_list_entries():
        entries = []
        for line in smiles_batch_widget.value.splitlines():
            line = line.strip()
            if not line or ';' not in line:
                continue
            parts = [p.strip() for p in line.split(';') if p.strip()]
            if len(parts) >= 2:
                name = parts[0]
                smi = parts[1]
                extras = {}
                for part in parts[2:]:
                    if '=' in part:
                        k, v = part.split('=', 1)
                        extras[k.strip()] = v.strip()
                entries.append((name, smi, extras))
        return entries

    def update_smiles_preview_label():
        entries = get_smiles_list_entries()
        total = len(entries)
        if total == 0:
            smiles_preview_label.value = '<span style="font-size:12px;">0 / 0</span>'
        else:
            current = state['smiles_preview_index'] + 1
            smiles_preview_label.value = f'<span style="font-size:12px;">{current} / {total}</span>'

    def preview_smiles_at_index(index):
        entries = get_smiles_list_entries()
        if not entries:
            with mol_output:
                clear_output()
                print('No valid SMILES entries in the batch list.')
            return
        if index < 0:
            index = 0
        if index >= len(entries):
            index = len(entries) - 1

        state['smiles_preview_index'] = index
        update_smiles_preview_label()
        name, smi, extras = entries[index]

        with mol_output:
            clear_output()
            print(f'Preview: {name}')
            print(f'SMILES: {smi}')
            if extras:
                print(f'Options: {extras}')
            print('Converting...')

        xyz_string, num_atoms, method, error = smiles_to_xyz(smi)
        with mol_output:
            clear_output()
            if error:
                print(f'Preview: {name}')
                print(f'SMILES: {smi}')
                if extras:
                    print(f'Options: {extras}')
                print(f'Error: {error}')
                return
            print(f'Preview: {name} ({method})')
            print(f'SMILES: {smi}')
            if extras:
                print(f'Options: {extras}')
            xyz_data = f'{num_atoms}\nPreview: {name}\n{xyz_string}'
            view = py3Dmol.view(width='100%', height=SUBMIT_MOL_HEIGHT)
            view.addModel(xyz_data, 'xyz')
            apply_molecule_view_style(view)
            view.show()

    def handle_smiles_prev(button):
        entries = get_smiles_list_entries()
        if not entries:
            with mol_output:
                clear_output()
                print('No valid SMILES entries in the batch list.')
            return
        new_index = state['smiles_preview_index'] - 1
        if new_index < 0:
            new_index = len(entries) - 1
        preview_smiles_at_index(new_index)

    def handle_smiles_next(button):
        entries = get_smiles_list_entries()
        if not entries:
            with mol_output:
                clear_output()
                print('No valid SMILES entries in the batch list.')
            return
        new_index = state['smiles_preview_index'] + 1
        if new_index >= len(entries):
            new_index = 0
        preview_smiles_at_index(new_index)

    def handle_submit_smiles_list(button):
        with smiles_batch_output:
            clear_output()
            job_prefix = job_name_widget.value.strip()
            if not job_prefix:
                print('Error: Job name cannot be empty!')
                return

            control_content_base = control_widget.value
            control_errors = validate_control_text(control_content_base)
            if control_errors:
                print('CONTROL.txt validation failed:')
                for err in control_errors:
                    print(f'- {err}')
                return

            time_limit = resolve_time_limit(job_type_widget, custom_time_widget, '48:00:00')

            entries = [l.strip() for l in smiles_batch_widget.value.splitlines() if l.strip()]
            if not entries:
                print('Error: SMILES list is empty.')
                return

            for idx, entry in enumerate(entries, 1):
                if ';' not in entry:
                    print(f'Line {idx}: Missing \';\' delimiter -> {entry}')
                    continue
                parts = [p.strip() for p in entry.split(';') if p.strip()]
                if len(parts) < 2:
                    print(f'Line {idx}: Missing name or SMILES -> {entry}')
                    continue

                name_raw = parts[0]
                smi = parts[1]
                extra_parts = parts[2:]

                if not name_raw or not smi:
                    print(f'Line {idx}: Missing name or SMILES -> {entry}')
                    continue

                extras = {}
                for part in extra_parts:
                    if '=' not in part:
                        print(f"Line {idx}: Invalid override '{part}' (expected key=value)")
                        continue
                    key, value = part.split('=', 1)
                    key, value = key.strip(), value.strip()
                    if not key:
                        print(f"Line {idx}: Invalid override '{part}' (empty key)")
                        continue
                    extras[key] = value

                safe_name = ''.join(c for c in name_raw if c.isalnum() or c in ('_', '-'))
                if not safe_name:
                    print(f'Line {idx}: Invalid name -> {name_raw}')
                    continue

                full_job_name = f'{job_prefix}_{safe_name}'
                safe_job_name = ''.join(c for c in full_job_name if c.isalnum() or c in ('_', '-'))
                if not safe_job_name:
                    print(f'Line {idx}: Invalid job name -> {full_job_name}')
                    continue

                xyz_string, num_atoms, method, error = smiles_to_xyz(smi)
                if error:
                    print(f'Line {idx}: {safe_name} - SMILES error: {error}')
                    continue

                job_dir = ctx.calc_dir / safe_job_name
                job_dir.mkdir(parents=True, exist_ok=True)

                control_content = re.sub(
                    r'(?m)^SMILES=.*$', f'SMILES={smi}', control_content_base,
                )
                for key, value in extras.items():
                    pattern = rf'(?m)^{re.escape(key)}\s*=.*$'
                    replacement = f'{key}={value}'
                    if re.search(pattern, control_content):
                        control_content = re.sub(pattern, replacement, control_content)
                    else:
                        control_content = control_content.rstrip() + f'\n{replacement}\n'

                (job_dir / 'CONTROL.txt').write_text(control_content)
                (job_dir / 'input.txt').write_text(xyz_string)

                pal, maxcore = parse_resource_settings(control_content)
                result = ctx.backend.submit_delfin(
                    job_dir=job_dir, job_name=safe_job_name, mode='delfin',
                    time_limit=time_limit, pal=pal or 40, maxcore=maxcore or 6000,
                )

                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    print(f'Submitted {safe_job_name} (ID: {job_id})')
                else:
                    print(f'Failed {safe_job_name}: {result.stderr or result.stdout}')

    def reset_form():
        job_name_widget.value = ''
        coords_widget.value = ''
        smiles_batch_widget.value = ''
        control_widget.value = ctx.default_control
        job_type_widget.value = '48h'
        custom_time_widget.value = 72
        only_goat_charge.value = 0
        only_goat_solvent.value = ''
        co2_species_delta.value = -2
        state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
        state['isomers'] = []
        state['isomer_index'] = 0
        isomer_nav_row.layout.display = 'none'
        with mol_output:
            clear_output()
            print('Please enter XYZ coordinates or SMILES.')

    def handle_submit(button):
        with output_area:
            clear_output()
            job_name = job_name_widget.value.strip()
            control_content = control_widget.value
            raw_input = coords_widget.value.strip()

            if not job_name:
                print('Error: Job name cannot be empty!')
                return
            if not raw_input:
                print('Error: Input (coordinates or SMILES) cannot be empty!')
                return

            control_errors = validate_control_text(control_content)
            if control_errors:
                print('CONTROL.txt validation failed:')
                for err in control_errors:
                    print(f'- {err}')
                return

            input_content, input_type = clean_input_data(raw_input)
            cache = state['converted_xyz_cache']
            if input_type == 'smiles' and cache.get('xyz'):
                input_content = cache['xyz']
                input_type = 'xyz (from SMILES)'

            if not input_content:
                print('Error: No valid input found!')
                return

            safe_job_name = ''.join(c for c in job_name if c.isalnum() or c in ('_', '-'))
            if not safe_job_name:
                print('Error: Job name contains only invalid characters!')
                return

            job_dir = ctx.calc_dir / safe_job_name
            time_limit = resolve_time_limit(job_type_widget, custom_time_widget, '48:00:00')

            try:
                job_dir.mkdir(parents=True, exist_ok=True)

                if cache.get('smiles'):
                    control_content = control_content.replace(
                        'SMILES=', f"SMILES={cache['smiles']}",
                    )

                (job_dir / 'CONTROL.txt').write_text(control_content)
                (job_dir / 'input.txt').write_text(input_content)

                pal, maxcore = parse_resource_settings(control_content)
                result = ctx.backend.submit_delfin(
                    job_dir=job_dir, job_name=safe_job_name, mode='delfin',
                    time_limit=time_limit, pal=pal or 40, maxcore=maxcore or 6000,
                )

                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    print('Job successfully submitted!')
                    print(f'Job ID: {job_id}')
                    print(f'Time Limit: {time_limit}')
                    print(f'Input Type: {input_type.upper()}')
                    print(f'Directory: {job_dir}')
                    print('')
                    print('Check status in Job Status tab')
                    reset_form()
                else:
                    print('Error submitting job:')
                    print(result.stderr)
            except Exception as e:
                print(f'Error creating job: {e}')

    def handle_validate_control(button):
        with validate_output:
            clear_output()
            errors = validate_control_text(control_widget.value)
            if errors:
                print('CONTROL.txt validation failed:')
                for err in errors:
                    print(f'- {err}')
            else:
                print('CONTROL.txt looks valid.')

    def handle_only_goat_submit(button):
        with only_goat_output:
            clear_output()
            job_name = job_name_widget.value.strip()
            raw_input = coords_widget.value.strip()
            charge_value = only_goat_charge.value
            solvent_value = only_goat_solvent.value.strip()

            if not job_name:
                print('Error: Job name cannot be empty!')
                return
            if not raw_input:
                print('Error: Input (coordinates or SMILES) cannot be empty!')
                return
            if not solvent_value:
                print('Error: Solvent cannot be empty!')
                return

            input_content, input_type = clean_input_data(raw_input)
            cache = state['converted_xyz_cache']
            if input_type == 'smiles' and cache.get('xyz'):
                input_content = cache['xyz']
                input_type = 'xyz (from SMILES)'

            if not input_content:
                print('Error: No valid input found!')
                return

            safe_job_name = ''.join(c for c in job_name if c.isalnum() or c in ('_', '-'))
            if not safe_job_name:
                print('Error: Job name contains only invalid characters!')
                return

            job_dir = ctx.calc_dir / safe_job_name
            time_limit = resolve_time_limit(job_type_widget, custom_time_widget, '48:00:00')

            try:
                job_dir.mkdir(parents=True, exist_ok=True)

                # Use template from file (BwUni) or inline (local)
                template = ctx.only_goat_template
                if ctx.only_goat_template_path and ctx.only_goat_template_path.exists():
                    template = ctx.only_goat_template_path.read_text()

                control_content = (
                    template
                    .replace('[CHARGE]', str(charge_value))
                    .replace('[SOLVENT]', solvent_value)
                )

                if cache.get('smiles'):
                    control_content = control_content.replace(
                        'SMILES=', f"SMILES={cache['smiles']}",
                    )

                control_errors = validate_control_text(control_content)
                if control_errors:
                    print('CONTROL.txt validation failed:')
                    for err in control_errors:
                        print(f'- {err}')
                    return

                (job_dir / 'CONTROL.txt').write_text(control_content)
                (job_dir / 'input.txt').write_text(input_content)

                pal, maxcore = parse_resource_settings(control_content)
                result = ctx.backend.submit_delfin(
                    job_dir=job_dir, job_name=safe_job_name, mode='delfin',
                    time_limit=time_limit, pal=pal or 40, maxcore=maxcore or 6000,
                )

                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    print('Only GOAT job successfully submitted!')
                    print(f'Job ID: {job_id}')
                    print(f'Time Limit: {time_limit}')
                    print(f'Input Type: {input_type.upper()}')
                    print(f'Charge: {charge_value}')
                    print(f'Solvent: {solvent_value}')
                    print(f'Directory: {job_dir}')
                    print('')
                    print('Check status in Job Status tab')
                    reset_form()
                else:
                    print('Error submitting job:')
                    print(result.stderr)
            except Exception as e:
                print(f'Error creating job: {e}')

    def handle_co2_chain_submit(button):
        with co2_output:
            clear_output()
            job_name = job_name_widget.value.strip()
            raw_input = coords_widget.value.strip()
            control_content = control_widget.value
            delta = co2_species_delta.value

            if not job_name:
                print('Error: Job name cannot be empty!')
                return
            if not raw_input:
                print('Error: Input (coordinates or SMILES) cannot be empty!')
                return

            control_errors = validate_control_text(control_content)
            if control_errors:
                print('CONTROL.txt validation failed:')
                for err in control_errors:
                    print(f'- {err}')
                return

            input_content, input_type = clean_input_data(raw_input)
            cache = state['converted_xyz_cache']
            if input_type == 'smiles' and cache.get('xyz'):
                input_content = cache['xyz']

            if not input_content:
                print('Error: No valid input found!')
                return

            safe_job_name = ''.join(c for c in job_name if c.isalnum() or c in ('_', '-'))
            if not safe_job_name:
                print('Error: Job name contains only invalid characters!')
                return

            job_dir = ctx.calc_dir / safe_job_name
            time_limit = resolve_time_limit(job_type_widget, custom_time_widget, '48:00:00')

            try:
                job_dir.mkdir(parents=True, exist_ok=True)

                if cache.get('smiles'):
                    control_content = control_content.replace(
                        'SMILES=', f"SMILES={cache['smiles']}",
                    )

                (job_dir / 'CONTROL.txt').write_text(control_content)
                (job_dir / 'input.txt').write_text(input_content)

                pal, maxcore = parse_resource_settings(control_content)
                result = ctx.backend.submit_delfin(
                    job_dir=job_dir, job_name=safe_job_name,
                    mode='delfin-co2-chain',
                    time_limit=time_limit, pal=pal or 40, maxcore=maxcore or 6000,
                    co2_species_delta=delta,
                )

                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    print('DELFIN + CO2 chain job submitted!')
                    print(f'Job ID: {job_id}')
                    print(f'Time Limit: {time_limit}')
                    print(f'CO2 Species Delta: {delta}')
                    print(f'Directory: {job_dir}')
                    print('')
                    print('Check status in Job Status tab')
                    reset_form()
                else:
                    print('Error submitting job:')
                    print(result.stderr)
            except Exception as e:
                print(f'Error creating job: {e}')

    # -- wiring ---------------------------------------------------------
    xyz_copy_btn.on_click(on_xyz_copy)
    coords_widget.observe(update_molecule_view, names='value')
    convert_smiles_button.on_click(handle_convert_smiles)
    convert_smiles_uff_button.on_click(handle_convert_smiles_uff)
    build_complex_button.on_click(handle_build_complex)
    guppy_submit_button.on_click(handle_guppy_submit)
    submit_smiles_list_button.on_click(handle_submit_smiles_list)
    smiles_prev_button.on_click(handle_smiles_prev)
    smiles_next_button.on_click(handle_smiles_next)
    isomer_prev_btn.on_click(handle_isomer_prev)
    isomer_next_btn.on_click(handle_isomer_next)
    only_goat_submit_button.on_click(handle_only_goat_submit)
    co2_submit_button.on_click(handle_co2_chain_submit)
    validate_button.on_click(handle_validate_control)
    submit_button.on_click(handle_submit)

    # -- layout ---------------------------------------------------------
    spacer = widgets.Label(value='', layout=widgets.Layout(height='10px'))
    spacer_large = widgets.Label(value='', layout=widgets.Layout(height='20px'))

    submit_left = widgets.VBox([
        job_name_widget, spacer,
        job_type_widget, custom_time_widget, spacer_large,
        widgets.HTML('<b>Input (XYZ or SMILES):</b>'), coords_widget, spacer,
        widgets.HBox([convert_smiles_button, convert_smiles_uff_button,
                      build_complex_button, guppy_submit_button],
                     layout=widgets.Layout(gap='10px', flex_wrap='wrap')),
        spacer_large,
        widgets.HTML('<b>Batch SMILES:</b>'), smiles_batch_widget, spacer,
        widgets.HBox(
            [smiles_prev_button, smiles_preview_label,
             smiles_next_button, submit_smiles_list_button],
            layout=widgets.Layout(gap='2px', align_items='center', flex_wrap='wrap'),
        ),
        smiles_batch_output,
        spacer_large,
        widgets.HTML('<b>CONTROL.txt:</b>'), control_widget, spacer,
        widgets.HBox([validate_button, submit_button]),
        output_area, validate_output,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='10px',
        box_sizing='border-box', overflow_x='hidden',
    ))

    submit_right = widgets.VBox([
        widgets.HTML('<b>Molecule Preview:</b>'), mol_output,
        isomer_nav_row,
        widgets.HBox([xyz_copy_btn, xyz_copy_status],
                     layout=widgets.Layout(gap='6px', align_items='center', flex_wrap='wrap')),
        spacer_large,
        widgets.HTML('<b>Only GOAT:</b>'),
        widgets.HBox(
            [only_goat_charge, only_goat_solvent, only_goat_submit_button],
            layout=widgets.Layout(gap='8px', flex_wrap='wrap', align_items='center'),
        ),
        only_goat_output,
        spacer_large,
        widgets.HTML('<b>CO2 Coordinator:</b>'),
        widgets.HBox(
            [co2_species_delta, co2_submit_button],
            layout=widgets.Layout(gap='8px', flex_wrap='wrap', align_items='center'),
        ),
        co2_output,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='10px',
        box_sizing='border-box', overflow_x='hidden',
    ))

    tab_widget = widgets.HBox(
        [submit_left, submit_right],
        layout=widgets.Layout(width='100%', align_items='stretch', overflow_x='hidden'),
    )
    submit_css = widgets.HTML(
        """
        <style>
        .submit-split-root, .submit-split-root * {
            box-sizing: border-box;
        }
        .submit-split-root {
            width: 100% !important;
            overflow-x: hidden !important;
        }
        .submit-split-pane {
            min-width: 0 !important;
            overflow-x: hidden !important;
        }
        .submit-split-pane .widget-box,
        .submit-split-pane .widget-hbox,
        .submit-split-pane .widget-vbox {
            max-width: 100% !important;
        }
        .submit-split-pane .widget-output,
        .submit-split-pane .output_area,
        .submit-split-pane .output_subarea,
        .submit-split-pane .output_wrapper,
        .submit-split-pane .jp-OutputArea,
        .submit-split-pane .jp-OutputArea-child,
        .submit-split-pane .jp-OutputArea-output {
            max-width: 100% !important;
            overflow-x: hidden !important;
        }
        </style>
        """
    )
    tab_widget.add_class('submit-split-root')
    submit_left.add_class('submit-split-pane')
    submit_right.add_class('submit-split-pane')
    mol_output.add_class('submit-split-pane')
    tab_widget = widgets.VBox([submit_css, tab_widget], layout=widgets.Layout(width='100%'))

    return tab_widget, {'reset_form': reset_form, 'mol_output': mol_output}
