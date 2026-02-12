"""Submit Job tab: main DELFIN job submission form."""

import re

import py3Dmol
import ipywidgets as widgets
from IPython.display import clear_output

from delfin.config import validate_control_text
from delfin.smiles_converter import contains_metal

from .constants import COMMON_LAYOUT, COMMON_STYLE
from .helpers import resolve_time_limit, create_time_limit_widgets, disable_spellcheck
from .input_processing import (
    smiles_to_xyz, is_smiles, clean_input_data, parse_resource_settings,
)


def create_tab(ctx):
    """Create the Submit Job tab.

    Returns ``(tab_widget, refs_dict)``.
    """
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
        description='Input:',
        layout=widgets.Layout(width='500px', height='200px'),
        style=COMMON_STYLE,
    )

    button_spacer = widgets.Label(value='', layout=widgets.Layout(height='10px'))

    convert_smiles_button = widgets.Button(
        description='CONVERT SMILES', button_style='info',
        layout=widgets.Layout(width='150px'),
    )

    build_complex_button = widgets.Button(
        description='BUILD COMPLEX', button_style='warning',
        layout=widgets.Layout(width='150px'),
    )

    smiles_batch_help = widgets.Label("Batch SMILES list: one per line 'Name;SMILES;key=value;key=value'")
    smiles_batch_widget = widgets.Textarea(
        value='',
        placeholder='name;SMILES;key=value;...\nNi_1;[Ni];charge=2;solvent=water\nCo_1;[Co];charge=3',
        description='SMILES List:',
        layout=widgets.Layout(width='500px', height='160px'),
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
        description='CONTROL.txt:',
        layout=widgets.Layout(width='500px', height='500px'),
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
        border='2px solid #1976d2', width='560px', height='420px', overflow='hidden',
    ))

    xyz_copy_btn = widgets.Button(
        description='\U0001f4cb Copy Coordinates', button_style='success',
        layout=widgets.Layout(width='150px'), disabled=True,
    )
    xyz_copy_status = widgets.HTML(value='', layout=widgets.Layout(margin='0 0 0 6px'))

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

    # -- state ----------------------------------------------------------
    state = {
        'converted_xyz_cache': {'smiles': None, 'xyz': None},
        'current_xyz_for_copy': {'content': None},
        'smiles_preview_index': 0,
    }

    # -- handlers -------------------------------------------------------
    def update_molecule_view(change=None):
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
                print("SMILES erkannt. Bitte 'CONVERT SMILES' klicken.")
                return

            state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
            coords = cleaned_data
            lines = [l for l in coords.split('\n') if l.strip()]
            num_atoms = len(lines)
            xyz_data = f'{num_atoms}\nGenerated by widget\n{coords}'
            state['current_xyz_for_copy'] = {'content': xyz_data}
            xyz_copy_btn.disabled = False
            xyz_copy_status.value = '<span style="color:#388e3c;">XYZ ready to copy</span>'
            view = py3Dmol.view(width=560, height=420)
            view.addModel(xyz_data, 'xyz')
            view.setStyle({}, {'stick': {'radius': 0.15}, 'sphere': {'scale': 0.22}})
            view.zoomTo()
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

    def handle_convert_smiles(button):
        raw_input = coords_widget.value.strip()
        if not raw_input:
            with mol_output:
                clear_output()
                print('Please enter SMILES in the input box.')
            return

        with mol_output:
            clear_output()
            print('Converting SMILES...')

        cleaned_data, input_type = clean_input_data(raw_input)
        if input_type != 'smiles':
            with mol_output:
                clear_output()
                print('Input is not a SMILES string.')
            return

        xyz_string, num_atoms, method, error = smiles_to_xyz(cleaned_data)
        if error:
            with mol_output:
                clear_output()
                print(f'SMILES: {cleaned_data}')
                print(f'Fehler: {error}')
            state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
            return

        state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': xyz_string}
        coords_widget.value = f'{num_atoms}\nConverted from SMILES ({method})\n{xyz_string}'

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
            view = py3Dmol.view(width=560, height=350)
            view.addModel(xyz_data, 'xyz')
            view.setStyle({}, {'stick': {'radius': 0.15}, 'sphere': {'scale': 0.22}})
            view.zoomTo()
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
        state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
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

    # -- wiring ---------------------------------------------------------
    xyz_copy_btn.on_click(on_xyz_copy)
    coords_widget.observe(update_molecule_view, names='value')
    convert_smiles_button.on_click(handle_convert_smiles)
    build_complex_button.on_click(handle_build_complex)
    submit_smiles_list_button.on_click(handle_submit_smiles_list)
    smiles_prev_button.on_click(handle_smiles_prev)
    smiles_next_button.on_click(handle_smiles_next)
    only_goat_submit_button.on_click(handle_only_goat_submit)
    validate_button.on_click(handle_validate_control)
    submit_button.on_click(handle_submit)

    # -- layout ---------------------------------------------------------
    spacer = widgets.Label(value='', layout=widgets.Layout(height='10px'))
    spacer_large = widgets.Label(value='', layout=widgets.Layout(height='20px'))

    submit_left = widgets.VBox([
        job_name_widget, spacer,
        job_type_widget, custom_time_widget, spacer_large,
        widgets.HTML('<b>Input (XYZ or SMILES):</b>'), coords_widget, spacer,
        widgets.HBox([convert_smiles_button, build_complex_button],
                     layout=widgets.Layout(gap='10px')),
        spacer_large,
        widgets.HTML('<b>Batch SMILES:</b>'), smiles_batch_widget, spacer,
        widgets.HBox(
            [smiles_prev_button, smiles_preview_label,
             smiles_next_button, submit_smiles_list_button],
            layout=widgets.Layout(gap='2px', align_items='center'),
        ),
        smiles_batch_output,
        spacer_large,
        widgets.HTML('<b>CONTROL.txt:</b>'), control_widget, spacer,
        widgets.HBox([validate_button, submit_button]),
        output_area, validate_output,
    ], layout=widgets.Layout(width='50%', padding='10px'))

    submit_right = widgets.VBox([
        widgets.HTML('<b>Molecule Preview:</b>'), mol_output,
        widgets.HBox([xyz_copy_btn, xyz_copy_status],
                     layout=widgets.Layout(gap='6px', align_items='center')),
        spacer_large,
        widgets.HTML('<b>Only GOAT:</b>'),
        widgets.HBox([only_goat_charge, only_goat_solvent, only_goat_submit_button]),
        only_goat_output,
    ], layout=widgets.Layout(width='50%', padding='10px'))

    tab_widget = widgets.HBox([submit_left, submit_right])

    return tab_widget, {'reset_form': reset_form, 'mol_output': mol_output}
