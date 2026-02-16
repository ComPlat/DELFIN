"""TURBOMOLE Builder tab: interactive define + job submission (SLURM only)."""

import json
import time as _time
from pathlib import Path

import py3Dmol
import ipywidgets as widgets
from IPython.display import clear_output

from .input_processing import smiles_to_xyz, is_smiles
from .molecule_viewer import coord_to_xyz, apply_molecule_view_style


TM_MODULE_OPTIONS = [
    'ridft', 'dscf', 'ricc2', 'escf', 'grad', 'egrad', 'rdgrad',
    'rigrad', 'jobex', 'jobex-noopt', 'aoforce', 'NumForce', 'freeh',
    'ccsdf12', 'pnoccsd',
]
TM_PARA_OPTIONS = ['SMP', 'MPI']


def _xyz_to_coord(xyz_text, smiles_cache):
    """Convert XYZ (or SMILES) text to TURBOMOLE coord format."""
    raw = xyz_text.strip()
    if not raw:
        smiles_cache.update(smiles=None, xyz=None, method=None)
        return raw, None

    lines = raw.split('\n')
    if lines[0].strip().startswith('$coord'):
        smiles_cache.update(smiles=None, xyz=None, method=None)
        return raw, None

    note = None
    try:
        if is_smiles(raw):
            xyz_string, num_atoms, method, error = smiles_to_xyz(raw)
            if error or not xyz_string:
                smiles_cache.update(smiles=None, xyz=None, method=None)
                return raw, f'SMILES conversion failed: {error}'
            raw = f'{num_atoms}\nConverted from SMILES ({method})\n{xyz_string}'
            lines = raw.split('\n')
            note = f'SMILES converted via {method}'
            smiles_cache.update(smiles=xyz_text.strip(), xyz=xyz_string, method=method)
        else:
            smiles_cache.update(smiles=None, xyz=None, method=None)
    except Exception as e:
        smiles_cache.update(smiles=None, xyz=None, method=None)
        return raw, f'SMILES conversion error: {e}'

    bohr = 1.8897259886
    coord_lines = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 4:
            try:
                elem = parts[0]
                if not elem[0].isalpha():
                    continue
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                coord_lines.append(
                    f'  {x * bohr:14.8f}  {y * bohr:14.8f}  {z * bohr:14.8f}  {elem.lower()}'
                )
            except (ValueError, IndexError):
                continue

    if not coord_lines:
        # Fallback: single-line SMILES that wasn't caught
        if len(lines) == 1 and ' ' not in raw and '\t' not in raw:
            try:
                xyz_string, num_atoms, method, error = smiles_to_xyz(raw)
                if error or not xyz_string:
                    smiles_cache.update(smiles=None, xyz=None, method=None)
                    return raw, f'SMILES conversion failed: {error}'
                coord_lines = []
                for line in xyz_string.split('\n'):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            elem = parts[0]
                            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                            coord_lines.append(
                                f'  {x * bohr:14.8f}  {y * bohr:14.8f}  {z * bohr:14.8f}  {elem.lower()}'
                            )
                        except (ValueError, IndexError):
                            continue
                if coord_lines:
                    smiles_cache.update(smiles=raw, xyz=xyz_string, method=method)
                    return '$coord\n' + '\n'.join(coord_lines) + '\n$end', f'SMILES converted via {method}'
            except Exception as e:
                smiles_cache.update(smiles=None, xyz=None, method=None)
                return raw, f'SMILES conversion error: {e}'
        smiles_cache.update(smiles=None, xyz=None, method=None)
        return raw, note

    return '$coord\n' + '\n'.join(coord_lines) + '\n$end', note


def create_tab(ctx):
    """Create the TURBOMOLE Builder tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    # -- widgets --------------------------------------------------------
    tm_job_name = widgets.Text(
        value='', placeholder='Enter job name (required)',
        description='Job Name:', layout=widgets.Layout(width='90%'),
    )
    tm_coords = widgets.Textarea(
        value='', placeholder='Paste XYZ or TURBOMOLE coord format',
        description='Coordinates:',
        layout=widgets.Layout(width='90%', height='180px', overflow_y='auto'),
    )
    tm_convert_smiles_btn = widgets.Button(
        description='CONVERT SMILES', button_style='info', icon='exchange',
        layout=widgets.Layout(width='150px'),
    )
    tm_module = widgets.Dropdown(
        options=TM_MODULE_OPTIONS, value='ridft',
        description='TM Module:', layout=widgets.Layout(width='200px'),
    )
    tm_para_arch = widgets.Dropdown(
        options=TM_PARA_OPTIONS, value='SMP',
        description='Parallel:', layout=widgets.Layout(width='150px'),
    )
    tm_nprocs = widgets.IntText(value=40, description='CPUs:',
                                layout=widgets.Layout(width='150px'))
    tm_memory = widgets.IntText(value=6000, description='Mem/CPU (MB):',
                                layout=widgets.Layout(width='180px'))
    tm_slurm_time = widgets.Text(value='48:00:00', description='Time Limit:',
                                 layout=widgets.Layout(width='180px'))

    tm_define_output = widgets.Output(
        layout=widgets.Layout(
            width='90%', height='400px', border='1px solid #333',
            overflow='auto', padding='5px', flex='0 0 320px',
        ),
    )
    tm_define_output.add_class('tm-define-output')

    tm_define_input = widgets.Textarea(
        value='', placeholder='Type response, press Enter',
        layout=widgets.Layout(width='90%', height='32px'),
        continuous_update=True,
    )
    tm_define_send = widgets.Button(
        description='Send', button_style='primary',
        layout=widgets.Layout(width='150px'),
    )
    tm_control_label = widgets.HTML(
        '<b>CONTROL (editable):</b>',
        layout=widgets.Layout(margin='8px 0 4px 0'),
    )
    tm_control_area = widgets.Textarea(
        value='', layout=widgets.Layout(width='100%', height='260px', display='none'),
    )
    tm_control_status = widgets.HTML(
        value='', layout=widgets.Layout(width='100%', overflow_x='hidden'),
    )
    tm_start_define_btn = widgets.Button(
        description='Start define', button_style='info', icon='terminal',
        layout=widgets.Layout(width='150px'),
    )
    tm_stop_define_btn = widgets.Button(
        description='Stop define', button_style='danger', icon='stop',
        layout=widgets.Layout(width='150px'), disabled=True,
    )
    tm_save_btn = widgets.Button(
        description='Save to Folder', button_style='warning', icon='save',
        layout=widgets.Layout(width='150px'),
    )
    tm_submit_btn = widgets.Button(
        description='SUBMIT JOB', button_style='success', icon='rocket',
        layout=widgets.Layout(width='150px'),
    )
    tm_status = widgets.HTML(value='')
    tm_load_path = widgets.Text(
        value='', placeholder='Job name or path', description='Load from:',
        layout=widgets.Layout(width='350px'),
    )
    tm_load_btn = widgets.Button(
        description='Load', button_style='info', icon='upload',
        layout=widgets.Layout(width='100px'),
    )
    tm_replay_btn = widgets.Button(
        description='Replay', button_style='warning', icon='play',
        layout=widgets.Layout(width='100px'),
    )
    tm_mol_output = widgets.Output(layout=widgets.Layout(
        width='100%', min_height='300px', border='2px solid #1976d2', overflow='hidden',
    ))

    tm_input_row = widgets.HBox([tm_define_input])

    # -- state ----------------------------------------------------------
    state = {
        'child': None,
        'job_dir': None,
        'running': False,
        'ready': False,
        'input_enabled': False,
        'sending': False,
        'last_edit_at': 0.0,
        'empty_block_until': 0.0,
        'smiles_cache': {'smiles': None, 'xyz': None, 'method': None},
        'pending_replay': [],
        'loaded_control': None,
    }

    # -- helpers --------------------------------------------------------
    def _append_to_json(cmd):
        job_dir = state['job_dir']
        if job_dir is None:
            return
        json_file = job_dir / 'delfin_import.json'
        try:
            if json_file.exists():
                data = json.loads(json_file.read_text())
            else:
                data = {'define_commands': [], 'settings': {}}
            data['define_commands'].append(cmd)
            json_file.write_text(json.dumps(data, indent=2))
        except Exception:
            pass

    def _scroll_to_bottom():
        ctx.run_js('''
            setTimeout(function() {
                var out = document.querySelector('.tm-define-output');
                if (out) { out.scrollTop = out.scrollHeight; }
            }, 100);
        ''')

    def _show_control_editor(job_dir):
        try:
            control_path = Path(job_dir) / 'control'
            tm_control_area.value = control_path.read_text() if control_path.exists() else ''
            tm_control_area.layout.display = 'block'
            tm_control_label.layout.display = 'block'
            tm_control_status.value = ''
            tm_input_row.layout.display = 'none'
            tm_start_define_btn.layout.display = 'none'
            tm_stop_define_btn.layout.display = 'none'
        except Exception as e:
            tm_control_status.value = f'<span style="color:red;">Control load error: {e}</span>'

    def _hide_control_editor():
        tm_control_area.layout.display = 'none'
        tm_control_label.layout.display = 'none'
        tm_control_status.value = ''
        tm_input_row.layout.display = 'block'
        tm_start_define_btn.layout.display = 'inline-flex'
        tm_stop_define_btn.layout.display = 'inline-flex'

    def _update_output():
        import re as _re
        child = state['child']
        if child is None:
            return
        updated = False
        try:
            import pexpect
            while True:
                try:
                    chunk = child.read_nonblocking(size=4096, timeout=0)
                except pexpect.TIMEOUT:
                    break
                except pexpect.EOF:
                    with tm_define_output:
                        print('\n[define finished]')
                    if state['job_dir'] and state['loaded_control']:
                        try:
                            (state['job_dir'] / 'control').write_text(state['loaded_control'])
                        except Exception:
                            pass
                    if state['job_dir']:
                        _show_control_editor(state['job_dir'])
                    state['running'] = False
                    state['ready'] = False
                    state['input_enabled'] = False
                    tm_start_define_btn.disabled = False
                    tm_stop_define_btn.disabled = True
                    tm_define_input.disabled = True
                    tm_define_send.disabled = True
                    break
                if not chunk:
                    break
                with tm_define_output:
                    clean = _re.sub(r'\x1b\[[0-9;]*[a-zA-Z]', '', chunk)
                    clean = clean.replace('\r\n', '\n').replace('\r', '')
                    print(clean, end='', flush=True)
                updated = True
                state['ready'] = True
                if not state['input_enabled']:
                    tm_define_input.disabled = False
                    tm_define_send.disabled = False
                    state['input_enabled'] = True
                    try:
                        ctx.run_js("""
                            setTimeout(function () {
                                var el = document.querySelector('textarea[placeholder="Type response, press Enter"]');
                                if (el) { el.focus(); }
                            }, 50);
                        """)
                    except Exception:
                        pass
            if child is not None and not child.isalive() and state['running']:
                with tm_define_output:
                    print('\n[define finished]')
                if state['job_dir']:
                    _show_control_editor(state['job_dir'])
                state['running'] = False
                state['ready'] = False
                state['input_enabled'] = False
                tm_start_define_btn.disabled = False
                tm_stop_define_btn.disabled = True
                tm_define_input.disabled = True
                tm_define_send.disabled = True
                updated = True
        except Exception as e:
            with tm_define_output:
                print(f'[Output error: {e}]')
        if updated:
            _scroll_to_bottom()

    # -- event handlers -------------------------------------------------
    def update_tm_molecule_view(change=None):
        with tm_mol_output:
            clear_output(wait=True)
            raw = tm_coords.value.strip()
            if not raw:
                print('No coordinates.')
                return
            xyz = coord_to_xyz(raw)
            if not xyz:
                lines = [l for l in raw.split('\n')
                         if len(l.split()) >= 4 and l.split()[0][0:1].isalpha()]
                xyz = f'{len(lines)}\nPreview\n' + '\n'.join(lines) if lines else None
            if not xyz:
                print('Could not parse.')
                return
            try:
                v = py3Dmol.view(width='100%', height=400)
                v.addModel(xyz, 'xyz')
                apply_molecule_view_style(v)
                v.show()
            except Exception as e:
                print(f'Error: {e}')

    def handle_tm_convert_smiles(button=None):
        raw = tm_coords.value.strip()
        if not raw:
            tm_status.value = '<span style="color:red;">Coordinates/SMILES required!</span>'
            return
        if not is_smiles(raw):
            tm_status.value = '<span style="color:orange;">Input is not a SMILES string.</span>'
            return
        tm_status.value = '<span style="color:orange;">Converting SMILES...</span>'
        xyz_string, num_atoms, method, error = smiles_to_xyz(raw)
        if error or not xyz_string:
            tm_status.value = f'<span style="color:red;">SMILES conversion failed: {error}</span>'
            return
        xyz_block = f'{num_atoms}\nConverted from SMILES ({method})\n{xyz_string}'
        coord_text, coord_note = _xyz_to_coord(xyz_block, state['smiles_cache'])
        tm_coords.value = coord_text
        if coord_note:
            tm_status.value = f'<span style="color:green;">{coord_note}</span>'
        else:
            tm_status.value = f'<span style="color:green;">SMILES converted via {method}</span>'

    def handle_start_define(button):
        job_name = tm_job_name.value.strip().lstrip('/')
        if not job_name:
            tm_status.value = '<span style="color:red;">Job name required!</span>'
            return
        if not tm_coords.value.strip():
            tm_status.value = '<span style="color:red;">Coordinates required!</span>'
            return

        _hide_control_editor()
        job_dir = ctx.calc_dir / job_name
        state['job_dir'] = job_dir
        job_dir.mkdir(parents=True, exist_ok=True)

        (job_dir / 'delfin_import.json').write_text(
            json.dumps({'define_commands': [], 'settings': {}}, indent=2)
        )

        coord_text, coord_note = _xyz_to_coord(
            tm_coords.value.strip(), state['smiles_cache'],
        )
        (job_dir / 'coord').write_text(coord_text)

        with tm_define_output:
            clear_output()
            print(f'Job directory: {job_dir}')
            print('Starting define... (this may take a moment)')
            print('=' * 50)

        tm_status.value = f'<span style="color:green;">Dir: {job_dir}</span>'
        if coord_note:
            tm_status.value = f'<span style="color:green;">Dir: {job_dir} \u00b7 {coord_note}</span>'
            with tm_define_output:
                print(coord_note)

        try:
            import pexpect
            cmd = 'source /etc/profile && module load chem/turbomole/7.9 2>&1 && define'
            child = pexpect.spawn(
                '/bin/bash', ['-lc', cmd],
                cwd=str(job_dir), encoding='utf-8', echo=False,
            )
            state['child'] = child
            state['running'] = True
            state['ready'] = False
            state['input_enabled'] = False
            state['empty_block_until'] = _time.time() + 2.0

            tm_start_define_btn.disabled = True
            tm_stop_define_btn.disabled = False
            tm_define_input.disabled = True
            tm_define_send.disabled = True

            for _ in range(20):
                _time.sleep(0.2)
                _update_output()

            def periodic_update():
                if state['running']:
                    _update_output()
                    import asyncio
                    try:
                        loop = asyncio.get_event_loop()
                        loop.call_later(0.2, periodic_update)
                    except Exception:
                        pass

            import asyncio
            try:
                loop = asyncio.get_event_loop()
                loop.call_later(0.5, periodic_update)
            except Exception:
                pass

        except Exception as e:
            tm_status.value = f'<span style="color:red;">Error: {e}</span>'
            with tm_define_output:
                print(f'Error starting define: {e}')
                import traceback
                traceback.print_exc()

    def handle_stop_define(button):
        state['running'] = False
        state['ready'] = False
        state['input_enabled'] = False
        state['empty_block_until'] = 0.0
        child = state['child']
        if child:
            try:
                child.close(force=True)
            except Exception:
                try:
                    child.terminate(force=True)
                except Exception:
                    pass
        state['child'] = None

        tm_start_define_btn.disabled = False
        tm_stop_define_btn.disabled = True
        tm_define_input.disabled = True
        tm_define_send.disabled = True

        with tm_define_output:
            print('\n[define stopped]')

    def handle_send_input(button=None, text_override=None):
        child = state['child']
        if child is None or (hasattr(child, 'isalive') and not child.isalive()):
            with tm_define_output:
                print('[Process not running]')
            return

        if not state['input_enabled']:
            for _ in range(20):
                _time.sleep(0.1)
                _update_output()
                if state['input_enabled']:
                    break

        if state['sending']:
            return
        state['sending'] = True

        if text_override is None and button is None:
            _time.sleep(0.05)
        inp = tm_define_input.value if text_override is None else text_override
        if text_override is None:
            tm_define_input.value = ''

        try:
            if state['child'] is not None:
                state['child'].sendline(inp)
                with tm_define_output:
                    print(f'>>> {inp}' if inp else '>>> [Enter]')
                _append_to_json(inp)
                for _ in range(40):
                    _time.sleep(0.1)
                    _update_output()
        except Exception as e:
            with tm_define_output:
                print(f'[Send error: {e}]')
        finally:
            state['sending'] = False

    def _on_textarea_change(change):
        if change.get('name') != 'value':
            return
        if state['sending'] or not state['input_enabled']:
            return
        new = change.get('new', '') or ''
        if '\n' not in new:
            return
        cmd = new.split('\n', 1)[0]
        tm_define_input.value = ''
        handle_send_input(text_override=cmd)

    def handle_tm_save(button):
        job_name = tm_job_name.value.strip().lstrip('/')
        if not job_name:
            tm_status.value = '<span style="color:red;">Job name required!</span>'
            return
        job_dir = ctx.calc_dir / job_name
        state['job_dir'] = job_dir
        if not (job_dir / 'control').exists():
            tm_status.value = '<span style="color:orange;">No control file. Run define first!</span>'
            return
        if tm_control_area.layout.display != 'none':
            try:
                (job_dir / 'control').write_text(tm_control_area.value)
                tm_control_status.value = '<span style="color:green;">CONTROL updated.</span>'
            except Exception as e:
                tm_control_status.value = f'<span style="color:red;">Control save error: {e}</span>'
        tm_status.value = f'<span style="color:green;">Saved: {job_dir}</span>'

    def handle_tm_submit(button):
        job_name = tm_job_name.value.strip().lstrip('/')
        if not job_name:
            tm_status.value = '<span style="color:red;">Job name required!</span>'
            return
        job_dir = ctx.calc_dir / job_name
        state['job_dir'] = job_dir
        if not (job_dir / 'control').exists():
            tm_status.value = '<span style="color:red;">No control file. Run define first!</span>'
            return
        if not (job_dir / 'coord').exists():
            tm_status.value = '<span style="color:red;">No coord file!</span>'
            return
        if tm_control_area.layout.display != 'none':
            try:
                (job_dir / 'control').write_text(tm_control_area.value)
                tm_control_status.value = '<span style="color:green;">CONTROL updated.</span>'
            except Exception as e:
                tm_control_status.value = f'<span style="color:red;">Control save error: {e}</span>'

        result = ctx.backend.submit_turbomole(
            job_dir=job_dir, job_name=job_name,
            module=tm_module.value,
            time_limit=tm_slurm_time.value,
            nprocs=tm_nprocs.value,
            mem_per_cpu=tm_memory.value,
            para_arch=tm_para_arch.value,
        )

        if result.returncode == 0:
            # Save settings
            json_file = job_dir / 'delfin_import.json'
            try:
                data = json.loads(json_file.read_text()) if json_file.exists() else {'define_commands': [], 'settings': {}}
                data['settings'] = {
                    'module': tm_module.value,
                    'para_arch': tm_para_arch.value,
                    'nprocs': tm_nprocs.value,
                    'memory_per_core': tm_memory.value,
                    'slurm_time': tm_slurm_time.value,
                }
                json_file.write_text(json.dumps(data, indent=2))
            except Exception:
                pass
            tm_status.value = f'<span style="color:green;">Submitted! {result.stdout.strip()}</span>'
            tm_job_name.value = ''
            tm_coords.value = ''
            with tm_define_output:
                clear_output()
            with tm_mol_output:
                clear_output()
            tm_define_input.value = ''
            tm_define_input.disabled = True
            tm_define_send.disabled = True
            _hide_control_editor()
        else:
            tm_status.value = f'<span style="color:red;">{result.stderr or result.stdout}</span>'

    def handle_tm_load(button):
        if not tm_job_name.value.strip():
            tm_status.value = '<span style="color:red;">Enter Job Name first</span>'
            return
        if not tm_coords.value.strip():
            tm_status.value = '<span style="color:red;">Enter Coordinates first</span>'
            return
        load_path = tm_load_path.value.strip()
        if not load_path:
            tm_status.value = '<span style="color:red;">Enter path to load from</span>'
            return
        load_dir = ctx.calc_dir / load_path if not load_path.startswith('/') else Path(load_path)
        json_file = load_dir / 'delfin_import.json'
        if not json_file.exists():
            tm_status.value = '<span style="color:red;">No delfin_import.json found</span>'
            return
        try:
            data = json.loads(json_file.read_text())
            s = data.get('settings', {})
            control_path = load_dir / 'control'
            state['loaded_control'] = control_path.read_text() if control_path.exists() else None
            if 'module' in s:
                tm_module.value = s['module']
            if 'para_arch' in s:
                tm_para_arch.value = s['para_arch']
            if 'nprocs' in s:
                tm_nprocs.value = s['nprocs']
            if 'memory_per_core' in s:
                tm_memory.value = s['memory_per_core']
            if 'slurm_time' in s:
                tm_slurm_time.value = s['slurm_time']
            state['pending_replay'] = data.get('define_commands', [])
            with tm_define_output:
                clear_output()
                print(f"Loaded {len(state['pending_replay'])} commands. Start Define, then click Replay.")
            tm_status.value = f'<span style="color:green;">Loaded from {load_dir.name}</span>'
            _hide_control_editor()
        except Exception as e:
            tm_status.value = f'<span style="color:red;">{e}</span>'

    def handle_tm_replay(button):
        if not state['pending_replay']:
            tm_status.value = '<span style="color:red;">No commands loaded</span>'
            return
        child = state['child']
        if child is None or not child.isalive():
            tm_status.value = '<span style="color:red;">Start Define first</span>'
            return
        if tm_control_area.layout.display != 'none' and state['job_dir']:
            try:
                (state['job_dir'] / 'control').write_text(tm_control_area.value)
                tm_control_status.value = '<span style="color:green;">CONTROL synced.</span>'
            except Exception as e:
                tm_control_status.value = f'<span style="color:red;">Control sync error: {e}</span>'
                return
        cmds = state['pending_replay'][:]
        state['pending_replay'] = []
        tm_status.value = f'<span style="color:orange;">Replaying {len(cmds)} commands...</span>'

        for cmd in cmds:
            try:
                state['child'].sendline(cmd)
                with tm_define_output:
                    print(f'>>> {cmd}')
                _append_to_json(cmd)
                _time.sleep(3.0)
                _update_output()
            except Exception:
                break
        tm_status.value = '<span style="color:green;">Replay done</span>'

    # -- wiring ---------------------------------------------------------
    import warnings
    warnings.filterwarnings('ignore', message='.*on_submit is deprecated.*',
                            category=DeprecationWarning)

    tm_start_define_btn.on_click(handle_start_define)
    tm_stop_define_btn.on_click(handle_stop_define)
    tm_convert_smiles_btn.on_click(handle_tm_convert_smiles)
    tm_define_send.on_click(handle_send_input)
    tm_define_input.observe(lambda c: state.update(last_edit_at=_time.time())
                            if c.get('name') == 'value' else None, names='value')
    tm_define_input.observe(_on_textarea_change, names='value')
    tm_save_btn.on_click(handle_tm_save)
    tm_submit_btn.on_click(handle_tm_submit)
    tm_load_btn.on_click(handle_tm_load)
    tm_replay_btn.on_click(handle_tm_replay)
    tm_coords.observe(update_tm_molecule_view, names='value')
    tm_define_input.disabled = True
    tm_define_send.disabled = True
    update_tm_molecule_view()

    # -- layout ---------------------------------------------------------
    tm_define_panel = widgets.VBox([
        widgets.HTML('<b>Interactive define Terminal:</b>'),
        tm_define_output,
        tm_control_label, tm_control_area, tm_control_status,
        tm_input_row,
        widgets.HBox([tm_start_define_btn, tm_stop_define_btn],
                     layout=widgets.Layout(gap='10px', margin='5px 0')),
    ], layout=widgets.Layout(flex='0 0 auto'))

    tm_left = widgets.VBox([
        tm_job_name, tm_coords,
        widgets.HBox([tm_convert_smiles_btn], layout=widgets.Layout(width='90%')),
        widgets.HTML('<b>Molecule Preview:</b>'),
        widgets.HBox([tm_mol_output],
                     layout=widgets.Layout(justify_content='center', width='100%')),
    ], layout=widgets.Layout(width='44%', padding='10px', overflow_y='hidden'))

    tm_right = widgets.VBox([
        tm_define_panel,
        widgets.HTML('<hr style="margin: 10px 0;">'),
        widgets.HTML('<b>Job Settings:</b>'),
        widgets.HBox([
            widgets.VBox([tm_module, tm_nprocs]),
            widgets.VBox([tm_para_arch, tm_memory]),
        ], layout=widgets.Layout(gap='10px')),
        widgets.HBox([tm_slurm_time]),
        widgets.HBox([tm_load_path, tm_load_btn, tm_replay_btn],
                     layout=widgets.Layout(gap='5px', margin='5px 0')),
        widgets.HBox([tm_save_btn, tm_submit_btn],
                     layout=widgets.Layout(gap='10px', margin='10px 0')),
        tm_status,
    ], layout=widgets.Layout(width='56%', padding='10px'))

    tab_widget = widgets.VBox([
        widgets.HTML('<h3>TURBOMOLE Builder</h3>'),
        widgets.HTML(
            '<a href="https://wiki.bwhpc.de/e/Turbomole" target="_blank" '
            'rel="noopener noreferrer">bwHPC Wiki</a> | '
            '<a href="https://www.turbomole.org/turbomole/turbomole-documentation/" '
            'target="_blank" rel="noopener noreferrer">TURBOMOLE Manuals</a>'
        ),
        widgets.HBox([tm_left, tm_right]),
    ], layout=widgets.Layout(padding='10px'))

    return tab_widget, {}
