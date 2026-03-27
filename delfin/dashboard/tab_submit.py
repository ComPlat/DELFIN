"""Submit Job tab: main DELFIN job submission form."""

import json
import re
import shlex
import threading

import py3Dmol
import ipywidgets as widgets
from IPython.display import clear_output

from delfin.config import parse_control_text, validate_control_text, get_esd_hints
from delfin.smiles_converter import contains_metal

from .constants import COMMON_LAYOUT, COMMON_STYLE
from .helpers import resolve_time_limit, create_time_limit_widgets, disable_spellcheck, parse_time_to_seconds
from .molecule_viewer import apply_molecule_view_style
from .input_processing import (
    smiles_to_xyz, smiles_to_xyz_quick, smiles_to_xyz_quick_with_previews,
    append_hapto_previews_to_isomers,
    smiles_to_xyz_isomers, is_smiles,
    clean_input_data, parse_resource_settings,
)


def _smiles_to_architector_input(smiles):
    """Decompose a metal-complex SMILES into an architector input dict.

    Uses RDKit to parse the SMILES, identify the metal center, split
    ligands, and track coordinating atom indices (coordList) for each
    ligand fragment.

    Returns ``{'core': {...}, 'ligands': [...], 'parameters': {}}``
    or ``None`` on failure.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return None

    from delfin.build_up_complex import _METALS

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None

    # ── Identify metal atoms ──
    metal_indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in _METALS:
            metal_indices.append(atom.GetIdx())
    if not metal_indices:
        return None

    # Use first metal (mononuclear)
    metal_idx = metal_indices[0]
    metal_atom = mol.GetAtomWithIdx(metal_idx)
    metal_symbol = metal_atom.GetSymbol()
    metal_charge = metal_atom.GetFormalCharge()

    # ── Find coordinating atoms and their bonds ──
    coord_atom_to_metal = {}  # old_idx -> metal_idx it was bonded to
    bonds_to_remove = []
    for bond in mol.GetBonds():
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        if a in metal_indices:
            coord_atom_to_metal[b] = a
            bonds_to_remove.append((a, b))
        elif b in metal_indices:
            coord_atom_to_metal[a] = b
            bonds_to_remove.append((a, b))

    # ── Remove metal bonds and atoms, track index mapping ──
    edit_mol = Chem.RWMol(mol)

    # Save original atom properties
    orig_props = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        orig_props[idx] = {
            'formal_charge': atom.GetFormalCharge(),
            'explicit_h': atom.GetNumExplicitHs(),
            'no_implicit': atom.GetNoImplicit(),
            'rad_e': atom.GetNumRadicalElectrons(),
        }

    for a, b in bonds_to_remove:
        edit_mol.RemoveBond(a, b)

    # Build old→new index map before removing atoms
    metal_set = set(metal_indices)
    idx_map = {}
    removed = 0
    for old_idx in range(mol.GetNumAtoms()):
        if old_idx in metal_set:
            removed += 1
        else:
            idx_map[old_idx] = old_idx - removed

    for mi in sorted(metal_indices, reverse=True):
        edit_mol.RemoveAtom(mi)

    mol_no_metal = edit_mol.GetMol()

    # Restore atom properties on non-metal atoms
    for old_idx, props in orig_props.items():
        if old_idx in metal_set:
            continue
        new_idx = idx_map.get(old_idx)
        if new_idx is None:
            continue
        try:
            a = mol_no_metal.GetAtomWithIdx(new_idx)
            a.SetNumRadicalElectrons(props['rad_e'])
            if old_idx in coord_atom_to_metal:
                # Distinguish covalent vs dative M-L bonds using bond
                # order sum to non-metal neighbors + explicit H vs
                # typical valence.
                # If (bo_sum + explicit_h) >= valence → dative donor
                #   (NH₃, H₂O, PPh₃ — valence fully satisfied)
                # If (bo_sum + explicit_h) < valence → covalent, atom
                #   lost an H to bond metal → anionic ([C-], [O-], etc.)
                # Bond ORDER (not count) is critical: aromatic C has 3
                # bonds but bo_sum ≈ 3.0 < typ_val 4 → correctly
                # detected as covalent.
                _TYPICAL_VALENCE = {
                    'C': 4, 'N': 3, 'O': 2, 'S': 2, 'Se': 2,
                    'P': 3, 'As': 3, 'Si': 4, 'B': 3, 'Te': 2,
                    'F': 1, 'Cl': 1, 'Br': 1, 'I': 1,
                }
                sym = a.GetSymbol()
                typ_val = _TYPICAL_VALENCE.get(sym, 2)
                # Sum bond orders to non-metal neighbors only
                bo_sum = 0.0
                for bond in mol.GetBonds():
                    ba, bb = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    if old_idx in (ba, bb):
                        neighbor = bb if ba == old_idx else ba
                        if neighbor not in metal_set:
                            bo_sum += bond.GetBondTypeAsDouble()
                saturated = (bo_sum + props['explicit_h']) >= typ_val

                if props['formal_charge'] != 0:
                    # Already charged (e.g. N+, O+) → keep
                    a.SetFormalCharge(props['formal_charge'])
                    a.SetNumExplicitHs(0)
                    a.SetNoImplicit(True)
                elif saturated:
                    # Dative donor: valence already satisfied → keep
                    # neutral and allow implicit H.
                    a.SetFormalCharge(0)
                else:
                    # Covalent: lost H → anionic
                    a.SetFormalCharge(-1)
                    a.SetNumExplicitHs(0)
                    a.SetNoImplicit(True)
            else:
                a.SetNumExplicitHs(props['explicit_h'])
                a.SetNoImplicit(props['no_implicit'])
        except Exception:
            pass
    try:
        mol_no_metal.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    # ── Split into fragments and track coordList per fragment ──
    frag_atom_lists = Chem.GetMolFrags(mol_no_metal, asMols=False)
    frag_mols = Chem.GetMolFrags(mol_no_metal, asMols=True, sanitizeFrags=False)

    # Map new_idx (in mol_no_metal) back to old_idx
    new_to_old = {v: k for k, v in idx_map.items()}

    ligands = []
    for frag_atoms, frag_mol in zip(frag_atom_lists, frag_mols):
        # frag_atoms: tuple of atom indices in mol_no_metal
        # Find which of these were coordinating atoms
        coord_positions = []
        for pos_in_frag, new_idx in enumerate(frag_atoms):
            old_idx = new_to_old.get(new_idx)
            if old_idx is not None and old_idx in coord_atom_to_metal:
                coord_positions.append(pos_in_frag)

        # Get canonical SMILES and map fragment-local coord positions
        # to atom indices in the canonical mol.
        # Tag each atom with its fragment-local index, then re-parse
        # the canonical SMILES to find where the tagged atoms ended up.
        for atom in frag_mol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)  # 1-based

        mapped_smiles = Chem.MolToSmiles(frag_mol, canonical=True)

        # Parse the mapped SMILES to get old→new index mapping
        mapped_mol = Chem.MolFromSmiles(mapped_smiles, sanitize=False)
        if mapped_mol is None:
            continue

        # Build mapping: fragment-local index → canonical index
        frag_to_canon = {}
        for atom in mapped_mol.GetAtoms():
            orig_frag_idx = atom.GetAtomMapNum() - 1
            frag_to_canon[orig_frag_idx] = atom.GetIdx()

        canon_coord = []
        for pos in coord_positions:
            mapped_idx = frag_to_canon.get(pos)
            if mapped_idx is not None:
                canon_coord.append(mapped_idx)

        if not canon_coord:
            continue

        # Clean SMILES: strip map numbers from the SAME canonical mol
        # to ensure atom ordering matches coordList.
        for atom in mapped_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        canon_smiles = Chem.MolToSmiles(mapped_mol, canonical=False)

        lig_dict = {'smiles': canon_smiles, 'coordList': sorted(canon_coord)}
        ligands.append(lig_dict)

    if not ligands:
        return None

    # ── Compute total charge ──
    total_charge = metal_charge
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in metal_set:
            total_charge += atom.GetFormalCharge()

    return {
        'core': {'metal': metal_symbol},
        'ligands': ligands,
        'parameters': {'full_charge': total_charge},
    }


def create_tab(ctx):
    """Create the Submit Job tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    SUBMIT_MOL_HEIGHT = 650

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
        layout=widgets.Layout(width='155px'),
        tooltip='Full isomer search',
    )
    convert_smiles_quick_button = widgets.Button(
        description='QUICK CONVERT SMILES', button_style='info',
        layout=widgets.Layout(width='185px'),
        tooltip='Fast single structure (no isomer search, no UFF)',
    )
    convert_smiles_uff_button = widgets.Button(
        description='CONVERT SMILES + UFF', button_style='info',
        layout=widgets.Layout(width='185px'),
    )

    build_complex_button = widgets.Button(
        description='BUILD COMPLEX', button_style='warning',
        layout=widgets.Layout(width='150px'),
    )

    architector_button = widgets.Button(
        description='ARCHITECTOR', button_style='warning',
        layout=widgets.Layout(width='150px'),
        tooltip='Convert metal-complex SMILES to 3D via architector',
    )

    smiles_batch_widget = widgets.Textarea(
        value='',
        placeholder=(
            "name;SMILES;key=value;...\n"
            "Ni_1;[Ni];charge=2;solvent=water\n"
            "Co_1;[Co];charge=3\n"
            "\n"
            "name1;source=20/build_complex.xyz;key=value;\n"
            "C  0.0  0.0  0.0\n"
            "H  0.0  0.0  1.0\n"
            "*\n"
            "name2;charge=1;\n"
            "XYZ\n"
            "C  0.0  0.0  0.0\n"
            "*"
        ),
        layout=widgets.Layout(width='100%', height='220px', box_sizing='border-box'),
        style=COMMON_STYLE,
    )

    submit_smiles_list_button = widgets.Button(
        description='SUBMIT BATCH LIST', button_style='success',
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

    guppy_pal = widgets.IntText(
        value=12, description='PAL:',
        style=COMMON_STYLE, layout=widgets.Layout(width='140px'),
    )
    guppy_goat_topk = widgets.Dropdown(
        options=[(str(i), i) for i in range(4)],
        value=0, description='GOAT:',
        style=COMMON_STYLE, layout=widgets.Layout(width='150px'),
    )
    guppy_timeout = widgets.Text(
        value='02:00:00', description='Timeout:',
        style=COMMON_STYLE, layout=widgets.Layout(width='200px'),
    )
    guppy_submit_button = widgets.Button(
        description='SUBMIT GUPPY', button_style='warning',
        layout=widgets.Layout(width='170px'),
    )
    guppy_output = widgets.Output()

    # -- state ----------------------------------------------------------
    state = {
        'converted_xyz_cache': {'smiles': None, 'xyz': None},
        'current_xyz_for_copy': {'content': None},
        'smiles_preview_index': 0,
        'isomers': [],
        'isomer_index': 0,
        'smiles_task_id': 0,
        'batch_preview_task_id': 0,
        'batch_preview_timer': None,
        'batch_preview_cache': {},
    }

    # -- handlers -------------------------------------------------------
    def _set_buttons_disabled(buttons, disabled):
        for button in buttons:
            button.disabled = disabled

    def _set_smiles_conversion_busy(is_busy):
        _set_buttons_disabled(
            [
                convert_smiles_button,
                convert_smiles_quick_button,
                convert_smiles_uff_button,
                isomer_prev_btn,
                isomer_next_btn,
            ],
            is_busy,
        )

    def _set_batch_preview_busy(is_busy):
        _set_buttons_disabled([smiles_prev_button, smiles_next_button], is_busy)

    def _cancel_batch_preview_timer():
        timer = state.get('batch_preview_timer')
        if timer is not None:
            timer.cancel()
            state['batch_preview_timer'] = None

    def update_molecule_view(change=None):
        state['smiles_task_id'] += 1
        _set_smiles_conversion_busy(False)
        # User manually edited coords -> clear isomer navigation and reset convert toggle
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
        text_payload = json.dumps(str(content))
        js_code = (
            "(function(){"
            f"const text={text_payload};"
            "function _manualPrompt(){"
            "try{window.prompt('Copy to clipboard (Cmd+C/Ctrl+C, Enter):', text);}catch(_e){}"
            "}"
            "function _legacyCopy(){"
            "try{"
            "const ta=document.createElement('textarea');"
            "ta.value=text;"
            "ta.setAttribute('readonly','readonly');"
            "ta.style.position='fixed';"
            "ta.style.top='-1000px';"
            "ta.style.left='-1000px';"
            "ta.style.opacity='0';"
            "document.body.appendChild(ta);"
            "ta.focus();"
            "ta.select();"
            "ta.setSelectionRange(0, ta.value.length);"
            "const ok=document.execCommand('copy');"
            "document.body.removeChild(ta);"
            "return !!ok;"
            "}catch(_e){return false;}"
            "}"
            "if(navigator.clipboard && navigator.clipboard.writeText){"
            "navigator.clipboard.writeText(text).catch(function(){"
            "if(!_legacyCopy()) _manualPrompt();"
            "});"
            "}else{"
            "if(!_legacyCopy()) _manualPrompt();"
            "}"
            "})();"
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

    def _start_smiles_conversion(*, apply_uff: bool, quick: bool):
        cached_smiles = state['converted_xyz_cache'].get('smiles') if quick else None
        raw_input = (cached_smiles or coords_widget.value).strip()
        if not raw_input:
            with mol_output:
                clear_output()
                print('Please enter SMILES in the input box.')
            return

        cleaned_data, input_type = clean_input_data(raw_input)
        if input_type != 'smiles':
            with mol_output:
                clear_output()
                print('Please enter SMILES in the input box.')
            return

        state['smiles_task_id'] += 1
        task_id = state['smiles_task_id']
        _set_smiles_conversion_busy(True)

        with mol_output:
            clear_output()
            if quick:
                print('Quick convert (single structure)...')
            elif apply_uff:
                print('Converting SMILES with UFF...')
            else:
                print('Converting SMILES (no UFF)...')

        def _worker():
            try:
                if quick:
                    xyz_string, num_atoms, _method, preview_items, error = (
                        smiles_to_xyz_quick_with_previews(cleaned_data)
                    )
                    result = {
                        'error': error,
                        'xyz_string': xyz_string,
                        'num_atoms': num_atoms,
                        'preview_items': preview_items,
                    }
                else:
                    isomers, error = smiles_to_xyz_isomers(
                        cleaned_data,
                        apply_uff=apply_uff,
                        collapse_label_variants=False,
                        include_binding_mode_isomers=True,
                    )
                    if not error and isomers:
                        isomers = append_hapto_previews_to_isomers(
                            isomers,
                            cleaned_data,
                            include_quick=apply_uff,
                        )
                    result = {'error': error, 'isomers': isomers}
            except Exception as exc:
                result = {'error': str(exc)}

            if task_id != state['smiles_task_id']:
                return

            _set_smiles_conversion_busy(False)
            if quick:
                xyz_string = result.get('xyz_string')
                preview_items = result.get('preview_items') or []
                error = result.get('error')
                if error or not xyz_string:
                    with mol_output:
                        clear_output()
                        print(f'Error: {error or "Conversion failed"}')
                    return
                state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': xyz_string}
                state['isomers'] = [(xyz_string, result['num_atoms'], 'quick')] + preview_items
                _show_isomer_at_index(0)
                return

            isomers = result.get('isomers') or []
            error = result.get('error')
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

        threading.Thread(target=_worker, daemon=True).start()

    def _convert_smiles(*, apply_uff: bool):
        _start_smiles_conversion(apply_uff=apply_uff, quick=False)

    def handle_convert_smiles(button):
        _convert_smiles(apply_uff=False)

    def handle_convert_smiles_quick(button):
        _start_smiles_conversion(apply_uff=False, quick=True)

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

    def handle_architector_convert(button):
        """Convert metal-complex SMILES to 3D XYZ via architector."""
        cached_smiles = state['converted_xyz_cache'].get('smiles')
        raw_input = cached_smiles or coords_widget.value.strip()
        if not raw_input:
            with mol_output:
                clear_output()
                print('Please enter a metal-complex SMILES in the input box.')
            return
        cleaned_data, input_type = clean_input_data(raw_input)
        if input_type != 'smiles':
            with mol_output:
                clear_output()
                print('Please enter a SMILES string in the input box.')
            return
        if not contains_metal(cleaned_data):
            with mol_output:
                clear_output()
                print('SMILES does not contain a metal atom.\n'
                      'Architector is for metal complexes — use CONVERT SMILES for organic molecules.')
            return
        with mol_output:
            clear_output()
            print('Running architector...')
        try:
            import importlib.util
            if importlib.util.find_spec('architector') is None:
                with mol_output:
                    clear_output()
                    print('architector is not installed.\n'
                          'Install via: pip install architector\n'
                          'Or use the Install button in Settings → AI Tools.')
                return

            from architector.complex_construction import build_complex_driver
            from architector.io_process_input import inparse
            from architector import io_ptable

            # ── Decompose SMILES into core + ligands via RDKit ──
            input_dict = _smiles_to_architector_input(cleaned_data)
            if input_dict is None:
                with mol_output:
                    clear_output()
                    print('Could not decompose SMILES into metal + ligands.\n'
                          'Try CONVERT SMILES or BUILD COMPLEX instead.')
                return

            with mol_output:
                clear_output()
                metal = input_dict['core'].get('metal', '?')
                n_lig = len(input_dict.get('ligands', []))
                print(f'Running architector: {metal} + {n_lig} ligands...')

            input_dict = inparse(input_dict)
            results = build_complex_driver(input_dict)

            # Retry with scaled radii when no results and has
            # multidentate ligands.
            real_keys = [k for k in results if '_init_only' not in k]
            ligands = input_dict.get('ligands', [])
            max_dent = max((len(l.get('coordList', [])) for l in ligands), default=0)
            if not real_keys and max_dent >= 2:
                for larger in (True, False):
                    scaled = io_ptable.map_metal_radii(
                        inparse(input_dict), larger=larger,
                    )
                    extra = build_complex_driver(scaled)
                    suffix = '_larger_scaled' if larger else '_smaller_scaled'
                    for k, v in extra.items():
                        results[k + suffix] = v
                    if any('_init_only' not in k for k in extra):
                        break

            if not results:
                with mol_output:
                    clear_output()
                    print('Architector returned no structures for this SMILES.')
                return

            isomers = []
            for i, (key, mol) in enumerate(results.items()):
                if '_init_only' in key:
                    continue
                atoms = mol.get('ase_atoms')
                if atoms is None:
                    continue

                xyz_lines = []
                for atom, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
                    xyz_lines.append(f'{atom}  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}')
                xyz_string = '\n'.join(xyz_lines)
                num_atoms = len(atoms)
                energy = mol.get('energy', None)
                e_str = f' ({energy:.4f} Ha)' if energy is not None else ''
                label = f'{key}{e_str}'
                isomers.append((xyz_string, num_atoms, label))

            if not isomers:
                with mol_output:
                    clear_output()
                    print('Architector could not produce valid 3D structures.')
                return

            state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': isomers[0][0]}
            state['isomers'] = isomers
            _show_isomer_at_index(0)
        except Exception as exc:
            with mol_output:
                clear_output()
                print(f'Architector error: {exc}')

    def handle_guppy_submit(button):
        """Submit GUPPY either for batch entries or a single SMILES input."""
        with guppy_output:
            clear_output()
            job_name = job_name_widget.value.strip()
            pal_value = int(guppy_pal.value or 0)
            goat_topk_value = int(guppy_goat_topk.value or 0)
            timeout_value = str(guppy_timeout.value or '').strip() or '02:00:00'

            if pal_value <= 0:
                print('Error: GUPPY PAL must be > 0.')
                return
            if goat_topk_value not in (0, 1, 2, 3):
                print('Error: GUPPY GOAT must be 0, 1, 2, or 3.')
                return
            try:
                parse_time_to_seconds(timeout_value)
            except Exception:
                print('Error: GUPPY timeout must use HH:MM:SS, e.g. 02:00:00.')
                return

            goat_template = ctx.only_goat_template
            if ctx.only_goat_template_path and ctx.only_goat_template_path.exists():
                goat_template = ctx.only_goat_template_path.read_text()
            guppy_maxcore = 500
            guppy_runs = 20
            guppy_parallel_jobs = 4
            guppy_goat_parallel_jobs = guppy_parallel_jobs
            guppy_cli_command = shlex.join([
                'python',
                '-m',
                'delfin.guppy_sampling',
                'input.txt',
                '--runs',
                str(guppy_runs),
                '--pal',
                str(pal_value),
                '--maxcore',
                str(guppy_maxcore),
                '--parallel-jobs',
                str(guppy_parallel_jobs),
                '--goat-topk',
                str(goat_topk_value),
                '--goat-parallel-jobs',
                str(guppy_goat_parallel_jobs),
                '--output',
                'GUPPY_try.xyz',
            ])
            guppy_env = {
                'GUPPY_RUNS': str(guppy_runs),
                'GUPPY_PAL': str(pal_value),
                'GUPPY_MAXCORE': str(guppy_maxcore),
                'GUPPY_PARALLEL_JOBS': str(guppy_parallel_jobs),
                'GUPPY_GOAT_TOPK': str(goat_topk_value),
                'GUPPY_GOAT_PARALLEL_JOBS': str(guppy_goat_parallel_jobs),
            }

            batch_text = smiles_batch_widget.value.strip()
            if batch_text:
                if not job_name:
                    print('Error: Job name is required for GUPPY batch submission.')
                    return
                entries, parse_errors = parse_batch_entries()
                for err in parse_errors:
                    print(err)
                if parse_errors:
                    return
                if not entries:
                    print('Error: No valid batch entries.')
                    return

                submitted = 0
                for entry in entries:
                    line_no = entry.get('line_no', '?')
                    if entry.get('input_kind') != 'smiles':
                        print(f"Line {line_no}: GUPPY batch supports only SMILES entries.")
                        continue

                    name_raw = entry.get('name', '').strip()
                    smiles_value = str(entry.get('input_raw') or '').strip()
                    safe_name = ''.join(c for c in name_raw if c.isalnum() or c in ('_', '-'))
                    if not safe_name:
                        print(f'Line {line_no}: Invalid name -> {name_raw}')
                        continue
                    if not smiles_value:
                        print(f'Line {line_no}: Missing SMILES payload for {safe_name}')
                        continue

                    full_job_name = f'{job_name}_{safe_name}'
                    safe_job_name = ''.join(c for c in full_job_name if c.isalnum() or c in ('_', '-'))
                    if not safe_job_name:
                        print(f'Line {line_no}: Invalid job name -> {full_job_name}')
                        continue

                    job_dir = ctx.calc_dir / safe_job_name
                    try:
                        job_dir.mkdir(parents=True, exist_ok=True)
                        (job_dir / 'input.txt').write_text(smiles_value + '\n')
                        (job_dir / 'guppy_settings.json').write_text(
                            json.dumps(
                                {
                                    'mode': 'guppy',
                                    'submit_mode': 'batch',
                                    'runs': guppy_runs,
                                    'pal': pal_value,
                                    'maxcore': guppy_maxcore,
                                    'parallel_jobs': guppy_parallel_jobs,
                                    'goat_topk': goat_topk_value,
                                    'goat_parallel_jobs': guppy_goat_parallel_jobs,
                                    'time_limit': timeout_value,
                                    'line_no': line_no,
                                    'batch_name': safe_name,
                                    'job_name': safe_job_name,
                                    'input_file': 'input.txt',
                                    'output_file': 'GUPPY_try.xyz',
                                    'cli_command': guppy_cli_command,
                                },
                                indent=2,
                            ),
                            encoding='utf-8',
                        )

                        result = ctx.backend.submit_delfin(
                            job_dir=job_dir,
                            job_name=safe_job_name,
                            mode='guppy',
                            time_limit=timeout_value,
                            pal=pal_value,
                            maxcore=guppy_maxcore,
                            extra_env=guppy_env,
                        )
                    except Exception as exc:
                        print(f'Failed {safe_job_name}: {exc}')
                        continue

                    if result.returncode == 0:
                        job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                        print(
                            f'Submitted {safe_job_name} [SMILES] '
                            f'(ID: {job_id}, PAL: {pal_value}, GOAT: {goat_topk_value}, Timeout: {timeout_value})'
                        )
                        submitted += 1
                    else:
                        print(f'Failed {safe_job_name}: {result.stderr or result.stdout}')

                if submitted:
                    print('')
                    print('Check status in Job Status tab')
                return

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

            try:
                # Match ONLY GOAT behavior: allow existing dir and reuse same naming flow.
                job_dir.mkdir(parents=True, exist_ok=True)

                # Required input for guppy mode: raw SMILES in input.txt
                (job_dir / 'input.txt').write_text(cleaned_data + '\n')
                (job_dir / 'guppy_settings.json').write_text(
                    json.dumps(
                        {
                            'mode': 'guppy',
                            'submit_mode': 'single',
                            'runs': guppy_runs,
                            'pal': pal_value,
                            'maxcore': guppy_maxcore,
                            'parallel_jobs': guppy_parallel_jobs,
                            'goat_topk': goat_topk_value,
                            'goat_parallel_jobs': guppy_goat_parallel_jobs,
                            'time_limit': timeout_value,
                            'job_name': safe_job_name,
                            'input_file': 'input.txt',
                            'output_file': 'GUPPY_try.xyz',
                            'cli_command': guppy_cli_command,
                        },
                        indent=2,
                    ),
                    encoding='utf-8',
                )

                result = ctx.backend.submit_delfin(
                    job_dir=job_dir,
                    job_name=safe_job_name,
                    mode='guppy',
                    time_limit=timeout_value,
                    pal=pal_value,
                    maxcore=guppy_maxcore,
                    extra_env=guppy_env,
                )

                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    print('GUPPY sampling job submitted!')
                    print(f'Job ID: {job_id}')
                    print(f'Time Limit: {timeout_value}')
                    print(f'PAL: {pal_value}')
                    print(f'GOAT top-k: {goat_topk_value}')
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

    def _clean_xyz_block(raw_xyz):
        text = (raw_xyz or '').strip()
        if not text:
            return ''
        lines = text.splitlines()
        if len(lines) >= 3:
            try:
                int(lines[0].strip())
                return '\n'.join(lines[2:]).strip()
            except ValueError:
                pass
        return '\n'.join(lines).strip()

    def _get_smiles_charge(smi):
        """Sum formal charges of all atoms in a SMILES string via RDKit."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            if mol is None:
                return 0
            return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        except Exception:
            return 0

    def _needs_smiles_charge(control_content, extras):
        """Return True if charge should be auto-extracted from SMILES.

        Triggers when:
        - extras['charge'] is empty or '[CHARGE]'
        - charge is not in extras AND CONTROL.txt has no charge, charge=, or charge=[CHARGE]
        """
        charge_extra = extras.get('charge', None)
        if charge_extra is not None:
            v = charge_extra.strip()
            return not v or v == '[CHARGE]'
        m = re.search(r'(?m)^charge\s*=\s*(.*)$', control_content, re.IGNORECASE)
        if not m:
            return True
        val = m.group(1).strip()
        return not val or val == '[CHARGE]'

    def _batch_token_looks_like_source(token):
        token = str(token or '').strip()
        if not token:
            return False
        token_lower = token.lower()
        return (
            token_lower.endswith('.xyz')
            or token.startswith('/')
            or token.startswith('./')
            or token.startswith('../')
            or token.startswith('~/')
            or bool(re.match(r'^[A-Za-z]:[\\/]', token))
        )

    def _set_control_value(control_content, key, value):
        pattern = rf'(?m)^{re.escape(key)}\s*=.*$'
        replacement = f'{key}={value}'
        if re.search(pattern, control_content):
            return re.sub(
                pattern,
                lambda _m, repl=replacement: repl,
                control_content,
            )
        return control_content.rstrip() + f'\n{replacement}\n'

    def _set_control_smiles(control_content, smiles):
        return _set_control_value(control_content, 'SMILES', smiles)

    def _prepare_delfin_submit_input(raw_input, cache):
        """Return the exact payload DELFIN should receive in input.txt.

        Dashboard conversions are previews only. DELFIN itself should decide
        whether to run quick SMILES conversion or GUPPY based on the original
        SMILES string that reached input.txt.
        """
        input_content, input_type = clean_input_data(raw_input)
        cached_smiles = str((cache or {}).get('smiles') or '').strip()
        cached_xyz = str((cache or {}).get('xyz') or '').strip()

        if input_type == 'smiles':
            submit_smiles = input_content.strip()
            return submit_smiles + '\n', 'smiles', submit_smiles

        cleaned_xyz = str(input_content or '').strip()
        if input_type == 'xyz' and cached_smiles and cached_xyz and cleaned_xyz == cached_xyz:
            return cached_smiles + '\n', 'smiles', cached_smiles

        return input_content, input_type, None

    def parse_batch_entries():
        """Parse mixed SMILES/XYZ batch textarea.

        Supported formats:
        1) SMILES line:
           Name;SMILES;key=value;...
        2) XYZ block:
           Name;source=20/build_complex.xyz;key=value;...
           XYZ
           <coordinates ...>    # with or without XYZ header lines
           *
           or short form:
           Name;20/build_complex.xyz;
           <coordinates ...>
           *
        """
        entries = []
        errors = []
        lines = smiles_batch_widget.value.splitlines()
        i = 0
        while i < len(lines):
            raw_line = lines[i]
            line_no = i + 1
            line = raw_line.strip()
            i += 1
            if not line:
                continue
            if ';' not in line:
                errors.append(f"Line {line_no}: Missing ';' delimiter -> {line}")
                continue

            parts = [p.strip() for p in line.split(';')]
            name = parts[0] if parts else ''
            header_tokens = parts[1:]
            if not name:
                errors.append(f'Line {line_no}: Missing name -> {line}')
                continue

            extras = {}
            source_path = ''
            header_free_tokens = []
            force_xyz = False
            for token in header_tokens:
                if not token:
                    continue
                if '=' in token:
                    key_candidate = token.split('=', 1)[0].strip()
                    # Valid key: starts with lowercase letter, only word chars
                    # (SMILES atoms like O, C, N, Cl are uppercase → treated as SMILES)
                    if (key_candidate[:1].islower()
                            and key_candidate.replace('_', '').isalnum()):
                        key, value = token.split('=', 1)
                        key = key.strip()
                        value = value.strip()
                        if key in {'source', 'path', 'origin'}:
                            source_path = value
                            force_xyz = True
                            continue
                        extras[key] = value
                        continue
                if token.upper() == 'XYZ':
                    force_xyz = True
                    continue
                header_free_tokens.append(token)

            smiles_payload = None
            if header_free_tokens:
                if not source_path and len(header_free_tokens) == 1 and _batch_token_looks_like_source(header_free_tokens[0]):
                    source_path = header_free_tokens[0]
                    force_xyz = True
                elif len(header_free_tokens) == 1 and not force_xyz:
                    smiles_payload = header_free_tokens[0]
                else:
                    invalid_tokens = ', '.join(repr(token) for token in header_free_tokens)
                    errors.append(f'Line {line_no}: Invalid header token(s): {invalid_tokens}')
                    continue

            if smiles_payload and not force_xyz:
                entries.append({
                    'line_no': line_no,
                    'name': name,
                    'input_kind': 'smiles',
                    'input_raw': smiles_payload,
                    'input_content': None,
                    'extras': extras,
                    'source_path': '',
                    'header_raw': raw_line.rstrip(),
                })
                continue

            # XYZ mode: optional explicit "XYZ" marker line, then block until "*"
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i < len(lines) and lines[i].strip().upper() == 'XYZ':
                i += 1

            xyz_lines = []
            while i < len(lines):
                if lines[i].strip() == '*':
                    i += 1
                    break
                xyz_lines.append(lines[i].rstrip())
                i += 1
            else:
                errors.append(f"Line {line_no}: XYZ block for '{name}' missing terminating '*'")
                continue

            xyz_raw = '\n'.join(xyz_lines).strip()
            xyz_content = _clean_xyz_block(xyz_raw)
            if not xyz_content:
                errors.append(f"Line {line_no}: XYZ block for '{name}' is empty")
                continue

            entries.append({
                'line_no': line_no,
                'name': name,
                'input_kind': 'xyz',
                'input_raw': xyz_raw,
                'input_content': xyz_content,
                'extras': extras,
                'source_path': source_path,
                'header_raw': raw_line.rstrip(),
            })
        return entries, errors

    def get_smiles_list_entries():
        entries, _errors = parse_batch_entries()
        return entries

    def update_smiles_preview_label():
        entries = get_smiles_list_entries()
        total = len(entries)
        if total == 0:
            smiles_preview_label.value = '<span style="font-size:12px;">0 / 0</span>'
        else:
            current = state['smiles_preview_index'] + 1
            smiles_preview_label.value = f'<span style="font-size:12px;">{current} / {total}</span>'

    def _batch_preview_key(entry):
        return (
            entry.get('input_kind'),
            entry.get('name'),
            entry.get('input_raw'),
            entry.get('input_content'),
            tuple(sorted((entry.get('extras') or {}).items())),
            str(entry.get('source_path') or ''),
        )

    def _render_batch_preview(entry, preview_payload):
        name = entry['name']
        extras = entry['extras']
        input_kind = entry['input_kind']
        source_path = str(entry.get('source_path') or '').strip()
        xyz_data = preview_payload.get('xyz_data')
        error = preview_payload.get('error')
        method = preview_payload.get('method')

        with mol_output:
            clear_output()
            if error:
                print(f'Preview: {name}')
                if input_kind == 'smiles':
                    print(f"SMILES: {entry['input_raw']}")
                elif source_path:
                    print(f'Source: {source_path}')
                if extras:
                    print(f'Options: {extras}')
                print(f'Error: {error}')
                return

            if input_kind == 'smiles':
                print(f'Preview: {name} ({method})')
                print(f"SMILES: {entry['input_raw']}")
            else:
                print(f'Preview: {name} (XYZ)')
                if source_path:
                    print(f'Source: {source_path}')
            if extras:
                print(f'Options: {extras}')
            view = py3Dmol.view(width='100%', height=SUBMIT_MOL_HEIGHT)
            view.addModel(xyz_data, 'xyz')
            apply_molecule_view_style(view)
            view.show()

    def preview_smiles_at_index(index, *, delay=0.0):
        entries = get_smiles_list_entries()
        if not entries:
            _cancel_batch_preview_timer()
            state['batch_preview_task_id'] += 1
            _set_batch_preview_busy(False)
            update_smiles_preview_label()
            with mol_output:
                clear_output()
                print('No valid batch entries.')
            return
        if index < 0:
            index = 0
        if index >= len(entries):
            index = len(entries) - 1

        state['smiles_preview_index'] = index
        update_smiles_preview_label()
        entry = entries[index]
        state['batch_preview_task_id'] += 1
        task_id = state['batch_preview_task_id']
        _cancel_batch_preview_timer()
        _set_batch_preview_busy(True)

        with mol_output:
            clear_output()
            print(f'Preview: {entry["name"]} [{entry["input_kind"].upper()}]')
            if entry['input_kind'] == 'smiles':
                print(f"SMILES: {entry['input_raw']}")
            else:
                print('XYZ coordinates')
                source_path = str(entry.get('source_path') or '').strip()
                if source_path:
                    print(f'Source: {source_path}')
            if entry['extras']:
                print(f'Options: {entry["extras"]}')
            print('Converting...')

        def _worker():
            preview_key = _batch_preview_key(entry)
            preview_payload = state['batch_preview_cache'].get(preview_key)
            if preview_payload is None:
                try:
                    if entry['input_kind'] == 'smiles':
                        smi = entry['input_raw']
                        xyz_string, num_atoms, method, error = smiles_to_xyz_quick(smi)
                        if error:
                            preview_payload = {'error': error}
                        else:
                            preview_payload = {
                                'xyz_data': f'{num_atoms}\nPreview: {entry["name"]}\n{xyz_string}',
                                'method': method,
                            }
                    else:
                        xyz_coords = entry['input_content']
                        coord_lines = [ln for ln in xyz_coords.splitlines() if ln.strip()]
                        if not coord_lines:
                            preview_payload = {'error': 'Empty XYZ coordinates'}
                        else:
                            num_atoms = len(coord_lines)
                            preview_payload = {
                                'xyz_data': f'{num_atoms}\nPreview: {entry["name"]}\n{xyz_coords}',
                                'method': 'XYZ',
                            }
                except Exception as exc:
                    preview_payload = {'error': str(exc)}
                state['batch_preview_cache'][preview_key] = preview_payload

            if task_id != state['batch_preview_task_id']:
                return

            _set_batch_preview_busy(False)
            _render_batch_preview(entry, preview_payload)

        def _launch_worker():
            state['batch_preview_timer'] = None
            threading.Thread(target=_worker, daemon=True).start()

        if delay > 0:
            timer = threading.Timer(delay, _launch_worker)
            timer.daemon = True
            state['batch_preview_timer'] = timer
            timer.start()
        else:
            _launch_worker()

    def handle_smiles_prev(button):
        entries = get_smiles_list_entries()
        if not entries:
            with mol_output:
                clear_output()
                print('No valid batch entries.')
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
                print('No valid batch entries.')
            return
        new_index = state['smiles_preview_index'] + 1
        if new_index >= len(entries):
            new_index = 0
        preview_smiles_at_index(new_index)

    def resolve_co2_submit_mode(control_content: str):
        """Return (mode, species_delta) derived from optional CONTROL CO2 keys."""
        try:
            parsed = parse_control_text(control_content, keep_steps_literal=True)
        except Exception:
            return 'delfin', None

        raw_flag = parsed.get('co2_coordination', 'off')
        flag = str(raw_flag).strip().lower()
        if flag not in {'on', 'yes', 'true', '1'}:
            return 'delfin', None

        raw_delta = parsed.get('co2_species_delta', 0)
        try:
            delta = int(raw_delta)
        except Exception:
            delta = 0
        return 'delfin-co2-chain', delta

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

            entries, parse_errors = parse_batch_entries()
            for err in parse_errors:
                print(err)
            if not entries:
                print('Error: No valid batch entries.')
                return

            for entry in entries:
                line_no = entry.get('line_no', '?')
                name_raw = entry.get('name', '').strip()
                extras = dict(entry.get('extras', {}))
                input_kind = entry.get('input_kind', 'smiles')
                source_path = str(entry.get('source_path') or '').strip()
                safe_name = ''.join(c for c in name_raw if c.isalnum() or c in ('_', '-'))
                if not safe_name:
                    print(f'Line {line_no}: Invalid name -> {name_raw}')
                    continue

                full_job_name = f'{job_prefix}_{safe_name}'
                safe_job_name = ''.join(c for c in full_job_name if c.isalnum() or c in ('_', '-'))
                if not safe_job_name:
                    print(f'Line {line_no}: Invalid job name -> {full_job_name}')
                    continue

                job_dir = ctx.calc_dir / safe_job_name
                job_dir.mkdir(parents=True, exist_ok=True)

                if input_kind == 'smiles':
                    smi = entry.get('input_raw', '').strip()
                    if not smi:
                        print(f'Line {line_no}: Missing SMILES payload for {safe_name}')
                        continue
                    input_content = smi + '\n'
                    control_content = _set_control_smiles(control_content_base, smi)
                else:
                    input_content = (entry.get('input_content') or '').strip()
                    if not input_content:
                        print(f'Line {line_no}: {safe_name} - empty XYZ payload')
                        continue
                    input_content = _clean_xyz_block(input_content)
                    if not input_content:
                        print(f'Line {line_no}: {safe_name} - empty XYZ payload')
                        continue
                    control_content = control_content_base

                # Auto-fill charge from SMILES if needed
                if input_kind == 'smiles' and _needs_smiles_charge(control_content, extras):
                    extras['charge'] = str(_get_smiles_charge(smi))

                for key, value in extras.items():
                    control_content = _set_control_value(control_content, key, value)

                if source_path:
                    control_content = _set_control_value(control_content, 'source', source_path)

                (job_dir / 'CONTROL.txt').write_text(control_content)
                (job_dir / 'input.txt').write_text(input_content)
                batch_metadata = {
                    'line_no': line_no,
                    'name': safe_name,
                    'input_kind': input_kind,
                    'source': source_path or None,
                    'extras': extras,
                    'header_raw': entry.get('header_raw', ''),
                }
                (job_dir / 'input_metadata.json').write_text(
                    json.dumps(batch_metadata, indent=2),
                    encoding='utf-8',
                )

                pal, maxcore = parse_resource_settings(control_content)
                mode, co2_delta = resolve_co2_submit_mode(control_content)
                result = ctx.backend.submit_delfin(
                    job_dir=job_dir, job_name=safe_job_name, mode=mode,
                    time_limit=time_limit, pal=pal or 40, maxcore=maxcore or 6000,
                    co2_species_delta=co2_delta,
                )

                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    if mode == 'delfin-co2-chain':
                        print(
                            f'Submitted {safe_job_name} '
                            f'[{input_kind.upper()}] (ID: {job_id}, CO2 delta: {co2_delta})'
                        )
                    else:
                        print(f'Submitted {safe_job_name} [{input_kind.upper()}] (ID: {job_id})')
                else:
                    print(f'Failed {safe_job_name}: {result.stderr or result.stdout}')

    def reset_form():
        _cancel_batch_preview_timer()
        job_name_widget.value = ''
        coords_widget.value = ''
        smiles_batch_widget.value = ''
        control_widget.value = ctx.default_control
        job_type_widget.value = '48h'
        custom_time_widget.value = 72
        only_goat_charge.value = 0
        only_goat_solvent.value = ''
        co2_species_delta.value = -2
        guppy_pal.value = 12
        guppy_goat_topk.value = 0
        guppy_timeout.value = '02:00:00'
        state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
        state['isomers'] = []
        state['isomer_index'] = 0
        state['smiles_task_id'] += 1
        state['batch_preview_task_id'] += 1
        state['batch_preview_cache'] = {}
        _set_smiles_conversion_busy(False)
        _set_batch_preview_busy(False)
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

            input_content, input_type, submit_smiles = _prepare_delfin_submit_input(
                raw_input,
                state['converted_xyz_cache'],
            )

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

                smiles_for_charge = submit_smiles
                if submit_smiles:
                    control_content = _set_control_smiles(control_content, submit_smiles)

                # Auto-fill charge from SMILES if CONTROL.txt has none/empty/[CHARGE]
                if smiles_for_charge and _needs_smiles_charge(control_content, {}):
                    auto_charge = _get_smiles_charge(smiles_for_charge)
                    charge_pat = r'(?m)^charge\s*=.*$'
                    if re.search(charge_pat, control_content, re.IGNORECASE):
                        control_content = re.sub(
                            charge_pat, f'charge={auto_charge}', control_content,
                            flags=re.IGNORECASE,
                        )
                    else:
                        control_content = control_content.rstrip() + f'\ncharge={auto_charge}\n'

                (job_dir / 'CONTROL.txt').write_text(control_content)
                (job_dir / 'input.txt').write_text(input_content)

                pal, maxcore = parse_resource_settings(control_content)
                mode, co2_delta = resolve_co2_submit_mode(control_content)
                result = ctx.backend.submit_delfin(
                    job_dir=job_dir, job_name=safe_job_name, mode=mode,
                    time_limit=time_limit, pal=pal or 40, maxcore=maxcore or 6000,
                    co2_species_delta=co2_delta,
                )

                if result.returncode == 0:
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    if mode == 'delfin-co2-chain':
                        print('DELFIN + CO2 chain job successfully submitted!')
                        print(f'CO2 Species Delta: {co2_delta}')
                    else:
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
            hints = get_esd_hints(control_widget.value)
            if hints:
                print('ESD hints (non-blocking):')
                for h in hints:
                    print(f'  ℹ {h}')

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

            input_content, input_type, submit_smiles = _prepare_delfin_submit_input(
                raw_input,
                state['converted_xyz_cache'],
            )

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

                if submit_smiles:
                    control_content = _set_control_smiles(control_content, submit_smiles)

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

            input_content, input_type, submit_smiles = _prepare_delfin_submit_input(
                raw_input,
                state['converted_xyz_cache'],
            )

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

                if submit_smiles:
                    control_content = _set_control_smiles(control_content, submit_smiles)

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
    convert_smiles_quick_button.on_click(handle_convert_smiles_quick)
    convert_smiles_uff_button.on_click(handle_convert_smiles_uff)
    build_complex_button.on_click(handle_build_complex)
    architector_button.on_click(handle_architector_convert)
    guppy_submit_button.on_click(handle_guppy_submit)
    submit_smiles_list_button.on_click(handle_submit_smiles_list)
    smiles_prev_button.on_click(handle_smiles_prev)
    smiles_next_button.on_click(handle_smiles_next)

    def on_batch_change(change):
        state['smiles_preview_index'] = 0
        state['batch_preview_cache'] = {}
        preview_smiles_at_index(0, delay=0.35)

    smiles_batch_widget.observe(on_batch_change, names='value')
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
                      convert_smiles_quick_button],
                     layout=widgets.Layout(gap='10px', flex_wrap='wrap')),
        widgets.HBox([build_complex_button, architector_button],
                     layout=widgets.Layout(gap='10px', flex_wrap='wrap')),
        spacer_large,
        widgets.HTML('<b>Batch SMILES/XYZ:</b>'),
        smiles_batch_widget, spacer,
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
        spacer_large,
        widgets.HTML('<b>GUPPY:</b>'),
        widgets.HBox(
            [guppy_pal, guppy_goat_topk, guppy_timeout, guppy_submit_button],
            layout=widgets.Layout(gap='8px', flex_wrap='wrap', align_items='center'),
        ),
        guppy_output,
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
        .submit-mol-output .output_subarea,
        .submit-mol-output .jp-OutputArea-output,
        .submit-mol-output .jp-OutputArea-child,
        .submit-mol-output .output_area {
            padding: 0 !important;
            margin: 0 !important;
        }
        </style>
        """
    )
    tab_widget.add_class('submit-split-root')
    submit_left.add_class('submit-split-pane')
    submit_right.add_class('submit-split-pane')
    mol_output.add_class('submit-split-pane')
    mol_output.add_class('submit-mol-output')
    tab_widget = widgets.VBox([submit_css, tab_widget], layout=widgets.Layout(width='100%'))

    return tab_widget, {'reset_form': reset_form, 'mol_output': mol_output}
