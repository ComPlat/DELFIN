"""Submit Job tab: main DELFIN job submission form."""

import base64
import html
import os
import importlib.resources
import json
import re
import shlex
import shutil
import threading
import time
from pathlib import Path

import ipywidgets as widgets
import py3Dmol
from IPython import get_ipython
from IPython.display import clear_output

from delfin.config import parse_control_text, validate_control_text, get_esd_hints
from delfin.smiles_converter import contains_metal

from .constants import COMMON_LAYOUT, COMMON_STYLE
from .helpers import resolve_time_limit, create_time_limit_widgets, disable_spellcheck, parse_time_to_seconds
from . import gfn_optimize as _gfn
from . import mopac_optimize as _mopac
from . import solvents as _solvents
from . import structure_editor as _structure_editor
from . import ketcher as _ketcher
from . import separate_systems as _separate
from .molecule_viewer import (
    apply_molecule_view_style, structure_viewer_fullscreen_bootstrap_js,
    structure_viewer_fullscreen_css, structure_viewer_fullscreen_kind_js,
    submit_manip_bootstrap_js, submit_manip_version,
)
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


# ---------------------------------------------------------------------------
# MANTA "best version" button: SHIP-31 champion construction config + GFN2
# single-point energy RANKING, then GFN2 geometry-OPT of the top-N ranked
# structures (parallel, laptop-bounded) -> complete manifold + best-possible
# top geometries in the viewer.
# ---------------------------------------------------------------------------
# 2026-07-04 DE-BLOAT (mirrors cli_manta._CHAMPION_FLAGS — CLI + Submit tab ship the identical config;
# see cli_manta for the forensic result: old 29-flag stack 13.7% vs de-bloated 28.6% topology, 2.1x,
# never-worse in every donor class). KAPPA4 folded into the tuple; CONF_ENERGY_RANK removed (net-neg).
# SINGLE SOURCE OF TRUTH (2026-07-21): the dashboard MUST build with the EXACT
# same champion construction as the CLI and the development loop -- "das dashboard
# konstruiert wie wir".  Import cli_manta's champion list instead of a stale
# hand-copied tuple that silently DRIFTED: it was missing DET_CLASSIFY (so the
# dashboard was not even byte-deterministic) + METALLOID_MD_LEN/JOINT_DECLASH/
# DECLASH_METALLOID_LIGAND/SP3C_TET_SEAT/METALLOID_MD_CLAMP(#17)/D8_SQ_ADD(#19)/
# CN6_OH_ADD(#20)/AROM_SEAT(#21).  Now the dashboard and CLI can never diverge.
# Top-N to GFN2-optimize + parallel worker count (laptop-bounded; tune here).








def create_tab(ctx):
    """Create the Submit Job tab.

    Returns ``(tab_widget, refs_dict)``.
    """
    SUBMIT_MOL_HEIGHT = 650
    SMILES_CONVERTER_PLACEHOLDER = '[QUICK|NORMAL|GUPPY|ARCHITECTOR]'
    main_io_loop = getattr(getattr(get_ipython(), 'kernel', None), 'io_loop', None)

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

    # Draw the structure instead of typing it.  Ketcher is a browser
    # application of thirty-odd megabytes, so it is fetched when it is first

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


    only_goat_label = widgets.Label('Only GOAT:')
    only_goat_charge = widgets.IntText(
        value=0, description='Charge:', style=COMMON_STYLE,
        layout=widgets.Layout(width='160px'),
    )
    only_goat_solvent = widgets.Text(
        value='', placeholder='e.g. water', description='Solvent:',
        style=COMMON_STYLE, layout=widgets.Layout(width='220px'),
    )
    only_goat_smiles_converter = widgets.Dropdown(
        options=['QUICK', 'NORMAL', 'GUPPY', 'ARCHITECTOR'],
        value='QUICK', description='Converter:', style=COMMON_STYLE,
        layout=widgets.Layout(width='220px'),
    )
    only_goat_pal = widgets.IntText(
        value=12, description='PAL:', style=COMMON_STYLE,
        layout=widgets.Layout(width='140px'),
    )
    only_goat_maxcore = widgets.IntText(
        value=600, description='MaxCore:', style=COMMON_STYLE,
        layout=widgets.Layout(width='170px'),
    )
    only_goat_time = widgets.Text(
        value='24:00:00', description='Time Limit:', style=COMMON_STYLE,
        layout=widgets.Layout(width='220px'),
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
        value='02:00:00', description='Time Limit:',
        style=COMMON_STYLE, layout=widgets.Layout(width='220px'),
    )
    guppy_submit_button = widgets.Button(
        description='SUBMIT GUPPY', button_style='warning',
        layout=widgets.Layout(width='170px'),
    )
    guppy_output = widgets.Output()

    # ── Fukui submit block ────────────────────────────────────────────
    # Input (SMILES or XYZ) reuses the shared coords_widget from above.
    # Method / basis / dispersion / solvation are read from the shared
    # CONTROL.txt textarea (control_widget) — same plumbing as initial.inp.
    # Pre-OPT is auto-derived from the input type (SMILES → yes, XYZ → no),
    # so the user doesn't need a separate toggle.
    fukui_skip_cubes = widgets.Checkbox(
        value=False, description='Skip cube generation',
        style=COMMON_STYLE, layout=widgets.Layout(width='240px'),
    )
    fukui_pal = widgets.IntText(
        value=4, description='PAL:',
        style=COMMON_STYLE, layout=widgets.Layout(width='140px'),
    )
    fukui_maxcore = widgets.IntText(
        value=3000, description='MaxCore:',
        style=COMMON_STYLE, layout=widgets.Layout(width='170px'),
    )
    fukui_time_limit = widgets.Text(
        value='04:00:00', description='Time Limit:',
        style=COMMON_STYLE, layout=widgets.Layout(width='220px'),
    )
    fukui_submit_button = widgets.Button(
        description='SUBMIT FUKUI', button_style='primary',
        layout=widgets.Layout(width='170px'),
    )
    fukui_output = widgets.Output()

    # -- state ----------------------------------------------------------
    state = {
        # Which room the editor is standing in, for a bug report written from
        # it: the same gesture behaves differently over a tab that keeps one
        # structure and a tab that keeps named blocks of them, so a report
        # that did not say which is a report that has to be guessed at.
        'editor_host': 'Submit',
        'converted_xyz_cache': {'smiles': None, 'xyz': None},
        'current_xyz_for_copy': {'content': None},
        'smiles_preview_index': 0,
        'isomers': [],
        'isomer_index': 0,
        'smiles_task_id': 0,
        'batch_preview_task_id': 0,
        'batch_preview_timer': None,
        'batch_preview_cache': {},
        'smiles_busy': False,
        'batch_preview_busy': False,
        'manip_inflight': False,
        'manip_bootstrap_done': False,
    }

    # -- handlers -------------------------------------------------------
    def _set_buttons_disabled(buttons, disabled):
        for button in buttons:
            button.disabled = disabled

    def _set_batch_preview_busy(is_busy):
        state['batch_preview_busy'] = bool(is_busy)
        _set_buttons_disabled([smiles_prev_button, smiles_next_button], is_busy)
        ctx.set_busy(state['smiles_busy'] or state['batch_preview_busy'])

    def _cancel_batch_preview_timer():
        timer = state.get('batch_preview_timer')
        if timer is not None:
            timer.cancel()
            state['batch_preview_timer'] = None

    def _schedule_ui_update(func, *args, **kwargs):
        if main_io_loop is not None:
            main_io_loop.add_callback(lambda: func(*args, **kwargs))
            return
        func(*args, **kwargs)
    # The structure editor: the viewer, its toolbar and everything that acts
    # on the structure. It is one part rather than this tab's own so that a
    # second tab can have the same editor over its own coordinates -- the ORCA
    # Builder over the block it is showing, and later the TURBOMOLE one.
    _editor = _structure_editor.build(
        ctx,
        state=state,
        coords_widget=coords_widget,
        viewer_height=SUBMIT_MOL_HEIGHT,
        schedule_ui_update=_schedule_ui_update,
        # Defined further down, so they are looked up when they are called.
        update_view=lambda *a, **k: update_molecule_view(*a, **k),
        get_smiles_charge=lambda *a, **k: _get_smiles_charge(*a, **k),
        set_buttons_disabled=lambda *a, **k: _set_buttons_disabled(*a, **k),
    )
    mol_status = _editor.mol_status
    mol_status_fs = _editor.mol_status_fs
    mol_output = _editor.mol_output
    submit_scope_id = _editor.submit_scope_id
    submit_manip_toolbar = _editor.submit_manip_toolbar
    submit_ff_notes = _editor.submit_ff_notes
    submit_scan_plot = _editor.submit_scan_plot
    submit_labels_btn = _editor.submit_labels_btn
    submit_label_size = _editor.submit_label_size
    _set_mol_status = _editor._set_mol_status
    _clear_mol_status = _editor._clear_mol_status
    _set_manip_toolbar_enabled = _editor._set_manip_toolbar_enabled
    _ensure_manip_bootstrap = _editor._ensure_manip_bootstrap
    _run_manip_js = _editor._run_manip_js
    _remember = _editor._remember
    _push_hand_bonds = _editor._push_hand_bonds
    _refresh_constraints = _editor._refresh_constraints
    convert_smiles_button = _editor.convert_smiles_button
    convert_smiles_quick_button = _editor.convert_smiles_quick_button
    convert_smiles_uff_button = _editor.convert_smiles_uff_button
    manta_button = _editor.manta_button
    manta_settings_row = _editor.manta_settings_row
    xyz_copy_btn = _editor.xyz_copy_btn
    xyz_copy_status = _editor.xyz_copy_status
    isomer_nav_row = _editor.isomer_nav_row
    _clear_mol_output = _editor._clear_mol_output
    _replace_mol_output_text = _editor._replace_mol_output_text
    _replace_mol_output_view = _editor._replace_mol_output_view
    _show_mol_busy = _editor._show_mol_busy
    _set_smiles_conversion_busy = _editor._set_smiles_conversion_busy
    _show_isomer_at_index = _editor._show_isomer_at_index
    _clean_xyz_block = _editor._clean_xyz_block
    _apply_smiles_conversion_result = _editor._apply_smiles_conversion_result
    _build_mol_output_bundle = _editor._build_mol_output_bundle
    submit_draw_open_btn = _editor.submit_draw_open_btn
    submit_draw_get_btn = _editor.submit_draw_get_btn
    submit_draw_update_btn = _editor.submit_draw_update_btn
    submit_draw_frame = _editor.submit_draw_frame
    submit_draw_sync = _editor.submit_draw_sync


    def _apply_batch_preview_result(task_id, entry, preview_payload):
        if task_id != state['batch_preview_task_id']:
            return
        _set_batch_preview_busy(False)
        _render_batch_preview(entry, preview_payload)

    def update_molecule_view(change=None):
        # The editor has put a structure in the box that is already on screen
        # -- an isomer it stepped to, a conversion it just took. Drawing it
        # again would throw the camera away and re-perceive its bonds.
        if state.get('editor_quiet'):
            return
        if state.get('manip_inflight'):
            # Change originated from JS-side atom manipulation: the 3Dmol viewer
            # has already reflected the new coordinates. Update state/copy only.
            raw = coords_widget.value.strip()
            if raw:
                cleaned_data, input_type = clean_input_data(raw)
                if input_type != 'smiles':
                    lines = [l for l in cleaned_data.split('\n') if l.strip()]
                    xyz_data = (
                        f'{len(lines)}\nEdited in DELFIN viewer\n{cleaned_data}'
                    )
                    state['current_xyz_for_copy'] = {'content': xyz_data}
                    xyz_copy_btn.disabled = False
                    xyz_copy_status.value = (
                        '<span style="color:#388e3c;">XYZ ready to copy</span>'
                    )
            state['manip_inflight'] = False
            return
        # A genuinely new structure invalidates the remembered bonding. The
        # element sequence alone cannot tell two constitutional isomers apart,
        # and this is the path every real change comes through -- a paste, a
        # conversion, an isomer step, an optimisation result. Drags do not
        # reach here: they take the manip_inflight branch above, which is
        # exactly what keeps a dragged double bond a double bond.
        state['perceived'] = None
        state['perceived_for'] = None
        # Every constraint names atoms by index, which says nothing about a
        # different molecule. An edit is not a different molecule though: it
        # renumbers the one being worked on, and _apply_structure carries
        # everything across that itself, so none of this may be thrown away
        # underneath it.
        if not state.get('structure_edit_inflight'):
            # A structure the editor has not seen: it starts on this one the
            # way it starts on any other.
            _editor.reset_controls()
            state['constraints'] = []
            # A scan armed on one molecule names atoms the next one may not
            # have.  Left armed it threw on the first click after loading a
            # smaller structure, and run it walked a coordinate that was not
            # there at all -- reported as a completed scan, because the run
            # reads only whether xtb answered.
            state['scan_legs'] = []
            state['bond_edits'] = {}
            state['hand_bonds'] = {}
            state['hyb_overrides'] = {}
            state['structure_undo'] = []
            # A different molecule starts a different history, and its first
            # entry is the structure as it arrived -- which is what Reset goes
            # back to and what Undo stops at.
            state['history'] = []
            state['history_seed_pending'] = True
            state['poly_applied'] = None
            state['poly_metal'] = None
            state['poly_assignment'] = None
            # What Reset goes back to.  This is the one moment a structure is
            # the one that arrived rather than one the viewer has been working
            # on, so it is the only place worth remembering.
            state['pristine_coords'] = coords_widget.value
        state['smiles_task_id'] += 1
        _set_smiles_conversion_busy(False)
        # User manually edited coords -> clear isomer navigation and reset convert toggle
        state['isomers'] = []
        state['isomer_index'] = 0
        isomer_nav_row.layout.display = 'none'
        raw_input = coords_widget.value.strip()

        if not raw_input:
            _replace_mol_output_text('Please enter XYZ coordinates or SMILES.')
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
            _replace_mol_output_text(
                "SMILES detected. Click 'CONVERT SMILES' or 'CONVERT SMILES + UFF'."
            )
            return

        state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
        coords = cleaned_data
        lines = [l for l in coords.split('\n') if l.strip()]
        num_atoms = len(lines)
        xyz_data = f'{num_atoms}\nGenerated by widget\n{coords}'
        state['current_xyz_for_copy'] = {'content': xyz_data}
        xyz_copy_btn.disabled = False
        xyz_copy_status.value = '<span style="color:#388e3c;">XYZ ready to copy</span>'
        _replace_mol_output_view(xyz_data)

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
            _replace_mol_output_text('Please enter a metal-complex SMILES in the input box.')
            return
        cleaned_data, input_type = clean_input_data(raw_input)
        if input_type != 'smiles':
            _replace_mol_output_text('Please enter a SMILES string in the input box.')
            return
        if not contains_metal(cleaned_data):
            _replace_mol_output_text(
                'SMILES does not contain a metal atom.',
                'Architector is for metal complexes — use CONVERT SMILES for organic molecules.',
            )
            return
        _clear_mol_output()
        _set_mol_status('Running architector...', spinner=True)
        try:
            import importlib.util
            if importlib.util.find_spec('architector') is None:
                _replace_mol_output_text(
                    'architector is not installed.',
                    'Install via: pip install architector',
                    'Or use the Install button in Settings → AI Tools.',
                )
                return

            from architector.complex_construction import build_complex_driver
            from architector.io_process_input import inparse
            from architector import io_ptable

            # ── Decompose SMILES into core + ligands via RDKit ──
            input_dict = _smiles_to_architector_input(cleaned_data)
            if input_dict is None:
                _replace_mol_output_text(
                    'Could not decompose SMILES into metal + ligands.',
                    'Try CONVERT SMILES or BUILD COMPLEX instead.',
                )
                return

            metal = input_dict['core'].get('metal', '?')
            n_lig = len(input_dict.get('ligands', []))
            _set_mol_status(f'Running architector: {metal} + {n_lig} ligands...', spinner=True)

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
                _replace_mol_output_text('Architector returned no structures for this SMILES.')
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
                _replace_mol_output_text('Architector could not produce valid 3D structures.')
                return

            state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': isomers[0][0]}
            state['isomers'] = isomers
            _show_isomer_at_index(0)
        except Exception as exc:
            _replace_mol_output_text(f'Architector error: {exc}')

    def handle_fukui_submit(button):
        """Submit a Fukui-indices job (3 ORCA SPs) via local or SLURM backend.

        Reads:
          - Input geometry from ``coords_widget`` (SMILES or XYZ, same field
            the rest of the submit tab uses).
          - ORCA method / basis / dispersion / solvation / metal-basis from
            ``control_widget`` (the dashboard's CONTROL.txt textarea).
        """
        with fukui_output:
            clear_output()
            job_name = job_name_widget.value.strip()
            raw_input_value = str(coords_widget.value or '').strip()
            control_text = str(control_widget.value or '').strip()
            time_limit_value = str(fukui_time_limit.value or '').strip() or '04:00:00'
            pal_value = int(fukui_pal.value or 0)
            maxcore_value = int(fukui_maxcore.value or 0)
            skip_cubes_value = bool(fukui_skip_cubes.value)

            if not job_name:
                print('Error: Job name is required!')
                return
            if not raw_input_value:
                print('Error: Paste a SMILES or XYZ block in the Input field (top of the tab).')
                return
            if pal_value <= 0:
                print('Error: Fukui PAL must be > 0.')
                return
            if maxcore_value <= 0:
                print('Error: Fukui MaxCore must be > 0.')
                return
            try:
                parse_time_to_seconds(time_limit_value)
            except Exception:
                print('Error: Fukui Time Limit must use HH:MM:SS, e.g. 04:00:00.')
                return

            # CONTROL.txt is the canonical source for charge / multiplicity /
            # functional / basis / dispersion / solvation / metal_basisset
            # — run the standard DELFIN validator before submitting so
            # the user sees the missing-key list in the dashboard instead
            # of a cryptic ORCA "INPUT ERROR" later.
            if not control_text:
                print(
                    'Error: CONTROL.txt is empty. Paste a CONTROL.txt block in '
                    'the field below the SMILES/XYZ input — Fukui needs charge, '
                    'functional, basis, dispersion, solvation etc. set there.'
                )
                return
            control_errors = validate_control_text(control_text)
            if control_errors:
                print('Error: CONTROL.txt is missing or invalid for Fukui submission:')
                for err in control_errors:
                    print(f'  • {err}')
                print(
                    '\nFix the values above (e.g. replace [SOLVENT] / [METAL] '
                    'placeholders with real values) and submit again.'
                )
                return

            safe_job_name = ''.join(c for c in job_name if c.isalnum() or c in ('_', '-'))
            if not safe_job_name:
                print('Error: Job name contains only invalid characters!')
                return

            cleaned_data, input_type = clean_input_data(raw_input_value)

            job_dir = ctx.calc_dir / safe_job_name
            try:
                job_dir.mkdir(parents=True, exist_ok=True)

                # DELFIN convention: start structure ALWAYS lives in
                # ``input.txt`` (SMILES line or XYZ block). The CLI then
                # calls ``normalize_input_file`` to produce start.txt and
                # invokes the canonical bang-line builder, exactly like
                # the classic pipeline does.
                if input_type == 'smiles':
                    (job_dir / 'input.txt').write_text(cleaned_data + '\n', encoding='utf-8')
                elif input_type == 'xyz':
                    (job_dir / 'input.txt').write_text(cleaned_data + '\n', encoding='utf-8')
                else:
                    print(f'Error: Could not classify input as SMILES or XYZ: {input_type}')
                    return

                # Pre-OPT auto-decided: SMILES inputs need optimisation,
                # XYZ inputs are assumed to be already optimised (either
                # by the user or via the QUICK CONVERT SMILES button which
                # left an XYZ in the field).
                pre_opt_value = (input_type == 'smiles')

                # CONTROL.txt drives functional / basis / dispersion /
                # solvation / metal_basisset / RI / relativity through
                # the canonical _build_bang_line path. Persisting it in
                # job_dir matches how every other DELFIN workflow runs.
                if control_text:
                    (job_dir / 'CONTROL.txt').write_text(control_text + '\n', encoding='utf-8')

                (job_dir / 'fukui_settings.json').write_text(
                    json.dumps(
                        {
                            'mode': 'fukui',
                            'input_type': input_type,
                            'pre_opt': pre_opt_value,
                            'skip_cubes': skip_cubes_value,
                            'pal': pal_value,
                            'maxcore': maxcore_value,
                            'time_limit': time_limit_value,
                        },
                        indent=2,
                    ),
                    encoding='utf-8',
                )

                # The CLI auto-detects input.txt + CONTROL.txt in the
                # job_dir, so no explicit --input is needed here.
                fukui_env = {
                    'DELFIN_FUKUI_PRE_OPT': '1' if pre_opt_value else '0',
                    'DELFIN_FUKUI_SKIP_CUBES': '1' if skip_cubes_value else '0',
                }

                result = ctx.backend.submit_delfin(
                    job_dir=job_dir,
                    job_name=safe_job_name,
                    mode='fukui',
                    time_limit=time_limit_value,
                    pal=pal_value,
                    maxcore=maxcore_value,
                    extra_env=fukui_env,
                )
            except Exception as exc:
                print(f'Failed to submit Fukui job: {exc}')
                return

            if result.returncode == 0:
                stdout = (result.stdout or '').strip()
                job_id = stdout.split()[-1] if stdout else '(unknown)'
                try:  # auto-watch (no-op unless job monitoring enabled)
                    from delfin.agent.job_monitor import auto_watch_if_enabled
                    auto_watch_if_enabled(job_id)
                except Exception:
                    pass
                print(f'Submitted Fukui job {safe_job_name} (ID: {job_id})')
                print(f'  Input type: {input_type}')
                print(f'  Pre-OPT:    {"yes" if pre_opt_value else "no"}')
                print(f'  PAL:        {pal_value}')
                print(f'  MaxCore:    {maxcore_value}')
                print(f'  Time Limit: {time_limit_value}')
                print('')
                print('Check status in Job Status tab.')
                # Clear input fields so the next submit starts fresh.
                # CONTROL.txt and resource sliders (PAL / MaxCore / Time
                # Limit / Pre-OPT / Skip-Cubes) stay as user preferences.
                job_name_widget.value = ''
                coords_widget.value = ''
            else:
                print(f'Submit failed: {result.stderr or result.stdout}')

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

            raw_input = coords_widget.value.strip()
            batch_text = smiles_batch_widget.value.strip()
            if not raw_input and batch_text:
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
                        try:  # auto-watch (no-op unless job monitoring enabled)
                            from delfin.agent.job_monitor import auto_watch_if_enabled
                            auto_watch_if_enabled(job_id)
                        except Exception:
                            pass
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
                    try:  # auto-watch (no-op unless job monitoring enabled)
                        from delfin.agent.job_monitor import auto_watch_if_enabled
                        auto_watch_if_enabled(job_id)
                    except Exception:
                        pass
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

    def _normalize_smiles_converter_value(value):
        text = str(value or '').strip()
        if not text or text == SMILES_CONVERTER_PLACEHOLDER:
            return ''
        normalized = text.upper()
        if normalized in {'QUICK', 'NORMAL', 'GUPPY', 'ARCHITECTOR'}:
            return normalized
        return ''

    def _get_smiles_converter_from_control(control_content):
        try:
            parsed = parse_control_text(control_content)
        except Exception:
            parsed = {}
        return _normalize_smiles_converter_value(parsed.get('smiles_converter', ''))

    def _validate_smiles_converter_requirement(control_content, raw_input, *, batch_has_smiles=False):
        converter = _get_smiles_converter_from_control(control_content)
        single_has_smiles = False

        raw_input = str(raw_input or '').strip()
        if raw_input:
            try:
                cleaned_data, input_type = clean_input_data(raw_input)
                single_has_smiles = input_type == 'smiles' and bool(cleaned_data.strip())
            except Exception:
                cleaned_data, input_type = clean_input_data(raw_input)
                single_has_smiles = input_type == 'smiles' and bool(cleaned_data.strip())

        if (single_has_smiles or batch_has_smiles) and not converter:
            return [
                'SMILES detected in Input or Batch list: set '
                'smiles_converter=QUICK, NORMAL, GUPPY, or ARCHITECTOR.'
            ]
        return []

    def _prepare_delfin_submit_input(raw_input, cache):
        """Return the exact payload DELFIN should receive in input.txt.

        Any XYZ visible in the textarea (preview isomer, manual paste, or
        Manipulate-edited geometry) is submitted verbatim. The cached SMILES
        from a previous conversion is still returned so CONTROL.txt can be
        annotated with it for charge auto-fill, but input.txt carries the
        XYZ coordinates the user is actually looking at.
        """
        input_content, input_type = clean_input_data(raw_input)
        cached_smiles = str((cache or {}).get('smiles') or '').strip() or None

        if input_type == 'smiles':
            submit_smiles = input_content.strip()
            return submit_smiles + '\n', 'smiles', submit_smiles

        return input_content, input_type, cached_smiles

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
                    # Distinguish key=value pairs from SMILES containing '='.
                    # Use RDKit: if the whole token parses as a molecule it is
                    # SMILES, not a key=value pair.
                    _is_smi = False
                    try:
                        from rdkit import Chem
                        _mol = Chem.MolFromSmiles(token, sanitize=False)
                        _is_smi = _mol is not None and _mol.GetNumAtoms() > 0
                    except Exception:
                        _is_smi = False
                    if not _is_smi:
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

            # --- Lookahead: does a *-terminated XYZ block follow? ---
            has_xyz_block = False
            j = i
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and lines[j].strip().upper() == 'XYZ':
                j += 1
            while j < len(lines):
                sj = lines[j].strip()
                if sj == '*':
                    has_xyz_block = True
                    break
                if ';' in sj:          # next header line → no block
                    break
                j += 1

            is_xyz = has_xyz_block or force_xyz

            # --- Resolve free tokens based on mode ---
            smiles_payload = None
            if header_free_tokens:
                if is_xyz:
                    # Pick the first source-path-looking token if we don't
                    # have one yet; silently ignore the rest.
                    if not source_path:
                        for ft in header_free_tokens:
                            if _batch_token_looks_like_source(ft):
                                source_path = ft
                                break
                else:
                    smiles_payload = header_free_tokens[0]

            # --- SMILES mode ---
            if not is_xyz:
                if not smiles_payload:
                    errors.append(f"Line {line_no}: No SMILES or XYZ coordinates for '{name}'")
                    continue
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

            # --- XYZ mode: read optional marker + block until "*" ---
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i < len(lines) and lines[i].strip().upper() == 'XYZ':
                i += 1

            xyz_lines = []
            if has_xyz_block:
                while i < len(lines):
                    if lines[i].strip() == '*':
                        i += 1
                        break
                    xyz_lines.append(lines[i].rstrip())
                    i += 1

            xyz_raw = '\n'.join(xyz_lines).strip()
            xyz_content = _clean_xyz_block(xyz_raw) if xyz_raw else ''

            if not xyz_content and not source_path:
                errors.append(f"Line {line_no}: XYZ block for '{name}' is empty and no source path given")
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

    def batch_has_smiles_entries():
        return any(entry.get('input_kind') == 'smiles' for entry in get_smiles_list_entries())

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

        if error:
            lines = [f'Preview: {name}']
            if input_kind == 'smiles':
                lines.append(f"SMILES: {entry['input_raw']}")
            elif source_path:
                lines.append(f'Source: {source_path}')
            if extras:
                lines.append(f'Options: {extras}')
            lines.append(f'Error: {error}')
            _replace_mol_output_text(*lines)
            return

        lines = []
        if input_kind == 'smiles':
            lines.append(f'Preview: {name} ({method})')
            lines.append(f"SMILES: {entry['input_raw']}")
        else:
            lines.append(f'Preview: {name} (XYZ)')
            if source_path:
                lines.append(f'Source: {source_path}')
        if extras:
            lines.append(f'Options: {extras}')
        _set_mol_status(*lines)
        mol_output.outputs = _build_mol_output_bundle(xyz_data)

    def preview_smiles_at_index(index, *, delay=0.0):
        entries = get_smiles_list_entries()
        if not entries:
            _cancel_batch_preview_timer()
            state['batch_preview_task_id'] += 1
            _set_batch_preview_busy(False)
            update_smiles_preview_label()
            _replace_mol_output_text('No valid batch entries.')
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

        status_lines = [f'Preview: {entry["name"]} [{entry["input_kind"].upper()}]']
        if entry['input_kind'] == 'smiles':
            status_lines.append(f"SMILES: {entry['input_raw']}")
        else:
            status_lines.append('XYZ coordinates')
            source_path = str(entry.get('source_path') or '').strip()
            if source_path:
                status_lines.append(f'Source: {source_path}')
        if entry['extras']:
            status_lines.append(f'Options: {entry["extras"]}')
        status_lines.append('Converting...')
        _clear_mol_output()
        _set_mol_status(*status_lines, spinner=True)

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

            _schedule_ui_update(_apply_batch_preview_result, task_id, entry, preview_payload)

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
            _replace_mol_output_text('No valid batch entries.')
            return
        new_index = state['smiles_preview_index'] - 1
        if new_index < 0:
            new_index = len(entries) - 1
        preview_smiles_at_index(new_index)

    def handle_smiles_next(button):
        entries = get_smiles_list_entries()
        if not entries:
            _replace_mol_output_text('No valid batch entries.')
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

    def _build_goat_control_content(
        *,
        template,
        charge_value,
        solvent_value,
        smiles_converter_value,
        pal_value,
        maxcore_value,
        submit_smiles=None,
        source_path='',
        extras=None,
    ):
        control_content = (
            template
            .replace('[CHARGE]', str(charge_value))
            .replace('[SOLVENT]', solvent_value)
        )
        control_content = _set_control_value(
            control_content,
            'smiles_converter',
            smiles_converter_value,
        )
        control_content = _set_control_value(
            control_content,
            'PAL',
            str(pal_value),
        )
        control_content = _set_control_value(
            control_content,
            'maxcore',
            str(maxcore_value),
        )
        if submit_smiles:
            control_content = _set_control_smiles(control_content, submit_smiles)
        for key, value in (extras or {}).items():
            control_content = _set_control_value(control_content, key, value)
        if source_path:
            control_content = _set_control_value(control_content, 'source', source_path)
        return control_content

    def handle_submit_smiles_list(button):
        with smiles_batch_output:
            clear_output()
            job_prefix = job_name_widget.value.strip()
            if not job_prefix:
                print('Error: Job name cannot be empty!')
                return

            control_content_base = control_widget.value
            control_errors = validate_control_text(control_content_base)
            control_errors.extend(
                _validate_smiles_converter_requirement(
                    control_content_base,
                    coords_widget.value,
                    batch_has_smiles=batch_has_smiles_entries(),
                )
            )
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

            any_success = False
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
                    any_success = True
                    job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                    try:  # auto-watch (no-op unless job monitoring enabled)
                        from delfin.agent.job_monitor import auto_watch_if_enabled
                        auto_watch_if_enabled(job_id)
                    except Exception:
                        pass
                    if mode == 'delfin-co2-chain':
                        print(
                            f'Submitted {safe_job_name} '
                            f'[{input_kind.upper()}] (ID: {job_id}, CO2 delta: {co2_delta})'
                        )
                    else:
                        print(f'Submitted {safe_job_name} [{input_kind.upper()}] (ID: {job_id})')
                else:
                    print(f'Failed {safe_job_name}: {result.stderr or result.stdout}')

            if any_success:
                reset_form()

    def reset_form():
        _cancel_batch_preview_timer()
        job_name_widget.value = ''
        coords_widget.value = ''
        smiles_batch_widget.value = ''
        control_widget.value = ctx.default_control
        job_type_widget.value = '48h'
        custom_time_widget.hours_widget.value = 72
        custom_time_widget.minutes_widget.value = 0
        custom_time_widget.seconds_widget.value = 0
        only_goat_charge.value = 0
        only_goat_solvent.value = ''
        only_goat_smiles_converter.value = 'QUICK'
        only_goat_pal.value = 12
        only_goat_maxcore.value = 600
        only_goat_time.value = '24:00:00'
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
        _replace_mol_output_text('Please enter XYZ coordinates or SMILES.')

    def handle_submit(button):
        raw_input = coords_widget.value.strip()
        batch_text = smiles_batch_widget.value.strip()
        if raw_input and batch_text:
            with output_area:
                clear_output(wait=True)
                print(
                    'Warning: both the single-job input and the batch list '
                    'are filled — submission cannot proceed. Please clear '
                    'one of the two before submitting.'
                )
            return
        if not raw_input and batch_text:
            handle_submit_smiles_list(button)
            return
        with output_area:
            clear_output()
            job_name = job_name_widget.value.strip()
            control_content = control_widget.value

            if not job_name:
                print('Error: Job name cannot be empty!')
                return
            if not raw_input:
                print('Error: Input (coordinates or SMILES) cannot be empty!')
                return

            control_errors = validate_control_text(control_content)
            control_errors.extend(
                _validate_smiles_converter_requirement(
                    control_content,
                    raw_input,
                )
            )
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
                    try:  # auto-watch (no-op unless job monitoring enabled)
                        from delfin.agent.job_monitor import auto_watch_if_enabled
                        auto_watch_if_enabled(job_id)
                    except Exception:
                        pass
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
            errors.extend(
                _validate_smiles_converter_requirement(
                    control_widget.value,
                    coords_widget.value,
                    batch_has_smiles=batch_has_smiles_entries(),
                )
            )
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
            smiles_converter_value = (
                _normalize_smiles_converter_value(only_goat_smiles_converter.value) or 'QUICK'
            )
            pal_value = int(only_goat_pal.value or 0)
            maxcore_value = int(only_goat_maxcore.value or 0)
            time_limit = str(only_goat_time.value or '').strip() or '48:00:00'

            if not job_name:
                print('Error: Job name cannot be empty!')
                return
            if not solvent_value:
                print('Error: Solvent cannot be empty!')
                return
            if pal_value <= 0:
                print('Error: PAL must be > 0.')
                return
            if maxcore_value <= 0:
                print('Error: MaxCore must be > 0.')
                return
            try:
                parse_time_to_seconds(time_limit)
            except Exception:
                print('Error: Time must use HH:MM:SS, e.g. 48:00:00.')
                return

            try:
                # Use template from file (BwUni) or inline (local)
                template = ctx.only_goat_template
                if ctx.only_goat_template_path and ctx.only_goat_template_path.exists():
                    template = ctx.only_goat_template_path.read_text()
                if raw_input:
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
                    job_dir.mkdir(parents=True, exist_ok=True)

                    control_content = _build_goat_control_content(
                        template=template,
                        charge_value=charge_value,
                        solvent_value=solvent_value,
                        smiles_converter_value=smiles_converter_value,
                        pal_value=pal_value,
                        maxcore_value=maxcore_value,
                        submit_smiles=submit_smiles,
                    )

                    control_errors = validate_control_text(control_content)
                    if control_errors:
                        print('CONTROL.txt validation failed:')
                        for err in control_errors:
                            print(f'- {err}')
                        return

                    (job_dir / 'CONTROL.txt').write_text(control_content)
                    (job_dir / 'input.txt').write_text(input_content)

                    result = ctx.backend.submit_delfin(
                        job_dir=job_dir, job_name=safe_job_name, mode='delfin',
                        time_limit=time_limit, pal=pal_value, maxcore=maxcore_value,
                    )

                    if result.returncode == 0:
                        job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                        try:  # auto-watch (no-op unless job monitoring enabled)
                            from delfin.agent.job_monitor import auto_watch_if_enabled
                            auto_watch_if_enabled(job_id)
                        except Exception:
                            pass
                        print('GOAT job successfully submitted!')
                        print(f'Job ID: {job_id}')
                        print(f'Time Limit: {time_limit}')
                        print(f'Input Type: {input_type.upper()}')
                        print(f'Charge: {charge_value}')
                        print(f'Solvent: {solvent_value}')
                        print(f'PAL: {pal_value}')
                        print(f'MaxCore: {maxcore_value}')
                        if input_type == 'smiles':
                            print(f'SMILES Converter: {smiles_converter_value}')
                        print(f'Directory: {job_dir}')
                        print('')
                        print('Check status in Job Status tab')
                        reset_form()
                    else:
                        print('Error submitting job:')
                        print(result.stderr)
                    return

                entries, parse_errors = parse_batch_entries()
                for err in parse_errors:
                    print(err)
                if not entries:
                    print('Error: Input is empty and no valid batch entries were found.')
                    return

                submitted = 0
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

                    full_job_name = f'{job_name}_{safe_name}'
                    safe_job_name = ''.join(c for c in full_job_name if c.isalnum() or c in ('_', '-'))
                    if not safe_job_name:
                        print(f'Line {line_no}: Invalid job name -> {full_job_name}')
                        continue

                    job_dir = ctx.calc_dir / safe_job_name
                    job_dir.mkdir(parents=True, exist_ok=True)

                    entry_charge = charge_value
                    submit_smiles = None
                    if input_kind == 'smiles':
                        smi = entry.get('input_raw', '').strip()
                        if not smi:
                            print(f'Line {line_no}: Missing SMILES payload for {safe_name}')
                            continue
                        input_content = smi + '\n'
                        submit_smiles = smi
                        if _needs_smiles_charge(template, extras):
                            entry_charge = _get_smiles_charge(smi)
                    else:
                        input_content = (entry.get('input_content') or '').strip()
                        if not input_content:
                            print(f'Line {line_no}: {safe_name} - empty XYZ payload')
                            continue
                        input_content = _clean_xyz_block(input_content)
                        if not input_content:
                            print(f'Line {line_no}: {safe_name} - empty XYZ payload')
                            continue
                        if 'charge' in extras:
                            try:
                                entry_charge = int(str(extras['charge']).strip())
                            except Exception:
                                entry_charge = extras['charge']

                    control_content = _build_goat_control_content(
                        template=template,
                        charge_value=entry_charge,
                        solvent_value=solvent_value,
                        smiles_converter_value=smiles_converter_value,
                        pal_value=pal_value,
                        maxcore_value=maxcore_value,
                        submit_smiles=submit_smiles,
                        source_path=source_path,
                        extras=extras,
                    )

                    control_errors = validate_control_text(control_content)
                    if control_errors:
                        print(f'Line {line_no}: CONTROL.txt validation failed for {safe_job_name}:')
                        for err in control_errors:
                            print(f'- {err}')
                        continue

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

                    result = ctx.backend.submit_delfin(
                        job_dir=job_dir, job_name=safe_job_name, mode='delfin',
                        time_limit=time_limit, pal=pal_value, maxcore=maxcore_value,
                    )
                    if result.returncode == 0:
                        job_id = result.stdout.strip().split()[-1] if result.stdout.strip() else '(unknown)'
                        try:  # auto-watch (no-op unless job monitoring enabled)
                            from delfin.agent.job_monitor import auto_watch_if_enabled
                            auto_watch_if_enabled(job_id)
                        except Exception:
                            pass
                        print(f'Submitted {safe_job_name} [GOAT {input_kind.upper()}] (ID: {job_id})')
                        submitted += 1
                    else:
                        print(f'Failed {safe_job_name}: {result.stderr or result.stdout}')

                if submitted:
                    print('')
                    print('Check status in Job Status tab')
                    reset_form()
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
                    try:  # auto-watch (no-op unless job monitoring enabled)
                        from delfin.agent.job_monitor import auto_watch_if_enabled
                        auto_watch_if_enabled(job_id)
                    except Exception:
                        pass
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
    coords_widget.observe(update_molecule_view, names='value')
    build_complex_button.on_click(handle_build_complex)
    architector_button.on_click(handle_architector_convert)
    guppy_submit_button.on_click(handle_guppy_submit)
    fukui_submit_button.on_click(handle_fukui_submit)
    smiles_prev_button.on_click(handle_smiles_prev)
    smiles_next_button.on_click(handle_smiles_next)

    def on_batch_change(change):
        state['smiles_preview_index'] = 0
        state['batch_preview_cache'] = {}
        preview_smiles_at_index(0, delay=0.35)

    smiles_batch_widget.observe(on_batch_change, names='value')
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
        # Drawing comes before converting, and the editor stands between them:
        # the order on screen is the order of the work -- draw it, hand it over
        # as a SMILES, then turn that into coordinates.  With the editor open,
        # the buttons that act on what it produced are below it, where the eye
        # arrives after the drawing rather than before it.
        widgets.HBox([submit_draw_open_btn, submit_draw_get_btn,
                      submit_draw_update_btn],
                     layout=widgets.Layout(gap='10px', flex_wrap='wrap')),
        submit_draw_frame, submit_draw_sync,
        widgets.HBox([convert_smiles_button, convert_smiles_uff_button,
                      convert_smiles_quick_button],
                     layout=widgets.Layout(gap='10px', flex_wrap='wrap')),
        widgets.HBox([build_complex_button, architector_button],
                     layout=widgets.Layout(gap='10px', flex_wrap='wrap')),
        manta_settings_row,
        manta_button,          # MANTA button sits directly UNDER its settings
        spacer_large,
        widgets.HTML('<b>Batch SMILES/XYZ:</b>'),
        smiles_batch_widget, spacer,
        widgets.HBox(
            [smiles_prev_button, smiles_preview_label, smiles_next_button],
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

    xyz_copy_row = widgets.HBox(
        [xyz_copy_btn, xyz_copy_status],
        layout=widgets.Layout(gap='6px', align_items='center', flex_wrap='wrap'),
    )
    # What travels into the overlay, and what it is once it is there. These are
    # the same names the ORCA Builder, the Calculations browser and the Archive
    # use: there is one fullscreen, and a fix to it is a fix in all of them.
    # The order in the overlay is the order here, which is the order of
    # submit_right's children -- toolbar, status, viewer, isomer nav, copy row.
    # The picture, with the status line lying along its bottom edge rather than
    # standing in a row above it. Above it, every message of a different length
    # moved the structure the user was aiming at; here it costs no layout and
    # cannot move anything. Both copies live in the stack -- the ordinary one
    # shows in this view, the other is taken by fullscreen and shows there.
    mol_viewer_stack = widgets.Box(
        # The profile of the last walk shares this box with the
        # structure and swaps with it: one or the other, in the same
        # place, so nothing above or below moves. In a row of its own
        # under the panel it was a second panel that pushed the page
        # down and landed off the bottom of the screen.
        [mol_output, _editor.submit_scan_plot, _editor.submit_view_panel,
         _editor.submit_manip_status, mol_status],
        layout=widgets.Layout(width='100%', min_width='0'),
    )
    mol_viewer_stack.add_class('delfin-structure-viewer-stack')
    # The stack is what fullscreen takes, so the line goes with the picture it
    # is about and comes back with it. That is what the second copy existed
    # for -- a line moved by hand did not come back -- and nothing is moved by
    # hand any more.
    mol_viewer_stack.add_class('delfin-structure-fs-viewer')

    # The toolbar and the status copy are marked by the editor itself, being
    # the editor's wherever it is built; these are this tab's own.
    for _member in (mol_viewer_stack, isomer_nav_row, xyz_copy_row):
        _member.add_class('delfin-structure-fs-member')
    # A strip that keeps its own height and takes the full width, which is what
    # the shared sheet calls a toolbar; the Builder's block stepper is one too.
    isomer_nav_row.add_class('delfin-structure-fs-toolbar')
    xyz_copy_row.add_class('delfin-structure-fs-toolbar')
    mol_output.add_class('delfin-structure-fs-viewer')

    submit_right = widgets.VBox([
        widgets.HTML('<b>Molecule Preview:</b>'),
        submit_manip_toolbar, mol_viewer_stack, isomer_nav_row, xyz_copy_row,
        submit_ff_notes,
        spacer_large,
        widgets.HTML('<b>GOAT:</b>'),
        widgets.VBox([
            widgets.HBox(
                [
                    only_goat_charge,
                    only_goat_solvent,
                    only_goat_smiles_converter,
                ],
                layout=widgets.Layout(gap='8px', flex_wrap='wrap', align_items='center'),
            ),
            widgets.HBox(
                [
                    only_goat_pal,
                    only_goat_maxcore,
                    only_goat_time,
                    only_goat_submit_button,
                ],
                layout=widgets.Layout(gap='8px', flex_wrap='wrap', align_items='center'),
            ),
        ], layout=widgets.Layout(gap='8px')),
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
        spacer_large,
        widgets.HTML('<b>Fukui Indices (atomic):</b>'),
        widgets.VBox([
            widgets.HBox(
                [fukui_skip_cubes],
                layout=widgets.Layout(gap='8px', flex_wrap='wrap', align_items='center'),
            ),
            widgets.HBox(
                [fukui_pal, fukui_maxcore, fukui_time_limit, fukui_submit_button],
                layout=widgets.Layout(gap='8px', flex_wrap='wrap', align_items='center'),
            ),
        ], layout=widgets.Layout(gap='8px')),
        fukui_output,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='10px',
        box_sizing='border-box', overflow_x='hidden',
    ))
    submit_right.add_class(submit_scope_id)
    # The box the shared fullscreen gathers its members from, and the word it
    # looks this tab up by. It is the whole right-hand column rather than a new
    # wrapper around the preview: a wrapper would be another flex level between
    # submit_right and its children, and every gap and width below it would be
    # laid out differently for no gain -- only the marked members travel, and
    # they are direct children either way.
    submit_right.add_class('delfin-structure-fs-module')
    submit_right.add_class('submit-structure-fs-module')

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
        .submit-mol-output,
        .submit-mol-output .output_subarea,
        .submit-mol-output .jp-OutputArea,
        .submit-mol-output .jp-OutputArea-output,
        .submit-mol-output .jp-OutputArea-child,
        .submit-mol-output .output_area {
            overflow: hidden !important;
        }
        .submit-mol-output [id^="3dmolviewer_"] {
            overflow: hidden !important;
            max-width: 100% !important;
            max-height: 100% !important;
        }
        /* The gutter a notebook keeps for "Out[7]:". There is no Out[7] here
           and never will be, but the column is reserved all the same -- and in
           fullscreen it is a white band down the left of the picture, inside
           the blue frame. It is taken away in the ordinary view as well: the
           viewer should start where its frame starts. The overlay's half of
           this lives in the shared sheet below. */
        .submit-mol-output .jp-OutputPrompt,
        .submit-mol-output .jp-OutputArea-prompt,
        .submit-mol-output .prompt {
            display: none !important;
            width: 0 !important;
            min-width: 0 !important;
            flex: 0 0 0 !important;
            padding: 0 !important;
            margin: 0 !important;
        }
        .submit-mol-output .jp-OutputArea-child {
            padding-left: 0 !important;
            margin-left: 0 !important;
        }
        """
        # This tab's overlay used to be described here, in rules of its own,
        # beside a second implementation that did the moving. Both are the
        # shared ones now -- the same sheet and the same script the ORCA
        # Builder, the Calculations browser and the Archive are laid out by.
        + structure_viewer_fullscreen_css()
        + """
        </style>
        """
    )
    tab_widget.add_class('submit-split-root')
    submit_left.add_class('submit-split-pane')
    submit_right.add_class('submit-split-pane')
    mol_output.add_class('submit-split-pane')
    mol_output.add_class('submit-mol-output')
    _replace_mol_output_text('Please enter XYZ coordinates or SMILES.')
    # The one fullscreen, and this tab saying which viewer is its own. Sent the
    # way the other three tabs send it, so whichever the user opens first
    # installs it and the rest find it already there.
    ctx.add_init_js(
        structure_viewer_fullscreen_bootstrap_js()
        + '\n'
        + structure_viewer_fullscreen_kind_js(
            'submit', 'submit-scope-', ['_submitMolViewerByScope']))
    tab_widget = widgets.VBox([submit_css, tab_widget], layout=widgets.Layout(width='100%'))

    return tab_widget, {
        'reset_form': reset_form,
        'mol_output': mol_output,
        'control_widget': control_widget,
        'coords_widget': coords_widget,
        'submit_button': submit_button,
        'job_name_widget': job_name_widget,
        'smiles_batch_widget': smiles_batch_widget,
        'job_type_widget': job_type_widget,
        'custom_time_widget': custom_time_widget,
        'handle_submit': handle_submit,
        'handle_validate_control': handle_validate_control,
        # Everything the editor is made of, under the names it uses.
        **_editor.exported,
        'submit_draw_open_btn': submit_draw_open_btn,
        'submit_draw_get_btn': submit_draw_get_btn,
        'submit_draw_update_btn': submit_draw_update_btn,
        'submit_draw_frame': submit_draw_frame,
        'submit_draw_sync': submit_draw_sync,
        # The two channels the page speaks through, and the line it is
        # answered on: a test that cannot use them has to drive the editor
        # through its internals instead of the way the browser drives it.
        'editor_state': state,
        'refresh_constraints': _refresh_constraints,
    }
