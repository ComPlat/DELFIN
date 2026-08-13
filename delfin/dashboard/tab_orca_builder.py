"""ORCA Input Builder tab: build, preview, and submit standalone ORCA jobs."""

import re
import shutil
from pathlib import Path
from collections import Counter
from itertools import permutations, product

import py3Dmol
import ipywidgets as widgets
import html
import json
import uuid
import numpy as np

from . import structure_editor as _structure_editor

from IPython import get_ipython
from IPython.display import clear_output, display, HTML

from delfin.common.control_validator import (
    ORCA_FUNCTIONALS, ORCA_BASIS_SETS, DISP_CORR_VALUES, _RI_JKX_KEYWORDS,
)
from delfin.user_settings import (
    load_orca_templates, save_orca_template, delete_orca_template,
)

from .molecule_viewer import (
    strip_xyz_header,
    DEFAULT_3DMOL_ZOOM,
    get_viewer_profile, viewer_disabled_html,
    molecule_view_style_js,
    patch_viewer_mouse_controls_js,
    structure_viewer_fullscreen_bootstrap_js,
    structure_viewer_fullscreen_css,
    structure_viewer_fullscreen_kind_js,
)
from .input_processing import (
    parse_inp_resources, sanitize_orca_input, clean_input_data,
    smiles_to_xyz_quick_with_previews,
)


_COORD_LINE_RE = re.compile(r'^\s*\*\s*(?:xyzfile|xyz|gzmt|internal)\b', re.IGNORECASE)

# On-screen size of the atom numbers in the molecule preview, as the factor the
# high-resolution label texture is down-scaled by.  Selectable in the preview
# toolbar; this is the size a fresh viewer starts with.


def strip_coord_block(text):
    """Return the INP text with a trailing ORCA coordinate block removed.

    Cuts everything from the first coordinate specifier line (``* xyz``,
    ``* xyzfile``, ...) onward, so the remaining ``inp_body`` holds only the
    settings and can be templated independently of any coordinates.
    """
    if not text:
        return ''
    lines = text.split('\n')
    for i, line in enumerate(lines):
        if _COORD_LINE_RE.match(line):
            return '\n'.join(lines[:i]).rstrip()
    return text.rstrip()


def _orca_kabsch_align(reference_coords, target_coords, mapping=None):
    """Superpose the reference onto the target and return (aligned, RMSD).

    The rotation that minimises the RMSD is ``R = V U^T`` from the SVD of the
    covariance ``ref^T target``, and it acts on *column* vectors.  Coordinates
    here are rows, so it has to be applied transposed -- ``rows @ R.T``.
    Applying it untransposed rotates the reference by the inverse and reports
    the RMSD of that, which is why an exactly rotated copy of a structure used
    to come out several Angstrom apart from itself.
    """
    ref = np.asarray(reference_coords, dtype=float)
    target = np.asarray(target_coords, dtype=float)
    if ref.shape != target.shape:
        raise ValueError(f'Coordinate shape mismatch: {ref.shape} vs {target.shape}.')
    if mapping is not None:
        map_idx = np.asarray(mapping, dtype=int)
        if map_idx.ndim != 1 or map_idx.size != ref.shape[0]:
            raise ValueError('Invalid mapping length for Kabsch alignment.')
        target = target[map_idx]
    ref_centroid = ref.mean(axis=0)
    target_centroid = target.mean(axis=0)
    ref_centered = ref - ref_centroid
    target_centered = target - target_centroid
    covariance = ref_centered.T @ target_centered
    u, _s, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1.0
        rotation = vt.T @ u.T
    aligned = ref_centered @ rotation.T + target_centroid
    diff = aligned - target
    rmsd = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))
    return aligned, rmsd


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
            'Paste XYZ coordinates, SMILES, or use named file blocks:\n\n'
            'SMILES example:  [Fe+2]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]\n'
            'Use CONVERT SMILES button to generate 3D coordinates.\n\n'
            'name1.xyz;\n6\nComment\nC  0.0  0.0  0.0\n...\n*\n\n'
            'name2.xyz;\nFe  0.0  0.0  0.0\nC   1.5  0.0  0.0\n*\n\n'
            'Named blocks write .xyz files to the job directory and use\n'
            '* xyzfile in the INP (header optional - auto-added if missing).\n'
            'Navigate all molecules with the ◀ ▶ buttons in the preview.\n\n'
            'Or plain XYZ (with or without header):\n'
            '6\nComment\nC  0.0  0.0  0.0\n...\nor just:\nC  0.0  0.0  0.0\n...'
        ),
        description='Coordinates:',
        layout=widgets.Layout(width='100%', height='520px', box_sizing='border-box'), style=ws,
        # Every keystroke used to rebuild the 3D preview from scratch, so typing
        # coordinates fought a viewer that was being torn down and recreated
        # between characters.  The browser now reports the box once typing
        # stops (see debounce_input), and the preview follows in one step.
        continuous_update=False,
    )
    orca_coords.add_class('delfin-debounced')
    orca_copy_coords_btn = widgets.Button(
        description='COPY COORDINATES',
        button_style='',
        layout=widgets.Layout(width='200px'),
    )
    orca_check_numbering_btn = widgets.Button(
        description='CHECK NUMBERING',
        button_style='warning',
        layout=widgets.Layout(width='200px'),
    )
    orca_apply_numbering_btn = widgets.Button(
        description='APPLY NUMBERING FIX',
        button_style='success',
        disabled=True,
        layout=widgets.Layout(width='220px'),
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
        value='acetonitrile', description='Solvent:',
        layout=widgets.Layout(width='250px'), style=ws,
    )

    orca_print_mos = widgets.Checkbox(value=False, description='Print MOs',
                                      layout=widgets.Layout(width='250px'), style=ws)
    orca_print_basis = widgets.Checkbox(value=False, description='Print Basis',
                                        layout=widgets.Layout(width='250px'), style=ws)

    orca_pal = widgets.IntText(value=12, description='PAL (cores):',
                               layout=widgets.Layout(width='250px'), style=ws)
    orca_maxcore = widgets.IntText(value=6000, description='MaxCore (MB):',
                                   layout=widgets.Layout(width='250px'), style=ws)
    orca_slurm_time = widgets.Combobox(
        value='06:00:00', placeholder='e.g. 02:00:00',
        options=['00:30:00', '01:00:00', '02:00:00', '06:00:00',
                 '12:00:00', '24:00:00', '48:00:00', '72:00:00'],
        ensure_option=False,  # allow arbitrary walltimes, not only the presets
        description='Time Limit:',
        layout=widgets.Layout(width='250px'), style=ws)

    orca_file_upload = widgets.FileUpload(
        accept='', multiple=True, description='',
        layout=widgets.Layout(width='1px', height='1px', overflow='hidden'),
    )
    orca_file_upload.add_class('orca-hidden-upload')
    orca_uploaded_files_label = widgets.HTML(
        value='<i>Drag & drop files here or click to browse (e.g. .gbw, .xyz, .hess)</i>',
        layout=widgets.Layout(width='100%'),
    )
    orca_drop_zone = widgets.HTML(
        value=(
            '<div class="orca-drop-zone" style="'
            'border:2px dashed #aaa; border-radius:8px; padding:18px 12px;'
            'text-align:center; cursor:pointer; color:#666;'
            'min-height:80px; display:flex; align-items:center; justify-content:center;'
            'transition: border-color 0.2s, background 0.2s;'
            '">'
            '<span style="font-size:14px;">📁 Drop files here or click to upload<br>'
            '<small style="color:#999;">.gbw, .xyz, .hess, etc.</small></span>'
            '</div>'
        ),
        layout=widgets.Layout(width='100%'),
    )
    orca_drop_zone.add_class('orca-drop-zone-wrap')
    orca_path_files = widgets.Textarea(
        value='',
        placeholder='Paste file paths (one per line):\n/path/to/file.gbw\n/path/to/file.xyz',
        description='File Paths:',
        layout=widgets.Layout(width='100%', height='120px', box_sizing='border-box'), style=ws,
    )

    orca_preview = widgets.Textarea(
        value='', description='INP Preview:',
        layout=widgets.Layout(width='100%', height='620px', box_sizing='border-box'), style=ws,
        disabled=False,
    )

    orca_save_btn = widgets.Button(description='SAVE', button_style='warning',
                                   layout=widgets.Layout(width='150px'))
    orca_submit_btn = widgets.Button(description='SUBMIT ORCA JOB', button_style='success',
                                     layout=widgets.Layout(width='150px'))
    orca_output = widgets.Output()

    orca_mol_output = widgets.Output(layout=widgets.Layout(
        border='2px solid #1976d2', width='100%', min_height='560px', height='560px',
        overflow='hidden', box_sizing='border-box', padding='0', margin='0',
    ))
    orca_mol_output.add_class('orca-mol-output')
    orca_mol_output.add_class('delfin-structure-fs-member')
    orca_mol_output.add_class('delfin-structure-fs-viewer')
    orca_scope_id = f'orca-scope-{abs(id(orca_mol_output))}'

    orca_mol_prev_btn = widgets.Button(
        description='◀', tooltip='Previous molecule',
        layout=widgets.Layout(width='36px', height='28px'),
    )
    orca_mol_next_btn = widgets.Button(
        description='▶', tooltip='Next molecule',
        layout=widgets.Layout(width='36px', height='28px'),
    )
    orca_mol_nav_label = widgets.HTML(value='')
    orca_mol_nav_row = widgets.HBox(
        [
            orca_mol_prev_btn,
            orca_mol_nav_label,
            orca_mol_next_btn,
        ],
        layout=widgets.Layout(display='none', align_items='center', gap='6px'),
    )
    orca_mol_nav_row.add_class('delfin-structure-fs-member')
    orca_mol_nav_row.add_class('delfin-structure-fs-toolbar')
    # -- state ----------------------------------------------------------
    state = {
        'extra_files': {},
        'last_auto_keywords': '',
        'is_resetting': False,
        'xyz_blocks': [],
        'xyz_view_idx': 0,
        'numbering_check_active': False,
        'numbering_check_results': {},
        'numbering_check_block_idx': 1,
        'numbering_view_step': 0,
        # the viewer currently on screen: its div, and whether it is still the
        # thing the Output widget holds (see _update_molecule_js)
        'viewer_div': None,
        'viewer_live': False,
    }

    # -- the structure editor -------------------------------------------
    # The same one the Submit tab holds, over the block this tab is showing.
    #
    # The editor writes a structure as plain XYZ into a coordinates box and
    # lets its tab take it from there. This tab's box holds named blocks --
    # "name.xyz;comment", the atoms, "*" -- which the editor knows nothing
    # about, so it is given a box of its own and what it writes is folded back
    # into the block on screen. That is the same rebuild Apply Numbering Fix
    # does, so Check Numbering and Apply Numbering Fix are untouched by this.
    orca_editor_coords = widgets.Textarea(
        value='', layout=widgets.Layout(display='none'))
    _main_io_loop = getattr(getattr(get_ipython(), 'kernel', None), 'io_loop', None)

    def _orca_schedule_ui_update(func, *args, **kwargs):
        if _main_io_loop is not None:
            _main_io_loop.add_callback(lambda: func(*args, **kwargs))
            return
        func(*args, **kwargs)

    def _write_orca_coords(text):
        """What was drawn, into the box this tab shows."""
        orca_coords.value = text

    orca_editor = _structure_editor.build(
        ctx,
        state=state,
        coords_widget=orca_editor_coords,
        viewer_height=560,
        schedule_ui_update=_orca_schedule_ui_update,
        # Defined further down, so they are looked up when they are called.
        update_view=lambda *a, **k: _orca_editor_wrote(*a, **k),
        # This tab takes its charge from its own box, not from a SMILES.
        get_smiles_charge=lambda *a, **k: None,
        # Several structures belong in several named blocks here.
        offer_structures=lambda *a, **k: _take_structures(*a, **k),
        # A SMILES is typed into the tab's own box, not the editor's.
        read_input=lambda: orca_coords.value,
        # A loader or a line of text goes where this tab's structures go, not
        # into the editor's own output widget -- which this tab never places.
        show_output=lambda items: _show_in_viewer(*items),
        # Every named block is a frame, so "all" reaches all of them.
        list_structures=lambda: [
            (xyz, len([r for r in strip_xyz_header(xyz).split('\n') if r.strip()]),
             name[:-4] if name.lower().endswith('.xyz') else name)
            for name, xyz in (state.get('xyz_blocks') or [])],
        write_input=_write_orca_coords,
    )
    orca_editor_scope = orca_editor.submit_scope_id
    # The editor's own fullscreen button, at the head of its toolbar where the
    # Submit tab has it, rather than a second one in a row of its own. There is
    # nothing to make over any more: the editor carries the shared button, and
    # its toolbar is a member of whatever overlay it finds itself in. This tab
    # only says which of the shared overlays is being opened.
    orca_mol_fullscreen_btn = orca_editor.submit_fullscreen_btn
    orca_mol_fullscreen_btn.add_class('orca-structure-fullscreen-btn')
    # The status line that travels is the copy, not the one the small view
    # keeps. Lending the real one out is a move ipywidgets knows nothing about,
    # and the small view did not get it back -- the Submit tab was burned by
    # that and grew a second widget for it; this tab holds the same second
    # widget rather than finding out again.
    orca_editor.mol_status_fs.layout.margin = '0 0 6px 0'
    # The force-field notes stay in the small view, as they do in the Submit
    # tab: they are several lines of prose about what had to be approximated,
    # and in fullscreen they take that space off the structure they describe.

    # -- helpers --------------------------------------------------------
    def _orca_parse_xyz_block_records(text):
        """Parse named XYZ blocks like ``name.xyz;`` or shorthand ``1;comment``."""
        records = []
        pattern = re.compile(
            r'^([^;\n]+?)\s*;[ \t]*([^\n]*)\n(.*?)^\s*\*\s*$',
            re.MULTILINE | re.DOTALL,
        )
        for m in pattern.finditer(text):
            raw_name = str(m.group(1) or '').strip()
            suffix_comment = str(m.group(2) or '').strip()
            content = m.group(3).strip()
            if not raw_name:
                continue
            safe_name = ''.join(ch for ch in raw_name if ch.isalnum() or ch in ('_', '-', '.')).strip('.')
            if not safe_name:
                continue
            if not safe_name.lower().endswith('.xyz'):
                safe_name += '.xyz'
            lines = [l for l in content.split('\n') if l.strip()]
            if not lines:
                continue
            try:
                int(lines[0].strip())
                full_xyz = content
            except ValueError:
                n = len(lines)
                comment = suffix_comment or safe_name
                full_xyz = f'{n}\n{comment}\n' + '\n'.join(lines)
            records.append(
                {
                    'filename': safe_name,
                    'raw_name': raw_name,
                    'suffix_comment': suffix_comment,
                    'full_xyz': full_xyz,
                }
            )
        return records

    def parse_xyz_blocks(text):
        """Parse named XYZ blocks like ``name.xyz;`` or shorthand ``1;comment``.

        Returns a list of ``(filename, full_xyz_str)`` tuples, or *None* if
        no named blocks are found.  The returned XYZ strings always include
        the standard two-line header (atom count + empty comment line).
        """
        blocks = []
        for record in _orca_parse_xyz_block_records(text):
            blocks.append((record['filename'], record['full_xyz']))
        return blocks if blocks else None

    def _orca_parse_xyz_symbols_coords(xyz_text):
        lines = [line.rstrip() for line in str(xyz_text or '').splitlines()]
        if len(lines) < 3:
            raise ValueError('XYZ block is too short.')
        try:
            n_atoms = int(lines[0].strip())
        except Exception as exc:
            raise ValueError('First XYZ line must contain the atom count.') from exc
        coord_lines = lines[2: 2 + n_atoms]
        if len(coord_lines) != n_atoms:
            raise ValueError('XYZ block does not contain the declared number of atoms.')
        symbols = []
        coords = []
        for line in coord_lines:
            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f'Invalid XYZ coordinate line: {line}')
            symbols.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        return symbols, np.asarray(coords, dtype=float)

    def _orca_build_xyz_from_symbols_coords(symbols, coords, comment=''):
        arr = np.asarray(coords, dtype=float)
        body = '\n'.join(
            f'{sym:<2} {xyz[0]: .8f} {xyz[1]: .8f} {xyz[2]: .8f}'
            for sym, xyz in zip(symbols, arr)
        )
        return f'{len(symbols)}\n{comment}\n{body}'

    def _orca_sq_distance_matrix(a, b):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        diff = a[:, None, :] - b[None, :, :]
        return np.einsum('ijk,ijk->ij', diff, diff)

    def _orca_element_assignment_for_rotation(ref_symbols_lc, target_symbols_lc, ref_rot_centered, target_centered):
        try:
            from scipy.optimize import linear_sum_assignment
        except Exception:
            return None

        mapping = np.full(len(ref_symbols_lc), -1, dtype=int)
        ref_groups = {}
        target_groups = {}
        for idx, symbol in enumerate(ref_symbols_lc):
            ref_groups.setdefault(symbol, []).append(idx)
        for idx, symbol in enumerate(target_symbols_lc):
            target_groups.setdefault(symbol, []).append(idx)
        if set(ref_groups) != set(target_groups):
            return None

        for symbol, ref_idx in ref_groups.items():
            target_idx = target_groups.get(symbol, [])
            if len(ref_idx) != len(target_idx):
                return None
            if len(ref_idx) == 1:
                mapping[ref_idx[0]] = target_idx[0]
                continue
            costs = _orca_sq_distance_matrix(
                ref_rot_centered[np.asarray(ref_idx, dtype=int)],
                target_centered[np.asarray(target_idx, dtype=int)],
            )
            row_ind, col_ind = linear_sum_assignment(costs)
            for row_pos, col_pos in zip(row_ind, col_ind):
                mapping[int(ref_idx[row_pos])] = int(target_idx[col_pos])
        if np.any(mapping < 0):
            return None
        return mapping

    def _orca_generate_proper_axis_rotations():
        mats = []
        for perm in permutations((0, 1, 2)):
            for signs in product((-1.0, 1.0), repeat=3):
                mat = np.zeros((3, 3), dtype=float)
                for new_axis, old_axis in enumerate(perm):
                    mat[old_axis, new_axis] = signs[new_axis]
                if np.linalg.det(mat) > 0.0:
                    mats.append(mat)
        return mats

    def _orca_topology_mapping_from_xyz(ref_symbols, ref_coords, target_symbols, target_coords):
        try:
            from rdkit import Chem
            from rdkit.Chem import rdDetermineBonds, rdMolAlign
        except Exception:
            return None

        ref_xyz = _orca_build_xyz_from_symbols_coords(ref_symbols, ref_coords, comment='Reference')
        target_xyz = _orca_build_xyz_from_symbols_coords(target_symbols, target_coords, comment='Target')
        ref_mol = Chem.MolFromXYZBlock(ref_xyz)
        target_mol = Chem.MolFromXYZBlock(target_xyz)
        if ref_mol is None or target_mol is None:
            return None
        try:
            rdDetermineBonds.DetermineConnectivity(ref_mol)
            rdDetermineBonds.DetermineConnectivity(target_mol)
            _rmsd, _transform, atom_map = rdMolAlign.GetBestAlignmentTransform(
                ref_mol, target_mol, maxMatches=20000, reflect=False
            )
        except Exception:
            return None
        if not atom_map:
            return None
        mapping = np.full(len(ref_symbols), -1, dtype=int)
        for probe_idx, ref_idx in atom_map:
            mapping[int(probe_idx)] = int(ref_idx)
        if np.any(mapping < 0) or np.unique(mapping).size != len(ref_symbols):
            return None
        ref_seq = [s.lower() for s in ref_symbols]
        target_seq = [s.lower() for s in target_symbols]
        if not all(ref_seq[i] == target_seq[mapping[i]] for i in range(len(ref_seq))):
            return None
        return mapping

    def _orca_check_numbering_pair(ref_symbols, ref_coords, target_symbols, target_coords):
        if len(ref_symbols) != len(target_symbols):
            raise ValueError(
                f'Atom count mismatch: ref {len(ref_symbols)} vs target {len(target_symbols)}.'
            )
        ref_seq = [s.lower() for s in ref_symbols]
        target_seq = [s.lower() for s in target_symbols]
        if Counter(ref_seq) != Counter(target_seq):
            raise ValueError('Element composition mismatch.')

        n_atoms = len(ref_symbols)
        identity = np.arange(n_atoms, dtype=int)
        direct_rmsd = None
        direct_aligned = None
        if ref_seq == target_seq:
            direct_aligned, direct_rmsd = _orca_kabsch_align(ref_coords, target_coords)

        # Primary strategy: DELFIN's universal atom mapper (metal- and
        # reaction-aware, returns a *verified* graph isomorphism).  The legacy
        # geometry strategies below are kept only as an automatic fallback for
        # the cases the mapper cannot resolve (too many simultaneous bond
        # changes, or a broken / ambiguous connectivity graph).
        best_mapping = None
        best_source = None
        best_aligned = None
        best_rmsd = float('inf')
        n_bond_edits = 0
        bond_edits = []
        used_fallback = False

        xyzmap_result = None
        try:
            from delfin.atom_mapping import map_atoms as _delfin_map_atoms
            xyzmap_result = _delfin_map_atoms(
                ref_symbols, ref_coords, target_symbols, target_coords
            )
        except ValueError:
            raise  # different molecular formula -> reported as not comparable
        except Exception:
            xyzmap_result = None  # any other failure -> legacy fallback below

        if xyzmap_result is not None and xyzmap_result.get('verified'):
            best_mapping = np.asarray(xyzmap_result['order'], dtype=int)
            best_source = xyzmap_result.get('method', 'xyzmap')
            best_aligned, best_rmsd = _orca_kabsch_align(
                ref_coords, target_coords, mapping=best_mapping
            )
            best_rmsd = float(best_rmsd)
            n_bond_edits = int(xyzmap_result.get('n_bond_edits', 0))
            bond_edits = (
                list(xyzmap_result.get('bond_edits_ref', []))
                + list(xyzmap_result.get('bond_edits_tgt', []))
            )
        else:
            used_fallback = True
            best_mapping = np.arange(n_atoms, dtype=int)
            best_source = 'direct'
            best_aligned = direct_aligned
            best_rmsd = float(direct_rmsd) if direct_rmsd is not None else float('inf')

            topo_mapping = _orca_topology_mapping_from_xyz(ref_symbols, ref_coords, target_symbols, target_coords)
            if topo_mapping is not None:
                aligned, rmsd = _orca_kabsch_align(ref_coords, target_coords, mapping=topo_mapping)
                if rmsd < best_rmsd - 1e-12:
                    best_mapping, best_source, best_aligned, best_rmsd = topo_mapping, 'rdkit-topology', aligned, float(rmsd)

            ref_centered = np.asarray(ref_coords, dtype=float) - np.asarray(ref_coords, dtype=float).mean(axis=0)
            target_centered = np.asarray(target_coords, dtype=float) - np.asarray(target_coords, dtype=float).mean(axis=0)
            for rot_guess in _orca_generate_proper_axis_rotations():
                mapping = _orca_element_assignment_for_rotation(
                    ref_seq,
                    target_seq,
                    ref_centered @ rot_guess,
                    target_centered,
                )
                if mapping is None:
                    continue
                aligned, rmsd = _orca_kabsch_align(ref_coords, target_coords, mapping=mapping)
                if rmsd < best_rmsd - 1e-12:
                    best_mapping, best_source, best_aligned, best_rmsd = mapping, 'global-permutation', aligned, float(rmsd)

        if n_bond_edits > 0:
            # Genuine reaction mapping (educt <-> product): not "wrong
            # numbering", but a reorder onto the reference numbering IS available.
            numbering_ok = False
            suspicious = False
        elif not used_fallback and not np.array_equal(best_mapping, identity):
            # The verified mapper found a real same-molecule reordering.
            numbering_ok = False
            suspicious = True
        else:
            numbering_ok = bool(
                np.array_equal(best_mapping, identity)
                or (
                    direct_rmsd is not None
                    and best_rmsd >= float(direct_rmsd) - 1e-4
                )
            )
            suspicious = bool(
                not numbering_ok
                and direct_rmsd is not None
                and best_rmsd + 0.10 < float(direct_rmsd)
                and best_rmsd <= 0.60
            )
            if direct_rmsd is None and not np.array_equal(best_mapping, identity) and best_rmsd <= 0.60:
                suspicious = True

        return {
            'direct_rmsd': None if direct_rmsd is None else float(direct_rmsd),
            'best_rmsd': float(best_rmsd),
            'best_mapping': [int(v) for v in np.asarray(best_mapping, dtype=int).tolist()],
            'best_source': best_source,
            'n_bond_edits': int(n_bond_edits),
            'bond_edits': [
                (int(i) + 1, str(ei), int(j) + 1, str(ej))
                for (i, j, ei, ej) in bond_edits
            ],
            'used_fallback': bool(used_fallback),
            'numbering_ok': numbering_ok,
            'suspicious': suspicious,
            'reordered_target_xyz': _orca_build_xyz_from_symbols_coords(
                [target_symbols[int(v)] for v in np.asarray(best_mapping, dtype=int).tolist()],
                np.asarray(target_coords, dtype=float)[np.asarray(best_mapping, dtype=int)],
                comment='Reordered target',
            ),
            'aligned_reference_xyz': _orca_build_xyz_from_symbols_coords(
                ref_symbols,
                best_aligned if best_aligned is not None else ref_coords,
                comment='Aligned reference',
            ),
        }

    def _update_numbering_fix_button():
        idx = int(state.get('numbering_check_block_idx', state.get('xyz_view_idx', 0)))
        result = (state.get('numbering_check_results') or {}).get(idx) or {}
        has_fix = bool(
            state.get('numbering_check_active')
            and idx > 0
            and result.get('best_mapping')
            and result.get('best_mapping') != list(range(len(result.get('best_mapping') or [])))
            and result.get('reordered_target_xyz')
        )
        orca_apply_numbering_btn.disabled = not has_fix

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

    def _build_coord_block():
        xyz_blocks = parse_xyz_blocks(orca_coords.value)
        if xyz_blocks:
            # Use external XYZ file – coordinates are written to the job dir
            return (
                f'* xyzfile {orca_charge.value} {orca_multiplicity.value}'
                f' {xyz_blocks[0][0]}'
            )
        coords = strip_xyz_header(orca_coords.value)
        return f'* xyz {orca_charge.value} {orca_multiplicity.value}\n{coords}\n*'

    #: What ORCA calls a held bond, angle and dihedral.
    _ORCA_CONSTRAINT = {'distance': ('B', 2), 'angle': ('A', 3),
                        'dihedral': ('D', 4)}

    def _input_structure_atom_count():
        """How many atoms the structure this input reads actually has."""
        blocks = state.get('xyz_blocks') or []
        text = blocks[0][1] if blocks else orca_coords.value
        return len([row for row in strip_xyz_header(text).split('\n')
                    if len(row.split()) >= 4])

    def _build_constraints_block():
        """The coordinates held in the editor, in ORCA's own syntax.

        Only the ones held exactly -- Hold in "fix" mode. A "pull" is a spring
        the browser relaxes against while a structure is being dragged; it has
        no counterpart in a geometry optimisation, and writing it as a
        constraint would claim something nobody asked for.

        The atom numbers are the ones on the atoms in the viewer: ORCA counts
        from zero, and so does the numbering.
        """
        held = [c for c in (state.get('constraints') or [])
                if c.get('mode') == 'fix']
        if not held:
            return ''
        # A held value names atoms by number, and the numbers belong to the
        # structure it was set on. Held between atoms 0 and 2 of a water and
        # then asked of a two-atom CO, "{ B 0 2 1.5000 C }" reaches ORCA, which
        # stops on it. Anything naming an atom this input does not have is left
        # out and said out loud rather than written down.
        reads = _input_structure_atom_count()
        dropped = 0
        if reads:
            keep = [c for c in held
                    if all(0 <= int(i) < reads for i in (c.get('atoms') or []))]
            dropped = len(held) - len(keep)
            held = keep
        if not held:
            return (f'# {dropped} held value(s) name atoms this structure does '
                    f'not have, and are left out.\n') if dropped else ''
        lines = []
        for entry in held:
            word, wanted = _ORCA_CONSTRAINT.get(entry.get('kind'), ('', 0))
            atoms = list(entry.get('atoms') or [])
            if not word or len(atoms) != wanted:
                continue
            lines.append('    { %s %s %.4f C }' % (
                word, ' '.join(str(int(i)) for i in atoms),
                float(entry.get('value') or 0.0)))
        if not lines:
            return ''
        note = ''
        blocks = state.get('xyz_blocks') or []
        shown = int(state.get('xyz_view_idx', 0))
        if len(blocks) > 1 and shown != 0:
            # The input reads the first block; the editor was working on
            # another one. Saying so is the only honest thing to do -- the
            # numbers below name atoms of a structure ORCA will not see.
            note = (f'# Held in the editor on {blocks[shown][0]}, while this\n'
                    f'# input reads {blocks[0][0]}. Check the atom numbers.\n')
        if dropped:
            note += (f'# {dropped} more held value(s) name atoms this structure '
                     f'does not have, and are left out.\n')
        return note + '%geom Constraints\n' + '\n'.join(lines) + '\n  end\nend'

    def generate_orca_input():
        keyword_line = _build_keyword_line()
        pal_block = f'%pal\n  nprocs {orca_pal.value}\nend'
        maxcore_line = f'%maxcore {orca_maxcore.value}'
        output_block = _build_output_block()
        constraints_block = _build_constraints_block()
        coord_block = _build_coord_block()
        inp = f'{keyword_line}\n\n{pal_block}\n\n{maxcore_line}\n'
        if output_block:
            inp += f'\n{output_block}\n'
        if constraints_block:
            inp += f'\n{constraints_block}\n'
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
            orca_solvent.value = 'acetonitrile'
            orca_print_mos.value = False
            orca_print_basis.value = False
            orca_pal.value = 12
            orca_maxcore.value = 6000
            orca_slurm_time.value = '06:00:00'
            orca_path_files.value = ''
            orca_preview.value = ''
            state['extra_files'].clear()
            state['last_auto_keywords'] = ''
            state['numbering_check_active'] = False
            state['numbering_check_results'] = {}
            state['numbering_check_block_idx'] = 1
            state['numbering_view_step'] = 0
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
    _VIEWER_JS_TMPL = (
        '<div id="__DIV__" style="width:100%;height:560px;position:relative;margin:0;padding:0;"></div>\n'
        '<script>\n'
        'if(typeof $3Dmol==="undefined"){\n'
        '  var _s=document.createElement("script");\n'
        '  _s.src="https://3Dmol.org/build/3Dmol-min.js";\n'
        '  document.head.appendChild(_s);\n'
        '}\n'
        '(function(){\n'
        '  var tries=0;\n'
        '  function init(){\n'
        '    var el=document.getElementById("__DIV__");\n'
        '    if(!el||typeof $3Dmol==="undefined"){\n'
        '      tries++;if(tries<400)setTimeout(init,tries<40?50:250);return;\n'
        '    }\n'
        '    __RESET__\n'
        '    window._orcaBuildViewState=window._orcaBuildViewState||null;\n'
        '    var prev=window._orcaBuildViewer||null;\n'
        '    if(!__RESETFLAG__&&prev&&typeof prev.getView==="function"){\n'
        '      try{window._orcaBuildViewState=prev.getView();}catch(_e){}\n'
        '    }\n'
        '    var saved=window._orcaBuildViewState;\n'
        '    if(prev&&window.__delfinDisposeViewer){window.__delfinDisposeViewer(prev);}\n'
        '    var viewer=window.__delfinCreateViewer(el,__VIEWER_CONFIG__);\n'
        '    __MOUSE__\n'
        '    viewer.addModel(__XYZ__,"xyz");\n'
        '    viewer.setStyle({},__STYLE__);\n'
        '    __LABELS__\n'
        '    if(saved&&typeof viewer.setView==="function"){\n'
        '      try{viewer.setView(saved);}catch(_e){\n'
        '        viewer.zoomTo();viewer.center();viewer.zoom(__ZOOM__);\n'
        '      }\n'
        '    }else{\n'
        '      viewer.zoomTo();viewer.center();viewer.zoom(__ZOOM__);\n'
        '    }\n'
        '    viewer.render();\n'
        '    window._orcaBuildViewer=viewer;\n'
        '    __REGISTER__\n'
        '  }\n'
        '  setTimeout(init,0);\n'
        '})();\n'
        '</script>\n'
    )

    def _register_with_editor_js(viewer='viewer', element='el'):
        """Hand the viewer to the structure editor working in this scope.

        This is what the Submit tab's viewer does when it appears, and it is
        all the toolbar needs: the editor addresses a viewer by the scope it
        belongs to, never by a name of its own.
        """
        scope = json.dumps(orca_editor_scope)
        # Waited for, not assumed. The editor is installed by the kernel, and a
        # viewer can appear before it: the first structure of a session that
        # began with a SMILES is exactly that, since a SMILES in the box has no
        # viewer and so nothing asked for the editor yet. The introduction was
        # made once, eighty milliseconds later, and if the editor was not there
        # by then it was never made at all -- the picture turned in space and
        # nothing else worked, because the editor had never heard of it.
        return (
            'try{'
            'window._submitMolViewerByScope=window._submitMolViewerByScope||{};'
            f'window._submitMolViewerByScope[{scope}]={viewer};'
            f'(function(el){{var tries=0;'
            'var meet=function(){'
            'if(window.__delfinSubmitManip){'
            f'try{{window.__delfinSubmitManip.onViewerReady({scope},el);}}catch(_e){{}}'
            'return;}'
            'if(++tries<200)setTimeout(meet,80);};'
            f'setTimeout(meet,80);}})({element});'
            '}catch(_e){}'
        )

    # JS installed once per viewer: labels sit at the exact atom centre and are
    # drawn on top (inFront:true), so a number never drifts off its atom and is
    # never hidden by its OWN sphere -- no matter the zoom (no parallax).  To
    # still hide numbers of atoms that are *behind other atoms*, an occlusion
    # test runs once after an interaction settles. During mouse movement labels
    # remain visible with their last settled visibility, while no expensive
    # all-label pass runs in the render hot path. A small projected-coordinate
    # grid limits the settled occlusion pass to nearby candidates instead of
    # comparing every atom with every other atom.
    def _label_scale():
        return _structure_editor.scale_for_px(orca_editor.submit_label_size.value)

    def _labels_js(var='viewer'):
        """Atom numbers, gated by the preview's labels on/off toggle.

        This was 140 lines here and reachable from nowhere else -- the size
        control resolved the viewer as ``window._orcaBuildViewer``. Numbering
        belongs to a viewer, so it lives beside the rest of the editor now and
        the Submit tab has it too. The coordinates are no longer handed over:
        the numbers are read off the atoms the viewer holds, which is also
        what keeps them on the cores when those move.
        """
        if not orca_editor.submit_labels_btn.value:
            return ''
        return _structure_editor.show_atom_numbers_js(
            var=var, scale=_label_scale())

    def _update_molecule_js(xyz_data, label_js=''):
        """Script that shows another molecule in the viewer that is already there.

        Building a viewer means disposing a WebGL context and creating a new one,
        which is what browsing molecules and switching the numbers on and off
        used to do on every click.  Swapping the model inside the living viewer
        costs a render, keeps the camera exactly where the user left it, and
        does not consume one of the handful of contexts a browser grants.

        Returns '' when there is nothing live to update, so the caller falls
        back to a full render.
        """
        div_id = state.get('viewer_div')
        if not div_id or not state.get('viewer_live'):
            return ''
        profile = get_viewer_profile()
        if not profile['enabled']:
            return ''
        return (
            '(function(){\n'
            '  var tries=0;\n'
            '  function go(){\n'
            f'    var el=document.getElementById({json.dumps(div_id)});\n'
            '    var viewer=window._orcaBuildViewer;\n'
            '    if(!el||!viewer||typeof viewer.addModel!=="function"){\n'
            # a viewer displayed moments ago may still be starting up
            '      if(++tries<40){setTimeout(go,50);}\n'
            '      return;\n'
            '    }\n'
            '    try{\n'
            '      if(window.__delfinAtomNumbers)window.__delfinAtomNumbers.clear(viewer);\n'
            '      viewer.removeAllModels();\n'
            f'      viewer.addModel({json.dumps(xyz_data)},"xyz");\n'
            f'      viewer.setStyle({{}},{profile["style_js"]});\n'
            f'      {label_js}\n'
            '      viewer.render();\n'
            f'      {_register_with_editor_js("viewer", "el")}\n'
            '      var hs=viewer.__delfinInteractionEndHandlers||[];\n'
            '      for(var i=0;i<hs.length;i++){try{hs[i]();}catch(_e){}}\n'
            '    }catch(_err){}\n'
            '  }\n'
            '  go();\n'
            '})();'
        )

    def _show_molecule_in_place(full_xyz):
        """Show *full_xyz* without rebuilding the viewer.  True if it happened."""
        script = _update_molecule_js(full_xyz, _labels_js())
        if not script:
            return False
        ctx.run_js(script)
        _hand_to_editor(full_xyz)
        return True

    def _hand_to_editor(full_xyz):
        """Tell the editor which structure is on screen.

        Quietly: this is the tab saying what it has just drawn, not the editor
        saying it has changed something, and the two travel down the same wire.

        Stepping to another block also puts aside what the editor knows about
        the one being left and hands back what it knew about the one arriving.
        Without that, coming back to a structure meant perceiving it again from
        its coordinates -- with an atom where it had been dragged to, and the
        bond to its neighbour simply gone.
        """
        blocks = state.get('xyz_blocks') or []
        index = int(state.get('xyz_view_idx', 0))
        here = blocks[index][0] if 0 <= index < len(blocks) else ''
        there = state.get('editor_block')
        memory = state.setdefault('editor_memory', {})
        if there != here and there is not None:
            memory[there] = orca_editor.remember_structure()
        # The coordinates first, then the memory: giving the memory back while
        # the box still held the block being left made the editor perceive that
        # one again under the new block's name, and the bonding that was being
        # handed back was overwritten in the same breath.
        state['editor_quiet'] = True
        try:
            orca_editor_coords.value = full_xyz or ''
        finally:
            state['editor_quiet'] = False
        if there != here:
            orca_editor.restore_structure(memory.get(here))
            state['editor_block'] = here
        if full_xyz:
            orca_editor._ensure_manip_bootstrap()
        orca_editor._set_manip_toolbar_enabled(bool(full_xyz))
        if there != here and full_xyz:
            # The block on screen changed, so a running force field has to be
            # worked out again for it.
            orca_editor.structure_changed()

    def _take_structures(isomers, quick=False):
        """Every structure a conversion produced, as named blocks.

        A conversion hands over one structure or several -- the isomers of a
        SMILES, MANTA's ranked candidates, whatever was drawn. The Submit tab
        shows the first and steps through the rest; this tab keeps them all at
        once, in the layout it reads: a name, the atoms, a closing star. So
        they are written that way and the block stepper walks them, which is
        the stepper this tab already had.

        Except for the quick conversion, which answers with a structure rather
        than a set to choose from: that one writes plain coordinates, the way
        it always has here.
        """
        if not isomers:
            return False
        if quick:
            xyz_string, num_atoms, _label = isomers[0]
            rows = [row for row in (strip_xyz_header(xyz_string)
                                    or xyz_string).split('\n') if row.strip()]
            orca_coords.value = '%d\nConverted from SMILES\n%s' % (
                num_atoms or len(rows), '\n'.join(rows))
            return True
        # The quick embedding rides along at the end of a conformer set as a
        # fallback to step to. It is not a conformer, and a block called
        # quick.xyz beside conf-1 and conf-2 says it is one.
        named = [entry for entry in isomers
                 if str(entry[2] or '').strip().lower() != 'quick']
        if not named:
            named = isomers
        blocks = []
        used = {}
        for xyz_string, num_atoms, label in named:
            name = re.sub(r'[^A-Za-z0-9_.-]+', '-', str(label or 'structure')).strip('-.')
            if not name:
                name = 'structure'
            if not name.lower().endswith('.xyz'):
                name += '.xyz'
            # Two isomers may carry the same label; a block needs its own name.
            used[name] = used.get(name, 0) + 1
            if used[name] > 1:
                name = '%s-%d.xyz' % (name[:-4], used[name])
            body = strip_xyz_header(xyz_string) or xyz_string
            rows = [row for row in body.split('\n') if row.strip()]
            # Name, nothing after the semicolon, then the structure itself --
            # the same header every other block in this box carries.
            blocks.append('%s;\n%d\n%s\n%s\n*' % (
                name, num_atoms or len(rows), label or '', '\n'.join(rows)))
        orca_coords.value = '\n\n'.join(blocks)
        return True

    def _blocks_with_edit(full_xyz):
        """The coordinates box with the shown block replaced by *full_xyz*.

        The same rebuild Apply Numbering Fix does -- name, comment, atoms, the
        closing star -- so an edit leaves every other block and every header
        exactly as it found them.
        """
        records = _orca_parse_xyz_block_records(orca_coords.value)
        idx = int(state.get('xyz_view_idx', 0))
        if not records or not 0 <= idx < len(records):
            # A box that was never written as named blocks holds one plain
            # XYZ, header and all. Handing back the bare atom lines took that
            # header away on every edit, and after an optimisation the box read
            # as a list of coordinates with no count and no comment.
            return full_xyz
        records[idx]['full_xyz'] = full_xyz
        rebuilt = []
        for record in records:
            label = record.get('raw_name') or record.get('filename') or 'block'
            suffix = record.get('suffix_comment', '')
            header = f'{label};{suffix}' if suffix else f'{label};'
            rebuilt.append(f'{header}\n{record["full_xyz"].strip()}\n*')
        return '\n\n'.join(rebuilt)

    def _orca_editor_wrote(change=None):
        """The editor has changed the structure. Put it in the block on screen.

        Rewriting the coordinates box would normally throw the preview away and
        step back to the first block, which is not what dragging an atom asks
        for -- so the box is written quietly and the block list corrected by
        hand.
        """
        if state.get('editor_quiet'):
            return
        if state.get('numbering_check_active'):
            # A comparison is on screen -- an overlay, the reference turned to
            # lie over the target, or a proposed renumbering. None of them is a
            # block, and a write that arrives while one is shown (a drag still
            # in flight when the check began) would land in the target.
            return
        text = orca_editor_coords.value
        atoms = strip_xyz_header(text)
        if not atoms.strip():
            return
        lines = [line for line in atoms.split('\n') if line.strip()]
        full_xyz = '%d\nEdited in DELFIN viewer\n%s' % (len(lines), atoms)
        # The browser moved the atoms and is already showing them; only the
        # coordinates have to catch up. Redrawing here would take the model
        # apart and perceive its bonds again, which is exactly what the editor
        # has been working against.
        drawn_already = bool(state.pop('manip_inflight', False))
        state['editor_quiet'] = True
        try:
            orca_coords.value = _blocks_with_edit(full_xyz)
            blocks = parse_xyz_blocks(orca_coords.value) or []
            if blocks:
                state['xyz_blocks'] = blocks
                state['xyz_view_idx'] = min(state.get('xyz_view_idx', 0),
                                            len(blocks) - 1)
        finally:
            state['editor_quiet'] = False
        if not drawn_already and not _show_molecule_in_place(full_xyz):
            _refresh_mol_view(reset_view=False)

    def _viewer_html(xyz_data, label_js='', reset_view=False):
        """Build a self-contained HTML block that renders xyz_data in a $3Dmol viewer.

        The viewer resets on new coordinates, but preserves the same camera view
        while switching between multiple blocks of the same system.
        """
        div_id = 'orca-mol-' + uuid.uuid4().hex[:10]
        profile = get_viewer_profile()
        if not profile['enabled']:
            state['viewer_live'] = False
            return viewer_disabled_html()
        state['viewer_div'] = div_id
        state['viewer_live'] = True
        mouse_js = patch_viewer_mouse_controls_js('viewer', 'el')
        zoom = str(DEFAULT_3DMOL_ZOOM if DEFAULT_3DMOL_ZOOM is not None else 0.9)
        reset_js = 'window._orcaBuildViewState=null;' if reset_view else ''
        html = (
            _VIEWER_JS_TMPL
            .replace('__DIV__', div_id)
            .replace('__RESET__', reset_js)
            .replace('__RESETFLAG__', 'true' if reset_view else 'false')
            .replace('__MOUSE__', mouse_js)
            .replace('__XYZ__', json.dumps(xyz_data))
            .replace('__STYLE__', profile['style_js'])
            .replace('__VIEWER_CONFIG__', profile['viewer_config_js'])
            .replace('__LABELS__', label_js)
            .replace('__REGISTER__', _register_with_editor_js())
            .replace('__ZOOM__', zoom)
        )
        return html

    def _overlay_viewer_html(reference_xyz, target_xyz, reset_view=False):
        profile = get_viewer_profile()
        # two models in one viewer: not the single molecule an in-place swap
        # would be updating
        state['viewer_live'] = False
        if not profile['enabled']:
            return viewer_disabled_html()
        div_id = 'orca-overlay-' + uuid.uuid4().hex[:10]
        mouse_js = patch_viewer_mouse_controls_js('viewer', 'el')
        zoom = str(DEFAULT_3DMOL_ZOOM if DEFAULT_3DMOL_ZOOM is not None else 0.9)
        reset_js = 'window._orcaBuildViewState=null;' if reset_view else ''
        target_style_js = molecule_view_style_js(profile['style'], color='#1f5fff')
        reference_style_js = molecule_view_style_js(profile['style'], color='#d32f2f')
        return (
            '<div id="' + div_id + '" style="width:100%;height:560px;position:relative;margin:0;padding:0;"></div>\n'
            '<script>\n'
            'if(typeof $3Dmol==="undefined"){\n'
            '  var _s=document.createElement("script");\n'
            '  _s.src="https://3Dmol.org/build/3Dmol-min.js";\n'
            '  document.head.appendChild(_s);\n'
            '}\n'
            '(function(){\n'
            '  var tries=0;\n'
            '  function init(){\n'
            f'    var el=document.getElementById("{div_id}");\n'
            '    if(!el||typeof $3Dmol==="undefined"){\n'
            '      tries++;if(tries<400)setTimeout(init,tries<40?50:250);return;\n'
            '    }\n'
            f'    {reset_js}\n'
            '    window._orcaBuildViewState=window._orcaBuildViewState||null;\n'
            '    var prev=window._orcaBuildViewer||null;\n'
            '    if(!' + ('true' if reset_view else 'false') + '&&prev&&typeof prev.getView==="function"){\n'
            '      try{window._orcaBuildViewState=prev.getView();}catch(_e){}\n'
            '    }\n'
            '    var saved=window._orcaBuildViewState;\n'
            '    if(prev&&window.__delfinDisposeViewer){window.__delfinDisposeViewer(prev);}\n'
            f'    var viewer=window.__delfinCreateViewer(el,{profile["viewer_config_js"]});\n'
            f'    {mouse_js}\n'
            f'    viewer.addModel({json.dumps(target_xyz)},"xyz");\n'
            f'    viewer.setStyle({{model:0}},{target_style_js});\n'
            f'    viewer.addModel({json.dumps(reference_xyz)},"xyz");\n'
            f'    viewer.setStyle({{model:1}},{reference_style_js});\n'
            '    if(saved&&typeof viewer.setView==="function"){\n'
            '      try{viewer.setView(saved);}catch(_e){viewer.zoomTo();viewer.center();viewer.zoom(' + zoom + ');}\n'
            '    }else{\n'
            '      viewer.zoomTo();viewer.center();viewer.zoom(' + zoom + ');\n'
            '    }\n'
            '    viewer.render();\n'
            '    window._orcaBuildViewer=viewer;\n'
            f'    {_register_with_editor_js()}\n'
            '  }\n'
            '  setTimeout(init,0);\n'
            '})();\n'
            '</script>\n'
        )

    def _numbering_check_view_html(reference_xyz, target_xyz, reordered_target_xyz, step, reset_view=False):
        step = int(step)
        if step == 1:
            return _viewer_html(
                reference_xyz,
                _labels_js(),
                reset_view=reset_view,
            )
        if step == 2:
            return _viewer_html(
                reordered_target_xyz,
                _labels_js(),
                reset_view=reset_view,
            )
        return _overlay_viewer_html(reference_xyz, target_xyz, reset_view=reset_view)

    def _update_nav_label():
        blocks = state['xyz_blocks']
        if state.get('numbering_check_active'):
            block_idx = int(state.get('numbering_check_block_idx', 1))
            step = int(state.get('numbering_view_step', 0))
            labels = [
                'Overlay',
                'Aligned reference',
                'Reordered target',
            ]
            block_name = blocks[block_idx][0] if 0 <= block_idx < len(blocks) else 'Comparison'
            orca_mol_nav_label.value = (
                f'<span style="font-size:12px;">'
                f'{step + 1}&thinsp;/&thinsp;3: {labels[step]} for {block_name}'
                f'</span>'
            )
            orca_mol_prev_btn.layout.display = ''
            orca_mol_next_btn.layout.display = ''
            orca_mol_nav_row.layout.display = ''
            return

        n = len(blocks)
        if n > 1:
            idx = state['xyz_view_idx']
            orca_mol_nav_label.value = (
                f'<span style="font-size:12px;">'
                f'{idx + 1}&thinsp;/&thinsp;{n}: {blocks[idx][0]}'
                f'</span>'
            )
            orca_mol_prev_btn.layout.display = ''
            orca_mol_next_btn.layout.display = ''
            orca_mol_nav_row.layout.display = ''
        else:
            # One structure, or one that never came in named blocks at all:
            # hide the stepper, keep the row. It carries the fullscreen button,
            # and a structure can be made large whether or not somebody wrote a
            # "name;" over it -- paste a bare XYZ and the row disappeared, and
            # with it the only way to enlarge the editor.
            orca_mol_nav_label.value = ''
            orca_mol_prev_btn.layout.display = 'none'
            orca_mol_next_btn.layout.display = 'none'
            # Only the stepper is left in this row -- the fullscreen button
            # sits in the toolbar now -- so with nothing to step to there is
            # nothing to show.
            orca_mol_nav_row.layout.display = 'none'

    def _show_in_viewer(*items):
        """Put *items* in the preview, by assignment rather than by capture.

        ``with output: display(...)`` only reaches the widget while the kernel
        is running the cell it belongs to. A conversion answers from a thread,
        through the interface loop, and there is no cell there: the coordinates
        arrived in the box, the preview was asked to draw them, and nothing
        appeared -- the viewer kept whatever it had, which after a SMILES was a
        one-atom model of the letters. Assigning the outputs works wherever the
        call comes from, which is how the Submit tab has always done it.
        """
        orca_mol_output.outputs = tuple(items)

    def _as_html(markup):
        return {'output_type': 'display_data',
                'data': {'text/html': markup}, 'metadata': {}}

    def _as_text(message):
        return {'output_type': 'stream', 'name': 'stdout',
                'text': str(message) + '\n'}

    def _refresh_mol_view(reset_view=False):
        """Re-render the molecule viewer, preserving orientation unless *reset_view*."""
        blocks = state['xyz_blocks']
        _update_nav_label()
        _update_numbering_fix_button()
        orca_mol_output.layout.height = '560px'
        orca_mol_output.layout.min_height = '560px'
        state['viewer_live'] = False   # nothing on screen until we draw it
        # Ask for the editor before drawing the thing it has to be told about.
        if blocks or strip_xyz_header(orca_coords.value).strip():
            orca_editor._ensure_manip_bootstrap()
        if blocks:
            idx = state['xyz_view_idx']
            _block_name, full_xyz = blocks[idx]
            try:
                overlay_idx = int(state.get('numbering_check_block_idx', 1))
                overlay_result = (state.get('numbering_check_results') or {}).get(overlay_idx)
                if (
                    state.get('numbering_check_active')
                    and overlay_idx > 0
                    and overlay_result
                    and overlay_result.get('aligned_reference_xyz')
                ):
                    target_xyz = blocks[overlay_idx][1]
                    reordered_target_xyz = overlay_result.get('reordered_target_xyz') or target_xyz
                    step = int(state.get('numbering_view_step', 0))
                    _show_in_viewer(_as_html(_numbering_check_view_html(
                        overlay_result['aligned_reference_xyz'],
                        target_xyz,
                        reordered_target_xyz,
                        step,
                        reset_view=reset_view,
                    )))
                    # All three are comparison pictures, and none of them is
                    # a block. The overlay is two structures at once; the
                    # aligned reference is the reference turned to lie over the
                    # target, and the reordered target is a proposal that has
                    # not been applied. Editing any of them wrote into the
                    # target block: drag a hydrogen in the aligned reference
                    # and the target came back holding the reference's
                    # geometry. So they are shown and numbered, and the toolbar
                    # waits until there is a structure to edit again.
                    _hand_to_editor('')
                    orca_editor._set_mol_status(*{
                        0: ('Overlay: reference in red, target in blue.',),
                        1: ('The reference, turned to lie over the target.',),
                        2: ('The target as it would be renumbered. Apply '
                            'Numbering Fix writes it into the block.',),
                    }.get(step, ()) + (
                        'A comparison, not a block -- editing waits until you '
                        'are back on the structure itself.',))
                else:
                    _show_in_viewer(_as_html(
                        _viewer_html(full_xyz, _labels_js(), reset_view=reset_view)))
                    _hand_to_editor(full_xyz)
            except Exception as e:
                _show_in_viewer(_as_text(f'Could not visualize: {e}'))
            return
        raw = orca_coords.value.strip()
        if not raw:
            _show_in_viewer(_as_text('Paste XYZ coordinates to see 3D preview.'))
            _hand_to_editor('')
            return
        # A SMILES is not coordinates. Read as some, "c1ccccc1" became a
        # one-atom model of its own first letter, and after a conversion had
        # plainly filled the box the viewer still showed that atom -- which is
        # what "it does not show it" looked like from the outside.
        if clean_input_data(raw)[1] == 'smiles':
            _show_in_viewer(_as_text(
                'SMILES detected. Use CONVERT SMILES, QUICK CONVERT or '
                'CONVERT SMILES + UFF to turn it into coordinates.'))
            _hand_to_editor('')
            return
        coords = strip_xyz_header(raw)
        if not coords:
            _show_in_viewer(_as_text('No valid coordinates.'))
            _hand_to_editor('')
            return
        try:
            atom_lines = [l for l in coords.split('\n') if l.strip()]
            xyz_data = f'{len(atom_lines)}\nORCA Builder Preview\n{coords}'
            _show_in_viewer(_as_html(
                _viewer_html(xyz_data, _labels_js(), reset_view=reset_view)))
            _hand_to_editor(xyz_data)
        except Exception as e:
            _show_in_viewer(_as_text(f'Could not visualize: {e}'))

    def update_orca_molecule_view(change=None):
        # An edit from the editor has already put itself in the box and told
        # the viewer; starting over here would step back to the first block
        # and reset the camera in the middle of a drag.
        if state.get('editor_quiet'):
            return
        # Coordinates the editor has not seen: everything it knew about the
        # last ones goes, not only the switches. A held value names atoms of
        # the structure it was set on -- kept across a paste, a bond held
        # between atoms 0 and 2 of a water was written into the input for a
        # two-atom CO as "{ B 0 2 1.5000 C }", and ORCA stops on it.
        orca_editor.restore_structure(None)
        state.pop('editor_memory', None)
        state.pop('editor_block', None)
        state['xyz_blocks'] = parse_xyz_blocks(orca_coords.value) or []
        state['xyz_view_idx'] = 0
        state['numbering_check_active'] = False
        state['numbering_check_results'] = {}
        state['numbering_check_block_idx'] = 1
        state['numbering_view_step'] = 0
        _refresh_mol_view(reset_view=True)  # new coords → reset camera

    def on_mol_prev(btn):
        blocks = state['xyz_blocks']
        if not blocks:
            return
        if state.get('numbering_check_active'):
            state['numbering_view_step'] = (int(state.get('numbering_view_step', 0)) - 1) % 3
            _update_nav_label()
            _refresh_mol_view(reset_view=False)
            return
        state['xyz_view_idx'] = (state['xyz_view_idx'] - 1) % len(blocks)
        _update_nav_label()
        _show_next_molecule()

    def on_mol_next(btn):
        blocks = state['xyz_blocks']
        if not blocks:
            return
        if state.get('numbering_check_active'):
            state['numbering_view_step'] = (int(state.get('numbering_view_step', 0)) + 1) % 3
            _update_nav_label()
            _refresh_mol_view(reset_view=False)
            return
        state['xyz_view_idx'] = (state['xyz_view_idx'] + 1) % len(blocks)
        _update_nav_label()
        _show_next_molecule()

    def _show_next_molecule():
        """Browse to the selected block, in place where the viewer allows it."""
        blocks = state['xyz_blocks']
        idx = state['xyz_view_idx']
        if 0 <= idx < len(blocks) and _show_molecule_in_place(blocks[idx][1]):
            return
        _refresh_mol_view(reset_view=False)  # keep orientation

    def _current_xyz_for_copy():
        xyz_blocks = parse_xyz_blocks(orca_coords.value)
        if xyz_blocks:
            idx = min(max(int(state.get('xyz_view_idx', 0)), 0), len(xyz_blocks) - 1)
            return xyz_blocks[idx][1]

        raw = orca_coords.value.strip()
        if not raw:
            return ''
        coords = strip_xyz_header(raw)
        if not coords:
            return ''
        atom_lines = [line for line in coords.split('\n') if line.strip()]
        if not atom_lines:
            return ''
        return f'{len(atom_lines)}\nORCA Builder Coordinates\n' + '\n'.join(atom_lines)

    def handle_orca_copy_coordinates(button):
        xyz_text = _current_xyz_for_copy()
        if not xyz_text:
            with orca_output:
                clear_output(wait=True)
                print('No coordinates available to copy.')
            return
        text_payload = json.dumps(str(xyz_text))
        ctx.run_js(
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
        with orca_output:
            clear_output(wait=True)
            print('Copied coordinates as XYZ to clipboard.')

    orca_copy_coords_btn.on_click(handle_orca_copy_coordinates)

    def handle_orca_check_numbering(button):
        xyz_blocks = parse_xyz_blocks(orca_coords.value) or []
        if len(xyz_blocks) < 2:
            with orca_output:
                clear_output(wait=True)
                print('Check Numbering needs at least two named XYZ blocks.')
            return

        # The comparison holds the kernel, so say so before it starts: the
        # message reaches the browser while the mapping is still running.
        with orca_output:
            clear_output(wait=True)
            print(f'Comparing {len(xyz_blocks) - 1} block(s) against the reference...')
        orca_check_numbering_btn.disabled = True
        try:
            _run_numbering_check(xyz_blocks)
        finally:
            orca_check_numbering_btn.disabled = False

    def _run_numbering_check(xyz_blocks):
        ref_name, ref_xyz = xyz_blocks[0]
        try:
            ref_symbols, ref_coords = _orca_parse_xyz_symbols_coords(ref_xyz)
        except Exception as exc:
            with orca_output:
                clear_output(wait=True)
                print(f'Reference block could not be parsed: {exc}')
            return

        results = {}
        lines = [f'Reference: {ref_name}']
        for idx in range(1, len(xyz_blocks)):
            name, xyz_text = xyz_blocks[idx]
            try:
                target_symbols, target_coords = _orca_parse_xyz_symbols_coords(xyz_text)
                result = _orca_check_numbering_pair(ref_symbols, ref_coords, target_symbols, target_coords)
                results[idx] = result
                n_edits = int(result.get('n_bond_edits', 0))
                mapping = list(result.get('best_mapping') or [])
                is_identity = mapping == list(range(len(mapping)))
                if n_edits > 0 and is_identity:
                    verdict = (
                        f'reaction mapping ({n_edits} bond edit(s)) - '
                        'numbering already matches'
                    )
                elif n_edits > 0:
                    verdict = f'reaction mapping ({n_edits} bond edit(s)) - reorder available'
                elif result.get('suspicious'):
                    verdict = 'numbering could be wrong'
                elif result.get('numbering_ok'):
                    verdict = 'numbering looks consistent'
                else:
                    verdict = 'mapping differs, please inspect overlay'
                direct_text = (
                    f"{result['direct_rmsd']:.4f} A"
                    if result.get('direct_rmsd') is not None
                    else 'n/a'
                )
                lines.append(
                    f'- {name}: {verdict} '
                    f'(direct={direct_text}, best={result["best_rmsd"]:.4f} A, '
                    f'method={result["best_source"]})'
                )
                for a_idx, a_el, b_idx, b_el in result.get('bond_edits', []):
                    lines.append(
                        f'    formed/broken bond: {a_el}{a_idx}-{b_el}{b_idx}'
                    )
            except Exception as exc:
                results[idx] = {'error': str(exc)}
                lines.append(f'- {name}: not comparable ({exc})')

        state['numbering_check_active'] = True
        state['numbering_check_results'] = results
        if len(xyz_blocks) > 1:
            state['xyz_view_idx'] = min(max(state.get('xyz_view_idx', 1), 1), len(xyz_blocks) - 1)
        state['numbering_check_block_idx'] = int(state.get('xyz_view_idx', 1)) if len(xyz_blocks) > 1 else 1
        state['numbering_view_step'] = 0
        _refresh_mol_view(reset_view=True)
        with orca_output:
            clear_output(wait=True)
            print('\n'.join(lines))
            print()
            print('Molecule Preview cycles through overlay, aligned reference, and reordered target for the selected comparison block.')

    orca_check_numbering_btn.on_click(handle_orca_check_numbering)

    def handle_orca_apply_numbering_fix(button):
        xyz_records = _orca_parse_xyz_block_records(orca_coords.value)
        idx = int(state.get('numbering_check_block_idx', state.get('xyz_view_idx', 0)))
        result = (state.get('numbering_check_results') or {}).get(idx) or {}
        if not xyz_records or idx <= 0 or idx >= len(xyz_records):
            with orca_output:
                clear_output(wait=True)
                print('No checked block is selected for numbering fix.')
            return
        reordered_xyz = str(result.get('reordered_target_xyz') or '').strip()
        if not reordered_xyz:
            with orca_output:
                clear_output(wait=True)
                print('No numbering fix is available for the selected block.')
            return

        xyz_records[idx]['full_xyz'] = reordered_xyz
        rebuilt = []
        for record in xyz_records:
            label = record.get('raw_name') or record.get('filename') or 'block'
            suffix = record.get('suffix_comment', '')
            header = f'{label};{suffix}' if suffix else f'{label};'
            rebuilt.append(f'{header}\n{record["full_xyz"].strip()}\n*')
        # Written quietly, so the comparison stays up: the reordered block is
        # the thing to look at after a fix, and starting the tab over would
        # drop it back to the first structure with nothing to compare against.
        state['editor_quiet'] = True
        try:
            orca_coords.value = '\n\n'.join(rebuilt)
            state['xyz_blocks'] = parse_xyz_blocks(orca_coords.value) or []
        finally:
            state['editor_quiet'] = False
        result['reordered_target_xyz'] = reordered_xyz
        _refresh_mol_view(reset_view=False)
        with orca_output:
            clear_output(wait=True)
            print(f'Applied numbering fix to block {idx + 1}.')
            print('Please inspect the Molecule Preview again to confirm the reordered block looks correct.')

    orca_apply_numbering_btn.on_click(handle_orca_apply_numbering_fix)

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
        # Match the entire %pal block in either form: inline (%pal nprocs N end)
        # or multi-line (%pal\n nprocs N\nend), with optional '=' syntax. (?s)
        # makes '.' span newlines so the multi-line form is also covered.
        text = re.sub(r'(?is)%pal\b.*?\bend\b', new_pal, text, count=1)
        text = re.sub(r'^%maxcore\s*=?\s*\d+', f'%maxcore {orca_maxcore.value}',
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
            text = re.sub(
                r'(\n\* xyz(?:file)? )',
                f'\n{new_output}\n\\1', text, count=1,
            )

        # What the editor holds, kept in step the same way %output is: the
        # preview is the user's to edit, so the block is replaced where it is
        # rather than the whole input written again.
        new_geom = _build_constraints_block()
        held = re.search(r'%geom\b.*?\n\s*end\s*\nend', text,
                         flags=re.DOTALL | re.IGNORECASE)
        if held:
            if new_geom:
                text = (text[:held.start()] + new_geom + text[held.end():])
            else:
                text = (text[:held.start()].rstrip('\n') + '\n\n'
                        + text[held.end():].lstrip('\n'))
        elif new_geom:
            text = re.sub(r'(\n\* xyz(?:file)? )',
                          f'\n{new_geom}\n\\1', text, count=1)

        xyz_blocks = parse_xyz_blocks(orca_coords.value)
        if xyz_blocks:
            new_coord = (
                f'* xyzfile {orca_charge.value} {orca_multiplicity.value}'
                f' {xyz_blocks[0][0]}'
            )
        else:
            coords = strip_xyz_header(orca_coords.value)
            new_coord = f'* xyz {orca_charge.value} {orca_multiplicity.value}\n{coords}\n*'
        # Replace existing coord block – handles both xyzfile (single line) and
        # inline xyz (multi-line) formats, including switches between the two.
        new_text = re.sub(
            r'\*\s*xyzfile\s+[-\d]+\s+\d+\s+\S+',
            new_coord, text, count=1, flags=re.MULTILINE | re.IGNORECASE,
        )
        if new_text == text:
            new_text = re.sub(
                r'\*\s*xyz\s+[-\d]+\s+\d+\s*\n.*?\n\*',
                new_coord, text, count=1, flags=re.DOTALL,
            )
        text = new_text
        orca_preview.value = sanitize_orca_input(text)

    def _save_orca_job():
        """Save .inp, XYZ files, and extra files to the job directory.

        Returns ``(job_dir, safe_job_name, inp_content, saved_files)`` on
        success, or ``None`` if validation failed (error already printed).
        """
        job_name = orca_job_name.value.strip()
        if not job_name:
            print('Error: Job name cannot be empty!')
            return None

        safe_job_name = ''.join(c for c in job_name if c.isalnum() or c in ('_', '-'))
        if not safe_job_name:
            print('Error: Job name contains only invalid characters!')
            return None

        preview_content = orca_preview.value.strip()
        if preview_content:
            inp_content = preview_content
        else:
            coords = strip_xyz_header(orca_coords.value)
            if not coords:
                print('Error: Coordinates or INP preview cannot be empty!')
                return None
            inp_content = generate_orca_input()

        inp_content = sanitize_orca_input(inp_content)

        job_dir = ctx.calc_dir / safe_job_name
        job_dir.mkdir(parents=True, exist_ok=True)

        inp_path = job_dir / f'{safe_job_name}.inp'
        inp_path.write_text(inp_content)

        # Write named XYZ files (from name.xyz;...* blocks) to job dir
        xyz_blocks = parse_xyz_blocks(orca_coords.value)
        if xyz_blocks:
            for filename, xyz_content in xyz_blocks:
                (job_dir / filename).write_text(xyz_content)

        saved_files = save_uploaded_files(job_dir)
        return job_dir, safe_job_name, inp_content, saved_files

    def handle_orca_save(button):
        with orca_output:
            clear_output()
            result = _save_orca_job()
            if result is None:
                return
            job_dir, safe_job_name, inp_content, saved_files = result
            print(f'Saved to {job_dir}')
            print(f'  {safe_job_name}.inp')
            if saved_files:
                for fn in saved_files:
                    print(f'  {fn}')

    def handle_orca_submit(button):
        with orca_output:
            clear_output()
            result = _save_orca_job()
            if result is None:
                return
            job_dir, safe_job_name, inp_content, saved_files = result

            pal_used, maxcore_used = parse_inp_resources(inp_content)
            if pal_used is None:
                pal_used = orca_pal.value
            if maxcore_used is None:
                maxcore_used = orca_maxcore.value

            submit_result = ctx.backend.submit_orca(
                job_dir=job_dir, job_name=safe_job_name,
                inp_file=f'{safe_job_name}.inp',
                time_limit=orca_slurm_time.value,
                pal=pal_used, maxcore=maxcore_used,
            )

            if submit_result.returncode == 0:
                job_id = submit_result.stdout.strip().split()[-1] if submit_result.stdout.strip() else '(unknown)'
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
                print(submit_result.stderr or submit_result.stdout)

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
    orca_save_btn.on_click(handle_orca_save)
    orca_submit_btn.on_click(handle_orca_submit)
    orca_file_upload.observe(update_uploaded_files_label, names='value')

    for w in [orca_method, orca_job_type, orca_basis, orca_dispersion, orca_ri,
              orca_aux_basis, orca_charge, orca_multiplicity, orca_pal, orca_maxcore,
              orca_solvation_type,
              orca_solvent, orca_print_mos, orca_print_basis, orca_autoaux]:
        w.observe(update_orca_preview, names='value')

    # The coordinates box, in this order: what the editor knew about the last
    # structure goes first, and only then is the input written. The other way
    # round, the input was built from held values that belonged to coordinates
    # that had just been replaced.
    orca_coords.observe(update_orca_molecule_view, names='value')
    orca_coords.observe(update_orca_preview, names='value')
    orca_mol_prev_btn.on_click(on_mol_prev)
    orca_mol_next_btn.on_click(on_mol_next)
    orca_editor_coords.observe(_orca_editor_wrote, names='value')
    # The list of held coordinates is rebuilt whenever one is added, retuned or
    # deleted, so this is where the input preview hears about it.
    orca_editor.submit_constraint_dd.observe(
        lambda _change: update_orca_preview(), names='options')
    update_orca_molecule_view()
    update_orca_preview()
    state['last_auto_keywords'] = _build_keyword_line()

    # -- config templates ----------------------------------------------
    # A named template stores the level of theory, keys and resources plus the
    # full editable INP preview body (see inp_body below), so a saved setup can
    # be re-applied to a new system. Charge, multiplicity and coordinates are
    # intentionally NOT stored (they are system specific). This list is the
    # single source for the structured collect/apply.
    _TEMPLATE_FIELDS = [
        ('method', orca_method), ('job_type', orca_job_type), ('basis', orca_basis),
        ('dispersion', orca_dispersion), ('ri', orca_ri), ('aux_basis', orca_aux_basis),
        ('autoaux', orca_autoaux), ('solvation_type', orca_solvation_type),
        ('solvent', orca_solvent), ('print_mos', orca_print_mos),
        ('print_basis', orca_print_basis), ('pal', orca_pal),
        ('maxcore', orca_maxcore), ('slurm_time', orca_slurm_time),
    ]

    orca_template_dd = widgets.Dropdown(
        options=[], description='Template:',
        layout=widgets.Layout(width='300px'), style=ws,
    )
    orca_template_load_btn = widgets.Button(
        description='LOAD', button_style='primary',
        layout=widgets.Layout(width='90px'),
    )
    orca_template_save_btn = widgets.Button(
        description='SAVE TEMPLATE', button_style='info',
        layout=widgets.Layout(width='150px'),
    )
    orca_template_delete_btn = widgets.Button(
        description='DELETE', button_style='danger',
        layout=widgets.Layout(width='90px'),
    )
    orca_template_status = widgets.HTML(value='')

    orca_template_name = widgets.Text(
        value='', placeholder='Template name', description='Name:',
        layout=widgets.Layout(width='320px'), style=ws,
    )
    orca_template_save_confirm_btn = widgets.Button(
        description='SAVE', button_style='success',
        layout=widgets.Layout(width='90px'),
    )
    orca_template_save_cancel_btn = widgets.Button(
        description='CANCEL', layout=widgets.Layout(width='90px'),
    )
    orca_template_save_dialog = widgets.VBox(
        [
            widgets.HTML('<b>Save current settings as template:</b>'),
            widgets.HBox(
                [orca_template_name, orca_template_save_confirm_btn,
                 orca_template_save_cancel_btn],
                layout=widgets.Layout(gap='6px', align_items='center', flex_wrap='wrap'),
            ),
        ],
        layout=widgets.Layout(
            display='none', flex_flow='column', border='1px solid #1976d2',
            border_radius='6px', padding='8px', margin='4px 0', gap='4px',
        ),
    )

    orca_template_delete_prompt = widgets.HTML(value='')
    orca_template_delete_confirm_btn = widgets.Button(
        description='DELETE', button_style='danger',
        layout=widgets.Layout(width='90px'),
    )
    orca_template_delete_cancel_btn = widgets.Button(
        description='CANCEL', layout=widgets.Layout(width='90px'),
    )
    orca_template_delete_dialog = widgets.VBox(
        [
            orca_template_delete_prompt,
            widgets.HBox(
                [orca_template_delete_confirm_btn, orca_template_delete_cancel_btn],
                layout=widgets.Layout(gap='6px', align_items='center'),
            ),
        ],
        layout=widgets.Layout(
            display='none', flex_flow='column', border='1px solid #d32f2f',
            border_radius='6px', padding='8px', margin='4px 0', gap='4px',
        ),
    )

    def _set_template_status(message='', color='#555'):
        orca_template_status.value = (
            f'<span style="color:{color};">{message}</span>' if message else ''
        )

    def _show_save_dialog(show):
        orca_template_save_dialog.layout.display = 'flex' if show else 'none'

    def _show_delete_dialog(show):
        orca_template_delete_dialog.layout.display = 'flex' if show else 'none'

    def _refresh_template_dd(select=None):
        names = sorted(load_orca_templates().keys())
        current = select if select is not None else orca_template_dd.value
        orca_template_dd.options = names
        if names:
            orca_template_dd.value = current if current in names else names[0]
        orca_template_load_btn.disabled = not names
        orca_template_delete_btn.disabled = not names

    def _collect_template_payload():
        payload = {key: widget.value for key, widget in _TEMPLATE_FIELDS}
        # Also capture whatever the user typed directly into the editable INP
        # preview (manual keyword suffixes, %-blocks, extra lines), minus the
        # coordinate block, so hand edits survive a save/load round-trip.
        payload['inp_body'] = strip_coord_block(orca_preview.value)
        return payload

    def _apply_template_payload(payload):
        skipped = []
        if not isinstance(payload, dict):
            return skipped
        state['is_resetting'] = True
        try:
            for key, widget in _TEMPLATE_FIELDS:
                if key not in payload:
                    continue
                val = payload[key]
                if isinstance(widget, widgets.Dropdown):
                    if val not in widget.options:
                        skipped.append(f'{key}={val!r}')
                        continue
                    widget.value = val
                elif isinstance(widget, widgets.Checkbox):
                    widget.value = bool(val)
                elif isinstance(widget, widgets.IntText):
                    try:
                        widget.value = int(val)
                    except (TypeError, ValueError):
                        skipped.append(f'{key}={val!r}')
                else:
                    widget.value = '' if val is None else str(val)
        finally:
            state['is_resetting'] = False
        body = payload.get('inp_body')
        if isinstance(body, str) and body.strip():
            # Restore the manually edited settings body verbatim and append a
            # fresh coordinate block from the CURRENT coordinates (empty until
            # the user pastes new ones). Set last_auto_keywords so subsequent
            # field edits keep detecting the manual keyword suffix correctly.
            state['last_auto_keywords'] = _build_keyword_line()
            orca_preview.value = f'{body.rstrip()}\n\n{_build_coord_block()}\n'
        else:
            # No stored body (older template) -> full regen from the fields.
            orca_preview.value = ''
            update_orca_preview()
        return skipped

    def _on_template_save_click(_btn=None):
        _show_delete_dialog(False)
        orca_template_name.value = orca_template_dd.value or ''
        _show_save_dialog(True)

    def _on_template_save_confirm(_btn=None):
        name = (orca_template_name.value or '').strip()
        if not name:
            _set_template_status('Template name must not be empty.', '#d32f2f')
            return
        try:
            save_orca_template(name, _collect_template_payload())
        except Exception as exc:
            _set_template_status(f'Could not save template: {html.escape(str(exc))}', '#d32f2f')
            return
        _show_save_dialog(False)
        _refresh_template_dd(select=name)
        _set_template_status(f"Template '{html.escape(name)}' saved.", '#2e7d32')

    def _on_template_save_cancel(_btn=None):
        _show_save_dialog(False)

    def _on_template_load_click(_btn=None):
        name = orca_template_dd.value
        if not name:
            _set_template_status('No template selected.', '#d32f2f')
            return
        templates = load_orca_templates()
        if name not in templates:
            _refresh_template_dd()
            _set_template_status(f"Template '{html.escape(str(name))}' no longer exists.", '#d32f2f')
            return
        skipped = _apply_template_payload(templates[name])
        if skipped:
            _set_template_status(
                f"Loaded template '{html.escape(name)}' "
                f"(skipped invalid: {html.escape(', '.join(skipped))}).",
                '#ef6c00',
            )
        else:
            _set_template_status(f"Loaded template '{html.escape(name)}'.", '#2e7d32')

    def _on_template_delete_click(_btn=None):
        name = orca_template_dd.value
        if not name:
            _set_template_status('No template selected.', '#d32f2f')
            return
        _show_save_dialog(False)
        orca_template_delete_prompt.value = (
            f"Delete template '<b>{html.escape(str(name))}</b>'? This cannot be undone."
        )
        _show_delete_dialog(True)

    def _on_template_delete_confirm(_btn=None):
        name = orca_template_dd.value
        _show_delete_dialog(False)
        if not name:
            return
        try:
            delete_orca_template(name)
        except Exception as exc:
            _set_template_status(f'Could not delete template: {html.escape(str(exc))}', '#d32f2f')
            return
        _refresh_template_dd()
        _set_template_status(f"Deleted template '{html.escape(str(name))}'.", '#2e7d32')

    def _on_template_delete_cancel(_btn=None):
        _show_delete_dialog(False)

    orca_template_save_btn.on_click(_on_template_save_click)
    orca_template_save_confirm_btn.on_click(_on_template_save_confirm)
    orca_template_save_cancel_btn.on_click(_on_template_save_cancel)
    orca_template_load_btn.on_click(_on_template_load_click)
    orca_template_delete_btn.on_click(_on_template_delete_click)
    orca_template_delete_confirm_btn.on_click(_on_template_delete_confirm)
    orca_template_delete_cancel_btn.on_click(_on_template_delete_cancel)
    _refresh_template_dd()

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

    orca_save_submit_row = _row([orca_save_btn, orca_submit_btn], wrap=False)
    orca_save_submit_row.layout.margin = '14px 0 0 0'

    orca_left = widgets.VBox([
        _row([orca_job_name], wrap=False),
        _row([orca_coords], wrap=False),
        # The conversions and the 2D editor are the editor's, the same ones
        # the Submit tab has.
        _row([orca_editor.submit_draw_open_btn,
              orca_editor.convert_smiles_button,
              orca_editor.convert_smiles_quick_button,
              orca_editor.convert_smiles_uff_button,
              orca_editor.manta_button]),
        orca_editor.manta_settings_row,
        _row([orca_editor.submit_draw_get_btn,
              orca_editor.submit_draw_update_btn]),
        orca_editor.submit_draw_frame,
        orca_editor.submit_draw_sync,
        _row([orca_copy_coords_btn, orca_check_numbering_btn, orca_apply_numbering_btn]),
        widgets.HTML('<b>Config Templates:</b>'),
        _row([orca_template_dd, orca_template_load_btn, orca_template_save_btn, orca_template_delete_btn]),
        orca_template_save_dialog,
        orca_template_delete_dialog,
        orca_template_status,
        _row([orca_charge, orca_multiplicity]),
        _row([orca_method, orca_job_type]),
        _row([orca_basis, orca_dispersion]),
        _row([orca_ri, orca_aux_basis]),
        _row([orca_autoaux]),
        _row([orca_solvation_type, orca_solvent]),
        _row([orca_print_mos, orca_print_basis]),
        _row([orca_pal, orca_maxcore]),
        _row([orca_slurm_time]),
        widgets.VBox([orca_drop_zone, orca_file_upload, orca_uploaded_files_label],
                     layout=widgets.Layout(width='100%', min_width='0', overflow='hidden', padding='0 8px 0 0')),
        _row([orca_path_files], wrap=False),
        orca_save_submit_row,
        orca_output,
    ], layout=widgets.Layout(
        flex='0 1 48%', max_width='48%', min_width='0', padding='8px', gap='6px',
        box_sizing='border-box', overflow_x='hidden',
    ))

    orca_mol_header = widgets.HTML(
        '<b>Molecule Preview:</b>',
        layout=widgets.Layout(margin='10px 0 0 0'),
    )
    orca_mol_header.add_class('delfin-structure-fs-member')
    orca_mol_header.add_class('delfin-structure-fs-header')
    orca_mol_module = widgets.VBox(
        [orca_mol_header, orca_mol_nav_row, orca_editor.submit_manip_toolbar,
         orca_editor.mol_status, orca_editor.mol_status_fs, orca_mol_output,
         orca_editor.submit_ff_notes],
        layout=widgets.Layout(width='100%', min_width='0', gap='6px'),
    )
    # The editor finds its own controls by this class, and only inside it.
    orca_mol_module.add_class(orca_editor_scope)
    orca_mol_module.add_class('delfin-structure-fs-module')
    orca_mol_module.add_class('orca-structure-fs-module')

    orca_right = widgets.VBox([
        widgets.HTML('<b>ORCA Input Preview (editable):</b>'),
        orca_preview,
        orca_mol_module,
    ], layout=widgets.Layout(
        flex='1 1 0', min_width='0', padding='8px', gap='6px',
        box_sizing='border-box', overflow_x='hidden',
    ))

    split = widgets.HBox(
        [orca_left, orca_right],
        layout=widgets.Layout(width='100%', align_items='stretch', overflow_x='hidden', gap='10px'),
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
        .orca-mol-output .widget-output,
        .orca-mol-output .output_area,
        .orca-mol-output .output_subarea,
        .orca-mol-output .output_wrapper,
        .orca-mol-output .jp-OutputArea,
        .orca-mol-output .jp-OutputArea-child,
        .orca-mol-output .jp-OutputArea-output {
            width: 100% !important;
            height: 100% !important;
            margin: 0 !important;
            padding: 0 !important;
            border: 0 !important;
            overflow: hidden !important;
        }
        .orca-mol-output canvas,
        .orca-mol-output [id^="orca-mol-"] {
            display: block !important;
            margin: 0 !important;
        }
        """
        + structure_viewer_fullscreen_css()
        + """
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
            'target="_blank" rel="noopener noreferrer" '
            'style="color:#1a73e8; text-decoration:underline; cursor:pointer;">'
            'ORCA User Manual</a>'
        ),
        orca_css,
        split,
    ], layout=widgets.Layout(width='100%', padding='8px', box_sizing='border-box'))
    tab_widget.add_class(orca_scope_id)

    # -- Drag-and-drop / click JS for the ORCA drop zone --------------------
    _orca_drop_js = r"""
    (function(){
        if (window._delfinOrcaDropReady) return;
        window._delfinOrcaDropReady = true;

        function injectFiles(uploadBtn, files) {
            if (!uploadBtn || !files || !files.length) return false;
            var dt = new DataTransfer();
            for (var i = 0; i < files.length; i++) {
                if (files[i]) dt.items.add(files[i]);
            }
            if (!dt.files.length) return false;
            var capturedInput = null;
            var origClick = HTMLInputElement.prototype.click;
            HTMLInputElement.prototype.click = function(){
                if (this.type === 'file') {
                    capturedInput = this;
                    return;
                }
                return origClick.call(this);
            };
            try { uploadBtn.click(); } finally {
                HTMLInputElement.prototype.click = origClick;
            }
            if (!capturedInput) return false;
            capturedInput.files = dt.files;
            capturedInput.dispatchEvent(new Event('change', { bubbles: true }));
            return true;
        }

        function install(dropZone) {
            if (!dropZone || dropZone._delfinDropReady) return;
            dropZone._delfinDropReady = true;
            var root = dropZone.closest('.jupyter-widgets') || dropZone.parentElement;
            function findUploadBtn() {
                var p = dropZone.closest('.widget-vbox, .widget-box') || (root && root.parentElement);
                return p ? (p.querySelector('.orca-hidden-upload') || document.querySelector('.orca-hidden-upload')) : null;
            }

            dropZone.addEventListener('click', function(e) {
                e.preventDefault();
                var btn = findUploadBtn();
                if (btn) btn.click();
            });

            dropZone.addEventListener('dragover', function(e) {
                e.preventDefault();
                e.stopPropagation();
                dropZone.style.borderColor = '#1a73e8';
                dropZone.style.background = '#e8f0fe';
                try { e.dataTransfer.dropEffect = 'copy'; } catch(_){}
            });
            dropZone.addEventListener('dragleave', function(e) {
                dropZone.style.borderColor = '#aaa';
                dropZone.style.background = '';
            });
            dropZone.addEventListener('drop', function(e) {
                e.preventDefault();
                e.stopPropagation();
                dropZone.style.borderColor = '#aaa';
                dropZone.style.background = '';
                var files = Array.from(e.dataTransfer.files || []);
                if (!files.length) return;
                var btn = findUploadBtn();
                if (btn) {
                    var ok = injectFiles(btn, files);
                    console.log('[DELFIN] ORCA drop inject:', ok, files.length, 'files');
                }
            });
        }

        function scan(root) {
            if (!root || !root.querySelectorAll) return;
            root.querySelectorAll('.orca-drop-zone').forEach(install);
        }
        scan(document.body);
        new MutationObserver(function(ms) {
            ms.forEach(function(m) {
                m.addedNodes.forEach(function(n) {
                    if (n && n.nodeType === 1) scan(n);
                });
            });
            scan(document.body);
        }).observe(document.body, { childList: true, subtree: true });

        function _hasExternalFiles(e) {
            var dt = e && e.dataTransfer;
            if (!dt) return false;
            if (dt.files && dt.files.length) return true;
            var types = Array.prototype.slice.call(dt.types || []);
            return types.indexOf('Files') >= 0;
        }
        document.addEventListener('dragenter', function(e) {
            if (_hasExternalFiles(e)) { e.preventDefault(); e.stopPropagation(); }
        }, true);
        document.addEventListener('dragover', function(e) {
            if (_hasExternalFiles(e)) {
                e.preventDefault();
                e.stopPropagation();
                try { e.dataTransfer.dropEffect = 'copy'; } catch(_){}
            }
        }, true);
        document.addEventListener('drop', function(e) {
            if (!_hasExternalFiles(e)) return;
            e.preventDefault();
            e.stopPropagation();
            var zone = document.querySelector('.orca-drop-zone');
            if (!zone) return;
            var files = Array.from(e.dataTransfer.files || []);
            if (!files.length) return;
            var btn = document.querySelector('.orca-hidden-upload');
            if (btn) {
                var ok = injectFiles(btn, files);
                console.log('[DELFIN] ORCA document drop inject:', ok, files.length, 'files');
            }
        }, true);
    })();
    """
    # Both startup scripts of this tab, registered together: the
    # structure viewer's fullscreen bootstrap arrived on one branch
    # while the drop handler had already moved to the context on
    # another, and a merge that keeps only one of the two loses a
    # feature without any test noticing which.
    ctx.add_init_js(structure_viewer_fullscreen_bootstrap_js()
                    + '\n' + structure_viewer_fullscreen_kind_js(
                        'orca', 'orca-scope-', ['_orcaBuildViewer'])
                    + '\n' + _orca_drop_js)

    return tab_widget, {
        'orca_pal': orca_pal,
        'orca_maxcore': orca_maxcore,
        'orca_coords': orca_coords,
        'orca_charge': orca_charge,
        'orca_multiplicity': orca_multiplicity,
        'orca_method': orca_method,
        'orca_job_type': orca_job_type,
        'orca_basis': orca_basis,
        'orca_dispersion': orca_dispersion,
        'orca_solvent': orca_solvent,
        'orca_template_dd': orca_template_dd,
        'orca_template_name': orca_template_name,
        'orca_template_save_btn': orca_template_save_btn,
        'orca_template_save_confirm_btn': orca_template_save_confirm_btn,
        'orca_template_load_btn': orca_template_load_btn,
        'orca_template_delete_btn': orca_template_delete_btn,
        'orca_template_delete_confirm_btn': orca_template_delete_confirm_btn,
        'orca_preview': orca_preview,
        'orca_submit_btn': orca_submit_btn,          # destructive: starts real ORCA job
        'orca_copy_coords_btn': orca_copy_coords_btn,
        'orca_check_numbering_btn': orca_check_numbering_btn,
        'orca_apply_numbering_btn': orca_apply_numbering_btn,
        'orca_save_btn': orca_save_btn,
        'orca_mol_prev_btn': orca_mol_prev_btn,
        'orca_mol_next_btn': orca_mol_next_btn,
        'orca_mol_fullscreen_btn': orca_mol_fullscreen_btn,
        'orca_mol_nav_row': orca_mol_nav_row,
        'orca_mol_output': orca_mol_output,
        'update_orca_preview': update_orca_preview,
        # The structure editor this tab holds, under the names it uses.
        **orca_editor.exported,
        'editor_state': state,
        'editor_coords': orca_editor_coords,
        'editor_scope': orca_editor_scope,
        'offer_structures': _take_structures,
        # The editor's own funnel, which is what a conversion goes through.
        'editor_offer_isomers': orca_editor._offer_isomers,
        'read_input': lambda: orca_coords.value,
        'list_structures': lambda: [
            (xyz, len([r for r in strip_xyz_header(xyz).split('\n') if r.strip()]),
             name[:-4] if name.lower().endswith('.xyz') else name)
            for name, xyz in (state.get('xyz_blocks') or [])],
    }
