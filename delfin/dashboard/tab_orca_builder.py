"""ORCA Input Builder tab: build, preview, and submit standalone ORCA jobs."""

import re
import shutil
from pathlib import Path
from collections import Counter
from itertools import permutations, product

import py3Dmol
import ipywidgets as widgets
import json
import uuid
import numpy as np

from IPython.display import clear_output, display, HTML

from delfin.common.control_validator import (
    ORCA_FUNCTIONALS, ORCA_BASIS_SETS, DISP_CORR_VALUES, _RI_JKX_KEYWORDS,
)

from .molecule_viewer import (
    strip_xyz_header,
    DEFAULT_3DMOL_STYLE_JS, DEFAULT_3DMOL_ZOOM,
    patch_viewer_mouse_controls_js,
)
from .input_processing import (
    parse_inp_resources, sanitize_orca_input, clean_input_data,
    smiles_to_xyz_quick_with_previews,
)


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
    )
    orca_convert_smiles_btn = widgets.Button(
        description='CONVERT SMILES', button_style='info',
        layout=widgets.Layout(width='200px'),
    )
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
    orca_additional = widgets.Text(value='', placeholder='e.g. FinalGrid6 NormalPrint',
                                   description='Additional:',
                                   layout=widgets.Layout(width='400px'), style=ws)

    orca_pal = widgets.IntText(value=12, description='PAL (cores):',
                               layout=widgets.Layout(width='250px'), style=ws)
    orca_maxcore = widgets.IntText(value=6000, description='MaxCore (MB):',
                                   layout=widgets.Layout(width='250px'), style=ws)
    orca_slurm_time = widgets.Text(value='12:00:00', placeholder='e.g. 02:00:00',
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
        [orca_mol_prev_btn, orca_mol_nav_label, orca_mol_next_btn],
        layout=widgets.Layout(display='none', align_items='center', gap='6px'),
    )
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
    }

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

    def _orca_kabsch_align(reference_coords, target_coords, mapping=None):
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
        aligned = ref_centered @ rotation + target_centroid
        diff = aligned - target
        rmsd = float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))
        return aligned, rmsd

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
        direct_rmsd = None
        direct_aligned = None
        if ref_seq == target_seq:
            direct_aligned, direct_rmsd = _orca_kabsch_align(ref_coords, target_coords)

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

        identity = np.arange(n_atoms, dtype=int)
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
        xyz_blocks = parse_xyz_blocks(orca_coords.value)
        if xyz_blocks:
            # Use external XYZ file – coordinates are written to the job dir
            coord_block = (
                f'* xyzfile {orca_charge.value} {orca_multiplicity.value}'
                f' {xyz_blocks[0][0]}'
            )
        else:
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
            orca_solvent.value = 'acetonitrile'
            orca_print_mos.value = False
            orca_print_basis.value = False
            orca_additional.value = ''
            orca_pal.value = 12
            orca_maxcore.value = 6000
            orca_slurm_time.value = '12:00:00'
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
        '      tries++;if(tries<80)setTimeout(init,50);return;\n'
        '    }\n'
        '    __RESET__\n'
        '    window._orcaBuildViewState=window._orcaBuildViewState||null;\n'
        '    var prev=window._orcaBuildViewer||null;\n'
        '    if(!__RESETFLAG__&&prev&&typeof prev.getView==="function"){\n'
        '      try{window._orcaBuildViewState=prev.getView();}catch(_e){}\n'
        '    }\n'
        '    var saved=window._orcaBuildViewState;\n'
        '    var viewer=$3Dmol.createViewer(el,{backgroundColor:"white"});\n'
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
        '  }\n'
        '  setTimeout(init,0);\n'
        '})();\n'
        '</script>\n'
    )

    def _atom_labels_js(full_xyz, var='viewer'):
        """Return JS fragment adding atom-index labels at atom centers."""
        lines = full_xyz.split('\n')
        try:
            n_atoms = int(lines[0].strip())
        except (ValueError, IndexError):
            return ''
        calls = []
        for i, line in enumerate(lines[2: 2 + n_atoms]):
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            except ValueError:
                continue
            calls.append(
                f'{var}.addLabel("{i + 1}",'
                f'{{position:{{x:{x:.6f},y:{y:.6f},z:{z:.6f}}},'
                f'fontSize:15,fontColor:"black",showBackground:false,inFront:true}});'
            )
        return '\n    '.join(calls)

    def _viewer_html(xyz_data, label_js='', reset_view=False):
        """Build a self-contained HTML block that renders xyz_data in a $3Dmol viewer.

        The viewer resets on new coordinates, but preserves the same camera view
        while switching between multiple blocks of the same system.
        """
        div_id = 'orca-mol-' + uuid.uuid4().hex[:10]
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
            .replace('__STYLE__', DEFAULT_3DMOL_STYLE_JS)
            .replace('__LABELS__', label_js)
            .replace('__ZOOM__', zoom)
        )
        return html

    def _overlay_viewer_html(reference_xyz, target_xyz, reset_view=False):
        div_id = 'orca-overlay-' + uuid.uuid4().hex[:10]
        mouse_js = patch_viewer_mouse_controls_js('viewer', 'el')
        zoom = str(DEFAULT_3DMOL_ZOOM if DEFAULT_3DMOL_ZOOM is not None else 0.9)
        reset_js = 'window._orcaBuildViewState=null;' if reset_view else ''
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
            '      tries++;if(tries<80)setTimeout(init,50);return;\n'
            '    }\n'
            f'    {reset_js}\n'
            '    window._orcaBuildViewState=window._orcaBuildViewState||null;\n'
            '    var prev=window._orcaBuildViewer||null;\n'
            '    if(!' + ('true' if reset_view else 'false') + '&&prev&&typeof prev.getView==="function"){\n'
            '      try{window._orcaBuildViewState=prev.getView();}catch(_e){}\n'
            '    }\n'
            '    var saved=window._orcaBuildViewState;\n'
            '    var viewer=$3Dmol.createViewer(el,{backgroundColor:"white"});\n'
            f'    {mouse_js}\n'
            f'    viewer.addModel({json.dumps(target_xyz)},"xyz");\n'
            '    viewer.setStyle({model:0},{stick:{radius:0.18,color:"#1f5fff"},sphere:{scale:0.24,color:"#1f5fff"}});\n'
            f'    viewer.addModel({json.dumps(reference_xyz)},"xyz");\n'
            '    viewer.setStyle({model:1},{stick:{radius:0.18,color:"#d32f2f"},sphere:{scale:0.24,color:"#d32f2f"}});\n'
            '    if(saved&&typeof viewer.setView==="function"){\n'
            '      try{viewer.setView(saved);}catch(_e){viewer.zoomTo();viewer.center();viewer.zoom(' + zoom + ');}\n'
            '    }else{\n'
            '      viewer.zoomTo();viewer.center();viewer.zoom(' + zoom + ');\n'
            '    }\n'
            '    viewer.render();\n'
            '    window._orcaBuildViewer=viewer;\n'
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
                _atom_labels_js(reference_xyz, var='viewer'),
                reset_view=reset_view,
            )
        if step == 2:
            return _viewer_html(
                reordered_target_xyz,
                _atom_labels_js(reordered_target_xyz, var='viewer'),
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
            orca_mol_nav_row.layout.display = ''
        else:
            orca_mol_nav_label.value = ''
            orca_mol_nav_row.layout.display = 'none'

    def _refresh_mol_view(reset_view=False):
        """Re-render the molecule viewer, preserving orientation unless *reset_view*."""
        blocks = state['xyz_blocks']
        _update_nav_label()
        _update_numbering_fix_button()
        if state.get('numbering_check_active'):
            orca_mol_output.layout.height = '560px'
            orca_mol_output.layout.min_height = '560px'
        else:
            orca_mol_output.layout.height = '560px'
            orca_mol_output.layout.min_height = '560px'
        with orca_mol_output:
            clear_output(wait=True)
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
                        display(
                            HTML(
                                _numbering_check_view_html(
                                    overlay_result['aligned_reference_xyz'],
                                    target_xyz,
                                    reordered_target_xyz,
                                    state.get('numbering_view_step', 0),
                                    reset_view=reset_view,
                                )
                            )
                        )
                    else:
                        label_js = _atom_labels_js(full_xyz)
                        display(HTML(_viewer_html(full_xyz, label_js, reset_view=reset_view)))
                except Exception as e:
                    print(f'Could not visualize: {e}')
            else:
                raw = orca_coords.value.strip()
                if not raw:
                    print('Paste XYZ coordinates to see 3D preview.')
                    return
                coords = strip_xyz_header(raw)
                if not coords:
                    print('No valid coordinates.')
                    return
                try:
                    atom_lines = [l for l in coords.split('\n') if l.strip()]
                    n = len(atom_lines)
                    xyz_data = f'{n}\nORCA Builder Preview\n{coords}'
                    label_js = _atom_labels_js(xyz_data)
                    display(HTML(_viewer_html(xyz_data, label_js, reset_view=reset_view)))
                except Exception as e:
                    print(f'Could not visualize: {e}')

    def update_orca_molecule_view(change=None):
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
        _refresh_mol_view(reset_view=False)  # keep orientation

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
        _refresh_mol_view(reset_view=False)  # keep orientation

    def handle_orca_convert_smiles(button):
        raw_input = orca_coords.value.strip()
        if not raw_input:
            with orca_output:
                clear_output()
                print('Please enter a SMILES string in the Coordinates box.')
            return
        cleaned_data, input_type = clean_input_data(raw_input)
        if input_type != 'smiles':
            with orca_output:
                clear_output()
                print('Input is not a SMILES string. Please enter a SMILES.')
            return
        with orca_output:
            clear_output()
            print('Converting SMILES to 3D coordinates...')
        xyz_string, num_atoms, _method, preview_items, error = smiles_to_xyz_quick_with_previews(
            cleaned_data
        )
        if error or not xyz_string:
            with orca_output:
                clear_output()
                print(f'Error: {error or "Conversion failed"}')
            return
        xyz_blocks = [('quick.xyz', xyz_string)]
        for _idx, (preview_xyz, _preview_num_atoms, label) in enumerate(preview_items, start=1):
            safe_label = re.sub(r'[^A-Za-z0-9_.-]+', '-', label or f'preview-{_idx}').strip('-')
            if not safe_label:
                safe_label = f'preview-{_idx}'
            xyz_blocks.append((f'{safe_label}.xyz', preview_xyz))
        if len(xyz_blocks) == 1:
            orca_coords.value = xyz_string
        else:
            block_text = []
            for filename, block_xyz in xyz_blocks:
                block_text.append(f'{filename};\n{block_xyz}\n*')
            orca_coords.value = '\n\n'.join(block_text)
        with orca_output:
            clear_output()
            print(f'Converted SMILES to {num_atoms} atoms.')

    orca_convert_smiles_btn.on_click(handle_orca_convert_smiles)

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
                clear_output()
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
            clear_output()
            print('Copied coordinates as XYZ to clipboard.')

    orca_copy_coords_btn.on_click(handle_orca_copy_coordinates)

    def handle_orca_check_numbering(button):
        xyz_blocks = parse_xyz_blocks(orca_coords.value) or []
        if len(xyz_blocks) < 2:
            with orca_output:
                clear_output()
                print('Check Numbering needs at least two named XYZ blocks.')
            return

        ref_name, ref_xyz = xyz_blocks[0]
        try:
            ref_symbols, ref_coords = _orca_parse_xyz_symbols_coords(ref_xyz)
        except Exception as exc:
            with orca_output:
                clear_output()
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
                if result.get('suspicious'):
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
            clear_output()
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
                clear_output()
                print('No checked block is selected for numbering fix.')
            return
        reordered_xyz = str(result.get('reordered_target_xyz') or '').strip()
        if not reordered_xyz:
            with orca_output:
                clear_output()
                print('No numbering fix is available for the selected block.')
            return

        xyz_records[idx]['full_xyz'] = reordered_xyz
        rebuilt = []
        for record in xyz_records:
            label = record.get('raw_name') or record.get('filename') or 'block'
            suffix = record.get('suffix_comment', '')
            header = f'{label};{suffix}' if suffix else f'{label};'
            rebuilt.append(f'{header}\n{record["full_xyz"].strip()}\n*')
        orca_coords.value = '\n\n'.join(rebuilt)
        with orca_output:
            clear_output()
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
            text = re.sub(
                r'(\n\* xyz(?:file)? )',
                f'\n{new_output}\n\\1', text, count=1,
            )

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
              orca_coords, orca_additional, orca_solvation_type, orca_solvent,
              orca_print_mos, orca_print_basis, orca_autoaux]:
        w.observe(update_orca_preview, names='value')

    orca_coords.observe(update_orca_molecule_view, names='value')
    orca_mol_prev_btn.on_click(on_mol_prev)
    orca_mol_next_btn.on_click(on_mol_next)
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

    orca_save_submit_row = _row([orca_save_btn, orca_submit_btn], wrap=False)
    orca_save_submit_row.layout.margin = '14px 0 0 0'

    orca_left = widgets.VBox([
        _row([orca_job_name], wrap=False),
        _row([orca_coords], wrap=False),
        _row([orca_convert_smiles_btn, orca_copy_coords_btn, orca_check_numbering_btn, orca_apply_numbering_btn]),
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
        widgets.VBox([orca_drop_zone, orca_file_upload, orca_uploaded_files_label],
                     layout=widgets.Layout(width='100%', min_width='0', overflow='hidden', padding='0 8px 0 0')),
        _row([orca_path_files], wrap=False),
        orca_save_submit_row,
        orca_output,
    ], layout=widgets.Layout(
        flex='0 1 48%', max_width='48%', min_width='0', padding='8px', gap='6px',
        box_sizing='border-box', overflow_x='hidden',
    ))

    orca_right = widgets.VBox([
        widgets.HTML('<b>ORCA Input Preview (editable):</b>'),
        orca_preview,
        widgets.HTML('<b>Molecule Preview:</b>', layout=widgets.Layout(margin='10px 0 0 0')),
        orca_mol_nav_row,
        orca_mol_output,
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
        'orca_preview': orca_preview,
        'orca_submit_btn': orca_submit_btn,          # destructive: starts real ORCA job
        'orca_convert_smiles_btn': orca_convert_smiles_btn,  # safe: SMILES→XYZ conversion
        'orca_copy_coords_btn': orca_copy_coords_btn,
        'orca_check_numbering_btn': orca_check_numbering_btn,
        'orca_apply_numbering_btn': orca_apply_numbering_btn,
        'orca_save_btn': orca_save_btn,
        'orca_mol_prev_btn': orca_mol_prev_btn,
        'orca_mol_next_btn': orca_mol_next_btn,
        'update_orca_preview': update_orca_preview,
        'init_js': _orca_drop_js,
    }
