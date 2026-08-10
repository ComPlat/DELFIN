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
from . import ketcher as _ketcher
from .molecule_viewer import apply_molecule_view_style, submit_manip_bootstrap_js
from .input_processing import (
    smiles_to_xyz, smiles_to_xyz_quick, smiles_to_xyz_quick_with_previews,
    append_hapto_previews_to_isomers,
    smiles_to_xyz_isomers, is_smiles,
    clean_input_data, parse_resource_settings,
)


_MANTA_GIF_DATA_URI_CACHE = None


def _manta_gif_data_uri():
    """Return the packaged MANTA loading animation as a data URI (cached).

    Mirrors the dashboard logo loader: read the gif shipped under
    ``delfin/logo`` and inline it as a base64 ``data:`` URI so it renders
    inside the molecule viewer without needing a served static path.
    Returns ``''`` when the asset is unavailable.
    """
    global _MANTA_GIF_DATA_URI_CACHE
    if _MANTA_GIF_DATA_URI_CACHE is not None:
        return _MANTA_GIF_DATA_URI_CACHE
    uri = ''
    try:
        gif = importlib.resources.files('delfin').joinpath(
            'logo', 'MANTA_readme_demo.gif')
        if gif.is_file():
            data = gif.read_bytes()
            if data:
                uri = 'data:image/gif;base64,' + base64.b64encode(data).decode('ascii')
    except Exception:
        uri = ''
    _MANTA_GIF_DATA_URI_CACHE = uri
    return uri


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
from delfin.cli_manta import _CHAMPION_FLAGS as _MANTA_CHAMPION_FLAGS
# Top-N to GFN2-optimize + parallel worker count (laptop-bounded; tune here).
_MANTA_OPT_TOPN = 10
_MANTA_OPT_WORKERS = 4


def _manta_best_env(charge, construction="champion", method="gfn2", rank=True):
    """Env for the chosen construction config + GFN2 energy ranking, GFN2 charge
    from the SMILES. construction: 'champion' (full SHIP-31 rich + KAPPA4 reach,
    DEFAULT/maximum richness) | 'builder' (lean core + reach) | 'default' (legacy)."""
    env = {"DELFIN_GFNFF_CHARGE": str(int(charge))}
    if rank:
        env["DELFIN_FFFREE_GFNFF_RANK"] = "1"
        env["DELFIN_CONF_RANK_METHOD"] = method
    if construction != "default":
        env["DELFIN_FFFREE_BUILDER"] = "1"
        env["DELFIN_FRAME_RANK_FIX"] = "1"
        env["DELFIN_CHIRAL_ENUM"] = "1"        # Λ/Δ enantiomer enumeration (>=2 chelate pairs)
        if construction == "champion":
            for _f in _MANTA_CHAMPION_FLAGS:   # de-bloated set (KAPPA4 included; CONF_ENERGY_RANK dropped)
                env["DELFIN_FFFREE_" + _f] = "1"
        else:  # builder = lean core + reach
            env["DELFIN_FFFREE_KAPPA4"] = "1"
            env["DELFIN_FFFREE_SIGMA_ENSEMBLE"] = "1"
            env["DELFIN_FFFREE_CONF_ENERGY_RANK"] = "1"
    return env


def _manta_rank_only(isomers, charge, method="gfn2", spin="auto"):
    """RANK the manifold by xtb SINGLE-POINT energy: reorder best (lowest-energy) first WITHOUT
    changing any geometry.  Each item is ``(xyz_string, num_atoms, label)``; the emitted structures
    stay byte-identical to construction — only their ORDER changes.  spin='auto' -> parity-correct
    uhf per structure (even electrons=singlet, odd=doublet); a fixed multiplicity sets uhf=mult-1.
    Best-effort: any structure whose energy eval fails sinks to the end keeping its geometry.
    Returns the list unchanged if xtb is unavailable or there is nothing to reorder."""
    if not isomers or len(isomers) < 2:
        return isomers
    try:
        from delfin.manta import _gfnff_rank as _gff
    except Exception:
        return isomers
    if not _gff.available():
        return isomers
    import concurrent.futures as _cf

    def _uhf_for(xyz):
        if str(spin) != "auto":
            return max(0, int(spin) - 1)
        try:
            return _gff._n_electrons(xyz, int(charge)) % 2   # parity-correct ground-state multiplicity
        except Exception:
            return 0

    def _energy_one(item):
        xyz = item[0]
        try:
            return _gff.gfnff_energy(xyz, charge=int(charge), uhf=_uhf_for(xyz), method=method)
        except Exception:
            return None
    _max_workers = max(1, min(len(isomers), (os.cpu_count() or 4)))
    try:
        with _cf.ThreadPoolExecutor(max_workers=_max_workers) as ex:
            energies = list(ex.map(_energy_one, isomers))
    except Exception:
        return isomers
    # Ascending by energy; failed evals (None) sink to the end preserving their relative order.
    order = sorted(range(len(isomers)),
                   key=lambda i: (energies[i] is None, energies[i] if energies[i] is not None else 0.0, i))
    return [isomers[i] for i in order]


def _manta_opt_top(isomers, charge, topn=None, method="gfn2", spin="auto"):
    """Geometry-optimize the top-N ranked isomers in parallel (laptop-bounded),
    replace their geometry + label, re-sort the optimized head by opt energy. The
    opt ``method`` FOLLOWS the Rank selection (gfn2/gfnff/gfn1/gfn0) so one switch
    controls both. Each item is ``(xyz_string, num_atoms, label)``. Best-effort:
    any structure whose optimization fails keeps its unrelaxed geometry.
    ``topn`` (user-settable): None -> _MANTA_OPT_TOPN; 0 -> ALL structures (optimise
    the complete ranked manifold, slowest/best); N>0 -> top-N; N<0 -> none."""
    if not isomers:
        return isomers
    if topn is None:
        _n = _MANTA_OPT_TOPN
    elif int(topn) == 0:
        _n = len(isomers)                  # 0 = ALL (optimise everything)
    elif int(topn) < 0:
        return isomers                     # negative = none
    else:
        _n = int(topn)
    import concurrent.futures as _cf
    try:
        from delfin.manta import _gfnff_rank as _gff
    except Exception:
        return isomers
    if not _gff.available():
        return isomers
    head = list(isomers[:_n])
    tail = list(isomers[_n:])

    def _opt_one(item):
        xyz, _na, label = item
        try:
            if str(spin) == "auto":
                # auto-spin: scan multiplicity -> GFN2 ground state (parity-correct)
                r = _gff.gfnff_optimize_autospin(xyz, charge=int(charge), method=method)
            else:
                # fixed multiplicity chosen by the user: uhf = multiplicity - 1
                _uhf = max(0, int(spin) - 1)
                r = _gff.gfnff_optimize(xyz, charge=int(charge), uhf=_uhf, method=method)
        except Exception:
            r = None
        if r and r[0]:
            opt_xyz, e = r
            na = len([ln for ln in opt_xyz.splitlines() if ln.strip()])
            _m = method.upper()
            tag = (" [%s-opt %.1f kcal]" % (_m, e)) if e is not None else " [%s-opt]" % _m
            return ((opt_xyz, na, (label or "isomer") + tag), e)
        return (item, None)

    workers = max(1, min(_MANTA_OPT_WORKERS, len(head)))
    try:
        with _cf.ThreadPoolExecutor(max_workers=workers) as ex:
            opted = list(ex.map(_opt_one, head))
    except Exception:
        return isomers
    # optimized structures first, sorted by GFN2-opt energy (failed/None last)
    opted.sort(key=lambda t: (t[1] is None, t[1] if t[1] is not None else 0.0))
    return [it for (it, _e) in opted] + tail


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

    # Draw the structure instead of typing it.  Ketcher is a browser
    # application of thirty-odd megabytes, so it is fetched when it is first
    # wanted rather than carried in the repository -- see dashboard/ketcher.py.
    submit_draw_open_btn = widgets.ToggleButton(
        value=False, description='DRAW', icon='pencil', button_style='info',
        tooltip=('Draw the structure in Ketcher and hand it back as a SMILES. '
                 'The editor is fetched the first time it is opened.'),
        layout=widgets.Layout(width='110px'),
    )
    submit_draw_get_btn = widgets.Button(
        description='TO SMILES', icon='arrow-down', button_style='success',
        tooltip='Put what is drawn into the input box above, as a SMILES.',
        layout=widgets.Layout(width='130px', display='none'),
    )
    submit_draw_update_btn = widgets.Button(
        description='Update', icon='refresh',
        tooltip='Fetch the newest published Ketcher.',
        layout=widgets.Layout(width='100px', display='none'),
    )
    submit_draw_frame = widgets.HTML(value='', layout=widgets.Layout(
        width='100%', display='none'))
    submit_draw_frame.add_class('submit-ketcher-frame')
    # What the editor hands back, on the same kind of channel the viewer uses:
    # a widget value is ordered and survives a background thread, where a
    # script sent through run_js can be replaced before the page has run it.
    submit_draw_sync = widgets.Textarea(
        value='', layout=widgets.Layout(display='none'))
    submit_draw_sync.add_class('submit-ketcher-sync')

    build_complex_button = widgets.Button(
        description='BUILD COMPLEX', button_style='warning',
        layout=widgets.Layout(width='150px'),
    )

    architector_button = widgets.Button(
        description='ARCHITECTOR', button_style='warning',
        layout=widgets.Layout(width='150px'),
        tooltip='Convert metal-complex SMILES to 3D via architector',
    )

    manta_button = widgets.Button(
        description='MANTA', button_style='',
        layout=widgets.Layout(width='150px'),
        tooltip='MANTA: build the complete coordination-isomer manifold from the '
                'SMILES and rank it by GFN2 energy (best isomer/conformer first; '
                'needs xtb). Results show in the viewer with isomer navigation.',
    )
    # MANTA-logo teal/cyan.
    manta_button.style.button_color = '#1FA9C0'
    manta_button.style.font_weight = 'bold'

    # --- MANTA settings (the 5 keys a user actually needs; MANTA button sits BELOW) ---
    # Power-user knobs (construction / seeds / confs-per-isomer / rank-method /
    # merge-variants) are CLI-only on purpose: each has one sensible value for the
    # dashboard (construction is ALWAYS champion = best; the rest are redundant with
    # Quality), so exposing them only confuses users.  Pinned here.
    _MANTA_DASH_DEFAULTS = dict(construction="champion", num_confs=None,
                               collapse=False)
    # Quality preset -> conformer-seed count.  Selecting a preset auto-fills the
    # Seeds field (transparent: extreme = 60), but the field stays editable so a
    # user can dial a custom seed count on top of the preset's cap/templates.
    _MANTA_PROFILE_SEEDS = {'fast': 12, 'normal': 20, 'max': 40, 'extreme': 60}
    manta_quality = widgets.Dropdown(
        options=['fast', 'normal', 'max', 'extreme'], value='extreme',
        description='Quality:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='160px'),
        tooltip='Conformer depth (the main completeness<->speed switch). extreme guarantees the '
                'GFN2 global minimum is in the manifold (convergence-proven); fast/normal miss it '
                'by ~2.5 kcal/mol on multi-isomer systems.')
    manta_seeds = widgets.IntText(
        value=_MANTA_PROFILE_SEEDS['extreme'], description='Seeds:',
        style={'description_width': 'initial'}, layout=widgets.Layout(width='130px'),
        tooltip='ETKDG conformer seeds. Auto-fills from Quality (extreme=60) for transparency; '
                'edit it to use a custom seed count (more seeds = more conformers, slower).')

    def _sync_seeds_to_quality(change):
        try:
            manta_seeds.value = _MANTA_PROFILE_SEEDS.get(change['new'], manta_seeds.value)
        except Exception:
            pass
    manta_quality.observe(_sync_seeds_to_quality, names='value')
    manta_max_iso = widgets.IntText(
        value=0, description='Max isomers (0=ALL):', style={'description_width': 'initial'},
        layout=widgets.Layout(width='190px'),
        tooltip='Cap emitted isomers. 0 = COMPLETE manifold, never cut off (recommended).')
    # DEFAULT = No / No -> users get ONLY the manifold (no post-processing).  Rank and Opt are
    # SEPARATE opt-ins: Rank reorders by xtb single-point energy (geometry UNCHANGED); Opt xtb-
    # geometry-optimises structures (changes geometry).  You can Rank without Opt, or Opt without Rank.
    manta_rank = widgets.Dropdown(
        options=['No', 'gfn2', 'gfnff', 'gfn1', 'gfn0'], value='No',
        description='Rank:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='150px'),
        tooltip='Energy-RANK the manifold so the global-minimum isomer/conformer is first — '
                'single-point xtb energy, geometry UNCHANGED (pure reordering, no post-processing of '
                'the structures). No = keep construction order (already best-first by realism). '
                'gfn2 = standard; gfnff = fast force-field.')
    manta_opt = widgets.Dropdown(
        options=['No', 'Top 5', 'Top 10', 'Top 20', 'All'], value='No',
        description='Opt:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='140px'),
        tooltip='Geometry-OPTimise structures to a clean final geometry (xtb; method FOLLOWS Rank, or '
                'gfn2 if Rank=No). No = keep the construction geometry (pure manifold, DEFAULT). '
                'Top-N = optimise the N best; All = the whole manifold (slowest/best). '
                'Independent of Rank.')
    manta_spin = widgets.Dropdown(
        options=['auto', '1', '2', '3', '4', '5', '6', '7'], value='auto',
        description='Spin:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='130px'),
        tooltip='Spin multiplicity for the GFN2 energy/rank. auto = scan parity/+2/+4 and take the '
                'GFN2 ground state (parity-correct for open-shell TM). Or FIX it (1=singlet, 2=doublet, '
                '3=triplet, ...) when you know the state. NOTE: GFN2 spin energetics for TM are '
                'approximate -> set it explicitly for accuracy.')
    manta_det = widgets.Dropdown(
        options=['On', 'Off'], value='On',
        description='Determinism:', style={'description_width': 'initial'},
        layout=widgets.Layout(width='160px'),
        tooltip='DETERMINISTIC construction (DEFAULT On) — byte-identical, reproducible manifold, '
                'IDENTICAL to what the CLI and the development loop build (ship = validate). '
                'Off = non-deterministic embedding (racy frame count; only for exploring extra '
                'ETKDG variety). Keep On for a stable, reproducible isomer x conformer manifold.')
    manta_settings_row = widgets.VBox([
        widgets.HTML("<b style='color:#1FA9C0'>MANTA settings</b> "
                     "<span style='color:#888;font-size:90%'>— complete coordination-isomer "
                     "&times; conformer manifold</span>"),
        widgets.HBox([manta_quality, manta_seeds, manta_max_iso, manta_rank, manta_opt, manta_spin, manta_det],
                     layout=widgets.Layout(gap='12px', flex_wrap='wrap', align_items='center')),
    ], layout=widgets.Layout(border='1px solid #d0e7ec', padding='8px', margin='4px 0'))

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

    mol_status = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', margin='0 0 6px 0'),
    )
    # A second one, for the overlay.  Not the same widget moved: fullscreen
    # relocates its members by hand and ipywidgets knows nothing about it, so
    # the line borrowed for the big view came back somewhere else and was lost
    # from the small one.  Two widgets, one text, nothing moved.
    mol_status_fs = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', margin='0 0 6px 0'),
    )
    mol_status_fs.add_class('submit-fs-member-status')
    # It lives in the ordinary layout so fullscreen can pick it up from there,
    # but it must not be seen next to the line it is a copy of -- the message
    # would simply be printed twice.  Hidden here, shown in the overlay.
    mol_status_fs.layout.display = 'none'
    mol_output = widgets.Output(layout=widgets.Layout(
        border='2px solid #1976d2', width='100%', height=f'{SUBMIT_MOL_HEIGHT}px',
        overflow='hidden', box_sizing='border-box',
    ))
    mol_output.add_class('submit-mol-output')

    submit_scope_id = f'submit-scope-{abs(id(coords_widget))}'

    submit_fullscreen_btn = widgets.Button(
        description='', icon='expand',
        tooltip='Toggle fullscreen (Esc to exit)',
        layout=widgets.Layout(width='40px', height='30px'),
        disabled=True,
    )
    submit_fullscreen_btn.add_class('submit-fullscreen-btn')

    submit_select_btn = widgets.ToggleButton(
        value=False, description='Select', icon='crosshairs',
        button_style='',
        tooltip=(
            'Click atoms to pick or unpick them. Hold Shift and drag for a '
            'rectangle; add Ctrl to keep the previous selection.'
        ),
        layout=widgets.Layout(width='96px', height='30px'),
        disabled=True,
    )
    submit_manip_btn = widgets.ToggleButton(
        value=False, description='Manipulate', icon='arrows',
        button_style='',
        tooltip=(
            'Grab any atom and drag it; grabbing a selected atom moves the '
            'whole selection. Drag empty space to turn the view. Right-click '
            'an atom to set the pivot, right-drag to rotate the selection '
            'about it.'
        ),
        layout=widgets.Layout(width='112px', height='30px'),
        disabled=True,
    )
    from .molecule_builder import DRAW_ELEMENTS

    submit_draw_btn = widgets.ToggleButton(
        value=False, description='Draw', icon='pencil',
        button_style='',
        tooltip=(
            'Click empty space to place an atom of the chosen element, drag '
            'from an atom to grow a new one bonded to it, drag onto another '
            'atom to bond them at the chosen order, tap an atom to change its '
            'element, and press Delete to remove the selection. Free valences '
            'are filled with hydrogen. The right button still turns the view.'
        ),
        layout=widgets.Layout(width='88px', height='30px'),
        disabled=True,
    )
    submit_element_dd = widgets.Dropdown(
        options=list(DRAW_ELEMENTS), value='C',
        layout=widgets.Layout(width='78px', display='none'),
        disabled=True,
    )
    submit_manip_clear_btn = widgets.Button(
        description='Clear', button_style='warning',
        tooltip='Clear selection & pivot',
        layout=widgets.Layout(width='66px', height='30px'),
        disabled=True,
    )
    submit_manip_undo_btn = widgets.Button(
        description='Undo', button_style='info', icon='undo',
        tooltip='Undo last move/rotate (Ctrl-Z)',
        layout=widgets.Layout(width='84px', height='30px'),
        disabled=True,
    )
    submit_relax_btn = widgets.ToggleButton(
        value=False, description='Dynamik Opt', icon='magic',
        button_style='',
        tooltip=(
            'Keep the force field running: the structure settles continuously, '
            'and an atom you grab drags the relaxing molecule along with it. '
            'Use the strength control to make it gentler. Toggle off to stop.'
        ),
        layout=widgets.Layout(width='132px', height='30px'),
        disabled=True,
    )
    submit_optimize_btn = widgets.ToggleButton(
        value=False,
        description='Optimize', button_style='success', icon='compress',
        tooltip=(
            'Minimise the frame on screen. A switch: press again to stop, and '
            'what is stopped is what you were looking at. Undo restores the '
            'geometry from before the run in one step.'
        ),
        layout=widgets.Layout(width='112px', height='30px'),
        disabled=True,
    )
    submit_optimize_all_btn = widgets.ToggleButton(
        value=False,
        description='all', button_style='success',
        tooltip=(
            'Minimise every frame that is loaded -- all isomers or batch '
            'entries, not only the one on screen. They run one after another, '
            'never side by side: a login node is shared, and a set of isomers '
            'would otherwise start one xtb per frame at once.'
        ),
        layout=widgets.Layout(width='52px', height='30px'),
        disabled=True,
    )
    submit_ff_dd = widgets.Dropdown(
        options=[('UFF', 'uff'), ('MMFF94', 'mmff94'),
                 ('GFN-FF', 'gfnff'), ('GFN2-xTB', 'gfn2'), ('g-xTB', 'gxtb')],
        value='uff',
        tooltip=(
            'What Optimise minimises with. UFF and MMFF94 run in the browser '
            'and also drive the live relaxation while you drag. GFN-FF, '
            'GFN2-xTB and g-xTB run xtb on the server, and they know about the '
            'metal where UFF guesses. g-xTB approximates wB97M-V/def2-TZVPPD '
            'and needs a build of its own; Install g-xTB fetches it.'
        ),
        layout=widgets.Layout(width='128px'),
        disabled=True,
    )
    # xtb needs both, and there is no honest default for a metal complex: the
    # wrong number of unpaired electrons gives a confident wrong answer rather
    # than an error.  Shown only when a GFN method is chosen.
    submit_gfn_charge = widgets.IntText(
        value=0, description='q', step=1,
        style={'description_width': '12px'},
        layout=widgets.Layout(width='72px', display='none'),
    )
    submit_gfn_mult = widgets.IntText(
        value=1, description='M', step=1,
        style={'description_width': '14px'},
        layout=widgets.Layout(width='72px', display='none'),
    )
    # The lines between the atoms are 3Dmol's own bond list, worked out once
    # when the model was built.  Moving atoms does not touch it, so a bond
    # pulled apart goes on being drawn -- the picture keeps the connectivity
    # the structure had rather than the one it has.  Off by default: it costs a
    # rebuild per frame, and in a crowded coordination sphere the perception is
    # at its limit and the lines flicker.
    submit_dyn_bonds_btn = widgets.ToggleButton(
        value=False, description='Dyn. bonds', icon='link',
        tooltip=(
            'Let the lines between the atoms follow the distances while the '
            'structure moves, instead of keeping the bonds it was drawn with. '
            'Only the picture changes -- what the calculation holds together '
            'is a separate question, and Bond and Unbond are how that is said.'
        ),
        layout=widgets.Layout(width='104px', height='30px'),
    )
    submit_gfn_solvent = widgets.Dropdown(
        options=[(label, name) for name, label in _gfn.SOLVENTS.items()],
        value='',
        tooltip=(
            'Optimise with an implicit solvent around the structure, using '
            'ALPB -- xtb\'s own recommendation, and the model that covers '
            'every solvent it is parametrised for. A geometry optimised in '
            'the gas phase and one optimised in water are different answers; '
            'the status line says which you got.'
        ),
        layout=widgets.Layout(width='140px', display='none'),
    )
    submit_gfn_autospin = widgets.Checkbox(
        value=False, description='auto M', indent=False,
        tooltip=(
            'Try the multiplicities the electron count allows and keep the '
            'one that comes out lowest. For an open-shell metal a fixed guess '
            'gives a confidently wrong energy and, through it, a wrong '
            'geometry -- but it costs three runs instead of one.'
        ),
        layout=widgets.Layout(width='86px', display='none'),
    )
    # Shown only when a GFN method is chosen and xtb cannot be found.  Two
    # presses, not one: the install fetches a conda environment of a few
    # hundred megabytes, and on a cluster the right answer is often "load the
    # module instead" -- so what it would run is shown before it runs.
    submit_xtb_install_btn = widgets.Button(
        description='Install xtb', icon='download',
        tooltip='Fetch xtb with DELFIN\'s own installer. It will say what it '
                'runs before it runs it.',
        layout=widgets.Layout(width='118px', height='30px', display='none'),
    )
    submit_xtb_confirm_btn = widgets.Button(
        description='Yes, install', button_style='warning',
        layout=widgets.Layout(width='108px', height='30px', display='none'),
    )
    submit_xtb_cancel_btn = widgets.Button(
        description='Cancel',
        layout=widgets.Layout(width='72px', height='30px', display='none'),
    )
    submit_internal_label = widgets.HTML(
        value=(
            '<span class="submit-internal-label" '
            'style="color:#888;font-size:0.9em;white-space:nowrap;">'
            'pick 2-4 atoms</span>'
        ),
        layout=widgets.Layout(margin='0 0 0 4px'),
    )
    submit_internal_value = widgets.FloatText(
        value=0.0, step=0.01,
        layout=widgets.Layout(width='92px', height='30px'),
        disabled=True,
    )
    submit_internal_value.add_class('submit-internal-value')
    submit_internal_btn = widgets.ToggleButton(
        value=False, description='Set', button_style='primary',
        tooltip=(
            'Turn the value by hand and watch it: while Set is on, the box '
            'drives the selection live -- the arrow keys step a bond by '
            '0.01 A, an angle by 0.1 and a dihedral by 0.5 degrees, and the '
            'fragment on the far side of the coordinate follows. Two atoms are a bond '
            'length, three an angle, four a dihedral. Hold is the other '
            'question: it keeps a value at its target while the field runs, '
            'with pull or fix.'
        ),
        layout=widgets.Layout(width='58px', height='30px'),
        disabled=True,
    )
    submit_manip_status = widgets.HTML(
        value='<span class="submit-manip-status" style="color:#888;font-size:0.9em;">— viewer empty —</span>',
        # Takes a share of the row when there is room -- fullscreen keeps the
        # toolbar on one line -- and wraps to its own line when there is not,
        # which is what happens inside the tab on a laptop.
        layout=widgets.Layout(
            flex='1 1 260px', min_width='0', overflow_x='hidden',
        ),
    )
    submit_manip_sync = widgets.Textarea(value='', layout=widgets.Layout(display='none'))
    submit_manip_sync.add_class('submit-manip-sync')
    # Coordinates from the kernel, one frame at a time.  Not through run_js:
    # that writes into a single Output and clears it first, so twenty scripts a
    # second overwrite each other before the page has rendered them -- the
    # relaxation ran, the last structure appeared, and nothing in between did.
    # A widget value is ordered, survives a background thread, and cannot be
    # clobbered by the next one.
    submit_gfn_frame = widgets.Textarea(value='', layout=widgets.Layout(display='none'))
    submit_gfn_frame.add_class('submit-gfn-frame')

    submit_strength_slider = widgets.IntSlider(
        value=20, min=1, max=200, step=1,
        description='Strength', continuous_update=False,
        readout=True, readout_format='d',
        style={'description_width': '58px'},
        layout=widgets.Layout(width='200px'),
        disabled=True,
    )
    submit_pick_sync = widgets.Text(value='', layout=widgets.Layout(display='none'))
    submit_pick_sync.add_class('submit-pick-sync')
    # Keyboard shortcuts for things Python owns. Unbond is not a picture edit:
    # it changes the topology the force field is built from, so the browser
    # cannot carry it out alone and has to ask through here.
    submit_cmd_sync = widgets.Text(value='', layout=widgets.Layout(display='none'))
    submit_cmd_sync.add_class('submit-cmd-sync')
    submit_poly_dd = widgets.Dropdown(
        options=[('— polyhedron —', '')], value='',
        layout=widgets.Layout(width='190px', display='none'),
        disabled=True,
    )
    submit_poly_turn_btn = widgets.Button(
        description='Turn', icon='refresh', button_style='',
        tooltip=(
            'Step to the next distinct arrangement on this polyhedron: which '
            'ligands take the axial or apical positions. Only shown where the '
            'vertices are not all alike -- an octahedron has nothing to turn, '
            'a trigonal bipyramid has ten arrangements.'
        ),
        layout=widgets.Layout(width='78px', height='30px', display='none'),
        disabled=True,
    )
    submit_hyb_dd = widgets.Dropdown(
        options=[('— hybridisation —', '')], value='',
        layout=widgets.Layout(width='190px', display='none'),
        disabled=True,
    )
    submit_hyb_auto_btn = widgets.Button(
        description='Types', icon='magic', button_style='',
        tooltip=(
            'Read each carbon\'s hybridisation off how many partners it is '
            'bonded to: four is tetrahedral, three trigonal planar, two '
            'linear. Stronger than perception, which goes through bond '
            'orders and misses a double bond it cannot see in the geometry. '
            'Applies to the selection, or to every carbon when nothing is '
            'selected.'
        ),
        layout=widgets.Layout(width='84px', height='30px'),
        disabled=True,
    )
    submit_settle_btn = widgets.ToggleButton(
        value=True, description='Settle', icon='level-down',
        button_style='info',
        tooltip=(
            'When you let go of an atom, let the structure relax around its '
            'new position instead of keeping the strain of the drag. Switch '
            'off to leave atoms exactly where you put them.'
        ),
        layout=widgets.Layout(width='92px', height='30px'),
        disabled=True,
    )
    submit_swap_btn = widgets.Button(
        description='Swap', button_style='', icon='exchange',
        tooltip=(
            'Exchange the two selected ligands on the polyhedron: they are '
            'pulled onto each other\'s vertex instead of back to their own.'
        ),
        layout=widgets.Layout(width='78px', height='30px', display='none'),
        disabled=True,
    )
    submit_bond_btn = widgets.Button(
        description='Bond', icon='link', button_style='',
        tooltip=(
            'Draw a bond between the two selected atoms. Distance-based '
            'perception is unreliable in a crowded coordination sphere, and '
            'the coordination number and the force field both follow from '
            'these bonds.'
        ),
        layout=widgets.Layout(width='74px', height='30px'),
        disabled=True,
    )
    submit_unbond_btn = widgets.Button(
        description='Unbond', icon='unlink', button_style='',
        tooltip='Remove the bond between the two selected atoms.',
        layout=widgets.Layout(width='90px', height='30px'),
        disabled=True,
    )
    submit_hold_mode = widgets.Dropdown(
        options=[('pull', 'pull'), ('fix', 'fix')], value='pull',
        layout=widgets.Layout(width='78px'),
        disabled=True,
    )
    submit_hold_btn = widgets.Button(
        description='Hold', button_style='warning', icon='thumb-tack',
        tooltip=(
            'Hold the value the selection describes while the field runs, '
            'instead of only setting it once. Held values appear in the list '
            'beside this button and can be dropped again there.'
        ),
        layout=widgets.Layout(width='72px', height='30px'),
        disabled=True,
    )
    submit_constraint_dd = widgets.Dropdown(
        options=[('no constraints', '')], value='',
        layout=widgets.Layout(width='210px', display='none'),
        disabled=True,
    )
    submit_reset_btn = widgets.Button(
        description='Reset', icon='undo', button_style='danger',
        tooltip=(
            'Back to the structure as it was loaded, and drop everything set '
            'on it since: held values, polyhedron, hand-made bonds, '
            'hybridisation overrides and the edit history.'
        ),
        layout=widgets.Layout(width='84px', height='30px'),
        disabled=True,
    )
    submit_constraint_del = widgets.Button(
        description='', icon='times', button_style='danger',
        tooltip='Drop the selected constraint',
        layout=widgets.Layout(width='40px', height='30px', display='none'),
        disabled=True,
    )
    submit_internal_group = widgets.HBox(
        [submit_internal_label, submit_internal_value,
         submit_internal_btn, submit_hold_btn, submit_hold_mode],
        layout=widgets.Layout(
            gap='6px', align_items='center', flex_flow='row nowrap',
            flex='0 0 auto', overflow='visible',
        ),
    )

    submit_manip_toolbar = widgets.HBox(
        [
            submit_fullscreen_btn,
            submit_select_btn, submit_manip_btn, submit_draw_btn,
            submit_element_dd,
            submit_manip_clear_btn, submit_manip_undo_btn, submit_reset_btn,
            submit_ff_dd, submit_gfn_charge, submit_gfn_mult,
            submit_gfn_autospin, submit_gfn_solvent,
            submit_xtb_install_btn, submit_xtb_confirm_btn,
            submit_xtb_cancel_btn,
            submit_strength_slider,
            submit_optimize_btn, submit_optimize_all_btn,
            submit_relax_btn, submit_settle_btn,
            submit_poly_dd, submit_poly_turn_btn,
            submit_hyb_dd, submit_hyb_auto_btn,
            submit_internal_group,
            submit_bond_btn, submit_unbond_btn, submit_dyn_bonds_btn,
            submit_swap_btn, submit_constraint_dd, submit_constraint_del,
            submit_pick_sync, submit_cmd_sync,
            submit_manip_status, submit_manip_sync, submit_gfn_frame,
        ],
        layout=widgets.Layout(
            display='none', gap='6px', align_items='center',
            width='100%', flex_flow='row wrap',
            margin='0 0 6px 0', overflow='visible',
        ),
    )

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

    def _set_smiles_conversion_busy(is_busy):
        state['smiles_busy'] = bool(is_busy)
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
        ctx.set_busy(state['smiles_busy'] or state['batch_preview_busy'])

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

    def _set_mol_status(*lines, spinner=False):
        # Both copies always say the same thing; which one is on screen is the
        # overlay's business, not the caller's.
        rendered = [html.escape(str(line)) for line in lines if line not in (None, '')]
        spinner_html = (
            "<span class='delfin-busy' style='margin-right:6px; vertical-align:middle;' "
            "title='Working'></span>"
            if spinner else ''
        )
        text_html = '<br>'.join(rendered)
        if not spinner_html and not text_html:
            mol_status.value = ''
            mol_status_fs.value = ''
            return
        if spinner_html and text_html:
            first, *rest = rendered
            body = spinner_html + first
            if rest:
                body += '<br>' + '<br>'.join(rest)
        else:
            body = spinner_html + text_html
        rendered_html = (
            "<div style='font-family: monospace; white-space: pre-wrap; "
            "font-size: 13px; line-height: 1.35;'>"
            f"{body}</div>"
        )
        mol_status.value = rendered_html
        # The fullscreen copy is there to report work -- a spinner, a
        # trajectory, a result.  Asking for coordinates is not work, and in
        # fullscreen there is a structure on screen to look at, so the prompt
        # would be a permanent line saying nothing.
        prompt = any('enter XYZ' in str(line) or 'Load a structure' in str(line)
                     for line in lines)
        mol_status_fs.value = '' if prompt else rendered_html

    def _clear_mol_status():
        mol_status.value = ''

    def _build_mol_output_bundle(xyz_data):
        view = py3Dmol.view(width='100%', height=SUBMIT_MOL_HEIGHT)
        view.addModel(xyz_data, 'xyz')
        apply_molecule_view_style(view)
        scope_key_js = json.dumps(submit_scope_id)
        registration = (
            '\n'
            + submit_manip_bootstrap_js()
            + '\n(function(){\n'
            '  try {\n'
            '    window._submitMolViewerByScope = window._submitMolViewerByScope || {};\n'
            # Every render of this tab creates a fresh WebGL viewer. Without
            # releasing the previous one its context, observers and window-level
            # mouse listeners stay alive; browsers cap contexts and start
            # killing the oldest, which blacks out the viewers in other tabs.
            f'    var prev = window._submitMolViewerByScope[{scope_key_js}];\n'
            # Keep the camera across a re-render. Optimising or stepping to
            # another isomer rebuilds the viewer, and the view otherwise snapped
            # back to the default orientation every time.
            '    window._submitMolViewByScope = window._submitMolViewByScope || {};\n'
            '    if (prev && prev !== viewer_UNIQUEID) {\n'
            '      try {\n'
            '        var prevModel = prev.getModel ? prev.getModel() : null;\n'
            '        var prevCount = prevModel && prevModel.selectedAtoms\n'
            '          ? prevModel.selectedAtoms({}).length : -1;\n'
            f'        window._submitMolViewByScope[{scope_key_js}] =\n'
            '          {view: prev.getView(), atoms: prevCount};\n'
            '      } catch(e) {}\n'
            '      if (window.__delfinDisposeViewer) window.__delfinDisposeViewer(prev);\n'
            '    }\n'
            '    try {\n'
            f'      var saved = window._submitMolViewByScope[{scope_key_js}];\n'
            '      var model = viewer_UNIQUEID.getModel();\n'
            '      var count = model && model.selectedAtoms\n'
            '        ? model.selectedAtoms({}).length : -1;\n'
            # Same structure, just moved -- keep the view. A different one
            # deserves a fresh look at it. An edit that adds or removes an
            # atom is neither: the count changes but it is still the molecule
            # being worked on, so the camera stays put rather than snapping
            # back to the default every time something is drawn.
            '      var edited = !!window.__delfinStructureEdit;\n'
            '      window.__delfinStructureEdit = false;\n'
            '      if (saved && saved.view && (edited || saved.atoms === count)) {\n'
            '        viewer_UNIQUEID.setView(saved.view);\n'
            '        viewer_UNIQUEID.__delfinUserInteracted = true;\n'
            '        viewer_UNIQUEID.render();\n'
            '      }\n'
            '    } catch(e) {}\n'
            f'    window._submitMolViewerByScope[{scope_key_js}] = viewer_UNIQUEID;\n'
            '    var el = document.getElementById("3dmolviewer_UNIQUEID");\n'
            '    var fire = function(){\n'
            '      if (window.__delfinSubmitManip) {\n'
            f'        window.__delfinSubmitManip.onViewerReady({scope_key_js}, el);\n'
            '      }\n'
            '    };\n'
            '    setTimeout(fire, 80);\n'
            '  } catch(e) {}\n'
            '})();\n'
        )
        if hasattr(view, 'startjs'):
            view.startjs += registration
        html_payload = view._make_html()
        return ({
            'output_type': 'display_data',
            'data': {
                'application/3dmoljs_load.v0': html_payload,
                'text/html': html_payload,
            },
            'metadata': {},
        },)

    def _clear_mol_output():
        mol_output.outputs = ()

    def _replace_mol_output_text(*lines):
        _set_mol_status(*lines)
        _clear_mol_output()
        _set_manip_toolbar_enabled(False)

    def _ensure_manip_bootstrap():
        if state['manip_bootstrap_done']:
            return
        try:
            ctx.run_js(submit_manip_bootstrap_js())
            state['manip_bootstrap_done'] = True
        except Exception:
            pass

    def _set_manip_toolbar_enabled(enabled):
        submit_fullscreen_btn.disabled = not enabled
        submit_select_btn.disabled = not enabled
        submit_manip_btn.disabled = not enabled
        submit_draw_btn.disabled = not enabled
        submit_element_dd.disabled = not enabled
        submit_manip_clear_btn.disabled = not enabled
        submit_relax_btn.disabled = not enabled
        submit_strength_slider.disabled = not enabled
        submit_settle_btn.disabled = not enabled
        submit_bond_btn.disabled = not enabled
        submit_unbond_btn.disabled = not enabled
        submit_hyb_auto_btn.disabled = not enabled
        submit_ff_dd.disabled = not enabled
        submit_optimize_btn.disabled = not enabled
        submit_optimize_all_btn.disabled = not enabled
        submit_internal_value.disabled = not enabled
        submit_internal_btn.disabled = not enabled
        submit_hold_btn.disabled = not enabled
        submit_hold_mode.disabled = not enabled
        submit_manip_undo_btn.disabled = not enabled
        submit_reset_btn.disabled = not enabled
        submit_manip_toolbar.layout.display = 'flex' if enabled else 'none'
        if not enabled:
            if submit_select_btn.value:
                submit_select_btn.value = False
            if submit_manip_btn.value:
                submit_manip_btn.value = False

    def _replace_mol_output_view(xyz_data):
        _clear_mol_status()
        _ensure_manip_bootstrap()
        mol_output.outputs = _build_mol_output_bundle(xyz_data)
        _set_manip_toolbar_enabled(True)

    def _show_mol_busy(message):
        """Render the animated MANTA loader centered in the viewer.

        While the manifold build / GFN2 ranking / optimization runs there is
        no structure to show, so fill the 3D viewer box with the floating
        MANTA mark plus a status caption. ``_replace_mol_output_view`` later
        swaps it out for the finished structure (same as before).
        """
        _clear_mol_status()
        safe_msg = html.escape(str(message))
        gif_uri = _manta_gif_data_uri()
        if gif_uri:
            visual_html = (
                f"<img src='{gif_uri}' alt='MANTA working' "
                "style='width:92%; max-width:1050px; height:auto; "
                "object-fit:contain; image-rendering:auto;'/>"
            )
        else:
            visual_html = (
                "<span class='delfin-busy' style='width:42px; height:42px; "
                "border-width:4px;'></span>"
            )
        payload = (
            "<div class='manta-busy-stage' style='display:flex; "
            "flex-direction:column; align-items:center; justify-content:center; "
            f"width:100%; min-height:{SUBMIT_MOL_HEIGHT - 6}px; gap:20px; "
            "background:#fcfcfc;'>"
            f"{visual_html}"
            "<div class='manta-busy-caption' style='display:flex; "
            "align-items:center; justify-content:center; gap:8px; max-width:84%; "
            "font-family:monospace; font-size:13px; line-height:1.4; "
            "color:#1976d2; text-align:center;'>"
            "<span class='delfin-busy' style='vertical-align:middle;'></span>"
            f"<span>{safe_msg}</span></div></div>"
        )
        mol_output.outputs = ({
            'output_type': 'display_data',
            'data': {'text/html': payload},
            'metadata': {},
        },)
        _set_manip_toolbar_enabled(False)

    def _apply_smiles_conversion_result(task_id, *, quick, cleaned_data, result):
        if task_id != state['smiles_task_id']:
            return

        _set_smiles_conversion_busy(False)
        if quick:
            xyz_string = result.get('xyz_string')
            preview_items = result.get('preview_items') or []
            error = result.get('error')
            if error or not xyz_string:
                _replace_mol_output_text(f'Error: {error or "Conversion failed"}')
                return
            state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': xyz_string}
            state['isomers'] = [(xyz_string, result['num_atoms'], 'quick')] + preview_items
            _show_isomer_at_index(0)
            return

        isomers = result.get('isomers') or []
        error = result.get('error')
        if error or not isomers:
            _replace_mol_output_text(
                f'SMILES: {cleaned_data}',
                f'Error: {error or "No isomers generated"}',
            )
            state['converted_xyz_cache'] = {'smiles': None, 'xyz': None}
            state['isomers'] = []
            isomer_nav_row.layout.display = 'none'
            return

        state['converted_xyz_cache'] = {'smiles': cleaned_data, 'xyz': isomers[0][0]}
        state['isomers'] = isomers
        _show_isomer_at_index(0)

    def _apply_batch_preview_result(task_id, entry, preview_payload):
        if task_id != state['batch_preview_task_id']:
            return
        _set_batch_preview_busy(False)
        _render_batch_preview(entry, preview_payload)

    def update_molecule_view(change=None):
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
            state['constraints'] = []
            state['bond_edits'] = {}
            state['hyb_overrides'] = {}
            state['structure_undo'] = []
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

        xyz_data = f'{num_atoms}\nIsomer: {label}\n{xyz_string}'
        _replace_mol_output_view(xyz_data)

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

    def _start_smiles_conversion(*, apply_uff: bool, quick: bool, rank: bool = False,
                                 quality_mode=None, seeds_override=None,
                                 max_isomers=None, opt_topn=None, construction=None,
                                 method="gfn2", num_confs=None, collapse=None, spin="auto",
                                 deterministic: bool = True):
        cached_smiles = state['converted_xyz_cache'].get('smiles') if quick else None
        raw_input = (cached_smiles or coords_widget.value).strip()
        if not raw_input:
            _replace_mol_output_text('Please enter SMILES in the input box.')
            return

        cleaned_data, input_type = clean_input_data(raw_input)
        if input_type != 'smiles':
            _replace_mol_output_text('Please enter SMILES in the input box.')
            return

        state['smiles_task_id'] += 1
        task_id = state['smiles_task_id']
        _set_smiles_conversion_busy(True)

        # The MANTA logo animation must ALWAYS play while MANTA builds the manifold — WITH or WITHOUT
        # post-processing.  A MANTA submit is (not quick) AND has a construction preset (the convert
        # buttons pass construction=None); Rank=No / Opt<0 = manifold-only, but the logo still shows.
        if (not quick) and (construction is not None):
            _method_lbl = str(method).upper()
            _topn = -1 if opt_topn is None else int(opt_topn)
            _opt_on = _topn >= 0
            if not rank and not _opt_on:
                _caption = ('MANTA: building the complete coordination-isomer × conformer manifold '
                            '(no post-processing)...')
            else:
                _parts = []
                if rank:
                    _parts.append(f'{_method_lbl} energy-ranking')
                if _opt_on:
                    _parts.append('optimizing all' if _topn == 0 else f'optimizing top {_topn}')
                _caption = ('MANTA: building manifold + ' + ' + '.join(_parts) +
                            ' (needs xtb; takes a bit)...')
            _show_mol_busy(_caption)
        else:
            _clear_mol_output()
            if quick:
                _set_mol_status('Quick convert (single structure)...', spinner=True)
            elif apply_uff:
                _set_mol_status('Converting SMILES with UFF...', spinner=True)
            else:
                _set_mol_status('Converting SMILES (no UFF)...', spinner=True)

        def _worker():
            import os
            # MANTA "best version": derive the GFN2 charge from the SMILES, then
            # apply the SHIP-31 champion construction + GFN2-rank env for this
            # build only (snapshot + restore so the global env isn't polluted).
            _chg = 0
            if rank or construction:
                try:
                    from rdkit import Chem as _Chem
                    _m = _Chem.MolFromSmiles(cleaned_data, sanitize=False)
                    if _m is not None:
                        _chg = _Chem.GetFormalCharge(_m)
                except Exception:
                    _chg = 0
            # construction env applies ONLY for MANTA (construction set); the plain
            # convert/build-complex buttons pass construction=None -> unchanged behaviour.
            _best_env = (_manta_best_env(_chg, construction=construction, method=method,
                                         rank=rank) if construction else {})
            _saved_env = {k: os.environ.get(k) for k in _best_env}
            os.environ.update(_best_env)
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
                    # Interactive metal-complex conversion should prioritize
                    # isomer diversity over strict reproducibility.
                    _iso_kwargs = dict(
                        apply_uff=apply_uff,
                        collapse_label_variants=(bool(collapse) if collapse is not None else False),
                        include_binding_mode_isomers=True,
                        deterministic=deterministic,
                    )
                    # user-exposed completeness/speed switches (MANTA settings row);
                    # None -> library default. max_isomers None/0 -> COMPLETE (no cut).
                    if quality_mode:
                        _iso_kwargs["quality_mode"] = quality_mode
                    if seeds_override:
                        _iso_kwargs["seeds_override"] = int(seeds_override)
                    if max_isomers:
                        _iso_kwargs["max_isomers"] = int(max_isomers)
                    if num_confs:
                        _iso_kwargs["num_confs"] = int(num_confs)
                    isomers, error = smiles_to_xyz_isomers(cleaned_data, **_iso_kwargs)
                    if not error and isomers:
                        isomers = append_hapto_previews_to_isomers(
                            isomers,
                            cleaned_data,
                            include_quick=apply_uff,
                        )
                    if rank and not error and isomers:
                        # RANK (opt-in): reorder best-first by xtb SINGLE-POINT energy.
                        # Geometry UNCHANGED — the emitted structures stay byte-identical to
                        # construction, only their order changes.
                        isomers = _manta_rank_only(isomers, _chg, method=method, spin=spin)
                    if (opt_topn is not None and int(opt_topn) >= 0
                            and not error and isomers):
                        # OPT (opt-in, independent of Rank): xtb geometry-optimise the top-N
                        # (0 = the whole manifold) for best-possible final geometry.
                        isomers = _manta_opt_top(isomers, _chg, topn=opt_topn, method=method, spin=spin)
                    result = {'error': error, 'isomers': isomers}
            except Exception as exc:
                result = {'error': str(exc)}
            finally:
                for _k, _v in _saved_env.items():
                    if _v is None:
                        os.environ.pop(_k, None)
                    else:
                        os.environ[_k] = _v

            _schedule_ui_update(
                _apply_smiles_conversion_result,
                task_id,
                quick=quick,
                cleaned_data=cleaned_data,
                result=result,
            )

        threading.Thread(target=_worker, daemon=True).start()

    def _convert_smiles(*, apply_uff: bool):
        _start_smiles_conversion(apply_uff=apply_uff, quick=False)

    def handle_convert_smiles(button):
        _convert_smiles(apply_uff=False)

    def handle_convert_smiles_quick(button):
        _start_smiles_conversion(apply_uff=False, quick=True)

    def handle_manta(button):
        # MANTA: full coordination-isomer manifold (with UFF cleanup) + GFN2
        # energy ranking, shown in the viewer with the existing isomer nav.
        # Read the exposed settings row (Quality/Seeds/Max-iso/Rank/Opt-top).
        _rank_sel = manta_rank.value
        # Opt dropdown -> top-N int: No = -1 (off, keep construction geometry); All = 0; Top-N = N.
        _opt_map = {'No': -1, 'Top 5': 5, 'Top 10': 10, 'Top 20': 20, 'All': 0}
        _start_smiles_conversion(
            apply_uff=True, quick=False,
            rank=(_rank_sel != 'No'),
            method=(_rank_sel if _rank_sel != 'No' else 'gfn2'),
            quality_mode=(manta_quality.value or None),
            # seeds field = preset value (transparent) unless the user edited it ->
            # custom seed count on top of the preset's cap/templates.
            seeds_override=(int(manta_seeds.value) or None),
            # 0 -> COMPLETE manifold (never cut off); else user cap
            max_isomers=(int(manta_max_iso.value) or 100000),
            opt_topn=_opt_map.get(manta_opt.value, -1),
            spin=str(manta_spin.value),     # 'auto' (scan) or fixed multiplicity (1/2/3/...)
            # Determinism toggle (default On) -> byte-identical, IDENTICAL to the CLI
            # and the development loop (ship = validate).  Off = non-deterministic embed.
            deterministic=(manta_det.value == 'On'),
            # construction always champion + power-user knobs CLI-only -> pinned here
            **_MANTA_DASH_DEFAULTS,
        )

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

    # -- atom-selection / manipulation handlers -------------------------
    def _run_manip_js(code):
        try:
            ctx.run_js(code)
        except Exception:
            pass

    def _apply_manip_mode_js(mode):
        _ensure_manip_bootstrap()
        _run_manip_js(
            f'if(window.__delfinSubmitManip) '
            f'window.__delfinSubmitManip.setMode({json.dumps(submit_scope_id)}, '
            f'{json.dumps(mode)});'
        )

    def _mode_buttons_mutex(keep):
        """Only one mode at a time; the others switch themselves off."""
        for button in (submit_select_btn, submit_manip_btn, submit_draw_btn):
            if button is not keep and button.value:
                button.value = False

    def _refresh_draw_controls():
        drawing = bool(submit_draw_btn.value)
        submit_element_dd.layout.display = '' if drawing else 'none'

    def on_submit_select_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_select_btn.value)
        submit_select_btn.button_style = 'info' if active else ''
        if active:
            _mode_buttons_mutex(submit_select_btn)
        _apply_manip_mode_js('select' if active else 'off')

    def on_submit_manip_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_manip_btn.value)
        submit_manip_btn.button_style = 'info' if active else ''
        if active:
            _mode_buttons_mutex(submit_manip_btn)
        _apply_manip_mode_js('manipulate' if active else 'off')

    def on_submit_draw_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_draw_btn.value)
        submit_draw_btn.button_style = 'info' if active else ''
        if active:
            _mode_buttons_mutex(submit_draw_btn)
        _refresh_draw_controls()
        _apply_manip_mode_js('draw' if active else 'off')
        if active:
            on_submit_draw_choice(None)

    def on_submit_draw_choice(_change=None):
        """Hand over the element the browser draws with.

        Not the bond order: a drawn bond is single, and what it should be is
        decided afterwards by tapping the stick, where it can be seen. Having
        to choose in advance was a setting that was almost always wrong.
        """
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setDrawElement('
            f'{json.dumps(submit_scope_id)},{json.dumps(submit_element_dd.value)});'
        )

    def _set_ff_notes(notes):
        """Show what the force field had to approximate, under the viewer."""
        rendered = [html.escape(str(note)) for note in notes if str(note).strip()]
        if rendered and _gfn.is_gfn_method(submit_ff_dd.value):
            # These come from the field that runs in the browser, which is UFF
            # whatever the method box says -- and read as though GFN had
            # produced them, which is how "GFN behaves exactly like UFF" gets
            # concluded from a panel that never claimed otherwise.
            label = _gfn.GFN_METHODS[str(submit_ff_dd.value)]['label']
            rendered.insert(0, html.escape(
                f'These notes are about the live field, which is UFF: it runs '
                f'in the browser so that dragging follows the mouse. '
                f'{label} runs when Optimise is pressed, and says so in the '
                f'status line when it has.'))
        if not rendered:
            submit_ff_notes.value = ''
            return
        items = ''.join(
            f'<li style="margin:0 0 2px 0;">{note}</li>' for note in rendered
        )
        submit_ff_notes.value = (
            "<div style='font-size:12px; line-height:1.4; color:#5a6570; "
            "background:#f6f7f9; border:1px solid #e0e4e8; border-radius:4px; "
            "padding:6px 10px;'>"
            "<b style='color:#455a64;'>Force field notes</b>"
            f"<ul style='margin:4px 0 0 16px; padding:0;'>{items}</ul>"
            "</div>"
        )

    def _ensure_ff_bootstrap():
        if state.get('ff_bootstrap_done'):
            return
        try:
            from .molecule_forcefield_js import molecule_ff_bootstrap_js
            ctx.run_js(molecule_ff_bootstrap_js())
            state['ff_bootstrap_done'] = True
        except Exception:
            pass

    def _structure_fingerprint(xyz):
        """Element column of an XYZ block -- what makes it the same molecule."""
        rows = [line.split() for line in (xyz or '').splitlines() if line.strip()]
        return tuple(r[0] for r in rows if len(r) >= 4)

    def _perception_for(xyz):
        """Perceive the bonding once per structure and keep it.

        Bond orders are read from the geometry, and a twisted double bond stops
        looking like one: turning a C=C by 30 degrees is enough for perception
        to call it a single bond, which drops the barrier holding the two
        halves coplanar from 19.5 to 1.1 kcal/mol. Everything downstream then
        lets the double bond rotate freely, and nothing brings it back.

        Editing moves atoms; it must not be able to change what the molecule
        is. So the bonding is perceived from the structure as it arrived and
        reused until a genuinely different one is loaded.
        """
        from .molecule_forcefield import perceive_molecule

        fingerprint = _structure_fingerprint(xyz)
        cached = state.get('perceived')
        if cached and state.get('perceived_for') == fingerprint:
            return cached
        perceived = perceive_molecule(xyz)
        _apply_bond_edits(perceived)
        # After the bond edits, never before: rebuilding the typing molecule
        # sanitizes it, and sanitisation re-perceives hybridisation.
        _apply_hyb_overrides(perceived)
        state['perceived'] = perceived
        state['perceived_for'] = fingerprint
        return perceived

    def _apply_bond_edits(perceived):
        """Lay the user's bond corrections over what perception found.

        The correction has to reach the molecules the force-field parameters
        are read from, not only the bond list -- otherwise a drawn bond keeps
        the length it was drawn at instead of contracting.
        """
        from .molecule_forcefield import apply_bond_edits

        apply_bond_edits(perceived, state.get('bond_edits') or {})

    def _apply_hyb_overrides(perceived):
        """Force the hybridisations the user has chosen by hand.

        Bond orders are perceived from the geometry, and a double bond that is
        not seen leaves its carbon typed sp3: angles at 109.5 degrees, and a
        centre that puckers where it should stay flat.
        """
        from .molecule_forcefield import apply_hybridisation_overrides

        apply_hybridisation_overrides(perceived, state.get('hyb_overrides') or {})

    def _install_gfn_frame_watcher():
        """Teach the page to play the trajectory, once.

        Two things were wrong with sending one frame at a time.  ipywidgets
        writes a new value into the DOM without firing an event, so there is
        nothing to listen for -- the value has to be read on a timer.  And the
        kernel produces frames faster than the reader can look, so all but the
        last of each burst were never seen: the structure jumped instead of
        moving.

        The field carries the whole trajectory instead, and the page keeps its
        place in it.  Nothing can be missed, because a frame that was written
        while the reader was busy is still there when it looks again.  Playback
        interpolates between frames -- twenty computed steps a second shown as
        continuous motion.  The positions in between are drawn, not calculated;
        every frame the structure actually passed through is one of the ends.
        """
        if state.get('gfn_watcher_installed'):
            return
        state['gfn_watcher_installed'] = True
        _emit_watcher_js(
            '(function(){\n'
            '  var scope=' + json.dumps(submit_scope_id) + ';\n'
            '  window.__delfinGfnPlay=window.__delfinGfnPlay||{};\n'
            '  if(window.__delfinGfnPlay[scope]) return;\n'
            '  var play={queue:[],at:0,started:0,last:null,seen:0,run:null};\n'
            '  window.__delfinGfnPlay[scope]=play;\n'
            '  var STEP_MS=55;\n'
            '  function stepMs(){\n'
            '    /* xtb computes faster than this plays: 75 frames arrive in\n'
            '       0.4 s and would take 4 s to show, so the picture trails the\n'
            '       calculation and keeps trailing it further.  A backlog is\n'
            '       played faster -- the whole path is still shown, in the time\n'
            '       the run actually takes. */\n'
            '    var n=play.queue.length;\n'
            '    if(n>60) return 8;\n'
            '    if(n>25) return 20;\n'
            '    if(n>10) return 35;\n'
            '    /* One answer at a time -- a hand being followed -- is drawn\n'
            '       over exactly the time the next answer takes to come, or the\n'
            '       picture sits still between them.  A relaxation arrives in\n'
            '       bursts instead, and those are paced by the rules above: a\n'
            '       burst of thirty drawn at a follow pace of a third of a\n'
            '       second each would put the picture ten seconds behind a\n'
            '       calculation that took one. */\n'
            '    if(play.follow&&play.gap&&n<=3) return play.gap;\n'
            '    return STEP_MS;\n'
            '  }\n'
            '  function read(arrivedAt){\n'
            '    /* Fullscreen moves the viewer into an overlay that carries\n'
            '       the same scope class, and the frame field is not one of the\n'
            '       things it takes -- so looking only inside the first element\n'
            '       with that class found the overlay and no field, and the\n'
            '       playback that worked in the small view showed nothing in\n'
            '       the big one.  Every element with the class is tried, and\n'
            '       the document stands behind them. */\n'
            '    var field=null;\n'
            '    var roots=document.querySelectorAll("."+scope);\n'
            '    for(var r=0;r<roots.length&&!field;r++){\n'
            '      field=roots[r].querySelector('
            '".submit-gfn-frame textarea, .submit-gfn-frame input");\n'
            '    }\n'
            '    if(!field) field=document.querySelector('
            '".submit-gfn-frame textarea, .submit-gfn-frame input");\n'
            '    if(!field) return;\n'
            '    var text=field.value||"";\n'
            '    if(!text){ play.queue=[]; play.seen=0; play.last=null;'
            ' play.run=null; return; }\n'
            '    var data=null;\n'
            '    try{ data=JSON.parse(text); }catch(e){ return; }\n'
            '    if(data&&data.halt){\n'
            '      /* The run was switched off.  Playing out the queue after\n'
            '         that is the picture carrying on without the thing it is\n'
            '         a picture of. */\n'
            '      play.queue=[]; play.seen=(data.frames||[]).length;\n'
            '      if(!play.toldStop){ play.toldStop=1;\n'
            '        /* Which frame is on screen.  Stopping keeps that one:\n'
            '           frames xtb had already computed but nobody had seen\n'
            '           are not what the user stopped at. */\n'
            '        say("stopped at frame "+(play.shown||0)); }\n'
            '      return;\n'
            '    }\n'
            '    var frames=(data&&data.frames)||[];\n'
            '    var run=(data&&data.run)||0;\n'
            '    /* Whether these frames belong to a molecule following a hand\n'
            '       rather than to a minimisation.  The two are told apart\n'
            '       because Optimise is not pressed during a follow, and the\n'
            '       check that abandons a queue when that switch goes up would\n'
            '       otherwise throw every followed frame away. */\n'
            '    play.follow=(data&&data.follow)?1:0;\n'
            '    if(run!==play.run){\n'
            '      /* A new run. Without this the count of frames already\n'
            '         played carried over, so a shorter run than the one\n'
            '         before it played nothing at all -- which is what made\n'
            '         the playback look like it worked only sometimes. */\n'
            '      if(play.queue.length){\n'
            '        /* Land the run that is ending on its last frame first.\n'
            '           Dropped there, the viewer keeps whatever it had drawn\n'
            '           while the kernel keeps the geometry it computed, and\n'
            '           the two drift apart: the next drag then hands over a\n'
            '           structure that is behind, and the relaxation nobody saw\n'
            '           is walked again -- which is every earlier drag being\n'
            '           made a second time. */\n'
            '        show(play.last,play.queue[play.queue.length-1],1);\n'
            '      }\n'
            '      play.run=run; play.seen=0; play.queue=[]; play.last=null;\n'
            '      play.shown=0; play.toldStop=0;\n'
            '    }\n'
            '    /* Where in the run these frames start.  A long run sends the\n'
            '       tail rather than the whole path -- every write is a message\n'
            '       and the whole path grows without end -- so counting from\n'
            '       the front of the message would replay frames already shown\n'
            '       and then stop showing new ones altogether. */\n'
            '    var from=(data&&data.from)||0;\n'
            '    if(play.follow&&play.seen===0&&frames.length>1){\n'
            '      /* A live run arriving from the start: only its newest frame\n'
            '         is worth anything.  The ones before it describe where the\n'
            '         structure was on the way here -- the drag that has just\n'
            '         happened, or a relaxation of it -- and playing those is\n'
            '         showing the user their own past.  What is wanted is the\n'
            '         frame that is current. */\n'
            '      from=from+frames.length-1;\n'
            '      frames=[frames[frames.length-1]];\n'
            '    }\n'
            '    if(from+frames.length>play.seen){\n'
            '      var start=Math.max(0,play.seen-from);\n'
            '      /* A frame arriving into an empty queue starts its own\n'
            '         interpolation.  Left over from the last one, the clock\n'
            '         reads far past a step and the frame is jumped to instead\n'
            '         of moved to -- which is every frame of a follow, where\n'
            '         they arrive further apart than a step. */\n'
            '      if(!play.queue.length) play.started=0;\n'
            '      /* How long the machine took over the last answer, which is\n'
            '         how long this one has to be drawn over. */\n'
            '      if(play.arrived){\n'
            '        var measured=Math.min(600,Math.max(24,\n'
            '          arrivedAt-play.arrived));\n'
            '        /* Averaged, not taken raw.  Answers do not arrive evenly,\n'
            '           and pacing each frame by the last interval alone makes\n'
            '           the drawing speed jump about as much as the arrivals\n'
            '           do -- which is felt as jerkiness even when nothing is\n'
            '           being missed. */\n'
            '        play.gap=play.gap?(play.gap*0.6+measured*0.4):measured;\n'
            '      }\n'
            '      play.arrived=arrivedAt;\n'
            '      for(var i=start;i<frames.length;i++) play.queue.push(frames[i]);\n'
            '      play.seen=from+frames.length;\n'
            '      say("received "+play.seen+" frames");\n'
            '    }\n'
            '    /* xtb produces frames faster than any frame rate can show\n'
            '       them, so a queue that is allowed to grow puts the picture\n'
            '       permanently behind the calculation.  Beyond a second of\n'
            '       backlog the older ones are skipped: what is on screen is\n'
            '       then always close to what xtb is doing, which is the point\n'
            '       of watching at all. */\n'
            '    if(play.queue.length>20){\n'
            '      play.queue=play.queue.slice(-20);\n'
            '    }\n'
            '  }\n'
            '  function say(text){ send("gfnplay", text); }\n'
            '  function send(verb,text){\n'
            '    /* The page reports back through the command bridge the editor\n'
            '       already has, so a playback that does not appear says why by\n'
            '       itself instead of being read out of a console. */\n'
            '    try{\n'
            '      var root=document.querySelector("."+scope);\n'
            '      var wrap=root&&root.querySelector(".submit-cmd-sync");\n'
            '      var input=wrap&&wrap.querySelector("input, textarea");\n'
            '      if(!input) return;\n'
            '      play.serial=(play.serial||0)+1;\n'
            '      var proto=(input.tagName==="TEXTAREA")\n'
            '        ? window.HTMLTextAreaElement.prototype\n'
            '        : window.HTMLInputElement.prototype;\n'
            '      var setter=Object.getOwnPropertyDescriptor(proto,"value");\n'
            '      var line=verb+":"+play.serial+":"+text;\n'
            '      if(setter&&setter.set) setter.set.call(input,line);\n'
            '      else input.value=line;\n'
            '      input.dispatchEvent(new Event("input",{bubbles:true}));\n'
            '      input.dispatchEvent(new Event("change",{bubbles:true}));\n'
            '    }catch(e){}\n'
            '  }\n'
            '  function show(a,b,t){\n'
            '    if(!window.__delfinSubmitManip||'
            '!window.__delfinSubmitManip.setPositions){\n'
            '      if(!play.toldNoApi){ play.toldNoApi=1;'
            ' say("no setPositions on the page"); }\n'
            '      return;\n'
            '    }\n'
            '    var out=new Array(b.length);\n'
            '    if(!a||a.length!==b.length){ out=b; }\n'
            '    else { for(var i=0;i<b.length;i++) out[i]=a[i]+(b[i]-a[i])*t; }\n'
            '    var ok=false;\n'
            '    try{ ok=window.__delfinSubmitManip.setPositions('
            'scope,out,heldSerials()); }\n'
            '    catch(e){ ok=false; }\n'
            '    play.drawn=(play.drawn||0)+(ok?1:0);\n'
            '    if(!ok&&!play.toldNoDraw){ play.toldNoDraw=1;'
            ' say("setPositions did not draw"); }\n'
            '    if(ok&&!play.toldDrawing){ play.toldDrawing=1;'
            ' say("drawing"); }\n'
            '  }\n'
            '  function grabbed(){\n'
            '    /* Whether an atom is being moved right now.  A playback that\n'
            '       keeps writing positions during a drag puts the atom back\n'
            '       where xtb had it once per animation frame, so it cannot be\n'
            '       moved at all -- the drag has to own the picture while it\n'
            '       lasts.  A rectangle being pulled over the molecule selects\n'
            '       and moves nothing, so it is not a grab. */\n'
            '    var held=(window._submitManipStateByScope||{})[scope];\n'
            '    var drag=held&&held.drag;\n'
            '    if(!drag) return false;\n'
            '    return drag.kind==="translate"||drag.kind==="rotate"'
            '||drag.kind==="draw";\n'
            '  }\n'
            '  function heldSerials(){\n'
            '    /* The atoms the hand is moving.  Coordinates that come back\n'
            '       describe where they were when they were sent, and the\n'
            '       cursor has moved on since -- so those are the ones the\n'
            '       playback must leave alone. */\n'
            '    var st=(window._submitManipStateByScope||{})[scope];\n'
            '    return (st&&st.drag&&st.drag.targets)||[];\n'
            '  }\n'
            '  function followIsOn(){\n'
            '    /* Relax, read off the page rather than asked of the kernel --\n'
            '       the same reason as the switch below. */\n'
            '    var holder=document.querySelector(".submit-gfn-follow");\n'
            '    if(!holder) return false;\n'
            '    var btn=(holder.tagName==="BUTTON")?holder'
            ':holder.querySelector("button");\n'
            '    if(!btn) return false;\n'
            '    return btn.classList.contains("mod-active");\n'
            '  }\n'
            '  function switchIsOn(){\n'
            '    /* ipywidgets marks a pressed toggle with mod-active.  Reading\n'
            '       it here is instant; asking the kernel costs a round trip,\n'
            '       and the picture ran on for the length of it. */\n'
            '    var holder=document.querySelector(".submit-optimize-switch");\n'
            '    if(!holder) return true;\n'
            '    var btn=(holder.tagName==="BUTTON")?holder'
            ':holder.querySelector("button");\n'
            '    if(!btn) return true;\n'
            '    return btn.classList.contains("mod-active");\n'
            '  }\n'
            '  function frame(now){\n'
            '    /* An atom picked up while xtb is running: the kernel is told\n'
            '       at the grab rather than at the release, because a GFN2 run\n'
            '       is thirteen seconds and every one of them would be spent\n'
            '       minimising a structure the user is in the middle of\n'
            '       changing.  The queue goes with it -- those frames belong to\n'
            '       the geometry that has just been altered. */\n'
            '    var held=grabbed();\n'
            '    if(held!==!!play.held){\n'
            '      play.held=held?1:0;\n'
            '      if(held){ play.queue=[]; play.last=null; }\n'
            '      else {\n'
            '        /* Let go of.  What is still queued describes the drag:\n'
            '           each of those frames carries the dragged atom where the\n'
            '           hand had it when the frame was computed, and while the\n'
            '           hand was down they were drawn around it.  Drawn now,\n'
            '           with nothing held any more, they walk that atom back\n'
            '           through the whole drag in front of the user.  Land on\n'
            '           the newest -- which is where the hand actually left it\n'
            '           -- and drop the rest. */\n'
            '        if(play.queue.length){\n'
            '          show(play.last,play.queue[play.queue.length-1],1);\n'
            '        }\n'
            '        play.queue=[]; play.last=null;\n'
            '        play.follow=0; play.pushed=0;\n'
            '      }\n'
            '      send(held?"gfngrab":"gfnfree","");\n'
            '    }\n'
            '    if(play.held&&!followIsOn()){\n'
            '      window.requestAnimationFrame(frame); return;\n'
            '    }\n'
            '    if(play.held){\n'
            '      /* Following: the geometry goes over while the mouse is\n'
            '         still down, and the pace is set by the machine rather\n'
            '         than by a clock.  The next one goes as soon as the last\n'
            '         answer has landed -- GFN-FF answers a small molecule in\n'
            '         under twenty milliseconds and a fixed fifth of a second\n'
            '         threw nine tenths of that away.  A floor keeps the\n'
            '         messages from becoming the bottleneck, and a ceiling\n'
            '         starts again if an answer never comes at all. */\n'
            '      /* The floor is the machine\'s own answer time, not a\n'
            '         constant: never ask more than twice as often as it can\n'
            '         answer, and never faster than an animation frame can show.\n'
            '         On a small molecule that is about twenty milliseconds; on\n'
            '         a hundred atoms it settles at sixty by itself. */\n'
            '      var floor=Math.max(16,Math.min(120,(play.gap||60)/2));\n'
            '      var since=now-(play.pushed||0);\n'
            '      var answered=play.seen>(play.pushedAt||0);\n'
            '      if(since>500||(answered&&since>floor)){\n'
            '        play.pushed=now; play.pushedAt=play.seen;\n'
            '        var api=window.__delfinSubmitManip;\n'
            '        if(!api||!api.pushXyz){\n'
            '          /* An editor from before the follow existed.  Swallowed,\n'
            '             this is a drag that does nothing and says nothing --\n'
            '             which is indistinguishable from a broken kernel. */\n'
            '          if(!play.toldNoPush){ play.toldNoPush=1;\n'
            '            say("this page has an editor without pushXyz; '
            'reload the page"); }\n'
            '        } else {\n'
            '          try{ api.pushXyz(scope,"drag-follow"); }\n'
            '          catch(e){ if(!play.toldNoPush){ play.toldNoPush=1;\n'
            '            say("pushXyz failed: "+e); } }\n'
            '        }\n'
            '      }\n'
            '    }\n'
            '    if(play.queue.length&&!play.follow&&!switchIsOn()){\n'
            '      /* Optimise going up abandons what it had computed but\n'
            '         nobody had seen.  A follow has no Optimise behind it --\n'
            '         checking that switch there threw away every frame the\n'
            '         drag produced, which is what "it does nothing" was. */\n'
            '      play.queue=[];\n'
            '      if(!play.toldStop){ play.toldStop=1;\n'
            '        say("stopped at frame "+(play.shown||0)); }\n'
            '    }\n'
            '    read(now);\n'
            '    if(play.queue.length){\n'
            '      if(!play.started) play.started=now;\n'
            '      var t=(now-play.started)/stepMs();\n'
            '      if(t>=1){\n'
            '        var next=play.queue.shift();\n'
            '        show(play.last,next,1);\n'
            '        play.last=next; play.started=now;\n'
            '        play.shown=(play.shown||0)+1;\n'
            '      } else if(play.last){\n'
            '        show(play.last,play.queue[0],t);\n'
            '      } else {\n'
            '        play.last=play.queue.shift(); show(null,play.last,1);\n'
            '        play.started=now;\n'
            '      }\n'
            '    }\n'
            '    window.requestAnimationFrame(frame);\n'
            '  }\n'
            '  window.requestAnimationFrame(frame);\n'
            '})();'
        )

    def _emit_watcher_js(script):
        """Send the player with the start-up scripts, not through run_js.

        run_js writes into a single Output and clears it first, so a script
        sent at click time can be replaced before the page has run it -- which
        is how the player came to be missing while everything it depends on
        was working.  add_init_js is the channel the explorer's own JS arrives
        on, and that one has never been in doubt.
        """
        try:
            ctx.add_init_js(script)
            return
        except Exception:
            pass
        _run_manip_js(script)

    #: A follow step is a whole xtb process, so it is a few cycles rather than
    #: a minimisation.  Measured on a 102-atom complex: one cycle 0.06 s, five
    #: 0.09 s, ten 0.12 s -- five is about ten answers a second, which reads as
    #: the molecule following the hand rather than catching up with it.
    _GFN_FOLLOW_CYCLES = 5

    def _begin_gfn_follow():
        """A drag has started and the molecule is to follow it."""
        if not (submit_relax_btn.value
                and _gfn.is_gfn_method(submit_ff_dd.value)):
            return False
        # The method on screen is the method that runs.  It used to be GFN-FF
        # whatever the box said, which is a picture of a calculation nobody
        # asked for -- and indistinguishable, from the outside, from the right
        # one.
        state['gfn_follow_method'] = str(submit_ff_dd.value)
        if state.get('gfn_follow'):
            return True     # already following; it keeps the run it began
        _gfn_new_generation()
        # What the molecule looked like before this drag: the bonding is read
        # from here, not from a frame that has already been pulled about.
        state['gfn_topology_source'] = (
            state.get('current_xyz_for_copy') or {}).get('content')
        state['gfn_follow'] = True
        state['gfn_follow_steps'] = 0
        state['gfn_follow_frames'] = []
        run = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = run
        state['gfn_follow_run'] = run
        return True

    def _end_gfn_follow():
        state['gfn_follow'] = False

    def _gfn_follow_step(xyz, holding=()):
        """Relax around the atom the hand is holding, and send that back.

        The dragged atom is not held by xtb: it cannot be.  ``$fix atoms:`` is
        broken in xtb 6.7.1 and ``$constrain atoms:`` naming one atom does
        nothing at all, both measured -- xtb holds internal coordinates, not
        positions.  The hold is the page's instead: the frames that come back
        are written to every atom *except* the ones being dragged, so the
        cursor keeps the one it is holding and xtb arranges the rest around it.

        Only one process at a time, and the newest geometry wins: a hand moves
        faster than xtb answers, and a queue of answers about where the atom
        used to be is worse than no answer at all.
        """
        if not state.get('gfn_follow') and not _begin_gfn_follow():
            return          # not following: Relax is up, or GFN is not chosen
        state['gfn_follow_xyz'] = (xyz, tuple(holding or ()))
        if state.get('gfn_follow_busy'):
            return
        state['gfn_follow_busy'] = True
        method = str(state.get('gfn_follow_method') or submit_ff_dd.value)
        label = _gfn.GFN_METHODS[method]['label']
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        constraints = list(state.get('constraints') or [])
        wet = str(submit_gfn_solvent.value or '') or None

        def _work():
            try:
                while state.get('gfn_follow'):
                    newest = state.pop('gfn_follow_xyz', None)
                    if newest is None:
                        return
                    current, holding = newest
                    began = time.perf_counter()
                    outcome = _gfn.relax_steps(
                        current, method=method, charge=charge, uhf=uhf,
                        cycles=_GFN_FOLLOW_CYCLES, timeout=30.0,
                        constraints=constraints, solvent=wet,
                        topology=_gfn_topology_dir(
                            len(_gfn.atom_lines(current))),
                    )
                    if not outcome.get('ok'):
                        note = str(outcome.get('status') or 'it did not run')
                        _schedule_ui_update(
                            _set_mol_status,
                            f'The molecule stopped following: {note}')
                        return
                    # Say it out loud, every step.  A follow that is working
                    # and a follow that is not both look like a molecule that
                    # is not moving much, and the difference is the whole
                    # question when something in the chain is broken.
                    steps = int(state.get('gfn_follow_steps') or 0) + 1
                    state['gfn_follow_steps'] = steps
                    # How many atoms the hand is on.  Grabbing an atom that is
                    # part of a selection drags the whole selection, so a
                    # selection left over from earlier makes every drag move
                    # everything -- which reads as the molecule fighting
                    # itself, and is invisible unless it is counted here.
                    many = (f' holding {len(holding)} atoms,'
                            if len(holding) > 1 else '')
                    said = (f'{label} is following the drag:{many} {steps} '
                            f'step(s), '
                            f'{(time.perf_counter() - began) * 1000:.0f} ms '
                            'each.')
                    state['gfn_last_status'] = said
                    _schedule_ui_update(_set_mol_status, said, spinner=True)
                    # The atoms under the cursor go back where they were sent.
                    # xtb pulls them most of the way home in five cycles, and
                    # this answer outlives the drag -- applied after the
                    # release it would take them with it, which is the spring
                    # back that looked like the drag being undone.
                    settled = _gfn.hold_atoms_at(
                        outcome['xyz'], current, holding)
                    frames = list(state.get('gfn_follow_frames') or [])
                    frames.append(_gfn.coordinates_of(settled))
                    state['gfn_follow_frames'] = frames
                    # The tail, not the whole drag: every write is a message,
                    # and a minute of dragging is three hundred frames of a
                    # hundred atoms in each of them.
                    trail = frames[-40:]
                    payload = json.dumps({'run': state.get('gfn_follow_run'),
                                          'from': len(frames) - len(trail),
                                          'follow': 1, 'frames': trail})
                    _schedule_ui_update(
                        lambda text=payload: setattr(
                            submit_gfn_frame, 'value', text))
            finally:
                state['gfn_follow_busy'] = False

        threading.Thread(target=_work, daemon=True).start()

    #: Letting go with Settle on.  More cycles than a follow step, because
    #: nothing is being held any more and the structure should reach somewhere
    #: it would stay rather than take one step towards it -- but still a bound,
    #: because this happens on every release.
    _GFN_SETTLE_CYCLES = 40

    #: How many rounds before it gives up on converging.  Held values that
    #: cannot all be met at once never converge -- measured on a propane with
    #: two fixed distances and an angle fighting each other, not converged in
    #: any round -- and a relaxation that will not end is a process per round
    #: for as long as the switch is down, and a structure that visibly jitters.
    _GFN_SETTLE_ROUNDS = 12
    #: And a round that moved nothing has settled, whatever xtb calls it.
    _GFN_SETTLE_STILL = 0.005

    def _gfn_topology_dir(atoms):
        """Where GFN-FF's perceived bonding is kept while a structure is worked
        on.

        GFN-FF reads its topology out of the geometry it is handed, once.  On a
        drag that is a fresh perception per step, and it has a cliff in it:
        measured on a propane, a C-C pulled to 1.87 A relaxes back to 1.49,
        and at 1.96 A the bond is not seen at all and the same relaxation
        pushes it to 2.80.  So the perception is made once and kept, the way
        the browser's field assigns its parameters once when the switch goes
        down -- at 2.33 A it then still pulls back, to 1.51.

        It belongs to one molecule: an atom added or taken away makes it a
        different one, and the count is what says so.
        """
        import tempfile

        kept = state.get('gfn_topology')
        if kept and kept.get('atoms') == atoms:
            return Path(kept['dir'])
        _drop_gfn_topology()
        folder = tempfile.mkdtemp(prefix='delfin-gfnff-topo-')
        state['gfn_topology'] = {'dir': folder, 'atoms': atoms}
        # Perceived here and now, from the structure as it stood before a hand
        # was laid on it.  Left to the first calculation that needs it, the
        # perception is made from a geometry that has already been pulled
        # about -- and if the drag has gone past where a bond is still
        # recognised, the topology that is then kept for the whole session is
        # one with the bond missing.  Measured: a propane whose C-C had been
        # pulled to 2.1 A perceived that way came back at 3.57 A.
        seed = (state.get('gfn_topology_source')
                or (state.get('current_xyz_for_copy') or {}).get('content'))
        if seed and len(_gfn.atom_lines(seed)) == atoms:
            try:
                _gfn.relax_steps(seed, cycles=1, timeout=30.0,
                                 topology=Path(folder))
            except Exception:
                pass
        return Path(folder)

    def _drop_gfn_topology():
        """Forget the bonding: the molecule is not the same one any more."""
        kept = state.pop('gfn_topology', None)
        if kept:
            shutil.rmtree(kept.get('dir') or '', ignore_errors=True)

    def _gfn_new_generation():
        """Everything in flight belongs to the structure it was started for.

        A drag, a Hold, a Set: each makes a new structure, and every timer and
        every worker started for the one before it is now about something that
        no longer exists.  They used to be left to finish -- writing their
        geometry over the new one, or holding a flag that made the next
        relaxation skip itself, which is why switching the toggle off and on
        again finished the job that letting go should have.  One counter
        settles all of it: whatever does not belong to the current generation
        does nothing at all.
        """
        state['gfn_generation'] = int(state.get('gfn_generation', 0)) + 1
        state['gfn_settle_forced'] = False
        state['gfn_settle_rounds'] = 0
        state['gfn_settle_again'] = False
        return state['gfn_generation']

    def _gfn_uhf_now():
        """How many unpaired electrons the live relaxation should assume.

        The box says a multiplicity and xtb counts unpaired electrons, so
        M = 2S+1 becomes M - 1.  With auto M on, the box is not the answer at
        all -- Optimise scans and keeps the lowest -- and a live relaxation
        running the box's fixed M while Optimise ran a scanned one is two
        answers about two different molecules, which is why pressing Optimise
        after it moved the structure again.  Scanning every round would cost
        three runs a round, so what the last scan settled on is used instead.
        """
        if submit_gfn_autospin.value and state.get('gfn_scanned_uhf') is not None:
            return int(state['gfn_scanned_uhf'])
        return max(0, int(submit_gfn_mult.value or 1) - 1)

    def _gfn_live_is_on():
        """Whether something on screen is meant to act on a change at once."""
        return (submit_relax_btn.value
                and _gfn.is_gfn_method(submit_ff_dd.value)
                and _gfn.find_binary(submit_ff_dd.value) is not None)

    def _arm_gfn_takeup(note=''):
        """Take up a change to what is held, straight away.

        Setting a value, watching nothing happen and having to press Optimise
        for it is the switch claiming to be live and not being it.  Under the
        browser's field a held value acts the moment it is set, because the
        field is already running; here nothing is running between drags, so
        the change is what starts it.
        """
        if not _gfn_live_is_on():
            return
        _gfn_new_generation()
        state['gfn_settle_note'] = note
        _arm_gfn_settle(forced=True)

    def _arm_gfn_settle(forced=False):
        """Tidy the structure after a release, with the method on screen.

        The same promise Settle makes under the browser's field: what reaches
        the coordinate box is a structure that has relaxed around where the
        atom was put, rather than wherever the cursor happened to stop --
        measured at 176 kcal/mol above a settled one on a real complex.

        Waits a moment first, because the coordinates of the release arrive on
        a message of their own and this has to start from them.
        """
        if forced:
            # Asked for by a hand -- a value held, a value set.  That is not a
            # tidy-up to be skipped when something else is in the air; it is
            # the answer to something the user just did, and it has to happen.
            state['gfn_settle_forced'] = True
        if not (_gfn.is_gfn_method(submit_ff_dd.value)
                and _gfn.find_binary(submit_ff_dd.value) is not None
                and (submit_settle_btn.value
                     or state.get('gfn_settle_forced'))):
            return
        if not forced and state.get('optimize_interrupted') is not None:
            # An optimisation was interrupted by this very drag and is coming
            # back.  A whole minimisation is more than a settle, not less.
            return
        # A round that follows another round waits only long enough not to
        # spin; a release waits for its coordinates to arrive first.
        state['gfn_settle_at'] = time.monotonic() + (
            0.05 if state.get('gfn_settle_again') else _GFN_RESTART_DELAY)
        if state.get('gfn_settle_armed'):
            return
        state['gfn_settle_armed'] = True
        generation = int(state.get('gfn_generation', 0))

        def _wait():
            while True:
                left = state.get('gfn_settle_at', 0.0) - time.monotonic()
                if left <= 0:
                    break
                time.sleep(min(left, 0.05))
            state['gfn_settle_armed'] = False
            if int(state.get('gfn_generation', 0)) != generation:
                return          # the structure it was armed for is history
            _schedule_ui_update(_gfn_settle_now)

        threading.Thread(target=_wait, daemon=True).start()

    def _gfn_settle_now():
        if state.get('gfn_follow'):
            return                      # a new drag started in the meantime
        if state.get('gfn_settle_busy'):
            # One is already running.  What has just been asked for is a
            # different question from the one it is answering, so it is asked
            # again the moment that one is done -- dropped here, a value held
            # while something was in flight did nothing at all until the
            # toggle was switched off and on again.
            state['gfn_settle_pending'] = True
            return
        if state.get('optimize_run') is not None:
            return                      # Optimise is doing this already
        generation = int(state.get('gfn_generation', 0))
        method = str(submit_ff_dd.value)
        xyz = (state.get('current_xyz_for_copy') or {}).get('content')
        if not xyz or not _gfn.is_gfn_method(method):
            return
        label = _gfn.GFN_METHODS[method]['label']
        state['gfn_settle_busy'] = True
        charge = int(submit_gfn_charge.value or 0)
        uhf = _gfn_uhf_now()
        constraints = list(state.get('constraints') or [])
        wet = str(submit_gfn_solvent.value or '') or None
        rounds = int(state.get('gfn_settle_rounds') or 0) + 1
        state['gfn_settle_rounds'] = rounds
        # One run number for the whole relaxation, not one per round.  A new
        # run number resets the player, which drops what it had not drawn yet
        # and applies the next frame outright instead of moving to it: with a
        # round every few tenths of a second that is a twitch per round, and
        # what it looks like is the structure jittering.
        if rounds == 1:
            state['gfn_run'] = int(state.get('gfn_run', 0)) + 1
            state['gfn_settle_offset'] = 0
        run = int(state.get('gfn_run', 0))
        offset = int(state.get('gfn_settle_offset') or 0)
        note = str(state.get('gfn_settle_note') or '')
        _set_mol_status(
            (f'{note}: {label} is moving the structure to it...' if note
             else f'{label} is settling the structure...')
            + (f' (round {rounds})' if rounds > 1 else ''), spinner=True)
        perceived = _gfn_topology_dir(len(_gfn.atom_lines(xyz)))
        # Auto M, and no scan has happened yet: this run does the scanning, so
        # that it and Optimise are asking about the same molecule.
        scanning = bool(submit_gfn_autospin.value
                        and state.get('gfn_scanned_uhf') is None)

        def _settle_stopped():
            # A hand on an atom takes over from a relaxation of the whole
            # thing; Optimise is the same calculation asked for by hand and
            # supersedes it outright -- two xtb runs writing the same
            # coordinate box is how the first press came to look like it had
            # only worked out an energy.  And switching off means what it says,
            # whichever of the two switches was keeping this alive.
            return bool(state.get('gfn_follow')
                        or state.get('optimize_run') is not None
                        or not (submit_settle_btn.value or _gfn_live_is_on()))

        def _push(frames):
            walked = list(frames)
            trail = walked[-60:]
            state['gfn_settle_walked'] = len(walked)
            _schedule_ui_update(
                lambda t=trail, f=offset + len(walked) - len(trail): setattr(
                    submit_gfn_frame, 'value',
                    json.dumps({'run': run, 'from': f, 'follow': 1,
                                'frames': t})))

        def _work():
            began = time.perf_counter()
            if scanning:
                # Auto M with nothing scanned yet.  Optimise would scan and
                # keep the lowest; running the box's M here instead makes the
                # two answer about different molecules, and pressing Optimise
                # afterwards then moves the structure for no visible reason.
                # It costs three runs, once -- after that the answer is known.
                outcome = _gfn.optimize_autospin(
                    xyz, method, charge=charge, constraints=constraints,
                    on_frames=_push, topology=perceived, timeout=None,
                    solvent=wet, should_stop=_settle_stopped,
                )
                if outcome.get('uhf') is not None:
                    state['gfn_scanned_uhf'] = int(outcome['uhf'])
            else:
                # No cycle cap and no clock: this is the ordinary
                # optimisation, run on the frame that is on screen now, on the
                # same terms Optimise runs on.  Chopping it into rounds bought
                # nothing -- the geometry between two rounds is not a place
                # anyone wants to stop at -- and cost a stutter at every
                # boundary and an early ending whenever a round moved little.
                outcome = _gfn.optimize_with_gfn(
                    xyz, method, charge=charge, uhf=uhf,
                    max_steps=None, timeout=None,
                    constraints=constraints, on_frames=_push, solvent=wet,
                    topology=perceived, should_stop=_settle_stopped,
                )

            def _done():
                state['gfn_settle_busy'] = False
                if state.pop('gfn_settle_pending', False):
                    # Something was asked for while this was running.
                    _schedule_ui_update(lambda: _arm_gfn_settle(forced=True))
                if (int(state.get('gfn_generation', 0)) != generation
                        or state.get('optimize_run') is not None):
                    # The structure moved on while this was running, or
                    # Optimise took over.  Either way this geometry is about a
                    # molecule that is no longer the current one, and writing
                    # it would be the past reaching into the present.
                    return
                state['gfn_settle_again'] = False
                # Where the next round's frames carry on from.
                state['gfn_settle_offset'] = offset + int(
                    state.pop('gfn_settle_walked', 0) or 0)
                if not outcome.get('ok'):
                    state['gfn_settle_forced'] = False
                    state['gfn_settle_rounds'] = 0
                    _set_mol_status(
                        f'{label} could not settle it: '
                        f'{outcome.get("status") or "it did not run"}')
                    return
                lines = [line for line in outcome['xyz'].splitlines()[2:]
                         if line.strip()]
                if lines:
                    # The playback has drawn this already; the box is what Copy
                    # and Submit read, and it has to be true whether or not a
                    # frame happened to land.
                    state['manip_inflight'] = True
                    coords_widget.value = (
                        f'{len(lines)}\nSettled with {label}\n'
                        + '\n'.join(lines))
                # Not converged and the switch is still down: keep going.  That
                # is what makes this a relaxation rather than a single push.
                # It ends three ways -- converged, standing still, or out of
                # rounds -- because held values that cannot all be met at once
                # never converge, and a relaxation that will not end is a
                # process per round for as long as the switch is down.
                moved = _gfn.largest_shift(xyz, outcome['xyz'])
                if (not outcome.get('converged') and _gfn_live_is_on()
                        and not state.get('gfn_follow')
                        and rounds < _GFN_SETTLE_ROUNDS
                        and moved > _GFN_SETTLE_STILL):
                    state['gfn_settle_again'] = True
                    state['gfn_settle_forced'] = True
                    _arm_gfn_settle()
                    return
                state['gfn_settle_forced'] = False
                state['gfn_settle_note'] = ''
                took = time.perf_counter() - began
                if outcome.get('converged'):
                    said = (f'{label} relaxed the structure: converged after '
                            f'{rounds} round(s).')
                elif not _gfn_live_is_on() or state.get('gfn_follow'):
                    said = f'{label} settled the structure in {took:.1f} s.'
                else:
                    # Said out loud, because a structure that has stopped
                    # improving and one that is finished look identical, and
                    # only one of them is worth pressing Optimise on.
                    why = ('it is no longer moving' if moved <= _GFN_SETTLE_STILL
                           else f'{rounds} rounds is as far as this goes')
                    said = (
                        f'{label} stopped without converging: {why}. '
                        + ('Held values that cannot all be met at once are the '
                           'usual reason -- the list beside the toolbar is '
                           'what it is trying to satisfy.'
                           if state.get('constraints') else
                           'Optimise will run it without a round limit.'))
                state['gfn_settle_rounds'] = 0
                state['gfn_last_status'] = said
                _set_mol_status(said)

            _schedule_ui_update(_done)

        threading.Thread(target=_work, daemon=True).start()

    def _enable_live_forcefield():
        """Assign UFF parameters for the geometry now in the viewer.

        Runs once, when the toggle is switched on -- never during a drag. The
        browser relaxes from these parameters alone; a round trip per frame
        would cap the drag at about 13 Hz.

        Refused outright while a GFN method is chosen.  A dozen handlers call
        this -- Hold, a polyhedron, a hybridisation, a bond edit -- and any one
        of them installing UFF parameters put a UFF relaxation under a molecule
        whose method box said GFN: it settled on release, and the geometry that
        reached the coordinate box was UFF's. The method that is chosen is the
        method that acts, and nothing else may touch the structure.
        """
        if _gfn.is_gfn_method(submit_ff_dd.value):
            _stop_browser_field()
            return
        xyz = (state.get('current_xyz_for_copy') or {}).get('content')
        if not xyz:
            _set_mol_status('Load a structure before enabling Relax.')
            submit_relax_btn.value = False
            return
        try:
            from .molecule_forcefield import export_forcefield_terms
            polyhedron = None
            if state.get('poly_applied') and state.get('poly_metal') is not None:
                polyhedron = {
                    'metal': state['poly_metal'],
                    'geometry': state['poly_applied'],
                    # None means: work it out from where the ligands are now.
                    'assignment': state.get('poly_assignment'),
                }
            payload = export_forcefield_terms(
                xyz,
                perceived=_perception_for(xyz),
                # The live field is what the browser relaxes with on every
                # frame, so it is always one of the two that live there.
                method=_live_ff_method(),
                polyhedron=polyhedron,
                restraints=[
                    c for c in (state.get('constraints') or [])
                    if c.get('mode', 'pull') == 'pull'
                ],
            )
        except Exception as exc:
            _set_mol_status(f'Force field unavailable: {exc}')
            submit_relax_btn.value = False
            return
        if not payload.get('ok'):
            _set_mol_status('Force field could not be assigned for this structure.')
            submit_relax_btn.value = False
            return
        _ensure_manip_bootstrap()
        _ensure_ff_bootstrap()
        _push_bond_orders()
        # The resume flag is set here, in the same script that hands over the
        # parameters, and not earlier: setting it before the re-render meant
        # it could be consumed against the viewer that was going away, and the
        # relaxation came back stuck until the toggle was cycled by hand.
        resume = 'true' if submit_relax_btn.value else 'false'
        _run_manip_js(
            f'window.__delfinResumeAutoOpt = {resume};'
            'if(window.__delfinSubmitManip){'
            'window.__delfinSubmitManip.setForceField('
            f'{json.dumps(submit_scope_id)},{json.dumps(payload)});'
            'window.__delfinSubmitManip.setOptimizerStrength('
            f'{json.dumps(submit_scope_id)},{int(submit_strength_slider.value)});'
            'window.__delfinSubmitManip.setSettleOnRelease('
            f'{json.dumps(submit_scope_id)},'
            f'{"true" if submit_settle_btn.value else "false"});'
            'window.__delfinSubmitManip.setFixedInternals('
            f'{json.dumps(submit_scope_id)},'
            + json.dumps([
                {'kind': c['kind'], 'atoms': c['atoms'], 'value': c['value']}
                for c in (state.get('constraints') or [])
                if c.get('mode') == 'fix'
            ])
            + ');'
            '}'
        )
        # Terms derived from the input geometry rather than real UFF typing --
        # the transition-metal case -- are worth saying out loud, and they
        # belong under the structure they describe rather than in the preview's
        # status line, which conversion messages keep overwriting.
        _set_ff_notes(payload.get('warnings') or [])

    def _stop_browser_field():
        """Take the browser's own field off the molecule.

        It has to be said to the page, not merely stopped being asked for: the
        relaxation runs on its own animation frames and keeps running until it
        is told, and the parameters stay installed until they are cleared --
        which is how a GFN method could be chosen while UFF went on relaxing
        underneath it.
        """
        _ensure_manip_bootstrap()
        _set_ff_notes([])
        _run_manip_js(
            'if(window.__delfinSubmitManip){'
            'window.__delfinSubmitManip.stopAutoOptimize('
            f'{json.dumps(submit_scope_id)});'
            'window.__delfinSubmitManip.setForceField('
            f'{json.dumps(submit_scope_id)},null);'
            '}'
        )

    def on_submit_relax_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_relax_btn.value)
        submit_relax_btn.button_style = 'info' if active else ''
        if _gfn.is_gfn_method(submit_ff_dd.value):
            # There is no GFN engine in the browser to run per frame, so this
            # switch means something else here: while it is on, the molecule
            # follows the atom being dragged -- one short xtb run per push,
            # and nothing at all when nothing is being dragged.  A continuous
            # loop of processes on a shared machine is what this deliberately
            # is not.
            _end_gfn_follow()
            # Whether it goes on or off, the browser's field has no business
            # here.  Left running from a UFF session it went on relaxing under
            # a molecule whose method box said GFN.
            _stop_browser_field()
            label = _gfn.GFN_METHODS[str(submit_ff_dd.value)]['label']
            if not active:
                state['gfn_settle_forced'] = False
                state['gfn_settle_rounds'] = 0
                _set_mol_status('The structure is no longer being relaxed.')
                return
            if _gfn.find_binary(submit_ff_dd.value) is None:
                _set_mol_status(f'{label} needs a program that was not found.')
                submit_relax_btn.value = False
                return
            if not submit_manip_btn.value:
                submit_manip_btn.value = True  # dragging is what it is for
            _ensure_manip_bootstrap()
            _install_gfn_frame_watcher()
            _set_mol_status(
                f'Drag an atom and the rest of the molecule follows it with '
                f'{label}. Letting go leaves it where you put it -- Optimise '
                'is what takes it downhill, and Settle asks for that on every '
                'release.'
                + ('' if submit_ff_dd.value == 'gfnff' else
                   f' {label} is the slow one -- if it drags heavily, '
                   'GFN-FF answers about twenty times faster.'))
            return
        if not active:
            _ensure_manip_bootstrap()
            _set_ff_notes([])
            _run_manip_js(
                'if(window.__delfinSubmitManip){'
                'window.__delfinSubmitManip.stopAutoOptimize('
                f'{json.dumps(submit_scope_id)});'
                'window.__delfinSubmitManip.setForceField('
                f'{json.dumps(submit_scope_id)},null);'
                '}'
            )
            return
        if not submit_manip_btn.value:
            submit_manip_btn.value = True   # dragging is what it is there for
        _enable_live_forcefield()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.startAutoOptimize('
            f'{json.dumps(submit_scope_id)});'
        )

    #: How long to wait after the last change before starting again.  Letting
    #: go of an atom arrives as a burst -- the release, then the coordinates,
    #: sometimes a settled version behind them -- and starting on the first of
    #: those would launch an xtb for each one.
    _GFN_RESTART_DELAY = 0.35

    def _interrupt_gfn():
        """End the running optimisation because the structure under it changed.

        xtb is minimising a geometry that stopped existing the moment an atom
        was moved, so the run is ended rather than raced.  It is not ended the
        way the switch ends it: nothing has been stopped from where the user
        is standing, and the frame they were shown is not a result to keep --
        the optimisation is about to start again from what they have made.
        """
        token = state.get('optimize_run')
        if token is None:
            return False
        state['optimize_run'] = None
        state['optimize_interrupted'] = token
        # No halt report: "stopped at frame 12" belongs to the switch.
        state['gfn_halt_sent'] = True
        # A run number the page has never seen, carrying nothing.  It resets
        # the player, so the frames of the abandoned run cannot play out over
        # the geometry the user has just made.
        blank = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = blank
        submit_gfn_frame.value = json.dumps({'run': blank, 'frames': []})
        return True

    def _arm_gfn_restart():
        """Start the optimisation again, once the changing has stopped."""
        if state.get('optimize_interrupted') is None:
            return
        state['gfn_restart_at'] = time.monotonic() + _GFN_RESTART_DELAY
        if state.get('gfn_restart_armed'):
            return                      # already waiting; it will wait longer
        state['gfn_restart_armed'] = True

        def _wait():
            while True:
                left = state.get('gfn_restart_at', 0.0) - time.monotonic()
                if left <= 0:
                    break
                time.sleep(min(left, 0.05))
            state['gfn_restart_armed'] = False
            _schedule_ui_update(_restart_gfn)

        threading.Thread(target=_wait, daemon=True).start()

    def _restart_gfn():
        if state.pop('optimize_interrupted', None) is None:
            return
        every_frame = bool(state.get('optimize_every_frame'))
        button = submit_optimize_all_btn if every_frame else submit_optimize_btn
        if not button.value:
            return                      # switched off while it was waiting
        state['gfn_restarting'] = True
        on_submit_optimize(None, every_frame=every_frame)

    def on_submit_optimize_all(change=None):
        on_submit_optimize(change, every_frame=True)

    def on_submit_optimize(change=None, every_frame=False):
        """A switch, not a push: on starts it, off stops it, and it turns
        itself off when the optimisation has converged or failed.

        *every_frame* is the difference between the two buttons: Optimize takes
        the frame on screen, all takes the whole set -- one after another,
        because a login node is shared.
        """
        button = submit_optimize_all_btn if every_frame else submit_optimize_btn
        if isinstance(change, dict) and change.get('name') == 'value':
            if not button.value:
                state['optimize_run'] = None      # off: the run ends itself
                return
        elif not button.value:
            return
        if submit_optimize_btn.value and submit_optimize_all_btn.value:
            # One run at a time, whichever was asked for second stands down.
            other = submit_optimize_btn if every_frame else submit_optimize_all_btn
            other.value = False
        """Minimise every frame that is loaded, not just the one on screen.

        The Submit tab can hold a whole set at once -- generated isomers, or
        the frames of a batch -- and any of them can end up submitted, so
        optimising only the visible one would leave the rest untouched.

        The geometries from before the run are kept so Undo can put them back:
        the browser's own undo stack cannot, because the results arrive from
        Python and re-render the viewer.
        """
        frames = list(state.get('isomers') or []) if every_frame else []
        single = (state.get('current_xyz_for_copy') or {}).get('content')
        if not frames and not single:
            _set_mol_status('Load a structure before optimising.')
            return
        method = submit_ff_dd.value
        gfn = _gfn.is_gfn_method(method)
        label = (_gfn.GFN_METHODS[method]['label'] if gfn else method.upper())
        charge = int(submit_gfn_charge.value or 0)
        # xtb counts unpaired electrons, not multiplicity: M = 2S+1.
        uhf = max(0, int(submit_gfn_mult.value or 1) - 1)
        autospin = bool(submit_gfn_autospin.value)
        count = len(frames) or 1
        # Which switch is running, so a restart presses the same one.
        state['optimize_every_frame'] = bool(every_frame)
        again = bool(state.pop('gfn_restarting', False))
        _set_mol_status(
            (f'Moved while it ran; {label} starts again from the structure '
             'you made...' if again else
             f'Optimising {count} frame(s) with {label}...'), spinner=True,
        )
        if gfn:
            # One call, not two: run_js clears its output before displaying,
            # so a bootstrap followed immediately by a watcher is a bootstrap
            # the page may never have run.
            _ensure_manip_bootstrap()
            _schedule_ui_update(_install_gfn_frame_watcher)
        played = [False]
        state['gfn_energy'] = None
        state['gfn_held'] = None
        state['gfn_halt_sent'] = False
        run_id = int(state.get('gfn_run', 0)) + 1
        state['gfn_run'] = run_id

        def _push_frames(frames):
            """Hand the path over while xtb is still walking it."""
            played[0] = True
            walked = list(frames)
            trail = walked[-400:]
            begins = len(walked) - len(trail)

            def _write(t=trail, first=begins):
                # A run that has been replaced does not draw.  An interrupted
                # one has frames in hand when it is told to stop, and writing
                # them afterwards played the abandoned path over the structure
                # the user had just made.
                if state.get('gfn_run') != run_id:
                    return
                submit_gfn_frame.value = json.dumps(
                    {'run': run_id, 'from': first, 'frames': t})

            _schedule_ui_update(_write)

        token = object()
        state['optimize_run'] = token

        def _stopped():
            halted = state.get('optimize_run') is not token
            if halted and not state.get('gfn_halt_sent'):
                state['gfn_halt_sent'] = True
                _schedule_ui_update(
                    lambda: setattr(submit_gfn_frame, 'value',
                                    json.dumps({'run': run_id, 'halt': 1,
                                                'frames': []})))
            return halted

        # The run this one replaces, taken before it is overwritten below.
        earlier = state.get('optimize_thread')
        # What the user is holding.  Set and Hold mean the same to xtb as they
        # mean to the browser's field -- a pull negotiates, a fix is met -- so
        # a value held on screen is held in the optimisation too, rather than
        # being quietly given up the moment GFN is chosen.
        held = list(state.get('constraints') or [])
        wet = str(submit_gfn_solvent.value or '') or None

        def _work():
            from .molecule_forcefield import relax_xyz
            # One xtb at a time.  An interrupted run is still shutting its
            # process down when the replacement starts, and a login node is
            # shared -- so the new one waits for the old one to be gone rather
            # than briefly doubling up.
            if earlier is not None and earlier.is_alive():
                earlier.join(timeout=15)
            results, failures = [], []
            targets = frames or [(single, None, None)]
            for position, item in enumerate(targets):
                xyz = item[0]
                try:
                    if _stopped():
                        failures.append(f'frame {position + 1}: stopped')
                        results.append(item)
                        continue
                    # The bonding the editor has been working with, so that
                    # pressing Optimise on a structure that has been pulled
                    # about does not re-perceive it and shove the molecule
                    # apart -- the same cliff the drag had.  Only for the one
                    # frame on screen: a set of isomers is a set of different
                    # molecules, and one perception cannot describe them all.
                    perceived = (None if every_frame
                                 else _gfn_topology_dir(
                                     len(_gfn.atom_lines(xyz))))
                    if gfn and autospin:
                        outcome = _gfn.optimize_autospin(
                            xyz, method, charge=charge, should_stop=_stopped,
                            timeout=None, on_frames=_push_frames,
                            constraints=held, topology=perceived, solvent=wet)
                    elif gfn:
                        outcome = _gfn.optimize_with_gfn(
                            xyz, method, charge=charge, uhf=uhf,
                            should_stop=_stopped, timeout=None,
                            constraints=held, topology=perceived, solvent=wet,
                            on_frames=_push_frames if position == 0 else None)
                    else:
                        outcome = relax_xyz(
                            xyz,
                            max_steps=500,
                            perceived=_perception_for(xyz),
                            method=method,
                        )
                except Exception as exc:
                    failures.append(f'frame {position + 1}: {exc}')
                    results.append(item)
                    continue
                if outcome.get('ok'):
                    if outcome.get('energy') is not None and position == 0:
                        state['gfn_energy'] = float(outcome['energy'])
                    if position == 0:
                        state['gfn_held'] = outcome.get('held')
                        if outcome.get('multiplicity'):
                            # What the scan settled on, so the live relaxation
                            # asks the same question this one did rather than
                            # a different one with the box's fixed M.
                            state['gfn_scanned_uhf'] = int(outcome['uhf'])
                    kept = outcome['xyz']
                    if _stopped() and outcome.get('frames'):
                        # Stop means the frame that was on screen.  xtb runs
                        # ahead of the picture, and geometries nobody saw are
                        # not what the user stopped at.
                        shown = state.get('gfn_shown_frame')
                        walked = outcome['frames']
                        if isinstance(shown, int) and 0 < shown <= len(walked):
                            symbols = [line.split()[0]
                                       for line in _gfn.atom_lines(xyz)]
                            frame = walked[shown - 1]
                            if len(symbols) * 3 == len(frame):
                                rows = [
                                    f'{symbols[i]} {frame[3*i]:.8f} '
                                    f'{frame[3*i+1]:.8f} {frame[3*i+2]:.8f}'
                                    for i in range(len(symbols))
                                ]
                                kept = (f'{len(rows)}\nstopped at the frame on '
                                        f'screen\n' + '\n'.join(rows) + '\n')
                    results.append((kept,) + tuple(item[1:]))
                    if gfn and outcome.get('frames') and position == 0:
                        played[0] = True
                        # xtb writes every cycle to xtbopt.log, so the path
                        # costs nothing extra -- one run, and the viewer plays
                        # what the optimiser really walked through.
                        _push_frames(outcome['frames'])
                    note = str(outcome.get('status') or '')
                    if 'before converging' in note:
                        # It came back with a geometry, but not a finished one.
                        failures.append(f'frame {position + 1}: {note}')
                else:
                    failures.append(
                        f"frame {position + 1}: {outcome.get('status') or 'failed'}"
                    )
                    results.append(item)

            def _apply():
                if state.get('optimize_interrupted') is token:
                    # The structure changed under this run.  What it reached is
                    # a minimum of a geometry that no longer exists, so neither
                    # the coordinates nor the switch are touched: the run that
                    # replaces it is already on its way.
                    return
                # Converged, failed or stopped -- the switch goes back up by
                # itself, so it never claims to be working when it is not.
                if state.get('optimize_run') is token:
                    state['optimize_run'] = None
                for switch in (submit_optimize_btn, submit_optimize_all_btn):
                    if switch.value:
                        switch.value = False
                    switch.disabled = False
                state['pre_optimize_frames'] = {
                    'isomers': frames,
                    'coords': coords_widget.value,
                }
                if frames:
                    state['isomers'] = results
                    if not played[0]:
                        # Showing the isomer rebuilds the viewer, which is the
                        # other way the playback was being torn down.
                        _show_isomer_at_index(state.get('isomer_index', 0))
                else:
                    lines = [
                        line for line in results[0][0].splitlines()[2:] if line.strip()
                    ]
                    if played[0]:
                        # The trajectory is playing, and its last frame is this
                        # very geometry.  Writing the coordinates the ordinary
                        # way rebuilds the viewer, which tore the playback down
                        # milliseconds after it started -- so only the end of
                        # the optimisation was ever seen.  The box is updated
                        # for Copy and Submit; the picture is already right.
                        state['manip_inflight'] = True
                    coords_widget.value = (
                        f"{len(lines)}\nOptimised in DELFIN viewer\n"
                        + '\n'.join(lines)
                    )
                done = count - len(failures)
                said = f'Optimised {done} of {count} frame(s) with {label}.'
                # The energy, the way the force field shows one.  xtb reports
                # it in hartree; kcal/mol is what the rest of the tab speaks,
                # and both are given because a total energy is compared
                # against other totals and a difference against chemistry.
                energy = state.get('gfn_energy')
                if energy is not None:
                    said += (f' E = {energy:.6f} Eh '
                             f'({energy * 627.5094740631:.2f} kcal/mol).')
                # What was held while it ran, and what became of it.  A value
                # the user is holding on screen that the optimisation quietly
                # ignored would make the result an answer to a question nobody
                # asked.
                said += _gfn.solvent_note(submit_gfn_solvent.value)
                said += _gfn.held_note(state.get('gfn_held') or {
                    'held': 0, 'dropped': [], 'mixed': False, 'force': None})
                state['gfn_last_status'] = said
                if gfn:
                    if autospin:
                        picked = results[0][0] if results else None
                        del picked
                        said += f' charge {charge}, multiplicity scanned.'
                    else:
                        said += f' charge {charge}, multiplicity {uhf + 1}.'
                # The page's report about the playback is kept: it arrives
                # while this message is being built and would otherwise be
                # wiped by it a moment later.
                note = state.get('gfn_play_note')
                extra = [f'Trajectory: {note}.'] if note else []
                _set_mol_status(said, *(list(failures[:2]) + extra))

            _schedule_ui_update(_apply)

        worker = threading.Thread(target=_work, daemon=True)
        state['optimize_thread'] = worker
        worker.start()

    def _clear_selection():
        """Drop the picks so the next constraint starts from a clean set."""
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.clearSelection('
            f'{json.dumps(submit_scope_id)});'
        )

    def _step_for_selection(indices):
        """How far one press of an arrow key moves the value.

        Three different quantities, three different steps.  A bond length is
        Angstrom, where a hundredth is fine.  An angle is degrees and a tenth
        is the useful step.  A dihedral is what gets turned through a whole
        rotation -- half a degree there, so holding the key sweeps it instead
        of creeping.
        """
        picked = len(indices or ())
        if picked == 2:
            submit_internal_value.step = 0.01     # bond length, in Angstrom
        elif picked == 4:
            submit_internal_value.step = 0.5      # dihedral, swept by hand
        else:
            submit_internal_value.step = 0.1      # angle

    def _apply_internal_now():
        """Put the selection at the value in the box, and leave it selected."""
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setInternal('
            + json.dumps(submit_scope_id) + ','
            + repr(float(submit_internal_value.value)) + ');'
        )
        # The browser moves the atoms; under GFN the rest of the molecule has
        # to be given the chance to arrange itself around where they landed.
        _arm_gfn_takeup('Set')

    def on_submit_set_internal(change=None):
        """Set is a mode: while it is on, the box turns the selection by hand.

        Switching it on puts the selection at what the box says, and every
        further change of the box does the same -- so the arrow keys turn a
        dihedral a tenth of a degree at a time and the structure follows.  The
        picks are kept, because letting go of them after every step is what
        made turning something by hand impossible.
        """
        if not submit_internal_btn.value:
            return
        _apply_internal_now()

    def on_submit_internal_value(change):
        """The box changed.  Who owns it depends on what is selected."""
        if change.get('name') != 'value' or state.get('internal_quiet'):
            return
        if _selected_constraint()[1] is not None:
            return          # a held value is being retuned, not the geometry
        if submit_internal_btn.value:
            _apply_internal_now()

    def on_submit_pick_sync(change):
        """Offer the coordination polyhedra of a metal the moment it is picked.

        Which ones are possible follows from its coordination number, and the
        tables are MANTA's own -- the same ideal donor vectors it builds
        complexes with.
        """
        if change.get('name') != 'value':
            return
        raw = (submit_pick_sync.value or '').strip()
        indices = [int(part) for part in raw.split(',') if part.strip().isdigit()]
        state['picked'] = indices
        _step_for_selection(indices)
        _refresh_swap(indices)
        xyz = (state.get('current_xyz_for_copy') or {}).get('content')
        options = None
        perceived = None
        if xyz and indices:
            # Perceive on demand. Waiting for the cache meant the offer only
            # appeared once the force field had been switched on at least
            # once -- tapping a metal before that did nothing at all, and said
            # nothing either.
            try:
                perceived = _perception_for(xyz)
            except Exception:
                perceived = None
        if perceived is not None and len(indices) == 1:
            try:
                from .molecule_forcefield import polyhedron_options
                options = polyhedron_options(perceived, indices[0])
            except Exception:
                options = None
        _refresh_hybridisation(indices, perceived)
        if not options:
            # Only the offer follows the selection. The applied polyhedron
            # stays: clearing it here meant that selecting three ligand atoms
            # to hold an angle silently threw the polyhedron away, and the very
            # next export went out without it.
            submit_poly_dd.layout.display = 'none'
            submit_poly_dd.disabled = True
            state['poly_offer_metal'] = None
            # Say why, when a single atom was picked and could have qualified.
            if perceived is not None and len(indices) == 1:
                index = indices[0]
                symbol = perceived.symbols[index]
                if index in set(perceived.metal_indices or ()):
                    donors = sorted(
                        j for pair in perceived.bonds for j in pair
                        if index in pair and j != index
                    )
                    _set_mol_status(
                        f'{symbol}: coordination number {len(donors)} — no '
                        'polyhedron table for that (2 to 9 are covered).'
                    )
            return

        coordination, current, choices = options
        state['poly_offer_metal'] = indices[0]
        symbol = perceived.symbols[indices[0]]
        entries = [(f'{symbol} · CN {coordination} · free', '')]
        for code, label in choices:
            mark = ' (current)' if code == current else ''
            entries.append((f'{label}{mark}', code))
        state['poly_quiet'] = True
        try:
            submit_poly_dd.options = entries
            # Only a code this metal actually offers. A polyhedron held on one
            # metal has nothing to say about the next one picked, and assigning
            # a value the options do not contain raises.
            applied = state.get('poly_applied')
            offered = {code for code, _label in choices}
            submit_poly_dd.value = applied if applied in offered else ''
        finally:
            state['poly_quiet'] = False
        submit_poly_dd.disabled = False
        submit_poly_dd.layout.display = ''

    def _refresh_hybridisation(indices, perceived):
        """Offer to overrule the hybridisation of the picked atoms.

        Any number of them: a double bond that went unperceived usually cost
        both of its carbons their type, and retyping a ring one atom at a time
        is busywork. Metals are dropped from the selection rather than
        blocking it -- RDKit's UFF has no types for one at all, so its bonds
        and angles come from the geometry either way.
        """
        metals = set(perceived.metal_indices or ()) if perceived else set()
        # An index the structure no longer has: the browser pushes its picks
        # after a re-render, and an edit that deleted atoms renumbers them.
        # Asking the perception about one is an IndexError, and this handler
        # runs on every click in the viewer.
        chosen = [
            i for i in indices
            if i not in metals and 0 <= i < len(perceived.symbols)
        ] if perceived else []
        if not chosen:
            submit_hyb_dd.layout.display = 'none'
            submit_hyb_dd.disabled = True
            state['hyb_offer_atoms'] = []
            return

        from .molecule_forcefield import (
            HYBRIDISATION_CHOICES, perceived_hybridisation_of,
        )

        overrides = state.get('hyb_overrides') or {}
        auto = {perceived_hybridisation_of(perceived, i) for i in chosen}
        if len(chosen) == 1:
            index = chosen[0]
            head = (f'{perceived.symbols[index]}{index} · '
                    f'{auto.pop() or "no type"} (automatic)')
        else:
            named = ''.join(sorted({perceived.symbols[i] for i in chosen}))
            found = auto.pop() if len(auto) == 1 else 'mixed'
            head = f'{len(chosen)} atoms ({named}) · {found or "no type"} (automatic)'
        entries = [(head, '')]
        for name in HYBRIDISATION_CHOICES:
            entries.append((f'force {name}', name))
        # Only a value they all already share can be shown as the current one.
        held = {overrides.get(i, '') for i in chosen}
        state['hyb_offer_atoms'] = chosen
        state['hyb_quiet'] = True
        try:
            submit_hyb_dd.options = entries
            submit_hyb_dd.value = held.pop() if len(held) == 1 else ''
        finally:
            state['hyb_quiet'] = False
        submit_hyb_dd.disabled = False
        submit_hyb_dd.layout.display = ''

    def on_submit_hyb_changed(change):
        if change.get('name') != 'value' or state.get('hyb_quiet'):
            return
        atoms = list(state.get('hyb_offer_atoms') or [])
        if not atoms:
            return
        index = atoms[0]
        overrides = dict(state.get('hyb_overrides') or {})
        chosen = submit_hyb_dd.value or ''
        for atom in atoms:
            if chosen:
                overrides[atom] = chosen
            else:
                overrides.pop(atom, None)
        state['hyb_overrides'] = overrides
        # Perception is cached by element sequence, which this does not
        # change, so the cache has to be dropped or the override never reaches
        # the force field.
        state['perceived'] = None
        state['perceived_for'] = None
        state['poly_assignment'] = None
        _enable_live_forcefield()
        perceived = state.get('perceived')
        if len(atoms) == 1:
            symbol = perceived.symbols[index] if perceived else '?'
            named = f'{symbol}{index}'
        else:
            named = f'{len(atoms)} atoms'
        if chosen:
            shape = {'sp': 'linear', 'sp2': 'trigonal planar',
                     'sp3': 'tetrahedral'}[chosen]
            _set_mol_status(f'{named} typed as {chosen}: {shape}.')
        else:
            _set_mol_status(f'{named} back to the perceived type.')
        _clear_selection()

    def _refresh_poly_turn():
        """Offer Turn only where the vertices are not all alike.

        An octahedron has nothing to turn -- every vertex is the same, and
        which ligand is trans to which is what Swap is for. A trigonal
        bipyramid has two kinds, so which pair is axial is a real choice.
        """
        geometry = state.get('poly_applied')
        metal = state.get('poly_metal')
        perceived = state.get('perceived')
        turnable = False
        if geometry and metal is not None and perceived is not None:
            try:
                from .molecule_forcefield import polyhedron_vertex_classes
                donors = len(perceived.neighbours()[int(metal)])
                grouped = polyhedron_vertex_classes(donors, geometry)
                turnable = bool(grouped) and len(set(grouped[0])) > 1
            except Exception:
                turnable = False
        submit_poly_turn_btn.layout.display = '' if turnable else 'none'
        submit_poly_turn_btn.disabled = not turnable
        if not turnable:
            state['poly_arrangements'] = []
            state['poly_arrangement_index'] = 0

    def on_submit_poly_turn(_button=None):
        """Step to the next way the ligands can sit on this polyhedron."""
        geometry = state.get('poly_applied')
        metal = state.get('poly_metal')
        xyz = (state.get('current_xyz_for_copy') or {}).get('content')
        if not geometry or metal is None or not xyz:
            _set_mol_status('Choose a polyhedron for a metal first.')
            return
        try:
            from .molecule_forcefield import (
                describe_polyhedron_arrangement, parse_xyz,
                polyhedron_arrangements,
            )
            perceived = _perception_for(xyz)
            # The coordinates as they are now, not as they were perceived: a
            # ligand that has been dragged has to be scored where it sits.
            parsed = parse_xyz(xyz)
            coords = perceived.coords
            if parsed is not None and list(parsed[0]) == list(perceived.symbols):
                coords = parsed[1]
            arrangements = polyhedron_arrangements(
                perceived, int(metal), geometry, coords)
        except Exception as exc:
            _set_mol_status(f'Could not work out the arrangements: {exc}')
            return
        if len(arrangements) < 2:
            _set_mol_status(
                'Every vertex of this polyhedron is the same, so there is '
                'nothing to turn. Swap exchanges two ligands.'
            )
            return
        position = (int(state.get('poly_arrangement_index') or 0) + 1) % len(arrangements)
        state['poly_arrangements'] = arrangements
        state['poly_arrangement_index'] = position
        state['poly_assignment'] = arrangements[position]
        _enable_live_forcefield()
        described = describe_polyhedron_arrangement(
            perceived, geometry, arrangements[position])
        _set_mol_status(
            f'Arrangement {position + 1} of {len(arrangements)} — {described}.'
        )

    def on_submit_hyb_auto(_button=None):
        """Derive the carbon types from the connectivity and hold them."""
        xyz = (state.get('current_xyz_for_copy') or {}).get('content')
        if not xyz:
            _set_mol_status('Load a structure first.')
            return
        try:
            from .molecule_forcefield import hybridisation_from_connectivity
            perceived = _perception_for(xyz)
            picked = list(state.get('picked') or [])
            derived = hybridisation_from_connectivity(perceived, picked or None)
        except Exception as exc:
            _set_mol_status(f'Could not read the types: {exc}')
            return
        if not derived:
            _set_mol_status(
                'No carbon in the selection — the count only fixes the shape '
                'for carbon, which has no lone pair.'
            )
            return
        changed = [
            i for i, name in derived.items()
            if (state.get('hyb_overrides') or {}).get(i) != name
        ]
        overrides = dict(state.get('hyb_overrides') or {})
        overrides.update(derived)
        state['hyb_overrides'] = overrides
        state['perceived'] = None
        state['perceived_for'] = None
        state['poly_assignment'] = None
        _enable_live_forcefield()
        where = f'{len(picked)} selected' if picked else 'the whole structure'
        counts = ', '.join(
            f'{sum(1 for v in derived.values() if v == name)}x {name}'
            for name in ('sp', 'sp2', 'sp3')
            if any(v == name for v in derived.values())
        )
        _set_mol_status(
            f'{len(derived)} carbons typed from their partners in {where} '
            f'({counts}); {len(changed)} changed.'
        )
        _clear_selection()

    #: How many structural edits can be taken back.  The browser keeps 50
    #: coordinate snapshots; there is no reason for this to be shorter.
    _STRUCTURE_UNDO_LIMIT = 50

    _CONSTRAINT_KINDS = {2: 'distance', 3: 'angle', 4: 'dihedral'}

    def _describe_constraint(entry):
        symbols = []
        perceived = state.get('perceived')
        for index in entry['atoms']:
            symbol = perceived.symbols[index] if perceived else '?'
            symbols.append(f'{symbol}{index}')
        unit = 'A' if entry['kind'] == 'distance' else 'deg'
        mode = entry.get('mode', 'pull')
        return f"{'-'.join(symbols)} = {entry['value']:.3g} {unit} ({mode})"

    def _selected_constraint():
        """The held entry the list is pointing at, or (None, None)."""
        key = submit_constraint_dd.value or ''
        if not key.startswith('c'):
            return None, None
        held = list(state.get('constraints') or [])
        position = int(key[1:])
        if not (0 <= position < len(held)):
            return None, None
        return position, held[position]

    def _sync_constraint_selection():
        """Mark the atoms of the selected entry and offer its value for editing.

        With nothing selected the picture keeps whatever the user picked and the
        value box goes back to serving the selection, which is what it does when
        no held value is being looked at.
        """
        _position, entry = _selected_constraint()
        if entry is None:
            return
        state['hold_mode_quiet'] = True
        try:
            submit_hold_mode.value = entry.get('mode', 'pull')
        finally:
            state['hold_mode_quiet'] = False
        kind = str(entry['kind'])
        unit = 'Å' if kind == 'distance' else '°'
        number = float(entry['value'])
        value = '{:.3f}'.format(number) if kind == 'distance' else '{:.1f}'.format(number)
        label = 'held <b>{}</b> ({})'.format(kind, unit)
        atoms = [int(i) for i in entry['atoms']]
        _run_manip_js(
            'if(window.__delfinSubmitManip&&window.__delfinSubmitManip.setPicks)'
            'window.__delfinSubmitManip.setPicks('
            + json.dumps(submit_scope_id) + ','
            + json.dumps(atoms) + ','
            + json.dumps(value) + ','
            + json.dumps(label) + ');'
        )

    def on_submit_constraint_selected(change):
        if change.get('name') != 'value' or state.get('constraint_quiet'):
            return
        if not (submit_constraint_dd.value or ''):
            _clear_selection()
            return
        _sync_constraint_selection()

    def on_submit_constraint_retune(change):
        """Editing the value while a held entry is selected retunes that entry."""
        if change.get('name') != 'value':
            return
        position, entry = _selected_constraint()
        if entry is None:
            return
        if abs(float(submit_internal_value.value) - float(entry['value'])) < 1e-9:
            return  # the box was filled with the held value, not edited
        held = list(state.get('constraints') or [])
        held[position] = dict(entry, value=float(submit_internal_value.value))
        state['constraints'] = held
        _refresh_constraints()
        _enable_live_forcefield()
        _set_mol_status(f'Holding {_describe_constraint(held[position])}.')

    def _refresh_constraints():
        """Show what the field is currently being held to."""
        entries = []
        if state.get('poly_applied') and state.get('poly_metal') is not None:
            entries.append(('poly', f"polyhedron: {state['poly_applied']}"))
        for position, entry in enumerate(state.get('constraints') or []):
            entries.append((f'c{position}', _describe_constraint(entry)))
        visible = bool(entries)
        submit_constraint_dd.layout.display = '' if visible else 'none'
        submit_constraint_del.layout.display = '' if visible else 'none'
        submit_constraint_dd.disabled = not visible
        submit_constraint_del.disabled = not visible
        if not visible:
            # Leave nothing behind: a hidden dropdown that still holds the
            # last constraint shows it again the moment another one is set.
            state['constraint_quiet'] = True
            try:
                submit_constraint_dd.options = [('no constraints', '')]
                submit_constraint_dd.value = ''
            except Exception:
                pass
            finally:
                state['constraint_quiet'] = False
            return
        # Nothing is selected to begin with.  Selecting an entry marks the atoms
        # it holds, so a preselected one would mean the picture always shows a
        # marked set nobody asked for.
        previous = submit_constraint_dd.value
        state['constraint_quiet'] = True
        try:
            submit_constraint_dd.options = (
                [(f'{len(entries)} held · show which', '')]
                + [(label, key) for key, label in entries]
            )
            submit_constraint_dd.value = (
                previous if previous in dict(submit_constraint_dd.options).values() else ''
            )
        finally:
            state['constraint_quiet'] = False
        _sync_constraint_selection()

    def _refresh_swap(indices):
        """Offer an exchange whenever two donors of one metal are selected.

        It does not need a polyhedron: exchanging two ligands is a move across
        a barrier, which is useful with or without a target shape.
        """
        perceived = state.get('perceived')
        metals = set(getattr(perceived, 'metal_indices', ()) or ())
        donors = set()
        if perceived is not None and len(metals) == 1:
            metal = next(iter(metals))
            donors = {
                j for pair in perceived.bonds for j in pair
                if metal in pair and j != metal
            }
        ready = len(indices) == 2 and donors and all(i in donors for i in indices)
        submit_swap_btn.layout.display = '' if ready else 'none'
        submit_swap_btn.disabled = not ready

    def _edit_bond(connect):
        """Draw or remove a bond, and remember it.

        Bond perception reads distances, and in a crowded coordination sphere
        that is simply not reliable: on a real Pt complex it counted two ipso
        carbons of a phosphine's phenyls as donors, giving CN 6 for a
        four-coordinate metal, while the viewer's own perception invented a
        Pt-H bond instead. Neither is trustworthy, so the correction has to be
        remembered and re-applied -- otherwise the next perception, which runs
        from the geometry, would quietly undo it.
        """
        indices = list(state.get('picked') or [])
        if len(indices) != 2:
            _set_mol_status('Select exactly two atoms to change a bond.')
            return
        pair = (min(indices), max(indices))
        edits = {tuple(k): v for k, v in (state.get('bond_edits') or {}).items()}
        edits[pair] = bool(connect)
        state['bond_edits'] = edits
        # The perception is cached by element sequence, which a bond edit does
        # not change -- so the cache has to be dropped explicitly, or the
        # correction would never reach the force field at all.
        state['perceived'] = None
        state['perceived_for'] = None
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.editBond('
            f'{json.dumps(submit_scope_id)},{pair[0]},{pair[1]},'
            f'{"true" if connect else "false"});'
        )
        # Re-assign the parameters straight away: the topology decides the
        # bonds, angles and torsions the field works with, and until it is
        # rebuilt the relaxation is still holding the bond that was just cut.
        state['poly_assignment'] = None
        _enable_live_forcefield()
        verb = 'Bonded' if connect else 'Unbonded'
        _set_mol_status(f'{verb} atoms {pair[0]} and {pair[1]}.')
        _clear_selection()

    def _apply_structure(structure, note):
        """Put an edited structure back and let the tab rebuild around it.

        The coordinate box is the tab's single source of truth, so writing to
        it re-renders the viewer and re-perceives everything. What does *not*
        survive that is the bond orders: perception reads them off the
        geometry and does not get them back -- ethene built here came back as
        a single bond at 1.514 A. So the edited topology is re-seeded as the
        hand corrections it is, which is exactly what it is: the user built
        this, and their bonds outrank a distance table until a different
        structure is loaded.
        """
        from .molecule_builder import to_xyz

        # Remember what it looked like first: a structural edit changes the
        # atom count, which the browser's coordinate snapshots cannot express.
        history = list(state.get('structure_undo') or [])
        history.append({
            'coords': coords_widget.value,
            'bond_edits': dict(state.get('bond_edits') or {}),
        })
        state['structure_undo'] = history[-_STRUCTURE_UNDO_LIMIT:]

        # An atom added or taken away is a different molecule from the one xtb
        # has in hand, down to the number of coordinates -- so a run under it
        # is ended here and starts again on the molecule that now exists, and
        # the bonding GFN-FF had perceived goes with it.
        _drop_gfn_topology()
        _interrupt_gfn()
        _arm_gfn_restart()

        xyz = to_xyz(structure, note)
        lines = [line for line in xyz.splitlines()[2:] if line.strip()]
        # The write below re-renders through update_molecule_view, which
        # clears the history a new structure invalidates. This is not a new
        # structure, it is a step in the one being edited.
        state['structure_edit_inflight'] = True
        _mark_structure_edit()
        try:
            coords_widget.value = f'{len(lines)}\n{note}\n' + '\n'.join(lines)
        finally:
            state['structure_edit_inflight'] = False
        # After update_molecule_view, which clears them.
        state['bond_edits'] = {
            (int(i), int(j)): int(order)
            for (i, j), order in structure.bonds.items()
        }
        state['perceived'] = None
        state['perceived_for'] = None
        # Types, derived and held, after every build step -- which is what
        # pressing the button by hand after one was doing. Perception reads
        # bond orders off the geometry, and a structure that has just been
        # built is exactly where that geometry is least settled: a centre
        # comes back sp3 and its angles at 109.5 where they should be 120, so
        # the field pulls the new part into the wrong shape. The number of
        # partners says it outright and does not depend on the geometry at
        # all, which is why doing it by hand fixed things.
        # Everything the tab holds by index has to follow the renumbering an
        # edit causes -- a deleted hydrogen moves every atom after it. A held
        # value that quietly pointed at different atoms, or vanished, is worse
        # than one that is dropped and said so.
        renumber = structure.renumbering()

        def _follow(indices):
            moved = [renumber.get(int(i)) for i in indices]
            return None if any(x is None for x in moved) else moved

        kept, lost = [], 0
        for entry in (state.get('constraints') or []):
            moved = _follow(entry.get('atoms') or [])
            if moved is None:
                lost += 1
                continue
            kept.append(dict(entry, atoms=moved))
        state['constraints'] = kept

        state['hyb_overrides'] = {
            renumber[i]: name
            for i, name in (state.get('hyb_overrides') or {}).items()
            if i in renumber
        }
        metal = state.get('poly_metal')
        if metal is not None:
            if metal in renumber:
                state['poly_metal'] = renumber[metal]
                assignment = state.get('poly_assignment') or {}
                followed = {
                    renumber[d]: v for d, v in assignment.items() if d in renumber
                }
                state['poly_assignment'] = (
                    followed if len(followed) == len(assignment) else None)
            else:
                state['poly_applied'] = None
                state['poly_metal'] = None
                state['poly_assignment'] = None

        # Straight off the structure that was just built, not off a fresh
        # perception of it: perception can fail or come back empty at exactly
        # this moment, and then nothing was derived at all -- which is why the
        # button still had to be pressed by hand. The builder knows every bond
        # it made, so the count is certain.
        by_count = {2: 'sp', 3: 'sp2', 4: 'sp3'}
        derived = {}
        for index, symbol in enumerate(structure.symbols):
            if symbol != 'C':
                continue
            partners = structure.neighbours(index)
            # A side-on alkene is the one case where the metal does not count.
            side_on = [
                m for m in partners
                if structure.symbols[m] not in ('H', 'C', 'N', 'O', 'S', 'P')
                and any(o in structure.neighbours(m) for o in partners
                        if structure.symbols[o] == 'C')
            ]
            name = by_count.get(len(partners) - len(side_on))
            if name:
                derived[index] = name
        if derived:
            overrides = dict(state.get('hyb_overrides') or {})
            overrides.update(derived)
            state['hyb_overrides'] = overrides
            state['perceived'] = None
            state['perceived_for'] = None
        if lost:
            note = f'{note} {lost} held value(s) lost their atoms and were dropped.'
        _set_mol_status(note)
        _refresh_constraints()
        _push_bond_orders(structure.bonds)
        if submit_relax_btn.value or state.get('ff_bootstrap_done'):
            _enable_live_forcefield()

    def _undo_structure():
        """Take back the last structural edit.

        Reached when the browser's own stack is empty, which after any
        structural edit it always is: the re-render clears it. So the two
        stacks stay in order without either side keeping a clock.
        """
        history = list(state.get('structure_undo') or [])
        if not history:
            _set_mol_status('Nothing left to undo.')
            return
        snapshot = history.pop()
        state['structure_edit_inflight'] = True
        _mark_structure_edit()
        try:
            coords_widget.value = snapshot['coords']
        finally:
            state['structure_edit_inflight'] = False
        state['structure_undo'] = history
        state['bond_edits'] = dict(snapshot.get('bond_edits') or {})
        state['perceived'] = None
        state['perceived_for'] = None
        _set_mol_status('Took back the last structural edit.')
        if submit_relax_btn.value or state.get('ff_bootstrap_done'):
            _enable_live_forcefield()

    def _push_bond_orders(bonds=None):
        """Let the picture show what the bonds are.

        A model read from an XYZ block has no orders in it -- the format
        carries none -- so every bond was drawn as one stick whatever it was.
        3Dmol draws a double as two cylinders and a triple as three once the
        model knows, so the orders are handed over after every render.
        """
        if bonds is None:
            xyz = (state.get('current_xyz_for_copy') or {}).get('content')
            if not xyz:
                return
            try:
                perceived = _perception_for(xyz)
                from .molecule_forcefield import _orders_from_mol
                known = _orders_from_mol(perceived.typing_mol)
                bonds = {
                    pair: int(known.get(pair, 1)) for pair in perceived.bonds
                }
            except Exception:
                return
        triples = [
            [int(i), int(j), int(order)]
            for (i, j), order in dict(bonds).items() if int(order) > 1
        ]
        if not triples:
            return
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setBondOrders('
            f'{json.dumps(submit_scope_id)},{json.dumps(triples)});'
        )

    def _mark_structure_edit():
        """Tell the re-render that this is an edit, not a different molecule.

        Two things follow from it: the camera stays where the user put it, and
        the continuous relaxation picks up again with the atom that was just
        drawn in it -- which is the point of being able to draw while it runs.
        """
        _ensure_manip_bootstrap()
        _run_manip_js('window.__delfinStructureEdit = true;')

    def _structure_now():
        from .molecule_builder import structure_from_xyz

        xyz = (state.get('current_xyz_for_copy') or {}).get('content')
        if not xyz:
            return None
        return structure_from_xyz(xyz, state.get('bond_edits') or {})

    def on_submit_cmd(change):
        """Carry out a gesture the browser cannot finish on its own.

        The value is ``verb:serial:payload``; the serial only exists to make
        the same command twice in a row read as two changes. Placing an atom
        is cheap, but how many hydrogens it needs and where they go is decided
        here, where RDKit's valences and covalent radii are.
        """
        if change.get('name') != 'value':
            return
        parts = (submit_cmd_sync.value or '').strip().split(':')
        if len(parts) != 3:
            return
        verb, payload = parts[0], parts[2]

        if verb == 'gfnplay':
            # What the playback is doing, said by the page.  Without this the
            # only way to tell an invisible trajectory from a missing one was
            # to read the browser's console.
            state['gfn_play_note'] = str(payload)
            if str(payload).startswith('stopped at frame '):
                try:
                    state['gfn_shown_frame'] = int(str(payload).rsplit(' ', 1)[1])
                except ValueError:
                    pass
            _set_mol_status(*[line for line in (
                state.get('gfn_last_status') or '', f'Trajectory: {payload}.'
            ) if line])
            return

        if verb == 'gfngrab':
            # An atom has been picked up while xtb was minimising.  Ending the
            # run at the grab rather than at the release is the whole point: a
            # GFN2 run is thirteen seconds, and all of them would otherwise be
            # spent on a structure the user is in the middle of changing.
            if _interrupt_gfn():
                _set_mol_status('Moved while it ran; the optimisation stops '
                                'there and starts again from what you make.',
                                spinner=True)
            _begin_gfn_follow()
            return

        if verb == 'gfnfree':
            _end_gfn_follow()
            # Let go of.  If nothing else has already asked for the restart --
            # a drag that moved something sends its coordinates first -- this
            # is what keeps the switch from being left on with nothing running.
            _arm_gfn_restart()
            # Letting go leaves the structure where the hand left it.  Dragging
            # moves along the potential surface; going down it is what Optimise
            # is for, and doing that unasked took the choice away -- the atom
            # you had just placed was pulled off the place you placed it.
            # Settle is the way to ask for it on every release.
            _arm_gfn_settle()
            return

        if verb == 'undo':
            _undo_structure()
            return

        if verb == 'unbond':
            indices = [int(p) for p in payload.split(',') if p.strip().isdigit()]
            if len(indices) == 2:
                state['picked'] = sorted(indices)
                _edit_bond(False)
            return

        from .molecule_builder import (
            delete_atoms, grow_from, normalise_element, place_atom,
            set_bond_order, set_element,
        )

        structure = _structure_now()
        if structure is None:
            _set_mol_status('Load a structure first.')
            return
        fields = payload.split(',')
        try:
            if verb == 'addatom' and len(fields) == 4:
                element = normalise_element(fields[0])
                if element is None:
                    _set_mol_status(f'{fields[0]} is not an element.')
                    return
                place_atom(structure, element,
                           [float(v) for v in fields[1:4]])
                _apply_structure(structure, f'Placed {element}.')
            elif verb == 'grow' and len(fields) == 6:
                element = normalise_element(fields[1])
                if element is None:
                    _set_mol_status(f'{fields[1]} is not an element.')
                    return
                grow_from(structure, int(fields[0]), element,
                          order=int(fields[2]),
                          direction=[float(v) for v in fields[3:6]])
                _apply_structure(structure, f'Grew {element}.')
            elif verb == 'setelement' and len(fields) == 2:
                element = normalise_element(fields[1])
                if element is None:
                    _set_mol_status(f'{fields[1]} is not an element.')
                    return
                index = int(fields[0])
                was = structure.symbols[index] if index < len(structure) else '?'
                if set_element(structure, index, element):
                    _apply_structure(structure, f'{was}{index} is now {element}.')
            elif verb == 'bondcycle' and len(fields) == 2:
                first, second = (int(v) for v in fields)
                ends = [structure.symbols[i] for i in (first, second)]
                current = structure.order(first, second)
                stepped = None
                for step in (1, 2, 3):
                    candidate = (current - 1 + step) % 3 + 1
                    if candidate == current:
                        break
                    if set_bond_order(structure, first, second, candidate):
                        stepped = candidate
                        break
                if stepped is None:
                    _set_mol_status(
                        f'{ends[0]}{first}-{ends[1]}{second} can only be '
                        'single: neither end has valence for more.'
                    )
                else:
                    named = {1: 'single', 2: 'double', 3: 'triple'}[stepped]
                    _apply_structure(
                        structure,
                        f'{ends[0]}{first}-{ends[1]}{second} is now {named}.')
            elif verb == 'bondorder' and len(fields) == 3:
                first, second, order = (int(v) for v in fields)
                named = {1: 'single', 2: 'double', 3: 'triple'}.get(order, '')
                ends = [structure.symbols[i] for i in (first, second)]
                if not structure.order(first, second):
                    # A bond that is not there yet is made the way the Bond
                    # button makes one: the topology changes and nothing else.
                    # Moving fragments and re-placing hydrogens to go with it
                    # was more than was asked for, and it wrecked the molecule
                    # on the way -- while Bond, which does none of that, has
                    # always worked.
                    state['picked'] = [first, second]
                    _edit_bond(True)
                elif set_bond_order(structure, first, second, order):
                    _apply_structure(
                        structure,
                        f'{ends[0]}{first}-{ends[1]}{second} is now {named}.')
                else:
                    _set_mol_status(
                        f'{ends[0]}{first}-{ends[1]}{second} cannot be '
                        f'{named}: one of them has no valence left for it.'
                    )
            elif verb == 'delatoms':
                doomed = [int(p) for p in fields if p.strip().lstrip('-').isdigit()]
                gone = delete_atoms(structure, doomed)
                if gone:
                    _apply_structure(structure, f'Deleted {gone} atom(s).')
        except Exception as exc:
            _set_mol_status(f'That edit did not work: {exc}')

    def on_submit_bond(_button=None):
        _edit_bond(True)

    def on_submit_unbond(_button=None):
        _edit_bond(False)

    def on_submit_swap(_button=None):
        """Exchange two ligands outright rather than dragging one at another.

        The two arrangements are separate minima and the relaxation only runs
        downhill, so it can never cross between them: a ligand dragged part of
        the way simply rolls back. The exchange is therefore performed in one
        step and the field is left to tidy up afterwards.
        """
        indices = list(state.get('picked') or [])
        metal = state.get('poly_metal')
        if metal is None:
            perceived = state.get('perceived')
            metals = list(getattr(perceived, 'metal_indices', ()) or ())
            metal = metals[0] if len(metals) == 1 else None
        if len(indices) != 2 or metal is None:
            _set_mol_status('Select two ligands of one metal to exchange them.')
            return
        state['poly_assignment'] = None
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.exchangeLigands('
            f'{json.dumps(submit_scope_id)},{int(metal)},'
            f'{int(indices[0])},{int(indices[1])});'
        )
        _set_mol_status(
            'Exchanged the two ligands. The field is settling the result; '
            'Undo puts them back.'
        )

    def on_submit_hold(_button=None):
        """Hold the value the selection describes while the field runs."""
        indices = list(state.get('picked') or [])
        kind = _CONSTRAINT_KINDS.get(len(indices))
        if not kind:
            _set_mol_status('Pick 2, 3 or 4 atoms before holding a value.')
            return
        entry = {
            'kind': kind,
            'atoms': indices,
            'value': float(submit_internal_value.value),
            # 'pull' negotiates with the chemistry and settles at a compromise;
            # 'fix' is restored after every relaxation step, so the value is met
            # exactly and the rest of the molecule arranges itself around it.
            'mode': submit_hold_mode.value,
        }
        held = list(state.get('constraints') or [])
        held = [c for c in held if c['atoms'] != indices]
        held.append(entry)
        state['constraints'] = held
        _refresh_constraints()
        _set_mol_status(f'Holding {_describe_constraint(entry)}.')
        _enable_live_forcefield()
        # Under GFN nothing is running between drags, so the change is what
        # starts it: a value set while the switch is down used to sit there
        # until Optimise was pressed.
        _arm_gfn_takeup(f'Holding {_describe_constraint(entry)}')
        # A fresh set for the next one: several values can then be held at
        # once, which is the whole point of a list.
        _clear_selection()

    def on_submit_hold_mode(change):
        """Retune the selected constraint, so a mode can be changed without
        having to select the atoms and set the value again."""
        if change.get('name') != 'value':
            return
        if state.get('hold_mode_quiet'):
            return
        position, entry = _selected_constraint()
        if entry is None:
            return
        held = list(state.get('constraints') or [])
        held[position] = dict(entry, mode=submit_hold_mode.value)
        state['constraints'] = held
        _refresh_constraints()
        _enable_live_forcefield()
        _arm_gfn_takeup(f'Holding {_describe_constraint(held[position])}')

    def on_submit_reset(_button=None):
        """Back to the structure that was loaded, with nothing set on it.

        Editing in the viewer is a one-way street otherwise: undo takes back
        one step at a time, and a structure that has been pulled apart over
        twenty of them has no way home short of pasting the coordinates again.
        """
        pristine = state.get('pristine_coords')
        if not pristine:
            _set_mol_status('Nothing to go back to yet.')
            return
        state['constraints'] = []
        state['bond_edits'] = {}
        state['hyb_overrides'] = {}
        state['structure_undo'] = []
        state['poly_applied'] = None
        state['poly_metal'] = None
        state['poly_assignment'] = None
        state['poly_arrangements'] = []
        state['poly_arrangement_index'] = 0
        submit_poly_turn_btn.layout.display = 'none'
        submit_poly_turn_btn.disabled = True
        state['poly_quiet'] = True
        try:
            submit_poly_dd.value = ''
        except Exception:
            pass
        finally:
            state['poly_quiet'] = False
        state['hold_mode_quiet'] = True
        try:
            submit_hold_mode.value = 'pull'
        except Exception:
            pass
        finally:
            state['hold_mode_quiet'] = False
        submit_internal_value.value = 0.0
        _clear_selection()
        # Writing the coordinates is what re-renders and re-perceives; it also
        # clears everything above a second time, which is the point.
        if coords_widget.value == pristine:
            # Same text, so the write would be a no-op and nothing would be
            # redrawn -- the whole reason the viewer looks destroyed is that
            # the *coordinates* changed underneath it.
            update_molecule_view()
        else:
            coords_widget.value = pristine
        _refresh_constraints()
        _set_mol_status('Back to the structure as it was loaded.')

    def on_submit_constraint_del(_button=None):
        key = submit_constraint_dd.value or ''
        if key == 'poly':
            state['poly_applied'] = None
            state['poly_metal'] = None
            state['poly_assignment'] = None
            state['poly_arrangements'] = []
            state['poly_arrangement_index'] = 0
            submit_poly_turn_btn.layout.display = 'none'
            submit_poly_turn_btn.disabled = True
            state['poly_quiet'] = True
            try:
                submit_poly_dd.value = ''
            except Exception:
                pass
            finally:
                state['poly_quiet'] = False
        elif key.startswith('c'):
            held = list(state.get('constraints') or [])
            position = int(key[1:])
            if 0 <= position < len(held):
                held.pop(position)
            state['constraints'] = held
        _refresh_constraints()
        _enable_live_forcefield()

    def on_submit_poly_changed(change):
        if change.get('name') != 'value' or state.get('poly_quiet'):
            return
        state['poly_applied'] = submit_poly_dd.value or None
        state['poly_assignment'] = None
        if state['poly_applied']:
            state['poly_metal'] = state.get('poly_offer_metal')
        else:
            state['poly_metal'] = None
        if state['poly_applied'] and state.get('poly_metal') is not None:
            try:
                from .molecule_forcefield import polyhedron_assignment
                state['poly_assignment'] = polyhedron_assignment(
                    state['perceived'], state['poly_metal'], state['poly_applied'],
                )
            except Exception:
                state['poly_assignment'] = None
        state['poly_arrangements'] = []
        state['poly_arrangement_index'] = 0
        _refresh_poly_turn()
        _refresh_constraints()
        # Re-assigning the parameters is what makes the pull start; with the
        # field running the complex visibly moves into the polyhedron.
        if submit_relax_btn.value or state.get('ff_bootstrap_done'):
            _enable_live_forcefield()
        # The metal has done its job once the polyhedron is on, exactly as
        # after Set or Hold: leaving it picked meant the next atom clicked
        # joined it, and the highlight sphere read as though something were
        # still waiting to be chosen. What is held stays in the list below.
        if state['poly_applied']:
            _set_mol_status(
                f'{submit_poly_dd.label}: the donors are pulled onto it.'
            )
        else:
            _set_mol_status('Polyhedron released.')
        _clear_selection()

    def on_submit_dyn_bonds(change):
        """Whether the drawn lines follow the distances."""
        if change.get('name') != 'value':
            return
        active = bool(submit_dyn_bonds_btn.value)
        submit_dyn_bonds_btn.button_style = 'info' if active else ''
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setDynamicBonds('
            f'{json.dumps(submit_scope_id)},{"true" if active else "false"});'
        )
        _set_mol_status(
            'The lines follow the distances now. Only the picture: what the '
            'calculation holds together is Bond and Unbond\'s business.'
            if active else
            'The lines keep the bonds the structure was drawn with.')

    def on_submit_settle_toggle(change):
        if change.get('name') != 'value':
            return
        active = bool(submit_settle_btn.value)
        submit_settle_btn.button_style = 'info' if active else ''
        if _gfn.is_gfn_method(submit_ff_dd.value):
            # Under GFN this is the server's job, not the browser's: a release
            # runs one short optimisation with the method on screen.  Telling
            # the browser to settle would be telling it to use a field that is
            # deliberately not installed.
            label = _gfn.GFN_METHODS[str(submit_ff_dd.value)]['label']
            _set_mol_status(
                f'Letting go of an atom now lets {label} tidy the structure '
                'around it.' if active else
                'Atoms will stay exactly where you put them.')
            return
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setSettleOnRelease('
            f'{json.dumps(submit_scope_id)},{"true" if active else "false"});'
        )

    def on_submit_strength_changed(change):
        if change.get('name') != 'value':
            return
        _ensure_manip_bootstrap()
        _run_manip_js(
            'if(window.__delfinSubmitManip)'
            'window.__delfinSubmitManip.setOptimizerStrength('
            f'{json.dumps(submit_scope_id)},{int(submit_strength_slider.value)});'
        )

    def on_submit_solvent(change):
        if change.get('name') != 'value' or state.get('solvent_quiet'):
            return
        wet = str(submit_gfn_solvent.value or '')
        _set_mol_status(
            f'Optimising in {_gfn.SOLVENTS[wet]} from now on (ALPB). The '
            'structure on screen was not, so it is worth optimising again.'
            if wet else
            'Optimising in the gas phase from now on.')

    def on_submit_autospin(change):
        if change.get('name') != 'value':
            return
        scanning = bool(submit_gfn_autospin.value)
        submit_gfn_mult.disabled = scanning
        if scanning:
            _set_mol_status(
                'The multiplicity is scanned: the ones the electron count '
                'allows are each optimised and the lowest is kept. Three runs '
                'instead of one.')

    def _fill_charge_from_smiles():
        """Take the charge off the SMILES the structure was built from.

        A SMILES states the formal charges outright, so asking the user for a
        number the input already carries is asking them to repeat themselves.
        Nothing else can be read that way -- a pasted XYZ says nothing about
        charge, and no input says anything about the spin -- so those stay the
        user's to set.  Returns the SMILES it read, or '' when there was none.
        """
        smiles = str((state.get('converted_xyz_cache') or {}).get('smiles') or '')
        if not smiles:
            return ''
        try:
            submit_gfn_charge.value = int(_get_smiles_charge(smiles))
        except Exception:
            return ''
        return smiles

    def _live_ff_method():
        """The method the browser relaxes with while an atom is dragged.

        GFN runs on the server; a round trip per frame would cap a drag at
        about 13 Hz, so dragging keeps the force field that lives in the
        browser.  Choosing GFN changes what Optimise does, not what a drag
        does, and the status line says so once rather than silently.
        """
        chosen = str(submit_ff_dd.value or 'uff')
        return 'uff' if _gfn.is_gfn_method(chosen) else chosen

    def _offer_xtb_install(offer):
        """Show the offer, or the two buttons that carry it out, or neither."""
        submit_xtb_install_btn.layout.display = '' if offer else 'none'
        submit_xtb_confirm_btn.layout.display = 'none'
        submit_xtb_cancel_btn.layout.display = 'none'

    def _missing_tool():
        """Which program the chosen method needs and has not got.

        g-xTB is not the xtb beside it: it ships as a build of its own, and an
        ordinary xtb accepts --gxtb and silently runs GFN2.  So what is missing
        is a question with two answers, and offering the wrong one would
        install something that changes nothing.
        """
        method = str(submit_ff_dd.value)
        if not _gfn.is_gfn_method(method):
            return None
        return None if _gfn.find_binary(method) else (
            'gxtb' if method == 'gxtb' else 'xtb')

    def _refresh_xtb_offer():
        """A missing program under a GFN method is a thing to fix, not a wall."""
        missing = _missing_tool()
        state['xtb_install_tool'] = missing
        submit_xtb_install_btn.description = (
            'Install g-xTB' if missing == 'gxtb' else 'Install xtb')
        _offer_xtb_install(
            missing is not None
            and _gfn.install_script() is not None
            and not state.get('xtb_installing'))

    def on_submit_xtb_install(_button=None):
        """The first press only says what the second one would do.

        Fetching a conda environment of a few hundred megabytes is not a thing
        to start on a single click, and on a cluster the right answer is often
        to load the module instead -- so the command is on screen before it
        runs, and cancelling is as easy as agreeing.
        """
        command = _gfn.install_command(state.get('xtb_install_tool') or 'xtb')
        if command is None:
            _set_mol_status('DELFIN\'s installer is not next to this copy of '
                            'the dashboard, so there is nothing to run.')
            return
        submit_xtb_install_btn.layout.display = 'none'
        submit_xtb_confirm_btn.layout.display = ''
        submit_xtb_cancel_btn.layout.display = ''
        _set_mol_status(
            'This will run:  ' + ' '.join(command),
            'It fetches xtb from conda-forge into DELFIN\'s own tool '
            'directory -- a few hundred megabytes, several minutes, and it '
            'needs the network. On a cluster, module load xtb may be the '
            'better answer.',
        )

    def on_submit_xtb_cancel(_button=None):
        _refresh_xtb_offer()
        _set_mol_status('Left alone. Optimise will say xtb is missing rather '
                        'than doing nothing.')

    def on_submit_xtb_confirm(_button=None):
        if state.get('xtb_installing'):
            return
        state['xtb_installing'] = True
        _offer_xtb_install(False)
        _set_mol_status(
            f'Installing {state.get("xtb_install_tool") or "xtb"}...',
            spinner=True)

        def _work():
            def _line(text):
                if text.strip():
                    _schedule_ui_update(
                        _set_mol_status,
                        f'Installing {state.get("xtb_install_tool") or "xtb"}'
                        '...', text.strip(), spinner=True)

            outcome = _gfn.install_xtb(
                on_line=_line, tool=state.get('xtb_install_tool') or 'xtb')

            def _done():
                state['xtb_installing'] = False
                _refresh_xtb_offer()
                if outcome.get('ok'):
                    # Asked again rather than assumed: the resolver is what
                    # every run will use, and a binary it cannot see is not
                    # installed as far as the dashboard is concerned.
                    _set_mol_status(outcome['status'],
                                    'Optimise will use it from now on.')
                else:
                    _set_mol_status(outcome.get('status') or 'The install '
                                    'did not finish.',
                                    *outcome.get('lines', [])[-2:])

            _schedule_ui_update(_done)

        threading.Thread(target=_work, daemon=True).start()

    def on_submit_ff_changed(change):
        if change.get('name') != 'value':
            return
        gfn = _gfn.is_gfn_method(submit_ff_dd.value)
        submit_gfn_charge.layout.display = '' if gfn else 'none'
        submit_gfn_mult.layout.display = '' if gfn else 'none'
        submit_gfn_autospin.layout.display = '' if gfn else 'none'
        # A method without solvation gets no solvent box: a control that can
        # only produce a refusal is worse than no control.
        solvated = gfn and _gfn.GFN_METHODS[
            str(submit_ff_dd.value)].get('solvation') is not False
        submit_gfn_solvent.layout.display = '' if solvated else 'none'
        if not solvated and submit_gfn_solvent.value:
            state['solvent_quiet'] = True
            try:
                submit_gfn_solvent.value = ''
            finally:
                state['solvent_quiet'] = False
        # Relax means the browser's own field running once per frame, and there
        # is no GFN engine in the browser to run.  Under GFN it means the other
        # half of the same idea: while an atom is being dragged, the rest of
        # the molecule follows it -- one short xtb run per push, and nothing at
        # all when nothing is being dragged.  It was switched off and hidden
        # here before there was anything for it to do.
        # A follow armed under the old method must not outlive the choice that
        # armed it: the toggle handler reads the method that is chosen *now*.
        _end_gfn_follow()
        if gfn and submit_settle_btn.value:
            # Under GFN, letting go leaves the structure where the hand left
            # it.  Tidying on every release is a thing to ask for, not the
            # default -- an atom pulled off the place it was just put is the
            # opposite of placing it.
            submit_settle_btn.value = False
        submit_relax_btn.disabled = False
        # The page reads this switch itself, and only follows when the class is
        # on it: asking the kernel whether to follow would cost a round trip
        # per push, and a UFF drag would be pushing for no reason.
        if gfn:
            submit_relax_btn.add_class('submit-gfn-follow')
        else:
            submit_relax_btn.remove_class('submit-gfn-follow')
        # Out of the way entirely, not merely greyed: a control that cannot do
        # anything under the chosen method is clutter that invites the question
        # of why it is dead.  Strength is how many steps the browser's field
        # takes per animation frame, and that field does not run here.  Settle
        # does have a meaning under GFN -- letting go and having the structure
        # tidied rather than left at the cursor -- so it stays, and the chosen
        # method is what tidies it.
        submit_strength_slider.layout.display = 'none' if gfn else ''
        # Dynamik Opt without an xtb behind it cannot do anything at all, so
        # it goes rather than sitting there being pressed to no effect.
        usable = not gfn or _gfn.find_binary(submit_ff_dd.value) is not None
        submit_relax_btn.layout.display = '' if usable else 'none'
        if not usable and submit_relax_btn.value:
            submit_relax_btn.value = False
        label = (_gfn.GFN_METHODS[str(submit_ff_dd.value)]['label']
                 if gfn else str(submit_ff_dd.value).upper())
        submit_relax_btn.tooltip = (
            f'While this is on, dragging an atom lets the rest of the molecule '
            f'follow it -- {label} on the server, a few cycles per push, and '
            f'each step says how long it took. Nothing runs while nothing is '
            f'being dragged.'
            if gfn else
            'Relax the structure continuously while you work on it.'
        )
        submit_settle_btn.tooltip = (
            f'When you let go of an atom, let {label} tidy the structure '
            f'around its new position instead of leaving it where the cursor '
            f'stopped. One short run on the server per release.'
            if gfn else
            'When you let go of an atom, let the structure relax around its '
            'new position instead of keeping the strain of the drag. Switch '
            'off to leave atoms exactly where you put them.'
        )
        # The engine follows the method, and the switch keeps its position:
        # switching back and forth is something a user does constantly, and it
        # must not cost a press each time -- nor leave the previous engine
        # running underneath the new choice.
        if gfn:
            _stop_browser_field()
        elif submit_relax_btn.value:
            _enable_live_forcefield()
            _run_manip_js(
                'if(window.__delfinSubmitManip)'
                'window.__delfinSubmitManip.startAutoOptimize('
                f'{json.dumps(submit_scope_id)});'
            )
        _refresh_xtb_offer()
        if gfn:
            label = _gfn.GFN_METHODS[str(submit_ff_dd.value)]['label']
            source = _fill_charge_from_smiles()
            if _gfn.find_binary(submit_ff_dd.value) is None:
                _set_mol_status(
                    f'{label} needs a program that was not found. The button '
                    'fetches it with DELFIN\'s own installer, and says what '
                    'it will run before it runs it.'
                    if _gfn.install_script() is not None else
                    f'{label} needs xtb on the PATH; it was not found. '
                    'Optimise will say so rather than doing nothing.')
            elif source:
                _set_mol_status(
                    f'Optimise now uses {label}. Charge {submit_gfn_charge.value} '
                    f'read from the SMILES; the multiplicity is yours to set. '
                    'Switch Relax on to have the molecule follow an atom you '
                    'drag.')
            else:
                _set_mol_status(
                    f'Optimise now uses {label}. Set the charge (q) and the '
                    'multiplicity (M): xtb needs both, and a wrong spin on a '
                    'metal gives a confident wrong answer rather than an error. '
                    'Switch Relax on to have the molecule follow an atom you '
                    'drag.')
        # Re-assign parameters under the newly chosen method, but only if the
        # live relaxation is actually switched on.
        if submit_relax_btn.value:
            _enable_live_forcefield()

    def on_submit_manip_clear(_button=None):
        _ensure_manip_bootstrap()
        _run_manip_js(
            f'if(window.__delfinSubmitManip) '
            f'window.__delfinSubmitManip.clear({json.dumps(submit_scope_id)});'
        )

    def on_submit_manip_undo(_button=None):
        snapshot = state.pop('pre_optimize_frames', None)
        if snapshot:
            if snapshot['isomers']:
                state['isomers'] = snapshot['isomers']
                _show_isomer_at_index(state.get('isomer_index', 0))
            else:
                coords_widget.value = snapshot['coords']
            _set_mol_status('Reverted to the geometries from before the optimisation.')
            return
        _ensure_manip_bootstrap()
        _run_manip_js(
            f'if(window.__delfinSubmitManip) '
            f'window.__delfinSubmitManip.undo({json.dumps(submit_scope_id)});'
        )

    def on_submit_manip_sync(change):
        if change.get('name') != 'value':
            return
        new_xyz = submit_manip_sync.value
        if not new_xyz or not new_xyz.strip():
            return
        # Extract only the new coordinate lines; drop JS-side count + comment.
        new_lines = new_xyz.splitlines()
        if len(new_lines) >= 2:
            try:
                int(new_lines[0].strip())
                coord_lines = new_lines[2:]
            except ValueError:
                coord_lines = new_lines
        else:
            coord_lines = new_lines
        coord_body = '\n'.join(line for line in coord_lines if line.strip())
        # Preserve the user's original header (atom count + comment line) if
        # present in the current coords_widget value.
        old_lines = coords_widget.value.splitlines()
        header = ''
        if len(old_lines) >= 2:
            try:
                int(old_lines[0].strip())
                header = f'{old_lines[0]}\n{old_lines[1]}\n'
            except ValueError:
                pass
        # A drag has just finished. If a polyhedron is being held, work out
        # again which donor is now nearest which vertex: dragging a ligand
        # towards another position and having the field haul it straight back
        # is not an exchange, it is a fight. Recomputing here means the
        # polyhedron accepts the ligand where it has been put and pulls it the
        # rest of the way onto the vertex it is now closest to.
        lines = new_xyz.splitlines()
        note = lines[1].strip() if len(lines) > 1 else ''
        drag_ended = note.startswith('DELFIN drag-end')
        # Sent while the mouse is still down, so the molecule can follow the
        # atom rather than wait for it to be let go.  The comment line names
        # the atoms the hand is on, so the answer can keep them there.
        if note.startswith('DELFIN drag-follow'):
            holding = []
            for word in note.split():
                if word.startswith('held='):
                    holding = [int(n) for n in word[5:].split(',')
                               if n.strip().lstrip('-').isdigit()]
            _gfn_follow_step(new_xyz, holding)
        if (drag_ended and state.get('poly_applied')
                and state.get('poly_metal') is not None):
            # Only a real end of a drag, not the twice-a-second heartbeat the
            # running optimiser sends: reassigning on every heartbeat reloaded
            # the whole field twice a second and never let a moved ligand
            # settle onto its new vertex.
            state['poly_assignment'] = None
            state['poly_recheck'] = True

        if drag_ended:
            # Set, Hold, a bond edit and a drag all arrive here.  Any of them
            # during an optimisation makes what xtb is doing about a structure
            # that is no longer on the screen, so it starts again from the one
            # that is -- and the geometry lands first, so it is the new one it
            # starts from.
            _interrupt_gfn()
            _arm_gfn_restart()

        payload = header + coord_body
        # The guard is cleared by update_molecule_view, which traitlets only
        # calls when the value actually changes. Dragging an atom out and back,
        # or any edit that lands on the same coordinates, would otherwise leave
        # the flag set for the rest of the session and swallow the user's next
        # genuine edit of the coordinate box.
        if coords_widget.value == payload:
            state['manip_inflight'] = False
            return
        state['manip_inflight'] = True
        coords_widget.value = payload
        if state.pop('poly_recheck', False):
            # After the coordinates have landed, so the assignment is worked
            # out from where the ligands actually are now.
            _schedule_ui_update(_enable_live_forcefield)

    # -- drawing the structure ------------------------------------------
    def _draw_frame_html(url):
        """The editor itself, in a frame of its own.

        A frame rather than the page: Ketcher is a React application that owns
        its own globals, and the dashboard is another one.  Same origin, so the
        page may reach in and ask it for the drawing -- across origins there
        would be nothing to ask with, because Ketcher speaks no messages.
        """
        return (
            "<iframe src='" + html.escape(url, quote=True) + "' "
            "style='width:100%; height:560px; border:1px solid #d0d0d0; "
            "border-radius:6px; background:#fff;' "
            "title='Ketcher'></iframe>"
        )

    def _refresh_draw_controls():
        drawn = bool(submit_draw_open_btn.value)
        ready = _ketcher.is_installed()
        submit_draw_frame.layout.display = '' if (drawn and ready) else 'none'
        submit_draw_get_btn.layout.display = '' if (drawn and ready) else 'none'
        submit_draw_update_btn.layout.display = '' if (drawn and ready) else 'none'

    def on_submit_draw_open(change):
        if change.get('name') != 'value':
            return
        if not submit_draw_open_btn.value:
            submit_draw_frame.value = ''      # stop the app when it is closed
            _refresh_draw_controls()
            return
        url = _ketcher.app_url()
        if url:
            version = _ketcher.installed_version()
            submit_draw_frame.value = _draw_frame_html(url)
            _refresh_draw_controls()
            _set_mol_status(
                f'Ketcher {version}: draw the structure, then press TO SMILES '
                'to put it in the input box.')
            return
        # Not there yet.  Offered rather than fetched: it is thirty-odd
        # megabytes, and on a machine without a network it is a wait that ends
        # in nothing.
        _refresh_draw_controls()
        if _ketcher.app_directory() is None:
            submit_draw_open_btn.value = False
            _set_mol_status(
                'The drawing editor needs a directory the browser can load it '
                'from, and this dashboard is not serving one.')
            return
        state['draw_installing'] = True
        _set_mol_status('Looking up the newest Ketcher...', spinner=True)

        def _work():
            newest = _ketcher.latest_release()

            def _ask():
                state['draw_installing'] = False
                if not newest['ok']:
                    submit_draw_open_btn.value = False
                    _set_mol_status(newest['status'])
                    return
                state['draw_offer'] = newest
                submit_draw_get_btn.layout.display = 'none'
                submit_draw_update_btn.description = 'Fetch it'
                submit_draw_update_btn.button_style = 'warning'
                submit_draw_update_btn.layout.display = ''
                _set_mol_status(
                    f'Ketcher {newest["version"]} is not here yet: '
                    f'{newest["size"] / 1e6:.0f} MB from '
                    'github.com/epam/ketcher, unpacked to about 32 MB beside '
                    'the dashboard. Press Fetch it to get it.')

            _schedule_ui_update(_ask)

        threading.Thread(target=_work, daemon=True).start()

    def on_submit_draw_update(_button=None):
        """Fetch the newest build -- the first time, or over an older one."""
        if state.get('draw_installing'):
            return
        state['draw_installing'] = True
        submit_draw_update_btn.layout.display = 'none'
        _set_mol_status('Fetching Ketcher...', spinner=True)

        def _work():
            def _line(text):
                _schedule_ui_update(_set_mol_status, f'Ketcher: {text}',
                                    spinner=True)

            outcome = _ketcher.install(on_line=_line)

            def _done():
                state['draw_installing'] = False
                state['draw_offer'] = None
                submit_draw_update_btn.description = 'Update'
                submit_draw_update_btn.button_style = ''
                if not outcome['ok']:
                    submit_draw_open_btn.value = False
                    _refresh_draw_controls()
                    _set_mol_status(outcome['status'])
                    return
                url = _ketcher.app_url()
                submit_draw_frame.value = _draw_frame_html(url) if url else ''
                if not submit_draw_open_btn.value:
                    submit_draw_open_btn.value = True
                _refresh_draw_controls()
                _set_mol_status(outcome['status'],
                                'Draw the structure, then press TO SMILES.')

            _schedule_ui_update(_done)

        threading.Thread(target=_work, daemon=True).start()

    def on_submit_draw_get(_button=None):
        """Ask the editor for what has been drawn.

        The molfile, not the SMILES the editor could write itself: everything
        downstream here reads structures with RDKit, and a SMILES RDKit wrote
        is one RDKit will certainly read back.
        """
        _set_mol_status('Reading the drawing...', spinner=True)
        ctx.run_js(
            "(function(){\n"
            "  var box=document.querySelector('.submit-ketcher-sync');\n"
            "  var input=box&&box.querySelector('textarea, input');\n"
            "  function hand(text){\n"
            "    if(!input) return;\n"
            "    var proto=(input.tagName==='TEXTAREA')\n"
            "      ? window.HTMLTextAreaElement.prototype\n"
            "      : window.HTMLInputElement.prototype;\n"
            "    var setter=Object.getOwnPropertyDescriptor(proto,'value');\n"
            "    /* A serial in front, so drawing the same thing twice reads\n"
            "       as two answers rather than as one that never came. */\n"
            "    var line=(Date.now())+'\\n'+text;\n"
            "    if(setter&&setter.set) setter.set.call(input,line);\n"
            "    else input.value=line;\n"
            "    input.dispatchEvent(new Event('input',{bubbles:true}));\n"
            "    input.dispatchEvent(new Event('change',{bubbles:true}));\n"
            "  }\n"
            "  var host=document.querySelector('.submit-ketcher-frame');\n"
            "  var frame=host&&host.querySelector('iframe');\n"
            "  var api=null;\n"
            "  try{ api=frame&&frame.contentWindow&&frame.contentWindow.ketcher; }\n"
            "  catch(e){ api=null; }\n"
            "  if(!api){ hand('!no-editor'); return; }\n"
            "  try{\n"
            "    Promise.resolve(api.getMolfile()).then(function(mol){\n"
            "      hand(mol||''); }, function(err){ hand('!'+err); });\n"
            "  }catch(e){ hand('!'+e); }\n"
            "})();"
        )

    def on_submit_draw_sync(change):
        """What the editor handed back."""
        if change.get('name') != 'value':
            return
        raw = submit_draw_sync.value or ''
        if '\n' not in raw:
            return
        molfile = raw.split('\n', 1)[1]
        if molfile.startswith('!'):
            trouble = molfile[1:]
            _set_mol_status(
                'The editor is not open yet, so there is nothing to read.'
                if trouble == 'no-editor' else
                f'The drawing could not be read: {trouble}')
            return
        outcome = _ketcher.smiles_from_molfile(molfile)
        if not outcome['ok']:
            _set_mol_status(outcome['status'])
            return
        # Into the box the rest of the tab reads, so Convert, Build and every
        # other button downstream sees it exactly as a typed SMILES.
        coords_widget.value = outcome['smiles']
        _set_mol_status(f'Drawn: {outcome["smiles"]}',
                        'It is in the input box -- Convert turns it into '
                        'coordinates.')

    # -- wiring ---------------------------------------------------------
    xyz_copy_btn.on_click(on_xyz_copy)
    coords_widget.observe(update_molecule_view, names='value')
    submit_select_btn.observe(on_submit_select_toggle, names='value')
    submit_manip_btn.observe(on_submit_manip_toggle, names='value')
    submit_manip_clear_btn.on_click(on_submit_manip_clear)
    submit_manip_undo_btn.on_click(on_submit_manip_undo)
    submit_relax_btn.observe(on_submit_relax_toggle, names='value')
    # On the page from the start: a player installed at click time races the
    # run it is meant to show.
    _install_gfn_frame_watcher()
    submit_ff_dd.observe(on_submit_ff_changed, names='value')
    submit_xtb_install_btn.on_click(on_submit_xtb_install)
    submit_xtb_confirm_btn.on_click(on_submit_xtb_confirm)
    submit_xtb_cancel_btn.on_click(on_submit_xtb_cancel)
    submit_gfn_autospin.observe(on_submit_autospin, names='value')
    submit_gfn_solvent.observe(on_submit_solvent, names='value')
    submit_strength_slider.observe(on_submit_strength_changed, names='value')
    submit_settle_btn.observe(on_submit_settle_toggle, names='value')
    submit_dyn_bonds_btn.observe(on_submit_dyn_bonds, names='value')
    submit_pick_sync.observe(on_submit_pick_sync, names='value')
    submit_poly_dd.observe(on_submit_poly_changed, names='value')
    submit_hyb_dd.observe(on_submit_hyb_changed, names='value')
    submit_hyb_auto_btn.on_click(on_submit_hyb_auto)
    submit_poly_turn_btn.on_click(on_submit_poly_turn)
    submit_cmd_sync.observe(on_submit_cmd, names='value')
    submit_draw_btn.observe(on_submit_draw_toggle, names='value')
    submit_element_dd.observe(on_submit_draw_choice, names='value')
    submit_hold_btn.on_click(on_submit_hold)
    submit_swap_btn.on_click(on_submit_swap)
    submit_hold_mode.observe(on_submit_hold_mode, names='value')
    submit_bond_btn.on_click(on_submit_bond)
    submit_unbond_btn.on_click(on_submit_unbond)
    submit_constraint_del.on_click(on_submit_constraint_del)
    submit_constraint_dd.observe(on_submit_constraint_selected, names='value')
    submit_internal_value.observe(on_submit_constraint_retune, names='value')
    submit_reset_btn.on_click(on_submit_reset)
    submit_internal_btn.observe(on_submit_set_internal, names='value')
    submit_internal_value.observe(on_submit_internal_value, names='value')
    # The page watches this button itself: waiting for the kernel to say the
    # switch went off costs a round trip, and the playback ran on for it.
    submit_optimize_btn.add_class('submit-optimize-switch')
    submit_optimize_btn.observe(on_submit_optimize, names='value')
    submit_optimize_all_btn.add_class('submit-optimize-switch')
    submit_optimize_all_btn.observe(on_submit_optimize_all, names='value')
    submit_manip_sync.observe(on_submit_manip_sync, names='value')
    submit_draw_open_btn.observe(on_submit_draw_open, names='value')
    submit_draw_get_btn.on_click(on_submit_draw_get)
    submit_draw_update_btn.on_click(on_submit_draw_update)
    submit_draw_sync.observe(on_submit_draw_sync, names='value')
    convert_smiles_button.on_click(handle_convert_smiles)
    convert_smiles_quick_button.on_click(handle_convert_smiles_quick)
    convert_smiles_uff_button.on_click(handle_convert_smiles_uff)
    build_complex_button.on_click(handle_build_complex)
    architector_button.on_click(handle_architector_convert)
    manta_button.on_click(handle_manta)
    guppy_submit_button.on_click(handle_guppy_submit)
    fukui_submit_button.on_click(handle_fukui_submit)
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
    xyz_copy_row.add_class('submit-fs-member-copyrow')
    # What the force field had to approximate belongs under the structure it
    # describes, not in the preview's status line where it competes with
    # conversion messages and scrolls away.
    submit_ff_notes = widgets.HTML(
        value='',
        layout=widgets.Layout(width='100%', margin='4px 0 0 0'),
    )
    submit_ff_notes.add_class('submit-ff-notes')
    submit_manip_toolbar.add_class('submit-fs-member-toolbar')
    mol_output.add_class('submit-fs-member-viewer')
    isomer_nav_row.add_class('submit-fs-member-isomer')

    submit_right = widgets.VBox([
        widgets.HTML('<b>Molecule Preview:</b>'), mol_status,
        submit_manip_toolbar, mol_status_fs, mol_output, isomer_nav_row, xyz_copy_row,
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
        .submit-fs-overlay {
            position: fixed !important;
            top: 0 !important;
            left: 0 !important;
            right: 0 !important;
            bottom: 0 !important;
            z-index: 9999 !important;
            background: #ffffff !important;
            padding: 8px !important;
            margin: 0 !important;
            overflow: hidden !important;
            box-sizing: border-box !important;
            display: flex !important;
            flex-direction: column !important;
            gap: 6px;
        }
        .submit-fs-overlay .submit-fs-member-status {
            display: block !important;
        }
        .submit-fs-overlay .submit-fs-member-viewer {
            flex: 1 1 auto !important;
            height: auto !important;
            min-height: 0 !important;
            width: 100% !important;
        }
        .submit-fs-overlay .submit-fs-member-viewer .output_area,
        .submit-fs-overlay .submit-fs-member-viewer .output_subarea,
        .submit-fs-overlay .submit-fs-member-viewer .jp-OutputArea,
        .submit-fs-overlay .submit-fs-member-viewer .jp-OutputArea-output,
        .submit-fs-overlay .submit-fs-member-viewer .jp-OutputArea-child {
            height: 100% !important;
            width: 100% !important;
        }
        .submit-fs-overlay .submit-fs-member-viewer [id^="3dmolviewer_"] {
            width: 100% !important;
            height: 100% !important;
            max-width: none !important;
            max-height: none !important;
        }
        .submit-fs-overlay .submit-fs-member-viewer [id^="3dmolviewer_"] canvas {
            width: 100% !important;
            height: 100% !important;
        }
        </style>
        """
    )
    tab_widget.add_class('submit-split-root')
    submit_left.add_class('submit-split-pane')
    submit_right.add_class('submit-split-pane')
    mol_output.add_class('submit-split-pane')
    mol_output.add_class('submit-mol-output')
    _replace_mol_output_text('Please enter XYZ coordinates or SMILES.')
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
        # editor state, for the held-value list and Reset
        'submit_constraint_dd': submit_constraint_dd,
        'submit_hold_btn': submit_hold_btn,
        'submit_hold_mode': submit_hold_mode,
        'submit_internal_value': submit_internal_value,
        'submit_internal_btn': submit_internal_btn,
        'submit_ff_dd': submit_ff_dd,
        'submit_gfn_charge': submit_gfn_charge,
        'submit_gfn_mult': submit_gfn_mult,
        'submit_gfn_autospin': submit_gfn_autospin,
        'submit_gfn_solvent': submit_gfn_solvent,
        'submit_dyn_bonds_btn': submit_dyn_bonds_btn,
        'submit_draw_open_btn': submit_draw_open_btn,
        'submit_draw_get_btn': submit_draw_get_btn,
        'submit_draw_update_btn': submit_draw_update_btn,
        'submit_draw_frame': submit_draw_frame,
        'submit_draw_sync': submit_draw_sync,
        'submit_xtb_install_btn': submit_xtb_install_btn,
        'submit_xtb_confirm_btn': submit_xtb_confirm_btn,
        'submit_xtb_cancel_btn': submit_xtb_cancel_btn,
        'submit_optimize_btn': submit_optimize_btn,
        'submit_optimize_all_btn': submit_optimize_all_btn,
        'submit_settle_btn': submit_settle_btn,
        'submit_strength_slider': submit_strength_slider,
        'submit_relax_btn': submit_relax_btn,
        'submit_gfn_frame': submit_gfn_frame,
        # The two channels the page speaks through, and the line it is
        # answered on: a test that cannot use them has to drive the editor
        # through its internals instead of the way the browser drives it.
        'submit_cmd_sync': submit_cmd_sync,
        'submit_manip_sync': submit_manip_sync,
        'mol_status': mol_status,
        'submit_pick_sync': submit_pick_sync,
        'submit_reset_btn': submit_reset_btn,
        'editor_state': state,
        'refresh_constraints': _refresh_constraints,
    }
