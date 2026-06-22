"""ChemDarwin tab – Custom Reaction SMARTS morphing for the DELFIN Dashboard.

This module is **optional and local-only**: it is listed in ``.gitignore`` so
that only developers who have the file locally will see the tab.  All
chemistry logic is self-contained (no imports from ChemDarwin2/).
"""

import io
import json
import re
import base64
from pathlib import Path

import ipywidgets as widgets
import numpy as np
from IPython.display import HTML, clear_output, display, Markdown
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import (
    AllChem,
    Descriptors,
    rdChemReactions,
    rdDepictor,
    rdFMCS,
)
from rdkit.Chem.Draw import MolsToGridImage, MolToImage

from .molecule_viewer import apply_molecule_view_style

RDLogger.DisableLog('rdApp.*')

try:
    from PIL import Image
except ImportError:
    Image = None

try:
    import py3Dmol
except ImportError:
    py3Dmol = None

try:
    import stk  # noqa: F401
    STK_AVAILABLE = True
except ImportError:
    STK_AVAILABLE = False

# ---------------------------------------------------------------------------
# Metal list
# ---------------------------------------------------------------------------
METALS = [
    'Li', 'Na', 'K', 'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
    'Ir', 'Pt', 'Au', 'Hg', 'Al', 'Ga', 'In', 'Tl', 'Sn', 'Pb', 'Bi',
]


# ===================================================================
# Helper utilities
# ===================================================================

def _count_aromatic_rings(mol):
    """Count aromatic rings (replaces m5.count_aromatic_rings)."""
    ri = mol.GetRingInfo()
    return sum(
        1 for ring in ri.BondRings()
        if all(mol.GetBondWithIdx(b).GetIsAromatic() for b in ring)
    )


def contains_metal(smiles):
    for metal in METALS:
        if re.search(rf'\[{metal}[+\-\d\]@H]', smiles, re.IGNORECASE):
            return True
        if re.search(rf'\[{metal}\]', smiles, re.IGNORECASE):
            return True
    return False


# ===================================================================
# Metal-complex support (dative bonds)
# ===================================================================

def aromatize_ligand_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return smiles
    try:
        Chem.SanitizeMol(mol, sanitizeOps=(
            Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        ))
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        return Chem.MolToSmiles(mol)
    except Exception:
        pass
    try:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C' and atom.IsInRing():
                has_double = any(
                    b.GetBondType() == Chem.BondType.DOUBLE for b in atom.GetBonds()
                )
                if has_double:
                    atom.SetIsAromatic(True)
        for bond in mol.GetBonds():
            if bond.GetBeginAtom().GetIsAromatic() and bond.GetEndAtom().GetIsAromatic():
                bond.SetIsAromatic(True)
                bond.SetBondType(Chem.BondType.AROMATIC)
        return Chem.MolToSmiles(mol)
    except Exception:
        return smiles


def extract_ligands_from_complex(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None, None, None, None, None, "Could not parse SMILES"

    metal_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in METALS]
    if not metal_atoms:
        return None, None, None, None, None, "No metal found in SMILES"

    metal_idx = metal_atoms[0]
    metal_atom = mol.GetAtomWithIdx(metal_idx)
    metal_symbol = metal_atom.GetSymbol()
    metal_charge = metal_atom.GetFormalCharge()

    binding_info = [nbr.GetIdx() for nbr in metal_atom.GetNeighbors()]

    orig_props = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        orig_props[idx] = {
            'formal_charge': atom.GetFormalCharge(),
            'explicit_h': atom.GetNumExplicitHs(),
            'no_implicit': atom.GetNoImplicit(),
            'rad_e': atom.GetNumRadicalElectrons(),
            'symbol': atom.GetSymbol(),
        }

    for i, idx in enumerate(binding_info):
        try:
            mol.GetAtomWithIdx(idx).SetAtomMapNum(901 + i)
        except Exception:
            pass

    for idx in binding_info:
        try:
            props = orig_props.get(idx)
            if props is not None:
                props['no_implicit'] = True
                props['explicit_h'] = 0
        except Exception:
            pass

    edit_mol = Chem.RWMol(mol)
    bonds_to_remove = []
    for bond in edit_mol.GetBonds():
        if bond.GetBeginAtomIdx() == metal_idx or bond.GetEndAtomIdx() == metal_idx:
            bonds_to_remove.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    for b in bonds_to_remove:
        edit_mol.RemoveBond(b[0], b[1])
    edit_mol.RemoveAtom(metal_idx)

    try:
        mol_no_metal = edit_mol.GetMol()
        for old_idx, props in orig_props.items():
            if old_idx == metal_idx:
                continue
            new_idx = old_idx - 1 if old_idx > metal_idx else old_idx
            try:
                a = mol_no_metal.GetAtomWithIdx(new_idx)
                a.SetFormalCharge(props['formal_charge'])
                a.SetNumExplicitHs(props['explicit_h'])
                a.SetNoImplicit(props['no_implicit'])
                a.SetNumRadicalElectrons(props['rad_e'])
            except Exception:
                pass
        try:
            mol_no_metal.UpdatePropertyCache(strict=False)
        except Exception:
            pass

        frags = Chem.GetMolFrags(mol_no_metal, asMols=True, sanitizeFrags=False)
        ligand_smiles = []
        ligand_binding_maps = []
        ligand_binding_symbols = []

        for frag_mol in frags:
            frag_map_nums = []
            frag_binding_syms = []
            for atom in frag_mol.GetAtoms():
                amap = atom.GetAtomMapNum()
                if amap >= 901:
                    frag_map_nums.append(amap)
                    frag_binding_syms.append(atom.GetSymbol())
            smi_with_maps = Chem.MolToSmiles(frag_mol, canonical=True)
            ligand_smiles.append(smi_with_maps)
            ligand_binding_maps.append(sorted(frag_map_nums))
            ligand_binding_symbols.append(frag_binding_syms)

        return metal_symbol, metal_charge, ligand_smiles, ligand_binding_maps, ligand_binding_symbols, None
    except Exception as e:
        return None, None, None, None, None, str(e)


def _strip_atom_maps(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return smiles
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, canonical=True)


def _binding_atom_signature(atom):
    return (atom.GetSymbol(), atom.GetIsAromatic())


def _find_binding_atoms_by_map(mol, map_nums):
    map_to_idx = {}
    for atom in mol.GetAtoms():
        amap = atom.GetAtomMapNum()
        if amap in map_nums:
            map_to_idx[amap] = atom.GetIdx()
    result = []
    for m in map_nums:
        if m in map_to_idx:
            result.append(map_to_idx[m])
        else:
            return None
    return result


def _map_binding_by_substructure(orig_mol, mut_mol, bind_indices):
    if orig_mol is None or mut_mol is None:
        return None
    try:
        match = mut_mol.GetSubstructMatch(orig_mol)
    except Exception:
        match = ()
    if match and len(match) == orig_mol.GetNumAtoms():
        return [match[b] for b in bind_indices if b < len(match)]
    try:
        o = Chem.Mol(orig_mol)
        m = Chem.Mol(mut_mol)
        Chem.Kekulize(o, clearAromaticFlags=True)
        Chem.Kekulize(m, clearAromaticFlags=True)
        match = m.GetSubstructMatch(o)
        if match and len(match) == o.GetNumAtoms():
            return [match[b] for b in bind_indices if b < len(match)]
    except Exception:
        pass
    try:
        res = rdFMCS.FindMCS(
            [orig_mol, mut_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            timeout=2,
        )
        if res and res.smartsString:
            mcs = Chem.MolFromSmarts(res.smartsString)
            if mcs is None:
                return None
            om = orig_mol.GetSubstructMatch(mcs)
            mm = mut_mol.GetSubstructMatch(mcs)
            if om and mm and len(om) == len(mm):
                om_index = {o_idx: i for i, o_idx in enumerate(om)}
                mapped = []
                for b in bind_indices:
                    if b in om_index:
                        mapped.append(mm[om_index[b]])
                return mapped if mapped else None
    except Exception:
        pass
    return None


def reassemble_complex_from_mols(
    metal_symbol, metal_charge, ligand_mols, ligand_binding_indices,
    expected_binding_symbols=None, use_dative_bonds=True,
):
    metal_atom = Chem.Atom(metal_symbol)
    metal_atom.SetFormalCharge(metal_charge)
    combo = Chem.RWMol()
    metal_idx = combo.AddAtom(metal_atom)

    for lig_i, (lig, bind_idxs) in enumerate(zip(ligand_mols, ligand_binding_indices)):
        if lig is None:
            return None
        if expected_binding_symbols and lig_i < len(expected_binding_symbols):
            expected_syms = expected_binding_symbols[lig_i]
            for j, bidx in enumerate(bind_idxs):
                if bidx >= lig.GetNumAtoms():
                    return None
                actual_sym = lig.GetAtomWithIdx(bidx).GetSymbol()
                if j < len(expected_syms) and actual_sym != expected_syms[j]:
                    return None

        offset = combo.GetNumAtoms()
        combo.InsertMol(lig)
        for bidx in bind_idxs:
            if use_dative_bonds:
                combo.AddBond(offset + bidx, metal_idx, Chem.BondType.DATIVE)
            else:
                combo.AddBond(metal_idx, offset + bidx, Chem.BondType.SINGLE)
                try:
                    a = combo.GetAtomWithIdx(offset + bidx)
                    a.SetNoImplicit(True)
                    a.SetNumExplicitHs(0)
                except Exception:
                    pass
    return combo.GetMol()


def _pick_mutation_smiles(mut):
    if isinstance(mut, str):
        return None, mut
    if not isinstance(mut, tuple) or len(mut) == 0:
        return None, None
    label = None
    for j, item in enumerate(mut[:3]):
        if isinstance(item, str):
            m = Chem.MolFromSmiles(item, sanitize=False)
            if m is not None and m.GetNumAtoms() > 0:
                return (None if j == 0 else mut[0]), item
            else:
                if j == 0:
                    label = item
    if len(mut) >= 3 and isinstance(mut[2], str):
        return (mut[0] if isinstance(mut[0], str) else None), mut[2]
    if len(mut) >= 2 and isinstance(mut[1], str):
        return (mut[0] if isinstance(mut[0], str) else None), mut[1]
    return (label, None)


def mutate_complex_ligands(complex_smiles, mutation_func, use_dative_bonds=True, **mutation_kwargs):
    result = extract_ligands_from_complex(complex_smiles)
    metal, charge, ligands, binding_maps, binding_symbols, error = result
    if error:
        raise ValueError(f"Could not extract ligands: {error}")
    if not ligands:
        raise ValueError("No ligands found in complex")

    orig_lig_mols = []
    orig_binding_indices = []
    for i, lig in enumerate(ligands):
        m = Chem.MolFromSmiles(lig, sanitize=False)
        if m is not None:
            try:
                m.UpdatePropertyCache(strict=False)
            except Exception:
                pass
        orig_lig_mols.append(m)
        if m is not None:
            bind_idx = _find_binding_atoms_by_map(m, binding_maps[i])
            orig_binding_indices.append(bind_idx if bind_idx else [])
        else:
            orig_binding_indices.append([])

    unique_ligands = {}
    for i, lig in enumerate(ligands):
        key = lig
        if key not in unique_ligands:
            unique_ligands[key] = []
        unique_ligands[key].append(i)

    results = []
    seen_complexes = set()

    for lig_smiles_with_maps, indices in unique_ligands.items():
        first_idx = indices[0]
        map_nums = binding_maps[first_idx]
        expected_symbols = binding_symbols[first_idx] if binding_symbols else None

        mol = Chem.MolFromSmiles(lig_smiles_with_maps, sanitize=False)
        if mol is None:
            continue
        try:
            mol.UpdatePropertyCache(strict=False)
        except Exception:
            pass
        if mol.GetNumHeavyAtoms() < 5:
            continue

        bind_indices = _find_binding_atoms_by_map(mol, map_nums)
        if not bind_indices:
            continue

        orig_signatures = {}
        for map_num, idx in zip(map_nums, bind_indices):
            orig_signatures[map_num] = _binding_atom_signature(mol.GetAtomWithIdx(idx))

        mapped_smiles = lig_smiles_with_maps

        mapped_smiles_arom = None
        try:
            mol_arom = Chem.Mol(mol)
            Chem.SanitizeMol(mol_arom, sanitizeOps=Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            mapped_smiles_arom = Chem.MolToSmiles(mol_arom, canonical=True)
        except Exception:
            pass

        mutations = None
        try:
            mutations = mutation_func(mapped_smiles, **mutation_kwargs)
        except Exception:
            mutations = None
        if not mutations and mapped_smiles_arom and mapped_smiles_arom != mapped_smiles:
            try:
                mutations = mutation_func(mapped_smiles_arom, **mutation_kwargs)
            except Exception:
                mutations = None
        if not mutations:
            try:
                mol_kek = Chem.MolFromSmiles(mapped_smiles_arom or mapped_smiles)
                if mol_kek:
                    Chem.Kekulize(mol_kek)
                    kek_smiles = Chem.MolToSmiles(mol_kek, kekuleSmiles=True)
                    mutations = mutation_func(kek_smiles, **mutation_kwargs)
            except Exception:
                pass
        if not mutations:
            continue

        for i, mut in enumerate(mutations):
            if isinstance(mut, (tuple, str)):
                lbl, smi = _pick_mutation_smiles(mut)
                if smi is None:
                    continue
                new_lig = smi
                label = lbl if lbl else f"mut_{i + 1}"
            else:
                continue

            mut_mol = Chem.MolFromSmiles(new_lig, sanitize=False)
            if mut_mol is None:
                continue
            try:
                mut_mol.UpdatePropertyCache(strict=False)
            except Exception:
                pass

            bind_indices_mut = _find_binding_atoms_by_map(mut_mol, map_nums)
            if bind_indices_mut is None or len(bind_indices_mut) != len(map_nums):
                mol_clean = Chem.Mol(mol)
                for a in mol_clean.GetAtoms():
                    a.SetAtomMapNum(0)
                mut_clean = Chem.Mol(mut_mol)
                for a in mut_clean.GetAtoms():
                    a.SetAtomMapNum(0)
                mapped = _map_binding_by_substructure(mol_clean, mut_clean, bind_indices)
                if mapped and len(mapped) == len(bind_indices):
                    bind_indices_mut = mapped
                    for map_num, idx in zip(map_nums, mapped):
                        mut_mol.GetAtomWithIdx(idx).SetAtomMapNum(map_num)
                else:
                    continue

            valid = True
            for map_num, new_idx in zip(map_nums, bind_indices_mut):
                new_sig = _binding_atom_signature(mut_mol.GetAtomWithIdx(new_idx))
                orig_sig = orig_signatures.get(map_num)
                if orig_sig and new_sig[0] != orig_sig[0]:
                    valid = False
                    break
            if not valid:
                continue

            clean_lig = Chem.Mol(mut_mol)
            for atom in clean_lig.GetAtoms():
                atom.SetAtomMapNum(0)

            new_lig_mols = []
            new_binding = []
            new_binding_symbols = []
            for j in range(len(ligands)):
                if j in indices:
                    new_lig_mols.append(Chem.Mol(clean_lig))
                    new_binding.append(list(bind_indices_mut))
                    new_binding_symbols.append(list(expected_symbols) if expected_symbols else [])
                else:
                    if orig_lig_mols[j] is not None:
                        clean_orig = Chem.Mol(orig_lig_mols[j])
                        for a in clean_orig.GetAtoms():
                            a.SetAtomMapNum(0)
                        new_lig_mols.append(clean_orig)
                    else:
                        new_lig_mols.append(None)
                    new_binding.append(
                        list(orig_binding_indices[j]) if orig_binding_indices[j] else []
                    )
                    new_binding_symbols.append(
                        list(binding_symbols[j]) if binding_symbols else []
                    )

            try:
                combo = reassemble_complex_from_mols(
                    metal, charge, new_lig_mols, new_binding,
                    expected_binding_symbols=new_binding_symbols,
                    use_dative_bonds=use_dative_bonds,
                )
                if combo is None:
                    continue
                canon_complex = Chem.MolToSmiles(combo)
                if canon_complex in seen_complexes:
                    continue
                seen_complexes.add(canon_complex)
                results.append((f"L:{label}", canon_complex))
            except Exception:
                continue
    return results


def get_complex_info(smiles):
    result = extract_ligands_from_complex(smiles)
    metal, charge, ligands, binding_maps, binding_symbols, error = result
    if error:
        return None
    clean_ligands = [_strip_atom_maps(lig) for lig in ligands]
    return {
        'metal': metal,
        'charge': charge,
        'ligands': clean_ligands,
        'binding_maps': binding_maps,
        'binding_symbols': binding_symbols,
        'num_ligands': len(ligands),
        'unique_ligands': len(set(clean_ligands)),
    }


# ===================================================================
# Drawing / visualisation helpers
# ===================================================================

def _precompute_2d_coords(mols):
    """Compute 2D coords for all mols in-place (call once after generation)."""
    for m in mols:
        if m is not None:
            try:
                rdDepictor.Compute2DCoords(m)
            except Exception:
                pass


def draw_grid(mols, legends=None, mols_per_row=4, size=(260, 260), coords_ready=False):
    valid = [m for m in mols if m is not None]
    if legends is not None:
        legends = [l for m, l in zip(mols, legends) if m is not None]
    if not coords_ready:
        for m in valid:
            try:
                rdDepictor.Compute2DCoords(m)
            except Exception:
                pass
    img = MolsToGridImage(
        valid, legends=legends, molsPerRow=mols_per_row,
        subImgSize=size, returnPNG=True,
    )
    return img


def _parse_complex_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return mol
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol, sanitizeOps=(
            Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
        ))
        return mol
    except Exception:
        pass
    try:
        mol.UpdatePropertyCache(strict=False)
        return mol
    except Exception:
        pass
    return Chem.MolFromSmiles(smiles, sanitize=False)


def smiles_to_3d_view(smiles, width=350, height=300):
    if py3Dmol is None:
        return None
    is_complex = contains_metal(smiles)
    mol = None
    if is_complex:
        mol = _parse_complex_mol(smiles)
        if mol is None:
            return None
        try:
            mol = Chem.AddHs(mol, addCoords=False)
        except Exception:
            pass
        try:
            params = AllChem.ETKDGv3()
            params.useRandomCoords = True
            params.randomSeed = 42
            params.maxIterations = 500
            conf_id = AllChem.EmbedMolecule(mol, params)
            if conf_id < 0:
                AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
        except Exception:
            try:
                rdDepictor.Compute2DCoords(mol)
                conf = mol.GetConformer()
                for i in range(mol.GetNumAtoms()):
                    pos = conf.GetAtomPosition(i)
                    conf.SetAtomPosition(i, (pos.x, pos.y, 0.0))
            except Exception:
                return None
    else:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                return None
        try:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            pass
    try:
        mol_block = Chem.MolToMolBlock(mol)
    except Exception:
        return None
    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(mol_block, 'mol')
    apply_molecule_view_style(viewer)
    return viewer


# ===================================================================
# Chemical-space plots
# ===================================================================

def _pca_2d(X):
    X = X.astype(float)
    X -= X.mean(axis=0, keepdims=True)
    U, S, _Vt = np.linalg.svd(X, full_matrices=False)
    return U[:, :2] * S[:2]


def _featurize_morgan(mols, n_bits=2048, radius=2):
    arr = np.zeros((len(mols), n_bits), dtype=float)
    for i, m in enumerate(mols):
        if m is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=n_bits)
        DataStructs.ConvertToNumpyArray(fp, arr[i])
    return arr


def _featurize_physchem(mols):
    feats = []
    for m in mols:
        if m is None:
            feats.append([0, 0, 0, 0, 0, 0])
            continue
        feats.append([
            Descriptors.MolWt(m),
            Descriptors.MolLogP(m),
            Descriptors.TPSA(m),
            Descriptors.NumHDonors(m),
            Descriptors.NumHAcceptors(m),
            Descriptors.NumRotatableBonds(m),
        ])
    return np.array(feats, dtype=float)


def _plot_space(coords, title, output_dir=None, filename=None):
    try:
        import matplotlib.pyplot as plt
    except Exception:
        display(Markdown('*Matplotlib not available — 2D space plot skipped.*'))
        return
    x, y = coords[:, 0], coords[:, 1]
    plt.figure(figsize=(5.2, 4.2))
    plt.scatter(x, y, s=18, alpha=0.7, edgecolors='none')
    plt.title(title)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.tight_layout()
    fig = plt.gcf()
    display(fig)
    if output_dir and filename:
        try:
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            fig.savefig(Path(output_dir) / filename, dpi=150)
        except Exception:
            pass
    plt.close()


def show_space_views(smiles_list, output_dir=None):
    import warnings

    mols = []
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m is None:
            m = Chem.MolFromSmiles(s, sanitize=False)
        mols.append(m)
    if not mols:
        return

    # -- sample for large datasets ----------------------------------------
    MAX_VIZ = 3000
    if len(mols) > MAX_VIZ:
        import random as _rng
        _rng.seed(42)
        idx = sorted(_rng.sample(range(len(mols)), MAX_VIZ))
        mols_viz = [mols[i] for i in idx]
        display(Markdown(f'*Sampling {MAX_VIZ} of {len(mols)} molecules for visualization.*'))
    else:
        mols_viz = mols

    # -- featurize upfront (shared by multiple plots) --------------------
    X_fp = _featurize_morgan(mols_viz)
    X_phys = _featurize_physchem(mols_viz)

    # -- fast PCA plots (sequential, cheap) ------------------------------
    display(Markdown('**Morgan FP PCA (2D)**'))
    _plot_space(_pca_2d(X_fp), 'Morgan FP PCA (2D)', output_dir, 'morgan_fp_pca_2d.png')
    display(Markdown('Projection of Morgan fingerprints; nearby points indicate similar substructure patterns.'))

    display(Markdown('**PhysChem PCA (2D)** – MolWt, LogP, TPSA, HBD/HBA, RotB.'))
    _plot_space(_pca_2d(X_phys), 'PhysChem PCA (2D)', output_dir, 'physchem_pca_2d.png')
    display(Markdown('Projection of basic physicochemical properties to show overall property diversity.'))

    # -- expensive computations: keep bounded so large spaces stay responsive
    n_mols = len(mols_viz)
    HEAVY_THRESHOLD = 100
    HEAVY_SKIP_THRESHOLD = 1500
    UMAP_MAX = 1200
    TSNE_MAX = 500
    SCAFFOLD_MAX = 1500
    SIMILARITY_MAX = 150

    if n_mols > HEAVY_SKIP_THRESHOLD:
        display(Markdown(
            f'*Large space detected ({n_mols} molecules). '
            'Skipping UMAP, t-SNE, scaffold PCA, and similarity heatmap to keep the tab responsive.*'
        ))

    def _compute_umap():
        if n_mols < HEAVY_THRESHOLD or n_mols > HEAVY_SKIP_THRESHOLD:
            return None
        try:
            import umap as _umap
            n = min(len(X_fp), UMAP_MAX)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', module='umap')
                reducer = _umap.UMAP(n_neighbors=15, min_dist=0.1, metric='jaccard', random_state=0)
                return reducer.fit_transform(X_fp[:n])
        except Exception:
            return None

    def _compute_tsne():
        if n_mols < HEAVY_THRESHOLD or n_mols > HEAVY_SKIP_THRESHOLD:
            return None
        try:
            from sklearn.manifold import TSNE
            n = min(len(X_fp), TSNE_MAX)
            perp = min(30, max(2, n - 1))
            tsne = TSNE(n_components=2, init='random', learning_rate='auto', perplexity=perp, random_state=0)
            return tsne.fit_transform(X_fp[:n])
        except Exception:
            return None

    def _compute_scaffold():
        if n_mols > HEAVY_SKIP_THRESHOLD:
            return None
        try:
            from rdkit.Chem.Scaffolds import MurckoScaffold
            scaff_mols = []
            for m in mols_viz[:SCAFFOLD_MAX]:
                if m is None:
                    scaff_mols.append(None)
                else:
                    try:
                        scaff_mols.append(MurckoScaffold.GetScaffoldForMol(m))
                    except Exception:
                        scaff_mols.append(None)
            return _pca_2d(_featurize_morgan(scaff_mols))
        except Exception:
            return None

    def _compute_similarity():
        try:
            if n_mols > HEAVY_SKIP_THRESHOLD:
                return None
            n = min(len(mols_viz), SIMILARITY_MAX)
            fps = [
                AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m is not None else None
                for m in mols_viz[:n]
            ]
            sim = np.zeros((n, n), dtype=float)
            for i in range(n):
                if fps[i] is None:
                    continue
                bulk = DataStructs.BulkTanimotoSimilarity(fps[i], [f for f in fps if f is not None])
                k = 0
                for j in range(n):
                    if fps[j] is not None:
                        sim[i, j] = bulk[k]
                        k += 1
            return sim
        except Exception:
            return None

    # This function already runs in a background thread from the UI. Running
    # another local thread pool here tends to oversubscribe CPU and makes
    # large result spaces look stuck, so keep the heavy work sequential.
    umap_coords = _compute_umap()
    if umap_coords is not None:
        display(Markdown('**UMAP (Morgan FP)**'))
        _plot_space(umap_coords, 'UMAP (Morgan FP)', output_dir, 'umap_morgan_fp.png')
        display(Markdown('Nonlinear embedding of fingerprints; clusters indicate structural families.'))

    tsne_coords = _compute_tsne()
    if tsne_coords is not None:
        display(Markdown('**t-SNE (Morgan FP)**'))
        _plot_space(tsne_coords, 't-SNE (Morgan FP)', output_dir, 'tsne_morgan_fp.png')
        display(Markdown('Nonlinear embedding emphasizing local neighborhoods in fingerprint space.'))

    scaff_coords = _compute_scaffold()
    if scaff_coords is not None:
        display(Markdown('**Scaffold Space (Bemis–Murcko)**'))
        _plot_space(scaff_coords, 'Scaffold PCA (2D)', output_dir, 'scaffold_pca_2d.png')
        display(Markdown('Scaffold-level projection to compare core frameworks.'))

    try:
        import matplotlib.pyplot as plt
        display(Markdown('**Property Plots**'))
        mw, logp, tpsa = X_phys[:, 0], X_phys[:, 1], X_phys[:, 2]
        hbd, hba, rot = X_phys[:, 3], X_phys[:, 4], X_phys[:, 5]
        for x, y, xl, yl in [(mw, logp, 'MolWt', 'LogP'), (tpsa, hbd, 'TPSA', 'HBD'), (hba, rot, 'HBA', 'RotB')]:
            plt.figure(figsize=(5.2, 4.2))
            plt.scatter(x, y, s=18, alpha=0.7, edgecolors='none')
            plt.xlabel(xl); plt.ylabel(yl)
            plt.title(f'{xl} vs {yl}')
            plt.tight_layout()
            fig = plt.gcf()
            display(fig)
            if output_dir:
                try:
                    Path(output_dir).mkdir(parents=True, exist_ok=True)
                    fig.savefig(Path(output_dir) / f'property_{xl.lower()}_{yl.lower()}.png', dpi=150)
                except Exception:
                    pass
            plt.close()
        display(Markdown('Pairwise property relationships to spot trends and outliers.'))
    except Exception:
        display(Markdown('*Property plots not available — skipped.*'))

    sim_matrix = _compute_similarity()
    if sim_matrix is not None:
        try:
            import matplotlib.pyplot as plt
            display(Markdown('**Similarity Heatmap**'))
            plt.figure(figsize=(5.4, 4.6))
            plt.imshow(sim_matrix, cmap='viridis', vmin=0.0, vmax=1.0)
            plt.colorbar(label='Tanimoto')
            plt.title('Similarity Heatmap (sample)')
            plt.tight_layout()
            fig = plt.gcf()
            display(fig)
            if output_dir:
                try:
                    Path(output_dir).mkdir(parents=True, exist_ok=True)
                    fig.savefig(Path(output_dir) / 'similarity_heatmap.png', dpi=150)
                except Exception:
                    pass
            plt.close()
            display(Markdown('Pairwise similarity matrix; warmer colors indicate more similar molecules.'))
        except Exception:
            display(Markdown('*Similarity heatmap not available — skipped.*'))


def _space_views_gallery_html(output_dir):
    out_dir = Path(output_dir)
    if not out_dir.exists():
        return '<div style="color:#777;">No chemical-space images were generated.</div>'

    descriptions = {
        'morgan_fp_pca_2d.png': (
            'Morgan fingerprint projection; nearby points indicate similar '
            'substructure patterns.'
        ),
        'physchem_pca_2d.png': (
            'Projection of basic physicochemical properties to show overall '
            'property diversity.'
        ),
        'umap_morgan_fp.png': (
            'Nonlinear embedding of fingerprints; clusters indicate structural '
            'families.'
        ),
        'tsne_morgan_fp.png': (
            'Nonlinear embedding emphasizing local neighborhoods in fingerprint '
            'space.'
        ),
        'scaffold_pca_2d.png': (
            'Scaffold-level projection to compare core frameworks.'
        ),
        'similarity_heatmap.png': (
            'Pairwise similarity matrix; warmer colors indicate more similar '
            'molecules.'
        ),
    }

    preferred = [
        'morgan_fp_pca_2d.png',
        'physchem_pca_2d.png',
        'umap_morgan_fp.png',
        'tsne_morgan_fp.png',
        'scaffold_pca_2d.png',
        'similarity_heatmap.png',
    ]
    files = []
    for name in preferred:
        path = out_dir / name
        if path.exists():
            files.append(path)
    files.extend(sorted(out_dir.glob('property_*.png')))

    if not files:
        return '<div style="color:#777;">No chemical-space images were generated.</div>'

    blocks = []
    for path in files:
        try:
            payload = base64.b64encode(path.read_bytes()).decode('ascii')
        except Exception:
            continue
        label = path.stem.replace('_', ' ')
        description = descriptions.get(path.name, '')
        if path.name.startswith('property_'):
            description = 'Pairwise property relationship plot to spot trends and outliers.'
        blocks.append(
            '<div style="margin:0 0 18px 0;">'
            f'<div style="font-weight:600; margin:0 0 6px 0;">{label}</div>'
            f'<div style="color:#555; margin:0 0 8px 0;">{description}</div>'
            f'<img src="data:image/png;base64,{payload}" '
            'style="max-width:100%; border:1px solid #ddd; display:block;">'
            '</div>'
        )
    return ''.join(blocks) if blocks else '<div style="color:#777;">No chemical-space images were generated.</div>'


# ===================================================================
# Reaction engine (Custom SMARTS)
# ===================================================================

def _parse_forbidden_line(line):
    line = line.strip()
    if not line or line.lower() in ('-', 'none', 'n/a', '.', '*'):
        return []
    patterns = [p.strip() for p in line.split(';') if p.strip()]
    return [p for p in patterns if p.lower() not in ('-', 'none', 'n/a', '.', '*')]


def _compile_smarts_list(smarts_list):
    """Pre-compile a list of SMARTS strings into RDKit pattern objects (once)."""
    compiled = []
    for smarts in smarts_list:
        smarts = smarts.strip()
        if not smarts:
            continue
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None:
            compiled.append((smarts, pattern))
    return compiled


def check_forbidden_patterns(seed_mol, product_mol, forbidden_smarts_list, debug=False, _compiled=None):
    if not forbidden_smarts_list and not _compiled:
        return True
    patterns = _compiled if _compiled is not None else _compile_smarts_list(forbidden_smarts_list)
    for smarts, pattern in patterns:
        try:
            count_before = len(seed_mol.GetSubstructMatches(pattern))
        except Exception:
            count_before = 0
        try:
            count_after = len(product_mol.GetSubstructMatches(pattern))
        except Exception:
            count_after = 0
        if count_after > count_before:
            return False
    return True


def check_protected_patterns(seed_mol, product_mol, protected_smarts_list, debug=False, _compiled=None):
    if not protected_smarts_list and not _compiled:
        return True
    patterns = _compiled if _compiled is not None else _compile_smarts_list(protected_smarts_list)
    for smarts, pattern in patterns:
        try:
            count_before = len(seed_mol.GetSubstructMatches(pattern))
        except Exception:
            count_before = 0
        try:
            count_after = len(product_mol.GetSubstructMatches(pattern))
        except Exception:
            count_after = 0
        if count_after < count_before:
            return False
    return True


def apply_custom_reaction_iter(
    seed_smiles, rxn_smarts, iterations=1, keep_rings=True,
    forbidden_patterns=None, protected_patterns=None, debug=False,
):
    base = Chem.MolFromSmiles(seed_smiles)
    if base is None:
        base = Chem.MolFromSmiles(seed_smiles, sanitize=False)
        if base is not None:
            try:
                base.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(base, sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                ))
            except Exception:
                try:
                    base.UpdatePropertyCache(strict=False)
                except Exception:
                    pass
    if base is None:
        raise ValueError('Invalid seed SMILES.')

    raw_lines = [s.strip() for s in rxn_smarts.splitlines() if s.strip()]
    if not raw_lines:
        raise ValueError('Could not read reaction SMARTS.')

    smarts_list = []
    for line in raw_lines:
        if '>>' not in line and ':' in line:
            line = line.split(':', 1)[1].strip()
        if '>>' in line:
            smarts_list.append(line)
    if not smarts_list:
        raise ValueError('Reaction SMARTS must be in A>>B format (one per line).')

    # Parse & pre-compile forbidden patterns (once, not per product)
    forbidden_compiled = []
    if forbidden_patterns:
        for line in forbidden_patterns.splitlines():
            forbidden_compiled.append(_compile_smarts_list(_parse_forbidden_line(line)))
    while len(forbidden_compiled) < len(smarts_list):
        forbidden_compiled.append([])

    protected_compiled = []
    if protected_patterns:
        for line in protected_patterns.splitlines():
            protected_compiled.append(_compile_smarts_list(_parse_forbidden_line(line)))
    while len(protected_compiled) < len(smarts_list):
        protected_compiled.append([])

    if debug:
        print(f"[DEBUG] Seed SMILES: {seed_smiles}")
        print(f"[DEBUG] Reaction SMARTS: {smarts_list}")

    rxns = []
    for smarts in smarts_list:
        try:
            rxn = rdChemReactions.ReactionFromSmarts(smarts)
        except Exception:
            rxn = None
        if rxn is not None:
            rxns.append(rxn)
    if not rxns:
        raise ValueError('Could not read reaction SMARTS.')

    try:
        input_rings = base.GetRingInfo().NumRings()
        input_aromatic = _count_aromatic_rings(base)
    except Exception:
        input_rings = 0
        input_aromatic = 0
        keep_rings = False

    seen = set()
    out = []
    filtered_count = 0
    protected_filtered_count = 0

    current = [base]
    for _iter_num in range(max(1, int(iterations))):
        next_round = []
        for mol in current:
            for rxn_idx, rxn in enumerate(rxns):
                forb_comp = forbidden_compiled[rxn_idx] if rxn_idx < len(forbidden_compiled) else []
                prot_comp = protected_compiled[rxn_idx] if rxn_idx < len(protected_compiled) else []
                try:
                    ps = rxn.RunReactants((Chem.Mol(mol),))
                except Exception:
                    continue
                for p in ps:
                    m = p[0]
                    try:
                        Chem.SanitizeMol(m)
                    except Exception:
                        continue
                    # Dedup early via canonical SMILES (much faster than InchiKey)
                    s = Chem.MolToSmiles(m, canonical=True, isomericSmiles=False)
                    if s in seen:
                        continue
                    if forb_comp:
                        if not check_forbidden_patterns(mol, m, None, debug=debug, _compiled=forb_comp):
                            filtered_count += 1
                            continue
                    if prot_comp:
                        if not check_protected_patterns(mol, m, None, debug=debug, _compiled=prot_comp):
                            protected_filtered_count += 1
                            continue
                    if keep_rings:
                        try:
                            if (m.GetRingInfo().NumRings() != input_rings
                                    or _count_aromatic_rings(m) != input_aromatic):
                                continue
                        except Exception:
                            pass
                    seen.add(s)
                    out.append((s, m))
                    next_round.append(m)
        current = next_round
        if not current:
            break

    if debug:
        print(f"\n[DEBUG] Filtered by forbidden patterns: {filtered_count}")
        print(f"[DEBUG] Filtered by protected patterns: {protected_filtered_count}")
        print(f"[DEBUG] Results: {len(out)}")

    if not out:
        if filtered_count > 0:
            raise ValueError(
                f'All {filtered_count} products were filtered by forbidden patterns.'
            )
        raise ValueError('Reaction SMARTS produced no products.')
    return out


# ===================================================================
# Tab factory
# ===================================================================

def create_tab(ctx):
    """Build the ChemDarwin tab and return ``(tab_widget, refs_dict)``."""

    VIZ_SIZE = 520
    IMG_SIZE = VIZ_SIZE
    PAGE_SIZE = 50

    # -- mutable state shared by closures --------------------------------
    state = {
        'last_mols': [],
        'last_legends': [],
        'page_idx': 0,
        'page_cache': {},        # page_idx -> rendered PIL image
        'coords_ready': False,   # True after _precompute_2d_coords
    }

    # -- widgets ---------------------------------------------------------
    job_name_input = widgets.Text(
        placeholder='Enter job name (required)',
        layout=widgets.Layout(width='95%'),
    )

    seed_input = widgets.Textarea(
        placeholder='Enter seed SMILES (organic or metal complex)',
        layout=widgets.Layout(width='95%', height='80px'),
    )

    complex_info_label = widgets.HTML(value='')
    seed_2d_out = widgets.Output(layout=widgets.Layout(width='100%', height='100%', margin='0', padding='0'))
    seed_2d_out.add_class('chemdarwin-seed-out')
    seed_3d_out = widgets.Output(layout=widgets.Layout(width='100%', height='100%', margin='0', padding='0'))
    seed_3d_out.add_class('chemdarwin-seed-out')

    custom_smarts = widgets.Textarea(
        placeholder='Reaction SMARTS (one per line)',
        layout=widgets.Layout(width='95%', height='240px'),
    )

    custom_forbidden = widgets.Textarea(
        placeholder=(
            'Forbidden patterns per line (separate multiple with ;)\n'
            'Example: [#8]~c~[#8];[N+]\n'
            'Only newly formed patterns are blocked.'
        ),
        layout=widgets.Layout(width='95%', height='240px'),
    )
    custom_forbidden_label = widgets.HTML(
        '<small style="color:#666;">One reaction SMARTS per line. '
        'Separate multiple patterns with ;. Only newly formed patterns are blocked.</small>'
    )

    custom_protected = widgets.Textarea(
        placeholder=(
            'Protected patterns per line (separate multiple with ;)\n'
            'These patterns must remain unchanged.'
        ),
        layout=widgets.Layout(width='95%', height='240px'),
    )
    custom_protected_label = widgets.HTML(
        '<small style="color:#666;">One reaction SMARTS per line. '
        'Separate multiple patterns with ;. Existing patterns must be preserved.</small>'
    )

    custom_keep_rings = widgets.Checkbox(value=False, description='Preserve ring count')
    custom_iters = widgets.IntText(value=1, description='Iterations', layout=widgets.Layout(width='30%'))
    custom_iters.add_class('chemdarwin-iter')
    custom_debug = widgets.Checkbox(value=False, description='Debug (show pattern matches)')
    custom_space_views = widgets.Checkbox(value=True, description='Generate space views (slower)')

    run_btn = widgets.Button(description='Run', button_style='primary')
    out = widgets.Output(layout=widgets.Layout(width='100%', overflow_x='hidden'))
    space_out = widgets.HTML(value='', layout=widgets.Layout(width='100%', overflow_x='hidden', display='none'))

    grid_out = widgets.Output(layout=widgets.Layout(width='100%', overflow_x='hidden'))
    page_label = widgets.HTML('')
    btn_prev_page = widgets.Button(description='Prev 50', button_style='')
    btn_next_page = widgets.Button(description='Next 50', button_style='')
    page_jump_input = widgets.IntText(
        value=1, min=1, description='Page:',
        layout=widgets.Layout(width='120px'),
    )
    page_jump_input.add_class('chemdarwin-iter')
    btn_jump = widgets.Button(description='Go', button_style='info', layout=widgets.Layout(width='40px'))

    def _safe_job_name(name_raw):
        return ''.join(c for c in name_raw if c.isalnum() or c in ('_', '-'))

    job_status_label = widgets.HTML(value='')

    def _load_existing_job(safe_name):
        """Load metadata from an existing job and populate widgets."""
        job_dir = Path(ctx.calc_dir) / safe_name
        meta_path = job_dir / 'run_metadata.json'
        if not meta_path.exists():
            return False
        try:
            meta = json.loads(meta_path.read_text(encoding='utf-8'))
        except Exception:
            return False
        seed_input.value = meta.get('seed_smiles', '')
        custom_smarts.value = meta.get('reaction_smarts', '')
        custom_forbidden.value = meta.get('forbidden_patterns', '')
        custom_protected.value = meta.get('protected_patterns', '')
        custom_iters.value = int(meta.get('iterations', 1))
        custom_keep_rings.value = bool(meta.get('preserve_ring_count', False))
        custom_debug.value = bool(meta.get('debug', False))
        result_count = meta.get('result_count', '?')
        job_status_label.value = (
            '<div style="background:#e8f5e9; padding:6px 10px; border-radius:4px; margin:4px 0;">'
            f'<b>Job loaded</b> – {result_count} previous results. '
            'Modify parameters and click <b>Run</b> to re-run.</div>'
        )
        return True

    def _update_run_enabled(_change=None):
        raw = (job_name_input.value or '').strip()
        safe = _safe_job_name(raw)
        run_btn.disabled = not bool(safe)
        if safe:
            loaded = _load_existing_job(safe)
            if not loaded:
                job_status_label.value = ''
        else:
            job_status_label.value = ''

    job_name_input.observe(_update_run_enabled, names='value')
    _update_run_enabled()

    # -- seed visualisation ----------------------------------------------
    def update_seed_visualization(_change=None):
        smiles = (seed_input.value or '').strip()
        is_complex = contains_metal(smiles) if smiles else False
        if is_complex:
            info = get_complex_info(smiles)
            if info:
                complex_info_label.value = (
                    '<div style="background:#e3f2fd; padding:8px; border-radius:4px; margin:5px 0;">'
                    f'<b>Metal complex detected:</b> [{info["metal"]}] with {info["num_ligands"]} ligands '
                    f'({info["unique_ligands"]} unique)<br>'
                    f'<small>Ligands: {", ".join(info["ligands"][:5])}'
                    f'{"..." if len(info["ligands"]) > 5 else ""}</small></div>'
                )
            else:
                complex_info_label.value = (
                    '<div style="background:#fff3e0; padding:8px; border-radius:4px; margin:5px 0;">'
                    '<b>Metal complex:</b> Could not extract ligands</div>'
                )
        else:
            complex_info_label.value = ''

        with seed_2d_out:
            clear_output(wait=True)
            if smiles:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    mol = Chem.MolFromSmiles(smiles, sanitize=False)
                if mol is not None:
                    try:
                        rdDepictor.Compute2DCoords(mol)
                    except Exception:
                        pass
                    img = MolToImage(mol, size=(IMG_SIZE, IMG_SIZE))
                    display(img)
                else:
                    display(Markdown('*Invalid SMILES*'))

        with seed_3d_out:
            clear_output(wait=True)
            if not smiles:
                display(HTML(
                    '<div style="width:100%;height:100%;display:flex;align-items:center;'
                    'justify-content:center;color:#777;">No 3D preview.</div>'
                ))
                return
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                try:
                    viewer = smiles_to_3d_view(smiles, width=VIZ_SIZE, height=VIZ_SIZE)
                    if viewer:
                        viewer.show()
                    else:
                        display(HTML(
                            '<div style="width:100%;height:100%;display:flex;align-items:center;'
                            'justify-content:center;color:#777;">3D preview unavailable.</div>'
                        ))
                except Exception as e:
                    display(Markdown(f'*3D error: {e}*'))
            else:
                display(Markdown('*Invalid SMILES*'))

    seed_input.observe(update_seed_visualization, names='value')

    seed_2d_wrapper = widgets.Box(
        [seed_2d_out],
        layout=widgets.Layout(
            width=f'{VIZ_SIZE}px', height=f'{VIZ_SIZE}px',
            max_width='100%',
            border='1px solid #ccc', overflow='hidden',
            margin='0', padding='0',
            display='flex', align_items='center', justify_content='center',
        ),
    )
    seed_2d_wrapper.add_class('chemdarwin-seed-wrapper')
    seed_3d_wrapper = widgets.Box(
        [seed_3d_out],
        layout=widgets.Layout(
            width=f'{VIZ_SIZE}px', height=f'{VIZ_SIZE}px',
            max_width='100%',
            border='1px solid #ccc', overflow='hidden',
            margin='0', padding='0',
            display='flex', align_items='center', justify_content='center',
        ),
    )
    seed_3d_wrapper.add_class('chemdarwin-seed-wrapper')
    seed_2d_box = widgets.VBox([widgets.HTML('<b>2D Structure</b>'), seed_2d_wrapper])
    seed_3d_box = widgets.VBox([widgets.HTML('<b>3D Structure</b>'), seed_3d_wrapper])
    seed_viz_row = widgets.HBox(
        [seed_2d_box, seed_3d_box],
        layout=widgets.Layout(gap='20px', overflow_x='hidden', flex_flow='row wrap'),
    )

    seed_section = widgets.VBox([
        widgets.HTML('<b>Seed SMILES</b>'),
        seed_input,
        complex_info_label,
        seed_viz_row,
    ], layout=widgets.Layout(overflow_x='hidden'))

    # -- grid pagination -------------------------------------------------
    def _render_grid_page():
        with grid_out:
            clear_output()
            mols = state['last_mols']
            legends = state['last_legends']
            if not mols:
                return
            pi = state['page_idx']
            start = pi * PAGE_SIZE
            end = start + PAGE_SIZE
            total = len(mols)
            max_pages = max(1, (total + PAGE_SIZE - 1) // PAGE_SIZE)

            if pi in state['page_cache']:
                img = state['page_cache'][pi]
            else:
                sub_mols = mols[start:end]
                sub_legends = legends[start:end] if legends else None
                img = ''
                if sub_mols:
                    img = draw_grid(
                        sub_mols, legends=sub_legends,
                        mols_per_row=4, size=(260, 260),
                        coords_ready=state['coords_ready'],
                    )
                    state['page_cache'][pi] = img

            display(img)

            page_label.value = (
                f'<b>Page {pi + 1} / {max_pages} '
                f'(showing {min(end, total)}/{total})</b>'
            )
            page_jump_input.max = max_pages
            page_jump_input.value = pi + 1

    def _on_prev_page(_):
        if state['page_idx'] > 0:
            state['page_idx'] -= 1
            _render_grid_page()

    def _on_next_page(_):
        if state['last_mols']:
            max_page = (len(state['last_mols']) - 1) // PAGE_SIZE
            if state['page_idx'] < max_page:
                state['page_idx'] += 1
                _render_grid_page()

    def _on_jump_page(_):
        if state['last_mols']:
            max_page = (len(state['last_mols']) - 1) // PAGE_SIZE
            target = max(0, min(page_jump_input.value - 1, max_page))
            if target != state['page_idx']:
                state['page_idx'] = target
                _render_grid_page()

    btn_prev_page.on_click(_on_prev_page)
    btn_next_page.on_click(_on_next_page)
    btn_jump.on_click(_on_jump_page)

    # -- run button handler ----------------------------------------------
    def on_run(_btn):
        with out:
            clear_output()
        with grid_out:
            clear_output()
        pagination_bar.layout.display = 'none'
        space_out.layout.display = 'none'
        space_out.value = ''
        with out:
            job_name_raw = (job_name_input.value or '').strip()
            safe_job_name = _safe_job_name(job_name_raw)
            if not safe_job_name:
                display(Markdown('**Error:** Please enter a valid job name.'))
                return

            seed = (seed_input.value or '').strip()
            if not seed:
                display(Markdown('**Error:** Please enter a seed SMILES.'))
                return

            rxn_smarts_text = custom_smarts.value.strip()
            if not rxn_smarts_text:
                display(Markdown('**Error:** Please enter reaction SMARTS.'))
                return

            forbidden_text = custom_forbidden.value.strip() or None
            protected_text = custom_protected.value.strip() or None
            debug_mode = custom_debug.value
            is_complex = contains_metal(seed)

            if is_complex:
                display(Markdown(
                    '**Complex mode:** Mutating ligands and rebuilding the complex...'
                ))

            results = []
            legends = []
            mols = []

            if is_complex:
                def _custom_wrapper(smi):
                    return apply_custom_reaction_iter(
                        smi, rxn_smarts_text,
                        iterations=custom_iters.value,
                        keep_rings=custom_keep_rings.value,
                        forbidden_patterns=forbidden_text,
                        protected_patterns=protected_text,
                        debug=debug_mode,
                    )
                try:
                    res = mutate_complex_ligands(seed, _custom_wrapper)
                    results = res
                    legends = [lbl for lbl, _ in res]
                    mols = [Chem.MolFromSmiles(s, sanitize=False) for _, s in res]
                except Exception as e:
                    display(Markdown(f'**Error:** {e}'))
                    return
            else:
                try:
                    res = apply_custom_reaction_iter(
                        seed, rxn_smarts_text,
                        iterations=custom_iters.value,
                        keep_rings=custom_keep_rings.value,
                        forbidden_patterns=forbidden_text,
                        protected_patterns=protected_text,
                        debug=debug_mode,
                    )
                    results = res
                    legends = [f'prod_{i + 1}' for i in range(len(res))]
                    mols = [m for _, m in res]
                except Exception as e:
                    display(Markdown(f'**Error:** {e}'))
                    return

            mols = [m for m in mols if m is not None]
            if not results:
                display(Markdown('**No results.**'))
                return

            display(Markdown(f'**Result count:** {len(results)}'))

            job_dir = Path(ctx.calc_dir) / safe_job_name
            job_dir.mkdir(parents=True, exist_ok=True)

            # Pre-compute 2D coords once for all molecules
            display(Markdown('*Computing 2D layouts...*'))
            _precompute_2d_coords(mols)

            state['last_mols'] = mols
            state['last_legends'] = legends
            state['page_idx'] = 0
            state['page_cache'] = {}
            state['coords_ready'] = True

            if mols:
                pagination_bar.layout.display = ''
                _render_grid_page()

            # Build CSV lines
            csv_lines = []
            counts = {}
            header_lines = ['Label;SMILES']
            if is_complex:
                for lbl, s in results:
                    counts[lbl] = counts.get(lbl, 0) + 1
                    entry = f'{lbl}_{counts[lbl]};{s}'
                    header_lines.append(entry)
                    csv_lines.append(entry)
            else:
                for s, _ in results:
                    key = 'prod'
                    counts[key] = counts.get(key, 0) + 1
                    entry = f'{key}_{counts[key]};{s}'
                    header_lines.append(entry)
                    csv_lines.append(entry)

            display(widgets.Textarea(
                value='\n'.join(header_lines),
                layout=widgets.Layout(width='95%', height='200px'),
            ))

            csv_path = job_dir / 'complete_mutation_space.csv'
            csv_path.write_text('\n'.join(header_lines), encoding='utf-8')

            metadata = {
                'job_name_raw': job_name_raw,
                'job_name_safe': safe_job_name,
                'seed_smiles': seed,
                'reaction_smarts': rxn_smarts_text,
                'forbidden_patterns': forbidden_text or '',
                'protected_patterns': protected_text or '',
                'iterations': custom_iters.value,
                'preserve_ring_count': custom_keep_rings.value,
                'debug': debug_mode,
                'is_complex': bool(is_complex),
                'result_count': len(results),
            }
            meta_path = job_dir / 'run_metadata.json'
            meta_path.write_text(json.dumps(metadata, indent=2), encoding='utf-8')
            meta_txt_lines = [
                f'Job Name: {job_name_raw}',
                f'Seed SMILES: {seed}',
                'Reaction SMARTS:',
                rxn_smarts_text,
                'Forbidden Patterns:',
                forbidden_text or '',
                'Protected Patterns:',
                protected_text or '',
                f'Iterations: {custom_iters.value}',
                f'Preserve Ring Count: {custom_keep_rings.value}',
                f'Debug: {debug_mode}',
                f'Complex Mode: {bool(is_complex)}',
                f'Result Count: {len(results)}',
            ]
            meta_txt_path = job_dir / 'run_metadata.txt'
            meta_txt_path.write_text('\n'.join(meta_txt_lines), encoding='utf-8')

            if custom_space_views.value:
                smiles_for_plots = [ln.split(';', 1)[1] for ln in csv_lines]

                # Run space views in background so UI stays responsive
                import threading
                space_out.layout.display = ''
                space_out.value = '<div style="color:#555;font-style:italic;">Generating space views in background...</div>'

                def _run_space_bg(smi_list, out_dir, output_widget):
                    try:
                        show_space_views(smi_list, output_dir=out_dir)
                        output_widget.value = _space_views_gallery_html(out_dir)
                    except Exception as e:
                        output_widget.value = (
                            f'<div style="color:#d32f2f;">Space plot error: {e}</div>'
                        )

                t = threading.Thread(
                    target=_run_space_bg,
                    args=(smiles_for_plots, job_dir, space_out),
                    daemon=True,
                )
                t.start()

    run_btn.on_click(on_run)

    # -- layout ----------------------------------------------------------
    param_box = widgets.VBox([
        widgets.HTML('<b>Reaction SMARTS</b>'),
        custom_smarts,
        widgets.HTML('<b>Forbidden patterns (block newly formed)</b>'),
        custom_forbidden,
        custom_forbidden_label,
        widgets.HTML('<b>Protected patterns (keep unchanged)</b>'),
        custom_protected,
        custom_protected_label,
        custom_keep_rings,
        custom_iters,
        custom_debug,
        custom_space_views,
    ], layout=widgets.Layout(overflow_x='hidden'))

    chemdarwin_css = widgets.HTML(
        '<style>'
        '.chemdarwin-tab, .chemdarwin-tab * { box-sizing: border-box; max-width: 100% !important; }'
        '.chemdarwin-tab { overflow-x: hidden !important; min-width: 0 !important; }'
        '.chemdarwin-tab .widget-text:not(.chemdarwin-iter), .chemdarwin-tab .widget-textarea { '
        'width: 95% !important; min-width: 0 !important; }'
        '.chemdarwin-tab .chemdarwin-iter { width: 30% !important; min-width: 0 !important; }'
        '.chemdarwin-tab .widget-text input, .chemdarwin-tab .widget-textarea textarea { '
        'width: 100% !important; max-width: 100% !important; }'
        '.chemdarwin-tab .widget-hbox, .chemdarwin-tab .widget-vbox, .chemdarwin-tab .widget-box { '
        'max-width: 100% !important; min-width: 0 !important; }'
        '.chemdarwin-tab .widget-hbox { flex-wrap: wrap !important; }'
        '.chemdarwin-tab .jupyter-widgets-output-area { max-width: 100% !important; min-width: 0 !important; }'
        '.chemdarwin-tab .chemdarwin-seed-wrapper { margin:0 !important; padding:0 !important; }'
        '.chemdarwin-tab .chemdarwin-seed-out { width:100% !important; height:100% !important;'
        ' margin:0 !important; padding:0 !important; overflow:hidden !important; }'
        '.chemdarwin-tab .chemdarwin-seed-out .output_area,'
        ' .chemdarwin-tab .chemdarwin-seed-out .output_subarea,'
        ' .chemdarwin-tab .chemdarwin-seed-out .output_wrapper,'
        ' .chemdarwin-tab .chemdarwin-seed-out .jp-OutputArea-child,'
        ' .chemdarwin-tab .chemdarwin-seed-out .jp-OutputArea-output {'
        ' margin:0 !important; padding:0 !important; max-width:100% !important;'
        ' width:100% !important; height:100% !important; overflow:hidden !important; }'
        '</style>'
    )

    pagination_bar = widgets.HBox([
        btn_prev_page, btn_next_page,
        page_jump_input, btn_jump, page_label,
    ], layout=widgets.Layout(align_items='center', gap='4px', display='none'))

    tab_body = widgets.VBox([
        chemdarwin_css,
        widgets.HTML('<h3 style="margin:0 0 8px 0;">ChemDarwin</h3>'),
        widgets.HTML('<b>Job Name</b>'),
        job_name_input,
        job_status_label,
        seed_section,
        param_box,
        run_btn,
        out,
        pagination_bar,
        grid_out,
        space_out,
    ], layout=widgets.Layout(overflow_x='hidden'))
    tab_body.add_class('chemdarwin-tab')

    refs = {
        'seed_input': seed_input,
        'custom_smarts': custom_smarts,
        'run_btn': run_btn,
    }

    return tab_body, refs

