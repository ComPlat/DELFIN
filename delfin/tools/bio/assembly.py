"""Explicit single-pose PDB assembly with a lossless ligand sidecar."""
from __future__ import annotations

from pathlib import Path
import re

from delfin.tools.bio import read_structure
from delfin.tools.bio.artifacts import file_hash, validate_artifact
from delfin.tools.bio.processing import load_ligands, select_ligand


def assemble_complex(receptor_file: Path, poses_file: Path, output_dir: Path, *,
                     record_id: int, shared_coordinate_frame: bool, ligand_chain: str = 'L') -> dict:
    from rdkit import Chem
    if shared_coordinate_frame is not True:
        raise ValueError('explicit shared_coordinate_frame=True is required; assembly does not align coordinates')
    if not isinstance(ligand_chain, str) or not re.fullmatch('[A-Za-z0-9]', ligand_chain):
        raise ValueError('ligand_chain must be a single alphanumeric PDB chain ID')
    validate_artifact(receptor_file, 'prepared_receptor')
    validate_artifact(poses_file, 'pose_set')
    receptor = read_structure(receptor_file)
    if receptor.source_format != 'pdb' or len(receptor.models) != 1:
        raise ValueError('assembly currently requires one prepared PDB receptor model')
    model = receptor.first_model
    if ligand_chain in {chain.chain_id for chain in model.chains}:
        raise ValueError('ligand_chain collides with an existing receptor chain')
    if any(residue.altlocs for residue in model.residues):
        raise ValueError('resolve receptor alternate locations explicitly before assembly')
    # Preserve residue-referencing metadata and coordinate columns; rewrite only
    # atom serial references. Unsupported serial-bearing records are refused.
    raw = Path(receptor_file).read_text().splitlines()
    if any(line[:6].strip() in {'ANISOU', 'SIGATM', 'SIGUIJ'} for line in raw):
        raise ValueError('anisotropic/uncertainty records require a dedicated preserving converter')
    serials, mapping = {}, []
    for line in raw:
        if line.startswith(('ATOM  ', 'HETATM', 'TER   ')):
            old = int(line[6:11])
            if old in serials:
                raise ValueError('duplicate receptor atom/TER serial')
            serials[old] = len(serials) + 1
            mapping.append({'source_serial': old, 'output_serial': serials[old], 'record': line[:6].strip()})
    selected_path = Path(output_dir) / 'selected_ligand.sdf'
    selected = select_ligand(poses_file, selected_path, record_id=record_id)
    mol = load_ligands(selected_path)[0]
    if any(len(f'{value:8.3f}') != 8 for row in mol.GetConformer().GetPositions() for value in row):
        raise ValueError('ligand coordinates exceed the PDB coordinate field width')
    if len(serials) + mol.GetNumAtoms() > 99999:
        raise ValueError('PDB atom serial limit exceeded; explicit mmCIF assembly is required')
    for i, atom in enumerate(mol.GetAtoms(), 1):
        if i > 999 or len(atom.GetSymbol()) > 2 or abs(atom.GetFormalCharge()) > 9:
            raise ValueError('ligand atom identity/charge cannot be represented by this PDB subset')
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f'L{i:03d}')
        info.SetResidueName('LIG')
        info.SetResidueNumber(1)
        info.SetChainId(ligand_chain)
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    ligand_block = Chem.MolToPDBBlock(mol).splitlines()
    lines, connections = [], []
    for line in raw:
        kind = line[:6].strip()
        if kind in {'END', 'ENDMDL', 'MODEL', 'MASTER', 'NUMMDL'}:
            continue
        if kind in {'ATOM', 'HETATM', 'TER'}:
            lines.append(line[:6] + f'{serials[int(line[6:11])]:5d}' + line[11:])
        elif kind == 'CONECT':
            entries = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
            if any(value not in serials for value in entries):
                raise ValueError('receptor CONECT refers to an absent serial')
            connections.append('CONECT' + ''.join(f'{serials[value]:5d}' for value in entries))
        else:
            lines.append(line)
    offset = len(serials)
    ligand_map = []
    for line in ligand_block:
        if line.startswith(('ATOM  ', 'HETATM')):
            old = int(line[6:11])
            lines.append(line[:6] + f'{old + offset:5d}' + line[11:])
            ligand_map.append({'source_atom_index': old - 1, 'output_serial': old + offset,
                               'atom_name': line[12:16].strip(), 'chain': ligand_chain, 'residue': 1})
        elif line.startswith('CONECT'):
            entries = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
            connections.append('CONECT' + ''.join(f'{value + offset:5d}' for value in entries))
    output = Path(output_dir) / 'complex.pdb'
    output.write_text('\n'.join(lines + connections + ['END']) + '\n')
    validate_artifact(output, 'complex_structure')
    actual = read_structure(output).first_model
    if len(actual.atoms) != len(model.atoms) + mol.GetNumAtoms():
        raise ValueError('assembled atom count mismatch')
    return {'schema_version': 1, 'receptor_sha256': file_hash(receptor_file), 'poses_sha256': file_hash(poses_file),
            'selection': selected, 'receptor_mapping': mapping, 'ligand_mapping': ligand_map,
            'complex_sha256': file_hash(output), 'ligand_sdf_sha256': file_hash(selected_path),
            'coordinate_frame': 'caller asserted shared frame; no alignment or minimization',
            'coordinate_precision': 'PDB rounds ligand coordinates to 0.001 angstrom; sidecar uses RDKit SDF serialization precision',
            'chemistry': 'selected_ligand.sdf is authoritative for ligand bond orders, charge, isotope and stereo; PDB alone is not lossless',
            'scientific_validation': 'not_performed'}
