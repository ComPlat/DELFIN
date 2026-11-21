"""Portable single-ligand force-field bundles and explicit charge contracts."""
from pathlib import Path
import math
import shutil
import xml.etree.ElementTree as ET

from delfin.tools.bio.artifacts import child_path, file_hash, read_json
from delfin.tools.bio.processing import load_ligands, chemical_identity


def _xml(path):
    try:
        return ET.parse(path).getroot()
    except ET.ParseError as exc:
        raise ValueError('invalid ligand force-field XML') from exc


def ligand_input(path):
    from rdkit import Chem
    path = Path(path)
    if path.suffix.lower() != '.sdf':
        raise ValueError('ligand parameterization requires an authoritative SDF')
    molecules = load_ligands(path)
    if len(molecules) != 1:
        raise ValueError('ligand parameterization requires exactly one molecule')
    mol = molecules[0]
    if len(Chem.GetMolFrags(mol)) != 1 or not 1 < mol.GetNumAtoms() <= 999:
        raise ValueError('ligand must be one connected molecule with 2..999 atoms')
    if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        raise ValueError('ligand requires one explicit 3D conformer')
    if any(a.GetNumImplicitHs() or a.GetNumExplicitHs() for a in mol.GetAtoms()):
        raise ValueError('all ligand hydrogens must be explicit atoms')
    if any(a.GetIsotope() or a.GetNumRadicalElectrons() or a.GetSymbol() not in
           {'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I'} for a in mol.GetAtoms()):
        raise ValueError('isotopes, radicals and non-organic elements require a separate parameterization policy')
    if any(s.specified != Chem.StereoSpecified.Specified for s in Chem.FindPotentialStereo(mol)):
        raise ValueError('undefined ligand stereochemistry must be resolved explicitly')
    return mol


def charge_input(path, ligand_path, mol):
    from rdkit import Chem
    value = read_json(path)
    required = {'schema_version', 'ligand_sha256', 'unit', 'method', 'source', 'values'}
    if not isinstance(value, dict) or set(value) != required or value['schema_version'] != 1 or value['unit'] != 'elementary_charge':
        raise ValueError('charges require schema 1, SDF hash, elementary_charge units, method, source and ordered values')
    if value['ligand_sha256'] != file_hash(ligand_path):
        raise ValueError('charges belong to a different SDF or atom ordering')
    if any(not isinstance(value[k], str) or not value[k].strip() for k in ('method', 'source')):
        raise ValueError('charge method and provenance source must be explicit')
    charges = value['values']
    if not isinstance(charges, list) or len(charges) != mol.GetNumAtoms() or any(type(q) not in (int, float) or not math.isfinite(q) for q in charges):
        raise ValueError('one finite partial charge is required per SDF atom')
    if abs(sum(charges) - Chem.GetFormalCharge(mol)) > 1e-4:
        raise ValueError('partial charges do not sum to the ligand formal charge')
    if all(abs(q) <= 1e-8 for q in charges):
        raise ValueError('all-zero charges are unsupported: the selected generator would infer another charge method')
    return value


def validate_offxml(path):
    root = _xml(path)
    allowed = {'Author', 'Date', 'Constraints', 'Bonds', 'Angles', 'ProperTorsions',
               'ImproperTorsions', 'vdW', 'Electrostatics', 'LibraryCharges', 'ToolkitAM1BCC'}
    if root.tag != 'SMIRNOFF' or any(child.tag not in allowed for child in root):
        raise ValueError('local SMIRNOFF file contains unsupported handlers or includes')
    if root.find('vdW') is None or root.find('Electrostatics') is None:
        raise ValueError('SMIRNOFF file requires explicit vdW and electrostatic handlers')


def validate_parameters_xml(path, mol, charges):
    root = _xml(path)
    allowed = {'AtomTypes', 'Residues', 'HarmonicBondForce', 'HarmonicAngleForce', 'PeriodicTorsionForce', 'NonbondedForce'}
    if root.tag != 'ForceField' or {child.tag for child in root} != allowed or len(root) != len(allowed):
        raise ValueError('ligand parameters must be self-contained without scripts, includes or virtual sites')
    residues = root.findall('./Residues/Residue')
    if len(residues) != 1 or any(child.tag not in {'Atom', 'Bond'} for child in residues[0]):
        raise ValueError('ligand parameters require one residue without external bonds')
    atoms = residues[0].findall('Atom')
    names = [f'L{i + 1:03d}' for i in range(mol.GetNumAtoms())]
    if [a.get('name') for a in atoms] != names:
        raise ValueError('ligand parameter atom ordering differs from SDF mapping')
    types = root.findall('./AtomTypes/Type')
    by_type = {t.get('name'): t for t in types}
    if len(by_type) != len(types):
        raise ValueError('duplicate ligand atom types')
    if len(types) != len(atoms) or any(t.get('class') != t.get('name') for t in types):
        raise ValueError('generated ligand parameters require one unique class/type per atom')
    for atom, expected, charge in zip(atoms, mol.GetAtoms(), charges):
        t = by_type.get(atom.get('type'))
        if t is None or t.get('element') != expected.GetSymbol():
            raise ValueError('ligand parameter element mismatch')
        mass = float(t.get('mass', 'nan'))
        if not math.isfinite(mass) or mass <= 0:
            raise ValueError('ligand atom type requires positive finite mass')
        if not math.isclose(float(atom.get('charge', 'nan')), charge, abs_tol=1e-7):
            raise ValueError('ligand parameters changed the supplied partial charges')
    bonds = [tuple(sorted((b.get('atomName1'), b.get('atomName2')))) for b in residues[0].findall('Bond')]
    wanted = {tuple(sorted((names[b.GetBeginAtomIdx()], names[b.GetEndAtomIdx()]))) for b in mol.GetBonds()}
    if len(set(bonds)) != len(bonds) or set(bonds) != wanted:
        raise ValueError('ligand parameter bonds differ from authoritative SDF')
    nb = root.findall('NonbondedForce')
    if len(nb) != 1 or nb[0].find("UseAttributeFromResidue[@name='charge']") is None:
        raise ValueError('ligand nonbonded parameters must use the checked residue charges')
    records = {a.get('class'): a for a in nb[0].findall('Atom')}
    if len(records) != len(nb[0].findall('Atom')) or set(records) != {a.get('type') for a in atoms}:
        raise ValueError('incomplete ligand nonbonded atom parameters')
    for record in records.values():
        if any(not math.isfinite(float(record.get(k, 'nan'))) or float(record.get(k, 'nan')) < 0 for k in ('sigma', 'epsilon')):
            raise ValueError('invalid ligand Lennard-Jones parameters')
    for setting in ('coulomb14scale', 'lj14scale'):
        value = float(nb[0].get(setting, 'nan'))
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('explicit finite 1-4 nonbonded scaling required')
    type_indices = {a.get('type'): i for i, a in enumerate(atoms)}
    from itertools import combinations
    expected_bonds = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()}
    expected_angles = {min((a.GetIdx(), center.GetIdx(), b.GetIdx()), (b.GetIdx(), center.GetIdx(), a.GetIdx()))
                       for center in mol.GetAtoms() for a, b in combinations(center.GetNeighbors(), 2)}
    for query, size, expected, length in (('./HarmonicBondForce/Bond', 2, expected_bonds, 'length'),
                                         ('./HarmonicAngleForce/Angle', 3, expected_angles, 'angle')):
        observed = []
        for record in root.findall(query):
            try:
                indices = tuple(type_indices[record.get(f'class{i + 1}')] for i in range(size))
            except KeyError as exc:
                raise ValueError('bonded parameter refers to an unknown atom class') from exc
            observed.append(min(indices, indices[::-1]))
            values = [float(record.get(k, 'nan')) for k in (length, 'k')]
            if not all(math.isfinite(v) for v in values) or values[0] <= 0 or values[1] < 0:
                raise ValueError('invalid bonded ligand parameters')
        if len(observed) != len(set(observed)) or set(observed) != expected:
            raise ValueError('ligand bond/angle parameters are incomplete or duplicated')
    for record in root.findall('./PeriodicTorsionForce/*'):
        if record.tag not in {'Proper', 'Improper'} or any(record.get(f'class{i}') not in type_indices for i in range(1, 5)):
            raise ValueError('invalid ligand torsion mapping')
        for key, raw in record.attrib.items():
            if not key.startswith('class') and not math.isfinite(float(raw)):
                raise ValueError('nonfinite ligand torsion parameter')
    return residues[0].get('name')


def validate_ligand_parameters(path):
    path = Path(path).resolve()
    bundle = read_json(path)
    if not isinstance(bundle, dict) or bundle.get('schema_version') != 1 or bundle.get('method') != 'smirnoff_supplied_charges':
        raise ValueError('invalid ligand parameter bundle schema or method')
    files = bundle.get('files')
    if not isinstance(files, dict) or set(files) != {'ligand', 'charges', 'offxml', 'parameters'}:
        raise ValueError('ligand bundle requires SDF, supplied charges, OFFXML and generated parameters')
    paths = {}
    for name, record in files.items():
        if not isinstance(record, dict) or set(record) != {'path', 'sha256'}:
            raise ValueError('invalid ligand parameter file record')
        p = child_path(path.parent, record['path'])
        if not p.is_file() or not p.stat().st_size or file_hash(p) != record['sha256']:
            raise ValueError('missing or changed ligand parameter artifact: ' + name)
        paths[name] = p
    if len(set(paths.values())) != len(paths):
        raise ValueError('ligand parameter artifacts must be distinct')
    mol = ligand_input(paths['ligand'])
    charges = charge_input(paths['charges'], paths['ligand'], mol)
    validate_offxml(paths['offxml'])
    template = validate_parameters_xml(paths['parameters'], mol, charges['values'])
    identity = chemical_identity(mol)
    if bundle.get('chemical_state_id') != identity['chemical_state_id'] or bundle.get('atom_count') != mol.GetNumAtoms():
        raise ValueError('ligand parameter manifest identity mismatch')
    return {'atoms': mol.GetNumAtoms(), 'chemical_state_id': identity['chemical_state_id'],
            'charge_e': sum(charges['values']), 'template': template}


def copy_ligand_parameters(source, destination):
    source, destination = Path(source).resolve(), Path(destination)
    validate_ligand_parameters(source)
    bundle = read_json(source)
    destination.mkdir(parents=True, exist_ok=False)
    for entry in bundle['files'].values():
        original = child_path(source.parent, entry['path'])
        target = child_path(destination, entry['path'])
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original, target)
    result = destination / 'ligand_parameters.json'
    shutil.copyfile(source, result)
    validate_ligand_parameters(result)
    return result


def complex_ligand_mapping(structure, bundle_path):
    """Check the explicit DELFIN L001... mapping; never infer chemistry from PDB."""
    path = Path(bundle_path).resolve()
    validate_ligand_parameters(path)
    bundle = read_json(path)
    mol = ligand_input(child_path(path.parent, bundle['files']['ligand']['path']))
    residues = [r for r in structure.first_model.residues if r.name == 'LIG']
    if len(residues) != 1 or residues[0].altlocs or residues[0].insertion_code:
        raise ValueError('complex requires exactly one explicitly mapped LIG residue')
    ligand = residues[0]
    if len(ligand.atoms) != mol.GetNumAtoms():
        raise ValueError('complex ligand atom count differs from parameter bundle')
    names = {f'L{i + 1:03d}': i for i in range(mol.GetNumAtoms())}
    actual = [a.name for a in ligand.atoms]
    if len(set(actual)) != len(actual) or set(actual) != set(names):
        raise ValueError('complex ligand requires unique L001... atom names')
    for atom in ligand.atoms:
        i = names[atom.name]
        if atom.element != mol.GetAtomWithIdx(i).GetSymbol() or any(abs(a - b) > 0.0006 for a, b in zip(atom.position, mol.GetConformer().GetAtomPosition(i))):
            raise ValueError('complex ligand identity/coordinates differ from authoritative SDF')
    identities = {id(a): names[a.name] for a in ligand.atoms}
    mapping = {index: identities[id(atom)] for index, atom in enumerate(structure.first_model.atoms) if id(atom) in identities}
    if len(mapping) != mol.GetNumAtoms():
        raise ValueError('complex ligand mapping is incomplete')
    return mapping, mol
