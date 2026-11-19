"""Checked OpenMM system bundles; validation does not import a simulation engine."""
from pathlib import Path
import math
import xml.etree.ElementTree as ET

from delfin.tools.bio.artifacts import child_path, file_hash, read_json


def validate_system(path):
    path = Path(path).resolve()
    bundle = read_json(path)
    if not isinstance(bundle, dict) or bundle.get('schema_version') != 1 or bundle.get('engine') != 'openmm':
        raise ValueError('invalid simulation system schema or engine')
    files = bundle.get('files')
    required = {'structure', 'system', 'topology', 'positions', 'forcefield'}
    if not isinstance(files, dict) or set(files) != required:
        raise ValueError('simulation system requires structure, system, topology, positions and forcefield')
    paths = {}
    for name, record in files.items():
        if not isinstance(record, dict) or set(record) != {'path', 'sha256'}:
            raise ValueError('invalid simulation artifact record')
        target = child_path(path.parent, record['path'])
        if not target.is_file() or not target.stat().st_size or file_hash(target) != record['sha256']:
            raise ValueError(f'simulation {name} missing or checksum mismatch')
        paths[name] = target
    if len(set(paths.values())) != len(paths):
        raise ValueError('simulation artifacts must be distinct files')
    topology, positions = read_json(paths['topology']), read_json(paths['positions'])
    if not isinstance(topology, dict) or not isinstance(positions, dict):
        raise ValueError('simulation topology and positions must be JSON objects')
    atoms, bonds = topology.get('atoms'), topology.get('bonds')
    count = bundle.get('particle_count')
    if type(count) is not int or count < 1 or not isinstance(atoms, list) or len(atoms) != count:
        raise ValueError('simulation particle/atom count mismatch')
    if any(not isinstance(a, dict) or type(a.get('index')) is not int or a['index'] != i for i, a in enumerate(atoms)):
        raise ValueError('simulation atoms require consecutive indices')
    if not isinstance(bonds, list):
        raise ValueError('simulation topology requires bonds')
    seen = set()
    for bond in bonds:
        if not isinstance(bond, list) or len(bond) != 2 or any(type(i) is not int or not 0 <= i < count for i in bond) or bond[0] == bond[1]:
            raise ValueError('invalid simulation bond indices')
        key = tuple(sorted(bond))
        if key in seen:
            raise ValueError('duplicate simulation bond')
        seen.add(key)
    coordinates = positions.get('values')
    if positions.get('unit') != 'nanometer' or not isinstance(coordinates, list) or len(coordinates) != count:
        raise ValueError('simulation positions require nanometer units and one coordinate per atom')
    if any(not isinstance(row, list) or len(row) != 3 or any(type(v) not in (int, float) or not math.isfinite(v) for v in row) for row in coordinates):
        raise ValueError('invalid simulation coordinates')
    from delfin.tools.bio import read_structure
    structure = read_structure(paths['structure'])
    source = [(r, a) for r in structure.first_model.residues for a in r.atoms]
    if len(structure.models) != 1 or len(source) != count:
        raise ValueError('simulation structure atom count mismatch')
    for (residue, atom), record, xyz in zip(source, atoms, coordinates):
        identity = dict(name=atom.name, element=atom.element, chain=residue.chain_id,
                        residue_id=str(residue.seq_id), residue_name=residue.name, insertion_code=residue.insertion_code)
        if any(record.get(k) != v for k, v in identity.items()):
            raise ValueError('simulation topology changed atom identity or ordering')
        if any(abs(a / 10 - b) > 0.000011 for a, b in zip(atom.position, xyz)):
            raise ValueError('simulation coordinates differ from input structure')
    try:
        system = ET.parse(paths['system']).getroot()
    except ET.ParseError as exc:
        raise ValueError('invalid OpenMM system XML') from exc
    particles = system.findall('./Particles/Particle')
    if system.tag != 'System' or len(particles) != count:
        raise ValueError('serialized system particle count mismatch')
    if any(len(particle) for particle in particles):
        raise ValueError('virtual sites are not supported by this system schema')
    if any(not math.isfinite(float(p.get('mass', 'nan'))) or float(p.get('mass', 'nan')) <= 0 for p in particles):
        raise ValueError('simulation particle masses must be positive and finite')
    forces = system.findall('./Forces/Force')
    nonbonded = [force for force in forces if force.get('type') == 'NonbondedForce']
    if len(nonbonded) != 1 or len(nonbonded[0].findall('./Particles/Particle')) != count:
        raise ValueError('system requires one complete NonbondedForce')
    charges = [float(p.get('q', 'nan')) for p in nonbonded[0].findall('./Particles/Particle')]
    expected = bundle.get('expected_charge')
    if type(expected) is not int or not all(math.isfinite(q) for q in charges) or abs(sum(charges) - expected) > 1e-4:
        raise ValueError('simulation charge differs from explicit expected charge')
    if bundle.get('environment') != 'vacuum' or bundle.get('nonbonded_method') != 'NoCutoff' or nonbonded[0].get('method') != '0':
        raise ValueError('only explicit vacuum NoCutoff systems are currently supported')
    if bundle.get('forcefield_name') != 'amber14/protein.ff14SB.xml':
        raise ValueError('unsupported force field')
    if bundle.get('constraints') not in ('none', 'HBonds'):
        raise ValueError('explicit constraints policy required')
    actual_constraints = set()
    for constraint in system.findall('./Constraints/Constraint'):
        pair = tuple(sorted((int(constraint.get('p1', '-1')), int(constraint.get('p2', '-1')))))
        distance = float(constraint.get('d', 'nan'))
        if pair in actual_constraints or pair not in seen or not math.isfinite(distance) or distance <= 0:
            raise ValueError('invalid or duplicated system constraint')
        actual_constraints.add(pair)
    expected_constraints = {pair for pair in seen if any(atoms[i]['element'] == 'H' for i in pair)} if bundle['constraints'] == 'HBonds' else set()
    if actual_constraints != expected_constraints:
        raise ValueError('system constraints differ from explicit policy and topology')
    try:
        forcefield = ET.parse(paths['forcefield']).getroot()
    except ET.ParseError as exc:
        raise ValueError('invalid force field XML') from exc
    if forcefield.tag != 'ForceField' or forcefield.find('Include') is not None:
        raise ValueError('force field must be a self-contained XML artifact')
    return {'particles': count, 'bonds': len(bonds), 'charge_e': sum(charges), 'engine': 'openmm'}
