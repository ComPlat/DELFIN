"""Checked OpenMM system bundles; validation does not import a simulation engine."""
from pathlib import Path
import math
import xml.etree.ElementTree as ET

from delfin.tools.bio.artifacts import child_path, file_hash, read_json


def validate_system(path):
    path = Path(path).resolve()
    bundle = read_json(path)
    if not isinstance(bundle, dict) or bundle.get('schema_version') not in (1, 2, 3, 4) or bundle.get('engine') != 'openmm':
        raise ValueError('invalid simulation system schema or engine')
    files = bundle.get('files')
    required = {'structure', 'system', 'topology', 'positions', 'forcefield'}
    periodic = bundle['schema_version'] in (2, 4)
    ligand = bundle['schema_version'] in (3, 4)
    if ligand:
        required.add('ligand_parameters')
    if periodic:
        required |= {'source', 'solvent_forcefield'}
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
    if ligand:
        from delfin.tools.bio.ligand_parameters import complex_ligand_mapping
        mapping, mol = complex_ligand_mapping(structure, paths['ligand_parameters'])
        parameters = read_json(paths['ligand_parameters'])
        values = read_json(child_path(paths['ligand_parameters'].parent, parameters['files']['charges']['path']))['values']
        if any(not math.isclose(charges[index], values[original], abs_tol=1e-7) for index, original in mapping.items()):
            raise ValueError('system ligand partial charges differ from supplied values')
        ligand_xml = ET.parse(child_path(paths['ligand_parameters'].parent, parameters['files']['parameters']['path'])).getroot()
        ligand_atoms = ligand_xml.findall('./Residues/Residue/Atom')
        ligand_types = {a.get('name'): a for a in ligand_xml.findall('./AtomTypes/Type')}
        ligand_nb = {a.get('class'): a for a in ligand_xml.findall('./NonbondedForce/Atom')}
        actual_nb = nonbonded[0].findall('./Particles/Particle')
        for index, original in mapping.items():
            atom_type = ligand_atoms[original].get('type')
            if not math.isclose(float(particles[index].get('mass')), float(ligand_types[atom_type].get('mass')), abs_tol=1e-7):
                raise ValueError('system ligand mass differs from generated parameters')
            for actual, wanted in (('sig', 'sigma'), ('eps', 'epsilon')):
                if not math.isclose(float(actual_nb[index].get(actual, 'nan')), float(ligand_nb[atom_type].get(wanted)), abs_tol=1e-7):
                    raise ValueError('system ligand Lennard-Jones parameters differ from parameter bundle')
        actual = set()
        for a, b in seen:
            if (a in mapping) != (b in mapping):
                raise ValueError('covalent ligand links require explicit boundary parameterization')
            if a in mapping:
                actual.add(tuple(sorted((mapping[a], mapping[b]))))
        wanted = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()}
        if actual != wanted:
            raise ValueError('system ligand bonds differ from authoritative SDF')
    if periodic:
        _validate_periodic(bundle, paths, system, nonbonded[0], source, coordinates, charges)
    elif bundle.get('environment') != 'vacuum' or bundle.get('nonbonded_method') != 'NoCutoff' or nonbonded[0].get('method') != '0':
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
    return {'particles': count, 'bonds': len(bonds), 'charge_e': sum(charges), 'engine': 'openmm', 'periodic': periodic, 'ligand': ligand}


def _validate_periodic(bundle, paths, system, nonbonded, prepared, coordinates, charges):
    from delfin.tools.bio import read_structure
    if bundle.get('environment') != 'tip3p_pme' or bundle.get('nonbonded_method') != 'PME' or nonbonded.get('method') != '4':
        raise ValueError('periodic system requires explicit TIP3P/PME policy')
    settings = bundle.get('solvation')
    validate_solvation(settings)
    if any(system.find(f'./PeriodicBoxVectors/{name}') is None for name in 'ABC'):
        raise ValueError('serialized periodic box is incomplete')
    vectors = [[float(system.find(f'./PeriodicBoxVectors/{name}').get(k, 'nan')) for k in 'xyz'] for name in 'ABC']
    expected = [[settings['box_size_nm'][i] if i == j else 0. for j in range(3)] for i in range(3)]
    if any(not math.isclose(a, b, abs_tol=0.00006) for row, target in zip(vectors, expected) for a, b in zip(row, target)):
        raise ValueError('serialized periodic box differs from requested orthorhombic box')
    if not math.isclose(float(nonbonded.get('cutoff', 'nan')), settings['cutoff_nm'], abs_tol=1e-10) or not math.isclose(float(nonbonded.get('ewaldTolerance', 'nan')), settings['ewald_tolerance'], abs_tol=1e-12):
        raise ValueError('PME cutoff/tolerance differs from requested policy')
    original = read_structure(paths['source'])
    original_atoms = [(r, a) for r in original.first_model.residues for a in r.atoms]
    count = len(original_atoms)
    if len(original.models) != 1 or len(prepared) <= count:
        raise ValueError('solvation must retain the source and add explicit solvent')
    for i, ((r, a), (new_r, new_a)) in enumerate(zip(original_atoms, prepared)):
        if (r.chain_id, r.seq_id, r.insertion_code, r.name, a.name, a.element) != (new_r.chain_id, new_r.seq_id, new_r.insertion_code, new_r.name, new_a.name, new_a.element):
            raise ValueError('solvation changed source atom identity or ordering')
        if any(abs(x / 10 - y) > 0.000011 for x, y in zip(a.position, coordinates[i])):
            raise ValueError('solvation moved source atoms')
    if abs(sum(charges[:count]) - settings['solute_charge']) > 1e-4:
        raise ValueError('solute parameter charge differs from requested charge')
    allowed_elements = {'O', 'H', 'Cl', 'Na' if settings['positive_ion'] == 'Na+' else 'K'}
    if any(r.name not in {'HOH', 'NA', 'CL', 'K'} or a.element not in allowed_elements for r, a in prepared[count:]):
        raise ValueError('unexpected added solvent/ion chemistry')
    if not any(r.name == 'HOH' for r, a in prepared[count:]):
        raise ValueError('solvated system contains no water')
    residues = {}
    for index, (r, a) in enumerate(prepared[count:], count):
        residues.setdefault((r.chain_id, r.seq_id, r.insertion_code, r.name), []).append((a, charges[index]))
    waters, positive, negative = 0, 0, 0
    for key, members in residues.items():
        if key[-1] == 'HOH':
            if sorted(a.element for a, _ in members) != ['H', 'H', 'O'] or abs(sum(q for _, q in members)) > 1e-5:
                raise ValueError('TIP3P water requires neutral H2O residues')
            waters += 1
        else:
            if len(members) != 1:
                raise ValueError('ion residue must contain exactly one atom')
            atom, charge = members[0]
            sign = -1 if atom.element == 'Cl' else 1
            if abs(charge - sign) > 1e-5:
                raise ValueError('ion charge differs from declared monovalent model')
            positive += sign == 1
            negative += sign == -1
    neutral = abs(settings['solute_charge']) if settings['neutralize'] else 0
    pairs = math.floor((waters + positive + negative - neutral) * settings['ionic_strength_m'] / 55.4 + 0.5)
    expected_positive = pairs + (max(-settings['solute_charge'], 0) if settings['neutralize'] else 0)
    expected_negative = pairs + (max(settings['solute_charge'], 0) if settings['neutralize'] else 0)
    if (positive, negative) != (expected_positive, expected_negative):
        raise ValueError('ion counts differ from explicit neutralization and concentration policy')
    if settings['neutralize'] and bundle['expected_charge'] != 0:
        raise ValueError('neutralized system requires zero final charge')
    if not settings['neutralize'] and bundle['expected_charge'] != settings['solute_charge']:
        raise ValueError('non-neutralized system must retain solute net charge')
    try:
        solvent = ET.parse(paths['solvent_forcefield']).getroot()
    except ET.ParseError as exc:
        raise ValueError('invalid solvent force field XML') from exc
    if solvent.tag != 'ForceField' or solvent.find('Include') is not None:
        raise ValueError('solvent force field must be self-contained')


def validate_solvation(settings):
    required = {'box_size_nm', 'ionic_strength_m', 'neutralize', 'positive_ion', 'negative_ion', 'seed', 'solute_charge', 'cutoff_nm', 'ewald_tolerance'}
    if not isinstance(settings, dict) or set(settings) != required:
        raise ValueError('solvation requires explicit box, ions, concentration, neutralization, seed, solute charge and PME settings')
    box = settings['box_size_nm']
    if not isinstance(box, list) or len(box) != 3 or any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0 for v in box):
        raise ValueError('box_size_nm requires three positive finite lengths')
    if settings['positive_ion'] not in ('Na+', 'K+') or settings['negative_ion'] != 'Cl-' or type(settings['neutralize']) is not bool:
        raise ValueError('explicit Na+/K+, Cl- and neutralization policy required')
    if type(settings['solute_charge']) is not int or type(settings['seed']) is not int or not 1 <= settings['seed'] < 2**31:
        raise ValueError('solute charge and nonzero solvation seed must be integers')
    for key in ('ionic_strength_m', 'cutoff_nm', 'ewald_tolerance'):
        if type(settings[key]) not in (int, float) or not math.isfinite(settings[key]) or settings[key] < 0:
            raise ValueError(f'invalid {key}')
    if not 0 < settings['cutoff_nm'] < min(box) / 2 or not 0 < settings['ewald_tolerance'] < 1:
        raise ValueError('PME requires positive cutoff below half the shortest box length and tolerance in (0,1)')
