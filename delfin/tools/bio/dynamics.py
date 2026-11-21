"""Portable simulation outputs with exact units, system identity and step coverage."""
from pathlib import Path
import json
import math
import shutil
import xml.etree.ElementTree as ET

from delfin.tools.bio.artifacts import child_path, file_hash, read_json, write_json
from delfin.tools.bio.simulation import validate_system

UNITS = {'positions': 'nanometer', 'time': 'picosecond', 'energy': 'kilojoule_per_mole'}


def copy_system(source, destination):
    source, destination = Path(source).resolve(), Path(destination)
    validate_system(source)
    bundle = read_json(source)
    destination.mkdir(parents=True, exist_ok=False)
    for name, entry in bundle['files'].items():
        if name == 'ligand_parameters':
            from delfin.tools.bio.ligand_parameters import copy_ligand_parameters
            target = copy_ligand_parameters(child_path(source.parent, entry['path']), destination / 'ligand_parameters')
            entry['path'] = str(target.relative_to(destination))
            continue
        target = destination / (name + child_path(source.parent, entry['path']).suffix)
        shutil.copyfile(child_path(source.parent, entry['path']), target)
        entry['path'] = target.name
    output = destination / 'simulation_system.json'
    write_json(output, bundle)
    validate_system(output)
    return output


def system_identity(path):
    validate_system(path)
    bundle = read_json(path)
    return {name: row['sha256'] for name, row in bundle['files'].items()}


def checked_file(root, record):
    if not isinstance(record, dict) or set(record) != {'path', 'sha256'}:
        raise ValueError('invalid simulation output file record')
    path = child_path(root, record['path'])
    if not path.is_file() or not path.stat().st_size or file_hash(path) != record['sha256']:
        raise ValueError('missing or changed simulation output: ' + record['path'])
    return path


def file_record(root, path):
    return {'path': str(Path(path).relative_to(root)), 'sha256': file_hash(path)}


def number(value):
    return type(value) in (int, float) and math.isfinite(value)


def system_box(path):
    path = Path(path).resolve()
    bundle = read_json(path)
    root = ET.parse(child_path(path.parent, bundle['files']['system']['path'])).getroot()
    vectors = []
    for name in 'ABC':
        element = root.find(f'./PeriodicBoxVectors/{name}')
        if element is None:
            raise ValueError('periodic system has no complete box vectors')
        vectors.append([float(element.get(k, 'nan')) for k in 'xyz'])
    return vectors


def same_box(actual, expected):
    return (isinstance(actual, list) and len(actual) == 3 and all(
        isinstance(row, list) and len(row) == 3 and all(number(a) and math.isclose(a, b, abs_tol=1e-7) for a, b in zip(row, target))
        for row, target in zip(actual, expected)))


def validate_frame(frame, particles):
    if not isinstance(frame, dict) or type(frame.get('step')) is not int or frame['step'] < 0:
        raise ValueError('invalid trajectory step')
    for key in ('time_ps', 'potential_kj_mol', 'kinetic_kj_mol'):
        if not number(frame.get(key)):
            raise ValueError('missing or nonfinite trajectory time/energy')
    if frame['time_ps'] < 0 or frame['kinetic_kj_mol'] < 0:
        raise ValueError('negative trajectory time or kinetic energy')
    coords = frame.get('positions_nm')
    if not isinstance(coords, list) or len(coords) != particles or any(
            not isinstance(row, list) or len(row) != 3 or not all(number(x) for x in row) for row in coords):
        raise ValueError('trajectory coordinates must contain one finite 3-vector per particle')


def iter_frames(path):
    with Path(path).open() as stream:
        for line in stream:
            if not line.strip():
                raise ValueError('empty trajectory frame')
            yield json.loads(line)


def frame_steps(start, steps, interval):
    yield start
    offset = interval
    while offset < steps:
        yield start + offset
        offset += interval
    if steps:
        yield start + steps


def validate_trajectory(path):
    path = Path(path).resolve()
    value = read_json(path)
    if not isinstance(value, dict) or value.get('schema_version') != 1 or value.get('format') != 'delfin_openmm_jsonl' or value.get('units') != UNITS:
        raise ValueError('unsupported trajectory schema, format or units')
    system = checked_file(path.parent, value.get('system'))
    details = validate_system(system)
    periodic = details.get('periodic', False)
    if value.get('periodic', False) is not periodic:
        raise ValueError('trajectory periodicity differs from system')
    box = system_box(system) if periodic else None
    if value.get('system_identity') != system_identity(system):
        raise ValueError('trajectory system identity mismatch')
    frames = checked_file(path.parent, value.get('frames'))
    for key in ('start_step', 'steps', 'report_interval'):
        if type(value.get(key)) is not int or value[key] < (1 if key == 'report_interval' else 0):
            raise ValueError('invalid trajectory step range')
    timestep = value.get('timestep_ps')
    start_time = value.get('start_time_ps')
    if not number(timestep) or timestep <= 0 or not number(start_time) or start_time < 0:
        raise ValueError('invalid trajectory time specification')
    expected = iter(frame_steps(value['start_step'], value['steps'], value['report_interval']))
    count, final = 0, None
    for frame in iter_frames(frames):
        validate_frame(frame, details['particles'])
        if periodic and not same_box(frame.get('box_vectors_nm'), box):
            raise ValueError('trajectory box differs from fixed-volume system')
        wanted = next(expected, None)
        if frame['step'] != wanted:
            raise ValueError('trajectory has missing, duplicate or out-of-order steps')
        wanted_time = start_time + (wanted - value['start_step']) * timestep
        if not math.isclose(frame['time_ps'], wanted_time, abs_tol=1e-7, rel_tol=1e-9):
            raise ValueError('trajectory time does not match timestep and step count')
        count += 1
        final = frame
    if next(expected, None) is not None or count != value.get('frame_count'):
        raise ValueError('trajectory is incomplete')
    return {'frames': count, 'particles': details['particles'], 'final_step': final['step'], 'final_time_ps': final['time_ps'], 'periodic': periodic}


def validate_checkpoint(path):
    path = Path(path).resolve()
    value = read_json(path)
    if not isinstance(value, dict) or value.get('schema_version') != 1 or value.get('engine') != 'openmm':
        raise ValueError('invalid checkpoint schema')
    system = checked_file(path.parent, value.get('system'))
    if value.get('system_identity') != system_identity(system):
        raise ValueError('checkpoint system identity mismatch')
    outputs = {key: checked_file(path.parent, value.get(key)) for key in ('binary', 'state', 'integrator')}
    if type(value.get('step')) is not int or value['step'] < 0 or not number(value.get('time_ps')) or value['time_ps'] < 0:
        raise ValueError('invalid checkpoint step/time')
    protocol = value.get('protocol')
    if not isinstance(protocol, dict) or not protocol.get('expected_version') or protocol.get('platform') not in ('CPU', 'Reference'):
        raise ValueError('checkpoint requires explicit runtime protocol')
    try:
        state, integrator = (ET.parse(outputs[key]).getroot() for key in ('state', 'integrator'))
    except ET.ParseError as exc:
        raise ValueError('invalid checkpoint state or integrator XML') from exc
    if state.tag != 'State' or int(state.get('stepCount', '-1')) != value['step'] or not math.isclose(float(state.get('time', 'nan')), value['time_ps'], abs_tol=1e-7):
        raise ValueError('checkpoint state step/time mismatch')
    particles = validate_system(system)['particles']
    for query in ('./Positions/Position', './Velocities/Velocity'):
        vectors = state.findall(query)
        if len(vectors) != particles or any(not math.isfinite(float(v.get(k, 'nan'))) for v in vectors for k in ('x', 'y', 'z')):
            raise ValueError('checkpoint state requires finite positions and velocities for every particle')
    if integrator.tag != 'Integrator' or integrator.get('type') != 'LangevinMiddleIntegrator' or protocol.get('integrator') != 'LangevinMiddleIntegrator':
        raise ValueError('checkpoint integrator mismatch')
    for attribute, setting in [('stepSize', 'timestep_ps'), ('temperature', 'temperature_k'), ('friction', 'friction_per_ps'), ('constraintTolerance', 'constraint_tolerance')]:
        if not number(protocol.get(setting)) or protocol[setting] <= 0 or not math.isclose(float(integrator.get(attribute, 'nan')), protocol[setting], abs_tol=1e-12):
            raise ValueError('serialized integrator differs from checkpoint protocol')
    if type(protocol.get('seed')) is not int or protocol['seed'] <= 0 or int(integrator.get('randomSeed', '-1')) != protocol['seed']:
        raise ValueError('checkpoint seed mismatch')
    if type(protocol.get('cores')) is not int or protocol['cores'] < 1:
        raise ValueError('checkpoint requires positive CPU count')
    trajectory_path = checked_file(path.parent, value.get('trajectory'))
    trajectory = read_json(trajectory_path)
    verified = validate_trajectory(trajectory_path)
    if trajectory['system_identity'] != value['system_identity'] or trajectory.get('protocol') != protocol or verified['final_step'] != value['step']:
        raise ValueError('checkpoint does not match its trajectory')
    final = None
    for final in iter_frames(checked_file(trajectory_path.parent, trajectory['frames'])):
        pass
    state_positions = [[float(v.get(k)) for k in ('x', 'y', 'z')] for v in state.findall('./Positions/Position')]
    if any(not math.isclose(a, b, abs_tol=1e-8) for row, expected in zip(state_positions, final['positions_nm']) for a, b in zip(row, expected)):
        raise ValueError('checkpoint positions differ from final trajectory frame')
    if not math.isclose(value['time_ps'], verified['final_time_ps'], abs_tol=1e-7):
        raise ValueError('checkpoint time differs from trajectory')
    if trajectory.get('periodic'):
        vectors = state.find('PeriodicBoxVectors')
        if vectors is None or any(vectors.find(k) is None for k in 'ABC'):
            raise ValueError('periodic checkpoint has no box vectors')
        state_box = [[float(vectors.find(name).get(k, 'nan')) for k in 'xyz'] for name in 'ABC']
        if not same_box(state_box, final['box_vectors_nm']):
            raise ValueError('checkpoint box differs from trajectory')
    return value
