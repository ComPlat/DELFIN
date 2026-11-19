"""Bounded OpenMM minimization/dynamics with immutable system identity and restart."""
from pathlib import Path
import math
import shutil
import subprocess
import sys
import time

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec, DataKeySpec
from delfin.tools._types import StepStatus, ErrorKind
from delfin.tools.bio.artifacts import child_path, read_json, write_json
from delfin.tools.bio.optional_runtime import probe_distribution, run_worker
from delfin.tools.bio.dynamics import (UNITS, copy_system, system_identity, file_record,
                                      validate_frame, validate_trajectory, validate_checkpoint, checked_file)
from delfin.tools.bio.simulation import validate_system

_WORKER = '''import json,sys,math
from pathlib import Path
import openmm as mm
from openmm import unit
options=json.load(open(sys.argv[1])); protocol=options['protocol']
bundle=json.load(open('system/simulation_system.json'))
root=Path('system')
system=mm.XmlSerializer.deserialize((root/bundle['files']['system']['path']).read_text())
integrator=mm.LangevinMiddleIntegrator(protocol['temperature_k']*unit.kelvin,
    protocol['friction_per_ps']/unit.picosecond,protocol['timestep_ps']*unit.picosecond)
integrator.setRandomNumberSeed(protocol['seed'])
integrator.setConstraintTolerance(protocol['constraint_tolerance'])
platform=mm.Platform.getPlatformByName(protocol['platform'])
properties={'Threads':str(protocol['cores'])} if protocol['platform']=='CPU' else {}
context=mm.Context(system,integrator,platform,properties)
minimization={'performed':False}
if options['restart']:
    context.loadCheckpoint(Path('restart.chk').read_bytes())
    if context.getStepCount()!=options['start_step'] or not math.isclose(context.getTime().value_in_unit(unit.picosecond),options['start_time_ps'],abs_tol=1e-7):
        raise ValueError('loaded checkpoint step/time differs from its manifest')
else:
    positions=json.load(open(root/bundle['files']['positions']['path']))['values']
    context.setPositions(positions*unit.nanometer)
    context.applyConstraints(protocol['constraint_tolerance'])
    if options['minimize']:
        before=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        mm.LocalEnergyMinimizer.minimize(context,options['min_tolerance_kj_mol_nm'],options['min_max_iterations'])
        after=context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        minimization={'performed':True,'before_kj_mol':before,'after_kj_mol':after,
                      'convergence':'not independently established; bounded OpenMM minimizer returned'}
    context.setVelocitiesToTemperature(protocol['temperature_k']*unit.kelvin,protocol['seed'])
    context.applyVelocityConstraints(protocol['constraint_tolerance'])
    context.setTime(0*unit.picosecond); context.setStepCount(0)
def frame():
    state=context.getState(getPositions=True,getEnergy=True)
    return dict(step=context.getStepCount(),time_ps=state.getTime().value_in_unit(unit.picosecond),
        potential_kj_mol=state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole),
        kinetic_kj_mol=state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole),
        positions_nm=[list(p) for p in state.getPositions().value_in_unit(unit.nanometer)])
with open('frames.jsonl','w') as stream:
    stream.write(json.dumps(frame(),allow_nan=False)+'\\n'); stream.flush()
    remaining=options['steps']
    while remaining:
        count=min(options['report_interval'],remaining)
        integrator.step(count); remaining-=count
        stream.write(json.dumps(frame(),allow_nan=False)+'\\n'); stream.flush()
state=context.getState(getPositions=True,getVelocities=True,getEnergy=True,getParameters=True,getIntegratorParameters=True)
Path('state.xml').write_text(mm.XmlSerializer.serialize(state))
Path('integrator.xml').write_text(mm.XmlSerializer.serialize(integrator))
Path('checkpoint.chk').write_bytes(context.createCheckpoint())
Path('worker_report.json').write_text(json.dumps(dict(minimization=minimization,final=frame(),
    platform=platform.getName(),platform_properties={name:platform.getPropertyValue(context,name) for name in platform.getPropertyNames()}),allow_nan=False))
'''


class OpenMMRun(StepAdapter):
    name = 'bio_openmm_run'
    category = 'structure'
    description = 'Minimize and/or run explicit Langevin dynamics from a checked OpenMM system; export a trajectory and exact-protocol checkpoint'
    produces_geometry = False
    consumes = ('simulation_system',)
    wires = (('simulation_system', 'system_file'),)
    produces = ('trajectory', 'checkpoint', 'json_report')
    required_artifacts = produces
    params = (ParamSpec('system_file', 'path', required=True), ParamSpec('expected_version', 'str', required=True),
              ParamSpec('platform', 'str', required=True, enum=('CPU', 'Reference')),
              ParamSpec('temperature_k', 'float', required=True), ParamSpec('friction_per_ps', 'float', required=True),
              ParamSpec('timestep_ps', 'float', required=True), ParamSpec('seed', 'int', required=True),
              ParamSpec('constraint_tolerance', 'float', required=True), ParamSpec('steps', 'int', required=True),
              ParamSpec('report_interval', 'int', required=True), ParamSpec('minimize', 'bool', required=True),
              ParamSpec('min_tolerance_kj_mol_nm', 'float', default=10.0), ParamSpec('min_max_iterations', 'int', default=1000),
              ParamSpec('checkpoint_file', 'path'), ParamSpec('max_frames', 'int', default=10000),
              ParamSpec('timeout', 'float', default=3600.0), ParamSpec('python', 'path', default=sys.executable))
    data_keys = (DataKeySpec('runtime', 'dict'), DataKeySpec('validation', 'dict'), DataKeySpec('minimization', 'dict'))
    example_kwargs = dict(system_file='simulation_system.json', expected_version='8.4.0', platform='CPU', temperature_k=300.,
                          friction_per_ps=1., timestep_ps=0.001, seed=7, constraint_tolerance=1e-5, steps=100,
                          report_interval=10, minimize=True)

    def validate_params(self, **kwargs):
        for p in self.params:
            if p.required and (kwargs.get(p.name) is None or kwargs.get(p.name) == ''):
                raise ValueError(f'{p.name} is required')
        if kwargs['platform'] not in ('CPU', 'Reference') or type(kwargs['minimize']) is not bool:
            raise ValueError('explicit platform and minimization policy required')
        for key, default in [('temperature_k', None), ('friction_per_ps', None), ('timestep_ps', None),
                             ('constraint_tolerance', None), ('min_tolerance_kj_mol_nm', 10.), ('timeout', 3600.)]:
            value = kwargs.get(key, default)
            if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
                raise ValueError(f'{key} must be positive and finite')
        if kwargs['constraint_tolerance'] >= 1:
            raise ValueError('constraint_tolerance must be less than one')
        for key, default, low in [('steps', None, 0), ('report_interval', None, 1), ('seed', None, 1),
                                  ('min_max_iterations', 1000, 1), ('max_frames', 10000, 1)]:
            value = kwargs.get(key, default)
            if type(value) is not int or not low <= value < 2**31:
                raise ValueError(f'invalid {key}')
        if kwargs['steps'] == 0 and not kwargs['minimize']:
            raise ValueError('zero steps requires explicit minimization')
        if kwargs.get('checkpoint_file') and kwargs['minimize']:
            raise ValueError('checkpoint continuation cannot minimize or reinitialize velocities')
        count = 1 + (kwargs['steps'] + kwargs['report_interval'] - 1) // kwargs['report_interval']
        if count > kwargs.get('max_frames', 10000):
            raise ValueError('requested trajectory exceeds max_frames')

    def runtime_preflight(self, **kwargs):
        return probe_distribution('openmm', python=kwargs.get('python', sys.executable), expected_version=kwargs['expected_version'],
                                  work_dir=Path.cwd(), allow_install=kwargs.get('_allow_install', False))

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start, root = time.monotonic(), Path(work_dir)
        try:
            self.validate_params(**kwargs)
            if type(cores) is not int or cores < 1 or (kwargs['platform'] == 'Reference' and cores != 1):
                raise ValueError('positive cores required; Reference platform requires cores=1')
            details = validate_system(kwargs['system_file'])
            identity = system_identity(kwargs['system_file'])
            protocol = {key: kwargs[key] for key in ('expected_version', 'platform', 'temperature_k', 'friction_per_ps', 'timestep_ps', 'seed', 'constraint_tolerance')}
            protocol.update(integrator='LangevinMiddleIntegrator', cores=cores)
            checkpoint = None
            if kwargs.get('checkpoint_file'):
                checkpoint = validate_checkpoint(kwargs['checkpoint_file'])
                if checkpoint['protocol'] != protocol or checkpoint['system_identity'] != identity:
                    raise ValueError('restart refused: system, runtime or integrator protocol changed')
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('simulation requires an empty step directory')
            system = copy_system(kwargs['system_file'], root / 'system')
            if checkpoint:
                shutil.copyfile(checked_file(Path(kwargs['checkpoint_file']).resolve().parent, checkpoint['binary']), root / 'restart.chk')
            start_step, start_time = (checkpoint['step'], checkpoint['time_ps']) if checkpoint else (0, 0.)
            options = dict(protocol=protocol, steps=kwargs['steps'], report_interval=kwargs['report_interval'],
                           minimize=kwargs['minimize'], min_tolerance_kj_mol_nm=kwargs.get('min_tolerance_kj_mol_nm', 10.),
                           min_max_iterations=kwargs.get('min_max_iterations', 1000), restart=bool(checkpoint), start_step=start_step, start_time_ps=start_time)
            runtime = run_worker(_WORKER, options, {'openmm': kwargs['expected_version']}, python=kwargs.get('python', sys.executable),
                                 work_dir=root, timeout=kwargs.get('timeout', 3600.))
            data = read_json(root / 'worker_report.json')
            validate_frame(data['final'], details['particles'])
            if data.get('platform') != kwargs['platform']:
                raise ValueError('engine platform differs from requested platform')
            minimized = data.get('minimization', {})
            if minimized.get('performed') is not kwargs['minimize']:
                raise ValueError('engine minimization policy mismatch')
            if kwargs['minimize']:
                before, after = minimized.get('before_kj_mol'), minimized.get('after_kj_mol')
                if any(type(x) not in (int, float) or not math.isfinite(x) for x in (before, after)) or after > before + max(1e-5, abs(before)*1e-8):
                    raise ValueError('minimization energy is invalid or increased')
            count = 1 + (kwargs['steps'] + kwargs['report_interval'] - 1) // kwargs['report_interval']
            trajectory = dict(schema_version=1, format='delfin_openmm_jsonl', units=UNITS, system=file_record(root, system),
                              system_identity=identity, frames=file_record(root, root / 'frames.jsonl'), frame_count=count,
                              start_step=start_step, start_time_ps=start_time, steps=kwargs['steps'], report_interval=kwargs['report_interval'],
                              timestep_ps=kwargs['timestep_ps'], protocol=protocol, runtime=runtime)
            write_json(root / 'trajectory.json', trajectory)
            validation = validate_trajectory(root / 'trajectory.json')
            if data['final']['step'] != validation['final_step'] or not math.isclose(data['final']['time_ps'], validation['final_time_ps'], abs_tol=1e-7):
                raise ValueError('engine final state differs from trajectory')
            restart = dict(schema_version=1, engine='openmm', system=file_record(root, system), system_identity=identity,
                           protocol=protocol, step=validation['final_step'], time_ps=validation['final_time_ps'],
                           binary=file_record(root, root / 'checkpoint.chk'), state=file_record(root, root / 'state.xml'),
                           integrator=file_record(root, root / 'integrator.xml'), trajectory=file_record(root, root / 'trajectory.json'), runtime=runtime,
                           compatibility='Exact protocol required; OpenMM may additionally reject different hardware or builds. No state-file fallback.')
            write_json(root / 'checkpoint.json', restart)
            validate_checkpoint(root / 'checkpoint.json')
            data.update(runtime=runtime, validation=validation, protocol=protocol, scientific_validation='not_performed')
            write_json(root / 'report.json', data)
            return self._make_result(self.name, StepStatus.SUCCESS, root, start, data=data,
                                     artifacts={key: root / filename for key, filename in dict(trajectory='trajectory.json', checkpoint='checkpoint.json', json_report='report.json').items()})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, TypeError, KeyError, OSError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, root, start, error=error, error_kind=kind)


register(OpenMMRun())
