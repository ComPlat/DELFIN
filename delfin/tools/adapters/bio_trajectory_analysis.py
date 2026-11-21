"""Explicit-selection MDAnalysis fitting and streaming trajectory statistics."""
from pathlib import Path
import math
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec, DataKeySpec
from delfin.tools._types import StepStatus, ErrorKind
from delfin.tools.bio.artifacts import child_path, read_json, write_json, file_hash
from delfin.tools.bio.dynamics import validate_trajectory, system_identity, checked_file, iter_frames, number
from delfin.tools.bio.optional_runtime import probe_distribution, run_worker

_WORKER = '''import json,sys
from pathlib import Path
import numpy as np
from MDAnalysis.analysis.align import rotation_matrix
options=json.load(open(sys.argv[1])); indices=options['atom_indices']
masses=np.asarray(options['masses'],dtype=float)
weights=masses if options['weighting']=='mass' else np.ones(len(indices))
reference=None; mean=None; m2=None; rows=[]; count=0
with open('frames.jsonl') as stream:
    for index,line in enumerate(stream):
        if index<options['discard_frames'] or (index-options['discard_frames'])%options['stride']:
            continue
        frame=json.loads(line); coordinates=np.asarray(frame['positions_nm'],dtype=float)[indices]
        center=np.average(coordinates,axis=0,weights=weights)
        centered=coordinates-center
        rg=float(np.sqrt(np.average(np.sum(centered**2,axis=1),weights=weights)))
        if options['alignment']=='fit':
            if np.linalg.matrix_rank(centered)<2:
                raise ValueError('alignment requires non-collinear selected coordinates in every frame')
            if reference is None:
                reference=centered.copy()
            rotation,_=rotation_matrix(centered,reference,weights=weights)
            observed=centered@rotation.T
        else:
            observed=coordinates
            if reference is None:
                reference=observed.copy()
        rmsd=float(np.sqrt(np.average(np.sum((observed-reference)**2,axis=1),weights=weights)))
        count+=1
        if mean is None:
            mean=observed.copy(); m2=np.zeros_like(observed)
        else:
            delta=observed-mean; mean+=delta/count; m2+=delta*(observed-mean)
        rows.append(dict(step=frame['step'],time_ps=frame['time_ps'],rmsd_nm=rmsd,radius_gyration_nm=rg,
                         potential_kj_mol=frame['potential_kj_mol'],kinetic_kj_mol=frame['kinetic_kj_mol']))
if count==0:
    raise ValueError('analysis selection retained no frames')
rmsf=np.sqrt(np.maximum(np.sum(m2/count,axis=1),0))
Path('analysis.json').write_text(json.dumps(dict(schema_version=1,frames=rows,
    atoms=[dict(index=i,rmsf_nm=float(value)) for i,value in zip(indices,rmsf)],
    frame_count=count,reference_step=rows[0]['step'],alignment=options['alignment'],weighting=options['weighting']),allow_nan=False))
'''


class TrajectoryAnalysis(StepAdapter):
    name = 'bio_mdanalysis'
    category = 'structure'
    description = 'Analyze a checked trajectory with explicit atom/frame selection, optional MDAnalysis fitting, RMSD/RMSF and radius of gyration'
    produces_geometry = False
    consumes = ('trajectory', 'simulation_system')
    wires = (('trajectory', 'trajectory_file'), ('simulation_system', 'system_file'))
    produces = ('trajectory_report', 'json_report')
    required_artifacts = produces
    params = (ParamSpec('trajectory_file', 'path', required=True), ParamSpec('system_file', 'path', required=True),
              ParamSpec('expected_version', 'str', required=True), ParamSpec('atom_indices', 'list', required=True),
              ParamSpec('alignment', 'str', required=True, enum=('none', 'fit')),
              ParamSpec('weighting', 'str', required=True, enum=('uniform', 'mass')),
              ParamSpec('discard_frames', 'int', required=True), ParamSpec('stride', 'int', required=True),
              ParamSpec('periodic_policy', 'str', default='reject', enum=('reject', 'as_stored')),
              ParamSpec('timeout', 'float', default=600.), ParamSpec('python', 'path', default=sys.executable))
    data_keys = (DataKeySpec('runtime', 'dict'), DataKeySpec('frame_count', 'int'), DataKeySpec('frames', 'list'))
    example_kwargs = dict(trajectory_file='trajectory.json', system_file='simulation_system.json', expected_version='2.9.0',
                          atom_indices=[0, 1, 2], alignment='fit', weighting='mass', discard_frames=0, stride=1)

    def validate_params(self, **kwargs):
        for p in self.params:
            if p.required and (kwargs.get(p.name) is None or kwargs.get(p.name) == ''):
                raise ValueError(f'{p.name} is required')
        indices = kwargs['atom_indices']
        if not isinstance(indices, list) or not indices or any(type(i) is not int or i < 0 for i in indices) or len(set(indices)) != len(indices):
            raise ValueError('atom_indices must be a unique nonempty list of zero-based indices')
        if kwargs['alignment'] not in ('none', 'fit') or kwargs['weighting'] not in ('uniform', 'mass'):
            raise ValueError('explicit alignment and weighting policies required')
        if kwargs.get('periodic_policy', 'reject') not in ('reject', 'as_stored'):
            raise ValueError('invalid periodic analysis policy')
        if kwargs['alignment'] == 'fit' and len(indices) < 3:
            raise ValueError('fitting requires at least three selected atoms')
        for name, minimum in [('discard_frames', 0), ('stride', 1)]:
            if type(kwargs[name]) is not int or kwargs[name] < minimum:
                raise ValueError(f'invalid {name}')
        if not number(kwargs.get('timeout', 600.)) or kwargs.get('timeout', 600.) <= 0:
            raise ValueError('timeout must be positive and finite')

    def runtime_preflight(self, **kwargs):
        return probe_distribution('mdanalysis', python=kwargs.get('python', sys.executable), expected_version=kwargs['expected_version'],
                                  work_dir=Path.cwd(), allow_install=kwargs.get('_allow_install', False))

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start, root = time.monotonic(), Path(work_dir)
        try:
            self.validate_params(**kwargs)
            trajectory_path = Path(kwargs['trajectory_file']).resolve()
            details = validate_trajectory(trajectory_path)
            if details['periodic'] and kwargs.get('periodic_policy', 'reject') != 'as_stored':
                raise ValueError('periodic trajectory analysis requires explicit periodic_policy=as_stored; no unwrapping is performed')
            trajectory = read_json(trajectory_path)
            if system_identity(kwargs['system_file']) != trajectory['system_identity']:
                raise ValueError('analysis system differs from trajectory system')
            if max(kwargs['atom_indices']) >= details['particles'] or kwargs['discard_frames'] >= details['frames']:
                raise ValueError('atom/frame selection outside trajectory')
            bundle_path = Path(kwargs['system_file']).resolve()
            bundle = read_json(bundle_path)
            particles = ET.parse(child_path(bundle_path.parent, bundle['files']['system']['path'])).getroot().findall('./Particles/Particle')
            masses = [float(particles[i].get('mass')) for i in kwargs['atom_indices']]
            frames = checked_file(trajectory_path.parent, trajectory['frames'])
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('analysis requires an empty step directory')
            shutil.copyfile(frames, root / 'frames.jsonl')
            options = {key: kwargs[key] for key in ('atom_indices', 'alignment', 'weighting', 'discard_frames', 'stride')}
            options['masses'] = masses
            options['periodic_policy'] = kwargs.get('periodic_policy', 'reject')
            runtime = run_worker(_WORKER, options, {'mdanalysis': kwargs['expected_version']}, python=kwargs.get('python', sys.executable),
                                 work_dir=root, timeout=kwargs.get('timeout', 600.))
            data = read_json(root / 'analysis.json')
            expected = [(f['step'], f['time_ps'], f['potential_kj_mol'], f['kinetic_kj_mol']) for i, f in enumerate(iter_frames(frames))
                        if i >= kwargs['discard_frames'] and (i-kwargs['discard_frames']) % kwargs['stride'] == 0]
            rows, atoms = data.get('frames'), data.get('atoms')
            if data.get('schema_version') != 1 or not isinstance(rows, list) or any(not isinstance(r, dict) for r in rows) or [(r.get('step'), r.get('time_ps'), r.get('potential_kj_mol'), r.get('kinetic_kj_mol')) for r in rows] != expected or data.get('frame_count') != len(expected):
                raise ValueError('analysis must account for every selected frame in order')
            if not isinstance(atoms, list) or any(not isinstance(a, dict) for a in atoms) or [a.get('index') for a in atoms] != kwargs['atom_indices']:
                raise ValueError('analysis atom mapping differs from selection')
            if any(not number(row.get(key)) or row[key] < 0 for row in rows for key in ('rmsd_nm', 'radius_gyration_nm')) or any(not number(a.get('rmsf_nm')) or a['rmsf_nm'] < 0 for a in atoms):
                raise ValueError('invalid trajectory statistics')
            if data.get('alignment') != kwargs['alignment'] or data.get('weighting') != kwargs['weighting'] or data.get('reference_step') != expected[0][0]:
                raise ValueError('analysis policy or reference mismatch')
            data.update(runtime=runtime, selection=options, trajectory_sha256=file_hash(trajectory_path), frames_sha256=file_hash(frames),
                        system_identity=trajectory['system_identity'], units=dict(distance='nanometer', time='picosecond', energy='kilojoule_per_mole'),
                        interpretation='Selected-atom statistics relative to first retained frame; no periodic unwrapping, automatic equilibration detection or scientific convergence claim')
            write_json(root / 'analysis.json', data)
            return self._make_result(self.name, StepStatus.SUCCESS, root, start, data=data,
                                     artifacts={'trajectory_report': root / 'analysis.json', 'json_report': root / 'analysis.json'})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, TypeError, KeyError, OSError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, root, start, error=error, error_kind=kind)


register(TrajectoryAnalysis())
