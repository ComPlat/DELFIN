"""ESMFold v1 inference against checksum-verified local model assets."""
from __future__ import annotations

import math
from pathlib import Path
import re
import sys
import time
import subprocess

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec, DataKeySpec
from delfin.tools._types import StepStatus, ErrorKind
from delfin.tools.bio import read_fasta, read_structure
from delfin.tools.bio.artifacts import read_json, write_json, file_hash, validate_artifact
from delfin.tools.bio.assets import validate_assets
from delfin.tools.bio.structures import import_predictions
from delfin.tools.bio.optional_runtime import run_worker, probe_distribution


_REQUIRED = ('checkpoints/esmfold_3B_v1.pt', 'checkpoints/esm2_t36_3B_UR50D.pt',
             'checkpoints/esm2_t36_3B_UR50D-contact-regression.pt')
_AMINO_ACIDS = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL UNK'.split(),
                        'ARNDCQEGHILKMFPSTWYVX'))
_WORKER = '''import json,sys,random
from pathlib import Path
import numpy as np
import torch
import esm
options = json.load(open(sys.argv[1]))
random.seed(options['seed']); np.random.seed(options['seed']); torch.manual_seed(options['seed'])
torch.set_num_threads(options['cores'])
device = torch.device(options['device'])
if device.type == 'cuda':
    if not torch.cuda.is_available() or device.index >= torch.cuda.device_count():
        raise ValueError('requested CUDA device unavailable')
    torch.cuda.manual_seed_all(options['seed'])
torch.hub.set_dir(options['assets_root'])
def refuse_download(*args, **kwargs):
    raise RuntimeError('model asset is missing; ESMFold automatic downloads are disabled')
torch.hub.download_url_to_file = refuse_download
model = esm.pretrained.esmfold_v1().eval()
model.set_chunk_size(options['chunk_size'])
if device.type == 'cpu':
    model.esm.float()
model = model.to(device)
with torch.no_grad():
    output = model.infer([options['sequence']], num_recycles=options['num_recycles'])
    pdbs = model.output_to_pdb(output)
if len(pdbs) != 1:
    raise ValueError('ESMFold did not return exactly one model')
Path('model.pdb').write_text(pdbs[0])
with open('confidence.json', 'w') as stream:
    json.dump({'engine': 'esmfold', 'mean_plddt': float(output['mean_plddt'][0]),
               'ptm': float(output['ptm'][0]), 'b_factor_interpretation': 'plddt',
               'device': str(device), 'torch_version': torch.__version__}, stream, allow_nan=False)
'''


class EsmFold(StepAdapter):
    name = 'esmfold_predict'
    category = 'structure'
    description = 'Run one ESMFold v1 protein prediction with exact runtime versions and verified offline model assets'
    produces_geometry = False
    consumes = ('esmfold_input',)
    wires = (('esmfold_input', 'request_file'),)
    produces = ('prediction_bundle', 'json_report')
    required_artifacts = produces
    params = (ParamSpec('request_file', 'path', required=True), ParamSpec('assets_file', 'path', required=True),
              ParamSpec('expected_version', 'str', required=True), ParamSpec('torch_version', 'str', required=True),
              ParamSpec('device', 'str', required=True), ParamSpec('seed', 'int', required=True),
              ParamSpec('chunk_size', 'int', default=128), ParamSpec('num_recycles', 'int', default=4),
              ParamSpec('timeout', 'float', default=3600.0), ParamSpec('python', 'path', default=sys.executable))
    data_keys = (DataKeySpec('assets', 'dict'), DataKeySpec('runtime', 'dict'))
    example_kwargs = dict(request_file='input.fasta', assets_file='assets.json', expected_version='2.0.0',
                          torch_version='2.2.1', device='cpu', seed=7)

    def validate_params(self, **kwargs):
        for p in self.params:
            if p.required and (kwargs.get(p.name) is None or kwargs.get(p.name) == ''):
                raise ValueError(f'{p.name} is required')
        if not isinstance(kwargs['device'], str) or not re.fullmatch(r'cpu|cuda:[0-9]+', kwargs['device']):
            raise ValueError('device must be cpu or an explicit cuda:N')
        for key, default, minimum in (('seed', None, 0), ('chunk_size', 128, 1), ('num_recycles', 4, 1)):
            value = kwargs.get(key, default)
            if type(value) is not int or not minimum <= value < 2**31:
                raise ValueError(f'invalid {key}')
        timeout = kwargs.get('timeout', 3600.0)
        if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError('timeout must be positive and finite')

    def runtime_preflight(self, **kwargs):
        assets = validate_assets(kwargs['assets_file'], engine='esmfold', required_paths=_REQUIRED)
        python = kwargs.get('python', sys.executable)
        versions = {name: probe_distribution(name, python=python, expected_version=version, work_dir=Path.cwd(), discover=False)
                    for name, version in (('fair-esm', kwargs['expected_version']), ('torch', kwargs['torch_version']))}
        return {'assets': assets, 'versions': versions}

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start = time.monotonic()
        root = Path(work_dir)
        try:
            self.validate_params(**kwargs)
            if type(cores) is not int or cores < 1:
                raise ValueError('cores must be positive')
            records = read_fasta(kwargs['request_file'])
            if len(records) != 1:
                raise ValueError('ESMFold adapter requires exactly one protein sequence')
            assets = validate_assets(kwargs['assets_file'], engine='esmfold', required_paths=_REQUIRED)
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('prediction requires an empty step directory')
            options = dict(sequence=records[0].sequence, assets_root=assets['root'], device=kwargs['device'],
                           seed=kwargs['seed'], cores=cores, chunk_size=kwargs.get('chunk_size', 128),
                           num_recycles=kwargs.get('num_recycles', 4))
            runtime = run_worker(_WORKER, options, {'fair-esm': kwargs['expected_version'], 'torch': kwargs['torch_version']},
                                 python=kwargs.get('python', sys.executable), work_dir=root, timeout=kwargs.get('timeout', 3600.0))
            model = root / 'model.pdb'
            validate_artifact(model, 'macromolecule')
            structure = read_structure(model)
            if len(structure.models) != 1 or len(structure.first_model.residues) != len(records[0].sequence):
                raise ValueError('predicted structure does not contain the expected single-model residue count')
            if len(structure.first_model.chains) != 1:
                raise ValueError('monomer prediction returned multiple chains')
            residues = structure.first_model.residues
            if ''.join(_AMINO_ACIDS.get(residue.name, '?') for residue in residues) != records[0].sequence:
                raise ValueError('predicted protein sequence differs from the input')
            if any(not {'N', 'CA', 'C', 'O'} <= {atom.name for atom in residue.atoms} for residue in residues):
                raise ValueError('predicted model has an incomplete backbone')
            confidence = read_json(root / 'confidence.json')
            for key, maximum in (('mean_plddt', 100), ('ptm', 1)):
                value = confidence.get(key)
                if isinstance(value, bool) or not isinstance(value, (int, float)) or not 0 <= value <= maximum:
                    raise ValueError(f'invalid ESMFold confidence: {key}')
            bundle_dir = root / 'bundle'
            bundle_dir.mkdir()
            bundle = import_predictions(root, [{'id': 'model_1', 'path': 'model.pdb', 'confidence': 'confidence.json'}],
                                        bundle_dir, engine='esmfold', engine_version=kwargs['expected_version'])
            value = read_json(bundle)
            value.update(provenance_origin='engine_run', runtime=runtime, assets=assets)
            write_json(bundle, value)
            data = {'runtime': runtime, 'assets': assets, 'sequence_sha256': file_hash(Path(kwargs['request_file'])),
                    'settings': options, 'confidence': confidence, 'scientific_validation': 'not_performed'}
            report = root / 'prediction_report.json'
            write_json(report, data)
            return self._make_result(self.name, StepStatus.SUCCESS, root, start, data=data,
                                     artifacts={'prediction_bundle': bundle, 'json_report': report, 'log': root / 'engine.log'})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, TypeError, KeyError, OSError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, root, start, error=error, error_kind=kind)


register(EsmFold())
