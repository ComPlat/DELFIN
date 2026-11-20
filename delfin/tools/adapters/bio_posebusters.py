"""Optional PoseBusters redocking checks, preserving every check result."""
from __future__ import annotations

import csv
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._spec import ParamSpec, DataKeySpec
from delfin.tools._types import ErrorKind, StepStatus
from delfin.tools.bio.artifacts import file_hash, validate_artifact, write_json
from delfin.tools.bio.optional_runtime import probe_distribution, run_console


class PoseBusters(StepAdapter):
    name = 'posebusters_validate'
    category = 'structure'
    description = 'Run optional PoseBusters redocking checks with an explicit reference; retain all raw checks without declaring affinity'
    produces_geometry = False
    produces = ('pose_validation', 'json_report')
    required_artifacts = produces
    consumes = ('pose_set', 'ligand', 'prepared_receptor')
    wires = (('pose_set', 'poses_file'), ('ligand', 'reference_file'), ('prepared_receptor', 'receptor_file'))
    requires_python = ('rdkit',)
    params = (ParamSpec('poses_file', 'path', required=True), ParamSpec('reference_file', 'path', required=True),
              ParamSpec('receptor_file', 'path', required=True), ParamSpec('expected_version', 'str', required=True),
              ParamSpec('python', 'path', default=sys.executable), ParamSpec('timeout', 'float', default=600.0))
    data_keys = (DataKeySpec('checks', 'list'), DataKeySpec('runtime', 'dict'))
    example_kwargs = dict(poses_file='poses.sdf', reference_file='reference.sdf', receptor_file='receptor.pdb', expected_version='0.5.1')

    def validate_params(self, **kwargs):
        for param in self.params:
            if param.required and not kwargs.get(param.name):
                raise ValueError(f'{param.name} is required')
        value = kwargs.get('timeout', 600.0)
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or value <= 0:
            raise ValueError('timeout must be positive and finite')

    def runtime_preflight(self, **kwargs):
        return probe_distribution('posebusters', python=kwargs.get('python', sys.executable),
                                  expected_version=kwargs['expected_version'], work_dir=Path.cwd())

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        from rdkit import Chem
        from delfin.tools.bio import read_structure
        from delfin.tools.bio.processing import load_ligands
        start = time.monotonic()
        root = Path(work_dir)
        try:
            self.validate_params(**kwargs)
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('PoseBusters requires an empty step directory')
            hashes = {}
            for key, kind in zip(('poses_file', 'reference_file', 'receptor_file'), self.consumes):
                hashes[key] = validate_artifact(Path(kwargs[key]), kind)['sha256']
            receptor = read_structure(kwargs['receptor_file'])
            if receptor.source_format != 'pdb' or len(receptor.models) != 1:
                raise ValueError('PoseBusters integration requires a single prepared PDB receptor')
            poses, reference = load_ligands(kwargs['poses_file']), load_ligands(kwargs['reference_file'])
            if len(reference) != 1:
                raise ValueError('PoseBusters redocking requires exactly one explicit reference ligand')
            if any(not mol.GetNumConformers() for mol in poses + reference):
                raise ValueError('PoseBusters requires existing coordinates in a shared reference frame')
            for name, molecules in (('poses.sdf', poses), ('reference.sdf', reference)):
                with Chem.SDWriter(str(root / name)) as writer:
                    for mol in molecules:
                        writer.write(mol)
            shutil.copyfile(kwargs['receptor_file'], root / 'receptor.pdb')
            argv = ['poses.sdf', '-l', 'reference.sdf', '-p', 'receptor.pdb', '--outfmt', 'csv', '--output', 'checks.csv', '--full-report']
            runtime = run_console('posebusters', 'bust', argv, python=kwargs.get('python', sys.executable),
                                  expected_version=kwargs['expected_version'], work_dir=root, timeout=kwargs.get('timeout', 600.0))
            output = root / 'checks.csv'
            validate_artifact(output, 'csv')
            with output.open(newline='') as stream:
                reader = csv.DictReader(stream)
                if not {'file', 'molecule', 'position'} <= set(reader.fieldnames or ()):
                    raise ValueError('unsupported PoseBusters CSV schema')
                rows = list(reader)
            if len(rows) != len(poses) or any(None in row or any(v is None for v in row.values()) for row in rows):
                raise ValueError('incomplete PoseBusters output: expected one complete row per pose')
            if any(len(row) <= 3 for row in rows):
                raise ValueError('PoseBusters output contains no checks')
            positions = []
            for row in rows:
                source = Path(row['file'])
                if not source.is_absolute():
                    source = root / source
                if source.resolve() != (root / 'poses.sdf').resolve():
                    raise ValueError('PoseBusters report refers to an unexpected pose file')
                if not row['position'].isdigit():
                    raise ValueError('PoseBusters report contains an invalid pose position')
                positions.append(int(row['position']))
            if sorted(positions) != list(range(len(poses))):
                raise ValueError('PoseBusters report must cover every pose position exactly once')
            data = dict(runtime=runtime, checks=rows, input_sha256=hashes, output_sha256=file_hash(output),
                        mode='redock', coordinate_frame='caller-supplied; no alignment',
                        interpretation='Execution success means a complete report, not that poses passed; raw checks retained, no affinity claim')
            report = root / 'validation_report.json'
            write_json(report, data)
            return self._make_result(self.name, StepStatus.SUCCESS, root, start, data=data,
                                     artifacts={'pose_validation': report, 'json_report': report, 'raw_checks': output, 'log': root / 'engine.log'})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, KeyError, TypeError, OSError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, root, start, error=error, error_kind=kind)


register(PoseBusters())
