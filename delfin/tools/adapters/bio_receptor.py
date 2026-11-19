"""Restricted PDBFixer/OpenMM preparation in an explicitly versioned environment."""
from __future__ import annotations

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
from delfin.tools.bio import read_structure
from delfin.tools.bio.artifacts import file_hash, read_json, validate_artifact, write_json
from delfin.tools.bio.optional_runtime import probe_distribution, run_worker


_WORKER = '''import json, random, sys
import numpy as np
from pdbfixer import PDBFixer
from openmm import Platform
from openmm.app import PDBFile
options = json.load(open(sys.argv[1]))
random.seed(options['seed'])
np.random.seed(options['seed'])
fixer = PDBFixer(filename='input.pdb', platform=Platform.getPlatformByName('Reference'))
fixer.findMissingResidues()
if fixer.missingResidues:
    raise ValueError('missing residues detected; explicit structural modelling is required')
fixer.findNonstandardResidues()
if fixer.nonstandardResidues:
    raise ValueError('nonstandard residues require explicit chemistry handling')
if options['waters'] == 'drop':
    fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
missing = [{'chain': r.chain.id, 'residue': r.id, 'name': r.name,
            'atoms': [a.name for a in atoms]} for r, atoms in fixer.missingAtoms.items()]
terminals = [{'chain': r.chain.id, 'residue': r.id, 'atoms': atoms}
             for r, atoms in fixer.missingTerminals.items()]
if (missing or terminals) and options['missing_atoms'] == 'reject':
    raise ValueError('missing heavy/terminal atoms detected; repair must be explicitly requested')
if options['missing_atoms'] == 'repair':
    fixer.addMissingAtoms(seed=options['seed'])
fixer.addMissingHydrogens(pH=options['ph'])
with open('receptor.pdb', 'w') as stream:
    PDBFile.writeFile(fixer.topology, fixer.positions, stream, keepIds=True)
with open('worker_report.json', 'w') as stream:
    json.dump({'missing_heavy_atoms': missing, 'missing_terminal_atoms': terminals,
               'platform': 'Reference', 'ph': options['ph'], 'seed': options['seed'],
               'sequence_records_present': bool(fixer.sequences)}, stream, allow_nan=False)
'''


def atom_records(structure):
    rows = {}
    for residue in structure.first_model.residues:
        for atom in residue.atoms:
            key = (residue.chain_id, residue.seq_id, residue.insertion_code, residue.name, atom.name)
            if key in rows:
                raise ValueError('duplicate residue/atom identity is ambiguous')
            rows[key] = (atom, residue)
    return rows


class ReceptorPreparation(StepAdapter):
    name = 'bio_receptor_prepare'
    description = 'Prepare one canonical protein PDB using version-pinned PDBFixer/OpenMM, explicit pH/water/repair policy and atom mapping'
    category = 'structure'
    produces_geometry = False
    consumes = ('macromolecule',)
    wires = (('macromolecule', 'structure_file'),)
    produces = ('prepared_receptor', 'json_report')
    required_artifacts = produces
    data_keys = (DataKeySpec('runtime', 'dict'), DataKeySpec('atom_mapping', 'list'))
    params = (ParamSpec('structure_file', 'path', required=True),
              ParamSpec('pdbfixer_version', 'str', required=True), ParamSpec('openmm_version', 'str', required=True),
              ParamSpec('ph', 'float', required=True), ParamSpec('seed', 'int', required=True),
              ParamSpec('waters', 'str', required=True, enum=('keep', 'drop')),
              ParamSpec('missing_atoms', 'str', required=True, enum=('reject', 'repair')),
              ParamSpec('python', 'path', default=sys.executable), ParamSpec('timeout', 'float', default=600.0))
    example_kwargs = dict(structure_file='protein.pdb', pdbfixer_version='1.12', openmm_version='8.4',
                          ph=7.4, seed=7, waters='drop', missing_atoms='reject')

    def validate_params(self, **kwargs):
        for param in self.params:
            if param.required and (kwargs.get(param.name) is None or kwargs.get(param.name) == ''):
                raise ValueError(f'{param.name} is required')
        for key, low, high in (('ph', 0, 14), ('timeout', 0, float('inf'))):
            value = kwargs.get(key, 600.0)
            if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not low <= value <= high:
                raise ValueError(f'invalid {key}')
        if kwargs.get('timeout', 600.0) <= 0:
            raise ValueError('timeout must be positive')
        if type(kwargs['seed']) is not int or not 0 <= kwargs['seed'] < 2**31:
            raise ValueError('seed must be an integer in 0..2**31-1')
        if kwargs['waters'] not in ('keep', 'drop') or kwargs['missing_atoms'] not in ('reject', 'repair'):
            raise ValueError('explicit waters and missing_atoms policies are required')

    def runtime_preflight(self, **kwargs):
        return {name: probe_distribution(name, python=kwargs.get('python', sys.executable),
                                        expected_version=kwargs[name + '_version'], work_dir=Path.cwd())
                for name in ('pdbfixer', 'openmm')}

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start = time.monotonic()
        root = Path(work_dir)
        try:
            self.validate_params(**kwargs)
            source = Path(kwargs['structure_file']).resolve()
            validate_artifact(source, 'macromolecule')
            original = read_structure(source)
            if original.source_format != 'pdb' or len(original.models) != 1:
                raise ValueError('preparation currently requires one explicitly selected PDB model')
            standard = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
            residues = original.first_model.residues
            if any(r.altlocs or r.is_modified or (r.name not in standard and not r.is_water) for r in residues):
                raise ValueError('resolve altlocs; modified residues, cofactors, metals and nonprotein residues require explicit preparation')
            before = atom_records(original)
            if any(atom.element in ('H', 'D', 'T') for atom, _ in before.values()):
                raise ValueError('this preparation subset requires a hydrogen-free input; existing protonation needs a separate policy')
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('preparation requires an empty step directory')
            shutil.copyfile(source, root / 'input.pdb')
            options = {key: kwargs[key] for key in ('ph', 'seed', 'waters', 'missing_atoms')}
            runtime = run_worker(_WORKER, options, {name: kwargs[name + '_version'] for name in ('pdbfixer', 'openmm')},
                                 python=kwargs.get('python', sys.executable), work_dir=root, timeout=kwargs.get('timeout', 600.0))
            output = root / 'receptor.pdb'
            validate_artifact(output, 'prepared_receptor')
            after = atom_records(read_structure(output))
            allowed_residues = {key[:4] for key, (_, residue) in before.items()
                                if not (kwargs['waters'] == 'drop' and residue.is_water)}
            if any(key[:4] not in allowed_residues for key in after):
                raise ValueError('preparation introduced an unexpected residue')
            mapping, removed = [], []
            for key, (atom, residue) in before.items():
                if kwargs['waters'] == 'drop' and residue.is_water:
                    if key in after:
                        raise ValueError('requested water removal was not performed')
                    removed.append({'identity': key, 'source_serial': atom.serial, 'reason': 'water policy'})
                    continue
                if key not in after or after[key][0].element != atom.element:
                    raise ValueError('preparation lost or renamed an original atom')
                new = after[key][0]
                if any(abs(a - b) > 0.0011 for a, b in zip(atom.position, new.position)):
                    raise ValueError('preparation moved an original heavy atom beyond PDB rounding tolerance')
                mapping.append({'identity': key, 'source_serial': atom.serial, 'output_serial': new.serial})
            added = [{'identity': key, 'output_serial': atom.serial, 'element': atom.element}
                     for key, (atom, _) in after.items() if key not in before]
            if not any(item['element'] == 'H' for item in added):
                raise ValueError('preparation did not produce hydrogen atoms')
            if kwargs['missing_atoms'] == 'reject' and any(item['element'] != 'H' for item in added):
                raise ValueError('heavy atoms added despite the reject policy')
            data = {'runtime': runtime, 'policies': options, 'source_sha256': file_hash(source),
                    'output_sha256': file_hash(output), 'atom_mapping': mapping, 'added_atoms': added, 'removed_atoms': removed,
                    'engine_report': read_json(root / 'worker_report.json'),
                    'method': 'PDBFixer repair and OpenMM default-residue pH hydrogen assignment',
                    'metadata': 'Output PDB contains the prepared topology and coordinates; original annotations remain in the copied input.pdb',
                    'limits': 'No pKa calculation, loop modelling, biological assembly inference, force-field parameterization or scientific validation. Missing-residue detection depends on sequence records.'}
            report = root / 'preparation_report.json'
            write_json(report, data)
            return self._make_result(self.name, StepStatus.SUCCESS, root, start, data=data,
                                     artifacts={'prepared_receptor': output, 'json_report': report, 'log': root / 'engine.log'})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, KeyError, TypeError, OSError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, root, start, error=error, error_kind=kind)


register(ReceptorPreparation())
