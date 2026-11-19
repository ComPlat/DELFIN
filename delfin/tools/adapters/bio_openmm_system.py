"""Explicit protein-only ff14SB parameterization into a checked OpenMM bundle."""
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
from delfin.tools.bio import read_structure
from delfin.tools.bio.artifacts import file_hash, write_json
from delfin.tools.bio.optional_runtime import probe_distribution, run_worker
from delfin.tools.bio.simulation import validate_system


_WORKER = '''import json,sys,shutil
from pathlib import Path
import xml.etree.ElementTree as ET
import openmm
from openmm import app, unit
options=json.load(open(sys.argv[1]))
pdb=app.PDBFile('input.pdb')
if pdb.topology.getPeriodicBoxVectors() is not None:
    raise ValueError('vacuum parameterization refuses a periodic input box')
source=Path(app.__file__).parent/'data'/'amber14'/'protein.ff14SB.xml'
if ET.parse(source).getroot().find('Include') is not None:
    raise ValueError('force field includes require an explicitly versioned dependency bundle')
shutil.copyfile(source,'forcefield.xml')
forcefield=app.ForceField('forcefield.xml')
system=forcefield.createSystem(pdb.topology,nonbondedMethod=app.NoCutoff,
    constraints=app.HBonds if options['constraints']=='HBonds' else None,
    rigidWater=False,removeCMMotion=False)
if system.getNumParticles()!=pdb.topology.getNumAtoms() or any(system.isVirtualSite(i) for i in range(system.getNumParticles())):
    raise ValueError('extra particles/virtual sites require explicit mapping support')
Path('system.xml').write_text(openmm.XmlSerializer.serialize(system))
# Round-trip the engine serialization before declaring the bundle usable.
restored=openmm.XmlSerializer.deserialize(Path('system.xml').read_text())
if restored.getNumParticles()!=system.getNumParticles():
    raise ValueError('system serialization did not round-trip')
atoms=[dict(index=a.index,name=a.name,element=a.element.symbol,chain=a.residue.chain.id,
            residue_id=a.residue.id,residue_name=a.residue.name,insertion_code=a.residue.insertionCode)
       for a in pdb.topology.atoms()]
bonds=[[a.index,b.index] for a,b in pdb.topology.bonds()]
for name,value in [('topology.json',{'atoms':atoms,'bonds':bonds}),
                   ('positions.json',{'unit':'nanometer','values':[list(p) for p in pdb.positions.value_in_unit(unit.nanometer)]})]:
    Path(name).write_text(json.dumps(value,allow_nan=False))
'''


class OpenMMSystem(StepAdapter):
    name = 'bio_openmm_parameterize'
    description = 'Parameterize an explicitly hydrogenated canonical protein in vacuum with ff14SB; export a hashed OpenMM system bundle'
    category = 'structure'
    produces_geometry = False
    consumes = ('complex_structure',)
    wires = (('complex_structure', 'structure_file'),)
    produces = ('simulation_system', 'json_report')
    required_artifacts = produces
    params = (ParamSpec('structure_file', 'path', required=True), ParamSpec('expected_version', 'str', required=True),
              ParamSpec('forcefield', 'str', required=True, enum=('amber14/protein.ff14SB.xml',)),
              ParamSpec('environment', 'str', required=True, enum=('vacuum',)),
              ParamSpec('constraints', 'str', required=True, enum=('none', 'HBonds')),
              ParamSpec('expected_charge', 'int', required=True),
              ParamSpec('python', 'path', default=sys.executable), ParamSpec('timeout', 'float', default=600.0))
    data_keys = (DataKeySpec('runtime', 'dict'), DataKeySpec('validation', 'dict'))
    example_kwargs = dict(structure_file='protein.pdb', expected_version='8.4.0', forcefield='amber14/protein.ff14SB.xml',
                          environment='vacuum', constraints='HBonds', expected_charge=0)

    def validate_params(self, **kwargs):
        for param in self.params:
            if param.required and (kwargs.get(param.name) is None or kwargs.get(param.name) == ''):
                raise ValueError(f'{param.name} is required')
        if kwargs['forcefield'] != 'amber14/protein.ff14SB.xml' or kwargs['environment'] != 'vacuum' or kwargs['constraints'] not in ('none', 'HBonds'):
            raise ValueError('explicit ff14SB/vacuum/constraints policies required')
        if type(kwargs['expected_charge']) is not int:
            raise ValueError('expected_charge must be an integer in elementary charge units')
        timeout = kwargs.get('timeout', 600.0)
        if type(timeout) not in (int, float) or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError('timeout must be positive and finite')

    def runtime_preflight(self, **kwargs):
        return probe_distribution('openmm', python=kwargs.get('python', sys.executable), expected_version=kwargs['expected_version'],
                                  work_dir=Path.cwd(), allow_install=kwargs.get('_allow_install', False))

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start, root = time.monotonic(), Path(work_dir)
        try:
            self.validate_params(**kwargs)
            structure = read_structure(kwargs['structure_file'])
            allowed = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
            if structure.source_format != 'pdb' or len(structure.models) != 1 or any(r.name not in allowed or r.altlocs or r.is_modified or r.is_hetero for r in structure.first_model.residues):
                raise ValueError('parameterization currently requires a single canonical protein PDB without solvent, ligands, modifications or altlocs')
            if not any(a.element == 'H' for a in structure.first_model.atoms):
                raise ValueError('protein must already contain explicit hydrogens; no protonation changes are made')
            if any(line.startswith('CRYST1') for line in Path(kwargs['structure_file']).read_text().splitlines()):
                raise ValueError('vacuum parameterization refuses a periodic input box')
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('parameterization requires an empty step directory')
            shutil.copyfile(kwargs['structure_file'], root / 'input.pdb')
            runtime = run_worker(_WORKER, {'constraints': kwargs['constraints']}, {'openmm': kwargs['expected_version']},
                                 python=kwargs.get('python', sys.executable), work_dir=root, timeout=kwargs.get('timeout', 600.0))
            files = {key: {'path': filename, 'sha256': file_hash(root / filename)} for key, filename in
                     dict(structure='input.pdb', system='system.xml', topology='topology.json', positions='positions.json', forcefield='forcefield.xml').items()}
            bundle = dict(schema_version=1, engine='openmm', files=files, particle_count=len(structure.first_model.atoms),
                          forcefield_name=kwargs['forcefield'], environment=kwargs['environment'], constraints=kwargs['constraints'],
                          nonbonded_method='NoCutoff', expected_charge=kwargs['expected_charge'], runtime=runtime,
                          provenance_origin='engine_run', scientific_validation='not_performed')
            output = root / 'simulation_system.json'
            write_json(output, bundle)
            validation = validate_system(output)
            report = dict(runtime=runtime, validation=validation, scope='Protein-only vacuum parameterization; no minimization, dynamics, solvent or ligand parameters')
            write_json(root / 'report.json', report)
            return self._make_result(self.name, StepStatus.SUCCESS, root, start, data=report,
                                     artifacts={'simulation_system': output, 'json_report': root / 'report.json'})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, TypeError, KeyError, OSError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, root, start, error=error, error_kind=kind)


register(OpenMMSystem())
