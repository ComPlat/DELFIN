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
from delfin.tools.bio.artifacts import file_hash, write_json, read_json, child_path
from delfin.tools.bio.ligand_parameters import copy_ligand_parameters, complex_ligand_mapping
from delfin.tools.bio.optional_runtime import probe_distribution, run_worker
from delfin.tools.bio.simulation import validate_system, validate_solvation


_WORKER = '''import json,sys,shutil,random
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
forcefield_files=['forcefield.xml']
if options['environment']=='tip3p_pme':
    solvent=Path(app.__file__).parent/'data'/'amber14'/'tip3p.xml'
    if ET.parse(solvent).getroot().find('Include') is not None:
        raise ValueError('solvent force field must be self-contained')
    shutil.copyfile(solvent,'solvent_forcefield.xml'); forcefield_files.append('solvent_forcefield.xml')
ligand_template=None
if options.get('ligand_forcefield'):
    ligand_template=ET.parse(options['ligand_forcefield']).getroot().find('./Residues/Residue').get('name')
    forcefield_files.append(options['ligand_forcefield'])
forcefield=app.ForceField(*forcefield_files)
def residue_templates(topology):
    return {r:ligand_template for r in topology.residues() if r.name=='LIG'} if ligand_template else {}
settings=options.get('solvation')
nonbonded={'nonbondedMethod':app.NoCutoff}
if settings:
    random.seed(settings['seed'])
    model=app.Modeller(pdb.topology,pdb.positions)
    model.addSolvent(forcefield,model='tip3p',boxSize=openmm.Vec3(*settings['box_size_nm'])*unit.nanometer,
        ionicStrength=settings['ionic_strength_m']*unit.molar,neutralize=settings['neutralize'],
        positiveIon=settings['positive_ion'],negativeIon=settings['negative_ion'],residueTemplates=residue_templates(model.topology))
    if model.topology.getNumAtoms()>99999 or any(len(c.id)!=1 for c in model.topology.chains()) or any(not r.id.isdigit() or not 1<=int(r.id)<=9999 for r in model.topology.residues()):
        raise ValueError('solvated system exceeds the supported PDB identity limits')
    with open('solvated.pdb','w') as stream:
        app.PDBFile.writeFile(model.topology,model.positions,stream,keepIds=True)
    # Use the stored coordinate precision consistently in all bundle members.
    pdb=app.PDBFile('solvated.pdb')
    nonbonded=dict(nonbondedMethod=app.PME,nonbondedCutoff=settings['cutoff_nm']*unit.nanometer,ewaldErrorTolerance=settings['ewald_tolerance'])
system=forcefield.createSystem(pdb.topology,**nonbonded,
    constraints=app.HBonds if options['constraints']=='HBonds' else None,
    rigidWater=False,removeCMMotion=False,residueTemplates=residue_templates(pdb.topology))
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
    description = 'Parameterize an explicitly hydrogenated canonical protein with ff14SB in vacuum or an explicit TIP3P/PME box; export a checked system bundle'
    category = 'structure'
    produces_geometry = False
    consumes = ('complex_structure',)
    wires = (('complex_structure', 'structure_file'),)
    produces = ('simulation_system', 'json_report')
    required_artifacts = produces
    params = (ParamSpec('structure_file', 'path', required=True), ParamSpec('expected_version', 'str', required=True),
              ParamSpec('forcefield', 'str', required=True, enum=('amber14/protein.ff14SB.xml',)),
              ParamSpec('environment', 'str', required=True, enum=('vacuum', 'tip3p_pme')),
              ParamSpec('solvation', 'dict'),
              ParamSpec('ligand_parameters_file', 'path'),
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
        if kwargs['forcefield'] != 'amber14/protein.ff14SB.xml' or kwargs['environment'] not in ('vacuum', 'tip3p_pme') or kwargs['constraints'] not in ('none', 'HBonds'):
            raise ValueError('explicit ff14SB/environment/constraints policies required')
        if kwargs['environment'] == 'tip3p_pme':
            validate_solvation(kwargs.get('solvation'))
            expected = 0 if kwargs['solvation']['neutralize'] else kwargs['solvation']['solute_charge']
            if kwargs['expected_charge'] != expected:
                raise ValueError('expected final charge contradicts solvation policy')
        elif kwargs.get('solvation') is not None:
            raise ValueError('vacuum construction does not accept solvation settings')
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
            ligand_path = kwargs.get('ligand_parameters_file')
            allowed = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
            protein = [r for r in structure.first_model.residues if not (ligand_path and r.name == 'LIG')]
            if structure.source_format != 'pdb' or len(structure.models) != 1 or not protein or any(r.name not in allowed or r.altlocs or r.is_modified or r.is_hetero for r in protein):
                raise ValueError('parameterization currently requires a single canonical protein PDB without solvent, ligands, modifications or altlocs')
            if ligand_path:
                complex_ligand_mapping(structure, ligand_path)
            if not any(a.element == 'H' for r in protein for a in r.atoms):
                raise ValueError('protein must already contain explicit hydrogens; no protonation changes are made')
            if any(line.startswith('CRYST1') for line in Path(kwargs['structure_file']).read_text().splitlines()):
                raise ValueError('parameterization requires an unboxed prepared protein; existing periodic boxes are not silently replaced')
            if kwargs['environment'] == 'tip3p_pme':
                for axis, length in enumerate(kwargs['solvation']['box_size_nm']):
                    values = [atom.position[axis] / 10 for atom in structure.first_model.atoms]
                    if max(values) - min(values) >= length:
                        raise ValueError('requested periodic box is smaller than the source protein extent')
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('parameterization requires an empty step directory')
            shutil.copyfile(kwargs['structure_file'], root / 'input.pdb')
            options = {'constraints': kwargs['constraints'], 'environment': kwargs['environment'], 'solvation': kwargs.get('solvation')}
            if ligand_path:
                ligand_copy = copy_ligand_parameters(ligand_path, root / 'ligand_parameters')
                record = read_json(ligand_copy)['files']['parameters']
                options['ligand_forcefield'] = str(child_path(ligand_copy.parent, record['path']).relative_to(root.resolve()))
            runtime = run_worker(_WORKER, options, {'openmm': kwargs['expected_version']},
                                 python=kwargs.get('python', sys.executable), work_dir=root, timeout=kwargs.get('timeout', 600.0))
            periodic = kwargs['environment'] == 'tip3p_pme'
            names = dict(structure='input.pdb', system='system.xml', topology='topology.json', positions='positions.json', forcefield='forcefield.xml')
            if periodic:
                names.update(structure='solvated.pdb', source='input.pdb', solvent_forcefield='solvent_forcefield.xml')
            if ligand_path:
                names['ligand_parameters'] = str(ligand_copy.relative_to(root))
            files = {key: {'path': filename, 'sha256': file_hash(root / filename)} for key, filename in names.items()}
            prepared = read_structure(root / names['structure'])
            bundle = dict(schema_version=(2 if periodic else 1) + (2 if ligand_path else 0), engine='openmm', files=files, particle_count=len(prepared.first_model.atoms),
                          forcefield_name=kwargs['forcefield'], environment=kwargs['environment'], constraints=kwargs['constraints'],
                          nonbonded_method='PME' if periodic else 'NoCutoff', expected_charge=kwargs['expected_charge'], runtime=runtime,
                          provenance_origin='engine_run', scientific_validation='not_performed')
            if periodic:
                bundle['solvation'] = kwargs['solvation']
            output = root / 'simulation_system.json'
            write_json(output, bundle)
            validation = validate_system(output)
            report = dict(runtime=runtime, validation=validation, scope='ff14SB protein with optional mapped ligand parameters; explicit vacuum or TIP3P/PME policy')
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


class OpenMMComplexSystem(OpenMMSystem):
    name = 'bio_openmm_complex_parameterize'
    description = 'Build a checked ff14SB protein plus one explicitly parameterized LIG system, in vacuum or TIP3P/PME'
    consumes = ('complex_structure', 'ligand_parameters')
    wires = (('complex_structure', 'structure_file'), ('ligand_parameters', 'ligand_parameters_file'))
    params = tuple(ParamSpec('ligand_parameters_file', 'path', required=True) if p.name == 'ligand_parameters_file' else p for p in OpenMMSystem.params)
    example_kwargs = {**OpenMMSystem.example_kwargs, 'structure_file': 'complex.pdb', 'ligand_parameters_file': 'ligand_parameters.json'}


register(OpenMMComplexSystem())
