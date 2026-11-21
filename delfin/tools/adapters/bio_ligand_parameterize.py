"""Optional OpenFF/SMIRNOFF parameter generation with supplied, mapped charges."""
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
from delfin.tools.bio.artifacts import file_hash, write_json
from delfin.tools.bio.ligand_parameters import ligand_input, charge_input, validate_offxml, validate_ligand_parameters
from delfin.tools.bio.processing import chemical_identity
from delfin.tools.bio.optional_runtime import probe_distribution, run_worker


_WORKER = '''import json,sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from openff.toolkit import Molecule
from openff.toolkit.utils import RDKitToolkitWrapper, ToolkitRegistry, toolkit_registry_manager
from openff.units import unit
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
options=json.load(open(sys.argv[1]))
source=Chem.SDMolSupplier('ligand.sdf',removeHs=False)[0]
def refuse(*args,**kwargs):
    raise ValueError('implicit charge, conformer or fractional-bond-order calculation is forbidden')
Molecule.assign_partial_charges=refuse
Molecule.generate_conformers=refuse
Molecule.assign_fractional_bond_orders=refuse
with toolkit_registry_manager(ToolkitRegistry([RDKitToolkitWrapper()])):
    molecule=Molecule.from_rdkit(source,hydrogens_are_explicit=True,allow_undefined_stereo=False)
    restored=molecule.to_rdkit()
    if Chem.MolToSmiles(source,isomericSmiles=True)!=Chem.MolToSmiles(restored,isomericSmiles=True):
        raise ValueError('OpenFF conversion changed ligand chemistry')
    if [(a.GetAtomicNum(),a.GetFormalCharge()) for a in source.GetAtoms()]!=[(a.GetAtomicNum(),a.GetFormalCharge()) for a in restored.GetAtoms()]:
        raise ValueError('OpenFF conversion changed atom ordering')
    def bonds(m):
        return sorted((min(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),max(b.GetBeginAtomIdx(),b.GetEndAtomIdx()),b.GetBondTypeAsDouble()) for b in m.GetBonds())
    if bonds(source)!=bonds(restored):
        raise ValueError('OpenFF conversion changed mapped bonds')
    for i,atom in enumerate(molecule.atoms):
        atom.name=f'L{i+1:03d}'
    molecule.partial_charges=np.asarray(json.load(open('charges.json'))['values'])*unit.elementary_charge
    generator=SMIRNOFFTemplateGenerator(molecules=[molecule],forcefield=str(Path('forcefield.offxml').resolve()))
    Path('parameters.xml').write_text(generator.generate_residue_template(molecule))
'''


class LigandParameterize(StepAdapter):
    name = 'bio_ligand_parameterize'
    category = 'structure'
    description = 'Generate an OpenFF ligand parameter bundle from explicit-H SDF, supplied mapped charges and a local SMIRNOFF file'
    produces_geometry = False
    requires_python = ('rdkit',)
    consumes = ('ligand',)
    wires = (('ligand', 'ligand_file'),)
    produces = ('ligand_parameters', 'json_report')
    required_artifacts = produces
    params = (ParamSpec('ligand_file', 'path', required=True), ParamSpec('charges_file', 'path', required=True),
              ParamSpec('forcefield_file', 'path', required=True), ParamSpec('expected_version', 'str', required=True),
              ParamSpec('toolkit_version', 'str', required=True), ParamSpec('openmm_version', 'str', required=True),
              ParamSpec('python', 'path', default=sys.executable), ParamSpec('timeout', 'float', default=600.0))
    data_keys = (DataKeySpec('runtime', 'dict'), DataKeySpec('validation', 'dict'))
    example_kwargs = dict(ligand_file='ligand.sdf', charges_file='charges.json', forcefield_file='forcefield.offxml',
                          expected_version='0.15.0', toolkit_version='0.17.0', openmm_version='8.4.0')

    def validate_params(self, **kwargs):
        for parameter in self.params:
            if parameter.required and not kwargs.get(parameter.name):
                raise ValueError(parameter.name + ' is required')
        timeout = kwargs.get('timeout', 600.)
        if type(timeout) not in (int, float) or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError('timeout must be positive and finite')

    def runtime_preflight(self, **kwargs):
        first = probe_distribution('openmmforcefields', python=kwargs.get('python', sys.executable),
                                   expected_version=kwargs['expected_version'], work_dir=Path.cwd(),
                                   allow_install=kwargs.get('_allow_install', False))
        result = {'openmmforcefields': first}
        for name, key in (('openff-toolkit', 'toolkit_version'), ('openmm', 'openmm_version')):
            result[name] = probe_distribution(name, python=first['python'], expected_version=kwargs[key],
                                               work_dir=Path.cwd(), discover=False, environment_module=first.get('environment_module'))
        return result

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        start, root = time.monotonic(), Path(work_dir)
        try:
            self.validate_params(**kwargs)
            mol = ligand_input(kwargs['ligand_file'])
            charge_input(kwargs['charges_file'], kwargs['ligand_file'], mol)
            validate_offxml(kwargs['forcefield_file'])
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('ligand parameterization requires an empty step directory')
            files = {'ligand': 'ligand.sdf', 'charges': 'charges.json', 'offxml': 'forcefield.offxml', 'parameters': 'parameters.xml'}
            for parameter, name in [('ligand_file', 'ligand'), ('charges_file', 'charges'), ('forcefield_file', 'offxml')]:
                shutil.copyfile(kwargs[parameter], root / files[name])
            versions = {'openmmforcefields': kwargs['expected_version'], 'openff-toolkit': kwargs['toolkit_version'], 'openmm': kwargs['openmm_version']}
            runtime = run_worker(_WORKER, {}, versions, python=kwargs.get('python', sys.executable),
                                 work_dir=root, timeout=kwargs.get('timeout', 600.))
            bundle = dict(schema_version=1, method='smirnoff_supplied_charges', atom_count=mol.GetNumAtoms(),
                          chemical_state_id=chemical_identity(mol)['chemical_state_id'], runtime=runtime,
                          files={k: {'path': v, 'sha256': file_hash(root / v)} for k, v in files.items()},
                          scientific_validation='not_performed')
            output = root / 'ligand_parameters.json'
            write_json(output, bundle)
            report = {'runtime': runtime, 'validation': validate_ligand_parameters(output),
                      'charge_policy': 'caller-supplied mapped values; no implicit charge or conformer calculation'}
            write_json(root / 'report.json', report)
            return self._make_result(self.name, StepStatus.SUCCESS, root, start, data=report,
                                     artifacts={'ligand_parameters': output, 'json_report': root / 'report.json'})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, TypeError, KeyError, OSError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, root, start, error=error, error_kind=kind)


register(LigandParameterize())
