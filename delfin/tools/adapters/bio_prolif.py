"""Explicit ProLIF interaction fingerprints with indexed molecule sidecars."""
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
from delfin.tools._types import StepStatus, ErrorKind
from delfin.tools.bio import read_structure
from delfin.tools.bio.artifacts import file_hash, read_json, validate_artifact, write_json
from delfin.tools.bio.optional_runtime import run_worker, probe_distribution
from delfin.tools.bio.processing import load_ligands, chemical_identity


def _validate_mapped_pose(source, exported, mapping):
    """Validate chemistry and geometry against the reported atom permutation."""
    indices = list(range(source.GetNumAtoms()))
    if (not isinstance(mapping, list) or exported.GetNumAtoms() != len(indices)
            or any(not isinstance(m, dict) or type(m.get('parent_index')) is not int
                   or type(m.get('source_atom_index')) is not int for m in mapping)
            or sorted(m['parent_index'] for m in mapping) != indices
            or sorted(m['source_atom_index'] for m in mapping) != indices):
        raise ValueError('ProLIF ligand atom mapping is incomplete')
    if chemical_identity(source)['chemical_state_id'] != chemical_identity(exported)['chemical_state_id']:
        raise ValueError('ProLIF export changed ligand chemistry')
    lookup = {m['parent_index']: m['source_atom_index'] for m in mapping}
    for parent, original in lookup.items():
        a, b = source.GetAtomWithIdx(original), exported.GetAtomWithIdx(parent)
        if (a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge()) != (b.GetAtomicNum(), b.GetIsotope(), b.GetFormalCharge()):
            raise ValueError('ProLIF atom mapping changed atom identity')
        p, q = source.GetConformer().GetAtomPosition(original), exported.GetConformer().GetAtomPosition(parent)
        if any(abs(p[i] - q[i]) > 0.0002 for i in range(3)):
            raise ValueError('ProLIF export changed pose coordinates')
    def bonds(mol, remap):
        return {(tuple(sorted((remap[b.GetBeginAtomIdx()], remap[b.GetEndAtomIdx()]))),
                 str(b.GetBondType()), b.GetIsAromatic()) for b in mol.GetBonds()}
    if bonds(source, dict(zip(indices, indices))) != bonds(exported, lookup):
        raise ValueError('ProLIF atom mapping changed bond connectivity')


def _validate_contacts(contacts, requested, ligand_atoms, receptor_atoms):
    if not isinstance(contacts, list):
        raise ValueError('invalid interaction rows')
    for contact in contacts:
        if not isinstance(contact, dict) or contact.get('interaction') not in requested:
            raise ValueError('unexpected interaction type in output')
        metadata = contact.get('metadata')
        parent = metadata.get('parent_indices') if isinstance(metadata, dict) else None
        if not isinstance(parent, dict):
            raise ValueError('interaction metadata lacks parent atom indices')
        for name, count in (('ligand', ligand_atoms), ('protein', receptor_atoms)):
            values = parent.get(name)
            if not isinstance(values, list) or not values or any(type(i) is not int or not 0 <= i < count for i in values):
                raise ValueError('interaction parent atom index is invalid')


_WORKER = '''import json,sys
from pathlib import Path
from rdkit import Chem
import prolif as plf
options=json.load(open(sys.argv[1]))
protein=Chem.MolFromPDBFile('receptor.pdb',removeHs=False,sanitize=True,proximityBonding=True)
if protein is None:
    raise ValueError('receptor could not be converted to an RDKit molecular graph')
if any(a.GetTotalNumHs(includeNeighbors=False) for a in protein.GetAtoms()):
    raise ValueError('receptor graph has implicit/unexpanded hydrogens')
prot=plf.Molecule(protein)
Chem.MolToMolFile(prot,'prolif_receptor.mol',forceV3000=True)
fp=plf.Fingerprint(interactions=options['interactions'],count=True)
rows=[]
for record,mol in enumerate(Chem.SDMolSupplier('poses.sdf',removeHs=False),1):
    if mol is None:
        raise ValueError('invalid pose; none may be silently skipped')
    if any(a.GetTotalNumHs(includeNeighbors=False) for a in mol.GetAtoms()):
        raise ValueError('pose graph has implicit/unexpanded hydrogens')
    for a in mol.GetAtoms():
        a.SetIntProp('DELFIN_SOURCE_INDEX',a.GetIdx())
    lig=plf.Molecule(mol)
    name=f'prolif_pose_{record}.mol'
    Chem.MolToMolFile(lig,name,forceV3000=True)
    mapping=[{'parent_index':a.GetIdx(),'source_atom_index':a.GetIntProp('DELFIN_SOURCE_INDEX')} for a in lig.GetAtoms()]
    contacts=[]
    for (ligand_residue,protein_residue),types in fp.generate(lig,prot,metadata=True).items():
        for interaction,occurrences in types.items():
            for metadata in occurrences:
                contacts.append({'ligand_residue':str(ligand_residue),'protein_residue':str(protein_residue),
                                 'interaction':interaction,'metadata':metadata})
    rows.append({'record_id':record,'molecule_file':name,'atom_mapping':mapping,'interactions':contacts})
def scalar(value):
    if hasattr(value,'tolist'): return value.tolist()
    if hasattr(value,'item'): return value.item()
    raise TypeError(type(value).__name__)
with open('interactions.json','w') as stream:
    json.dump({'schema_version':1,'records':rows,'receptor_molecule':'prolif_receptor.mol',
               'index_convention':'parent_indices refer to the exported ProLIF molecule sidecars'},stream,default=scalar,allow_nan=False)
'''


class Prolif(StepAdapter):
    name = 'prolif_analyze'
    description = 'Compute explicit ProLIF interactions for every prepared pose, retaining metadata and indexed molecule sidecars'
    category = 'structure'
    produces_geometry = False
    consumes = ('pose_set', 'prepared_receptor')
    wires = (('pose_set', 'poses_file'), ('prepared_receptor', 'receptor_file'))
    produces = ('interaction_report', 'json_report')
    required_artifacts = produces
    requires_python = ('rdkit',)
    params = (ParamSpec('poses_file', 'path', required=True), ParamSpec('receptor_file', 'path', required=True),
              ParamSpec('interactions', 'list', required=True), ParamSpec('expected_version', 'str', required=True),
              ParamSpec('bond_policy', 'str', required=True, enum=('rdkit_proximity',)),
              ParamSpec('shared_coordinate_frame', 'bool', required=True),
              ParamSpec('python', 'path', default=sys.executable), ParamSpec('timeout', 'float', default=600.0))
    data_keys = (DataKeySpec('runtime', 'dict'), DataKeySpec('records', 'list'))
    example_kwargs = dict(poses_file='poses.sdf', receptor_file='receptor.pdb', interactions=['Hydrophobic', 'HBDonor', 'HBAcceptor'],
                          expected_version='2.2.3', bond_policy='rdkit_proximity', shared_coordinate_frame=True)

    def validate_params(self, **kwargs):
        for param in self.params:
            if param.required and kwargs.get(param.name) is None:
                raise ValueError(f'{param.name} is required')
        values = kwargs['interactions']
        if not isinstance(kwargs['expected_version'], str) or not kwargs['expected_version'].strip():
            raise ValueError('expected_version must be nonempty')
        if not isinstance(values, list) or not values or any(not isinstance(x, str) or not x.isidentifier() for x in values) or len(values) != len(set(values)):
            raise ValueError('interactions must be a nonempty unique list of named ProLIF interactions')
        if kwargs['bond_policy'] != 'rdkit_proximity' or kwargs['shared_coordinate_frame'] is not True:
            raise ValueError('explicit rdkit_proximity bond policy and shared coordinate frame are required')
        timeout = kwargs.get('timeout', 600.0)
        if isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or not math.isfinite(timeout) or timeout <= 0:
            raise ValueError('timeout must be positive and finite')

    def runtime_preflight(self, **kwargs):
        return probe_distribution('prolif', python=kwargs.get('python', sys.executable), expected_version=kwargs['expected_version'],
                                  work_dir=Path.cwd(), allow_install=kwargs.get('_allow_install', False))

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        from rdkit import Chem
        from delfin.tools.bio.artifacts import child_path
        start = time.monotonic()
        root = Path(work_dir)
        try:
            self.validate_params(**kwargs)
            validate_artifact(Path(kwargs['poses_file']), 'pose_set')
            validate_artifact(Path(kwargs['receptor_file']), 'prepared_receptor')
            structure = read_structure(kwargs['receptor_file'])
            standard = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())
            if structure.source_format != 'pdb' or len(structure.models) != 1 or any(r.name not in standard or r.altlocs or r.is_modified for r in structure.first_model.residues):
                raise ValueError('ProLIF RDKit subset requires one canonical protein PDB without waters, altlocs or cofactors')
            if not any(a.element == 'H' for a in structure.first_model.atoms):
                raise ValueError('ProLIF requires explicit receptor hydrogens')
            poses = load_ligands(kwargs['poses_file'])
            if any(a.GetTotalNumHs(includeNeighbors=False) for mol in poses for a in mol.GetAtoms()):
                raise ValueError('ProLIF requires explicit ligand hydrogens')
            root.mkdir(parents=True, exist_ok=True)
            if any(root.iterdir()):
                raise ValueError('analysis requires an empty step directory')
            shutil.copyfile(kwargs['receptor_file'], root / 'receptor.pdb')
            with Chem.SDWriter(str(root / 'poses.sdf')) as writer:
                for mol in poses:
                    writer.write(mol)
            runtime = run_worker(_WORKER, {'interactions': kwargs['interactions']}, {'prolif': kwargs['expected_version']},
                                 python=kwargs.get('python', sys.executable), work_dir=root, timeout=kwargs.get('timeout', 600.0))
            data = read_json(root / 'interactions.json')
            rows = data.get('records')
            if not isinstance(rows, list) or any(not isinstance(r, dict) for r in rows) or [r.get('record_id') for r in rows] != list(range(1, len(poses) + 1)):
                raise ValueError('interaction report must contain exactly one ordered record per pose')
            validate_artifact(child_path(root, data['receptor_molecule']), 'ligand')
            protein = load_ligands(child_path(root, data['receptor_molecule']))[0]
            original_atoms = sorted((a.element.upper(), *(round(v, 3) for v in a.position)) for a in structure.first_model.atoms)
            exported_atoms = sorted((a.GetSymbol().upper(), *(round(v, 3) for v in protein.GetConformer().GetAtomPosition(a.GetIdx()))) for a in protein.GetAtoms())
            if original_atoms != exported_atoms:
                raise ValueError('ProLIF export changed receptor atoms or coordinates')
            for row, pose in zip(rows, poses):
                validate_artifact(child_path(root, row['molecule_file']), 'ligand')
                exported = load_ligands(child_path(root, row['molecule_file']))[0]
                _validate_mapped_pose(pose, exported, row['atom_mapping'])
                _validate_contacts(row.get('interactions'), kwargs['interactions'], exported.GetNumAtoms(), protein.GetNumAtoms())
            data.update(runtime=runtime, requested_interactions=kwargs['interactions'], bond_policy=kwargs['bond_policy'],
                        receptor_sha256=file_hash(Path(kwargs['receptor_file'])), poses_sha256=file_hash(Path(kwargs['poses_file'])),
                        interpretation='Interaction geometry, not affinity; receptor bond orders inferred by RDKit, no alignment performed')
            write_json(root / 'interactions.json', data)
            return self._make_result(self.name, StepStatus.SUCCESS, root, start, data=data,
                                     artifacts={'interaction_report': root / 'interactions.json', 'json_report': root / 'interactions.json'})
        except subprocess.TimeoutExpired as exc:
            kind, error = ErrorKind.TIMEOUT, str(exc)
        except ImportError as exc:
            kind, error = ErrorKind.BINARY_NOT_FOUND, str(exc)
        except (ValueError, TypeError, KeyError, OSError) as exc:
            kind, error = ErrorKind.UNSUITABLE_INPUT, str(exc)
        except RuntimeError as exc:
            kind, error = ErrorKind.TOOL_FAILED, str(exc)
        return self._make_result(self.name, StepStatus.FAILED, root, start, error=error, error_kind=kind)


register(Prolif())
