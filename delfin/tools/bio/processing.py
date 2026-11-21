"""RDKit-backed processing primitives that preserve molecular chemistry."""
from __future__ import annotations

import gzip
import hashlib
import io
import math
from pathlib import Path


def load_ligands(path: Path) -> list:
    from rdkit import Chem
    path = Path(path)
    raw = gzip.decompress(path.read_bytes()) if path.suffix.lower() == ".gz" else path.read_bytes()
    suffix = Path(path.stem).suffix.lower() if path.suffix.lower() == ".gz" else path.suffix.lower()
    if suffix == ".sdf":
        records = list(Chem.ForwardSDMolSupplier(io.BytesIO(raw), removeHs=False))
    elif suffix == ".mol":
        records = [Chem.MolFromMolBlock(raw.decode(), removeHs=False)]
    elif suffix == ".mol2":
        if raw.count(b"@<TRIPOS>MOLECULE") != 1:
            raise ValueError("MOL2 processing currently requires exactly one record")
        records = [Chem.MolFromMol2Block(raw.decode(), removeHs=False)]
    elif suffix in {".smi", ".smiles"}:
        records = []
        for line in raw.decode().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split(maxsplit=1)
            mol = Chem.MolFromSmiles(fields[0])
            if mol is not None and len(fields) > 1:
                mol.SetProp("_Name", fields[1])
            records.append(mol)
    else:
        raise ValueError("supported ligand formats: SDF, MOL, MOL2, SMI/SMILES (optionally gzip)")
    if not records:
        raise ValueError("ligand file contains no records")
    for index, molecule in enumerate(records):
        if molecule is None or molecule.GetNumAtoms() == 0:
            raise ValueError(f"invalid molecule record {index + 1}; no records silently skipped")
        for conf in molecule.GetConformers():
            if not all(math.isfinite(v) for row in conf.GetPositions() for v in row):
                raise ValueError(f"non-finite coordinates in record {index + 1}")
    return records


def chemical_identity(mol) -> dict:
    from rdkit import Chem, __version__
    copy = Chem.RemoveHs(Chem.Mol(mol))
    for atom in copy.GetAtoms():
        atom.SetAtomMapNum(0)
    smiles = Chem.MolToSmiles(copy, isomericSmiles=True)
    return {"chemical_state_id": hashlib.sha256(smiles.encode()).hexdigest(),
            "canonical_isomeric_smiles": smiles, "formal_charge": Chem.GetFormalCharge(mol),
            "rdkit_version": __version__}


def ligand_properties(path: Path) -> dict:
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from rdkit import Chem
    rows = []
    for index, mol in enumerate(load_ligands(path)):
        rows.append({"record_id": index + 1, "name": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
                     **chemical_identity(mol), "molecular_weight": Descriptors.MolWt(mol),
                     "molecular_weight_unit": "g/mol", "clogp": Descriptors.MolLogP(mol),
                     "tpsa_angstrom2": rdMolDescriptors.CalcTPSA(mol), "hbd": Lipinski.NumHDonors(mol),
                     "hba": Lipinski.NumHAcceptors(mol), "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
                     "heavy_atoms": mol.GetNumHeavyAtoms(),
                     "scaffold_smiles": Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))})
    return {"schema_version": 1, "kind": "computed_molecular_descriptors", "records": rows,
            "note": "Calculated descriptors, not measured ADMET or biological activity"}


def enumerate_states(path: Path, output: Path, *, max_states: int = 32) -> dict:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    if type(max_states) is not int or not 1 <= max_states <= 1024:
        raise ValueError("max_states must be an integer in 1..1024")
    sources = load_ligands(path)
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(max_states)
    generated, report = [], []
    for record, mol in enumerate(sources, 1):
        parent = chemical_identity(mol)
        source = Chem.RemoveHs(mol)
        source.RemoveAllConformers()
        for atom in source.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        states = enumerator.Enumerate(source)
        if str(states.status) != "Completed":
            raise ValueError(f"record {record}: tautomer enumeration incomplete ({states.status}); increase explicit limits")
        seen = set()
        for mol_state in states:
            identity = chemical_identity(mol_state)
            if identity["chemical_state_id"] in seen:
                continue
            seen.add(identity["chemical_state_id"])
            mol_state.SetProp("DELFIN_PARENT_STATE_ID", parent["chemical_state_id"])
            mol_state.SetProp("DELFIN_STATE_ID", identity["chemical_state_id"])
            mol_state.SetProp("DELFIN_SOURCE_RECORD", str(record))
            generated.append(mol_state)
            report.append({"source_record": record, "output_record": len(generated), **identity})
    with Chem.SDWriter(str(output)) as writer:
        for mol in generated:
            writer.write(mol)
    return {"schema_version": 1, "method": "RDKit tautomer enumeration", "states": report,
            "protonation": "no pH model", "coordinates": "must regenerate 3D coordinates per state",
            "atom_mapping": "atom-map numbers refer to heavy atoms of each input record"}


def ligand_grid_svg(path: Path, *, max_molecules: int = 64) -> str:
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor
    molecules = load_ligands(path)
    if len(molecules) > max_molecules:
        raise ValueError(f"{len(molecules)} molecules exceed explicit depiction limit {max_molecules}")
    copies = []
    for mol in molecules:
        copy = Chem.RemoveHs(Chem.Mol(mol))
        rdDepictor.Compute2DCoords(copy)
        copies.append(copy)
    return Draw.MolsToGridImage(copies, molsPerRow=4, subImgSize=(260, 220),
                               legends=[f"Record {i + 1}" for i in range(len(copies))], useSVG=True)


def select_ligand(path: Path, output: Path, *, record_id: int) -> dict:
    from rdkit import Chem
    from delfin.tools.bio.artifacts import file_hash
    molecules = load_ligands(path)
    if type(record_id) is not int or not 1 <= record_id <= len(molecules):
        raise ValueError(f'record_id must be an explicit 1-based record in 1..{len(molecules)}')
    selected = molecules[record_id - 1]
    identity = chemical_identity(selected)
    with Chem.SDWriter(str(output)) as writer:
        writer.write(selected)
    if chemical_identity(load_ligands(output)[0])['chemical_state_id'] != identity['chemical_state_id']:
        raise ValueError('SDF export changed the selected chemical state')
    return {'schema_version': 1, 'source_sha256': file_hash(path), 'source_record': record_id,
            'source_count': len(molecules), 'output_record': 1, **identity,
            'selection': 'explicit record, not automatic score ranking',
            'coordinates': 'retained when present; no embedding or alignment'}


def box_from_ligand(path: Path, *, padding: float = 4.0) -> dict:
    from delfin.tools.bio import BindingBox
    if isinstance(padding, bool) or not isinstance(padding, (float, int)) or not math.isfinite(padding) or padding <= 0:
        raise ValueError('padding must be positive and finite, in angstrom')
    molecules = load_ligands(path)
    if len(molecules) != 1 or not molecules[0].GetNumConformers() or not molecules[0].GetConformer().Is3D():
        raise ValueError('box requires exactly one ligand with existing 3D coordinates')
    mol = molecules[0]
    conf = mol.GetConformer()
    coords = [conf.GetAtomPosition(atom.GetIdx()) for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if not coords:
        raise ValueError('box requires heavy atoms')
    low = [min(point[i] for point in coords) for i in range(3)]
    high = [max(point[i] for point in coords) for i in range(3)]
    return BindingBox(center=tuple((a + b) / 2 for a, b in zip(low, high)),
                      size=tuple(b - a + 2 * padding for a, b in zip(low, high))).to_dict()
