"""Ligand preparation: a SMILES or an SDF/MOL2 ligand becomes a 3D ligand
with explicit hydrogens and Gasteiger charges.

RDKit is already a DELFIN dependency (BSD-3 licence), so the layer may lean
on it where a format reader must not.  The embedding method is named
(``ETDGv3`` — the third generation of RDKit's distance-geometry embedding
with experimental torsion terms, so the coordinates are usable before any
force field), the seed is fixed, and the protonation state is taken as
given — no pKa model is applied, and this file says so.

A salt, a counter-ion, a co-crystallised water: an input that holds more
than one fragment is refused by name unless the caller asks for the largest
one to be kept, in which case the removed fragments are named in what the
function returns — a ligand file that silently sheds a component answers a
different question.

Reactions to a file already carrying 3D coordinates (an SDF/MOL2 ligand):
the coordinates are kept as they are — hydrogens are added onto them, and
the embedding step is skipped; embedding is only for a molecule without
coordinates (a SMILES).
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import List, Optional, Tuple, Union

from delfin.tools.bio._errors import StructureFormatError
from delfin.tools.bio._model import LigandAtom, SmallMolecule

# RDKit is imported lazily so a machine without it still imports this
# module; the adapters declare the requirement in their contract.
_RDKIT = None


def _rdkit():
    """Import RDKit once; ``Chem`` is imported alongside ``rdkit`` because
    ``import rdkit`` alone does not attach its submodules."""
    global _RDKIT
    if _RDKIT is None:
        try:
            import rdkit as module
            import rdkit.Chem  # noqa: F401
            import rdkit.Chem.AllChem  # noqa: F401 — embedding lives here
        except ImportError as exc:
            raise StructureFormatError(
                "ligand preparation needs RDKit, which is not importable "
                f"({exc}); install the 'rdkit' package (BSD-3 licence)"
            ) from None
        _RDKIT = module
    return _RDKIT


# The embedding seed.  Fixed so the same input gives the same coordinates on
# every run and every machine; ETKDG is a distance-geometry method whose
# random phase depends on it.
EMBED_SEED = 0xF00D
EMBED_METHOD = "ETDGv3"

_PathLike = Union[str, Path]


@dataclass(frozen=True)
class PreparedLigand:
    """The prepared ligand plus what the preparation decided."""

    molecule: SmallMolecule
    method: str
    seed: int
    removed_fragments: Tuple[str, ...] = ()

    @property
    def n_atoms(self) -> int:
        return self.molecule.n_atoms

    @property
    def has_hydrogens(self) -> bool:
        return self.molecule.has_hydrogens

    @property
    def charges(self) -> Tuple[float, ...]:
        return tuple(a.charge for a in self.molecule.atoms)


def prepare_ligand(
    *,
    smiles: Optional[str] = None,
    ligand_file: Optional[_PathLike] = None,
    keep_largest: bool = False,
    seed: int = EMBED_SEED,
) -> PreparedLigand:
    """Prepare one ligand from a SMILES string or a small-molecule file.

    Exactly one of ``smiles`` / ``ligand_file`` must be given.  A SMILES is
    embedded with :attr:`EMBED_METHOD` and ``seed``; a file's 3D coordinates
    are kept and only completed with hydrogens.  The result carries explicit
    hydrogens and Gasteiger charges on every atom.  Raises
    :class:`StructureFormatError` for a refused input (two or no sources, an
    unparseable SMILES, an unreadable file, more than one fragment without
    ``keep_largest``).
    """
    if (smiles is None) == (ligand_file is None):
        raise StructureFormatError(
            "prepare_ligand takes exactly one of 'smiles' or 'ligand_file'"
        )
    if smiles is not None:
        mol = _from_smiles(smiles)
        source_name = "<smiles>"
        embed = True
    else:
        from delfin.tools.bio._read import read_small_molecule
        smalls = read_small_molecule(ligand_file)
        if len(smalls) != 1:
            file_name = Path(ligand_file).name
            raise StructureFormatError(
                f"{file_name}: the file holds {len(smalls)} records, ligand "
                "preparation takes exactly one"
            )
        small = smalls[0]
        source_name = small.source_name
        mol = _from_small_molecule(small)
        embed = not mol.GetConformer().Is3D()

    for atom in mol.GetAtoms():
        atom.SetIntProp("_delfin_input_index", atom.GetIdx())
    if embed:
        mol = _embed(mol, seed=seed)
    mol = _add_hydrogens(mol, embed=embed, seed=seed)

    fragments = _fragments(mol)
    if len(fragments) > 1 and not keep_largest:
        sizes = ", ".join(str(len(f.GetAtoms())) for f in fragments)
        raise StructureFormatError(
            f"{source_name}: the molecule is {len(fragments)} fragments "
            f"({sizes} atoms); pass keep_largest=True to keep the largest "
            "one and remove the others"
        )
    removed: Tuple[str, ...] = ()
    if len(fragments) > 1:
        largest = max(range(len(fragments)),
                      key=lambda i: fragments[i].GetNumAtoms())
        dropped = [f for i, f in enumerate(fragments) if i != largest]
        removed = tuple(
            f"{_formula(d)} ({d.GetNumAtoms()} atoms)" for d in dropped
        )
        mol = fragments[largest]

    m = _to_small_molecule(mol, name=_source_label(smiles, source_name))
    return PreparedLigand(molecule=m, method=EMBED_METHOD if embed else
                          "coordinates-kept", seed=seed,
                          removed_fragments=removed)


def _from_smiles(smiles: str):
    """Parse one SMILES string; an empty or unparseable one is refused."""
    rdkit = _rdkit()
    text = str(smiles).strip()
    if not text:
        raise StructureFormatError("the SMILES string is empty")
    mol = rdkit.Chem.MolFromSmiles(text)
    if mol is None:
        raise StructureFormatError(f"the SMILES {text!r} could not be parsed")
    return mol


def _from_small_molecule(small: SmallMolecule):
    """Build an RDKit molecule, bonds included, from a read record.

    Bond order 4 in an SDF is aromatic; the reader of
    :mod:`delfin.tools.bio._ligand` passes it through as it stands.
    """
    rdkit = _rdkit()
    if small.molblock:
        mol = rdkit.Chem.MolFromMolBlock(small.molblock, removeHs=False)
        if mol is None or mol.GetNumAtoms() != small.n_atoms:
            raise StructureFormatError(f"{small.source_name}: invalid molecular graph")
        conformer = mol.GetConformer()
        conformer.Set3D(True)
        for k, atom in enumerate(small.atoms):
            conformer.SetAtomPosition(k, atom.position)
        return mol
    rw = rdkit.Chem.RWMol()
    for atom in small.atoms:
        try:
            rw.AddAtom(rdkit.Chem.Atom(atom.element))
        except (AttributeError, RuntimeError, ValueError) as exc:
            raise StructureFormatError(
                f"{small.source_name}: atom {atom.name!r} "
                f"({atom.element}) cannot become an RDKit atom ({exc})"
            ) from None
    for i, j, order in small.bonds:
        if not (0 <= i < small.n_atoms and 0 <= j < small.n_atoms):
            raise StructureFormatError(
                f"{small.source_name}: bond ({i}, {j}) names an atom outside "
                f"the record's {small.n_atoms} atoms"
            )
        bond_type = {1: rdkit.Chem.BondType.SINGLE,
                     2: rdkit.Chem.BondType.DOUBLE,
                     3: rdkit.Chem.BondType.TRIPLE,
                     4: rdkit.Chem.BondType.AROMATIC}.get(int(order))
        if bond_type is None:
            raise StructureFormatError(
                f"{small.source_name}: bond order {order} between atoms "
                f"{i} and {j} is not 1, 2, 3 or 4"
            )
        rw.AddBond(int(i), int(j), bond_type)
    mol = rw.GetMol()
    conformer = rdkit.Chem.Conformer(mol.GetNumAtoms())
    conformer.Set3D(True)
    for k, atom in enumerate(small.atoms):
        conformer.SetAtomPosition(k, (atom.x, atom.y, atom.z))
    mol.AddConformer(conformer, assignId=True)
    rdkit.Chem.SanitizeMol(mol)
    return mol


def _embed(mol, *, seed: int):
    """Embed a molecule without coordinates using ETKDGv3 with *seed*."""
    rdkit = _rdkit()
    params = rdkit.Chem.AllChem.ETKDGv3()
    params.randomSeed = seed
    if rdkit.Chem.AllChem.EmbedMolecule(mol, params) != 0:
        raise StructureFormatError(
            "the molecule could not be embedded in 3D with "
            f"{EMBED_METHOD} (seed {seed}); it may be too flexible or its "
            "SMILES may be chemically impossible"
        )
    return mol


def _add_hydrogens(mol, *, embed: bool, seed: int):
    """Add explicit hydrogens; a freshly embedded molecule is re-embedded
    with them, so the coordinates include the hydrogen positions."""
    rdkit = _rdkit()
    work = rdkit.Chem.AddHs(mol, addCoords=not embed)
    if embed:
        params = rdkit.Chem.AllChem.ETKDGv3()
        params.randomSeed = seed
        if rdkit.Chem.AllChem.EmbedMolecule(work, params) != 0:
            raise StructureFormatError(
                "the molecule could not be embedded with its hydrogens "
                f"using {EMBED_METHOD} (seed {seed})"
            )
    return work


def _fragments(mol) -> List["object"]:
    """The connected components of *mol* as separate RDKit molecules."""
    rdkit = _rdkit()
    parts = rdkit.Chem.rdmolops.GetMolFrags(mol, asMols=True,
                                            sanitizeFrags=True)
    return list(parts)


def _formula(mol) -> str:
    counts: dict = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        counts[sym] = counts.get(sym, 0) + 1
    return "".join(
        f"{sym}{n if n != 1 else ''}"
        for sym, n in sorted(counts.items(), key=lambda kv: (kv[0] != "C",
                                                             kv[0]))
    )


def _source_label(smiles, source_name: str) -> str:
    if smiles is not None:
        return str(smiles).strip()[:60]
    return source_name


def _to_small_molecule(mol, *, name: str) -> SmallMolecule:
    """An RDKit molecule back into a :class:`SmallMolecule`.

    Gasteiger partial charges are computed here and land on every atom as
    its ``charge``; the same numbers go into the ``GASTEIGER_CHARGES`` SD
    property so a file reader sees them too.  RDKit writes NaN for an atom
    it cannot parametrise; a NaN is not a number a consumer can use, so
    such an atom's charge is refused rather than written.
    """
    rdkit = _rdkit()
    rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(mol)
    charges: List[float] = []
    for atom in mol.GetAtoms():
        value = float(atom.GetDoubleProp("_GasteigerCharge"))
        if not math.isfinite(value):
            raise StructureFormatError(
                f"atom {atom.GetIdx()} ({atom.GetSymbol()}) has no Gasteiger "
                "charge; the element may be unsupported"
            )
        charges.append(value)
    conformer = mol.GetConformer()
    atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conformer.GetAtomPosition(i)
        atoms.append(LigandAtom(
            name=f"{atom.GetSymbol()}{i + 1}",
            element=atom.GetSymbol(),
            x=float(pos.x), y=float(pos.y), z=float(pos.z),
            charge=charges[i],
        ))
    bonds = []
    for bond in mol.GetBonds():
        if bond.GetIsAromatic():
            order = 4
        else:
            order = {rdkit.Chem.BondType.SINGLE: 1,
                     rdkit.Chem.BondType.DOUBLE: 2,
                     rdkit.Chem.BondType.TRIPLE: 3}.get(bond.GetBondType(), 1)
        bonds.append((int(bond.GetBeginAtomIdx()),
                      int(bond.GetEndAtomIdx()), order))
    canonical = rdkit.Chem.RemoveHs(rdkit.Chem.Mol(mol))
    for atom in canonical.GetAtoms():
        atom.SetAtomMapNum(0)
    smiles = rdkit.Chem.MolToSmiles(canonical, isomericSmiles=True)
    identity = {
        "schema_version": 1,
        "canonical_isomeric_smiles": smiles,
        "chemical_state_id": hashlib.sha256(smiles.encode()).hexdigest(),
        "identity_method": "RDKit canonical isomeric SMILES; no tautomer/protonation normalization",
        "rdkit_version": rdkit.__version__,
        "formal_charge": rdkit.Chem.GetFormalCharge(mol),
        "atom_mapping": [{"output_index": atom.GetIdx(),
                          "input_index": atom.GetIntProp("_delfin_input_index") if atom.HasProp("_delfin_input_index") else None,
                          "element": atom.GetSymbol(),
                          "origin": "input" if atom.HasProp("_delfin_input_index") else "added_hydrogen"}
                         for atom in mol.GetAtoms()],
    }
    return SmallMolecule(
        name=name,
        atoms=tuple(atoms),
        bonds=tuple(bonds),
        properties={"DELFIN_IDENTITY": json.dumps(identity, sort_keys=True), "GASTEIGER_CHARGES":
                    " ".join(f"{c:.6f}" for c in charges)},
        source_format="prepared",
        source_name=name,
        molblock=rdkit.Chem.MolToMolBlock(mol),
    )


def write_sdf(ligand: PreparedLigand, path: _PathLike) -> None:
    """Write the prepared ligand as one SDF (V2000) record.

    Formal charges and stereochemistry travel in the molecular graph;
    Gasteiger partial charges remain a separate SD property.
    """
    molecule = ligand.molecule
    if molecule.molblock:
        mol = _from_small_molecule(molecule)
        mol.SetProp("_Name", molecule.name)
        for key, value in molecule.properties.items():
            mol.SetProp(key, value)
        with _rdkit().Chem.SDWriter(str(path)) as writer:
            writer.write(mol)
        return
    lines: List[str] = [molecule.name or "ligand",
                        "  DELFIN  prepared ligand", ""]
    n_atoms = molecule.n_atoms
    n_bonds = len(molecule.bonds)
    lines.append(f"{n_atoms:3d}{n_bonds:3d}  0  0  0  0  0  0  0  0999 V2000")
    for atom in molecule.atoms:
        lines.append(
            f"{atom.x:10.4f}{atom.y:10.4f}{atom.z:10.4f} "
            f"{atom.element:<3}  0  0  0  0  0  0  0  0  0  0  0  0"
        )
    for i, j, order in molecule.bonds:
        lines.append(f"{i + 1:3d}{j + 1:3d}{order:3d}  0  0  0  0")
    lines.append("M  END")
    lines.append(">  <GASTEIGER_CHARGES>")
    lines.append(molecule.properties.get("GASTEIGER_CHARGES", ""))
    lines.append("")
    lines.append("$$$$")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


__all__ = ["PreparedLigand", "prepare_ligand", "write_sdf",
           "EMBED_SEED", "EMBED_METHOD"]
