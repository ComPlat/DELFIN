"""The shared structure layer for biomolecular systems.

Proteins and nucleic acids with their bound small molecules need the same
thing the quantum-chemistry blocks of :mod:`delfin.tools` already have: one
place where the objects they hand each other are defined, so a block that
prepares a receptor, one that places a ligand and one that reports a
composition all mean the same thing by "chain", "residue" and "binding site".

What is here
------------
``read_structure`` / ``read_pdb`` / ``read_mmcif``
    A macromolecular structure from PDB or mmCIF → a
    :class:`MacromolecularStructure` with models, chains, residues (insertion
    codes and alternate locations kept), hetero components, waters, the
    missing residues the file itself reports, and its connection records.
``read_small_molecule`` / ``read_sdf`` / ``read_mol2``
    A ligand record from SDF/MOL or MOL2 → :class:`SmallMolecule`.
``BindingBox`` / ``ResidueSelection``
    The two typed ways of naming a binding site, both JSON round-trippable.

What is deliberately **not** here: anything that runs.  No docking, no
preparation, no structure prediction, no download and no force field.  This
module is the vocabulary those blocks share; each of them stays its own
building block with its own contract.

Dependencies
------------
None.  The readers are plain Python over the file formats, so the layer
works in any environment DELFIN itself runs in, and a block built on it is
never unavailable because an optional package is missing.  A heavier
structure library may be added later behind a feature check as an accelerator
for large assemblies — it must never become the only way to read a file,
because the formats are the contract, not the library.

Everything is immutable: a structure is passed between blocks and must not
be edited in place by one of them.  A block that changes a system produces a
new file and a new structure, which is also what makes the change reviewable.
"""

from delfin.tools.bio._errors import (StructureFormatError,
                                      StructureTooLargeError)
from delfin.tools.bio._ligand import parse_mol2, parse_sdf
from delfin.tools.bio._mmcif import parse_mmcif
from delfin.tools.bio._model import (
    WATER_NAMES,
    AtomRef,
    BioAtom,
    BondedLigandLink,
    Chain,
    LigandAtom,
    MacromolecularStructure,
    MissingResidue,
    Residue,
    SmallMolecule,
    StructureComposition,
    StructureModel,
    normalise_element,
)
from delfin.tools.bio._pdb import parse_pdb
from delfin.tools.bio._read import (
    MMCIF_SUFFIXES,
    MOL2_SUFFIXES,
    PDB_SUFFIXES,
    SDF_SUFFIXES,
    read_mmcif,
    read_mol2,
    read_pdb,
    read_sdf,
    read_small_molecule,
    read_structure,
)
from delfin.tools.bio._site import (
    BindingBox,
    BindingSite,
    ResidueRef,
    ResidueSelection,
    binding_site_from_dict,
    residues_around,
)

__all__ = [
    # failure
    "StructureFormatError",
    "StructureTooLargeError",
    # macromolecular model
    "BioAtom",
    "Residue",
    "Chain",
    "StructureModel",
    "MacromolecularStructure",
    "StructureComposition",
    "MissingResidue",
    "AtomRef",
    "BondedLigandLink",
    "WATER_NAMES",
    "normalise_element",
    # small molecules
    "LigandAtom",
    "SmallMolecule",
    # readers (paths)
    "read_structure",
    "read_pdb",
    "read_mmcif",
    "read_small_molecule",
    "read_sdf",
    "read_mol2",
    "PDB_SUFFIXES",
    "MMCIF_SUFFIXES",
    "SDF_SUFFIXES",
    "MOL2_SUFFIXES",
    # readers (text)
    "parse_pdb",
    "parse_mmcif",
    "parse_sdf",
    "parse_mol2",
    # binding sites
    "BindingBox",
    "BindingSite",
    "ResidueRef",
    "ResidueSelection",
    "binding_site_from_dict",
    "residues_around",
]
