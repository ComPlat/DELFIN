"""Typed, dependency-free description of a biomolecular system.

Input: nothing is read here — the readers in :mod:`delfin.tools.bio._pdb`,
:mod:`delfin.tools.bio._mmcif` and :mod:`delfin.tools.bio._ligand` build these
objects.  Output: frozen dataclasses that every later block (docking,
receptor preparation, structure prediction, QM/MM region selection) can agree
on without agreeing on a parser.

Semantics, and the four decisions that make this layer honest about the two
formats it reads:

* **A residue is identified by chain, sequence number, insertion code, name
  and record type.**  Not by the sequence number: a bound ligand written as
  HETATM may carry the number of a polymer residue in the same chain, and an
  insertion code makes two polymer residues share one number.  Keying on the
  number alone merges residues that the file keeps apart.
* **Alternate locations stay in the residue.**  An altloc is two coordinate
  sets for one atom, so it changes neither the residue count nor the atom
  identity.  :attr:`Residue.atoms` holds every record; selecting one
  conformation is :meth:`Residue.atoms_for_altloc`, an explicit act.
* **Polymer / hetero / water follows the file's own record type and residue
  name**, never a built-in list of standard residue names.  A modified
  residue (MSE, phosphoserine) is written as HETATM and is reported as a
  hetero component, which is what the file asserts; a block that needs a
  chemistry-level polymer definition applies its own rule on top.
* **Element symbols are normalised to the canonical spelling** ("ZN" → "Zn"),
  so a comparison against an element never depends on which format the system
  came from.

Rejected alternatives: a residue keyed by its sequence number only (merges a
ligand onto a polymer residue, the single most common silent error in this
format family); dropping all but the first altloc at read time (changes the
occupancy-weighted composition without saying so); mutable objects (a
structure passed between blocks must not be edited in place by one of them).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Dict, Iterator, Optional, Tuple

# Residue names that denote solvent water.  Counted apart from the other
# hetero components because every downstream block treats them separately
# (kept, stripped, or re-placed) and because a water count is the cheapest
# signal of whether a structure has been prepared at all.
WATER_NAMES = frozenset({"HOH", "DOD", "WAT", "H2O", "TIP", "TIP3", "SOL"})

_HYDROGEN_SYMBOLS = frozenset({"H", "D", "T"})


def normalise_element(symbol: str) -> str:
    """Return *symbol* in the canonical element spelling ("ZN" → "Zn").

    Input: an element symbol as a file writes it (PDB columns 77-78 are
    upper case, mmCIF ``type_symbol`` and SYBYL types are mixed case).
    Output: first letter upper, remainder lower; the empty string stays empty.
    """
    s = str(symbol).strip()
    return s[:1].upper() + s[1:].lower() if s else ""


@dataclass(frozen=True)
class BioAtom:
    """One ATOM / HETATM record.

    ``altloc`` is the alternate-location indicator ("" when the record carries
    none); ``is_hetero`` is the record type the file used, not a judgement
    about the chemistry.
    """

    serial: int
    name: str
    element: str
    x: float
    y: float
    z: float
    altloc: str = ""
    occupancy: float = 1.0
    b_factor: float = 0.0
    is_hetero: bool = False

    @property
    def position(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_hydrogen(self) -> bool:
        return self.element in _HYDROGEN_SYMBOLS


@dataclass(frozen=True)
class Residue:
    """A residue, a bound ligand, an ion or a water — one grouped record set.

    ``is_hetero`` stays the record type the file used (``HETATM`` vs
    ``ATOM``), not a judgement about the chemistry.  A residue the file
    *declares* modified — PDB ``MODRES`` or mmCIF
    ``_pdbx_struct_mod_residue`` — keeps that record type here, and the
    classification (polymer / hetero / water) applies the declaration on
    top: ``modified_residues``, ``parent_name`` and the polymer counting in
    :class:`Chain` and :meth:`MacromolecularStructure.composition` read the
    declaration, because the file itself says this residue is a polymer
    residue with a chemical modification.  ``name`` is then the parent's
    (MSE → MET); the name the records carry is kept in ``original_name``.
    """

    chain_id: str
    seq_id: int
    insertion_code: str
    name: str
    is_hetero: bool
    atoms: Tuple[BioAtom, ...] = ()
    original_name: str = ""
    parent_name: str = ""

    @property
    def is_modified(self) -> bool:
        """The file declares this residue a modification of a parent."""
        return bool(self.parent_name)
    @property
    def key(self) -> Tuple[str, int, str, str, bool]:
        """The identity this layer groups records by (see the module docstring)."""
        return (self.chain_id, self.seq_id, self.insertion_code, self.name,
                self.is_hetero)

    @property
    def label(self) -> str:
        """Short human-readable identity, e.g. ``A/52A/SER``."""
        return f"{self.chain_id}/{self.seq_id}{self.insertion_code}/{self.name}"

    @property
    def is_water(self) -> bool:
        return self.name.upper() in WATER_NAMES

    @property
    def altlocs(self) -> Tuple[str, ...]:
        """Sorted alternate-location indicators present in this residue."""
        return tuple(sorted({a.altloc for a in self.atoms if a.altloc}))

    @property
    def has_hydrogens(self) -> bool:
        return any(a.is_hydrogen for a in self.atoms)

    def atoms_for_altloc(self, code: str) -> Tuple[BioAtom, ...]:
        """Atoms of one conformation: those with no altloc plus those with *code*.

        A record with a blank altloc belongs to every conformation, which is
        what the format means by leaving the column empty.
        """
        return tuple(a for a in self.atoms if not a.altloc or a.altloc == code)

    @property
    def primary_atoms(self) -> Tuple[BioAtom, ...]:
        """One atom per position: the first alternate location, or all atoms.

        "First" is the lowest indicator (A before B), the convention the
        formats themselves use for the preferred conformation.  Occupancy is
        deliberately not consulted: two conformations at 0.50 would make the
        choice depend on file order.
        """
        alt = self.altlocs
        return self.atoms_for_altloc(alt[0]) if alt else self.atoms


@dataclass(frozen=True)
class Chain:
    """All residues written under one chain identifier, in file order."""

    chain_id: str
    residues: Tuple[Residue, ...] = ()

    @property
    def polymer_residues(self) -> Tuple[Residue, ...]:
        return tuple(r for r in self.residues
                     if (not r.is_hetero or r.is_modified) and not r.is_water)

    @property
    def modified_residues(self) -> Tuple[Residue, ...]:
        """Residues the file declares modified, polymer ones only.

        A MODRES entry whose residue is not in the coordinate records is not
        held by this chain (nothing modelled, nothing to report)."""
        return tuple(r for r in self.residues if r.is_modified)

    @property
    def hetero_residues(self) -> Tuple[Residue, ...]:
        """Hetero components that are not solvent water and not declared
        modifications of a polymer residue."""
        return tuple(r for r in self.residues
                     if r.is_hetero and not r.is_water and not r.is_modified)

    @property
    def waters(self) -> Tuple[Residue, ...]:
        return tuple(r for r in self.residues if r.is_water)

    @property
    def atoms(self) -> Tuple[BioAtom, ...]:
        return tuple(a for r in self.residues for a in r.atoms)


@dataclass(frozen=True)
class StructureModel:
    """One model of a structure (an NMR ensemble member, a docking pose, …)."""

    model_id: int
    chains: Tuple[Chain, ...] = ()

    @property
    def residues(self) -> Tuple[Residue, ...]:
        return tuple(r for c in self.chains for r in c.residues)

    @property
    def atoms(self) -> Tuple[BioAtom, ...]:
        return tuple(a for c in self.chains for a in c.atoms)

    def chain(self, chain_id: str) -> Optional[Chain]:
        for c in self.chains:
            if c.chain_id == chain_id:
                return c
        return None


@dataclass(frozen=True)
class MissingResidue:
    """A residue the file itself reports as not located in the experiment.

    Read from PDB ``REMARK 465`` and from the mmCIF category
    ``_pdbx_unobs_or_zero_occ_residues``.  It is what the depositor states,
    not a gap this layer infers from the numbering — an inferred gap cannot
    distinguish a disordered loop from a construct that never had those
    residues.
    """

    name: str
    chain_id: str
    seq_id: int
    insertion_code: str = ""
    model_id: int = 1

    @property
    def label(self) -> str:
        return f"{self.chain_id}/{self.seq_id}{self.insertion_code}/{self.name}"


@dataclass(frozen=True)
class AtomRef:
    """One end of a connection record: which atom of which residue."""

    chain_id: str
    residue_name: str
    seq_id: int
    insertion_code: str = ""
    atom_name: str = ""
    altloc: str = ""

    @property
    def label(self) -> str:
        return (f"{self.chain_id}/{self.seq_id}{self.insertion_code}/"
                f"{self.residue_name}/{self.atom_name}")


@dataclass(frozen=True)
class BondedLigandLink:
    """A connection the file records between two residues.

    ``kind`` is the file's own word: ``link`` for a PDB LINK record,
    ``covale`` / ``metalc`` / ``disulf`` for an mmCIF ``_struct_conn`` entry.
    Non-bonded connection types (hydrogen bonds, salt bridges) are not kept:
    a block that must know whether a ligand is covalently attached would
    otherwise read a hydrogen bond as an attachment.
    """

    atom1: AtomRef
    atom2: AtomRef
    distance: Optional[float] = None
    kind: str = "link"


@dataclass(frozen=True)
class StructureComposition:
    """What one model of a structure holds, as comparable, printable numbers.

    Every field is a tuple or a scalar so two compositions compare with ``==``
    and serialise to JSON without a custom encoder.  Mappings are tuples of
    pairs for the same reason.
    """

    n_models: int
    model_id: int
    chain_ids: Tuple[str, ...]
    polymer_residues: Tuple[Tuple[str, int], ...]      # (chain_id, count)
    n_polymer_residues: int
    hetero_components: Tuple[Tuple[str, int], ...]     # (residue name, count)
    n_hetero_residues: int
    n_waters: int
    n_atoms: int
    n_atoms_primary: int
    n_hydrogens: int
    has_hydrogens: bool
    n_missing_residues: int
    n_bonded_ligand_links: int
    n_modified_residues: int
    modified_residues: Tuple[Tuple[str, int, str, str, str], ...]  # (chain, seq, icode, name, parent)
    altloc_residues: Tuple[str, ...]

    def to_dict(self) -> Dict[str, object]:
        """JSON-ready view; the pair tuples become objects keyed by name."""
        return {
            "n_models": self.n_models,
            "model_id": self.model_id,
            "chain_ids": list(self.chain_ids),
            "polymer_residues": {c: n for c, n in self.polymer_residues},
            "n_polymer_residues": self.n_polymer_residues,
            "hetero_components": {c: n for c, n in self.hetero_components},
            "n_hetero_residues": self.n_hetero_residues,
            "n_waters": self.n_waters,
            "n_atoms": self.n_atoms,
            "n_atoms_primary": self.n_atoms_primary,
            "n_hydrogens": self.n_hydrogens,
            "has_hydrogens": self.has_hydrogens,
            "n_missing_residues": self.n_missing_residues,
            "n_bonded_ligand_links": self.n_bonded_ligand_links,
            "n_modified_residues": self.n_modified_residues,
            "modified_residues": [
                {"chain_id": c, "seq_id": s, "insertion_code": i,
                 "name": n, "parent_name": p}
                for c, s, i, n, p in self.modified_residues
            ],
            "altloc_residues": list(self.altloc_residues),
        }


@dataclass(frozen=True)
class MacromolecularStructure:
    """A macromolecular system as one file described it.

    ``source_name`` is the file's base name only: a structure travels between
    blocks and into reports, and the directory it was read from is not part of
    what the file says.
    """

    models: Tuple[StructureModel, ...]
    missing_residues: Tuple[MissingResidue, ...] = ()
    bonded_ligands: Tuple[BondedLigandLink, ...] = ()
    entry_id: str = ""
    source_format: str = ""
    source_name: str = ""

    @property
    def first_model(self) -> StructureModel:
        return self.models[0]

    def model(self, model_id: int) -> StructureModel:
        """The model with this id.  Raises ``KeyError`` when the file has none."""
        for m in self.models:
            if m.model_id == model_id:
                return m
        raise KeyError(
            f"no model {model_id} in {self.source_name or 'structure'} "
            f"(models: {[m.model_id for m in self.models]})"
        )

    def residues(self, model_id: Optional[int] = None) -> Iterator[Residue]:
        model = self.first_model if model_id is None else self.model(model_id)
        for chain in model.chains:
            for residue in chain.residues:
                yield residue

    def composition(self, model_id: Optional[int] = None) -> StructureComposition:
        """Composition of one model (the first one unless *model_id* says otherwise).

        The structure-level counts (models, missing residues, connection
        records) are carried along so one object answers "what is in this
        file" without a second call.
        """
        model = self.first_model if model_id is None else self.model(model_id)
        polymer_counts: Dict[str, int] = {}
        hetero_counts: Dict[str, int] = {}
        modified: List[Tuple[str, int, str, str, str]] = []
        n_waters = n_atoms = n_primary = n_hydrogens = n_hetero = 0
        altloc_labels = []
        for chain in model.chains:
            for residue in chain.residues:
                n_atoms += len(residue.atoms)
                n_primary += len(residue.primary_atoms)
                n_hydrogens += sum(1 for a in residue.atoms if a.is_hydrogen)
                if residue.altlocs:
                    altloc_labels.append(residue.label)
                if residue.is_water:
                    n_waters += 1
                elif residue.is_modified:
                    modified.append((chain.chain_id, residue.seq_id,
                                     residue.insertion_code, residue.original_name,
                                     residue.parent_name))
                    polymer_counts[chain.chain_id] = (
                        polymer_counts.get(chain.chain_id, 0) + 1)
                elif residue.is_hetero:
                    hetero_counts[residue.name] = hetero_counts.get(residue.name, 0) + 1
                    n_hetero += 1
                else:
                    polymer_counts[chain.chain_id] = (
                        polymer_counts.get(chain.chain_id, 0) + 1)
        return StructureComposition(
            n_models=len(self.models),
            model_id=model.model_id,
            chain_ids=tuple(c.chain_id for c in model.chains),
            polymer_residues=tuple(sorted(polymer_counts.items())),
            n_polymer_residues=sum(polymer_counts.values()),
            hetero_components=tuple(sorted(hetero_counts.items())),
            n_hetero_residues=n_hetero,
            n_waters=n_waters,
            n_atoms=n_atoms,
            n_atoms_primary=n_primary,
            n_hydrogens=n_hydrogens,
            has_hydrogens=n_hydrogens > 0,
            n_missing_residues=len(self.missing_residues),
            n_bonded_ligand_links=len(self.bonded_ligands),
            n_modified_residues=len(modified),
            modified_residues=tuple(sorted(modified)),
            altloc_residues=tuple(altloc_labels),
        )


def _apply_modifications(model: StructureModel,
                         modifications: Dict[Tuple[str, int, str], Tuple[str, str]],
                         ) -> StructureModel:
    """Rename declared modified residues to their parent's name.

    A MODRES entry maps (chain, seq, icode) → (resName, stdRes).  A residue
    the coordinates hold under that key and whose record name matches the
    declaration's ``resName`` becomes a polymer residue named ``stdRes``
    (MSE → MET), keeping its original name.  A declaration whose residue is
    not in the coordinates, or whose record name differs, changes nothing:
    the file declares, the coordinates hold, and only their agreement
    counts as a modification this layer can see.
    """
    if not modifications:
        return model
    chains = []
    for chain in model.chains:
        residues = []
        for residue in chain.residues:
            declared = modifications.get(
                (chain.chain_id, residue.seq_id, residue.insertion_code))
            if declared and residue.name == declared[0]:
                residues.append(replace(
                    residue,
                    name=declared[1],
                    original_name=residue.name,
                    parent_name=declared[1],
                ))
            else:
                residues.append(residue)
        chains.append(replace(chain, residues=tuple(residues)))
    return replace(model, chains=tuple(chains))


# --- small molecules -------------------------------------------------------


@dataclass(frozen=True)
class LigandAtom:
    """One atom of a small molecule (SDF / MOL2)."""

    name: str
    element: str
    x: float
    y: float
    z: float
    charge: float = 0.0

    @property
    def position(self) -> Tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def is_hydrogen(self) -> bool:
        return self.element in _HYDROGEN_SYMBOLS


@dataclass(frozen=True)
class SmallMolecule:
    """A small molecule as one SDF or MOL2 record described it.

    ``bonds`` are ``(i, j, order)`` with **zero-based** atom indices into
    :attr:`atoms`; both formats number atoms from one, and keeping their
    numbering here would make every consumer subtract one.
    """

    name: str
    atoms: Tuple[LigandAtom, ...] = ()
    bonds: Tuple[Tuple[int, int, int], ...] = ()
    properties: Dict[str, str] = field(default_factory=dict)
    source_format: str = ""
    source_name: str = ""

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def has_hydrogens(self) -> bool:
        return any(a.is_hydrogen for a in self.atoms)

    @property
    def positions(self) -> Tuple[Tuple[float, float, float], ...]:
        return tuple(a.position for a in self.atoms)

    @property
    def formula(self) -> str:
        """Molecular formula in Hill order (C, H, then the rest alphabetically)."""
        counts: Dict[str, int] = {}
        for a in self.atoms:
            if a.element:
                counts[a.element] = counts.get(a.element, 0) + 1

        def piece(sym: str) -> str:
            n = counts[sym]
            return sym if n == 1 else f"{sym}{n}"

        if "C" in counts:
            order = ["C"] + (["H"] if "H" in counts else [])
            order += sorted(s for s in counts if s not in ("C", "H"))
        else:
            order = sorted(counts)
        return "".join(piece(s) for s in order)

    @property
    def centroid(self) -> Tuple[float, float, float]:
        if not self.atoms:
            raise ValueError("a molecule with no atoms has no centroid")
        n = float(len(self.atoms))
        return (sum(a.x for a in self.atoms) / n,
                sum(a.y for a in self.atoms) / n,
                sum(a.z for a in self.atoms) / n)


__all__ = [
    "WATER_NAMES",
    "normalise_element",
    "BioAtom",
    "Residue",
    "Chain",
    "StructureModel",
    "MissingResidue",
    "AtomRef",
    "BondedLigandLink",
    "StructureComposition",
    "MacromolecularStructure",
    "LigandAtom",
    "SmallMolecule",
]
