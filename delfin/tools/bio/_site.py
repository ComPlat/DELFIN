"""Typed descriptors for a binding site, and the two ways of naming one.

Input: coordinates (from a bound ligand, a set of atoms, or numbers a user
typed) or a structure plus a point.  Output: :class:`BindingBox` — a centre
and an edge length in Angstrom — or :class:`ResidueSelection` — the residues
that define the site.  Both carry a ``kind`` through :meth:`to_dict`, so a
site survives a JSON round trip and can travel between blocks as an artifact.

Semantics: the two descriptors are not interchangeable and neither is derived
from the other here.  A box is what a grid-based program needs and is stated
in Angstrom, never in grid points, because a grid spacing belongs to the
program and not to the site.  A residue selection is what a flexible-residue
setup, a QM region or a pocket comparison needs, and it keeps the residue
identity (chain, number, insertion code, name) rather than atom indices,
which change the moment hydrogens are added or waters are stripped.

Rejected alternatives: one descriptor with optional fields for both forms (a
consumer then has to check which half is filled, and both halves can be
filled inconsistently); a box stored as minimum/maximum corners (every
program that takes a box takes a centre and a size, and the corner form
invites a silently inverted axis); a selection stored as atom indices (see
above).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from delfin.tools.bio._model import MacromolecularStructure, SmallMolecule

_Point = Tuple[float, float, float]


def _as_point(value: Sequence[float], what: str) -> _Point:
    values = tuple(float(v) for v in value)
    if len(values) != 3:
        raise ValueError(f"{what} needs three numbers, got {len(values)}")
    return values  # type: ignore[return-value]


@dataclass(frozen=True)
class BindingBox:
    """An axis-aligned search box: a centre and an edge length, in Angstrom."""

    center: _Point
    size: _Point

    def __post_init__(self) -> None:
        object.__setattr__(self, "center", _as_point(self.center, "box centre"))
        size = _as_point(self.size, "box size")
        if any(s <= 0.0 for s in size):
            raise ValueError(
                f"box edge lengths must be positive Angstrom values, got {size}"
            )
        object.__setattr__(self, "size", size)

    @property
    def kind(self) -> str:
        return "binding_box"

    @property
    def volume(self) -> float:
        """Box volume in cubic Angstrom."""
        return self.size[0] * self.size[1] * self.size[2]

    @classmethod
    def from_points(cls, points: Iterable[Sequence[float]],
                    padding: float = 4.0) -> "BindingBox":
        """Smallest box covering *points*, grown by *padding* on every side.

        The padding is added to both faces of each axis (a 2 A padding makes
        the box 4 A longer), which is how the clearance around a ligand is
        meant.  A flat set of points (a planar ligand, a single atom) still
        yields a usable box because the padding alone is positive.
        """
        values = [tuple(float(c) for c in p) for p in points]
        if not values:
            raise ValueError("a binding box needs at least one point")
        if padding <= 0.0:
            raise ValueError(f"padding must be a positive distance, got {padding}")
        mins = [min(v[i] for v in values) for i in range(3)]
        maxs = [max(v[i] for v in values) for i in range(3)]
        return cls(
            center=tuple((mins[i] + maxs[i]) / 2.0 for i in range(3)),
            size=tuple((maxs[i] - mins[i]) + 2.0 * padding for i in range(3)),
        )

    @classmethod
    def from_molecule(cls, molecule: SmallMolecule,
                      padding: float = 4.0) -> "BindingBox":
        """Box around every atom of *molecule* (the usual "box around the ligand")."""
        return cls.from_points(molecule.positions, padding=padding)

    def contains(self, point: Sequence[float]) -> bool:
        """Whether *point* lies inside the box (its faces count as inside)."""
        p = _as_point(point, "point")
        return all(abs(p[i] - self.center[i]) <= self.size[i] / 2.0 + 1e-9
                   for i in range(3))

    def to_dict(self) -> Dict[str, object]:
        return {"kind": self.kind, "center": list(self.center),
                "size": list(self.size), "unit": "angstrom"}

    @classmethod
    def from_dict(cls, data: Dict[str, object]) -> "BindingBox":
        if data.get("kind") not in (None, "binding_box"):
            raise ValueError(f"not a binding box: kind={data.get('kind')!r}")
        return cls(center=data["center"], size=data["size"])  # type: ignore[arg-type]


@dataclass(frozen=True)
class ResidueRef:
    """A residue named the way a file names it, without its coordinates."""

    chain_id: str
    seq_id: int
    insertion_code: str = ""
    residue_name: str = ""

    @property
    def label(self) -> str:
        return f"{self.chain_id}/{self.seq_id}{self.insertion_code}/{self.residue_name}"

    def to_dict(self) -> Dict[str, object]:
        return {"chain_id": self.chain_id, "seq_id": self.seq_id,
                "insertion_code": self.insertion_code,
                "residue_name": self.residue_name}

    @classmethod
    def from_dict(cls, data: Dict[str, object]) -> "ResidueRef":
        return cls(chain_id=str(data.get("chain_id", "")),
                   seq_id=int(data["seq_id"]),       # type: ignore[arg-type]
                   insertion_code=str(data.get("insertion_code", "")),
                   residue_name=str(data.get("residue_name", "")))


@dataclass(frozen=True)
class ResidueSelection:
    """A binding site given as the residues that line it."""

    residues: Tuple[ResidueRef, ...] = ()

    @property
    def kind(self) -> str:
        return "residue_selection"

    @property
    def labels(self) -> Tuple[str, ...]:
        return tuple(r.label for r in self.residues)

    def to_dict(self) -> Dict[str, object]:
        return {"kind": self.kind,
                "residues": [r.to_dict() for r in self.residues]}

    @classmethod
    def from_dict(cls, data: Dict[str, object]) -> "ResidueSelection":
        if data.get("kind") not in (None, "residue_selection"):
            raise ValueError(f"not a residue selection: kind={data.get('kind')!r}")
        rows = data.get("residues") or []
        return cls(tuple(ResidueRef.from_dict(r) for r in rows))  # type: ignore[arg-type]


BindingSite = Union[BindingBox, ResidueSelection]


def binding_site_from_dict(data: Dict[str, object]) -> BindingSite:
    """Rebuild whichever descriptor *data* describes, by its ``kind``."""
    kind = data.get("kind")
    if kind == "binding_box":
        return BindingBox.from_dict(data)
    if kind == "residue_selection":
        return ResidueSelection.from_dict(data)
    raise ValueError(f"unknown binding-site kind: {kind!r}")


def residues_around(
    structure: MacromolecularStructure,
    point: Sequence[float],
    radius: float,
    *,
    model_id: Optional[int] = None,
    include_waters: bool = True,
) -> ResidueSelection:
    """Residues with at least one atom within *radius* Angstrom of *point*.

    Input: a parsed structure, a point in the structure's own coordinates and
    a radius in Angstrom.  Output: a :class:`ResidueSelection` in file order.

    Every atom counts, alternate locations included: a residue that reaches
    the site in one of its conformations lines the site.  Waters are included
    by default because a selection reports what is in the sphere; a caller
    that treats them as solvent says so with ``include_waters=False``.
    """
    if radius <= 0.0:
        raise ValueError(f"radius must be a positive distance, got {radius}")
    cx, cy, cz = _as_point(point, "point")
    limit = radius * radius
    found: List[ResidueRef] = []
    for residue in structure.residues(model_id):
        if residue.is_water and not include_waters:
            continue
        for atom in residue.atoms:
            dx, dy, dz = atom.x - cx, atom.y - cy, atom.z - cz
            if dx * dx + dy * dy + dz * dz <= limit:
                found.append(ResidueRef(residue.chain_id, residue.seq_id,
                                        residue.insertion_code, residue.name))
                break
    return ResidueSelection(tuple(found))


__all__ = [
    "BindingBox",
    "BindingSite",
    "ResidueRef",
    "ResidueSelection",
    "binding_site_from_dict",
    "residues_around",
]
