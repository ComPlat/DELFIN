"""Coordinate manipulation on instruction, so an agent can do what a chemist
asks with a structure.

Public surface (every edit returns the *new* structure plus a short record of
what changed, and raises :class:`GeometryError` for impossible requests --
atoms that do not exist, a ring-central dihedral that would break ring
closure -- giving the reason):

 Read / write
   parse_xyz, to_xyz            bare xyz text (count/comment header optional)
   parse_orca_block, to_orca_block   ORCA ``* xyz {charge} {mult}`` blocks
   build_from_smiles            embedded 3D structure via RDKit
 Merge / align
   merge_fragments, align_rmsd
 Internal-coordinate edits (move the SMALLER fragment)
   set_bond, set_angle, set_dihedral
 Whole-fragment motions
   translate, rotate, mirror
 Atoms
   add_atom, remove_atoms, replace_ligand
 Constrained optimisation
   freeze_atoms, constraints_block
 Measurement (also the reference the tests compare against)
   measured_length, measured_angle, measured_dihedral

Conventions: atom indices are 0-based and follow the element order of the
input block; coordinates are in angstrom; angles and dihedrals in degrees
(signed, (-180, 180] for dihedrals); ``change`` records ``op``/``atoms_moved``
and the old/new internal coordinate where one exists.

Reuses by import delfin/xyz_io (element labels, atom parsing) and RDKit
(covalent radii, SMILES embedding); owns nothing, edits nothing of theirs.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:  # numpy is present in both CI interpreters
    import numpy as np
except Exception:  # pragma: no cover - not exercised on either CI interpreter
    np = None  # type: ignore[assignment]

from delfin.xyz_io import _elem_from_label, _parse_xyz_atoms  # reuse by import


class GeometryError(ValueError):
    """An edit was impossible: a missing atom, a ring-breaking dihedral, ..."""


@dataclass
class Geom:
    """A structure plus the metadata an edit has to carry along.

    ``coords`` is a list of 3-tuples in angstrom (tuples so plain tuple
    equality works on them); ``bonds`` is an optional {(lo, hi): order} map
    that operations preserve and extending edits update.
    """

    symbols: List[str]
    coords: List[Tuple[float, float, float]]
    charge: int = 0
    multiplicity: int = 1
    comment: str = ""
    bonds: Dict[Tuple[int, int], int] = field(default_factory=dict)
    frozen: set = field(default_factory=set)  # 0-based atom indices

    def __len__(self) -> int:
        return len(self.symbols)

    def _copy(self) -> "Geom":
        return replace(
            self,
            symbols=list(self.symbols),
            coords=[tuple(c) for c in self.coords],
            bonds=dict(self.bonds),
            frozen=set(self.frozen),
        )


# ---------------------------------------------------------------------------
# small vector helpers (numpy-backed)
# ---------------------------------------------------------------------------

def _arr(p) -> "np.ndarray":
    return np.asarray(p, dtype=float)


def _norm(v) -> "np.ndarray":
    n = np.linalg.norm(v)
    if n == 0:
        raise GeometryError("undefined: zero-length direction vector")
    return v / n


def _sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _dot(a, b) -> float:
    return sum(x * y for x, y in zip(a, b))


def _dist(a, b) -> float:
    return math.dist(a, b)


# ---------------------------------------------------------------------------
# measured internal coordinates (reference used by the tests)
# ---------------------------------------------------------------------------

def measured_length(g: Geom, i: int, j: int) -> float:
    _require_atoms(g, i, j)
    return _dist(g.coords[i], g.coords[j])


def measured_angle(g: Geom, i: int, j: int, k: int) -> float:
    """Angle vertex ``j`` between atoms ``i`` and ``k``, degrees in [0, 180]."""
    _require_atoms(g, i, j, k)
    a = _arr(_sub(g.coords[i], g.coords[j]))
    b = _arr(_sub(g.coords[k], g.coords[j]))
    c = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def measured_dihedral(g: Geom, i: int, j: int, k: int, l: int) -> float:
    """Signed dihedral i-j-k-l in degrees, (-180, 180].

    Matches the sign of ASE ``get_dihedrals`` (verified numerically; a
    clockwise turn about the j->k axis viewed from j increases the value).
    """
    _require_atoms(g, i, j, k, l)
    p0, p1, p2, p3 = (_arr(g.coords[x]) for x in (i, j, k, l))
    b2 = p2 - p1
    n2 = _norm(np.cross(b2, p3 - p2)) if np.linalg.norm(p3 - p2) else b2
    b1 = p1 - p0
    n1 = _norm(np.cross(b1, b2)) if np.linalg.norm(b1) else b2
    m1 = _norm(np.cross(n1, b2))
    x = float(np.dot(n1, n2))
    y = float(np.dot(m1, n2))
    deg = float(np.degrees(np.arctan2(y, x)))
    # arctan2 gives (-180, 180]; keep 180 as +180
    if deg <= -180.0:
        deg = 180.0
    return deg


def _require_atoms(g: Geom, *idx: int) -> None:
    for a in idx:
        if not (isinstance(a, int) and not isinstance(a, bool) and 0 <= a < len(g)):
            raise GeometryError(
                f"atom index {a} does not exist: structure has {len(g)} atoms "
                f"(0..{len(g) - 1})"
            )


# ---------------------------------------------------------------------------
# read / write
# ---------------------------------------------------------------------------

def parse_xyz(text: str, *, charge: int = 0, multiplicity: int = 1,
              comment: str = "") -> Geom:
    """Parse an xyz block (count+comment header optional, trailing ``*`` ok)."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise GeometryError("empty structure: no atoms to parse")
    start, body = 0, lines
    n = None
    if len(lines) >= 2 and lines[0].strip().isdigit():
        n = int(lines[0].strip())
        start, body = 1, lines[1:]
    hdr = lines[1] if start else ""
    atoms = _parse_xyz_atoms(body)
    if n is not None and len(atoms) != n:
        raise GeometryError(
            f"header says {n} atoms but the block holds {len(atoms)}"
        )
    if not atoms:
        raise GeometryError("no atoms found in structure block")
    sym = [_elem_from_label(a["elem"]) for a in atoms]
    coords = [(float(a["x"]), float(a["y"]), float(a["z"])) for a in atoms]
    comment = comment or (hdr.strip() if not hdr.lstrip().isdigit() else "")
    return Geom(symbols=sym, coords=coords, charge=charge,
                multiplicity=multiplicity, comment=comment)


def to_xyz(g: Geom, *, comment: Optional[str] = None) -> str:
    """Serialise to xyz text with a count/comment header."""
    c = g.comment if comment is None else comment
    lines = [str(len(g)), c or ""]
    for sym, (x, y, z) in zip(g.symbols, g.coords):
        lines.append(f"{sym} {x:.8f} {y:.8f} {z:.8f}")
    return "\n".join(lines) + "\n"


def parse_orca_block(text: str) -> Geom:
    """Parse an ORCA ``* xyz {charge} {mult}`` .. ``*`` coordinate block."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise GeometryError("empty ORCA block")
    head = lines[0]
    m = re.match(r"^\*\s*(?:xyz|xyzfile|int|frac)\s+(.+?)\s*$", head)
    if not m:
        raise GeometryError(
            f"not an ORCA coordinate block (expected '* xyz charge mult'): {head!r}"
        )
    spec = m.group(1).split()
    try:
        charge, mult = int(spec[0]), int(spec[1])
    except (ValueError, IndexError) as exc:  # pragma: no cover - defensive
        raise GeometryError(f"malformed charge/multiplicity in {head!r}") from exc
    body = [ln for ln in lines[1:] if ln != "*"]
    return parse_xyz("\n".join(body), charge=charge, multiplicity=mult)


def to_orca_block(g: Geom, *, charge: Optional[int] = None,
                  multiplicity: Optional[int] = None) -> str:
    """Serialise as an ORCA ``* xyz {charge} {mult}`` .. ``*`` block."""
    ch = g.charge if charge is None else charge
    mu = g.multiplicity if multiplicity is None else multiplicity
    lines = [f"* xyz {ch} {mu}"]
    for sym, (x, y, z) in zip(g.symbols, g.coords):
        lines.append(f"{sym}   {x:.8f}   {y:.8f}   {z:.8f}")
    lines.append("*")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# covalent-radius connectivity used to split fragments and detect rings
# ---------------------------------------------------------------------------

def _covalent_radii(g: Geom) -> List[float]:
    """Angstrom covalent radii per atom (RDKit), fallback table if RDKit is not
    importable so the module still works on bare structures."""
    try:
        from rdkit import Chem
    except Exception:  # pragma: no cover - RDKit present on both CI interpreters
        Chem = None  # type: ignore[assignment]
    if Chem is not None:
        pt = Chem.GetPeriodicTable()
    radii = []
    for sym in g.symbols:
        r = None
        if Chem is not None:
            try:
                r = pt.GetRcovalent(int(pt.GetAtomicNumber(sym)))
            except Exception:
                r = None
        radii.append(float(r) if r else _FALLBACK_RADII.get(sym, 1.0))
    return radii


_FALLBACK_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05,
    "P": 1.07, "Cl": 1.02, "Br": 1.20, "I": 1.39, "B": 0.84,
}


def _bond_graph(g: Geom, tolerance: float = 0.45) -> Dict[int, List[int]]:
    """Neighbourhood map from covalent radii + a small tolerance."""
    radii = _covalent_radii(g)
    nbrs: Dict[int, List[int]] = {a: [] for a in range(len(g))}
    for i in range(len(g)):
        for j in range(i + 1, len(g)):
            if _dist(g.coords[i], g.coords[j]) <= radii[i] + radii[j] + tolerance:
                nbrs[i].append(j)
                nbrs[j].append(i)
    return nbrs


def _connected_component(
    g: Geom, start: int, excluded: set, nbrs: Optional[Dict[int, List[int]]] = None
) -> set:
    """Atoms reachable from ``start`` without entering ``excluded``."""
    nbrs = nbrs or _bond_graph(g)
    seen = {start}
    stack = [start]
    while stack:
        a = stack.pop()
        for b in nbrs.get(a, []):
            if b not in excluded and b not in seen:
                seen.add(b)
                stack.append(b)
    return seen


def _cutted_connected(g: Geom, a: int, b: int, nbrs) -> Tuple[set, set]:
    """The two connected components of atom ``a`` and atom ``b`` after slicing
    the edge a-b (a,b themselves excluded on each other's side)."""
    side_a = _connected_component(g, a, excluded={b}, nbrs=nbrs)
    side_b = _connected_component(g, b, excluded={a}, nbrs=nbrs)
    return side_a, side_b


def _edge_in_cycle(g: Geom, i: int, j: int, nbrs) -> bool:
    """True if there is a path from ``i`` to ``j`` that does not use the direct
    edge (i, j) -- i.e. that edge lies on a ring.

    OLD approach (removed): ``_connected_component(g, i, excluded={j})`` pruned
    node *j* entirely, so ``j`` could never be reached and the test was always
    False for any edge. Here we drop only the (i, j) step and allow ``j`` as a
    reachable target via the rest of the ring.
    """
    backrefs = {}
    for a, ns in nbrs.items():
        for b in ns:
            if (a, b) == (i, j) or (a, b) == (j, i):
                continue
            backrefs.setdefault(a, []).append(b)
    seen = {i}
    stack = [i]
    while stack:
        a = stack.pop()
        for b in backrefs.get(a, []):
            if b == j:
                return True
            if b not in seen:
                seen.add(b)
                stack.append(b)
    return False


# ---------------------------------------------------------------------------
# rigid rotations
# ---------------------------------------------------------------------------

def _rotate_points(points: Iterable[Tuple[float, float, float]],
                   axis, angle_deg: float, center) -> List[Tuple[float, float, float]]:
    """Rotate points (angstrom) about an axis through ``center`` by ``angle``
    degrees right-handed about the given axis vector (Rodrigues)."""
    ax = _norm(_arr(axis))
    a = math.radians(angle_deg)
    out = []
    for p in points:
        v = _arr(p) - _arr(center)
        vn = (v * math.cos(a)
              + np.cross(ax, v) * math.sin(a)
              + ax * np.dot(ax, v) * (1 - math.cos(a)))
        r = vn + _arr(center)
        out.append((float(r[0]), float(r[1]), float(r[2])))
    return out


def _translate_points(points: Iterable[Tuple[float, float, float]],
                      vector) -> List[Tuple[float, float, float]]:
    dv = _arr(vector)
    return [(float(p[0] + dv[0]), float(p[1] + dv[1]), float(p[2] + dv[2]))
            for p in points]


# ---------------------------------------------------------------------------
# edits -- each returns (Geom, change)
# ---------------------------------------------------------------------------

def _apply_to(g: Geom, new_coords: Dict[int, Tuple[float, float, float]]) -> Geom:
    out = g._copy()
    for idx, c in new_coords.items():
        out.coords[idx] = (float(c[0]), float(c[1]), float(c[2]))
    return out


def set_bond(g: Geom, i: int, j: int, length: float) -> Tuple[Geom, dict]:
    """Move the smaller fragment so the i-j distance becomes ``length`` (A)."""
    _require_atoms(g, i, j)
    if length <= 0:
        raise GeometryError(f"bond length must be positive, got {length}")
    old = _dist(g.coords[i], g.coords[j])
    if old == 0:
        raise GeometryError("atoms coincide; cannot set a bond length on them")
    nbrs = _bond_graph(g)
    side_i, side_j = _cutted_connected(g, i, j, nbrs)
    # choose the smaller fragment to move; tie -> the j-side
    move = (side_j if len(side_j) <= len(side_i) else side_i)
    axis = _norm(_arr(g.coords[j]) - _arr(g.coords[i]))
    scale = (length / old - 1.0)
    out = g._copy()
    moved = []
    for a in move:
        v = _arr(out.coords[a]) - _arr(out.coords[i])
        nv = v + axis * (_dot(v, axis) * scale)
        r = nv + _arr(out.coords[i])
        out.coords[a] = (float(r[0]), float(r[1]), float(r[2]))
        moved.append(a)
    new = _dist(out.coords[i], out.coords[j])
    return out, {
        "op": "set_bond", "atoms_moved": moved, "old": old, "new": new,
        "unit": "angstrom", "type": "bond",
    }


def set_angle(g: Geom, i: int, j: int, k: int, degrees: float) -> Tuple[Geom, dict]:
    """Move the smaller fragment so the angle i-j-k (vertex j) = ``degrees``."""
    _require_atoms(g, i, j, k)
    old = measured_angle(g, i, j, k)
    if not (0.0 < degrees < 180.0):
        raise GeometryError(f"angle must be in (0,180), got {degrees}")
    nbrs = _bond_graph(g)
    # the two sides swing about the vertex j
    side_i = _connected_component(g, i, excluded={j}, nbrs=nbrs)
    side_k = _connected_component(g, k, excluded={j}, nbrs=nbrs)
    move = side_i if len(side_i) <= len(side_k) else side_k
    move = side_i  # tie-break: move the i side (deterministic)
    anchor = _arr(g.coords[j])
    a_vec = _arr(g.coords[i]) - anchor       # i points away from vertex
    b_vec = _arr(g.coords[k]) - anchor
    axis = _norm(np.cross(a_vec, b_vec))       # rotation axis through j
    if np.linalg.norm(np.cross(a_vec, b_vec)) < 1e-9:
        raise GeometryError("angle atoms are collinear; cannot set an angle")
    delta = degrees - old
    # which ray is being moved decides the rotation hand:
    #   i-side: rotate by -delta about n=cross(a,b) (verified: +delta swings i
    #           toward k and lowers the angle, -delta raises it)
    #   k-side: the mirror image, so rotate by +delta instead
    move_is_i = (move is side_i)
    move_by = (-delta if move_is_i else delta)
    newcoords: Dict[int, Tuple[float, float, float]] = {}
    for a in move:
        if a == j:
            continue
        p = _rotate_points([g.coords[a]], axis, move_by, tuple(anchor))[0]
        newcoords[a] = p
    out = _apply_to(g, newcoords)
    moved = [a for a in sorted(move) if a != j]
    return out, {
        "op": "set_angle", "atoms_moved": moved, "old": old, "new": degrees,
        "unit": "deg", "type": "angle",
    }


def set_dihedral(g: Geom, i: int, j: int, k: int, l: int,
                 degrees: float) -> Tuple[Geom, dict]:
    """Move the smaller fragment about the j-k bond so the dihedral = degrees."""
    _require_atoms(g, i, j, k, l)
    if len({i, j, k, l}) < 4:
        raise GeometryError("a dihedral needs four distinct atoms")
    nbrs = _bond_graph(g)
    if _edge_in_cycle(g, j, k, nbrs):
        raise GeometryError(
            f"dihedral across the ring bond ({j}-{k}) would break ring closure; "
            "refusing to rotate one side of a ring"
        )
    old = measured_dihedral(g, i, j, k, l)
    # the two halves lie on the two sides of the j-k bond
    side_a = _connected_component(g, i, excluded={k}, nbrs=nbrs)
    side_b = _connected_component(g, l, excluded={j}, nbrs=nbrs)
    move = side_a if len(side_a) <= len(side_b) else side_b
    axis = _arr(g.coords[k]) - _arr(g.coords[j])     # j -> k
    pivot = tuple(g.coords[j])
    delta = degrees - old
    newcoords: Dict[int, Tuple[float, float, float]] = {}
    for a in move:
        if a in (j, k):
            continue
        p = _rotate_points([g.coords[a]], axis, delta, pivot)[0]
        newcoords[a] = p
    out = _apply_to(g, newcoords)
    moved = [a for a in sorted(move) if a not in (j, k)]
    return out, {
        "op": "set_dihedral", "atoms_moved": moved, "old": old, "new": degrees,
        "unit": "deg", "type": "dihedral",
    }


def translate(g: Geom, indices: Sequence[int], vector) -> Tuple[Geom, dict]:
    idx = _as_indices(g, indices, "translate")
    newcoords = {
        a: tuple(_arr(g.coords[a]) + _arr(vector)) for a in idx
    }
    return _apply_to(g, newcoords), {"op": "translate", "atoms_moved": list(idx)}


def rotate(g: Geom, indices: Sequence[int], axis, angle_deg: float,
           center=None) -> Tuple[Geom, dict]:
    """Rotate ``indices`` about ``axis`` (a 3-vector) by ``angle`` degrees.
    ``center`` defaults to the centroid of the moved atoms."""
    idx = _as_indices(g, indices, "rotate")
    if axis is None or len(list(axis)) != 3 or math.sqrt(sum(x * x for x in axis)) == 0:
        raise GeometryError("rotation axis must be a non-zero 3-vector")
    c = tuple(center) if center is not None else _centroid(g, idx)
    newcoords = {a: p for a, p in zip(
        idx, _rotate_points([g.coords[a] for a in idx], axis, angle_deg, c))}
    return _apply_to(g, newcoords), {
        "op": "rotate", "atoms_moved": list(idx), "angle_deg": float(angle_deg),
    }


def _centroid(g: Geom, idx: Iterable[int]) -> Tuple[float, float, float]:
    cs = [_arr(g.coords[a]) for a in idx]
    return tuple(float(x) for x in np.mean(cs, axis=0))


def _as_indices(g: Geom, indices, op: str) -> List[int]:
    if isinstance(indices, (int,)) and not isinstance(indices, bool):
        indices = [indices]
    if isinstance(indices, str):
        raise GeometryError(f"{op}: integer indices expected, got a string")
    idx = list(indices)
    _require_atoms(g, *idx)
    return sorted(set(idx))


# ---------------------------------------------------------------------------
# atoms
# ---------------------------------------------------------------------------

def add_atom(g: Geom, element: str, bonded_to: int, bond_length: float) -> Tuple[Geom, dict]:
    _require_atoms(g, bonded_to)
    sym = _elem_from_label(str(element))
    if bond_length <= 0:
        raise GeometryError(f"bond_length must be positive, got {bond_length}")
    anchor = _arr(g.coords[bonded_to])
    # place along a random-but-deterministic non-collinear direction from the
    # anchor, using the anchor's existing neighbours to pick a sensible vector
    nbrs = sorted(_bond_graph(g).get(bonded_to, []))
    if nbrs:
        ref = _arr(g.coords[nbrs[0]]) - anchor
    else:
        ref = _arr((1.0, 0.0, 0.0))
    cand = np.cross(ref, _arr((0.0, 0.0, 1.0)))
    if np.linalg.norm(cand) < 1e-9:
        cand = np.cross(ref, _arr((1.0, 0.0, 0.0)))
    pos = anchor + _norm(cand) * bond_length
    out = g._copy()
    new_idx = len(out)
    out.symbols.append(sym)
    out.coords.append((float(pos[0]), float(pos[1]), float(pos[2])))
    out.bonds[(min(new_idx, bonded_to), max(new_idx, bonded_to))] = 1
    return out, {
        "op": "add_atom", "atoms_moved": [new_idx], "added": [new_idx],
        "element": sym,
    }


def remove_atoms(g: Geom, indices: Sequence[int]) -> Tuple[Geom, dict]:
    idx = _as_indices(g, indices, "remove_atoms")
    keep = [a for a in range(len(g)) if a not in set(idx)]
    out = Geom(
        symbols=[g.symbols[a] for a in keep],
        coords=[g.coords[a] for a in keep],
        charge=g.charge, multiplicity=g.multiplicity, comment=g.comment,
    )
    # rebuild frozen as the kept (renumbered) set
    out.frozen = {k for k, a in enumerate(keep) if a in g.frozen}
    # remap bonds
    remap = {a: n for n, a in enumerate(keep)}
    out.bonds = {tuple(sorted((remap[x], remap[y]))): o
                 for (x, y), o in g.bonds.items()
                 if x in remap and y in remap}
    return out, {"op": "remove_atoms", "atoms_moved": list(idx), "removed": list(idx)}


# ---------------------------------------------------------------------------
# ligand swap
# ---------------------------------------------------------------------------

def _ligand_component(g: Geom, connect: int, neighbour: int,
                      nbrs) -> set:
    """The peripheral fragment of ``neighbour`` that hangs off ``connect``."""
    return _connected_component(g, neighbour, excluded={connect}, nbrs=nbrs)


def replace_ligand(g: Geom, connect: int, ligand) -> Tuple[Geom, dict]:
    """Replace the ligand on ``connect``.

    ``ligand`` is an element symbol (a single atom) -- the common case -- or a
    SMILES string, which is embedded with RDKit and attached at ``connect``.
    The replaced fragment is dropped; the new ligand sits along the bond that
    used to hold the old one.
    """
    _require_atoms(g, connect)
    nbrs = _bond_graph(g)
    bonded = nbrs.get(connect, [])
    # choose the dangling neighbour whose component is a pendant (not the rest
    # of the molecule): prefer a non-carbon heteroatom, else the smallest side
    def comp_size(n):
        return len(_ligand_component(g, connect, n, nbrs))

    choice = None
    for n in bonded:
        if g.symbols[n] != "C":
            choice = n
            break
    if choice is None and bonded:
        choice = min(bonded, key=comp_size)
    if choice is None:
        raise GeometryError(
            f"atom {connect} has no bonded neighbour to treat as a ligand"
        )
    old_side = sorted(_ligand_component(g, connect, choice, nbrs))
    old_vec = _arr(g.coords[choice]) - _arr(g.coords[connect])  # ligand direction

    # build the replacement coordinates
    if isinstance(ligand, str) and _single_symbol(ligand):
        sym = _elem_from_label(ligand)
        radii = _covalent_radii(g)
        bond_len = radii[connect] + _FALLBACK_RADII.get(sym, 1.0)
        pos = _arr(g.coords[connect]) + _norm(old_vec) * bond_len
        add_syms, add_coords = [sym], [tuple(pos)]
    else:  # SMILES -> embed, then place the centroid on the ligand direction
        base = build_from_smiles(str(ligand))
        add_syms, add_coords = list(base.symbols), list(base.coords)
        cent = _centroid(base, range(len(base)))
        bond_len = _dist(g.coords[connect], tuple(cent))
        shift = _arr(g.coords[connect]) + _norm(old_vec) * max(bond_len, 1.5) - cent
        add_coords = [_translate_points([c], shift)[0] for c in add_coords]

    out = g._copy()
    # drop old ligand atoms
    drop = set(old_side)
    keep = [a for a in range(len(out)) if a not in drop]
    out.symbols = [out.symbols[a] for a in keep]
    out.coords = [out.coords[a] for a in keep]
    remap = {a: n for n, a in enumerate(keep)}
    connect_new = remap[connect]
    out.bonds = {tuple(sorted((remap[x], remap[y]))): o
                 for (x, y), o in out.bonds.items() if x in remap and y in remap}
    # append the new ligand atoms
    added = []
    for s, c in zip(add_syms, add_coords):
        idx = len(out)
        out.symbols.append(s)
        out.coords.append(c)
        if idx != connect_new:
            out.bonds[tuple(sorted((connect_new, idx)))] = 1
        added.append(idx)
    return out, {
        "op": "replace_ligand", "atoms_moved": added, "removed": old_side,
        "added": added, "connect": connect_new,
    }


def _single_symbol(s: str) -> bool:
    return bool(re.fullmatch(r"[A-Z][a-z]?(?:\([^)]+\))?", s.strip()))


# ---------------------------------------------------------------------------
# mirror
# ---------------------------------------------------------------------------

def mirror(g: Geom, plane: str = "yz") -> Tuple[Geom, dict]:
    """Mirror all atoms across one of the coordinate planes."""
    p = (plane or "yz").strip().lower()
    if p not in ("xy", "xz", "yz"):
        raise GeometryError(f"unknown mirror plane {p!r}; use 'xy','xz' or 'yz'")
    flip = {"xy": 2, "xz": 1, "yz": 0}[p]
    out = g._copy()
    for a in range(len(out)):
        c = list(out.coords[a])
        c[flip] = -c[flip]
        out.coords[a] = (c[0], c[1], c[2])
    # mirroring turns right-handed into left-handed: drop any absolute
    # stereochemical expectations by not carrying per-atom handedness here
    return out, {"op": "mirror", "plane": p, "atoms_moved": list(range(len(g)))}


# ---------------------------------------------------------------------------
# SMILES / merge / align
# ---------------------------------------------------------------------------

def build_from_smiles(smiles: str) -> Geom:
    """Embed a molecule from SMILES with RDKit (UFF force field for a sensible
    starting geometry). Raises GeometryError when RDKit cannot parse/embed."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:
        raise GeometryError(f"RDKit is required to build from SMILES: {exc}") from exc
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        raise GeometryError(f"RDKit could not parse SMILES {smiles!r}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=0) != 0:
        AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=0)
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except Exception:
        pass  # embed only is fine; coordinates must stay finite
    conf = mol.GetConformer()
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = []
    for a in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(a)
        coords.append((float(p.x), float(p.y), float(p.z)))
    # derive connectivity/bonds from the RDKit molecule
    bonds = {}
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonds[(min(i, j), max(i, j))] = int(round(b.GetBondTypeAsDouble()))
    return Geom(symbols=syms, coords=coords, bonds=bonds)


def merge_fragments(a: Geom, b: Geom, *rest: Geom) -> Tuple[Geom, dict]:
    """Concatenate two (or more) structures into one, atoms of ``a`` first."""
    parts = [a, b, *rest]
    off = 0
    symbols: List[str] = []
    coords = []
    bonds: Dict[Tuple[int, int], int] = {}
    for p in parts:
        for sym in p.symbols:
            symbols.append(sym)
        for c in p.coords:
            coords.append(tuple(c))
        for (x, y), o in p.bonds.items():
            bonds[(x + off, y + off)] = o
        off += len(p)
    out = Geom(symbols=symbols, coords=coords,
               charge=a.charge, multiplicity=a.multiplicity,
               comment=a.comment, bonds=bonds)
    return out, {"op": "merge_fragments", "atoms_moved": [], "added": range(len(a), len(out))}


def align_rmsd(mobile: Geom, reference: Geom):
    """Align ``mobile`` onto ``reference`` (Kabsch) when both have the same
    atom order. Returns ``(aligned_geom, rmsd_angstrom)``."""
    if len(mobile) != len(reference):
        raise GeometryError(
            f"align needs matching atom counts: mobile has {len(mobile)}, "
            f"reference {len(reference)}"
        )
    P = np.asarray(mobile.coords, dtype=float)
    Q = np.asarray(reference.coords, dtype=float)
    mu_p = P.mean(axis=0)
    mu_q = Q.mean(axis=0)
    Pc = P - mu_p
    Qc = Q - mu_q
    cov = Qc.T @ Pc
    U, _S, Vt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = U @ D @ Vt
    rotated = (R @ Pc.T).T + mu_q
    aligned = Geom(
        symbols=list(mobile.symbols),
        coords=[(float(c[0]), float(c[1]), float(c[2])) for c in rotated],
        charge=mobile.charge, multiplicity=mobile.multiplicity,
        comment=mobile.comment, bonds=dict(mobile.bonds),
    )
    diff = rotated - Q
    rmsd = float(math.sqrt((diff * diff).sum() / len(diff)))
    return aligned, rmsd


# ---------------------------------------------------------------------------
# freeze atoms -> constrained optimisation spec
# ---------------------------------------------------------------------------

def freeze_atoms(g: Geom, indices: Sequence[int]) -> Tuple[Geom, dict]:
    """Mark atoms as frozen for a following constrained optimisation. Moving
    nothing now, it records the constraints so ``constraints_block`` can emit
    an ORCA ``%geom Constraints`` block."""
    idx = _as_indices(g, indices, "freeze_atoms")
    out = g._copy()
    out.frozen = set(idx)
    return out, {"op": "freeze_atoms", "atoms_moved": [], "frozen": sorted(idx)}


def constraints_block(g: Geom) -> str:
    """ORCA ``%geom Constraints`` block for the frozen atoms (0-based Cartesian
    constraints), e.g::

        %geom
          Constraints
            { C 0 C }
          end
        end
    """
    frozen = sorted(g.frozen)
    if not frozen:
        return ""
    lines = ["%geom", "  Constraints"]
    for idx in frozen:
        lines.append(f"    {{ C {idx} C }}")
    lines.append("  end")
    lines.append("end")
    return "\n".join(lines) + "\n"
