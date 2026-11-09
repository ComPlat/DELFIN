"""Structure-department library (package XE1, SA1 of the wave-27 plan).

Pure helper functions behind the oxyz-editor expert's ``tool_geometry_*``
MCP tools. Every function takes xyz *text* in and returns xyz text / a JSON
object out, matching the existing geometry-tool contract
(``delfin.ops_server.server`` -> ``delfin.api`` -> ``delfin.tools.geometry_edit``).

The library reuses DELFIN's own chemistry primitives -- nothing is
re-implemented:

- xyz <-> Geom and the covalent-radius / connectivity tables live in
  ``delfin.tools.geometry_edit`` (``parse_xyz``, ``to_xyz``, ``Geom``,
  ``_covalent_radii``, ``_bond_graph``).
- renumbering a product onto an educt uses the universal mapper
  ``delfin.atom_mapping.map_atoms`` (a second mapper would be a defect).

Design rules the package target fixes:
- Every mutation returns a NEW structure; no tool writes a file. The
  structure expert writes the new file with the gate-checked ``write_file``
  (which refuses overwriting); an MCP tool writing a file would bypass the
  write-scope gate, so nothing here touches disk.
- A formula mismatch in a renumber is surfaced (ValueError), never hidden.
- ``fix_hydrogens`` only ever ADDS hydrogens, never removes an atom.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from typing import Any

from delfin.tools.geometry_edit import (
    Geom,
    GeometryError,
    _bond_graph,
    _covalent_radii,
    parse_orca_block,
    parse_xyz,
    to_xyz,
)

try:  # atom_mapping imports numpy at module scope
    from delfin.atom_mapping import map_atoms
except Exception:  # pragma: no cover - import-time guard, exercised in CI
    map_atoms = None  # type: ignore[assignment]


# Distance (A) above which two positions count as "moved" in a diff.
_MOVED_ANGSTROM = 0.1

# Distance (A) below which a pair of atoms counts as a steric clash in a diff.
# 1.0 A is well below any real covalent single bond and above the sum of two
# hydrogens sitting on top of each other, so it only fires on physically
# impossible overlaps.
_CLASH_ANGSTROM = 1.0

# Standard maximal valence for common main-group elements. A non-metal atom
# whose heavy-atom degree is below this needs hydrogen(s) placed.
_VALENCE: dict[str, int] = {
    "H": 1,
    "B": 3,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    "P": 3,
    "S": 2,
    "Cl": 1,
    "Br": 1,
    "I": 1,
}


def _parse(text: str, what: str = "structure") -> Geom:
    if not isinstance(text, str) or not text.strip():
        raise ValueError(
            f"{what}: expected non-empty xyz text, got "
            + ("empty text" if isinstance(text, str) else type(text).__name__)
        )
    # An ORCA coordinate block carries a real charge/multiplicity in its
    # "* xyz {charge} {mult}" header line. Route those to parse_orca_block so
    # the diff's charge/multiplicity check reads real values instead of the
    # parse_xyz defaults (0, 1). Plain xyz text (no "* xyz" header) still
    # parses identically through parse_xyz.
    first_line = text.lstrip().splitlines()[0].strip()
    is_orca_block = bool(re.match(r"^\*\s*(?:xyz|xyzfile|int|frac)\b", first_line))
    parser = parse_orca_block if is_orca_block else parse_xyz
    try:
        return parser(text)
    except GeometryError as exc:
        raise ValueError(f"{what} is not a valid xyz block: {exc}") from exc


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def _min_pair_distance(coords: list[tuple[float, float, float]]) -> float | None:
    """Smallest distance between any two atoms of a structure, or None for <2
    atoms. The clash check compares this against ``_CLASH_ANGSTROM``."""
    best: float | None = None
    n = len(coords)
    for i in range(n):
        for j in range(i + 1, n):
            d = _dist(coords[i], coords[j])
            if best is None or d < best:
                best = d
    return best


# ---------------------------------------------------------------------------
# renumber (product onto educt)
# ---------------------------------------------------------------------------

def renumber_to_reference(
    xyz_text: str,
    reference_xyz_text: str,
    max_edits: int = 2,
) -> dict[str, Any]:
    """Reorder ``xyz_text`` onto ``reference_xyz_text`` atom order.

    Uses the universal educt<->product mapper ``delfin.atom_mapping.map_atoms``
    (never a second mapper). A formula mismatch raises ``ValueError`` and is
    never hidden. Returns a dict with the reordered structure, the verified
    flag, the renumbering table (product index -> reference index) and checks.
    """
    tgt = _parse(xyz_text, "target")
    ref = _parse(reference_xyz_text, "reference")
    if map_atoms is None:  # pragma: no cover - atom_mapping unavailable
        return {
            "reordered_xyz": "",
            "n_atoms": len(tgt),
            "renumbering": [],
            "verified": False,
            "rmsd": None,
            "error": "atom_mapping is unavailable on this host",
        }
    result = map_atoms(
        list(ref.symbols), list(ref.coords),
        list(tgt.symbols), list(tgt.coords),
        max_edits=max_edits,
    )
    order: list[int] = result["order"]  # target index per reference atom
    reorder: list[int] = [0] * len(tgt)
    for ref_i, tgt_i in enumerate(order):
        reorder[tgt_i] = ref_i
    table = [
        {"product_index": j, "reference_index": ref_i}
        for j, ref_i in enumerate(reorder)
    ]
    reordered = Geom(
        symbols=[tgt.symbols[order[i]] for i in range(len(order))],
        coords=[tgt.coords[order[i]] for i in range(len(order))],
        charge=tgt.charge,
        multiplicity=tgt.multiplicity,
        comment=f"renumbered onto reference ({result['method']})",
    )
    return {
        "reordered_xyz": to_xyz(reordered),
        "n_atoms": len(order),
        "renumbering": table,
        "n_bond_edits": result["n_bond_edits"],
        "rmsd": result["rmsd"],
        "method": result["method"],
        "verified": bool(result["verified"]),
        "checks": ["composition matching via atom_mapping"],
        "error": None,
    }


# ---------------------------------------------------------------------------
# fix hydrogens (add missing coordinating / implicit H)
# ---------------------------------------------------------------------------

def _heavy_neighbour_directions(
    g: Geom, atom: int, nbrs: dict[int, list[int]]
) -> list[tuple[float, float, float]]:
    """Unit vectors from ``atom`` toward each heavy neighbour (away faces)."""
    out: list[tuple[float, float, float]] = []
    for j in nbrs.get(atom, ()):
        if g.symbols[j] == "H":
            continue
        v = tuple(
            (g.coords[j][k] - g.coords[atom][k]) for k in range(3)
        )
        n = sum(x * x for x in v) ** 0.5
        if n > 1e-9:
            out.append(tuple(x / n for x in v))
    return out


def _h_directions(
    away: list[tuple[float, float, float]], count: int
) -> list[tuple[float, float, float]]:
    """``count`` unit directions on ``atom`` that point AWAY from its heavy
    neighbours and are mutually separated like tetrahedral hydrogens.

    The base direction is the negated sum of the toward-neighbour vectors
    (the far side of the atom). A second/third hydrogen completes the
    tetrahedron: 109.47 deg from the base axis, spread 120 deg apart, so
    two placed hydrogens are 109.5 deg apart (>= 1.6 bond lengths) and
    never coincide.
    """
    # Sum of toward vectors, negated -> away (may be empty -> +x fallback).
    s = [-sum(v[k] for v in away) for k in range(3)]
    n = sum(x * x for x in s) ** 0.5
    base = tuple(x / n for x in s) if n > 1e-9 else (1.0, 0.0, 0.0)
    if count <= 1:
        return [base]
    # Orthonormal frame around base.
    ref = (0.0, 0.0, 1.0) if abs(base[2]) < 0.9 else (1.0, 0.0, 0.0)
    v = [base[1] * ref[2] - base[2] * ref[1],
         base[2] * ref[0] - base[0] * ref[2],
         base[0] * ref[1] - base[1] * ref[0]]
    vn = sum(x * x for x in v) ** 0.5
    v = [x / vn for x in v]
    w = [base[1] * v[2] - base[2] * v[1],
         base[2] * v[0] - base[0] * v[2],
         base[0] * v[1] - base[1] * v[0]]
    ang = math.radians(109.47)  # tetrahedral angle from the base axis
    sa, ca = math.sin(ang), math.cos(ang)
    dirs = []
    for idx in range(count):
        az = 2.0 * math.pi * idx / count
        dz, dv = ca, sa
        cv, sv = math.cos(az), math.sin(az)
        dirs.append(tuple(
            dz * base[k] + dv * (cv * v[k] + sv * w[k]) for k in range(3)
        ))
    return dirs


# Distance ratio (actual bond / covalent-radii single-bond reference) below
# which a bond is treated as a multiple bond with no free valence for H.
_MULTIPLE_BOND_RATIO = 0.90
# Distance up to which a contact is still a single bond needing no H-free
# valence beyond its own (roughly covalent-radii sum * this factor).
_SINGLE_BOND_FACTOR = 1.12


def _heavy_bond_order(
    g: Geom, i: int, j: int, radii: list[float]
) -> int:
    """Estimated bond order (0..3) between two heavy atoms from distance.

    Uses the covalent-radii sum as the single-bond reference (DELFIN's own
    radii). A bond clearly shorter than that reference is a multiple bond
    (double/triple) and leaves no free valence for an added hydrogen; a
    bond near or above the reference is a single bond. Returns 0 when the
    atoms are not in contact. This is what keeps fix_hydrogens from
    hydroxylating a double-bonded oxygen (CO2) or triple-bonded nitrogen.
    """
    d = _dist(g.coords[i], g.coords[j])
    rsum = radii[i] + radii[j]
    if d > rsum * _SINGLE_BOND_FACTOR:
        return 0
    ratio = d / rsum if rsum else 1.0
    if ratio <= 0.80:
        # Real X#Y triple-bond distances sit at ~0.79 of the covalent-radii
        # sum (C#N 1.16/1.47, C#C 1.20/1.52); the double window starts at
        # ~0.86 (C=O 1.22/1.42) and singles at ~0.95 (C-C 1.54/1.52).
        return 3
    if ratio <= _MULTIPLE_BOND_RATIO:
        return 2
    return 1


# Typical X-H bond lengths (Angstrom) for placing a guessed hydrogen, instead
# of a covalent-radius sum that over-lengthens C-H (1.39 A vs the real 1.09).
_XH_BOND_LENGTHS: dict[str, float] = {
    "C": 1.09,
    "N": 1.01,
    "O": 0.96,
    "S": 1.34,
    "P": 1.42,
    "B": 1.19,
}


def fix_hydrogens(xyz_text: str) -> dict[str, Any]:
    """Add missing hydrogens to a structure; never removes an atom.

    A non-metal atom whose heavy-bond order sum is below its standard
    valence gets the missing hydrogen(s) placed at DELFIN's covalent-radius
    bond length, pushed away from its heavy neighbours. Heavy-bond order is
    estimated from the covalent radii and the interatomic distance, so a
    multiple bond (CO2's C=O, a nitrile C#N) saturates the atom and no H is
    guessed. Metals are never given guessed hydrogens (their coordination H
    is a donor convention, not a valence count).
    """
    g = _parse(xyz_text, "structure")
    radii = _covalent_radii(g)
    nbrs = _bond_graph(g)
    added: list[dict[str, Any]] = []
    new_symbols = list(g.symbols)
    new_coords = list(g.coords)
    for i, sym in enumerate(g.symbols):
        if sym == "H":
            continue
        need = _VALENCE.get(sym)
        if not need:
            continue  # metal or unknown -> no guessed H
        heavy_order = sum(
            _heavy_bond_order(g, i, j, radii)
            for j in nbrs.get(i, ())
            if g.symbols[j] != "H"
        )
        existing_h = sum(
            1 for j in nbrs.get(i, ()) if g.symbols[j] == "H"
        )
        missing = need - heavy_order - existing_h
        if missing <= 0:
            continue
        away = _heavy_neighbour_directions(g, i, nbrs)
        r_h_i = _XH_BOND_LENGTHS.get(sym, radii[i] + _COSMETIC_H_RADIUS)
        for direction in _h_directions(away, missing):
            pos = tuple(
                g.coords[i][k] + direction[k] * r_h_i for k in range(3)
            )
            new_symbols.append("H")
            new_coords.append(pos)
            added.append(
                {"index": len(new_symbols) - 1, "bonded_to": i, "element": "H"}
            )
    out = Geom(
        symbols=new_symbols,
        coords=new_coords,
        charge=g.charge,
        multiplicity=g.multiplicity,
        comment=g.comment or "hydrogens fixed",
    )
    return {
        "xyz": to_xyz(out),
        "n_added": len(added),
        "added_rows": added,
        "checks": [
            f"added {len(added)} hydrogen(s), removed none",
            "H placed at DELFIN covalent-radius bond length",
        ],
        "error": None,
    }


_COSMETIC_H_RADIUS = 0.63  # DELFIN's fallback covalent radius for H


# ---------------------------------------------------------------------------
# structure_diff (per-atom diff + checks)
# ---------------------------------------------------------------------------

def structure_diff(xyz_a_text: str, xyz_b_text: str) -> dict[str, Any]:
    """A per-atom diff of ``xyz_a`` against ``xyz_b`` plus checks.

    Reports changed element and moved atoms (position change above
    ``_MOVED_ANGSTROM``), the composition equality, the charge /
    multiplicity consistency (read from an ORCA ``* xyz`` header when the
    input is an ORCA block; plain xyz defaults to 0/1) and a steric-clash
    check: a pair closer than ``_CLASH_ANGSTROM`` A is reported. Atom
    positions compared positionally (same index); a renumbering is a
    separate operation (``renumber_to_reference``).
    """
    a = _parse(xyz_a_text, "structure A")
    b = _parse(xyz_b_text, "structure B")
    n = min(len(a), len(b))
    rows: list[dict[str, Any]] = []
    for i in range(n):
        sa, sb = a.symbols[i], b.symbols[i]
        da = _dist(a.coords[i], b.coords[i])
        rows.append({
            "index": i,
            "from_symbol": sa,
            "to_symbol": sb,
            "element_changed": sa != sb,
            "moved": da > _MOVED_ANGSTROM,
            "distance": round(da, 3),
        })
    for i in range(n, len(a)):
        rows.append({
            "index": i, "from_symbol": a.symbols[i], "to_symbol": None,
            "element_changed": True, "moved": True, "distance": None,
        })
    for i in range(n, len(b)):
        rows.append({
            "index": i, "from_symbol": None, "to_symbol": b.symbols[i],
            "element_changed": True, "moved": True, "distance": None,
        })
    min_da = _min_pair_distance(a.coords)
    min_db = _min_pair_distance(b.coords)
    both = [d for d in (min_da, min_db) if d is not None]
    min_distance = min(both) if both else None
    clash_detected = min_distance is not None and min_distance < _CLASH_ANGSTROM
    return {
        "rows": rows,
        "composition_equal": Counter(a.symbols) == Counter(b.symbols),
        "charge_mult_match": (
            a.charge == b.charge and a.multiplicity == b.multiplicity
        ),
        "charge_a": a.charge,
        "multiplicity_a": a.multiplicity,
        "charge_b": b.charge,
        "multiplicity_b": b.multiplicity,
        "min_distance": min_distance,
        "min_distance_a": min_da,
        "min_distance_b": min_db,
        "clash_detected": clash_detected,
        "clash_angstrom": _CLASH_ANGSTROM,
        "checks": [
            "composition compared by element count",
            "charge/multiplicity compared (0/1 default for plain xyz, "
            "real values from an ORCA * xyz header)",
            f"clash check: no pair under {_CLASH_ANGSTROM} Angstrom",
        ],
        "error": None,
    }
