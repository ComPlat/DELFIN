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
from delfin.xyz_io import _elem_from_label


def _xyz_coord_like(text: str) -> bool:
    """True when ``text`` parses as a float (an xyz coordinate field)."""
    try:
        float(text)
        return True
    except ValueError:
        return False

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


# Real X-X-H bond angle (degrees) for an atom with ONE heavy neighbour
# whose H placement is sp3-like but element-specific — linear 180 is
# wrong for these (operator round-7 rerun: methanol C-O-H 180.0,
# real ~108.5; methanethiol C-S-H 180.0, real ~96). Absent entries
# (C, N) keep the sp3/sp2 machinery.
_BOND_ANGLE_ONE_NEIGHBOUR: dict[str, float] = {
    "O": 108.5,
    "S": 96.0,
}


def _filter_graph(nbrs, g, radii):
    """Drop contacts above 1.05 x the covalent-radii sum. Real covalent
    bonds sit at <= ~1.01 x (C-C 1.54/1.52, C-S 1.82/1.81; multiple and
    aromatic bonds well below); ``_bond_graph``'s 0.45-tolerance cutoff
    is purely geometric and admits phantom bonds: a benzene C-C stretched
    to 1.75 A passed, giving that carbon 2 H instead of 1 (probe R7-F1
    v2: 9 H on benzene)."""
    out: dict[int, list[int]] = {}
    for i, js in nbrs.items():
        keep = []
        for j in js:
            d = math.dist(g.coords[i], g.coords[j])
            if d <= 1.05 * (radii[i] + radii[j]):
                keep.append(j)
        if keep:
            out[i] = keep
    return out


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
    away: list[tuple[float, float, float]], count: int,
    heavy_order: int = 1,
    polar_override: float | None = None,
) -> list[tuple[float, float, float]]:
    """``count`` unit directions for hydrogens on ``atom`` so that the
    placed angles match real hybridisation geometry, measured AT THE
    VERTEX (operator real-run round 4):

    - lone atom (no heavy neighbours): a regular tetrahedron for 4 H
      (CH4, H-C-H 109.47), a 104.5-deg pair for 2 H (H2O), a 106.7-deg
      ring for 3 H (NH3), one direction for 1 H.
    - one heavy neighbour: H sit 109.47 deg from the bond (sp3; the
      previous 70.53-deg cone with opposite azimuths produced
      H-C-H 141.1), or 120 deg from a double bond (order 2 via
      ``half_deg`` override) -> C-X-H ~121.5 for sp2 centres.
    - two heavy neighbours (sp2, e.g. =CH2): H in the bond plane at
      121.5 deg from EACH bond -> C-C-H 121.5, H-C-H 117, coplanar
      (the previous cone put them out of plane).
    """
    if not away:
        base = (0.0, 0.0, 1.0)
    elif len(away) == 1:
        # One heavy neighbour. The sp2 case is detected by the heavy
        # bond ORDER (=CH2 sits on a double bond, order 2), NOT by the
        # missing-H count: an amine N (valence 3, single bond, 2
        # missing H) has the same count but is sp3/pyramidal, and the
        # count-based detection made methylamine's NH2 planar 120/120
        # (operator real-run round 4 rerun).
        base = tuple(-x for x in away[0])
        n = sum(x * x for x in base) ** 0.5
        base = tuple(x / n for x in base) if n > 1e-9 else (0.0, 0.0, 1.0)
        # sp2 CH2= centre: two H at 120 deg from the bond, 120 apart
        # (real Csp2-H2: C-C-H ~121.5 with small pyramidalisation; 120 is
        # the planar ideal and keeps H-C-H 120 too).
        if heavy_order >= 2 and count == 2:
            # =CH2-type sp2 centre: two H at 120 deg from the double bond,
            # 120 apart (planar ideal; real C-C-H ~121.5 carries small
            # pyramidalisation).
            return _ring(base, 60.0, 2, offset=180.0)
        if heavy_order >= 2 and count == 1:
            # sp carbon on a triple bond (order 3, e.g. acetylene): its
            # ONE H continues the bond axis -> X-C-H 180 (sp). The imine
            # ring at 120 deg trans-bent it instead (operator round 5
            # rerun). An imine N keeps the 120-deg in-plane H (order 2).
            if heavy_order == 3:
                return [base]
            return _ring(base, 60.0, 2, offset=180.0)[:1]
        if heavy_order == 1 and polar_override is not None and count == 1:
            # O/S with one heavy neighbour is sp3-like: the H sits at the
            # element's real bond angle from the parent bond (O ~108.5,
            # S ~96), NOT linear (operator round 7 rerun). polar_override
            # is the X-X-H angle; the away-axis polar angle is 180-angle.
            return _ring(base, 180.0 - polar_override, 1)
        return _tetrahedral_from_axis(base, count)
    else:
        # Two heavy neighbours. For ONE H the operator's rule (round 5
        # rerun) is universal and right: direction = -(sum of unit
        # vectors to the neighbours), normalised — the H points AWAY
        # from all of them (formamide's C-H sat 0.265 A from N in the
        # 121.5-from-bond-1 placement). For TWO H (=CH2-type sp2
        # centres) the in-plane mirrored 120-deg placement stays.
        if count == 1:
            s = [0.0, 0.0, 0.0]
            for a in away:
                u = _unit(a)
                for k in range(3):
                    s[k] -= u[k]
            return [_unit(tuple(s))]
        # =CH2-type sp2: H in the plane spanned by the two bonds,
        # mirrored across the C-C axis, so X-X-H and H-X-H both come
        # out at real sp2 values (ethylene: 121.5 / 117).
        b1 = _unit(away[0])
        b2 = _unit(away[1])
        # Internal bisector (points between the bonds) and the angle
        # between the bonds.
        dot = max(-1.0, min(1.0, sum(b1[k] * b2[k] for k in range(3))))
        theta = math.degrees(math.acos(dot))  # e.g. 0 for ethylene
        # Half the heavy-heavy angle; H must sit at 180 - half - 121.5
        # from the bisector measured toward each bond side... For ethylene
        # the bonds are collinear (theta=0): H at 121.5 from each bond on
        # opposite sides of the C-C axis, in-plane.
        out: list[tuple[float, float, float]] = []
        # In-plane frame: axis = b1 (toward bond 1), perp = component of
        # b2 orthogonal to b1 (zero for collinear -> pick any perp).
        perp = [b2[k] - dot * b1[k] for k in range(3)]
        pn = sum(x * x for x in perp) ** 0.5
        if pn > 1e-9:
            perp = [x / pn for x in perp]
        else:
            # Collinear bonds (ethylene): any plane contains the axis;
            # choose one deterministic perp.
            perp = _any_perp(b1)
        for idx in range(count):
            side = 1.0 if idx % 2 == 0 else -1.0
            # Angle from bond-1 axis toward +/-perp: 121.5 deg.
            a1 = math.radians(121.5)
            ca, sa = math.cos(a1), math.sin(a1)
            d = [ca * b1[k] + side * sa * perp[k] for k in range(3)]
            out.append(tuple(d))
        return out
    # Lone atom: exact small-n geometries around +z.
    if count == 1:
        return [base]
    if count == 2:
        return _ring(base, 52.25, 2, offset=180.0)  # H-O-H 104.5
    if count == 3:
        return _ring(base, 67.98, 3)  # H-N-H 106.7 ring
    return _tetrahedron()  # CH4: 109.47 / H-H 1.78


def _unit(v):
    n = sum(x * x for x in v) ** 0.5
    return tuple(x / n for x in v) if n > 1e-9 else v


def _any_perp(v):
    ref = (0.0, 0.0, 1.0) if abs(v[2]) < 0.9 else (1.0, 0.0, 0.0)
    w = [v[1] * ref[2] - v[2] * ref[1],
         v[2] * ref[0] - v[0] * ref[2],
         v[0] * ref[1] - v[1] * ref[0]]
    return _unit(w)


def _ring(base, half_deg, count, offset=0.0):
    """``count`` points on a cone of half-angle ``half_deg`` around base."""
    a = math.radians(half_deg)
    sa, ca = math.sin(a), math.cos(a)
    v = _any_perp(base)
    w = [base[1] * v[2] - base[2] * v[1],
         base[2] * v[0] - base[0] * v[2],
         base[0] * v[1] - base[1] * v[0]]
    dirs = []
    for idx in range(count):
        az = 2.0 * math.pi * idx / count + math.radians(offset)
        cv, sv = math.cos(az), math.sin(az)
        dirs.append(tuple(
            ca * base[k] + sa * (cv * v[k] + sv * w[k]) for k in range(3)
        ))
    return dirs


def _tetrahedron():
    """Regular tetrahedron around the origin (pairwise 109.47 deg).

    Corners are NORMALISED: raw (+-1,+-1,+-1) corners have magnitude
    sqrt(3), and the caller multiplies by the X-H bond length, which gave
    C-H 1.09*sqrt(3) = 1.888 A (operator real-run round 4).
    """
    return [tuple(x / math.sqrt(3.0) for x in corner) for corner in (
        (1.0, 1.0, 1.0),
        (1.0, -1.0, -1.0),
        (-1.0, 1.0, -1.0),
        (-1.0, -1.0, 1.0),
    )]


def _tetrahedral_from_axis(base, count):
    """sp3 H directions: each 109.47 deg from the bond axis (-base)."""
    if count == 1:
        return [base]
    if count == 2:
        # Two of the three 120-deg-apart ring points: each at 70.53 from
        # the away axis = 109.47 from the bond, and the pair angle
        # between two ring points 120 deg apart in azimuth is 109.47.
        return _ring(base, 70.53, 3)[:2]
    # 3 H (CH3): classic 109.47 staggered ring: each at 70.53 from base
    # (= 109.47 from the bond), azimuths 120 apart -> pair angles 109.47.
    return _ring(base, 70.53, 3)


# Distance ratio (actual bond / covalent-radii single-bond reference) below
# which a bond is treated as a multiple bond with no free valence for H.
_MULTIPLE_BOND_RATIO = 0.90
# Distance up to which a contact is still a single bond needing no H-free
# valence beyond its own (roughly covalent-radii sum * this factor).
_SINGLE_BOND_FACTOR = 1.12


def _heavy_bond_order(
    g: Geom, i: int, j: int, radii: list[float]
) -> float:
    """Estimated bond order (0..3, or 1.5 for aromatic) between heavy atoms.

    Uses the covalent-radii sum as the single-bond reference (DELFIN's own
    radii). A bond clearly shorter than that reference is a multiple bond
    (double/triple) and leaves no free valence for an added hydrogen; a
    bond near or above the reference is a single bond. Returns 0 when the
    atoms are not in contact. This is what keeps fix_hydrogens from
    hydroxylating a double-bonded oxygen (CO2) or triple-bonded nitrogen.

    Aromatic bonds (operator round 7): a C-C pair at 0.895-0.945 of the
    radii sum in which BOTH partners have at least two heavy neighbours
    (ring context) is order 1.5, so a benzene ring C carries valence
    1.5+1.5+1 = 4 -> exactly one H, and a pyridine N carries
    1.5+1.5 = 3 -> no H.
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
    # Aromatic window (operator round 7): benzene C-C 1.39/1.52 = 0.914,
    # C-N aromatic ~1.34/1.42 = 0.943. Both partners have >= 2 heavy
    # contacts INCLUDING each other (ring context): a benzene ring C has
    # exactly its two ring partners in contact. An isolated single C-C
    # at the same ratio would have its partner + at most substituents —
    # the distance window is the primary discriminator.
    if 0.89 <= ratio <= 0.95 and g.symbols[i] in _SP2_ELEMENTS:
        ring_i = sum(
            1 for k in range(len(g.symbols))
            if k != i
            and g.symbols[k] != "H"
            and _dist(g.coords[i], g.coords[k]) <= (
                radii[i] + radii[k]) * _SINGLE_BOND_FACTOR
        ) >= 2
        ring_j = sum(
            1 for k in range(len(g.symbols))
            if k != j
            and g.symbols[k] != "H"
            and _dist(g.coords[j], g.coords[k]) <= (
                radii[j] + radii[k]) * _SINGLE_BOND_FACTOR
        ) >= 2
        if ring_i and ring_j:
            return 1.5
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

# Elements that can take an aromatic (order 1.5) bond in the ring window.
_SP2_ELEMENTS = {"C", "N"}


def _on_aromatic_cycle(g, i, nbrs, radii) -> bool:
    """True when heavy atom i lies on a ring cycle that carries an
    aromatic (order 1.5) bond.

    Detection: with exactly two heavy neighbours a cycle through i
    exists iff a and b stay connected without passing through i (BFS);
    a carbon with 0, 1 or 3+ heavy neighbours is either not on a cycle
    this way or is a fused junction with no free valence. Only cycle
    edges are tested for order 1.5, so a CH2 bridge between two
    separate aromatic rings (diphenylmethane) is NOT clamped -- its
    neighbours are in different rings and do not connect around it.
    """
    heavy = [j for j in nbrs.get(i, ()) if g.symbols[j] != "H"]
    if len(heavy) != 2:
        return False
    a, b = heavy
    seen = {i, a}
    prev: dict[int, int | None] = {a: None}
    queue = [a]
    while queue:
        cur = queue.pop(0)
        for j in nbrs.get(cur, ()):
            if j in seen or g.symbols[j] == "H":
                continue
            seen.add(j)
            prev[j] = cur
            queue.append(j)
    if b not in prev:
        return False
    edges = [(i, a), (i, b)]
    cur = b
    while prev[cur] is not None:
        edges.append((prev[cur], cur))
        cur = prev[cur]
    return any(
        _heavy_bond_order(g, x, y, radii) == 1.5 for x, y in edges
    )


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
    # parse_xyz strips blank lines BEFORE reading the comment header, so a
    # blank comment line makes it read the first atom line as the comment
    # and the written block repeats the atom line as its header (operator
    # real-run round 3). Contained here: an empty comment, or one whose
    # first word parses as an element followed by coordinates (an atom
    # line), is replaced by a neutral header. Real text comments survive.
    parts = (g.comment or "").split()
    atom_like = bool(
        parts
        and len(parts) >= 4
        and _elem_from_label(parts[0]) is not None
        and all(_xyz_coord_like(p) for p in parts[1:4])
    )
    if not (g.comment or "").strip() or atom_like:
        g.comment = "hydrogens fixed"
    radii = _covalent_radii(g)
    nbrs = _filter_graph(_bond_graph(g), g, radii)
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
        # Aromatic / partial bond orders make the sum fractional (a
        # benzene carbon: 3 x 1.5 = 4.5, missing 1.5). The hydrogen
        # COUNT is an integer: round to nearest, minimum 1 when the
        # fractional part says at least one H is needed (reviewer
        # R7-F1: count=1.5 crashed _h_directions with a float slice).
        # QS-XE1-F1: on a ring cycle that carries an aromatic bond the
        # count is clamped to 1 -- a fractional valence sum there means
        # "at least one", never "round up to a CH2", which would make
        # the carbon sp3 and break the pi sextet. Non-ring contexts
        # keep int(round()).
        if _on_aromatic_cycle(g, i, nbrs, radii):
            # Cap at 1; zero or negative passes through (pyridine N:
            # 1.5+1.5 = 3.0, missing 0 -> no H).
            missing = min(1, int(round(missing)))
        else:
            missing = int(round(missing))
        if missing <= 0:
            continue
        away = _heavy_neighbour_directions(g, i, nbrs)
        # Planar amide N (operator round 7): an N whose ONE heavy
        # neighbour is a carbon that itself carries a double bond
        # (N-C(=O), e.g. formamide) is planar — its H pair sits at
        # 120 deg in-plane, not pyramidal 109.5. Implemented as an
        # sp2 hint for the DIRECTION choice only (missing was already
        # computed from the real orders).
        dir_order = heavy_order
        if (
            sym == "N"
            and sum(1 for j in nbrs.get(i, ()) if g.symbols[j] != "H") == 1
        ):
            n_carbon = next(
                (j for j in nbrs.get(i, ()) if g.symbols[j] != "H"), None
            )
            if n_carbon is not None and g.symbols[n_carbon] == "C":
                carbon_order = sum(
                    _heavy_bond_order(g, n_carbon, k, radii)
                    for k in nbrs.get(n_carbon, ())
                    if k != i and g.symbols[k] != "H"
                )
                if carbon_order >= 2:
                    dir_order = 2
        r_h_i = _XH_BOND_LENGTHS.get(sym, radii[i] + _COSMETIC_H_RADIUS)
        # Element-real bond angle for an O/S with one heavy neighbour
        # (operator round-7 rerun: linear 180 was wrong; O ~108.5,
        # S ~96). C and N keep the sp3/sp2 machinery.
        polar_override = _BOND_ANGLE_ONE_NEIGHBOUR.get(sym)
        for direction in _h_directions(away, missing, dir_order,
                                       polar_override=polar_override):
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
    checks = [
        f"added {len(added)} hydrogen(s), removed none",
        "H placed at real X-H bond lengths (C-H 1.09, N-H 1.01, O-H 0.96 A)",
    ]
    # Self-check: flag (do not silently ship) placed H that clash with
    # anything already in the structure (min interatomic pair distance).
    clash = _clash_report(new_symbols, new_coords)
    if clash:
        checks.append("WARNING: " + clash)
    return {
        "xyz": to_xyz(out),
        "n_added": len(added),
        "added_rows": added,
        "checks": checks,
        "error": None,
    }


# Hydrogen's parent (the heavy atom it is bonded to): the NEAREST heavy
# atom. In a fixed-up structure a placed H sits at its real X-H bond
# length (~0.96-1.09 A) while any other heavy atom is farther away, so
# "nearest heavy" identifies the parent without a bond table.
def _parent_index(symbols, coords, h_index: int) -> int | None:
    parent = None
    best = None
    for b, sym in enumerate(symbols):
        if sym == "H" or b == h_index:
            continue
        d = math.dist(coords[h_index], coords[b])
        if best is None or d < best:
            best = d
            parent = b
    return parent


def _covalent_radii_from_symbols(symbols) -> list[float]:
    """Covalent radii (A) per element symbol, from DELFIN's RDKit-backed
    table (``_covalent_radii``) with the same fallbacks, so the clash
    check can run on a plain (symbol, coordinate) pair without a Geom."""
    from delfin.tools.geometry_edit import Geom as _G

    g = _G(symbols=list(symbols), coords=[(0.0, 0.0, 0.0)] * len(symbols),
           charge=0, multiplicity=1, comment="radii only")
    return _covalent_radii(g)


def _clash_report(symbols, coords) -> str | None:
    """First clash between a placed H and the rest of the structure.

    Rules (operator rounds 4 and 7 rerun):
    - H-H pair below 1.4 A (two hydrogens in a real structure are never
      closer than ~1.4 A even gauche/ortho);
    - any pair below 0.8 A (physical overlap);
    - a placed H closer than 0.9 x the covalent-radii sum to a HEAVY
      atom that is NOT its parent (ethanol round-7 rerun: a CH2 hydrogen
      landed 0.42 A from the O and the check stayed silent). The parent
      is excluded because its X-H bond length (e.g. C-H 1.09 < 0.9 x
      1.07) is by definition a bond, not a clash.
    """
    limit_same = 1.4  # two H cannot be this close in a real structure
    limit_any = 0.8
    radii = _covalent_radii_from_symbols(symbols)
    n = len(symbols)
    for a in range(n):
        if symbols[a] != "H":
            continue
        parent = _parent_index(symbols, coords, a)
        for b in range(n):
            if b == a:
                continue
            d = math.dist(coords[a], coords[b])
            if symbols[b] == "H" and d < limit_same:
                return (
                    f"H-H clash {d:.3f} A (< {limit_same} A) between "
                    f"atoms {a} and {b}"
                )
            if d < limit_any:
                return (
                    f"clash {d:.3f} A (< {limit_any} A) between atoms "
                    f"{a} ({symbols[a]}) and {b} ({symbols[b]})"
                )
            # H..heavy to a NON-parent heavy atom below 0.9 x cov sum.
            if symbols[b] != "H" and parent is not None and b != parent:
                cov = 0.9 * (radii[a] + radii[b])
                if d < cov:
                    return (
                        f"H..{symbols[b]} clash {d:.3f} A "
                        f"(< 0.9 x covalent sum {cov:.3f} A) between "
                        f"atoms {a} and {b}"
                    )
    return None


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
