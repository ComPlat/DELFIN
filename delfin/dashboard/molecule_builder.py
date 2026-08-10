"""Structural edits to a coordinate block: place, delete and retype atoms.

Everything here works on an XYZ block and hands one back, because that is what
the Submit tab keeps, what the viewer renders and what the force field is
assigned from.  Nothing in this module knows about widgets.

Why the chemistry lives on the Python side at all: placing an atom is cheap,
but deciding how many hydrogens it needs and where they go is not something a
browser can be trusted with.  RDKit's periodic table has the valences and the
covalent radii; the geometry of a free valence follows from the directions
already taken at that centre.

Hydrogens are adjusted for the atoms an edit touches, and only those.  Filling
every unsatisfied valence in the molecule would quietly hydrogenate a radical
or a coordinating donor the user put there on purpose.
"""
from __future__ import annotations

import logging
import math
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    RDKIT_AVAILABLE = True
except Exception:  # pragma: no cover - exercised only without RDKit
    Chem = None
    RDKIT_AVAILABLE = False


#: Elements offered in the editor, in the order chemists reach for them.
#: Metals are included because this dashboard exists for coordination
#: chemistry; they take no hydrogens (see :data:`_NO_HYDROGENS`).
DRAW_ELEMENTS: Tuple[str, ...] = (
    'C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si',
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ru', 'Rh', 'Pd', 'Ag', 'Ir', 'Pt', 'Au',
)

#: Elements whose valence the editor does not try to satisfy.  A metal's
#: coordination number is the user's decision, not a table's, and filling it
#: with hydrides would be wrong in almost every case.
_NO_HYDROGENS = frozenset({
    'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ru', 'Rh', 'Pd', 'Ag', 'Ir', 'Pt', 'Au',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Cd', 'Hf',
    'Ta', 'W', 'Re', 'Os', 'Hg', 'La', 'Ce', 'Eu', 'Gd', 'Yb', 'Lu', 'U',
})

#: Fallback valences for when RDKit is unavailable or has no default.
_VALENCE = {
    'H': 1, 'B': 3, 'C': 4, 'N': 3, 'O': 2, 'F': 1,
    'Si': 4, 'P': 3, 'S': 2, 'Cl': 1, 'Br': 1, 'I': 1,
}

#: Lone pairs on a neutral atom of each element.  The shape at a centre is
#: set by its steric number -- the number of *partners* plus the lone pairs --
#: and not by the sum of its bond orders: the carbons of ethene have three
#: partners and no lone pair, so they are trigonal planar, while the double
#: bond counts once.  Using the valence sum instead made ethene tetrahedral.
_LONE_PAIRS = {
    'H': 0, 'B': 0, 'C': 0, 'Si': 0,
    'N': 1, 'P': 1,
    'O': 2, 'S': 2,
    'F': 3, 'Cl': 3, 'Br': 3, 'I': 3,
}

#: Covalent radii in Angstrom, used when RDKit cannot supply one.
_RADII = {
    'H': 0.31, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57,
    'Si': 1.11, 'P': 1.07, 'S': 1.05, 'Cl': 1.02, 'Br': 1.20, 'I': 1.39,
}

_DEFAULT_RADIUS = 1.30


def normalise_element(raw: Any) -> Optional[str]:
    """Return a capitalised element symbol, or None if it is not one."""
    text = str(raw or '').strip()
    if not text or not text.isalpha() or len(text) > 2:
        return None
    symbol = text[0].upper() + text[1:].lower()
    if RDKIT_AVAILABLE:
        # Quietly: RDKit's periodic table does not merely raise for a symbol
        # it does not know, it prints a post-condition violation and a C++
        # stack trace to stderr first, which in a notebook lands in the user's
        # face for what is only a typo in a dropdown.
        from .molecule_forcefield import _RDKitQuiet

        try:
            with _RDKitQuiet():
                if Chem.GetPeriodicTable().GetAtomicNumber(symbol) > 0:
                    return symbol
        except Exception:
            return None
        return None
    return symbol if symbol in _VALENCE or symbol in _NO_HYDROGENS else None


def covalent_radius(symbol: str) -> float:
    """Covalent radius in Angstrom."""
    if RDKIT_AVAILABLE:
        try:
            table = Chem.GetPeriodicTable()
            value = float(table.GetRcovalent(table.GetAtomicNumber(symbol)))
            if value > 0.1:
                return value
        except Exception:
            pass
    return _RADII.get(symbol, _DEFAULT_RADIUS)


def bond_length(first: str, second: str, order: int = 1) -> float:
    """Where a bond between two elements wants to sit.

    The covalent radii are for single bonds; a double bond is shorter by
    about 0.1 A and a triple by about 0.2, which is close enough for placing
    an atom the force field will relax anyway.
    """
    base = covalent_radius(first) + covalent_radius(second)
    return max(0.7, base - 0.10 * max(0, int(order) - 1))


def default_valence(symbol: str) -> Optional[int]:
    """How many bonds this element wants, or None when it is not our business."""
    if symbol in _NO_HYDROGENS:
        return None
    if RDKIT_AVAILABLE:
        try:
            table = Chem.GetPeriodicTable()
            value = int(table.GetDefaultValence(table.GetAtomicNumber(symbol)))
            if value > 0:
                return value
        except Exception:
            pass
    return _VALENCE.get(symbol)


# --------------------------------------------------------------------------
# vector helpers -- deliberately plain, so this module needs no numpy
# --------------------------------------------------------------------------

def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _scale(a, k):
    return (a[0] * k, a[1] * k, a[2] * k)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _norm(a):
    return math.sqrt(_dot(a, a))


def _unit(a, fallback=(0.0, 0.0, 1.0)):
    length = _norm(a)
    if length < 1e-9:
        return fallback
    return _scale(a, 1.0 / length)


def _perpendicular(a):
    """Any unit vector at right angles to ``a``."""
    trial = (0.0, 0.0, 1.0) if abs(a[2]) < 0.9 else (1.0, 0.0, 0.0)
    return _unit(_cross(a, trial))


def free_directions(taken: Sequence[Sequence[float]], count: int,
                    total: int) -> List[Tuple[float, float, float]]:
    """Unit vectors completing a coordination sphere.

    Args:
        taken: The directions already used at this centre, as unit vectors.
        count: How many more are needed.
        total: The coordination number the centre will end up with, which is
            what decides the shape -- 2 linear, 3 trigonal, 4 tetrahedral.

    The construction is geometric rather than tabulated: with nothing taken it
    lays out the ideal shape, and with some directions taken it puts the rest
    as far from them as it can.  Either way the force field tidies up
    afterwards; this only has to be close enough not to look wrong and not to
    start inside another atom.
    """
    taken = [_unit(v) for v in taken]
    out: List[Tuple[float, float, float]] = []
    angle = {1: math.pi, 2: math.pi, 3: math.radians(120.0),
             4: math.radians(109.471)}.get(total, math.radians(109.471))

    for _ in range(max(0, count)):
        here = taken + out
        if not here:
            out.append((0.0, 0.0, 1.0))
            continue
        if len(here) == 1:
            axis = here[0]
            if total <= 2:
                out.append(_scale(axis, -1.0))
                continue
            # Tilt away by the ideal angle, in an arbitrary but fixed plane.
            side = _perpendicular(axis)
            out.append(_unit(_add(_scale(axis, math.cos(angle)),
                                  _scale(side, math.sin(angle)))))
            continue
        if len(here) == 2:
            bisector = _unit(_add(here[0], here[1]), fallback=_perpendicular(here[0]))
            if total <= 3:
                out.append(_scale(bisector, -1.0))
                continue
            normal = _cross(here[0], here[1])
            if _norm(normal) < 1e-6:
                normal = _perpendicular(here[0])
            normal = _unit(normal)
            # Half of the tetrahedral pair straddling the bisector's opposite.
            back = _scale(bisector, -1.0)
            out.append(_unit(_add(_scale(back, math.cos(math.radians(54.75))),
                                  _scale(normal, math.sin(math.radians(54.75))))))
            continue
        # Three or more: opposite the sum, which is where the gap is.
        total_vec = (sum(v[0] for v in here),
                     sum(v[1] for v in here),
                     sum(v[2] for v in here))
        if _norm(total_vec) < 1e-6:
            out.append(_perpendicular(here[0]))
        else:
            out.append(_scale(_unit(total_vec), -1.0))
    return out


# --------------------------------------------------------------------------
# the edits
# --------------------------------------------------------------------------

class Structure:
    """Symbols, coordinates and the bonds an edit has to keep consistent.

    Bond orders are carried alongside because an XYZ block cannot hold them
    and the force field needs them: a bond drawn as a double bond has to stay
    one through the next perception, which reads orders off the geometry and
    would call a long C=C single.
    """

    def __init__(self, symbols: Sequence[str],
                 coords: Sequence[Sequence[float]],
                 bonds: Optional[Dict[Tuple[int, int], int]] = None):
        self.symbols: List[str] = [str(s) for s in symbols]
        self.coords: List[Tuple[float, float, float]] = [
            (float(c[0]), float(c[1]), float(c[2])) for c in coords
        ]
        self.bonds: Dict[Tuple[int, int], int] = dict(bonds or {})
        # Where each atom came from, so anything the tab holds by index --
        # a constraint, a forced type, the metal a polyhedron is on -- can
        # follow an edit that renumbers the structure. -1 means the atom was
        # created here and had no index before.
        self.origin: List[int] = list(range(len(self.symbols)))

    # -- queries ---------------------------------------------------------
    def __len__(self) -> int:
        return len(self.symbols)

    def neighbours(self, index: int) -> List[int]:
        return sorted(
            (j if i == index else i)
            for i, j in self.bonds
            if index in (i, j)
        )

    def order(self, first: int, second: int) -> int:
        return int(self.bonds.get((min(first, second), max(first, second)), 0))

    def used_valence(self, index: int) -> int:
        return sum(
            int(order) for pair, order in self.bonds.items() if index in pair
        )

    def hydrogens_on(self, index: int) -> List[int]:
        return [j for j in self.neighbours(index) if self.symbols[j] == 'H']

    # -- edits -----------------------------------------------------------
    def add_atom(self, element: str, position: Sequence[float]) -> int:
        self.symbols.append(element)
        self.coords.append((float(position[0]), float(position[1]),
                            float(position[2])))
        self.origin.append(-1)
        return len(self.symbols) - 1

    def set_bond(self, first: int, second: int, order: int) -> None:
        if first == second:
            return
        key = (min(first, second), max(first, second))
        if order <= 0:
            self.bonds.pop(key, None)
        else:
            self.bonds[key] = int(order)

    def remove_atoms(self, indices: Iterable[int]) -> Dict[int, int]:
        """Delete atoms and return the old-index to new-index map."""
        drop = {int(i) for i in indices if 0 <= int(i) < len(self.symbols)}
        if not drop:
            return {i: i for i in range(len(self.symbols))}
        mapping: Dict[int, int] = {}
        symbols, coords, origin = [], [], []
        for old in range(len(self.symbols)):
            if old in drop:
                continue
            mapping[old] = len(symbols)
            symbols.append(self.symbols[old])
            coords.append(self.coords[old])
            origin.append(self.origin[old])
        bonds = {}
        for (i, j), order in self.bonds.items():
            if i in drop or j in drop:
                continue
            bonds[(mapping[i], mapping[j])] = order
        self.symbols, self.coords, self.bonds = symbols, coords, bonds
        self.origin = origin
        return mapping

    def renumbering(self) -> Dict[int, int]:
        """Where each atom of the structure this was built from ended up."""
        return {
            was: now for now, was in enumerate(self.origin) if was >= 0
        }


def structure_from_xyz(xyz_text: str,
                       bond_orders: Optional[Dict[Tuple[int, int], int]] = None
                       ) -> Optional[Structure]:
    """Parse a coordinate block and perceive its bonds.

    The perception is the force field's own, so the editor and the force field
    never disagree about what is bonded -- including the corrections the user
    has made by hand, which the caller passes in as ``bond_orders``.
    """
    from .molecule_forcefield import parse_xyz, perceive_molecule

    parsed = parse_xyz(xyz_text)
    if parsed is None:
        return None
    symbols, coords, _had_header = parsed
    bonds: Dict[Tuple[int, int], int] = {}
    perceived = perceive_molecule(xyz_text)
    if perceived is not None:
        from .molecule_forcefield import _orders_from_mol
        known = _orders_from_mol(perceived.typing_mol)
        for pair in perceived.bonds:
            key = (min(pair), max(pair))
            bonds[key] = int(known.get(key, 1))
    for key, order in (bond_orders or {}).items():
        pair = (min(key), max(key))
        if max(pair) < len(symbols):
            if order:
                bonds[pair] = int(order)
            else:
                bonds.pop(pair, None)
    return Structure(symbols, coords, bonds)


def to_xyz(structure: Structure, comment: str = 'edited') -> str:
    """Render a structure back to a standard XYZ block."""
    rows = [f'{len(structure)}', comment]
    for symbol, (x, y, z) in zip(structure.symbols, structure.coords):
        rows.append(f'{symbol} {x:.6f} {y:.6f} {z:.6f}')
    return '\n'.join(rows) + '\n'


def adjust_hydrogens(structure: Structure,
                     touched: Iterable[int]) -> Dict[int, int]:
    """Fill or trim the hydrogens on the atoms an edit touched.

    Only those atoms: filling every unsatisfied valence in the molecule would
    quietly hydrogenate a radical or a donor the user put there on purpose.

    Every removal is collected and carried out in one pass at the end.
    Deleting as we went renumbered the structure underneath the loop, so the
    atoms still queued -- and the index the caller was holding on to -- became
    something else: growing a chain with a double bond returned the wrong
    atom and the bond it had just made looked absent.

    Returns:
        The old-index to new-index map, empty when nothing was removed.
    """
    doomed: set = set()
    for index in sorted({int(i) for i in touched}):
        if not 0 <= index < len(structure):
            continue
        symbol = structure.symbols[index]
        if symbol == 'H':
            continue
        wanted = default_valence(symbol)
        if wanted is None:
            continue
        # Hydrogens already condemned no longer count towards the valence.
        used = structure.used_valence(index) - sum(
            structure.order(index, j) for j in structure.hydrogens_on(index)
            if j in doomed
        )
        deficit = wanted - used
        if deficit > 0:
            centre = structure.coords[index]
            taken = [
                _sub(structure.coords[j], centre)
                for j in structure.neighbours(index) if j not in doomed
            ]
            # Steric number again: the lone pairs take room too, which is why
            # ammonia is pyramidal and water bent. Without them both came out
            # of here flat -- three hydrogens at 120 degrees, and two at 180.
            total = len(taken) + deficit + _LONE_PAIRS.get(symbol, 0)
            distance = bond_length(symbol, 'H')
            for direction in free_directions(taken, deficit, total):
                new = structure.add_atom(
                    'H', _add(centre, _scale(direction, distance)))
                structure.set_bond(index, new, 1)
        elif deficit < 0:
            spare = [j for j in structure.hydrogens_on(index) if j not in doomed]
            # Farthest first: the one most likely to be the odd hydrogen out.
            spare.sort(key=lambda j: -_norm(
                _sub(structure.coords[j], structure.coords[index])))
            doomed.update(spare[:min(len(spare), -deficit)])
    if not doomed:
        return {}
    return structure.remove_atoms(doomed)


def place_atom(structure: Structure, element: str, position: Sequence[float],
               bonded_to: Optional[int] = None, order: int = 1,
               adjust_h: bool = True) -> int:
    """Put a new atom down, optionally bonded to an existing one.

    *adjust_h* is the toolbar's Adjust H switch.  With it off the hydrogens
    are left exactly as they are, which is what somebody wants who is drawing
    a radical, an open coordination site, or a fragment that is going to be
    joined to something else in a moment.
    """
    index = structure.add_atom(element, position)
    touched = [index]
    if bonded_to is not None and 0 <= int(bonded_to) < len(structure) - 1:
        structure.set_bond(int(bonded_to), index, max(1, int(order)))
        touched.append(int(bonded_to))
    mapping = adjust_hydrogens(structure, touched) if adjust_h else {}
    if mapping:
        touched = [mapping.get(i, i) for i in touched if i in mapping]
        index = mapping.get(index, index)
    for atom in touched:
        rearrange_hydrogens(structure, atom)
    return index


def grow_from(structure: Structure, anchor: int, element: str,
              order: int = 1,
              direction: Optional[Sequence[float]] = None,
              adjust_h: bool = True) -> Optional[int]:
    """Hang a new atom off ``anchor`` at a sensible length and angle.

    A hydrogen already sitting in the way is consumed rather than pushed
    aside, which is what makes repeated growing build a chain instead of a
    thicket -- the same thing Avogadro does.
    """
    if not 0 <= int(anchor) < len(structure):
        return None
    anchor = int(anchor)
    centre = structure.coords[anchor]
    spare = structure.hydrogens_on(anchor)
    if spare:
        # Reuse the hydrogen's own direction when the caller gave none.
        if direction is None:
            direction = _sub(structure.coords[spare[0]], centre)
        structure.remove_atoms([spare[0]])
        if spare[0] < anchor:
            anchor -= 1
        centre = structure.coords[anchor]
    if direction is None or _norm(direction) < 1e-6:
        taken = [_sub(structure.coords[j], centre)
                 for j in structure.neighbours(anchor)]
        wanted = default_valence(structure.symbols[anchor]) or (len(taken) + 1)
        direction = free_directions(taken, 1, max(wanted, len(taken) + 1))[0]
    distance = bond_length(structure.symbols[anchor], element, order)
    position = _add(centre, _scale(_unit(direction), distance))
    index = structure.add_atom(element, position)
    structure.set_bond(anchor, index, max(1, int(order)))
    mapping = adjust_hydrogens(structure, [anchor, index]) if adjust_h else {}
    if mapping:
        anchor, index = mapping.get(anchor, anchor), mapping.get(index, index)
    # Both ends, or the hydrogens that were already on the anchor stay where
    # the old shape put them and the new atom lands among them: growing
    # straight at one of them left two atoms 0.45 A apart, which the renderer
    # then drew as a second bond to that hydrogen.
    rearrange_hydrogens(structure, anchor)
    rearrange_hydrogens(structure, index)
    return index


def set_element(structure: Structure, index: int, element: str,
                adjust_h: bool = True) -> bool:
    """Change one atom's element and re-satisfy its valence.

    *adjust_h* is the toolbar's Adjust H switch.  With it off the hydrogens
    are left exactly as they are, which is what somebody wants who is drawing
    a radical, an open coordination site, or a fragment that is going to be
    joined to something else in a moment.
    """
    if not 0 <= int(index) < len(structure):
        return False
    index = int(index)
    # The same element is not a no-op: a tap on an atom is also how a valence
    # that is short gets filled, so it falls through to the adjustment below.
    # Returning early here meant tapping a carbon with C selected did nothing
    # at all, and a half-built centre stayed half-built.
    same = structure.symbols[index] == element
    structure.symbols[index] = element
    # The bond it sits on wants a different length now.
    for neighbour in structure.neighbours(index):
        order = structure.order(index, neighbour)
        wanted = bond_length(element, structure.symbols[neighbour], order)
        offset = _sub(structure.coords[neighbour], structure.coords[index])
        length = _norm(offset)
        if length > 1e-6 and structure.symbols[neighbour] == 'H':
            structure.coords[neighbour] = _add(
                structure.coords[index], _scale(offset, wanted / length))
    before = len(structure)
    mapping = adjust_hydrogens(structure, [index]) if adjust_h else {}
    moved = mapping.get(index, index) if mapping else index
    rearrange_hydrogens(structure, moved)
    if same and len(structure) == before:
        return False        # nothing was missing and nothing changed
    return True


def delete_atoms(structure: Structure, indices: Iterable[int],
                 adjust_h: bool = True) -> int:
    """Remove atoms with their bonds, then tidy the valences left behind.

    *adjust_h* is the toolbar's Adjust H switch.  With it off the hydrogens
    are left exactly as they are, which is what somebody wants who is drawing
    a radical, an open coordination site, or a fragment that is going to be
    joined to something else in a moment.
    """
    drop = {int(i) for i in indices if 0 <= int(i) < len(structure)}
    if not drop:
        return 0
    # A hydrogen that hangs off nothing but a deleted atom goes with it.
    # Leaving it behind turned deleting one carbon of ethane into a carbon
    # with seven hydrogens: three orphans, plus the ones grown to fill the
    # valence the deletion had just freed.
    for index in list(drop):
        for neighbour in structure.neighbours(index):
            if structure.symbols[neighbour] != 'H':
                continue
            if all(j in drop for j in structure.neighbours(neighbour)):
                drop.add(neighbour)
    # A hydrogen the user deleted stays deleted. Refilling the valence it
    # freed grew it straight back, so removing one was a no-op: methane came
    # back as methane. Deleting a heavy atom is different -- there the
    # neighbour is left with a hole nobody asked for, and it gets a hydrogen.
    kept_open = {
        n for i in drop if structure.symbols[i] == 'H'
        for n in structure.neighbours(i)
    }
    neighbours = set()
    for index in drop:
        if structure.symbols[index] == 'H':
            continue
        neighbours.update(structure.neighbours(index))
    neighbours -= drop | kept_open
    mapping = structure.remove_atoms(drop)
    if adjust_h:
        adjust_hydrogens(
            structure, [mapping[n] for n in neighbours if n in mapping])
    return len(drop)


def rearrange_hydrogens(structure: Structure, index: int) -> int:
    """Move the hydrogens on one atom to the shape its bonds now imply.

    Raising a bond to double does not only change how many hydrogens fit, it
    changes where they belong: the carbons of ethene are trigonal planar, not
    two tetrahedra sharing an edge. Leaving them where they were made a
    molecule that looked like a twisted ethane and perceived like one -- a
    double bond read back off that geometry came out single.

    Only hydrogens move. A heavy substituent carries a whole branch with it,
    and shifting that to satisfy one centre would wreck the rest.
    """
    if not 0 <= int(index) < len(structure):
        return 0
    index = int(index)
    symbol = structure.symbols[index]
    if symbol == 'H' or default_valence(symbol) is None:
        return 0
    centre = structure.coords[index]
    hydrogens = structure.hydrogens_on(index)
    if not hydrogens:
        return 0
    heavy = [j for j in structure.neighbours(index) if structure.symbols[j] != 'H']
    taken = [_sub(structure.coords[j], centre) for j in heavy]
    # The steric number: partners plus lone pairs. Counting the sum of the
    # bond orders instead made ethene's carbons tetrahedral, because the
    # double bond was counted twice -- they have three partners, so they are
    # trigonal planar and the hydrogens belong at 120 degrees in that plane.
    total = len(heavy) + len(hydrogens) + _LONE_PAIRS.get(symbol, 0)
    distance = bond_length(symbol, 'H')
    directions = free_directions(taken, len(hydrogens), total)
    for hydrogen, direction in zip(hydrogens, directions):
        structure.coords[hydrogen] = _add(centre, _scale(_unit(direction), distance))
    return len(hydrogens)


def displaced_hydrogens(structure: Structure, first: int, second: int,
                        order: int = 1) -> List[int]:
    """Which hydrogens a new bond between these two atoms takes the place of.

    One per unit of valence the bond costs an end, and the ones it takes are
    the ones pointing most nearly at the new partner -- the hydrogens that
    were standing where the bond now is.

    An end with no valence of ours -- a metal -- gives up nothing, and neither
    does the other end of a bond to one: ammonia coordinating to platinum
    keeps all three of its hydrogens, and so does the platinum keep whatever
    it was given.
    """
    ends = (int(first), int(second))
    if not all(0 <= i < len(structure) for i in ends) or ends[0] == ends[1]:
        return []
    if any(default_valence(structure.symbols[i]) is None for i in ends):
        return []
    doomed: List[int] = []
    for index, partner in (ends, ends[::-1]):
        symbol = structure.symbols[index]
        if symbol == 'H':
            continue
        wanted = default_valence(symbol)
        excess = structure.used_valence(index) + max(1, int(order)) - wanted
        if excess <= 0:
            continue
        centre = structure.coords[index]
        towards = _unit(_sub(structure.coords[partner], centre))
        spare = [j for j in structure.hydrogens_on(index)
                 if j not in doomed and j not in ends]
        spare.sort(key=lambda j: -_dot(
            _unit(_sub(structure.coords[j], centre)), towards))
        doomed.extend(spare[:min(len(spare), excess)])
    return doomed


def connect_atoms(structure: Structure, first: int, second: int,
                  order: int = 1, adjust_h: bool = True) -> Dict[int, int]:
    """Bond two atoms that are both already there.

    Drawing a bond where there was none is the one edit that has to make room
    for itself: two methanes bonded at their carbons are not C2H8.  So with
    *adjust_h* on each end gives up the hydrogen the bond replaces, and with it
    off nothing goes -- which is what somebody wants who is correcting a
    perception rather than building.

    Nothing moves.  Only the topology changes and the displaced hydrogens go;
    pulling the two ends to a bond length and re-placing the hydrogens around
    them drags whole fragments through each other in a coordination sphere,
    which is why the Bond button never did it and still does not.

    Returns:
        The old-index to new-index map, empty when no hydrogen was displaced.
    """
    if not all(0 <= int(i) < len(structure) for i in (first, second)):
        return {}
    first, second = int(first), int(second)
    if first == second or structure.order(first, second):
        return {}
    doomed = displaced_hydrogens(structure, first, second, order) if adjust_h else []
    structure.set_bond(first, second, max(1, min(3, int(order))))
    return structure.remove_atoms(doomed) if doomed else {}


def set_bond_order(structure: Structure, first: int, second: int,
                   order: int, adjust_h: bool = True) -> bool:
    """Retype a bond and re-satisfy both ends.

    Setting a bond to double is not just a label: ethane becomes ethene, and
    the two hydrogens that no longer fit have to go. Without that the carbons
    are five-valent, the typing molecule refuses the order and the force field
    quietly hands back a single bond -- measured on benzene, where a drawn
    double bond came out at 1.514 A, exactly the single it had been asked to
    stop being.

    Passing 0 removes the bond and fills the valences it frees.
    """
    if not (0 <= int(first) < len(structure) and 0 <= int(second) < len(structure)):
        return False
    first, second = int(first), int(second)
    if first == second:
        return False
    order = max(0, min(3, int(order)))
    # No end can carry more than its own valence: a hydrogen has one, so a
    # tap asking for a double bond to it has to be refused rather than
    # producing a C=H that nothing downstream can type.
    if order:
        for index in (first, second):
            limit = default_valence(structure.symbols[index])
            if limit is not None and order > limit:
                return False
    if structure.order(first, second) == order:
        return False
    structure.set_bond(first, second, order)
    if order:
        wanted = bond_length(structure.symbols[first], structure.symbols[second], order)
        offset = _sub(structure.coords[second], structure.coords[first])
        length = _norm(offset)
        if length > 1e-6:
            structure.coords[second] = _add(
                structure.coords[first], _scale(offset, wanted / length))
    mapping = adjust_hydrogens(structure, [first, second]) if adjust_h else {}
    for index in (first, second):
        rearrange_hydrogens(structure, mapping.get(index, index) if mapping else index)
    return True
