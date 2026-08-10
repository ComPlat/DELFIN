"""A SMILES that describes two systems should come back as two systems.

``A.B`` in a SMILES means two molecules that are not bonded to each other -- a
complex and its counter-ion, a host and its guest, a substrate beside a
catalyst.  Handed to a converter whole, they are embedded into one coordinate
set and land on top of one another: measured on
``hexaphenylbenzene.benzene``, the benzene came out *inside* the other
molecule, closest approach 0.877 A and the two centres 2.22 A apart.  Nothing
can be done with that in a viewer -- the atoms cannot be told apart by eye, and
dragging one drags through the other.

So each part is built on its own and they are then set beside each other, far
enough apart to be two things and near enough to be one picture.  Which
converter builds them is the caller's business; this only splits, places and
joins.
"""

from __future__ import annotations

import math
from typing import Any, Callable, Dict, List, Optional, Tuple

__all__ = ['split_smiles', 'has_separate_systems', 'place_beside',
           'convert_each', 'combine_isomers', 'closest_approach']

#: How much clear space to leave between two systems, in Angstrom, measured
#: between the spheres that just contain them.  Three is wide enough that no
#: pair of atoms is within bonding distance and narrow enough that both are on
#: screen at a useful size.
GAP = 3.0


def split_smiles(smiles: Any) -> List[str]:
    """The parts of a SMILES that are not bonded to each other.

    Only a dot at the top level separates molecules.  One inside brackets is
    part of an atom -- a ``[C.]`` radical dot in some dialects -- and one
    inside parentheses is still the same molecule, so the nesting is counted
    rather than the string split.
    """
    text = str(smiles or '').strip()
    if not text:
        return []
    parts, current = [], []
    brackets = rounds = 0
    for letter in text:
        if letter == '[':
            brackets += 1
        elif letter == ']':
            brackets = max(0, brackets - 1)
        elif letter == '(':
            rounds += 1
        elif letter == ')':
            rounds = max(0, rounds - 1)
        if letter == '.' and not brackets and not rounds:
            parts.append(''.join(current))
            current = []
            continue
        current.append(letter)
    parts.append(''.join(current))
    return [p.strip() for p in parts if p.strip()]


def has_separate_systems(smiles: Any) -> bool:
    return len(split_smiles(smiles)) > 1


def _atoms(xyz_text: Any) -> List[Tuple[str, float, float, float]]:
    """The coordinate lines, whether or not a header is in front of them."""
    lines = [line for line in str(xyz_text or '').splitlines() if line.strip()]
    if not lines:
        return []
    start = 0
    try:
        int(lines[0].split()[0])
        start = 2                      # a count and its comment line
    except (ValueError, IndexError):
        start = 0
    out = []
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            out.append((parts[0], float(parts[1]), float(parts[2]),
                        float(parts[3])))
        except ValueError:
            continue
    return out


def _centre(atoms: List[Tuple[str, float, float, float]]) -> Tuple[float, float, float]:
    n = len(atoms)
    return tuple(sum(a[i] for a in atoms) / n for i in (1, 2, 3))


def _radius(atoms: List[Tuple[str, float, float, float]],
            centre: Tuple[float, float, float]) -> float:
    """How far the outermost atom sits from the middle."""
    return max((math.dist(a[1:], centre) for a in atoms), default=0.0)


def place_beside(blocks: List[str], gap: float = GAP,
                 comment: str = 'built separately, set side by side') -> str:
    """Put each block next to the last one instead of on top of it.

    Along one axis, each system centred on that axis and pushed out by its own
    reach plus the one before it: two spheres that just contain them, set apart
    by *gap*.  A row rather than a cloud, because a row is what a viewer shows
    without having to be turned.
    """
    groups = [_atoms(block) for block in blocks]
    groups = [g for g in groups if g]
    if not groups:
        return ''
    placed: List[str] = []
    edge = 0.0                          # where the last system stopped
    for atoms in groups:
        centre = _centre(atoms)
        reach = _radius(atoms, centre)
        offset = edge + reach if placed else 0.0
        for symbol, x, y, z in atoms:
            placed.append(
                f'{symbol} {x - centre[0] + offset:.8f} '
                f'{y - centre[1]:.8f} {z - centre[2]:.8f}')
        edge = offset + reach + gap
    return f'{len(placed)}\n{comment}\n' + '\n'.join(placed) + '\n'


def closest_approach(blocks: List[str]) -> float:
    """The nearest any two atoms of different systems come, in Angstrom.

    What "they do not clash" means, said as a number a test can check.
    """
    groups = [_atoms(block) for block in blocks]
    groups = [g for g in groups if g]
    worst = float('inf')
    for i, first in enumerate(groups):
        for second in groups[i + 1:]:
            for a in first:
                for b in second:
                    worst = min(worst, math.dist(a[1:], b[1:]))
    return worst


def convert_each(
    smiles: Any,
    convert: Callable[[str], Tuple[Any, int, Any, Any]],
    gap: float = GAP,
) -> Dict[str, Any]:
    """Build every part of *smiles* on its own and set them side by side.

    *convert* is whichever conversion the user pressed -- it is handed one
    SMILES at a time and is expected to answer the way they all do,
    ``(xyz, atoms, method, error)``.

    Returns ``{'ok', 'xyz', 'atoms', 'parts', 'method', 'error'}``.  A single
    system goes through untouched: there is nothing to place beside anything,
    and rebuilding it would only move it.
    """
    parts = split_smiles(smiles)
    if len(parts) < 2:
        xyz, atoms, method, error = convert(str(smiles or '').strip())
        return {'ok': error is None, 'xyz': xyz, 'atoms': atoms,
                'parts': 1, 'method': method, 'error': error}

    blocks, methods = [], []
    for position, part in enumerate(parts, start=1):
        xyz, _atom_count, method, error = convert(part)
        if error or not xyz:
            return {
                'ok': False, 'xyz': None, 'atoms': 0, 'parts': len(parts),
                'method': method,
                'error': (f'part {position} of {len(parts)} ({part[:40]}) '
                          f'could not be built: {error or "nothing came back"}'),
            }
        blocks.append(xyz)
        methods.append(method)
    joined = place_beside(
        blocks, gap=gap,
        comment=f'{len(parts)} systems, built separately and set side by side')
    return {
        'ok': True, 'xyz': joined, 'atoms': len(_atoms(joined)),
        'parts': len(parts),
        'method': methods[0] if len(set(map(str, methods))) == 1
        else ', '.join(str(m) for m in methods),
        'error': None,
    }


def combine_isomers(per_part: List[List[Any]], gap: float = GAP) -> List[Any]:
    """One row of systems per frame, from a list of isomers for each part.

    A metal complex may come back as twelve arrangements and the counter-ion
    beside it as one.  Every pairing of the two would be twelve frames of the
    same counter-ion, which is what a combinatorial product buys here: nothing.
    So the part with the most arrangements drives the navigation, and the
    others stand at whichever of theirs exists -- their last, once they have
    run out.

    Each isomer is whatever the converter returned, a tuple whose first element
    is the coordinates; the rest is carried through from the part that has the
    most to say.
    """
    lists = [list(part) for part in per_part if part]
    if not lists:
        return []
    frames = max(len(part) for part in lists)
    leader = max(range(len(lists)), key=lambda i: len(lists[i]))
    out: List[Any] = []
    for index in range(frames):
        blocks = []
        for part in lists:
            item = part[min(index, len(part) - 1)]
            blocks.append(item[0] if isinstance(item, (tuple, list)) else item)
        joined = place_beside(
            blocks, gap=gap,
            comment=f'{len(lists)} systems, built separately '
                    f'(arrangement {index + 1} of {frames})')
        head = lists[leader][min(index, len(lists[leader]) - 1)]
        rest = tuple(head[1:]) if isinstance(head, (tuple, list)) else ()
        out.append((joined,) + rest)
    return out
