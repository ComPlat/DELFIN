"""Force-field terms and clean-up relaxation for interactive molecule editing.

This module is the *Python* half of the Submit-tab live manipulation feature
(grab an atom, drag it, the rest of the molecule follows).  The force field
itself runs in the browser: a kernel round trip costs ~45 ms, which collapses
a continuous drag to ~13 Hz, so no per-frame Python is possible.  Python is
therefore called exactly twice per manipulation session:

* :func:`export_forcefield_terms` -- once, when manipulate mode is switched
  on.  It returns a JSON-serialisable dict of force-field terms that the
  browser engine evaluates itself, frame after frame, with zero further
  kernel traffic.
* :func:`relax_xyz` -- once, on mouse release, for a proper constrained
  minimisation of the dragged structure (the browser engine is deliberately
  cheap and only does a handful of steepest-descent steps per frame).

Payload contract
================

``export_forcefield_terms`` returns::

    {
      'ok':       bool,            # False -> the other fields are empty/meaningless
      'source':   str,             # 'rdkit-uff' | 'openbabel-uff' | 'geometric-fallback'
      'n_atoms':  int,
      'elements': ['C', 'H', ...], # length n_atoms, IN INPUT XYZ ORDER
      'bonds':    [{'i': int, 'j': int, 'k': float, 'r0': float}, ...],
      'angles':   [{'i': int, 'j': int, 'k': int, 'kt': float,
                    'theta0': float, 'n': int}, ...],
      'torsions': [{'i': int, 'j': int, 'k': int, 'l': int,
                    'v': float, 'n': int, 'phi0': float}, ...],
      'vdw':      [{'x': float, 'd': float}, ...],   # length n_atoms, per atom
      'metals':   [{'index': int, 'element': str}, ...],
      'warnings': [str, ...],
      'version':  int,
    }

**Indices are 0-based and index the input XYZ atom order exactly.**  Nothing
in this module ever reorders, adds or drops atoms, so entry ``i`` of
``elements``/``vdw`` is atom ``i`` of the XYZ block the browser already holds
(and therefore 3Dmol atom ``i``).  Element symbols are returned canonically
capitalised (``NI`` in, ``Ni`` out); only the spelling changes, never the
order.

Recommended call sequence::

    perceived = perceive_molecule(xyz)                       # switch-on
    payload   = export_forcefield_terms(xyz, perceived=perceived)
    ...                                                      # drag: browser only
    result    = relax_xyz(dragged_xyz, [dragged_index], perceived=perceived)

Passing the switch-on ``perceived`` into :func:`relax_xyz` is worth doing:
it skips a second perception and, more importantly, pins the topology, so a
drag can never silently re-bond the molecule between frames.

Units are UFF's: kcal/mol, Angstrom, and *degrees* for ``theta0``/``phi0``
(the consumer converts).  The energy expressions the terms belong to are the
UFF ones as implemented by RDKit, verified numerically against
``UFFGetMoleculeForceField(...).CalcEnergy()``:

``bonds``
    ``E = 0.5 * k * (r - r0)**2``
``angles``
    ``n == 0``  ->  general cosine expansion, with ``t0 = radians(theta0)``::

        C2 = 1 / (4 * sin(t0)**2)
        C1 = -4 * C2 * cos(t0)
        C0 = C2 * (2 * cos(t0)**2 + 1)
        E  = kt * (C0 + C1*cos(theta) + C2*cos(2*theta))

    ``n > 0``  ->  periodic form used by UFF for the linear / trigonal /
    square-planar-octahedral ideal angles::

        E = (kt / n**2) * (1 - cos(n * theta0) * cos(n * theta))

    ``n`` is derived from ``theta0`` (180 deg -> 1, 120 deg -> 3, 90 deg -> 4,
    anything else -> 0) and is exported only so the browser cannot get that
    switch wrong.
``torsions``
    ``E = v * (1 - cos(n*phi0) * cos(n*phi))``.  ``v`` is *ready to use*: it
    already contains UFF's factor 1/2 and the division by the number of
    torsions about the central bond, both of which RDKit applies internally
    but does *not* include in ``GetUFFTorsionParams``.
``vdw``
    Per-atom entries combined with the plain geometric mean (verified against
    RDKit's pair values)::

        x_ij = sqrt(x_i * x_j)
        D_ij = sqrt(d_i * d_j)
        E    = D_ij * ((x_ij/r)**12 - 2 * (x_ij/r)**6)

    Per-atom export is deliberate: N entries (~10 KiB at 400 atoms) instead
    of the O(N^2) pair table (2.85 MiB at 392 atoms).  **Exclusions:** pairs
    separated by fewer than three bonds carry no vdW term (1-2 and 1-3
    excluded, 1-4 included) -- the consumer derives that from ``bonds``.

Deliberate omissions the consumer should know about:

* No UFF *inversion* (out-of-plane) terms.  Planarity of sp2 centres is held
  by the torsion terms alone during a drag; :func:`relax_xyz` restores it
  exactly on mouse release because it runs the full force field.
* No electrostatics.  UFF's QEq charges are not part of RDKit's default UFF
  either, so this matches the reference implementation.
* Torsions whose central bond touches a metal are dropped (see below).
"""

from __future__ import annotations

import glob
import logging
import math
import os
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

try:  # pragma: no cover - import guard
    from rdkit import Chem
    from rdkit import rdBase
    from rdkit.Chem import rdDetermineBonds
    from rdkit.Chem import rdForceFieldHelpers as _uff_helpers
    from rdkit.Geometry import Point3D
    RDKIT_AVAILABLE = True
except ImportError:  # pragma: no cover - depends on environment
    RDKIT_AVAILABLE = False

try:  # pragma: no cover - import guard
    from openbabel import pybel
    OPENBABEL_AVAILABLE = True
    try:
        pybel.ob.obErrorLog.SetOutputLevel(pybel.ob.obError)
    except Exception:
        pass
except ImportError:  # pragma: no cover - depends on environment
    OPENBABEL_AVAILABLE = False


#: Payload schema version.  Bump when a field changes meaning so the browser
#: engine can refuse a payload it does not understand.
PAYLOAD_VERSION = 1

SOURCE_RDKIT = 'rdkit-uff'
SOURCE_OPENBABEL = 'openbabel-uff'
SOURCE_GEOMETRIC = 'geometric-fallback'

#: Above this atom count the export is noticeably slow (~84 ms at 392 atoms)
#: and the browser engine can no longer hold 30 fps; the payload is still
#: produced, but a warning is attached so the UI can say so.
LARGE_MOLECULE_ATOMS = 400

#: Force constants used only when no UFF parameters at all can be obtained
#: for an atom pair/triple (unknown element, missing parameter file).
_FALLBACK_BOND_K = 300.0     # kcal/mol/A^2
_FALLBACK_ANGLE_K = 50.0     # kcal/mol/rad^2

#: UFF bond force constant prefactor, k = 664.12 * Zi * Zj / r0^3.
_UFF_K_PREFACTOR = 664.12

#: Group-16 elements: UFF gives sp3-sp3 bonds between them a two-fold,
#: 90 degree torsion instead of the usual three-fold, 180 degree one.
_GROUP16 = {8, 16, 34, 52, 84}

#: Angle theta0 values (degrees) that switch UFF to the periodic angle form,
#: mapped to their periodicity.
_ANGLE_PERIODICITY = ((180.0, 1), (120.0, 3), (90.0, 4))

#: Beyond this deviation from an ideal angle the periodic form is not used.
_ANGLE_IDEAL_TOL = 0.5

#: A metal-centred angle this close to linear cannot use the general cosine
#: expansion (C2 diverges as sin(theta0) -> 0); it gets the n = 1 form.
_LINEAR_ANGLE_CUTOFF = 175.0

#: Highest coordination number a light element can reach.  Used only to prune
#: bonds that purely geometric perception invented on a strained structure.
_MAX_COORDINATION = {
    'H': 1, 'B': 4, 'C': 4, 'N': 4, 'O': 2, 'F': 1,
    'Si': 4, 'P': 5, 'S': 6, 'Cl': 1, 'Br': 1, 'I': 1,
}

#: Iteration bound for RDKit's combinatorial bond-order search, which is only
#: ever reached when Open Babel is unavailable.  Without a bound it can run
#: for minutes on a strained geometry.
_BOND_ORDER_MAX_ITERATIONS = 2000


# --------------------------------------------------------------------------
# small utilities
# --------------------------------------------------------------------------

class _RDKitQuiet:
    """Silence RDKit's C++ log while querying parameters for odd molecules.

    Typing a metal complex makes RDKit print one ``UFFTYPER: Unrecognized
    atom type`` line per attempt, which would flood the notebook.  The block
    is released again when the context exits.
    """

    def __enter__(self):
        self._blocker = None
        if RDKIT_AVAILABLE:
            try:
                self._blocker = rdBase.BlockLogs()
            except Exception:
                self._blocker = None
        return self

    def __exit__(self, exc_type, exc, tb):
        self._blocker = None
        return False


def _normalise_symbol(raw: str) -> str:
    """Return an element symbol in canonical capitalisation ('NI' -> 'Ni')."""
    token = ''.join(ch for ch in str(raw).strip() if ch.isalpha())
    if not token:
        return ''
    return token[0].upper() + token[1:].lower()


def _is_finite(value: float) -> bool:
    return not (math.isnan(value) or math.isinf(value))


def parse_xyz(xyz_text: str) -> Optional[Tuple[List[str], List[Tuple[float, float, float]], bool]]:
    """Parse a coordinate block into ``(symbols, coords, had_header)``.

    Both flavours used across the code base are accepted:

    * standard XYZ -- atom count on line 1, comment on line 2, then one
      ``symbol x y z`` line per atom;
    * the DELFIN block -- ``symbol x y z`` lines only, no header.

    Trailing junk after the declared atom count is ignored.  Returns ``None``
    instead of raising when the text is not a coordinate block at all.
    """
    if not xyz_text or not str(xyz_text).strip():
        return None

    raw_lines = [ln for ln in str(xyz_text).splitlines()]
    stripped = [ln.strip() for ln in raw_lines]
    non_empty = [ln for ln in stripped if ln]
    if not non_empty:
        return None

    had_header = False
    body = non_empty
    declared = None
    try:
        declared = int(non_empty[0].split()[0])
    except (ValueError, IndexError):
        declared = None
    if declared is not None and declared > 0:
        # Standard XYZ: the comment line may legitimately be blank, so work
        # from the raw (unfiltered) lines here.
        if len(raw_lines) >= 2:
            candidate = [ln.strip() for ln in raw_lines[2:] if ln.strip()]
            if len(candidate) >= declared:
                body = candidate[:declared]
                had_header = True

    symbols: List[str] = []
    coords: List[Tuple[float, float, float]] = []
    for line in body:
        parts = line.split()
        if len(parts) < 4:
            return None
        symbol = _normalise_symbol(parts[0])
        if not symbol:
            return None
        try:
            xyz = (float(parts[1]), float(parts[2]), float(parts[3]))
        except ValueError:
            return None
        if not all(_is_finite(v) for v in xyz):
            return None
        symbols.append(symbol)
        coords.append(xyz)

    if not symbols:
        return None
    return symbols, coords, had_header


def format_xyz(
    symbols: Sequence[str],
    coords: Sequence[Sequence[float]],
    *,
    header: bool = True,
    comment: str = '',
) -> str:
    """Render coordinates back to text (standard XYZ or the DELFIN block)."""
    lines: List[str] = []
    if header:
        lines.append(str(len(symbols)))
        lines.append(comment)
    for symbol, xyz in zip(symbols, coords):
        lines.append(
            f"{symbol:4s} {float(xyz[0]):12.6f} {float(xyz[1]):12.6f} {float(xyz[2]):12.6f}"
        )
    return '\n'.join(lines) + '\n'


def _metal_symbols() -> frozenset:
    """Return the repo's metal element set.

    ``delfin.smiles_converter._METAL_SET`` is the single list the rest of
    DELFIN classifies metals with (it also backs ``contains_metal``), so it is
    reused here rather than duplicated.  Only if that import is unavailable
    (stand-alone use of this module) is a rule-based substitute derived from
    atomic numbers.
    """
    try:
        from delfin.smiles_converter import _METAL_SET  # noqa: WPS433 (lazy on purpose)
        return frozenset(_METAL_SET)
    except Exception:  # pragma: no cover - only without the full package
        logger.debug('smiles_converter._METAL_SET unavailable; using periodic-table rule')
        if not RDKIT_AVAILABLE:
            return frozenset()
        table = Chem.GetPeriodicTable()
        metals = set()
        for number in range(3, 95):
            try:
                symbol = table.GetElementSymbol(number)
            except Exception:
                continue
            # Everything that is not H/He, a noble gas, a halogen or one of
            # the classic non-metals/metalloids counts as a metal here.
            if symbol in {
                'B', 'C', 'N', 'O', 'F', 'Ne', 'Si', 'P', 'S', 'Cl', 'Ar',
                'Ge', 'As', 'Se', 'Br', 'Kr', 'Sb', 'Te', 'I', 'Xe',
                'At', 'Rn', 'He',
            }:
                continue
            metals.add(symbol)
        return frozenset(metals)


# --------------------------------------------------------------------------
# UFF parameter table (Open Babel's UFF.prm)
# --------------------------------------------------------------------------

_UFF_PRM_CACHE: Optional[Dict[str, Dict[str, float]]] = None
_UFF_PRM_BY_ELEMENT: Optional[Dict[str, List[Tuple[str, Dict[str, float]]]]] = None

_PRM_COLUMNS = (
    'r1', 'theta0', 'x1', 'D1', 'zeta', 'Z1', 'Vi', 'Uj', 'Xi', 'hard', 'radius',
)

#: Coordination number implied by the geometry digit of a UFF atom type
#: ('Zn3+2' -> 3 -> tetrahedral -> CN 4).
_UFF_GEOMETRY_CN = {'1': 2, '2': 3, '3': 4, '4': 4, '5': 5, '6': 6}


def _uff_prm_candidates() -> List[str]:
    """Return plausible filesystem locations of Open Babel's ``UFF.prm``."""
    patterns: List[str] = []
    data_dir = os.environ.get('BABEL_DATADIR')
    if data_dir:
        patterns.append(os.path.join(data_dir, 'UFF.prm'))
    if OPENBABEL_AVAILABLE:
        try:
            import openbabel as _ob_pkg
            root = os.path.dirname(os.path.abspath(_ob_pkg.__file__))
            patterns.append(os.path.join(root, 'share', 'openbabel', '*', 'UFF.prm'))
        except Exception:
            pass
    patterns.append(os.path.join(sys.prefix, 'share', 'openbabel', '*', 'UFF.prm'))
    patterns.append('/usr/share/openbabel/*/UFF.prm')
    patterns.append('/usr/local/share/openbabel/*/UFF.prm')

    found: List[str] = []
    for pattern in patterns:
        try:
            found.extend(sorted(glob.glob(pattern)) if '*' in pattern else
                         ([pattern] if os.path.isfile(pattern) else []))
        except Exception:
            continue
    return found


def _element_of_uff_type(uff_type: str) -> str:
    """Extract the element symbol from a UFF atom type ('Fe6+2' -> 'Fe')."""
    if not uff_type:
        return ''
    if len(uff_type) == 1:
        return uff_type.upper()
    second = uff_type[1]
    if second == '_':
        return uff_type[0].upper()
    if second.isalpha() and second.islower():
        return uff_type[0].upper() + second
    return uff_type[0].upper()


def _load_uff_prm() -> Dict[str, Dict[str, float]]:
    """Parse Open Babel's ``UFF.prm`` once and cache it.

    This is the published UFF parameter table (Rappe et al., JACS 1992,
    114, 10024) shipped as data with Open Babel 3.x.  It is the source of
    truth for the elements RDKit's UFF implementation refuses to type --
    every transition metal, in particular -- so the fallback path can use
    real UFF numbers instead of invented ones.  An empty dict is returned
    (and the caller degrades further) when the file cannot be found.
    """
    global _UFF_PRM_CACHE, _UFF_PRM_BY_ELEMENT
    if _UFF_PRM_CACHE is not None:
        return _UFF_PRM_CACHE

    table: Dict[str, Dict[str, float]] = {}
    for path in _uff_prm_candidates():
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as handle:
                for line in handle:
                    if not line.startswith('param'):
                        continue
                    parts = line.split()
                    if len(parts) < 2 + len(_PRM_COLUMNS) - 1:
                        continue
                    name = parts[1]
                    values = parts[2:]
                    entry: Dict[str, float] = {}
                    for key, raw in zip(_PRM_COLUMNS, values):
                        try:
                            entry[key] = float(raw)
                        except ValueError:
                            entry[key] = 0.0
                    if len(entry) >= 6:
                        table.setdefault(name, entry)
        except Exception as exc:  # pragma: no cover - filesystem dependent
            logger.debug('Could not read UFF parameter file %s: %s', path, exc)
            continue
        if table:
            logger.debug('Loaded %d UFF atom types from %s', len(table), path)
            break

    by_element: Dict[str, List[Tuple[str, Dict[str, float]]]] = {}
    for name, entry in table.items():
        by_element.setdefault(_element_of_uff_type(name), []).append((name, entry))

    _UFF_PRM_CACHE = table
    _UFF_PRM_BY_ELEMENT = by_element
    return table


def uff_atom_parameters(symbol: str, coordination: int = 0) -> Optional[Dict[str, float]]:
    """Return UFF parameters for ``symbol`` at the given coordination number.

    When an element has several UFF types (``Fe3+2`` tetrahedral vs
    ``Fe6+2`` octahedral) the one whose implied coordination number is
    closest to the observed one is chosen, which also picks the right
    ``theta0``.  Returns ``None`` for elements absent from the table.
    """
    _load_uff_prm()
    if not _UFF_PRM_BY_ELEMENT:
        return None
    candidates = _UFF_PRM_BY_ELEMENT.get(_normalise_symbol(symbol))
    if not candidates:
        return None
    if len(candidates) == 1 or coordination <= 0:
        return candidates[0][1]

    def _score(item: Tuple[str, Dict[str, float]]) -> Tuple[int, int]:
        name = item[0]
        digit = next((ch for ch in name[len(_element_of_uff_type(name)):] if ch.isdigit()), None)
        implied = _UFF_GEOMETRY_CN.get(digit, 0) if digit else 0
        if implied == 0:
            return (99, 0)
        return (abs(implied - coordination), 0)

    return min(candidates, key=_score)[1]


def _uff_bond_force_constant(z_i: float, z_j: float, r0: float) -> float:
    """UFF bond force constant ``k = 664.12 * Zi * Zj / r0**3`` (kcal/mol/A^2).

    Verified against RDKit: C(sp3)-C(sp3) with Z* = 1.912 and r0 = 1.514 A
    gives 699.6, RDKit's ``GetUFFBondStretchParams`` returns 699.59.
    """
    if r0 <= 0.0:
        return _FALLBACK_BOND_K
    return _UFF_K_PREFACTOR * z_i * z_j / (r0 ** 3)


def _uff_angle_force_constant(
    z_i: float, z_k: float, r_ij: float, r_jk: float, theta0_rad: float
) -> float:
    """UFF angle force constant (kcal/mol/rad^2) for i-j-k.

    ``K = 664.12 * Zi*Zk/r_ik^5 * (3*r_ij*r_jk*(1-cos^2 t0) - r_ik^2 * cos t0)``
    with ``r_ik`` from the law of cosines.  Verified against RDKit: H-C-C
    gives 117.4 against RDKit's 117.32 (the difference is only the rounding
    of the input bond lengths).
    """
    cos_t0 = math.cos(theta0_rad)
    r_ik_sq = r_ij * r_ij + r_jk * r_jk - 2.0 * r_ij * r_jk * cos_t0
    if r_ik_sq <= 1e-8 or r_ij <= 0.0 or r_jk <= 0.0:
        return _FALLBACK_ANGLE_K
    r_ik = math.sqrt(r_ik_sq)
    bracket = 3.0 * r_ij * r_jk * (1.0 - cos_t0 * cos_t0) - r_ik_sq * cos_t0
    value = _UFF_K_PREFACTOR * z_i * z_k / (r_ik ** 5) * bracket
    if not _is_finite(value) or value <= 0.0:
        return _FALLBACK_ANGLE_K
    return value


# --------------------------------------------------------------------------
# perception
# --------------------------------------------------------------------------

@dataclass
class PerceivedMolecule:
    """Everything downstream needs about one coordinate block.

    Attributes:
        symbols: Element symbols, input XYZ order.
        coords: Cartesian coordinates in Angstrom, input XYZ order.
        n_atoms: Number of atoms.
        mol: RDKit molecule carrying the perceived connectivity *including*
            metal-ligand bonds, with one conformer holding ``coords``.
        typing_mol: Same atoms and indices, but with metal-ligand bonds
            removed and bond orders assigned, sanitized so RDKit can assign
            UFF atom types to the organic part.  ``None`` when sanitisation
            was impossible.
        bonds: Perceived bonds as sorted ``(i, j)`` index pairs.
        metal_indices: Indices of metal atoms (repo-wide metal definition).
        has_metal: True when ``metal_indices`` is non-empty.
        bond_orders_known: False when bond-order perception failed and every
            bond is treated as single (hybridisation is then unreliable).
        had_header: True when the input carried the standard XYZ count line.
        warnings: Human-readable notes about approximations made.
        auto_hybridisation: What perception said about an atom, recorded before
            an override replaced it, so the offer can still name it.
        forced_hybridisation: Atom index to the hybridisation the user forced.
            The angles at such an atom are built from it directly.
    """

    symbols: List[str]
    coords: List[Tuple[float, float, float]]
    n_atoms: int
    mol: Any = None
    typing_mol: Any = None
    bonds: List[Tuple[int, int]] = field(default_factory=list)
    metal_indices: List[int] = field(default_factory=list)
    has_metal: bool = False
    bond_orders_known: bool = True
    had_header: bool = True
    warnings: List[str] = field(default_factory=list)
    auto_hybridisation: Dict[int, str] = field(default_factory=dict)
    forced_hybridisation: Dict[int, str] = field(default_factory=dict)

    def neighbours(self) -> List[List[int]]:
        """Return the adjacency list implied by :attr:`bonds`."""
        adjacency: List[List[int]] = [[] for _ in range(self.n_atoms)]
        for i, j in self.bonds:
            adjacency[i].append(j)
            adjacency[j].append(i)
        return adjacency


def _openbabel_bond_orders(xyz_text: str) -> Optional[Dict[Tuple[int, int], int]]:
    """Perceive Kekule bond orders with Open Babel.

    Open Babel's XYZ reader preserves atom order, so the returned keys are
    0-based indices into the same atom order.  Open Babel is used because it
    perceives bond orders for metal complexes, which RDKit's
    ``DetermineBondOrders`` refuses to do (it raises for any element outside
    its valence table).
    """
    if not OPENBABEL_AVAILABLE:
        return None
    try:
        ob_mol = pybel.readstring('xyz', xyz_text)
    except Exception as exc:
        logger.debug('Open Babel could not read the XYZ block: %s', exc)
        return None
    orders: Dict[Tuple[int, int], int] = {}
    try:
        for bond in pybel.ob.OBMolBondIter(ob_mol.OBMol):
            i = bond.GetBeginAtomIdx() - 1
            j = bond.GetEndAtomIdx() - 1
            orders[(min(i, j), max(i, j))] = int(bond.GetBondOrder())
    except Exception as exc:
        logger.debug('Open Babel bond iteration failed: %s', exc)
        return None
    return orders


def _build_typing_mol(
    mol: Any,
    metal_indices: Sequence[int],
    orders: Optional[Dict[Tuple[int, int], int]],
    warnings: List[str],
) -> Tuple[Any, bool]:
    """Return a sanitized copy of ``mol`` that RDKit can assign UFF types to.

    Atom indices are preserved exactly (bonds are removed, never atoms).
    Metal-ligand bonds are stripped because RDKit's UFF typer cannot type a
    metal at all and an attached metal poisons the typing of its donors;
    donor atoms get their implicit-hydrogen ban lifted so the free ligand
    perceives as the neutral, aromatic species it normally is.

    Assemblies are tried in order of how much they give up.  A strained
    geometry -- which is what a drag hands over -- and a bond the user drew
    can both over-fill an atom's valence; keeping such a molecule would only
    defer the exception to the moment the force field is built.  So the
    perceived Kekule pattern is first repaired *locally*, by lowering the
    order of a multiple bond at the offending atom, and only if that still
    fails are the bond orders dropped wholesale.  The difference matters: a
    single drawn bond used to cost every other atom in the molecule its
    hybridisation.
    """
    bond_types = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
    }
    metals = set(int(m) for m in metal_indices)

    def _prune_overbonded(editable) -> int:
        """Drop the longest bonds of atoms that exceed their maximum coordination.

        ``DetermineConnectivity`` is purely geometric, so a strained drag can
        give a carbon five neighbours -- a valence RDKit rejects outright and
        which would leave the whole molecule untypeable.  The surplus bonds
        are removed from the *typing* copy only; the exported bond list still
        contains them and they simply fall back to geometric parameters.
        """
        try:
            positions = mol.GetConformer().GetPositions()
        except Exception:
            return 0
        removed = 0
        for atom in editable.GetAtoms():
            limit = _MAX_COORDINATION.get(atom.GetSymbol())
            if limit is None:
                continue
            neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
            if len(neighbours) <= limit:
                continue
            idx = atom.GetIdx()
            ranked = sorted(
                neighbours,
                key=lambda other: -float(
                    ((positions[idx] - positions[other]) ** 2).sum()
                ),
            )
            for other in ranked[:len(neighbours) - limit]:
                editable.RemoveBond(idx, other)
                removed += 1
        return removed

    def _reduce_overvalence(editable) -> int:
        """Lower a multiple bond at every atom that carries too much valence.

        Drawing a bond onto an aromatic carbon takes it to five bonds, which
        RDKit rejects.  The chemically honest reading is that the carbon gave
        up its double bond, so that is what happens here: the longest multiple
        bond at the atom drops one order, the partner keeps a radical electron
        and stays sp2, and every ring the edit did not touch keeps its own
        double bonds.  Bonds are only lowered, never removed -- the exported
        topology is unchanged.
        """
        try:
            positions = mol.GetConformer().GetPositions()
        except Exception:
            positions = None
        table = Chem.GetPeriodicTable()
        lowered = 0
        # Atoms left one bond short by a lowering.  Only these may be paired
        # up again below: every other under-valent atom was already like that
        # in the perceived molecule and is not this function's business.
        orphaned: set = set()

        def _capacity(atom) -> int:
            """Valence still free at this atom; negative when it is over-full."""
            try:
                default = table.GetDefaultValence(atom.GetAtomicNum())
            except Exception:
                return 0
            if default <= 0:
                return 0
            used = sum(
                int(round(b.GetBondTypeAsDouble())) for b in atom.GetBonds()
            ) + atom.GetNumExplicitHs()
            return default - used

        # Lowering a bond can push its partner over the edge in turn, so this
        # sweeps until the molecule stops changing.  The bound is a guard
        # against a pathological structure, not an expected exit.
        for _sweep in range(8):
            busy = False
            for atom in editable.GetAtoms():
                surplus = -_capacity(atom)
                if surplus <= 0:
                    continue
                bonds = list(atom.GetBonds())
                multiples = [b for b in bonds if b.GetBondTypeAsDouble() > 1.5]
                if not multiples:
                    continue  # too many single bonds: _prune_overbonded's job
                if positions is not None:
                    multiples.sort(key=lambda b: -float(
                        ((positions[b.GetBeginAtomIdx()]
                          - positions[b.GetEndAtomIdx()]) ** 2).sum()
                    ))
                victim = multiples[0]
                order = int(round(victim.GetBondTypeAsDouble()))
                victim.SetBondType(bond_types[max(1, order - surplus)])
                idx = atom.GetIdx()
                partner = (victim.GetEndAtomIdx() if victim.GetBeginAtomIdx() == idx
                           else victim.GetBeginAtomIdx())
                orphaned.add(partner)
                lowered += 1
                busy = True
            if not busy:
                break

        # Two neighbours that both lost a partner pair up again, because that
        # is what the chemistry does: joining two carbons across a benzene
        # ring gives Dewar benzene, whose remaining double bonds are real, not
        # a pair of radicals sitting next to each other.
        for a in sorted(orphaned):
            atom_a = editable.GetAtomWithIdx(a)
            if _capacity(atom_a) < 1:
                continue
            for neighbour in atom_a.GetNeighbors():
                b = neighbour.GetIdx()
                if b not in orphaned or _capacity(neighbour) < 1:
                    continue
                bond = editable.GetBondBetweenAtoms(a, b)
                order = int(round(bond.GetBondTypeAsDouble())) if bond else 0
                if not 1 <= order <= 2:
                    continue
                bond.SetBondType(bond_types[order + 1])
                break
        return lowered

    def _repair_charges(editable) -> None:
        """Give over-valent main-group atoms the formal charge they imply.

        A dragged carbon that ends up bonded to an amine nitrogen makes that
        nitrogen four-valent, which RDKit only accepts as N+.  Assigning the
        charge keeps the molecule sanitizable; it changes no UFF type, since
        those follow hybridisation.
        """
        table = Chem.GetPeriodicTable()
        for atom in editable.GetAtoms():
            if atom.GetSymbol() not in ('N', 'O', 'P', 'S'):
                continue
            valence = sum(
                int(round(bond.GetBondTypeAsDouble())) for bond in atom.GetBonds()
            ) + atom.GetNumExplicitHs()
            try:
                default = table.GetDefaultValence(atom.GetAtomicNum())
            except Exception:
                continue
            if 0 < default < valence:
                atom.SetFormalCharge(valence - default)

    def _assemble(with_orders: bool, prune: bool = False, reduce_valence: bool = False):
        editable = Chem.RWMol(mol)
        if prune:
            _prune_overbonded(editable)
            _repair_charges(editable)
        donors: List[int] = []
        for bond in list(editable.GetBonds()):
            i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if i in metals or j in metals:
                donors.append(j if i in metals else i)
                continue
            if with_orders and orders:
                order = orders.get((min(i, j), max(i, j)))
                if order:
                    bond.SetBondType(bond_types.get(order, Chem.BondType.SINGLE))
            else:
                bond.SetBondType(Chem.BondType.SINGLE)
        for i, j in [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in list(editable.GetBonds())]:
            if i in metals or j in metals:
                editable.RemoveBond(i, j)
        for idx in donors:
            atom = editable.GetAtomWithIdx(idx)
            atom.SetNoImplicit(False)
            atom.SetNumRadicalElectrons(0)
        # After the metal bonds are gone, so a donor is not counted as
        # over-valent for coordinating.
        if reduce_valence:
            _reduce_overvalence(editable)
        return editable.GetMol()

    attempts = (
        ('perceived bond orders', True, False, False, True),
        ('perceived bond orders, over-valence lowered', True, False, True, True),
        ('single bonds', False, False, False, True),
        ('single bonds, over-bonded atoms pruned and charges repaired', False, True, False, True),
        ('single bonds, lenient', False, True, False, False),
    )
    for label, with_orders, prune, reduce_valence, strict in attempts:
        candidate = Chem.Mol(_assemble(with_orders, prune, reduce_valence))
        try:
            candidate.UpdatePropertyCache(strict=False)
            if strict:
                Chem.SanitizeMol(candidate)
            else:
                Chem.SanitizeMol(
                    candidate,
                    Chem.SanitizeFlags.SANITIZE_ALL
                    ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES
                    ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                )
            if reduce_valence:
                warnings.append(
                    'An atom carried more bonds than its valence allows; a '
                    'multiple bond at each such atom was lowered by one order '
                    'when assigning force-field types, which is what drawing '
                    'a bond onto it implies.'
                )
            if not with_orders and orders:
                warnings.append(
                    'The perceived bond orders gave an impossible valence '
                    '(the geometry is strained); every bond was treated as '
                    'single, so hybridisation-dependent parameters are '
                    'approximate.'
                )
            if prune:
                warnings.append(
                    'The geometry implies impossible valences; surplus '
                    'contacts were ignored and formal charges adjusted when '
                    'assigning force-field types.'
                )
            return candidate, True
        except Exception as exc:
            logger.debug('Sanitisation (%s) of the typing molecule failed: %s', label, exc)

    warnings.append(
        'The molecular graph could not be sanitized; force-field atom types '
        'were derived from the geometry alone.'
    )
    return None, False


def _drop_impossible_metal_contacts(
    mol: Any,
    metal_indices: Sequence[int],
    symbols: Sequence[str],
    coords: Sequence[Sequence[float]],
    warnings: List[str],
) -> Any:
    """Return ``mol`` with any metal contact that gave a carbon a fifth bond cut.

    ``DetermineConnectivity`` is purely geometric, and a metal's covalent
    radius is large: on a scandium complex the cutoff reaches past 2.9 A, so
    the backbone carbons of a triazacyclononane -- held near the metal by
    their own chelate rings, bonded to nothing but their nitrogen -- were
    counted as donors.  The coordination number came out as 9 instead of 6,
    and only a nine-vertex polyhedron was offered for what is an octahedron.

    Carbon is the one element where this is unambiguous.  Every real carbon
    donor stays within four bonds: an alkyl has M and three substituents, an
    N-heterocyclic carbene M and two nitrogens, a side-on alkene carbon M, its
    partner and two hydrogens.  Five means the geometry was over-read.  The
    longest offending contact goes first, and only far enough to bring the
    carbon back to four.

    Deliberately carbon only.  Oxygen would be the obvious next candidate and
    is exactly wrong: a coordinated ether or a bridging alkoxide is
    three-connected on purpose.
    """
    metals = set(int(m) for m in metal_indices)
    if not metals:
        return mol
    removed = 0
    editable = Chem.RWMol(mol)
    for atom in list(editable.GetAtoms()):
        index = atom.GetIdx()
        if symbols[index] != 'C':
            continue
        neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
        surplus = len(neighbours) - 4
        if surplus <= 0:
            continue
        attached = [n for n in neighbours if n in metals]
        if not attached:
            continue
        attached.sort(key=lambda m: -_distance(coords[index], coords[m]))
        for metal in attached[:surplus]:
            editable.RemoveBond(index, metal)
            removed += 1
    if not removed:
        return mol
    warnings.append(
        f'{removed} metal-carbon contact(s) were dropped: they would have '
        'given a carbon five bonds, which the geometric perception reaches '
        'only because a metal has a large covalent radius. Use Bond if one '
        'of them was real.'
    )
    return editable.GetMol()


def perceive_molecule(xyz_text: str) -> Optional[PerceivedMolecule]:
    """Parse an XYZ block and perceive connectivity, bond orders and metals.

    The pipeline is: ``Chem.MolFromXYZBlock`` ->
    ``rdDetermineBonds.DetermineConnectivity`` (geometry-based connectivity,
    never raises for exotic elements) -> bond orders, first with RDKit's
    ``DetermineBondOrders`` and, when that refuses (metals, unusual charge
    states), with Open Babel's perception transferred onto the very same
    bonds.  Atom order and count are preserved throughout.

    Args:
        xyz_text: Standard XYZ text or a bare DELFIN ``symbol x y z`` block.

    Returns:
        A :class:`PerceivedMolecule`, or ``None`` when the input is not a
        usable coordinate block or RDKit is unavailable.  This function does
        not raise for odd input.
    """
    parsed = parse_xyz(xyz_text)
    if parsed is None:
        logger.debug('perceive_molecule: input is not a coordinate block')
        return None
    symbols, coords, had_header = parsed

    if not RDKIT_AVAILABLE:
        logger.debug('perceive_molecule: RDKit unavailable')
        return None

    warnings: List[str] = []
    metal_set = _metal_symbols()
    metal_indices = [i for i, s in enumerate(symbols) if s in metal_set]

    standard_xyz = format_xyz(symbols, coords, header=True, comment='perceive')
    try:
        with _RDKitQuiet():
            mol = Chem.MolFromXYZBlock(standard_xyz)
    except Exception as exc:
        logger.debug('MolFromXYZBlock failed: %s', exc)
        mol = None
    if mol is None or mol.GetNumAtoms() != len(symbols):
        logger.debug('perceive_molecule: RDKit did not accept the block')
        return None

    try:
        with _RDKitQuiet():
            rdDetermineBonds.DetermineConnectivity(mol)
    except Exception as exc:
        logger.debug('DetermineConnectivity failed: %s', exc)
        warnings.append('Connectivity perception failed; no bonded terms could be built.')

    mol = _drop_impossible_metal_contacts(
        mol, metal_indices, symbols, coords, warnings
    )

    bonds = sorted(
        (min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        for b in mol.GetBonds()
    )

    # Bond orders.  Open Babel is asked first, deliberately: its perception is
    # a linear-time geometric one that works for metals and, more importantly,
    # cannot run away.  RDKit's DetermineBondOrders searches valence/charge
    # assignments combinatorially -- on a *strained* geometry (exactly what a
    # drag produces) it has been measured to run for minutes on a 103-atom
    # peptide, which would freeze the notebook on mouse release.  It is kept
    # only as a fallback and is then bounded by ``maxIterations``.
    # Kekule orders are transferred, not aromatic flags: RDKit re-perceives
    # aromaticity itself during sanitisation and would otherwise refuse to
    # kekulize rings it did not choose.
    orders: Optional[Dict[Tuple[int, int], int]] = _openbabel_bond_orders(standard_xyz)
    ordered_mol = None
    # An empty mapping is a valid answer (a molecule without bonds); only
    # ``None`` means the perception itself failed.
    bond_orders_known = orders is not None

    if not bond_orders_known and not metal_indices:
        candidate = Chem.Mol(mol)
        try:
            with _RDKitQuiet():
                rdDetermineBonds.DetermineBondOrders(
                    candidate, charge=0, maxIterations=_BOND_ORDER_MAX_ITERATIONS
                )
                Chem.SanitizeMol(candidate)
            if candidate.GetNumAtoms() == len(symbols):
                ordered_mol = candidate
                bond_orders_known = True
        except Exception as exc:
            logger.debug('DetermineBondOrders failed: %s', exc)

    if not bond_orders_known:
        warnings.append(
            'Bond orders could not be perceived; every bond is treated as '
            'single, so hybridisation-dependent parameters are approximate.'
        )

    conformer = Chem.Conformer(mol.GetNumAtoms())
    for idx, (x, y, z) in enumerate(coords):
        conformer.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))
    mol.RemoveAllConformers()
    mol.AddConformer(conformer, assignId=True)

    with _RDKitQuiet():
        if ordered_mol is not None and not metal_indices:
            typing_mol = ordered_mol
        else:
            typing_mol, _ = _build_typing_mol(mol, metal_indices, orders, warnings)

    return PerceivedMolecule(
        symbols=symbols,
        coords=coords,
        n_atoms=len(symbols),
        mol=mol,
        typing_mol=typing_mol,
        bonds=bonds,
        metal_indices=metal_indices,
        has_metal=bool(metal_indices),
        bond_orders_known=bond_orders_known,
        had_header=had_header,
        warnings=warnings,
    )


def _orders_from_mol(mol: Any) -> Dict[Tuple[int, int], int]:
    """Read integer Kekule bond orders out of a molecule, keyed by index pair.

    Aromatic flags are resolved first: RDKit re-perceives aromaticity when it
    sanitizes, and would refuse to kekulize a ring it did not choose itself,
    so only whole orders are ever passed on.
    """
    if mol is None:
        return {}
    try:
        kekulized = Chem.Mol(mol)
        Chem.Kekulize(kekulized, clearAromaticFlags=True)
    except Exception as exc:
        logger.debug('Kekulisation for order transfer failed: %s', exc)
        kekulized = mol
    orders: Dict[Tuple[int, int], int] = {}
    for bond in kekulized.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = int(round(bond.GetBondTypeAsDouble()))
        if order >= 1:
            orders[(min(i, j), max(i, j))] = order
    return orders


def apply_bond_edits(perceived: Any, edits: Any) -> bool:
    """Lay hand-drawn bond corrections onto a perceived molecule.

    ``edits`` maps an ``(i, j)`` atom pair to True (this bond exists) or False
    (it does not).  Distance-based perception is not reliable in a crowded
    coordination sphere, so the user can overrule it -- and the correction has
    to reach the *parameters*, not only the bond list.

    Correcting the list alone left the typing molecule without the drawn bond.
    RDKit then had no entry for it and the exporter fell back to its geometric
    estimate, whose equilibrium is whatever distance the bond happened to be
    drawn at: joining two carbons across a benzene ring gave r0 = 2.798 A with
    k = 111 instead of 1.514 A with k = 700, so the bond never contracted and
    neither carbon changed hybridisation.  Both molecules are therefore
    rebuilt here, and RDKit re-perceives the chemistry from the corrected
    connectivity.

    Args:
        perceived: A :class:`PerceivedMolecule`, mutated in place.
        edits: Mapping of atom-index pairs to whether they are bonded.

    Returns:
        True when the bond list changed.
    """
    if perceived is None or not edits:
        return False

    wanted: Dict[Tuple[int, int], bool] = {}
    for pair, connect in dict(edits).items():
        try:
            i, j = (int(x) for x in pair)
        except Exception:
            continue
        if i == j or min(i, j) < 0 or max(i, j) >= perceived.n_atoms:
            continue
        wanted[(min(i, j), max(i, j))] = bool(connect)

    current = {(min(i, j), max(i, j)) for i, j in perceived.bonds}
    target = set(current)
    for key, connect in wanted.items():
        if connect:
            target.add(key)
        else:
            target.discard(key)
    if target == current:
        return False

    perceived.bonds = sorted(target)
    if not RDKIT_AVAILABLE or perceived.mol is None:
        return True

    added = sorted(target - current)
    removed = sorted(current - target)
    orders = _orders_from_mol(perceived.typing_mol)
    for key in removed:
        orders.pop(key, None)
    for key in added:
        # A drawn bond is a single bond.  Nothing in a click says otherwise,
        # and the relaxation is free to shorten it from there.
        orders[key] = 1

    try:
        with _RDKitQuiet():
            editable = Chem.RWMol(perceived.mol)
            for i, j in removed:
                if editable.GetBondBetweenAtoms(i, j) is not None:
                    editable.RemoveBond(i, j)
            for i, j in added:
                if editable.GetBondBetweenAtoms(i, j) is None:
                    editable.AddBond(i, j, Chem.BondType.SINGLE)
            rebuilt = editable.GetMol()
            rebuilt.UpdatePropertyCache(strict=False)
            fresh: List[str] = []
            typing_mol, _ok = _build_typing_mol(
                rebuilt, perceived.metal_indices, orders, fresh
            )
    except Exception as exc:
        logger.debug('Rebuilding the molecule after a bond edit failed: %s', exc)
        perceived.warnings.append(
            'The edited bond could not be applied to the force-field types; '
            'it is parameterised from the geometry instead.'
        )
        return True

    perceived.mol = rebuilt
    if typing_mol is None and perceived.typing_mol is not None:
        # Keeping the molecule that was typed before the edit is worth more
        # than the edit itself: without one, *every* atom falls back to
        # geometric parameters, while the stale one still types everything
        # the edit did not touch.
        perceived.warnings.append(
            'The edited connectivity could not be sanitized; the changed '
            'bonds are parameterised from the geometry.'
        )
    else:
        perceived.typing_mol = typing_mol
        for message in fresh:
            if message not in perceived.warnings:
                perceived.warnings.append(message)
    return True


#: The hybridisations a user may force on an atom.
HYBRIDISATION_CHOICES: Tuple[str, ...] = ('sp', 'sp2', 'sp3')

#: The angle each one means at that centre.  Needed because setting the
#: hybridisation is not always enough on its own: RDKit picks the UFF type of
#: phosphorus and silicon from valence and charge rather than hybridisation
#: (P stays at 93.8 degrees and Si at 109.5 whatever it is set to), and no
#: angle that touches a metal is typed by RDKit at all -- which is every angle
#: that orients a ligand against its metal.
HYBRIDISATION_ANGLES: Dict[str, float] = {
    'sp': 180.0, 'sp2': 120.0, 'sp3': 109.471,
}

_HYBRIDISATION_TYPES: Dict[str, Any] = {}


def _hybridisation_type(name: str) -> Any:
    """Map ``'sp2'`` to RDKit's enum member, once."""
    if not _HYBRIDISATION_TYPES and RDKIT_AVAILABLE:
        _HYBRIDISATION_TYPES.update({
            'sp': Chem.HybridizationType.SP,
            'sp2': Chem.HybridizationType.SP2,
            'sp3': Chem.HybridizationType.SP3,
        })
    return _HYBRIDISATION_TYPES.get(name)


def perceived_hybridisation_of(perceived: Any, index: int) -> Optional[str]:
    """What perception made of an atom, whether or not an override replaced it.

    The offer has to be able to say what ``automatic`` would mean, and once an
    override is in place the molecule itself no longer knows.
    """
    if perceived is None:
        return None
    recorded = (perceived.auto_hybridisation or {}).get(int(index))
    if recorded is not None:
        return recorded
    return _hybridisation(perceived.typing_mol, index)


def hybridisation_from_connectivity(
    perceived: Any, indices: Optional[Iterable[int]] = None
) -> Dict[int, str]:
    """Read each carbon's hybridisation off how many partners it is bonded to.

    A carbon has no lone pair, so its sigma count fixes its shape outright:
    four partners is tetrahedral, three is trigonal planar, two is linear.
    That is a stronger statement than perception can make, because perception
    goes through bond *orders* -- and those are read from the geometry, so a
    double bond twisted out of plane, or one at an unusual length, is simply
    not seen and its carbon comes back sp3.

    Deliberately carbon only.  Everywhere else a lone pair decides and the
    count cannot see it: nitrogen with three partners is pyramidal in an amine
    and planar in an amide, oxygen with two is bent either way but at quite
    different angles.  Guessing there would trade one wrong answer for another.

    **A coordinated atom raises the obvious question -- does the metal count
    as a partner?**  For carbon it almost always does, and for a reason that
    does not hold for the other donors: carbon has no lone pair to give away
    unless it is a carbene, so an M-C bond is a real sigma bond.  A methyl
    ligand is C plus three substituents and the metal, four, tetrahedral.  An
    N-heterocyclic carbene donates from an sp2 orbital and is two nitrogens
    plus the metal, three, trigonal planar -- both right only *because* the
    metal is counted.  That is exactly the judgement that cannot be made for
    N, O or P, where a dative bond and a covalent one look the same to a
    counter and the lone pair changes the answer; those are left alone.

    The one genuine exception is a side-on alkene, where both carbons are
    bonded to the same metal.  That is a three-membered M-C-C ring and is
    visible without touching bond orders, so it is checked for: the metal is
    then *not* counted, and the carbons come back sp2.  Reality sits between
    the free alkene (sp2) and the metallacyclopropane (sp3), and for the
    weakly bound alkenes this normally means, sp2 is the closer of the two.

    Args:
        perceived: A :class:`PerceivedMolecule`.
        indices: Restrict to these atoms; all of them when omitted.

    Returns:
        Atom index to ``'sp'``/``'sp2'``/``'sp3'``, ready for
        :func:`apply_hybridisation_overrides`.
    """
    if perceived is None:
        return {}
    by_count = {2: 'sp', 3: 'sp2', 4: 'sp3'}
    adjacency = perceived.neighbours()
    wanted = (range(perceived.n_atoms) if indices is None
              else [int(i) for i in indices])
    metals = set(perceived.metal_indices or ())
    derived: Dict[int, str] = {}
    for index in wanted:
        if not 0 <= index < perceived.n_atoms:
            continue
        if index in metals or perceived.symbols[index] != 'C':
            continue
        partners = list(adjacency[index])
        side_on = [
            m for m in partners if m in metals
            and any(other in adjacency[m] for other in partners
                    if other not in metals and perceived.symbols[other] == 'C')
        ]
        count = len(partners) - len(side_on)
        name = by_count.get(count)
        if name is not None:
            derived[index] = name
    return derived


def apply_hybridisation_overrides(perceived: Any, overrides: Any) -> int:
    """Force the hybridisation of individual atoms.

    Hybridisation is perceived from the bond orders, and those are read from
    the geometry, which gets it wrong often enough to matter: a carbon whose
    double bond went unperceived comes out sp3, so its three angles are typed
    at 109.5 degrees and the centre puckers where it should stay flat.

    RDKit's UFF typer reads the atom's hybridisation for most main-group
    elements, so setting it is usually enough -- forcing a carbon to sp2 types
    it ``C_2`` and its angles come back at 120 degrees.  Three angles of 120
    degrees at a three-coordinate centre *is* trigonal planar, which is what
    makes this work without an inversion term.

    Two places where setting it is *not* enough, so the choice is recorded in
    ``forced_hybridisation`` and the exporter builds the angles at that centre
    from it directly:

    * RDKit picks the UFF type of phosphorus and silicon from valence and
      charge, not hybridisation.  Measured: P stays at 93.8 degrees and Si at
      109.5 whatever it is set to, while C, N, O and S all follow.
    * No angle that touches a metal is typed by RDKit at all -- those are
      restrained to the input geometry.  At a donor atom that is every angle
      holding the ligand against its metal, so the shape would not move.

    Must run after :func:`apply_bond_edits`: rebuilding the typing molecule
    sanitizes it, and sanitisation re-perceives hybridisation.

    Args:
        perceived: A :class:`PerceivedMolecule`, mutated in place.
        overrides: Mapping of atom index to ``'sp'``, ``'sp2'`` or ``'sp3'``.

    Returns:
        How many atoms were actually changed.
    """
    if perceived is None or not overrides or not RDKIT_AVAILABLE:
        return 0
    applied = 0
    for raw_index, raw_name in dict(overrides).items():
        try:
            index = int(raw_index)
        except Exception:
            continue
        name = str(raw_name or '').strip().lower()
        if not 0 <= index < perceived.n_atoms or name not in HYBRIDISATION_CHOICES:
            continue
        target = _hybridisation_type(name)
        if target is None:
            continue
        # Recorded whether or not a molecule can carry it: the exporter builds
        # the angles at this centre from the choice directly, and that is the
        # half that works for phosphorus, silicon and every angle at a metal.
        if index not in perceived.auto_hybridisation:
            was = _hybridisation(perceived.typing_mol, index)
            if was is not None:
                perceived.auto_hybridisation[index] = was
        perceived.forced_hybridisation[index] = name
        applied += 1
        for molecule in (perceived.typing_mol, perceived.mol):
            if molecule is None:
                continue
            try:
                atom = molecule.GetAtomWithIdx(index)
            except Exception:
                continue
            # An aromatic flag reads as sp2 wherever it is checked, so a ring
            # carbon forced to sp3 or sp would otherwise ignore the override.
            if name != 'sp2':
                atom.SetIsAromatic(False)
            atom.SetHybridization(target)
    return applied


# --------------------------------------------------------------------------
# UFF torsion periodicity / phase
# --------------------------------------------------------------------------

def _hybridisation(mol: Any, idx: int) -> Optional[str]:
    """Return ``'sp'``/``'sp2'``/``'sp3'`` for an atom, or ``None``."""
    if mol is None:
        return None
    try:
        atom = mol.GetAtomWithIdx(int(idx))
    except Exception:
        return None
    if atom.GetIsAromatic():
        return 'sp2'
    hybrid = atom.GetHybridization()
    name = str(hybrid)
    if name in ('SP3', 'SP3D', 'SP3D2'):
        return 'sp3'
    if name in ('SP2', 'SP2D'):
        return 'sp2'
    if name == 'SP':
        return 'sp'
    return None


def uff_torsion_form(mol: Any, i: int, j: int, k: int, l: int) -> Optional[Tuple[int, float]]:
    """Derive UFF's torsion periodicity ``n`` and phase ``phi0`` for i-j-k-l.

    RDKit hands out only the barrier ``V`` (``GetUFFTorsionParams``), so the
    rest of the functional form has to be reconstructed from UFF's rules
    (Rappe et al. 1992, eq. 16 and the paragraph following it).  ``j-k`` is
    the central bond:

    =====================================  ====  =======  ====================
    case                                   n     phi0     note
    =====================================  ====  =======  ====================
    sp3-sp3                                3     180 deg  general single bond
    sp3-sp3, both group 16, order 1        2      90 deg  the O-O/S-S rule
    sp2-sp2 (incl. aromatic/resonant)      2     180 deg  double/resonant bond
    sp2-sp3, order 1, sp3 is group 16      2      90 deg  e.g. anisole's O
    sp2-sp3, order 1, sp2 end atom is sp2  3     180 deg  the propene rule
    sp2-sp3 otherwise                      6       0 deg  e.g. an isolated CH3
    anything involving sp                  --    --       no torsion term
    =====================================  ====  =======  ====================

    Approximations and subtleties, honestly stated:

    * "sp2" also covers UFF's resonant (aromatic) type.  UFF gives both the
      same ``n``/``phi0`` and differs only in ``V``, which RDKit supplies.
    * The conjugation ("propene") exception is decided on the *end* atom
      attached to the sp2 centre, not on the sp2 centre's neighbourhood in
      general.  That is what RDKit does, and the difference is observable:
      in acetone RDKit returns V = 2 kcal/mol for H-C-C=O and V = 1 kcal/mol
      for H-C-C-C about the very same central bond.  A per-bond rule would
      therefore pair the wrong periodicity with one of those barriers.
    * Hybridisation comes from RDKit's perception on the metal-stripped
      molecule.  If bond-order perception failed, everything looks sp3 and
      the exceptions above cannot trigger.
    * Bonds to metals never get a torsion term (see
      :func:`export_forcefield_terms`).

    Returns:
        ``(n, phi0_degrees)``, or ``None`` when UFF assigns no torsion here.
    """
    hyb_j = _hybridisation(mol, j)
    hyb_k = _hybridisation(mol, k)
    if hyb_j is None or hyb_k is None:
        return None
    if hyb_j == 'sp' or hyb_k == 'sp':
        return None

    bond = mol.GetBondBetweenAtoms(int(j), int(k))
    if bond is None:
        return None
    try:
        order = bond.GetBondTypeAsDouble()
    except Exception:
        order = 1.0
    is_single = abs(order - 1.0) < 1e-6

    atom_j = mol.GetAtomWithIdx(int(j))
    atom_k = mol.GetAtomWithIdx(int(k))

    if hyb_j == 'sp3' and hyb_k == 'sp3':
        if is_single and atom_j.GetAtomicNum() in _GROUP16 and atom_k.GetAtomicNum() in _GROUP16:
            return (2, 90.0)
        return (3, 180.0)

    if hyb_j == 'sp2' and hyb_k == 'sp2':
        return (2, 180.0)

    # Mixed sp2 / sp3: the end atom bonded to the sp2 centre decides.
    if hyb_j == 'sp3':
        sp3_atom, end_at_sp2 = atom_j, l
    else:
        sp3_atom, end_at_sp2 = atom_k, i
    if is_single:
        if sp3_atom.GetAtomicNum() in _GROUP16:
            return (2, 90.0)
        if _hybridisation(mol, end_at_sp2) == 'sp2':
            return (3, 180.0)
    return (6, 0.0)


def _angle_periodicity(theta0_deg: float) -> int:
    """Return UFF's angle periodicity for an ideal ``theta0`` (0 = general)."""
    for ideal, order in _ANGLE_PERIODICITY:
        if abs(theta0_deg - ideal) < _ANGLE_IDEAL_TOL:
            return order
    return 0


# --------------------------------------------------------------------------
# term export
# --------------------------------------------------------------------------

def _distance(a: Sequence[float], b: Sequence[float]) -> float:
    return math.sqrt(
        (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
    )


def _angle_degrees(a: Sequence[float], b: Sequence[float], c: Sequence[float]) -> float:
    v1 = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    v2 = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
    n1 = math.sqrt(sum(v * v for v in v1))
    n2 = math.sqrt(sum(v * v for v in v2))
    if n1 <= 1e-9 or n2 <= 1e-9:
        return 109.47
    cos_theta = sum(x * y for x, y in zip(v1, v2)) / (n1 * n2)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_theta))))


def _empty_payload(source: str, warnings: Sequence[str]) -> Dict[str, Any]:
    return {
        'ok': False,
        'source': source,
        'n_atoms': 0,
        'elements': [],
        'bonds': [],
        'angles': [],
        'torsions': [],
        'vdw': [],
        'metals': [],
        'warnings': list(warnings),
        'version': PAYLOAD_VERSION,
    }


METHODS = ('uff', 'mmff94')


def polyhedron_options(perceived, metal_index):
    """Which coordination polyhedra the donors around this metal could take.

    Returns ``(coordination_number, current_geometry, [(code, label), ...])``,
    or ``None`` when the atom is not a metal or its coordination number is
    outside the catalogue. The geometry tables come from
    :mod:`delfin.manta._polyhedron_targets`, the same ones MANTA builds
    complexes with.
    """
    try:
        from delfin.manta import _polyhedron_targets as targets
    except Exception:
        return None
    if perceived is None:
        return None
    metal_index = int(metal_index)
    if metal_index not in set(perceived.metal_indices or ()):
        return None
    donors = sorted(
        j for pair in perceived.bonds for j in pair
        if metal_index in pair and j != metal_index
    )
    cn = len(donors)
    if not 2 <= cn <= 9:
        return None

    available = []
    for key, vectors in targets._IDEAL_VECTORS.items():
        if vectors.shape[0] != cn:
            continue
        if key not in available:
            available.append(key)
    if not available:
        return None

    # Which polyhedron the complex actually sits closest to. The CN-based
    # classifier answers what a given coordination number *usually* is, which
    # for CN=4 is always tetrahedral -- it would have labelled a clearly
    # square-planar Ni complex 'tetrahedral'. Measure instead.
    current = _closest_polyhedron(perceived, donors, metal_index, available)
    labels = {
        'linear_2': 'linear', 'bent_2': 'bent',
        'trigonal_planar': 'trigonal planar',
        'Td': 'tetrahedral', 'sqp_4': 'square planar', 'see_saw': 'see-saw',
        'tbp': 'trigonal bipyramidal', 'sqp_5': 'square pyramidal',
        'Oh': 'octahedral', 'trig_prism': 'trigonal prismatic',
        'pbp': 'pentagonal bipyramidal', 'capped_oct': 'capped octahedral',
        'sq_antiprism': 'square antiprismatic', 'cube': 'cubic',
        'dodecahedron': 'dodecahedral',
        'bicapped_trig_antiprism': 'bicapped trigonal antiprismatic',
        'capped_sap': 'capped square antiprismatic',
        'tricapped_tp': 'tricapped trigonal prismatic',
    }
    options = sorted(
        ((key, labels.get(key, key)) for key in available),
        key=lambda pair: pair[1],
    )
    return cn, current, options


def _pairwise_angles(vectors):
    """Sorted L-M-L angles of a set of unit vectors, in degrees.

    Rotation invariant, which is the whole point: a polyhedron has to be
    recognised whatever way round the molecule happens to lie.
    """
    import numpy as np

    out = []
    for a in range(len(vectors)):
        for b in range(a + 1, len(vectors)):
            cosine = float(np.clip(np.dot(vectors[a], vectors[b]), -1.0, 1.0))
            out.append(math.degrees(math.acos(cosine)))
    return sorted(out)


def _unit_donor_vectors(perceived, donors, metal_index, coords=None):
    import numpy as np

    coords = coords if coords is not None else perceived.coords
    centre = np.asarray(coords[metal_index], dtype=float)
    vectors = np.asarray([coords[j] for j in donors], dtype=float) - centre
    lengths = np.linalg.norm(vectors, axis=1)
    if not len(lengths) or float(lengths.min()) < 1e-9:
        return None
    return vectors / lengths[:, None]


def _closest_polyhedron(perceived, donors, metal_index, candidates):
    """Name the candidate polyhedron the donors currently sit closest to.

    Compared through the sorted set of L-M-L angles, so the answer does not
    depend on how the molecule happens to be oriented.
    """
    import numpy as np

    from delfin.manta import _polyhedron_targets as targets

    observed = _unit_donor_vectors(perceived, donors, metal_index)
    if observed is None:
        return None
    measured = _pairwise_angles(observed)

    best, best_cost = None, None
    for key in candidates:
        try:
            ideal = np.asarray(
                targets.get_ideal_donor_vectors(len(donors), key), dtype=float,
            )
        except Exception:
            continue
        wanted = _pairwise_angles(ideal)
        cost = float(np.mean(np.abs(np.array(measured) - np.array(wanted))))
        if best_cost is None or cost < best_cost:
            best, best_cost = key, cost
    return best


def polyhedron_assignment(perceived, metal_index, geometry, coords=None):
    """Which donor currently sits on which vertex of the chosen polyhedron.

    Returned as ``{donor_index: vertex_number}``. Swapping two entries and
    handing the result back as ``polyhedron['assignment']`` is what lets a user
    exchange two ligands: without it the assignment is recomputed on every
    export and always picks the nearest match, which is exactly the
    arrangement they are trying to leave.
    """
    _targets, assigned = _assign_donors_to_vertices(
        perceived, coords if coords is not None else perceived.coords,
        metal_index, geometry, None,
    )
    return assigned


def polyhedron_vertex_classes(coordination, geometry):
    """Group a polyhedron's vertices into the kinds that are not equivalent.

    A trigonal bipyramid has two: three equatorial vertices at 120 degrees to
    each other, and two axial ones facing each other at 180.  Which ligand
    ends up in which kind is a real chemical choice -- ``PF5`` is one molecule
    but Berry pseudorotation moves ligands between the two -- and the
    assignment made by matching the polyhedron onto the geometry as it stands
    only ever finds the nearest one.

    Vertices are compared through the sorted angles from each one to all the
    others, which does not depend on how the polyhedron is oriented.

    Returns:
        ``(classes, labels)``: a class index per vertex, and a readable name
        per class ordered by how many vertices it holds.  ``None`` when the
        geometry is not in the catalogue.
    """
    import numpy as np

    try:
        from delfin.manta import _polyhedron_targets as targets
        ideal = np.asarray(
            targets.get_ideal_donor_vectors(int(coordination), geometry),
            dtype=float,
        )
    except Exception:
        return None
    count = ideal.shape[0]
    signatures = []
    for i in range(count):
        angles = []
        for j in range(count):
            if i == j:
                continue
            cosine = float(np.clip(float(np.dot(ideal[i], ideal[j])), -1.0, 1.0))
            angles.append(round(math.degrees(math.acos(cosine)), 1))
        signatures.append(tuple(sorted(angles)))

    order: List[Any] = []
    for signature in signatures:
        if signature not in order:
            order.append(signature)
    # Smallest group first, because that is the distinguished one people name:
    # the two axial vertices of a bipyramid, the single apex of a pyramid.
    order.sort(key=lambda s: (signatures.count(s), s))
    classes = [order.index(s) for s in signatures]

    sizes = [signatures.count(s) for s in order]
    if len(order) == 1:
        labels = ['all equivalent']
    elif len(order) == 2 and sorted(sizes) == [1, count - 1]:
        labels = ['apical', 'basal'] if sizes[0] == 1 else ['basal', 'apical']
    elif len(order) == 2:
        labels = ['axial', 'equatorial'] if sizes[0] < sizes[1] \
            else ['equatorial', 'axial']
    else:
        labels = [f'set {n + 1}' for n in range(len(order))]
    return classes, labels


def polyhedron_arrangements(perceived, metal_index, geometry, coords=None):
    """Every distinct way of spreading the donors over the vertex kinds.

    A trigonal bipyramid with five different ligands can be built ten ways --
    which two are axial -- and they are different molecules, not different
    views of one.  This enumerates them, best fit first, each as an
    ``{donor: vertex}`` mapping ready to hand back as the polyhedron's
    ``assignment``.

    Only the split between *kinds* is enumerated.  Which of the three
    equatorial vertices a given ligand takes is left to the usual nearest
    match, so choosing an arrangement never moves a ligand further than it
    has to.  On a polyhedron whose vertices are all equivalent -- an
    octahedron, a tetrahedron -- there is one arrangement and turning it
    changes nothing; exchanging two ligands there is what Swap is for.
    """
    import numpy as np
    from itertools import combinations

    from delfin.manta import _polyhedron_targets as targets

    if perceived is None:
        return []
    metal_index = int(metal_index)
    donors = sorted(
        j for pair in perceived.bonds for j in pair
        if metal_index in pair and j != metal_index
    )
    cn = len(donors)
    grouped = polyhedron_vertex_classes(cn, geometry)
    if not grouped or cn < 2:
        return []
    classes, labels = grouped
    if len(set(classes)) < 2:
        return [{}]

    coords = coords if coords is not None else perceived.coords
    observed = _unit_donor_vectors(perceived, donors, metal_index, coords)
    if observed is None:
        return []
    try:
        ideal = np.asarray(
            targets.get_ideal_donor_vectors(cn, geometry), dtype=float)
    except Exception:
        return []

    try:
        from scipy.optimize import linear_sum_assignment
    except Exception:
        linear_sum_assignment = None

    by_class: Dict[int, List[int]] = {}
    for vertex, klass in enumerate(classes):
        by_class.setdefault(klass, []).append(vertex)
    ordered_classes = sorted(by_class)

    def _best(fixed_pairs, frame):
        """Kabsch the polyhedron onto this split, then score it."""
        rows = [donors.index(d) for d, _v in fixed_pairs]
        cols = [v for _d, v in fixed_pairs]
        for _pass in range(3):
            matched = frame[np.asarray(cols)]
            target = observed[np.asarray(rows)]
            u, _sv, vt = np.linalg.svd(matched.T @ target)
            rotation = u @ vt
            if np.linalg.det(rotation) < 0:
                u[:, -1] *= -1.0
                rotation = u @ vt
            frame = frame @ rotation
            # Re-match inside each kind only: the split itself is the choice.
            cols = list(cols)
            for klass in ordered_classes:
                members = by_class[klass]
                mine = [n for n, v in enumerate(cols) if v in members]
                if len(mine) < 2:
                    continue
                sub = 1.0 - observed[np.asarray([rows[n] for n in mine])] \
                    @ frame[np.asarray(members)].T
                if linear_sum_assignment is not None:
                    r, c = linear_sum_assignment(sub)
                else:
                    r, c = range(len(mine)), range(len(mine))
                for a, b in zip(r, c):
                    cols[mine[int(a)]] = members[int(b)]
        score = float(np.sum(1.0 - np.einsum(
            'ij,ij->i', observed[np.asarray(rows)], frame[np.asarray(cols)])))
        return score, {donors[rows[n]]: int(cols[n]) for n in range(cn)}

    # Enumerate the splits: choose the members of each kind in turn.
    splits: List[List[Tuple[int, int]]] = []

    def _walk(remaining, index, chosen):
        if index == len(ordered_classes):
            splits.append(list(chosen))
            return
        klass = ordered_classes[index]
        members = by_class[klass]
        for pick in combinations(remaining, len(members)):
            rest = [d for d in remaining if d not in pick]
            _walk(rest, index + 1,
                  chosen + [(d, v) for d, v in zip(pick, members)])

    _walk(donors, 0, [])
    scored = []
    for split in splits:
        try:
            scored.append(_best(split, np.array(ideal, copy=True)))
        except Exception:
            continue
    scored.sort(key=lambda pair: pair[0])
    out, seen = [], set()
    for _score, mapping in scored:
        key = tuple(sorted((d, classes[v]) for d, v in mapping.items()))
        if key in seen:
            continue
        seen.add(key)
        out.append(mapping)
    return out


def describe_polyhedron_arrangement(perceived, geometry, assignment):
    """Name an arrangement by which donors sit on the distinguished vertices."""
    if not assignment:
        return ''
    grouped = polyhedron_vertex_classes(len(assignment), geometry)
    if not grouped:
        return ''
    classes, labels = grouped
    if len(set(classes)) < 2:
        return 'all vertices equivalent'
    members: Dict[int, List[int]] = {}
    for donor, vertex in assignment.items():
        members.setdefault(classes[int(vertex)], []).append(int(donor))
    klass = min(members)
    named = ', '.join(
        f'{perceived.symbols[d]}{d}' for d in sorted(members[klass])
    )
    return f'{labels[klass]}: {named}'


def _assign_donors_to_vertices(perceived, coords, metal_index, geometry, fixed):
    """Match donors to polyhedron vertices, honouring a caller's swap."""
    import numpy as np

    from delfin.manta import _polyhedron_targets as targets

    donors = sorted(
        j for pair in perceived.bonds for j in pair
        if metal_index in pair and j != metal_index
    )
    cn = len(donors)
    ideal = np.asarray(
        targets.get_ideal_donor_vectors(cn, geometry), dtype=float,
    )
    observed = _unit_donor_vectors(perceived, donors, metal_index, coords)
    if observed is None:
        return {}, {}

    if fixed:
        # The caller says which donor belongs on which vertex; only the
        # orientation of the polyhedron still has to be found.
        pairs = [
            (donors.index(int(d)), int(v))
            for d, v in fixed.items()
            if int(d) in donors and 0 <= int(v) < cn
        ]
        if len(pairs) == cn:
            rows = [r for r, _c in pairs]
            cols = [c for _r, c in pairs]
        else:
            fixed = None
    if not fixed:
        try:
            from scipy.optimize import linear_sum_assignment
        except Exception:
            linear_sum_assignment = None
        rows, cols = list(range(cn)), list(range(cn))
        for _pass in range(3):
            cost = 1.0 - observed @ ideal.T
            if linear_sum_assignment is not None:
                rows, cols = linear_sum_assignment(cost)
            matched = ideal[np.asarray(cols)]
            target = observed[np.asarray(rows)]
            u, _sv, vt = np.linalg.svd(matched.T @ target)
            rotation = u @ vt
            if np.linalg.det(rotation) < 0:
                u[:, -1] *= -1.0
                rotation = u @ vt
            ideal = ideal @ rotation

    assigned_vec = {}
    assigned_vertex = {}
    for row, col in zip(rows, cols):
        assigned_vec[donors[int(row)]] = np.asarray(ideal[int(col)], dtype=float)
        assigned_vertex[donors[int(row)]] = int(col)

    wanted = {}
    for a in range(cn):
        for b in range(a + 1, cn):
            i, k = donors[a], donors[b]
            if i not in assigned_vec or k not in assigned_vec:
                continue
            cosine = float(np.clip(np.dot(assigned_vec[i], assigned_vec[k]), -1.0, 1.0))
            wanted[(min(i, k), max(i, k))] = math.degrees(math.acos(cosine))
    return wanted, assigned_vertex


def _polyhedron_target_angles(perceived, coords, metal_index, geometry, fixed=None):
    """Ideal L-M-L angles for one metal forced onto a polyhedron.

    The donors are matched to the polyhedron's vertices -- by Hungarian
    assignment after the polyhedron has been rotated onto them, so the complex
    is drawn onto the target the shortest way round, or by the caller's own
    mapping when two ligands are being exchanged.

    These become the equilibrium values of ordinary harmonic angle terms with
    UFF force constants -- a pull, not a constraint. A chelate whose bite angle
    cannot reach the ideal vertex separation settles at a compromise between
    ligand and polyhedron rather than being torn apart.
    """
    wanted, _assigned = _assign_donors_to_vertices(
        perceived, coords, metal_index, geometry, fixed,
    )
    return wanted


def normalise_method(method: Optional[str]) -> str:
    """Validate a force-field choice, falling back to UFF."""
    name = str(method or 'uff').strip().lower().replace('-', '').replace('_', '')
    if name in ('mmff', 'mmff94', 'mmff94s'):
        return 'mmff94'
    return 'uff'


#: Force constants for user restraints, chosen to sit in the same range as the
#: real UFF terms (a C-C stretch is around 700 kcal/mol/A^2, an angle bend
#: 100-300 kcal/mol/rad^2).  A restraint that dominated everything would hold
#: its value by tearing the rest of the structure, which is the opposite of
#: what a constraint is for.
RESTRAINT_FORCE_CONSTANTS = {
    'distance': 500.0,
    'angle': 200.0,
    'dihedral': 50.0,
}


def build_restraints(entries) -> List[Dict[str, Any]]:
    """Turn the user's held values into force-field restraint terms."""
    out: List[Dict[str, Any]] = []
    for entry in entries or ():
        kind = str(entry.get('kind') or '')
        if kind not in RESTRAINT_FORCE_CONSTANTS:
            continue
        atoms = [int(a) for a in entry.get('atoms') or ()]
        needed = {'distance': 2, 'angle': 3, 'dihedral': 4}[kind]
        if len(atoms) != needed or len(set(atoms)) != needed:
            continue
        try:
            value = float(entry.get('value'))
        except (TypeError, ValueError):
            continue
        out.append({
            'kind': kind,
            'atoms': atoms,
            'value': value,
            'k': float(entry.get('k') or RESTRAINT_FORCE_CONSTANTS[kind]),
        })
    return out


def export_forcefield_terms(
    xyz_text: str,
    *,
    perceived: Optional[PerceivedMolecule] = None,
    method: Optional[str] = 'uff',
    polyhedron: Optional[Dict[str, Any]] = None,
    restraints=None,
) -> Dict[str, Any]:
    """Build the JSON payload of force-field terms for the browser engine.

    Called once, when manipulate mode is switched on.  See the module
    docstring for the exact payload shape, units and energy expressions --
    the browser engine is written against that contract.

    Parameter sources, in the order they are tried per term:

    1. **RDKit UFF** (``source = 'rdkit-uff'``) whenever RDKit can type every
       atom involved.  These are the real published UFF parameters.
    2. **UFF table + input geometry** (``source = 'geometric-fallback'``) for
       every term RDKit refuses.  In practice that means *metals*: RDKit's
       UFF silently ignores transition metals -- ``UFFHasAllMoleculeParams``
       is False, ``GetUFFBondStretchParams(Zn, N)`` is ``None``, and the
       metal ends up in no bonded and no van-der-Waals term at all, so
       moving it costs exactly zero energy.  That trap is closed here: the
       metal's bonds and angles are restrained at their **observed** length
       and angle with force constants from UFF's own formulas (using the
       effective charges Z* from Open Babel's ``UFF.prm``), and its
       van-der-Waals parameters are read from the same table.  A warning
       naming the metal is always attached.

    Args:
        xyz_text: Standard XYZ text or a bare DELFIN coordinate block.
        perceived: Optional pre-computed :func:`perceive_molecule` result, to
            avoid perceiving the same structure twice.  When given,
            ``xyz_text`` is not read at all.

    Returns:
        The payload dict.  On failure ``{'ok': False, ...}`` with empty term
        lists and an explanatory warning -- never an exception.
    """
    if perceived is None:
        perceived = perceive_molecule(xyz_text)
    if perceived is None:
        reason = ('RDKit is not available.' if not RDKIT_AVAILABLE
                  else 'The coordinates could not be parsed.')
        return _empty_payload(SOURCE_GEOMETRIC, [reason])

    warnings: List[str] = list(perceived.warnings)
    symbols = perceived.symbols
    coords = perceived.coords
    # A perception handed in by the caller is cached deliberately -- the
    # bonding must not be re-read from a geometry the user has been dragging,
    # or a twisted double bond stops being one. Its *coordinates* are another
    # matter: taking those from the cache too meant every geometry-derived
    # value was computed from the structure as it was when the field was first
    # switched on. A ligand dragged to another vertex was therefore assigned
    # the vertex it used to be nearest, which is why exchanging two of them by
    # dragging never took.
    parsed = parse_xyz(xyz_text)
    if parsed is not None:
        fresh_symbols, fresh_coords, _had_header = parsed
        if list(fresh_symbols) == list(symbols):
            coords = fresh_coords
    n_atoms = perceived.n_atoms
    metals = set(perceived.metal_indices)
    typing_mol = perceived.typing_mol
    adjacency = perceived.neighbours()

    # An atom whose hybridisation the user forced takes its angles from that
    # choice, not from RDKit.  Setting the hybridisation alone reaches neither
    # phosphorus nor silicon -- RDKit types those from valence and charge, so
    # P stays at 93.8 degrees and Si at 109.5 whatever it is set to -- nor any
    # angle that touches a metal, which is never typed at all.  At a donor
    # atom that is every angle holding the ligand against its metal, so
    # forcing a hybridisation there did nothing whatsoever.
    forced_hyb = dict(getattr(perceived, 'forced_hybridisation', None) or {})
    if forced_hyb:
        listed = ', '.join(
            f'{symbols[i]}{i} as {forced_hyb[i]}' for i in sorted(forced_hyb)
        )
        warnings.append(
            f'Hybridisation forced by hand: {listed}. Every angle at those '
            'atoms is built from that choice rather than from perception, so '
            'it also holds where the atom is bonded to a metal.'
        )

    if perceived.has_metal:
        named = ', '.join(
            f'{symbols[i]} (atom {i})' for i in sorted(metals)
        )
        warnings.append(
            f'RDKit UFF has no parameters for {named}; the bonds and angles at '
            'the metal are restrained to the input geometry with UFF force '
            'constants, and its van-der-Waals parameters come from the UFF '
            'table shipped with Open Babel.'
        )
    if n_atoms > LARGE_MOLECULE_ATOMS:
        warnings.append(
            f'{n_atoms} atoms: the browser force field will not hold 30 fps '
            'at this size; expect fewer relaxation steps per frame.'
        )

    # Per-atom UFF parameter lookups from the fallback table, cached by
    # (element, coordination number).
    prm_cache: Dict[Tuple[str, int], Optional[Dict[str, float]]] = {}

    def _prm(index: int) -> Optional[Dict[str, float]]:
        key = (symbols[index], len(adjacency[index]))
        if key not in prm_cache:
            prm_cache[key] = uff_atom_parameters(key[0], key[1])
        return prm_cache[key]

    def _effective_charge(index: int) -> float:
        entry = _prm(index)
        if entry and entry.get('Z1'):
            return float(entry['Z1'])
        return 2.0  # a mid-range Z*, only reached without the UFF table

    used_fallback = False
    fallback_elements: set = set()
    missing_vdw: List[str] = []

    bonds: List[Dict[str, Any]] = []
    angles: List[Dict[str, Any]] = []
    torsions: List[Dict[str, Any]] = []
    vdw: List[Dict[str, float]] = []

    with _RDKitQuiet():
        # ---- bonds ----------------------------------------------------
        for i, j in perceived.bonds:
            params = None
            if typing_mol is not None and i not in metals and j not in metals:
                try:
                    params = _uff_helpers.GetUFFBondStretchParams(typing_mol, i, j)
                except Exception:
                    params = None
            if params is not None:
                k_bond, r0 = float(params[0]), float(params[1])
            else:
                used_fallback = True
                if i not in metals and j not in metals:
                    fallback_elements.update((symbols[i], symbols[j]))
                r0 = _distance(coords[i], coords[j])
                k_bond = _uff_bond_force_constant(
                    _effective_charge(i), _effective_charge(j), r0
                )
            bonds.append({'i': i, 'j': j, 'k': round(k_bond, 4), 'r0': round(r0, 5)})

        # ---- angles ---------------------------------------------------
        # A polyhedron the user asked for replaces the observed angles at that
        # one metal with the ideal ones. They stay ordinary harmonic terms with
        # UFF force constants, so the complex is pulled towards the shape and a
        # chelate that cannot reach it settles at a compromise.
        forced_angles = {}
        forced_centre = None
        if polyhedron:
            try:
                forced_centre = int(polyhedron.get('metal'))
                forced_angles = _polyhedron_target_angles(
                    perceived, coords, forced_centre,
                    str(polyhedron.get('geometry')),
                    polyhedron.get('assignment'),
                )
            except Exception as exc:
                warnings.append(f'Could not apply the polyhedron: {exc}')
                forced_angles, forced_centre = {}, None

        for centre in range(n_atoms):
            neighbours = adjacency[centre]
            forced_theta = HYBRIDISATION_ANGLES.get(forced_hyb.get(centre, ''))
            for a in range(len(neighbours)):
                for b in range(a + 1, len(neighbours)):
                    i, k = neighbours[a], neighbours[b]
                    params = None
                    if (forced_theta is None and typing_mol is not None
                            and centre not in metals and i not in metals and k not in metals):
                        try:
                            params = _uff_helpers.GetUFFAngleBendParams(typing_mol, i, centre, k)
                        except Exception:
                            params = None
                    if params is not None:
                        k_theta, theta0 = float(params[0]), float(params[1])
                        periodicity = _angle_periodicity(theta0)
                    else:
                        if forced_theta is None:
                            used_fallback = True
                            if not ({i, centre, k} & metals):
                                fallback_elements.update(
                                    (symbols[i], symbols[centre], symbols[k])
                                )
                        theta0 = (
                            forced_theta if forced_theta is not None
                            else _angle_degrees(coords[i], coords[centre], coords[k])
                        )
                        if centre == forced_centre:
                            theta0 = forced_angles.get(
                                (min(i, k), max(i, k)), theta0,
                            )
                        r_ij = _distance(coords[i], coords[centre])
                        r_jk = _distance(coords[k], coords[centre])
                        k_theta = _uff_angle_force_constant(
                            _effective_charge(i), _effective_charge(k),
                            r_ij, r_jk, math.radians(theta0),
                        )
                        # The general cosine expansion diverges as theta0
                        # approaches 180 deg (C2 = 1/(4 sin^2 theta0)), which
                        # is exactly the trans-donor case in an octahedron:
                        # use UFF's linear form there instead.
                        periodicity = 1 if theta0 >= _LINEAR_ANGLE_CUTOFF else 0
                    angles.append({
                        'i': i, 'j': centre, 'k': k,
                        'kt': round(k_theta, 4),
                        'theta0': round(theta0, 4),
                        'n': periodicity,
                    })

        # ---- torsions -------------------------------------------------
        # UFF assigns no torsional barrier to coordination bonds, and RDKit
        # cannot type a metal in any case, so every torsion that touches a
        # metal is dropped.  The metal's orientation is held by its angle
        # terms instead.
        if typing_mol is not None:
            for j, k in perceived.bonds:
                if j in metals or k in metals:
                    continue
                nbrs_j = [x for x in adjacency[j] if x != k and x not in metals]
                nbrs_k = [x for x in adjacency[k] if x != j and x not in metals]
                if not nbrs_j or not nbrs_k:
                    continue
                n_torsions = len(nbrs_j) * len(nbrs_k)
                for i in nbrs_j:
                    for l in nbrs_k:
                        if i == l:
                            continue  # three-membered ring
                        form = uff_torsion_form(typing_mol, i, j, k, l)
                        if form is None:
                            continue
                        periodicity, phi0 = form
                        try:
                            barrier = _uff_helpers.GetUFFTorsionParams(typing_mol, i, j, k, l)
                        except Exception:
                            barrier = None
                        if not barrier:
                            continue
                        # RDKit's V is the raw barrier: the factor 1/2 and the
                        # division by the number of torsions about the central
                        # bond are applied inside its contrib, so they are
                        # folded into the exported value here.
                        value = 0.5 * float(barrier) / n_torsions
                        if abs(value) < 1e-9:
                            continue
                        torsions.append({
                            'i': i, 'j': j, 'k': k, 'l': l,
                            'v': round(value, 6),
                            'n': periodicity,
                            'phi0': phi0,
                        })

        # ---- van der Waals (per atom, geometric-mean combining) --------
        for index in range(n_atoms):
            params = None
            if typing_mol is not None and index not in metals:
                try:
                    params = _uff_helpers.GetUFFVdWParams(typing_mol, index, index)
                except Exception:
                    params = None
            if params is not None:
                x_i, d_i = float(params[0]), float(params[1])
            else:
                used_fallback = True
                fallback_elements.add(symbols[index])
                entry = _prm(index)
                if entry and entry.get('x1'):
                    x_i, d_i = float(entry['x1']), float(entry.get('D1', 0.05))
                else:
                    # Last resort: a soft, weakly attractive sphere sized from
                    # the covalent radius.  Only reached when neither RDKit nor
                    # the UFF table knows the element.
                    missing_vdw.append(symbols[index])
                    radius = 1.5
                    if RDKIT_AVAILABLE:
                        try:
                            radius = Chem.GetPeriodicTable().GetRcovalent(symbols[index])
                        except Exception:
                            radius = 1.5
                    x_i, d_i = 2.0 * max(0.6, radius), 0.05
            vdw.append({'x': round(x_i, 4), 'd': round(d_i, 5)})

    if missing_vdw:
        warnings.append(
            'No UFF van-der-Waals parameters for '
            + ', '.join(sorted(set(missing_vdw)))
            + '; a covalent-radius estimate was used for those atoms.'
        )

    other_fallback = sorted(fallback_elements - {symbols[i] for i in metals})
    if other_fallback:
        warnings.append(
            'RDKit UFF could not type ' + ', '.join(other_fallback)
            + '; the terms involving those atoms were built from the input '
            'geometry with parameters from the UFF table.'
        )

    source = SOURCE_GEOMETRIC if (used_fallback or typing_mol is None) else SOURCE_RDKIT
    if not bonds and n_atoms > 1:
        warnings.append('No bonds were perceived; only van-der-Waals terms are exported.')

    if polyhedron and forced_angles:
        warnings.append(
            'Atom {} ({}) is being pulled towards a {} polyhedron: the angle '
            'terms at the metal carry its ideal values instead of the measured '
            'ones. They stay ordinary harmonic terms, so a ligand that cannot '
            'reach the ideal settles at a compromise rather than being '
            'strained apart.'.format(
                forced_centre, symbols[forced_centre],
                str(polyhedron.get('geometry')),
            )
        )

    chosen = normalise_method(method)
    if chosen == 'mmff94':
        # MMFF94 is not a drop-in for UFF here. Its bond stretch is quartic,
        # its angle bend cubic-corrected, its torsion a three-term Fourier
        # series and its van der Waals a buffered 14-7 -- none of which are the
        # functional forms the browser engine evaluates. Silently relabelling
        # UFF terms as MMFF94 would be a fabrication, so the live drag keeps
        # UFF-shaped terms and MMFF94 governs the clean-up minimisation on
        # release, where RDKit evaluates the real thing.
        warnings.append(
            'MMFF94 selected: the interactive drag uses UFF-shaped terms, and '
            'MMFF94 is applied by the relaxation that runs when you let go.'
        )

    payload: Dict[str, Any] = {
        'ok': True,
        'source': source,
        'method': chosen,
        'restraints': build_restraints(restraints),
        'polyhedron': (
            {'metal': forced_centre, 'geometry': str(polyhedron.get('geometry'))}
            if (polyhedron and forced_angles) else None
        ),
        'n_atoms': n_atoms,
        'elements': list(symbols),
        'bonds': bonds,
        'angles': angles,
        'torsions': torsions,
        'vdw': vdw,
        'metals': [{'index': i, 'element': symbols[i]} for i in sorted(metals)],
        'warnings': warnings,
        'version': PAYLOAD_VERSION,
    }

    # Contract check: the browser maps these indices straight onto 3Dmol
    # atoms, so an out-of-range index would be a silent, hard-to-debug bug.
    if not _payload_indices_valid(payload):
        payload['ok'] = False
        payload['warnings'].append(
            'Internal consistency check failed: term indices do not match the '
            'atom order; the payload was rejected.'
        )
    return payload


def _payload_indices_valid(payload: Dict[str, Any]) -> bool:
    """Return True when every index in ``payload`` is a valid atom index."""
    n = int(payload.get('n_atoms', 0))
    if len(payload.get('elements', [])) != n or len(payload.get('vdw', [])) != n:
        return False
    for term, keys in (
        ('bonds', ('i', 'j')),
        ('angles', ('i', 'j', 'k')),
        ('torsions', ('i', 'j', 'k', 'l')),
    ):
        for entry in payload.get(term, []):
            for key in keys:
                index = entry.get(key)
                if not isinstance(index, int) or index < 0 or index >= n:
                    return False
    return True


# --------------------------------------------------------------------------
# clean-up relaxation
# --------------------------------------------------------------------------

def _pairs_closer_than(
    positions: Sequence[Sequence[float]], cutoff: float
) -> List[Tuple[int, int]]:
    """Return all atom pairs closer than ``cutoff``, via a cell list.

    A quadratic scan would cost 77k distance evaluations at 392 atoms on
    every mouse release; bucketing into ``cutoff``-sized cells makes it
    linear in the atom count.
    """
    if cutoff <= 0.0:
        return []
    cells: Dict[Tuple[int, int, int], List[int]] = {}
    inverse = 1.0 / cutoff
    for index, point in enumerate(positions):
        key = (
            int(math.floor(point[0] * inverse)),
            int(math.floor(point[1] * inverse)),
            int(math.floor(point[2] * inverse)),
        )
        cells.setdefault(key, []).append(index)

    pairs: List[Tuple[int, int]] = []
    for (cx, cy, cz), members in cells.items():
        candidates: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    candidates.extend(cells.get((cx + dx, cy + dy, cz + dz), ()))
        for i in members:
            for j in candidates:
                if j <= i:
                    continue
                if _distance(positions[i], positions[j]) < cutoff:
                    pairs.append((i, j))
    return pairs


def _relieve_hard_clashes(
    coords: List[Tuple[float, float, float]],
    fixed: Sequence[int],
    threshold: float = 0.7,
    target: float = 0.9,
) -> Tuple[List[Tuple[float, float, float]], int]:
    """Push apart atom pairs that are absurdly close, before minimising.

    A drag can leave two atoms almost on top of each other.  Both engines
    cope badly with that: UFF's r^-12 term is then of order 1e12 kcal/mol,
    and the repo's Open Babel wrapper treats such an energy as the marker of
    an unparameterised metal and returns the input geometry untouched.  A
    single geometric push (fixed atoms never move) costs nothing and lets the
    real minimisation do its job.

    Returns the (possibly corrected) coordinates and the number of pairs
    that had to be separated.
    """
    frozen = set(int(i) for i in fixed)
    positions = [list(map(float, xyz)) for xyz in coords]
    relieved = 0
    for _ in range(4):
        pairs = _pairs_closer_than(positions, threshold)
        if not pairs:
            break
        for i, j in pairs:
            if i in frozen and j in frozen:
                continue
            d = _distance(positions[i], positions[j])
            if d >= threshold:
                continue
            relieved += 1
            if d < 1e-6:
                axis = [1.0, 0.0, 0.0]
                d = 1e-6
            else:
                axis = [(positions[j][k] - positions[i][k]) / d for k in range(3)]
            push = target - d
            if i in frozen:
                share_i, share_j = 0.0, 1.0
            elif j in frozen:
                share_i, share_j = 1.0, 0.0
            else:
                share_i = share_j = 0.5
            for k in range(3):
                positions[i][k] -= axis[k] * push * share_i
                positions[j][k] += axis[k] * push * share_j
    if not relieved:
        return list(coords), 0
    return [(p[0], p[1], p[2]) for p in positions], relieved


def _relax_with_rdkit(
    perceived: PerceivedMolecule,
    coords: Sequence[Sequence[float]],
    fixed: Sequence[int],
    max_steps: int,
    method: str = 'uff',
) -> Optional[Tuple[List[Tuple[float, float, float]], str]]:
    """Constrained RDKit UFF minimisation with ``fixed`` atoms pinned.

    ``perceived`` supplies the topology only; ``coords`` are the coordinates
    to relax, so a perception made before the drag can be reused afterwards.

    The force field is built, used and dropped inside this function on
    purpose: an RDKit ``ForceField`` holds raw pointers into the conformer of
    its owning molecule, so a cached force field whose molecule has been
    garbage-collected segfaults the kernel.  A cold rebuild costs only
    0.1-10 ms for 18-392 atoms, so there is nothing to gain from caching one.
    """
    mol = perceived.typing_mol if perceived.typing_mol is not None else perceived.mol
    if mol is None or mol.GetNumAtoms() != len(coords):
        return None
    # The typing molecule has no conformer of its own; give it the input
    # geometry (same atom order, same count).
    work = Chem.Mol(mol)
    conformer = Chem.Conformer(work.GetNumAtoms())
    for idx, (x, y, z) in enumerate(coords):
        conformer.SetAtomPosition(idx, Point3D(float(x), float(y), float(z)))
    work.RemoveAllConformers()
    work.AddConformer(conformer, assignId=True)

    force_field = None
    try:
        with _RDKitQuiet():
            # A strained geometry can leave an atom over-bonded (a dragged
            # carbon that came within bonding distance of an amine nitrogen,
            # say).  Refreshing the property cache non-strictly accepts that
            # valence instead of letting the force-field constructor raise.
            work.UpdatePropertyCache(strict=False)
            engine = 'UFF'
            force_field = None
            if normalise_method(method) == 'mmff94':
                # MMFF94 covers organics well and is the better clean-up there,
                # but it has no transition-metal parameters at all, so fall
                # back rather than silently produce nothing.
                if _uff_helpers.MMFFHasAllMoleculeParams(work):
                    props = _uff_helpers.MMFFGetMoleculeProperties(work)
                    if props is not None:
                        force_field = _uff_helpers.MMFFGetMoleculeForceField(
                            work, props,
                        )
                        if force_field is not None:
                            engine = 'MMFF94'
            if force_field is None:
                if not _uff_helpers.UFFHasAllMoleculeParams(work):
                    return None
                force_field = _uff_helpers.UFFGetMoleculeForceField(work)
                engine = 'UFF'
            if force_field is None:
                return None
            for index in fixed:
                force_field.AddFixedPoint(int(index))
            force_field.Initialize()
            converged = force_field.Minimize(maxIts=int(max_steps))
            energy = force_field.CalcEnergy()
        positions = work.GetConformer().GetPositions()
        result = [(float(p[0]), float(p[1]), float(p[2])) for p in positions]
    finally:
        # Drop the force field before the molecule goes out of scope.
        force_field = None

    status = (
        f'RDKit {engine} relaxed {work.GetNumAtoms()} atoms '
        f'({len(fixed)} fixed), E = {energy:.2f} kcal/mol'
        + ('' if converged == 0 else f'; not converged in {max_steps} steps')
    )
    return result, status


def _relax_with_openbabel(
    perceived: PerceivedMolecule,
    coords: Sequence[Sequence[float]],
    fixed: Sequence[int],
    max_steps: int,
) -> Optional[Tuple[List[Tuple[float, float, float]], str]]:
    """Constrained Open Babel UFF minimisation via the repo's hardened wrapper.

    ``_optimize_xyz_openbabel_safe`` is used rather than Open Babel's
    ``FindForceField`` directly: the force-field object is process-global and
    keeps the constraint set of the previous call (there is no public "clear
    constraints" entry point), which the wrapper documents and handles.  Open
    Babel is also the only one of the two engines whose UFF actually has
    transition-metal parameters.
    """
    try:
        from delfin.smiles_converter import _optimize_xyz_openbabel_safe
    except Exception as exc:
        logger.debug('Open Babel relax unavailable: %s', exc)
        return None

    # Open Babel perceives its own bonds from the coordinates and regularly
    # misses part of a coordination sphere (on Zn(NH3)4 it finds two of the
    # four Zn-N bonds), which leaves ligands floating free of the metal.  The
    # metal-donor distances of the reference geometry are therefore pinned
    # explicitly -- the same restraint idea the exported payload uses, and
    # what the repo's own metal path does through
    # ``_build_coordination_constraints_from_xyz``.
    metal_pins: List[Tuple[int, int, float]] = []
    if perceived.metal_indices:
        metals = set(perceived.metal_indices)
        reference = perceived.coords
        for i, j in perceived.bonds:
            if (i in metals) == (j in metals):
                continue
            metal_pins.append((i, j, _distance(reference[i], reference[j])))

    block = format_xyz(perceived.symbols, coords, header=False)
    constraints = {
        'fix_atoms': [int(i) for i in fixed],
        'distances': metal_pins,
        'angles': [],
        'torsions': [],
    }
    try:
        relaxed = _optimize_xyz_openbabel_safe(
            block, steps=int(max_steps), coord_constraints=constraints
        )
    except Exception as exc:
        logger.debug('Open Babel relaxation failed: %s', exc)
        return None

    parsed = parse_xyz(relaxed)
    if parsed is None or len(parsed[0]) != perceived.n_atoms:
        return None
    relaxed_coords = parsed[1]
    unchanged = all(
        _distance(a, b) < 1e-9 for a, b in zip(relaxed_coords, coords)
    )
    status = (
        f'Open Babel UFF relaxed {perceived.n_atoms} atoms ({len(fixed)} fixed)'
        + (' - geometry unchanged (the wrapper rejected the result or it was '
           'already minimal)' if unchanged else '')
    )
    return relaxed_coords, status


def relax_xyz(
    xyz_text: str,
    fixed_indices: Optional[Iterable[int]] = None,
    max_steps: int = 200,
    *,
    perceived: Optional[PerceivedMolecule] = None,
    method: Optional[str] = 'uff',
) -> Dict[str, Any]:
    """Run the mouse-release clean-up minimisation.

    The browser's per-frame force field is deliberately cheap (a handful of
    steepest-descent steps per frame, no inversion terms).  When the drag
    ends, this restores a properly minimised geometry while holding the
    atoms the user positioned.

    Engine choice mirrors :func:`export_forcefield_terms`: Open Babel UFF for
    anything containing a metal (RDKit's UFF ignores transition metals
    entirely and would happily leave the metal where the drag put it), RDKit
    UFF with ``AddFixedPoint`` otherwise.  The force field is built and
    discarded inside this call; it is never stored.

    Args:
        xyz_text: Standard XYZ text or a bare DELFIN coordinate block.
        fixed_indices: 0-based indices of atoms to hold in place, in input
            XYZ order.  Out-of-range and duplicate entries are dropped.
        max_steps: Maximum minimisation steps.
        perceived: Optional pre-computed :func:`perceive_molecule` result,
            used for its **topology only** -- the coordinates always come
            from ``xyz_text``.  Passing the perception made when manipulate
            mode was switched on is the recommended call pattern: it saves
            the perception cost and, more importantly, guarantees that the
            drag cannot change the perceived bonding underneath the user.

    Returns:
        ``{'ok': bool, 'xyz': str, 'source': str, 'status': str,
        'warnings': [str, ...]}``.  ``xyz`` is the relaxed structure in the
        same format the input used (header or bare block) and always has the
        same atoms in the same order; on failure it is the unchanged input
        and ``ok`` is False.
    """
    parsed = parse_xyz(xyz_text)
    if parsed is None:
        return {
            'ok': False,
            'xyz': xyz_text,
            'source': SOURCE_GEOMETRIC,
            'status': 'Coordinates could not be parsed; nothing was relaxed.',
            'warnings': ['Coordinates could not be parsed.'],
        }
    symbols, coords, had_header = parsed

    if perceived is not None and (
        perceived.n_atoms != len(symbols) or list(perceived.symbols) != list(symbols)
    ):
        logger.debug('Supplied perception does not match the coordinates; re-perceiving')
        perceived = None
    if perceived is None:
        perceived = perceive_molecule(xyz_text)
    if perceived is None or perceived.n_atoms != len(symbols):
        return {
            'ok': False,
            'xyz': xyz_text,
            'source': SOURCE_GEOMETRIC,
            'status': 'Perception failed; the geometry was left untouched.',
            'warnings': ['The molecule could not be perceived.'],
        }

    fixed: List[int] = []
    for raw in (fixed_indices or []):
        try:
            index = int(raw)
        except (TypeError, ValueError):
            continue
        if 0 <= index < perceived.n_atoms and index not in fixed:
            fixed.append(index)

    warnings = list(perceived.warnings)
    chosen = normalise_method(method)
    if len(fixed) == perceived.n_atoms and perceived.n_atoms:
        return {
            'ok': True,
            'xyz': xyz_text,
            'source': SOURCE_GEOMETRIC,
            'status': 'Every atom is fixed; nothing to relax.',
            'warnings': warnings,
        }

    coords, relieved = _relieve_hard_clashes(coords, fixed)
    if relieved:
        warnings.append(
            f'{relieved} atom pair(s) were closer than 0.7 A and had to be '
            'pushed apart before the minimisation could run.'
        )

    attempts: List[Tuple[str, Any]] = []
    if perceived.has_metal:
        attempts.append((SOURCE_OPENBABEL, _relax_with_openbabel))
        attempts.append((SOURCE_RDKIT, _relax_with_rdkit))
        warnings.append(
            'Metal present: relaxed with Open Babel UFF, which has '
            'transition-metal parameters; RDKit UFF does not.'
        )
        if chosen == 'mmff94':
            warnings.append(
                'MMFF94 has no transition-metal parameters at all, so this '
                'structure was relaxed with UFF regardless of the choice.'
            )
    else:
        attempts.append((SOURCE_RDKIT, _relax_with_rdkit))
        attempts.append((SOURCE_OPENBABEL, _relax_with_openbabel))

    for source, runner in attempts:
        try:
            if runner is _relax_with_rdkit:
                outcome = runner(perceived, coords, fixed, max_steps, chosen)
            else:
                outcome = runner(perceived, coords, fixed, max_steps)
        except Exception as exc:  # defensive: never break the drag UI
            logger.debug('Relaxation with %s failed: %s', source, exc)
            outcome = None
        if outcome is None:
            continue
        relaxed_coords, status = outcome
        if len(relaxed_coords) != perceived.n_atoms:
            continue
        relaxed_coords = list(relaxed_coords)
        # Fixed atoms must come back exactly where the user put them; both
        # engines honour the constraint, this only removes numerical drift.
        for index in fixed:
            relaxed_coords[index] = coords[index]
        return {
            'ok': True,
            'xyz': format_xyz(
                symbols, relaxed_coords, header=had_header, comment='relaxed',
            ),
            'source': source,
            'status': status,
            'warnings': warnings,
        }

    warnings.append('No force field could be set up for this structure.')
    return {
        'ok': False,
        'xyz': xyz_text,
        'source': SOURCE_GEOMETRIC,
        'status': 'No usable force field; the geometry was left untouched.',
        'warnings': warnings,
    }
