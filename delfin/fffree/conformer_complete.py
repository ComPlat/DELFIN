"""conformer_complete.py — UNIVERSAL conformer-completeness pass (FF-free).

ROOT-CAUSE upgrade to the conformer/rotamer engine.  Operates on the FULL
emitted ensemble (a list of ``(xyz, label)`` frames sharing atom order per
structure) at the public ``smiles_to_xyz_isomers`` chokepoint, AFTER the legacy
+ native build and the coordination-integrity / topology gates.

The three problems this cures (user eye-validation, root cause)
---------------------------------------------------------------
1.  INCOMPLETE coverage — many rotatable single bonds in ligands are NEVER
    rotated.  The pre-existing :mod:`delfin._rotamer_diversity` picks the
    *cheaper* (fewer-heavy) side of each rotatable bond, so a tBu / large
    pendant substituent (whose heavy side is the bulky one) only ever spins its
    terminal methyls and the bulky group itself stays frozen (ABUSEY "too few",
    ATIDEM "tBu doesn't rotate").  Here we enumerate EVERY rotatable single bond
    that moves >= 1 heavy atom and rotate the HEAVY-bearing subtree on a
    deterministic discrete grid (hindered, not continuous).

2.  TOPOLOGY-BREAKING free rotation — a free 360 deg rotation can swing a group
    so it (a) tears an M-D coordinate bond, (b) puts an H onto the M-C axis /
    into the metal (AYOZIX), or (c) collides into another ligand forging a
    spurious non-bonded contact at bonding distance.  Every emitted rotamer is
    passed through a HARD TOPOLOGY-PRESERVING gate: the metal + its first-shell
    donors stay FROZEN (rotations act on peripheral subtrees only); a rotamer is
    rejected if any M-D distance leaves its base coordinated range, a ring bond
    breaks / over-stretches, or any previously-non-bonded heavy pair falls to
    bonding distance (spurious bond).

3.  DUPLICATE bloat — symmetry-equivalent / near-identical rotamers (ATOSAE,
    indistinguishable atoms) are not deduped by the combo-key logic upstream.
    Here the resulting ensemble is RMSD-deduped (heavy-atom Kabsch, threshold
    ~0.5 Angstrom, keep first per cluster) so duplicates collapse.  A
    deterministic per-structure cap bounds bloat.

Design constraints
-------------------
* UNIVERSAL / graph-only.  DOFs come from perceived connectivity (Open Babel),
  never from SMILES strings, refcodes or named-ligand patterns.
* FF-FREE.  Pure geometry: Rodrigues rotation on a fixed staggered grid +
  geometric topology gate + Kabsch RMSD.  No force field is invoked.
* DETERMINISTIC.  Fixed grid, fixed ordering, no RNG / time / hash dependence.
* env-gated, default-OFF, BYTE-IDENTICAL when off.  Master flag
  ``DELFIN_FFFREE_CONF_COMPLETE`` (default ``0``).  With it unset, the public
  helpers are identity (return the input list untouched).
* NEVER raises.  Any failure falls back to the input frames.

Env-flags (all read once per call)
----------------------------------
``DELFIN_FFFREE_CONF_COMPLETE``      (default ``0``) — master switch.
``DELFIN_FFFREE_CONF_STATES``        (default ``3``) — staggered rotor states
                                       per DOF (3 = 120 deg gauche+/anti/gauche-).
``DELFIN_FFFREE_CONF_MAX_DOFS``      (default ``5``) — DOF cap per frame; the
                                       most-backbone-moving rotors win the grid.
``DELFIN_FFFREE_CONF_GRID_CAP``      (default ``96``) — combinatorial grid cap.
``DELFIN_FFFREE_CONF_RMSD``          (default ``0.5``) — Angstrom heavy-RMSD
                                       dedup threshold.
``DELFIN_FFFREE_CONF_MAX_PER_STRUCT`` (default ``24``) — emitted-frame cap per
                                       input structure (incl. its base frame).
``DELFIN_FFFREE_CONF_MD_TOL``        (default ``0.45``) — Angstrom M-D distance
                                       slack vs the base frame (coordinated
                                       range; crystals pack tight but stay
                                       coordinated).
``DELFIN_FFFREE_CONF_COUPLE``        (default ``1``) — enable coupled-grid search
                                       (rotate a neighbour rotor too when the
                                       single rotation clashes).

Spec lives in this docstring (re-buildable from it).
"""

from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

# Reuse the proven graph / DOF / rotation machinery — graph-only, FF-free.
from delfin._rotamer_diversity import (
    _build_ob_mol_from_xyz,
    _graph_from_ob,
    _parse_delfin_xyz,
    _format_delfin_xyz,
    _fragment_atoms_on_side,
    _rodrigues_rotate,
)

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Covalent radii (Cordero 2008, subset) — for the geometric topology gate.
# ---------------------------------------------------------------------------

_COV: Dict[int, float] = {
    1: 0.31, 5: 0.84, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 14: 1.11,
    15: 1.07, 16: 1.05, 17: 1.02, 33: 1.21, 34: 1.20, 35: 1.20, 52: 1.38,
    53: 1.39,
}
_COV_DEFAULT = 1.25


def _cov(z: int) -> float:
    return _COV.get(int(z), _COV_DEFAULT)


# ---------------------------------------------------------------------------
# Environment configuration helpers (mirror _rotamer_diversity conventions)
# ---------------------------------------------------------------------------


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in ("1", "true", "yes", "on")


def _env_int(name: str, default: int, lo: int = 1, hi: int = 4096) -> int:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        out = int(raw)
    except (TypeError, ValueError):
        return default
    return max(lo, min(hi, out))


def _env_float(name: str, default: float, lo: float = 0.0, hi: float = 10.0) -> float:
    raw = os.environ.get(name)
    if raw is None:
        return default
    try:
        out = float(raw)
    except (TypeError, ValueError):
        return default
    return max(lo, min(hi, out))


def _is_enabled() -> bool:
    return _env_bool("DELFIN_FFFREE_CONF_COMPLETE", False)


# ---------------------------------------------------------------------------
# DOF identification — COMPLETE coverage (rotate the heavy-bearing subtree).
# ---------------------------------------------------------------------------


def identify_complete_dofs(graph: Dict, max_dofs: int = 5) -> List[Dict]:
    """Identify EVERY rotatable single bond that moves >= 1 heavy atom.

    Differs from :func:`delfin._rotamer_diversity.identify_rotamer_dofs` in two
    ways that are the whole point of the completeness pass:

    * We rotate the side that ACTUALLY MOVES A HEAVY ATOM.  The upstream picker
      rotates the *cheaper* side, which for a tBu / biaryl / pincer-arm bond is
      the terminal-methyl side, leaving the bulky group frozen.  Here, of the two
      BFS subtrees on a bond, we pick the side that (a) does not contain a metal
      and (b) moves the FEWER atoms but still moves >= 1 heavy atom; ties broken
      deterministically.  This guarantees a genuine backbone/substituent torsion.

    * A bond is INCLUDED iff rotating the chosen side moves >= 1 heavy atom
      (``moved_heavy >= 1``).  Pure terminal-methyl / terminal-OH / terminal-NH
      rotors (moved_heavy == 0) are SKIPPED — they only spin H and are removed by
      the RMSD-dedup anyway, so emitting them is pure bloat.

    Rejected always:
      * coordination bonds (either endpoint a metal),
      * aromatic / ring / multiple bonds (rings are pucker-preserving, never
        freely rotated open),
      * bonds whose movable subtree wraps back onto a metal (multi-metal bridge).

    Returns a list of DOF dicts: ``pivot``, ``anchor`` (axis = anchor->pivot),
    ``rotating`` (atom indices to rotate, EXCLUDING anchor/metal),
    ``moved_heavy``.  Ranked by heavy atoms moved descending (bulky / backbone
    rotors first), deterministic tie-break by (pivot, anchor).
    """
    n = int(graph.get("n_atoms", 0))
    if n == 0:
        return []
    atomic_nums = graph["atomic_nums"]
    neighbours = graph["neighbours"]
    is_metal = graph["is_metal"]
    bonds = graph["bonds"]

    dofs: List[Dict] = []
    seen: set = set()
    for a, b, order, aromatic, ring in bonds:
        if order != 1 or aromatic or ring:
            continue
        if is_metal[a] or is_metal[b]:
            continue
        if atomic_nums[a] == 1 or atomic_nums[b] == 1:
            continue
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)

        # BFS the two subtrees (each excludes the partner endpoint).
        side_a = _fragment_atoms_on_side(neighbours, a, b)  # includes a
        side_b = _fragment_atoms_on_side(neighbours, b, a)  # includes b

        # A subtree that wraps onto a metal (bridging non-coord bond) is never
        # rotatable — never spin through a metal.
        a_has_metal = any(is_metal[i] for i in side_a)
        b_has_metal = any(is_metal[i] for i in side_b)

        def _moved(side: Sequence[int], root: int) -> Tuple[int, int]:
            # heavy / total atoms moved when rotating `side` around the bond,
            # excluding the root pivot atom itself? No -- the pivot atom moves
            # too (it is off-axis relative to the anchor unless collinear), so
            # count it.  But the "backbone leverage" is heavy atoms beyond the
            # root: a methyl root moves 0 heavy beyond itself.
            heavy = sum(1 for i in side if atomic_nums[i] > 1)
            beyond = sum(1 for i in side if i != root and atomic_nums[i] > 1)
            return heavy, beyond

        a_heavy, a_beyond = _moved(side_a, a)
        b_heavy, b_beyond = _moved(side_b, b)

        # Candidate sides: not wrapping a metal, and moving >= 1 heavy atom
        # beyond the rotor root (a genuine torsion that displaces structure).
        # (beyond==0 means only the root C + its own H's move = methyl-like.)
        cand = []
        if not a_has_metal and a_beyond >= 1:
            cand.append(("a", side_a, a, b, a_heavy, a_beyond, len(side_a)))
        if not b_has_metal and b_beyond >= 1:
            cand.append(("b", side_b, b, a, b_heavy, b_beyond, len(side_b)))
        if not cand:
            continue

        # Pick the CHEAPER side that still moves heavy structure: fewer total
        # atoms wins (cheap to rotate), tie-break fewer-heavy, then lower root
        # index for determinism.  Both sides are valid torsions; rotating the
        # smaller subtree is geometrically equivalent (rigid-body) and faster.
        cand.sort(key=lambda c: (c[6], c[4], c[2]))
        _tag, side, root, anchor, heavy, beyond, total = cand[0]

        rotating = [i for i in side if i != anchor]  # root included, anchor not
        dofs.append({
            "pivot": root,
            "anchor": anchor,
            "rotating": rotating,
            "moved_heavy": heavy,
            "moved_beyond": beyond,
        })

    # Rank: rotors that displace the MOST heavy atoms first (bulky tBu / biaryl /
    # pincer arms occupy the capped grid), deterministic tie-break.
    dofs.sort(key=lambda d: (-d["moved_beyond"], -d["moved_heavy"],
                             d["pivot"], d["anchor"]))
    if len(dofs) > max_dofs:
        dofs = dofs[:max_dofs]
    return dofs


# ---------------------------------------------------------------------------
# Geometric topology-preservation gate (FF-free, per rotamer).
# ---------------------------------------------------------------------------


def _base_bond_set(
    graph: Dict,
    coords: Sequence[Tuple[float, float, float]],
    bond_tol: float = 1.30,
) -> set:
    """Set of (i,j) heavy/H pairs that are bonded in the BASE geometry.

    Bonded iff dist <= bond_tol * (cov_i + cov_j).  Computed once on the base
    frame; the gate then forbids ANY new heavy-heavy bond not in this set
    (spurious bond) and requires every base ring / backbone bond to survive.
    """
    n = graph.get("n_atoms", 0)
    nums = graph["atomic_nums"]
    out = set()
    for i in range(n):
        ri = _cov(nums[i])
        xi, yi, zi = coords[i]
        for j in range(i + 1, n):
            t = bond_tol * (ri + _cov(nums[j]))
            dx = xi - coords[j][0]
            dy = yi - coords[j][1]
            dz = zi - coords[j][2]
            if dx * dx + dy * dy + dz * dz <= t * t:
                out.add((i, j))
    return out


def _md_pairs(graph: Dict) -> List[Tuple[int, int]]:
    """Metal -> first-shell heavy donor pairs from the graph (M first)."""
    n = graph.get("n_atoms", 0)
    nums = graph["atomic_nums"]
    nbrs = graph["neighbours"]
    is_metal = graph["is_metal"]
    out = []
    for m in range(n):
        if not is_metal[m]:
            continue
        for d in nbrs[m]:
            if nums[d] != 1:
                out.append((m, d))
    return out


def _ring_bonds(graph: Dict) -> List[Tuple[int, int]]:
    """In-ring bonds from the graph (each as a sorted (i,j) pair)."""
    out = []
    for a, b, _order, _arom, ring in graph["bonds"]:
        if ring:
            out.append((min(a, b), max(a, b)))
    return out


def _dist(c, i, j) -> float:
    dx = c[i][0] - c[j][0]
    dy = c[i][1] - c[j][1]
    dz = c[i][2] - c[j][2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def topology_preserved(
    graph: Dict,
    base_coords: Sequence[Tuple[float, float, float]],
    new_coords: Sequence[Tuple[float, float, float]],
    base_bonds: set,
    md_pairs: Sequence[Tuple[int, int]],
    ring_bonds: Sequence[Tuple[int, int]],
    md_tol: float = 0.45,
    ring_tol: float = 0.40,
    spurious_tol: float = 1.10,
    h_to_metal_min: float = 1.70,
) -> bool:
    """Return True iff *new_coords* preserve the base topology.

    Hard checks (any failure -> reject the rotamer):

    1.  M-D INTACT: every metal-donor distance stays within ``md_tol`` Angstrom
        of its base value (the donor stays coordinated; crystals pack tight but
        never decoordinate, so this is generous but absolute).
    2.  RING INTACT: every in-ring bond stays within ``ring_tol`` Angstrom of its
        base length (pucker is allowed; ring-open / over-stretch is not).
    3.  NO SPURIOUS BOND: no heavy-heavy pair that was NON-bonded in the base
        falls below ``spurious_tol * (cov_i + cov_j)`` (a forged bond / a group
        crashing into another ligand).
    4.  NO H ON / IN THE METAL: a hydrogen that is NOT bonded to a coordinating
        donor atom (i.e. not a legitimate alpha / agostic H) must not invade the
        metal coordination sphere -- it must stay no closer than the metal's
        nearest-heavy-donor distance MEASURED ON THE BASE FRAME (minus a tiny
        epsilon), and never closer than an absolute floor.  This catches the
        AYOZIX artifact where an SiMe3 methyl H swings onto the M-C axis / into
        the metal.  alpha-H's (bonded directly to a donor carbon, e.g. the
        Y-CH(SiMe3)2 methine H) are exempt because they ARE legitimately near
        the metal in the base structure.  Geometry-only.

    Pure geometry, deterministic.
    """
    n = graph.get("n_atoms", 0)
    nums = graph["atomic_nums"]
    is_metal = graph["is_metal"]
    nbrs = graph.get("neighbours", [[] for _ in range(n)])

    # 1. M-D intact
    for (m, d) in md_pairs:
        if abs(_dist(new_coords, m, d) - _dist(base_coords, m, d)) > md_tol:
            return False

    # 2. ring intact
    for (i, j) in ring_bonds:
        if abs(_dist(new_coords, i, j) - _dist(base_coords, i, j)) > ring_tol:
            return False

    # 4. no peripheral H invading a metal coordination sphere (AYOZIX).
    metals = [m for m in range(n) if is_metal[m]]
    if metals:
        # donor atoms per metal + nearest-donor distance on the BASE frame (the
        # reference shell the rotation must not penetrate).
        donors_of = {m: [] for m in metals}
        for (m, d) in md_pairs:
            if m in donors_of:
                donors_of[m].append(d)
        base_donor_min = {}
        for m in metals:
            base_donor_min[m] = min(
                (_dist(base_coords, m, d) for d in donors_of[m]),
                default=h_to_metal_min,
            )
        # an H is "exempt" if it is bonded to ANY coordinating donor atom (alpha
        # / agostic H -- legitimately close to the metal already in the base).
        donor_atoms = {d for ds in donors_of.values() for d in ds}
        exempt_h = set()
        for d in donor_atoms:
            for h in nbrs[d]:
                if 0 <= h < n and nums[h] == 1:
                    exempt_h.add(h)
        for m in metals:
            shell = base_donor_min.get(m, h_to_metal_min)
            for h in range(n):
                if nums[h] != 1 or h in exempt_h:
                    continue
                dh = _dist(new_coords, m, h)
                if dh < h_to_metal_min or dh < shell - 1e-3:
                    return False

    # 3. no spurious heavy-heavy bond (group crashed into another ligand).
    heavy = [k for k in range(n) if nums[k] > 1 and not is_metal[k]]
    for a in range(len(heavy)):
        i = heavy[a]
        ri = _cov(nums[i])
        xi, yi, zi = new_coords[i]
        for b in range(a + 1, len(heavy)):
            j = heavy[b]
            if (i, j) in base_bonds:
                continue
            t = spurious_tol * (ri + _cov(nums[j]))
            dx = xi - new_coords[j][0]
            dy = yi - new_coords[j][1]
            dz = zi - new_coords[j][2]
            if dx * dx + dy * dy + dz * dz < t * t:
                return False
    return True


# ---------------------------------------------------------------------------
# Deterministic combinatorial grid (capped, identity first).
# ---------------------------------------------------------------------------


def _grid_iter(n_dofs: int, n_states: int, cap: int):
    """Yield up to *cap* state-index tuples over n_dofs DOFs, identity first.

    Deterministic radix order.  When the full product exceeds *cap*, the identity
    is emitted first, then evenly-spaced codes, so coverage is spread.
    """
    total = n_states ** n_dofs
    if total <= cap:
        for code in range(total):
            digits = []
            v = code
            for _ in range(n_dofs):
                digits.append(v % n_states)
                v //= n_states
            yield tuple(digits)
        return
    yield tuple([0] * n_dofs)
    step = max(1, total // (cap - 1))
    emitted = 1
    for code in range(step, total, step):
        digits = []
        v = code
        for _ in range(n_dofs):
            digits.append(v % n_states)
            v //= n_states
        yield tuple(digits)
        emitted += 1
        if emitted >= cap:
            return


# ---------------------------------------------------------------------------
# Heavy-atom Kabsch RMSD (FF-free dedup).
# ---------------------------------------------------------------------------


def _kabsch_rmsd_heavy(
    syms: Sequence[str],
    A: Sequence[Tuple[float, float, float]],
    B: Sequence[Tuple[float, float, float]],
) -> float:
    """Heavy-atom RMSD between two SAME-ORDER coordinate sets after a proper-
    rotation Kabsch superposition.  Returns a large number if degenerate.

    NumPy when available; a pure-Python 3x3 SVD-free fallback otherwise (never
    raises).  Deterministic.
    """
    heavy = [k for k in range(len(syms)) if syms[k] != "H"]
    if len(heavy) < 3:
        # Too few heavy atoms to align — compare raw centroid-shifted RMSD.
        heavy = list(range(len(syms)))
        if len(heavy) < 1:
            return 1e9
    try:
        import numpy as np
        P = np.array([A[k] for k in heavy], float)
        Q = np.array([B[k] for k in heavy], float)
        if P.shape != Q.shape or P.shape[0] < 1:
            return 1e9
        Pc = P - P.mean(0)
        Qc = Q - Q.mean(0)
        if P.shape[0] >= 3:
            H = Pc.T @ Qc
            U, _S, Vt = np.linalg.svd(H)
            d = np.sign(np.linalg.det(Vt.T @ U.T))
            D = np.diag([1.0, 1.0, d])
            R = Vt.T @ D @ U.T
            diff = (R @ Pc.T).T - Qc
        else:
            diff = Pc - Qc
        return float(math.sqrt((diff * diff).sum() / len(heavy)))
    except Exception:
        # Pure-Python centroid-aligned (no rotation) RMSD fallback.
        try:
            ax = sum(A[k][0] for k in heavy) / len(heavy)
            ay = sum(A[k][1] for k in heavy) / len(heavy)
            az = sum(A[k][2] for k in heavy) / len(heavy)
            bx = sum(B[k][0] for k in heavy) / len(heavy)
            by = sum(B[k][1] for k in heavy) / len(heavy)
            bz = sum(B[k][2] for k in heavy) / len(heavy)
            s = 0.0
            for k in heavy:
                dx = (A[k][0] - ax) - (B[k][0] - bx)
                dy = (A[k][1] - ay) - (B[k][1] - by)
                dz = (A[k][2] - az) - (B[k][2] - bz)
                s += dx * dx + dy * dy + dz * dz
            return math.sqrt(s / len(heavy))
        except Exception:
            return 1e9


# ---------------------------------------------------------------------------
# Core: complete the conformer manifold for ONE frame.
# ---------------------------------------------------------------------------


def complete_frame(
    xyz: str,
    n_states: int = 3,
    max_dofs: int = 5,
    grid_cap: int = 96,
    rmsd_dedup: float = 0.5,
    max_emit: int = 24,
    md_tol: float = 0.45,
    couple: bool = True,
) -> List[str]:
    """Return ``[base_xyz, conf_1, ...]`` for one input frame.

    Enumerates all heavy-moving rotatable axes on a discrete staggered grid,
    hard-gates every rotamer for topology preservation, RMSD-dedups the result,
    and caps the count.  ``base_xyz`` is always element 0 (byte-identical to the
    input string).  On any failure / no DOFs, returns ``[xyz]``.
    """
    ob_mol = _build_ob_mol_from_xyz(xyz)
    if ob_mol is None:
        return [xyz]
    graph = _graph_from_ob(ob_mol)
    if not graph:
        return [xyz]
    try:
        symbols, base_coords = _parse_delfin_xyz(xyz)
    except Exception:
        return [xyz]
    if len(symbols) != graph["n_atoms"]:
        return [xyz]

    dofs = identify_complete_dofs(graph, max_dofs=max_dofs)
    if not dofs:
        return [xyz]

    base_bonds = _base_bond_set(graph, base_coords)
    md_pairs = _md_pairs(graph)
    ring_bonds = _ring_bonds(graph)

    step_rad = 2.0 * math.pi / float(max(2, n_states))

    # ----- enumerate the grid, gate, collect accepted coordinate sets -----
    accepted: List[Tuple[Tuple[int, ...], List[Tuple[float, float, float]]]] = []
    # base is always kept first
    accepted.append((tuple([0] * len(dofs)), list(base_coords)))

    for combo in _grid_iter(len(dofs), n_states, grid_cap):
        if all(s == 0 for s in combo):
            continue
        coords = list(base_coords)
        for di, state in enumerate(combo):
            if state == 0:
                continue
            dof = dofs[di]
            anchor = dof["anchor"]
            pivot = dof["pivot"]
            ax, ay, az = coords[anchor]
            px, py, pz = coords[pivot]
            axis = (px - ax, py - ay, pz - az)
            coords = _rodrigues_rotate(
                coords, (ax, ay, az), axis, step_rad * state, dof["rotating"]
            )
        if topology_preserved(
            graph, base_coords, coords, base_bonds, md_pairs, ring_bonds,
            md_tol=md_tol,
        ):
            accepted.append((combo, coords))

    # ----- COUPLED rescue: for single-DOF rotations that all failed the gate
    # (a rotor that clashes unless a neighbour also moves), search the small
    # 2-DOF coupled grid over the top pair of rotors.  Deterministic. ----------
    if couple and len(dofs) >= 2:
        # which single-DOF axes produced NO accepted non-identity frame?
        moved_single = {di for (combo, _c) in accepted
                        for di in range(len(dofs)) if combo[di] != 0
                        and sum(1 for s in combo if s != 0) == 1}
        stuck = [di for di in range(len(dofs)) if di not in moved_single]
        # couple each stuck DOF with the next DOF in rank order
        for di in stuck[:max_dofs]:
            dj = (di + 1) % len(dofs)
            if dj == di:
                continue
            for si in range(1, n_states):
                for sj in range(0, n_states):
                    coords = list(base_coords)
                    for (dd, st) in ((di, si), (dj, sj)):
                        if st == 0:
                            continue
                        dof = dofs[dd]
                        ax, ay, az = coords[dof["anchor"]]
                        px, py, pz = coords[dof["pivot"]]
                        axis = (px - ax, py - ay, pz - az)
                        coords = _rodrigues_rotate(
                            coords, (ax, ay, az), axis,
                            step_rad * st, dof["rotating"]
                        )
                    if topology_preserved(
                        graph, base_coords, coords, base_bonds, md_pairs,
                        ring_bonds, md_tol=md_tol,
                    ):
                        combo = [0] * len(dofs)
                        combo[di] = si
                        combo[dj] = sj
                        accepted.append((tuple(combo), coords))

    # ----- RMSD-dedup (keep first per cluster, deterministic order) -----
    kept_coords: List[List[Tuple[float, float, float]]] = []
    out_xyz: List[str] = []
    # deterministic ordering: identity first, then by combo radix
    accepted.sort(key=lambda ac: (sum(1 for s in ac[0] if s != 0), ac[0]))
    for (_combo, coords) in accepted:
        dup = False
        for kc in kept_coords:
            if _kabsch_rmsd_heavy(symbols, kc, coords) < rmsd_dedup:
                dup = True
                break
        if dup:
            continue
        kept_coords.append(coords)
        if not out_xyz:
            # element 0 must be byte-identical to the input string
            out_xyz.append(xyz)
        else:
            out_xyz.append(_format_delfin_xyz(symbols, coords))
        if len(out_xyz) >= max_emit:
            logger.debug("conf-complete cap %d hit", max_emit)
            break

    if not out_xyz:
        return [xyz]
    # guarantee element 0 is the verbatim input
    if out_xyz[0] != xyz:
        out_xyz.insert(0, xyz)
    return out_xyz


# ---------------------------------------------------------------------------
# Public ensemble-level wire-in.
# ---------------------------------------------------------------------------


def apply_to_ensemble(isomers):
    """Expand each emitted ``(xyz, label)`` frame with completeness conformers.

    Identity (returns *isomers* untouched) unless
    ``DELFIN_FFFREE_CONF_COMPLETE=1`` -> output BYTE-IDENTICAL when off.

    Each input frame ``(xyz, label)`` becomes ``(base_xyz, label)`` followed by
    ``(conf_xyz, "<label>_conf-N")`` for each accepted, deduped rotamer.
    Deterministic; never raises (falls back to *isomers* on error).
    """
    if not isomers or not _is_enabled():
        return isomers
    try:
        n_states = _env_int("DELFIN_FFFREE_CONF_STATES", 3, lo=2, hi=12)
        max_dofs = _env_int("DELFIN_FFFREE_CONF_MAX_DOFS", 5, lo=1, hi=32)
        grid_cap = _env_int("DELFIN_FFFREE_CONF_GRID_CAP", 96, lo=2, hi=8192)
        rmsd_dedup = _env_float("DELFIN_FFFREE_CONF_RMSD", 0.5, lo=0.0, hi=5.0)
        max_per = _env_int("DELFIN_FFFREE_CONF_MAX_PER_STRUCT", 24, lo=1, hi=512)
        md_tol = _env_float("DELFIN_FFFREE_CONF_MD_TOL", 0.45, lo=0.0, hi=3.0)
        couple = _env_bool("DELFIN_FFFREE_CONF_COUPLE", True)
    except Exception:
        return isomers

    out = []
    try:
        for item in isomers:
            if isinstance(item, (tuple, list)) and len(item) >= 2:
                xyz, label = item[0], item[1]
                rest = tuple(item[2:])
            else:
                xyz, label, rest = item, "", ()
            frames = complete_frame(
                xyz, n_states=n_states, max_dofs=max_dofs, grid_cap=grid_cap,
                rmsd_dedup=rmsd_dedup, max_emit=max_per, md_tol=md_tol,
                couple=couple,
            )
            if not frames:
                out.append(item)
                continue
            # frame 0 keeps the original label verbatim
            out.append((frames[0], label) + rest if rest else (frames[0], label))
            for k, fxyz in enumerate(frames[1:], start=1):
                clbl = f"{label}_conf-{k}" if label else f"conf-{k}"
                out.append((fxyz, clbl) + rest if rest else (fxyz, clbl))
    except Exception as exc:  # pragma: no cover — safety net
        logger.debug("conf-complete ensemble pass failed: %s", exc)
        return isomers
    return out or isomers
