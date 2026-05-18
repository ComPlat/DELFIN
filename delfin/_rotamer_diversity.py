"""Universal post-UFF rotamer-diversity generator (Welle-5l Track-6).

Goal
----
For every isomer XYZ produced by the SMILES → 3D pipeline, generate K
additional rotamer-frames so that downstream DFT/ML refinement can pick the
best relative orientation of bulky substituents (tBu, PMe3, iPr, NMe2, …).

Pattern addressed (D-SOWKAQ Test-File 5)
----------------------------------------
PMe3 ligands carry three rotatable C–P–C torsions. UFF after ETKDG often
parks two methyl groups in an eclipsed/near-clashing arrangement
(H-H @ 2.08 Å) because the local minimum it found is not the global
minimum. By sampling 3 staggered states per rotatable DOF and re-minimising,
we recover the globally best rotamer per isomer.

Universal features only
-----------------------
DOFs are detected from the **molecular graph** (atomic number, bond order,
ring membership, neighbour count). Never via SMILES strings, refcodes or
named-ligand patterns. The detector picks every single bond between two
heavy atoms that is

* not part of a ring,
* not a coordination (M–donor) bond,
* not a double/triple/aromatic bond,
* not terminal-degenerate (both endpoints have a non-H neighbour or are
  themselves an sp3 carbon with 3 H — i.e. methyl rotors),

and where rotation would actually move atoms (the rotating fragment
contains at least one heavy atom OR ≥ 3 hydrogens).

Wire-in
-------
``delfin.smiles_converter`` adds ``_apply_5l_t6_rotamer_diversity_if_enabled``
inside the per-isomer emission loop. Default OFF (env-flag gated).

Env-flags
---------
``DELFIN_5L_T6_ROTAMER_DIVERSITY`` (default ``0``) — master switch
``DELFIN_5L_T6_ROTAMER_K``         (default ``3``) — K extra frames / isomer
``DELFIN_5L_T6_ROTAMER_STATES``    (default ``3``) — staggered states / DOF
``DELFIN_5L_T6_ROTAMER_MAX_DOFS``  (default ``6``) — DOF cap per isomer
``DELFIN_5L_T6_ROTAMER_GRID_CAP``  (default ``64``) — combinatorial grid cap
``DELFIN_5L_T6_ROTAMER_MD_TOL``    (default ``0.05``) — Å, M–D invariant

Welle-5p-A hard-gate (default ON when ``DELFIN_5P_A_TOPOLOGY_HARDGATE=1``)
Every candidate frame is additionally passed through
:func:`delfin._topology_hash.topology_preserved` which rejects any
rotation that
    * changes the bond multiset,
    * moves an M–D edge,
    * flips an amine-H (or P-H / As-H) onto the metal side
      (∠(M-D-H) < 60° or H · · · M < 2.30 Å).
This catches the X10-ALEQEO Fe(CO)2(NH2-CH2-CH2-S)2 amine-H umbrella
inversion (User-Direktive 2026-05-18).
"""

from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.common.logging import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Environment configuration helpers
# ---------------------------------------------------------------------------


def _env_bool(name: str, default: bool = False) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    val = raw.strip().lower()
    return val in ("1", "true", "yes", "on")


def _env_int(name: str, default: int, lo: int = 1, hi: int = 1024) -> int:
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
    return _env_bool("DELFIN_5L_T6_ROTAMER_DIVERSITY", False)


def _topo_hardgate_enabled() -> bool:
    """Welle-5p-A master flag — read here so the wire-in cost is one
    `os.environ.get` per candidate when the gate is disabled.

    Default reverted 1 -> 0 on 2026-05-18 (Iter-17) — voll-pool b5defcd
    showed pool-wide sigma -2783 isomere over-rejection.  Per-class
    adaptive thresholds deferred to Iter-18.  Env-flag opt-in preserved.
    """
    return _env_bool("DELFIN_5P_A_TOPOLOGY_HARDGATE", False)


_METAL_ATOMIC_NUMBERS = frozenset(
    # alkali + alkaline earth
    [3, 11, 19, 37, 55, 4, 12, 20, 38, 56]
    # 3d transition (Sc–Zn)
    + list(range(21, 31))
    # 4d transition (Y–Cd)
    + list(range(39, 49))
    # lanthanides (La–Lu)
    + list(range(57, 72))
    # 5d transition (Hf–Hg)
    + list(range(72, 81))
    # main-group metals
    + [13, 31, 49, 50, 81, 82, 83]
    # actinides up to Pu
    + list(range(89, 95))
)


def _is_metal(atomic_num: int) -> bool:
    return int(atomic_num) in _METAL_ATOMIC_NUMBERS


# ---------------------------------------------------------------------------
# XYZ parsing / writing (OpenBabel-free path for default-off safety)
# ---------------------------------------------------------------------------


def _parse_delfin_xyz(xyz: str) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """Parse DELFIN-format XYZ (no atom-count header).

    Returns (symbols, coords). Raises ``ValueError`` on malformed lines.
    """
    symbols: List[str] = []
    coords: List[Tuple[float, float, float]] = []
    for line in xyz.splitlines():
        s = line.strip()
        if not s:
            continue
        parts = s.split()
        if len(parts) < 4:
            raise ValueError(f"malformed XYZ line: {line!r}")
        try:
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
        except ValueError as exc:
            raise ValueError(f"bad coord in line {line!r}") from exc
        symbols.append(parts[0])
        coords.append((x, y, z))
    if not symbols:
        raise ValueError("empty XYZ")
    return symbols, coords


def _format_delfin_xyz(
    symbols: Sequence[str], coords: Sequence[Tuple[float, float, float]]
) -> str:
    lines = []
    for sym, (x, y, z) in zip(symbols, coords):
        lines.append(f"{sym:4s} {x:12.6f} {y:12.6f} {z:12.6f}")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Graph construction via Open Babel (perceived connectivity)
# ---------------------------------------------------------------------------


def _build_ob_mol_from_xyz(xyz: str):
    """Return an OBMol for the given DELFIN XYZ string, or ``None`` on failure.

    We rely on Open Babel's bond perception. It is the same engine the
    surrounding pipeline uses for round-trip checks, which keeps DOF detection
    consistent with topology checks downstream.
    """
    try:
        from openbabel import pybel  # noqa: WPS433 – lazy import
    except Exception:
        return None
    try:
        symbols, coords = _parse_delfin_xyz(xyz)
    except Exception:
        return None
    n_atoms = len(symbols)
    if n_atoms == 0:
        return None
    std = f"{n_atoms}\n\n"
    for sym, (x, y, z) in zip(symbols, coords):
        std += f"{sym}  {x}  {y}  {z}\n"
    try:
        return pybel.readstring("xyz", std)
    except Exception:
        return None


def _graph_from_ob(ob_mol) -> Dict:
    """Return a plain Python graph dict from an Open Babel molecule.

    The graph captures everything we need: atomic numbers, neighbours per
    atom, single-bond marker, ring marker, and metal flag. We deliberately
    avoid keeping the OB handles alive past this function so the rotamer
    sampler operates on pure-Python objects.
    """
    try:
        from openbabel import pybel  # noqa: WPS433
    except Exception:
        return {}
    obmol = ob_mol.OBMol
    n = obmol.NumAtoms()
    atomic_nums: List[int] = []
    is_metal: List[bool] = []
    neighbours: List[List[int]] = [[] for _ in range(n)]
    for i in range(1, n + 1):
        atom = obmol.GetAtom(i)
        z = atom.GetAtomicNum()
        atomic_nums.append(z)
        is_metal.append(_is_metal(z))

    bonds: List[Tuple[int, int, int, bool, bool]] = []
    # (a, b, order, is_aromatic, is_ring)
    for bond in pybel.ob.OBMolBondIter(obmol):
        a = bond.GetBeginAtomIdx() - 1
        b = bond.GetEndAtomIdx() - 1
        order = bond.GetBondOrder()
        aromatic = bool(bond.IsAromatic())
        ring = bool(bond.IsInRing())
        bonds.append((a, b, order, aromatic, ring))
        neighbours[a].append(b)
        neighbours[b].append(a)

    return {
        "n_atoms": n,
        "atomic_nums": atomic_nums,
        "is_metal": is_metal,
        "neighbours": neighbours,
        "bonds": bonds,
    }


# ---------------------------------------------------------------------------
# DOF identification — universal graph features only
# ---------------------------------------------------------------------------


def _fragment_atoms_on_side(
    neighbours: Sequence[Sequence[int]],
    pivot: int,
    other: int,
) -> List[int]:
    """BFS from *pivot* without crossing the bond pivot–other.

    Returns the list of atom indices on the *pivot* side, including pivot.
    """
    visited = {pivot, other}
    out = [pivot]
    stack = [pivot]
    while stack:
        cur = stack.pop()
        for nbr in neighbours[cur]:
            if nbr in visited:
                continue
            visited.add(nbr)
            out.append(nbr)
            stack.append(nbr)
    return out


def identify_rotamer_dofs(graph: Dict, max_dofs: int = 6) -> List[Dict]:
    """Identify rotatable DOFs from a Python graph dict.

    A rotatable bond is a single, non-aromatic, non-ring bond between two
    heavy atoms where

    1. Neither endpoint is a metal (M-X stays fixed),
    2. The rotating fragment has at least one heavy atom OR ≥ 3 hydrogens
       (so a non-trivial group actually moves),
    3. The opposite endpoint has at least one heavy neighbour other than the
       pivot (otherwise the rotor sees an atom with only Hs and pivots
       cylindrical symmetry-equivalent — terminal C–H, terminal O–H, …).

    The "rotating fragment" is the side of the bond chosen to minimise the
    atom count so we get the cheap end. DOFs are ranked by **bulkiness**:
    heavy-atom count on the rotating side, descending — bulky rotors first.

    Returns a list of dicts with keys:
        ``"pivot"``, ``"anchor"`` — atom indices defining the rotation axis
        ``"rotating"``           — atom indices to rotate
        ``"score"``              — bulkiness score (heavy atoms moved)
        ``"is_methyl"``          — True for terminal CH3 rotors
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
            # never touch coordination bonds — rotamer pass is not coord redo
            continue
        if atomic_nums[a] == 1 or atomic_nums[b] == 1:
            continue
        key = (min(a, b), max(a, b))
        if key in seen:
            continue
        seen.add(key)

        # heavy-neighbour counts at each endpoint (excluding the partner)
        def _hcount(atom_idx: int, other_idx: int) -> Tuple[int, int]:
            heavy = 0
            hyd = 0
            for nbr in neighbours[atom_idx]:
                if nbr == other_idx:
                    continue
                if atomic_nums[nbr] == 1:
                    hyd += 1
                else:
                    heavy += 1
            return heavy, hyd

        a_heavy, a_h = _hcount(a, b)
        b_heavy, b_h = _hcount(b, a)

        # Avoid bonds with no movable atom on *both* sides (would be a lone-atom
        # to lone-atom case — impossible for heavy bond but stay safe)
        if a_heavy + a_h == 0 and b_heavy + b_h == 0:
            continue

        # Pick the *cheaper* side to rotate: fewer heavy atoms wins. This is
        # the chemically meaningful choice — a methyl on a phosphine should
        # rotate the methyl (3 Hs, 0 heavies), not the phosphine (3 heavies).
        # Tie-break by total atom count (heavy + H), then by lower pivot index
        # for determinism.
        a_total = a_heavy + a_h
        b_total = b_heavy + b_h
        side_a_cheaper = (
            a_heavy < b_heavy
            or (a_heavy == b_heavy and a_total < b_total)
            or (a_heavy == b_heavy and a_total == b_total and a < b)
        )
        if side_a_cheaper:
            rotating_root, anchor = a, b
            rot_heavy, rot_h = a_heavy, a_h
            anc_heavy = b_heavy
        else:
            rotating_root, anchor = b, a
            rot_heavy, rot_h = b_heavy, b_h
            anc_heavy = a_heavy

        # Anchor must have at least one OTHER heavy neighbour, else rotation
        # is degenerate (cylindrical symmetry: rotating a methyl on a terminal
        # C-H carbon still produces a unique geometry, but rotating around an
        # axis whose anchor has only H neighbours has no preferred reference —
        # exclude to keep the grid focussed on chemically relevant rotors).
        if anc_heavy < 1:
            continue

        # Movement filter: skip if the rotating fragment has neither a heavy
        # atom nor at least 3 hydrogens. Single -CH2- or single -OH bonds with
        # one H are skipped (covered by VSEPR-H realism pass elsewhere).
        if rot_heavy == 0 and rot_h < 3:
            continue

        rotating_atoms = _fragment_atoms_on_side(neighbours, rotating_root, anchor)
        # If the BFS wandered back to the metal (multi-metal complex with
        # bridging non-coord bond), reject — never rotate through a metal.
        if any(is_metal[idx] for idx in rotating_atoms):
            continue

        # bulkiness score: heavy atoms moved, then total atoms moved
        moved_heavy = sum(1 for idx in rotating_atoms if atomic_nums[idx] > 1)
        moved_total = len(rotating_atoms)

        dofs.append({
            "pivot": rotating_root,
            "anchor": anchor,
            "rotating": rotating_atoms,
            "score": moved_heavy * 100 + moved_total,
            "is_methyl": (rot_heavy == 0 and rot_h >= 3),
            "heavy_moved": moved_heavy,
            "total_moved": moved_total,
        })

    # Rank: bulky DOFs first (tBu / iPr / NMe2 / PMe3), then methyls
    dofs.sort(key=lambda d: (-d["score"], d["pivot"], d["anchor"]))
    if len(dofs) > max_dofs:
        dofs = dofs[:max_dofs]
    return dofs


# ---------------------------------------------------------------------------
# Cartesian rotation around an axis
# ---------------------------------------------------------------------------


def _rodrigues_rotate(
    coords: List[Tuple[float, float, float]],
    axis_origin: Tuple[float, float, float],
    axis_dir: Tuple[float, float, float],
    angle_rad: float,
    atom_indices: Sequence[int],
) -> List[Tuple[float, float, float]]:
    """Return a new coord list where *atom_indices* are rotated by angle.

    Rodrigues formula. ``axis_dir`` is a unit-vector along the rotation axis;
    if not unit, it is normalised here.
    """
    ax, ay, az = axis_dir
    norm = math.sqrt(ax * ax + ay * ay + az * az)
    if norm < 1e-9:
        return list(coords)
    ax /= norm
    ay /= norm
    az /= norm
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    one_minus_cos = 1.0 - cos_a
    ox, oy, oz = axis_origin
    out = list(coords)
    for idx in atom_indices:
        x, y, z = out[idx]
        dx = x - ox
        dy = y - oy
        dz = z - oz
        dot = ax * dx + ay * dy + az * dz
        # rotated = d*cos + (axis × d)*sin + axis*(axis·d)*(1-cos)
        cx = ay * dz - az * dy
        cy = az * dx - ax * dz
        cz = ax * dy - ay * dx
        rx = dx * cos_a + cx * sin_a + ax * dot * one_minus_cos
        ry = dy * cos_a + cy * sin_a + ay * dot * one_minus_cos
        rz = dz * cos_a + cz * sin_a + az * dot * one_minus_cos
        out[idx] = (ox + rx, oy + ry, oz + rz)
    return out


# ---------------------------------------------------------------------------
# Combinatorial grid (capped)
# ---------------------------------------------------------------------------


def _grid_iter(n_dofs: int, n_states: int, cap: int):
    """Yield up to *cap* combinations of state indices over n_dofs DOFs.

    Yields the identity combination (all 0) first so it can be skipped if
    it equals the baseline. The remaining combinations are produced in a
    deterministic radix-style order.
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
    # Cap exceeded: emit identity, then evenly spaced samples
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
# Energy + clash scoring via OB UFF
# ---------------------------------------------------------------------------


def _evaluate_xyz(xyz: str) -> Optional[float]:
    """Return a UFF energy estimate (lower is better) or ``None`` on failure.

    The scorer adds a heavy clash penalty so geometries with atom overlap are
    pushed out of the top-K. Bonded pairs are excluded — the UFF energy
    handles bond-length excursions internally.
    """
    try:
        from openbabel import pybel  # noqa: WPS433
    except Exception:
        return None
    try:
        symbols, coords = _parse_delfin_xyz(xyz)
    except Exception:
        return None
    n_atoms = len(symbols)
    if n_atoms == 0:
        return None
    std = f"{n_atoms}\n\n"
    for sym, (x, y, z) in zip(symbols, coords):
        std += f"{sym}  {x}  {y}  {z}\n"
    try:
        obmol = pybel.readstring("xyz", std)
    except Exception:
        return None
    try:
        ff = pybel._forcefields.get("uff")
        if ff is None or not ff.Setup(obmol.OBMol):
            return None
        energy = float(ff.Energy())
    except Exception:
        return None
    if not math.isfinite(energy):
        return None

    # Clash penalty: penalise non-bonded contacts that fall below VSEPR-safe
    # thresholds. Heavy-heavy < 1.5 Å and H-H < 1.9 Å are both penalised so
    # the ranker prefers rotamers that relieve PMe3-type H-H crowding.
    bonded: set = set()
    bonded_via_neighbour: set = set()
    try:
        # bonds (direct 1,2 pairs)
        for bond in pybel.ob.OBMolBondIter(obmol.OBMol):
            a = bond.GetBeginAtomIdx() - 1
            b = bond.GetEndAtomIdx() - 1
            bonded.add((min(a, b), max(a, b)))
        # 1,3 pairs: two atoms bonded to a common neighbour. These are normal
        # VSEPR-determined contacts (e.g. geminal H-H on the same C); ignoring
        # them prevents the penalty from punishing geometrically required
        # neighbours.
        atom_bonds: Dict[int, List[int]] = {i: [] for i in range(len(symbols))}
        for (i, j) in bonded:
            atom_bonds[i].append(j)
            atom_bonds[j].append(i)
        for mid in range(len(symbols)):
            nbrs = atom_bonds.get(mid, [])
            for ii in range(len(nbrs)):
                for jj in range(ii + 1, len(nbrs)):
                    a, b = nbrs[ii], nbrs[jj]
                    bonded_via_neighbour.add((min(a, b), max(a, b)))
    except Exception:
        pass
    penalty = 0.0
    n_at = len(symbols)
    for i in range(n_at):
        xi, yi, zi = coords[i]
        si = symbols[i]
        for j in range(i + 1, n_at):
            pair = (i, j)
            if pair in bonded or pair in bonded_via_neighbour:
                continue
            xj, yj, zj = coords[j]
            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            d2 = dx * dx + dy * dy + dz * dz
            sj = symbols[j]
            if si == "H" and sj == "H":
                # H-H non-bonded contact: penalise below 1.9 Å
                cut = 1.9
                weight = 200.0
            elif si == "H" or sj == "H":
                # X-H contact: tighter threshold 1.7 Å
                cut = 1.7
                weight = 200.0
            else:
                # heavy-heavy: 1.5 Å
                cut = 1.5
                weight = 1000.0
            if d2 < cut * cut:
                penalty += weight * (cut * cut - d2)
    return energy + penalty


# ---------------------------------------------------------------------------
# Topology / M-D invariant guard
# ---------------------------------------------------------------------------


def _coord_bond_invariant_holds(
    graph: Dict,
    base_coords: Sequence[Tuple[float, float, float]],
    new_coords: Sequence[Tuple[float, float, float]],
    tol: float = 0.05,
) -> bool:
    """Verify M–D bond distances are preserved within *tol* Å.

    We pull the metal indices and their first-shell heavy neighbours and
    re-measure distances on the rotated coords. Any deviation > tol fails.
    """
    n = graph.get("n_atoms", 0)
    if n == 0 or n != len(new_coords):
        return False
    atomic_nums = graph["atomic_nums"]
    neighbours = graph["neighbours"]
    is_metal = graph["is_metal"]
    for m_idx in range(n):
        if not is_metal[m_idx]:
            continue
        for d_idx in neighbours[m_idx]:
            if atomic_nums[d_idx] == 1:
                continue
            mx, my, mz = new_coords[m_idx]
            dx, dy, dz = new_coords[d_idx]
            d_new = math.sqrt(
                (mx - dx) ** 2 + (my - dy) ** 2 + (mz - dz) ** 2
            )
            bmx, bmy, bmz = base_coords[m_idx]
            bdx, bdy, bdz = base_coords[d_idx]
            d_base = math.sqrt(
                (bmx - bdx) ** 2 + (bmy - bdy) ** 2 + (bmz - bdz) ** 2
            )
            if abs(d_new - d_base) > tol:
                return False
    return True


def _topology_hash(graph: Dict) -> Tuple:
    """Return a tuple summarising the bond graph (order-independent)."""
    bonds = sorted(
        (min(a, b), max(a, b), order)
        for (a, b, order, _arom, _ring) in graph.get("bonds", [])
    )
    return tuple(bonds)


# ---------------------------------------------------------------------------
# Main entry — apply rotamer diversity to a single XYZ
# ---------------------------------------------------------------------------


def apply(
    xyz: str,
    n_per_isomer: int = 3,
    n_states: int = 3,
    max_dofs: int = 6,
    grid_cap: int = 64,
    md_tol: float = 0.05,
) -> List[str]:
    """Return [base_xyz, rotamer_xyz_1, …, rotamer_xyz_K] (K = n_per_isomer).

    On failure or when no DOFs are found, returns ``[xyz]`` unchanged.
    """
    if n_per_isomer < 1:
        return [xyz]
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
        # OB perceived different atom count than XYZ has — bail out
        return [xyz]

    dofs = identify_rotamer_dofs(graph, max_dofs=max_dofs)
    if not dofs:
        return [xyz]

    base_topo = _topology_hash(graph)
    base_energy = _evaluate_xyz(xyz)

    # Build candidate set (sorted by energy ascending)
    candidates: List[Tuple[float, str, Tuple[int, ...]]] = []
    if base_energy is not None:
        candidates.append((base_energy, xyz, tuple([0] * len(dofs))))

    step_rad = 2.0 * math.pi / float(max(2, n_states))

    for combo in _grid_iter(len(dofs), n_states, grid_cap):
        if all(s == 0 for s in combo):
            # identity == base, already in candidates
            continue
        # Apply each DOF rotation in sequence to the base coords
        coords = list(base_coords)
        for dof_idx, state in enumerate(combo):
            if state == 0:
                continue
            dof = dofs[dof_idx]
            anchor = dof["anchor"]
            pivot = dof["pivot"]
            ax, ay, az = coords[anchor]
            px, py, pz = coords[pivot]
            axis = (px - ax, py - ay, pz - az)
            angle = step_rad * state
            coords = _rodrigues_rotate(
                coords, (ax, ay, az), axis, angle, dof["rotating"]
            )

        # Hard-guard: M–D bond lengths unchanged
        if not _coord_bond_invariant_holds(
            graph, base_coords, coords, tol=md_tol
        ):
            continue

        cand_xyz = _format_delfin_xyz(symbols, coords)
        # topology preservation (cheap: rebuild graph + compare bond multiset)
        cand_mol = _build_ob_mol_from_xyz(cand_xyz)
        if cand_mol is None:
            continue
        cand_graph = _graph_from_ob(cand_mol)
        if _topology_hash(cand_graph) != base_topo:
            continue

        # Welle-5p-A universal topology hard-gate (env-flag gated, default OFF).
        # Catches amine-H umbrella inversion + M-D set / bond-multiset changes
        # that the cheap OB topology hash above can miss when OB perceives
        # a different bond network on the modified geometry.
        if _topo_hardgate_enabled():
            try:
                from delfin import _topology_hash as _th  # local import
                _gate_res = _th.topology_preserved(
                    symbols,
                    base_coords,
                    coords,
                    md_tol=md_tol,
                    amine_h_min_deg=_env_float(
                        "DELFIN_5P_A_AMINE_H_MIN_DEG", 60.0, lo=0.0, hi=180.0
                    ),
                    h_to_metal_min_dist=_env_float(
                        "DELFIN_5P_A_CLASH_HM_MIN", 2.30, lo=0.0, hi=10.0
                    ),
                    bond_tol=_env_float(
                        "DELFIN_5P_A_BOND_TOL", 1.30, lo=0.5, hi=3.0
                    ),
                )
                if not _gate_res.passed:
                    logger.debug(
                        "5p-A rotamer rejected: %s", _gate_res.violations[:3]
                    )
                    continue
            except Exception as _gate_exc:  # pragma: no cover — safety
                logger.debug("5p-A rotamer gate failed open: %s", _gate_exc)

        energy = _evaluate_xyz(cand_xyz)
        if energy is None:
            continue
        candidates.append((energy, cand_xyz, combo))

    # Deduplicate by combo & rank by energy
    seen_combos: set = set()
    ranked: List[Tuple[float, str, Tuple[int, ...]]] = []
    candidates.sort(key=lambda x: x[0])
    for e, c_xyz, combo in candidates:
        if combo in seen_combos:
            continue
        seen_combos.add(combo)
        ranked.append((e, c_xyz, combo))

    if not ranked:
        return [xyz]

    # Always emit base XYZ first to preserve byte-identical default semantics
    # when n_per_isomer is treated as "extra". The wire-in callsite decides how
    # to label/merge these.
    out = [xyz]
    emitted = 0
    for e, c_xyz, combo in ranked:
        if c_xyz == xyz:
            continue
        out.append(c_xyz)
        emitted += 1
        if emitted >= n_per_isomer:
            break
    return out


# ---------------------------------------------------------------------------
# Public wire-in helper used by smiles_converter
# ---------------------------------------------------------------------------


def apply_if_enabled(xyz: str) -> List[str]:
    """Wire-in entry: returns ``[xyz, rotamer_1, …]`` when the env flag is set.

    Default-OFF semantics: if the master flag is unset, returns ``[xyz]``
    unchanged. The wire-in caller iterates the returned list and appends
    each extra frame as a labelled isomer.
    """
    if not _is_enabled():
        return [xyz]
    k = _env_int("DELFIN_5L_T6_ROTAMER_K", 3, lo=1, hi=32)
    n_states = _env_int("DELFIN_5L_T6_ROTAMER_STATES", 3, lo=2, hi=12)
    max_dofs = _env_int("DELFIN_5L_T6_ROTAMER_MAX_DOFS", 6, lo=1, hi=32)
    grid_cap = _env_int("DELFIN_5L_T6_ROTAMER_GRID_CAP", 64, lo=2, hi=4096)
    md_tol = _env_float("DELFIN_5L_T6_ROTAMER_MD_TOL", 0.05, lo=0.0, hi=2.0)
    try:
        return apply(
            xyz,
            n_per_isomer=k,
            n_states=n_states,
            max_dofs=max_dofs,
            grid_cap=grid_cap,
            md_tol=md_tol,
        )
    except Exception as exc:  # pragma: no cover - safety net
        logger.debug("5l-T6 rotamer-diversity failed: %s", exc)
        return [xyz]
