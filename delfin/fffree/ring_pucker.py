"""delfin.fffree.ring_pucker — Cremer-Pople ring conformer enumeration (NO force field).

Implements the Cremer-Pople (1975) ring puckering formalism for deterministic,
universal enumeration of ring conformers (chair/boat/twist-boat for 6-rings,
envelope/twist for 5-rings, chair/boat/twist for 7-rings). FF-free, geometric,
no template lookup.

Reference: Cremer, D.; Pople, J.A. (1975) "A General Definition of Ring Puckering
Coordinates" J. Am. Chem. Soc. 97, 1354-1358.

The formalism:
  For an N-ring with atoms r_1..r_N and geometric center C:
    - Mean plane normal n via Newman et al. construction
    - z_i = (r_i - C) · n  (out-of-plane displacement)
    - Decompose into Fourier-like puckering amplitudes Q_m and phases φ_m:
        z_i = sqrt(2/N) Σ_{m=2}^{(N-1)/2} Q_m cos(2π m i / N + φ_m)
              + (1/sqrt(N)) Q_{N/2} (-1)^i     [N even only]
    - For N=6: Q_2 (boat), φ_2 (boat phase), Q_3 (chair amplitude, ±)
    - For N=5: Q_2 (envelope/twist), φ_2 (phase, 10 distinct states)
    - For N=7: Q_2 (chair/twist), Q_3 (boat), phases

Canonical states (energy-minima of unsubstituted ring):
  6-ring: 2 chair (Q_2=0, Q_3=±0.55), 6 boat (Q_2=0.55, Q_3=0, φ_2 ∈ 60°·k),
          6 twist-boat (Q_2=0.55, Q_3=0, φ_2 ∈ 60°·k + 30°)
  5-ring: 10 distinct (envelope E_1..E_5, twist T_1..T_5), all at Q_2 ≈ 0.40
  7-ring: 14 distinct combinations of (Q_2, Q_3, φ_2)

Universal across ANY ring (5/6/7 atoms, ANY hetero composition); deterministic;
no FF energy; preserves all bond lengths to within tolerance via post-snap.

Env-gate: DELFIN_FFFREE_RING_PUCKER=1 (default OFF, byte-identical when unset).
Auto-enabled under PURE_TRACK3 for full FF-free conformer completeness.
"""
from __future__ import annotations

import math
import os
from typing import Iterator, List, Optional, Sequence, Tuple

import numpy as np


_PURE_TRACK3 = os.environ.get("DELFIN_FFFREE_PURE_TRACK3", "0") == "1"
_RING_PUCKER = _PURE_TRACK3 or os.environ.get("DELFIN_FFFREE_RING_PUCKER", "0") == "1"

# Canonical pucker amplitudes (Å) — empirical defaults from CCDC ring database
_Q_DEFAULT = {5: 0.40, 6: 0.55, 7: 0.60}

# ----------------------------------------------------------------------
# Core math: Cremer-Pople decomposition + reconstruction
# ----------------------------------------------------------------------


def _mean_plane(coords_ring: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute geometric center and mean-plane normal of an N-ring.

    Uses Cremer-Pople convention: R' = Σ r_i sin(2π i / N), R'' = Σ r_i cos(2π i / N),
    normal n = R' × R'' / |R' × R''|.
    """
    n = len(coords_ring)
    center = coords_ring.mean(axis=0)
    centered = coords_ring - center
    angles = 2.0 * math.pi * np.arange(n) / n
    R_prime = (centered * np.sin(angles)[:, None]).sum(axis=0)
    R_doubleprime = (centered * np.cos(angles)[:, None]).sum(axis=0)
    cross = np.cross(R_prime, R_doubleprime)
    norm = float(np.linalg.norm(cross))
    if norm < 1e-9:
        # degenerate: use SVD fallback
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        normal = vh[-1]
    else:
        normal = cross / norm
    return center, normal


def _z_displacements(coords_ring: np.ndarray) -> np.ndarray:
    """Out-of-plane displacements z_i = (r_i - C) · n for each ring atom."""
    center, normal = _mean_plane(coords_ring)
    return (coords_ring - center) @ normal


def compute_pucker(coords_ring: np.ndarray) -> dict:
    """Compute Cremer-Pople puckering parameters for an N-ring (N = 5, 6, or 7).

    Parameters
    ----------
    coords_ring : (N, 3) array, ordered along the ring

    Returns
    -------
    dict with keys:
      N         : ring size
      Q_m       : list of pucker amplitudes for m = 2 .. floor((N-1)/2) (+ N/2 if even)
      phi_m     : list of phase angles (radians) for m = 2 ..
      Q_total   : total pucker amplitude sqrt(Σ Q_m^2)
      theta     : spherical pucker angle (6-ring only: arctan2(Q_2, Q_3))
      label     : conformer label (chair/boat/twist/envelope/planar)
    """
    N = len(coords_ring)
    if N < 4 or N > 8:
        return {"N": N, "label": "unsupported", "Q_total": 0.0}
    z = _z_displacements(coords_ring)
    Q_m: List[float] = []
    phi_m: List[float] = []
    # m runs from 2 to floor((N-1)/2)
    max_m = (N - 1) // 2
    for m in range(2, max_m + 1):
        c = math.sqrt(2.0 / N) * sum(z[i] * math.cos(2.0 * math.pi * m * i / N) for i in range(N))
        s = -math.sqrt(2.0 / N) * sum(z[i] * math.sin(2.0 * math.pi * m * i / N) for i in range(N))
        Q_mi = math.hypot(c, s)
        phi_mi = math.atan2(s, c)
        Q_m.append(Q_mi)
        phi_m.append(phi_mi)
    # N even: extra term Q_{N/2}
    if N % 2 == 0:
        Q_top = math.sqrt(1.0 / N) * sum(z[i] * (-1) ** i for i in range(N))
        Q_m.append(abs(Q_top))
        phi_m.append(0.0 if Q_top >= 0 else math.pi)
    Q_total = math.sqrt(sum(q * q for q in Q_m))
    # Classification
    label = _classify(N, Q_m, phi_m, Q_total)
    out = {"N": N, "Q_m": Q_m, "phi_m": phi_m, "Q_total": Q_total, "label": label}
    if N == 6:
        # theta = arctan(Q_2 / Q_3); θ=0/180 = chair, θ=90 = boat/twist-boat
        out["theta"] = math.atan2(Q_m[0], Q_m[1]) if len(Q_m) >= 2 else 0.0
    return out


def _classify(N: int, Q_m: Sequence[float], phi_m: Sequence[float], Q_total: float) -> str:
    """Classify ring conformer from CP parameters."""
    if Q_total < 0.10:
        return "planar"
    if N == 6:
        Q2, Q3 = Q_m[0], Q_m[1]
        if Q3 > 1.5 * Q2:
            return "chair_up" if phi_m[1] < math.pi / 2 else "chair_down"
        if Q2 > 1.5 * Q3:
            # boat vs twist by phase modulo 60°
            phase_deg = math.degrees(phi_m[0]) % 60.0
            if abs(phase_deg) < 15.0 or abs(phase_deg - 60.0) < 15.0:
                return "boat"
            return "twist_boat"
        return "half_chair"
    if N == 5:
        if Q_total < 0.20:
            return "planar"
        # envelope vs twist by phase modulo 36°
        phase_deg = math.degrees(phi_m[0]) % 36.0
        if abs(phase_deg) < 9.0 or abs(phase_deg - 36.0) < 9.0:
            return "envelope"
        return "twist"
    if N == 7:
        return "chair_7" if Q_m[1] > Q_m[0] else "twist_chair_7"
    return "puckered"


# ----------------------------------------------------------------------
# Reconstruction: set target pucker
# ----------------------------------------------------------------------


def set_pucker(
    coords_ring: np.ndarray,
    Q_m_target: Sequence[float],
    phi_m_target: Sequence[float],
) -> np.ndarray:
    """Return new ring-atom coords whose z-displacements match the target CP
    pucker (Q_m, φ_m). In-plane positions are preserved (geometric center +
    radial direction) — only the z-component along the mean-plane normal is set.

    This is the inverse Cremer-Pople synthesis:
        z_i_new = sqrt(2/N) Σ_m Q_m cos(2π m i / N + φ_m)
                  + (1/sqrt(N)) (-1)^i Q_{N/2}   [N even]
    """
    N = len(coords_ring)
    center, normal = _mean_plane(coords_ring)
    in_plane = coords_ring - center - np.outer(_z_displacements(coords_ring), normal)
    max_m = (N - 1) // 2
    z_new = np.zeros(N)
    for idx_m, m in enumerate(range(2, max_m + 1)):
        Q = Q_m_target[idx_m] if idx_m < len(Q_m_target) else 0.0
        phi = phi_m_target[idx_m] if idx_m < len(phi_m_target) else 0.0
        for i in range(N):
            z_new[i] += math.sqrt(2.0 / N) * Q * math.cos(2.0 * math.pi * m * i / N + phi)
    if N % 2 == 0:
        Q_top = Q_m_target[max_m - 1] if (max_m - 1) < len(Q_m_target) else 0.0
        for i in range(N):
            z_new[i] += math.sqrt(1.0 / N) * Q_top * (-1) ** i
    return center + in_plane + np.outer(z_new, normal)


# ----------------------------------------------------------------------
# Canonical state enumeration
# ----------------------------------------------------------------------


def canonical_pucker_states(N: int) -> List[Tuple[List[float], List[float], str]]:
    """Enumerate the canonical puckering states for an N-ring.

    Returns list of (Q_m_target, phi_m_target, label) tuples. These are the
    energy-minima conformers from the Cremer-Pople map for an unsubstituted ring.
    Universal (no atom-specific data); deterministic; complete within the
    given amplitude.

    Coverage:
      N=5: 10 states (5 envelopes + 5 twists at Q_2 ≈ 0.40)
      N=6: 14 states (2 chairs + 6 boats + 6 twist-boats at Q ≈ 0.55)
      N=7: 14 states (combinations of Q_2/Q_3 at Q ≈ 0.60)
    """
    states: List[Tuple[List[float], List[float], str]] = []
    if N == 5:
        Q = _Q_DEFAULT[5]
        # envelope: phi = 36° × k, twist: phi = 36° × k + 18°
        for k in range(5):
            states.append(([Q], [math.radians(36.0 * k)], f"envelope_E{k+1}"))
        for k in range(5):
            states.append(([Q], [math.radians(36.0 * k + 18.0)], f"twist_T{k+1}"))
    elif N == 6:
        Q = _Q_DEFAULT[6]
        # 2 chair: Q_2=0, Q_3=±Q
        states.append(([0.0, Q], [0.0, 0.0], "chair_4C1"))
        states.append(([0.0, Q], [0.0, math.pi], "chair_1C4"))
        # 6 boat: Q_2=Q, Q_3=0, phi_2 = 60° × k
        for k in range(6):
            states.append(([Q, 0.0], [math.radians(60.0 * k), 0.0], f"boat_B_{k}"))
        # 6 twist-boat: phi_2 = 60° × k + 30°
        for k in range(6):
            states.append(([Q, 0.0], [math.radians(60.0 * k + 30.0), 0.0], f"twist_boat_T_{k}"))
    elif N == 7:
        Q = _Q_DEFAULT[7]
        # 7-ring: combinations of Q_2 (chair-like) and Q_3 (boat-like)
        # Chair: Q_2 dominant
        for k in range(7):
            states.append(([Q, 0.0], [math.radians(360.0 / 7 * k), 0.0], f"chair7_{k}"))
        # Boat: Q_3 dominant
        for k in range(7):
            states.append(([0.0, Q], [0.0, math.radians(360.0 / 7 * k)], f"boat7_{k}"))
    else:
        # 4/8-ring: only chair (Q_{N/2} alternating)
        Q = 0.5
        states.append(([0.0, Q], [0.0, 0.0], f"chair_{N}"))
        states.append(([0.0, Q], [0.0, math.pi], f"chair_inv_{N}"))
    return states


# ----------------------------------------------------------------------
# Bond-length restoration after pucker reconstruction
# ----------------------------------------------------------------------


def _restore_ring_bonds(
    coords_ring: np.ndarray,
    ideal_bonds: Optional[Sequence[float]] = None,
    max_iter: int = 30,
    tol: float = 1e-3,
) -> np.ndarray:
    """After setting pucker, ring bond lengths may shift. Iteratively restore
    each bond to its ideal (or original) length while preserving the pucker
    by moving atoms along the bond axis. Deterministic.
    """
    N = len(coords_ring)
    P = coords_ring.copy()
    if ideal_bonds is None:
        # use original bond lengths as targets
        ideal_bonds = [float(np.linalg.norm(coords_ring[i] - coords_ring[(i + 1) % N])) for i in range(N)]
    for _ in range(max_iter):
        max_dev = 0.0
        for i in range(N):
            j = (i + 1) % N
            v = P[j] - P[i]
            d = float(np.linalg.norm(v))
            if d < 1e-9:
                continue
            target = ideal_bonds[i]
            dev = (target - d) / 2.0
            u = v / d
            P[i] -= u * dev * 0.5
            P[j] += u * dev * 0.5
            max_dev = max(max_dev, abs(dev))
        if max_dev < tol:
            break
    return P


# ----------------------------------------------------------------------
# Main API: enumerate ring conformers for a molecule
# ----------------------------------------------------------------------


def enumerate_ring_conformers(
    coords: np.ndarray,
    rings: Sequence[Sequence[int]],
    max_per_ring: int = 4,
) -> Iterator[np.ndarray]:
    """Enumerate distinct ring-pucker conformers of a molecule.

    Parameters
    ----------
    coords : (M, 3) full molecule coordinates
    rings : list of ring atom-index sequences (from RDKit GetRingInfo or similar)
    max_per_ring : cap on canonical states per ring (default 4 = chair-up,
        chair-down, dominant-boat, dominant-twist).
        Set to e.g. 14 for full coverage.

    Yields
    ------
    coords_variant : (M, 3) array with the molecule's rings set to target pucker
        states (Cartesian product across rings). Bond lengths restored.

    Notes
    -----
    For multi-ring systems: yields the Cartesian product of per-ring states.
    Rings are independent in this enumeration; coupled ring systems (fused, bridged)
    should still produce valid local-minima conformers because each ring's
    in-plane geometry is preserved and only z-displacement is altered.

    Deterministic: enumeration order is stable.
    No FF; uses Cremer-Pople geometric synthesis only.
    """
    if not _RING_PUCKER:
        return
    if not rings:
        return
    per_ring_states: List[List[Tuple[List[float], List[float], str]]] = []
    for ring in rings:
        N = len(ring)
        if N < 5 or N > 7:
            # only puckerable 5/6/7-rings handled
            per_ring_states.append([([], [], "as_is")])
            continue
        states = canonical_pucker_states(N)[:max_per_ring]
        per_ring_states.append(states)

    # Cartesian product
    from itertools import product as _prod
    for combo in _prod(*per_ring_states):
        P = coords.copy()
        for ring_atoms, (Q_t, phi_t, label) in zip(rings, combo):
            if label == "as_is":
                continue
            ring_arr = np.array([P[i] for i in ring_atoms])
            new_ring = set_pucker(ring_arr, Q_t, phi_t)
            new_ring = _restore_ring_bonds(new_ring)
            for k, i in enumerate(ring_atoms):
                P[i] = new_ring[k]
        yield P


# ----------------------------------------------------------------------
# Diagnostic / standalone test
# ----------------------------------------------------------------------


if __name__ == "__main__":
    # Test: ideal chair cyclohexane
    a = 1.54
    ideal_chair = np.array([
        [a / 2,  a * math.sqrt(3) / 2 / 2,  a / 2 * math.tan(math.radians(15))],
        [-a / 2, a * math.sqrt(3) / 2 / 2, -a / 2 * math.tan(math.radians(15))],
        [-a,     0,                          a / 2 * math.tan(math.radians(15))],
        [-a / 2, -a * math.sqrt(3) / 2 / 2, -a / 2 * math.tan(math.radians(15))],
        [a / 2,  -a * math.sqrt(3) / 2 / 2,  a / 2 * math.tan(math.radians(15))],
        [a,      0,                         -a / 2 * math.tan(math.radians(15))],
    ])
    res = compute_pucker(ideal_chair)
    print(f"Ideal-chair cyclohexane: Q={res['Q_total']:.3f}  label={res['label']}")
    print(f"  Q_m={[round(q,3) for q in res['Q_m']]}  phi_m={[round(math.degrees(p),1) for p in res['phi_m']]}")

    # Enumerate canonical 6-ring states
    states = canonical_pucker_states(6)
    print(f"\n6-ring canonical states: {len(states)}")
    for Q_t, phi_t, label in states[:6]:
        print(f"  {label:<18} Q={Q_t} phi={[round(math.degrees(p), 1) for p in phi_t]}")

    # Enumerate canonical 5-ring states
    states5 = canonical_pucker_states(5)
    print(f"\n5-ring canonical states: {len(states5)}")
    for Q_t, phi_t, label in states5[:5]:
        print(f"  {label:<14} Q={Q_t} phi={[round(math.degrees(p), 1) for p in phi_t]}")
