"""GRIP -- Donor-Cone Degree of Freedom (Phase 4).

User direction (2026-06-05, follow-up to ebb8cc3 angle-to-metal layer):

    "winkel zu TM können ja auch optimiert werden oder? nicht
    bindungslängen und wo donor sitzt könnte das das problem beheben?"

Translation: optimise the angles to the metal but also "where the donor
sits".  The angle-to-metal layer (commit ebb8cc3) addresses the polar
M-D-X cone angle (how far X tilts off the M-D axis), but it does NOT
explicitly enumerate / optimise the AZIMUTHAL rotation of the donor's
substituent plane around the M-D bond axis.  That azimuth is one extra
degree of freedom per donor; rotating it does NOT change M-D length, the
donor position, or the M-D-X polar angle, but it DOES change inter-
ligand contacts and the relative orientation of substituents on adjacent
donors.

This module adds that DoF as a JOINT variable in the GRIP L-BFGS state.

Mechanism
---------
For each donor d that has >= 2 heavy non-metal, non-donor neighbours we
classify "cone X atoms" as the donor's heavy first-shell ligand atoms
(skipping the metal, other donors, and hapto-pi atoms).  We then expose a
single scalar `theta_d` (radians) that controls a rigid rotation around
the M-d axis applied to the entire cone X atoms (the rotation is anchored
at the donor position, axis = unit vector from M to D).

In the L-BFGS objective we augment the optimisation variable

    x = (P_flat, theta_1, ..., theta_K)

where K is the number of donors with a non-trivial cone.  At each loss
evaluation we reconstruct the effective coordinates

    P_eval[i] = R(theta_d, axis_d) @ (x_pos[i] - donor_d) + donor_d
                                  for i in cone_atoms(d)
    P_eval[i] = x_pos[i]                                      otherwise

and compute the loss + gradient on P_eval.  The chain rule yields the
gradient w.r.t. each theta_d as

    dL/dtheta_d = sum_{i in cone(d)} (dL/dP_eval[i]) . (axis_d x (P_eval[i] - donor_d))

and the gradient w.r.t. x_pos[i] for cone atoms is

    dL/dx_pos[i] = R(theta_d)^T @ (dL/dP_eval[i])

The whole construction is differentiable in theta_d, and the inter-
ligand clash floor terms see the rotated X positions so the cone
rotation is genuinely steered by the loss surface (not a heuristic
search).

Env-flag (default OFF, byte-identical to the legacy polish path)
----------------------------------------------------------------
    DELFIN_FFFREE_GRIP_DONOR_CONE=1            -- enable the cone DoF
    DELFIN_FFFREE_GRIP_DONOR_CONE_INCLUDE_H=1  -- include H atoms in the
                                                  cone subtree (default
                                                  heavy-only; H rides along
                                                  via existing C-H bond
                                                  pull, but the cone
                                                  rotation IS rigid so
                                                  including H costs
                                                  nothing geometrically
                                                  -- default OFF preserves
                                                  byte-identity with the
                                                  legacy heavy-only sphere
                                                  the angle-to-metal layer
                                                  uses).

Determinism
-----------
* Donor iteration order: sorted ascending.
* Cone-atom iteration order: sorted ascending.
* axis construction: M -> D direction, deterministic.
* No RNG.

Universality
------------
* Cone classification uses graph topology + (metal, donors, hapto_atoms)
  only.  No SMILES patterns, no class-specific weights, no force-field
  constants.  FF-free contract preserved.
"""
from __future__ import annotations

import math
import os
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

__all__ = [
    "DONOR_CONE_ENV",
    "DONOR_CONE_INCLUDE_H_ENV",
    "donor_cone_active",
    "donor_cone_include_h",
    "DonorCone",
    "build_donor_cones",
    "apply_donor_cone_rotations",
    "augmented_loss_and_grad",
]


# ---------------------------------------------------------------------------
# Env-flag plumbing
# ---------------------------------------------------------------------------
DONOR_CONE_ENV = "DELFIN_FFFREE_GRIP_DONOR_CONE"
DONOR_CONE_INCLUDE_H_ENV = "DELFIN_FFFREE_GRIP_DONOR_CONE_INCLUDE_H"


def donor_cone_active() -> bool:
    """``True`` iff the donor-cone DoF is enabled (default OFF)."""
    raw = os.environ.get(DONOR_CONE_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


def donor_cone_include_h() -> bool:
    """``True`` iff H atoms should be included in the cone subtree
    (default OFF)."""
    raw = os.environ.get(DONOR_CONE_INCLUDE_H_ENV, "").strip().lower()
    if not raw:
        return False
    return raw in ("1", "true", "yes", "on")


# ---------------------------------------------------------------------------
# Cone container
# ---------------------------------------------------------------------------
class DonorCone:
    """One cone DoF per donor.

    Attributes
    ----------
    donor : int
        Atom index of the donor.
    axis : ndarray (3,)
        Unit vector pointing from the metal to the donor.  Frozen for the
        duration of the polish (M and D are M-D-rigid).
    anchor : ndarray (3,)
        Coordinates of the donor (the rotation anchor).
    atoms : tuple of int
        Sorted indices of the cone X atoms whose positions are rotated by
        ``theta``.
    """

    __slots__ = ("donor", "axis", "anchor", "atoms")

    def __init__(self, donor: int, axis: np.ndarray, anchor: np.ndarray,
                 atoms: Sequence[int]):
        self.donor = int(donor)
        self.axis = np.asarray(axis, dtype=np.float64).reshape(3)
        self.anchor = np.asarray(anchor, dtype=np.float64).reshape(3)
        self.atoms = tuple(int(a) for a in atoms)


# ---------------------------------------------------------------------------
# Cone subtree construction
# ---------------------------------------------------------------------------
def _heavy_first_shell(mol, donor: int, exclude: Set[int]) -> List[int]:
    """Heavy (non-H) first-shell neighbours of ``donor`` not in ``exclude``."""
    out: List[int] = []
    try:
        atom = mol.GetAtomWithIdx(int(donor))
    except Exception:
        return out
    for nb in atom.GetNeighbors():
        idx = int(nb.GetIdx())
        if idx in exclude:
            continue
        try:
            sym = nb.GetSymbol()
        except Exception:
            sym = "*"
        if sym == "H":
            continue
        out.append(idx)
    return sorted(out)


def _include_h_neighbours(mol, donor: int, exclude: Set[int]) -> List[int]:
    """H first-shell neighbours of ``donor`` not in ``exclude``."""
    out: List[int] = []
    try:
        atom = mol.GetAtomWithIdx(int(donor))
    except Exception:
        return out
    for nb in atom.GetNeighbors():
        idx = int(nb.GetIdx())
        if idx in exclude:
            continue
        try:
            sym = nb.GetSymbol()
        except Exception:
            sym = "*"
        if sym != "H":
            continue
        out.append(idx)
    return sorted(out)


def _ligand_subtree(
    mol,
    donor: int,
    first_shell: Sequence[int],
    exclude: Set[int],
) -> List[int]:
    """Return the FULL ligand subtree rooted at each X in first_shell.

    Walk the bond graph BFS starting from each X in ``first_shell`` but
    NOT entering ``exclude`` (which contains the donor + metal + other
    donors + hapto-pi atoms).  This gives the entire downstream ligand
    fragment hanging off the donor on the X side.

    Including the entire subtree is necessary because the cone rotation
    is a RIGID rotation around the M-D axis: rotating only the first
    shell would stretch every bond from the first shell to its downstream
    atoms.  Rotating the whole subtree preserves all internal distances
    and angles by construction.

    A bridging ligand (X atom that is also bonded to a second donor)
    can cause the subtree to include the OTHER donor's downstream atoms,
    which would conflict with the other donor's cone DoF.  When the BFS
    encounters another donor (in ``exclude``) we stop the walk on that
    branch.

    Deterministic: sorted-output traversal.
    """
    # Start the BFS frontier with first-shell atoms.
    visited: Set[int] = set(int(x) for x in first_shell)
    # Atoms NOT to enter as walk targets: exclude (donor + metal + other
    # donors + hapto-pi) and the donor itself.
    barrier = set(exclude) | {int(donor)}
    frontier = sorted(visited)
    while frontier:
        cur = frontier.pop(0)
        try:
            atom = mol.GetAtomWithIdx(int(cur))
        except Exception:
            continue
        for nb in atom.GetNeighbors():
            idx = int(nb.GetIdx())
            if idx in barrier:
                continue
            if idx in visited:
                continue
            visited.add(idx)
            frontier.append(idx)
            frontier.sort()
    return sorted(visited)


def build_donor_cones(
    mol,
    metal: int,
    donors: Sequence[int],
    P: np.ndarray,
    hapto_atoms: Optional[Iterable[int]] = None,
    *,
    include_h: bool = False,
    min_x_count: int = 1,
) -> List[DonorCone]:
    """Construct one :class:`DonorCone` per donor that has a non-trivial cone.

    Parameters
    ----------
    mol : RDKit-Mol-like
        Molecular graph with bonds populated.
    metal : int
        Metal atom index.
    donors : sequence of int
        Donor atom indices.
    P : ndarray (N, 3)
        Initial coordinates -- used only to extract the M->D axis
        direction and the donor anchor position.  The axis is frozen for
        the duration of a polish (M-D rigidity guarantees the axis does
        not move).
    hapto_atoms : iterable of int, optional
        Atoms in the hapto-pi system; both the cone atoms and donors that
        are themselves hapto are excluded from cone construction.
    include_h : bool
        If True, H atoms in the donor's first shell join the cone subtree.
        Default False (heavy only).
    min_x_count : int
        Minimum number of cone atoms required to build a cone.  Default 1
        -- even a single X atom benefits from azimuth optimisation when
        another donor's cone overlaps it (inter-ligand clash).

    Returns
    -------
    list of DonorCone
        Deterministic: cones are returned sorted by donor index.
    """
    out: List[DonorCone] = []
    if metal is None or metal < 0:
        return out
    if mol is None:
        return out
    P = np.asarray(P, dtype=np.float64)
    if P.ndim == 1:
        P = P.reshape(-1, 3)
    try:
        n_atoms = int(mol.GetNumAtoms())
    except Exception:
        n_atoms = P.shape[0]
    if P.shape[0] < n_atoms:
        return out

    hapto_set: Set[int] = set(int(i) for i in (hapto_atoms or ()))
    donor_seq = sorted(int(d) for d in donors)
    # Frozen-or-donor set: never put metal or another donor into the cone
    # (other donors have their own cone; the metal is M-D-rigid).
    frozen_or_donor: Set[int] = {int(metal)} | set(donor_seq)
    exclude_for_cone: Set[int] = frozen_or_donor | hapto_set

    metal_pos = P[int(metal)]

    for d in donor_seq:
        if d in hapto_set:
            # Hapto-pi donors: piano-stool placement geometry, not a
            # discrete cone -- skip.
            continue
        try:
            donor_pos = P[int(d)]
        except Exception:
            continue
        delta = donor_pos - metal_pos
        n_axis = float(np.linalg.norm(delta))
        if n_axis < 1e-9:
            # Coincident M and D -- degenerate, skip.
            continue
        axis = delta / n_axis
        first_shell = _heavy_first_shell(mol, int(d), exclude_for_cone)
        if include_h:
            first_shell = sorted(
                set(first_shell)
                | set(_include_h_neighbours(mol, int(d), exclude_for_cone))
            )
        if len(first_shell) < int(min_x_count):
            continue
        # Expand to full ligand subtree (BFS from first_shell, barrier =
        # exclude + donor itself).  Without this the cone rotation
        # would rigid-rotate only the first shell and break every X-Y
        # bond from X into the donor's downstream fragment.
        subtree = _ligand_subtree(mol, int(d), first_shell, exclude_for_cone)
        # Drop atoms whose position is too close to the donor (1e-6 Å)
        # to make the cross product numerically meaningful.
        keep: List[int] = []
        for a in subtree:
            r = P[int(a)] - donor_pos
            if float(np.linalg.norm(r)) >= 1e-6:
                keep.append(int(a))
        if len(keep) < int(min_x_count):
            continue
        out.append(DonorCone(donor=int(d), axis=axis,
                             anchor=donor_pos.copy(), atoms=tuple(keep)))
    return out


# ---------------------------------------------------------------------------
# Rodrigues rotation -- vectorised
# ---------------------------------------------------------------------------
def _rodrigues_rotmat(axis: np.ndarray, theta: float) -> np.ndarray:
    """3x3 rotation matrix for ``theta`` (radians) around unit ``axis``.

    Uses Rodrigues' formula; ``axis`` is assumed unit-length.
    """
    c = math.cos(float(theta))
    s = math.sin(float(theta))
    one_minus_c = 1.0 - c
    x, y, z = float(axis[0]), float(axis[1]), float(axis[2])
    R = np.array([
        [c + x * x * one_minus_c, x * y * one_minus_c - z * s, x * z * one_minus_c + y * s],
        [y * x * one_minus_c + z * s, c + y * y * one_minus_c, y * z * one_minus_c - x * s],
        [z * x * one_minus_c - y * s, z * y * one_minus_c + x * s, c + z * z * one_minus_c],
    ], dtype=np.float64)
    return R


def apply_donor_cone_rotations(
    P_flat: np.ndarray,
    cones: Sequence[DonorCone],
    thetas: Sequence[float],
) -> np.ndarray:
    """Return ``P_eval`` obtained by rotating each cone's atoms.

    The rotation is anchored at ``cone.anchor`` (donor position at the
    start of the polish) and is around ``cone.axis``.  Atoms that are
    not in any cone get their ``P_flat`` coordinates unchanged.

    ``P_flat`` may be a (N, 3) array or a (3N,) flat array.  Returns a
    (N, 3) array always (caller reshapes as needed).
    """
    P_in = np.asarray(P_flat, dtype=np.float64)
    if P_in.ndim == 1:
        P_in = P_in.reshape(-1, 3)
    P_eval = P_in.copy()
    if not cones:
        return P_eval
    for cone, theta in zip(cones, thetas):
        if cone is None or not cone.atoms:
            continue
        th = float(theta)
        if th == 0.0:
            # Byte-identical short-circuit -- no rotation applied.
            continue
        R = _rodrigues_rotmat(cone.axis, th)
        for a in cone.atoms:
            r = P_in[a] - cone.anchor
            P_eval[a] = cone.anchor + R @ r
    return P_eval


# ---------------------------------------------------------------------------
# Augmented loss-and-gradient
# ---------------------------------------------------------------------------
def augmented_loss_and_grad(
    loss_and_grad_pos,
    cones: Sequence[DonorCone],
    n_atoms: int,
):
    """Wrap a positional ``loss_and_grad`` to optimise (positions, thetas).

    Parameters
    ----------
    loss_and_grad_pos : callable
        Takes a (3*N,) flat positional array and returns ``(L, G_flat)``
        where ``G_flat`` has shape ``(3*N,)``.
    cones : sequence of DonorCone
        Cones whose theta values are appended to the variable vector.
    n_atoms : int
        Number of atoms (used to split the augmented vector).

    Returns
    -------
    callable
        A function with signature ``f(x_aug) -> (L, grad_aug)``
        compatible with scipy.minimize(method="L-BFGS-B").

    Notes
    -----
    The convention is that ``x_aug = concat(P_flat, thetas)`` so the
    initial theta values default to 0 (no rotation).  The optimisation
    then explores both translational adjustments AND the azimuthal cone
    rotations jointly.
    """
    K = len(cones)
    n_pos = 3 * int(n_atoms)

    def f(x_aug: np.ndarray) -> Tuple[float, np.ndarray]:
        x_aug = np.asarray(x_aug, dtype=np.float64)
        P_flat = x_aug[:n_pos]
        thetas = x_aug[n_pos:n_pos + K]
        # Reshape and apply cone rotations to obtain effective positions.
        P_eval = apply_donor_cone_rotations(P_flat, cones, thetas)
        # Evaluate the position-based loss + gradient on the effective
        # coordinates (this is the bridge to the existing GRIP machinery).
        L, G_flat = loss_and_grad_pos(P_eval.reshape(-1))
        G_pos = np.asarray(G_flat, dtype=np.float64).reshape(n_atoms, 3)
        # Chain rule -- distribute G_pos back into x_aug.
        G_xpos = G_pos.copy()
        G_theta = np.zeros(K, dtype=np.float64)
        for k, cone in enumerate(cones):
            if cone is None or not cone.atoms:
                continue
            th = float(thetas[k])
            R = _rodrigues_rotmat(cone.axis, th)
            # Replace x_pos gradient for cone atoms by R^T @ G_pos so
            # the dx_pos gradient is in the un-rotated frame.
            RT = R.T
            for a in cone.atoms:
                G_xpos[a] = RT @ G_pos[a]
            # theta gradient = sum over cone atoms of dot(G_pos[a],
            #                  axis x (P_eval[a] - anchor))
            axis = cone.axis
            theta_grad = 0.0
            for a in cone.atoms:
                r_eval = P_eval[a] - cone.anchor
                tangent = np.cross(axis, r_eval)
                theta_grad += float(np.dot(G_pos[a], tangent))
            G_theta[k] = theta_grad
        grad_aug = np.empty(n_pos + K, dtype=np.float64)
        grad_aug[:n_pos] = G_xpos.reshape(-1)
        grad_aug[n_pos:n_pos + K] = G_theta
        return float(L), grad_aug

    return f
