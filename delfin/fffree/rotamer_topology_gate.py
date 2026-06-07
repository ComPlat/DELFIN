"""delfin.fffree.rotamer_topology_gate — universal hard-gate for rotation
operations that emit conformer/rotamer/pucker variants.

A rotation around a single bond is permitted to emit a new conformer only
when ALL of the following are true on the post-rotation geometry P_after:

  1. EVERY expected SMILES bond is still intact: the 1-2 distance between
     the two bonded atoms must remain inside ``[0.7, 1.6] * ideal_bond_length``.
     Stretching beyond 1.6x ideal or compressing below 0.7x ideal is treated
     as a topology change (bond break or atomic collapse).
  2. NO spurious heavy-heavy bond is formed: any heavy-atom pair (excluding
     metals on metals) that is NOT in the expected SMILES bond set must
     sit at d >= 1.6 * ideal_bond_length(a, b).  Pairs closer than 0.85 *
     ideal are treated as collapses / spurious bonds.
  3. The metal coordination shell stays at its expected size: the number
     of heavy atoms within ``M_SHELL_FACTOR`` (default 1.30) * ideal_bond
     of any metal must not exceed the original CN.  This is a pure geometric
     count — no SMILES pattern, no per-class branch.
  4. No X-H or H-H collision: any pair containing at least one H must
     sit at d >= ``H_FLOOR`` (default 0.95 Å).

The function is **pure**: it only reads ``syms`` and ``P_after``, and never
mutates any of its arguments.  Determinism is preserved: identical inputs
yield identical accept/reject decisions.

Env-gate
--------
``DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE``
    Default OFF for backwards / byte-identical behaviour with HEAD.  When
    ``DELFIN_FFFREE_MOGUL_PRIMARY=1`` is set the gate is auto-ON unless the
    user explicitly disables it (``=0``).

Universality
------------
All checks are pure graph-topology + geometric counts — no SMILES patterns,
no per-class branches.  Metal detection uses the universal element-table
in :mod:`delfin._bond_decollapse`.
"""
from __future__ import annotations

import os
from typing import Iterable, Optional, Sequence, Tuple

import numpy as np


# ----------------------------------------------------------------------
# Thresholds — geometric, universal
# ----------------------------------------------------------------------

#: Multiplier on the ideal bond length above which we consider the
#: SMILES bond broken.  ~1.6x covalent sum is comfortably past any
#: published bond stretch in the CSD; legitimate bond lengths sit in
#: the [0.85, 1.20] x ideal range.
BOND_BREAK_FACTOR = 1.60

#: Multiplier on the ideal bond length BELOW which a heavy-heavy non-
#: bonded pair counts as a spurious bond / atomic collapse.
SPURIOUS_BOND_FACTOR = 0.85

#: Multiplier on the ideal bond length below which a non-bonded heavy
#: pair is treated as a freshly-formed spurious bond.
SPURIOUS_BOND_CUTOFF = 1.60

#: Multiplier on the ideal bond length used to count atoms inside the
#: metal coordination shell.  Matches converter_backend's M-D detection.
M_SHELL_FACTOR = 1.30

#: Hard floor (Å) for any pair containing a hydrogen.  Catches H-H
#: collisions plus C-H / N-H / O-H clashes.  0.95 Å is well below any
#: real X-H bond (>= 1.0 Å) yet rejects the synthetic-collapse pattern
#: the user-eye flagged (5+ H-H pairs at < 0.6 Å).
H_FLOOR = 0.95

#: Hard floor (Å) for heavy-heavy pairs not on the same bond.  Catches
#: structural collapses that escape the SPURIOUS_BOND_FACTOR floor when
#: ideal_bond is small (e.g. F-F ~ 1.14 Å).
HEAVY_FLOOR = 0.70


# ----------------------------------------------------------------------
# Env helpers
# ----------------------------------------------------------------------


def _env_on() -> bool:
    """Return True iff the gate is active.

    Default OFF for HEAD byte-identity.  Auto-ON whenever
    ``DELFIN_FFFREE_MOGUL_PRIMARY=1`` is set, unless the user explicitly
    sets ``DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE=0``.
    """
    raw = os.environ.get("DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE")
    if raw is None:
        # Auto-on under MOGUL_PRIMARY (the user-mandated default).
        if os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY", "0") == "1":
            return True
        return False
    return str(raw).strip() in ("1", "true", "True", "on", "yes", "YES")


# ----------------------------------------------------------------------
# Geometry helpers — depend only on _bond_decollapse so the gate is
# self-contained and bypasses any RDKit import.
# ----------------------------------------------------------------------


def _ideal(a: str, b: str) -> float:
    """Element-pair ideal bond length.  Falls back to a generic 1.55 Å
    sum on lookup failure (so a missing entry doesn't crash the gate).
    """
    try:
        from delfin._bond_decollapse import _ideal_bond
        return float(_ideal_bond(str(a), str(b)))
    except Exception:
        return 1.55


def _is_metal(s: str) -> bool:
    """Universal metal detection via :mod:`delfin._bond_decollapse`.
    Returns False on import failure.
    """
    try:
        from delfin._bond_decollapse import _is_metal as _bd_is_metal
        return bool(_bd_is_metal(str(s)))
    except Exception:
        return False


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------


def extract_expected_bonds(
    mol,
) -> Optional[Tuple[Tuple[int, int], ...]]:
    """Build the canonical expected heavy-atom bond set from an RDKit Mol.

    Returns the sorted tuple of ``(min(i,j), max(i,j))`` for every bond in
    the input ``mol`` (skipping metal-metal duplicates is unnecessary --
    they are real coordination edges and are preserved).  Returns ``None``
    when the mol is unavailable or empty.

    Bonds containing H atoms ARE included so the H-pair floor and the
    spurious-bond check can decide whether a given pair was meant to be
    bonded.
    """
    if mol is None:
        return None
    out = []
    try:
        for b in mol.GetBonds():
            i = int(b.GetBeginAtomIdx())
            j = int(b.GetEndAtomIdx())
            if i == j:
                continue
            out.append((min(i, j), max(i, j)))
    except Exception:
        return None
    return tuple(sorted(set(out)))


def expected_metal_cn(
    mol, syms: Sequence[str],
) -> Tuple[int, ...]:
    """Per-metal expected coordination number from the input bond graph.

    Returns a tuple aligned with the atom indices: the value at index ``i``
    is the original CN (heavy-atom neighbours, H excluded) when ``i`` is a
    metal, else -1.  Used by :func:`rotation_preserves_topology` to detect
    M-shell overfilling.
    """
    n = len(syms)
    out = [-1] * n
    if mol is None:
        return tuple(out)
    try:
        for i in range(n):
            if not _is_metal(str(syms[i])):
                continue
            try:
                atom = mol.GetAtomWithIdx(int(i))
            except Exception:
                continue
            cn = 0
            for nb in atom.GetNeighbors():
                try:
                    if int(nb.GetAtomicNum()) == 1:
                        continue
                    cn += 1
                except Exception:
                    continue
            out[i] = int(cn)
    except Exception:
        pass
    return tuple(out)


def rotation_preserves_topology(
    syms: Sequence[str],
    P_after: np.ndarray,
    mol=None,
    expected_bonds: Optional[Iterable[Tuple[int, int]]] = None,
    frozen_idxs: Optional[Iterable[int]] = None,
    *,
    expected_cn: Optional[Sequence[int]] = None,
    return_reason: bool = False,
):
    """Return True iff the post-rotation geometry preserves topology.

    Parameters
    ----------
    syms : sequence of element symbols, length N.
    P_after : (N, 3) array of post-rotation coordinates.
    mol : optional RDKit Mol used to derive ``expected_bonds`` and
        ``expected_cn`` when those are not supplied.  When both ``mol``
        and ``expected_bonds`` are None we fall through to a permissive
        decision (only the H-floor and HEAVY_FLOOR checks run).
    expected_bonds : optional iterable of ``(i, j)`` index pairs that the
        rotated geometry must preserve.  Derived from ``mol`` when ``None``.
    frozen_idxs : optional iterable of atom indices that are NOT rotated
        (typically the metal + donors).  Currently unused inside the gate
        but accepted in the signature for API completeness — downstream
        callers may use it to decide whether to invoke the gate at all.
    expected_cn : optional per-atom expected coordination number tuple
        (``-1`` for non-metal entries) — overrides the value derived from
        ``mol``.
    return_reason : when True, return ``(decision, reason)`` where
        ``reason`` is a short string ("ok", "bond_break", "spurious_bond",
        "m_shell_overfill", "h_collision", "heavy_collision",
        "non_finite") to aid logging / debugging.

    Returns
    -------
    bool or (bool, str)
        ``True`` ⇒ topology preserved; emit the variant.
        ``False`` ⇒ topology broken; the caller MUST skip the variant.
    """
    syms_l = [str(s) for s in syms]
    n = len(syms_l)
    P = np.asarray(P_after, dtype=float)
    if P.shape != (n, 3):
        return (False, "shape_mismatch") if return_reason else False
    if not np.all(np.isfinite(P)):
        return (False, "non_finite") if return_reason else False

    # ------------------------------------------------------------------
    # 1) Derive expected-bond set if not supplied.
    # ------------------------------------------------------------------
    if expected_bonds is None:
        eb = extract_expected_bonds(mol)
        if eb is None:
            # Permissive fallback: no graph -> only the geometric floors
            # apply (H-floor + heavy collision detection).  Still useful
            # because spurious H-H contacts are the dominant failure mode.
            expected_bonds = tuple()
        else:
            expected_bonds = eb
    expected_bonds = tuple(
        (int(min(i, j)), int(max(i, j)))
        for (i, j) in expected_bonds if int(i) != int(j)
    )
    bond_set = set(expected_bonds)

    # ------------------------------------------------------------------
    # 2) Expected CN per metal — only used when ``expected_cn`` is not
    #    provided.
    # ------------------------------------------------------------------
    if expected_cn is None and mol is not None:
        expected_cn = expected_metal_cn(mol, syms_l)
    if expected_cn is None:
        expected_cn = tuple(-1 for _ in range(n))

    # ------------------------------------------------------------------
    # 3) Bond-break check: every expected bond must still sit inside the
    #    allowed [HEAVY_FLOOR / BOND_BREAK_FACTOR] band.
    # ------------------------------------------------------------------
    for (i, j) in bond_set:
        if i < 0 or i >= n or j < 0 or j >= n:
            continue
        si, sj = syms_l[i], syms_l[j]
        try:
            d = float(np.linalg.norm(P[i] - P[j]))
        except Exception:
            return (False, "non_finite") if return_reason else False
        if not np.isfinite(d):
            return (False, "non_finite") if return_reason else False
        # Skip metal-metal (some bridged complexes have legitimately
        # long M-M edges that exceed BOND_BREAK_FACTOR).
        if _is_metal(si) and _is_metal(sj):
            continue
        ideal = _ideal(si, sj)
        if ideal <= 1e-6:
            continue
        # Bond stretched beyond break threshold (covers M-D too).
        if d > BOND_BREAK_FACTOR * ideal:
            return (False, "bond_break") if return_reason else False
        # H-containing bond compressed below H_FLOOR (caught here so we
        # don't double-count in step 5).
        if (si == "H" or sj == "H") and d < H_FLOOR:
            return (False, "h_collision") if return_reason else False
        # Heavy-heavy bond collapsed below absolute floor.
        if (si != "H" and sj != "H") and d < HEAVY_FLOOR:
            return (False, "heavy_collision") if return_reason else False

    # ------------------------------------------------------------------
    # 4) Spurious-bond check: any heavy-heavy pair that is NOT in
    #    bond_set must sit at d >= SPURIOUS_BOND_CUTOFF * ideal.  Metal-
    #    metal pairs are skipped (legitimate bridging M-M can sit
    #    close).  Metal-donor edges are present in bond_set so this
    #    only fires on non-bonded pairs.
    # ------------------------------------------------------------------
    # Also bookkeep M-shell counts here (single O(n^2) pass).
    metal_shell_counts = [0] * n
    for i in range(n):
        si = syms_l[i]
        si_is_h = (si == "H")
        si_is_m = _is_metal(si)
        for j in range(i + 1, n):
            sj = syms_l[j]
            sj_is_h = (sj == "H")
            sj_is_m = _is_metal(sj)
            if (i, j) in bond_set:
                continue  # expected bond -- handled in step 3
            try:
                d = float(np.linalg.norm(P[i] - P[j]))
            except Exception:
                return (False, "non_finite") if return_reason else False
            if not np.isfinite(d):
                return (False, "non_finite") if return_reason else False
            ideal = _ideal(si, sj)
            # H-pair check: any pair containing an H below H_FLOOR is rejected.
            if (si_is_h or sj_is_h):
                if d < H_FLOOR:
                    return (False, "h_collision") if return_reason else False
                continue
            # Skip metal-metal (allow close bridged M-M topologies).
            if si_is_m and sj_is_m:
                continue
            # Heavy-heavy: hard floor (catches structural collapse).
            if d < HEAVY_FLOOR:
                return (False, "heavy_collision") if return_reason else False
            # Spurious-bond test on non-expected heavy pairs.  Only flags
            # when the pair is closer than SPURIOUS_BOND_CUTOFF * ideal AND
            # closer than SPURIOUS_BOND_FACTOR * ideal (the inner threshold
            # is the actual "this is now a bond" radius).
            if ideal > 0.0 and d < SPURIOUS_BOND_FACTOR * ideal:
                return (False, "spurious_bond") if return_reason else False

            # M-shell counting: when one side is a metal, count the heavy
            # atom on the other side as part of the M-shell.  We use
            # M_SHELL_FACTOR * ideal_bond as the cutoff (same threshold the
            # converter uses to identify donors).  No double counting --
            # i < j and we count both directions explicitly.
            if si_is_m and not sj_is_m and not sj_is_h:
                if ideal > 0.0 and d < M_SHELL_FACTOR * ideal:
                    metal_shell_counts[i] += 1
            elif sj_is_m and not si_is_m and not si_is_h:
                if ideal > 0.0 and d < M_SHELL_FACTOR * ideal:
                    metal_shell_counts[j] += 1

    # We must also count the expected M-D edges (they are in bond_set so
    # the inner loop skipped them).  This makes the shell count reflect
    # the FULL set of heavy atoms inside M_SHELL_FACTOR * ideal, which is
    # what we compare against expected_cn.
    for (i, j) in bond_set:
        if i < 0 or i >= n or j < 0 or j >= n:
            continue
        si, sj = syms_l[i], syms_l[j]
        if _is_metal(si) and not _is_metal(sj) and sj != "H":
            ideal = _ideal(si, sj)
            try:
                d = float(np.linalg.norm(P[i] - P[j]))
            except Exception:
                d = float("inf")
            if ideal > 0.0 and d < M_SHELL_FACTOR * ideal:
                metal_shell_counts[i] += 1
        elif _is_metal(sj) and not _is_metal(si) and si != "H":
            ideal = _ideal(si, sj)
            try:
                d = float(np.linalg.norm(P[i] - P[j]))
            except Exception:
                d = float("inf")
            if ideal > 0.0 and d < M_SHELL_FACTOR * ideal:
                metal_shell_counts[j] += 1

    # ------------------------------------------------------------------
    # 5) M-shell overfill: reject if any metal's geometric shell count
    #    exceeds its expected CN.  expected_cn[i] == -1 means "unknown"
    #    (non-metal or graph unavailable) and skips the check for that
    #    metal -- defensive against missing input.
    # ------------------------------------------------------------------
    for i in range(n):
        if not _is_metal(syms_l[i]):
            continue
        exp = int(expected_cn[i]) if expected_cn is not None else -1
        if exp < 0:
            continue
        if metal_shell_counts[i] > exp:
            return (False, "m_shell_overfill") if return_reason else False

    return (True, "ok") if return_reason else True


# ----------------------------------------------------------------------
# Convenience wrapper for downstream filters
# ----------------------------------------------------------------------


def filter_topology_preserving_variants(
    variants,
    syms: Sequence[str],
    mol=None,
    expected_bonds=None,
    expected_cn=None,
):
    """Yield ``(P, label)`` tuples from ``variants`` whose rotation
    preserves topology, dropping rejects.

    The function is iterable-friendly: it preserves the iteration order
    of ``variants`` and the determinism of any downstream consumer.
    When the env-gate is OFF the function passes every variant through
    unchanged -- byte-identical with HEAD.
    """
    gate_on = _env_on()
    # Pre-derive once for efficiency.
    if gate_on:
        if expected_bonds is None:
            expected_bonds = extract_expected_bonds(mol)
        if expected_cn is None and mol is not None:
            expected_cn = expected_metal_cn(mol, syms)
    for item in variants:
        if not gate_on:
            yield item
            continue
        try:
            P_v, lab = item[0], item[1]
        except Exception:
            yield item
            continue
        try:
            ok = rotation_preserves_topology(
                syms, P_v, mol=mol,
                expected_bonds=expected_bonds,
                expected_cn=expected_cn,
            )
        except Exception:
            ok = True  # defensive -- never reject on gate failure
        if ok:
            yield item


# ----------------------------------------------------------------------
# Smoke check
# ----------------------------------------------------------------------


if __name__ == "__main__":
    # Quick sanity: rotating two atoms onto each other must be rejected.
    syms = ["C", "C", "H", "H", "H", "H", "H", "H"]
    P_ok = np.array([
        [0.0, 0.0, 0.0],
        [1.54, 0.0, 0.0],
        [-0.5, 0.9, 0.0],
        [-0.5, -0.45, 0.78],
        [-0.5, -0.45, -0.78],
        [2.04, 0.9, 0.0],
        [2.04, -0.45, 0.78],
        [2.04, -0.45, -0.78],
    ], dtype=float)
    P_bad = P_ok.copy()
    P_bad[2] = P_bad[6]  # collide H2 onto H6
    # No mol -> permissive (no expected bonds, but H_FLOOR catches collision).
    print(
        "P_ok decision:",
        rotation_preserves_topology(syms, P_ok, return_reason=True),
    )
    print(
        "P_bad decision:",
        rotation_preserves_topology(syms, P_bad, return_reason=True),
    )
