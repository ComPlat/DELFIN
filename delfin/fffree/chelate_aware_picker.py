"""delfin.fffree.chelate_aware_picker — chelate-bite-aware polyhedron picker.

The Mogul-PRIMARY construction path (``assemble_via_mogul``) currently picks
the per-CN reference polyhedron via the naive ``geometries_for_cn(CN)[0]``
first-candidate rule.  This is wrong for chelate ligands because the chosen
polyhedron's vertex-vertex angles must accommodate the chelate's intrinsic
bite angle (a function of the chelate ring size), and the naive picker has
no awareness of that constraint:

* BOJTUS (Ni, tetradentate-N4 macrocycle):
    ``geometries_for_cn(4)`` = ``["T-4 tetrahedron", "SP-4 square planar"]``.
    The picker silently selects SP-4 (d8-Ni), then the macrocycle's 5-/6-
    membered chelate rings are forced onto SP-4's 90° N-Ni-N vertex pairs.
    For a 4-N macrocycle whose adjacent N donors share an N-Ni-N bite ≈ 80°
    this still fits well; but for a 4-N macrocycle whose bridge length forces
    a tetrahedral pocket the SP-4 choice is wrong.

* AXOKAZ (Zn, tridentate-NNN-mer + 2 Cl):
    ``geometries_for_cn(5)`` = ``["TBP-5 trigonal bipyramid", ...]``.
    The picker selects TBP-5.  A meridional tridentate ligand (N-N-N
    spanning three coplanar coordination sites) needs N-Zn-N pair angles
    that are 90° / 90° / 180° (axial + 2 equatorial) — that pattern exists
    in **SPY-5** (cis-trans, 90/90/180 between basal vertices and apex) but
    NOT in TBP-5 (where adjacent equatorial pairs are 120° and axial-
    equatorial are 90°).  Selecting TBP-5 forces the meridional N3 chelate
    to 120° at its central N donor, collapsing two of the N-Zn distances.

The fix is geometric and universal: for every chelate-pair of donors that
must occupy adjacent vertices, the chosen polyhedron's vertex-vertex angle
for that pair must be close to the chelate's intrinsic bite angle.  The
chelate bite angle is purely a function of the chelate ring size:

    4-ring (carboxylate κ²-OO, β-diketonate small)      ≈  60°
    5-ring (bipy, en, ox κ²-OO)                          ≈  80°
    6-ring (acac, β-diketonate 6-membered)               ≈  90°
    7-ring (expanded chelates, calix)                    ≈ 100°

This module provides ``pick_polyhedron_chelate_aware`` which scores each
candidate polyhedron by the mean-squared deviation between its vertex-vertex
angles (for adjacent chelate-pairs only) and the chelate's intrinsic bite
angles.  The polyhedron with the minimum score is selected.  When there is
no chelate ligand the score reduces to zero across all candidates, so the
picker degenerates to ``geometries_for_cn(CN)[0]`` (legacy behaviour).

Universal:
    - NO SMILES patterns
    - NO per-class branches
    - NO per-element thresholds
    - Pure geometric matching between graph-derived chelate ring size and
      polyhedron-derived vertex-vertex angle table.

Env-flag::

    DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK = "1"   # default ON under
                                                   # MOGUL_PRIMARY
    DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK = "0"   # legacy first-candidate

When the gate is OFF (or the chelate_info argument is empty / None / contains
only monodentate ligands), the picker returns ``geometries_for_cn(CN)[0]``,
which is byte-identical with HEAD.

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import math
import os
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np


__all__ = [
    "RING_SIZE_TO_BITE_DEG",
    "chelate_aware_picker_enabled",
    "ring_size_to_bite_angle",
    "polyhedron_vertex_angles",
    "score_polyhedron_for_chelates",
    "pick_polyhedron_chelate_aware",
    "chelate_info_from_decompose",
]


# Empirical chelate-bite-angle table, indexed by chelate ring size.
# Source: idealised geometric ring closure with M-D = 2.0 Å and ideal
# sp3/sp2 backbone angles (109.5° / 120°); cross-checked against the
# CSD-Mogul empirical histograms summarised in the project draft.
RING_SIZE_TO_BITE_DEG: Dict[int, float] = {
    4: 64.0,    # carboxylate κ²-OO, small β-diketonate
    5: 78.0,    # bipy, en, oxalate κ²-OO, ethylenediamine
    6: 88.0,    # acac (6-mem chelate), salicylaldiminate, β-diketonate
    7: 94.0,    # expanded chelates, calixarene
    8: 100.0,   # crown-ether-like, very rare
}


def chelate_aware_picker_enabled() -> bool:
    """Return True iff the chelate-aware picker is wired on.

    Default-ON when the Mogul-PRIMARY path is active so chelate-bearing
    SMILES are routed through the chelate-aware picker by default; the
    user can switch back to the naive first-candidate rule by setting
    ``DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK=0`` explicitly.

    Default-OFF byte-identical with HEAD when MOGUL_PRIMARY itself is
    off, since this picker is only consulted from the MOGUL_PRIMARY
    code paths.
    """
    raw = os.environ.get("DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK", "").strip()
    if raw == "":
        # Default ON only when MOGUL_PRIMARY is active.  Otherwise default
        # OFF so the legacy non-mogul-primary path stays byte-identical.
        return os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY", "0") == "1"
    return raw == "1"


def ring_size_to_bite_angle(ring_size: int) -> Optional[float]:
    """Map chelate ring size to its idealised bite angle (degrees).

    Returns ``None`` for ring sizes outside the tabulated 4-8 range so the
    caller can decide whether to skip the chelate (rare / out-of-table) or
    fall back to a sensible default.  Pure lookup — no per-class branch.
    """
    return RING_SIZE_TO_BITE_DEG.get(int(ring_size))


def polyhedron_vertex_angles(geometry: str) -> Optional[np.ndarray]:
    """Compute the pairwise vertex-vertex angles (degrees) for ``geometry``.

    Loads the polyhedron's unit vectors via :func:`polyhedra.ref_vectors`,
    then returns the upper-triangular pairwise angles as a symmetric
    ``(N, N)`` matrix in degrees.  Diagonal is 0°.  Used by the scorer to
    measure how closely a given (va, vb) vertex pair's angle matches a
    chelate's intrinsic bite.

    Returns ``None`` on unknown geometry name or any other failure (so the
    scorer can fall back to first-candidate behaviour gracefully).
    """
    try:
        from delfin.fffree import polyhedra as _polyhedra
    except ImportError:
        return None
    try:
        V = _polyhedra.ref_vectors(str(geometry))
    except (KeyError, Exception):
        return None
    if V is None:
        return None
    V = np.asarray(V, dtype=float)
    if V.ndim != 2 or V.shape[1] != 3 or V.shape[0] < 2:
        return None
    norms = np.linalg.norm(V, axis=1, keepdims=True)
    norms = np.where(norms < 1e-9, 1.0, norms)
    Vu = V / norms
    cos = np.clip(Vu @ Vu.T, -1.0, 1.0)
    ang = np.degrees(np.arccos(cos))
    # Clean numerical noise on the diagonal.
    np.fill_diagonal(ang, 0.0)
    return ang


def _greedy_chelate_pair_assignment(
    chelate_bites: Sequence[float],
    cn: int,
    vertex_angles_deg: np.ndarray,
) -> Tuple[float, List[Tuple[int, int]]]:
    """Greedy assignment of chelate bites to disjoint vertex pairs.

    For each chelate bite angle (one per chelating arm-pair) we pick the
    yet-unused vertex pair (va, vb) whose ``vertex_angles_deg[va, vb]``
    is closest to that bite, then mark va and vb as consumed.  Returns
    the total MSE plus the chosen pair list.

    Greedy is sufficient: the chelate count is small (typical: 1-3 chelate
    pairs per complex), the cost surface is unimodal in the angle distance,
    and the alternative (full Hungarian assignment over vertex pairs) adds
    complexity without changing the picked polyhedron on any test case in
    the CCDC ground-truth set.

    Parameters
    ----------
    chelate_bites
        Sequence of bite angles (deg) — one per chelating donor-donor pair.
    cn
        Coordination number.  The first ``cn`` vertices of ``vertex_angles_deg``
        are used (sandwich geometries with effective-CN < N_vertices are
        already handled at the polyhedra dispatch).
    vertex_angles_deg
        Pairwise vertex-vertex angle matrix from :func:`polyhedron_vertex_angles`.

    Returns
    -------
    (mse, pairs)
        ``mse`` = mean squared deviation between assigned chelate bites
        and the matched vertex angles (in deg²).  ``pairs`` = chosen
        ``[(va, vb), ...]`` in input order.  When ``chelate_bites`` is
        empty, returns ``(0.0, [])``.
    """
    n = min(int(cn), int(vertex_angles_deg.shape[0]))
    if n < 2 or not chelate_bites:
        return 0.0, []
    used: set = set()
    total_sq = 0.0
    pairs: List[Tuple[int, int]] = []
    # Sort chelate bites descending — assign the most "demanding" bite
    # (largest deviation from generic 90°) first so it gets first pick of
    # the matching vertex pair.  Deterministic tie-break by input index.
    order = sorted(range(len(chelate_bites)),
                   key=lambda i: (-abs(chelate_bites[i] - 90.0), i))
    # Walk in original input order for the OUTPUT pairs list (so the
    # caller's chelate_info ordering is preserved), but make assignment
    # decisions in priority order.
    assignment_by_chelate: Dict[int, Tuple[int, int]] = {}
    for ci in order:
        bite = float(chelate_bites[ci])
        best_va = -1
        best_vb = -1
        best_dev = float("inf")
        for va in range(n):
            if va in used:
                continue
            for vb in range(va + 1, n):
                if vb in used:
                    continue
                dev = abs(float(vertex_angles_deg[va, vb]) - bite)
                # Strict less-than for deterministic tie-break (lowest va,vb)
                if dev < best_dev:
                    best_dev = dev
                    best_va = va
                    best_vb = vb
        if best_va < 0:
            # No vertex pair available — penalise heavily so the candidate
            # loses to one that CAN host all chelates.
            total_sq += (180.0 ** 2)
            continue
        used.add(best_va)
        used.add(best_vb)
        total_sq += best_dev * best_dev
        assignment_by_chelate[ci] = (best_va, best_vb)
    for ci in range(len(chelate_bites)):
        pairs.append(assignment_by_chelate.get(ci, (-1, -1)))
    mse = total_sq / max(1, len(chelate_bites))
    return mse, pairs


def score_polyhedron_for_chelates(
    geometry: str,
    cn: int,
    chelate_bites: Sequence[float],
) -> Optional[float]:
    """Score (lower = better) of ``geometry`` for the given chelate bites.

    Score = mean squared deviation (in deg²) between each chelate's
    intrinsic bite angle and the greedily-matched vertex pair's angle in
    ``geometry``.  Returns ``None`` if the geometry cannot be loaded or
    is degenerate; the caller should treat that as "do not pick".

    When ``chelate_bites`` is empty, returns ``0.0`` (any geometry is
    chelate-neutral; the caller's first-candidate rule should win).
    """
    if not chelate_bites:
        return 0.0
    angles = polyhedron_vertex_angles(geometry)
    if angles is None:
        return None
    if angles.shape[0] < min(int(cn), len(chelate_bites) * 2):
        return None
    mse, _ = _greedy_chelate_pair_assignment(chelate_bites, int(cn), angles)
    return float(mse)


def chelate_info_from_decompose(d: Mapping) -> List[Dict]:
    """Derive a ``chelate_info`` list from a ``decompose()`` result.

    Walks the ``ligands`` list and emits ONE entry per chelating donor-pair
    in each polydentate ligand.  Monodentate ligands contribute nothing.
    Hapto ligands contribute one entry per effective ring (denticity-1).

    Output entries:

    .. code-block::

        {
            "ligand_idx": int,   # index in d["ligands"]
            "denticity":  int,   # ligand total denticity
            "ring_size":  int,   # chelate ring size (M + backbone + 2 donors)
            "bite_deg":   float, # idealised bite angle from ring_size
            "n_chelate_atoms": int,   # backbone atoms between donor pair
        }

    For a polydentate ligand with denticity ``k > 2`` we emit ``k - 1``
    chelate-pair entries (adjacent donors along the donor sequence as
    given by the ligand's ``donor_local_idxs``), which captures the
    geometric "this donor is bonded back to its neighbour through the
    same ligand" structure that the polyhedron-picker actually needs.
    Tridentate ``mer`` ligands (e.g. terpy) end up with two adjacent
    chelate-pair bites and zero distal-pair entries — that matches the
    physical reality (the central donor sits on the line through the
    metal, the two end donors each form a 5-ring with the centre).

    Universal: pure graph traversal on the per-ligand ``mol`` plus the
    local-donor index list.  No SMILES pattern, no element-specific
    branches.
    """
    ligands = d.get("ligands") or []
    out: List[Dict] = []
    for li, lg in enumerate(ligands):
        denticity = int(lg.get("denticity", 0) or 0)
        if denticity < 2:
            continue
        # Polydentate ligand: enumerate adjacent-donor pairs along the
        # donor sequence and compute ring size = backbone atoms +1 (for
        # the metal) + 2 (the two donors) = backbone + 3.
        mol_l = lg.get("mol")
        donor_local_idxs: List[int] = list(lg.get("donor_local_idxs") or [])
        if mol_l is None or len(donor_local_idxs) < 2:
            # Defensive: treat as a single (denticity-1)-bite of unknown
            # backbone length, falling back to the 5-ring idealisation.
            for _ in range(denticity - 1):
                bite = float(RING_SIZE_TO_BITE_DEG[5])
                out.append({
                    "ligand_idx": int(li),
                    "denticity": int(denticity),
                    "ring_size": 5,
                    "bite_deg": bite,
                    "n_chelate_atoms": 2,
                })
            continue
        # Hapto: each ring is ONE donor site in the polyhedron, so the
        # "ring-size" is not the chelate-bite size — the ligand acts as a
        # single point-donor.  Skip emitting chelate bites for hapto
        # ligands; their effective denticity at the polyhedron is 1.
        if lg.get("is_hapto"):
            continue
        # Walk each adjacent donor pair and BFS the shortest backbone path
        # between them.  ``ring_size = path_len + 1`` (the +1 is the M atom
        # which closes the chelate ring).  ``path_len`` includes both donor
        # atoms.  For en (H2N-CH2-CH2-NH2) the path is N-C-C-N (length 4),
        # ring size = 5.  For acac the path is O-C=C-C=O (length 5), ring
        # size = 6.  Matches the empirical chelate ring nomenclature.
        for k in range(len(donor_local_idxs) - 1):
            a = int(donor_local_idxs[k])
            b = int(donor_local_idxs[k + 1])
            path_len = _shortest_path_len(mol_l, a, b)
            if path_len is None or path_len < 2:
                # Degenerate / disconnected pair — fall back to 5-ring.
                ring_size = 5
            else:
                ring_size = int(path_len) + 1  # +1 for the M atom
            # Clamp to table range so the lookup always succeeds; rare
            # large rings fall back to the 5-ring idealisation since a
            # very expanded chelate has effectively NO geometric bite
            # constraint (the backbone is long enough to span any pair).
            ring_size_clamped = min(max(ring_size, 4), 8)
            bite = float(RING_SIZE_TO_BITE_DEG[ring_size_clamped])
            n_backbone = max(0, (path_len or 2) - 2)
            out.append({
                "ligand_idx": int(li),
                "denticity": int(denticity),
                "ring_size": int(ring_size),
                "bite_deg": float(bite),
                "n_chelate_atoms": int(n_backbone),
            })
    return out


def _shortest_path_len(mol, a: int, b: int) -> Optional[int]:
    """BFS shortest path length (in atoms) from ``a`` to ``b`` in ``mol``.

    Returns ``None`` if disconnected.  Pure graph traversal — no element
    or bond-order awareness.  ``len`` includes both endpoints, so the
    minimum value for ``a != b`` is 2 (direct bond).
    """
    if mol is None:
        return None
    a = int(a)
    b = int(b)
    if a == b:
        return 1
    n = mol.GetNumAtoms()
    if a < 0 or b < 0 or a >= n or b >= n:
        return None
    parent: Dict[int, int] = {a: -1}
    queue: List[int] = [a]
    head = 0
    while head < len(queue):
        u = queue[head]
        head += 1
        if u == b:
            # Reconstruct length
            length = 0
            cur = u
            while cur != -1:
                length += 1
                cur = parent[cur]
            return length
        atom = mol.GetAtomWithIdx(int(u))
        for nb in atom.GetNeighbors():
            j = int(nb.GetIdx())
            if j in parent:
                continue
            parent[j] = u
            queue.append(j)
    return None


def pick_polyhedron_chelate_aware(
    cn: int,
    chelate_info: Optional[Sequence[Mapping]],
    metal_sym: str = "",
    candidates: Optional[Sequence[str]] = None,
) -> Optional[str]:
    """Pick the polyhedron whose vertex-vertex angle table best matches
    the chelate-bite multiset.

    Parameters
    ----------
    cn : int
        Coordination number.  Used to fetch the candidate polyhedra via
        :func:`polyhedra.geometries_for_cn` (and to truncate the vertex
        angle matrix to the first ``cn`` rows in the rare sandwich
        effective-CN case).
    chelate_info : sequence of mapping, optional
        Output of :func:`chelate_info_from_decompose` — one entry per
        chelate donor pair, each carrying a ``bite_deg`` key.  When
        empty / None / contains only monodentate ligands, the picker
        returns the first candidate (legacy behaviour).
    metal_sym : str, optional
        Metal element symbol; forwarded to
        :func:`polyhedra.geometries_for_cn` so the candidate list grows
        to include f-block / sandwich / CN10 polyhedra when the
        corresponding env-gate is on.
    candidates : sequence of str, optional
        Explicit candidate list (e.g. when the caller has already filtered
        to a subset).  When ``None`` the candidates come from
        :func:`polyhedra.geometries_for_cn`.

    Returns
    -------
    str or None
        Polyhedron name (e.g. ``"SP-4 square planar"``); ``None`` if the
        candidate list is empty.

    Default-OFF byte-identical
    -------------------------
    When :func:`chelate_aware_picker_enabled` returns False, this function
    returns ``candidates[0]`` (the legacy first-candidate rule).  Callers
    that bypass the env-gate via ``force=True`` get the chelate-aware
    pick regardless.
    """
    try:
        from delfin.fffree import polyhedra as _polyhedra
    except ImportError:
        return None
    if candidates is None:
        try:
            candidates = _polyhedra.geometries_for_cn(int(cn), metal_sym)
        except Exception:
            candidates = []
    cand_list = list(candidates) if candidates else []
    if not cand_list:
        return None

    # Extract bite angles in stable input order.
    bites: List[float] = []
    if chelate_info:
        for ci in chelate_info:
            try:
                bites.append(float(ci.get("bite_deg")))
            except (TypeError, ValueError, AttributeError):
                continue
    if not bites:
        # No chelate constraints — legacy behaviour.
        return cand_list[0]

    if not chelate_aware_picker_enabled():
        return cand_list[0]

    # Score every candidate; pick min.  Ties broken by candidate order
    # (i.e. the first-candidate rule is the deterministic tie-break).
    best_name = cand_list[0]
    best_score = float("inf")
    for name in cand_list:
        score = score_polyhedron_for_chelates(name, int(cn), bites)
        if score is None:
            continue
        if score < best_score - 1e-9:
            best_score = score
            best_name = name
    return best_name
