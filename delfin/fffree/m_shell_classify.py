"""delfin.fffree.m_shell_classify — Universal classifier for extra atoms
inside a metal's geometric coordination shell.

When the rotamer-topology gate flags a metal atom whose geometric shell
contains MORE heavy atoms than the SMILES-graph CN (``m_shell_overfill``),
that excess is one of two things:

  A) Genuine alternative donor (CHEMICALLY VALID, KEEP)
     -- a heteroatom (or carbene-C) carrying at least one lone pair that
        sits within bite-compatible range of an already-designated donor
        on the same ligand fragment.  These are emergent κⁿ states:
        κ²-acetate emerging from a κ¹ SMILES, hemilabile donor closing,
        agostic α/β CH brushing the shell, etc.  These MUST NOT be
        rejected — they are real chemistry that the SMILES under-described.

  B) Spurious overfill (BUG, REJECT)
     -- an atom WITHOUT a free lone pair (sp3-C backbone, plain H, etc.)
        OR an atom outside any bite-compatible range of the designated
        donor set on its fragment.  This is the soft-DG-drift signature:
        a ligand backbone carbon randomly sliding into the inner sphere,
        a hydrogen "tunnelling" past the M-D edge, a non-coordinating
        fragment that happens to lie near M.

The classifier is **pure**: it reads the mol graph, element symbols and
3-D coordinates only.  No SMARTS, no per-anion table, no SMILES pattern.
The lone-pair test is the textbook Lewis octet count (re-used from
``ambidentate_kappa_enum``).  The bite-compatibility test is graph-bond
distance ≤ ``MAX_GRAPH_HOPS_TO_DONOR`` AND through-space distance to the
nearest designated donor ≤ ``MAX_GEOM_DIST_TO_DONOR``.

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import os
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Thresholds — universal, geometry-only
# ---------------------------------------------------------------------------

#: Multiplier on the ideal bond length used to define the metal's
#: geometric shell.  Matches ``rotamer_topology_gate.M_SHELL_FACTOR``.
M_SHELL_FACTOR = 1.30

#: Hard ceiling on the M-X distance (Å) above which we never even consider
#: the atom as being in the shell.  Matches the user-spec 2.6 Å in the
#: mission brief for the brute-force "atoms within ..." pass.
M_SHELL_CEILING = 2.6

#: Graph-bond distance (hops on the molecular graph EXCLUDING the metal)
#: from a candidate extra atom to ANY designated donor on the same ligand.
#: 4 hops covers the same envelope as the κⁿ enumerator (carboxylate-O,
#: nitrate-O, salicylate-O, phenanthroline-N, bipyridine-N, hemilabile
#: ether-O 3 hops from a phosphine seed).
MAX_GRAPH_HOPS_TO_DONOR = 4

#: Through-space distance ceiling (Å) between a candidate extra atom and
#: the nearest designated donor on the same ligand.  3.5 Å is the upper
#: bite envelope used by ``ambidentate_kappa_enum`` (CCDC manifold,
#: 2026-06-06 Mogul-DG Phase A).  Beyond this distance the "chelate
#: extension" interpretation breaks down.
MAX_GEOM_DIST_TO_DONOR = 3.5


# ---------------------------------------------------------------------------
# Env-gate
# ---------------------------------------------------------------------------
def _env_on() -> bool:
    """Return True iff the classifier should soften the m-shell overfill
    decision.  Default ON when ``DELFIN_FFFREE_MOGUL_PRIMARY=1`` is set.
    Explicit ``DELFIN_FFFREE_M_SHELL_CLASSIFY=0`` disables.
    """
    raw = os.environ.get("DELFIN_FFFREE_M_SHELL_CLASSIFY")
    if raw is None:
        if os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY", "0") == "1":
            return True
        return False
    return str(raw).strip() in ("1", "true", "True", "on", "yes", "YES")


# ---------------------------------------------------------------------------
# Lone-pair test (re-exported from ambidentate_kappa_enum for a single
# source of truth on Lewis-octet semantics).
# ---------------------------------------------------------------------------
def _has_lone_pair(atom) -> bool:
    """Universal Lewis lone-pair test.  Delegates to the canonical
    implementation in :mod:`delfin.fffree.ambidentate_kappa_enum` so the
    classifier shares one source of truth.
    """
    try:
        from delfin.fffree.ambidentate_kappa_enum import _has_lone_pair as _akp
        return bool(_akp(atom))
    except Exception:
        # Conservative fallback: only common donor heteroatoms accepted
        # without the full Lewis count.  Returning False here would
        # over-reject (every extra atom becomes "spurious").  Returning
        # True here would over-accept.  We prefer the conservative
        # over-reject path so the rotamer-gate stays strict on import
        # failure.
        return False


def _is_metal(s: str) -> bool:
    try:
        from delfin._bond_decollapse import _is_metal as _bd_is_metal
        return bool(_bd_is_metal(str(s)))
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Graph distance helpers
# ---------------------------------------------------------------------------
def _graph_distance_excluding_metal(
    mol, src: int, dst: int, metal_idx: int,
) -> int:
    """Shortest-path length (in bonds) from ``src`` to ``dst`` in the
    molecular graph WITHOUT traversing the metal atom.  Returns -1 when no
    path exists in the metal-removed subgraph (different ligand fragments)
    or when either endpoint coincides with the metal.
    """
    if src == metal_idx or dst == metal_idx:
        return -1
    if src == dst:
        return 0
    try:
        n = mol.GetNumAtoms()
    except Exception:
        return -1
    # BFS, skipping the metal node entirely.
    visited: Set[int] = {int(src)}
    frontier: List[int] = [int(src)]
    depth = 0
    while frontier:
        depth += 1
        nxt: List[int] = []
        for u in frontier:
            try:
                atom = mol.GetAtomWithIdx(int(u))
            except Exception:
                continue
            for nb in atom.GetNeighbors():
                v = int(nb.GetIdx())
                if v == metal_idx or v in visited:
                    continue
                if v == dst:
                    return depth
                visited.add(v)
                nxt.append(v)
        frontier = nxt
        if depth > 64:
            break
    return -1


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def atoms_in_metal_shell(
    syms: Sequence[str],
    P: np.ndarray,
    metal_idx: int,
    *,
    factor: float = M_SHELL_FACTOR,
    ceiling: float = M_SHELL_CEILING,
) -> List[int]:
    """Return the indices of all heavy (non-H, non-metal-self) atoms whose
    Euclidean distance to ``P[metal_idx]`` is below
    ``max(factor * ideal_bond(M, X), ceiling)``.

    The classifier shell is INTENTIONALLY broader than the strict
    ``rotamer_topology_gate.M_SHELL_FACTOR * ideal`` cutoff: we want to
    catch every heavy atom the user-eye would consider "in the inner
    sphere" (per spec: within ~2.6 Å of the metal), including atoms that
    sit just outside the strict M-D ideal but well inside the soft
    shell envelope of the user's mission brief.

    Pure geometry, deterministic, no graph required.
    """
    try:
        from delfin._bond_decollapse import _ideal_bond
    except Exception:
        def _ideal_bond(a: str, b: str) -> float:  # noqa: E306
            return 1.55
    P_arr = np.asarray(P, dtype=float)
    n = len(syms)
    out: List[int] = []
    if metal_idx < 0 or metal_idx >= n:
        return out
    pm = P_arr[metal_idx]
    sm = str(syms[metal_idx])
    for j in range(n):
        if j == metal_idx:
            continue
        sj = str(syms[j])
        if sj == "H":
            continue
        try:
            d = float(np.linalg.norm(P_arr[j] - pm))
        except Exception:
            continue
        if not np.isfinite(d):
            continue
        ideal = 1.55
        try:
            ideal = float(_ideal_bond(sm, sj))
        except Exception:
            pass
        strict = factor * ideal if ideal > 0.0 else ceiling
        cutoff = max(strict, ceiling)
        if d < cutoff:
            out.append(int(j))
    return sorted(out)


def classify_m_shell_extras(
    syms: Sequence[str],
    P: np.ndarray,
    metal_idx: int,
    designated_donors: Iterable[int],
    mol=None,
) -> Dict[str, List[int]]:
    """Classify the *extra* heavy atoms inside ``metal_idx``'s coordination
    shell into ``valid_alternatives`` and ``spurious``.

    Parameters
    ----------
    syms : sequence of element symbols, length N.
    P : (N, 3) coordinate array.
    metal_idx : index of the metal whose shell we are auditing.
    designated_donors : iterable of atom indices that the SMILES graph
        marks as the metal's primary donors (i.e. ``M-D`` edges in the
        expected-bond set).
    mol : optional RDKit Mol.  Required for graph-distance / lone-pair
        checks.  When ``None`` the classifier returns every extra in the
        ``spurious`` bucket (safe — we cannot prove emergence without the
        graph).

    Returns
    -------
    dict with three lists of atom indices (all sorted ascending):
        ``valid_alternatives`` -- extras that pass the lone-pair AND
            bite-compatibility test (legitimate κⁿ-emergence candidates,
            chelate extensions, hemilabile arms).
        ``spurious`` -- extras that fail either test (random soft-DG
            drift, sp3-backbone, far-out atoms).
        ``shell`` -- the full set of heavy atoms inside the geometric
            shell (designated_donors UNION extras).

    Determinism
    -----------
    All returned lists are sorted, so the same inputs always produce the
    same output ordering.  No randomness anywhere in the path.
    """
    designated_set: Set[int] = {int(d) for d in designated_donors}
    shell = atoms_in_metal_shell(syms, P, int(metal_idx))
    extras = [a for a in shell if a not in designated_set]
    valid_alternatives: List[int] = []
    spurious: List[int] = []

    if mol is None or not extras:
        return {
            "valid_alternatives": valid_alternatives,
            "spurious": sorted(extras),
            "shell": shell,
        }

    P_arr = np.asarray(P, dtype=float)
    for ex in extras:
        try:
            atom = mol.GetAtomWithIdx(int(ex))
        except Exception:
            spurious.append(ex)
            continue

        # Step 1: lone-pair gate.  No lone pair -> spurious (sp3-C
        # backbone, plain H already excluded above).
        try:
            has_lp = _has_lone_pair(atom)
        except Exception:
            has_lp = False
        if not has_lp:
            spurious.append(ex)
            continue

        # Step 2: bite-compatibility — must be graph-close AND space-close
        # to at least one designated donor on the SAME ligand fragment.
        # "Same ligand fragment" is enforced by the metal-excluded BFS
        # in _graph_distance_excluding_metal — a finite hop count means
        # the path exists without going through the metal, i.e. both
        # atoms belong to the same ligand.
        bite_ok = False
        for dd in designated_set:
            hops = _graph_distance_excluding_metal(
                mol, int(ex), int(dd), int(metal_idx),
            )
            if hops < 0 or hops > MAX_GRAPH_HOPS_TO_DONOR:
                continue
            try:
                d_geom = float(np.linalg.norm(P_arr[ex] - P_arr[dd]))
            except Exception:
                continue
            if not np.isfinite(d_geom):
                continue
            if d_geom <= MAX_GEOM_DIST_TO_DONOR:
                bite_ok = True
                break

        if bite_ok:
            valid_alternatives.append(ex)
        else:
            spurious.append(ex)

    return {
        "valid_alternatives": sorted(valid_alternatives),
        "spurious": sorted(spurious),
        "shell": shell,
    }


def designated_donors_for_metal(
    mol, metal_idx: int,
) -> List[int]:
    """Return the heavy-atom neighbours of the metal in the graph.

    These are the SMILES-declared M-D edges.  Used by callers that do
    not pass the designated_donors set explicitly.
    """
    if mol is None:
        return []
    try:
        atom = mol.GetAtomWithIdx(int(metal_idx))
    except Exception:
        return []
    out: List[int] = []
    for nb in atom.GetNeighbors():
        try:
            if int(nb.GetAtomicNum()) == 1:
                continue
            out.append(int(nb.GetIdx()))
        except Exception:
            continue
    return sorted(out)


def shell_audit(
    syms: Sequence[str],
    P: np.ndarray,
    mol,
    *,
    expected_cn: Optional[Sequence[int]] = None,
) -> List[Dict]:
    """Per-metal audit of every metal's coordination shell.

    Returns a list (one dict per metal atom in the structure) with
    fields:
        ``metal_idx`` : index
        ``metal_sym`` : element
        ``expected_cn`` : -1 if unknown
        ``designated_donors`` : list of indices (SMILES M-D edges)
        ``shell_atoms`` : list of indices in geometric shell
        ``valid_alternatives`` : list of emergent-κⁿ candidates
        ``spurious`` : list of drift atoms
        ``verdict`` : one of "ok" / "emergent_kappa" / "spurious_overfill"
            * "ok"            : len(shell) <= expected_cn
            * "emergent_kappa": shell > expected_cn but ALL extras are
                                valid alternatives -> KEEP, tag as
                                "emergent_kappa".
            * "spurious_overfill" : at least one extra is spurious ->
                                TRIGGER ROLLBACK.

    Determinism: list is sorted by metal_idx, every sub-list is sorted.
    """
    out: List[Dict] = []
    n = len(syms)
    for i in range(n):
        if not _is_metal(str(syms[i])):
            continue
        dds = designated_donors_for_metal(mol, i)
        cls = classify_m_shell_extras(syms, P, i, dds, mol=mol)
        shell = cls["shell"]
        valid = cls["valid_alternatives"]
        spur = cls["spurious"]
        exp = -1
        if expected_cn is not None and i < len(expected_cn):
            try:
                exp = int(expected_cn[i])
            except Exception:
                exp = -1
        if exp < 0 or len(shell) <= exp:
            verdict = "ok"
        elif spur:
            verdict = "spurious_overfill"
        else:
            verdict = "emergent_kappa"
        out.append({
            "metal_idx": int(i),
            "metal_sym": str(syms[i]),
            "expected_cn": exp,
            "designated_donors": list(dds),
            "shell_atoms": list(shell),
            "valid_alternatives": list(valid),
            "spurious": list(spur),
            "verdict": verdict,
        })
    return out


# ---------------------------------------------------------------------------
# Smoke check
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Quick sanity: a κ²-acetate emergence on a κ¹ SMILES should land in
    # valid_alternatives; a random sp3-C drifting into the shell should
    # land in spurious.
    from rdkit import Chem
    # Cu-OAc with the second O drifting into the shell.
    m = Chem.MolFromSmiles("[Cu]OC(=O)C")  # κ¹ in SMILES
    m = Chem.AddHs(m)
    # Place atoms manually for the sanity test.
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    P = np.zeros((m.GetNumAtoms(), 3), dtype=float)
    idx = {a.GetIdx(): a.GetSymbol() for a in m.GetAtoms()}
    # Cu at origin; O1 (carboxylate-O bound) at +x; C(=O) at +x+1; O2 (κ²
    # emergence) at +x+1, +y; CH3 backbone behind.
    cu = next(i for i, s in idx.items() if s == "Cu")
    P[cu] = [0.0, 0.0, 0.0]
    # Designated donor = O bonded to Cu
    cu_atom = m.GetAtomWithIdx(cu)
    o1 = next(int(nb.GetIdx()) for nb in cu_atom.GetNeighbors()
              if nb.GetSymbol() == "O")
    P[o1] = [2.0, 0.0, 0.0]
    o1_atom = m.GetAtomWithIdx(o1)
    c_central = next(int(nb.GetIdx()) for nb in o1_atom.GetNeighbors()
                     if nb.GetSymbol() == "C")
    P[c_central] = [3.0, 0.6, 0.0]
    c_atom = m.GetAtomWithIdx(c_central)
    other_O = [int(nb.GetIdx()) for nb in c_atom.GetNeighbors()
               if nb.GetSymbol() == "O" and int(nb.GetIdx()) != o1]
    if other_O:
        # Put the κ² O at 2.3 Å from Cu (inside the shell).
        P[other_O[0]] = [2.1, 1.0, 0.0]
    methyl_c = [int(nb.GetIdx()) for nb in c_atom.GetNeighbors()
                if nb.GetSymbol() == "C"]
    if methyl_c:
        P[methyl_c[0]] = [4.0, 0.6, 0.0]
    res = classify_m_shell_extras(
        syms, P, cu, designated_donors=[o1], mol=m,
    )
    print("Cu-OAc κ² emergence test:", res)
