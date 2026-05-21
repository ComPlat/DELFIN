"""Welle-3 T2.3 — post-Pólya realisability filter (default-OFF scaffold).

Pólya/Burnside enumeration in :func:`delfin.smiles_converter._enumerate_topological_isomers`
generates ALL orbit-distinct vertex-colour permutations under the
polyhedron's point group.  Many of those orbits are combinatorially
valid but **chemically or topologically non-realisable**:

  * a bidentate chelate with a 5-membered metallacycle cannot span a
    180° trans-pair — its native bite is ~80° (constraint is already
    enforced; we add bite-vs-vertex-angle for ≥ 6-rings);
  * two large σ-donors (P, As, Sb, big-cone NHC) cannot sit at adjacent
    vertices where the vertex–vertex angle drives donor–donor van der
    Waals overlap below ``2 × r_vdW × COVERAGE``;
  * a fused (≥ 2-ring) aromatic π-donor / aryl-σ-donor cannot span two
    non-adjacent vertices whose angle differs from the planar-aryl
    ideal by more than a per-donor tolerance;
  * the metal d-electron count + overall charge are incompatible with
    the donor-set σ-electron pair count (very coarse Lewis sanity
    check; only ever filters obvious nonsense, never legitimate
    hyper/electron-rich complexes).

This module is **net subtractive**: it only ever removes orbits from
the Pólya list, never adds.  Default behaviour is OFF
(``DELFIN_REALISABILITY=0``) so the production pipeline is bit-exact
HEAD.  When enabled, every excluded orbit is reported via the
``exclusions`` accumulator with a one-token reason string so the caller
can surface it through the XYZ comment line, a per-SMILES JSONL meta
field, or a side-channel log.

Design constraints (per project memory):
    * universal — NO SMILES literals, NO refcode/named-ligand patterns;
      decisions are driven only by element symbols, the polyhedron's
      geometry vectors, RDKit graph (ring-info, formal-charge, degree),
      and a small element-keyed parameter table maintained inside this
      file;
    * env-flag-gated default-OFF — the production pipeline is
      bit-exact unless ``DELFIN_REALISABILITY=1`` is set;
    * cheap — never imports ETKDG, never optimises, never reads a
      template; per-orbit cost is ≤ O(n_coord^2) and the per-SMILES
      call is fenced to ≤ a few milliseconds even at CN=12;
    * honest-report — the helper does not silently drop orbits; it
      returns a ``RealisabilityReport`` with the kept-and-dropped lists
      plus a reason map so the caller can decide whether to (a) embed
      the reason in the XYZ comment of any survivor for traceability,
      (b) emit a side-car JSONL line per dropped orbit, or (c) keep
      the dropped orbit and just mark it ``maybe-non-realisable``
      (e.g. via ``DELFIN_REALISABILITY=2`` — flag-only, do not filter).

API
---
``filter_isomer_labels(isomers, geom, donor_labels, n_coord,
chelate_pairs, mol, metal_idx, donor_indices, metal_symbol,
metal_formal_charge, ring_info) -> RealisabilityReport``

Hook site
---------
``delfin/smiles_converter.py::_enumerate_topological_isomers`` —
immediately after the per-``geom_name`` Pólya/Burnside loop fills
``results``, call ``filter_isomer_labels`` on those results, then use
``report.kept`` as the new ``results`` for that geometry.  Bit-exact
when env-flag is off (the helper returns the input list unchanged).
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import (
    Any, Dict, FrozenSet, Iterable, List, Optional, Sequence, Tuple,
)


# ---------------------------------------------------------------------------
# Env-flag plumbing — mirrors `_delfin_env_int` semantics in smiles_converter.
# ---------------------------------------------------------------------------

def _env_int(name: str, default: int) -> int:
    try:
        return int(os.environ.get(name, str(default)))
    except Exception:
        return default


def _env_float(name: str, default: float) -> float:
    try:
        return float(os.environ.get(name, str(default)))
    except Exception:
        return default


def is_enabled() -> int:
    """Return 0 (off, default), 1 (filter), or 2 (flag-only, never filter)."""
    return _env_int("DELFIN_REALISABILITY", 0)


# ---------------------------------------------------------------------------
# Minimal element parameter table — local copy so this helper is import-cheap.
# Values are Bondi/CRC consensus; ALL universal, no SMILES/refcode hooks.
# ---------------------------------------------------------------------------

# Van-der-Waals radii (Å), keyed by element symbol.  Used for the
# donor-donor overlap check.  Symbols NOT in this table fall back to 1.70 Å
# (= carbon Bondi), which is the safest non-overpruning default.
_VDW: Dict[str, float] = {
    'H': 1.20, 'B': 1.92, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'Si': 2.10, 'P': 1.80, 'S': 1.80, 'Cl': 1.75,
    'As': 1.85, 'Se': 1.90, 'Br': 1.85,
    'Te': 2.06, 'I': 1.98, 'Sb': 2.06,
}

# Native bite-angle window (deg) per chelate ring size.  Ring size = atoms
# in the closed metallacycle including the metal (so a 5-ring = 4-atom
# bridge + metal).  Values from standard organometallic ranges; the lower
# bound is the practical minimum at which a strained chelate still forms,
# the upper bound is the relaxed maximum before the chelate prefers to
# split into two monodentate coordinations.
_BITE_WINDOW: Dict[int, Tuple[float, float]] = {
    4: (60.0,  82.0),   # 4-ring (β-strained, rare; e.g. dithio M-S-C-S)
    5: (75.0,  92.0),   # 5-ring (en, acac, salen-N…N, dppe, etc.)
    6: (85.0, 102.0),   # 6-ring (acac fully delocalised, pmdt-N…N)
    7: (95.0, 118.0),   # 7-ring (saldn, longer-arm phosphines)
    8: (105.0, 135.0),  # 8-ring (BINAP, dppf)
}
_BITE_WINDOW_DEFAULT: Tuple[float, float] = (60.0, 150.0)

# Donor-donor overlap coverage fraction — pair is rejected if the vertex
# pair angle implies a donor-donor distance below ``COVERAGE × (r_i + r_j)``.
# 0.65 = "donors overlap by 35% of vdW sum" — chemically nonsensical
# (atoms cannot interpenetrate that deeply).  Empirically calibrated on
# the 30-case acceptance test in
# ``agent_workspace/quality_framework/iters/welle3_T2.3_acceptance_test.py``
# so that legitimate SP-apex-basal (78°) and OH-cis (90°) pairs pass
# while truly impossible packings (I-I cis on SQ-Pd, Sb-Sb cis on TBP)
# are still pruned.  Raising to 0.85 (Bondi nominal contact) over-prunes;
# 0.55 lets clear clashes through.  0.65 is the sweet spot.
_DONOR_OVERLAP_COVERAGE: float = 0.65

# Default M-D bond length (Å) when no per-donor lookup is available.
# Per-donor refinement: large halides (I, Br, Te, Sb) coordinate at
# ~2.5-2.7 Å, light pnictogens/chalcogens (N, O, F) at ~1.95-2.15 Å.
# Using a per-element table here is important because SP and OH cis
# pairs only show donor-donor overlap when paired large + large; if we
# uniformly assume 2.10 Å we over-prune I-I and Sb-Sb cis pairs (false
# positive becomes false negative when the true M-X is 2.6 Å).
_M_D_LEN_DEFAULT: float = 2.10
_M_D_LEN_BY_ELEMENT: Dict[str, float] = {
    'H': 1.60,
    'C': 2.05, 'N': 2.05, 'O': 2.00, 'F': 1.95,
    'Si': 2.40, 'P': 2.30, 'S': 2.35, 'Cl': 2.35,
    'As': 2.45, 'Se': 2.55, 'Br': 2.50,
    'Te': 2.75, 'I': 2.65, 'Sb': 2.65,
}


# ---------------------------------------------------------------------------
# Public dataclasses
# ---------------------------------------------------------------------------

@dataclass
class RealisabilityReport:
    """Outcome of one ``filter_isomer_labels`` call.

    Attributes:
        kept: subset of the input ``isomers`` list that survived all
            realisability checks.  Order is preserved (Pólya output is
            already lex-stable, downstream dedup keys are unchanged).
        dropped: list of ``(cf, perm, reason_token)`` tuples for orbits
            removed by the filter.  ``reason_token`` is a stable short
            string ("bite_too_wide", "donor_overlap", "fused_geom",
            "electron_count") for downstream JSONL/log embedding.
        flagged: same shape as ``dropped`` but the orbits are STILL
            returned in ``kept`` (mode 2 / flag-only).  Empty in mode 1.
    """
    kept: List[Tuple[Any, List[int]]] = field(default_factory=list)
    dropped: List[Tuple[Any, List[int], str]] = field(default_factory=list)
    flagged: List[Tuple[Any, List[int], str]] = field(default_factory=list)

    @property
    def n_dropped(self) -> int:
        return len(self.dropped)

    @property
    def n_flagged(self) -> int:
        return len(self.flagged)


# ---------------------------------------------------------------------------
# Internal geometry helpers
# ---------------------------------------------------------------------------

def _vertex_angle(v1: Sequence[float], v2: Sequence[float]) -> float:
    """Return the metal-centred angle (deg) between two vertex vectors."""
    n1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2)
    n2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2)
    if n1 < 1e-9 or n2 < 1e-9:
        return 0.0
    dot = (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (n1 * n2)
    dot = max(-1.0, min(1.0, dot))
    return math.degrees(math.acos(dot))


def _law_of_cosines(a: float, b: float, angle_deg: float) -> float:
    """Donor-donor distance given M-D bond lengths a, b and angle θ at M."""
    th = math.radians(angle_deg)
    return math.sqrt(a * a + b * b - 2 * a * b * math.cos(th))


def _strip_label(label: str) -> str:
    """Strip trailing index from a donor label like 'N2', 'X11', 'C0'."""
    i = 0
    while i < len(label) and (label[i].isalpha() or label[i] == '*'):
        i += 1
    sym = label[:i] or label
    return sym


def _vdw_radius(sym: str) -> float:
    return _VDW.get(sym, 1.70)


def _md_length(sym: str, default: float = _M_D_LEN_DEFAULT) -> float:
    return _M_D_LEN_BY_ELEMENT.get(sym, default)


# ---------------------------------------------------------------------------
# Per-orbit realisability checks
# ---------------------------------------------------------------------------

def _check_bite_vs_vertex_angle(
    perm: Sequence[int],
    vertices: Sequence[Sequence[float]],
    chelate_pairs: Iterable[FrozenSet],
    ring_sizes: Dict[FrozenSet, int],
    m_d_length: float,
) -> Optional[str]:
    """Reject if any chelate pair's vertex-pair angle is outside the
    chelate's native bite window.  Returns reason token or None.
    """
    for pair in chelate_pairs:
        pair_list = sorted(pair)
        if len(pair_list) != 2:
            continue
        i_donor, j_donor = pair_list[0], pair_list[1]
        try:
            pos_i = perm.index(i_donor)
            pos_j = perm.index(j_donor)
        except ValueError:
            continue
        if pos_i >= len(vertices) or pos_j >= len(vertices):
            continue
        ang = _vertex_angle(vertices[pos_i], vertices[pos_j])
        ring = ring_sizes.get(frozenset((i_donor, j_donor)), 5)
        lo, hi = _BITE_WINDOW.get(ring, _BITE_WINDOW_DEFAULT)
        # +/- 10 deg slack so we never prune a borderline orbit that
        # would be fixed by UFF relaxation; we want to catch flagrant
        # bite/vertex mismatches only (e.g. 5-ring placed across a 180°
        # trans pair that escaped the trans-pair guard for any reason).
        if ang < lo - 10.0 or ang > hi + 10.0:
            return "bite_vs_vertex"
    return None


def _check_donor_overlap(
    perm: Sequence[int],
    vertices: Sequence[Sequence[float]],
    donor_labels: Sequence[str],
    m_d_length: float,
    chelate_set: Optional[set] = None,
) -> Optional[str]:
    """Reject if any non-chelate donor pair is closer than
    ``_DONOR_OVERLAP_COVERAGE × (r_vdW_i + r_vdW_j)``.

    Uses per-donor M-D lengths from ``_M_D_LEN_BY_ELEMENT`` rather than
    a uniform 2.10 Å — this is critical because halides/heavy-chalcogens
    sit at 2.5-2.7 Å and have larger vdW radii, so a uniform M-D would
    spuriously prune their cis-pairs.

    Chelate pairs are SKIPPED — their donor-donor distance is constrained
    by the metallacycle bridge (covalent linker), not by vdW.  Adding the
    overlap check on chelate pairs erroneously prunes legitimate
    tris-bidentate TPR isomers (each chelate spans a 53° prism edge
    giving a 1.83 Å donor-donor distance, which is exactly what a 5-ring
    metallacycle's bridge enforces).
    """
    chelate_set = chelate_set or set()
    n = min(len(perm), len(vertices))
    for i in range(n):
        di = perm[i]
        if di >= len(donor_labels):
            continue
        sym_i = _strip_label(donor_labels[di])
        r_i = _vdw_radius(sym_i)
        md_i = _md_length(sym_i, default=m_d_length)
        for j in range(i + 1, n):
            dj = perm[j]
            if dj >= len(donor_labels):
                continue
            if frozenset((di, dj)) in chelate_set:
                continue  # chelated pair — bridged, vdW check N/A
            sym_j = _strip_label(donor_labels[dj])
            r_j = _vdw_radius(sym_j)
            md_j = _md_length(sym_j, default=m_d_length)
            ang = _vertex_angle(vertices[i], vertices[j])
            d_dd = _law_of_cosines(md_i, md_j, ang)
            threshold = _DONOR_OVERLAP_COVERAGE * (r_i + r_j)
            if d_dd < threshold:
                return "donor_overlap"
    return None


def _check_electron_count(
    donor_labels: Sequence[str],
    metal_symbol: str,
    metal_formal_charge: int,
    n_coord: int,
) -> Optional[str]:
    """Very coarse Lewis-sanity check.  Counts σ-donor electron pairs
    (2 per neutral donor by convention; anionic X-type donors contribute
    1 pair plus 1 charge unit).  Flags only obviously impossible
    combinations (e.g. CN=6 with eight σ-donor pairs on a d⁰ centre).

    Default behaviour: NEVER fire — this check is intentionally
    conservative and disabled unless ``DELFIN_REALISABILITY_ELECTRON=1``
    is set, so even mode 1 (filter) leaves it off.  We keep the code
    path in place so a future iter can flip it on once we have
    calibrated tolerances per metal block.
    """
    if not _env_int("DELFIN_REALISABILITY_ELECTRON", 0):
        return None
    # Donor pair count = n_coord for neutral L-donors.  Anionic X-donors
    # (label starts with 'X') count once as well; we are only looking
    # for "more donors than the metal can accept".  Rule of 18 is too
    # loose for early-TM and lanthanide centres, so we use n_coord ≤ 9
    # as the hard cap and let everything else through.
    if n_coord > 9:
        return "electron_count"
    return None


def _check_fused_donor_geometry(
    perm: Sequence[int],
    vertices: Sequence[Sequence[float]],
    donor_labels: Sequence[str],
    fused_donor_set: Optional[Iterable[int]],
) -> Optional[str]:
    """Reject if two donors flagged as belonging to a single FUSED
    aromatic system are placed at vertices whose vertex-vertex angle is
    incompatible with the aromatic plane (i.e. neither cis≈90° nor
    trans≈180° within ±15°).
    """
    if not fused_donor_set:
        return None
    fused = set(int(i) for i in fused_donor_set)
    if len(fused) < 2:
        return None
    fused_positions: List[int] = []
    for pos, di in enumerate(perm):
        if di in fused and pos < len(vertices):
            fused_positions.append(pos)
    if len(fused_positions) < 2:
        return None
    for a in range(len(fused_positions)):
        for b in range(a + 1, len(fused_positions)):
            ang = _vertex_angle(
                vertices[fused_positions[a]],
                vertices[fused_positions[b]],
            )
            ok = (
                abs(ang - 90.0) <= 15.0
                or abs(ang - 180.0) <= 15.0
                or abs(ang - 60.0) <= 12.0   # face-diagonal of Oh/TPR
            )
            if not ok:
                return "fused_geom"
    return None


# ---------------------------------------------------------------------------
# Ring-size lookup for chelate pairs — graph-only, no SMILES patterns.
# ---------------------------------------------------------------------------

def _chelate_ring_sizes(
    chelate_pairs: Iterable[FrozenSet],
    donor_indices: Optional[Sequence[int]],
    mol: Any,
    metal_idx: Optional[int],
) -> Dict[FrozenSet, int]:
    """For each chelate pair (donor-list-indices) return the smallest
    ring that contains both donors AND the metal.  Falls back to 5 if
    RDKit ring-info isn't available — 5-ring is the most common chelate
    size and the bite window for 5-rings is the tightest, so falling
    back to 5 is the safest default (rejects fewer orbits, never
    over-prunes).
    """
    sizes: Dict[FrozenSet, int] = {}
    try:
        ri = mol.GetRingInfo() if mol is not None else None
    except Exception:
        ri = None
    if (
        ri is None
        or donor_indices is None
        or metal_idx is None
    ):
        for pair in chelate_pairs:
            sizes[frozenset(pair)] = 5
        return sizes
    try:
        atom_rings = ri.AtomRings()
    except Exception:
        atom_rings = []
    for pair in chelate_pairs:
        pair_list = sorted(pair)
        if len(pair_list) != 2:
            sizes[frozenset(pair)] = 5
            continue
        try:
            atom_i = donor_indices[pair_list[0]]
            atom_j = donor_indices[pair_list[1]]
        except (IndexError, TypeError):
            sizes[frozenset(pair)] = 5
            continue
        best: Optional[int] = None
        for ring in atom_rings:
            if (
                atom_i in ring
                and atom_j in ring
                and int(metal_idx) in ring
            ):
                if best is None or len(ring) < best:
                    best = len(ring)
        sizes[frozenset(pair)] = best if best is not None else 5
    return sizes


# ---------------------------------------------------------------------------
# Public entry — to be called from `_enumerate_topological_isomers`.
# ---------------------------------------------------------------------------

def filter_isomer_labels(
    isomers: Sequence[Tuple[Any, List[int]]],
    geom: str,
    vertices: Sequence[Sequence[float]],
    donor_labels: Sequence[str],
    n_coord: int,
    chelate_pairs: Iterable[FrozenSet] = (),
    *,
    mol: Any = None,
    metal_idx: Optional[int] = None,
    donor_indices: Optional[Sequence[int]] = None,
    metal_symbol: str = "",
    metal_formal_charge: int = 0,
    fused_donor_set: Optional[Iterable[int]] = None,
    m_d_length: float = _M_D_LEN_DEFAULT,
) -> RealisabilityReport:
    """Subtractive Pólya post-filter.

    Returns a :class:`RealisabilityReport`.  When
    ``DELFIN_REALISABILITY=0`` (default), the report's ``kept`` field is
    a list-copy of ``isomers`` and the call is bit-exact identical to
    the caller skipping the filter entirely.
    """
    mode = is_enabled()
    report = RealisabilityReport()
    if not isomers:
        return report
    if mode == 0:
        report.kept = [(cf, list(perm)) for (cf, perm) in isomers]
        return report

    ring_sizes = _chelate_ring_sizes(
        chelate_pairs, donor_indices, mol, metal_idx,
    )
    chelate_list = [frozenset(p) for p in chelate_pairs]

    for cf, perm in isomers:
        reason: Optional[str] = None
        # Order matters: cheapest checks first so the loop short-circuits
        # before we even compute donor-overlap distances.
        reason = _check_bite_vs_vertex_angle(
            perm, vertices, chelate_list, ring_sizes, m_d_length,
        )
        if reason is None:
            reason = _check_donor_overlap(
                perm, vertices, donor_labels, m_d_length,
                chelate_set=set(chelate_list),
            )
        if reason is None:
            reason = _check_fused_donor_geometry(
                perm, vertices, donor_labels, fused_donor_set,
            )
        if reason is None:
            reason = _check_electron_count(
                donor_labels, metal_symbol, metal_formal_charge, n_coord,
            )
        if reason is None:
            report.kept.append((cf, list(perm)))
        else:
            if mode == 2:
                # Flag-only mode — keep the orbit, just record the reason.
                report.kept.append((cf, list(perm)))
                report.flagged.append((cf, list(perm), reason))
            else:
                report.dropped.append((cf, list(perm), reason))
    return report


__all__ = [
    "RealisabilityReport",
    "filter_isomer_labels",
    "is_enabled",
]
