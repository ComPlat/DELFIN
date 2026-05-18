"""Per-system treatment-matrix dispatcher (Welle-5m-Z).

Universal-fundamental architecture for replacing global env-flag default-flips
with per-system, cell-conditional activation. The dispatcher classifies each
molecule (via :mod:`delfin._system_classifier`, falling back to a local stub
when the dedicated classifier module is unavailable) and looks up matching
treatment patterns. Each pattern maps a chemistry feature predicate to a set
of env-flag overrides that should be active for that conversion only.

Master gate
-----------
The whole dispatcher is gated by the environment variable
``DELFIN_5M_TREATMENT_MATRIX`` (``0`` / unset disables — bit-exact OFF).
When unset, :func:`dispatch_treatment` returns an empty mapping and
:func:`treatment_scope` is a no-op context manager.

Patterns
--------
Treatment patterns are encoded in :data:`TREATMENT_MATRIX` as tuples of
``(predicate, flags)``. Patterns are evaluated in declaration order; flag
dicts are merged left-to-right, so later patterns can override earlier ones
when both fire.

The current 6 cells (A..F) cover the highest-priority Welle-5l findings:

==  ==========================================  ===========================
ID  Predicate                                    Treatment flags
==  ==========================================  ===========================
A   sigma + CN=6 + maximal-asym + bid+mono      BURNSIDE_FULL + THEOREM_D
B   tris-bidentate-asym                         THEOREM_D
C   |q| >= 4                                    EXTREME_CHARGE_FALLBACK
D   hapto + max hapticity == 5                  CP_PIANO_STOOL
E   bulky alkyl-rich (tBu/iPr/NMe2/PMe3 ...)    ROTAMER_DIVERSITY (K=3)
F   class == "multi-hapto"                      MM_RIGID_DRAG
==  ==========================================  ===========================

The predicates read **only** from the classifier output dict — never from
SMILES strings or refcodes — so the matrix is universal-fundamental
(:doc:`feedback_universal_fundamental_doctrine`).
"""

from __future__ import annotations

import contextlib
import os
from typing import Any, Callable, Dict, List, Mapping, Optional, Tuple

# ---------------------------------------------------------------------------
# Classifier import: prefer Agent Y's dedicated module, fall back to stub
# ---------------------------------------------------------------------------
try:
    # Agent Y's deliverable lands here. If present, it must expose the same
    # function name + return contract.
    from delfin._system_classifier import (  # type: ignore[import-not-found]
        classify_complex_system,
    )
    _CLASSIFIER_SOURCE = "system_classifier"
except Exception:  # pragma: no cover - exercised when Y has not landed yet
    classify_complex_system = None  # type: ignore[assignment]
    _CLASSIFIER_SOURCE = "stub"


# ---------------------------------------------------------------------------
# Stub classifier (used when delfin._system_classifier is unavailable)
# ---------------------------------------------------------------------------
# These element sets are universal — they cover *all* possible TM centers and
# bulky donor patterns. No refcode / element allowlist.
_METAL_ATOMIC_NUMS = frozenset(
    list(range(21, 31))  # Sc..Zn
    + list(range(39, 49))  # Y..Cd
    + list(range(57, 81))  # La..Hg
    + list(range(89, 113))  # Ac..Cn
)

# Bulky alkyl SMARTS-like detection via degree counting (universal).
# A "bulky" carbon = sp3 carbon with degree>=4 (quaternary) bonded to >=3 heavy
# neighbours, OR sp3 nitrogen/phosphorus with degree>=3 (NMe2/PMe3 patterns).


def _stub_classify(mol: Any) -> Dict[str, Any]:
    """Minimal feature extractor used when the dedicated classifier is absent.

    Returned dict has the same keys the TREATMENT_MATRIX predicates expect
    (subset of Agent Y's contract — enough for the 6 cells):

    - ``class``: ``"sigma" | "hapto" | "multi-hapto" | "none"``
    - ``CN``: int, coordination number of primary metal (0 if no metal)
    - ``donor_heterogeneity``: ``"homo" | "bi" | "tri" | "maximal-asym"``
    - ``chelate_pattern``: ``"all-mono" | "mono-plus-bid" | "bis-bid" |
      "tris-bid"``
    - ``q_magnitude``: int, ``|formal charge|`` of primary metal
    - ``hapticity_max``: int, max contiguous metal-bound-C block size
    - ``bulky_flag``: bool

    The stub is *deliberately* conservative — it returns ``"none"`` whenever
    RDKit cannot parse the molecule, ensuring no patterns fire and bit-exact
    OFF behaviour is preserved.
    """
    result: Dict[str, Any] = {
        "class": "none",
        "CN": 0,
        "donor_heterogeneity": "homo",
        "chelate_pattern": "mono-only",
        "q_magnitude": 0,
        "hapticity_max": 0,
        "bulky_flag": False,
    }
    if mol is None:
        return result

    # Find primary metal
    metal_atom = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in _METAL_ATOMIC_NUMS:
            metal_atom = atom
            break
    if metal_atom is None:
        return result

    result["q_magnitude"] = abs(int(metal_atom.GetFormalCharge() or 0))

    # CN = number of heavy neighbours of primary metal
    neighbours = list(metal_atom.GetNeighbors())
    result["CN"] = len(neighbours)

    # donor-element multiset (matches Agent Y's "homo/bi/tri/maximal-asym" labels)
    donor_elements = tuple(sorted(n.GetSymbol() for n in neighbours))
    unique_elems = set(donor_elements)
    n_classes = len(unique_elems)
    if n_classes <= 1:
        result["donor_heterogeneity"] = "homo"
    elif n_classes == 2:
        result["donor_heterogeneity"] = "bi"
    elif n_classes == 3:
        result["donor_heterogeneity"] = "tri"
    else:
        result["donor_heterogeneity"] = "maximal-asym"

    # Hapto detection: contiguous metal-bound C block size
    c_neighbour_idxs = {n.GetIdx() for n in neighbours if n.GetAtomicNum() == 6}
    if c_neighbour_idxs:
        # Find largest contiguous block of metal-bound carbons connected via
        # ring/chain bonds among themselves.
        max_block = 1
        seen: set = set()
        for start in c_neighbour_idxs:
            if start in seen:
                continue
            block = [start]
            stack = [start]
            seen.add(start)
            while stack:
                cur_idx = stack.pop()
                cur = mol.GetAtomWithIdx(cur_idx)
                for nbr in cur.GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in c_neighbour_idxs and ni not in seen:
                        seen.add(ni)
                        block.append(ni)
                        stack.append(ni)
            if len(block) > max_block:
                max_block = len(block)
        result["hapticity_max"] = max_block
    else:
        result["hapticity_max"] = 0

    # Class assignment (universal: based on hapticity_max + metal neighbour set)
    n_hapto_blocks = 0
    if c_neighbour_idxs:
        # Re-walk to count blocks of size >= 2 (eta-2 or higher)
        seen2: set = set()
        for start in c_neighbour_idxs:
            if start in seen2:
                continue
            blk = [start]
            stack = [start]
            seen2.add(start)
            while stack:
                cur = mol.GetAtomWithIdx(stack.pop())
                for nbr in cur.GetNeighbors():
                    ni = nbr.GetIdx()
                    if ni in c_neighbour_idxs and ni not in seen2:
                        seen2.add(ni)
                        blk.append(ni)
                        stack.append(ni)
            if len(blk) >= 2:
                n_hapto_blocks += 1

    if n_hapto_blocks >= 2:
        result["class"] = "multi-hapto"
    elif result["hapticity_max"] >= 2:
        result["class"] = "hapto"
    elif result["CN"] >= 1:
        result["class"] = "sigma"

    # Chelate pattern: count rings containing the metal of size 5 or 6 with two
    # distinct metal-donor positions = bidentate ring. The stub uses RDKit ring
    # info if available.
    bidentate_count = 0
    try:
        ring_info = mol.GetRingInfo()
        metal_idx = metal_atom.GetIdx()
        for ring in ring_info.AtomRings():
            if metal_idx in ring and 4 <= len(ring) <= 7:
                # Count metal-neighbour atoms in this ring (excluding metal)
                in_ring_donors = sum(
                    1 for nidx in ring
                    if nidx != metal_idx
                    and mol.GetAtomWithIdx(nidx).GetIdx() in {n.GetIdx() for n in neighbours}
                )
                if in_ring_donors >= 2:
                    bidentate_count += 1
    except Exception:
        bidentate_count = 0

    # Chelate-pattern labels match Agent Y's contract.
    if bidentate_count >= 3:
        result["chelate_pattern"] = "tris-bid"
    elif bidentate_count == 2:
        result["chelate_pattern"] = "bis-bid"
    elif bidentate_count == 1:
        result["chelate_pattern"] = "mono-plus-bid"
    else:
        result["chelate_pattern"] = "all-mono"

    # Bulky flag: any sp3-C with >=3 carbon neighbours and degree 4
    # (tBu-like quaternary) OR sp3-N/P with >=3 methyl-equivalent neighbours.
    bulky = False
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z == 6 and atom.GetDegree() == 4:
            c_nbr = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
            if c_nbr >= 3:
                bulky = True
                break
        if z in (7, 15) and atom.GetDegree() >= 3:
            c_nbr = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
            if c_nbr >= 2:
                # Heuristic: NMe2-like (2 methyls) or PMe3-like
                methyl_like = sum(
                    1 for n in atom.GetNeighbors()
                    if n.GetAtomicNum() == 6 and n.GetDegree() == 1
                )
                if methyl_like >= 2:
                    bulky = True
                    break
    result["bulky_flag"] = bulky

    return result


def _classify(mol: Any) -> Dict[str, Any]:
    """Internal dispatch helper. Prefers Agent Y's classifier; falls back to
    the local stub. Returns the same dict contract either way."""
    if classify_complex_system is not None:
        try:
            out = classify_complex_system(mol)
            if isinstance(out, Mapping):
                return dict(out)
        except Exception:
            pass
    return _stub_classify(mol)


# ---------------------------------------------------------------------------
# Treatment matrix
# ---------------------------------------------------------------------------
# Predicate signature: ``Callable[[Mapping[str, Any]], bool]``.
# Treatment signature: ``Mapping[str, str]`` — string values so they can be
# round-tripped through ``os.environ``.
PredicateFn = Callable[[Mapping[str, Any]], bool]
TreatmentFlags = Mapping[str, str]


def _pattern_a(c: Mapping[str, Any]) -> bool:
    """sigma + CN6 + maximal-asym donor set with at least one chelate ring.

    Uses Agent Y's classifier contract: ``chelate_pattern`` values are
    ``all-mono / mono-plus-bid / bis-bid / tris-bid / macrocyclic / pincer``.
    "At least one chelate" = anything other than ``all-mono``.
    """
    return (
        c.get("class") == "sigma"
        and int(c.get("CN", 0)) == 6
        and c.get("donor_heterogeneity") == "maximal-asym"
        and c.get("chelate_pattern")
        in ("mono-plus-bid", "bis-bid", "tris-bid", "pincer", "macrocyclic")
    )


def _pattern_b(c: Mapping[str, Any]) -> bool:
    """Tris-bidentate (asymmetric or otherwise) — Theorem-D enumeration."""
    return c.get("chelate_pattern") == "tris-bid"


def _pattern_c(c: Mapping[str, Any]) -> bool:
    """Extreme formal charge magnitude (|q| >= 4)."""
    return int(c.get("q_magnitude", 0)) >= 4


def _pattern_d(c: Mapping[str, Any]) -> bool:
    """Classic Cp piano-stool: hapto class with eta-5 ligand."""
    return c.get("class") == "hapto" and int(c.get("hapticity_max", 0)) == 5


def _pattern_e(c: Mapping[str, Any]) -> bool:
    """Bulky alkyl-rich periphery requires rotamer diversity."""
    return bool(c.get("bulky_flag", False))


def _pattern_f(c: Mapping[str, Any]) -> bool:
    """Multi-hapto cluster (>=2 contiguous metal-bound-C blocks of size >=2)."""
    return c.get("class") == "multi-hapto"


TREATMENT_MATRIX: List[Tuple[PredicateFn, TreatmentFlags]] = [
    # Pattern A — sigma+CN6+maximal-asym+(mono+bid)
    (
        _pattern_a,
        {
            "DELFIN_BURNSIDE_FULL": "1",
            "DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE": "1",
        },
    ),
    # Pattern B — tris-bidentate-asym
    (
        _pattern_b,
        {"DELFIN_5L_T62_THEOREM_D_ASYM_BIDENTATE": "1"},
    ),
    # Pattern C — extreme charge magnitude
    (
        _pattern_c,
        {"DELFIN_EXTREME_CHARGE_FALLBACK": "1"},
    ),
    # Pattern D — hapto + eta-5 → Cp piano-stool
    (
        _pattern_d,
        {"DELFIN_5J_A_CP_PIANO_STOOL": "1"},
    ),
    # Pattern E — bulky alkyl-rich
    (
        _pattern_e,
        {
            "DELFIN_5L_T6_ROTAMER_DIVERSITY": "1",
            "DELFIN_5L_T6_ROTAMER_K": "3",
        },
    ),
    # Pattern F — multi-hapto clusters
    (
        _pattern_f,
        {"DELFIN_5J_G_MM_RIGID_DRAG": "1"},
    ),
]


# ---------------------------------------------------------------------------
# Master gate + public API
# ---------------------------------------------------------------------------
def _matrix_enabled() -> bool:
    """Master gate. Returns True only when DELFIN_5M_TREATMENT_MATRIX is set
    to a truthy value. Default OFF guarantees bit-exact reproducibility."""
    val = os.environ.get("DELFIN_5M_TREATMENT_MATRIX", "").strip().lower()
    return val in {"1", "true", "yes", "on"}


def dispatch_treatment(mol: Any) -> Dict[str, str]:
    """Return env-flag overrides for ``mol``.

    When the master gate is OFF the result is an empty dict — this preserves
    bit-exact default behaviour. When ON, the function classifies the
    molecule and merges flag dicts from every matching pattern (in matrix
    declaration order; later patterns override earlier ones).
    """
    if not _matrix_enabled():
        return {}
    if mol is None:
        return {}
    features = _classify(mol)
    flags: Dict[str, str] = {}
    for predicate, treatment in TREATMENT_MATRIX:
        try:
            if predicate(features):
                flags.update(treatment)
        except Exception:
            # A faulty predicate must never break the conversion pipeline.
            continue
    return flags


@contextlib.contextmanager
def treatment_scope(mol: Any):
    """Context manager: apply per-molecule env-flag overrides, restore on exit.

    Usage::

        with treatment_scope(mol):
            # pipeline calls here see the per-molecule env flags
            ...

    When the master gate is OFF the context is a no-op (no env reads,
    no env writes) and the runtime overhead is a single dict lookup.
    """
    overrides = dispatch_treatment(mol)
    if not overrides:
        # Fast no-op path — preserves bit-exact OFF behaviour.
        yield {}
        return

    previous: Dict[str, Optional[str]] = {}
    try:
        for key, value in overrides.items():
            previous[key] = os.environ.get(key)
            os.environ[key] = value
        yield dict(overrides)
    finally:
        for key, prev_val in previous.items():
            if prev_val is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = prev_val


# ---------------------------------------------------------------------------
# Introspection helpers (used by tests and the dashboard)
# ---------------------------------------------------------------------------
def matrix_summary() -> List[Dict[str, Any]]:
    """Return a JSON-serialisable summary of all matrix patterns for the
    dashboard / debugging tools."""
    return [
        {
            "predicate_name": predicate.__name__,
            "predicate_doc": (predicate.__doc__ or "").strip(),
            "flags": dict(flags),
        }
        for predicate, flags in TREATMENT_MATRIX
    ]


def classifier_source() -> str:
    """Return ``"system_classifier"`` if Agent Y's module is loaded, else
    ``"stub"``."""
    return _CLASSIFIER_SOURCE
