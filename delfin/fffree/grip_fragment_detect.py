"""GRIP — Fragment Detector (Phase 2, v1).

Enumerates the bond / angle / improper fragment instances in a molecule and
emits a deterministic, sorted list of :mod:`grip_loss_terms` terms whose
``(mu, sigma)`` come from the :mod:`grip_mogul_lookup` library.

Rules (SPEC §3.4 + §11):

* skip every fragment that touches a ``frozen_atom`` (default: metal + its
  donor atoms — coordination sphere is M-D-rigid, handled by other modules)
* (Heal 1, 2026-06-01, mddir): when ``donors`` is provided, also protect
  the donor's heavy first-shell neighbours from being repositioned as the
  *centre* of an angle or improper.  The Mahalanobis pull on those terms
  rotates the M-D-X axis around the donor, which keeps the M-D distance
  exactly but flips the lone-pair direction — exactly what mddir flags
  (``feedback_phase5_grip_smoke500_verdict``).  Bond terms involving the
  shell-1 atoms are kept (so the X-X' bond lengths still pull toward
  Mogul), only the central-atom rotations are dropped.
* (Heal 1b, 2026-06-01, amine-H regression): Heal-1 was too strict — it
  dropped EVERY angle/improper centred on a donor's shell-1 atom, which
  included C-N-H angles and let amine H drift to unrealistic orientations
  (smoke-50 surfaced +100% amine_h_pct_files_viol).  Heal-1b refines
  the rule: skip the shell-1-centred term ONLY when none of the non-centre
  participants is hydrogen.  Concretely: C-N-H, H-N-H angles centred on
  the shell-1 C are kept; C-N-C, C-N-Cα heavy-only triples still drop.
  Same for impropers (kept iff at least one neighbour is hydrogen).  This
  preserves the M-D-X heavy-rotation guard while restoring the H
  orientation priors.
* (Option B, 2026-06-01, adaptive coverage-aware protection):
  Heal-1/1b's unconditional drop of heavy-only shell-1-centred angles
  helped mddir but stripped CCDC-Mahalanobis pressure off the funcgrp
  internals, surfacing +126% F23_funcgrp_geom regression in the
  Heal-Stack smoke-500.  Option B refines the rule once more: when a
  GripLibrary is available, the heavy-only shell-1 skip activates ONLY
  if the library coverage for the SPECIFIC fragment class is thin
  (lookup returns None with ``min_n=adaptive_min_n``, default 5).  When
  the coverage is good (lookup hits with n ≥ 5), the term is KEPT so the
  funcgrp prior keeps pulling.  When the lookup is unavailable
  (``library=None``) the legacy Heal-1/1b behaviour is reproduced
  byte-exactly so existing tests stay green.  Universal: no
  SMILES-specific gating — every decision is driven by the
  graph-only key and the library's CCDC-Mogul statistics.
* (Hapto-class protection, 2026-06-02, voll-pool regression fix):
  The voll-pool 11458 SMILES Heal-Stack run surfaced 4 severe hapto-
  class regressions (hapto_geom_per_frame +182 %, hapto_geom_pct +104 %,
  haptounit_pct +44 %, hapto_count_pct +38 %).  Mechanism: when a metal
  is η-bound to a π-system (ferrocene Cp, η⁶-arene, η³-allyl, ...) the
  "donors" are π-carbons whose bond/angle priors live in a fundamentally
  different distribution from the CCDC Mogul library (which is dominated
  by discrete σ-donor chemistry).  GRIP's Mahalanobis pull then forces
  the hapto-π atoms toward standard sp2-aromatic geometries and breaks
  the piano-stool placement we worked so hard to construct.  The fix:
  callers pass the set of hapto-π atom indices via ``hapto_atoms`` and
  detect_fragments skips EVERY bond / angle / improper term that
  touches any of them.  Env-flag ``DELFIN_FFFREE_GRIP_SKIP_HAPTO`` (read
  per-call, default ``"1"``) lets the operator disable the skip for
  forensic A/B without a code edit; when ``hapto_atoms`` is empty or
  ``None`` the function is byte-identical to the pre-fix path.
  Universal: hapto atom-set is graph-derived (ring membership +
  metal-edge contiguity), never SMILES-pattern-matched.
* deterministic iteration order — bonds/angles/impropers are sorted by
  atom-index tuple lex before terms are constructed
* fragments without a library match (n < 5 across the whole fallback chain)
  are silently dropped — they have no statistical prior
* float64 throughout
"""
from __future__ import annotations

import os
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np

from .grip_loss_terms import (
    AngleTerm,
    BondTerm,
    ImproperTerm,
    sparse_downweight,
    default_weights,
)
from .grip_mogul_lookup import (
    GripLibrary,
    get_default_library,
    GRIP_LOOKUP_MIN_N,
)

__all__ = [
    "detect_fragments",
    "FragmentDetectionResult",
    "ANGLE_TO_METAL_ENV",
    "angle_to_metal_active",
    "build_mdx_angle_terms",
    "DMD_POLY_ENV",
    "dmd_polyhedron_active",
    "build_dmd_polyhedron_terms",
]


# ---------------------------------------------------------------------------
# Angle-to-Metal extension (2026-06-05, User direction).
#
# The legacy detector skips EVERY angle whose any atom is in the M-D frozen
# set (metal + donors).  That was deliberate -- M-D bond lengths and the
# M-position are rigid by construction; CCDC priors for X-Y-Z angles
# around the metal are dominated by sigma-donor chemistry.
#
# But this leaves TWO classes of geometric pressure on the table:
#
#   (1) M-D-X angles (centre = donor, neighbours = metal + ligand-X):
#       these decide WHERE the ligand-X atoms sit around the donor (the
#       angle the donor cone makes with the M-D axis).  In the legacy
#       detector this angle gets ZERO pressure.  The donor stays frozen
#       and the metal stays frozen, so the gradient acts ONLY on X --
#       which is exactly what the user asked for: optimise the position
#       of atoms around the donor without touching the M-D bond length
#       or the M-D-D' polyhedron arrangement.
#
#   (2) D-M-D' angles (centre = metal):  these encode the coordination
#       polyhedron.  Both endpoints (donors) and the centre (metal) are
#       frozen so gradient is identically zero.  Useful only as a
#       severity diagnostic; not actionable for L-BFGS without unfreezing
#       the polyhedron -- which the M-D constraint would then roll back.
#
# This module implements class (1) (and exposes class (2) for diagnostic /
# severity gating).  Default OFF -- when the env-flag is unset the angle
# stack is byte-identical with the legacy path (no new terms emitted,
# legacy frozen-skip preserved).
#
# Universal: angle classification uses only graph topology + the (metal,
# donors) tuple; no SMILES patterns, no class-specific weights.
# ---------------------------------------------------------------------------
ANGLE_TO_METAL_ENV = "DELFIN_FFFREE_GRIP_ANGLE_TO_METAL"
DMD_POLY_ENV = "DELFIN_FFFREE_GRIP_DMD_POLYHEDRON"
ANGLE_TO_METAL_INCLUDE_H_ENV = "DELFIN_FFFREE_GRIP_ANGLE_TO_METAL_INCLUDE_H"

# ---------------------------------------------------------------------------
# Mission D2 (2026-06-05): per-class weight env-overrides.
#
# The legacy default-weights table is (bond=1.0, angle=0.5, improper=2.0,
# torsion=0.2).  These overrides bump bond/angle weights so the Mahalanobis
# pull on the most-actionable internal geometry classes outranks the rest.
# Default-OFF byte-identical: when env unset, the legacy table is returned
# unchanged.
# ---------------------------------------------------------------------------
WEIGHT_BOND_ENV = "DELFIN_FFFREE_GRIP_WEIGHT_BOND"
WEIGHT_ANGLE_ENV = "DELFIN_FFFREE_GRIP_WEIGHT_ANGLE"
WEIGHT_IMPROPER_ENV = "DELFIN_FFFREE_GRIP_WEIGHT_IMPROPER"
WEIGHT_TORSION_ENV = "DELFIN_FFFREE_GRIP_WEIGHT_TORSION"


def _resolve_weight_env(env_name: str, default: float) -> float:
    """Read a float-valued env-var; return ``default`` when unset / invalid."""
    raw = os.environ.get(env_name, "").strip()
    if not raw:
        return float(default)
    try:
        v = float(raw)
        if np.isfinite(v) and v >= 0:
            return v
    except (TypeError, ValueError):
        pass
    return float(default)


def _apply_weight_env_overrides(w: dict) -> dict:
    """Apply Mission D2 env-flag weight overrides to ``w`` (mutates + returns).

    Default-OFF byte-identical: when no env-flags are set, ``w`` is unchanged.
    """
    w["bond"] = _resolve_weight_env(WEIGHT_BOND_ENV, w.get("bond", 1.0))
    w["angle"] = _resolve_weight_env(WEIGHT_ANGLE_ENV, w.get("angle", 0.5))
    w["improper"] = _resolve_weight_env(WEIGHT_IMPROPER_ENV, w.get("improper", 2.0))
    w["torsion"] = _resolve_weight_env(WEIGHT_TORSION_ENV, w.get("torsion", 0.2))
    return w


def _angle_to_metal_include_h() -> bool:
    """``True`` iff M-D-H terms should be emitted as well (default OFF).

    The H atoms are 1.0-1.1 Å from the donor so a small angular pull on
    them translates to a large *relative* bond-length change.  This
    triggers the topology safety rollback far more often than it helps,
    so the include-H path is gated by its own env-flag.
    """
    raw = os.environ.get(ANGLE_TO_METAL_INCLUDE_H_ENV, "").strip().lower()
    return raw in ("1", "true", "yes", "on")

# Per-hybridisation VSEPR fallback when the library has no TM-aware angle
# entry for the (M, donor, X) triple.  Sigma is generous (10°) so the soft
# pull does not fight the existing CCDC X-D-X' pressure -- it only adds
# direction to the otherwise-zero gradient on X.
#
# The hybridisation here is the donor's UNDERLYING hyb without the M bond
# (we count heavy non-metal neighbours + lone pairs to derive sp/sp2/sp3
# correctly).  This avoids RDKit's tendency to mark donor C/N/O as sp3d
# just because the metal adds a 5th neighbour.
_VSEPR_FALLBACK = {
    "sp":   (180.0, 18.0),
    "sp2":  (120.0, 15.0),
    "sp3":  (109.5, 15.0),
    "sp3d": (90.0,  20.0),
    "sp3d2":(90.0,  20.0),
}
_VSEPR_DEFAULT = (109.5, 18.0)


def _donor_underlying_hyb(mol, donor_idx: int, metal_idx: int) -> str:
    """Infer the donor's underlying hybridisation by counting NON-metal
    heavy neighbours.

    Rule (graph-only, universal):
      * deg_no_M = number of bonded heavy atoms minus the metal
      * sp3 for tetrahedral (deg_no_M >= 3 if z in {N,O,S,P,As,Se,C})
      * sp2 for trigonal-planar (deg_no_M == 2 with sp2 RDKit hint OR
        donor in aromatic ring)
      * sp for linear (deg_no_M == 1 and RDKit says sp, e.g. C of -C#N-M)
      * default sp3
    """
    try:
        atom = mol.GetAtomWithIdx(int(donor_idx))
    except Exception:
        return "sp3"
    # Underlying RDKit hyb (may be sp3d due to metal bond) — use as a hint.
    try:
        rdkit_hyb = atom.GetHybridization().name.upper()
    except Exception:
        rdkit_hyb = ""
    deg_no_M = 0
    for nb in atom.GetNeighbors():
        try:
            if int(nb.GetIdx()) == int(metal_idx):
                continue
            deg_no_M += 1
        except Exception:
            pass
    # Aromatic ring centre -> sp2.
    try:
        if atom.GetIsAromatic():
            return "sp2"
    except Exception:
        pass
    z = atom.GetSymbol()
    if z == "C":
        # Linear C donor (terminal -C#N coordinated through C is rare; the
        # cyanide N-coordinated form is more common, this is the C-end).
        if rdkit_hyb in ("SP",) or deg_no_M == 1:
            return "sp"
        if rdkit_hyb in ("SP2",) or deg_no_M == 2:
            return "sp2"
        return "sp3"
    if z in ("N", "P", "As"):
        if rdkit_hyb in ("SP",) or deg_no_M == 1:
            return "sp"
        if rdkit_hyb in ("SP2",) or deg_no_M == 2:
            return "sp2"
        return "sp3"
    if z in ("O", "S", "Se"):
        # O/S donors are almost always sp3 (water, alkoxide, thiolate).
        # An aromatic O (furan-O) would have been caught by the aromatic
        # branch above.
        return "sp3"
    # Halides + others: treat as sp3 (single lone pair acceptor).
    return "sp3"

# Weight multiplier applied to the M-D-X terms.  The default angle weight
# (0.5 from grip_loss_terms.default_weights) already balances bond/angle
# pressure; the M-D-X terms are a NEW pressure that did not exist in the
# legacy operator, so we down-weight them by default (0.1) so the
# Mahalanobis pull on X is gentle and does not over-stretch X-Y bonds
# beyond the topology safety multiplier (180 % stretch when the layer is
# on, see grip_polish._topo_mult_eff override).
# Override via env var ``DELFIN_FFFREE_GRIP_ANGLE_TO_METAL_WEIGHT``.
# Calibration: smoke 50 on 2792332-VOLLPOOL @ wt=0.1 gives delta_RMSE
# (vs OFF) = -0.60° M-D-X with 5/46 polishes accepted; @ wt=0.25 gives
# -0.14° / 4 accepted; @ wt=0.5 gives -0.07° / 3 accepted.  Lower weight
# = smaller per-polish pull but higher accept rate -> larger aggregate
# improvement.  0.1 is the empirical sweet spot.
_ANGLE_TO_METAL_WEIGHT_ENV = "DELFIN_FFFREE_GRIP_ANGLE_TO_METAL_WEIGHT"
_DEFAULT_ANGLE_TO_METAL_WEIGHT = 0.1


def angle_to_metal_active() -> bool:
    """``True`` iff the M-D-X angle layer is enabled.

    Pure env-flag check; default OFF preserves byte-identity with the
    legacy detector.  Read-per-call so subprocess pools inherit the
    parent's setting.
    """
    raw = os.environ.get(ANGLE_TO_METAL_ENV, "").strip().lower()
    return raw in ("1", "true", "yes", "on")


def dmd_polyhedron_active() -> bool:
    """``True`` iff the D-M-D polyhedron severity terms are emitted.

    Default OFF.  D-M-D terms have zero gradient (M and donors are both
    frozen) so they only add to the accept-if-better severity metric --
    enabling this flag without ``ANGLE_TO_METAL_ENV`` makes the gate
    polyhedron-aware without changing the L-BFGS trajectory.
    """
    raw = os.environ.get(DMD_POLY_ENV, "").strip().lower()
    return raw in ("1", "true", "yes", "on")


def _resolve_angle_to_metal_weight() -> float:
    raw = os.environ.get(_ANGLE_TO_METAL_WEIGHT_ENV, "").strip()
    if not raw:
        return float(_DEFAULT_ANGLE_TO_METAL_WEIGHT)
    try:
        v = float(raw)
        if np.isfinite(v) and v >= 0:
            return float(v)
    except (TypeError, ValueError):
        pass
    return float(_DEFAULT_ANGLE_TO_METAL_WEIGHT)


def build_mdx_angle_terms(
    mol,
    metal: int,
    donors: Sequence[int],
    hapto_atoms: Optional[Iterable[int]] = None,
    library: Optional[GripLibrary] = None,
    min_n: int = GRIP_LOOKUP_MIN_N,
    *,
    include_hydrogens: bool = False,
) -> List[AngleTerm]:
    """Build M-D-X angle terms (centre = donor, neighbours = metal + X).

    For each (donor d, non-metal/donor neighbour x of d) emit an
    :class:`AngleTerm` with ``a=metal, b=d, c=x``.  Mu/sigma come from
    the library (TM-aware fallback wins when v4 / pair-resolved data is
    available), otherwise from the VSEPR-template fallback keyed by the
    donor's hybridisation.

    Atoms in ``hapto_atoms`` (the π-system carbons of η³-η⁸ coordination)
    are SKIPPED -- the M-π-X angle is the piano-stool placement axis and
    has a fundamentally different distribution from σ M-D-X.

    Deterministic: sorted (donor, x) iteration, lex output order.

    Returns
    -------
    list of AngleTerm
        Possibly empty (e.g. when every donor is hapto or no neighbours).
    """
    out: List[AngleTerm] = []
    if metal is None or metal < 0:
        return out
    try:
        mol_metal = mol.GetAtomWithIdx(int(metal))  # noqa: F841 -- existence probe
    except Exception:
        return out
    donor_seq = sorted(int(d) for d in donors)
    if not donor_seq:
        return out
    hapto_set: Set[int] = set(int(i) for i in (hapto_atoms or ()))
    frozen_or_donor: Set[int] = {int(metal), *donor_seq}
    weight_scale = _resolve_angle_to_metal_weight()

    # Cache the metal symbol for library lookups.
    try:
        zM = str(mol.GetAtomWithIdx(int(metal)).GetSymbol())
    except Exception:
        zM = "*"

    lib = library  # may be None -- VSEPR fallback covers that case

    for d in donor_seq:
        if d in hapto_set:
            # η-carbons: M-C-C angles within the π-system are placement
            # geometry, not σ-donor chemistry.
            continue
        try:
            atom_d = mol.GetAtomWithIdx(int(d))
        except Exception:
            continue
        zD = atom_d.GetSymbol()
        # Underlying donor hyb (without the metal bond) -- used both for
        # the library lookup AND for the VSEPR fallback.  Avoids the
        # rdkit-sp3d-from-metal-coordination defect that maps tetrahedral
        # donor C/N to a 90° prior instead of 109.5°.
        hyb_d = _donor_underlying_hyb(mol, int(d), int(metal))
        vsepr_mu, vsepr_sigma = _VSEPR_FALLBACK.get(hyb_d, _VSEPR_DEFAULT)

        # Gather X neighbours (any neighbour of d that is not the metal
        # and not another donor in the same coordination sphere).
        neighbours = sorted(int(n.GetIdx()) for n in atom_d.GetNeighbors())
        for x in neighbours:
            if x in frozen_or_donor:
                # M-D-D' = D-M-D' relative to M and is handled by the
                # polyhedron layer; M-D-M doesn't exist in mono-metal
                # complexes; skip.
                continue
            if x in hapto_set:
                continue
            try:
                atom_x = mol.GetAtomWithIdx(int(x))
            except Exception:
                continue
            zX = atom_x.GetSymbol()
            hyb_x = _hyb_str(atom_x)
            # Skip hydrogens by default: H atoms are 1.0-1.1 Å from the
            # donor, so a small angular pull translates into a large
            # relative bond-length change.  This caused the topology
            # constraint to roll back the polish on every smoke-50 frame
            # in the 2792332 pool (C-H stretched > 1.5×).  Heavy-atom X
            # carries the bulk of the donor cone geometry; H is then
            # carried along by the existing C-H bond term.  Caller can
            # opt back into H via include_hydrogens=True.
            if not include_hydrogens and zX == "H":
                continue

            # Try the library first (TM-aware fallback is active per the
            # env var DELFIN_FFFREE_GRIP_TM_LOOKUP_FALLBACK; when v4 is
            # loaded this hits the triple_angle pair-resolved table).
            mu, sigma, n_lib = vsepr_mu, vsepr_sigma, 0
            source = "vsepr"
            if lib is not None:
                try:
                    hit = lib.lookup_angle(
                        zM, zD, hyb_d, zX,
                        ring_size_min=-1, in_aromatic=False,
                        hyb1="*", hyb3=hyb_x,
                        min_n=min_n,
                    )
                except Exception:
                    hit = None
                if hit is not None:
                    mu, sigma, n_lib = float(hit[0]), float(hit[1]), int(hit[2])
                    source = "library"
            # Sanity: angle priors must be in (0°, 180°) with sigma > 0.
            if not (0.0 < mu < 180.0 and sigma > 0.0):
                continue
            # Apply the weight scaling.  Library-sourced terms also get
            # the legacy sparse_downweight so under-sampled bins stay
            # gentle; VSEPR-sourced terms get a flat down-weight (n=0
            # would zero them via sparse_downweight, which is the
            # opposite of what we want for the VSEPR fallback).
            if source == "library":
                w = sparse_downweight(0.5, n_lib, n_floor=min_n) * weight_scale
            else:
                w = 0.5 * weight_scale
            if w <= 0.0:
                continue
            out.append(AngleTerm(
                a=int(metal), b=int(d), c=int(x),
                mu=mu, sigma=sigma, weight=float(w),
            ))
    # Lex-deterministic ordering (atom_indices triple).
    out.sort(key=lambda t: t.atom_indices)
    return out


def build_dmd_polyhedron_terms(
    mol,                        # noqa: ARG001 -- API parity (graph not needed)
    metal: int,
    donors: Sequence[int],
    P: np.ndarray,
    sigma_dmd: float = 5.0,
    weight: float = 0.25,
) -> List[AngleTerm]:
    """Build D-M-D' angle terms keyed off the OBSERVED initial geometry.

    The mu of each (Di, M, Dj) angle is the angle measured in the input
    coordinates ``P`` -- this freezes the polyhedron arrangement we
    constructed.  Sigma is a soft tolerance (default 5°) so the term
    contributes to severity but does not pull the (frozen) donors.

    Returns terms whose gradient acts only on the metal + donor positions,
    all of which the polish freezes -- so these terms are ZERO-GRADIENT
    by construction.  They only matter for the accept-if-better severity
    score: a polished geometry that distorts the polyhedron loses ground
    to one that preserves it.

    Pure-functional; deterministic ordering (sorted donor index pairs).
    """
    out: List[AngleTerm] = []
    if metal is None or metal < 0:
        return out
    P = np.asarray(P, dtype=np.float64)
    ds = sorted(int(d) for d in donors)
    if len(ds) < 2:
        return out
    rM = P[int(metal)]
    for i in range(len(ds)):
        for j in range(i + 1, len(ds)):
            di, dj = ds[i], ds[j]
            try:
                u = P[di] - rM
                v = P[dj] - rM
                nu = float(np.linalg.norm(u))
                nv = float(np.linalg.norm(v))
                if nu < 1e-9 or nv < 1e-9:
                    continue
                cos_t = float(np.dot(u, v) / (nu * nv))
                cos_t = max(-1.0, min(1.0, cos_t))
                import math as _math
                theta_deg = _math.degrees(_math.acos(cos_t))
            except Exception:
                continue
            if not (0.0 < theta_deg < 180.0):
                continue
            out.append(AngleTerm(
                a=int(di), b=int(metal), c=int(dj),
                mu=theta_deg, sigma=float(sigma_dmd), weight=float(weight),
            ))
    out.sort(key=lambda t: t.atom_indices)
    return out

# RDKit is imported lazily so the module also works in test paths that monkey
# patch / mock RDKit (and so the module imports even on installations without
# RDKit). detect_fragments raises a clear error if RDKit is unavailable.
def _import_rdkit():
    try:
        from rdkit import Chem  # noqa: F401
        from rdkit.Chem import rdchem  # noqa: F401
    except Exception as exc:
        raise ImportError(
            "RDKit is required by grip_fragment_detect.detect_fragments; "
            "install rdkit-pypi or activate the project's micromamba env."
        ) from exc
    return Chem


# Map RDKit hybridisation enum to the library's string convention.
_HYB_MAP = {
    "S": "sp",         # special-case unusual codes
    "SP": "sp",
    "SP2": "sp2",
    "SP3": "sp3",
    "SP3D": "sp3d",
    "SP3D2": "sp3d2",
    "UNSPECIFIED": "*",
    "OTHER": "*",
    "UNKNOWN": "*",
}


def _hyb_str(atom) -> str:
    """Return the library string for an RDKit atom's hybridisation."""
    try:
        name = atom.GetHybridization().name
    except Exception:
        return "*"
    return _HYB_MAP.get(name.upper(), "*")


def _ring_size_min(mol, atom_idx: int) -> int:
    """Smallest ring containing ``atom_idx``; -1 if not in any ring."""
    ri = mol.GetRingInfo()
    if not ri.NumAtomRings(atom_idx):
        return -1
    rings = ri.AtomRings()
    sizes = [len(r) for r in rings if atom_idx in r]
    if not sizes:
        return -1
    return int(min(sizes))


def _bond_ring_min(mol, a: int, b: int) -> int:
    """Smallest ring containing both atoms; -1 if none."""
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    sizes = [len(r) for r in rings if (a in r and b in r)]
    if not sizes:
        return -1
    return int(min(sizes))


def _atom_is_aromatic(atom) -> bool:
    try:
        return bool(atom.GetIsAromatic())
    except Exception:
        return False


class FragmentDetectionResult:
    """Container with detection diagnostics for tests / debugging.

    Holds the produced terms plus per-class match/miss counters so callers
    can reason about library coverage in equal-n smokes.
    """

    __slots__ = (
        "bond_terms", "angle_terms", "improper_terms",
        "n_bond_candidates", "n_bond_matched",
        "n_angle_candidates", "n_angle_matched",
        "n_improper_candidates", "n_improper_matched",
    )

    def __init__(self):
        self.bond_terms: List[BondTerm] = []
        self.angle_terms: List[AngleTerm] = []
        self.improper_terms: List[ImproperTerm] = []
        self.n_bond_candidates = 0
        self.n_bond_matched = 0
        self.n_angle_candidates = 0
        self.n_angle_matched = 0
        self.n_improper_candidates = 0
        self.n_improper_matched = 0

    def all_terms(self) -> List:
        """Return all terms in a single deterministic-order list."""
        out = list(self.bond_terms) + list(self.angle_terms) + list(self.improper_terms)
        out.sort(key=lambda t: tuple(t.atom_indices))
        return out


def detect_fragments(
    mol,
    P: np.ndarray,
    frozen_atoms: Optional[Iterable[int]] = None,
    *,
    donors: Optional[Iterable[int]] = None,
    hapto_atoms: Optional[Iterable[int]] = None,
    library: Optional[GripLibrary] = None,
    min_n: int = GRIP_LOOKUP_MIN_N,
    weights: Optional[dict] = None,
    return_result: bool = False,
    adaptive_shell1: bool = True,
    adaptive_min_n: int = GRIP_LOOKUP_MIN_N,
    metal: Optional[int] = None,
):
    """Enumerate bond/angle/improper fragments in ``mol`` and build GRIP terms.

    Parameters
    ----------
    mol : RDKit ``Mol`` / ``RWMol``
        Molecule with bonds, atoms, and hybridisation populated.
    P : ndarray (n_atoms, 3)
        Atom coordinates.  Not used at detection time (kept for API symmetry
        with the spec) — values are observed at :meth:`value_and_grad` time.
    frozen_atoms : iterable of int, optional
        Atom indices to treat as rigid: any fragment that touches at least
        one of them is skipped (SPEC §3.4 — M-D sphere is handled elsewhere).
        Typical caller composition: ``{metal, *donors}``.
    donors : iterable of int, optional
        Donor atom indices.  When supplied (Heal-1, mddir-fix 2026-06-01),
        the donor's first-shell heavy neighbours are additionally protected
        from being the *central* atom of an angle or improper term.  This
        prevents the Mahalanobis pull on those terms from rotating the
        M-D-X axis around the (frozen) donor — a subtle leak that lets
        mddir flag the polished geometry even when the M-D distance is
        exactly preserved.  Bonds and end-position angle terms involving
        the shell-1 atoms are still emitted so the X-X' bond/angle priors
        keep working.  If ``donors`` is ``None`` the legacy behaviour is
        recovered byte-exactly.
    hapto_atoms : iterable of int, optional
        Atom indices that belong to a hapto-coordinated π-system
        (η³-allyl, η⁴-diene, η⁵-Cp, η⁶-arene, η⁷-cycloheptatrienyl,
        η⁸-COT).  Every bond, angle, and improper term that touches at
        least one of these atoms is skipped — the CCDC Mogul priors are
        dominated by non-hapto chemistry and pull the η-coordinated
        carbons toward standard sp²/sp³ donor geometries, which
        destroys the piano-stool placement that fffree constructs and
        produced the +182 % hapto_geom regression in the 2026-06-02
        voll-pool.  When ``hapto_atoms`` is ``None`` / empty OR the
        environment variable ``DELFIN_FFFREE_GRIP_SKIP_HAPTO`` is
        ``"0"`` / ``"false"`` / ``"no"``, the function is byte-identical
        to the pre-fix path (Option-B + Heal-1/1b + frozen behaviour
        unchanged).  Universal: the caller is responsible for deriving
        the set from the molecule graph (see :mod:`grip_polish` for the
        canonical detector).
    library : GripLibrary, optional
        Pre-loaded library instance; uses :func:`get_default_library` if None.
    min_n : int
        Minimum sample size for a library entry to be trusted (SPEC §2.4).
    weights : dict, optional
        Per-class weight overrides.  Defaults to :func:`default_weights`.
    return_result : bool
        If True, returns a :class:`FragmentDetectionResult` with diagnostics;
        otherwise returns the flat sorted list of terms.
    adaptive_shell1 : bool
        Option-B switch (2026-06-01).  When True (default) AND a library
        is available AND ``donors`` is provided, the heavy-only shell-1
        skip (Heal-1) is performed ONLY for fragment classes whose
        library coverage is thin (the angle / improper lookup returns
        ``None`` or ``n < adaptive_min_n``).  When the coverage is good
        the term is kept so the CCDC-Mogul Mahalanobis pressure on the
        funcgrp internals is preserved.  Set False to recover the pure
        Heal-1/1b behaviour.  When ``library`` is None this flag is
        effectively a no-op (legacy fallback path is taken).
    adaptive_min_n : int
        Sample-size threshold used by the adaptive coverage check
        (defaults to :data:`GRIP_LOOKUP_MIN_N`, i.e. 5).  A library
        lookup that hits with ``n >= adaptive_min_n`` is considered
        "good coverage" and the shell-1 term is kept.

    Returns
    -------
    list of terms (default) or :class:`FragmentDetectionResult`
    """
    _import_rdkit()

    frozen: Set[int] = set(int(i) for i in (frozen_atoms or ()))
    donor_set: Set[int] = set(int(i) for i in (donors or ()))
    # Heal-1: union the donor set into the frozen membership test so the
    # caller does not need to remember to put donors into BOTH frozen and
    # donors.  This is the lex of "donors are M-D-rigid" already enforced
    # in grip_polish (frozen={metal, *donors}); duplicating it here makes
    # the function safe to call standalone with donors only.
    frozen_or_donor: Set[int] = frozen | donor_set

    # Hapto-class protection (2026-06-02): build the active hapto set.
    # When the env flag DELFIN_FFFREE_GRIP_SKIP_HAPTO is "0"/"false"/"no"
    # the operator wants to disable the skip (forensic A/B comparison),
    # so the active set is forced empty.  When unset (default) the skip
    # is ON.  Reading the env per-call ensures fork-based parallel pools
    # honour overrides written into the child environment.
    _skip_hapto_env = os.environ.get(
        "DELFIN_FFFREE_GRIP_SKIP_HAPTO", "1"
    ).strip().lower()
    _skip_hapto_on = _skip_hapto_env not in ("0", "false", "no", "")
    if _skip_hapto_on:
        hapto_set: Set[int] = set(int(i) for i in (hapto_atoms or ()))
    else:
        hapto_set = set()

    # Heal-1: build the donor first-shell protection set.  These atoms must
    # not become the *central* atom of an angle or improper, because such
    # terms tug them off the M-D ray (their motion direction is tangent to
    # the M-D sphere — the donor stays frozen but the X position rotates
    # around it, which is what mddir's sp3_anti and donor_linearized modes
    # detect).  Bond terms touching these atoms are KEPT so the X-X' bond
    # priors still pull toward Mogul medians (no F3_bond regression).
    donor_shell1: Set[int] = set()
    if donor_set:
        for di in donor_set:
            try:
                atom = mol.GetAtomWithIdx(int(di))
            except Exception:
                continue
            for nb in atom.GetNeighbors():
                try:
                    if nb.GetSymbol() == "H":
                        continue
                    j = int(nb.GetIdx())
                except Exception:
                    continue
                if j in frozen_or_donor:
                    continue
                donor_shell1.add(j)

    w = dict(default_weights())
    if weights:
        w.update(weights)
    # Mission D2 (2026-06-05): env-overrides for bond/angle/improper/torsion
    # global weights.  Default-OFF byte-identical (env unset -> no change).
    w = _apply_weight_env_overrides(w)

    lib = library if library is not None else get_default_library()
    P = np.asarray(P, dtype=np.float64)

    result = FragmentDetectionResult()

    # ---- 1) Bonds ----------------------------------------------------------
    # Deterministic order: sort by (min(a,b), max(a,b)).
    bond_pairs: List[Tuple[int, int]] = []
    for bond in mol.GetBonds():
        a = int(bond.GetBeginAtomIdx())
        b = int(bond.GetEndAtomIdx())
        if a == b:
            continue
        if a > b:
            a, b = b, a
        bond_pairs.append((a, b))
    bond_pairs = sorted(set(bond_pairs))

    for (a, b) in bond_pairs:
        # Heal-1: skip if either endpoint is metal or donor (frozen_or_donor).
        # Shell-1 atoms are NOT excluded here — we want their X-X' bond
        # length prior to keep working (preserves F3_bond/funcgrp gains).
        if a in frozen_or_donor or b in frozen_or_donor:
            continue
        # Hapto-class protection: skip every bond that touches a hapto-π
        # atom — the CCDC bond-length prior pulls the η-carbons toward
        # standard sp2-aromatic geometries that are inconsistent with the
        # piano-stool placement (voll-pool +182 % hapto_geom regression).
        if hapto_set and (a in hapto_set or b in hapto_set):
            continue
        result.n_bond_candidates += 1
        atom_a = mol.GetAtomWithIdx(a)
        atom_b = mol.GetAtomWithIdx(b)
        z1 = atom_a.GetSymbol()
        z2 = atom_b.GetSymbol()
        hyb1 = _hyb_str(atom_a)
        hyb2 = _hyb_str(atom_b)
        rmin = _bond_ring_min(mol, a, b)
        aromatic = _atom_is_aromatic(atom_a) and _atom_is_aromatic(atom_b)

        hit = lib.lookup_bond(z1, hyb1, z2, hyb2, rmin, aromatic, min_n=min_n)
        if hit is None:
            continue
        mu, sigma, n = hit
        weight = sparse_downweight(w["bond"], n, n_floor=min_n)
        result.bond_terms.append(BondTerm(a=a, b=b, mu=mu, sigma=sigma, weight=weight))
        result.n_bond_matched += 1

    # ---- 2) Angles ---------------------------------------------------------
    # For each atom b with degree >= 2, every ordered (a,c) pair with a<c.
    angle_triples: List[Tuple[int, int, int]] = []
    for atom in mol.GetAtoms():
        b = int(atom.GetIdx())
        nbrs = sorted(int(n.GetIdx()) for n in atom.GetNeighbors())
        if len(nbrs) < 2:
            continue
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                a, c = nbrs[i], nbrs[j]
                angle_triples.append((a, b, c))
    angle_triples = sorted(set(angle_triples))

    for (a, b, c) in angle_triples:
        # Heal-1: skip if ANY position is metal/donor (legacy frozen check
        # via union), AND skip when the centre b is a donor's first-shell
        # neighbour.  The central-atom gradient is the one that rotates b
        # around the donor (Newton's 3rd: g_b = -(g_a+g_c), pointing along
        # the angle bisector tangent) — exactly the motion that mddir
        # detects as a flipped lone-pair axis.  End-position (a/c) angle
        # terms involving shell-1 atoms are still emitted so they keep
        # contributing to F3_angle / funcgrp via the symmetric pair where
        # they ARE the centre of another triple.
        if a in frozen_or_donor or b in frozen_or_donor or c in frozen_or_donor:
            continue
        # Hapto-class protection: skip every angle that touches a hapto-π
        # atom in ANY position.  The CCDC priors for these angles are
        # dominated by non-hapto chemistry and the Mahalanobis pull
        # destroys η-coordination geometry.
        if hapto_set and (a in hapto_set or b in hapto_set or c in hapto_set):
            continue
        # Pre-compute geometry/hyb that we need either for the lookup or
        # for the Option-B coverage check.  Doing it once here keeps the
        # determinism / FP-order guarantees identical to the legacy path.
        atom_a = mol.GetAtomWithIdx(a)
        atom_b = mol.GetAtomWithIdx(b)
        atom_c = mol.GetAtomWithIdx(c)
        z1, z2, z3 = atom_a.GetSymbol(), atom_b.GetSymbol(), atom_c.GetSymbol()
        hyb1 = _hyb_str(atom_a)
        hyb2 = _hyb_str(atom_b)
        hyb3 = _hyb_str(atom_c)
        rmin = _ring_size_min(mol, b)
        aromatic = _atom_is_aromatic(atom_b)

        if b in donor_shell1:
            # Heal-1b (amine-H regression fix): retain the central-skip
            # only for heavy-only triples — when at least one endpoint is
            # hydrogen, the term is the C-N-H / H-X-H priors that anchor
            # H orientation.  Letting those drop is what produced the
            # +100% amine_h_pct_files_viol regression in Heal-1's smoke-50.
            # Gradient note: H endpoints carry the same Newton-3 g_b
            # term, but H mass is light and the loss-induced rotation of
            # the heavy-only triangle (mddir failure mode) is dominated
            # by heavy endpoints; H-endpoint triples just pin H direction.
            ai_a = atom_a.GetAtomicNum()
            ai_c = atom_c.GetAtomicNum()
            if ai_a != 1 and ai_c != 1:
                # Option-B (2026-06-01): adaptive coverage-aware protection.
                # When a library is available, query the angle distribution
                # for THIS fragment class.  If the lookup finds an entry with
                # n >= adaptive_min_n the CCDC pressure on the funcgrp
                # internals is informative — keep the term.  If the lookup
                # returns None or hits a sparse entry, fall back to the
                # Heal-1/1b heavy-only skip so we don't pollute mddir.
                # When library is None we strictly reproduce Heal-1/1b
                # (the gating below is False, so we drop unconditionally).
                keep_for_coverage = False
                if adaptive_shell1 and lib is not None:
                    try:
                        cov_hit = lib.lookup_angle(
                            z1, z2, hyb2, z3, rmin, aromatic,
                            hyb1=hyb1, hyb3=hyb3,
                            min_n=int(adaptive_min_n),
                        )
                    except Exception:
                        cov_hit = None
                    if cov_hit is not None and int(cov_hit[2]) >= int(adaptive_min_n):
                        keep_for_coverage = True
                if not keep_for_coverage:
                    continue
        result.n_angle_candidates += 1

        hit = lib.lookup_angle(
            z1, z2, hyb2, z3, rmin, aromatic,
            hyb1=hyb1, hyb3=hyb3, min_n=min_n,
        )
        if hit is None:
            continue
        mu, sigma, n = hit
        weight = sparse_downweight(w["angle"], n, n_floor=min_n)
        result.angle_terms.append(
            AngleTerm(a=a, b=b, c=c, mu=mu, sigma=sigma, weight=weight)
        )
        result.n_angle_matched += 1

    # ---- 3) Impropers ------------------------------------------------------
    # Center atom with exactly 3 heavy neighbours and a hybridisation that
    # is sp2 (true planar centres) — sp3 centres legitimately have non-zero
    # OOP and are also represented in the library, so we include them too,
    # but skip sp / sp3d / sp3d2 / unspecified.
    # We carry the "drop-pending" flag for shell-1 centres into the second
    # loop so the Option-B coverage lookup can override the unconditional
    # heavy-only Heal-1/1b drop with a kept-term when the library shows
    # good coverage for this specific improper class.
    improper_specs: List[Tuple[int, Tuple[int, int, int], bool]] = []
    for atom in mol.GetAtoms():
        c_idx = int(atom.GetIdx())
        # Heal-1: skip improper centre if it is metal/donor (legacy) OR a
        # donor's first-shell neighbour (the improper-OOP gradient on c
        # tilts it out of its current plane relative to the donor, which
        # is what re-orients the M-D-X vector and triggers mddir).
        if c_idx in frozen_or_donor:
            continue
        # Hapto-class protection: skip if the centre is a hapto-π atom.
        # (Neighbour-touching impropers are filtered in the lookup loop
        # below, mirroring the existing frozen_or_donor neighbour check.)
        if hapto_set and c_idx in hapto_set:
            continue
        nbrs = sorted(int(n.GetIdx()) for n in atom.GetNeighbors())
        if len(nbrs) != 3:
            continue
        heavy_only_shell1 = False
        if c_idx in donor_shell1:
            # Heal-1b (amine-H regression fix, improper symmetry case):
            # only drop the heavy-only impropers centred on shell-1.  If
            # ANY of the 3 neighbours is hydrogen, keep the term — it is
            # the planarity / pyramidalisation prior that holds the H
            # in its sp2/sp3-correct orientation (e.g. sp2-CH2 OOP).
            has_h_neighbor = any(
                mol.GetAtomWithIdx(int(ni)).GetAtomicNum() == 1
                for ni in nbrs
            )
            if not has_h_neighbor:
                # Mark for heavy-only shell-1 drop.  Option-B may revoke
                # the drop in the lookup loop below when the library has
                # good coverage for this specific improper class.
                heavy_only_shell1 = True
        hyb_c = _hyb_str(atom)
        if hyb_c not in ("sp2", "sp3"):
            continue
        improper_specs.append((c_idx, (nbrs[0], nbrs[1], nbrs[2]), heavy_only_shell1))
    improper_specs = sorted(set(improper_specs), key=lambda t: (t[0], t[1]))

    for (c_idx, (n1, n2, n3), heavy_only_shell1) in improper_specs:
        if n1 in frozen_or_donor or n2 in frozen_or_donor or n3 in frozen_or_donor:
            continue
        # Hapto-class protection: drop the improper if ANY of the three
        # neighbours is a hapto-π atom (the OOP gradient would tilt the
        # central atom relative to the η-ring, breaking the piano-stool).
        if hapto_set and (n1 in hapto_set or n2 in hapto_set or n3 in hapto_set):
            continue
        atom = mol.GetAtomWithIdx(c_idx)
        zc = atom.GetSymbol()
        hyb_c = _hyb_str(atom)
        # Library keys store neighbours in a CANONICAL order. For impropers
        # we sort by (element, hyb) tuple so the lookup is stable regardless
        # of which neighbour iteration order RDKit returned.
        nbrs_data = []
        for ni in (n1, n2, n3):
            a = mol.GetAtomWithIdx(ni)
            nbrs_data.append((a.GetSymbol(), _hyb_str(a), int(ni)))
        nbrs_data.sort(key=lambda t: (t[0], t[1]))
        zs_sorted = [t[0] for t in nbrs_data]
        hybs_sorted = [t[1] for t in nbrs_data]
        nbr_indices_sorted: Tuple[int, int, int] = tuple(t[2] for t in nbrs_data)  # type: ignore[assignment]
        rmin = _ring_size_min(mol, c_idx)

        if heavy_only_shell1:
            # Option-B coverage gate (2026-06-01).  Drop the heavy-only
            # shell-1 improper unless the library has good coverage
            # (n ≥ adaptive_min_n) for this specific improper class.
            keep_for_coverage = False
            if adaptive_shell1 and lib is not None:
                try:
                    cov_hit = lib.lookup_improper(
                        zc, hyb_c, zs_sorted,
                        ring_size_min=rmin,
                        neighbor_hybs_sorted=hybs_sorted,
                        min_n=int(adaptive_min_n),
                    )
                except Exception:
                    cov_hit = None
                if cov_hit is not None and int(cov_hit[2]) >= int(adaptive_min_n):
                    keep_for_coverage = True
            if not keep_for_coverage:
                continue
        result.n_improper_candidates += 1

        hit = lib.lookup_improper(
            zc, hyb_c, zs_sorted,
            ring_size_min=rmin,
            neighbor_hybs_sorted=hybs_sorted,
            min_n=min_n,
        )
        if hit is None:
            continue
        mu, sigma, n = hit
        weight = sparse_downweight(w["improper"], n, n_floor=min_n)
        result.improper_terms.append(
            ImproperTerm(
                center=c_idx,
                neighbors_3=nbr_indices_sorted,
                mu=mu, sigma=sigma, weight=weight,
            )
        )
        result.n_improper_matched += 1

    # ---- 4) Angle-to-Metal extension (2026-06-05, User direction). --------
    #
    # Default OFF: when neither ``ANGLE_TO_METAL_ENV`` nor ``DMD_POLY_ENV``
    # is set the loop is a pure no-op and the term list is byte-identical
    # with the legacy detector.
    #
    # When ANGLE_TO_METAL is ON:  add M-D-X terms (centre = donor) so the
    #   ligand-X atoms get a Mogul/VSEPR-grounded pull around the M-D axis.
    #   Gradient on M and D is zeroed by grip_polish; only X moves.
    #
    # When DMD_POLY is ON:  add D-M-D' polyhedron-anchored terms (centre =
    #   metal).  Gradient is identically zero (both endpoints + centre are
    #   frozen) but the term contributes to the accept-if-better severity
    #   so polishes that distort the polyhedron lose ground.
    if metal is not None and metal >= 0 and donor_set:
        try:
            if angle_to_metal_active():
                mdx_terms = build_mdx_angle_terms(
                    mol, int(metal), tuple(sorted(int(d) for d in donor_set)),
                    hapto_atoms=hapto_set, library=lib, min_n=min_n,
                    include_hydrogens=_angle_to_metal_include_h(),
                )
                for t in mdx_terms:
                    result.angle_terms.append(t)
                    result.n_angle_matched += 1
                    result.n_angle_candidates += 1
            if dmd_polyhedron_active():
                dmd_terms = build_dmd_polyhedron_terms(
                    mol, int(metal), tuple(sorted(int(d) for d in donor_set)),
                    P=P,
                )
                for t in dmd_terms:
                    result.angle_terms.append(t)
                    result.n_angle_matched += 1
                    result.n_angle_candidates += 1
        except Exception:
            # Defence in depth: a misconfigured fffree mol must not crash
            # the detector.  Both env flags default OFF so any failure
            # here simply drops the new terms and yields the legacy result.
            pass

    if return_result:
        return result
    return result.all_terms()
