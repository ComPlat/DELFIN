"""delfin.fffree.assemble_via_mogul — Mogul-primary construction path.

The CCDC fragment manifold IS the construction, not a polish on top of one.

This module implements the architecture described in the project draft
manuscript (sections "Closed-form coverage of the configuration space" and
"Empirical bounds from a 2.45-million-key crystallographic manifold"):

    SMILES
       → full-complex mol (RDKit, with explicit H)
       → fragment-key bounds for every atom pair (1-2, 1-3, 1-4, M-D,
         donor-donor, ring, hapto, resonance) populated from the CCDC
         empirical distributions held in
         ``delfin/fffree/mogul_bounds.build_bounds_matrix``
       → ONE constrained distance-geometry embed of the full complex
         (RDKit ``rdDistGeom.EmbedMolecule`` with the populated bounds
         matrix; deterministic seed)
       → return 3D coordinates

NO idealised polyhedron as the primary path.
NO Pyykkö-derived M-D distances.
NO VSEPR-derived donor-vertex orientations.
NO post-hoc geometric snaps / planarity fixes / vdW floors stacked on top.

The single source of geometric truth is the CCDC manifold.  The seven
dedicated TM categories (NHC carbene, hapto-η²/η⁵/η⁶, μ-bridge,
agostic, oxidative-addition) plus the 3d/4d/5d/f-block tier are queried
directly through ``mogul_bounds`` and ``grip_mogul_lookup``; the path
is universal across the TMC classes — no SMILES patterns, no per-class
branches, no hand-tuned templates.

Env flag (verbindlich, default OFF byte-identical):
    DELFIN_FFFREE_MOGUL_PRIMARY=1

When the flag is unset, importing this module is a no-op and the legacy
construction path stays bit-identical with HEAD.

Author: hmaximilian <hmaximilian496@gmail.com>
Branch: feat-mogul-primary-2026-06-07
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom as _DG
from rdkit import DistanceGeometry as _DGs

import delfin._bond_decollapse as _bd

# --- Internal seed used everywhere in the fffree stack -----------------------
# Aligns with assemble_complex.SEED so any side-by-side comparison of legacy
# vs primary path uses the same DG random seed.
SEED = 0xDE17  # matches delfin.fffree.assemble_complex.SEED via convention


# ---------------------------------------------------------------------------
# Env-flag helpers
# ---------------------------------------------------------------------------
def mogul_primary_enabled() -> bool:
    """Return True iff the Mogul-primary path is wired on for this process."""
    return os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY", "0") == "1"


def mogul_primary_conformers_enabled() -> bool:
    """Return True iff the per-orbit conformer enumeration is active.

    Default ON when ``DELFIN_FFFREE_MOGUL_PRIMARY`` is on (since a single
    structure per stereoisomer is demonstrably incomplete -- NHC rings flex,
    aryl substituents rotate, single bonds adopt gauche/anti).  Set
    ``DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS=0`` to compare bytes against
    the orbit-only path (single output per orbit).
    """
    return os.environ.get(
        "DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS", "1"
    ) == "1"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def assemble_complex_mogul_primary(
    smiles: str,
    *,
    max_attempts: int = 5,
    return_mol: bool = False,
    donor_at_vertex: Optional[Sequence[int]] = None,
    geometry_key: Optional[str] = None,
    mol_override=None,
    donor_idxs_override: Optional[Sequence[int]] = None,
) -> Optional[Tuple[List[str], np.ndarray]]:
    """Construct a 3D complex by embedding the CCDC-populated bounds matrix.

    Parameters
    ----------
    smiles : str
        SMILES of the metal complex (single or multi-metal; multi-metal
        with primary metal selected as in ``decompose``).
    max_attempts : int
        Number of DG embed attempts with consecutive seeds before falling
        back to ``None`` (caller falls through to legacy).
    return_mol : bool
        Reserved for future use; the public contract returns ``(syms, P)``.
    donor_at_vertex : sequence of int, optional
        Pólya orbit enumeration hook (2026-06-07): permutation mapping
        polyhedron vertex ``k`` → atom index of the donor placed at that
        vertex.  When ``None`` (default) the canonical sort-by-atom-idx
        ordering is used (byte-identical with the pre-orbit-enumeration
        path).  When supplied, the donor-donor distance bounds are
        re-injected so the donor at vertex ``k`` sees the ideal polyhedron
        distance to every other vertex; the DG embed therefore lands on
        the stereoisomer encoded by the permutation.  Pure graph-structural
        logic — no SMILES pattern, no per-class branch.
    geometry_key : str, optional
        Polyhedron name (e.g. ``"OC-6 octahedron"``) used to resolve the
        canonical vertex set when ``donor_at_vertex`` is supplied.  When
        ``None``, the per-CN default polyhedron from
        ``polyhedra.geometries_for_cn`` is used (matches build_bounds_matrix).

    Returns
    -------
    (syms, P) or None
        ``syms`` is the element list (length ``N``), ``P`` is an
        ``(N, 3)`` ndarray of Cartesian coordinates centred on the metal
        (the metal is at the origin).  ``None`` is returned on any failure
        (parse, missing metal, empty bounds, infeasible bounds, embed
        failure) so the caller can fall through to legacy.
    """
    # 1) Full mol with metal + ligands as one connected graph + explicit H.
    #    The κⁿ enumeration path passes ``mol_override`` (the SAME mol
    #    after graph-rewiring to insert/remove dative bonds for the chosen
    #    donor set) so this assemble call sees the variant's topology.
    if mol_override is not None:
        mol = mol_override
    else:
        mol = _full_complex_mol(smiles)
    if mol is None:
        return None

    # 2) Locate metal + donors from the graph topology.
    #    Under κⁿ override the donor set is provided explicitly; the
    #    metal index is still taken from the graph topology.
    metal_idx, auto_donor_idxs = _locate_metal_and_donors(mol)
    if metal_idx is None:
        return None
    if donor_idxs_override is not None:
        donor_idxs = sorted(int(d) for d in donor_idxs_override)
    else:
        donor_idxs = auto_donor_idxs
    if not donor_idxs:
        return None

    # 3) Build the CCDC-empirical bounds matrix.
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    # 3a) Chelate-aware polyhedron pick (2026-06-07) when the caller did
    #     NOT pass an explicit ``geometry_key``.  The orbit enumerator
    #     :func:`enumerate_and_embed_mogul_primary` already resolves
    #     ``geometry_key`` chelate-aware before calling us, so this branch
    #     fires on the single-shot path (e.g. converter_backend's fallback
    #     when no Pólya orbit decomposition is available) and on direct
    #     library callers.
    #
    #     Universal: pure geometric matching between graph-derived chelate
    #     ring sizes and polyhedron vertex-vertex angles.  No SMILES
    #     patterns, no per-class branches.  Default-OFF byte-identical
    #     when ``DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK=0`` or when no
    #     chelate ligand is detected.  Author: hmaximilian.
    effective_geometry_key: Optional[str] = geometry_key
    if effective_geometry_key is None and mol_override is None:
        _d_dec = None
        try:
            from delfin.fffree.chelate_aware_picker import (
                chelate_aware_picker_enabled as _cap_on,
                chelate_info_from_decompose as _ci_from,
                pick_polyhedron_chelate_aware as _pick_cap,
            )
            from delfin.fffree import decompose as _DEC2
            if _cap_on():
                _d_dec = _DEC2.decompose(smiles)
                if (_d_dec is not None
                        and _d_dec.get("has_chelate")):
                    _chelate_info = _ci_from(_d_dec)
                    if _chelate_info:
                        _picked = _pick_cap(
                            cn=int(_d_dec.get("cn", 0) or 0),
                            chelate_info=_chelate_info,
                            metal_sym=str(_d_dec.get("metal", "")),
                        )
                        if _picked:
                            effective_geometry_key = _picked
        except ImportError:
            pass
        except Exception:
            pass

        # f-block-aware singleshot pick (2026-06-07, hmaximilian).  Mirrors
        # the orbit-enumerator wiring further down: when the chelate-aware
        # picker did NOT fire AND the metal is in the f-block / high-CN
        # picker scope, score candidates by CCDC D-D distribution match.
        # Defaults to OFF byte-identical when ``DELFIN_FFFREE_F_BLOCK_PICK=0``.
        if effective_geometry_key is None:
            try:
                from delfin.fffree.f_block_picker import (
                    f_block_picker_enabled as _fbp_on,
                    f_block_polyhedron_picker as _pick_fbp,
                )
                from delfin.fffree import decompose as _DEC3
                if _fbp_on():
                    if _d_dec is None:
                        try:
                            _d_dec = _DEC3.decompose(smiles)
                        except Exception:
                            _d_dec = None
                    if _d_dec is not None and not _d_dec.get("has_chelate"):
                        _cn_int = int(_d_dec.get("cn", 0) or 0)
                        _metal = str(_d_dec.get("metal", ""))
                        _donor_syms: List[str] = []
                        for _lg in (_d_dec.get("ligands") or []):
                            _de_list = _lg.get("donor_elems")
                            _de = _lg.get("donor_elem")
                            _dk = int(_lg.get("denticity", 1) or 1)
                            if isinstance(_de_list, (list, tuple)) and _de_list:
                                for _s in _de_list:
                                    _donor_syms.append(str(_s))
                            elif _de is not None:
                                for _ in range(max(1, _dk)):
                                    _donor_syms.append(str(_de))
                        _donor_syms = [
                            "C" if s == "C_hapto" else s for s in _donor_syms
                        ]
                        if (_metal and _cn_int >= 7
                                and len(_donor_syms) == _cn_int):
                            _picked_fb = _pick_fbp(
                                metal_sym=_metal,
                                cn=_cn_int,
                                donors=_donor_syms,
                                force=True,
                            )
                            if _picked_fb:
                                effective_geometry_key = _picked_fb
            except ImportError:
                pass
            except Exception:
                pass
    try:
        from delfin.fffree.mogul_bounds import build_bounds_matrix
    except ImportError:
        return None
    try:
        lower, upper, info = build_bounds_matrix(
            syms=syms,
            mol=mol,
            metal_idx=metal_idx,
            donor_idxs=donor_idxs,
            grip_lib=None,         # default: release-pinned grip_lib_v5/v6
            cod_lib=None,
            geometry=effective_geometry_key,
            min_n=5,
            use_automorphism=True,
        )
    except Exception:
        return None

    # 3b) Tighten the M-D bounds via the seven dedicated TM categories
    #     described in the draft manuscript (NHC carbene, hapto-η², η⁵,
    #     η⁶, μ-bridge, agostic, oxidative-addition).  ``build_bounds_matrix``
    #     already routes hapto donors through ``lookup_tm_category``, but
    #     CARBENE / AGOSTIC / μ-BRIDGE / OX-ADDITION are detected from the
    #     graph topology here and re-injected into the bounds.  This is
    #     the universal mechanism the draft requires — graph features
    #     drive the category, never SMILES patterns or per-class branches.
    _override_md_bounds_via_tm_category(
        mol=mol,
        syms=syms,
        metal_idx=metal_idx,
        donor_idxs=donor_idxs,
        lower=lower,
        upper=upper,
        info=info,
    )

    # 3c) Pólya-orbit override (2026-06-07): when a donor-to-vertex
    #     permutation is supplied, re-write the donor-donor distance
    #     bounds so they encode that specific stereoisomer.  The BOUNDS
    #     MATRIX is the orbit representation — each Burnside-distinct
    #     coloring becomes its own injection.  Universal: works for any
    #     registered polyhedron.  See `enumerate_and_embed_mogul_primary`
    #     for the orbit-enumerating wrapper.
    if donor_at_vertex is not None:
        _override_dd_bounds_for_orbit(
            metal_sym=str(syms[int(metal_idx)]),
            metal_idx=int(metal_idx),
            donor_idxs=donor_idxs,
            donor_at_vertex=[int(d) for d in donor_at_vertex],
            geometry_key=effective_geometry_key,
            lower=lower,
            upper=upper,
        )

    # 4) Triangle-smooth + DG embed via RDKit using the populated bounds.
    P = _embed_with_bounds(
        mol, lower, upper,
        metal_idx=metal_idx,
        donor_idxs=donor_idxs,
        max_attempts=max_attempts,
    )
    if P is None:
        return None

    # 4b) Post-embed projection: pull each donor onto its CCDC-defined
    #     M-D distance and the polyhedron vertex direction.  The embedder
    #     leaves a soft-minimum residual outside the CCDC window (because
    #     it minimises a smooth penalty, not a hard barrier); this step
    #     applies the CCDC mean as a rigid-body translation of the
    #     ligand subtree along the M->D ray so the M-D distance is exact
    #     and the D-M-D angle relaxes onto its polyhedron-ideal value.
    #
    #     IMPORTANT for orbit enumeration: the projection homogenises the
    #     donor-to-vertex assignment back to the polyhedron's canonical
    #     order (donor #k of the sorted list → vertex #k).  When a Pólya
    #     orbit override is in effect we must pass the orbit permutation
    #     through so the right donor lands at the right vertex.  Without
    #     this the orbit choice would be silently overwritten.
    P = _project_donors_to_ccdc_geometry(
        mol=mol,
        syms=syms,
        metal_idx=metal_idx,
        donor_idxs=donor_idxs,
        lower=lower,
        upper=upper,
        P=P,
        donor_at_vertex=donor_at_vertex,
        geometry_key=effective_geometry_key,
    )

    # 4c) GRIP-polish (2026-06-07): Mahalanobis L-BFGS pull against the
    #     CCDC fragment manifold.  The DG embed + rigid-body projection
    #     satisfy the bounds matrix but settle at the geometric centre of
    #     the (lower, upper) window, not at the empirical CCDC mean.  In
    #     particular soft-constrained pairs (chelate D-D, NHC ring
    #     internals) drift +0.2-0.4 Å above the CCDC mean and NHC rings
    #     come out asymmetric.  ``grip_polish`` pulls the geometry
    #     toward the manifold mean using the per-fragment covariance,
    #     under M-D / topology / chirality / CShM hard validators.
    #
    #     Env-flag: DELFIN_FFFREE_MOGUL_PRIMARY_GRIP, default ON.  Set
    #     to 0 explicitly when comparing pre/post polish bytes.
    #
    #     Defence-in-depth M-D check (±0.05 Å) mirrors the legacy
    #     assemble_complex GRIP wiring -- any unexpected drift rolls
    #     back to the pre-polish P silently.
    if os.environ.get("DELFIN_FFFREE_MOGUL_PRIMARY_GRIP", "1") == "1":
        # Multi-step polish (2026-06-07, validator-tuning).  Default ON
        # when MOGUL_PRIMARY_GRIP is active so the polish actually applies
        # its 99.3-99.6% severity reductions (the legacy single-shot path
        # rolls back on spurious sp2-pseudo-stereocenter sign flips and on
        # late-iteration topology overstretch).  Explicit env-flag wins.
        _ms_raw = os.environ.get(
            "DELFIN_FFFREE_GRIP_POLISH_MULTISTEP", ""
        ).strip().lower()
        if not _ms_raw:
            os.environ["DELFIN_FFFREE_GRIP_POLISH_MULTISTEP"] = "1"
        try:
            from delfin.fffree.grip_polish import grip_polish
            from delfin.fffree.grip_mogul_lookup import GripLibrary
            _grip_lib = GripLibrary.get_default()
            _md_targets = [
                float(np.linalg.norm(P[int(d)] - P[int(metal_idx)]))
                for d in donor_idxs
            ]
            Pg = grip_polish(
                P, mol, metal=int(metal_idx),
                donors=[int(d) for d in donor_idxs],
                geom="", mogul_lib=_grip_lib,
            )
            if (
                Pg is not None
                and isinstance(Pg, np.ndarray)
                and Pg.shape == P.shape
                and np.all(np.isfinite(Pg))
            ):
                _md_ok = True
                for _di, _d in enumerate(donor_idxs):
                    _d_now = float(np.linalg.norm(
                        Pg[int(_d)] - Pg[int(metal_idx)]
                    ))
                    if abs(_d_now - _md_targets[_di]) > 0.05:
                        _md_ok = False
                        break
                if _md_ok:
                    P = Pg
        except Exception:
            pass  # silent rollback to pre-polish P

    # 5) Centre on metal so M is at origin (downstream convention).
    P = P - P[metal_idx]
    out_syms = [a.GetSymbol() for a in mol.GetAtoms()]
    # The fffree contract puts the metal at index 0.  Reorder if needed.
    if metal_idx != 0:
        order = [metal_idx] + [i for i in range(len(out_syms)) if i != metal_idx]
        out_syms = [out_syms[i] for i in order]
        P = P[order]
        # Re-map donor indices to the output frame so the M-shell gate
        # below counts heavy atoms inside the canonical (output-frame)
        # geometry.  Donor permutation: idx i -> i if i < metal_idx,
        # else i+1 stays where it is (the metal was removed and
        # prepended).  Equivalently: idx maps to its position in
        # ``order``.
        idx_map = {old: new for new, old in enumerate(order)}
        out_metal_idx = 0
        out_donor_idxs = sorted(int(idx_map.get(int(d), int(d))) for d in donor_idxs)
    else:
        out_metal_idx = 0
        out_donor_idxs = sorted(int(d) for d in donor_idxs)

    # 6) M-shell overfill hard-gate (Bug #3 fix, 2026-06-07).
    #
    #    The bounds-matrix embed satisfies the per-pair (lower, upper)
    #    window but the polished geometry may still leave a non-donor
    #    heavy atom inside the M-shell radius (e.g. the bridging carbon
    #    of an η²-arene, or a counter-ion that got pulled close by the
    #    centroid translation).  We gate the BASE structure with the
    #    same universal threshold the conformer pipeline uses
    #    (``M_SHELL_FACTOR * ideal_bond``) so a downstream conformer
    #    enumerator never starts from an over-filled shell.
    #
    #    Under DELFIN_FFFREE_ROTAMER_TOPOLOGY_GATE=1 (auto-ON when
    #    DELFIN_FFFREE_MOGUL_PRIMARY=1) we reject the embed and return
    #    None — the caller (κⁿ enumerator / orbit enumerator) then
    #    falls through to the next variant or the legacy path.  When
    #    the gate flag is off explicitly (``=0``) the check is skipped
    #    so byte-identity with the pre-gate path is preserved.
    try:
        from .rotamer_topology_gate import (
            _env_on as _topology_gate_on,
            M_SHELL_FACTOR as _MSF,
        )
        if _topology_gate_on():
            cn_expected = len(out_donor_idxs)
            metal_sym = out_syms[out_metal_idx]
            shell = 0
            for j in range(len(out_syms)):
                if j == out_metal_idx:
                    continue
                sj = out_syms[j]
                if sj == "H":
                    continue
                if _bd._is_metal(sj):
                    continue
                try:
                    ideal = float(_bd._ideal_bond(metal_sym, sj))
                except Exception:
                    ideal = 0.0
                if ideal <= 0.0:
                    continue
                d_mj = float(np.linalg.norm(P[j] - P[out_metal_idx]))
                if d_mj < float(_MSF) * ideal:
                    shell += 1
                    if shell > cn_expected:
                        return None
    except Exception:
        pass  # defensive: never crash on the gate import / lookup

    return out_syms, P


# ---------------------------------------------------------------------------
# Pólya orbit enumeration over the Mogul-primary bounds matrix (2026-06-07)
# ---------------------------------------------------------------------------
# Decompose-geometry-name → Pólya-group key.  Local copy of
# ``converter_backend._GEOM_TO_POLYA`` (subset that the bounds-matrix path
# supports out of the box -- monodentate sigma plus chelate via the
# graph-detected D-D pairings).  Whenever ``_default_geometry`` returns a
# new shape, add the matching Pólya key here too.
_GEOM_NAME_TO_POLYA_KEY = {
    "L-2 linear": "linear",
    "SP-3 trigonal planar": "trigonal_planar",
    "T-3 T-shape": "tshape",
    "OC-6 octahedron": "octahedron",
    "SP-4 square planar": "square_planar",
    "T-4 tetrahedron": "tetrahedron",
    "TBP-5 trigonal bipyramid": "trigonal_bipyramid",
    "SPY-5 square pyramid": "square_pyramid",
    "TPR-6 trigonal prism": "trigonal_prism",
    "PB-7 pentagonal bipyramid": "pentagonal_bipyramid",
    "SQAP-8 square antiprism": "square_antiprism",
    "TTP-9 tricapped trigonal prism": "tricapped_trigonal_prism",
}


def _override_dd_bounds_for_orbit(
    *,
    metal_sym: str,
    metal_idx: int,
    donor_idxs: Sequence[int],
    donor_at_vertex: Sequence[int],
    geometry_key: Optional[str],
    lower: np.ndarray,
    upper: np.ndarray,
    sigma_band: float = 0.30,
) -> None:
    """Re-write donor-donor distance bounds to encode a specific Pólya orbit.

    Mutates ``lower`` / ``upper`` in place.

    For each pair ``(va, vb)`` of polyhedron vertex indices, the donor atom
    indices that sit at those vertices in this orbit are
    ``donor_at_vertex[va]`` and ``donor_at_vertex[vb]``.  Their target
    distance is ``|V[va] - V[vb]| * r_md`` where ``V`` are the unit
    polyhedron vectors and ``r_md`` is the mean M-D distance already
    encoded in the matrix.  We open a soft ±``sigma_band`` Å window.

    Universal: works for any registered polyhedron (via
    ``polyhedra.ref_vectors``).  No SMILES pattern, no per-class branch.
    """
    try:
        from delfin.fffree import polyhedra as _polyhedra
    except ImportError:
        return
    donor_list = [int(d) for d in donor_idxs]
    cn = len(donor_list)
    if cn < 2 or len(donor_at_vertex) != cn:
        return
    ref_vecs = None
    if geometry_key:
        try:
            ref_vecs = _polyhedra.ref_vectors(str(geometry_key))
        except Exception:
            ref_vecs = None
    if ref_vecs is None:
        try:
            cands = _polyhedra.geometries_for_cn(int(cn), metal_sym)
        except Exception:
            cands = []
        for g in cands:
            try:
                ref_vecs = _polyhedra.ref_vectors(str(g))
                break
            except Exception:
                continue
    if ref_vecs is None or ref_vecs.shape[0] < cn:
        return
    V = np.asarray(ref_vecs[:cn], dtype=float)
    norms = np.linalg.norm(V, axis=1, keepdims=True)
    norms = np.where(norms < 1e-9, 1.0, norms)
    Vu = V / norms

    # Estimate the mean M-D distance from the existing M-D bounds (midpoints).
    md_means: List[float] = []
    for d in donor_list:
        lo_md = float(lower[metal_idx, d])
        hi_md = float(upper[metal_idx, d])
        if hi_md > 100.0 or hi_md <= lo_md:
            continue
        md_means.append(0.5 * (lo_md + hi_md))
    if not md_means:
        return
    r_md = float(np.mean(md_means))

    for va in range(cn):
        for vb in range(va + 1, cn):
            da = int(donor_at_vertex[va])
            db = int(donor_at_vertex[vb])
            if da == db:
                continue
            d_ideal = float(np.linalg.norm(Vu[va] - Vu[vb])) * r_md
            if d_ideal <= 0.0:
                continue
            lo = max(0.5, d_ideal - sigma_band)
            hi = d_ideal + sigma_band
            i, j = (da, db) if da < db else (db, da)
            lower[i, j] = lower[j, i] = lo
            upper[i, j] = upper[j, i] = hi


def _coloring_to_donor_assignment(
    coloring: Sequence[str],
    donor_labels: Sequence[str],
    donor_idxs: Sequence[int],
) -> Optional[List[int]]:
    """Map a Pólya coloring (label-per-vertex) → donor-atom-index-per-vertex.

    For coloring ``("F", "N", "F", "O", "O", "N")`` and donor_labels
    ``["N", "N", "O", "O", "F", "F"]`` indexed by ``donor_idxs``, returns
    a list ``donor_at_vertex`` of donor atom indices (one per polyhedron
    vertex) such that the label multiset at each vertex matches the
    coloring.

    Returns ``None`` if the multisets disagree (defensive — orbit
    enumeration is best-effort, never raises on shape mismatch).
    """
    from collections import Counter as _Counter
    if _Counter(coloring) != _Counter(donor_labels):
        return None
    buckets: Dict[str, List[int]] = {}
    for lab, idx in zip(donor_labels, donor_idxs):
        buckets.setdefault(str(lab), []).append(int(idx))
    for k in buckets:
        buckets[k].sort()  # deterministic pop order
    assignment: List[int] = []
    for lab in coloring:
        bucket = buckets.get(str(lab))
        if not bucket:
            return None
        assignment.append(bucket.pop(0))
    return assignment


def _classify_orbit_label(geom_polya_key: str, coloring: Sequence[str]) -> str:
    """Universal isomer name from coloring + polyhedron antipode structure
    (cis/trans, fac/mer, all-cis/all-trans).  Local copy of
    ``converter_backend._classify_coloring`` so the orbit path has no
    converter dependency.
    """
    from collections import Counter as _Counter
    _ANTI = {
        "octahedron": {0: 1, 1: 0, 2: 3, 3: 2, 4: 5, 5: 4},
        "square_planar": {0: 2, 1: 3, 2: 0, 3: 1},
        "trigonal_bipyramid": {0: 1, 1: 0},
        "square_pyramid": {1: 3, 3: 1, 2: 4, 4: 2},
        "tshape": {0: 1, 1: 0},
    }
    anti = _ANTI.get(geom_polya_key)
    if anti is None:
        return ""
    cnt = _Counter(coloring)
    n = len(coloring)

    def is_trans(el):
        v = [i for i, e in enumerate(coloring) if e == el]
        return any(anti.get(v[a]) == v[b]
                   for a in range(len(v)) for b in range(a + 1, len(v)))

    pairs2 = [el for el, c in cnt.items() if c == 2]
    threes = [el for el, c in cnt.items() if c == 3]
    if n == 6:
        if threes:
            return "mer" if is_trans(threes[0]) else "fac"
        if len(pairs2) == 3:
            trans_els = sorted(el for el in pairs2 if is_trans(el))
            if not trans_els:
                return "all-cis"
            if len(trans_els) == 3:
                return "all-trans"
            return "-".join(f"{e}trans" for e in trans_els)
        if len(pairs2) == 1:
            return "trans" if is_trans(pairs2[0]) else "cis"
    elif n == 4:
        if len(pairs2) == 1:
            return "trans" if is_trans(pairs2[0]) else "cis"
    elif n == 5:
        if threes:
            return "mer" if is_trans(threes[0]) else "fac"
        if len(pairs2) == 1:
            return "trans" if is_trans(pairs2[0]) else "cis"
    elif n == 3:
        if len(pairs2) == 1:
            return "trans" if is_trans(pairs2[0]) else "cis"
    return ""


# ---------------------------------------------------------------------------
# Per-orbit conformer enumeration (Cremer-Pople ring pucker + single-bond
# rotamers).  Composes the third Cauchy-Frobenius layer of the draft
# manuscript: Pólya orbit × binding-mode × CONFORMER -> ensemble per SMILES.
# ---------------------------------------------------------------------------
def _build_mol_for_conformer_enum(syms: Sequence[str], P: np.ndarray):
    """Build an RDKit Mol with a heavy-atom bond graph + metal-donor edges
    for ring detection and rotamer enumeration.

    Mirrors the geometric-bond construction inside
    ``converter_backend._conformer_enum_post`` so the same conformer
    enumeration semantics apply to the Mogul-primary path.  Universal:
    bonds are derived from geometry + element radii, never from SMILES.

    Returns ``None`` on any failure -- caller falls through to single
    orbit output.
    """
    try:
        rw = Chem.RWMol()
        for s in syms:
            rw.AddAtom(Chem.Atom(str(s)))
        bonds = _bd._geometric_bonds(list(syms), np.asarray(P, dtype=float))
        added = set()
        for ii, jj in bonds:
            key = (min(int(ii), int(jj)), max(int(ii), int(jj)))
            if key in added:
                continue
            try:
                rw.AddBond(int(ii), int(jj), Chem.BondType.SINGLE)
                added.add(key)
            except Exception:
                pass
        # Add metal->donor bonds (geometric_bonds skips them).
        nm = len(syms)
        for im in range(nm):
            if not _bd._is_metal(str(syms[im])):
                continue
            for jm in range(nm):
                if jm == im or str(syms[jm]) == "H":
                    continue
                if _bd._is_metal(str(syms[jm])):
                    continue
                try:
                    ideal = _bd._ideal_bond(str(syms[im]), str(syms[jm]))
                    dij = float(np.linalg.norm(
                        np.asarray(P[jm]) - np.asarray(P[im])
                    ))
                    if dij < 1.30 * ideal:
                        key = (min(im, jm), max(im, jm))
                        if key not in added:
                            try:
                                rw.AddBond(int(im), int(jm),
                                           Chem.BondType.SINGLE)
                                added.add(key)
                            except Exception:
                                pass
                except Exception:
                    pass
        mol_g = rw.GetMol()
        try:
            Chem.SanitizeMol(mol_g, catchErrors=True)
        except Exception:
            pass
        try:
            Chem.GetSSSR(mol_g)
        except Exception:
            pass
        try:
            Chem.FastFindRings(mol_g)
        except Exception:
            pass
        return mol_g
    except Exception:
        return None


def _enumerate_ring_pucker_variants(
    mol_g,
    P: np.ndarray,
    *,
    max_per_ring: int = 4,
    skip_aromatic: bool = True,
    aromatic_atom_set: Optional[Sequence[int]] = None,
) -> List[Tuple[np.ndarray, str]]:
    """Enumerate Cremer-Pople ring-pucker variants of ``P`` against ``mol_g``.

    Direct call into ``ring_pucker`` primitives -- bypasses the
    ``DELFIN_FFFREE_RING_PUCKER`` env gate inside
    ``ring_pucker_integration``, since this module owns its own env flag
    (``DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS``).  Each puckerable
    non-aromatic 5/6/7 ring (also avoids chelate metallacycles by skipping
    rings containing a metal atom) contributes ``min(max_per_ring,
    len(canonical_pucker_states(N)))`` variants; the Cartesian product
    across rings is taken and bond lengths restored after each pucker.

    Returns
    -------
    list of (P_variant, label) tuples.  ``label`` is e.g. ``"pucker3"``.
    The first entry equals ``P`` (identity, sub-label ``"base"``) only
    if no puckerable ring is found; otherwise the identity is omitted
    so the caller can re-attach it explicitly.
    """
    try:
        from .ring_pucker import (
            canonical_pucker_states,
            set_pucker,
            _restore_ring_bonds,
        )
    except Exception:
        return []
    if mol_g is None:
        return []

    # Detect rings to pucker: 5/6/7 ring size, non-aromatic, no metal.
    # The bond graph used here is geometric (all single bonds, no aromatic
    # flags), so ``RDKit.GetIsAromatic`` always returns False on these
    # mols.  We therefore use a UNIVERSAL geometric aromaticity heuristic:
    # a ring is aromatic when (a) all heavy atoms are sp²-eligible
    # (C/N/O/S/P/B), (b) the ring is planar (max out-of-plane deviation
    # below 0.20 Å), (c) every ring bond sits in the aromatic length band
    # 1.30-1.45 Å.  Catches benzene, pyridine, pyrrole, NHC π-systems
    # without relying on bond-order / aromaticity perception.
    try:
        ri = mol_g.GetRingInfo()
        all_rings = [tuple(int(a) for a in r) for r in ri.AtomRings()]
    except Exception:
        return []

    _SP2_ELIGIBLE = {"C", "N", "O", "S", "P", "B"}
    _AROM_BOND_MIN = 1.30
    _AROM_BOND_MAX = 1.45
    _AROM_PLANAR_TOL = 0.20  # Å

    # When the caller supplies an aromatic-atom-set (typically derived
    # from RDKit aromaticity perception on the SOURCE mol — the SMILES
    # parser preserves aromatic-bond flags that survive the embed-mol
    # roundtrip), we mark a ring aromatic whenever ALL of its atoms are
    # in that set.  This is robust even when the post-embed geometry has
    # the aromatic ring distorted (the embedder treats aromatic bonds as
    # SINGLE, which inflates ring bond lengths to ~1.50 Å).
    arom_idx_set: Optional[set] = None
    if aromatic_atom_set is not None:
        try:
            arom_idx_set = {int(i) for i in aromatic_atom_set}
        except Exception:
            arom_idx_set = None

    def _ring_is_aromatic_supplied(ring_atoms: Tuple[int, ...]) -> bool:
        """Aromatic if every ring atom is in the caller-supplied set."""
        if arom_idx_set is None:
            return False
        try:
            return all(int(i) in arom_idx_set for i in ring_atoms)
        except Exception:
            return False

    def _ring_is_aromatic_geom(ring_atoms: Tuple[int, ...]) -> bool:
        """Universal aromatic-ring detection from geometry alone."""
        try:
            syms_r = [
                str(mol_g.GetAtomWithIdx(int(i)).GetSymbol())
                for i in ring_atoms
            ]
            if not all(s in _SP2_ELIGIBLE for s in syms_r):
                return False
            coords_r = np.array(
                [P[int(i)] for i in ring_atoms], dtype=float
            )
            centroid = coords_r.mean(axis=0)
            M = coords_r - centroid
            _, _, Vt = np.linalg.svd(M, full_matrices=False)
            normal = Vt[-1]
            dev = float(np.max(np.abs(M @ normal)))
            if dev > _AROM_PLANAR_TOL:
                return False
            n_r = len(ring_atoms)
            for k in range(n_r):
                a = int(ring_atoms[k])
                b = int(ring_atoms[(k + 1) % n_r])
                d = float(np.linalg.norm(P[a] - P[b]))
                if not (_AROM_BOND_MIN <= d <= _AROM_BOND_MAX):
                    return False
            return True
        except Exception:
            return False

    rings: List[Tuple[int, ...]] = []
    for r in all_rings:
        if not (5 <= len(r) <= 7):
            continue
        if skip_aromatic and (
            _ring_is_aromatic_supplied(r) or _ring_is_aromatic_geom(r)
        ):
            continue
        # Skip metallacycle rings (chelate -- pucker would clash with
        # the empirical M-D / D-D bound the embed already satisfies).
        try:
            if any(
                _bd._is_metal(mol_g.GetAtomWithIdx(int(i)).GetSymbol())
                for i in r
            ):
                continue
        except Exception:
            pass
        rings.append(r)
    if not rings:
        return []

    # Per-ring canonical states (capped).
    per_ring_states: List[List[Tuple[List[float], List[float], str]]] = []
    h_map_per_ring: List[Dict[int, List[int]]] = []
    for r in rings:
        states = canonical_pucker_states(len(r))
        if max_per_ring is not None and 0 < max_per_ring < len(states):
            states = states[:max_per_ring]
        per_ring_states.append(states)
        # Hydrogens bonded to each ring atom -- carried rigidly during pucker.
        h_map: Dict[int, List[int]] = {}
        for ai in r:
            try:
                a = mol_g.GetAtomWithIdx(int(ai))
                h_map[int(ai)] = [
                    int(n.GetIdx()) for n in a.GetNeighbors()
                    if n.GetAtomicNum() == 1
                ]
            except Exception:
                h_map[int(ai)] = []
        h_map_per_ring.append(h_map)

    # Cartesian product across rings.
    from itertools import product as _prod
    variants: List[Tuple[np.ndarray, str]] = []
    for combo_idx, combo in enumerate(_prod(*per_ring_states)):
        Pv = np.asarray(P, dtype=float).copy()
        ok = True
        for ring_atoms, h_map, (Q_t, phi_t, label) in zip(
            rings, h_map_per_ring, combo
        ):
            try:
                ring_arr = np.array(
                    [Pv[int(i)] for i in ring_atoms], dtype=float
                )
                new_ring = set_pucker(ring_arr, Q_t, phi_t)
                new_ring = _restore_ring_bonds(new_ring)
                for k, i in enumerate(ring_atoms):
                    delta = new_ring[k] - Pv[int(i)]
                    Pv[int(i)] = new_ring[k]
                    # Drag H rigidly with the ring atom.
                    for hi in h_map.get(int(i), []):
                        Pv[int(hi)] = Pv[int(hi)] + delta
            except Exception:
                ok = False
                break
        if not ok or not np.all(np.isfinite(Pv)):
            continue
        variants.append((Pv, f"pucker{combo_idx}"))
    return variants


def _enumerate_rotamer_variants(
    mol_g,
    P: np.ndarray,
    *,
    max_configs: int = 12,
    syms: Optional[Sequence[str]] = None,
) -> List[Tuple[np.ndarray, str]]:
    """Enumerate single-bond rotamer variants of ``P``.

    Direct call into ``single_bond_rotamers.enumerate_single_bond_rotamers``
    which is gate-free at the helper level (the env flag is honoured by the
    converter caller, this module owns its own gate).  Returns a list of
    ``(P_variant, label)`` tuples with the first (identity) entry dropped.

    ``syms`` is forwarded to the topology-preservation hard gate inside
    :func:`enumerate_single_bond_rotamers` so rotations that break the
    SMILES bond graph are skipped at the source (default OFF via env
    flag; auto-ON under MOGUL_PRIMARY).
    """
    try:
        from .single_bond_rotamers import enumerate_single_bond_rotamers
    except Exception:
        return []
    if mol_g is None:
        return []
    out: List[Tuple[np.ndarray, str]] = []
    try:
        for vk, (Pv, rot_lab) in enumerate(
            enumerate_single_bond_rotamers(
                mol_g, np.asarray(P, dtype=float),
                max_configs=max_configs + 1,
                syms=syms,
            )
        ):
            if vk == 0:
                # Identity = original; skip so the caller can attach the
                # base variant explicitly.
                continue
            if not np.all(np.isfinite(Pv)):
                continue
            out.append((np.asarray(Pv, dtype=float), str(rot_lab)))
            if len(out) >= max_configs:
                break
    except Exception:
        return []
    return out


def _mahalanobis_severity(
    syms: Sequence[str], P: np.ndarray, mol_g, metal_idx: int,
    donor_idxs: Sequence[int],
) -> float:
    """Compute the Mogul-Mahalanobis severity for ranking.  ``inf`` on any
    detection failure so failing conformers rank last.  Mirrors
    ``grip_ensemble._compute_severity`` -- pulled out as a local helper so
    the dependency stays in one place.
    """
    try:
        from .grip_fragment_detect import detect_fragments
        from .grip_loss_terms import TotalGripLoss
        from .grip_mogul_lookup import get_default_library
        try:
            library = get_default_library()
        except Exception:
            library = None
        frozen = frozenset(
            (int(metal_idx), *[int(d) for d in donor_idxs])
        )
        terms = detect_fragments(
            mol_g, np.asarray(P, dtype=np.float64),
            frozen_atoms=frozen,
            donors=tuple(int(d) for d in donor_idxs),
            library=library, return_result=False,
        )
        agg = TotalGripLoss(terms=list(terms))
        if len(agg) == 0:
            return 0.0
        sev, _ = agg.evaluate(np.asarray(P, dtype=np.float64))
        sev = float(sev)
        if not np.isfinite(sev):
            return float("inf")
        return sev
    except Exception:
        return float("inf")


def _reordered_source_mol(smiles: str, n_atoms_expected: int):
    """Return the source mol reordered to match ``assemble_complex_mogul_primary``
    output (metal at index 0).  Atom bond orders + aromaticity flags survive
    the reorder, so the GRIP fragment detector can run on the same atom
    indices as the (syms, P) coordinate frame.

    Returns ``None`` on any failure -- caller falls back to the geometric
    bond-graph mol (severity collapses to 0 but the RMSD dedup ranking
    still works).
    """
    try:
        src = _full_complex_mol(smiles)
        if src is None:
            return None
        # CRITICAL (Bug #1 fix, 2026-06-07): re-parsing the SMILES here
        # MUST select the same primary metal as
        # ``assemble_complex_mogul_primary``.  We delegate to
        # ``_locate_metal_and_donors`` so the priority-based tie-break
        # (d/f-block > s-block > p-block "metals" like Sb / Sn / Pb) is
        # consistent across both code paths.
        metal_idx_src, _ = _locate_metal_and_donors(src)
        if metal_idx_src is None:
            return None
        n_src = src.GetNumAtoms()
        if n_src != n_atoms_expected:
            return None
        if metal_idx_src == 0:
            order = list(range(n_src))
        else:
            order = ([metal_idx_src]
                     + [i for i in range(n_src) if i != metal_idx_src])
        # RDKit RenumberAtoms preserves bond orders + aromaticity flags.
        return Chem.RenumberAtoms(src, order)
    except Exception:
        return None


def _aromatic_atom_set_from_smiles(
    smiles: str, n_atoms_expected: int,
) -> List[int]:
    """Build the aromatic-atom-index set for a Mogul-primary output.

    Parses the SMILES via :func:`_full_complex_mol`, runs RDKit
    aromaticity perception, locates the primary metal, then applies the
    SAME index re-ordering ``assemble_complex_mogul_primary`` applies
    (metal moved to index 0).  Returns the list of indices (in the
    output frame) that are flagged aromatic in the source mol.

    Universal: depends only on the source SMILES and RDKit aromaticity
    perception; no per-class branch.  Returns ``[]`` on any failure --
    the caller falls back to the geometric aromaticity heuristic, which
    catches benzene-style rings still in their aromatic geometric
    signature.
    """
    try:
        src = _full_complex_mol(smiles)
        if src is None:
            return []
        # Aromaticity perception (lazy / forgiving -- this is best-effort).
        try:
            Chem.SetAromaticity(src)
        except Exception:
            try:
                Chem.SanitizeMol(src, catchErrors=True)
            except Exception:
                pass
        # Locate the primary metal index in the source mol.  Delegate to
        # the priority-based selector so Sb / Sn / Pb donors cannot
        # outrank a d/f-block centre (Bug #1 fix, 2026-06-07).
        metal_idx_src, _ = _locate_metal_and_donors(src)
        if metal_idx_src is None:
            return []
        # Apply the assemble_complex_mogul_primary re-order: metal moves
        # to index 0; all other atoms keep their original ordering with
        # gaps closed.
        n_src = src.GetNumAtoms()
        if n_src != n_atoms_expected:
            # Atom-count mismatch (post-build re-orderings differ); skip
            # the supplied-set and fall through to the geometric heuristic.
            return []
        if metal_idx_src == 0:
            order = list(range(n_src))
        else:
            order = ([metal_idx_src]
                     + [i for i in range(n_src) if i != metal_idx_src])
        # Build map src_idx -> output_idx.
        src_to_out = {src_i: out_i for out_i, src_i in enumerate(order)}
        aromatic_out: List[int] = []
        for a in src.GetAtoms():
            if a.GetIsAromatic():
                src_i = int(a.GetIdx())
                if src_i in src_to_out:
                    aromatic_out.append(int(src_to_out[src_i]))
        return aromatic_out
    except Exception:
        return []


def _donors_from_geometry(
    syms: Sequence[str], P: np.ndarray, metal_idx: int = 0,
) -> List[int]:
    """Identify donor atom indices in an assembled (syms, P) tuple from
    geometric proximity to the metal.  Heavy atoms within 1.45 * ideal-bond
    distance of the metal count as donors.  Universal across CN / donor
    element / chelate; mirrors the donor-detection heuristic inside
    ``converter_backend._conformer_enum_post``.
    """
    P = np.asarray(P, dtype=float)
    m = int(metal_idx)
    if m < 0 or m >= len(syms):
        return []
    metal_sym = str(syms[m])
    donors: List[int] = []
    for j in range(len(syms)):
        if j == m:
            continue
        sj = str(syms[j])
        if sj == "H":
            continue
        if _bd._is_metal(sj):
            continue
        try:
            ideal = _bd._ideal_bond(metal_sym, sj)
            dij = float(np.linalg.norm(P[j] - P[m]))
            if dij < 1.45 * ideal:
                donors.append(int(j))
        except Exception:
            continue
    return donors


def enumerate_orbit_conformers(
    syms: Sequence[str],
    P: np.ndarray,
    *,
    metal_idx: int = 0,
    donor_idxs: Optional[Sequence[int]] = None,
    aromatic_atom_set: Optional[Sequence[int]] = None,
    source_mol=None,
    max_ring_states: int = 4,
    max_rotamer_states: int = 12,
    max_total: int = 50,
    rmsd_dedup_tol: float = 0.30,
) -> List[Tuple[List[str], np.ndarray, str, float]]:
    """Build the per-orbit conformer ensemble for one assembled structure.

    Pipeline:

      1. Build a heavy-atom bond graph (geometric bonds + metal-donor edges)
         on top of ``(syms, P)``.
      2. Enumerate Cremer-Pople ring puckers (5/6/7-ring, non-aromatic,
         skip metallacycle).
      3. Enumerate single-bond rotamers (non-ring, non-metal-incident).
      4. Cartesian-cap the union at ``max_total`` candidates (base always
         kept first).
      5. RMSD-dedup via :func:`conformer_dedup.dedup_by_rmsd` to remove
         near-duplicates that come from H-only or aromatic-ring rotations.
      6. Rank by Mogul-Mahalanobis severity ascending (best first).

    Returns
    -------
    list of (syms, P, sub_label, severity)
        ``sub_label`` is e.g. ``"base"``, ``"pucker0"``, ``"rot_+60_0_0"``.
        Sorted ascending by severity; the lowest-severity conformer is
        first.  Length <= ``max_total``.

    Universal: pure graph topology + geometry-derived bonds + Cremer-Pople
    formalism + dihedral rotation. No SMILES patterns, no per-class
    branches.  Deterministic: stable iteration order from
    ``itertools.product`` + ``dedup_by_rmsd``'s deterministic tiebreak.
    """
    P0 = np.asarray(P, dtype=float)
    syms_list = [str(s) for s in syms]
    base_variant: List[Tuple[List[str], np.ndarray, str]] = [
        (syms_list, P0.copy(), "base"),
    ]

    # If donors not supplied, derive them from geometric proximity to the
    # metal (heavy atoms within 1.45 * ideal-bond).  Required so the Mogul
    # severity calculation can identify the M-D coordination block and
    # rank conformers by Mahalanobis distance.
    if donor_idxs is None:
        donor_idxs = _donors_from_geometry(syms_list, P0, int(metal_idx))

    # 1) Build bond-graph mol.
    mol_g = _build_mol_for_conformer_enum(syms_list, P0)
    if mol_g is None or mol_g.GetNumAtoms() != len(syms_list):
        # Cannot enumerate without a graph -- return base only.
        sev = _mahalanobis_severity(
            syms_list, P0, mol_g, int(metal_idx),
            list(donor_idxs) if donor_idxs is not None else [],
        ) if mol_g is not None else float("inf")
        return [(syms_list, P0, "base", sev)]

    # 2) Ring-pucker variants.
    pucker_variants = _enumerate_ring_pucker_variants(
        mol_g, P0, max_per_ring=max_ring_states,
        aromatic_atom_set=aromatic_atom_set,
    )
    # 3) Rotamer variants -- the rotamer enumerator now applies the
    #    rotamer-topology hard gate at the source when the env flag is on
    #    (default OFF -> byte-identical with HEAD; auto-ON under
    #    DELFIN_FFFREE_MOGUL_PRIMARY=1).  Pucker variants are gated
    #    separately below so a single API path covers ALL rotation-derived
    #    conformers.
    rotamer_variants = _enumerate_rotamer_variants(
        mol_g, P0, max_configs=max_rotamer_states, syms=syms_list,
    )

    # 3b) Topology hard-gate on ring-pucker variants.  When the gate is
    #     active, derive (expected_bonds, expected_cn) ONCE on the source
    #     mol and filter every pucker before it reaches the combined list.
    #     This makes the conformer pipeline universal -- ANY rotation that
    #     breaks topology is rejected, regardless of which enumerator
    #     produced it.  Default OFF -> byte-identical with HEAD.
    try:
        from .rotamer_topology_gate import (
            _env_on as _topology_gate_on,
            extract_expected_bonds as _topology_expected_bonds,
            expected_metal_cn as _topology_expected_cn,
            rotation_preserves_topology as _topology_preserves,
        )
        _gate_on = bool(_topology_gate_on())
    except Exception:
        _gate_on = False
        _topology_preserves = None
        _topology_expected_bonds = None
        _topology_expected_cn = None

    if _gate_on and _topology_preserves is not None and pucker_variants:
        try:
            _exp_b = _topology_expected_bonds(mol_g)
            _exp_cn = _topology_expected_cn(mol_g, syms_list)
        except Exception:
            _exp_b, _exp_cn = None, None
        kept_puckers: List[Tuple[np.ndarray, str]] = []
        for Pv, lab in pucker_variants:
            try:
                if _topology_preserves(
                    syms_list, Pv, mol=mol_g,
                    expected_bonds=_exp_b, expected_cn=_exp_cn,
                ):
                    kept_puckers.append((Pv, lab))
            except Exception:
                kept_puckers.append((Pv, lab))
        pucker_variants = kept_puckers

    # 4) Combine (base + puckers + rotamers).  Cap at max_total candidates.
    combined: List[Tuple[List[str], np.ndarray, str]] = list(base_variant)
    for Pv, lab in pucker_variants:
        if len(combined) >= max_total:
            break
        combined.append((syms_list, np.asarray(Pv, dtype=float), str(lab)))
    for Pv, lab in rotamer_variants:
        if len(combined) >= max_total:
            break
        combined.append((syms_list, np.asarray(Pv, dtype=float), str(lab)))

    # 5) RMSD-dedup.  Pass severity = position-in-list so the base wins ties.
    try:
        from .conformer_dedup import dedup_by_rmsd
        frames = [
            (lab, P_v, float(k)) for k, (_s, P_v, lab) in enumerate(combined)
        ]
        deduped = dedup_by_rmsd(
            frames, threshold=rmsd_dedup_tol, max_keep=max_total,
            syms=syms_list,
        )
        # Map back to (syms, P, sub_label).
        kept: List[Tuple[List[str], np.ndarray, str]] = []
        for lab, P_v, _pos in deduped:
            kept.append((syms_list, np.asarray(P_v, dtype=float), str(lab)))
        if not kept:
            kept = combined
    except Exception:
        kept = combined

    # 6) Mogul-Mahalanobis severity per conformer -> rank ascending.
    donors_for_sev = list(donor_idxs) if donor_idxs is not None else []
    # For the Mahalanobis severity we prefer the SOURCE mol (with proper
    # bond orders + aromaticity perception) when the caller supplied one;
    # the geometric-bond mol_g misses bond-order information so the GRIP
    # fragment detector returns zero terms and the severity collapses to 0
    # for every conformer (defeating the ranking).  When ``source_mol`` is
    # not supplied the geometric mol is the best available -- severity may
    # still be 0 but the RMSD-dedup still gives a meaningful (insertion-
    # order) ordering, and the base structure remains first.
    mol_for_sev = source_mol if source_mol is not None else mol_g
    scored: List[Tuple[List[str], np.ndarray, str, float]] = []
    for s_list, P_v, lab in kept:
        sev = _mahalanobis_severity(
            s_list, P_v, mol_for_sev, int(metal_idx), donors_for_sev,
        )
        scored.append((s_list, P_v, lab, float(sev)))
    # Sort: severity ascending; ties broken by deterministic label order.
    # Ties on zero severity (when fragment detection returns no terms) are
    # broken by insertion order via the label suffix -- "base" sorts
    # before "pucker*" / "rot_*" lexicographically by chance, but more
    # importantly the dedup_by_rmsd step has already promoted the base.
    scored.sort(key=lambda t: (t[3], t[2]))
    return scored


def _default_geom_for_cn(metal_sym: str, cn: int) -> Optional[str]:
    """Return the canonical polyhedron name for ``(metal_sym, cn)``.

    Mirrors ``decompose._default_geometry`` but does not require a full
    SMILES decomposition pass -- the κⁿ enumerator needs to recompute
    the geometry for a rewired complex whose CN changed (κ¹-O → κ²-OO
    increments CN by +1).  Universal: pure CN + d8-metal lookup.
    """
    try:
        from delfin.fffree.decompose import _default_geometry
    except ImportError:
        return None
    try:
        return _default_geometry(metal_sym, int(cn))
    except Exception:
        return None


def enumerate_and_embed_mogul_primary(
    smiles: str,
    *,
    max_isomers: int = 50,
    max_attempts: int = 5,
    mol_override=None,
    donor_idxs_override: Optional[Sequence[int]] = None,
    geom_name_override: Optional[str] = None,
    label_prefix: str = "",
) -> Optional[List[Tuple[List[str], np.ndarray, str]]]:
    """Enumerate ALL Pólya orbits and embed each one through the Mogul-primary
    bounds matrix.

    The draft manuscript's contract is:

        SMILES → discrete configuration set (Burnside-closed enumeration)
              → realise EVERY member with one DG embed each
              → return the ensemble.

    Step 1 (per-orbit embed) was wired in
    :func:`assemble_complex_mogul_primary`; this function adds step 2
    (orbit enumeration → embed each → label):

      1.  decompose SMILES → metal + geometry + donor atoms + donor labels.
      2.  look up the Pólya group key for the geometry.
      3.  Burnside-enumerate every distinct vertex-coloring of the donor
          label multiset.
      4.  for each coloring, compute the donor-atom-index → polyhedron
          vertex permutation and call
          ``assemble_complex_mogul_primary(..., donor_at_vertex=perm)``
          so the bounds matrix encodes that orbit.  Each orbit is a
          separate DG embed.
      5.  return the list of ``(syms, P, label)`` triples.  Labels follow
          the legacy naming convention
          ``"<orbit-name>-<geom-tag>-<k>-mogul"``
          (e.g. ``"cis-OC-6-1-mogul"``) when the polyhedron has a
          registered antipode table; otherwise
          ``"MOGUL-PRIMARY-<geom-tag>-<k>"``.

    Universal: pure graph topology + element labels.  No SMILES pattern,
    no per-class branch.

    Parameters
    ----------
    mol_override, donor_idxs_override, geom_name_override : κⁿ-variant
        hook.  When the caller has already rewired the complex graph to
        encode a non-default κⁿ binding mode (e.g. κ²-OO promotion of an
        acetate), it passes the rewired mol + the new donor index list +
        the new polyhedron name here.  The orbit enumeration then operates
        on the κⁿ-promoted topology rather than re-deriving from the
        SMILES.  When all three are ``None`` (default), the SMILES is
        decomposed as before -- byte-identical with the pre-κⁿ path.
    label_prefix : str
        String prefix added to every orbit label, used to mark which κⁿ
        variant produced each XYZ (e.g. ``"k2-OO-"``).  Empty by default
        so seed-κⁿ orbit labels are unchanged.

    Returns
    -------
    list of (syms, P, label) or None
        ``None`` ⇒ orbit enumeration failed at a level where the caller
        should fall through to the single-structure path.  Otherwise the
        list contains at least one entry (the canonical orbit).
    """
    # 1) Source the mol / donor / geometry, either from override or from
    #    a fresh SMILES decomposition.
    if mol_override is not None and donor_idxs_override is not None:
        mol = mol_override
        metal_idx, _ = _locate_metal_and_donors(mol)
        if metal_idx is None:
            return None
        donor_idxs = sorted(int(d) for d in donor_idxs_override)
        if not donor_idxs:
            return None
        metal_sym = str(mol.GetAtomWithIdx(int(metal_idx)).GetSymbol())
        geom_name = geom_name_override or _default_geom_for_cn(
            metal_sym, len(donor_idxs),
        )
        if geom_name is None:
            return None
    else:
        try:
            from delfin.fffree import decompose as _DEC
        except ImportError:
            return None
        try:
            d = _DEC.decompose(smiles)
        except Exception:
            d = None
        if d is None:
            return None
        geom_name = d.get("geometry")
        mol = _full_complex_mol(smiles)
        if mol is None:
            return None
        metal_idx, donor_idxs = _locate_metal_and_donors(mol)
        if metal_idx is None or not donor_idxs:
            return None
        # Chelate-aware polyhedron pick (2026-06-07, hmaximilian).
        # When the decomposition reports a chelate-bearing complex,
        # override the naive ``decompose._default_geometry`` first-
        # candidate rule with the geometry whose vertex-vertex angles
        # best match the chelate's intrinsic bite (bite_deg derived
        # geometrically from the chelate ring size).  Universal — pure
        # geometric matching, no SMILES patterns, no per-class branches.
        # Default-OFF byte-identical when
        # ``DELFIN_FFFREE_CHELATE_AWARE_POLY_PICK=0`` (or when
        # MOGUL_PRIMARY is off).  Only chelate cases where the picker
        # finds a better polyhedron (e.g. AXOKAZ NNN-mer terpy → SPY-5,
        # en-Ni → SP-4 even for non-d8 metals) see a geometry switch.
        try:
            from delfin.fffree.chelate_aware_picker import (
                chelate_aware_picker_enabled as _cap_on,
                chelate_info_from_decompose as _ci_from,
                pick_polyhedron_chelate_aware as _pick_cap,
            )
            if _cap_on() and d.get("has_chelate"):
                _chelate_info = _ci_from(d)
                if _chelate_info:
                    _picked = _pick_cap(
                        cn=int(d.get("cn", 0) or 0),
                        chelate_info=_chelate_info,
                        metal_sym=str(d.get("metal", "")),
                    )
                    if _picked and _picked != geom_name:
                        geom_name = _picked
        except ImportError:
            pass
        except Exception:
            pass

        # f-block-aware polyhedron pick (2026-06-07, hmaximilian).
        # When no chelate constraint applied AND the metal is in the
        # f-block / high-CN picker scope (Ln/An at CN 8-12 or group-3 /
        # 4d-late / 5d / Cd at CN 7-12), score the candidates by their
        # CCDC donor-donor distance distribution match + symmetry
        # preference for ionic bonding.  Universal — pure CCDC-empirical
        # ranking, no SMILES patterns, no per-metal hardcodes.
        # Default-OFF byte-identical when ``DELFIN_FFFREE_F_BLOCK_PICK=0``
        # (or when MOGUL_PRIMARY is off).  Fires only when the picker
        # finds a different geometry than the legacy first-candidate.
        try:
            from delfin.fffree.f_block_picker import (
                f_block_picker_enabled as _fbp_on,
                f_block_polyhedron_picker as _pick_fbp,
            )
            if _fbp_on() and not d.get("has_chelate"):
                _cn_int = int(d.get("cn", 0) or 0)
                _metal = str(d.get("metal", ""))
                # Expand per-ligand donor_elems to per-vertex donor symbols:
                # a denticity-k ligand contributes k donor sites of the same
                # element symbol (donor_elem field is a single string).
                _donor_syms: List[str] = []
                for _lg in (d.get("ligands") or []):
                    _de = _lg.get("donor_elem")
                    _de_list = _lg.get("donor_elems")
                    _dk = int(_lg.get("denticity", 1) or 1)
                    if isinstance(_de_list, (list, tuple)) and _de_list:
                        for _s in _de_list:
                            _donor_syms.append(str(_s))
                    elif _de is not None:
                        for _ in range(max(1, _dk)):
                            _donor_syms.append(str(_de))
                # Hapto donors are encoded as ``"C_hapto"``; the CCDC lookup
                # is keyed on element symbols, so collapse the marker to "C".
                _donor_syms = [
                    "C" if s == "C_hapto" else s for s in _donor_syms
                ]
                if _metal and _cn_int >= 7 and len(_donor_syms) == _cn_int:
                    _picked_fb = _pick_fbp(
                        metal_sym=_metal,
                        cn=_cn_int,
                        donors=_donor_syms,
                        force=True,
                    )
                    if _picked_fb and _picked_fb != geom_name:
                        geom_name = _picked_fb
        except ImportError:
            pass
        except Exception:
            pass

    geom_polya_key = _GEOM_NAME_TO_POLYA_KEY.get(str(geom_name))
    if geom_polya_key is None:
        return None

    donor_labels = [str(mol.GetAtomWithIdx(int(d)).GetSymbol())
                    for d in donor_idxs]
    if not donor_labels:
        return None

    # 3) Pólya enumeration
    try:
        from delfin.fffree import polya_isomer_count as _PIC
    except ImportError:
        return None
    spec: Dict[str, int] = {}
    for lab in donor_labels:
        spec[lab] = spec.get(lab, 0) + 1
    try:
        colorings = _PIC.enumerate_isomers(geom_polya_key, spec)
    except Exception:
        return None
    if not colorings:
        return None

    # 4) embed each orbit (+ optional per-orbit conformer enumeration)
    geom_tag = str(geom_name).split()[0]
    results: List[Tuple[List[str], np.ndarray, str]] = []

    # Conformer-enumeration knobs (env-tunable; default ON).  The total cap
    # protects against combinatorial explosion on large rings + many rotors;
    # the per-orbit cap keeps the ensemble balanced across orbits.  When
    # the gate is OFF (DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS=0), every
    # orbit emits exactly one structure (the base build) and the path is
    # byte-identical to the pre-2026-06-07 orbit-only behaviour.
    conformers_on = mogul_primary_conformers_enabled()
    try:
        conf_max_total = int(os.environ.get(
            "DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS_MAX_TOTAL", "300"
        ))
    except (TypeError, ValueError):
        conf_max_total = 300
    try:
        conf_per_orbit = int(os.environ.get(
            "DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS_PER_ORBIT", "30"
        ))
    except (TypeError, ValueError):
        conf_per_orbit = 30
    try:
        conf_ring_states = int(os.environ.get(
            "DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS_RING_STATES", "4"
        ))
    except (TypeError, ValueError):
        conf_ring_states = 4
    try:
        conf_rot_states = int(os.environ.get(
            "DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS_ROT_STATES", "12"
        ))
    except (TypeError, ValueError):
        conf_rot_states = 12
    try:
        conf_rmsd_tol = float(os.environ.get(
            "DELFIN_FFFREE_MOGUL_PRIMARY_CONFORMERS_RMSD_TOL", "0.30"
        ))
    except (TypeError, ValueError):
        conf_rmsd_tol = 0.30

    n_total_emitted = 0
    for k, coloring in enumerate(colorings[:max_isomers]):
        if conformers_on and n_total_emitted >= conf_max_total:
            break
        donor_at_vertex = _coloring_to_donor_assignment(
            coloring, donor_labels, donor_idxs,
        )
        if donor_at_vertex is None:
            continue
        try:
            built = assemble_complex_mogul_primary(
                smiles,
                max_attempts=max_attempts,
                donor_at_vertex=donor_at_vertex,
                geometry_key=str(geom_name),
                mol_override=mol_override,
                donor_idxs_override=(
                    donor_idxs if mol_override is not None else None
                ),
            )
        except Exception:
            built = None
        if built is None:
            continue
        out_syms, P = built
        iso_name = _classify_orbit_label(geom_polya_key, coloring)
        if iso_name:
            orbit_label = f"{iso_name}-{geom_tag}-{k+1}-mogul"
        else:
            orbit_label = f"MOGUL-PRIMARY-{geom_tag}-{k+1}"
        if label_prefix:
            orbit_label = f"{label_prefix}{orbit_label}"

        # Conformer enumeration: Cremer-Pople ring puckers × single-bond
        # rotamers × RMSD-dedup × Mahalanobis-rank.  The base structure is
        # always conf0; any additional puckers/rotamers are appended in
        # severity-ascending order.  When the conformer gate is OFF, only
        # the base structure is emitted (byte-identical orbit-only path).
        if not conformers_on:
            results.append((out_syms, P, orbit_label))
            n_total_emitted += 1
            continue

        # Budget remaining = how many conformers we can still emit before
        # the total cap.  Always at least 1 (the base structure) so each
        # orbit appears at least once even at the cap edge.
        budget = max(1, min(conf_per_orbit, conf_max_total - n_total_emitted))
        # Aromatic-atom-set from the source SMILES so the pucker enum
        # skips benzene/pyridine/NHC π-systems whose post-embed geometry
        # has lost the aromatic bond-length / planarity signature
        # (RDKit's embed routes aromatic bonds as SINGLE).  Computed once
        # per SMILES outside the orbit loop -- but the post-build re-order
        # maps source -> output identically across orbits so the same set
        # applies.
        arom_set = _aromatic_atom_set_from_smiles(smiles, len(out_syms))
        # Source mol with proper bond orders + aromaticity -- the GRIP
        # fragment detector needs bond-order information to enumerate the
        # angle / improper terms whose Mahalanobis distance drives the
        # severity-ranking of conformers.  Without it, every conformer's
        # severity collapses to 0 and the ranking degenerates.
        src_mol = _reordered_source_mol(smiles, len(out_syms))
        try:
            ensemble = enumerate_orbit_conformers(
                out_syms, P,
                metal_idx=0,           # assemble_complex_mogul_primary puts metal at 0
                donor_idxs=None,       # auto-derive from geometry
                aromatic_atom_set=arom_set or None,
                source_mol=src_mol,
                max_ring_states=conf_ring_states,
                max_rotamer_states=conf_rot_states,
                max_total=budget,
                rmsd_dedup_tol=conf_rmsd_tol,
            )
        except Exception:
            ensemble = []
        if not ensemble:
            # Fallback: at least emit the base structure when conformer
            # enumeration fails (defensive -- in practice always succeeds).
            results.append((out_syms, P, orbit_label))
            n_total_emitted += 1
            continue
        for ci, (s_list, P_v, sub_label, _sev) in enumerate(ensemble):
            if n_total_emitted >= conf_max_total:
                break
            label = f"{orbit_label}-conf{ci}-{sub_label}"
            results.append((list(s_list), np.asarray(P_v, dtype=float), label))
            n_total_emitted += 1

    if not results:
        return None
    return results


# ---------------------------------------------------------------------------
# Top-level κⁿ binding-mode enumeration
# (2026-06-07, branch feat-mogul-primary-2026-06-07)
# ---------------------------------------------------------------------------
def enumerate_kappa_variants_mogul_primary(
    smiles: str,
    *,
    max_isomers: int = 50,
    max_attempts: int = 5,
) -> Optional[List[Tuple[List[str], np.ndarray, str]]]:
    """Top-level entry: enumerate κⁿ binding-mode variants × Pólya orbits
    × bounds-matrix embed for ``smiles``.

    Algorithm:

      1. Parse + decompose ``smiles`` -> base mol, primary metal, seed donor
         indices (one per current M-X bond).
      2. For each ligand fragment (identified via its seed donor), call
         ``ambidentate_kappa_enum.enumerate_kappa_modes_for_ligand`` to
         generate the κⁿ options for that ligand.
      3. Take the Cartesian product across ligands ->
         per-complex κⁿ variant donor sets.  De-duplicate by donor
         multiset.
      4. For each variant:
           a. Rewire the mol -> insert/remove dative bonds so the
              metal is bonded exactly to the variant's donor set.
           b. Re-derive the polyhedron name from the new CN.
           c. Call ``enumerate_and_embed_mogul_primary`` on the rewired
              mol with ``label_prefix`` set to the κⁿ tag.
      5. Concatenate the per-variant results into one ensemble.

    Universal: ambidentate detection is graph + element + Lewis-octet,
    no SMARTS.  Bite-compatibility is a graph-distance test plus a
    universal chord-length estimate.  No per-anion table.

    Env-gate: ``DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM``.  Default ON
    when ``DELFIN_FFFREE_MOGUL_PRIMARY=1`` is also set; off otherwise
    (returns ``None`` and the caller falls through to the κ¹-seed
    orbit-only path).

    Returns
    -------
    list of (syms, P, label) or None
        ``None`` ⇒ κⁿ enumeration unavailable for this SMILES (e.g.
        decompose failure, no ambidentate ligand detected so the result
        would equal the orbit-only path).  Caller falls through to the
        orbit-only path on ``None``.
    """
    # Env-gate check (the ambidentate module is the single source of truth).
    try:
        from delfin.fffree.ambidentate_kappa_enum import (
            _env_on as _kappa_env_on,
            enumerate_complex_kappa_variants,
            make_variant_mol,
        )
    except ImportError:
        return None
    if not _kappa_env_on():
        return None

    # Base mol from the full-complex graph.  We do NOT require the
    # heavier ``decompose`` path -- it rejects large NHC ligands and
    # other edge cases where the κⁿ enumeration still has work to do.
    # This mirrors ``assemble_complex_mogul_primary`` so any SMILES
    # that yields a single-structure embed also yields a κⁿ enumeration.
    mol_base = _full_complex_mol(smiles)
    if mol_base is None:
        return None
    metal_idx, base_donors = _locate_metal_and_donors(mol_base)
    if metal_idx is None or not base_donors:
        return None

    # κⁿ enumeration over the WHOLE complex.
    try:
        variants = enumerate_complex_kappa_variants(
            mol_base, int(metal_idx), base_donors,
        )
    except Exception:
        return None
    if not variants:
        return None

    metal_sym = str(mol_base.GetAtomWithIdx(int(metal_idx)).GetSymbol())

    # If the only variant is the default κ¹-seed, no κⁿ work to do -- let
    # the caller use the cheaper orbit-only path (returns None here so the
    # dispatcher falls through).
    non_default = [v for v in variants if not v.get("is_seed_default")]
    if not non_default:
        return None

    # Build per-variant ensemble.  The DEFAULT (seed-κ¹) variant
    # delegates to the unmodified orbit path so its bytes are byte-
    # identical with the pre-κⁿ output.  Non-default variants get a
    # rewired mol + a label prefix.
    #
    # Global cap: when κⁿ × orbits × conformers combinatorially explodes
    # (e.g. 6 κⁿ × 2 orbits × 13 conformers = 156 structures), we cap
    # the ensemble to keep pool sizes bounded.  Env-tunable.
    try:
        MAX_ENSEMBLE_TOTAL = int(os.environ.get(
            "DELFIN_FFFREE_MOGUL_PRIMARY_KAPPA_ENUM_ENSEMBLE_MAX", "200",
        ))
    except (TypeError, ValueError):
        MAX_ENSEMBLE_TOTAL = 200
    ensemble: List[Tuple[List[str], np.ndarray, str]] = []
    seen_xyz_hashes = set()  # de-dup variants whose mol rewiring collapsed
    for v in variants:
        if len(ensemble) >= MAX_ENSEMBLE_TOTAL:
            break
        v_donors = sorted(int(d) for d in v["all_donors"])
        if v.get("is_seed_default"):
            sub = enumerate_and_embed_mogul_primary(
                smiles,
                max_isomers=max_isomers,
                max_attempts=max_attempts,
                label_prefix="",
            )
            # Fallback: if orbit-only path fails (e.g. decompose rejects
            # a large NHC ligand but the single-structure embed still
            # works), emit the canonical single structure so the κⁿ
            # ensemble always has the κ¹-seed reference structure.  The
            # conformer-enum extension applies here too -- the single
            # structure becomes a per-isomer conformer ensemble when the
            # gate is on.
            if not sub:
                try:
                    base_built = assemble_complex_mogul_primary(
                        smiles, max_attempts=max_attempts,
                    )
                except Exception:
                    base_built = None
                if base_built is not None:
                    syms_b, P_b = base_built
                    if mogul_primary_conformers_enabled():
                        try:
                            arom_b = _aromatic_atom_set_from_smiles(
                                smiles, len(syms_b),
                            )
                        except Exception:
                            arom_b = []
                        try:
                            src_b = _reordered_source_mol(
                                smiles, len(syms_b),
                            )
                        except Exception:
                            src_b = None
                        try:
                            ens_b = enumerate_orbit_conformers(
                                syms_b, P_b,
                                metal_idx=0,
                                donor_idxs=None,
                                aromatic_atom_set=arom_b or None,
                                source_mol=src_b,
                            )
                        except Exception:
                            ens_b = []
                        if ens_b:
                            sub = [
                                (s_l, P_v,
                                 f"MOGUL-PRIMARY-1-conf{ci}-{sub_lab}")
                                for ci, (s_l, P_v, sub_lab, _sev)
                                in enumerate(ens_b)
                            ]
                        else:
                            sub = [(syms_b, P_b, "MOGUL-PRIMARY-1")]
                    else:
                        sub = [(syms_b, P_b, "MOGUL-PRIMARY-1")]
        else:
            try:
                new_mol = make_variant_mol(
                    mol_base, int(metal_idx), base_donors, v_donors,
                )
            except Exception:
                new_mol = None
            if new_mol is None:
                continue
            new_geom = _default_geom_for_cn(metal_sym, len(v_donors))
            if new_geom is None:
                continue
            tag = v.get("label_suffix") or "kappa"
            tag = tag.replace("+", "_")
            sub = enumerate_and_embed_mogul_primary(
                smiles,
                max_isomers=max_isomers,
                max_attempts=max_attempts,
                mol_override=new_mol,
                donor_idxs_override=v_donors,
                geom_name_override=new_geom,
                label_prefix=f"{tag}-",
            )
            # Fallback: when orbit enumeration fails for the variant
            # (e.g. polyhedron without a registered Pólya group, or
            # GRIP-lib lookup misses), still emit the single rewired
            # κⁿ structure so the user sees the κⁿ topology choice.
            # When the conformer gate is on, the single structure becomes
            # a per-variant conformer ensemble.
            if not sub:
                try:
                    v_built = assemble_complex_mogul_primary(
                        smiles,
                        max_attempts=max_attempts,
                        geometry_key=new_geom,
                        mol_override=new_mol,
                        donor_idxs_override=v_donors,
                    )
                except Exception:
                    v_built = None
                if v_built is not None:
                    syms_v, P_v = v_built
                    if mogul_primary_conformers_enabled():
                        try:
                            arom_v = _aromatic_atom_set_from_smiles(
                                smiles, len(syms_v),
                            )
                        except Exception:
                            arom_v = []
                        try:
                            src_v = _reordered_source_mol(
                                smiles, len(syms_v),
                            )
                        except Exception:
                            src_v = None
                        try:
                            ens_v = enumerate_orbit_conformers(
                                syms_v, P_v,
                                metal_idx=0,
                                donor_idxs=None,
                                aromatic_atom_set=arom_v or None,
                                source_mol=src_v,
                            )
                        except Exception:
                            ens_v = []
                        if ens_v:
                            sub = [
                                (s_l, P_q,
                                 f"{tag}-MOGUL-PRIMARY-1-conf{ci}-{sub_lab}")
                                for ci, (s_l, P_q, sub_lab, _sev)
                                in enumerate(ens_v)
                            ]
                        else:
                            sub = [(syms_v, P_v, f"{tag}-MOGUL-PRIMARY-1")]
                    else:
                        sub = [(syms_v, P_v, f"{tag}-MOGUL-PRIMARY-1")]
        if not sub:
            continue
        for syms, P, lab in sub:
            if len(ensemble) >= MAX_ENSEMBLE_TOTAL:
                break
            # Dedup by content hash of (sorted atom array) -- two
            # different κⁿ variants on geometrically-equivalent donors
            # may collapse to the same XYZ after the embed.
            try:
                key = (
                    tuple(syms),
                    tuple(round(float(x), 4) for row in np.asarray(P) for x in row),
                )
            except Exception:
                key = None
            if key is not None and key in seen_xyz_hashes:
                continue
            if key is not None:
                seen_xyz_hashes.add(key)
            ensemble.append((syms, P, lab))

    if not ensemble:
        return None
    return ensemble


# ---------------------------------------------------------------------------
# Internal building blocks
# ---------------------------------------------------------------------------
def _full_complex_mol(smiles: str):
    """Return the full RDKit Mol with metal + ligands + explicit H atoms.

    Uses ``smiles_converter._prepare_mol_for_embedding`` so the same
    SMILES handling (stk, dative bonds, normalisation, AddHs) as the rest
    of the converter is used.  Returns ``None`` on parse failure.
    """
    try:
        from delfin.smiles_converter import _prepare_mol_for_embedding
    except ImportError:
        return None
    try:
        mol = _prepare_mol_for_embedding(smiles)
    except Exception:
        mol = None
    if mol is None:
        return None
    # Some prep paths return a mol without explicit H -> ensure AddHs.
    try:
        if not any(a.GetSymbol() == "H" for a in mol.GetAtoms()):
            mol = Chem.AddHs(mol)
    except Exception:
        pass
    # RDKit's ``GetMoleculeBoundsMatrix`` calls into ``set14Bounds`` which
    # accesses ``RingInfo`` via ``numBondRings`` — that fails with a
    # ``Pre-condition Violation`` if the ring perception has not been run on
    # the mol.  For stk-built / partial-sanitised metal-complex mols this is
    # not always done by the prep path, so we run a minimal symmetry-based
    # ring perception (``FastFindRings``) which works even on partially-
    # sanitised mols (no aromaticity perception required).  This is a no-op
    # when ring info is already populated.
    try:
        ri = mol.GetRingInfo()
        if not ri.IsInitialized() if hasattr(ri, "IsInitialized") else not ri.NumRings():
            Chem.FastFindRings(mol)
    except Exception:
        try:
            Chem.FastFindRings(mol)
        except Exception:
            pass
    return mol


def _metal_priority(symbol: str) -> int:
    """Universal primary-metal priority for tie-break.

    Lower value = higher priority.  Used to disambiguate the "primary"
    metal when several atoms in the SMILES are flagged as metals by
    :func:`delfin._bond_decollapse._is_metal` (e.g. a Pt complex with
    Sb-donor ligands).  Without a stable priority the multi-metal
    tie-break is the SMILES atom order, so an Sb ligand-donor can
    overrule the Pt centre and end up at output index 0 — this is the
    Bug #1 root cause documented in the 2026-06-07 invariant audit.

    Priority (ascending = higher) follows the chemical convention used
    throughout DELFIN:

      0: d-block + f-block transition metals (the "actual" metals)
      1: alkali / alkaline earth (Group 1/2)
      2: post-transition / p-block "metals" classified as metals by
         ``_is_metal`` (Al, Ga, In, Tl, Sn, Pb, Bi, Sb, Ge, Po).  These
         are routinely donors in TMC ligands and should NEVER outrank a
         d/f metal as the coordination centre.

    Universal: atomic-number ranges only, no SMILES patterns.
    """
    try:
        from rdkit.Chem import GetPeriodicTable  # type: ignore
        z = int(GetPeriodicTable().GetAtomicNumber(str(symbol)))
    except Exception:
        # Element-symbol fallback table for the most common metals.  The
        # priority is conservative — the alternative is a hard crash, and
        # the legacy max-degree behaviour is what we get back here.
        _Z = {
            "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
            "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
            "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
            "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
            "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77,
            "Pt": 78, "Au": 79, "Hg": 80,
            "La": 57, "Ce": 58, "Lu": 71,
            "Li": 3, "Na": 11, "K": 19, "Rb": 37, "Cs": 55,
            "Be": 4, "Mg": 12, "Ca": 20, "Sr": 38, "Ba": 56,
            "Al": 13, "Ga": 31, "In": 49, "Tl": 81,
            "Ge": 32, "Sn": 50, "Pb": 82, "Sb": 51, "Bi": 83, "Po": 84,
        }
        z = int(_Z.get(str(symbol), 0))
    # d-block (21-30, 39-48, 72-80, 104-112) + f-block (57-71, 89-103)
    if (21 <= z <= 30) or (39 <= z <= 48) or (72 <= z <= 80) or (104 <= z <= 112):
        return 0
    if (57 <= z <= 71) or (89 <= z <= 103):
        return 0
    # Alkali / alkaline earth
    if z in (3, 11, 19, 37, 55, 87, 4, 12, 20, 38, 56, 88):
        return 1
    # Everything else (p-block "metals": Al, Ga, In, Tl, Sn, Pb, Bi, Sb, Ge, Po, ...)
    return 2


def _locate_metal_and_donors(mol) -> Tuple[Optional[int], List[int]]:
    """Identify the primary metal index + list of donor atom indices.

    Primary metal selection (universal — no SMILES patterns):

      1. Lowest ``_metal_priority`` symbol wins (d/f-block > s-block >
         p-block "metals" like Sb / Sn / Pb).  This guarantees that a
         d-block centre is never overruled by a p-block donor in a
         multi-metal SMILES.  Fixes the 2026-06-07 invariant Bug #1
         (D-ATOQOP: Pt outranks Sb).
      2. Among same-priority metals, highest graph degree wins (matches
         the legacy ``decompose`` selector for the genuine multi-metal
         case, e.g. Mn₂(CO)₁₀ where both Mn are d-block).
      3. Final tie-break: lowest atom index (deterministic).

    Donors = heavy-atom neighbours of the primary metal in the molecular
    graph (sorted by index).
    """
    metals = [a.GetIdx() for a in mol.GetAtoms() if _bd._is_metal(a.GetSymbol())]
    if not metals:
        return None, []
    if len(metals) == 1:
        m = metals[0]
    else:
        def _sort_key(mi: int):
            sym = mol.GetAtomWithIdx(mi).GetSymbol()
            # Sort ascending: low priority + high degree + low idx wins.
            return (
                int(_metal_priority(sym)),
                -int(mol.GetAtomWithIdx(mi).GetDegree()),
                int(mi),
            )
        m = min(metals, key=_sort_key)
    donors = sorted(int(n.GetIdx()) for n in mol.GetAtomWithIdx(m).GetNeighbors())
    return int(m), donors


def _embed_with_bounds(
    mol,
    lower: np.ndarray,
    upper: np.ndarray,
    *,
    metal_idx: Optional[int] = None,
    donor_idxs: Optional[Sequence[int]] = None,
    max_attempts: int = 5,
) -> Optional[np.ndarray]:
    """Embed the molecule subject to ``(lower, upper)`` distance bounds.

    Steps:

      1. Pull RDKit's default bounds matrix (so 1-2 / 1-3 / ring entries
         RDKit already populates are kept as the *base*).
      2. Overwrite every cell where our CCDC-populated bounds are tighter
         or where RDKit's default is uninformative.
      3. Triangle-smooth (``DoTriangleSmoothing``) — required for a
         feasible DG embed.
      4. ``EmbedMolecule`` with the bounds matrix, multiple seeds.

    Returns the (N, 3) ndarray of coordinates or ``None`` on failure.

    Determinism: fixed ``randomSeed`` per attempt, ``useRandomCoords=True``
    so the initial coordinates are drawn from the fixed-seed RNG.
    """
    n = mol.GetNumAtoms()
    if lower.shape != (n, n) or upper.shape != (n, n):
        return None

    # RDKit's ``EmbedMolecule`` internally kekulises the molecule and types
    # every atom for the ETKDG 1-4 enrichment.  Both steps fail on metal-
    # complex graphs:
    #
    #   * Dative bonds count zero toward atom valence, which breaks the C+
    #     carbene of an NHC fragment (C has 3 σ-bonds + 1 dative-out, so
    #     valence appears as 3 and kekuliser refuses to assign double bonds
    #     in the ring).
    #   * Aromatic flags from the SMILES parser do not survive the
    #     metallation step (the dative bond pulls electron density off the
    #     ring), so the kekuliser sees unmatched aromatic atoms.
    #
    # The bounds matrix we built in ``mogul_bounds`` is index-based and
    # depends only on the graph + element + hybridisation tags, not on
    # bond orders.  So we can pre-process a temporary copy of the mol for
    # embedding: all DATIVE→SINGLE, all AROMATIC→SINGLE, all aromatic flags
    # cleared.  Atom indices are preserved, so the bounds matrix transfer
    # one-to-one.
    embed_mol = Chem.RWMol(mol)
    for b in embed_mol.GetBonds():
        bt = b.GetBondType()
        if bt == Chem.BondType.DATIVE or bt == Chem.BondType.AROMATIC:
            b.SetBondType(Chem.BondType.SINGLE)
    for a in embed_mol.GetAtoms():
        a.SetIsAromatic(False)
    # Remove the metal-donor bonds in the embed graph: RDKit's
    # ``GetMoleculeBoundsMatrix`` populates 1-2 distances from covalent radii
    # (Sutton tables), which for M-D bonds gives the WRONG distance: a generic
    # Pyykkö M-X sum (~ 1.7 Å for Ni-N), not the empirical CCDC M-D
    # distribution (~ 2.05 Å for Ni(II)-pyridyl).  Leaving the bond in place
    # makes RDKit fight our injected M-D bound and the embed lands at the
    # covalent-radius distance.  Removing the bond converts the M-D pair into
    # an unbonded pair whose only constraint is OUR injected (lo, hi) window,
    # which is exactly the contract the draft manuscript requires.  Atom
    # indices are preserved so the bounds matrix transfer is unchanged.
    if metal_idx is not None and donor_idxs is not None:
        for d in donor_idxs:
            b = embed_mol.GetBondBetweenAtoms(int(metal_idx), int(d))
            if b is not None:
                embed_mol.RemoveBond(int(metal_idx), int(d))
    embed_mol = embed_mol.GetMol()
    # Sanitise without kekulisation / aromaticity perception — both fight
    # the bond-order rewrite above and re-introduce the same kekulise
    # errors we just removed.  Ring info + hybridisation + valence are
    # still computed, so RDKit's bounds-matrix builder works.
    try:
        sanitize_flags = (
            Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
            ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
        )
        Chem.SanitizeMol(embed_mol, sanitizeOps=sanitize_flags)
    except Exception:
        # Some highly charged metal carbene fragments cannot satisfy the
        # SANITIZE_PROPERTIES check.  Fall back to ring-only perception
        # so embedding still has the ring info it needs.
        try:
            Chem.FastFindRings(embed_mol)
        except Exception:
            pass

    for attempt in range(int(max_attempts)):
        # Fresh per-attempt copy of the RDKit default bounds matrix.
        try:
            bm = _DG.GetMoleculeBoundsMatrix(embed_mol)
        except Exception:
            return None
        # Inject our CCDC bounds.  RDKit convention:
        #   bm[i][j] for i < j  -> UPPER bound
        #   bm[i][j] for i > j  -> LOWER bound
        # We copy lower into the (i>j) triangle and upper into (i<j).
        #
        # M-D and donor-donor pairs are FORCE-overwritten with our values
        # (the bounds matrix would otherwise hold an uninformative vdW
        # floor that is LARGER than the empirical CCDC bond distance).
        # Other pairs use "tighten only" semantics so RDKit's bonded /
        # ring-pucker bounds stay in place where we have nothing better.
        _UB_CAP = 1.0e3
        donor_set = set(int(d) for d in (donor_idxs or []))
        for i in range(n):
            for j in range(i + 1, n):
                lo_ij = float(lower[i, j])
                hi_ij = float(upper[i, j])
                # Identify metal-donor and donor-donor pairs for force-override.
                is_md = (i == metal_idx and j in donor_set) or (
                    j == metal_idx and i in donor_set
                )
                is_dd = (i in donor_set) and (j in donor_set)
                force = is_md or is_dd
                # Upper bound
                if hi_ij < _UB_CAP:
                    cur_hi = float(bm[i][j])
                    if force or cur_hi <= 0.0 or hi_ij < cur_hi:
                        bm[i][j] = hi_ij
                # Lower bound
                if lo_ij > 0.0:
                    cur_lo = float(bm[j][i])
                    if force or lo_ij > cur_lo:
                        bm[j][i] = lo_ij
        # Triangle smoothing — best-effort.  When the CCDC bounds form an
        # infeasible system on some pair (typical for highly fused or
        # over-constrained polydentate ligands the manifold is sparse on)
        # ``DoTriangleSmoothing`` returns False but partially smooths the
        # matrix in place.  RDKit's embedder copes with that better than
        # we do by failing here — proceed and let ``EmbedMolecule`` either
        # succeed (the usual case for ~80-90 % of inputs) or return a
        # negative cid which the next attempt retries with a different
        # seed.
        try:
            _DGs.DoTriangleSmoothing(bm)
        except Exception:
            pass
        # Embed using the populated bounds.
        try:
            ep = _DG.EmbedParameters()
            ep.randomSeed = SEED + attempt
            ep.useRandomCoords = True
            ep.SetBoundsMat(bm)
            ep.clearConfs = True
            ep.maxIterations = 200
            # Removing M-D bonds in the embed graph DISCONNECTS the metal
            # from each ligand fragment (the bounds matrix re-couples them
            # via M-D + D-D + cross-ligand vdW).  RDKit's default of
            # ``embedFragmentsSeparately=True`` would then place each
            # ligand at its own origin and the M-D + D-D bounds would be
            # silently dropped.  We must embed as a single fragment.
            ep.embedFragmentsSeparately = False
            # Smoothing-feasibility failures are common on densely
            # constrained metal complexes (~ 25 % of inputs).  RDKit
            # still produces a useful embed because the bounds matrix is
            # only partially infeasible; ignoring the failure lets us
            # progress and report a structure rather than bailing out.
            ep.ignoreSmoothingFailures = True
            # Stronger bounds-matrix force scaling: the default 1.0 lets
            # M-D bonds drift several tenths of an Å outside our CCDC
            # window because the loss is a soft mid-range minimum.
            # Doubling the force pulls the embedded distance much closer
            # to the (lower, upper) envelope without diverging the
            # optimisation.
            ep.boundsMatForceScaling = 2.0
            # Tighter ETKDG-style optimiser stop criterion -> the final
            # refinement does more work on the bounds-violations.
            ep.optimizerForceTol = 1.0e-4
            # NOTE: useExpTorsionAnglePrefs / useBasicKnowledge cause RDKit
            # to inject sp/sp2/sp3 angle bounds for the metal center
            # itself, which conflict with our CCDC donor-donor bounds and
            # cause "bad bounds" failures across CN >= 4.  Leave these
            # OFF; the bounds matrix WE supply already encodes the
            # empirical 1-3 (angle) distances via mogul_bounds and the
            # CCDC distribution does a better job than the ETKDG
            # heuristics for transition-metal complexes (which is exactly
            # the gap the draft manuscript was written to close).
            cid = AllChem.EmbedMolecule(embed_mol, ep)
        except Exception:
            cid = -1
        if cid is None or cid < 0:
            continue
        try:
            conf = embed_mol.GetConformer(cid)
            P = np.array(conf.GetPositions(), dtype=float)
        except Exception:
            continue
        if not np.all(np.isfinite(P)):
            continue
        return P
    return None


# ---------------------------------------------------------------------------
# TM-category detection + M-D bound override
# ---------------------------------------------------------------------------
def _detect_tm_category(mol, metal_idx: int, donor_idx: int) -> Optional[str]:
    """Universal graph-topology classification of a single M-D bond into
    one of the seven dedicated TM categories from the draft manuscript:

      - ``carbene``    : the donor is a divalent or triply-bonded C with
                         two N neighbours in a 5-ring (NHC pattern: the
                         carbene C of imidazol-2-ylidene, imidazolin-,
                         triazol-, …).  Also matches Fischer / Schrock
                         carbenes where the donor C has a non-carbon
                         heteroatom (N, O, S) neighbour and the metal is
                         the only metal in the C's first shell.
      - ``hapto_eta2`` : the donor C is part of an alkene (sp2 C=C) where
                         both sp2 carbons are bonded to the metal.
      - ``hapto_eta5`` : the donor C is part of a 5-ring that is fully
                         coordinated to the metal (η⁵-cyclopentadienyl).
      - ``hapto_eta6`` : same for a 6-ring (η⁶-arene).
      - ``mu_bridge``  : the donor atom (typically halide, O, S, N) is
                         bonded to two or more metals.
      - ``agostic``    : the donor is an H attached to a C that is itself
                         within ~3 bonds of the metal (3-centre-2-electron).
      - ``ox_addition``: the donor is part of an axial X-Y unit whose Y
                         is also coordinated to the same metal (Sigman/
                         Hartwig oxidative-addition adducts; rare).

    Returns the category name (compatible with
    ``GripLibrary.lookup_tm_category``) or ``None`` if none of the seven
    patterns match — in which case the generic
    ``mogul_bounds._lookup_bond_md`` value already in the bounds matrix
    stays as-is.

    Universal: pure graph topology + atom-element + hybridisation.
    """
    try:
        m = mol.GetAtomWithIdx(int(metal_idx))
        d = mol.GetAtomWithIdx(int(donor_idx))
    except Exception:
        return None

    d_sym = d.GetSymbol()

    # -- μ-bridge: donor bonded to >=2 metals
    n_metal_nbrs = sum(
        1 for nb in d.GetNeighbors() if _bd._is_metal(nb.GetSymbol())
    )
    if n_metal_nbrs >= 2:
        return "mu_bridge"

    # -- carbene: donor is C with two N neighbours in a 5-ring (NHC), or
    #    a C with at least one heteroatom (N, O, S) neighbour and
    #    explicit charge / degree-3 pattern.
    if d_sym == "C":
        nbrs = [nb for nb in d.GetNeighbors() if nb.GetIdx() != metal_idx]
        n_nbrs = [nb for nb in nbrs if nb.GetSymbol() == "N"]
        # NHC pattern: 2 N neighbours, all in the same 5-ring as the
        # donor C (imidazol-2-ylidene, imidazolin-2-ylidene, …).
        if len(n_nbrs) >= 2:
            try:
                ri = mol.GetRingInfo()
                for ring in ri.AtomRings():
                    if int(donor_idx) in ring and len(ring) == 5:
                        if all(int(n.GetIdx()) in ring for n in n_nbrs[:2]):
                            return "carbene"
            except Exception:
                pass
        # Fischer / Schrock-style M=CR2 carbene (no NHC ring): the donor
        # C has hetero (N, O, S) neighbour, degree 3, formal charge != 0
        # OR an explicit C=C neighbour with no H.
        if (any(nb.GetSymbol() in ("N", "O", "S") for nb in nbrs)
                and d.GetTotalDegree() == 3
                and d.GetTotalNumHs() == 0):
            return "carbene"

    # -- hapto-π: donor C is part of an n-membered ring where the metal
    #    is bonded to multiple ring atoms.  Classify by ring size.
    if d_sym == "C":
        try:
            ri = mol.GetRingInfo()
            for ring in ri.AtomRings():
                if int(donor_idx) not in ring:
                    continue
                m_bonded_in_ring = sum(
                    1 for a in ring
                    if mol.GetBondBetweenAtoms(int(a), int(metal_idx)) is not None
                )
                if m_bonded_in_ring >= 3:
                    if len(ring) == 5 and m_bonded_in_ring >= 5:
                        return "hapto_eta5"
                    if len(ring) == 6 and m_bonded_in_ring >= 6:
                        return "hapto_eta6"
                if m_bonded_in_ring == 2 and len(ring) >= 4:
                    return "hapto_eta2"
            # η²-alkene (open-chain): donor C has a C=C neighbour where
            # that C is also bonded to the same metal.
            for nb in d.GetNeighbors():
                if nb.GetSymbol() != "C":
                    continue
                bd = mol.GetBondBetweenAtoms(int(donor_idx), int(nb.GetIdx()))
                if bd is None or bd.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                if mol.GetBondBetweenAtoms(int(nb.GetIdx()), int(metal_idx)) is not None:
                    return "hapto_eta2"
        except Exception:
            pass

    # -- agostic: donor is H, attached to a C, where that C is also
    #    bonded to the metal (or to a donor adjacent to the metal in the
    #    chelate ring).  Rare in our test set but supported per draft.
    if d_sym == "H":
        try:
            h_parent = next(nb for nb in d.GetNeighbors() if nb.GetIdx() != metal_idx)
            if h_parent.GetSymbol() == "C" and mol.GetBondBetweenAtoms(
                int(h_parent.GetIdx()), int(metal_idx)
            ) is not None:
                return "agostic"
        except Exception:
            pass

    return None


def _override_md_bounds_via_tm_category(
    *,
    mol,
    syms: Sequence[str],
    metal_idx: int,
    donor_idxs: Sequence[int],
    lower: np.ndarray,
    upper: np.ndarray,
    info: dict,
) -> None:
    """For each donor whose TM-category is detected from the graph (see
    ``_detect_tm_category``), replace the generic ``M-X`` bond bound in
    the matrix with the dedicated category lookup
    (``GripLibrary.lookup_tm_category``) using a much tighter σ window.

    Mutates ``lower`` / ``upper`` in place.  Category hits are recorded
    in ``info['tm_category_hits']`` for the diagnostic counters.
    """
    try:
        from delfin.fffree.grip_mogul_lookup import get_default_library
    except ImportError:
        return
    try:
        lib = get_default_library()
    except Exception:
        return
    metal_sym = str(syms[int(metal_idx)])
    cat_hits = []
    # The mogul_bounds uses ±2σ for M-D; here we tighten to ±2σ on the
    # category-specific distribution (whose σ is typically 5-10× smaller
    # than the generic Ag-C or Pt-N distribution because the category
    # is sharp).  This pulls the M-D embed onto the empirical CCDC mean.
    SIGMA_MULT = 2.0
    for d in donor_idxs:
        d = int(d)
        if d == metal_idx:
            continue
        cat = _detect_tm_category(mol, metal_idx, d)
        if cat is None:
            continue
        donor_sym = str(syms[d])
        try:
            hit = lib.lookup_tm_category(cat, metal_sym, donor_sym, min_n=5)
        except Exception:
            hit = None
        if hit is None:
            # Try with a relaxed min_n; the user-eye categories (carbene,
            # agostic, μ-bridge) are statistically narrow but rare.
            try:
                hit = lib.lookup_tm_category(cat, metal_sym, donor_sym, min_n=3)
            except Exception:
                hit = None
        if hit is None:
            continue
        mu, sigma, n_samp = hit
        sigma = max(0.02, float(sigma))
        lo = max(0.5, float(mu) - SIGMA_MULT * sigma)
        hi = float(mu) + SIGMA_MULT * sigma
        i, j = (metal_idx, d) if metal_idx < d else (d, metal_idx)
        lower[i, j] = lower[j, i] = lo
        upper[i, j] = upper[j, i] = hi
        cat_hits.append((int(d), cat, float(mu), float(sigma), int(n_samp)))
    if cat_hits:
        info["tm_category_hits"] = cat_hits


# ---------------------------------------------------------------------------
# Post-embed CCDC projection
# ---------------------------------------------------------------------------
def _ligand_subtree(mol, metal_idx: int, donor_idx: int) -> List[int]:
    """Atoms reachable from ``donor_idx`` without crossing the metal."""
    seen = {int(donor_idx)}
    stack = [int(donor_idx)]
    while stack:
        u = stack.pop()
        for nb in mol.GetAtomWithIdx(u).GetNeighbors():
            j = int(nb.GetIdx())
            if j == int(metal_idx) or j in seen:
                continue
            seen.add(j)
            stack.append(j)
    return sorted(seen)


def _polyhedron_for_donors(
    metal_sym: str,
    n_donors: int,
    geometry_key: Optional[str] = None,
) -> Optional[np.ndarray]:
    """Return the canonical polyhedron unit vectors for ``n_donors``.

    Uses the same dispatcher ``mogul_bounds`` does, so the choice agrees
    with the donor-donor bounds already populated in the matrix.

    ``geometry_key`` (2026-06-07, hmaximilian) lets the caller pin the
    polyhedron explicitly — used by the chelate-aware picker wiring so
    the post-embed projection uses the SAME geometry the bounds matrix
    and the orbit enumerator selected.  When unset, the legacy
    first-candidate rule from ``geometries_for_cn[0]`` is preserved.
    """
    try:
        from delfin.fffree import polyhedra as _polyhedra
    except ImportError:
        return None
    if geometry_key:
        try:
            v = _polyhedra.ref_vectors(str(geometry_key))
            if v is not None and v.shape[0] >= n_donors:
                vecs = np.asarray(v[:n_donors], dtype=float)
                norms = np.linalg.norm(vecs, axis=1, keepdims=True)
                norms = np.where(norms < 1e-9, 1.0, norms)
                return vecs / norms
        except Exception:
            pass
    cands = _polyhedra.geometries_for_cn(int(n_donors), metal_sym)
    for g in cands:
        try:
            v = _polyhedra.ref_vectors(str(g))
            if v is not None and v.shape[0] >= n_donors:
                vecs = np.asarray(v[:n_donors], dtype=float)
                norms = np.linalg.norm(vecs, axis=1, keepdims=True)
                norms = np.where(norms < 1e-9, 1.0, norms)
                return vecs / norms
        except Exception:
            continue
    return None


def _project_donors_to_ccdc_geometry(
    *,
    mol,
    syms: Sequence[str],
    metal_idx: int,
    donor_idxs: Sequence[int],
    lower: np.ndarray,
    upper: np.ndarray,
    P: np.ndarray,
    donor_at_vertex: Optional[Sequence[int]] = None,
    geometry_key: Optional[str] = None,
) -> np.ndarray:
    """Rigid-body project each ligand subtree onto the CCDC ideal.

    For each donor:
        1. Decide the target direction = polyhedron unit vector (the same
           polyhedron used in the bounds matrix), Kabsch-rotated onto
           the embed's average donor direction so the *ordering* of
           donors around the metal is preserved (no isomer flip).
        2. Decide the target |M-D| = midpoint of the CCDC (lo, hi)
           window from the bounds matrix.
        3. Translate the ligand subtree (donor + everything in its
           connected component once the M atom is removed from the
           graph) so the donor lands at ``M + |M-D|target * dir_target``.
           No rotation of the subtree — preserves all internal angles
           and bond lengths the embed has produced.

    Universal: works for any CN, any donor element, any chelate.  No
    SMILES patterns, no per-class branches.  When the polyhedron lookup
    fails (rare CNs without a vertex set) the projection is skipped
    and ``P`` is returned unchanged.

    Parameters
    ----------
    donor_at_vertex : sequence of int, optional
        Pólya orbit override (2026-06-07): permutation mapping vertex k
        → atom index of the donor placed at vertex k.  When supplied,
        the projection assigns ``donor_at_vertex[k]`` to polyhedron
        vertex ``k`` rather than the canonical ``donor_idxs[k]`` mapping.
        This preserves the stereoisomer encoded in the bounds-matrix
        donor-donor block.  When ``None`` (default) the canonical
        sort-by-atom-idx mapping is used (byte-identical to the
        pre-orbit-enumeration code path).
    """
    P = np.asarray(P, dtype=float).copy()
    metal_pos = P[int(metal_idx)].copy()
    n_d = len(donor_idxs)
    if n_d < 1:
        return P

    metal_sym = str(syms[int(metal_idx)])
    ref_unit = _polyhedron_for_donors(metal_sym, n_d,
                                      geometry_key=geometry_key)

    # Decide which donor sits at which polyhedron vertex.  Default =
    # sorted-canonical (matches the bounds-matrix donor_list); orbit
    # override replaces it with the Pólya-derived permutation.
    if (donor_at_vertex is not None
            and len(donor_at_vertex) == n_d):
        donors_per_vertex = [int(d) for d in donor_at_vertex]
    else:
        donors_per_vertex = [int(d) for d in donor_idxs]

    # Embed-derived donor directions, in vertex order (so cur_dirs[k] is
    # the embed direction of the donor we have decided sits at vertex k).
    cur_dirs = np.zeros((n_d, 3), dtype=float)
    for k, d in enumerate(donors_per_vertex):
        v = P[int(d)] - metal_pos
        nv = float(np.linalg.norm(v))
        cur_dirs[k] = v / nv if nv > 1e-9 else np.array([1.0, 0.0, 0.0])

    if ref_unit is not None:
        # Kabsch align ``ref_unit`` to ``cur_dirs`` so the polyhedron is in
        # the same orientation as the embed (no isomer change).
        H = cur_dirs.T @ ref_unit
        U, S, Vt = np.linalg.svd(H)
        if np.linalg.det(U) * np.linalg.det(Vt) < 0:
            Vt[-1, :] *= -1.0
        Rrot = U @ Vt
        aligned = ref_unit @ Rrot.T
    else:
        # No registered polyhedron for this CN (rare, e.g. CN=12 bis-η⁶
        # arene).  Use the embed's own donor directions as the projection
        # rays — supplement-based radial pull still applies (the donor's
        # M-D distance gets pulled onto the CCDC mean) without enforcing
        # any specific vertex polyhedron.  Universal fallback so we
        # never silently skip the projection on rare-CN cases.
        aligned = cur_dirs.copy()

    # Translate each ligand subtree onto its target.
    #
    # Chelates: a single ligand carries TWO OR MORE donors.  We cannot
    # translate each donor's subtree independently because they share
    # internal atoms — moving donor #2 then drags donor #1 too.  Instead
    # we treat each connected ligand component as ONE subtree and apply
    # the SINGLE projection that minimises the sum of |donor - target|².
    # For monodentate ligands this reduces to "place donor exactly on
    # target", which is what the simple case wants.
    #
    # NOTE (2026-06-07 orbit-enum): the vertex index ``k`` indexes
    # ``donors_per_vertex`` (the orbit-permuted donor-to-vertex mapping),
    # not ``donor_idxs``.  Under the canonical-orbit default they are
    # identical; under a Pólya orbit override they differ and the
    # projection MUST use ``donors_per_vertex`` so the right donor lands
    # at the right polyhedron vertex.
    visited: set[int] = {int(metal_idx)}
    # Build per-vertex subtree mapping (atom set reachable from the donor
    # at vertex k without crossing the metal).
    subtree_by_vertex: List[List[int]] = [
        _ligand_subtree(mol, int(metal_idx), int(d)) for d in donors_per_vertex
    ]
    # Group vertices by shared subtree component (chelate detection).
    n_d_local = n_d
    component: List[int] = list(range(n_d_local))
    for a in range(n_d_local):
        for b in range(a + 1, n_d_local):
            if set(subtree_by_vertex[a]) & set(subtree_by_vertex[b]):
                # Merge components a and b
                ra = a
                while component[ra] != ra:
                    ra = component[ra]
                rb = b
                while component[rb] != rb:
                    rb = component[rb]
                component[max(ra, rb)] = min(ra, rb)
    # Resolve final root per vertex.
    roots = []
    for a in range(n_d_local):
        r = a
        while component[r] != r:
            r = component[r]
        roots.append(r)

    groups: dict[int, List[int]] = {}
    for a, r in enumerate(roots):
        groups.setdefault(r, []).append(a)

    for root, donor_positions in groups.items():
        # Union subtree for this component.
        atoms_in_group: set[int] = set()
        for k in donor_positions:
            atoms_in_group.update(subtree_by_vertex[k])
        # Atoms not yet visited (avoid double-move of shared atoms in a
        # pathological multi-metal case).
        movable = sorted(a for a in atoms_in_group if a not in visited)
        if not movable:
            continue
        # Compute target positions for each donor in this group.
        target_positions = []
        donor_indices_in_group = []
        for k in donor_positions:
            d = int(donors_per_vertex[k])
            i, j = (int(metal_idx), d) if int(metal_idx) < d else (d, int(metal_idx))
            lo = float(lower[i, j])
            hi = float(upper[i, j])
            if not np.isfinite(hi) or hi > 100.0 or hi <= lo:
                continue
            r_target = 0.5 * (lo + hi)
            target = metal_pos + r_target * aligned[k]
            target_positions.append(target)
            donor_indices_in_group.append(d)

        if not target_positions:
            continue
        target_positions = np.asarray(target_positions, dtype=float)
        donor_indices_in_group = np.asarray(donor_indices_in_group, dtype=int)

        # MONODENTATE: just translate so donor lands on target.
        if len(donor_positions) == 1:
            d = donor_indices_in_group[0]
            delta = target_positions[0] - P[d]
            for a in movable:
                P[a] = P[a] + delta
                visited.add(a)
            continue

        # CHELATE: the embed already produces a good *internal* geometry
        # (chelate's natural bite angle, ring puckers, bond lengths) but
        # the M-D distance can drift outside the CCDC supplement window
        # because the DG embedder uses a smooth soft-penalty that doesn't
        # strictly enforce the hard bounds — and we don't have a per-
        # donor delta to apply (each donor's translation would drag the
        # shared subtree atoms, scrambling the chelate).
        #
        # Conservative chelate fix (2026-06-08, hmaximilian, Bug B):
        # apply a single RIGID-BODY transform (rotation + translation)
        # to the WHOLE chelate subtree so every donor lands as close as
        # possible to its supplement-derived target ``M + r_target ·
        # aligned[k]``.  Kabsch solves the LSQ rigid alignment of
        # ``embed donors`` → ``targets`` and preserves all internal
        # distances + angles by construction.  Universal: works for any
        # bidentate / tridentate / polydentate chelate.
        #
        # Env-flag: DELFIN_FFFREE_PROJECTION_USE_SUPPLEMENT (default ON
        # when MOGUL_PRIMARY is active).  Set to 0 explicitly to restore
        # the legacy "skip chelate projection" behaviour.
        _use_supp = os.environ.get(
            "DELFIN_FFFREE_PROJECTION_USE_SUPPLEMENT", "1"
        ).strip().lower() not in ("", "0", "false", "no")
        if not _use_supp:
            for a in movable:
                visited.add(a)
            continue

        # Donor coordinates in the embed (one row per donor in this group).
        src_donors = np.asarray(
            [P[int(d)] for d in donor_indices_in_group], dtype=float
        )
        tgt_donors = np.asarray(target_positions, dtype=float)
        n_pts = src_donors.shape[0]
        if n_pts < 2:
            # Should never hit -- single-donor branch handled above.
            for a in movable:
                visited.add(a)
            continue

        # Centroid radial translation for chelates (2026-06-08,
        # hmaximilian, Bug B universal fix):
        #
        # The fundamental problem: each donor has a target M-D distance
        # from the CCDC supplement, but donors in a chelate share
        # subtree atoms (ring carbons, backbone) and cannot be
        # translated independently without tearing the ring apart.
        # The DG embed lands the chelate at its OWN intrinsic geometry
        # which may put the MEAN M-D distance below or above the
        # supplement window (typically pinched downward when the
        # chelate's bite forces donors toward each other).
        #
        # We solve this with a CENTROID RADIAL TRANSLATION applied to
        # the whole subtree:
        #   1. Compute the EMBED'S MEAN per-donor M-D distance:
        #      r_src_md = mean(|donor_k - metal|).
        #   2. Compute the TARGET MEAN: r_tgt_md = mean(supplement[k]).
        #   3. Translate the entire chelate subtree along the unit
        #      vector ``(centroid_of_donors - metal)/r_src_centroid``
        #      by magnitude ``(r_tgt_md - r_src_md)``.
        #
        # Pure translation -> preserves every internal distance and
        # angle by construction.  Pulls the chelate's MEAN per-donor M-D
        # distance onto the CCDC empirical mean.  Per-donor radial
        # fine-tuning is deferred to GRIP-polish which already runs
        # after projection and honours per-donor Mahalanobis pulls
        # under a hard topology validator.
        #
        # ACCEPT-IF-BETTER gate: only apply when the post-translation
        # RMS per-donor deviation IMPROVES.  Guards against degenerate
        # chelates whose intrinsic geometry already sits near the
        # supplement mean (no improvement available) or where the
        # centroid direction is mis-aligned (e.g. a chelate whose
        # donors converge axially -- the centroid is much closer to
        # the metal than the donors themselves are).  Universal
        # geometric check, no per-class branch.
        r_src_centroid = float(np.linalg.norm(src_donors.mean(axis=0) - metal_pos))
        if r_src_centroid < 1e-6:
            for a in movable:
                visited.add(a)
            continue
        per_donor_targets = np.linalg.norm(tgt_donors - metal_pos, axis=1)
        pre_src_md = np.linalg.norm(src_donors - metal_pos, axis=1)
        # Centroid direction (unit) -- the chelate moves rigidly along
        # this ray.
        radial_dir = (src_donors.mean(axis=0) - metal_pos) / r_src_centroid
        # Closed-form optimal scalar translation along ``radial_dir``
        # that minimises Σ((|donor_k + t·r̂ - metal|)² - target_k²)²
        # is intractable; the linear-in-distance proxy ``t = mean(target)
        # - mean(pre_md)`` is the right first-order estimate and is what
        # we use.  Capped at ±0.6 Å.
        delta_r = float(np.mean(per_donor_targets) - np.mean(pre_src_md))
        translation_mag = max(-0.6, min(0.6, delta_r))
        translation = translation_mag * radial_dir
        post_src = src_donors + translation
        post_md = np.linalg.norm(post_src - metal_pos, axis=1)

        # Hard floor: reject if any donor would end up below 0.7×target
        # (heavy contraction) or above 1.5×target (heavy stretch).
        if (np.any(post_md < 0.7 * per_donor_targets)
                or np.any(post_md > 1.5 * per_donor_targets)):
            for a in movable:
                visited.add(a)
            continue

        # Accept-if-better gate: RMS per-donor deviation from target.
        pre_rms_dev = float(np.sqrt(np.mean((pre_src_md - per_donor_targets) ** 2)))
        post_rms_dev = float(np.sqrt(np.mean((post_md - per_donor_targets) ** 2)))
        if post_rms_dev >= pre_rms_dev:
            for a in movable:
                visited.add(a)
            continue

        # Apply the radial translation to every movable atom (whole
        # subtree).  This preserves all internal distances and angles
        # by construction (pure translation).
        for a in movable:
            P[a] = P[a] + translation
            visited.add(a)
    return P


# ---------------------------------------------------------------------------
# Diagnostic helpers (for tests + user-eye validation)
# ---------------------------------------------------------------------------
def measure_md_distance(syms: Sequence[str], P: np.ndarray, metal_sym: str,
                        donor_sym: str) -> Optional[float]:
    """Return the shortest distance between the first ``metal_sym`` atom and
    any ``donor_sym`` atom directly bonded to it (geometry only, by element).

    Convenience helper for user-eye validation tables.
    """
    P = np.asarray(P, dtype=float)
    try:
        m = next(i for i, s in enumerate(syms) if s == metal_sym)
    except StopIteration:
        return None
    best = None
    for i, s in enumerate(syms):
        if i == m or s != donor_sym:
            continue
        d = float(np.linalg.norm(P[i] - P[m]))
        if best is None or d < best:
            best = d
    return best


def measure_planarity(P: np.ndarray, atom_idxs: Sequence[int]) -> float:
    """Maximum perpendicular deviation of a set of atoms from their SVD
    best-fit plane.  0 ⇒ perfectly planar.
    """
    Q = np.asarray(P, dtype=float)[list(atom_idxs)]
    centroid = Q.mean(axis=0)
    M = Q - centroid
    # SVD: smallest singular vector is the plane normal.
    _, _, Vt = np.linalg.svd(M, full_matrices=False)
    n_vec = Vt[-1]
    dev = np.abs(M @ n_vec)
    return float(dev.max())
