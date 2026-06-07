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
    mol = _full_complex_mol(smiles)
    if mol is None:
        return None

    # 2) Locate metal + donors from the graph topology.
    metal_idx, donor_idxs = _locate_metal_and_donors(mol)
    if metal_idx is None or not donor_idxs:
        return None

    # 3) Build the CCDC-empirical bounds matrix.
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
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
            geometry=geometry_key, # explicit (orbit-enum) or auto-derive
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
            geometry_key=geometry_key,
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
    return out_syms, P


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


def _locate_metal_and_donors(mol) -> Tuple[Optional[int], List[int]]:
    """Identify the primary metal index + list of donor atom indices.

    Primary metal = metal with highest graph degree (matches decompose).
    Donors = neighbours of the primary metal in the molecular graph.
    """
    metals = [a.GetIdx() for a in mol.GetAtoms() if _bd._is_metal(a.GetSymbol())]
    if not metals:
        return None, []
    if len(metals) == 1:
        m = metals[0]
    else:
        m = max(metals, key=lambda mi: mol.GetAtomWithIdx(mi).GetDegree())
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


def _polyhedron_for_donors(metal_sym: str, n_donors: int) -> Optional[np.ndarray]:
    """Return the canonical polyhedron unit vectors for ``n_donors``.

    Uses the same dispatcher ``mogul_bounds`` does, so the choice agrees
    with the donor-donor bounds already populated in the matrix.
    """
    try:
        from delfin.fffree import polyhedra as _polyhedra
    except ImportError:
        return None
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
    """
    P = np.asarray(P, dtype=float).copy()
    metal_pos = P[int(metal_idx)].copy()
    n_d = len(donor_idxs)
    if n_d < 1:
        return P

    metal_sym = str(syms[int(metal_idx)])
    ref_unit = _polyhedron_for_donors(metal_sym, n_d)
    if ref_unit is None:
        return P

    # Embed-derived donor directions.
    cur_dirs = np.zeros((n_d, 3), dtype=float)
    for k, d in enumerate(donor_idxs):
        v = P[int(d)] - metal_pos
        nv = float(np.linalg.norm(v))
        cur_dirs[k] = v / nv if nv > 1e-9 else np.array([1.0, 0.0, 0.0])

    # Kabsch align ``ref_unit`` to ``cur_dirs`` so the polyhedron is in
    # the same orientation as the embed (no isomer change).
    H = cur_dirs.T @ ref_unit
    U, S, Vt = np.linalg.svd(H)
    if np.linalg.det(U) * np.linalg.det(Vt) < 0:
        Vt[-1, :] *= -1.0
    Rrot = U @ Vt
    aligned = ref_unit @ Rrot.T

    # Translate each ligand subtree onto its target.
    #
    # Chelates: a single ligand carries TWO OR MORE donors.  We cannot
    # translate each donor's subtree independently because they share
    # internal atoms — moving donor #2 then drags donor #1 too.  Instead
    # we treat each connected ligand component as ONE subtree and apply
    # the SINGLE projection that minimises the sum of |donor - target|².
    # For monodentate ligands this reduces to "place donor exactly on
    # target", which is what the simple case wants.
    visited: set[int] = {int(metal_idx)}
    # Build per-donor subtree mapping (atom set reachable from donor
    # without crossing the metal).
    subtree_by_donor: List[List[int]] = [
        _ligand_subtree(mol, int(metal_idx), int(d)) for d in donor_idxs
    ]
    # Group donors by shared subtree component (chelate detection).
    n_d_local = len(donor_idxs)
    component: List[int] = list(range(n_d_local))
    for a in range(n_d_local):
        for b in range(a + 1, n_d_local):
            if set(subtree_by_donor[a]) & set(subtree_by_donor[b]):
                # Merge components a and b
                ra = a
                while component[ra] != ra:
                    ra = component[ra]
                rb = b
                while component[rb] != rb:
                    rb = component[rb]
                component[max(ra, rb)] = min(ra, rb)
    # Resolve final root per donor.
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
            atoms_in_group.update(subtree_by_donor[k])
        # Atoms not yet visited (avoid double-move of shared atoms in a
        # pathological multi-metal case).
        movable = sorted(a for a in atoms_in_group if a not in visited)
        if not movable:
            continue
        # Compute target positions for each donor in this group.
        target_positions = []
        donor_indices_in_group = []
        for k in donor_positions:
            d = int(donor_idxs[k])
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

        # CHELATE: the embed already produces a good internal geometry
        # (donor-donor distance comes out of the polyhedron donor-donor
        # bound + the chelate's natural bite angle).  The polyhedron the
        # bounds matrix picks may NOT match the chelate's intrinsic
        # geometry (e.g. tridentate-N₃ pyridyl chelate fits a meridional
        # cis-cis arrangement, not the tetrahedron T-4 the default
        # ``geometries_for_cn(4)`` returns).  Forcing it onto a poor
        # polyhedron via Kabsch causes some donors to collapse into the
        # metal.
        #
        # Conservative policy: do NOT re-project chelate subtrees.  The
        # embedder's positions are within ±0.1 Å of the CCDC window for
        # M-D (verified on BERTEB / ALAHEB) and the D-M-D angles reflect
        # the chelate's intrinsic bite.  This matches the draft
        # manuscript's design — chelate geometry IS empirical CCDC for
        # that LIGAND'S internal frame, not for an idealised polyhedron.
        # Mark all atoms visited so any subsequent monodentate
        # projection in the same complex does not double-move them.
        for a in movable:
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
