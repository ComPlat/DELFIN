"""Geometric inter-ligand clash relief, template bond orders and the Open Babel UFF optimisation of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import math
import os
from typing import Dict, Optional

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _xyz_to_rdkit_conformer,
)
from delfin.manta.converter_flags import (
    _apply_uff_jitter,
    _class_conditional_flag,
    _delfin_env_int,
)
from delfin.manta.geometry_quality import (
    _geometry_quality_score,
    _has_bad_geometry,
)
from delfin.manta.ligand_placement import (
    _snap_aromatic_rings_in_xyz,
)
from delfin.manta.ml_tables import (
    Chem,
    OPENBABEL_AVAILABLE,
    RDKIT_AVAILABLE,
    _COVALENT_RADII,
    _METAL_SET,
    _get_ml_bond_length,
    pybel,
)
from delfin.manta.pre_uff_snap import (
    _clamp_metalloid_md_xyz,
    _flatten_d8_sq_planar_xyz,
)
from delfin.manta.topology_checks import (
    _flatten_sp2_atoms_xyz,
    _fragment_topology_ok,
    _has_severe_covalent_distortion,
    _no_spurious_bonds,
    _roundtrip_ring_count_ok,
    _verify_metal_connectivity,
)
from delfin.manta.uff_constraints import (
    _VDW_RADII_CLASH,
    _build_coordination_constraints_from_xyz,
    _build_uff_constraints_from_template,
)

logger = get_logger("delfin.smiles_converter")


def _geometric_inter_clash_relief(
    xyz_delfin: str,
    frozen_idx,
    threshold: float = 0.85,
    max_passes: int = 40,
    step_frac: float = 0.5,
):
    """Deterministic, FF-free inter-ligand clash relief.

    Replaces the (garbage-gradient, non-deterministic) OB-UFF
    ConjugateGradients step for UFF-unparameterised-metal complexes.  Mirrors
    the inter-ligand-clash detector's geometry contract exactly so the moves it
    makes target the metric it is meant to improve:

      1. infer bonds from covalent radii (metal_extra_tol 0.45, organic 0.40);
      2. cluster atoms into ligand fragments AFTER removing every
         metal-touching bond (so two ligands = two fragments);
      3. for every pair of atoms in DISTINCT fragments closer than
         ``threshold * (vdW_i + vdW_j)``, push the two atoms apart symmetrically
         along their separation axis.

    Determinism is guaranteed by construction: no force field, no RNG, fixed
    iteration order by sorted atom index, and a deterministic separation axis
    even for near-degenerate (overlapping) atoms (the project's
    sorted-index-seeded canonical-axis idiom).  Frozen atoms (the metal + its
    pinned donors + hapto ring atoms) never move; only one atom of a clashing
    pair is displaced when its partner is frozen (full correction onto the free
    atom), otherwise both move half-and-half.  Bonded atoms and same-fragment
    atoms are never pushed apart, so covalent connectivity and chelate-ring
    geometry are preserved.

    Args:
        xyz_delfin: DELFIN-format coordinate block (``symbol x y z`` per line).
        frozen_idx: iterable of 0-based atom indices that must NOT move
            (metal + donors + hapto atoms = ``constraints['fix_atoms']``).
        threshold: vdW-overlap fraction below which a pair is a clash (0.85,
            matching the detector default).
        max_passes: bounded relaxation passes.
        step_frac: fraction of the missing separation removed per pass
            (0.5 → half the deficit per pass; converges without overshoot).

    Returns:
        Relaxed DELFIN-format XYZ string (or the input unchanged on any error).
    """
    try:
        import numpy as _np

        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        n = len(lines)
        if n == 0:
            return xyz_delfin
        syms: list = []
        pos = _np.zeros((n, 3), dtype=float)
        for i, ln in enumerate(lines):
            p = ln.split()
            syms.append(p[0])
            pos[i, 0] = float(p[1])
            pos[i, 1] = float(p[2])
            pos[i, 2] = float(p[3])

        frozen = set(int(x) for x in (frozen_idx or []))
        is_metal = [s in _METAL_SET for s in syms]
        cov = _np.array([_COVALENT_RADII.get(s, 1.5) for s in syms], dtype=float)
        vdw = _np.array(
            [_VDW_RADII_CLASH.get(s, 1.70) for s in syms], dtype=float
        )

        # --- (1) bond inference (mirror derive_bonds) -------------------
        # bonded[i][j] True if i,j within (cov_i+cov_j+tol); metal-touching
        # bonds use the wider 0.45 tol, organic pairs 0.40.
        bonded = [set() for _ in range(n)]
        d_all = _np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
        for i in range(n):
            for j in range(i + 1, n):
                tol = 0.45 if (is_metal[i] or is_metal[j]) else 0.40
                if d_all[i, j] <= (cov[i] + cov[j] + tol):
                    bonded[i].add(j)
                    bonded[j].add(i)

        # --- (2) ligand fragments after removing metal-touching bonds ---
        metal_set = {i for i in range(n) if is_metal[i]}
        adj = [[] for _ in range(n)]
        for i in range(n):
            if i in metal_set:
                continue
            for j in bonded[i]:
                if j in metal_set or j <= i:
                    continue
                adj[i].append(j)
                adj[j].append(i)
        frag = [-1] * n
        cid = 0
        for start in range(n):
            if start in metal_set or frag[start] != -1:
                continue
            frag[start] = cid
            stack = [start]
            while stack:
                u = stack.pop()
                for v in adj[u]:
                    if frag[v] == -1:
                        frag[v] = cid
                        stack.append(v)
            cid += 1

        # --- (3) bounded symmetric push-apart of inter-fragment clashes -
        moved_any = False
        for _ in range(max_passes):
            # Collect all current inter-fragment clashes (deterministic order
            # by (i, j) ascending), then apply displacements in that order.
            disp = _np.zeros((n, 3), dtype=float)
            n_clash = 0
            d_cur = _np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
            for i in range(n):
                if i in metal_set or frag[i] < 0:
                    continue
                for j in range(i + 1, n):
                    if j in metal_set or frag[j] < 0:
                        continue
                    if frag[i] == frag[j]:
                        continue  # same ligand → intra, not our concern
                    if j in bonded[i]:
                        continue  # directly bonded across (rare) → leave
                    floor = threshold * (vdw[i] + vdw[j])
                    d = float(d_cur[i, j])
                    if d >= floor:
                        continue
                    n_clash += 1
                    deficit = floor - d
                    # Separation axis; deterministic fallback for the
                    # (near-)degenerate overlap case: a unit axis seeded by the
                    # sorted atom indices (no RNG), the project idiom.
                    axis = pos[i] - pos[j]
                    na = float(_np.linalg.norm(axis))
                    if na < 1e-9:
                        # canonical deterministic axis from index parity
                        k = (i * 131 + j) % 3
                        axis = _np.zeros(3)
                        axis[k] = 1.0
                        na = 1.0
                    axis = axis / na
                    push = step_frac * deficit * axis
                    fi = i in frozen
                    fj = j in frozen
                    if fi and fj:
                        continue  # both pinned: cannot relieve, skip
                    elif fi:
                        disp[j] -= push  # all correction onto free atom j
                    elif fj:
                        disp[i] += push  # all correction onto free atom i
                    else:
                        disp[i] += 0.5 * push
                        disp[j] -= 0.5 * push
            if n_clash == 0:
                break
            # Apply (frozen atoms forced to zero displacement for safety).
            for f in frozen:
                if 0 <= f < n:
                    disp[f] = 0.0
            pos = pos + disp
            moved_any = True

        if not moved_any:
            return xyz_delfin
        out = []
        for i in range(n):
            out.append(
                f"{syms[i]:4s} {pos[i,0]:12.6f} {pos[i,1]:12.6f} {pos[i,2]:12.6f}"
            )
        return "\n".join(out) + "\n"
    except Exception as _exc:
        logger.debug("Geometric inter-clash relief failed (%s); input kept", _exc)
        return xyz_delfin


def _apply_template_bond_orders(ob_mol, mol_template) -> int:
    """Write bond orders from the RDKit template into the OB molecule.

    WHY.  ``pybel.readstring("xyz", ...)`` hands over a BARE coordinate list.
    OpenBabel has to guess bonds AND bond orders from it itself
    (ConnectTheDots + PerceiveBondOrders), and UFF chooses its atom types based on
    this estimate.  But the bond order is NOT a geometric quantity -- it
    stands in the SMILES, i.e. in the template, and is thrown away here and guessed afterwards.
    ⛔ CAUTION, 14.08.2026: THIS RATIONALE IS NOT CONFIRMED.  Work item 2.9 attributed
    the unphysical C-S distance (measured 1.37 / 1.41 / 1.46 A against a shortest
    real C=S of 1.55) to this missing bond order.  Re-checked the same day:
    OB perceives C=S on an isolated thioketone at 1.61, 1.50, 1.433 AND 1.35 A
    every time CORRECTLY as order 2.  So the order is not lost, at least
    not without a metal nearby.  The DEFECT is real, the CAUSE is open.

    This switch therefore stays OFF and is NOT a repair of 2.9, but a
    hypothesis that can still fail.  Before a measurement, first reproduce the error on one of
    the three named frames -- in the metal complex, not on the model molecule.
    What is nevertheless right here: a quantity that STANDS in the template should not be
    left to guessing.

    ⚠ STRICT MAPPING.  Transfer happens ONLY if atom count AND symbol sequence match
    exactly.  The emitted XYZ need not follow the template order
    (the same concern is stated verbatim in _bond_decollapse), and writing a bond order onto
    the wrong atom pair would be worse than guessing it.  On any
    deviation: leave unchanged and let OB perceive as before.

    Return: number of bonds set (0 = nothing touched).
    """
    if mol_template is None or not RDKIT_AVAILABLE:
        return 0
    try:
        n_t = mol_template.GetNumAtoms()
        if ob_mol.NumAtoms() != n_t:
            return 0
        for i in range(n_t):
            ob_a = ob_mol.GetAtom(i + 1)            # OB zaehlt ab 1
            if pybel.ob.GetSymbol(ob_a.GetAtomicNum()) != \
                    mol_template.GetAtomWithIdx(i).GetSymbol():
                return 0                            # Reihenfolge weicht ab -> Finger weg
        _ORDER = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2,
                  Chem.BondType.TRIPLE: 3, Chem.BondType.AROMATIC: 5}
        n_set = 0
        for b in mol_template.GetBonds():
            o = _ORDER.get(b.GetBondType())
            if not o:
                continue
            ob_b = ob_mol.GetBond(b.GetBeginAtomIdx() + 1, b.GetEndAtomIdx() + 1)
            if ob_b is None:
                continue                            # OB did not perceive the bond
            if o == 5:
                ob_b.SetAromatic(True); ob_b.SetBondOrder(1)
            else:
                ob_b.SetBondOrder(o)
            n_set += 1
        return n_set
    except Exception:
        return 0


def _optimize_xyz_openbabel(
    xyz_delfin: str,
    steps: int = 500,
    constraints: Optional[Dict] = None,
    return_energy: bool = False,
    mol_template=None,
):
    """Optimize a DELFIN-format XYZ string using Open Babel's UFF force field.

    Open Babel's UFF implementation includes full parameters for transition
    metals (Fe, Co, Ni, Ti, etc.) that RDKit lacks.  This makes it ideal as
    a post-ETKDG refinement step for metal complexes.

    Args:
        xyz_delfin: Coordinate block in DELFIN format (``symbol x y z`` per line,
            no atom-count header).
        steps: Maximum number of conjugate-gradient steps.
        constraints: Optional dict with keys:
            - ``'fix_atoms'``: list of 0-based atom indices to freeze in place
            - ``'distances'``: list of ``(idx_a, idx_b, target_dist)``
            - ``'angles'``: list of ``(idx_a, idx_b, idx_c, target_angle_deg)``
            - ``'torsions'``: list of ``(idx_a, idx_b, idx_c, idx_d, target_deg)``

    Returns:
        Optimized DELFIN-format XYZ string.  If ``return_energy`` is True,
        returns ``(xyz_str, energy_kcal_mol)`` tuple; energy is ``None``
        when UFF setup fails.
    """
    if not OPENBABEL_AVAILABLE:
        return (xyz_delfin, None) if return_energy else xyz_delfin

    # Iter-8.7 (2026-05-12): hash-seeded micro-jitter on frozen atoms.
    # The topology builder places M + monodentate donors at the exact
    # symmetric _TOPO_GEOMETRY_VECTORS positions (Td/Oh/TPR/SAP/ICOS/
    # CUBO) and OB-UFF then freezes them via constraints["fix_atoms"].
    # A frozen polyhedron is a stationary point of UFF, so the template-
    # perfect bond lengths (FUPVOB Zn-Br all 2.3500 Å exact, CESNIZ Ru-N
    # 2.060/2.400 Å exact) survive UFF byte-identical and the surrounding
    # ligand frame can collapse symmetrically (ZOQTEE 49/73 atoms at
    # x ≈ 0).  A small ±0.05 Å offset on the frozen subset breaks the
    # symmetric stationary point while staying within the M-D tolerance
    # window of downstream _verify_metal_connectivity.  XYZ-hash-seeded
    # → deterministic same-input → same-jitter.
    # Phase 4D default-flip 2026-05-12: 0 → 1 — smoke500-verified +1.18pp sigma.
    # Wave-4 full-pool 11363 (A+B+I+K): jitter causes M-D break +7.71pp,
    # h-clash files% +11.44pp, poly_mean_dev_Oh +4.17°, F3-bond +7.8%, F19-Td
    # +10.3% across all classes.  Net loss far exceeds sigma gain.
    # Reverted default 1 → 0 (Wave-4 Bundle 1).  Env override still works.
    if (
        constraints
        and constraints.get("fix_atoms")
        and _delfin_env_int("DELFIN_UFF_JITTER", 0)
    ):
        try:
            _jitter_mag_ppm = _delfin_env_int("DELFIN_UFF_JITTER_PPM", 500)
            _jitter_mag = max(0.0, _jitter_mag_ppm) / 10000.0
            xyz_delfin = _apply_uff_jitter(
                xyz_delfin,
                constraints["fix_atoms"],
                xyz_delfin,  # XYZ used as deterministic seed source
                magnitude=_jitter_mag,
            )
        except Exception as exc:
            logger.debug("UFF jitter failed, continuing unjittered: %s", exc)

    try:
        # Parse DELFIN lines (skip blank lines)
        lines = [l for l in xyz_delfin.strip().splitlines() if l.strip()]
        n_atoms = len(lines)
        if n_atoms == 0:
            return (xyz_delfin, None) if return_energy else xyz_delfin

        # Build standard XYZ string (atom count + comment + coordinates)
        std_xyz = f"{n_atoms}\n\n"
        for line in lines:
            parts = line.split()
            # DELFIN format: "symbol x y z"
            std_xyz += f"{parts[0]}  {parts[1]}  {parts[2]}  {parts[3]}\n"

        # Read into Open Babel
        ob_mol = pybel.readstring("xyz", std_xyz)

        # BOND ORDERS FROM THE TEMPLATE INSTEAD OF FROM THE GEOMETRY (2026-08-14).
        # Default OFF -> byte-identical; see _apply_template_bond_orders for the why
        # (C-S at 1.433 A, 249 systems, legacy).
        if os.environ.get("DELFIN_FFFREE_OB_BOND_ORDERS", "0") == "1":
            _n_bo = _apply_template_bond_orders(ob_mol.OBMol, mol_template)
            if _n_bo:
                logger.debug("OB bond orders taken from template: %d bonds", _n_bo)

        # Run UFF conjugate-gradient optimization.
        #
        # DETERMINISM CONTRACT: ``pybel._forcefields["uff"]`` is a single
        # process-global force-field object.  Open Babel's force field
        # retains the constraint set from the previous ``SetConstraints``
        # call — there is no public "clear constraints" entry point — so
        # reusing the shared object means a constrained metal-complex
        # optimisation leaves its M-D distance pins / fixed atoms behind,
        # and the *next* call (even an unrelated non-metal molecule)
        # silently inherits them.  That makes a SMILES' output depend on
        # whatever was optimised before it in the same process.
        #
        # Fix: obtain a FRESH force-field instance per call via
        # ``OBForceField.FindForceField`` (which constructs a new object
        # rather than handing back the cached singleton).  Each
        # optimisation starts from a clean, constraint-free force-field
        # state, so the result depends only on this call's own input.
        # Falls back to the shared singleton only if the fresh-instance
        # lookup is somehow unavailable.
        ff = None
        try:
            ff = pybel.ob.OBForceField.FindForceField("uff")
        except Exception:
            ff = None
        if ff is None:
            ff = pybel._forcefields["uff"]
        # Welle-5i Agent A (2026-05-17): extreme-formal-charge OB-UFF fallback.
        # OB-UFF has no parameters for metals whose OB-perceived |formal charge|
        # exceeds normal coordination chemistry (e.g. user SMILES ``[Mo-4]``,
        # ``[Hg-6]``, ``[Re-3..-5]``, ``[Rh-3]``).  Such atoms hit the
        # unparam-TM uninitialised-memory bug and zero frames are emitted.
        # When enabled, the metal's perceived formal charge is scrubbed to 0
        # BEFORE Setup() so OB-UFF picks a parametrised atom type.  Read-back
        # at the end of the routine uses ``GetAtomicNum`` only, so this does
        # not propagate to downstream chemistry.  Universal: gated on
        # ``abs(q)`` magnitude, never on metal element or SMILES pattern.
        if _delfin_env_int("DELFIN_EXTREME_CHARGE_FALLBACK", 0):
            _ec_thresh = max(1, _delfin_env_int("DELFIN_EXTREME_CHARGE_THRESHOLD", 3))
            try:
                for _ob_atom in pybel.ob.OBMolAtomIter(ob_mol.OBMol):
                    _sym = pybel.ob.GetSymbol(_ob_atom.GetAtomicNum())
                    if _sym in _METAL_SET and abs(int(_ob_atom.GetFormalCharge())) >= _ec_thresh:
                        _ob_atom.SetFormalCharge(0)
            except Exception as _ec_exc:
                logger.debug("Extreme-charge fallback charge-scrub failed: %s", _ec_exc)
        if not ff.Setup(ob_mol.OBMol):
            logger.debug("Open Babel UFF setup failed, returning unoptimized geometry")
            return (xyz_delfin, None) if return_energy else xyz_delfin

        # DETERMINISM GATE: probe whether OB UFF actually parameterised
        # every atom type in this molecule.  When OB's UFFTYPER cannot
        # match a perceived atom type (typical for transition metals with
        # explicit formal charges — ``Pt+2``, ``Fe+2``, ``Pd+2``, ``Ru+3``
        # etc.), OB silently falls back to a default/uninitialised
        # parameter slot.  The energy stays at a marker value
        # (~5.9e11 kcal/mol — a recognisable uninitialised-memory
        # pattern) BUT the per-atom gradient on the affected atoms is
        # read from uninitialised heap memory and differs across process
        # invocations.  Calls that downstream FREEZE every donor +
        # metal via ``AddAtomConstraint`` are unaffected because the bad
        # gradients can't move fixed atoms.  Calls that swap a donor's
        # ``AddAtomConstraint`` for an ``AddDistanceConstraint`` (the
        # soft-donor mode added by Wave-1-7) leave the donor free, the
        # uninitialised gradient drags it across the workspace, and the
        # output geometry depends on whatever the freshly-allocated
        # heap page happened to contain.  The same input then produces
        # a different output every run — the bug that blocked v2-final-
        # prime even after Determinism Fix #1 + #2.
        #
        # Detection is purely runtime: an unconstrained ``ff.Energy()``
        # immediately after Setup returns the parameterisation marker.
        # Anything above 1e9 kcal/mol is the uninitialised pattern (a
        # physically-meaningful 5-atom complex has E < 1e5).  No metal
        # allowlist, no SMILES patterns.  When this flag is set, the
        # soft-donor block below skips the distance-pin path and falls
        # back to ``AddAtomConstraint`` on every donor — the legacy
        # behaviour that is empirically deterministic even with
        # unparameterised metals.
        _uff_param_unsafe = False
        try:
            _e_marker = float(ff.Energy())
            if not math.isfinite(_e_marker) or abs(_e_marker) > 1.0e9:
                _uff_param_unsafe = True
                logger.debug(
                    "UFF parameterisation marker = %.3e — soft-donor "
                    "distance pins disabled to keep this call "
                    "bit-deterministic",
                    _e_marker,
                )
        except Exception as _e_exc:
            # Energy probe failed → assume unsafe (conservative).
            _uff_param_unsafe = True
            logger.debug(
                "UFF parameterisation probe raised (%s); soft-donor "
                "distance pins disabled for safety", _e_exc,
            )

        # Apply constraints (if provided) to preserve coordination geometry
        if constraints:
            try:
                ob_constraints = pybel.ob.OBFFConstraints()

                # Baustein-5+6 Phase 3: opt-in soft-donor mode.
                # When DELFIN_UFF_SOFT_DONORS=1 AND constraints carry the
                # soft-donor meta block AND the complex class allows soft
                # donors (sigma / multi_sigma), the monodentate-donor
                # FixAtom entries are *replaced* by M-D distance pins so
                # the donor can REORIENT during UFF (substituents find a
                # tetrahedral configuration) while the M-D bond length
                # stays at the lookup-table ideal.  Metal remains FixAtom.
                # Hapto / multi_hapto / no_metal classes always fall back
                # to legacy FixAtom behaviour (helper enforces this).
                _soft_skip_donors: set = set()
                _soft_meta = constraints.get("_soft_donor_meta") if isinstance(constraints, dict) else None
                # Phase 3B per-class override (analogue of DELFIN_SIGMA_*_CLASSES):
                #   export DELFIN_UFF_SOFT_DONORS_CLASSES="sigma"
                #     → enable only for sigma (drop multi_sigma) when pool-verdict
                #     shows UFF-soft only helps sigma.
                # Empty _CLASSES env (default) → fall back to scalar
                # DELFIN_UFF_SOFT_DONORS.
                #
                # DETERMINISM FIX: this call previously passed a bare
                # ``mol`` — but ``_optimize_xyz_openbabel`` has no ``mol``
                # parameter and no module-level ``mol`` exists, so every
                # invocation raised ``NameError: name 'mol' is not
                # defined``.  That exception was swallowed by the broad
                # ``except`` around the whole constraint block, which
                # meant ``ff.SetConstraints`` was NEVER reached: every
                # constrained UFF call silently ran UNCONSTRAINED.  An
                # unconstrained UFF relax on a metal complex (metals OB
                # cannot fully parameterise) wanders to a different
                # geometry on every call, so the topology-enumerator
                # output drifted run-to-run.  ``_class_conditional_flag``
                # is fail-safe with ``mol=None`` (the ``_classify_complex_class``
                # call is wrapped in try/except, and with no ``_CLASSES``
                # env set it never touches ``mol`` at all).  The per-call
                # class guard below still applies via
                # ``_soft_meta["class_label"]`` + ``should_use_soft_donor``.
                _soft_enabled = _class_conditional_flag("DELFIN_UFF_SOFT_DONORS", None, default=1)
                # Determinism gate: when OB UFF could not parameterise
                # the metal (see ``_uff_param_unsafe`` probe above), the
                # distance-pin path is non-deterministic.  Bypass the
                # soft-donor block AND explicitly promote every donor
                # from the soft-meta into ``AddAtomConstraint`` (HARD
                # FixAtom).  Upstream ``_build_coordination_constraints_from_xyz``
                # may have dropped the donors from ``fix_atoms`` in
                # anticipation of soft-mode adding distance pins; if we
                # do not re-promote here, the donors stay FREE and the
                # uninitialised-gradient bug drags them across the
                # workspace just as if soft-mode had run.
                # Track explicit promotions so we never add the same
                # atom twice to the OB constraint set (idempotent).
                _gated_hard_added: set = set()
                # T4.1 audit 2026-05-15: class-agnostic HARD-promotion gate.
                # Default OFF — the legacy ``_sus()``-only guard remains
                # active so behaviour is byte-identical when the new env
                # is unset.  When ``DELFIN_UFF_PROBE_HARD_ALL_CLASSES=1``,
                # any soft-meta block (including hapto / multi_hapto)
                # gets HARD-promoted on _uff_param_unsafe, closing the
                # determinism gap exposed when a caller manually built
                # ``constraints`` with a hapto-class label and an empty
                # ``fix_atoms`` list (5/5 unique digests in audit Case D).
                _hard_all_classes = bool(
                    _delfin_env_int("DELFIN_UFF_PROBE_HARD_ALL_CLASSES", 0)
                )
                if _soft_meta and _uff_param_unsafe and not _soft_enabled:
                    # Soft path disabled outright; nothing to gate.
                    pass
                if _soft_enabled and _soft_meta and _uff_param_unsafe:
                    try:
                        from delfin.manta._uff_soft_donor import (
                            should_use_soft_donor as _sus,
                        )
                        _gated_cls = _soft_meta.get("class_label", "no_metal")
                        if _sus(_gated_cls) or _hard_all_classes:
                            # Promote every donor + metal recorded in
                            # the soft-meta block to HARD FixAtom.
                            # Upstream ``_build_coordination_constraints_from_xyz``
                            # already lists the metal in ``fix_atoms`` and
                            # — when ``DELFIN_UFF_RELAX_DONORS`` is the
                            # default 0 — the donors as well; we do NOT
                            # touch ``_soft_skip_donors`` here so those
                            # entries get added through the regular
                            # ``fix_atoms`` loop below.  This block
                            # covers the (rare) case where a caller
                            # built ``constraints`` manually and did not
                            # list the donors / metal: the resulting
                            # OB constraint set still pins every donor.
                            for _d in _soft_meta.get("donor_indices", []):
                                _di = int(_d)
                                if _di in _gated_hard_added:
                                    continue
                                try:
                                    ob_constraints.AddAtomConstraint(_di + 1)
                                    _gated_hard_added.add(_di)
                                except Exception:
                                    pass
                            for _m in _soft_meta.get("metal_indices", []):
                                _mi = int(_m)
                                if _mi in _gated_hard_added:
                                    continue
                                try:
                                    ob_constraints.AddAtomConstraint(_mi + 1)
                                    _gated_hard_added.add(_mi)
                                except Exception:
                                    pass
                    except Exception as _gex:
                        logger.debug(
                            "Soft-donor HARD-fallback promotion raised "
                            "(%s); relying on upstream fix_atoms list",
                            _gex,
                        )
                    logger.debug(
                        "DELFIN_UFF_SOFT_DONORS gated OFF for class=%s: "
                        "OB UFF cannot fully parameterise this metal "
                        "(soft-mode would be non-deterministic).  "
                        "Falling back to FixAtom on every donor.",
                        _soft_meta.get("class_label", "?"),
                    )
                if _soft_enabled and _soft_meta and not _uff_param_unsafe:
                    try:
                        from delfin.manta._uff_soft_donor import should_use_soft_donor  # lazy
                        _cls = _soft_meta.get("class_label", "no_metal")
                        if should_use_soft_donor(_cls):
                            _donor_set = set(int(d) for d in _soft_meta.get("donor_indices", []))
                            _pair_keys = {
                                tuple(sorted((int(m), int(d))))
                                for (m, d) in _soft_meta.get("pairs", [])
                            }
                            # Force-constant for the M-D distance pin
                            # (used only when OB version accepts the 4-arg
                            # AddDistanceConstraint signature; the 3-arg
                            # OBFFConstraints API has no force-constant).
                            _force_const = 10000.0
                            # Phase 3C per-donor-element gate
                            # (ITER-uffsoft 2026-05-13): M-C donors lose
                            # orientation under UFF because UFF has no
                            # transition-metal-bonded parameters; SOFT
                            # mode then dissociates 76% of M-C bonds in
                            # the smoke pool.  N/O/P/S/halide donors
                            # retain SOFT mode.  Carbon falls back to
                            # legacy FixAtom unless
                            # DELFIN_UFF_SOFT_DONORS_CARBON=1.
                            from delfin.manta._uff_soft_donor import (
                                should_soften_donor,  # lazy import
                            )
                            _allow_carbon_soft = bool(
                                _delfin_env_int(
                                    "DELFIN_UFF_SOFT_DONORS_CARBON", 0
                                )
                            )
                            # Replace donor FixAtom with M-D distance pin
                            # only for donor elements that pass the gate;
                            # everything else stays HARD (handled below
                            # by the regular fix_atoms loop).
                            _soft_skip_donors = set()
                            for (m_idx, d_idx) in _soft_meta.get("pairs", []):
                                try:
                                    # Look up element symbols via OB atomic
                                    # number to avoid UFF-typed strings.
                                    _m_atom = ob_mol.OBMol.GetAtom(int(m_idx) + 1)
                                    _d_atom = ob_mol.OBMol.GetAtom(int(d_idx) + 1)
                                    m_sym = pybel.ob.GetSymbol(_m_atom.GetAtomicNum())
                                    d_sym = pybel.ob.GetSymbol(_d_atom.GetAtomicNum())
                                    if not should_soften_donor(
                                        d_sym,
                                        allow_carbon=_allow_carbon_soft,
                                    ):
                                        # HARD branch: leave donor for
                                        # legacy FixAtom loop below.
                                        logger.debug(
                                            "Soft-donor gate: donor %s "
                                            "(idx=%s) treated as HARD",
                                            d_sym, d_idx,
                                        )
                                        continue
                                    d_ideal = float(_get_ml_bond_length(m_sym, d_sym))
                                    ob_constraints.AddDistanceConstraint(
                                        int(m_idx) + 1, int(d_idx) + 1, d_ideal
                                    )
                                    _soft_skip_donors.add(int(d_idx))
                                except Exception as _exc:
                                    # On any per-pair failure, fall back to
                                    # FixAtom for that donor.
                                    try:
                                        ob_constraints.AddAtomConstraint(int(d_idx) + 1)
                                    except Exception:
                                        pass
                                    _soft_skip_donors.discard(int(d_idx))
                                    logger.debug(
                                        "Soft-donor pin failed for M=%s D=%s: %s; "
                                        "falling back to FixAtom",
                                        m_idx, d_idx, _exc,
                                    )
                            # Quiet self-test/import sanity log.
                            logger.debug(
                                "DELFIN_UFF_SOFT_DONORS active: class=%s, "
                                "n_donor_pins=%d (carbon_soft=%s)",
                                _cls, len(_soft_skip_donors),
                                _allow_carbon_soft,
                            )
                    except Exception as _exc:
                        logger.debug(
                            "Soft-donor setup failed (%s); using legacy "
                            "FixAtom on donors", _exc,
                        )
                        _soft_skip_donors = set()

                # Fix atom positions (1-based indices for OB).
                # Skip donors that received an M-D distance pin above
                # (soft-mode active) or that were already promoted to
                # FixAtom by the determinism gate's HARD fallback above.
                for idx in constraints.get('fix_atoms', []):
                    if int(idx) in _soft_skip_donors:
                        continue
                    if int(idx) in _gated_hard_added:
                        continue
                    ob_constraints.AddAtomConstraint(idx + 1)
                # Distance constraints
                for idx_a, idx_b, target in constraints.get('distances', []):
                    ob_constraints.AddDistanceConstraint(
                        idx_a + 1, idx_b + 1, target
                    )
                # Angle constraints
                for idx_a, idx_b, idx_c, target in constraints.get('angles', []):
                    ob_constraints.AddAngleConstraint(
                        idx_a + 1, idx_b + 1, idx_c + 1, target
                    )
                # Torsion constraints
                for idx_a, idx_b, idx_c, idx_d, target in constraints.get('torsions', []):
                    ob_constraints.AddTorsionConstraint(
                        idx_a + 1, idx_b + 1, idx_c + 1, idx_d + 1, target
                    )
                ff.SetConstraints(ob_constraints)
            except Exception as e:
                logger.debug("OB UFF constraint setup failed, running unconstrained: %s", e)
        else:
            # DETERMINISM: explicitly install an EMPTY constraint set on
            # the (possibly fallback-shared) force field for unconstrained
            # calls.  Open Babel keeps the last ``SetConstraints`` block
            # alive at the force-field level; without this reset an
            # unconstrained optimisation that follows a constrained one in
            # the same process would silently inherit the previous call's
            # fixed atoms / distance pins, making the result depend on
            # call order.  An empty ``OBFFConstraints`` is a no-op for the
            # geometry but guarantees a clean, order-independent state.
            try:
                ff.SetConstraints(pybel.ob.OBFFConstraints())
            except Exception as e:
                logger.debug("OB UFF empty-constraint reset failed: %s", e)

        # DETERMINISM GATE (CG step): when OB-UFF could not parameterise a
        # metal in this molecule (``_uff_param_unsafe`` — the > 1e9 kcal/mol
        # energy marker above), every per-atom UFF gradient is read from
        # uninitialised heap memory and differs across process invocations.
        # The constraint set (built by ``_build_coordination_constraints_from_xyz``)
        # freezes the metal + its monodentate donor atoms, but leaves
        # *non-donor* atoms free (e.g. the carbonyl O of M-C≡O in W(CO)6 /
        # Cr(CO)6 / Mo(CO)6 — C is the donor, O is free).  ``ConjugateGradients``
        # then drags those free atoms along the garbage gradient, so the same
        # SMILES produces a different geometry every run.  This is the homoleptic
        # metal-carbonyl / unparameterised-TM determinism hole.
        #
        # Fix: when the parameterisation is unsafe the CG step provides no real
        # optimisation (the energy/gradients are meaningless), so skipping it
        # loses nothing physical and removes the only source of run-to-run
        # variation.  The geometry returned is the template-built, physically
        # sensible polyhedron the caller passed in.  Universal — gated purely on
        # the runtime energy marker, no metal allowlist, no SMILES patterns.
        # Default ON (deterministic is the correct default); revert with
        # ``DELFIN_UFF_UNSAFE_SKIP_CG=0`` to restore the legacy (random) CG run.
        _skip_unsafe_cg = (
            _uff_param_unsafe
            and _delfin_env_int("DELFIN_UFF_UNSAFE_SKIP_CG", 1)
        )
        if _skip_unsafe_cg:
            # DETERMINISM <-> CLASH HEAL (env DELFIN_UFF_FROZEN_METAL_CG,
            # default 0 → byte-identical to the b25c8b4 skip behaviour).
            #
            # The plain skip above kept the run deterministic but threw away the
            # incidental ligand relaxation that the (random) CG run used to do,
            # which had been relieving inter-ligand clashes → a confirmed
            # inter-clash regression vs golden.  When the flag is ON we instead
            # apply a clash-relief step that keeps BOTH properties:
            #
            #   mode 1 ("cg")  — metal-frozen OB-UFF: pin the unparameterised
            #       metal atom(s) via AddAtomConstraint and run CG on the rest.
            #       The garbage metal gradient can never move a fixed atom, but
            #       the surrounding ligand atoms still relax.  RISK: the M-L UFF
            #       energy terms are themselves unparameterised, so the LIGAND
            #       gradients can still read uninitialised heap → must be proven
            #       byte-identical empirically before trusting it.
            #
            #   mode 2 ("geom", DEFAULT when ON) — pure-geometric inter-fragment
            #       push-apart (``_geometric_inter_clash_relief``): deterministic
            #       BY CONSTRUCTION (no FF, no RNG, sorted-index order,
            #       canonical degenerate axis).  Mirrors the inter-ligand-clash
            #       detector's fragment/vdW contract so it improves exactly the
            #       measured metric.  Shipped default because empirical testing
            #       showed the metal-frozen CG (mode 1) can still leak
            #       non-determinism through the garbage M-L terms.
            #
            # Gated purely on the runtime energy marker — no metal allowlist,
            # no SMILES patterns.  Frozen-set = constraints['fix_atoms'] (metal +
            # pinned donors + hapto atoms), so chelate / hapto geometry is held.
            _heal_on = bool(_delfin_env_int("DELFIN_UFF_FROZEN_METAL_CG", 0))
            if _heal_on:
                _heal_mode = (os.environ.get(
                    "DELFIN_UFF_FROZEN_METAL_CG_MODE", "geom"
                ) or "geom").strip().lower()
                _frozen = []
                if isinstance(constraints, dict):
                    _frozen = list(constraints.get("fix_atoms", []) or [])
                if _heal_mode == "cg":
                    # mode 1: metal-frozen OB-UFF.  Ensure every metal atom is
                    # an AddAtomConstraint (the constraint set already freezes
                    # metal+donors via fix_atoms, but force it explicitly so the
                    # garbage metal gradient can never move the metal even if a
                    # caller passed bare constraints).
                    try:
                        _mfc = pybel.ob.OBFFConstraints()
                        _mfc_added = set()
                        for _fa in _frozen:
                            try:
                                _mfc.AddAtomConstraint(int(_fa) + 1)
                                _mfc_added.add(int(_fa))
                            except Exception:
                                pass
                        for _oba in pybel.ob.OBMolAtomIter(ob_mol.OBMol):
                            if pybel.ob.GetSymbol(_oba.GetAtomicNum()) in _METAL_SET:
                                _mi0 = _oba.GetIdx() - 1
                                if _mi0 not in _mfc_added:
                                    try:
                                        _mfc.AddAtomConstraint(_oba.GetIdx())
                                        _mfc_added.add(_mi0)
                                    except Exception:
                                        pass
                        # carry over the geometry constraints (distances/angles/
                        # torsions) so M-D bonds + CO linearity are still pinned.
                        if isinstance(constraints, dict):
                            for _a, _b, _t in constraints.get("distances", []):
                                try:
                                    _mfc.AddDistanceConstraint(_a + 1, _b + 1, _t)
                                except Exception:
                                    pass
                            for _a, _b, _c, _t in constraints.get("angles", []):
                                try:
                                    _mfc.AddAngleConstraint(
                                        _a + 1, _b + 1, _c + 1, _t
                                    )
                                except Exception:
                                    pass
                        ff.SetConstraints(_mfc)
                        ff.ConjugateGradients(steps)
                        ff.GetCoordinates(ob_mol.OBMol)
                        _hl_lines = []
                        for _oba in pybel.ob.OBMolAtomIter(ob_mol.OBMol):
                            _sym = pybel.ob.GetSymbol(_oba.GetAtomicNum())
                            _hl_lines.append(
                                f"{_sym:4s} {_oba.GetX():12.6f} "
                                f"{_oba.GetY():12.6f} {_oba.GetZ():12.6f}"
                            )
                        _xyz_heal = "\n".join(_hl_lines) + "\n"
                        logger.debug(
                            "Determinism+clash heal (mode=cg): metal-frozen "
                            "OB-UFF clash relief applied."
                        )
                        return (
                            (_xyz_heal, None) if return_energy else _xyz_heal
                        )
                    except Exception as _hl_exc:
                        logger.debug(
                            "Metal-frozen CG heal failed (%s); falling back to "
                            "geometric clash relief", _hl_exc,
                        )
                        # fall through to geometric mode
                # mode 2 (default): deterministic geometric clash relief.
                _xyz_heal = _geometric_inter_clash_relief(xyz_delfin, _frozen)
                logger.debug(
                    "Determinism+clash heal (mode=geom): deterministic "
                    "geometric inter-clash relief applied."
                )
                return (_xyz_heal, None) if return_energy else _xyz_heal
            logger.debug(
                "Skipping OB-UFF ConjugateGradients: metal unparameterised "
                "(energy marker indicates uninitialised gradients); returning "
                "input geometry to keep this call bit-deterministic."
            )
            return (xyz_delfin, None) if return_energy else xyz_delfin

        ff.ConjugateGradients(steps)
        ff.GetCoordinates(ob_mol.OBMol)

        # Capture final UFF energy (kcal/mol)
        energy: Optional[float] = None
        try:
            energy = float(ff.Energy())
        except Exception:
            energy = None

        # Convert back to DELFIN format
        opt_lines = []
        for ob_atom in pybel.ob.OBMolAtomIter(ob_mol.OBMol):
            symbol = pybel.ob.GetSymbol(ob_atom.GetAtomicNum())
            x, y, z = ob_atom.GetX(), ob_atom.GetY(), ob_atom.GetZ()
            opt_lines.append(f"{symbol:4s} {x:12.6f} {y:12.6f} {z:12.6f}")

        logger.debug("Open Babel UFF optimization completed (%d steps max)", steps)
        xyz_out = '\n'.join(opt_lines) + '\n'
        return (xyz_out, energy) if return_energy else xyz_out

    except Exception as e:
        logger.debug("Open Babel UFF optimization failed: %s", e)
        return (xyz_delfin, None) if return_energy else xyz_delfin


def _optimize_xyz_openbabel_safe(
    xyz_delfin: str,
    mol_template=None,
    smiles: Optional[str] = None,
    steps: int = 500,
    apply_template_constraints: bool = False,
    coord_constraints: Optional[Dict] = None,
) -> str:
    """Run OB-UFF, but keep original XYZ if optimization breaks topology.

    OB force fields operate on perceived connectivity from XYZ and can
    occasionally stretch/break covalent bonds for charged metal complexes.
    This wrapper accepts the optimized geometry only if it passes structural
    sanity checks against the original molecular graph.

    When ``coord_constraints`` is provided, those constraints are used
    directly (topology-enumerator path with M-D + L-M-L pinning).
    Otherwise, when ``apply_template_constraints`` is True and a
    ``mol_template`` is available, mild OCO/aromatic planarity constraints
    are passed to OB-UFF.
    """
    # d8 SQUARE-PLANAR seating (root fix): a d8 CN4 centre (Pd/Pt/Ni/Au/Rh/Ir) built by ETKDG comes
    # out TETRAHEDRAL (UFF's default for CN4 has no square-planar knowledge).  Rearrange the 4
    # monodentate donors to square-planar HERE -- BEFORE the constraint build below freezes the
    # monodentate donors at their positions -- so the freeze locks in the SQUARE and UFF preserves it.
    # (A post-UFF rearrange is undone: the freeze has already pinned the tetrahedron.)  Env-gated,
    # default OFF -> byte-identical; chelate d8 centres are skipped (backbone-coupling guard inside).
    if (mol_template is not None and RDKIT_AVAILABLE
            and os.environ.get("DELFIN_FFFREE_D8_SQ_FLATTEN", "0") == "1"):
        xyz_delfin = _flatten_d8_sq_planar_xyz(xyz_delfin, mol_template)

    constraints = coord_constraints
    if constraints is None and mol_template is not None:
        # Universal principle: EVERY UFF call on a metal complex gets
        # coordination constraints (M-D distances + L-M-L angles).
        # This prevents UFF from inventing distances for metals it has
        # no parameters for (Sc, Cd, lanthanides, etc.).
        _has_metal_in_template = any(
            a.GetSymbol() in _METAL_SET for a in mol_template.GetAtoms()
        ) if RDKIT_AVAILABLE else False
        if _has_metal_in_template:
            try:
                constraints = _build_coordination_constraints_from_xyz(
                    mol_template, xyz_delfin
                )
            except Exception as exc:
                logger.debug("Coordination constraint generation failed: %s", exc)
                constraints = None
        elif apply_template_constraints:
            try:
                constraints = _build_uff_constraints_from_template(
                    mol_template, xyz_delfin=xyz_delfin
                )
            except Exception as exc:
                logger.debug("Template constraint generation failed: %s", exc)
                constraints = None

    xyz_opt = _optimize_xyz_openbabel(xyz_delfin, steps=steps, constraints=constraints,
                                      mol_template=mol_template)
    if not xyz_opt or xyz_opt == xyz_delfin:
        return xyz_delfin

    # Post-UFF polish: project every sp2 3-coordinate atom onto the plane
    # of its three heavy non-metal neighbours.  OB-UFF torsion constraints
    # keep such atoms near-planar but do not fully enforce planarity; this
    # geometric step removes residual pyramidalisation at ring-junction
    # atoms (e.g. the fused-bicycle C in triazolothiadiazine ligands).
    if mol_template is not None and RDKIT_AVAILABLE:
        xyz_opt = _flatten_sp2_atoms_xyz(xyz_opt, mol_template)

    # POST-UFF metalloid M-D clamp (root fix): OB-UFF collapses soft, unparameterised M-metalloid
    # bonds (Sb/As/Bi/Te/Se/Ge/Sn/Pb) -- QIGFOF Ag-Sb -> 2.01 A vs ideal 2.65.  Re-snap them to ideal
    # here, BEFORE the connectivity/quality checks below so those judge the corrected geometry.
    # Env-gated, default OFF -> byte-identical.  N/O/P/S bonds (UFF handles them) are untouched.
    if (mol_template is not None and RDKIT_AVAILABLE
            and os.environ.get("DELFIN_FFFREE_METALLOID_MD_CLAMP", "0") == "1"):
        xyz_opt = _clamp_metalloid_md_xyz(xyz_opt, mol_template)

    # Fundamental check: UFF must not change the metal-donor connectivity.
    # If a donor drifted away or a non-donor collapsed onto the metal,
    # discard the UFF result and keep the original XYZ.
    if mol_template is not None and RDKIT_AVAILABLE:
        if not _verify_metal_connectivity(xyz_opt, mol_template):
            logger.debug("Discarding UFF geometry: metal connectivity changed.")
            return xyz_delfin

    if mol_template is not None and RDKIT_AVAILABLE:
        try:
            # Build baseline quality from the pre-UFF geometry so we can keep
            # constrained UFF structures that are improved even if still not
            # fully "good" by strict thresholds.
            orig_bad = False
            orig_score = float("inf")
            mol_orig = Chem.RWMol(mol_template)
            mol_orig.RemoveAllConformers()
            conf_orig = _xyz_to_rdkit_conformer(mol_orig.GetMol(), xyz_delfin)
            if conf_orig is not None:
                cid_orig = mol_orig.AddConformer(conf_orig, assignId=True)
                orig_bad = _has_bad_geometry(mol_orig.GetMol(), cid_orig)
                try:
                    orig_score = _geometry_quality_score(mol_orig.GetMol(), cid_orig)
                except Exception:
                    orig_score = float("inf")

            mol_tmp = Chem.RWMol(mol_template)
            mol_tmp.RemoveAllConformers()
            conf = _xyz_to_rdkit_conformer(mol_tmp.GetMol(), xyz_opt)
            if conf is None:
                logger.debug("Discarding UFF geometry: atom mapping failed.")
                return xyz_delfin
            cid = mol_tmp.AddConformer(conf, assignId=True)
            if _has_severe_covalent_distortion(mol_tmp.GetMol(), cid):
                logger.debug("Discarding UFF geometry: severe covalent distortion.")
                return xyz_delfin
            opt_bad = _has_bad_geometry(mol_tmp.GetMol(), cid)
            if opt_bad:
                # Keep constrained UFF if it improves relative to the original
                # geometry; otherwise fall back to the unoptimized structure.
                try:
                    opt_score = _geometry_quality_score(mol_tmp.GetMol(), cid)
                except Exception:
                    opt_score = float("inf")
                if (not orig_bad) or (opt_score >= orig_score - 1e-6):
                    logger.debug(
                        "Discarding UFF geometry: no bad-geometry improvement "
                        "(orig_bad=%s orig_score=%.3f opt_score=%.3f).",
                        orig_bad, orig_score, opt_score,
                    )
                    return xyz_delfin
        except Exception as e:
            logger.debug("Discarding UFF geometry: post-check failed (%s).", e)
            return xyz_delfin

    # For constrained, template-guided UFF in the isomer path we defer
    # topology/spurious checks to the caller (which applies the same checks
    # afterwards) to avoid double-rejecting all optimized geometries.
    if smiles and not apply_template_constraints and coord_constraints is None:
        try:
            if not _roundtrip_ring_count_ok(xyz_opt, smiles):
                logger.debug("Discarding UFF geometry: ring-count mismatch.")
                return xyz_delfin
            if not _no_spurious_bonds(xyz_opt, smiles):
                logger.debug("Discarding UFF geometry: spurious bonds.")
                return xyz_delfin
            if not _fragment_topology_ok(xyz_opt, smiles):
                logger.debug("Discarding UFF geometry: fragment topology mismatch.")
                return xyz_delfin
        except Exception:
            pass

    if mol_template is not None and RDKIT_AVAILABLE:
        try:
            xyz_opt = _snap_aromatic_rings_in_xyz(
                xyz_opt, mol_template, rms_threshold=0.05
            )
        except Exception as exc:
            logger.debug("Post-UFF aromatic-plane snap skipped: %s", exc)

    return xyz_opt
