"""Single-structure SMILES to XYZ conversion of the MANTA constructor: quick hapto previews, the sanitised path and the legacy path.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import os
import re
import threading
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_io import (
    _inject_openbabel_conformers_into_mol,
    _mol_to_xyz,
    _openbabel_generate_conformer_xyz,
)
from delfin.manta.converter_flags import (
    _EMBED_TIMEOUT,
    _PIPELINE_SEEDS,
    _TOP_LEVEL_SEEDS,
    _class_conditional_flag,
    _delfin_env_int,
    _deterministic_embed_seed,
    _deterministic_mode,
    _multihapto_etkdg_fallback_enabled,
    _resolve_top_level_seed_count,
)
from delfin.manta.embed_strategies import (
    _manual_metal_embed,
    _smiles_to_xyz_unsanitized_fallback,
    _try_multiple_strategies,
)
from delfin.manta.embed_timeout import (
    _embed_with_timeout,
    _make_random_embed_params,
)
from delfin.manta.geometry_quality import (
    _geometry_quality_score,
    _has_bad_geometry,
)
from delfin.manta.hapto_candidates import (
    _hapto_candidate_topology_ok,
    _select_best_hapto_candidate,
)
from delfin.manta.hapto_detect import (
    _apply_hapto_approximation,
    _find_hapto_groups,
    _hapto_approx_enabled,
    _hapto_failfast_error,
    _probe_hapto_groups_from_smiles,
    contains_metal,
    mol_from_smiles_rdkit,
)
from delfin.manta.hapto_scaffold import (
    _build_hapto_scaffold,
    _correct_hapto_geometry,
    _propagate_non_hapto_atoms,
)
from delfin.manta.hybrid_assembly import (
    _build_hybrid_hapto_complex,
    _enforce_donor_pi_coplanarity,
    _final_clash_resolution,
    _fix_secondary_metal_distances,
    _secondary_metal_variant_plans,
)
from delfin.manta.ligand_placement import (
    _scale_aromatic_rings_in_xyz_from_smiles,
    _scale_aromatic_rings_to_ideal_cc,
    _snap_aromatic_rings_in_xyz,
)
from delfin.manta.metal_smiles import (
    _convert_metal_bonds_to_dative,
    _denormalize_metal_smiles,
    _fix_hapto_donor_h,
    _fix_organometallic_carbon_h,
    _normalize_metal_smiles,
    _strip_h_on_metal_halogen,
)
from delfin.manta.ml_tables import (
    AllChem,
    Chem,
    OPENBABEL_AVAILABLE,
    Point3D,
    RDKIT_AVAILABLE,
    STK_AVAILABLE,
    _METAL_SET,
    _get_ml_bond_length,
    _is_metal_nitrogen_complex,
    _is_simple_organometallic,
    stk,
)
from delfin.manta.mol_prep import (
    _dearomatized_embedding_copy,
)
from delfin.manta.openbabel_optimize import (
    _optimize_xyz_openbabel_safe,
)
from delfin.manta.topology_checks import (
    _flatten_sp2_atoms_xyz,
    _has_atom_clash,
)

logger = get_logger("delfin.smiles_converter")


_HAPTO_QUICK_PREVIEW_CACHE: Dict[str, List[Tuple[str, str]]] = {}


_SMILES_TOKEN_RE = re.compile(
    r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\*|\(|\)|\.|=|#|-|\+|\\|/|:|~|@|\?"
    r"|>|<|\$|%[0-9]{2}|[0-9])")


_SMILES_ATOM_TOKEN_RE = re.compile(r"^(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\*)$")


def _pin_metal_neighbour_hydrogens(smiles):
    """Write the RDKit hydrogen count of every UNBRACKETED metal neighbour into the SMILES.

    ``[Pt](Cl)(Cl)(N)N``: RDKit gives each N two implicit hydrogens (valence 3, one
    bond to Pt).  Implicit hydrogens are not stored on the atom, they are recomputed
    from its valence -- and preparation paths mark atoms bonded to a metal
    ``NoImplicit`` (so a dative donor does not gain a spurious H) or cleave the M-L
    bond, both of which change what "implicit" comes to.  The unbracketed N lost both
    hydrogens on some paths and kept them on others, so one manifold held frames of 5
    and of 9 atoms, and with a counter-ion (``.Cl``) of 7 and of 11.

    The count is fixed where the user wrote it: each such atom token becomes the
    bracket atom RDKit itself reads it as (``N`` -> ``[NH2]``), so from here on the
    hydrogens are explicit and no downstream path can drop them.  Only unbracketed
    atoms that border a metal AND carry at least one implicit H are rewritten; any
    other SMILES (a bracket atom already states its H count) is returned as the
    identical string.  Never raises; on any doubt the input is returned unchanged.
    """
    if not smiles or not isinstance(smiles, str) or not RDKIT_AVAILABLE:
        return smiles
    try:
        tokens = _SMILES_TOKEN_RE.findall(smiles)
        if "".join(tokens) != smiles:
            return smiles
        atom_pos = [i for i, t in enumerate(tokens) if _SMILES_ATOM_TOKEN_RE.match(t)]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                return smiles
            mol.UpdatePropertyCache(strict=False)
        if mol.GetNumAtoms() != len(atom_pos):
            return smiles
        changed = False
        for atom in mol.GetAtoms():
            pos = atom_pos[atom.GetIdx()]
            tok = tokens[pos]
            if tok.startswith("[") or tok == "*":
                continue
            if not any(n.GetSymbol() in _METAL_SET for n in atom.GetNeighbors()):
                continue
            n_h = int(atom.GetNumImplicitHs())
            if n_h <= 0:
                continue
            tokens[pos] = f"[{tok}H{n_h if n_h > 1 else ''}]"
            changed = True
        return "".join(tokens) if changed else smiles
    except Exception:
        return smiles


def smiles_to_xyz_quick_hapto_previews(
    smiles: str,
    hapto_approx: Optional[bool] = None,
) -> List[Tuple[str, str]]:
    """Return extra quick-convert preview structures for hapto-specific builders."""
    smiles = _pin_metal_neighbour_hydrogens(smiles)   # same key as smiles_to_xyz_quick
    _ = _hapto_approx_enabled(hapto_approx)
    return list(_HAPTO_QUICK_PREVIEW_CACHE.get(smiles, []))


def smiles_to_xyz_quick(
    smiles: str,
    hapto_approx: Optional[bool] = None,
) -> Tuple[Optional[str], Optional[str]]:
    """Fast single-conformer SMILES → XYZ (DELFIN format, no header).

    Uses the same strategy as the pre-Avogadro smiles_to_xyz: single
    ETKDGv3 embedding per strategy (no OB restarts, no multi-seed loop).
    For metal complexes delegates to _try_multiple_strategies which tries
    stk → RDKit → unsanitized → no-valence-check → OB in order and
    returns as soon as one succeeds.  Returns ``(xyz, error)``.
    """
    # Unbracketed metal neighbours keep the H RDKit gives them (see
    # _pin_metal_neighbour_hydrogens); identical string for every other SMILES.
    smiles = _pin_metal_neighbour_hydrogens(smiles)
    if not RDKIT_AVAILABLE:
        return None, "RDKit not available"

    has_metal = contains_metal(smiles)
    hapto_mode = _hapto_approx_enabled(hapto_approx)

    if has_metal:
        hapto_groups = _probe_hapto_groups_from_smiles(smiles)
        if hapto_groups:
            if not hapto_mode:
                return None, _hapto_failfast_error(hapto_groups)
            # Experimental eta/hapto fallback: reuse the main converter with
            # approximation enabled and no UFF to keep quick-mode behavior.
            return smiles_to_xyz(smiles, apply_uff=False, hapto_approx=True)
        return _try_multiple_strategies(smiles)

    # Non-metal: single ETKDG embedding
    mol, note = mol_from_smiles_rdkit(smiles, allow_metal=False)
    if mol is None:
        return None, f"Could not parse SMILES: {note}"
    try:
        mol = Chem.AddHs(mol)
    except Exception:
        pass
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    result = AllChem.EmbedMolecule(mol, params)
    if result != 0:
        result = AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
    if result != 0:
        return None, f"Quick embedding failed for SMILES: {smiles[:60]}"
    return _mol_to_xyz(mol), None


def _topology_hard_gate_check(
    xyz: Optional[str], smiles: str,
) -> Tuple[Optional[str], Optional[str]]:
    """Env-flag-gated runtime topology-invariant validation.

    Modes (DELFIN_TOPOLOGY_HARD_GATE):
        0 — gate off, return (xyz, None) unchanged (default).
        1 — log-only: warn on violation but still return xyz.
        2+ — reject: return (None, "topology_hard_gate: ...") on violation.

    Pillar 1 of the Hybrid-Path. Per nature_project/15_HYBRID_PATH_FINAL.md.
    """
    if xyz is None:
        return None, None
    # Universal aryl-ring-size finalizer (env-gated, default OFF, byte-identical
    # when unset).  Runs here because this is the single post-processor every
    # meaningful build path returns through, including the metal / multi-strategy
    # builders that bypass _snap_aromatic_rings_in_xyz.  Strict no-op unless the
    # SMILES maps 1:1 onto the geometry.
    xyz = _scale_aromatic_rings_in_xyz_from_smiles(xyz, smiles)
    mode = _delfin_env_int("DELFIN_TOPOLOGY_HARD_GATE", 0)
    if mode == 0:
        return xyz, None
    try:
        from delfin.topology_hard_gate import validate_topology_invariant
        result = validate_topology_invariant(xyz, smiles)
        if not result.passed:
            reasons = ",".join(sorted({v.kind for v in result.violations}))
            if mode == 1:
                logger.warning(
                    "topology_hard_gate (log-only mode=1): %s [smiles=%s]",
                    reasons, smiles[:80],
                )
                return xyz, None
            return None, f"topology_hard_gate: {reasons}"
    except Exception as exc:
        logger.debug("topology_hard_gate skipped: %s", exc)
    return xyz, None


def _try_ensemble_router(smiles: str) -> Optional[Tuple[Optional[str], Optional[str]]]:
    """Phase 2: try class-aware ensemble router if DELFIN_ENSEMBLE_ROUTER>0.

    Returns:
        None — router disabled or no specialist routed (caller should
               fall through to legacy conversion path).
        (xyz, err) — specialist produced a result; caller should return this.

    Modes:
        0 (default): no-op, returns None.
        1: try specialist, fall through if no match or specialist fails.
        2: strict — specialist required, returns error tuple if no match.
    """
    mode = _delfin_env_int("DELFIN_ENSEMBLE_ROUTER", 0)
    if mode == 0:
        return None
    try:
        from delfin.class_modules.registry import register_all_specialists
        from delfin.class_modules.smiles_name_lookup import get_name_for_smiles
        from delfin.ensemble_router import route
        register_all_specialists()
        specialist = route(smiles)
        if specialist is None:
            if mode >= 2:
                from delfin.classify import classify_smiles
                f = classify_smiles(smiles)
                return (None,
                    f"ensemble_router strict: no specialist for "
                    f"({f.coord_class}, {f.metal_block})")
            return None  # fall through to legacy
        # D2: pass name so archive_readback can short-circuit to cached
        # champion XYZ when SMILES is in a known pool.
        smi_name = get_name_for_smiles(smiles)
        xyz, err = specialist.convert(smiles, name=smi_name)
        if xyz is None and mode < 2:
            return None  # let legacy retry
        # Apply topology gate to specialist output
        gated_xyz, gate_err = _topology_hard_gate_check(xyz, smiles)
        if gate_err and mode < 2:
            return None  # let legacy retry
        return (gated_xyz, gate_err or err)
    except Exception as exc:
        logger.debug("ensemble_router skipped due to exception: %s", exc)
        return None  # fall through


def _uff_seat_aromatic_bonds(mol, max_iters: int = 200,
                             force_constant: float = 1.0e4) -> None:
    """UFF-minimise ``mol`` with every AROMATIC ring bond harmonically pinned to
    its delocalised CCDC target length (DELFIN_FFFREE_AROM_SEAT root seat).

    This is the metal-FREE lever: metal-free molecules return from the isomer
    pool before the post-hoc aromatic corrector ever runs, so the only place a
    free organic aromatic ring can be seated at its delocalised length is the UFF
    starting-structure minimisation.  A strong distance constraint holds each
    aromatic bond at its target while UFF relaxes the ring angles and the
    substituents coherently (so junction angles co-adapt rather than the ring
    staying at the covalent-sum length).  Falls back to a plain UFF minimisation
    when the force field cannot be built.  Mutates ``mol``'s conformer in place,
    exactly like ``AllChem.UFFOptimizeMolecule``."""
    from delfin.fffree.aromatic_bond_targets import aromatic_target
    ff = AllChem.UFFGetMoleculeForceField(mol)
    if ff is None:                                   # unparametrised → plain UFF
        AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
        return
    for bond in mol.GetBonds():
        if not bond.GetIsAromatic():
            continue
        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()
        tgt = aromatic_target(a.GetSymbol(), b.GetSymbol())
        if tgt is None:                              # untabulated pair → leave to UFF
            continue
        ff.AddDistanceConstraint(a.GetIdx(), b.GetIdx(), tgt, tgt, force_constant)
    ff.Initialize()
    ff.Minimize(maxIts=max_iters)


def smiles_to_xyz(
    smiles: str,
    output_path: Optional[str] = None,
    apply_uff: bool = True,
    hapto_approx: Optional[bool] = None,
    deterministic: bool = True,
) -> Tuple[Optional[str], Optional[str]]:
    """Convert SMILES string to XYZ coordinates using RDKit.

    Uses RDKit's ETKDG (Experimental Torsion Knowledge Distance Geometry) method
    for generating reasonable 3D conformers. This is suitable for initial geometries
    that will be further optimized with GOAT/xTB/ORCA.

    Supports coordination bonds (>) in SMILES for metal complexes.

    Args:
        smiles: SMILES string to convert
        output_path: Optional path to write XYZ file
        apply_uff: If True, apply UFF refinement (RDKit/Open Babel) where available.
        deterministic: If True (default), the metal conformer pool skips the
            Open Babel ``make3D``/rotor injection.  Open Babel's coordinate
            builder is seeded from the wall clock with no Python-side seed
            hook, so including it makes the output non-reproducible run to
            run.  Set to ``False`` only when extra geometric diversity is
            wanted and reproducibility is not required (mirrors the
            ``deterministic`` flag of :func:`smiles_to_xyz_isomers`).

    Returns:
        Tuple of (xyz_content, error_message)
        - xyz_content: XYZ format string if successful, None on error
        - error_message: Error description if failed, None on success
    """
    # Unbracketed metal neighbours keep the H RDKit gives them (see
    # _pin_metal_neighbour_hydrogens); identical string for every other SMILES.
    smiles = _pin_metal_neighbour_hydrogens(smiles)
    if not RDKIT_AVAILABLE:
        error = "RDKit is not installed. Install with: pip install rdkit"
        logger.error(error)
        return None, error

    # Phase 2: try class-aware ensemble router (env-flag-gated, default OFF)
    _router_result = _try_ensemble_router(smiles)
    if _router_result is not None:
        if output_path and _router_result[0]:
            Path(output_path).write_text(_router_result[0], encoding='utf-8')
        return _router_result

    import numpy as np  # local import: keeps module import cheap

    try:
        has_metal = contains_metal(smiles)
        hapto_mode = _hapto_approx_enabled(hapto_approx)
        hapto_groups = _probe_hapto_groups_from_smiles(smiles) if has_metal else []
        _HAPTO_QUICK_PREVIEW_CACHE[smiles] = []
        legacy_hapto_xyz: Optional[str] = None
        if hapto_groups and not hapto_mode:
            error = _hapto_failfast_error(hapto_groups)
            logger.warning(error)
            return None, error
        if has_metal and hapto_mode and hapto_groups:
            # For single-hapto systems, try legacy path as backup.
            # For multi-hapto (>=2 groups), skip the slow legacy path.
            if len(hapto_groups) <= 1:
                legacy_xyz, legacy_err = _try_multiple_strategies(
                    smiles, deterministic=deterministic
                )
                if legacy_err is None and legacy_xyz:
                    legacy_hapto_xyz = legacy_xyz
                    logger.info(
                        "Hapto detected: collected legacy multi-strategy candidate "
                        "for universal candidate selection."
                    )
        method = None

        # Cluster compounds (borane/carborane cages, etc.): extremely dense
        # ring-closure topologies cause ETKDG to hang.  Route directly to
        # the manual metal embed fallback which places atoms geometrically.
        if has_metal and not hapto_groups:
            _cluster_mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if _cluster_mol is not None:
                _nb = _cluster_mol.GetNumBonds()
                _na = _cluster_mol.GetNumAtoms()
                if _na > 0 and _nb > 2 * _na:
                    logger.info(
                        "Detected cluster topology (%d bonds / %d atoms = %.1f), "
                        "using manual metal embed", _nb, _na, _nb / _na)
                    xyz_content, manual_err = _manual_metal_embed(smiles)
                    if xyz_content:
                        if apply_uff:
                            xyz_content = _optimize_xyz_openbabel_safe(
                                xyz_content, mol_template=_cluster_mol
                            )
                        if output_path:
                            Path(output_path).write_text(xyz_content, encoding='utf-8')
                        return xyz_content, None

        # For metal-nitrogen coordination complexes (both neutral and charged notation),
        # use the multi-strategy approach that tries multiple parsing methods
        if _is_metal_nitrogen_complex(smiles) and not hapto_groups:
            logger.info("Detected metal-nitrogen complex, using multi-strategy approach")
            return _try_multiple_strategies(
                smiles, output_path, deterministic=deterministic
            )

        # Parse SMILES - try stk first for metal complexes, then RDKit.
        # IMPORTANT: Prefer the original SMILES before charge-normalized
        # variants to avoid introducing artificial oxidation states
        # (e.g., neutral Cu complexes rewritten as Cu+2).
        mol = None
        normalized_smiles = _normalize_metal_smiles(smiles)
        if has_metal and STK_AVAILABLE and len(hapto_groups) <= 1:
            # Skip STK for multi-hapto systems (slow and usually fails).
            # Guard with timeout: stk internally runs ETKDG which can hang.
            _stk_mol_holder = [None]

            def _stk_parse():
                try:
                    bb = stk.BuildingBlock(smiles)
                    _stk_mol_holder[0] = bb.to_rdkit_mol()
                except Exception:
                    _stk_mol_holder[0] = None

            _stk_t = threading.Thread(target=_stk_parse, daemon=True)
            _stk_t.start()
            _stk_t.join(timeout=_EMBED_TIMEOUT)
            if _stk_t.is_alive():
                logger.info("stk conversion timed out after %.1fs, falling back to RDKit", _EMBED_TIMEOUT)
                mol = None
                method = None
            elif _stk_mol_holder[0] is not None:
                mol = _stk_mol_holder[0]
                method = "stk"
            else:
                logger.info("stk conversion failed, falling back to RDKit")
                mol = None
                method = None

        if mol is None:
            mol, rdkit_note = mol_from_smiles_rdkit(smiles, allow_metal=has_metal)
            method = "RDKit" if (mol is not None and rdkit_note is None) else (
                f"RDKit ({rdkit_note})" if mol is not None else None
            )
            if mol is None:
                # Try normalized charged form for neutral metal SMILES
                if normalized_smiles:
                    mol2, rdkit_note2 = mol_from_smiles_rdkit(normalized_smiles, allow_metal=True)
                    if mol2 is not None:
                        mol = mol2
                        method = "RDKit (normalized metal SMILES)"
                        rdkit_note = None
                    else:
                        rdkit_note = rdkit_note2 or rdkit_note

                # Try denormalized (neutral) SMILES as fallback
                if mol is None:
                    denormalized_smiles = _denormalize_metal_smiles(smiles)
                    if denormalized_smiles:
                        mol3, rdkit_note3 = mol_from_smiles_rdkit(denormalized_smiles, allow_metal=True)
                        if mol3 is not None:
                            mol = mol3
                            method = "RDKit (denormalized)"
                            rdkit_note = None

                if rdkit_note and ("Explicit valence" in rdkit_note or "kekulize" in rdkit_note):
                    # For hapto complexes: try unsanitized mol (no valence check)
                    # so we can still go through the hapto correction path.
                    if hapto_groups and hapto_mode:
                        try:
                            mol_unsan = Chem.MolFromSmiles(smiles, sanitize=False)
                            if mol_unsan is not None:
                                mol_unsan.UpdatePropertyCache(strict=False)
                                mol = mol_unsan
                                method = "RDKit (unsanitized for hapto)"
                                rdkit_note = None
                                logger.info(
                                    "Using unsanitized mol for hapto path "
                                    "(valence error bypassed)")
                        except Exception:
                            pass  # fall through to legacy

                    if mol is None:
                        legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                        if legacy_err is None and legacy_xyz:
                            if hapto_groups and hapto_mode:
                                # Save as fallback but don't return yet —
                                # try hapto path first via unsanitized mol
                                legacy_hapto_xyz = legacy_xyz
                            else:
                                if output_path:
                                    Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                                return legacy_xyz, None
                # Last resort: try multi-strategy approach for metal complexes
                if has_metal and mol is None:
                    if hapto_groups and hapto_mode:
                        error = (
                            "Failed to parse hapto complex for experimental approximation. "
                            "Try simplifying the SMILES around eta-coordination."
                        )
                        logger.error(error)
                        return None, error
                    logger.info("Trying multi-strategy fallback for unparseable metal SMILES")
                    return _try_multiple_strategies(
                        smiles, output_path, deterministic=deterministic
                    )
                error = f"Failed to parse SMILES string: {rdkit_note}"
                logger.error(error)
                return None, error

        if has_metal and hapto_mode:
            try:
                parsed_hapto = _find_hapto_groups(mol)
                if parsed_hapto:
                    mol, n_removed = _apply_hapto_approximation(mol, parsed_hapto)
                    logger.info(
                        "Applied experimental hapto approximation: "
                        "%d group(s), %d bond(s) converted to dative.",
                        len(parsed_hapto), n_removed,
                    )
            except Exception as e:
                logger.debug("Hapto approximation failed, continuing with original graph: %s", e)

        # For metal complexes: convert bonds to dative and recalculate hydrogens
        # This fixes the issue where metal coordination bonds are counted towards ligand valence
        # Each step is independent so failure of one (e.g., RemoveHs on tetracoordinate B)
        # doesn't block dative conversion or AddHs.
        if has_metal:
            # Step 1: Remove existing explicit H atoms (stk may have added incorrect ones)
            try:
                mol = Chem.RemoveHs(mol)
            except Exception as e:
                logger.debug(f"RemoveHs skipped (non-standard valence, e.g. tetracoordinate B): {e}")

            # Step 2: Convert single bonds to metals to dative bonds
            try:
                mol = _convert_metal_bonds_to_dative(mol)
                mol.UpdatePropertyCache(strict=False)
            except Exception as e:
                logger.debug(f"Dative conversion skipped: {e}")

            # Step 3: Add hydrogens with correct valence calculation
            try:
                if _is_simple_organometallic(smiles):
                    mol = Chem.AddHs(mol, addCoords=True)
                    mol = _strip_h_on_metal_halogen(mol)
                    mol = _fix_organometallic_carbon_h(mol)
                else:
                    mol = Chem.AddHs(mol, addCoords=True)
                if hapto_mode:
                    mol = _fix_hapto_donor_h(mol)
            except Exception as e:
                logger.warning(f"Could not add hydrogens to metal complex: {e}")
                # Fallback for valence errors: try unsanitized path
                if "Explicit valence" in str(e):
                    legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                    if legacy_err is None and legacy_xyz:
                        if output_path:
                            Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                            logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                        return legacy_xyz, None

            # Step 4: Fix atoms with non-standard valence (e.g., tetracoordinate B
            # in pyrazolylborate/scorpionate ligands) to prevent embedding failures.
            # B with 4 bonds is chemically B- (borate) - set charge so RDKit
            # accepts valence 4 during embedding.
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'B' and atom.GetDegree() >= 4:
                    atom.SetNoImplicit(True)
                    if atom.GetFormalCharge() == 0:
                        atom.SetFormalCharge(-1)
        else:
            # Non-metal molecules: just add hydrogens normally
            try:
                mol = Chem.AddHs(mol, addCoords=True)
            except Exception as e:
                if "Explicit valence" in str(e):
                    legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
                    if legacy_err is None and legacy_xyz:
                        if output_path:
                            Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                            logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                        return legacy_xyz, None
                pass

        # Hapto path: ETKDG for topology preservation + sphere correction
        # for correct hapto geometry.  Sphere scaffold only as last fallback.
        if has_metal and hapto_mode and hapto_groups:
            # Re-identify hapto groups on the current mol (atom indices may have
            # shifted due to RemoveHs/AddHs/dative bond conversion above).
            hapto_groups = _find_hapto_groups(mol)
            if not hapto_groups:
                hapto_groups = _probe_hapto_groups_from_smiles(smiles)
            preview_candidates: List[Tuple[str, str]] = []
            hapto_candidate_mols: List[Tuple[str, object]] = []

            # Primary hybrid path: analytical hapto scaffold plus rigidly
            # aligned ligand fragments. This avoids global ETKDG on the full
            # metal graph where multi-hapto systems are most brittle.
            try:
                # ---- Wave-5 MULTIHAPTO_SIMPLE_PATH (BEGIN) -----------------
                # Wave-4/Agent-3 root-cause: ``_build_multimetal_hapto_sequential``
                # (called via ``_build_hybrid_hapto_complex``) produces hapto
                # geometry with spurious extra-bonds on the secondary metal's
                # coordination sphere -- Step 3 (``donor_indices``) skips all
                # metals so SMILES-declared M-M sigma bonds (e.g. Sn-Ir) are
                # never enforced.  Downstream ``_select_best_hapto_candidate``
                # picks those broken seq-builder candidates over ETKDG-seed
                # and sphere-scaffold candidates because they have correct
                # eta-distance, masking the M-M topology breakage.
                #
                # SIMPLE_PATH bypasses the variant-plan loop +
                # ``_build_hybrid_hapto_complex`` entirely for the
                # ``multi_hapto`` class.  ETKDG-seed embedding (below) and
                # sphere-scaffold fallback still run and feed
                # ``hapto_candidate_mols``; ``_select_best_hapto_candidate``
                # then picks the best from those.
                #
                # Phase-1.5 wire-in (2026-05-15): class-conditional default-ON
                # for ``multi_hapto`` only.  The class is 0/22 broken on the
                # v2-final deterministic baseline; SIMPLE_PATH is the
                # single-flag rescue and its scope is class-gated so
                # mono-hapto (729 SMILES, 95.6% pass) stays bit-identical.
                #
                # Resolution precedence (see ``_class_conditional_flag``):
                #   1. DELFIN_MULTIHAPTO_SIMPLE_PATH_CLASSES set → use that list.
                #   2. DELFIN_MULTIHAPTO_SIMPLE_PATH=0 → fully disabled
                #      (per-class rollback), restores pre-wire-in behaviour.
                #   3. DELFIN_MULTIHAPTO_SIMPLE_PATH=1 → enabled on every class
                #      (operator escape hatch for cross-class experiments).
                #   4. env unset → enabled iff class is in ``default_classes``,
                #      i.e. only for ``multi_hapto``.
                _simple_path_active = _class_conditional_flag(
                    "DELFIN_MULTIHAPTO_SIMPLE_PATH", mol, default=0,
                    default_classes=("multi_hapto",),
                )
                if _simple_path_active:
                    logger.info(
                        "multi-hapto SIMPLE_PATH active: skipping "
                        "_secondary_metal_variant_plans + "
                        "_build_hybrid_hapto_complex (ETKDG-seeds + "
                        "scaffold-fallback still run)"
                    )
                    hybrid_variant_plans: List[Tuple[str, Dict[int, int]]] = []
                else:
                    hybrid_variant_plans = _secondary_metal_variant_plans(mol, hapto_groups)
                # ---- Wave-5 MULTIHAPTO_SIMPLE_PATH (END) -------------------
                seen_hybrid_xyz: set = set()
                for plan_idx, (hybrid_label, variant_plan) in enumerate(hybrid_variant_plans):
                    hybrid_mol = _build_hybrid_hapto_complex(
                        mol,
                        hapto_groups,
                        preview_store=preview_candidates if plan_idx == 0 else None,
                        secondary_variant_plan=variant_plan,
                    )
                    if hybrid_mol is None:
                        continue
                    try:
                        hybrid_xyz_key = "\n".join(
                            line.strip()
                            for line in _mol_to_xyz(hybrid_mol).splitlines()
                            if line.strip()
                        )
                    except Exception:
                        hybrid_xyz_key = ""
                    if hybrid_xyz_key and hybrid_xyz_key in seen_hybrid_xyz:
                        continue
                    if hybrid_xyz_key:
                        seen_hybrid_xyz.add(hybrid_xyz_key)
                    hapto_candidate_mols.append((hybrid_label, hybrid_mol))
                if seen_hybrid_xyz:
                    logger.info(
                        "Hybrid hapto scaffold/fragment builder produced %d ranked candidate(s)",
                        len(seen_hybrid_xyz),
                    )
            except Exception as e:
                logger.debug("Hybrid hapto scaffold/fragment builder failed: %s", e)
            _HAPTO_QUICK_PREVIEW_CACHE[smiles] = list(preview_candidates)

            # Fix O- atoms with double bonds that cause valence errors
            # during ETKDG embedding (e.g. C=[O-] in metal chelates).
            # Temporarily neutralize for embedding; formal charges don't
            # affect distance geometry.
            _fixed_o_atoms = []
            for _oa in mol.GetAtoms():
                if (_oa.GetSymbol() == 'O'
                        and _oa.GetFormalCharge() == -1
                        and any(b.GetBondType() == Chem.BondType.DOUBLE
                                for b in _oa.GetBonds())):
                    _oa.SetFormalCharge(0)
                    _fixed_o_atoms.append(_oa.GetIdx())
            if _fixed_o_atoms:
                try:
                    mol.UpdatePropertyCache(strict=False)
                except Exception:
                    pass

            # SECONDARY: Two-phase ETKDG embedding
            # Phase 1: ETKDG + analytical hapto correction → fix hapto atoms
            # Phase 2: candidate selection via topology/geometry filter.
            etkdg_seeds = _PIPELINE_SEEDS[:8]
            for seed in etkdg_seeds:
                try:
                    mol_try = Chem.Mol(mol)
                    params = AllChem.ETKDGv3()
                    params.randomSeed = seed
                    params.useRandomCoords = True
                    params.enforceChirality = False
                    r = _embed_with_timeout(mol_try, params)
                    if r != 0:
                        r = _embed_with_timeout(
                            mol_try,
                            _make_random_embed_params(seed),
                        )
                    if r != 0:
                        continue
                    cid = int(r)

                    # Phase 1: Analytical hapto correction
                    try:
                        _correct_hapto_geometry(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    # Phase 1.5: Fix sigma metal-ligand distances.
                    # ETKDG treats M-L bonds like organic bonds (~1.5Å)
                    # which is too short.  Scale entire ligand fragments
                    # rigidly to correct M-L distances.
                    try:
                        from collections import deque as _deque
                        _conf = mol_try.GetConformer(cid)
                        hapto_atom_set = set()
                        metal_set_idx = set()
                        for _mi, catoms in hapto_groups:
                            metal_set_idx.add(_mi)
                            for _ca in catoms:
                                hapto_atom_set.add(_ca)
                        for _ai in range(mol_try.GetNumAtoms()):
                            if mol_try.GetAtomWithIdx(_ai).GetSymbol() in _METAL_SET:
                                metal_set_idx.add(_ai)

                        for _ai in metal_set_idx:
                            _atom = mol_try.GetAtomWithIdx(_ai)
                            m_sym = _atom.GetSymbol()
                            mp = _conf.GetAtomPosition(_ai)
                            m_pos = np.array([mp.x, mp.y, mp.z])
                            for _nbr in _atom.GetNeighbors():
                                _ni = _nbr.GetIdx()
                                if _ni in hapto_atom_set:
                                    continue  # handled by Phase 1
                                if _ni in metal_set_idx:
                                    continue  # M-M bonds
                                l_sym = _nbr.GetSymbol()
                                lp = _conf.GetAtomPosition(_ni)
                                l_pos = np.array([lp.x, lp.y, lp.z])
                                cur_d = float(np.linalg.norm(l_pos - m_pos))
                                if cur_d < 1e-8:
                                    continue
                                target_d = float(
                                    _get_ml_bond_length(m_sym, l_sym))
                                if abs(cur_d - target_d) / target_d < 0.15:
                                    continue
                                # BFS to find entire fragment attached to
                                # this donor (not crossing metals)
                                frag = set()
                                q = _deque([_ni])
                                frag.add(_ni)
                                while q:
                                    cur = q.popleft()
                                    for fn in mol_try.GetAtomWithIdx(cur).GetNeighbors():
                                        fi = fn.GetIdx()
                                        if (fi not in frag
                                                and fi not in metal_set_idx
                                                and fi not in hapto_atom_set):
                                            frag.add(fi)
                                            q.append(fi)
                                # Rigid translation of entire fragment
                                delta = (target_d - cur_d) / cur_d * (
                                    l_pos - m_pos)
                                for fi in frag:
                                    fp = _conf.GetAtomPosition(fi)
                                    fv = np.array([fp.x, fp.y, fp.z])
                                    new_fv = fv + delta
                                    _conf.SetAtomPosition(
                                        fi, Point3D(float(new_fv[0]),
                                                    float(new_fv[1]),
                                                    float(new_fv[2])))
                    except Exception:
                        pass

                    # Phase 1.6: BFS propagation for non-hapto atoms
                    try:
                        _propagate_non_hapto_atoms(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    try:
                        _enforce_donor_pi_coplanarity(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    # Fix secondary metal positions AFTER coplanarity
                    try:
                        _fix_secondary_metal_distances(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    try:
                        _final_clash_resolution(mol_try, cid, hapto_groups)
                    except Exception:
                        pass

                    hapto_candidate_mols.append((f"etkdg-seed-{seed}", mol_try))
                except Exception:
                    continue

            # FALLBACK: Sphere-based constructive scaffold builder
            try:
                mol_scaffold = Chem.Mol(mol)
                if _build_hapto_scaffold(mol_scaffold, hapto_groups):
                    try:
                        _enforce_donor_pi_coplanarity(mol_scaffold, 0, hapto_groups)
                    except Exception:
                        pass
                    try:
                        _final_clash_resolution(mol_scaffold, 0, hapto_groups)
                    except Exception:
                        pass
                    hapto_candidate_mols.append(("scaffold", mol_scaffold))
            except Exception as e:
                logger.debug("Sphere scaffold failed: %s", e)

            best_mol = _select_best_hapto_candidate(
                smiles,
                hapto_groups,
                hapto_candidate_mols,
                apply_uff=apply_uff,
            )

            if best_mol is None:
                if legacy_hapto_xyz:
                    legacy_hapto_xyz = _scale_aromatic_rings_in_xyz_from_smiles(
                        legacy_hapto_xyz, smiles
                    )
                    if output_path:
                        Path(output_path).write_text(legacy_hapto_xyz, encoding='utf-8')
                    return legacy_hapto_xyz, None
                return None, "Hapto-approx embedding failed"

            # ---- Iter-23 MULTIHAPTO_ETKDG_FALLBACK (BEGIN) -----------------
            # WAVE7_Q: when the analytical scaffold yields only topology-broken
            # candidates, re-enable the ETKDG-seed42 fallback-as-feature that
            # 81f8a1f had by accident (killed by fdeb9cb).  Zero-downside: only
            # swap when the fallback XYZ is *strictly* topology-OK; otherwise
            # keep the scaffold result unchanged.  ``mol`` here is still the
            # original parsed mol (reassigned to best_mol below) so class
            # detection is reliable.
            try:
                _i23_enabled = _multihapto_etkdg_fallback_enabled(mol)
            except Exception:
                _i23_enabled = False
            if _i23_enabled:
                try:
                    _i23_best_xyz = _mol_to_xyz(best_mol)
                except Exception:
                    _i23_best_xyz = None
                _i23_best_ok = bool(_i23_best_xyz) and _hapto_candidate_topology_ok(
                    _i23_best_xyz, smiles, mol=best_mol, conf_id=0,
                    hapto_groups=hapto_groups,
                )
                if not _i23_best_ok:
                    try:
                        _i23_fb_xyz, _ = _try_multiple_strategies(smiles)
                    except Exception:
                        _i23_fb_xyz = None
                    if _i23_fb_xyz and _hapto_candidate_topology_ok(
                        _i23_fb_xyz, smiles, mol=None, conf_id=0,
                        hapto_groups=hapto_groups,
                    ):
                        logger.info(
                            "Iter-23 multihapto ETKDG-seed42 fallback recovered "
                            "topology (%s)", smiles[:48],
                        )
                        _i23_fb_xyz = _scale_aromatic_rings_in_xyz_from_smiles(
                            _i23_fb_xyz, smiles
                        )
                        if output_path:
                            Path(output_path).write_text(_i23_fb_xyz, encoding='utf-8')
                        return _i23_fb_xyz, None
            # ---- Iter-23 MULTIHAPTO_ETKDG_FALLBACK (END) -------------------

            mol = best_mol
            # Universal aryl-ring-size finalizer (env-gated, default OFF).  The
            # hapto / multi-hapto emit returns here without passing through
            # _topology_hard_gate_check.  Operate on best_mol's conformer so we
            # use its reliable bond topology (atoms match the geometry 1:1),
            # then emit.  M-D safe: rings reaching the metal are frozen.
            if _delfin_env_int("DELFIN_FFFREE_ARYL_RING_SIZE", 0):
                try:
                    if mol.GetNumConformers() > 0:
                        _scale_aromatic_rings_to_ideal_cc(mol, 0)
                except Exception as _exc:
                    logger.debug("aryl-ring-size (hapto emit) skipped: %s", _exc)
            xyz_content = _mol_to_xyz(mol)
            if output_path:
                Path(output_path).write_text(xyz_content, encoding='utf-8')
            return xyz_content, None

        # Generate 3D coordinates using a hybrid OB+RDKit conformer pool.
        # For metal complexes: OB WeightedRotorSearch in 3 independent restarts
        # (non-deterministic diversity) + RDKit ETKDG with
        # 12 diverse fixed seeds → up to ~500 conformers total.
        # Best geometry is selected by _geometry_quality_score.
        #
        # For large molecules (>50 heavy atoms): skip the expensive OB/ETKDG
        # multi-conformer pipeline.  OB's SystematicRotorSearch is
        # combinatorially explosive and holds the GIL (no thread-based timeout).
        # ETKDG can also hang on complex metal ring systems.  These are initial
        # geometries destined for GOAT/xTB anyway.
        _n_heavy_atoms = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
        _skip_expensive_pipeline = _n_heavy_atoms > 50
        result = -1
        if has_metal and not _skip_expensive_pipeline:
            all_conf_ids: List[int] = []

            # --- OB conformers (Avogadro-equivalent pipeline) ---
            # Skipped in deterministic mode: OpenBabel's
            # SystematicRotorSearch holds the GIL and can hang
            # indefinitely on heavily-decorated complexes (Pt-phosphine-
            # tBu, V-phenolate-THF, Mn-phosphazene-Te, etc.) where the
            # rotor count is small enough to bypass the rotor-cap guard
            # in _openbabel_generate_conformer_xyz but the search tree
            # is still combinatorially explosive.  RDKit ETKDG below
            # provides deterministic conformer diversity without this
            # hazard.
            ob_injection_ok = False
            if OPENBABEL_AVAILABLE and not deterministic:
                try:
                    _n_ob_r = 3
                    _per_ob_r = max(10, 200 // _n_ob_r)
                    ob_xyz_blocks = []
                    _ob_seen_s: set = set()
                    ob_err: Optional[str] = None
                    for _ri in range(_n_ob_r):
                        _bl, _er = _openbabel_generate_conformer_xyz(
                            smiles, num_confs=_per_ob_r, deterministic=False
                        )
                        if _bl:
                            for _b in _bl:
                                _k = "\n".join(
                                    l.strip() for l in _b.splitlines() if l.strip()
                                )
                                if _k not in _ob_seen_s:
                                    _ob_seen_s.add(_k)
                                    ob_xyz_blocks.append(_b)
                        elif _er and not ob_xyz_blocks:
                            ob_err = _er
                    if ob_xyz_blocks:
                        ob_ids = _inject_openbabel_conformers_into_mol(mol, ob_xyz_blocks)
                        if ob_ids:
                            all_conf_ids.extend(ob_ids)
                            ob_injection_ok = True
                            logger.debug(
                                "OB conformer injection: %d conformers (3 restarts) for %s",
                                len(ob_ids), smiles[:40],
                            )
                        else:
                            logger.debug("OB atom-order mismatch; using RDKit-only pool")
                            mol.RemoveAllConformers()
                    elif ob_err:
                        logger.debug("OB conformer generation: %s", ob_err)
                except Exception as ob_exc:
                    logger.debug("OB conformer generation exception: %s", ob_exc)
                    mol.RemoveAllConformers()

            # --- RDKit ETKDG with 12 diverse fixed seeds (~17 confs/seed) ---
            # Use a thread-pool with a global timeout to prevent hangs on
            # complex metal ring systems where ETKDG can stall.
            #
            # When ``DELFIN_CLASS_AWARE_SEEDS=1`` the seed count is replaced
            # by the per-class value (``_resolve_top_level_seed_count``);
            # default OFF preserves bit-exact pre-patch behaviour by falling
            # through to ``_TOP_LEVEL_SEEDS``.
            _seed_count = _resolve_top_level_seed_count(mol)
            if _seed_count != len(_TOP_LEVEL_SEEDS):
                seeds = list(_PIPELINE_SEEDS[:max(1, _seed_count)])
            else:
                seeds = list(_TOP_LEVEL_SEEDS)
            per_seed = max(1, 200 // len(seeds))
            _etkdg_deadline = _EMBED_TIMEOUT * 3  # total budget for all seeds
            _etkdg_start = __import__('time').monotonic()
            # Determinism: the master switch disables the multi-seed wall-clock
            # budget (and the per-seed join timeout below) so the full seed sweep
            # always runs and the result never depends on timing/CPU load.
            _det_multi = _deterministic_mode()
            try:
                for seed in seeds:
                    if (
                        not _det_multi
                        and __import__('time').monotonic() - _etkdg_start > _etkdg_deadline
                    ):
                        logger.debug("ETKDG multi-seed budget exhausted after %.1fs", _etkdg_deadline)
                        break
                    params_multi = AllChem.ETKDGv3()
                    params_multi.randomSeed = seed
                    params_multi.useRandomCoords = True
                    params_multi.enforceChirality = False
                    try:
                        params_multi.clearConfs = False  # append to OB conformers
                    except Exception:
                        pass
                    # Run with timeout to avoid hanging on single seed
                    _multi_ids = [None]
                    def _do_multi(m=mol, n=per_seed, p=params_multi):
                        try:
                            _multi_ids[0] = list(AllChem.EmbedMultipleConfs(m, numConfs=n, params=p))
                        except Exception:
                            # Kekulize / embed failure on a charged aromatic metal
                            # complex (Zn/Cu/Co-porphyrin, salen: RDKit "Can't
                            # kekulize") -> embed on a DE-AROMATISED copy (same atom
                            # order) and transfer the coordinates back, so the
                            # metallo-macrocycle still builds instead of dropping.
                            try:
                                _dc = _dearomatized_embedding_copy(m)
                                if _dc is not None:
                                    _kids = list(AllChem.EmbedMultipleConfs(_dc, numConfs=n, params=p))
                                    _multi_ids[0] = [
                                        m.AddConformer(Chem.Conformer(_dc.GetConformer(int(_c))), assignId=True)
                                        for _c in _kids
                                    ]
                            except Exception:
                                _multi_ids[0] = None
                    _mt = threading.Thread(target=_do_multi, daemon=True)
                    _mt.start()
                    _mt.join(timeout=None if _det_multi else _EMBED_TIMEOUT)
                    if _mt.is_alive():
                        logger.debug("EmbedMultipleConfs timed out for seed %d", seed)
                        continue
                    if _multi_ids[0]:
                        all_conf_ids.extend(_multi_ids[0])
            except Exception:
                pass

            if all_conf_ids:
                best_conf = None
                best_score = float("inf")
                for cid in all_conf_ids:
                    try:
                        if _has_atom_clash(mol, cid):
                            continue
                        if _has_bad_geometry(mol, cid):
                            continue
                        score = _geometry_quality_score(mol, cid)
                        if score < best_score:
                            best_score = score
                            best_conf = cid
                    except Exception:
                        continue

                if best_conf is None:
                    best_conf = all_conf_ids[0]

                # Keep only the selected conformer so downstream conversion
                # can continue using the default conformer accessor.
                selected = Chem.Conformer(mol.GetConformer(int(best_conf)))
                mol.RemoveAllConformers()
                mol.AddConformer(selected, assignId=True)
                result = 0

        if result != 0:
            if _skip_expensive_pipeline:
                # For large molecules: go straight to permissive embed.
                # ETKDG with standard knowledge can hang on complex metal
                # ring systems (GIL-holding C code, thread timeout ineffective).
                logger.info(
                    "Large molecule (%d heavy atoms): using permissive embed",
                    _n_heavy_atoms,
                )
                result = AllChem.EmbedMolecule(mol, _make_random_embed_params(42))
            else:
                params = AllChem.ETKDGv3()
                params.randomSeed = 42  # For reproducibility

                result = _embed_with_timeout(mol, params)

                if result != 0:
                    # Try with random coordinates if ETKDG fails.
                    # Use an explicit fixed seed: a bare AllChem.ETKDG() keeps
                    # RDKit's wall-clock default (randomSeed = -1) and would
                    # make this fallback path non-deterministic.
                    logger.warning("ETKDG embedding failed, trying random coordinates")
                    _etkdg_rand = AllChem.ETKDG()
                    _etkdg_rand.randomSeed = _deterministic_embed_seed(smiles)
                    result = _embed_with_timeout(mol, _etkdg_rand)

                if result != 0:
                    # Last resort: permissive embed without ETKDG knowledge
                    logger.warning("ETKDG timed out or failed, trying permissive embed")
                    result = AllChem.EmbedMolecule(mol, _make_random_embed_params(42))

        if result != 0:
            # Fallback to legacy embedding logic (used in older dashboards)
            logger.warning("ETKDG embedding failed, trying legacy SMILES embedding fallback")
            legacy_xyz, legacy_err = _smiles_to_xyz_legacy(smiles, has_metal=has_metal)
            if legacy_err is None and legacy_xyz:
                if output_path:
                    Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                    logger.info(f"Converted SMILES to XYZ using legacy fallback: {output_path}")
                return legacy_xyz, None
            # Try multi-strategy approach before giving up
            if has_metal:
                logger.info("Trying multi-strategy fallback after embedding failure")
                multi_xyz, multi_err = _try_multiple_strategies(
                    smiles, output_path, deterministic=deterministic
                )
                if multi_xyz:
                    return multi_xyz, None
                # Last resort: manual coordinate construction for metal complexes
                logger.info("Trying manual metal coordinate construction")
                manual_xyz, manual_err = _manual_metal_embed(smiles)
                if manual_xyz:
                    if output_path:
                        Path(output_path).write_text(manual_xyz, encoding='utf-8')
                        logger.info(f"Converted SMILES to XYZ using manual metal embed: {output_path}")
                    logger.warning(
                        "Used manual coordinate construction - geometry is very rough! "
                        "GOAT/xTB optimization is essential before any calculations."
                    )
                    return manual_xyz, None
            error = f"Failed to generate 3D coordinates for SMILES: {smiles}"
            if legacy_err:
                error = f"{error} (legacy fallback failed: {legacy_err})"
            logger.error(error)
            return None, error

        # Optional UFF refinement for better starting structures
        if apply_uff and not has_metal:
            if os.environ.get("DELFIN_FFFREE_AROM_SEAT", "0") == "1":
                # ROOT seat (metal-free path): constrain aromatic ring bonds to
                # their delocalised CCDC targets during the UFF minimisation so
                # free organic aromatics are BUILT at the right length.  Default
                # OFF → the else-branch below is the unchanged original path.
                try:
                    _uff_seat_aromatic_bonds(mol, max_iters=200)
                    logger.debug("RDKit UFF (aromatic-seat) optimization successful")
                except Exception as e:
                    logger.info(f"RDKit UFF (aromatic-seat) optimization skipped: {e}")
            else:
                try:
                    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                    logger.debug("RDKit UFF optimization successful")
                except Exception as e:
                    logger.info(f"RDKit UFF optimization skipped: {e}")

        # Convert to XYZ format
        xyz_content = _mol_to_xyz(mol)

        # Post-UFF universal sp2 flatten: RDKit UFF leaves residual
        # pyramidalisation at every 3-coordinate sp2 atom (carbonyl-C
        # of DMF/amides, carbamate-C, enamine-N, oxime-N, ring-junction
        # carbons of fused aromatics).  The same helper was already
        # applied in the metal UFF path via _optimize_xyz_openbabel_safe;
        # extending it here keeps every sp2 center planar regardless of
        # which UFF path produced the XYZ — one helper, one rule, no
        # per-system tuning.
        if apply_uff and not has_metal and RDKIT_AVAILABLE:
            try:
                xyz_content = _flatten_sp2_atoms_xyz(xyz_content, mol)
            except Exception as exc:
                logger.debug("Non-metal sp2 flatten skipped: %s", exc)

        # For metal complexes: UFF with universal coordination constraints.
        if apply_uff and has_metal:
            xyz_content = _optimize_xyz_openbabel_safe(
                xyz_content, mol_template=mol
            )
        # For non-metal molecules RDKit UFF above often leaves fused
        # aromatic cores with residual buckle (perylenediimide, carbazole,
        # acridine, …).  Snap each π-system onto its best-fit plane so
        # downstream conformer-pool de-duplication compares planar cores
        # and Avogadro renders them as users expect.
        elif apply_uff and RDKIT_AVAILABLE:
            try:
                xyz_content = _snap_aromatic_rings_in_xyz(
                    xyz_content, mol, rms_threshold=0.05,
                )
            except Exception as exc:
                logger.debug("Non-metal aromatic-plane snap skipped: %s", exc)

        # Check if molecule contains metals and warn about geometry quality
        has_metals = any(atom.GetAtomicNum() in range(21, 31) or  # 3d metals: Sc-Zn
                        atom.GetAtomicNum() in range(39, 49) or  # 4d metals: Y-Cd
                        atom.GetAtomicNum() in range(57, 81)     # Lanthanides + 5d metals
                        for atom in mol.GetAtoms())

        if has_metals:
            logger.warning(
                "SMILES contains metal atoms - generated geometry may be unrealistic! "
                "Coordination geometries from RDKit are rough approximations. "
                "STRONGLY RECOMMENDED: Use GOAT or manual geometry optimization before ORCA calculations."
            )

        if output_path:
            Path(output_path).write_text(xyz_content, encoding='utf-8')
            logger.info(f"Converted SMILES to XYZ using {method}: {output_path}")

        # Pillar 1: runtime topology-invariant gate (env-flag-gated, default OFF)
        return _topology_hard_gate_check(xyz_content, smiles)

    except Exception as e:
        msg = str(e)
        if "kekulize" in msg or "Explicit valence" in msg:
            legacy_xyz, legacy_err = _smiles_to_xyz_unsanitized_fallback(smiles)
            if legacy_err is None and legacy_xyz:
                if output_path:
                    Path(output_path).write_text(legacy_xyz, encoding='utf-8')
                    logger.info(f"Converted SMILES to XYZ using unsanitized fallback: {output_path}")
                return _topology_hard_gate_check(legacy_xyz, smiles)
        # Last resort: multi-strategy fallback for metal complexes
        if has_metal:
            logger.info("Trying multi-strategy fallback after exception: %s", e)
            multi_xyz, multi_err = _try_multiple_strategies(
                smiles, output_path, deterministic=deterministic
            )
            if multi_xyz:
                return _topology_hard_gate_check(multi_xyz, smiles)
        error = f"Error converting SMILES to XYZ: {e}"
        logger.error(error, exc_info=True)
        return None, error


def _smiles_to_xyz_legacy(smiles: str, has_metal: bool) -> Tuple[Optional[str], Optional[str]]:
    """Legacy embedding fallback (matches older dashboard behavior)."""
    if not RDKIT_AVAILABLE:
        return None, "RDKit is not installed"

    stk_error = None
    mol = None

    # Try stk first for metal complexes
    if has_metal and STK_AVAILABLE:
        try:
            bb = stk.BuildingBlock(smiles)
            mol = bb.to_rdkit_mol()
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol, randomSeed=42, useRandomCoords=True)
        except Exception as e:
            stk_error = str(e)
            mol = None

    # RDKit fallback
    if mol is None:
        mol, rdkit_note = mol_from_smiles_rdkit(smiles, allow_metal=has_metal)
        if mol is None:
            if stk_error:
                return None, f"RDKit: {rdkit_note}; stk: {stk_error}"
            return None, rdkit_note

        try:
            mol = Chem.AddHs(mol)
        except Exception:
            pass

        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useRandomCoords = True
            try:
                params.maxAttempts = 100
            except Exception:
                pass

            try:
                result = AllChem.EmbedMolecule(mol, params, maxAttempts=50)
            except TypeError:
                result = AllChem.EmbedMolecule(mol, params)

            if result == -1:
                params2 = AllChem.ETKDGv3()
                params2.randomSeed = 42
                params2.useRandomCoords = True
                try:
                    params2.maxAttempts = 200
                except Exception:
                    pass
                try:
                    result = AllChem.EmbedMolecule(mol, params2, maxAttempts=200)
                except TypeError:
                    result = AllChem.EmbedMolecule(mol, params2)
                if result == -1:
                    return None, "Legacy embed failed to generate 3D structure"
        except Exception as e:
            return None, f"Legacy RDKit error: {e}"

    try:
        xyz_content = _mol_to_xyz(mol)
        return xyz_content, None
    except Exception as e:
        return None, f"Legacy coordinate error: {e}"
