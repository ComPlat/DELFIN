"""Linkage isomers and alternative binding-site exploration of the MANTA constructor.

Moved verbatim from delfin/smiles_converter.py (MANTA split, 2026-10); every
statement is the original text, only the imports below are new.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

from delfin.common.logging import get_logger
from delfin.manta.conformer_pools import (
    _rank_template_conformers,
)
from delfin.manta.converter_flags import (
    _deterministic_mode,
    _multi_sigma_v2_active,
    _multi_sigma_v2_budget,
)
from delfin.manta.hapto_detect import (
    _classify_complex_class,
)
from delfin.manta.isomer_labels import (
    _is_viable_donor,
    _ligand_fragments,
)
from delfin.manta.ligand_placement import (
    _snap_aromatic_rings_in_xyz,
)
from delfin.manta.ml_tables import (
    Chem,
    _METAL_SET,
)
from delfin.manta.openbabel_optimize import (
    _optimize_xyz_openbabel_safe,
)
from delfin.manta.topo_isomers import (
    _build_topology_xyz,
)
from delfin.manta.topology_checks import (
    _fragment_topology_ok,
)

logger = get_logger("delfin.smiles_converter")


# ---------------------------------------------------------------------------
# Linkage isomers (Feature 2)
# ---------------------------------------------------------------------------

def _find_linkage_alternatives(mol) -> List[Tuple[int, int, int, str]]:
    """Detect ambidentate ligands and return alternative coordination modes.

    Supported patterns:
      - NO2⁻ (N-donor → O-donor, label 'nitrito')
      - SCN⁻ (S-donor → N-donor, label 'isothiocyanato-N')
      - CN⁻  (C-donor → N-donor, label 'isocyano')

    Returns list of (metal_idx, current_donor_idx, alt_donor_idx, label).
    """
    alternatives: List[Tuple[int, int, int, str]] = []
    metal_neighbor_sets: Dict[int, set] = {}

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        metal_idx = atom.GetIdx()
        bonded = {nbr.GetIdx() for nbr in atom.GetNeighbors()}
        metal_neighbor_sets[metal_idx] = bonded

        for donor in atom.GetNeighbors():
            donor_idx = donor.GetIdx()
            donor_sym = donor.GetSymbol()

            # --- NO2 → nitrito (N-donor to O-donor) ---
            if donor_sym == 'N':
                o_nbrs = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O' and n.GetIdx() not in bonded
                ]
                if len(o_nbrs) >= 2:
                    alt_o = o_nbrs[0]
                    alternatives.append((metal_idx, donor_idx, alt_o.GetIdx(), 'nitrito'))

            # --- SCN → isothiocyanato-N (S-donor to N-donor) ---
            if donor_sym == 'S':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for n_nbr in c_nbr.GetNeighbors():
                            if (n_nbr.GetSymbol() == 'N'
                                    and n_nbr.GetIdx() != donor_idx
                                    and n_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, n_nbr.GetIdx(), 'isothiocyanato-N')
                                )

            # --- CN → isocyano (C-donor to N-donor via triple bond) ---
            if donor_sym == 'C':
                for n_nbr in donor.GetNeighbors():
                    if (n_nbr.GetSymbol() == 'N'
                            and n_nbr.GetIdx() not in bonded):
                        bond = mol.GetBondBetweenAtoms(donor_idx, n_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 2.5:
                            alternatives.append(
                                (metal_idx, donor_idx, n_nbr.GetIdx(), 'isocyano')
                            )

            # --- NCS (N-donor → S-donor via C): thiocyanato-S ---
            if donor_sym == 'N':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for s_nbr in c_nbr.GetNeighbors():
                            if (s_nbr.GetSymbol() == 'S'
                                    and s_nbr.GetIdx() != donor_idx
                                    and s_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, s_nbr.GetIdx(), 'thiocyanato-S')
                                )

            # --- Carboxylate (O1 → O2): alternative carboxylate oxygen ---
            if donor_sym == 'O':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C':
                        o_nbrs = [
                            n for n in c_nbr.GetNeighbors()
                            if n.GetSymbol() == 'O'
                            and n.GetIdx() != donor_idx
                            and n.GetIdx() not in bonded
                        ]
                        if o_nbrs:
                            alternatives.append(
                                (metal_idx, donor_idx, o_nbrs[0].GetIdx(), 'carboxylato-alt')
                            )

            # --- Sulfoxide (S-donor → O-donor via double bond) ---
            if donor_sym == 'S':
                o_nbrs = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O'
                    and n.GetIdx() not in bonded
                    and mol.GetBondBetweenAtoms(donor_idx, n.GetIdx()) is not None
                    and mol.GetBondBetweenAtoms(donor_idx, n.GetIdx()).GetBondTypeAsDouble() >= 1.5
                ]
                if o_nbrs:
                    alternatives.append(
                        (metal_idx, donor_idx, o_nbrs[0].GetIdx(), 'sulfoxide-O')
                    )

            # --- Sulfite / Sulfonate (S-donor → O-donor, S with ≥3 O) ---
            if donor_sym == 'S':
                o_nbrs_all = [
                    n for n in donor.GetNeighbors()
                    if n.GetSymbol() == 'O' and n.GetIdx() not in bonded
                ]
                if len(o_nbrs_all) >= 2:
                    for alt_o in o_nbrs_all:
                        alternatives.append(
                            (metal_idx, donor_idx, alt_o.GetIdx(), 'sulfito-O')
                        )

            # --- Nitrosyl (N=O, N-donor → O-donor) ---
            if donor_sym == 'N':
                for o_nbr in donor.GetNeighbors():
                    if (o_nbr.GetSymbol() == 'O'
                            and o_nbr.GetIdx() not in bonded
                            and len(list(o_nbr.GetNeighbors())) == 1):
                        bond = mol.GetBondBetweenAtoms(donor_idx, o_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                            alternatives.append(
                                (metal_idx, donor_idx, o_nbr.GetIdx(), 'isonitrosyl-O')
                            )

            # --- Nitrosyl O-bound → N-bound ---
            if donor_sym == 'O':
                for n_nbr in donor.GetNeighbors():
                    if (n_nbr.GetSymbol() == 'N'
                            and n_nbr.GetIdx() not in bonded
                            and len(list(donor.GetNeighbors())) == 1):
                        bond = mol.GetBondBetweenAtoms(donor_idx, n_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                            alternatives.append(
                                (metal_idx, donor_idx, n_nbr.GetIdx(), 'nitrosyl-N')
                            )

            # --- Selenocyanate SeCN (Se-donor → N-donor via C) ---
            if donor_sym == 'Se':
                for c_nbr in donor.GetNeighbors():
                    if c_nbr.GetSymbol() == 'C' and c_nbr.GetIdx() not in bonded:
                        for n_nbr in c_nbr.GetNeighbors():
                            if (n_nbr.GetSymbol() == 'N'
                                    and n_nbr.GetIdx() != donor_idx
                                    and n_nbr.GetIdx() not in bonded):
                                alternatives.append(
                                    (metal_idx, donor_idx, n_nbr.GetIdx(), 'isoselenocyanato-N')
                                )

            # --- Pyrazole / Imidazole / Triazole (N1 → N2 in same ring) ---
            if donor_sym == 'N':
                try:
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if donor_idx in ring:
                            other_n = [
                                idx for idx in ring
                                if mol.GetAtomWithIdx(idx).GetSymbol() == 'N'
                                and idx != donor_idx
                                and idx not in bonded
                            ]
                            for alt_n in other_n:
                                alternatives.append(
                                    (metal_idx, donor_idx, alt_n, 'N-isomer')
                                )
                except Exception:
                    pass

    return alternatives


def _rewire_linkage(mol, metal_idx: int, old_donor: int, new_donor: int) -> Optional[object]:
    """Return a new mol with the M–old_donor bond replaced by M–new_donor.

    Returns None if the rewiring would duplicate an existing bond or fails.
    """
    try:
        rw = Chem.RWMol(mol)
        # Abort if new_donor is already bonded to metal
        if rw.GetBondBetweenAtoms(metal_idx, new_donor) is not None:
            return None
        rw.RemoveBond(metal_idx, old_donor)
        rw.AddBond(metal_idx, new_donor, Chem.BondType.SINGLE)
        rw.UpdatePropertyCache(strict=False)
        return rw.GetMol()
    except Exception as e:
        logger.debug("_rewire_linkage failed: %s", e)
        return None


def _generate_linkage_isomers(
    mol,
    smiles: str,
    apply_uff: bool = True,
    max_template_tries: int = 8,
) -> List[Tuple[str, str]]:
    """Build linkage isomers by rewiring ambidentate ligands and generating XYZ.

    For each alternative coordination mode, rewires the metal bond and builds
    an idealized topology structure (OB UFF optimized).  The builder is tried
    against the top ``max_template_tries`` template conformers (ranked by
    :func:`_rank_template_conformers`) so one bad conformer does not silently
    discard valid linkage isomers.

    Returns [(xyz_string, label), …] where label is e.g. 'nitrito'.
    """
    results: List[Tuple[str, str]] = []
    alternatives = _find_linkage_alternatives(mol)

    for metal_idx, old_donor, new_donor, type_label in alternatives:
        alt_mol = _rewire_linkage(mol, metal_idx, old_donor, new_donor)
        if alt_mol is None:
            continue
        try:
            # Find metal in alt_mol (same index)
            metal_atom = alt_mol.GetAtomWithIdx(metal_idx)
            donor_indices = [nbr.GetIdx() for nbr in metal_atom.GetNeighbors()]
            n_coord = len(donor_indices)
            geom_map = {2: 'LIN', 3: 'TP', 4: 'SQ', 5: 'TBP', 6: 'OH', 7: 'PBP'}
            if n_coord not in geom_map:
                continue
            geom = geom_map[n_coord]
            perm = list(range(n_coord))  # identity: donor i → position i

            candidate_cids = _rank_template_conformers(
                alt_mol, top_k=max_template_tries
            )
            if not candidate_cids:
                candidate_cids = [None]

            for cid in candidate_cids:
                # Build WITHOUT UFF first (fast), check topology, THEN UFF.
                xyz = _build_topology_xyz(
                    alt_mol, metal_idx, donor_indices, perm, geom, False,
                    conf_id=cid,
                )
                if xyz is None:
                    continue
                if not _fragment_topology_ok(xyz, smiles):
                    logger.debug(
                        "linkage %s via template cid=%s: "
                        "fragment topology mismatch, trying next template",
                        type_label, cid,
                    )
                    continue
                if apply_uff:
                    xyz = _optimize_xyz_openbabel_safe(
                        xyz, mol_template=alt_mol
                    )
                    xyz = _snap_aromatic_rings_in_xyz(xyz, alt_mol)
                results.append((xyz, type_label))
                break
        except Exception as e:
            logger.debug("Linkage isomer (%s) failed: %s", type_label, e)

    return results


def _generate_alternative_binding_modes(
    mol,
    smiles: str,
    apply_uff: bool = True,
    max_alternatives: int = 20,
    max_template_tries: int = 8,
) -> List[Tuple[str, str]]:
    """Generate isomers by swapping donor atoms with alternative binding sites.

    Only considers viable donors (atoms with available lone pairs) on the
    SAME ligand fragment as the current donor.  This avoids generating
    nonsensical structures from swapping donors across unrelated ligands
    or using C=O carbonyl oxygens as coordination donors.

    For each rewired candidate the builder is tried against the top
    ``max_template_tries`` template conformers (ranked by geometry quality)
    and the first result that passes ``_fragment_topology_ok`` is kept.  This
    decouples constitutional-isomer generation from a single, possibly poor,
    template conformer (see Codex findings for the reasoning).

    Returns [(xyz_string, label), ...].
    """
    results: List[Tuple[str, str]] = []

    # Iter-8.6 multi-σ timeout-mitigation: env-gated wall-clock budget.
    # Multi-metal fragments fire 8-seed × 6 s ETKDG sweeps per rewire
    # candidate; with ~12 alt donors × 2 metals × top_k=8 templates the
    # cumulative cost exceeds the 600 s per-SMILES driver budget, so the
    # outer pool_evaluator marks the SMILES as timeout and the entire
    # multi-σ result set is lost.  When DELFIN_ALT_MODE_BUDGET_S is set
    # to a positive integer, the candidate loop breaks once that many
    # wall-clock seconds have elapsed.  Partial results are returned.
    # Default 0 → disabled, bit-exact baseline behaviour.
    import time as _time_mod
    # Wave-4 (A+D forensics): Phase 4D's universal 60s cap killed large
    # bimetallic multi-σ conversions (Sn-Ir, Sn-Rh, Os-Sn, Fe-Fe μ-O,
    # In-Ru). class-conditional: 60s for sigma/hapto (safe), 0s/unlimited
    # for multi_sigma/multi_hapto (needs 100-500s).  Env override still
    # respected for testing.
    try:
        _env_budget_raw = os.environ.get("DELFIN_ALT_MODE_BUDGET_S")
        if _env_budget_raw is not None and _env_budget_raw.strip() != "":
            _alt_budget_s = float(_env_budget_raw)
        else:
            try:
                _cls = _classify_complex_class(mol)
            except Exception:
                _cls = "sigma"
            # Multi-sigma V2: opt out of the historical "unlimited" budget
            # for multi-metal sigma — that branch was the second biggest
            # contributor to the 600 s pool-evaluator timeouts (forensics
            # 2026-05-13).  Cap to the heavy-atom-scaled wall-clock so
            # alt-binding-mode exploration cannot starve the rest of the
            # pipeline.  Other classes keep their pre-patch budget.
            _v2_alt_cap: Optional[float] = None
            try:
                if _cls == "multi_sigma" and _multi_sigma_v2_active(mol):
                    _heavy_n_alt = sum(
                        1 for _a in mol.GetAtoms() if _a.GetAtomicNum() > 1
                    )
                    # Re-use the multi-metal augmentation budget — these
                    # two stages have similar per-call cost characteristics.
                    _v2_alt_cap = float(
                        _multi_sigma_v2_budget(_heavy_n_alt)["mm_walltime"]
                    )
            except Exception:
                _v2_alt_cap = None
            if _v2_alt_cap is not None:
                _alt_budget_s = _v2_alt_cap
            else:
                _alt_budget_s = (
                    0.0 if _cls in ("multi_sigma", "multi_hapto") else 60.0
                )
    except Exception:
        _alt_budget_s = 60.0
    # DETERMINISM (2026-07-28): this was the ONLY isomer-generating path whose wall-clock budget was
    # not gated on _deterministic_mode() -- every sibling already is (2-metal :27381, N-metal :27605,
    # MM-augmentation, ETKDG multi-seed, embed joins).  Consequence: alternative binding modes are
    # DROPPED under load, so the ISOMER SET itself became load-dependent.  Measured on LATTUW: 14
    # frames under load vs 16 on a free box -- the two missing ones are alt-bind-C isomers.  It passes
    # within-run determinism (both replica builds run under the same load) and only differs ACROSS
    # runs, which is why the byte-determinism gate never caught it.  That violates both "completeness is
    # sacred" and the deterministic-manifold claim.  0.0 disables the wall clock; termination falls back
    # to the deterministic caps, exactly like the siblings.
    if _deterministic_mode():
        _alt_budget_s = 0.0
    _alt_t0 = _time_mod.monotonic() if _alt_budget_s > 0 else None

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in _METAL_SET:
            continue
        if _alt_t0 is not None and (_time_mod.monotonic() - _alt_t0) > _alt_budget_s:
            logger.debug(
                "alt-binding budget %.1fs exceeded — stopping at metal loop",
                _alt_budget_s,
            )
            return results
        metal_idx = atom.GetIdx()
        bonded = {nbr.GetIdx() for nbr in atom.GetNeighbors()}
        current_donors = list(bonded)

        # Build ligand fragments to ensure we only swap within the same ligand
        fragments = _ligand_fragments(mol, metal_idx)
        atom_to_frag: Dict[int, int] = {}
        for fi, frag in enumerate(fragments):
            for aidx in frag:
                atom_to_frag[aidx] = fi

        # For each current donor, find viable alternative donors on the same fragment
        for current_d in current_donors:
            current_frag = atom_to_frag.get(current_d, -1)
            if current_frag < 0:
                continue
            current_sym = mol.GetAtomWithIdx(current_d).GetSymbol()

            for alt_d in fragments[current_frag]:
                if alt_d == current_d or alt_d in bonded:
                    continue
                if not _is_viable_donor(mol, alt_d, bonded):
                    continue
                if len(results) >= max_alternatives:
                    return results
                if _alt_t0 is not None and (_time_mod.monotonic() - _alt_t0) > _alt_budget_s:
                    logger.debug(
                        "alt-binding budget %.1fs exceeded — stopping at donor loop",
                        _alt_budget_s,
                    )
                    return results

                alt_sym = mol.GetAtomWithIdx(alt_d).GetSymbol()
                alt_mol = _rewire_linkage(mol, metal_idx, current_d, alt_d)
                if alt_mol is None:
                    continue
                try:
                    new_metal = alt_mol.GetAtomWithIdx(metal_idx)
                    donor_indices = [nbr.GetIdx() for nbr in new_metal.GetNeighbors()]
                    n_coord = len(donor_indices)
                    geom_map = {2: 'LIN', 3: 'TP', 4: 'SQ', 5: 'TBP', 6: 'OH', 7: 'PBP'}
                    if n_coord not in geom_map:
                        continue
                    geom = geom_map[n_coord]
                    perm = list(range(n_coord))

                    if current_sym == alt_sym:
                        label = f'alt-{alt_sym}-isomer'
                    else:
                        label = f'alt-bind-{alt_sym}'

                    # Multi-template trial: iterate over the best-scored
                    # conformer IDs and keep the first XYZ that validates.
                    candidate_cids = _rank_template_conformers(
                        alt_mol, top_k=max_template_tries
                    )
                    if not candidate_cids:
                        # No usable template conformer — try de-novo once.
                        candidate_cids = [None]

                    accepted = False
                    for cid in candidate_cids:
                        # Multi-sigma V2: per-template budget check.  The
                        # outer donor-loop check (above) can still overshoot
                        # by 30-60 s on multi-metal complexes where each
                        # template alignment + UFF + topology gate runs
                        # serially.  This inner check bounds the overshoot.
                        if (
                            _alt_t0 is not None
                            and (_time_mod.monotonic() - _alt_t0) > _alt_budget_s
                        ):
                            logger.debug(
                                "alt-binding budget %.1fs exceeded — "
                                "stopping at template loop",
                                _alt_budget_s,
                            )
                            return results
                        # Build WITHOUT UFF first, check topology, THEN UFF.
                        xyz = _build_topology_xyz(
                            alt_mol, metal_idx, donor_indices, perm,
                            geom, False, conf_id=cid,
                        )
                        if xyz is None:
                            continue
                        if not _fragment_topology_ok(xyz, smiles):
                            logger.debug(
                                "alt-mode %s via template cid=%s: "
                                "fragment topology mismatch, trying next template",
                                label, cid,
                            )
                            continue
                        if apply_uff:
                            xyz = _optimize_xyz_openbabel_safe(
                                xyz, mol_template=alt_mol
                            )
                        results.append((xyz, label))
                        accepted = True
                        break

                    if not accepted:
                        logger.debug(
                            "alt-mode %s: no template conformer produced a "
                            "topology-consistent XYZ (tried %d)",
                            label, len(candidate_cids),
                        )
                except Exception as e:
                    logger.debug("Alternative binding mode failed: %s", e)
                    continue

    return results
