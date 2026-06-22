"""Class-aware ensemble router — Phase 2 of Hybrid-Path.

Routes SMILES to the right class-specialist based on empirically validated
champion-routing table from quality_framework/results/champion_routing.md.

Default behaviour (env DELFIN_ENSEMBLE_ROUTER unset/0): no-op, returns None
to signal caller should fall back to legacy smiles_converter.smiles_to_xyz.
This makes the router fully opt-in for safe incremental rollout.

Per nature_project/15_HYBRID_PATH_FINAL.md.
"""
from __future__ import annotations

import logging
import os
from typing import Optional, Tuple

from delfin.class_modules import Specialist, get_specialist
from delfin.classify import (
    BLOCK_3D, BLOCK_4D, BLOCK_5D, BLOCK_LN, BLOCK_AN, BLOCK_S, BLOCK_P,
    COORD_SIGMA, COORD_HAPTO, COORD_MULTI_SIGMA, COORD_MULTI_HAPTO,
    COORD_NO_METAL, COORD_BORANE, COORD_CLUSTER,
    ClassFeatures, classify_smiles,
)

logger = logging.getLogger(__name__)


# Empirical routing-table from quality_framework/results/champion_routing.md
# Format: (coord_class, metal_block) → specialist_id
# Specialist-ID maps to a class_modules/<name>_specialist.py via the registry.
ROUTING_TABLE: dict[tuple[str, str], str] = {
    # multi-sigma — 123a130 dominates all metal blocks
    (COORD_MULTI_SIGMA, BLOCK_3D): "multi_sigma_123a130",
    (COORD_MULTI_SIGMA, BLOCK_4D): "multi_sigma_123a130",
    (COORD_MULTI_SIGMA, BLOCK_5D): "multi_sigma_123a130",
    (COORD_MULTI_SIGMA, BLOCK_P):  "multi_sigma_123a130",
    # sigma — re-routed 2026-05-12 after smoke500 verdict:
    # 1e7eefe-prev-wt-control had only ~5 files in sample → small-n misranking.
    # noHfix archive has full pool + 1.50° SAFLOA cis_rmsd (champion).
    (COORD_SIGMA, BLOCK_3D): "sigma_1e7eefe_noHfix",
    (COORD_SIGMA, BLOCK_4D): "sigma_229e5dc",
    (COORD_SIGMA, BLOCK_5D): "sigma_1e7eefe_noHfix",
    (COORD_SIGMA, BLOCK_LN): "sigma_1e7eefe_noHfix",
    # hapto — re-routed same reason: prefer noHfix for 3d/4d
    (COORD_HAPTO, BLOCK_3D): "hapto_1e7eefe_noHfix",
    (COORD_HAPTO, BLOCK_4D): "hapto_1e7eefe_noHfix",
    (COORD_HAPTO, BLOCK_5D): "hapto_3eb9aaa",
    (COORD_HAPTO, BLOCK_LN): "hapto_1e7eefe_noHfix",
    (COORD_HAPTO, BLOCK_P):  "hapto_5b3e0d2",
    # multi-hapto — noHfix 3d, c0b0451 5d, 29c2398 p
    (COORD_MULTI_HAPTO, BLOCK_3D): "multi_hapto_1e7eefe_noHfix",
    (COORD_MULTI_HAPTO, BLOCK_5D): "multi_hapto_c0b0451",
    (COORD_MULTI_HAPTO, BLOCK_P):  "multi_hapto_29c2398",
    # no-metal organic — fallback path (will be added in Phase 2)
    (COORD_NO_METAL, "no-metal"): "organic_fallback",
    # cluster / borane — Phase 3+ targets
    (COORD_CLUSTER, BLOCK_3D): "cluster_specialist",  # TBD
    (COORD_BORANE, "no-metal"): "borane_specialist",  # TBD
}


def lookup_specialist_id(features: ClassFeatures) -> Optional[str]:
    """Return specialist-ID for a given ClassFeatures, or None if no route."""
    return ROUTING_TABLE.get((features.coord_class, features.metal_block))


def route(smiles: str) -> Optional[Specialist]:
    """Classify SMILES and return the routed specialist, or None.

    Returns None if:
      - No routing rule matches this (coord_class, metal_block) pair
      - Routing rule exists but specialist not yet registered

    Pure: no env-flag check (the env-flag gate is in route_or_fallback).
    """
    features = classify_smiles(smiles)
    spec_id = lookup_specialist_id(features)
    if spec_id is None:
        return None
    return get_specialist(spec_id)


def route_or_fallback(
    smiles: str,
    fallback_fn,
    **kwargs,
) -> Tuple[Optional[str], Optional[str]]:
    """Env-flag-gated routing entry-point.

    Env DELFIN_ENSEMBLE_ROUTER (default 0):
      0 — router disabled, always delegates to fallback_fn (legacy)
      1 — router active, falls back to fallback_fn when no specialist
      2 — router strict, returns error when no specialist (no fallback)

    Args:
        smiles: input SMILES
        fallback_fn: callable(smiles, **kwargs) -> Tuple[xyz, error]
                     used when router disabled or no specialist matches
        **kwargs: passed to specialist.convert or fallback_fn

    Returns: (xyz, error) tuple matching smiles_to_xyz API.
    """
    try:
        mode = int(os.environ.get("DELFIN_ENSEMBLE_ROUTER", "0"))
    except ValueError:
        mode = 0

    if mode == 0:
        return fallback_fn(smiles, **kwargs)

    # Lazy-register Champion specialists (only when router actually used)
    try:
        from delfin.class_modules.registry import register_all_specialists
        register_all_specialists()
    except Exception as exc:
        logger.debug("specialist registry init failed: %s", exc)

    specialist = route(smiles)
    if specialist is None:
        if mode >= 2:
            features = classify_smiles(smiles)
            return None, (
                f"ensemble_router: no specialist for "
                f"({features.coord_class}, {features.metal_block})"
            )
        logger.debug("ensemble_router: no specialist, falling back")
        return fallback_fn(smiles, **kwargs)

    try:
        return specialist.convert(smiles, **kwargs)
    except Exception as exc:
        if mode >= 2:
            return None, f"ensemble_router: {specialist.id} failed: {exc}"
        logger.warning(
            "ensemble_router: %s failed (%s), falling back to legacy",
            specialist.id, exc,
        )
        return fallback_fn(smiles, **kwargs)


def routing_summary() -> dict:
    """Return diagnostic info about routing table + registered specialists.

    Useful for debugging which (class, block) cells have specialists.
    """
    from delfin.class_modules import list_registered
    registered = set(list_registered())
    summary = {
        "n_routes_defined": len(ROUTING_TABLE),
        "n_specialists_registered": len(registered),
        "routes": [],
    }
    for (coord, block), spec_id in sorted(ROUTING_TABLE.items()):
        summary["routes"].append({
            "coord_class": coord,
            "metal_block": block,
            "specialist_id": spec_id,
            "registered": spec_id in registered,
        })
    return summary
