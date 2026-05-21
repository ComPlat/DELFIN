"""Specialist registry — wires empirical champion commits to specialist IDs.

Each ID matches an entry in delfin/ensemble_router.py ROUTING_TABLE.
Champion selection backed by quality_framework/results/champion_routing.md.

This module is imported by delfin.ensemble_router when DELFIN_ENSEMBLE_ROUTER>0
to populate the registry.  Pure side-effect: register specialists.
"""
from __future__ import annotations

from delfin.class_modules import register
from delfin.class_modules.champion_specialist import ChampionSpecialist


# Per quality_framework/results/champion_routing.md (Phase 1 score-matrix)
_SPECIALIST_SPECS = [
    # (specialist_id,                  source_commit,        archive_subdir)
    ("multi_sigma_123a130",            "123a130",            "123a130"),
    ("sigma_1e7eefe_prev",             "1e7eefe",            "1e7eefe-prev-wt-control"),
    ("sigma_229e5dc",                  "229e5dc",            "229e5dc"),
    ("sigma_1e7eefe_noHfix",           "1e7eefe",            "1e7eefe-noHfix"),
    ("hapto_1e7eefe_prev",             "1e7eefe",            "1e7eefe-prev-wt-control"),
    ("hapto_3eb9aaa",                  "3eb9aaa",            "3eb9aaa"),
    ("hapto_1e7eefe_noHfix",           "1e7eefe",            "1e7eefe-noHfix"),
    ("hapto_5b3e0d2",                  "5b3e0d2",            "5b3e0d2_smoke500"),
    ("multi_hapto_1e7eefe_noHfix",     "1e7eefe",            "1e7eefe-noHfix"),
    ("multi_hapto_c0b0451",            "c0b0451",            "HEAD-c0b0451"),
    ("multi_hapto_29c2398",            "29c2398",            "29c2398"),
]


_INITIALIZED = False


def register_all_specialists() -> int:
    """Register every Champion specialist.  Idempotent.

    Returns: number of newly-registered specialists.
    """
    global _INITIALIZED
    if _INITIALIZED:
        return 0
    n_registered = 0
    for spec_id, commit, archive_sub in _SPECIALIST_SPECS:
        spec = ChampionSpecialist(
            id=spec_id,
            source_commit=commit,
            archive_subdir=archive_sub,
        )
        try:
            register(spec)
            n_registered += 1
        except ValueError:
            # Already registered
            pass
    _INITIALIZED = True
    return n_registered


def specialists_available() -> dict:
    """Return per-specialist availability diagnostic."""
    register_all_specialists()
    from delfin.class_modules import _REGISTRY
    return {
        spec_id: spec.is_available()
        for spec_id, spec in _REGISTRY.items()
        if isinstance(spec, ChampionSpecialist)
    }
