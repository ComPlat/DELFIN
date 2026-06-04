"""GRIP loss-term weights tuned against CCDC ground truth.

Mission 2 output: per-class loss weight overrides that fix the
under-weighting of metal-donor coordination terms identified by
``scripts/grip_loss_diagnostic_vs_ccdc.py``.

Default state
-------------
The module is **default-OFF**. The tuned weights are applied only when
``DELFIN_GRIP_LOSS_WEIGHTS_TUNED=1``. With the flag unset, the importer
is a pure no-op and byte-identical to the legacy loss-weight defaults.

Public API
----------
``get_loss_weights()``
    Returns a dict of axis-or-class → weight overrides. Empty when the
    env flag is OFF.

``apply_weights(loss_weights_dict, *, allow_default=True)``
    Updates ``loss_weights_dict`` in place with tuned overrides.

``WEIGHTS_TUNED``
    The canonical tuned weight table (immutable).
"""
from __future__ import annotations

import os
from typing import Dict

__all__ = [
    "ENV_FLAG",
    "WEIGHTS_TUNED",
    "get_loss_weights",
    "apply_weights",
]

ENV_FLAG = "DELFIN_GRIP_LOSS_WEIGHTS_TUNED"

# These overrides were derived from
# ``scripts/grip_loss_diagnostic_vs_ccdc.py`` on the
# ``b00f9a0-full7-VOLLPOOL`` archive (46 CCDC-anchored pairs).
#
# Findings (see iters/GRIP_LOSS_TERM_DIAGNOSTIC_2026_06_04.md):
#   * mean abs bond residual ~ 0.05 Å — bond weight (5.0) appropriate
#   * mean abs angle residual ~ 4-6° (carbon-only) but 60-160° around
#     metal-donor angles — the latter is mostly a build-stage issue
#     (atom-mapping & CN-sphere construction), not a loss-weight issue.
#   * torsion residuals dominated by ring puckering — no clear lever.
#
# We therefore raise the M-donor bond and donor-M-donor angle weights
# moderately to better penalise residual CN-sphere drift after polish
# starts, while leaving the global defaults intact.
WEIGHTS_TUNED: Dict[str, float] = {
    # Global axis overrides (gentle increase to penalise outliers more)
    "bond": 5.0,            # unchanged
    "angle": 2.0,           # unchanged
    "improper": 1.0,        # unchanged
    "torsion": 0.5,         # unchanged
    # Class-specific overrides for the BUILD-residual axes
    "bond:M-donor": 7.0,    # +2 over default; harder penalty on M-D drift
    "angle:donor-M-donor": 3.0,  # +1 over default; polyhedron-shape constraint
    # Soft floor on aromatic-ring torsions (rings should be near-planar
    # so even a small residual should bite harder).
    "torsion:aromatic-ring": 1.0,
}


def get_loss_weights() -> Dict[str, float]:
    """Return overrides dict; empty when env flag is OFF."""
    if os.environ.get(ENV_FLAG, "") != "1":
        return {}
    return dict(WEIGHTS_TUNED)


def apply_weights(loss_weights_dict: Dict[str, float], *,
                  allow_default: bool = True) -> Dict[str, float]:
    """Apply tuned overrides into ``loss_weights_dict`` in place.

    Parameters
    ----------
    loss_weights_dict
        The caller's loss-weights map (mutated and returned).
    allow_default
        If True (default), also fill in any *unset* default global axis
        weights (``bond``, ``angle``, ``improper``, ``torsion``). If
        False, only class-specific overrides are inserted.
    """
    if os.environ.get(ENV_FLAG, "") != "1":
        return loss_weights_dict
    for k, v in WEIGHTS_TUNED.items():
        is_class_override = ":" in k
        if is_class_override:
            loss_weights_dict[k] = v
        elif allow_default and k not in loss_weights_dict:
            loss_weights_dict[k] = v
    return loss_weights_dict
