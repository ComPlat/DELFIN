"""Declarative capability-wiring rules for the Pipeline.

Historically the pipeline auto-injected one artifact: an upstream ORCA ``.gbw``
became the ``moread`` kwarg of the next ``orca_*`` step (an inline
``if name.startswith("orca_")`` in :meth:`Pipeline.run`).  This module lifts
that single rule into a small data table so new artifact hand-offs (hessian,
molden, …) are a one-line :class:`WiringRule` plus a ``consumes`` tag on the
consuming adapter — no edits to the executor.

Backwards compatibility: :data:`DEFAULT_RULES` contains exactly the historical
gbw→moread rule, and :func:`autowire_kwargs` reproduces the old behaviour
byte-for-byte when the consumer declares no ``consumes`` tags.  Capability-based
firing (via ``consumes``) is strictly additive on top.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Tuple


@dataclass(frozen=True)
class WiringRule:
    """Map an upstream-produced capability onto a kwarg of a downstream step.

    A rule fires for a consumer step when the artifact is available upstream,
    the consumer has not set the kwarg explicitly, and *either* the name
    predicate matches *or* the consumer declares it ``consumes`` the capability.
    """

    capability: str                       # upstream capability tag, e.g. "gbw"
    inject_kwarg: str                     # kwarg to set on the consumer, e.g. "moread"
    applies_to: Callable[[str], bool]     # predicate on the consumer step name
    artifact_key: str = ""                # key in current_artifacts (defaults to capability)


def _is_orca(name: str) -> bool:
    return name.startswith("orca_")


# The historical default: upstream GBW → next ORCA step's MOREAD.  Keep this
# first and unchanged so existing pipelines behave identically.
DEFAULT_RULES: Tuple[WiringRule, ...] = (
    WiringRule(
        capability="gbw",
        inject_kwarg="moread",
        applies_to=_is_orca,
        artifact_key="gbw",
    ),
    # Future opt-in rules (inactive until a consumer declares `consumes`):
    # WiringRule("hessian", "hess_file", lambda n: n == "imag_fix", "hess"),
)


def autowire_kwargs(
    step_name: str,
    merged_kwargs: Dict[str, Any],
    current_artifacts: Dict[str, Any],
    *,
    consumes: frozenset = frozenset(),
    rules: Tuple[WiringRule, ...] = DEFAULT_RULES,
) -> Dict[str, Any]:
    """Inject artifact-backed kwargs into *merged_kwargs* in place and return it.

    For each rule: skip if the consumer already set the kwarg (explicit wins),
    skip if the artifact is not available upstream, otherwise inject it when the
    rule's name predicate matches or the consumer ``consumes`` the capability.
    """
    for rule in rules:
        if rule.inject_kwarg in merged_kwargs:
            continue
        key = rule.artifact_key or rule.capability
        if key not in current_artifacts:
            continue
        if rule.applies_to(step_name) or (rule.capability in consumes):
            merged_kwargs[rule.inject_kwarg] = str(current_artifacts[key])
    return merged_kwargs


def can_autowire(
    param_name: str,
    avail_caps: frozenset,
    *,
    rules: Tuple[WiringRule, ...] = DEFAULT_RULES,
) -> bool:
    """Whether *param_name* could be satisfied by a wiring rule.

    Used by the validator so it agrees with what :func:`autowire_kwargs` will
    actually do at run time — e.g. ``moread`` is not "missing" when a ``gbw``
    capability is available upstream.
    """
    return any(
        r.inject_kwarg == param_name and r.capability in avail_caps for r in rules
    )


__all__ = ["WiringRule", "DEFAULT_RULES", "autowire_kwargs", "can_autowire"]
