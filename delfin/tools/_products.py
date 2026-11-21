"""Declared products checked against delivered artifacts.

A step adapter declares the capability tags it ``produces``.  The wiring layer
trusts that declaration: an upstream step's artifacts are offered downstream
under the declared keys.  When an adapter mis-keys — declares ``("X",)`` but
puts its file under ``Y`` — the failure surfaces only at the consumer, far
from the producer, as a missing input capability with no hint that the
producer is the defect.

:meth:`check_declared_products` compares the declaration with what a
successful run actually left behind and returns the keys that are missing.
Two keys are satisfied outside the artifacts dict by design and are never
reported: ``geometry`` lives on ``StepResult.geometry`` and ``qc_output`` on
``StepResult.output_file`` (that mapping is what
``delfin.tools.pipeline._result_artifacts`` already relies on).

Legacy declarations remain report-only. Adapters can additionally declare
``required_artifacts``: missing, empty or unreadable files then fail the step.
Parsing and chemical validity remain the producer adapter's responsibility.
"""

from __future__ import annotations

from typing import Any, Dict, Tuple
from pathlib import Path

# Declared keys satisfied by StepResult fields rather than the artifacts dict.
# geometry: the adapter's output XYZ; qc_output: the engine's main log.
_SATISFIED_BY_RESULT_FIELD = {
    "geometry": ("geometry",),
    "qc_output": ("output_file",),
}

# The report's key inside StepResult.data.
DATA_KEY = "_missing_declared_products"


def check_declared_products(produced, result) -> Tuple[str, ...]:
    """Declared capability keys the run did not deliver, best-first order.

    *produced* is the adapter's declared ``produces`` set; *result* the
    successful :class:`~delfin.tools._types.StepResult`.
    """
    artifacts: Dict[str, Any] = result.artifacts or {}
    missing = []
    for key in sorted(produced):
        if key in artifacts:
            continue
        satisfied_by = _SATISFIED_BY_RESULT_FIELD.get(key, ())
        if any(getattr(result, field, None) for field in satisfied_by):
            continue
        missing.append(key)
    return tuple(missing)


def report_result(result, step_name: str, missing: Tuple[str, ...]) -> None:
    """Record a missing-products report on *result* and say it in the log."""
    if not missing:
        return
    result.data[DATA_KEY] = missing
    import logging
    logging.getLogger(__name__).warning(
        "step '%s' declares product(s) %s its result does not carry",
        step_name, sorted(missing),
    )


def invalid_required_artifacts(required, result) -> dict:
    """Opt-in physical output contract; content semantics belong to adapters."""
    invalid = {}
    for key in sorted(required):
        value = (result.artifacts or {}).get(key)
        if value is None:
            fields = _SATISFIED_BY_RESULT_FIELD.get(key, ())
            value = next((getattr(result, field, None) for field in fields
                          if getattr(result, field, None) is not None), None)
        try:
            if value is None:
                raise ValueError("missing artifact")
            path = Path(value)
            if not path.is_absolute() and result.work_dir is not None:
                path = Path(result.work_dir) / path
            if not path.is_file():
                raise ValueError("artifact is not a regular file")
            with path.open("rb") as handle:
                if not handle.read(1):
                    raise ValueError("artifact is empty")
        except (OSError, ValueError, TypeError) as exc:
            invalid[key] = str(exc)
    return invalid


__all__ = ["DATA_KEY", "check_declared_products", "report_result"]
