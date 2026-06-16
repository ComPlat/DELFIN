"""Error classification + recovery policies.

Maps a failure to *why* it failed (an :class:`~delfin.tools._types.ErrorKind`)
and *what to do about it* (retry, retry-with-modified-params, or give up), so
recovery is a reusable platform policy instead of hand-coded per caller.  Works
even when adapters do not set ``error_kind`` themselves, by classifying the
free-text error message.

    result = run_step_with_recovery("orca_opt", geometry=geom, charge=0)
    # on an SCF convergence failure this automatically retries with more SCF
    # iterations + SOSCF, before giving up.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Optional

from delfin.tools._types import ErrorKind

# Ordered (specific first) — first match wins.
_PATTERNS = [
    (ErrorKind.CONVERGENCE, re.compile(
        r"scf not converged|did not converge|not converged|convergence (?:failed|problem)",
        re.I)),
    (ErrorKind.TIMEOUT, re.compile(r"timed out|timeout", re.I)),
    (ErrorKind.BINARY_NOT_FOUND, re.compile(
        r"command not found|not found on path|executable not found|no such file or directory",
        re.I)),
    (ErrorKind.MISSING_PARAM, re.compile(r"parameter is required|is required\b", re.I)),
    (ErrorKind.PARSE, re.compile(r"could not parse|parse error|failed to parse", re.I)),
]


def classify_error(text: str) -> ErrorKind:
    """Best-effort classification of a free-text error message."""
    if not text:
        return ErrorKind.NONE
    for kind, pat in _PATTERNS:
        if pat.search(text):
            return kind
    return ErrorKind.TOOL_FAILED


@dataclass
class RecoveryAction:
    action: str                                    # "retry" | "modify" | "give_up"
    reason: str = ""
    kwargs_overrides: Dict[str, Any] = field(default_factory=dict)


def suggest_recovery(kind: ErrorKind, kwargs: Optional[Dict[str, Any]] = None) -> RecoveryAction:
    """Recommend a recovery action for a failure of category *kind*."""
    kwargs = kwargs or {}
    if kind is ErrorKind.CONVERGENCE:
        scf_maxiter = max(int(kwargs.get("scf_maxiter", 0) or 0), 300)
        extra = str(kwargs.get("scf_extra", "") or "")
        if "soscf" not in extra.lower():
            extra = (extra + " SOSCF true").strip()
        return RecoveryAction(
            "modify",
            reason="SCF did not converge — raising SCF iterations and enabling SOSCF",
            kwargs_overrides={"scf_maxiter": scf_maxiter, "scf_extra": extra},
        )
    if kind is ErrorKind.TIMEOUT:
        return RecoveryAction("retry", reason="transient timeout — retrying")
    if kind in (ErrorKind.BINARY_NOT_FOUND,):
        return RecoveryAction(
            "give_up",
            reason="a required tool is not installed — see platform.install_plan()",
        )
    if kind in (ErrorKind.MISSING_PARAM, ErrorKind.MISSING_INPUT, ErrorKind.PARSE,
                ErrorKind.INTERNAL):
        return RecoveryAction("give_up", reason=f"non-recoverable ({kind.value})")
    # TOOL_FAILED / NONE / unknown — one plain retry is reasonable.
    return RecoveryAction("retry", reason="generic failure — retrying once")


def run_step_with_recovery(
    step_name: str,
    *,
    geometry=None,
    cores: int = 1,
    work_dir=None,
    max_attempts: int = 3,
    **kwargs: Any,
):
    """Run a step, applying the recovery policy on failure up to *max_attempts*.

    Returns the final :class:`StepResult`; its ``data["_recovery"]`` carries the
    per-attempt log (error kind, chosen action, reason).
    """
    from delfin.tools._runner import run_step

    current = dict(kwargs)
    attempts = []
    result = None
    for i in range(max_attempts):
        wd = (Path(work_dir) / f"attempt_{i}") if work_dir else None
        result = run_step(step_name, geometry=geometry, cores=cores, work_dir=wd, **current)
        entry: Dict[str, Any] = {"attempt": i, "ok": result.ok, "error": result.error}
        if result.ok:
            attempts.append(entry)
            result.data["_recovery"] = attempts
            return result

        kind = result.error_kind if result.error_kind is not ErrorKind.NONE \
            else classify_error(result.error or "")
        action = suggest_recovery(kind, current)
        entry.update({"error_kind": kind.value, "action": action.action,
                      "reason": action.reason})
        attempts.append(entry)

        if action.action == "give_up" or i == max_attempts - 1:
            break
        if action.action == "modify":
            current.update(action.kwargs_overrides)
        # "retry" → loop with unchanged params

    if result is not None:
        result.data["_recovery"] = attempts
    return result


__all__ = [
    "classify_error",
    "RecoveryAction",
    "suggest_recovery",
    "run_step_with_recovery",
]
