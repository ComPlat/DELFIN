"""Self-healing for ORCA OPI steps: classify a failure and escalate the retry.

The block classifies an ORCA failure from the ``.out`` text an OPI step wrote
(SCF not converged, geometry not converged, frequency failure, disk/memory)
and proposes an **escalation ladder** — a progressively stronger, but
chemically identical, retry.  It never changes method, basis, charge,
multiplicity or electronic state: an escalation *cannot* carry such a key, and
:func:`guard_unchanged_chemistry` refuses one that tries.  It is wired into
the ``orca_sp_opi`` / ``orca_opt_opi`` / ``orca_opt_freq_opi`` adapters as an
opt-in knob (``self_heal``, off by default).

The ladder mirrors ``delfin.orca_recovery.RecoveryStrategy``'s own escalation
(SlowConv -> VerySlowConv + damping -> SOSCF; fixed trust radius -> RecalcHess
-> looser trust), re-stated for the OPI block surface: every emitted keyword is
checked against ``delfin.orca_knowledge.known_keyword`` and every block option
against ``known_block_option`` before it is proposed.

Manual citations (ORCA 6.1, plain-text extraction under a read-only dir):
- ``! SlowConv`` / ``! VerySlowConv`` / damping: *essentialelements/scf.txt*
  2.6.2 ("many things can be done to improve convergence"; the damping block
  ``DampFac`` default 0.7, ``DampErr`` default 0.1, and the value 0.001 for
  difficult cases).
- ``! KDIIS``: *essentialelements/scf.txt* 2.6.5 ("Kolmar's DIIS … superior to
  the standard DIIS" for hard convergence).
- ``CNVSOSCF``: *essentialelements/scf.txt* 2.6.6 (Approximate Second Order
  SCF; the ``%scf`` switch disables DIIS and runs SOSCF when the DIIS error
  drops below ``SOSCFStart``).
- Geometry ``%geom Trust`` (negative = fixed trust radius), ``MaxIter``
  (default 200) and ``RecalcHess`` : *structurereactivity/optimizations.txt*
  (the geometry block reference table).
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence

#: Chemistry-surface keys an escalation must never touch (the task's "never
#: change method, basis, charge, multiplicity or electronic state").
_CHEMISTRY_KEYS = frozenset({
    "functional", "basis", "ri", "aux_basis", "dispersion",
    "method_keywords", "charge", "mult", "state",
})


class SelfHealFailure(Enum):
    """The failure classes the block recognises and can escalate."""

    SCF_NOT_CONVERGED = "scf_not_converged"
    GEOMETRY_NOT_CONVERGED = "geometry_not_converged"
    FREQUENCY_FAILURE = "frequency_failure"
    MEMORY_ERROR = "memory_error"
    NONE = "none"


# Marker texts ORCA prints, ordered specific-first so a mixed tail names the
# *aborting* cause, not a symptom.  Memory is checked before SCF because ORCA's
# out-of-memory abort mentions the SCF module.
_MARKERS: List[tuple[SelfHealFailure, Sequence[str]]] = [
    (
        SelfHealFailure.MEMORY_ERROR,
        (
            "not enough memory",
            "increase maxcore",
            "out-of-memory",
            "memory exhaustion",
        ),
    ),
    (
        SelfHealFailure.GEOMETRY_NOT_CONVERGED,
        (
            "maximum number of optimization cycles reached",
            "convergence failure",
            "geometry has not yet converged",
            "optimization has not converged",
        ),
    ),
    (
        SelfHealFailure.SCF_NOT_CONVERGED,
        (
            "scf not converged",
            "scf has not converged",
            "scf convergence failure",
        ),
    ),
    (
        SelfHealFailure.FREQUENCY_FAILURE,
        (
            "frequency calculation aborted",
            "aborted in module orca_freq",
            "failed to converged in the cphf/leancf",
        ),
    ),
]
# FLAG/BASIS group for editing chemistry is never produced here; kept as a
# sentinel so guard_unchanged_chemistry can distinguish absent from present.

#: Maximum ladder rungs (beyond this no escalation is proposed, only a halt).
MAX_ESCALATION_STEPS = 3
#: Job kinds that optimise a geometry (they may carry a %geom block).
_OPT_JOB_KINDS = frozenset({"OPT", "OPT_FREQ", "FREQ"})


@dataclass
class Escalation:
    """One rung of the retry ladder: the (smallest) change to retry with.

    ``scf_options`` / ``geom_options`` / ``keywords`` are *additions* the
    caller composes onto the original step options via
    :func:`compose_overrides`.  They never contain a chemistry key; that is
    enforced by :func:`plan_escalation` (the ladders emit none) and re-enforced
    by :func:`compose_overrides`.
    """

    failure: SelfHealFailure
    step: int
    scf_options: dict = field(default_factory=dict)
    geom_options: dict = field(default_factory=dict)
    keywords: list = field(default_factory=list)
    maxcore: Optional[int] = None
    rationale: str = ""


def classify_failure(output_text: str) -> SelfHealFailure:
    """Classify an ORCA output tail into a :class:`SelfHealFailure`.

    ``output_text`` may be freely formatted; matching is on lowercase marker
    substrings ORCA prints around a failure.  An empty or clean text returns
    :attr:`SelfHealFailure.NONE`.
    """
    if not output_text:
        return SelfHealFailure.NONE
    lowered = output_text.lower()
    for failure, markers in _MARKERS:
        if any(m in lowered for m in markers):
            return failure
    return SelfHealFailure.NONE


# ---------------------------------------------------------------------------
# Escalation ladders
# ---------------------------------------------------------------------------


def _scf_ladder(step: int) -> Escalation:
    # Manual *essentialelements/scf.txt* 2.6.2/2.6.5/2.6.6.  Every keyword and
    # %scf option is asserted catalog-known in the tests.
    tables: Dict[int, Dict[str, Any]] = {
        1: {
            "keywords": ["SlowConv"],
            "scf": {"maxiter": 300},
            "why": "SlowConv (manual scf 2.6.2) raises SCF iterations and "
                   "relaxes the damping; the first, mildest rung.",
        },
        2: {
            "keywords": ["VerySlowConv", "KDIIS"],
            "scf": {"maxiter": 500, "dampfac": 0.9, "damperr": 0.02},
            "why": "VerySlowConv + Kolmar DIIS (manual 2.6.5) plus high damping "
                   "(dampfac 0.9, manual 2.6.2) for a stubborn SCF.",
        },
        3: {
            "keywords": ["VerySlowConv"],
            "scf": {"maxiter": 800, "cnvsoscf": "true",
                    "dampfac": 0.95, "damperr": 0.001},
            "why": "Second-order SCF via CNVSOSCF (manual 2.6.6) with "
                   "very high damping (dampfac 0.95, damperr 0.001) — "
                   "the last SCF rung.",
        },
    }
    cfg = tables[step]
    return Escalation(
        failure=SelfHealFailure.SCF_NOT_CONVERGED,
        step=step,
        scf_options=dict(cfg["scf"]),
        keywords=list(cfg["keywords"]),
        rationale=cfg["why"],
    )


def _geometry_ladder(step: int) -> Escalation:
    # Manual *structurereactivity/optimizations.txt*: %geom Trust (negative =
    # fixed radius), MaxIter (default 200), RecalcHess.
    tables: Dict[int, Dict[str, Any]] = {
        1: {
            "geom": {"trust": -0.1, "maxiter": 250},
            "why": "Fixed small trust radius (trust -0.1, negative = fixed, "
                   "manual optimizations.txt) to stop the geometry overshooting "
                   "an oscillating energy; more cycles (maxiter 250, default 200).",
        },
        2: {
            "geom": {"trust": -0.05, "recalc_hess": 10, "maxiter": 300},
            "why": "Smaller fixed trust radius plus Hessian recalculation every "
                   "10 steps (recalc_hess, manual) for a difficult PES.",
        },
        3: {
            "geom": {"trust": -0.05, "recalc_hess": 5, "maxiter": 400},
            "why": "Frequent Hessian updates (recalc_hess 5) and the widest cycle "
                   "budget (maxiter 400); the last geometry rung.",
        },
    }
    cfg = tables[step]
    return Escalation(
        failure=SelfHealFailure.GEOMETRY_NOT_CONVERGED,
        step=step,
        geom_options=dict(cfg["geom"]),
        rationale=cfg["why"],
    )


def _frequency_ladder(step: int) -> Escalation:
    # The analytic frequency's embedded CP-SCF / LEANSCF fails; escalate by
    # tightening that SCF (mirrors orca_recovery._leanscf_convergence_fixes).
    tables: Dict[int, Dict[str, Any]] = {
        1: {
            "keywords": ["SlowConv"],
            "scf": {"maxiter": 400},
            "why": "Frequency's coupled-perturbed SCF did not converge; "
                   "tighten it with SlowConv and more iterations.",
        },
        2: {
            "keywords": ["VerySlowConv", "KDIIS"],
            "scf": {"maxiter": 600, "dampfac": 0.9, "damperr": 0.02},
            "why": "VerySlowConv + Kolmar DIIS + high damping for the "
                   "frequency's SCF.",
        },
        3: {
            "keywords": ["VerySlowConv"],
            "scf": {"maxiter": 1000, "cnvsoscf": "true",
                    "dampfac": 0.95, "damperr": 0.001},
            "why": "Second-order SCF with very heavy damping; the last "
                   "frequency SCF rung.",
        },
    }
    cfg = tables[step]
    return Escalation(
        failure=SelfHealFailure.FREQUENCY_FAILURE,
        step=step,
        scf_options=dict(cfg["scf"]),
        keywords=list(cfg["keywords"]),
        rationale=cfg["why"],
    )


def _memory_ladder(step: int, maxcore: Optional[int] = None) -> Escalation:
    # Memory is not fixed by a block option but by giving the job more MaxCore
    # (Oracle's out-of-memory message names the figure needed).  Halving the
    # core count is *not* done here: it can cross a PAL dependency the user set.
    base = maxcore or 1000
    return Escalation(
        failure=SelfHealFailure.MEMORY_ERROR,
        step=step,
        maxcore=base * (2**step),
        rationale="ORCA ran out of memory; retry with more MaxCore "
                  "(doubled per rung). Chemistry untouched.",
    )


_LADDERS = {
    SelfHealFailure.SCF_NOT_CONVERGED: _scf_ladder,
    SelfHealFailure.GEOMETRY_NOT_CONVERGED: _geometry_ladder,
    SelfHealFailure.FREQUENCY_FAILURE: _frequency_ladder,
    SelfHealFailure.MEMORY_ERROR: _memory_ladder,
}


def plan_escalation(
    failure: SelfHealFailure,
    step: int,
    *,
    job_kind: str = "SP",
    scf_options: Optional[Mapping[str, Any]] = None,
    geom_options: Optional[Mapping[str, Any]] = None,
    maxcore: Optional[int] = None,
) -> Escalation:
    """The escalation for *failure* on ladder *step* (1-based, <= 3).

    Raises ``ValueError`` for a step outside ``1..MAX_ESCALATION_STEPS`` and
    for a geometry escalation on a job that does not optimise a geometry
    (``job_kind`` an SP).  The returned :class:`Escalation` carries only
    additive, catalog-known changes; it never mentions method/basis/charge/
    mult/state.
    """
    if not isinstance(step, int) or not 1 <= step <= MAX_ESCALATION_STEPS:
        raise ValueError(
            f"self-heal escalation step must be 1..{MAX_ESCALATION_STEPS}, got {step!r}")
    if failure is SelfHealFailure.NONE:
        raise ValueError("no escalation exists for a clean / unclassified run")
    if (failure is SelfHealFailure.GEOMETRY_NOT_CONVERGED
            and job_kind not in _OPT_JOB_KINDS):
        raise ValueError(
            "geometry self-heal escalation requires an optimizing job "
            f"(got job_kind={job_kind!r}); a single point has no geometry to fix")
    builder = _LADDERS[failure]
    if failure is SelfHealFailure.MEMORY_ERROR:
        return builder(step, maxcore=maxcore)
    return builder(step)


def compose_overrides(
    user_options: Optional[Mapping[str, Any]],
    additions: Optional[Mapping[str, Any]],
) -> dict:
    """Merge *additions* onto *user_options* into a new dict, escalation wins.

    The input dicts are never mutated.  A chemistry key in *additions*
    (functional / basis / charge / ...) is deliberately dropped — the block
    must not change chemistry, even if a rogue caller passes one.
    """
    out = dict(user_options or {})
    for key, value in (additions or {}).items():
        if str(key).strip().lower() in _CHEMISTRY_KEYS:
            continue
        out[key] = copy.deepcopy(value)
    return out


def guard_unchanged_chemistry(base: Mapping[str, Any], **escalated: Any) -> None:
    """Refuse an escalation that would change method / basis / charge / mult.

    ``base`` is the step's original chemistry surface (``functional``,
    ``basis``, ``ri``, ``aux_basis``, ``dispersion``, ``method_keywords``,
    ``charge``, ``mult``, ``state``).  Each key in ``escalated`` that names the
    chemistry surface must equal the base value; a differing one raises
    ``ValueError`` naming the key.  ``scf_options`` / ``geom_options`` /
    ``keywords`` are never chemistry and are ignored.  Used by the retry loop
    before every re-launch so no escalation can silently alter the chemistry.
    """
    for key, value in escalated.items():
        if str(key).strip().lower() not in _CHEMISTRY_KEYS:
            continue
        k = str(key).strip().lower()
        if base.get(k, base.get(key)) != value:
            raise ValueError(
                f"self-heal escalation would change {k} "
                f"(base {base.get(k, base.get(key))!r} -> {value!r}); refused. "
                "Self-healing never alters method, basis, charge, "
                "multiplicity or electronic state.")


def build_attempt_log(
    plans: Sequence[Escalation],
) -> List[Dict[str, Any]]:
    """Compose the per-attempt retry log an adapter surfaces in its result.

    Each entry is a typed dict: the failure class string, the ladder step, the
    additive keyword / ``%scf`` / ``%geom`` changes and the manual rationale.
    """
    out: List[Dict[str, Any]] = []
    for plan in plans:
        out.append({
            "failure": plan.failure.value,
            "step": plan.step,
            "keywords": list(plan.keywords),
            "scf_options": dict(plan.scf_options),
            "geom_options": dict(plan.geom_options),
            "maxcore": plan.maxcore,
            "rationale": plan.rationale,
        })
    return out


def run_with_self_heal(
    *,
    run_once,
    params: Mapping[str, Any],
    job_kind: str,
    max_attempts: int = 3,
) -> tuple[Any, Optional[List[Dict[str, Any]]]]:
    """Re-run a failed OPI step, escalating per the ladder, chemistry intact.

    ``run_once(**params)`` executes one attempt and returns a result with an
    ``.ok`` flag and (on failure) an ``.output_file`` the classifier reads.
    Attempt 0 uses ``params`` verbatim; each later attempt composes the ladder
    rung's ``scf_options`` / ``geom_options`` / ``keywords`` / ``maxcore`` onto
    the originals and re-runs, after :func:`guard_unchanged_chemistry`
    confirms the chemistry surface is untouched.  Returns
    ``(final_result, attempt_log)``; ``attempt_log`` is a list of typed
    per-retry records (empty when no retry ran), for the adapter to surface in
    its ``StepResult.data["_self_heal"]``.
    """
    current: dict = dict(params)
    plans: List[Escalation] = []
    result = run_once(**current)
    if result.ok or max_attempts <= 1:
        return result, None
    for step in range(1, max_attempts + 1):
        failure = classify_failure(_read_tail(result))
        if failure is SelfHealFailure.NONE:
            break
        try:
            plan = plan_escalation(failure, step, job_kind=job_kind,
                                   maxcore=current.get("maxcore"))
        except ValueError:
            break
        # The next attempt must keep the chemistry surface identical.
        guard_unchanged_chemistry(
            base=current,
            scf_options=plan.scf_options,
            geom_options=plan.geom_options,
            keywords=plan.keywords,
        )
        plans.append(plan)
        if plan.keywords:
            current["extra_keywords"] = list(
                current.get("extra_keywords") or ()) + list(plan.keywords)
        if plan.scf_options:
            current["scf_options"] = compose_overrides(
                current.get("scf_options"), plan.scf_options)
        if plan.geom_options:
            current["geom_options"] = compose_overrides(
                current.get("geom_options"), plan.geom_options)
        if plan.maxcore is not None:
            current["maxcore"] = plan.maxcore
        result = run_once(**current)
        if result.ok:
            break
    return result, (build_attempt_log(plans) if plans else None)


def _read_tail(result: Any) -> str:
    """The ``.out`` text of a result's output file, or the empty string."""
    out = getattr(result, "output_file", None)
    if not out:
        return ""
    try:
        return Path(out).read_text(errors="replace")
    except OSError:
        return ""


__all__ = [
    "MAX_ESCALATION_STEPS",
    "SelfHealFailure",
    "Escalation",
    "classify_failure",
    "plan_escalation",
    "compose_overrides",
    "guard_unchanged_chemistry",
    "build_attempt_log",
    "run_with_self_heal",
]
