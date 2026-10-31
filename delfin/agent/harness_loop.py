"""Harness^2 improvement loop (package HX1, wave 25).

One improvement round for the agent's own harness, after the method of
Harness^2: Recursive Agent Harnessing for an Open World (Xu et al.,
arXiv:2610.14224, 2026, doi:10.48550/arXiv.2610.14224):

* §2.2.1 "A Contrastive Improvement Step": propose candidate edits to
  the harness, probe them under MATCHED conditions (same tasks, same
  seeds, same model), compose only where the behavioural contrast
  supports the edit.  The paper stresses that no ground-truth verifier
  is needed — the CONTRAST between executions decides.
* §2.3 "Harnessing-Level Recursion": the improvement experience is
  recorded as evidence-gated lessons (here: KN1-style records shaped
  like :class:`delfin.agent.knowledge.Entry` — situation / lesson / why
  with a provenance evidence pointer), so later rounds reuse earlier
  findings instead of re-deriving them.
* §3.2: weaker executors gain most from harness edits; gains must be
  separated from run-to-run noise, which DELFIN measured directly
  (bench_watch: identical runs flip 6/48 tasks, quality stdev ~20).

Safety invariants (package HX1 brief):
* Candidates are TEXT diffs to a skill file or tool-guidance text
  only.  A candidate touching code, a protected path, or the
  pack/knowledge tree is refused before it is ever probed.
* NOTHING is applied automatically.  Adoption is a PR the person
  merges, accepted only if the family rises and no other family falls
  (holdout check happens at adoption time, outside this module).
* Every decision in the report cites an observation (task id + what
  went wrong) from a real reference transcript.
* Real mode states its cost up front and refuses above the budget;
  dry-run mode makes no model calls and replays recorded transcripts
  from fixtures.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

__all__ = [
    "pick_weakest_family", "DEFAULT_CEILING",
    "Candidate", "HarnessCandidateError", "validate_candidate",
    "ALLOWED_SUFFIXES", "PROTECTED_MARKERS",
    "contrast_supported",
    "RoundConfig", "RoundReport", "run_round",
]

#: Classic workflow families scored 91-95 in the 31.10. E2 baseline
#: (round 2); a family at or above this ceiling has no headroom and
#: must not be picked (brief, target item 1).
DEFAULT_CEILING = 90.0


def pick_weakest_family(
        scores: Mapping[str, float],
        *,
        ceiling: float = DEFAULT_CEILING) -> str:
    """Return the id of the weakest family that still has headroom.

    ``scores`` maps family id to a 0-100 quality score.  Families at or
    above ``ceiling`` are skipped (a family at the ceiling cannot gain
    from a harness edit and probing it would only burn budget).  When
    NO family has headroom, :class:`LookupError` is raised — the caller
    records that as a valid outcome of the round, not as an error.
    """
    if not scores:
        raise ValueError("no family scores given")
    # Reviewer F4: a NaN/None score is UNMEASURED data, not a family
    # "at the ceiling" — refusing beats silently skipping (the skip
    # manufactured a false LookupError claim about the data) and beats
    # a TypeError from the comparison.
    for fid, score in scores.items():
        try:
            bad = score is None or not math.isfinite(float(score))
        except (TypeError, ValueError):
            bad = True
        if bad:
            raise ValueError(
                "family %r has a non-numeric score %r; refusing to pick "
                "a target from unmeasured data" % (fid, score))
    headroom = [(score, fid) for fid, score in scores.items()
                if score < ceiling]
    if not headroom:
        raise LookupError(
            "every family is at or above the ceiling %s; no target "
            "with headroom" % ceiling)
    # min() over (score, id): ties fall back to the first family id in
    # iteration order, deterministic for a dict built from a list.
    return min(headroom)[1]


# ---------------------------------------------------------------------------
# Candidate edits: TEXT diffs to guidance/skill text only (brief item 3)
# ---------------------------------------------------------------------------

#: A candidate may only touch a text file that carries guidance: skills
#: and tool-guidance markdown.  Everything else — code, data, configs —
#: is refused: the paper's edits are harness TEXT, not code (§2.2.1).
ALLOWED_SUFFIXES = (".md", ".txt", ".rst")

#: The enumerated ALLOWED ROOTS (reviewer F1/F2, round 1: the deny-list
#: of suffix+markers let traversal and unrelated .md files through; an
#: allow-list is the fix, not a longer deny-list).  A candidate target
#: must be repo-relative, contain no "..", and live under one of these
#: guidance directories.
ALLOWED_ROOTS = (
    "delfin/agent/pack/skills/",
    "delfin/tools/guidance/",
)

#: Protected markers stay as a second, independent check on top of the
#: allow-list (defence in depth; the roots above are already outside
#: every protected tree, but a future root edit must not silently
#: widen into one).
PROTECTED_MARKERS = (
    "sandbox", "permission", "approval", "api_client", "engine.py",
    "benchmark_runner", "knowledge", "skill_safety", "skill_share",
    "mcp_", "web_tools", "text_sandbox", "tests/", "fixtures/",
)

#: The two candidate kinds the paper composes from (brief item 3):
#: one invariant (planning/tool use/context handling), one domain
#: (the family's procedure).
CANDIDATE_KINDS = ("invariant", "domain")


class HarnessCandidateError(ValueError):
    """A candidate edit is not a text-diff to an allowed target."""


@dataclass(frozen=True)
class Candidate:
    """One proposed harness edit: a TEXT diff plus its observation."""

    kind: str              # "invariant" or "domain"
    target: str            # repo-relative path of the guidance file
    old_text: str
    new_text: str
    observation_task_id: str
    observation: str       # what went wrong in the reference transcript


def validate_candidate(candidate: Candidate) -> None:
    """Refuse every candidate the brief forbids, BEFORE any probing.

    Refusals (all :class:`HarnessCandidateError`):
    * target not a guidance text file (code, data, notebooks);
    * target not repo-relative, containing "..", absolute, or outside
      the enumerated allowed guidance roots (allow-list — reviewer F1/F2);
    * target under a protected path (independent second check);
    * a no-op diff: empty or unchanged new_text (F3 — a zero-diff probe
      burns budget for a guaranteed no-contrast);
    * no observation (the report must cite one for every decision);
    * kind not in ("invariant", "domain").
    """
    if not candidate.target or not candidate.target.endswith(
            ALLOWED_SUFFIXES):
        raise HarnessCandidateError(
            "candidate target %r is not a guidance text file "
            "(allowed suffixes: %s)" % (candidate.target, ALLOWED_SUFFIXES))
    if candidate.target.startswith("/") or ".." in candidate.target.split("/"):
        raise HarnessCandidateError(
            "candidate target %r must be a repo-relative path without "
            "traversal" % candidate.target)
    if not any(candidate.target.startswith(root)
               for root in ALLOWED_ROOTS):
        raise HarnessCandidateError(
            "candidate target %r is outside the allowed guidance roots "
            "%s; the loop edits only enumerated guidance files"
            % (candidate.target, ALLOWED_ROOTS))
    lowered = candidate.target.lower()
    for marker in PROTECTED_MARKERS:
        if marker in lowered:
            raise HarnessCandidateError(
                "candidate target %r hits the protected marker %r; "
                "the safety frame is not editable by the loop"
                % (candidate.target, marker))
    if not (candidate.new_text or "").strip():
        raise HarnessCandidateError(
            "candidate new_text is empty; a diff that changes nothing is "
            "not an edit")
    if candidate.new_text == candidate.old_text:
        raise HarnessCandidateError(
            "candidate is a no-op (old_text == new_text); probing a "
            "zero-diff burns budget for a guaranteed no-contrast")
    if not (candidate.observation_task_id or "").strip():
        raise HarnessCandidateError(
            "candidate has no observation task id; every edit must cite "
            "what went wrong in a reference transcript")
    if not (candidate.observation or "").strip():
        raise HarnessCandidateError(
            "candidate has no observation text; every edit must cite "
            "what went wrong in a reference transcript")
    if candidate.kind not in CANDIDATE_KINDS:
        raise HarnessCandidateError(
            "candidate kind %r is not one of %s (one invariant, one "
            "domain per round)" % (candidate.kind, CANDIDATE_KINDS))


# ---------------------------------------------------------------------------
# Noise-aware contrast (brief item 4; QS requirement 1)
# ---------------------------------------------------------------------------

def contrast_supported(
        task_ids: list[str],
        reference: list[float],
        probe: list[float],
        *,
        spread: float) -> dict:
    """Pair reference and probe qualities and test against the noise band.

    The band is DELFIN's own measured threshold, taken from
    :mod:`bench_watch` — the same ``max(8, 1.5*spread)`` the nightly
    watcher requires before it calls a quality drop a regression
    (bench_watch.py:50 QUALITY_DROP_FLOOR, bench_watch.py:57
    SPREAD_FACTOR, bench_watch.py:70 suspect_threshold).  Identical
    runs flip 6/48 tasks with a quality stdev of ~20 (bench_watch.py:7),
    so a raw point estimate can never decide.

    Returns ``{"supported": bool, "band": float, "mean_gain": float,
    "per_task": [(task_id, ref, probe, gain), ...]}``.  ``supported``
    is True only when the mean gain is STRICTLY above the band.  The
    band and the mean gain travel with the verdict so the report can
    cite them (a verdict that hid its band would be a point estimate
    in disguise).

    Raises :class:`ValueError` when the paired lists are empty or do
    not match — unmatched samples give no contrast, and inventing one
    would be fabrication.
    """
    if not task_ids or not reference:
        raise ValueError("no paired reference samples")
    if len(task_ids) != len(reference) or len(reference) != len(probe):
        raise ValueError(
            "reference and probe must be paired per task "
            "(got %d tasks, %d reference, %d probe values)"
            % (len(task_ids), len(reference), len(probe)))
    # Reviewer F5: a non-finite sample is an UNMEASURED value — letting
    # it through would manufacture a clean-looking "no contrast" from
    # nothing.  Refuse loudly.
    for where, values in (("reference", reference), ("probe", probe)):
        for value in values:
            if not math.isfinite(float(value)):
                raise ValueError(
                    "non-finite %s quality %r; refusing to build a "
                    "contrast from unmeasured data" % (where, value))
    per_task = [(tid, ref, pr, pr - ref)
                for tid, ref, pr in zip(task_ids, reference, probe)]
    mean_gain = sum(gain for *_, gain in per_task) / len(per_task)
    from .bench_watch import suspect_threshold
    band = suspect_threshold(spread)
    return {
        "supported": bool(mean_gain > band),
        "band": band,
        "mean_gain": mean_gain,
        "per_task": per_task,
    }


# ---------------------------------------------------------------------------
# The improvement round (brief items 5 and 6)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RoundConfig:
    """One improvement round's inputs.

    In dry-run mode the paired quality samples come from RECORDED
    transcripts (fixtures); no runner is imported and no model is
    called.  In real mode the caller has already run the reference and
    probe rounds through bench_nightly's path and hands the paired
    qualities in — this module measures and composes, it never drives
    the model stack itself (one module, one responsibility; the driver
    is the operator's nightly orchestration).
    """

    family_id: str
    candidate_edits: dict  # candidate id -> Candidate
    paired_tasks: list  # task ids in matched order
    reference: list  # recorded reference qualities
    probe_quality: dict  # candidate id -> recorded probe qualities
    spread: float = 0.0
    estimated_cost_usd: float = 0.0
    budget_usd: float = float("inf")


class RoundReport:
    """Namespace for the round's outcome record and its budget error.

    The report is a plain, JSON-serialisable dict (the caller writes
    it wherever the review flow reads it); nothing here opens a file,
    touches a model, or applies a change.
    """

    class BudgetRefused(Exception):
        """The round's estimated cost exceeds the stated budget."""


def run_round(config: RoundConfig, *, dry_run: bool) -> dict:
    """Run one Propose-Probe-Compose round and return the report.

    Every candidate is probed against the SAME paired task ids and the
    SAME reference qualities (matched conditions, arXiv:2610.14224
    §2.2.1).  Only candidates whose contrast is supported are composed;
    a "no contrast" is recorded as a valid, named outcome — never
    hidden, never counted as a gain.  Each composed candidate becomes
    one KN1-style lesson record (fields shaped like
    :class:`delfin.agent.knowledge.Entry`: situation / lesson / why /
    provenance with an evidence pointer, kinds per knowledge.py:96).

    The report is a PROPOSAL.  This function applies nothing: adoption
    is a PR the person merges, gated at adoption time on the family
    rising and no other family falling (holdout check, outside this
    module).
    """
    if config.estimated_cost_usd > config.budget_usd:
        raise RoundReport.BudgetRefused(
            "round estimated at %.2f USD exceeds the stated budget of "
            "%.2f USD; refusing to start" % (config.estimated_cost_usd,
                                             config.budget_usd))
    # Reviewer F6: a round with no candidates is not a measured round —
    # an empty report would present noise_band=None as if something was
    # probed.  Refuse instead.
    if not config.candidate_edits:
        raise ValueError(
            "no candidate edits given; a round without candidates is "
            "not a measured round (report would carry noise_band=None)")
    composed = []
    no_contrast = []
    lessons = []
    band = None
    for cid, candidate in config.candidate_edits.items():
        validate_candidate(candidate)
        if cid not in config.probe_quality:
            raise KeyError(
                "no probe qualities recorded for candidate %r; the probe "
                "must run before composition" % cid)
        verdict = contrast_supported(
            config.paired_tasks, config.reference,
            list(config.probe_quality[cid]), spread=config.spread)
        band = verdict["band"]
        decision = {
            "id": cid,
            "kind": candidate.kind,
            "target": candidate.target,
            "observation_task_id": candidate.observation_task_id,
            "observation": candidate.observation,
            "mean_gain": verdict["mean_gain"],
            "band": verdict["band"],
            "supported": verdict["supported"],
        }
        if verdict["supported"]:
            composed.append(decision)
            lessons.append(_lesson_for(candidate, cid, verdict))
        else:
            no_contrast.append(decision)
    return {
        "mode": "dry-run" if dry_run else "real",
        "family_id": config.family_id,
        "noise_band": band,
        "spread": config.spread,
        "composed": composed,
        "no_contrast": no_contrast,
        "lessons": lessons,
        "adoption": "PR for the person to review; nothing is applied "
                    "automatically",
    }


def _lesson_for(candidate: Candidate, cid: str, verdict: dict) -> dict:
    """KN1-style lesson record for one composed candidate."""
    return {
        "topic": "harness-improvement",
        "situation": (
            "Workflow family tasks like %s: %s"
            % (candidate.observation_task_id, candidate.observation)),
        "lesson": (
            "Apply this harness text edit to %s: replace %r with %r "
            "(candidate %s)" % (candidate.target,
                                candidate.old_text,
                                candidate.new_text, cid)),
        "why": (
            "Matched-condition probe over the same instances gave a mean "
            "gain of %.1f points against a noise band of %.1f "
            "(stdev %.1f) — outside the band"
            % (verdict["mean_gain"], verdict["band"], verdict["band"] / 1.5
               if verdict["band"] > 8.0 else 0.0)),
        "provenance": {
            "evidence": "probe run of candidate %s on tasks %s"
                        % (cid, ",".join(map(str, verdict["per_task"]))),
            "evidence_kind": "test",
        },
    }
