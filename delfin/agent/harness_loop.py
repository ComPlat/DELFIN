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

from typing import Mapping

__all__ = ["pick_weakest_family", "DEFAULT_CEILING"]

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
    headroom = [(score, fid) for fid, score in scores.items()
                if score < ceiling]
    if not headroom:
        raise LookupError(
            "every family is at or above the ceiling %s; no target "
            "with headroom" % ceiling)
    # min() over (score, id): ties fall back to the first family id in
    # iteration order, deterministic for a dict built from a list.
    return min(headroom)[1]
