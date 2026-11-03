"""Report writer for the nightly holdout runner (PL4 phase 2).

Consumes the RunRecord that :func:`bench_nightly.run` produces plus the
history of earlier rounds' records, and emits the round report:

* score per tool family (task_class without the correctness suffix) and
  per correctness level (consistency, reproduction, reference, trap —
  CB0B's task families carry the level as the task_class suffix),
* cost next to every score (core-hours, EF1 accounting),
* the weakest tasks as a queue for the next round, weakest first,
* the stop verdict: two rounds without gain, or the daily budget
  exhausted.

The report is JSON text and a pure function of its arguments: the same
record and history render the same bytes (no timestamps from ``write``
itself; the record carries its own clock stamps). No secret ever enters
it — the seed never lands in the record, so it cannot land here.
"""

from __future__ import annotations

import json

from .bench_nightly import InstanceRecord, RunRecord

CORRECTNESS_LEVELS = ("consistency", "reproduction", "reference", "trap")
#: A task goes back into the next round's queue below this quality.
WEAKEST_QUALITY_THRESHOLD = 60
#: Stop after this many consecutive rounds without gain.
FLAT_ROUNDS_TO_STOP = 2


def _split_family_level(task_class: str) -> tuple[str, str]:
    """``family-a_consistency`` -> (``family-a``, ``consistency``).

    The family is everything before the last ``_`` WHEN that last part
    is a correctness level; anything else is an all-family class with
    no level.
    """
    if "_" in task_class:
        head, _, tail = task_class.rpartition("_")
        if tail in CORRECTNESS_LEVELS:
            return head, tail
    return task_class, ""


def _score(instances: list[InstanceRecord]) -> dict:
    """Mean quality / success over the instances that ran and graded.

    Unmeasured instances (the run never happened: engine-init failure,
    HB2 D2) are excluded — a dead run is not a quality-0 failure of
    the model, and scoring one fabricates a measurement.
    """
    graded = [i for i in instances
              if i.ran and not i.unmeasured
              and i.success is not None and i.quality is not None]
    if not graded:
        return {"quality": None, "success": None, "n": 0,
                "core_hours": 0.0}
    return {
        "quality": round(sum(i.quality for i in graded) / len(graded)),
        "success": all(i.success for i in graded),
        "n": len(graded),
        "core_hours": round(sum(i.core_hours for i in graded), 3),
    }


def _stop_verdict(record: RunRecord, history: list[RunRecord]) -> dict:
    """Two consecutive rounds without gain, or the budget exhausted.

    HB2: a history record that measured nothing (every started
    instance unmeasured — e.g. the 31.10. engine-init baseline) is
    INVALID as evidence and is skipped: its flat/failed outcome must
    not push the stop verdict. The mark is computed, not stored, so
    an old broken record is recognised without anyone editing it.
    """
    gains = list(record.previous_rounds_gain)
    if history:
        last = history[-1]
        if _all_unmeasured(last):
            gains = []  # the "last round" measured nothing; no gains to chain
        else:
            gains = list(last.previous_rounds_gain) + gains
    flat_tail = 0
    for gain in reversed(gains):
        if gain and gain > 0:
            break
        flat_tail += 1
    if flat_tail >= FLAT_ROUNDS_TO_STOP:
        return {"should_stop": True,
                "reason": "no gain in %d consecutive rounds" % flat_tail}
    if record.spent_core_hours >= record.daily_budget_core_hours:
        return {"should_stop": True,
                "reason": "budget exhausted: %.3f of %.3f core-hours"
                          % (record.spent_core_hours,
                             record.daily_budget_core_hours)}
    return {"should_stop": False, "reason": ""}


def _all_unmeasured(record: RunRecord) -> bool:
    """True when the round started but every started instance is
    marked unmeasured (HB2): the round measured nothing.

    Deliberately ONLY the explicit flag: a pre-HB2 record (all
    success=False, empty error) is indistinguishable from a real
    failed round that measured quality 0 — the flag did not exist
    then, and guessing marked real failed rounds invalid and broke
    the flat-round stop verdict (caught by test_pl4_report). The
    31.10. round is identified by the person, not by a heuristic.
    """
    started = record.started()
    if not started:
        return False
    return all(getattr(i, "unmeasured", False) for i in started)


def write(record: RunRecord, history: list[RunRecord] | None = None) -> str:
    """The round report as JSON text: a pure function of its arguments."""
    history = list(history or [])
    ran = [i for i in record.instances if i.ran]
    # HB2 D2: unmeasured instances are named, never scored.
    unmeasured = [i for i in ran if i.unmeasured]

    families: dict[str, dict] = {}
    levels: dict[str, dict] = {}
    for inst in ran:
        if inst.unmeasured:
            continue
        family, level = _split_family_level(inst.task_class)
        families.setdefault(family, []).append(inst)
        if level:
            levels.setdefault(level, []).append(inst)

    weakest = sorted(
        (i for i in ran
         if i.quality is not None and i.quality < WEAKEST_QUALITY_THRESHOLD),
        key=lambda i: i.quality)
    report = {
        "model": record.model,
        "mode": record.mode,
        "holdout_ran": record.holdout_ran,
        "seed_log": record.seed_log,
        "refused_over_budget": record.refused_over_budget,
        "unmeasured": [
            {"task_id": i.task_id, "task_class": i.task_class,
             "error": i.error}
            for i in unmeasured],
        "families": {f: _score(insts) for f, insts in sorted(families.items())},
        "correctness_levels": {lv: _score(insts)
                               for lv, insts in sorted(levels.items())},
        "cost": {
            "spent_core_hours": record.spent_core_hours,
            "daily_budget_core_hours": record.daily_budget_core_hours,
        },
        "weakest_queue": [
            {"task_id": i.task_id, "task_class": i.task_class,
             "quality": i.quality, "success": i.success}
            for i in weakest],
        "stop": _stop_verdict(record, history),
        "rounds_included": len(history) + 1,
    }
    # HB2: an invalid history record (measured nothing) is marked,
    # never deleted — a later reader sees why it was skipped.
    if history and _all_unmeasured(history[-1]):
        report["history_invalid"] = (
            "the previous round measured nothing (every started "
            "instance unmeasured or a pre-HB2 all-failed/empty-error "
            "record); its gains are excluded from the stop verdict")
    return json.dumps(report, sort_keys=True)
