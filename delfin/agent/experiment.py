"""Experiments the agent cannot fool itself with.

A candidate change to the agent sits behind ONE switch, default off.  To
claim it worked, an agent has to run an *experiment*: a pre-registered,
stamped, reach-checked A-vs-B comparison whose every blocker is classified
and that lands only through an explicit human approval record.

The whole point is that the agent cannot quietly redefine "worked".  The
refusals below are the mechanism: measuring without pre-registration is an
error; comparing two results measured under different instruments is an
error; a switch-off that is not byte-identical to baseline aborts; landing
without a human approval record is impossible.

This module is generic -- no paths, no server names, no run names, no data
sets.  It ships three layers, one per concern:

* :class:`Experiment`  -- the record: hypothesis, why-chain, switch,
  pre-registered expectation and reading, pool size, status.
* instrument stamping -- a content hash of the code/judge/tool files and
  the relevant environment, taken once, that refuses cross-stamp compares.
* verdict and landing   -- a noise-gated verdict that reads the effect
  against package D's statistics (``wilson_interval``, ``compare_runs``)
  with mandatory blocker classification, and a state machine whose only
  path to ``approved`` is a human approval record.

The statistics live in ``delfin.agent.benchmark`` (package D, mainline --
stable).  Every call from this module into ``benchmark`` goes through
``_noise_gate`` so that if the consumed interface ever changes at
integration there is exactly one place to rewire.
"""

from __future__ import annotations

import hashlib
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence


__all__ = [
    "Experiment",
    "ExperimentError",
    "InstrumentStamp",
    "Measurement",
    "assert_same_stamp",
    "check_reach",
    "instrument_stamp",
    "pre_register",
    "record_measurement",
    "require_switch_off_identical",
    "status_of",
]


class ExperimentError(Exception):
    """A refusal with a reason.  The message is the reason -- it says what
    is missing or wrong and what the valid state is, never bare a code."""


# Allowed pool sizes, exactly as the requirement names them.
_POOL_SIZES = ("small", "large")

# States of the experiment life-cycle, in order.  Phase 2 fixes the first
# three; later phases add the verdict/landing states on top of these.
_STATUSES = ("draft", "registered", "recording", "measured",
             "verdict", "human_review", "approved", "landed")

# States in which a new measurement may be recorded.
_RECORDING_STATES = ("registered", "recording")


@dataclass
class Measurement:
    """One observed outcome of an experiment on one case.  The ``outcome``
    dict is opaque to this module: the caller decides its keys (wall-clock,
    token counts, success flag, ...)."""

    case: str
    outcome: dict[str, Any] = field(default_factory=dict)


@dataclass
class Experiment:
    """The complete written-down definition of one experiment.

    ``why_chain`` is the reason chain from the mechanism to the observed
    outcome (each element one step).  ``switch`` is the single knob behind
    which the candidate change sits; it defaults off and the additivity
    check later proves that switch-off output is byte-identical to baseline.
    ``expectation`` and ``reading`` are the pre-registered expected effect
    and how each outcome is read -- set only by :func:`pre_register`.

    ``status`` is the state-machine field; its legal values are
    ``_STATUSES`` and it only moves forward.
    """

    id: str
    hypothesis: str
    why_chain: list[str] = field(default_factory=list)
    switch: str = ""
    expectation: Optional[str] = None
    reading: Optional[str] = None
    pool_size: str = "small"
    status: str = "draft"
    measurements: list[Measurement] = field(default_factory=list)

    def _move(self, status: str) -> None:
        self.status = status


def _require(condition: bool, message: str) -> None:
    """Raise :class:`ExperimentError` with ``message`` unless *condition*.

    One refusal channel for every guard in this module -- a helper so the
    user-facing reason is one string, not an exception type per case.
    """
    if not condition:
        raise ExperimentError(message)


def status_of(exp: Experiment) -> str:
    """The experiment's current state.  Read-only shorthand for callers."""
    return exp.status


def _validate_draft(exp: Experiment) -> None:
    """The invariants a draft must satisfy before it may be registered."""
    _require(bool(exp.id and exp.id.strip()), "experiment needs a non-empty id")
    _require(bool(exp.hypothesis and exp.hypothesis.strip()),
             "experiment needs a hypothesis")
    _require(bool(exp.switch and exp.switch.strip()),
             "experiment needs a switch behind which the change sits")
    _require(bool(exp.why_chain),
             "experiment needs a why-chain: the reason steps to the effect")
    _require(all(bool(step.strip()) for step in exp.why_chain),
             "why-chain steps must be non-empty")


def pre_register(
    exp: Experiment,
    *,
    expectation: str,
    reading: str,
    pool_size: str,
) -> Experiment:
    """Pre-register *exp*: write down the expected effect and how each
    outcome will be read BEFORE any measurement happens.

    Refuses (raises :class:`ExperimentError`) when the draft is incomplete
    -- no hypothesis, no why-chain, no switch, no expectation, no reading, or
    an unknown pool size -- or when the experiment has already been
    registered or started measuring.  On success the status moves
    ``draft -> registered`` and the expectation/reading are stored.
    """
    _require(exp.status == "draft", f"cannot pre-register an experiment in state {exp.status!r}")
    _validate_draft(exp)
    _require(bool(expectation.strip()), "pre-registration needs an expected effect")
    _require(bool(reading.strip()), "pre-registration needs a reading (how each outcome is judged)")
    _require(pool_size in _POOL_SIZES,
             f"pool_size must be one of {_POOL_SIZES!r}, got {pool_size!r}")
    exp.expectation = expectation
    exp.reading = reading
    exp.pool_size = pool_size
    exp._move("registered")
    return exp


def record_measurement(
    exp: Experiment,
    *,
    case: str,
    outcome: dict[str, Any],
) -> Measurement:
    """Record one observed outcome on *case*.

    Refuses an experiment that has NOT been pre-registered and refuses an
    experiment in a state where measuring is no longer allowed.  This is the
    gate the package exists for: a measurement taken without pre-registration
    is exactly the self-bias it stops, so it raises instead of recording.
    """
    _require(bool(case.strip()), "a measurement needs a case name")
    _require(exp.status in _RECORDING_STATES,
             f"cannot measure an experiment in state {exp.status!r}; "
             "pre-register it (expectation + reading written down) first")
    rec = Measurement(case=case, outcome=dict(outcome))
    exp.measurements.append(rec)
    if exp.status == "registered":
        exp._move("recording")
    return rec


# ---------------------------------------------------------------------------
# Phase 3 -- instrument stamp
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class InstrumentStamp:
    """A content stamp of the instrument a measurement was taken with.

    Hashes the DECLARED code/judge/tool files' contents, the DECLARED
    environment variables' values, and the judge version, all taken once at
    stamp time.  Frozen by construction: an instrument cannot be quietly
    changed under a measurement, and two results measured under different
    stamps must not be compared (see :func:`assert_same_stamp`).

    ``stamp_id`` is the single hash over every component -- content, env,
    judge -- so equality of ids is the equality of the whole instrument.
    """

    content_hash: str
    env_hash: str
    judge: str
    files: tuple[str, ...]
    env_keys: tuple[str, ...]
    taken_at: float

    @property
    def stamp_id(self) -> str:
        h = hashlib.sha256()
        h.update(self.content_hash.encode("utf-8", "replace"))
        h.update(self.env_hash.encode("utf-8", "replace"))
        h.update(str(self.judge).encode("utf-8", "replace"))
        return h.hexdigest()[:16]


def _file_content_hash(files: Sequence[str]) -> str:
    """sha256 over the concatenation of every file's bytes, in declared order.

    Content is read exactly as stored -- an empty file is "empty content",
    not a skip -- so a changed file always changes the hash.
    """
    _require(bool(files), "instrument stamp needs at least one code/judge/tool file")
    seen: set[str] = set()
    h = hashlib.sha256()
    for raw in files:
        path = Path(raw)
        _require(bool(str(path).strip()), "instrument stamp file path must be non-empty")
        key = os.path.abspath(str(path))
        _require(path.is_file(), f"instrument stamp cannot read a missing file: {key}")
        try:
            data = path.read_bytes()
        except OSError as exc:
            raise ExperimentError(
                f"instrument stamp cannot read {key}: {exc}") from exc
        h.update(len(data).to_bytes(8, "big"))
        h.update(data)
        h.update(b"\x00")
        seen.add(key)
    _require(len(seen) == len(files),
             "instrument stamp refuses the same file listed twice (ambiguous instrument)")
    return h.hexdigest()


def _env_hash(env_keys: Sequence[str]) -> str:
    """sha256 over the values of the DECLARED environment variables.

    Only the declared keys are read -- an undeclared environment variable is
    out of scope BY CONSTRUCTION and cannot change the stamp.  A declared key
    that is missing from the environment is a refusal: the instrument is not
    reproducible if a component it names does not exist.
    """
    _require(bool(env_keys), "instrument stamp needs at least one relevant environment key")
    h = hashlib.sha256()
    for key in env_keys:
        _require(bool(str(key).strip()), "environment keys must be non-empty")
        _require(key in os.environ,
                 f"instrument stamp declares environment variable {key!r} which is not set")
        h.update(str(key).encode("utf-8", "replace"))
        h.update(b"\x00")
        h.update(os.environ[key].encode("utf-8", "replace"))
    return h.hexdigest()


def instrument_stamp(
    *,
    files: Sequence[str],
    env_keys: Sequence[str],
    judge: str = "",
) -> InstrumentStamp:
    """Take the instrument stamp for *files* + *env_keys* + *judge*.

    Raises :class:`ExperimentError` on an unusable instrument: no files, a
    missing or unreadable file, a file listed twice, no environment keys, a
    declared environment variable that is not set.  The result is frozen and
    stamped with its creation time, so it can be stored as the record of
    exactly which instrument a comparison was measured under.
    """
    content_hash = _file_content_hash(files)
    env = _env_hash(env_keys)
    return InstrumentStamp(
        content_hash=content_hash,
        env_hash=env,
        judge=str(judge),
        files=tuple(files),
        env_keys=tuple(env_keys),
        taken_at=time.time(),
    )


def assert_same_stamp(a: InstrumentStamp, b: InstrumentStamp) -> None:
    """Refuse to compare two results measured under different instruments.

    Two results are comparable only when BOTH the content and the
    environment and the judge match.  A difference in any one is a refusal
    with a message naming the differing component, so an agent that moved a
    knob between the arms of an experiment sees exactly what changed.
    """
    _require(a.content_hash == b.content_hash,
             "refused: results were measured under different instrument content "
             "(the code, judge or tool files differ)")
    _require(a.env_hash == b.env_hash,
             "refused: results were measured under a different relevant environment "
             "(a declared environment variable changed)")
    _require(str(a.judge) == str(b.judge),
             f"refused: results were judged by different judge versions "
             f"({a.judge!r} vs {b.judge!r})")


# ---------------------------------------------------------------------------
# Phase 4 -- reach and additivity
# ---------------------------------------------------------------------------


def check_reach(
    *,
    switch_read: bool,
    switch_name: str,
    targeted: Sequence[str],
    effective: Sequence[str],
) -> None:
    """Prove the switch is READ and is EFFECTIVE on every case it targets.

    "Ran" is not "hit": before a measurement means anything, the switch must
    have been read (the change was actually in force) and must have changed
    something on exactly the cases it targets.

    Refuses (raises :class:`ExperimentError`) when:
    * the switch was never read -- the experiment would measure a no-op;
    * the switch changed zero targeted cases -- it did nothing at all;
    * the switch hit only a SUBSET of its targets -- that is a partial reach,
      the silent self-deception this package exists to stop;
    * the switch changed a case it did NOT target -- it reaches beyond its
      declared surface, which is a defect, not a hit.
    """
    _require(bool(switch_name.strip()), "reach check needs a switch name")
    _require(bool(targeted), "reach check needs the targeted cases")
    _require(switch_read,
             f"reach check: switch {switch_name!r} was never read; the "
             "'ran' is not a 'hit' -- measure nothing")
    _require(bool(effective),
             f"reach check: switch {switch_name!r} changed nothing on its "
             "targeted cases; measure nothing")
    targeted_set = set(targeted)
    effective_set = set(effective)
    missing = targeted_set - effective_set
    _require(not missing,
             f"reach check: switch {switch_name!r} hit only a subset of its "
             f"targets; missing {sorted(missing)!r} -- a partial reach is not a hit")
    strays = effective_set - targeted_set
    _require(not strays,
             f"reach check: switch {switch_name!r} changed untargeted cases "
             f"{sorted(strays)!r}; it reaches beyond its declared surface")


def require_switch_off_identical(
    *,
    switch_off: bytes,
    baseline: bytes,
    label: str,
) -> None:
    """ABORT unless the switch-off output is byte-identical to baseline.

    This is the additivity check: a candidate change behind one switch,
    default off, must produce exactly baseline when the switch is off.
    Byte-identity is STRICT -- no stripping, no normalising, no smoothing of
    timestamps or counters.  A switch-off output that differs by a single
    byte ABORTS (raises :class:`ExperimentError`) and the caller must not
    record the measurement: the change is not additive by construction.
    """
    _require(bool(label.strip()), "additivity check needs a label for the case")
    _require(
        isinstance(switch_off, bytes) and isinstance(baseline, bytes),
        "additivity check compares byte content only")
    _require(
        switch_off == baseline,
        f"additivity abort: switch-off output for {label!r} is NOT "
        "byte-identical to baseline; the change is not additive by "
        "construction -- abort, do not normalise")


