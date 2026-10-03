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

from dataclasses import dataclass, field
from typing import Any, Optional


__all__ = [
    "Experiment",
    "ExperimentError",
    "pre_register",
    "record_measurement",
    "status_of",
    "Measurement",
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
