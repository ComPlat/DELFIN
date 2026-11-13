"""Crash-safe chain state for calculation chains (JC1 target 1).

A calculation chain (e.g. GOAT -> relaxed scan -> DLPNO single points) runs
for hours to days across several SLURM jobs.  This module persists the chain
as a state machine so that a crashed agent session, a rebooted node or a
dead SLURM job can be picked up again without rerunning finished steps or
overwriting calculations.

Design rules this module encodes (see tests/test_jc1_red_chain_state.py and
the reviewer's tests/test_jc1_adversarial_checkpoint.py):

- The state file is written ATOMICALLY (write-temp-``os.replace``) next to
  the chain's work dir.  A crash between any two writes leaves the previous
  state readable; a torn file is never reported as a finished step.
- A step may only be marked ``finished`` with real EVIDENCE (normal ORCA
  termination AND converged SCF, as read by :mod:`delfin.orca_read`), never
  by guess.  Invalid state transitions are refused.

The consumer of this file is :mod:`delfin.tools.pipeline` (one mechanism,
not two — JC1 target 4) and :mod:`delfin.cli_chain` (the ``delfin chain``
status/resume entry point).
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

#: Per-step states.  ``finished`` requires evidence; ``failed`` requires a
#: reason.  ``pending -> submitted -> running -> finished | failed``.
STATES = ("pending", "submitted", "running", "finished", "failed")

#: Allowed transitions.  Anything else is refused with ChainStateError.
_TRANSITIONS: Dict[str, tuple] = {
    "pending": ("submitted",),
    "submitted": ("running", "failed"),
    "running": ("finished", "failed"),
    "finished": (),
    "failed": (),
}

#: Evidence keys a ``finished`` step must carry.  They come from
#: delfin.orca_read's normal-termination / SCF-convergence readers — the
#: only accepted proof that a step ended well.
REQUIRED_EVIDENCE = ("normal_termination", "scf_converged")


class ChainStateError(Exception):
    """An invalid chain-state operation (bad transition, missing evidence)."""


@dataclass
class ChainStep:
    """One step of a chain, persisted with the facts that ended it."""

    name: str
    tool: str
    state: str = "pending"
    job_id: Optional[str] = None
    work_dir: Optional[str] = None
    attempt: int = 0
    error: Optional[str] = None
    #: Evidence dict for ``finished`` (normal_termination, scf_converged).
    evidence: Optional[Dict[str, Any]] = None
    #: Hash of the step's inputs (geometry + parameters) at submission time,
    #: so a resume can detect that a step's inputs changed since it ran.
    inputs_hash: str = ""
    #: Recoverable-error classification from delfin.tools._recovery.
    error_kind: Optional[str] = None
    attempts: List[Dict[str, Any]] = field(default_factory=list)

    # -- transitions -------------------------------------------------------

    def _require(self, state: str, msg: str) -> None:
        if self.state != state:
            raise ChainStateError(
                f"step {self.name!r}: {msg} is not allowed from state "
                f"{self.state!r}"
            )

    def mark_submitted(self, job_id: str, inputs_hash: Optional[str] = None) -> None:
        self._require("pending", "mark_submitted")
        self.state = "submitted"
        self.job_id = job_id
        self.attempt += 1
        if inputs_hash is not None:
            self.inputs_hash = inputs_hash

    def mark_running(self) -> None:
        self._require("submitted", "mark_running")
        self.state = "running"

    def mark_finished(self, evidence: Optional[Dict[str, Any]]) -> None:
        self._require("running", "mark_finished")
        if not evidence or not isinstance(evidence, dict):
            raise ChainStateError(
                f"step {self.name!r}: finished requires evidence, not a guess"
            )
        missing = [k for k in REQUIRED_EVIDENCE if k not in evidence]
        if missing:
            raise ChainStateError(
                f"step {self.name!r}: finished evidence missing {missing} — "
                "evidence comes from delfin.orca_read (normal termination, "
                "converged SCF), never a guess"
            )
        self.state = "finished"
        self.evidence = dict(evidence)

    def mark_failed(self, error: str, error_kind: Optional[str] = None) -> None:
        if self.state not in ("submitted", "running"):
            raise ChainStateError(
                f"step {self.name!r}: mark_failed is not allowed from state "
                f"{self.state!r}"
            )
        self.state = "failed"
        self.error = error
        self.error_kind = error_kind

    def reset_for_retry(self, attempt_dir: Optional[str] = None) -> None:
        """Return a ``failed`` step to ``pending`` for a bounded resubmission.

        The OLD attempt's work dir is kept untouched (recalc path: never
        overwrite a calculation); the new attempt directory is recorded.
        """
        self._require("failed", "reset_for_retry")
        self.state = "pending"
        self.job_id = None
        self.error = None
        self.error_kind = None
        self.evidence = None
        if attempt_dir is not None:
            self.work_dir = attempt_dir

    # -- persistence -------------------------------------------------------

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "tool": self.tool,
            "state": self.state,
            "job_id": self.job_id,
            "work_dir": self.work_dir,
            "attempt": self.attempt,
            "error": self.error,
            "error_kind": self.error_kind,
            "evidence": self.evidence,
            "inputs_hash": self.inputs_hash,
            "attempts": list(self.attempts),
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ChainStep":
        return cls(
            name=data["name"],
            tool=data["tool"],
            state=data.get("state", "pending"),
            job_id=data.get("job_id"),
            work_dir=data.get("work_dir"),
            attempt=data.get("attempt", 0),
            error=data.get("error"),
            error_kind=data.get("error_kind"),
            evidence=data.get("evidence"),
            inputs_hash=data.get("inputs_hash", ""),
            attempts=list(data.get("attempts", [])),
        )


def compute_inputs_hash(geometry: Optional[str] = None, **params: Any) -> str:
    """Stable hash over a step's inputs, so a resume can detect drift."""
    payload = json.dumps(
        {"geometry": str(geometry) if geometry else None, "params": params},
        sort_keys=True,
        default=str,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


@dataclass
class ChainState:
    """The persisted state machine of one calculation chain."""

    name: str
    work_dir: Path
    steps: List[ChainStep] = field(default_factory=list)
    created_step_count: int = 0

    # -- construction ------------------------------------------------------

    @classmethod
    def create(cls, name: str, work_dir: str | Path) -> "ChainState":
        """Create a new chain; the state file lives NEXT TO the work dir."""
        wd = Path(work_dir)
        wd.mkdir(parents=True, exist_ok=True)
        state = cls(name=name, work_dir=wd)
        state.save()
        return state

    @classmethod
    def load(cls, state_path: str | Path) -> "ChainState":
        """Load a chain state file.

        A torn (half-written) file loads with steps degraded to ``pending``
        — a crash must never let a half-written step be reported finished.
        """
        path = Path(state_path)
        try:
            data = json.loads(path.read_text())
        except json.JSONDecodeError:
            logger.warning(
                "chain state file %s is torn (half-written); loading with "
                "no steps — nothing will be reported finished",
                path,
            )
            data = {"name": "torn-chain", "work_dir": str(path.parent), "steps": []}
        chain = cls(
            name=data["name"],
            work_dir=Path(data["work_dir"]),
        )
        for sd in data.get("steps", []):
            if sd.get("state") == "finished" and (
                not isinstance(sd.get("evidence"), dict)
                or any(k not in sd["evidence"] for k in REQUIRED_EVIDENCE)
            ):
                logger.warning(
                    "chain %s step %s: torn/unfinished state downgraded to "
                    "pending",
                    data["name"],
                    sd.get("name"),
                )
                sd = dict(sd, state="pending", evidence=None)
            chain.steps.append(ChainStep.from_dict(sd))
        return chain

    # -- steps -------------------------------------------------------------

    def add_step(self, name: str, tool: str, inputs_hash: str = "") -> ChainStep:
        step = ChainStep(name=name, tool=tool, inputs_hash=inputs_hash)
        self.steps.append(step)
        return step

    def step(self, name: str) -> ChainStep:
        for s in self.steps:
            if s.name == name:
                return s
        raise ChainStateError(f"chain {self.name!r}: no step named {name!r}")

    # -- persistence (atomic) ----------------------------------------------

    @property
    def state_path(self) -> Path:
        return self.work_dir / "chain_state.json"

    def save(self) -> None:
        """Write the state atomically: temp file in the same dir, then
        ``os.replace``.  A crash before the rename leaves the old state
        intact and readable; the temp file is cleaned up either way."""
        self.work_dir.mkdir(parents=True, exist_ok=True)
        data = {
            "name": self.name,
            "work_dir": str(self.work_dir),
            "steps": [s.to_dict() for s in self.steps],
        }
        fd, tmp_name = tempfile.mkstemp(
            prefix=".chain_state.", suffix=".tmp", dir=str(self.work_dir)
        )
        tmp_path = Path(tmp_name)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2)
                fh.flush()
                os.fsync(fh.fileno())
            os.replace(tmp_path, self.state_path)
        except BaseException:
            tmp_path.unlink(missing_ok=True)
            raise
