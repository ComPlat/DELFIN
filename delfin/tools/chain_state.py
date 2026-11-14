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
import re
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
    #: The step's parameters exactly as Pipeline passes them to run_step,
    #: so a resume can rebuild a REAL run (L1: charge/method/basis/mult,
    #: ...) and the inputs_hash has inputs to hash (L3).
    params: Dict[str, Any] = field(default_factory=dict)
    attempts: List[Dict[str, Any]] = field(default_factory=list)

    # -- transitions -------------------------------------------------------

    def _require(self, state: str, msg: str) -> None:
        if self.state != state:
            raise ChainStateError(
                f"step {self.name!r}: {msg} is not allowed from state "
                f"{self.state!r}"
            )

    def mark_submitted(self, job_id: str, inputs_hash: Optional[str] = None,
                       attempt_no: Optional[int] = None) -> None:
        self._require("pending", "mark_submitted")
        self.state = "submitted"
        self.job_id = job_id
        # The attempt number comes from whoever owns the attempt lifecycle:
        # a fresh step counts up on its first submission (attempt 1); a
        # retry's number was already advanced by reset_for_retry (L5/L6 —
        # double-counting mislabeled the ledger's entry as attempt 3).
        if attempt_no is not None:
            self.attempt = attempt_no
        elif self.attempt == 0:
            self.attempt = 1
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
        falsy = [k for k, v in evidence.items()
                 if (k in REQUIRED_EVIDENCE or k == "tool_reported_success")
                 and not v]
        if falsy:
            raise ChainStateError(
                f"step {self.name!r}: finished evidence {falsy} is falsy — "
                "each piece of evidence must be a truthy proof, not the "
                "recorded absence of one"
            )
        has_orca_proof = all(k in evidence for k in REQUIRED_EVIDENCE)
        # A non-ORCA step (xtb, crest, ...) may be finished on the word of
        # the run that observed its success; its output is re-verified on
        # resume (verify_finished_step).  Partial ORCA keys alone are a
        # guess for any tool (JC1 target 1 / finding F-C).
        has_run_report = bool(evidence.get("tool_reported_success"))
        if not (has_orca_proof or has_run_report):
            raise ChainStateError(
                f"step {self.name!r}: finished evidence missing "
                f"{[k for k in REQUIRED_EVIDENCE if k not in evidence]} — "
                "evidence comes from delfin.orca_read (normal termination, "
                "converged SCF) or from the run that observed the step, "
                "never a guess"
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

    def reset_for_retry(self, attempt_dir: Optional[str] = None,
                        end: Optional[str] = None) -> None:
        """Return a ``failed`` step to ``pending`` for a bounded resubmission.

        The OLD attempt's work dir is kept untouched (recalc path: never
        overwrite a calculation) and recorded in ``attempts[]`` with its
        end reason — L5: the killed attempt must remain visible in the
        state, not be overwritten by the retry.
        """
        self._require("failed", "reset_for_retry")
        self.attempts.append({
            "attempt": self.attempt,
            # L7: only JSON types go into the step — the work_dir on the
            # live path is a pathlib.Path (set from the attempt-dir
            # chain), which json.dump refused and the successful retry's
            # recording crashed on save.
            "work_dir": str(self.work_dir) if self.work_dir is not None else None,
            "end": end or self.error or "failed",
        })
        # L5: the retry is a NEW attempt — the counter must move, else the
        # attempt budget and the ledger disagree with the real runs.
        self.attempt += 1
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
            "params": self.params,
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
            params=dict(data.get("params", {})),
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
            if sd.get("state") == "finished":
                # Torn/unfinished check: a stored `finished` step must
                # carry a full evidence dict — either the ORCA markers or
                # the run's own success report (F-B/F-C); anything else is
                # a half-written state and is downgraded to pending.
                ev = sd.get("evidence")
                ok_evidence = isinstance(ev, dict) and (
                    all(k in ev for k in REQUIRED_EVIDENCE)
                    or bool(ev.get("tool_reported_success"))
                )
                if not ok_evidence:
                    logger.warning(
                        "chain %s step %s: torn/unfinished state downgraded "
                        "to pending",
                        data["name"],
                        sd.get("name"),
                    )
                    sd = dict(sd, state="pending", evidence=None)
            chain.steps.append(ChainStep.from_dict(sd))
        return chain

    # -- steps -------------------------------------------------------------

    def add_step(
        self,
        name: str,
        tool: str,
        inputs_hash: str = "",
        params: Optional[Dict[str, Any]] = None,
    ) -> ChainStep:
        if params is None:
            params = {}
        # L3: the inputs hash comes from the step's actual inputs, never ""
        if not inputs_hash and params:
            inputs_hash = compute_inputs_hash(**params)
        step = ChainStep(
            name=name, tool=tool, inputs_hash=inputs_hash,
            params=dict(params),
        )
        self.steps.append(step)
        # Persist immediately: a chain step must never be lost to a crash
        # between add_step and the next explicit save (the state file is
        # written by create() BEFORE the first add_step, so without this
        # the on-disk state would stay empty and a resume would see a
        # zero-step chain).
        self.save()
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
                # L7: convert centrally — any non-JSON value anywhere in
                # the payload (a stray Path in evidence/params/attempts)
                # is str()'d here, so a call site can never crash the
                # save.  JSON-safe values pass through untouched.
                json.dump(data, fh, indent=2, default=str)
                fh.flush()
                os.fsync(fh.fileno())
            os.replace(tmp_path, self.state_path)
        except BaseException:
            tmp_path.unlink(missing_ok=True)
            raise


# ---------------------------------------------------------------------------
# Resume decides from facts (JC1 target 2)
# ---------------------------------------------------------------------------


@dataclass
class FinishVerdict:
    """What re-verification of a step marked ``finished`` concluded."""

    #: The output still carries the proof the step ended well.
    verified: bool
    #: Where the step went when not verified (``pending``).
    downgraded_to: Optional[str] = None
    reason: str = ""


@dataclass
class RunningVerdict:
    """What a scheduler-backed classification of a ``running`` step said."""

    #: "wait" | "classify" | "mark_failed"
    action: str
    error_kind: Optional[str] = None
    reason: str = ""


# -- JSON type safety -----------------------------------------------------

def json_safe(value: Any) -> Any:
    """Convert a value to JSON-safe types (recursively).

    pathlib.Path -> str; dicts/lists recurse; anything JSON-safe passes
    through.  Used when the chain state is written so no call site can
    put a non-JSON type into the persisted step.
    """
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def _read_step_evidence(step: ChainStep, chain_work_dir=None) -> Dict[str, Any]:
    """Re-read a step's own output and return the evidence it carries.

    Uses the same readers the run-time path uses: normal termination via
    ``delfin.orca_read._terminated_normally_text`` (the ``ORCA TERMINATED
    NORMALLY`` marker near the file end) and SCF convergence via
    ``_scf_converged_text``.  Never a guess — the absence of a file or a
    marker is reported, not assumed away.
    """
    from delfin.orca_read import _scf_converged_text, _terminated_normally_text

    out_path = _step_output_path(step, chain_work_dir)
    if out_path is None or not out_path.is_file():
        return {}
    evidence: Dict[str, Any] = {}
    evidence["normal_termination"] = _terminated_normally_text(out_path)
    try:
        text = out_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return {}
    evidence["scf_converged"] = _scf_converged_text(text)
    return evidence


def _step_output_path(step: ChainStep, chain_work_dir=None) -> Optional[Path]:
    """The step's ``<base>.out`` inside its recorded work dir.

    Falls back to the pipeline's own per-step directory layout
    ``<NN>_<step-name>`` under the chain's work dir (pipeline.py writes step
    ``i`` into ``base / f"{i:02d}_{label}"``) when the step carries no
    work_dir of its own.
    """
    candidates: list = []
    if step.work_dir:
        wd = Path(step.work_dir)
        if wd.is_dir():
            candidates = sorted(wd.glob("*.out"))
    if not candidates and chain_work_dir is not None:
        base = Path(chain_work_dir)
        subs = [base / step.name]
        subs += sorted(
            p for p in base.glob(f"*_{step.name}")
            if p.is_dir()
        )
        for sub in subs:
            if sub.is_dir():
                candidates = sorted(sub.glob("*.out"))
                if candidates:
                    break
    return candidates[-1] if candidates else None


def _is_orca_tool(tool: Optional[str]) -> bool:
    """ORCA-family steps are verified by the ORCA markers; every other
    tool (xtb, crest, ...) by its own non-empty output (F-B: demanding
    ORCA markers of an xtb step would rerun finished work forever)."""
    return bool(tool) and str(tool).lower().startswith("orca")


def verify_finished_step(
    step: ChainStep, chain_work_dir: Optional[str | Path] = None
) -> FinishVerdict:
    """Re-verify a step marked ``finished`` from its output before resume
    skips it (JC1 target 2).  On any missing proof the step is DOWNGRADED
    to ``pending`` in place — a half-written or vanished step is never
    reported finished.

    The proof is the step's OWN tool's success evidence: ORCA-family
    steps show ``ORCA TERMINATED NORMALLY`` and a converged SCF (read
    through :mod:`delfin.orca_read`); other tools (xtb, crest, ...) end
    with their own output file, which must still exist and be non-empty.
    """
    if step.state != "finished":
        return FinishVerdict(verified=False, downgraded_to=None,
                             reason=f"step {step.name!r} is not finished")
    if not _is_orca_tool(step.tool):
        out_path = _step_output_path(step, chain_work_dir)
        if out_path is None or not out_path.is_file():
            step.state = "pending"
            step.evidence = None
            return FinishVerdict(
                verified=False,
                downgraded_to="pending",
                reason=f"the step's output ({step.tool}) is gone",
            )
        try:
            if out_path.stat().st_size == 0:
                step.state = "pending"
                step.evidence = None
                return FinishVerdict(
                    verified=False,
                    downgraded_to="pending",
                    reason=f"the step's output ({step.tool}) is empty",
                )
            stored_hash = (step.evidence or {}).get("output_sha256")
            if stored_hash:
                actual = hashlib.sha256(
                    out_path.read_bytes()
                ).hexdigest()
                if actual != stored_hash:
                    step.state = "pending"
                    step.evidence = None
                    return FinishVerdict(
                        verified=False,
                        downgraded_to="pending",
                        reason="the step's output changed since the step "
                        "finished (crash mid-write)",
                    )
        except OSError as exc:
            step.state = "pending"
            step.evidence = None
            return FinishVerdict(
                verified=False,
                downgraded_to="pending",
                reason=f"the step's output is unreadable: {exc}",
            )
        return FinishVerdict(verified=True)
    evidence = _read_step_evidence(step, chain_work_dir)
    if not evidence.get("normal_termination"):
        step.state = "pending"
        step.evidence = None
        return FinishVerdict(
            verified=False,
            downgraded_to="pending",
            reason="no ORCA TERMINATED NORMALLY marker in the step's output",
        )
    if not evidence.get("scf_converged"):
        step.state = "pending"
        step.evidence = None
        return FinishVerdict(
            verified=False,
            downgraded_to="pending",
            reason="output does not show SCF convergence",
        )
    return FinishVerdict(verified=True)


def _job_state(job_id: str, job_state) -> object:
    if job_state is not None:
        return job_state(job_id)
    from delfin.scheduler_client import query

    # squeue -j <id>: non-empty stdout means the job is queued or running.
    out = query("squeue", ["-j", str(job_id), "-h"])
    if out is delfin_scheduler_THROTTLED():
        return out
    if out is None:
        return None
    return str(out).strip() or None


def delfin_scheduler_THROTTLED():
    from delfin.scheduler_client import THROTTLED

    return THROTTLED


def classify_running_step(
    step: ChainStep, job_state=None, chain_work_dir=None
) -> RunningVerdict:
    """Check a step marked ``running`` against the scheduler (JC1 target 2).

    Still queued/running -> wait.  Gone from the scheduler -> the end is
    classified from the step's own output through
    ``delfin.tools._recovery.classify_error`` (never a guess): a job whose
    output shows no normal termination failed.  A throttled scheduler
    answer is NOT a fact — wait and re-ask later.
    """
    if step.state != "running":
        return RunningVerdict(action="wait", reason=f"step {step.name!r} is not running")
    state = _job_state(step.job_id, job_state)
    if state is delfin_scheduler_THROTTLED():
        return RunningVerdict(
            action="wait",
            reason="scheduler answer throttled — not a fact, re-ask later",
        )
    if state is not None:
        return RunningVerdict(
            action="wait", reason=f"job {step.job_id} is {state} in the scheduler"
        )
    # Job gone: classify the end from the step's own output.
    out_path = _step_output_path(step, chain_work_dir)
    evidence = _read_step_evidence(step, chain_work_dir)
    from delfin.tools._recovery import classify_error

    if evidence.get("normal_termination") and evidence.get("scf_converged"):
        step.mark_finished(
            evidence={
                "normal_termination": True,
                "scf_converged": True,
            }
        )
        return RunningVerdict(
            action="mark_finished",
            reason="job gone but the output proves normal termination + "
            "converged SCF",
        )
    if out_path is None:
        # Nothing to read: the classification rests on the job's absence
        # alone (recovery vocabulary, not a guess about the output).
        reason = (
            f"job {step.job_id} no longer in the scheduler and no step "
            "output exists"
        )
        kind = classify_error(reason)
        step.mark_failed(reason, error_kind=kind.value)
        return RunningVerdict(
            action="classify", error_kind=kind.value, reason=reason
        )
    reason = step.error or (
        f"job {step.job_id} no longer in the scheduler and the output shows "
        "no normal termination"
    )
    kind = classify_error(reason)
    step.mark_failed(reason, error_kind=kind.value)
    return RunningVerdict(action="mark_failed", error_kind=kind.value, reason=reason)


def build_resume_plan(
    chain: "ChainState",
    job_state=None,
    work_dir: Optional[str | Path] = None,
) -> Dict[str, Any]:
    """Produce the resume decision from FACTS (JC1 target 2).

    - finished steps are re-verified from their output; only verified ones
      are skipped;
    - the first running step is checked against the scheduler: live -> the
      plan says wait; gone -> classified to finished (with proof) or failed;
    - the next step is the first one that still needs to run.
    """
    steps_to_skip = []
    for step in chain.steps:
        if step.state == "finished":
            verdict = verify_finished_step(step, chain_work_dir=work_dir)
            if verdict.verified:
                steps_to_skip.append(step)
                continue
            # Downgraded to pending: it will run again.
            return _plan_after(chain, chain.steps.index(step), steps_to_skip, wait=False)
        if step.state == "running":
            verdict = classify_running_step(
                step, job_state=job_state, chain_work_dir=work_dir
            )
            if verdict.action == "wait":
                return {
                    "steps_to_skip": steps_to_skip,
                    "next_step": None,
                    "wait": True,
                    "reason": verdict.reason,
                }
            if verdict.action == "classify":
                # Job vanished without output: this step must be
                # resubmitted, it does not unblock the chain.
                return _plan_after(
                    chain, chain.steps.index(step), steps_to_skip, wait=False
                )
            # mark_finished or mark_failed: finished-with-proof is skipped
            # by the caller (record it here), failed becomes the next step.
            if verdict.action == "mark_finished":
                steps_to_skip.append(step)
                continue
            # mark_failed -> the chain must resubmit this step (phase 3).
            return _plan_after(
                chain, chain.steps.index(step), steps_to_skip, wait=False
            )
        if step.state == "failed":
            # A failed step CAN be the resume target (target 3): it is
            # resubmitted in a new attempt dir — but only if the recovery
            # policy allows it; a non-recoverable cause blocks the chain.
            from delfin.tools import chain_state as _cs

            decision = _cs.resubmission_decision(step)
            if decision.action == "stop":
                return {
                    "steps_to_skip": steps_to_skip,
                    "next_step": None,
                    "wait": False,
                    "blocked": True,
                    "reason": decision.reason,
                }
            return _plan_after(chain, chain.steps.index(step), steps_to_skip, wait=False)
        if step.state == "submitted":
            # JC1 F-A: a submitted step whose job has vanished from the
            # scheduler must be RESUBMITTED, not walked past — it never
            # unblocks the chain.  The step's own output decides whether
            # the job actually finished before disappearing (proof ->
            # finished) or not (mark_failed -> resubmit via the plan).
            from delfin.tools import chain_state as _cs

            state = _cs._job_state(step.job_id, job_state)
            if state is _cs.delfin_scheduler_THROTTLED():
                return {
                    "steps_to_skip": steps_to_skip,
                    "next_step": None,
                    "wait": True,
                    "reason": "scheduler answer throttled — not a fact, "
                              "re-ask later",
                }
            if state is not None:
                return {
                    "steps_to_skip": steps_to_skip,
                    "next_step": None,
                    "wait": True,
                    "reason": f"job {step.job_id} is {state} in the scheduler",
                }
            # Job gone: the step's own output decides the end (target 2:
            # never a guess from the state file alone).  The proof is the
            # step's OWN tool's success evidence (same rule as
            # verify_finished_step): ORCA-family steps need the ORCA
            # markers, other tools their own non-empty output.
            out_path = _cs._step_output_path(step, work_dir)
            evidence = _cs._read_step_evidence(step, work_dir) if out_path else {}
            if _cs._is_orca_tool(step.tool):
                finished = (
                    evidence.get("normal_termination")
                    and evidence.get("scf_converged")
                )
            else:
                finished = bool(out_path) and out_path.is_file() \
                    and out_path.stat().st_size > 0
            if finished:
                step.mark_running()
                step.mark_finished(evidence={
                    "tool_reported_success": True,
                    **({"normal_termination": True,
                        "scf_converged": True} if _cs._is_orca_tool(step.tool)
                       else {}),
                })
                steps_to_skip.append(step)
                continue
            step.mark_failed(
                "job vanished from the scheduler without completing "
                f"(job {step.job_id})"
            )
            return _plan_after(chain, chain.steps.index(step), steps_to_skip,
                               wait=False)
        if step.state == "pending":
            return _plan_after(chain, chain.steps.index(step), steps_to_skip, wait=False)
    return {
        "steps_to_skip": steps_to_skip,
        "next_step": None,
        "wait": False,
        "complete": True,
        "reason": "every step finished with verified evidence",
    }


def _plan_after(chain, idx, steps_to_skip, wait):
    remaining = chain.steps[idx:]
    return {
        "steps_to_skip": steps_to_skip,
        "next_step": remaining[0] if remaining else None,
        "wait": wait,
    }


# ---------------------------------------------------------------------------
# Bounded resubmission (JC1 target 3)
# ---------------------------------------------------------------------------


@dataclass
class ResubmissionDecision:
    """What to do with a failed step: retry in a NEW attempt dir, or stop."""

    #: "retry" | "stop"
    action: str
    reason: str = ""
    #: Empty by design: PAL, maxcore and time are NEVER changed to make a
    #: job start sooner (JC1 target 3); only _recovery's convergence kwargs
    #: (scf_maxiter / scf_extra) may come through here.
    kwargs_overrides: Dict[str, Any] = field(default_factory=dict)


def next_attempt_dir(step: ChainStep, chain_work_dir: str | Path) -> Path:
    """A NEW attempt directory for a failed step, keeping the old one.

    The recalc path never overwrites a calculation (delfin/cli_recalc.py,
    recalc_control.record_completed_run): the resubmission goes into
    ``attempt_<N+1>`` next to the old attempt dir.
    """
    base = Path(chain_work_dir) / step.name
    used = sorted(
        int(p.name.rsplit("_", 1)[-1])
        for p in base.glob("attempt_*")
        if p.name.rsplit("_", 1)[-1].isdigit()
    )
    if step.work_dir:
        old = Path(step.work_dir)
        name = old.name
        # The first attempt's dir is the step's numbered dir NN_<name>
        # (e.g. 01_orca_opt): the number LEADS, it IS attempt 1, and it
        # must be counted — else the second attempt is mislabeled
        # attempt_1 (live run #3).  A previous retry dir attempt_<n> is
        # counted through its trailing number.
        m = re.match(r"^(\d{2})_(.+)$", name)
        if m and m.group(2) == step.name:
            used.append(1)
        elif (tail := name.rsplit("_", 1)[-1]).isdigit() \
                and old.parent.name == step.name:
            used.append(int(tail))
    n = (used[-1] + 1) if used else 1
    return base / f"attempt_{n}"


def resubmission_decision(
    step: ChainStep, max_attempts: int = 3
) -> ResubmissionDecision:
    """Decide whether a failed step may be resubmitted, and how often.

    The decision comes from the recovery policy, not from hope: a cause
    ``delfin.tools._recovery.suggest_recovery`` calls non-recoverable stops
    the chain with its reason; a recoverable cause is retried at most
    ``max_attempts`` times, then the chain stops with the reason.
    """
    from delfin.tools._recovery import suggest_recovery

    suggestion = suggest_recovery(
        _error_kind_from_value(step.error_kind), kwargs={}
    )
    if suggestion.action == "give_up":
        return ResubmissionDecision(
            action="stop",
            reason=f"non-recoverable cause ({step.error_kind}): "
            f"{step.error or suggestion.reason}",
        )
    if step.attempt >= max_attempts:
        return ResubmissionDecision(
            action="stop",
            reason=f"attempt budget exhausted ({step.attempt}/{max_attempts} "
            f"attempts): {step.error or 'unknown error'}",
        )
    return ResubmissionDecision(
        action="retry",
        reason=f"recoverable cause ({step.error_kind}): {step.error}",
        kwargs_overrides=dict(suggestion.kwargs_overrides),
    )


def _error_kind_from_value(value: Optional[str]):
    """Back to the ErrorKind enum member for a stored value string."""
    if value is None:
        from delfin.tools._recovery import ErrorKind

        return ErrorKind.NONE
    from delfin.tools._recovery import ErrorKind

    for member in ErrorKind:
        if member.value == value:
            return member
    return ErrorKind.TOOL_FAILED
