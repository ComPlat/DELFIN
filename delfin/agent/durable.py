"""Durable sessions: a worker killed mid-turn comes back to its turn.

This is the decision layer under DELFIN's existing per-turn checkpoint
(``session_store.save_turn_checkpoint`` over ``<sid>.turn.json``) and the
append-only step journal (``delfin.agent.step_journal``). The journal already
records every tool call twice -- a ``start`` before it runs and a ``finish``
after -- and a ``start`` with no matching ``finish`` is an OPEN step that
"may or may not have run". That existing machinery is PRODUCED by the real
tool loop (``api_client`` brackets every call with
``step_journal.auto_start`` / ``auto_finish``).

What this module adds is the decision a restart needs but the note alone
does not give: for each OPEN side-effecting step, is its side effect (the
file, the commit, the job) already present, or must a person decide? The
rule is unconditional -- an OPEN side-effecting step is NEVER run again
automatically. The resumed session first checks whether the step's effect
happened; what a check cannot confirm, the session asks about. A finished
step is never repeated at all.

The guard the supervisor consults is :meth:`RecoveryPlan.would_block_rerun`:
before the tool loop re-dispatches a call after a crash, it asks whether that
call is a step that must not be auto-re-run. The answer is True for every
finished step and every OPEN step that is not a verified read-only step.

Records are resolved PER CALL against ``root`` (default
``session_store._SESSIONS_DIR``) so a test can redirect a store without
touching a real home -- the same rule ``step_journal.journal_path`` already
follows.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

from delfin.agent import session_store, step_journal

#: The possible answer to "did this OPEN step's side effect happen?".
CheckerResult = str  # "present" | "absent" | "unknown"


class RecoveryVerdict(enum.Enum):
    """What a restart may do with a journaled step.

    * ``VERIFIED_DONE`` -- an OPEN step whose side effect is confirmed
      present on disk: do NOT re-run, the outcome is known.
    * ``MUST_ASK`` -- an OPEN step whose side effect could not be confirmed
      (absent or unverifiable): never auto-rerun, a person decides.
    * ``FINISHED`` -- the journal recorded a finish: never repeated.
    * ``RE_RUN_OK`` -- an OPEN read-only step (no side effect): safe to
      re-run, it changes nothing.
    """

    VERIFIED_DONE = "verified_done"
    MUST_ASK = "must_ask"
    FINISHED = "finished"
    RE_RUN_OK = "re_run_ok"


@dataclass(frozen=True)
class StepDecision:
    """One journaled step's recovery status."""

    turn: int
    step: int
    tool: str
    side_effecting: bool
    verdict: RecoveryVerdict
    #: Whether the tool loop may re-dispatch this call after a crash.
    auto_rerun: bool
    #: The raw journal record this decision was derived from.
    record: dict[str, Any]


@dataclass
class RecoveryPlan:
    """What a restart should do with one session's in-flight turn."""

    session: str
    #: The turn number to resume (from the checkpoint when present, else the
    #: highest turn the journal knows). 0 when the session was idle.
    resume_turn: int
    #: True when the session died mid-turn (checkpoint or open step present).
    in_flight: bool
    decisions: list[StepDecision] = field(default_factory=list)

    def decision_for_open_step(self, step: int, turn: int | None = None,
                                ) -> StepDecision | None:
        """The OPEN recovery decision for one step, if it is still open.

        A registered FINISHED step is answered None: it is not open, it must
        never be repeated, and there is no recovery decision about re-running
        it -- the journal already closed it.
        """
        for d in self.decisions:
            if d.step == step and (turn is None or d.turn == turn) \
                    and d.verdict is not RecoveryVerdict.FINISHED:
                return d
        return None

    def would_block_rerun(self, session: str, *, tool: str,
                          step: int, turn: int) -> bool:
        """May the tool loop NOT re-dispatch this call after a crash?

        True for every step that must not run again automatically: any
        finished step, and any OPEN step that is not a read-only step.
        (Read-only OPEN steps may be re-run; they change nothing.)
        """
        if session != self.session:
            return False
        dec = self.decision_for_open_step(step=step, turn=turn)
        if dec is None:
            # Not a journaled step of this turn: normal execution. (A finished
            # step of an earlier turn carries no open decision.)
            return False
        # Only a read-only OPEN step may be re-drawn after a crash. Every
        # side-effecting OPEN step -- verified done OR must-ask -- is blocked
        # from an automatic re-run.
        return dec.verdict is not RecoveryVerdict.RE_RUN_OK


def _step_finished(session: str, *, root: "str | Path | None",
                   step: int, turn: int) -> bool:
    for entry in step_journal.read_steps(session, root=root):
        if entry.get("kind") == "finish" and entry.get("step") == step \
                and entry.get("turn") == turn:
            return True
    return False


def _default_checker(_step: dict[str, Any]) -> CheckerResult:
    """No checker supplied: the effect is unverifiable -> MUST_ASK."""
    return "unknown"


def decide_recovery(session: str, *, root: "str | Path | None" = None,
                    checkers: "dict[str, Callable[[dict[str, Any]], CheckerResult]] | None" = None,
                    ) -> RecoveryPlan:
    """Resolve how a restart resumes one session.

    Reads the REAL checkpoint and journal (the ones the tool loop wrote) and
    turns each OPEN step into a verdict. ``checkers`` maps a tool name (or
    ``"*"``) to a callable that inspects the step record and returns
    ``"present"`` / ``"absent"`` / ``"unknown"`` -- the actual side-effect
    check. Without a checker the verdict is ``MUST_ASK``.
    """
    checkers = checkers or {}
    base = Path(root) if root is not None else Path(session_store._SESSIONS_DIR)

    checkpoint = session_store.load_turn_checkpoint(session)
    cp_turn = checkpoint.get("turn") if checkpoint else None
    open_steps, _ = step_journal.open_steps(session, root=base)

    # Highest turn the journal knows (open or finished), for the resume point.
    journal_turn = 0
    all_steps, _ = step_journal.read_steps(session, root=base)
    for entry in all_steps:
        if not isinstance(entry, dict):
            continue
        t = entry.get("turn")
        if isinstance(t, int):
            journal_turn = max(journal_turn, t)

    resume_turn = cp_turn if cp_turn is not None else (journal_turn or 0)
    in_flight = checkpoint is not None or bool(open_steps)

    decisions: list[StepDecision] = []
    for st in open_steps:
        tool = str(st.get("tool") or "")
        side_effecting = bool(st.get("side_effect", step_journal.side_effecting(tool)))
        checker = checkers.get(tool) or checkers.get("*") or _default_checker
        result = checker(st) if side_effecting else "absent"
        if not side_effecting:
            verdict = RecoveryVerdict.RE_RUN_OK  # read-only: safe to re-run
            auto = True
        elif result == "present":
            verdict = RecoveryVerdict.VERIFIED_DONE
            auto = False
        else:  # "absent" or "unknown"
            verdict = RecoveryVerdict.MUST_ASK
            auto = False
        decisions.append(StepDecision(
            turn=int(st.get("turn") or 0),
            step=int(st.get("step") or 0),
            tool=tool,
            side_effecting=side_effecting,
            verdict=verdict,
            auto_rerun=auto,
            record=st,
        ))

    # Finished steps are never repeated: record them as FINISHED.
    all_steps, _ = step_journal.read_steps(session, root=base)
    for entry in all_steps:
        if not isinstance(entry, dict) or entry.get("kind") != "finish":
            continue
        decisions.append(StepDecision(
            turn=int(entry.get("turn") or 0),
            step=int(entry.get("step") or 0),
            tool=str(entry.get("tool") or ""),
            side_effecting=True,
            verdict=RecoveryVerdict.FINISHED,
            auto_rerun=False,
            record=entry,
        ))

    return RecoveryPlan(session=session, resume_turn=resume_turn,
                        in_flight=in_flight, decisions=decisions)
