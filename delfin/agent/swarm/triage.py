"""swarm triage: the user's approval queue for the swarm (package SW1).

The rule is deliberately narrow: only PLAIN partner messages between
members of the same swarm, and a member reading its own test output,
may be auto-answered. Everything else -- push requests, permission
requests, anything touching the safety frame -- escalates to the user.
Auto-answers are logged, so what the machine said on the user's behalf
is always auditable.

The rule file (``.swarm/triage_rules.json``) may NARROW this further,
never widen it: a rule can only demote an auto-answerable item to
escalation, never promote an escalated one.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

TRIAGE_LOG = ".swarm/triage_log.jsonl"

#: Words that make an item escalate no matter its kind. Permission and
#: push requests are the user's to answer: answering one automatically
#: would widen what a member may do without a person deciding.
ESCALATE_MARKERS = (
    "push", "merge to main", "permission", "allow list", "rm -rf",
    "secret", "credential", "token", "sudo", "force",
)


class TriageError(RuntimeError):
    """A triage action cannot be taken, with the reason in prose."""


def classify(item: dict[str, Any]) -> dict[str, Any]:
    """Decide what may happen to one queued item. Read-only."""
    kind = item.get("kind", "")
    sender = item.get("from", "")
    member = item.get("member", "")
    text = (item.get("text") or "").lower()
    if any(marker in text for marker in ESCALATE_MARKERS):
        return {"action": "escalate",
                "reason": "touches something only the user may grant"}
    if kind == "session_message" and sender and sender != member:
        return {"action": "auto_answer",
                "reason": "plain partner message between swarm members"}
    if kind == "test_output" and sender == member:
        return {"action": "auto_answer",
                "reason": "a member reading its own test output"}
    return {"action": "escalate", "reason": "not on the auto-answer list"}


def triage_list(queue: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Sort the queue: escalations first, then auto-answerables."""
    enriched = [dict(item, **classify(item)) for item in queue]
    escalations = [i for i in enriched if i["action"] == "escalate"]
    autos = [i for i in enriched if i["action"] == "auto_answer"]
    return escalations + autos


def triage_answer(item: dict[str, Any], verdict: dict[str, Any],
                  *, answer: str) -> dict[str, Any]:
    """Answer one item on the user's behalf; log it.

    Only items classified ``auto_answer`` may be answered here; an
    escalated item goes to the user, and answering it from code is
    refused (the log is the receipt, not a loophole).
    """
    if verdict.get("action") != "auto_answer":
        raise TriageError(
            f"item {item.get('id')!r} is {verdict.get('action')!r} -- "
            f"it goes to the user, not to an auto-answer")
    log_path = Path.cwd() / TRIAGE_LOG
    log_path.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "ts": time.time(),
        "item_id": item.get("id"),
        "member": item.get("member"),
        "kind": item.get("kind"),
        "action": "auto_answer",
        "reason": verdict.get("reason"),
        "answer": answer,
    }
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, sort_keys=True) + "\n")
    return {"answered": True, "item_id": item.get("id"),
            "logged": str(log_path)}
