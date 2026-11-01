"""swarm triage: the user's approval queue for the swarm (package SW1).

The rule is an ALLOW-list, deliberately narrow. An item may be
auto-answered only when BOTH hold:

- its kind is a plain ``session_message`` between two identities the
  queue's roster knows as members of THIS swarm (roster supplied by the
  caller; never inferred from the item's own text), or its kind is
  ``test_output`` from the member's own identity; and
- no narrowing rule in the queue's ``rules`` demotes it.

EVERYTHING ELSE escalates to the user -- a push request, a permission
request, a message from an unknown identity, anything whose wording
asks for an act only a person may grant. Nothing in the item's own data
can widen this: rules may only narrow, sender strings are not trusted,
and the verdict is re-derived inside ``triage_answer`` from the ITEM
(a caller-passed verdict dict is never trusted).

Auto-answers are logged to the swarm's STATE directory (outside the
workspace, resolved per call like DELFIN's other per-workspace state),
so what the machine said on the user's behalf is auditable and nothing
is written into the user's tree.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any


class TriageError(RuntimeError):
    """A triage action cannot be taken, with the reason in prose."""


def _log_dir() -> Path:
    """The log lives in per-workspace state, NOT in the workspace:
    resolved per call so a redirected home moves it, and never written
    inside a checkout (finding T4/D14: a probe wrote
    .swarm/triage_log.jsonl into a user tree)."""
    from .. import state_paths

    return state_paths.workspace_state_dir(Path.cwd())


def _log_path() -> Path:
    return _log_dir() / "triage_log.jsonl"


def classify(item: dict[str, Any]) -> dict[str, Any]:
    """Decide what may happen to one queued item. Read-only.

    ALLOW-list: every failure to match a permitted shape escalates.
    The queue may carry ``roster`` (the member identities of this
    swarm) and ``rules`` (narrowing only); neither is trusted beyond
    narrowing.
    """
    roster = frozenset(item.get("roster") or ())
    kind = item.get("kind", "")
    sender = item.get("from", "")
    member = item.get("member", "")

    def _escalate(reason: str) -> dict[str, Any]:
        return {"action": "escalate", "reason": reason}

    # Both ends of a partner message must be known members of THIS
    # swarm (T2: 'from' is self-claimed item data; an unknown identity
    # is not a partner). Unknown member or sender -> escalate.
    if not member or member not in roster:
        return _escalate("recipient is not a member of this swarm")
    if kind == "session_message":
        if not sender or sender not in roster:
            return _escalate("sender is not a member of this swarm")
        return {"action": "auto_answer",
                "reason": "plain partner message between swarm members"}
    if kind == "test_output" and sender == member:
        return {"action": "auto_answer",
                "reason": "a member reading its own test output"}
    return _escalate("not on the auto-answer list")


def _apply_rules(verdict: dict[str, Any], item: dict[str, Any]) -> dict:
    """The queue's rule file may NARROW only (T3): a rule demotes an
    auto-answerable item to escalation. A rule that tries to promote an
    escalation is ignored. Shape: rules["never_auto_answer"] is a list
    of "kind:identity" strings."""
    if verdict.get("action") != "auto_answer":
        return verdict
    rules = item.get("rules") or {}
    never = rules.get("never_auto_answer") or []
    key = f"{item.get('kind', '')}:{item.get('from', '')}"
    if key in never:
        return {"action": "escalate",
                "reason": "narrowed by the swarm's triage rules"}
    return verdict


def triage_list(queue: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Sort the queue: escalations first, then auto-answerables."""
    enriched = [dict(item, **_apply_rules(classify(item), item))
                for item in queue]
    escalations = [i for i in enriched if i["action"] == "escalate"]
    autos = [i for i in enriched if i["action"] == "auto_answer"]
    return escalations + autos


def triage_answer(item: dict[str, Any], verdict: dict[str, Any],
                  *, answer: str) -> dict[str, Any]:
    """Answer one item on the user's behalf; log it.

    The verdict is RE-DERIVED from the item here (T4/D13: a caller's
    verdict dict is untrusted data); only an item that classifies as
    auto-answerable -- after the narrowing rules -- may be answered.
    """
    derived = _apply_rules(classify(item), item)
    if derived.get("action") != "auto_answer":
        raise TriageError(
            f"item {item.get('id')!r} is {derived.get('action')!r} "
            f"({derived.get('reason')}) -- it goes to the user, not to "
            f"an auto-answer")
    log_path = _log_path()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    entry = {
        "ts": time.time(),
        "item_id": item.get("id"),
        "member": item.get("member"),
        "kind": item.get("kind"),
        "action": "auto_answer",
        "reason": derived.get("reason"),
        "answer": answer,
    }
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(entry, sort_keys=True) + "\n")
    return {"answered": True, "item_id": item.get("id"),
            "logged": str(log_path)}
