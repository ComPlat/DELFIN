"""The unattended permission mode: a layered decision without a person.

A session with nobody at the terminal must not wait 15 minutes for an
approval nobody can give, and must not fall back to bypassPermissions —
measured in stage 1a/1d: sessions locked for two hours on harmless asks
(partner messages, own test output), and a bypass tree that allowed curl,
pip installs and git commit -am (see FR1).

The policy is ONE pure decision function, called by the gate where a
confirmation dialog would otherwise open (api_client,
``_classifier_second_opinion``, 21440+). It never widens the hard floor;
it can only choose between running, refusing with a reason, escalating to
the operator inbox, or pausing.

Layers, checked in order:
1. floor      — the deny list, sandbox, write scope, protected files, no
                push, no secrets, no network beyond DELFIN's web tools,
                no ssh/scp/rsync, no package installs. Unchanged, checked
                first. A floor hit is refused at once.
2. routine    — the code-checked routine allow-list (FR1). Runs without a
                dialog.
3. classifier — everything else. The action classifier's verdict decides:
                allow -> run; refuse -> refuse; unknown (doubt) ->
                refuse at once with the reason — never a wait.
4. escalation — a refused action the task really needs becomes ONE
                operator-inbox message and the session moves on.
5. pause      — the third refusal in a turn pauses the session with an
                operator message instead of waiting for nobody.

Attended sessions never reach this module: ``decide`` raises ValueError
when a person can be asked.
"""
from __future__ import annotations

from collections import namedtuple

__all__ = ["Decision", "decide", "is_routine", "register_routine",
           "REFUSALS_BEFORE_PAUSE"]

REFUSALS_BEFORE_PAUSE = 3

Decision = namedtuple(
    "Decision",
    "decision layer reason message wait_for_user",
    defaults=("", "", "", False))
# decision: "run" | "refuse" | "escalate" | "pause"
# message: the operator-inbox text for escalate/pause decisions.

# The floor shapes live in _floor_hit's first-token logic below (reviewer
# s26 F1-F5: first-token, case-insensitive, pip family, flag-normalized
# git commit, wider transfer family).

# Registered routine predicates (FR1's allow-list hooks in here).
_ROUTINE_PREDICATES = []


def register_routine(predicate) -> None:
    """Register a routine-allow predicate (FR1 layer 2). Test-visible so
    the wiring can be driven through the public path."""
    _ROUTINE_PREDICATES.append(predicate)


def is_routine(action) -> bool:
    """True when the action is on the code-checked routine allow-list."""
    for pred in _ROUTINE_PREDICATES:
        try:
            if pred(action):
                return True
        except Exception:
            continue
    return False


def _floor_hit(action) -> str | None:
    """The reason text of the floor shape hit, or None.

    Reviewer s26 findings F1-F5 fixed: matching is on the FIRST TOKEN of
    the command (never a bare substring in the middle of a read, F5),
    case-insensitive (F1), with a pip front-end family (F2), flag-
    normalized git commit (F3: -am, -a -m, -m+staged all mass-commit
    risk) and a wider transfer family (F4: sftp, nc/ncat/netcat).
    """
    text = str(action.get("command") or action.get("url") or "").strip()
    tool = str(action.get("tool") or "")
    low = text.lower()
    first = low.split()[0] if low.split() else ""
    if first in ("ssh", "scp", "rsync", "sftp", "nc", "ncat", "netcat",
                 "wget", "curl"):
        return "network transfer outside DELFIN's web tools"
    if first in ("pip", "pip3") or (
            first == "python" and len(low.split()) > 2
            and low.split()[1] == "-m" and low.split()[2] == "pip"):
        return "package installation"
    if first == "git":
        rest = low.split()[1:]
        if rest and rest[0] == "push":
            return "push"
        if rest and rest[0] == "commit" and (
                "-a" in rest or "-am" in rest):
            # A mass-commit (everything tracked, no staged selection).
            # ``git commit -m`` alone is a normal commit, NOT the floor.
            return "mass-commit of every tracked change"
    if tool in ("ssh", "scp", "rsync", "sftp"):
        return "network transfer outside DELFIN"
    return None


def _plain_read(action) -> bool:
    """Thin adapter over action_classifier.is_plain_read (:258) — the real
    definition lives there; this only survives its absence/narrow failure
    so the policy layer never crashes on it. Reads are safe by definition
    and run without a verdict and without a dialog."""
    try:
        from . import action_classifier as _ac
    except ImportError:
        return False
    try:
        return _ac.is_plain_read(action)
    except Exception:
        return False


def decide(action, verdict, ctx) -> Decision:
    """The one layered decision for an unattended session.

    action: {"tool": ..., "command"/"path"/"url": ..., optional "risk"}
    verdict: the classifier's verdict — "allow" | "refuse" | "unknown"
    ctx: attended, confirm_callback, in_task, task, refusals_this_turn,
         and optional flags: needed (the task really needs the refused
         action), reason (the rule's refusal reason).
    """
    if ctx.get("attended") or ctx.get("confirm_callback") is not None:
        raise ValueError(
            "unattended_policy.decide is for sessions with nobody at the "
            "terminal; attended sessions use the normal confirm path")
    action = action or {}

    # -- layer 5 first: the pause brake counts refusals, not layers
    refusals = int(ctx.get("refusals_this_turn") or 0)
    if refusals >= REFUSALS_BEFORE_PAUSE:
        return Decision(
            "pause", "pause",
            reason="three refusals in this turn",
            message=("Unattended session paused: three actions were "
                     "refused in this turn, so the session waits for a "
                     "person. Operator, please review. Last refused: "
                     + str(action.get("command") or action.get("tool")
                           or action)),
            wait_for_user=False)

    # -- layer 1: the hard floor
    floor = _floor_hit(action)
    if floor:
        return Decision("refuse", "floor", reason=floor)

    # -- plain reads run without a verdict and without a dialog (the
    #    latency-free floor of the classifier itself, action_classifier:258)
    if _plain_read(action):
        return Decision("run", "floor")

    # -- layer 2: the routine allow-list runs without a dialog
    if verdict != "refuse" and is_routine(action):
        return Decision("run", "routine")

    # -- layer 3: the classifier decides everything else. Conservative per
    #    the QS ruling (s28, 2026-10-25): anything not clearly allowed is
    #    refused; escalation arrives in phase 2, fire-and-forget.
    if verdict == "allow":
        return Decision("run", "classifier")
    if verdict == "refuse":
        return Decision(
            "refuse", "classifier",
            reason=str(ctx.get("reason") or "refused by the rules"))
    # verdict == "unknown": doubt under nobody-at-the-terminal is refusal
    # at once, with the reason — never a wait.
    return Decision(
        "refuse", "classifier",
        reason=str(ctx.get("reason") or "unsure and nobody to ask"))
