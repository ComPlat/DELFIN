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

    Reviewer s26 round-2 (G1-G10) rebuild, QS taste ruling taken: the
    floor is ALLOW-FIRST at the token level. The command text is
    tokenised (whole line, never a bare substring: F5/G10), each line of
    a multi-line script is checked (G4), env prefixes are peeled
    (FOO=bar, env FOO=bar: G1-G2), command substitution voids the token
    picture so the line is refused outright ($(...) hides what runs:
    G3), path-qualified binaries resolve to their basename (G5-G7), and
    git commit with ANY flag of the -a/--all family mass-commits (G8-G9).
    The url field is floor-checked like a binary (G10).
    """
    tool = str(action.get("tool") or "")
    url = str(action.get("url") or "")
    # G10, allow-first: a url is only legal on DELFIN's own web-fetch
    # tool; on any other tool (bash carrying a url) it is an outward
    # target the tool's own contract never checks.
    if url and tool not in ("web_fetch", "web"):
        return "url field on a tool that does not define it"
    commands = str(action.get("command") or "")
    if not commands and not url:
        return None
    for line in commands.splitlines() or [""]:
        reason = _floor_hit_one_line(line, tool)
        if reason:
            return reason
    return None


def _floor_hit_one_line(line: str, tool: str) -> str | None:
    tokens = _tokens(line)
    if tokens is None:
        return "command substitution hides what runs; write it as a visible command"
    if not tokens:
        return None
    # G1: peel VAR=value prefixes; G2: peel env(1) prefixes. Loop until
    # stable: "env FOO=bar curl" peels env, then the next VAR=value.
    while True:
        if tokens and "=" in tokens[0] and tokens[0].split("=")[0]:
            tokens = tokens[1:]
            continue
        if tokens and tokens[0] == "env":
            tokens = tokens[1:]
            continue
        break
    if not tokens:
        return None
    binary = tokens[0].rsplit("/", 1)[-1] if "/" in tokens[0] else tokens[0]
    # H-class (reviewer s26 round 3), QS ruling: wrappers must not hide
    # the refused binary. bash/sh -c recurse into the wrapped command;
    # nohup/timeout/xargs drop their own arguments and recurse. sudo is
    # refused outright in unattended mode — no shape where it runs.
    if binary == "sudo":
        return "sudo needs the person's word; unattended sessions do not raise their own rights"
    if binary == "xargs":
        return "xargs feeds its input to a command the line does not name; write the command out"
    if binary in ("bash", "sh"):
        wrapped = [t for t in tokens[1:] if t != "-c"]
        if wrapped:
            inner = wrapped[-1] if len(wrapped) == 1 else " ".join(wrapped)
            inner_hit = _floor_hit_one_line(inner, tool)
            if inner_hit:
                return inner_hit
            return ("bash -c of a command the floor cannot fully see; "
                    "write the command directly")
        return None
    if binary == "nohup":
        return _floor_hit_one_line(" ".join(tokens[1:]), tool)
    if binary == "timeout":
        rest = tokens[1:]
        if rest and rest[0].replace(".", "").isdigit():
            rest = rest[1:]
        return _floor_hit_one_line(" ".join(rest), tool)
    if binary in _FLOOR_REASONS:
        if binary == "curl" and len(tokens) > 1:
            # A "$(echo curl)" first token resolves to "curl" only by the
            # value of the substitution, which _tokens refuses: it cannot
            # appear here. A real multi-word argument list after a plain
            # curl binary is still curl.
            pass
        return _FLOOR_REASONS[binary]
    if binary in ("pip", "pip3") or (
            len(tokens) > 2 and tokens[:3] == ["python", "-m", "pip"]):
        return "package installation"
    if binary == "git":
        rest = tokens[1:]
        if rest and rest[0] == "push":
            return "push"
        if rest and rest[0] == "commit" and any(
                f == "-a" or f == "-am" or f == "--all" or f.startswith("-am")
                or f.startswith("--all=") for f in rest):
            # A mass-commit (everything tracked, no staged selection).
            # ``git commit -m`` alone is a normal commit, NOT the floor.
            return "mass-commit of every tracked change"
    if tool in ("ssh", "scp", "rsync", "sftp"):
        return "network transfer outside DELFIN"
    return None


_FLOOR_REASONS = {
    "ssh": "ssh/scp/rsync to another host",
    "scp": "ssh/scp/rsync to another host",
    "rsync": "ssh/scp/rsync to another host",
    "sftp": "ssh/scp/rsync to another host",
    "nc": "network transfer outside DELFIN",
    "ncat": "network transfer outside DELFIN",
    "netcat": "network transfer outside DELFIN",
    "wget": "network fetch outside DELFIN's web tools",
    "curl": "network fetch outside DELFIN's web tools",
}


def _tokens(line: str) -> list[str] | None:
    """Tokenise one shell line, lowercased; None when the argv picture is
    false (command substitution: what runs is not what the line shows,
    and a false argv picture must never read as an allow)."""
    stripped = line.strip()
    if not stripped:
        return []
    if ("$(" in stripped or "`" in stripped or "<(" in stripped
            or "(" in stripped):
        # H-class: <( ) / >( ) process substitution runs the command
        # inside the parentheses, so the argv picture is false too.
        return None
    try:
        import shlex
        return [t.lower() for t in shlex.split(stripped)]
    except ValueError:
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
