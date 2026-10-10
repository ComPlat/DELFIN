"""U6 — a refused shell command is answered with the form the gate accepts.

The gate refuses a command with a one-line message that names the deny pattern
but not the way forward. This module gives the wiring a concrete, gate-legal
replacement for the shapes that dominate the wave-13 refusal log
(``.gate/refusals_w13.txt``), so "command rejected by deny-pattern '\\|\\|' …"
turns into "… try ``test -f <file>; echo rc=$?``" instead of a dead end for the
model.

Contract
--------
``suggest(command, reason)`` returns one of

* a short note whose concrete commands (backticked) are gate-legal — the gate
  would not refuse them again,
* ``""`` when no allowed form exists (the refusal stands as-is).

The commands a suggestion names are re-checked through the real gate in
``tests/test_u6_refusal_hints.py`` (each backticked command runs through
``KitToolPermissions.matches_bash_deny``); a hint that names a command the
gate would refuse again is a bug.
"""

from __future__ import annotations

import re
from pathlib import Path

#: Shell commands found in suggestions (a backticked runnable form).
_BACKTICK_CMD = re.compile(r"`([^`]+)`")


def _gate_bin() -> str:
    """Absolute path of the session's sanctioned test runner at the MAIN
    DELFIN checkout: ``<repo>/.delfin/gate-tools/gate``.

    The module runs from a git worktree in production (``<repo>/.delfin/
    worktrees/<name>/delfin/agent/refusal_hints.py``), whose own ``.delfin``
    has no ``gate-tools`` — the gate lives only at the main repo root. Walk up
    from this module to the first ancestor whose ``.delfin/gate-tools/gate``
    exists, which lands on the main checkout in both the main-tree and
    worktree cases. Ever derived at runtime, never hardcoded: the checkout
    lives under a user home, and a home path must not ship in source.
    """
    cur = Path(__file__).resolve().parent  # .../delfin/agent
    for _ in range(16):
        gate = cur / ".delfin" / "gate-tools" / "gate"
        if gate.is_file():
            return str(gate)
        if cur.parent == cur:  # reached the filesystem root, stop
            break
        cur = cur.parent
    return str(cur / ".delfin" / "gate-tools" / "gate")


def _concrete(command: str, *, drop_first: bool = True) -> str:
    """A concrete suggested command keeps only the backticked tail command.

    Suggestions are prose plus one backticked runnable form; this returns the
    runnable form so the tests can hand it to the real gate. For the four
    command-shaped returns that's the whole backticked snippet; for pure tool
    notes there is none and this returns "".
    """
    m = _BACKTICK_CMD.search(command)
    return m.group(1).strip() if m else ""


def _first_operand(command: str) -> str:
    """The first bare path-like operand past the program and any flags."""
    # Split not only on whitespace but also on the joiners (in case they are
    # glued to a token, e.g. ``cat file||echo x``) so the operand is a clean
    # file path, never ``file||echo``.
    words = [
        part
        for tok in command.split()
        for part in re.split(r"\s*\|\|\s*|\s*&&\s*|\s*;\s*", tok)
        if part not in ("&&", "||", "|", ";", "")
    ]
    rest = words[1:]                      # past the program itself
    for tok in rest:
        if not tok.startswith("-"):
            return tok
    return ""


# --- shape helpers ---------------------------------------------------------


def _hint_double_pipe(command: str) -> str:
    target = _first_operand(command)
    if target:
        return (
            "the gate refuses the 'or' joiner. Read the file with the "
            "read_file tool, or probe with "
            f"`test -f {target}; echo rc=$?`"
        )
    return (
        "the gate refuses the 'or' joiner. Run the two sides as two separate "
        "bash calls — the fallback as its own command (a failed probe is told "
        "by its exit code, read from the second call), never chained by the "
        "or-joiner."
    )


def _hint_read_before_write(command: str) -> str:
    target = _first_operand(command)
    if target:
        return (
            "write_file refuses to overwrite without a read. Call read_file on "
            f"{target!r} first (so the change is made against what is actually "
            "on disk), then write_file. There is no shell form — it is a "
            "read-then-write tool contract."
        )
    return (
        "write_file refuses to overwrite without a read: call read_file on the "
        "file you mean to change first, then write_file against it. There is "
        "no shell form — it is a read-then-write tool contract."
    )


def _hint_pytest(command: str) -> str:
    """Answer the pytest refusal (bare ``pytest X`` or ``python -m pytest X``)
    with the gate tool, never ``python -m pytest``: the operator refuses
    pytest outside the gate, so suggesting it would send the model straight
    into another refusal."""
    m = re.search(r"(?:-m\s+)?(?:pytest|py\.test)\b(.*)$", command)
    targets = [
        t for t in (m.group(1) or "").split() if not t.startswith("-") and t
    ] if m else []
    target = targets[0] if targets else ""
    joiner = (
        " and the '||' or-joiner with it"
        if re.search(r"\|\|", command)
        else ""
    )
    head = (
        "pytest outside the gate is refused" + joiner
        + "; the session's tests run through the gate tool, not by invoking "
        "pytest on the command line. "
    )
    if target:
        return head + f"From the worktree root run `{_gate_bin()} {target}`."
    return (
        head + "From the worktree root run the gate with a tests/*.py file, "
        f"e.g. `{_gate_bin()} tests/test_a.py`."
    )


def _hint_write_scope(command: str) -> str:
    return (
        "the write scope limits you to your owned files and the scratch dir "
        "(.gate, tests/test_zz_*). A write outside those goes to the confirm "
        "gate: name the path and what you need it for, and the user can allow "
        "it — or move the target inside the scope."
    )


def _hint_delete(command: str) -> str:
    target = _first_operand(command)
    if target:
        return (
            "deleting is never the way to get rid of a file. Move it aside "
            f"instead: `mv {target} .gate/` (the scratch dir) — nothing is "
            "destroyed, and the user can still see it there."
        )
    return (
        "deleting is never the way to get rid of a file: move it aside with "
        "`mv <file> .gate/` (the scratch dir) instead — nothing is destroyed."
    )


def _hint_remote_fetch(_command: str) -> str:
    return (
        "fetching from a remote over the shell is hard-denied in every mode "
        "(the deny list wins over the allow list, even when the verb looks "
        "read-only). There is no allowed shell form: tell the user what you "
        "need from the remote and let them fetch. Do not try to reach it "
        "another way."
    )


def _hint_scheduler(_command: str) -> str:
    return (
        "the job-queue tools are not driven from the shell here — the DELFIN "
        "dashboard owns them. Drive the job through the calculator's own "
        "actions (submit, recalc) instead of the raw shell verb."
    )


def _hint_stderr(command: str) -> str:
    """A command carrying a stderr redirect the gate refuses; the redirect
    adds nothing for reading the output, so the plain read is the allowed
    form."""
    reduced = re.sub(r"\s*2(?:>&?\d*|>/[^\s|&;]+)\s*", " ", command)
    reduced = re.sub(r"\s+", " ", reduced).strip()
    return (
        "the stderr redirect is what the gate refused — it adds nothing for "
        f"reading the output. Run `{reduced}` instead."
    )


def _hint_semicolon_chain(_command: str) -> str:
    """A ';'-joined chain of read steps, often with echo separators; each
    side is allowed on its own, so split the chain into separate calls."""
    return (
        "the ';'-joined chain is what the gate refused, not any one step. "
        "Run the two sides as two separate bash calls; a separator echo is "
        "unnecessary — read each command's own output."
    )


# --- public API ------------------------------------------------------------


def suggest(command: str, reason: str = "") -> str:
    """The gate-legal way to do what the refused ``command`` wanted.

    ``command`` is the whole refused command line; ``reason`` is the deny
    pattern (or refusal fragment) the gate reported — used to disambiguate
    shapes that share a command spelling. Returns a note naming the allowed
    form (backticked commands are gate-legal), or ``""`` when no allowed form
    exists.
    """
    if not command or not command.strip():
        return ""

    cmd = command.strip()
    reason = reason or ""   # reason may be None from a bare call
    joined = f"{cmd}\n{reason}"

    # Harmful shapes first: a write that destroys has no safe fallback, and a
    # remote/scheduler verb is decided by role, not by shape — say "no allowed
    # form" rather than invent a route around a hard deny.
    if re.search(r"\brm\s+", cmd):
        return _hint_delete(cmd)

    if re.search(r"(?:^|[;&|]\s*|-m\s+)(?:pytest|py\.test)\b", joined):
        return _hint_pytest(cmd)

    if re.search(r"\|\|", cmd):
        return _hint_double_pipe(cmd)

    if re.search(r"\bgit\s+fetch\b", joined) or "fetch" in reason.lower():
        return _hint_remote_fetch(cmd)
    if re.search(r"\b(?:squeue|sbatch|scontrol)\b", joined):
        return _hint_scheduler(cmd)

    # Read-only compound shapes with a plain allowed form: drop the stderr
    # redirect, or split a ';'-chain into separate calls.
    if re.search(r"\s*2(?:>&?\d*|>/[^\s|&;]+)\s*", cmd):
        return _hint_stderr(cmd)

    if re.search(r"\s*;\s*", cmd):
        return _hint_semicolon_chain(cmd)

    # Writes that must not be silent.
    low = joined.lower()
    if ("write scope" in reason.lower() or "outside what you may modify" in low
            or "may modify" in low):
        return _hint_write_scope(cmd)
    if ("overwrite" in low or "read_file first" in reason.lower()
            or "call read_file first" in reason.lower()):
        return _hint_read_before_write(cmd)

    # A known shape with no confident allowed form: stay silent rather than
    # advise something the gate may refuse again.
    return ""
