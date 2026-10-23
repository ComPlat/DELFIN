"""Per-session step journal: what already ran before a crash.

The mid-turn checkpoint (``session_store.save_turn_checkpoint``) records ONE
coarse snapshot per turn: the task, the partial answer and how many tool
calls had run. It cannot answer "step N of that turn — did it already
happen?". When a session dies mid-turn (node restart, OOM, a killed pane)
the restart begins the turn again, and a tool call that already ran — a
commit, a job submit, a message — is either repeated or silently lost.

This module keeps an APPEND-ONLY journal of steps: one line per tool call,
written in two halves.

* :func:`start_step` — appended BEFORE the call runs: turn, step number,
  tool name, a digest of the arguments, and whether the tool can change
  anything outside itself (side effect).
* :func:`finish_step` — appended AFTER: the step number, the result digest
  and whether the call succeeded. A step whose start exists and whose
  finish does not is an OPEN step: it "may or may not have run".

Storage: one JSONL file per session, ``<session>.steps.jsonl``, under the
SAME directory the session store uses (``session_store._SESSIONS_DIR``), so
it sits next to the session it belongs to and moves when the tests redirect
that directory. The file is created 0600 and appended through
``state_paths.open_append``, the same primitive the tool trace uses.

Nothing secret is ever written: arguments and results go in as digests
built by :mod:`delfin.agent.tool_digest` (secret-scrubbed, bounded), and a
short redacted head of the arguments for the resume line. Every write is
best-effort — a journal that cannot be written must never break the tool
call it is recording.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

# Read-only tools: calls whose ONLY effect is the value they return. This
# mirrors the plan-mode read-only surface (api_client._PLAN_READONLY_TOOLS,
# api_client.py:7119) minus the entries that belong to the agent's own
# planning machinery (task boards, plan exit, verdicts) — those change
# nothing on the machine either, but they are per-session state a resume
# rebuilds anyway, and marking them side-effecting errs on the safe side:
# an open step of this class is shown, never auto-re-run.
READONLY_TOOLS: frozenset[str] = frozenset({
    "read_file", "grep_file", "list_files", "find_definition",
    "find_references", "notebook_read", "view_image", "read_document",
    "compare_tables", "sum_column", "list_docs", "list_sections",
    "read_section", "search_docs", "search_calcs", "get_calc_info",
    "calc_summary", "check_environment", "list_changes_made",
    "history_search", "history_get", "task_list", "task_get",
    "session_search", "bash_status", "bash_output", "subagent_result",
    "cron_list", "extract_orbital_energies", "extract_excited_states",
    "extract_dipole", "extract_mulliken_charges", "extract_loewdin_charges",
    "extract_imaginary_frequencies", "extract_vibrational_modes",
    "extract_thermochem", "extract_optimization_trajectory",
    "extract_scf_convergence", "parse_orca_output", "extract_energy_table",
    "extract_delfin_json", "extract_calc_summary_table",
    "extract_excited_states_table", "find_orca_errors",
})

_MAX_ARGS_HEAD = 160
_MAX_NOTE_LINES = 40   # a read looks at the last lines only
#: Bytes read from the end of the journal per read. Every tool call reads
#: it (auto_start), so a read must not grow with the session: 40 lines of
#: at most ~1 KB each fit with room to spare.
_TAIL_BYTES = 256 * 1024


def _args_head(arguments: Any) -> str:
    """A short, redacted head of the arguments for the resume line.

    Goes through the secret scrubber (output_guard.scrub_secrets, the same
    redactor tool_digest uses), so what reaches disk is bounded and
    secret-free even when the digest itself is too coarse to show.
    """
    from .output_guard import scrub_secrets
    try:
        if isinstance(arguments, dict):
            parts = [f"{k}={arguments[k]!r}" for k in sorted(arguments)]
            text = " ".join(parts)
        else:
            text = repr(arguments)
    except Exception:
        text = ""
    text = " ".join(scrub_secrets(text).split())
    if len(text) > _MAX_ARGS_HEAD:
        text = text[: _MAX_ARGS_HEAD - 1] + "…"
    return text


def side_effecting(tool: str) -> bool:
    """True when the tool may change anything outside its return value.

    Everything not in :data:`READONLY_TOOLS` counts as side-effecting —
    the classification errs towards "do not re-run": bash, writes,
    submits, messages, MCP calls of every kind.
    """
    return tool not in READONLY_TOOLS


# ---------------------------------------------------------------------------
# Storage
# ---------------------------------------------------------------------------

def journal_path(session: str, *, root: "str | Path | None" = None) -> Path:
    """The journal file of one session.

    Resolved PER CALL from ``session_store._SESSIONS_DIR`` so the state
    redirect (state_paths.USER_STATE_SINKS, which owns that attribute)
    moves this journal with the rest of the per-run state — the wave-15b
    rule that no store is resolved at import time. An explicit ``root``
    (tests, callers that keep their own dir) wins.
    """
    from . import session_store as _ss
    base = Path(root) if root is not None else Path(_ss._SESSIONS_DIR)
    safe = "".join(
        c if (c.isalnum() or c in "-_.") else "_" for c in (session or "default")
    )[:80] or "default"
    return base / f"{safe}.steps.jsonl"


def _append(path: Path, entry: dict[str, Any]) -> bool:
    """Append one JSON line, 0600, best-effort. True when written.

    A torn append (the process died mid-``write``) leaves its last line
    without the trailing newline; a bare append after that would merge
    the next record INTO the torn line and lose both. The guard newline
    first costs one byte and keeps every later record parseable on its
    own: the torn line is still corrupt (skipped with a note), but the
    records after it survive whole.
    """
    from .state_paths import open_append, ensure_dir
    try:
        ensure_dir(path.parent)
        line = json.dumps(entry, ensure_ascii=False, sort_keys=True)
        torn = False
        try:
            size = path.stat().st_size
            if size:
                with path.open("rb") as fh:
                    fh.seek(size - 1)
                    torn = fh.read(1) != b"\n"
        except OSError:
            torn = False
        with open_append(path) as fh:
            if torn:
                fh.write("\n")
            fh.write(line + "\n")
        return True
    except Exception:
        return False


def start_step(session: str, *, turn: int, step: int, tool: str,
               arguments: Any = None) -> bool:
    """Record that a tool call is about to run. Returns True when written.

    Called BEFORE the call executes (the api_client hook, phase 3): the
    write landing is what makes "it started" true on disk, and the call's
    own crash cannot take the record with it.
    """
    return _append(journal_path(session), {
        "kind": "start",
        "ts": time.time(),
        "turn": int(turn),
        "step": int(step),
        "tool": str(tool),
        "side_effect": side_effecting(tool),
        "args_digest": _args_head(arguments),
    })


def finish_step(session: str, *, turn: int, step: int, tool: str,
                ok: bool, result: Any = "") -> bool:
    """Record that the call returned. The result goes in as a digest.

    ``turn`` travels with the finish record: a finish matches the start
    it closes only on (turn, step). Step numbers restart each turn — a
    finish without the turn would let a stale finish from an older turn
    mask a genuinely open step of the current one. ``tool`` travels too:
    a replayed result is only useful with the tool it came from.
    """
    from . import tool_digest
    content = result if isinstance(result, str) else str(result or "")
    return _append(journal_path(session), {
        "kind": "finish",
        "ts": time.time(),
        "turn": int(turn),
        "step": int(step),
        "tool": str(tool),
        "ok": bool(ok),
        "result_digest": tool_digest.digest("", {}, content),
    })


# ---------------------------------------------------------------------------
# The tool-loop bracket (phase 3): derive turn/step from the journal tail
# ---------------------------------------------------------------------------

def auto_start(session: str, *, tool: str, arguments: Any = None,
               turn: int | None = None,
               ) -> dict[str, Any]:
    """Derive (turn, step) from the journal and record the call start.

    The protected api_client hook (phase 3) has no turn/step counter of
    its own — this one is derived from the tail: step = last recorded
    step + 1 within the current turn, turn rolling over when the caller
    announces a new turn (``turn`` given) or when the previous open
    bracket ended a turn. Returns the record dict with its derived
    turn/step, which auto_finish consumes; on a write failure returns a
    record marked ``recorded: False`` so the finish is skipped rather
    than paired to a line that never landed.
    """
    record = {"recorded": False, "turn": 0, "step": 0}
    if turn is None:
        turn = _current_turn(session)
    step = _next_step(session, turn)
    ok = start_step(session, turn=turn, step=step, tool=tool,
                    arguments=arguments)
    record.update({"recorded": ok, "turn": turn, "step": step,
                   "tool": tool})
    return record


def auto_finish(session: str, record: dict[str, Any] | None, *,
                ok: bool, result: Any = "") -> bool:
    """Close the bracket auto_start opened. No-op without a real start."""
    if not isinstance(record, dict) or not record.get("recorded"):
        return False
    try:
        return finish_step(session, turn=record["turn"],
                           step=record["step"], tool=record["tool"],
                           ok=ok, result=result)
    except Exception:
        # Never let the journal take a tool call down with it: the hook
        # wraps real dispatch, and a failed close must not surface as a
        # tool error. An open step reads as "may or may not have run",
        # which is the honest state when the close never landed.
        return False


def _current_turn(session: str) -> int:
    """The turn the journal is currently in; 1 for an empty journal.

    Consecutive calls continue the same turn — steps must increment, not
    restart — so the fallback keeps the highest turn already recorded. A
    NEW turn is announced explicitly by the caller (auto_start(turn=...)),
    which also restarts the step numbering via _next_step.
    """
    steps, _ = read_steps(session)
    turns = [int(e["turn"]) for e in steps
             if isinstance(e.get("turn"), int) or _intable(e.get("turn"))]
    return max(turns) if turns else 1


def _intable(v: Any) -> bool:
    try:
        int(v)
        return True
    except (TypeError, ValueError):
        return False


def _next_step(session: str, turn: int) -> int:
    """1 + the highest step recorded in ``turn`` (steps restart per turn).

    A non-int turn (a corrupt tail line) is simply not in ``turn`` —
    the helper's no-raise contract (leaned on by the api_client hook,
    api_client.py:16660-16662) holds even against a broken journal.
    """
    steps, _ = read_steps(session)
    in_turn = [int(e["step"]) for e in steps
               if e.get("turn") == turn]
    return max(in_turn, default=0) + 1


# ---------------------------------------------------------------------------
# Reading: what the resume path sees
# ---------------------------------------------------------------------------

def read_steps(session: str, *, root: "str | Path | None" = None,
               max_lines: int | None = None,
              ) -> tuple[list[dict[str, Any]], list[str]]:
    """Every parseable journal line of a session, oldest first.

    Returns ``(steps, notes)``. A line that is not valid JSON, or JSON of
    the wrong shape, is SKIPPED and produces one note naming the line
    number — a corrupt or half-written line (a crash can tear an append)
    costs that one record, never the session: the resume must always come
    up with something, and never crash. ``notes`` also carries the
    self-trim report when the file has grown past its bound.

    The bound is read from the module global at call time (an ``int | None``
    default, not a baked-in constant) so an operator retuning the bound —
    or a test exercising it — actually changes the behaviour.
    """
    if max_lines is None:
        max_lines = _MAX_NOTE_LINES
    p = journal_path(session, root=root)
    if not p.is_file():
        if p.exists():
            # A failed close can leave a directory (or anything else
            # that is not a regular file) behind: the reader must say
            # so, not read it as silently "no journal" — a step whose
            # finish_step returned False would otherwise vanish from
            # every resume surface.
            return [], [f"step journal unreadable: {p.name} "
                        "(exists, not a regular file)"]
        return [], []
    notes: list[str] = []
    steps: list[dict[str, Any]] = []
    try:
        with p.open("rb") as fh:
            size = fh.seek(0, 2)
            start = max(0, size - _TAIL_BYTES)
            fh.seek(start)
            data = fh.read()
    except OSError:
        return [], ["step journal unreadable: " + p.name]
    lines = data.decode("utf-8", errors="replace").splitlines()
    if start and lines:
        lines = lines[1:]          # the first line is cut at the seek point
    if len(lines) > max_lines or start:
        lines = lines[-max_lines:]
        notes.append(
            f"step journal trimmed to the last {max_lines} lines")
    for i, raw in enumerate(lines, 1):
        raw = raw.strip()
        if not raw:
            continue
        try:
            entry = json.loads(raw)
        except (json.JSONDecodeError, ValueError):
            notes.append(f"step journal line {i} is corrupt — skipped")
            continue
        if not isinstance(entry, dict) or entry.get("kind") not in ("start", "finish"):
            notes.append(f"step journal line {i} has an unknown shape — skipped")
            continue
        steps.append(entry)
    return steps, notes


def open_steps(session: str, *, root: "str | Path | None" = None,
              ) -> tuple[list[dict[str, Any]], list[str]]:
    """Steps that started and never finished, oldest first.

    A start with no matching finish after it is an OPEN step — on resume it
    reads as "may or may not have run", no matter how plausible the rest
    of the journal looks: the crash landed somewhere inside the call, and
    only the call's own evidence (the workspace, a job list, git log) can
    tell whether it completed.
    """
    steps, notes = read_steps(session, root=root)
    finished: set[tuple[int, int]] = set()
    for entry in steps:
        if entry.get("kind") == "finish":
            key = _turn_step_key(entry)
            if key is not None:
                finished.add(key)
    open_out = [e for e in steps
                if e.get("kind") == "start"
                and _turn_step_key(e) not in finished]
    return open_out, notes


def _turn_step_key(entry: dict[str, Any]) -> "tuple[int, int] | None":
    """(turn, step) of a journal record, or None when unusable.

    Step numbers restart each turn, so a finish matches its start only
    on BOTH fields — a stale finish of step 1 from turn 7 must not mask
    the step-1 start of turn 8.
    """
    try:
        return int(entry.get("turn", -1)), int(entry.get("step", -1))
    except (TypeError, ValueError):
        return None


def finished_result(session: str, *, step: int, turn: int | None = None,
                    root: "str | Path | None" = None) -> dict[str, Any] | None:
    """The recorded finish of one step, or None when it never finished.

    With ``turn`` given (the resume path knows the turn it is resuming),
    only the finish of that (turn, step) matches. Without it, the newest
    finish of that step number is returned — the pre-turn-field journal
    shape, kept for older files.
    """
    steps, notes = read_steps(session, root=root)
    for entry in reversed(steps):
        if entry.get("kind") != "finish":
            continue
        try:
            if int(entry.get("step", -1)) != int(step):
                continue
        except (TypeError, ValueError):
            continue
        if turn is not None:
            try:
                if int(entry.get("turn", -1)) != int(turn):
                    continue
            except (TypeError, ValueError):
                continue
        return entry
    return None


def last_finished_step(session: str, *,
                       root: "str | Path | None" = None,
                       side_effect_only: bool = False,
                       ) -> dict[str, Any] | None:
    """The newest finished step — the number compaction keeps.

    ``working_state`` carries this into the compaction block (with
    ``side_effect_only``: the newest step that could change something) so
    a later context rebuild knows how far the turn had actually got, not
    just how many tool calls a coarse checkpoint had counted.
    """
    steps, _ = read_steps(session, root=root)
    for entry in reversed(steps):
        if entry.get("kind") != "finish":
            continue
        if side_effect_only and not side_effecting(str(entry.get("tool"))):
            continue
        return entry
    return None


def render_resume_note(session: str, *, root: "str | Path | None" = None,
                       max_open: int = 5) -> str | None:
    """The step-resume note, or None when every step finished.

    One line per open step. A side-effecting step names the tool and the
    arguments head and says CHECK BEFORE REPEATING — the person's words
    for this rule; the note never instructs a re-run, it only says what
    is uncertain. A read-only open step says the same thing with a milder
    tail: it MAY be re-run, but its outcome should be checked too.
    Corrupt-line notes are appended, so a torn journal is reported where
    it is read, never swallowed.
    """
    open_out, notes = open_steps(session, root=root)
    # That a read looked at the tail only is how every long journal is
    # read, not something the resumed session has to know.
    notes = [n for n in notes if not n.startswith("step journal trimmed")]
    if not open_out and not notes:
        return None
    lines: list[str] = []
    if open_out:
        lines.append(
            "[recovered steps] The previous process died inside the turn. "
            "These steps started but never recorded a result:")
        for e in open_out[:max_open]:
            step_no = e.get("step", "?")
            tool = e.get("tool", "?")
            args_head = e.get("args_digest") or "(no arguments recorded)"
            if e.get("side_effect"):
                tail = ("This tool can change things outside itself — "
                        "step may or may not have run: CHECK BEFORE "
                        "REPEATING it; verify the workspace state first.")
            else:
                tail = ("step may or may not have run: check before "
                        "repeating; a re-run is allowed but verify its "
                        "outcome.")
            lines.append(
                f"  step {step_no} ({tool}, args: {args_head}) — {tail}")
        if len(open_out) > max_open:
            lines.append(f"  … and {len(open_out) - max_open} more open steps")
        finished = last_finished_step(session, root=root)
        if finished:
            lines.append(
                f"  Last recorded finished step: {finished.get('step', '?')}. "
                "Finished steps must NOT be repeated — reuse their results.")
    lines.extend(f"  ({n})" for n in notes)
    return "\n".join(lines)
