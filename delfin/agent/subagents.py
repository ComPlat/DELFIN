"""Sub-agent delegation (.delfin Agent tool surface).

A sub-agent runs an isolated tool-calling loop: separate message
history, separate token budget, separate (usually narrower) tool
set. The parent agent invokes one with a single ``subagent`` tool
call and gets the final assistant message back as the tool result.

Use cases:

  - parallel research that would otherwise pollute the parent's
    context window (the parent only sees the summary)
  - read-only auditing of a candidate change
  - planning runs that should not be allowed to edit anything

Built-in subagent types::

    explore           read-only filesystem + grep + web search
    plan              plan mode, no writes/bash
    code-reviewer     read-only, geared toward review checklists
    general-purpose   full tool set (parent's permissions)

Each type maps to:

  - a system prompt (``SUBAGENT_PRESETS[type].system``)
  - a permission tweak (``SUBAGENT_PRESETS[type].mode``); the
    sub-agent's permissions are derived from the parent's by cloning
    and overriding ``mode``. The sandbox / deny-list / self-mod-guard
    inherit from the parent untouched.

Hard limits:

  - max 40 tool calls per sub-agent      (``_MAX_TOOL_CALLS``)
  - max 900 seconds wall-clock per run   (``_MAX_WALL_S``)
  - max 16000 tokens output              (``_MAX_OUTPUT_TOKENS``)

  All three are overridable via ``settings["agent"]["subagents"]``.

Failures don't propagate: a crash inside the sub-agent returns an
error string, never an exception.

Report verification (what the delegate SAYS vs. what it DID)
------------------------------------------------------------

A sub-agent's final message reaches the parent as a tool result and is
read as fact — the parent applies its own grounding guards to its own
answers, but had nothing to hold a delegate to. ``verify_subagent_report``
closes that gap: it classifies the claims in a report (test results,
file/location references, quantities, functional and completion claims)
and cross-checks each against the ONE record the delegate cannot
retroactively invent — the tool calls it actually made in that run. The
verdict rides along with the result (``attach_verification``), so the
parent reads it before the report body.

What this CANNOT catch — stated plainly, because the check is easy to
over-trust: it compares a report against its own trace, not against
reality. A delegate that runs the tests and then misreports the outcome,
reads a file and describes it wrongly, edits the wrong file, or draws a
false conclusion from real tool calls passes every rule here. A clean
verdict means "this report is consistent with what this run did", never
"this report is correct". Claims that matter still need the parent's own
verification.
"""

from __future__ import annotations

import dataclasses
import json
import os
import re
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING

from . import german as _de

# Structured returns — the JSON-Schema subset validator, JSON extraction
# and the schema-contract instruction moved to ``structured_output`` as a
# reusable harness service (any caller can request schema-validated JSON,
# not just this module). Re-imported (and re-exported via __all__) here so
# existing import paths stay valid.
from .structured_output import (
    extract_json_object,
    schema_instruction as _schema_instruction,
    validate_json_schema,
)

if TYPE_CHECKING:  # pragma: no cover
    from .api_client import KitToolPermissions, OpenAIClient


# Per-run hard caps. These are FALLBACK defaults; settings
# ["agent"]["subagents"] overrides them (see _subagent_limits). Wall-clock
# raised 60→300s: 60s truncated real exploration/research runs (especially on
# slower KIT/Qwen models) before they could report back — subagents were too
# short-leashed to be useful for anything but trivial lookups.
_MAX_TOOL_CALLS = 40
# Wall-clock, not call count, is the binding constraint for a delegated
# build task. Measured on a real delegation round (2026-07-29, KIT
# endpoint): the two runs that died at the 300 s cap had made only 10 and
# 3 tool calls — roughly 30 s per call, dominated by time-to-first-token
# on a large prompt. The four that finished needed 15-77 s; none was ever
# truncated on output. So the cap is raised and the other two budgets
# stay: they were never the limit.
_MAX_WALL_S = 900.0
_MAX_OUTPUT_TOKENS = 16000


def _subagent_limits() -> dict:
    """Resolve per-run subagent limits from settings, with fallbacks.

    ``settings["agent"]["subagents"]`` may set ``max_tool_calls``,
    ``max_wall_s`` and ``max_output_tokens`` — so Jerome can tune how deep and
    how long a delegated subagent may run from the dashboard, without code
    changes. Missing/zero values fall back to the module defaults above.
    """
    cfg: dict = {}
    try:
        from delfin.user_settings import load_settings
        cfg = (((load_settings() or {}).get("agent") or {})
               .get("subagents") or {})
    except Exception:
        cfg = {}

    def _num(key, default, cast):
        try:
            v = cfg.get(key)
            return cast(v) if v not in (None, "", 0) else default
        except Exception:
            return default

    return {
        "max_tool_calls": _num("max_tool_calls", _MAX_TOOL_CALLS, int),
        "max_wall_s": _num("max_wall_s", _MAX_WALL_S, float),
        "max_output_tokens": _num("max_output_tokens", _MAX_OUTPUT_TOKENS, int),
    }

_TELEMETRY_PATH = Path.home() / ".delfin" / "subagent_telemetry.jsonl"
_TELEMETRY_MAX_LINES = 5000

# Context-inheritance budgets. A spawned child used to start with a bare
# preset prompt — no DELFIN.MD/AGENTS.md project rules, no typed memories —
# so it could cheerfully violate standing project conventions the parent
# knows about. Both blocks are appended to the child's system prompt (unless
# ``inherit_context=False``) and are hard-capped so they can't crowd out the
# actual task briefing.
_INHERIT_RULES_MAX_CHARS = 2000
_INHERIT_MEMORY_MAX_CHARS = 1500
# Only these typed-memory kinds are inherited: they encode durable
# user/project preferences. Task-scoped or auto-distilled kinds stay out.
_INHERIT_MEMORY_TYPES = ("feedback", "user")


def _inherited_context(workspace, prompt: str) -> str:
    """Project rules + standing memories for a child's system prompt.

    (a) Project rules: ``load_project_memory`` walking up from the child's
        workspace (DELFIN.MD / AGENTS.md), capped at
        ``_INHERIT_RULES_MAX_CHARS``.
    (b) Standing memory: feedback/user typed memories from the child
        workspace's store, BM25-ranked against the subagent prompt so the
        most task-relevant ones survive the ``_INHERIT_MEMORY_MAX_CHARS``
        cap (ties broken by recency).

    Every step is best-effort: a memory failure must never break a spawn,
    so ANY exception collapses to omitting that block. Returns "" when the
    workspace is unknown or nothing is found.
    """
    if workspace is None:
        return ""
    parts: list[str] = []
    try:
        from .project_memory import load_project_memory
        rules = load_project_memory(
            cwd=workspace, max_chars=_INHERIT_RULES_MAX_CHARS)
        if rules and rules.strip():
            parts.append(
                "## Project rules (inherited from the parent workspace)\n"
                + rules.strip())
    except Exception:
        pass
    try:
        from . import memory_store as _ms
        recs = [r for r in _ms.list_typed_memories(workspace)
                if (r.get("type") or "") in _INHERIT_MEMORY_TYPES]
        if recs:
            texts = [
                " ".join(x for x in (r.get("name"), r.get("description"),
                                     r.get("body")) if x)
                for r in recs
            ]
            try:
                scores = _ms.bm25_scores(prompt or "", texts)
            except Exception:
                scores = [0.0] * len(recs)
            order = sorted(
                range(len(recs)),
                key=lambda i: (-scores[i],
                               -int(recs[i].get("updated_at") or 0)))
            lines: list[str] = []
            used = 0
            for i in order:
                r = recs[i]
                body = " ".join(
                    str(r.get("body") or r.get("description") or "").split())
                line = f"- [{r.get('type', 'user')}] {r.get('name', '')}: {body[:300]}"
                if used + len(line) > _INHERIT_MEMORY_MAX_CHARS and lines:
                    break
                lines.append(line)
                used += len(line) + 1
            if lines:
                parts.append("## Standing memory (inherited)\n"
                             + "\n".join(lines))
    except Exception:
        pass
    if not parts:
        return ""
    return "\n\n" + "\n\n".join(parts)


# ---------------------------------------------------------------------------
# Structured returns: validate_json_schema / extract_json_object /
# _schema_instruction moved to ``structured_output`` (imported at the top,
# still re-exported via __all__).
# ---------------------------------------------------------------------------


def _stream_text_only(sub_client, messages: list[dict], system_prompt: str,
                      max_output_tokens: int, deadline: float):
    """One extra text-only model round (used for the schema correction).

    Ignores tool activity — the correction round only reformats an answer
    the child already produced. Returns ``(text, in_tokens, out_tokens,
    error)`` and never raises.
    """
    parts: list[str] = []
    in_tok = out_tok = 0
    err = ""
    try:
        for event in sub_client.stream_message(
            messages=messages,
            system=system_prompt,
            max_tokens=max_output_tokens,
        ):
            if time.monotonic() > deadline:
                err = "wall-clock budget exhausted during schema correction"
                break
            if event.type == "text_delta" and event.text:
                parts.append(event.text)
            elif event.type == "message_delta":
                in_tok = max(in_tok, event.input_tokens)
                out_tok = max(out_tok, event.output_tokens)
    except Exception as exc:
        err = f"schema correction round raised: {exc}"
    return "".join(parts).strip(), in_tok, out_tok, err


# One file PER running subagent (not a shared dict file): 6+ parallel
# subagents update their live status on every tool call, and a shared
# read-modify-write would race-drop each other's entries (the panel would
# flicker subagents in and out). One owner per file → no races. The normal
# finally-cleanup removes a subagent's file when it finishes or errors.
_RUNNING_DIR = Path.home() / ".delfin" / "subagent_running"


def _format_action(name: str, tool_input) -> str:
    """A short 'tool: target' line for the live drill-down panel."""
    # Strip the MCP namespace so the panel shows `bash`, not the noisy
    # `mcp__kit-coding__bash` (user 2026-06-26: "die mcp_ Sachen … zu viel").
    name = (name or "").split("__")[-1]
    try:
        args = json.loads(tool_input) if isinstance(tool_input, str) else (tool_input or {})
    except Exception:
        args = {}
    if not isinstance(args, dict):
        args = {}
    if name in ("write_file", "edit_file", "multi_edit", "read_file", "Write", "Edit", "Read"):
        tgt = str(args.get("path") or args.get("file_path") or "")
        tgt = tgt.rsplit("/", 1)[-1] if tgt else ""
        return f"{name} {tgt}".strip()
    if name in ("bash", "Bash"):
        cmd = str(args.get("command") or "").strip().replace("\n", " ")
        return f"bash: {cmd[:60]}"
    if name in ("search_docs", "grep", "Grep", "glob", "Glob"):
        q = str(args.get("query") or args.get("pattern") or "")
        return f"{name} {q[:40]}".strip()
    return name


def _owner_stamp() -> dict:
    """Identify the process that owns a live registry entry.

    The entry file is removed in a ``finally``, which covers every way the
    RUN can end and none of the ways the PROCESS can end. Without an owner,
    a killed dashboard left the file behind and the subagent showed as
    working forever -- in the panel, and to a parent agent polling for a
    report that could no longer arrive.

    The start ticks come along because a pid alone is not an identity: the
    number is recycled, and a recycled pid would keep a dead entry alive.
    Same stamp the background-job registry already uses."""
    from .bash_jobs import _proc_start_ticks
    pid = os.getpid()
    return {"owner_pid": pid, "owner_start_ticks": _proc_start_ticks(pid)}


def _entry_owner_alive(entry: dict) -> bool:
    """Whether the process that wrote this entry is still around.

    An entry with no owner recorded was written before the field existed;
    an unknown owner is not evidence of a dead one, so it counts as alive.
    Those disappear on their own as soon as the run that owns them ends."""
    try:
        pid = int((entry or {}).get("owner_pid") or 0)
    except (TypeError, ValueError):
        return True
    if not pid:
        return True
    from .bash_jobs import _pid_alive
    return _pid_alive(pid, (entry or {}).get("owner_start_ticks"))


def _entry_owned_by_us(entry: dict) -> bool:
    """Whether THIS process wrote the marker.

    ``_entry_owner_alive`` asks a different question -- is the writer still
    around -- and the pending drain used only that one. So a second DELFIN
    session on the same machine, sharing ``~/.delfin/subagent_pending``,
    passed the liveness check on somebody else's marker and claimed it.
    Measured: a marker stamped with another live pid was drained here, its
    report handed to a parent that never started the delegate, and the
    marker unlinked -- so the exactly-once contract held perfectly and
    delivered to the wrong session, while the right one waited for a
    report that could no longer arrive.

    A marker with NO owner recorded predates the field. Nobody could match
    it, and leaving it to the TTL reaper would strand a report that used to
    be delivered, so those stay claimable by any live session -- exactly
    the old behaviour, kept on purpose.
    """
    try:
        pid = int((entry or {}).get("owner_pid") or 0)
    except (TypeError, ValueError):
        return True
    if not pid:
        return True
    if pid != os.getpid():
        return False
    # Same pid is not yet the same process: after a reboot or a pid wrap
    # the number gets reissued. The start-time stamp is what tells the two
    # apart, and it is already recorded.
    theirs = (entry or {}).get("owner_start_ticks")
    ours = _owner_stamp().get("owner_start_ticks")
    if theirs and ours and theirs != ours:
        return False
    return True


def _running_update(sa_id: str, entry: dict | None) -> None:
    """Maintain the live per-subagent status file (entry=None removes).

    File-based so the dashboard can render a live drill-down
    (name · task · steps · status) without sharing memory with the worker
    thread."""
    try:
        from .state_paths import ensure_dir, write_text
        ensure_dir(_RUNNING_DIR)
        f = _RUNNING_DIR / f"{sa_id}.json"
        if entry is None:
            try:
                f.unlink()
            except FileNotFoundError:
                pass
        else:
            if not f.exists():
                # Starting an entry is the one moment guaranteed to happen
                # and cheap enough to carry the cleanup of what an earlier
                # crash left behind.
                reap_dead_running()
            write_text(f, json.dumps({**entry, **_owner_stamp()}))
    except Exception:
        pass


def reserve_running(sa_id: str, *, subagent_type: str = "",
                    description: str = "") -> None:
    """Announce a subagent id BEFORE the work behind it starts.

    A background run hands its id to the parent immediately, but the live
    entry used to be written inside the run itself. Anything that failed in
    between -- the runner raising on the way in, a runner that stores
    nothing -- left that id resolving to "no such id", which the model
    reads as having invented the id rather than as a run that crashed.
    Reserving the entry with the id keeps the two answers distinct.

    Reserving is also what marks the run as one the parent is OWED a
    report for: it happens only on the background path, where the tool
    result carries an id instead of the work. ``_note_pending_report``
    records that debt so the finished report can be pushed into a turn
    rather than waiting for a poll that may never come."""
    _running_update(sa_id, {
        "type": subagent_type,
        "description": (description or "")[:120],
        "started_at": time.time(),
        "actions": [],
        "last_action": "",
        "transcript": [],
        "reserved": True,
    })
    _note_pending_report(sa_id, subagent_type=subagent_type,
                         description=description)


def mark_running_died(sa_id: str, error: str = "") -> None:
    """Record that a subagent ended without leaving a report.

    Not a removal: the id was handed to a parent, so it has to keep
    answering -- with what happened, so the parent can start the work
    again instead of polling for a report that can no longer arrive."""
    entry = dict(read_running(include_dead=True).get(sa_id) or {})
    entry.pop("dead", None)
    entry.setdefault("type", "")
    entry.setdefault("description", "")
    entry.setdefault("started_at", time.time())
    entry["died"] = True
    entry["error"] = (error or "").strip()[:400]
    _running_update(sa_id, entry)


def read_running(*, include_dead: bool = False) -> dict:
    """Live registry: {id: {type, description, started_at, actions, last_action}}.

    Entries whose owning process is gone are left out: they are leftovers
    from a killed or crashed session, not work in flight. So are entries
    whose run ended without a report (``died``) -- the process is fine,
    the run is not. Pass ``include_dead=True`` to see both, marked with
    ``dead: True``."""
    out: dict = {}
    try:
        for f in _RUNNING_DIR.glob("*.json"):
            try:
                entry = json.loads(f.read_text(encoding="utf-8"))
            except Exception:
                continue
            if entry.get("died"):
                if include_dead:
                    out[f.stem] = {**entry, "dead": True}
            elif _entry_owner_alive(entry):
                out[f.stem] = entry
            elif include_dead:
                out[f.stem] = {**entry, "dead": True}
    except Exception:
        pass
    return out


# How long a "this run died" record is kept answering for its id. Long
# enough that the parent that was given the id still gets the truthful
# answer; short enough that the records cannot pile up unbounded, since
# their owning process is alive and nothing else would ever remove them.
_DIED_ENTRY_TTL_S = 24 * 3600


def reap_dead_running() -> list[str]:
    """Delete registry entries whose owning process is gone. Returns the ids.

    Without this the files accumulate: nothing else ever removes an entry
    whose ``finally`` did not run, so every crash leaves one behind for
    good. Aged-out ``died`` records go the same way."""
    removed: list[str] = []
    now = time.time()
    try:
        for f in sorted(_RUNNING_DIR.glob("*.json")):
            try:
                entry = json.loads(f.read_text(encoding="utf-8"))
            except Exception:
                continue
            if entry.get("died"):
                try:
                    age = now - float(entry.get("started_at") or 0.0)
                except (TypeError, ValueError):
                    age = 0.0
                if age <= _DIED_ENTRY_TTL_S:
                    continue
            elif _entry_owner_alive(entry):
                continue
            try:
                f.unlink()
                removed.append(f.stem)
            except OSError:
                continue
    except Exception:
        pass
    return removed


# ---------------------------------------------------------------------------
# Background reports the parent has not been given yet (exactly-once)
# ---------------------------------------------------------------------------
#
# A backgrounded delegation returns an id instead of a report, and nothing
# ever pushed the finished report anywhere: the only route back was the
# parent spending a whole tool round on ``subagent_result(sa_id)``, against
# round budgets as small as ten. A parent that ended its turn first -- the
# normal "I started it, I'll wait" ending -- never took another turn, so the
# report was not late, it was gone.
#
# One file per outstanding report, claimed by exactly one drain. Same
# exactly-once shape as the bash-job registry, for the same reason: a
# completion announced twice is its own defect, and the claim has to
# survive a restart. The claim IS the unlink -- on POSIX only one caller
# can remove a given file -- so two drains racing cannot both report it.
_PENDING_DIR = Path.home() / ".delfin" / "subagent_pending"
# A marker whose owning process is gone belongs to a session that can no
# longer be told anything; kept only long enough to be reaped, like the
# died-run records above.
_PENDING_TTL_S = 24 * 3600
# How many finished reports one drain hands over. The rest stay claimed by
# nobody and arrive on the next drain -- dropping them here would destroy
# exactly what this exists to deliver.
_PENDING_DRAIN_LIMIT = 3


def _pending_path(sa_id: str) -> Path:
    return _PENDING_DIR / f"{(sa_id or '').strip()}.json"


def _atomic_write_json(path: Path, data: dict) -> None:
    """Write a marker so no reader can ever see it half-written.

    ``bash_jobs._atomic_write_json`` is the pattern this codebase already
    uses for exactly this — temp file beside the target, then
    ``os.replace``, which is atomic within a directory."""
    from .bash_jobs import _atomic_write_json as _write
    _write(path, data)


def _note_pending_report(sa_id: str, *, subagent_type: str = "",
                         description: str = "") -> None:
    """Record that this process owes its parent agent a report for ``sa_id``.

    Best-effort: a delegation must never fail because bookkeeping could
    not be written.

    The owner stamp is what keeps two concurrent sessions on one machine
    from draining each other's reports — and for a long time that sentence
    was true of nothing. The drain matched on ``_entry_owner_alive``, which
    asks whether the WRITER is still around, a question another live
    session passes. It matches the stamp against this process now, in
    ``_entry_owned_by_us``.

    Written atomically. ``write_text`` opens with ``"w"``, which truncates
    before it writes, so a marker existed as an empty file for as long as
    the write took — and the reaper that runs on every reservation reads
    exactly this directory. Measured over 24 concurrent reservations, two
    live markers per trial were destroyed at birth: the parent was never
    told the sub-agent had finished, and the report was not late but
    gone."""
    sa_id = (sa_id or "").strip()
    if not sa_id:
        return
    try:
        from .state_paths import ensure_dir, write_text
        ensure_dir(_PENDING_DIR)
        # Starting one is the moment to clear what an earlier crash left
        # behind — the same place the live registry does its reaping.
        reap_pending_reports()
        # Atomic AND owner-only: mkstemp creates the temporary file 0600
        # and os.replace carries that mode over, so this satisfies both the
        # torn-read fix and the state-file permission rule.
        _atomic_write_json(_pending_path(sa_id), {
            "sa_id": sa_id,
            "type": subagent_type or "",
            "description": (description or "")[:120],
            "started_at": time.time(),
            **_owner_stamp(),
        })
    except Exception:
        pass


def _claim_pending_report(sa_id: str) -> bool:
    """Take ownership of an outstanding report. True iff THIS call took it.

    Removing the file is the claim, so the exactly-once contract needs no
    lock and no flag: whoever the unlink succeeds for is the one caller
    that may announce this completion."""
    try:
        _pending_path(sa_id).unlink()
        return True
    except (FileNotFoundError, NotADirectoryError):
        return False
    except Exception:
        return False


def reap_pending_reports() -> list[str]:
    """Drop markers nobody can be told about any more. Returns the ids.

    A marker outlives its process only when that process died; nothing
    else would ever remove it, so without this they accumulate for good.

    A record that cannot be read is skipped, exactly as the sibling reaper
    over the running registry skips one. Substituting an empty record
    instead read as "owner pid 0" — which counts as alive — with a
    ``started_at`` of 0.0, so its age came out around 1.8e9 seconds, past
    every TTL, and the marker was deleted. The unreadable case is a marker
    being written this instant, so the one it deleted was always a live
    one."""
    removed: list[str] = []
    now = time.time()
    try:
        for f in sorted(_PENDING_DIR.glob("*.json")):
            try:
                rec = json.loads(f.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not isinstance(rec, dict):
                continue
            if _entry_owner_alive(rec):
                try:
                    age = now - float(rec.get("started_at") or 0.0)
                except (TypeError, ValueError):
                    age = 0.0
                if age <= _PENDING_TTL_S:
                    continue
            try:
                f.unlink()
                removed.append(f.stem)
            except OSError:
                continue
    except Exception:
        pass
    return removed


def drain_finished_subagents(limit: int = _PENDING_DRAIN_LIMIT) -> list[dict]:
    """Background delegations that ENDED since the last drain — exactly once.

    Returns the collected result of each run this process backgrounded and
    has not yet handed to its parent, in the shape ``subagent_result``
    would have returned: ``sa_id``, ``subagent_type``, ``description``,
    ``status`` (``finished`` / ``died``), ``final_text``, ``error``, the
    verification verdict, and the worktree summary when the run had one.

    A run still in flight is left alone, so a marker is consumed only when
    there is something to say. Never raises: this runs on the parent's turn
    and must not be able to end it."""
    out: list[dict] = []
    try:
        files = sorted(_PENDING_DIR.glob("*.json"))
    except Exception:
        return out
    if not files:
        return out
    try:
        live = read_running()
        dead = read_running(include_dead=True)
    except Exception:
        live, dead = {}, {}
    for f in files:
        if len(out) >= max(1, int(limit or 1)):
            break
        sa_id = f.stem
        try:
            rec = json.loads(f.read_text(encoding="utf-8"))
        except Exception:
            continue                    # unreadable this instant, not gone
        if not isinstance(rec, dict) or not _entry_owner_alive(rec):
            continue                    # a dead session's leftover; reaped
        if not _entry_owned_by_us(rec):
            continue                    # a LIVE session's report, and theirs
                                        # to deliver -- claiming it here is
                                        # how the right parent stops getting
                                        # one at all
        # Terminality is decided BEFORE the claim, and without
        # get_subagent_result -- that call acknowledges the report itself
        # (an explicit collection is a delivery), so asking it first would
        # consume the marker and leave this drain with nothing to hand on.
        # Same precedence it uses: a live entry outranks a stored report,
        # so a resumed run is not reported as finished mid-flight.
        if sa_id in live:
            continue                    # still in flight; say nothing
        if not (load_subagent_session(sa_id) or dead.get(sa_id)):
            continue                    # nothing recorded yet; not an ending
        if not _claim_pending_report(sa_id):
            continue                    # a concurrent drain announced it
        try:
            result = get_subagent_result(sa_id)
        except Exception:
            result = {}
        if str((result or {}).get("status") or "") in ("finished", "died"):
            out.append(result)
        else:
            # The run went back to running (a resume) or the store could
            # not be read. Put the debt back: a claim that hands nothing on
            # is how an event gets destroyed rather than delivered.
            _note_pending_report(sa_id, subagent_type=str(rec.get("type") or ""),
                                 description=str(rec.get("description") or ""))
    return out


# Finished-subagent sessions (resume-by-id): each run
# persists its conversation so a later ``resume_id`` call can continue
# the same subagent with its context intact.
_SESSIONS_DIR = Path.home() / ".delfin" / "subagent_sessions"
_SESSIONS_LIST_LIMIT = 20


def _trim_for_store(text: str, limit: int = 2000) -> str:
    """Head+tail trim for stored message content (sliding-window style)."""
    if not isinstance(text, str) or len(text) <= limit:
        return text
    return (
        text[:1200]
        + f"\n... [trimmed for subagent session store, "
        + f"{len(text) - 1600} chars dropped] ...\n"
        + text[-400:]
    )


_MAX_STORED_INTERACTIONS = 60

# The final report is stored separately from the trimmed conversation.
# The trim (head 1200 + tail 400) is right for the resume recap it was
# written for, and wrong for the report: the session store is the ONLY
# route by which a BACKGROUND report reaches its parent, so a long one
# came back as a stub with a trim marker where the findings had been.
# Bounded, because a runaway report must not be able to fill the disk --
# but at a size a report actually fits in, and it says so when it cuts.
_MAX_STORED_REPORT_CHARS = 200_000


def collectable_sa_id(*, reserved: str, resume_from: str) -> str:
    """The id a finished run will be stored under.

    Backgrounding reserves an id up front and tells the parent to collect
    with it, but a resumed run keeps the id it is resuming -- correctly,
    so one session accumulates under one id. The two disagreed, and the
    parent was handed the reserved one: it polled an id that never existed
    and was told "unknown" forever, for work that had finished."""
    return (resume_from or "").strip() or (reserved or "").strip()


def _save_subagent_session(
    sa_id: str,
    *,
    subagent_type: str,
    description: str,
    messages: list[dict],
    interactions: list[dict],
    error: str = "",
    worktree: dict | None = None,
) -> None:
    """Persist a finished subagent conversation for later resumption.

    Stores the logical user/assistant conversation AND the tool
    interactions WITH their (trimmed) outputs — so a resumed subagent
    sees what it actually read, not just its own conclusions. Best-effort,
    never raises.

    ``worktree`` is the framework's own account of the isolated tree the
    run used. A FOREGROUND run returns it in the payload; a background one
    returns an id, so without storing it here the account of the tree
    reached nobody — including the one case that has to be said out loud,
    a tree kept because the run left processes in it.
    """
    try:
        record = {
            "sa_id": sa_id,
            "subagent_type": subagent_type,
            "description": (description or "")[:200],
            "finished_at": time.time(),
            "messages": [
                {"role": m.get("role", ""),
                 "content": _trim_for_store(str(m.get("content", "")))}
                for m in messages
                if m.get("role") in ("user", "assistant")
            ],
            # Tool calls + the outputs they returned (the fidelity that the
            # old text-only replay dropped). Keep the most recent ones.
            "interactions": [
                {"name": it.get("name"),
                 "input": str(it.get("input", ""))[:200],
                 "output": _trim_for_store(str(it.get("output", "")))}
                for it in (interactions or [])[-_MAX_STORED_INTERACTIONS:]
            ],
            "error": error or "",
            "worktree": dict(worktree or {}),
        }
        # The report, whole. Taken before the trim above touches it.
        for m in reversed(messages or []):
            if m.get("role") == "assistant" and str(m.get("content", "")):
                full = str(m["content"])
                if len(full) > _MAX_STORED_REPORT_CHARS:
                    full = (full[:_MAX_STORED_REPORT_CHARS]
                            + "\n\n[report truncated at "
                            + f"{_MAX_STORED_REPORT_CHARS} characters]")
                record["final_report"] = full
                break
        # Owner-only from creation. This record is the densest single file
        # the agent writes: up to 60 tool outputs at 2000 chars each plus
        # the delegate's whole report, and it had no mode of its own at all
        # — all 108 on disk were observed group-readable.
        from .state_paths import ensure_dir, write_text
        ensure_dir(_SESSIONS_DIR)
        write_text(_SESSIONS_DIR / f"{sa_id}.json",
                   json.dumps(record, ensure_ascii=False))
    except Exception:
        pass


def _render_resume_recap(prior: dict) -> str:
    """Build a faithful, model-robust recap of a prior subagent run.

    Includes the earlier conversation AND the tool calls with their actual
    (trimmed) outputs, rendered as one readable context block — this is
    what closes the resume-fidelity gap without juggling tool_call_id
    protocol details across backends.
    """
    msgs = prior.get("messages") or []
    inter = prior.get("interactions") or []
    parts = ["[Resuming this subagent — your earlier context follows. "
             "Trust it as your own prior work.]"]
    if msgs:
        parts.append("\n## Earlier conversation")
        for m in msgs:
            role = str(m.get("role", "")).upper()
            content = _trim_for_store(str(m.get("content", "")), 800)
            if content.strip():
                parts.append(f"{role}: {content}")
    if inter:
        parts.append("\n## Tools you already ran and what they returned")
        for it in inter[-_MAX_STORED_INTERACTIONS:]:
            name = it.get("name", "?")
            inp = str(it.get("input", ""))[:200]
            out = _trim_for_store(str(it.get("output", "")), 800)
            parts.append(f"- {name}({inp}) -> {out}")
    parts.append("\nContinue from here with the new request below.")
    return "\n".join(parts)


def load_subagent_session(sa_id: str) -> dict | None:
    """Load one stored subagent session; ``None`` when unknown/corrupt."""
    try:
        path = _SESSIONS_DIR / f"{(sa_id or '').strip()}.json"
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def list_finished(last_n: int = _SESSIONS_LIST_LIMIT) -> list[dict]:
    """Most recently finished subagent sessions, newest first.

    The limit applies at READ time only — older session files stay on
    disk untouched.
    """
    try:
        files = sorted(_SESSIONS_DIR.glob("*.json"),
                       key=lambda p: p.stat().st_mtime, reverse=True)
    except Exception:
        return []
    out: list[dict] = []
    for p in files[: max(1, int(last_n or _SESSIONS_LIST_LIMIT))]:
        try:
            rec = json.loads(p.read_text(encoding="utf-8"))
            out.append({
                "sa_id": rec.get("sa_id", p.stem),
                "subagent_type": rec.get("subagent_type", ""),
                "description": rec.get("description", ""),
                "finished_at": rec.get("finished_at", 0),
                "error": rec.get("error", ""),
            })
        except Exception:
            continue
    return out


def _write_telemetry(record: dict) -> None:
    """Append a one-line JSON record to ``~/.delfin/subagent_telemetry.jsonl``.

    Best-effort: silently skip if the home dir isn't writable. Trims to
    the most recent ``_TELEMETRY_MAX_LINES`` entries to keep the file
    cheap to read for ``/agents stats``.
    """
    try:
        from .state_paths import ensure_dir, open_append, secure_file
        from .state_paths import write_text as _write_secure
        path = _TELEMETRY_PATH
        ensure_dir(path.parent)
        with open_append(path) as f:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")
        secure_file(path)
        # Best-effort trim — only every 50 writes to avoid I/O thrash
        try:
            stat = path.stat()
            if stat.st_size > 1_000_000:  # ~5k records of 200 bytes
                lines = path.read_text(encoding="utf-8").splitlines()
                if len(lines) > _TELEMETRY_MAX_LINES:
                    tail = lines[-_TELEMETRY_MAX_LINES:]
                    _write_secure(path, "\n".join(tail) + "\n")
        except Exception:
            pass
    except Exception:
        pass


def read_telemetry(path: Path | None = None, *, last_n: int | None = None) -> list[dict]:
    """Load telemetry records as a list of dicts. Newest last."""
    p = path or _TELEMETRY_PATH
    if not p.exists():
        return []
    out: list[dict] = []
    try:
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    except OSError:
        return []
    if last_n is not None and last_n > 0:
        return out[-last_n:]
    return out


def summarize_telemetry(records: list[dict]) -> dict:
    """Group telemetry records by subagent_type and return per-type
    aggregates (count, total elapsed_s, total in/out tokens, error rate)."""
    by_type: dict[str, dict] = {}
    for r in records:
        t = r.get("subagent_type") or "?"
        b = by_type.setdefault(t, {
            "count": 0, "elapsed_s_total": 0.0,
            "input_tokens_total": 0, "output_tokens_total": 0,
            "errors": 0, "truncated": 0,
        })
        b["count"] += 1
        b["elapsed_s_total"] += float(r.get("elapsed_s") or 0)
        b["input_tokens_total"] += int(r.get("input_tokens") or 0)
        b["output_tokens_total"] += int(r.get("output_tokens") or 0)
        if r.get("error"):
            b["errors"] += 1
        if r.get("truncated"):
            b["truncated"] += 1
    return by_type


@dataclass(frozen=True)
class SubagentPreset:
    name: str
    description: str
    system_prompt: str
    mode: str = "default"          # "plan" / "default" / "acceptEdits" / "bypassPermissions"


_BUILTIN_PRESETS: dict[str, SubagentPreset] = {
    "explore": SubagentPreset(
        name="explore",
        description=(
            "Read-only research agent for locating code. Use for "
            "'where is X defined?', 'find all callers of Y', 'how is "
            "Z structured?'."
        ),
        system_prompt=(
            "You are a read-only Explore sub-agent. Your job is to "
            "find code, summarise findings, and return concise "
            "results. You may use read_file, grep_file, list_files, "
            "web_search, web_fetch, notebook_read. You may NOT edit, "
            "write, or run bash. Report back in 2-5 short bullets, "
            "each with a path:line reference where possible. Do not "
            "speculate; if a thing is not found, say so."
        ),
        mode="plan",
    ),
    "plan": SubagentPreset(
        name="plan",
        description=(
            "Planning agent that drafts an implementation strategy "
            "without touching code."
        ),
        system_prompt=(
            "You are a Plan sub-agent. Read the relevant files, then "
            "draft a step-by-step implementation plan. Do NOT edit, "
            "write, or run bash. Return a numbered Markdown plan "
            "with: (1) files to change, (2) order of changes, "
            "(3) anything risky or irreversible. Be terse."
        ),
        mode="plan",
    ),
    "code-reviewer": SubagentPreset(
        name="code-reviewer",
        description=(
            "Read-only reviewer that audits a diff or file against "
            "common-sense correctness, security, and style criteria."
        ),
        system_prompt=(
            "You are a code-review sub-agent. Read the targeted "
            "files / diff and report concrete issues: bugs, security, "
            "missing validation, dead code. Do NOT edit. Format your "
            "answer as a checklist; lead with the highest-severity "
            "issue. If nothing concerning is found, say so explicitly."
        ),
        mode="plan",
    ),
    "general-purpose": SubagentPreset(
        name="general-purpose",
        description=(
            "Full-tool sub-agent that can read, edit, run bash, etc. "
            "Use sparingly — most tasks fit a narrower preset."
        ),
        system_prompt=(
            "You are a general-purpose sub-agent. Complete the "
            "assigned task end-to-end and return a brief summary of "
            "what you did and the final state."
        ),
        mode="default",
    ),
}


def _md_preset_search_dirs() -> list[Path]:
    """Directories scanned for ``*_subagent.md`` user-extensible presets."""
    return [
        Path(__file__).resolve().parent / "pack" / "agents",
        Path.home() / ".delfin" / "subagents",
    ]


def _load_md_presets() -> dict[str, SubagentPreset]:
    """Discover ``*_subagent.md`` presets with YAML frontmatter.

    Frontmatter fields accepted: ``name`` (required, kebab-case),
    ``description``, ``mode``. The body of the markdown file becomes
    the sub-agent's system_prompt. Built-in presets keep precedence —
    user-supplied files can extend but not silently override unless the
    user explicitly drops a file under ``~/.delfin/subagents/``.

    The frontmatter parser is reused from ``delfin.agent.skills``.
    """
    try:
        from .skills import _parse_frontmatter
    except Exception:
        return {}
    discovered: dict[str, SubagentPreset] = {}
    for d in _md_preset_search_dirs():
        if not d.is_dir():
            continue
        try:
            paths = sorted(d.glob("*_subagent.md"))
        except OSError:
            continue
        for p in paths:
            try:
                text = p.read_text(encoding="utf-8")
            except OSError:
                continue
            meta, body = _parse_frontmatter(text)
            name = (meta.get("name") or p.stem.replace("_subagent", "")).strip()
            if not name:
                continue
            description = (meta.get("description") or "").strip()
            mode = (meta.get("mode") or "default").strip() or "default"
            system_prompt = body.strip() or description or name
            discovered[name] = SubagentPreset(
                name=name,
                description=description or f"Custom subagent '{name}'",
                system_prompt=system_prompt,
                mode=mode,
            )
    return discovered


def _build_preset_registry() -> dict[str, SubagentPreset]:
    """Compose the live preset registry: built-ins first, then md overrides.

    MD files in ``~/.delfin/subagents/`` may override built-ins; MD files
    in the packaged ``delfin/agent/pack/agents/`` directory cannot —
    they only contribute new types. This keeps user-local extensions
    powerful while preventing accidental pack-shipped overrides.
    """
    registry = dict(_BUILTIN_PRESETS)
    pack_dir = Path(__file__).resolve().parent / "pack" / "agents"
    user_dir = Path.home() / ".delfin" / "subagents"
    for d, allow_override in ((pack_dir, False), (user_dir, True)):
        if not d.is_dir():
            continue
        try:
            paths = sorted(d.glob("*_subagent.md"))
        except OSError:
            continue
        try:
            from .skills import _parse_frontmatter
        except Exception:
            return registry
        for p in paths:
            try:
                text = p.read_text(encoding="utf-8")
            except OSError:
                continue
            meta, body = _parse_frontmatter(text)
            name = (meta.get("name") or p.stem.replace("_subagent", "")).strip()
            if not name:
                continue
            if name in _BUILTIN_PRESETS and not allow_override:
                continue
            description = (meta.get("description") or "").strip()
            mode = (meta.get("mode") or "default").strip() or "default"
            system_prompt = body.strip() or description or name
            registry[name] = SubagentPreset(
                name=name,
                description=description or f"Custom subagent '{name}'",
                system_prompt=system_prompt,
                mode=mode,
            )
    return registry


# Composed at import time: built-ins + any user/pack md files.
SUBAGENT_PRESETS: dict[str, SubagentPreset] = _build_preset_registry()


def reload_subagent_presets() -> dict[str, SubagentPreset]:
    """Re-scan markdown preset directories. Useful in tests + after the
    user drops a new file into ``~/.delfin/subagents/``."""
    SUBAGENT_PRESETS.clear()
    SUBAGENT_PRESETS.update(_build_preset_registry())
    return SUBAGENT_PRESETS


def subagent_type_names() -> list[str]:
    """Return current preset names — used by api_client to build the
    runtime tool-schema enum."""
    return list(SUBAGENT_PRESETS.keys())


def list_subagents() -> list[dict]:
    return [
        {
            "name": p.name,
            "subagent_type": p.name,
            "description": p.description,
            "mode": p.mode,
        }
        for p in SUBAGENT_PRESETS.values()
    ]


def mark_delegate_text(text: str) -> str:
    """Mark a delegate's own prose as untrusted before the parent reads it.

    A sub-agent reads whatever its task points it at — a checked-out
    repository's README, a fetched page, a tool result from an MCP server
    someone else configured. Its report is a MODEL's summary of that
    material, so any instruction inside the material can reach the parent
    through it. The same boundary is already drawn for web_search,
    web_fetch and MCP results in ``api_client._wrap_untrusted``; the
    delegation path was the one consumer that never called it, and it is
    the one whose text arrives with a colleague's authority rather than a
    stranger's.

    Only the PROSE is marked. The envelope around it — the id, the token
    counts, the tool-call names, the verification verdict — is built by
    this file from its own records, and marking that too would say the
    harness does not trust itself.

    The remaining surface, named rather than implied: ``structured_output``
    is not marked. Its shape is constrained by a schema the PARENT wrote,
    which is narrower, and stringifying it would destroy the one thing it
    exists to provide. A free-text field inside such a schema is still a
    route, and this is where that is written down.
    """
    body = (text or "").strip()
    if not body:
        return text or ""
    try:
        from .api_client import _wrap_untrusted
    except Exception:
        return body
    return _wrap_untrusted(body)


@dataclass
class SubagentResult:
    subagent_type: str
    description: str
    final_text: str
    tool_calls: list[dict] = field(default_factory=list)
    elapsed_s: float = 0.0
    input_tokens: int = 0
    output_tokens: int = 0
    truncated: bool = False
    error: str = ""
    worktree: dict = field(default_factory=dict)
    sa_id: str = ""
    model: str = ""          # model the run actually used ("" = unknown)
    model_tier: str = ""     # "cheap" when routed to the cheap tier, else "parent"
    # Workspace the child ran in. Not part of the payload — used as the repo
    # root when cross-checking file citations in the report (see
    # ``verify_subagent_report``).
    workspace: str = ""
    # Schema-validated returns (only meaningful when run_subagent was given
    # an output_schema): the parsed+validated dict, or None with
    # schema_error explaining why validation failed after the one
    # correction round. final_text stays available either way.
    structured_output: Optional[dict] = None
    schema_error: str = ""

    def to_payload(self) -> dict:
        payload = {
            "subagent_type": self.subagent_type,
            "description": self.description,
            "sa_id": self.sa_id,
            "model": self.model,
            "result": self.final_text,
            "tool_calls": [
                {"name": tc.get("name"), "input": str(tc.get("input"))[:200]}
                for tc in self.tool_calls
            ],
            "elapsed_s": round(self.elapsed_s, 2),
            "input_tokens": self.input_tokens,
            "output_tokens": self.output_tokens,
            "truncated": self.truncated,
            "error": self.error,
            "worktree": self.worktree or {},
        }
        if self.structured_output is not None:
            payload["structured_output"] = self.structured_output
        if self.schema_error:
            payload["schema_error"] = self.schema_error
        # Cross-check the report against this run's own tool trace BEFORE the
        # parent ever reads it. The full trace (with tool outputs) is passed
        # explicitly — the payload copy above carries names/inputs only.
        # Defined below; never raises, and leaves ``result`` untouched.
        payload = attach_verification(
            payload,
            tool_calls=self.tool_calls,
            repo_root=self.workspace or None,
        )
        # Marked LAST, and the order is load-bearing. The verifier scans
        # this field for claims and matches them against the run's own
        # tool trace; marking first would have it reading the wrapper's
        # own words as part of the delegate's report. So: verify the text
        # the delegate wrote, then hand the parent a marked copy of it.
        payload["result"] = mark_delegate_text(payload.get("result", ""))
        return payload


# ---------------------------------------------------------------------------
# Delegate-report verification — claims vs. the delegate's own trace
# ---------------------------------------------------------------------------
#
# Everything below is pure: dict in, dict out, no I/O, no model call, and no
# exception ever leaves these functions. The claim scanners for code
# locations, path citations, physical quantities and functional assertions
# live in ``verify_guard`` (the parent already runs them over its OWN
# answers) and are REUSED here — only the two claim classes the parent never
# needed for itself are new: test-result claims and completion claims.
#
# Evidence is whatever the delegate's own run recorded: the tools it called,
# the inputs it passed them, the commands it executed, and — when available —
# what those calls returned. See the module docstring for the limits.

_VERIFY_MAX_TEXT = 20_000        # report chars scanned (bounded work)
_VERIFY_MAX_PER_KIND = 4         # flags taken from any single scanner
_VERIFY_MAX_UNSUPPORTED = 8
_VERIFY_MAX_SUPPORTED = 6
_VERIFY_MAX_NOTICE_ITEMS = 4     # unsupported claims named in the notice line
_VERIFY_CLAIM_CHARS = 120
_VERIFY_MAX_OUTPUT_CHARS = 20_000
_VERIFY_MAX_CALLS = 200

_VERIFY_LIMITS_NOTE = (
    "consistency check only: it compares the report against this run's own "
    "tool trace, never against reality — a false conclusion drawn from real "
    "tool calls is not detectable here"
)

# Tools whose call means the delegate MUTATED something.
_WRITE_TOOL_NAMES = frozenset({
    "write_file", "edit_file", "multi_edit", "apply_patch", "notebook_edit",
    "write", "edit", "patch", "create_file", "undo_changes",
})
# Tools whose call means the delegate LOOKED AT the file it names. Mirrors
# the parent's own observation set; a tool outside it (bash above all) names
# its files only inside a command string, where naming and reading are
# indistinguishable.
_READ_TOOL_NAMES = frozenset({
    "read_file", "grep_file", "list_files", "notebook_read", "view_image",
    "read_document", "read", "read_section",
})
# Tools whose call means the delegate RAN tests (name-based; command-based
# detection reuses verify_guard's test-command pattern).
_TEST_TOOL_NAMES = frozenset({"run_tests", "pytest"})

# Tool-input keys that carry a file target.
_PATH_INPUT_KEYS = (
    "path", "file_path", "filepath", "file", "paths", "files",
    "notebook_path", "target", "directory", "dir",
)

# Loose path-ish token used to harvest file names out of command strings and
# git --stat output. Deliberately generous: over-collecting only makes the
# cross-check more lenient (fewer false alarms), never stricter.
_PATH_TOKEN_RE = re.compile(
    r"(?:[\w\-.]+/)+[\w\-.]+|\b[\w\-]+\.[A-Za-z][A-Za-z0-9]{0,5}\b")

# Shell shapes that write something (heredocs, redirects, in-place edits).
_BASH_WRITE_RE = re.compile(
    r"(?:>>?\s*[\w./\-]|\btee\b|\bsed\s+-i\b|\bpatch\s+-p\d|\bgit\s+apply\b|"
    r"\bgit\s+(?:commit|checkout|restore|revert)\b|\bmv\b|\bcp\b|\btouch\b|"
    r"\bmkdir\b|\brm\b|\bcat\s+>|\btruncate\b)")

# A redirect whose target is a discard sink changes nothing. The write
# pattern above matches any redirect, so `pytest -q > /dev/null` -- a
# command whose entire purpose is to throw the output away -- registered as
# file-writing evidence and backed a claim of having changed something.
_EPHEMERAL_SINK_RE = re.compile(r">>?\s*/dev/(?:null|stdout|stderr)\b")


def _writes_anything(cmd: str) -> bool:
    """Whether a shell command is write evidence at all.

    Only the redirect is discounted, and only when every redirect in the
    command goes to a discard sink: `echo a > /dev/null; echo b > real.txt`
    still writes."""
    if not _BASH_WRITE_RE.search(cmd or ""):
        return False
    stripped = _EPHEMERAL_SINK_RE.sub(" ", cmd or "")
    return bool(_BASH_WRITE_RE.search(stripped))


# Test output that contradicts a success claim. Anchored on the counted
# forms a runner prints, so "0 failed" and prose containing the word are
# not mistaken for a failing run.
_TEST_FAILURE_RE = re.compile(
    r"(?:^|[\s,=])(?!0\s)(\d+)\s+(?:failed|failures?|errors?)\b"
    r"|^FAILED\b|\bFAILED\s+\S+::"
    r"|(?:^|[\s,=])(?!0\s)(\d+)\s+error\b",
    re.IGNORECASE | re.MULTILINE)


def _output_shows_failures(output: str) -> bool:
    """Whether recorded test output contradicts a claim that all passed."""
    return bool(_TEST_FAILURE_RE.search(output or ""))


def _readable_output(out) -> str:
    """The part of a recorded tool result that is the run's own output.

    The shell tool returns a JSON envelope (command, cwd, exit_code,
    elapsed_s, stdout, stderr). Judging a test claim against the envelope
    let the COMMAND STRING stand in for the output it discarded:
    ``pytest -q > /dev/null`` leaves a result that mentions pytest, shows
    no failures, and can even supply a claimed count out of ``elapsed_s``.
    The quoting hid real failures too — ``"stdout": "3 failed, ...`` puts
    a quote where the failure pattern expects a separator.

    Anything that is not such an envelope is returned unchanged.
    """
    text = out if isinstance(out, str) else str(out)
    stripped = text.strip()
    if not (stripped.startswith("{") and '"stdout"' in stripped):
        return text
    try:
        data = json.loads(stripped)
    except (json.JSONDecodeError, ValueError):
        return text          # truncated envelope: the raw text is all there is
    if not isinstance(data, dict) or "stdout" not in data:
        return text
    return "\n".join(
        str(data.get(key) or "") for key in ("stdout", "stderr")
        if str(data.get(key) or "").strip()
    )


# A read that returned an error, or nothing, is not an observation. Paths
# are taken from a call's INPUT, so a refused or empty read still put the
# file in the observed set -- grounding the delegate's own citations and,
# because delegate evidence is merged upward, telling the PARENT it had
# seen a file nobody read.
_FAILED_READ_MARKERS = (
    "(file not found)", "(empty file)", "(no matches)", "(binary file)",
    "file not found", "no such file", "permission denied",
)


def _output_saw_content(output) -> bool:
    """Whether a call's recorded output shows it actually read something.

    ``None`` means no output was recorded — the payload trace carries
    inputs without outputs — and returns True. Missing bookkeeping is never
    treated as evidence of fabrication; concluding failure from an absent
    record would discard every path on that whole path."""
    if output is None:
        return True
    text = output if isinstance(output, str) else str(output)
    low = text.strip().lower()
    if not low:
        return False
    if low.startswith("{") and '"error"' in low:
        return False
    return not any(low.startswith(m) or low == m
                   for m in _FAILED_READ_MARKERS)

# --- new claim classes ------------------------------------------------------
#
# Test-result claims: a report asserting the suite is fine. Every pattern
# needs an explicit success predicate, so "run the tests" or "2 failed" never
# match. Bilingual because reports are written in the user's language.
_TEST_CLAIM_PATTERNS: tuple[tuple[str, "re.Pattern[str]"], ...] = (
    ("suite", re.compile(
        r"(?i)\b(?:all\s+|alle\s+)?(?:\d+\s+)?(?:unit[\s-]?)?tests?\b"
        r"[^\n]{0,32}?\b(?:pass(?:ed|es|ing)?|green|gr(?:ü|ue)n|"
        r"erfolgreich|bestanden|durchgelaufen)\b")),
    ("suite", re.compile(
        r"(?i)\b(?:test[\s-]?suite|testsuite|test\s+run)\b[^\n]{0,32}?"
        r"\b(?:pass(?:ed|es|ing)?|green|gr(?:ü|ue)n|clean|erfolgreich|"
        r"sauber)\b")),
    ("count", re.compile(r"(?i)\b(\d{1,5})\s+(?:tests?\s+)?"
                         r"(?:passed|bestanden)\b")),
    ("nofail", re.compile(
        r"(?i)\b(?:no|zero)\s+(?:test\s+)?(?:failures|failing\s+tests)\b|"
        r"\bkeine\s+(?:test)?(?:fehlschl(?:ä|ae)ge|fehler)\b")),
    ("runner", re.compile(
        r"(?i)\bpytest\b[^\n]{0,32}?\b(?:pass(?:ed|es|ing)?|green|"
        r"gr(?:ü|ue)n|clean|erfolgreich|durch)\b")),
)

# Completion claims, by family. The families differ in what evidence would
# back them, so they are judged by different rules (see _judge_completion).
_COMPLETION_FAMILIES: tuple[tuple[str, "re.Pattern[str]"], ...] = (
    # The ordinary way a delegate says it changed a file — in either
    # language. Measured on the shipped pattern, none of "geschrieben",
    # "erstellt", "eingetragen", "übertragen", "angelegt", "ergänzt",
    # "gespeichert" produced a claim; only "aktualisiert" did. The same
    # gap ran the other way: the list carried "rewrote" and not "wrote",
    # and no "saved" at all, so "I wrote the file" was no claim either.
    # The vocabulary is shared with the write-verb matchers.
    ("mutation", re.compile(
        r"(?i)\b(?:implemented|fixed|patched|refactored|migrated|"
        r"wired(?:\s+up)?|added|created|removed|deleted|renamed|updated|"
        r"rewrote|" + "|".join(_de.ENGLISH_WRITE_PARTICIPLES) + r"|"
        r"implementiert|umgesetzt|behoben|gefixt|"
        r"hinzugef(?:ü|ue)gt|entfernt|umbenannt|aktualisiert|eingebaut|"
        + "|".join(_de.GERMAN_WRITE_PARTICIPLES) + r")\b")),
    ("verification", re.compile(
        r"(?i)\b(?:verified|validated|confirmed|double[\s-]?checked|"
        r"cross[\s-]?checked|bit[\s-]?identical|identical|unchanged|"
        r"verifiziert|validiert|best(?:ä|ae)tigt|(?:ü|ue)berpr(?:ü|ue)ft|"
        r"gepr(?:ü|ue)ft|unver(?:ä|ae)ndert|identisch)\b|"
        r"\bno\s+(?:functional\s+)?(?:changes?|differences?)\b|"
        r"\bkeine\s+(?:funktionalen\s+)?"
        r"(?:(?:ä|ae)nderungen|unterschiede)\b")),
    ("generic", re.compile(
        r"(?i)\b(?:done|completed?|finished|all\s+set|erledigt|fertig|"
        r"abgeschlossen)\b")),
)

# First-person marker: a mutation verb only counts as the delegate's own
# completion claim when the delegate says IT did the thing (or writes it as a
# bare report bullet). "an enum was added to the parameter" is a FINDING
# about the code, not a claim that this run changed anything.
_FIRST_PERSON_RE = re.compile(r"(?i)\b(?:i|we|ich|wir)\b")
_BULLET_PREFIX = "-*•>#0123456789.)( \t"


@dataclass(frozen=True)
class ReportClaim:
    """One claim found in a delegate's report."""

    kind: str            # "test_result" | "completion"
    text: str            # normalized claim excerpt
    subject: str = ""    # count for test claims, "" otherwise
    family: str = ""     # completion family: mutation|verification|generic


def _vg():
    """The shared claim scanners, or None when unavailable.

    Import failures degrade the check to the two local claim classes rather
    than breaking a delegation.
    """
    try:
        from . import verify_guard as _vgmod
        return _vgmod
    except Exception:
        return None


def _scrubbed(text: str) -> str:
    """Blank non-claim regions (code fences, quotes) keeping offsets."""
    vg = _vg()
    if vg is not None:
        try:
            return vg._scrub_keep_offsets(text)
        except Exception:
            pass
    return text


def _hedged(scrubbed: str, start: int, end: int) -> bool:
    """True when the claim's line discloses uncertainty (never a target)."""
    vg = _vg()
    if vg is None:
        return False
    try:
        return bool(vg._is_hedged(scrubbed, start, end))
    except Exception:
        return False


def _sentence_of(scrubbed: str, start: int, end: int) -> str:
    """The sentence containing a match (falls back to its line)."""
    vg = _vg()
    if vg is not None:
        try:
            lo, hi = vg._sentence_span(scrubbed, start, end)
            return scrubbed[lo:hi]
        except Exception:
            pass
    lo = scrubbed.rfind("\n", 0, start) + 1
    hi = scrubbed.find("\n", end)
    return scrubbed[lo:hi if hi != -1 else len(scrubbed)]


def _line_of(scrubbed: str, start: int) -> tuple[str, int]:
    """``(line, offset_of_start_within_line)``."""
    lo = scrubbed.rfind("\n", 0, start) + 1
    hi = scrubbed.find("\n", start)
    return scrubbed[lo:hi if hi != -1 else len(scrubbed)], start - lo


def _is_agentive(line: str, verb_offset: int) -> bool:
    """True when the delegate presents the verb as its OWN action."""
    if _FIRST_PERSON_RE.search(line):
        return True
    return line[:verb_offset].strip(_BULLET_PREFIX).strip() == ""


def _norm_tool_name(name) -> str:
    """Bare tool name (``mcp__ns__bash`` -> ``bash``), lower-cased."""
    try:
        return str(name or "").strip().split("__")[-1].strip().lower()
    except Exception:
        return ""


def _tool_args(raw) -> dict:
    """Decode a tool input (dict or JSON string) to a dict; {} on anything else."""
    if isinstance(raw, dict):
        return raw
    if isinstance(raw, str) and raw.strip().startswith("{"):
        try:
            obj = json.loads(raw)
            return obj if isinstance(obj, dict) else {}
        except (json.JSONDecodeError, ValueError):
            return {}
    return {}


def _paths_in_grep_hits(output: Any) -> list[str]:
    """Files a repo-wide grep demonstrably SHOWED the delegate.

    Shares the parent's extractor rather than a second copy of the same
    regex: two implementations of one rule drift, and a grounding rule
    that means different things on the two sides is exactly the defect
    this closes. Imported inside the function — ``api_client`` imports
    this module, so a top-level import would be a cycle.

    An error result is already rejected by ``_output_saw_content`` before
    this is reached, so there is no second check for it here: a mutation
    test showed the one I first wrote could be deleted without a single
    assertion noticing, which is the definition of a guard that guards
    nothing. Two copies of one rule drift apart; the caller owns it.

    A non-string output does have to be handled — the trace carries dicts
    — and answers nothing.
    """
    text = output if isinstance(output, str) else ""
    if not text:
        return []
    try:
        from .api_client import _paths_in_grep_output
    except Exception:
        return []
    return sorted(_paths_in_grep_output(text))


def _declared_paths(args: dict) -> list[str]:
    """File paths a tool call NAMED in a recognised path argument.

    The narrow half of ``_paths_from_call``: a path here was handed to the
    tool as its target, so the tool acted on that file. Nothing is inferred
    from prose or from a command string.
    """
    out: list[str] = []
    for key in _PATH_INPUT_KEYS:
        val = args.get(key)
        if isinstance(val, str) and val.strip():
            out.append(val.strip())
        elif isinstance(val, (list, tuple)):
            out.extend(str(v).strip() for v in val
                       if isinstance(v, str) and v.strip())
    edits = args.get("edits")
    if isinstance(edits, (list, tuple)):
        for e in edits:
            if not isinstance(e, dict):
                continue
            for key in ("path", "file_path", "file"):
                val = e.get(key)
                if isinstance(val, str) and val.strip():
                    out.append(val.strip())
    return out[:60]


def _paths_from_call(args: dict, raw) -> list[str]:
    """File paths a tool call plausibly touched: the declared ones plus
    path-shaped tokens in the raw input (a bash command names its files
    too).

    Deliberately generous, and therefore usable in ONE direction only: a
    bigger set can suppress a false alarm about the delegate's own report,
    it can never ground a citation. ``echo "checked engine.py"`` names a
    file it did not open, so grounding reads ``_declared_paths`` instead.
    """
    out: list[str] = _declared_paths(args)
    text = raw if isinstance(raw, str) else ""
    if not text and args:
        text = " ".join(str(v) for v in args.values() if isinstance(v, str))
    for m in _PATH_TOKEN_RE.finditer(str(text)[:2000]):
        tok = m.group(0).strip(".,;:()[]\"'")
        if len(tok) > 3:
            out.append(tok)
    return out[:60]


def collect_report_evidence(payload, *, tool_calls=None) -> dict:
    """Reduce a sub-agent's recorded trace to the facts a report can be
    checked against.

    Evidence source order: an explicitly passed ``tool_calls`` list (the full
    in-process trace, including what each call returned), else the payload's
    own ``tool_calls`` entry. ``trace_available`` is False only when NO trace
    was recorded at all — missing bookkeeping is never treated as evidence of
    fabrication, so nothing is flagged in that case.

    Two path sets, because they are read in opposite directions:

    ``files``      — generous (command strings scraped for path-shaped
                     tokens, inputs of calls whose output was never
                     recorded). Only ever SUPPRESSES a flag about this
                     delegate's own report, so over-collecting is safe.
    ``files_read`` — a read or write tool, pointed at that path, that
                     returned something. This is the set published to the
                     parent as ``files_touched``, and the parent's honesty
                     guard is disarmed for every path in it, so a path a
                     delegate merely TYPED must never appear here.

    Returns ``{"trace_available", "calls", "tools", "files", "files_read",
    "commands", "test_runs", "writes", "write_paths", "test_output",
    "test_output_recorded"}``. Never raises.
    """
    ev: dict = {
        "trace_available": False,
        "calls": 0,
        "tools": [],
        "files": [],
        "files_read": [],
        "commands": [],
        "test_runs": [],
        "writes": [],
        "write_paths": [],
        "test_output": "",
        "test_output_recorded": False,
    }
    try:
        calls = tool_calls
        if calls is not None:
            ev["trace_available"] = True
        elif isinstance(payload, dict) and "tool_calls" in payload:
            ev["trace_available"] = True
            calls = payload.get("tool_calls")
        if not isinstance(calls, (list, tuple)):
            calls = []
        vg = _vg()
        test_cmd_re = getattr(vg, "_TEST_CMD_RE", None) if vg else None
        outs: list[str] = []
        for call in list(calls)[:_VERIFY_MAX_CALLS]:
            if not isinstance(call, dict):
                continue
            ev["calls"] += 1
            name = _norm_tool_name(call.get("name"))
            if name:
                ev["tools"].append(name)
            raw_in = call.get("input")
            args = _tool_args(raw_in)
            out_raw = call.get("output")
            if _output_saw_content(out_raw):
                for p in _paths_from_call(args, raw_in):
                    if p not in ev["files"]:
                        ev["files"].append(p)
                # Grounding needs the call to have RETURNED something: an
                # input without a recorded result is what a wall-clock cut
                # leaves behind, and counting it made a truncated run
                # ground more than a clean one, where a refused read is
                # dropped.
                if (out_raw is not None
                        and (name in _READ_TOOL_NAMES
                             or name in _WRITE_TOOL_NAMES)):
                    declared = _declared_paths(args)
                    for p in declared:
                        if p not in ev["files_read"]:
                            ev["files_read"].append(p)
                    # A repo-wide grep declares no path at all — its hits
                    # are in the OUTPUT, and the delegate read them there.
                    # The parent's own harvest has had this branch since
                    # grep was given the same treatment as the code-nav
                    # tools; the child's never got it, so a delegate sent
                    # out precisely to FIND something grounded the parent
                    # with nothing. Measured over 99 recorded delegate
                    # sessions: 107 grep_file calls, 22 of them with no
                    # path (21%), spread over 14 sessions. Every one of
                    # those made a relayed finding read as unsourced, so
                    # delegating cost a caveat that doing the work
                    # yourself did not.
                    #
                    # From the RESULT, never from the arguments: a path
                    # the delegate merely typed still grounds nothing,
                    # which is the rule #105 exists for.
                    if not declared and _norm_tool_name(name) == "grep_file":
                        for p in _paths_in_grep_hits(out_raw):
                            if p not in ev["files_read"]:
                                ev["files_read"].append(p)
            cmd = ""
            if vg is not None:
                try:
                    cmd = vg.extract_exec_command(call.get("name"), raw_in)
                except Exception:
                    cmd = ""
            if cmd:
                ev["commands"].append(cmd)
            ran_tests = name in _TEST_TOOL_NAMES or bool(
                cmd and test_cmd_re is not None and test_cmd_re.search(cmd))
            if ran_tests:
                ev["test_runs"].append(cmd or name)
                if out_raw is not None:
                    # Recorded, even when it is blank: an empty result is a
                    # fact about the run ("the output went to /dev/null"),
                    # unlike a trace that carries no results at all.
                    ev["test_output_recorded"] = True
                    outs.append(_readable_output(out_raw))
            if name in _WRITE_TOOL_NAMES or (cmd and _writes_anything(cmd)):
                ev["writes"].append(cmd or name)
                # Which files the write plausibly hit, so a mutation claim
                # can be paired with the path it names. Generous on purpose
                # (this direction only makes the pairing easier to
                # satisfy): a `sed -i` names its target in the command.
                for p in (_declared_paths(args) if name in _WRITE_TOOL_NAMES
                          else _paths_from_call(args, raw_in)):
                    if p not in ev["write_paths"]:
                        ev["write_paths"].append(p)
        # An isolated worktree reports WHAT it changed; those paths are write
        # evidence even though the edits happened inside the child.
        if isinstance(payload, dict):
            wt = payload.get("worktree")
            if isinstance(wt, dict):
                stat = str(wt.get("diff_summary") or "")[:4000]
                for m in _PATH_TOKEN_RE.finditer(stat):
                    tok = m.group(0).strip(".,;:()[]\"'|")
                    if len(tok) <= 3:
                        continue
                    # git --stat is the framework's own record of what the
                    # child changed, not the child's account of it.
                    for key in ("files", "files_read", "write_paths"):
                        if tok not in ev[key]:
                            ev[key].append(tok)
                    if tok not in ev["writes"]:
                        ev["writes"].append(tok)
        ev["test_output"] = "\n".join(outs)[:_VERIFY_MAX_OUTPUT_CHARS]
    except Exception:
        return ev
    return ev


def scan_report_test_claims(
    text: str, *, max_flags: int = _VERIFY_MAX_PER_KIND,
) -> list[ReportClaim]:
    """Find test-result claims ("all tests pass", "alle Tests grün", "60
    passed"). Hedged claims are skipped. Deterministic, de-duplicated,
    capped. Never raises."""
    claims: list[ReportClaim] = []
    try:
        if not text or not text.strip() or max_flags <= 0:
            return []
        scrubbed = _scrubbed(text[:_VERIFY_MAX_TEXT])
        found: list[tuple[int, str, str]] = []
        for shape, rx in _TEST_CLAIM_PATTERNS:
            for m in rx.finditer(scrubbed):
                if _hedged(scrubbed, m.start(), m.end()):
                    continue
                subject = ""
                if shape == "count":
                    try:
                        subject = str(m.group(1) or "")
                    except IndexError:
                        subject = ""
                found.append((m.start(), " ".join(m.group(0).split()), subject))
        found.sort(key=lambda t: t[0])
        seen: set[str] = set()
        for _pos, claim, subject in found:
            key = claim.lower()
            if key in seen:
                continue
            seen.add(key)
            claims.append(ReportClaim(kind="test_result",
                                      text=claim[:_VERIFY_CLAIM_CHARS],
                                      subject=subject))
            if len(claims) >= max_flags:
                break
    except Exception:
        return claims
    return claims


def scan_report_completion_claims(
    text: str, *, max_flags: int = _VERIFY_MAX_PER_KIND,
) -> list[ReportClaim]:
    """Find completion claims and tag their family.

    ``mutation``     — the delegate says it changed something ("implemented",
                       "fixed"); only counted in a first-person or bare-bullet
                       frame, so a finding about existing code is not mistaken
                       for a claim of authorship.
    ``verification`` — the delegate says it checked something or that
                       something is unchanged/identical.
    ``generic``      — "done", "erledigt", "finished".

    Hedged, negated, conditional and user-attributed sentences are skipped
    (reusing the shared claim filters). Never raises.
    """
    claims: list[ReportClaim] = []
    try:
        if not text or not text.strip() or max_flags <= 0:
            return []
        vg = _vg()
        scrubbed = _scrubbed(text[:_VERIFY_MAX_TEXT])
        found: list[tuple[int, str, str]] = []
        for family, rx in _COMPLETION_FAMILIES:
            for m in rx.finditer(scrubbed):
                if _hedged(scrubbed, m.start(), m.end()):
                    continue
                sentence = _sentence_of(scrubbed, m.start(), m.end())
                if not sentence.strip():
                    continue
                if vg is not None:
                    try:
                        if (vg._FUNC_DISCLOSURE_RE.search(sentence)
                                or vg._FUNC_NEGATION_RE.search(sentence)
                                or vg._FUNC_CONDITIONAL_RE.search(sentence)
                                or vg._FUNC_USER_SOURCE_RE.search(sentence)):
                            continue
                    except Exception:
                        pass
                if family == "mutation":
                    line, off = _line_of(scrubbed, m.start())
                    if not _is_agentive(line, off):
                        continue
                found.append((m.start(),
                              " ".join(sentence.split())[:_VERIFY_CLAIM_CHARS],
                              family))
        found.sort(key=lambda t: t[0])
        seen: set[str] = set()
        for _pos, claim, family in found:
            key = f"{family}:{claim.lower()}"
            if key in seen:
                continue
            seen.add(key)
            claims.append(ReportClaim(kind="completion", text=claim,
                                      family=family))
            if len(claims) >= max_flags:
                break
    except Exception:
        return claims
    return claims


def _report_text(payload) -> str:
    """The delegate's prose, wherever the payload carries it.

    ``result`` is the field the parent reads from a finished ``subagent``
    call; ``final_text`` is the equivalent in a collected background result.
    """
    if isinstance(payload, str):
        return payload
    if not isinstance(payload, dict):
        return ""
    for key in ("result", "final_text"):
        val = payload.get(key)
        if isinstance(val, str) and val.strip():
            return val
    return ""


def _unwritten_citation(claim: str, write_paths) -> str:
    """The file a mutation claim names when the run wrote none of them.

    ``ev["writes"]`` is a flat list of write-ish calls, so ANY write backed
    ANY mutation claim: ``mkdir -p /tmp/scratch`` corroborated "I
    implemented the fix in engine.py". Pairing the claim with the file it
    names closes that. A claim naming no file is unaffected, and one named
    file actually written is enough (a sentence may cite the test it fixed
    as well as the source it changed).

    Returns "" when there is nothing to flag. Never raises: a failure here
    means no extra flag, not a wrong one.
    """
    vg = _vg()
    if vg is None:
        return ""
    try:
        written = frozenset(str(p) for p in (write_paths or ()))
        cited: list[str] = []
        for path, _line in vg._iter_path_citations(claim):
            if vg._is_observed(path, written):
                return ""
            cited.append(path)
        return cited[0] if cited else ""
    except Exception:
        return ""


def _did_work(ev: dict) -> bool:
    """True when the trace contains a call that acted on the workspace.

    The pairing a generic "done" needs. ``ev["calls"] > 0`` was the whole
    test, so a run whose single tool call was ``task_update``,
    ``bash_status`` or a question to the user corroborated "erledigt" —
    any tool call at all counted as having finished the work. What can
    back a completion claim is a write, a command that ran, or a read
    that returned something; ``files_read`` rather than ``files`` on
    purpose, because the generous set includes paths merely typed into a
    command line.
    """
    return bool(ev.get("writes") or ev.get("commands")
                or ev.get("files_read") or ev.get("test_runs"))


def _inspected(ev: dict) -> bool:
    """True when the run looked at anything at all (read/grep/search/exec)."""
    if ev.get("files") or ev.get("commands"):
        return True
    vg = _vg()
    if vg is None:
        return bool(ev.get("tools"))
    try:
        return bool(vg._used_evidence_tool(frozenset(ev.get("tools") or ())))
    except Exception:
        return bool(ev.get("tools"))


def _verification_notice(verdict: dict) -> str:
    """One explicit line for the parent, naming what to re-check.

    Kept short and put at the FRONT of the tool result: the parent is a
    language model, so the verdict has to be in the text it reads, not only
    in a field it may skip.
    """
    status = verdict.get("status", "")
    prefix = ""
    if verdict.get("run_incomplete"):
        prefix = ("[subagent-verify] this run did not finish cleanly "
                  f"({verdict.get('run_incomplete')}) — the report describes "
                  "an incomplete run. ")
    if status == "no_report":
        return prefix.strip()
    if status == "no_trace":
        return (prefix + "[subagent-verify] no tool trace was recorded for "
                "this sub-agent, so none of its claims could be cross-checked "
                "— verify anything you rely on.")
    unsupported = verdict.get("unsupported") or []
    checked = int(verdict.get("claims_checked") or 0)
    if not unsupported:
        if not checked:
            # Silence read as approval. Saying nothing after a check that
            # ruled on nothing left the parent to assume the report had been
            # cleared, which is the opposite of what happened.
            return (prefix + "[subagent-verify] this report contains nothing "
                    "that can be cross-checked against the sub-agent's trace "
                    "— no claim was ruled on either way.")
        return (prefix + f"[subagent-verify] {checked} claim(s) cross-checked "
                f"against this sub-agent's own tool trace "
                f"({verdict.get('evidence', {}).get('tool_calls', 0)} tool "
                f"call(s)); all are backed by its own evidence. Limit: "
                f"{_VERIFY_LIMITS_NOTE}.")
    # Claims that fail for the SAME reason are listed together — the parent
    # needs the claims named and one instruction per reason, not the same
    # sentence repeated.
    groups: list[list] = []
    for u in unsupported:
        why = str(u.get("why", ""))
        recheck = str(u.get("recheck", ""))
        for g in groups:
            if g[0] == why and g[1] == recheck:
                g[2].append(str(u.get("claim", "")))
                break
        else:
            groups.append([why, recheck, [str(u.get("claim", ""))]])
    items = []
    shown = 0
    for i, (why, recheck, claims) in enumerate(
            groups[:_VERIFY_MAX_NOTICE_ITEMS], start=1):
        named = claims[:_VERIFY_MAX_NOTICE_ITEMS]
        shown += len(named)
        quoted = ", ".join(f'"{c}"' for c in named)
        if len(claims) > len(named):
            quoted += f" (+{len(claims) - len(named)} more)"
        items.append(f"({i}) {quoted} — {why}; re-check: {recheck}")
    more = ""
    if shown < len(unsupported):
        more = (f" (+{len(unsupported) - shown} further unsupported claim(s) "
                "in the 'verification' field)")
    return (prefix + f"[subagent-verify] {len(unsupported)} of {checked} "
            "claim(s) in this report are NOT backed by this sub-agent's own "
            "tool trace — do not pass them on unverified: "
            + " ".join(items) + more + ".")


def verify_subagent_report(payload, *, tool_calls=None, repo_root=None) -> dict:
    """Cross-check a delegate's report against the delegate's own evidence.

    ``payload`` is a sub-agent result payload (or the bare report text).
    ``tool_calls`` overrides the payload's trace with the full in-process one
    (tool outputs included). ``repo_root`` lets path citations be checked for
    existence as well as for having been opened.

    Claim classes and the rule each is judged by:

      test_result  — no test-running tool call in the trace -> unsupported;
                     recorded output reporting failures -> unsupported; a
                     test run whose recorded output is empty -> unsupported
                     (a discarded output is no evidence, not a pass); a
                     claimed count absent from the recorded test output ->
                     unsupported.
      file_reference / location / quantity / functional — delegated to the
                     shared scanners in ``verify_guard``, evaluated against
                     the files THIS sub-agent opened and the commands it ran.
      completion   — zero tool calls -> unsupported; a mutation claim with no
                     write in the trace, or naming only files this run never
                     wrote -> unsupported; a verification claim with nothing
                     read, grepped or executed -> unsupported; a generic
                     "done" with no write, no command and no read that
                     returned anything -> unsupported.

    Returns a compact, JSON-serializable verdict::

        {"status", "claims_checked", "supported", "unsupported",
         "evidence", "limits", "notice", "run_incomplete"?}

    Never raises: any internal failure degrades to
    ``{"status": "unavailable", ...}`` — the check must never cost the parent
    a delegate's work.
    """
    verdict: dict = {
        "status": "unavailable",
        "claims_checked": 0,
        "supported": [],
        "unsupported": [],
        "evidence": {},
        "limits": _VERIFY_LIMITS_NOTE,
        "notice": "",
    }
    try:
        # One cap for every scanner. The local ones already stopped at
        # _VERIFY_MAX_TEXT while the shared ones got the whole report, and
        # one of those is quadratic in the length of a single LINE: a report
        # arriving as one long unbroken line (a pasted blob, a minified
        # dump) took 7s at 20 kB and 47s at 50 kB, on the parent's turn.
        # Ordinary multi-line text of the same size costs 16ms.
        text = _report_text(payload)[:_VERIFY_MAX_TEXT]
        if isinstance(payload, dict):
            err = str(payload.get("error") or "").strip()
            if err:
                verdict["run_incomplete"] = err[:120]
            elif payload.get("truncated"):
                verdict["run_incomplete"] = "run was truncated (budget)"
        ev = collect_report_evidence(payload, tool_calls=tool_calls)
        verdict["evidence"] = {
            "tool_calls": ev["calls"],
            "tools": sorted(set(ev["tools"]))[:12],
            # The parent merges this into its own evidence ledger, so it
            # carries only what a read or write tool actually returned --
            # never a path scraped out of a command the delegate typed.
            "files_touched": ev["files_read"][:8],
            "test_runs": len(ev["test_runs"]),
            "writes": len(ev["writes"]),
        }
        if not text.strip():
            verdict["status"] = "no_report"
            verdict["notice"] = _verification_notice(verdict)
            return verdict
        if not ev["trace_available"]:
            verdict["status"] = "no_trace"
            verdict["notice"] = _verification_notice(verdict)
            return verdict

        verdict["status"] = "checked"
        supported: list[dict] = []
        unsupported: list[dict] = []

        def _add(ok: bool, kind: str, claim: str,
                 why: str = "", recheck: str = "") -> None:
            entry = {"kind": kind, "claim": str(claim)[:_VERIFY_CLAIM_CHARS]}
            if ok:
                if len(supported) < _VERIFY_MAX_SUPPORTED:
                    supported.append(entry)
                return
            entry["why"] = why
            entry["recheck"] = recheck
            if len(unsupported) < _VERIFY_MAX_UNSUPPORTED:
                unsupported.append(entry)

        # --- test-result claims -----------------------------------------
        for claim in scan_report_test_claims(text):
            if not ev["test_runs"]:
                _add(False, "test_result", claim.text,
                     "no test-running tool call in this sub-agent's trace",
                     "run the tests yourself before relying on this")
                continue
            if _output_shows_failures(ev["test_output"]):
                # The trace was consulted only for a claimed NUMBER, so a
                # delegate that ran the suite, saw failures and reported
                # success came back supported -- the check endorsing the
                # exact claim it exists to catch.
                _add(False, "test_result", claim.text,
                     "the test output recorded in this run reports failures",
                     "read the actual test output yourself")
                continue
            if ev["test_output_recorded"] and not ev["test_output"].strip():
                # A run that discarded its output is not a passing run; it
                # is no run anyone can check. Every rule below reads the
                # output, so an empty one used to fall through to
                # "supported" -- `pytest -q > /dev/null` certified as green.
                _add(False, "test_result", claim.text,
                     "the test run in this sub-agent's trace recorded no "
                     "output, so nothing in the run says how it ended",
                     "run the tests yourself and read the output")
                continue
            if (claim.subject and ev["test_output"]
                    and not re.search(rf"\b{re.escape(claim.subject)}\b",
                                      ev["test_output"])):
                _add(False, "test_result", claim.text,
                     "the test output recorded in this run does not contain "
                     f"the number {claim.subject}",
                     "read the actual test output yourself")
                continue
            _add(True, "test_result", claim.text)

        # --- completion claims ------------------------------------------
        for claim in scan_report_completion_claims(text):
            unwritten = (_unwritten_citation(claim.text, ev["write_paths"])
                         if claim.family == "mutation" else "")
            if ev["calls"] == 0:
                _add(False, "completion", claim.text,
                     "this sub-agent made no tool calls at all, so nothing in "
                     "the run backs this",
                     "treat the report as unverified and inspect the "
                     "workspace yourself")
            elif claim.family == "mutation" and not ev["writes"]:
                _add(False, "completion", claim.text,
                     "no file-writing tool call in this sub-agent's trace",
                     "diff the files before accepting the change")
            elif unwritten:
                _add(False, "completion", claim.text,
                     f"nothing in this sub-agent's trace wrote {unwritten}; a "
                     "write somewhere else does not make this change happen",
                     f"diff {unwritten} before accepting the change")
            elif claim.family == "verification" and not _inspected(ev):
                _add(False, "completion", claim.text,
                     "nothing was read, grepped or executed in this "
                     "sub-agent's trace",
                     "repeat the comparison yourself")
            elif claim.family == "generic" and not _did_work(ev):
                _add(False, "completion", claim.text,
                     "no tool call in this sub-agent's trace wrote, ran or "
                     "read anything, so nothing in the run shows work being "
                     "finished",
                     "check what was actually done before accepting this")
            else:
                _add(True, "completion", claim.text)

        # --- shared scanners (verify_guard) ------------------------------
        vg = _vg()
        observed = frozenset(str(p) for p in ev["files"])
        flagged_paths: set[str] = set()
        if vg is not None:
            try:
                for f in vg.scan_for_ungrounded_code_claims(
                        text, repo_root=repo_root, observed_files=observed,
                        max_flags=_VERIFY_MAX_PER_KIND):
                    ref = f"{f.path}:{f.line}" if f.line else f.path
                    flagged_paths.add(f.path.lower())
                    why = "this sub-agent never read or wrote that file"
                    if f.kind == "nonexistent":
                        why += " and it does not resolve as a file"
                    _add(False, "file_reference", ref, why,
                         f"open {f.path} yourself")
            except Exception:
                pass
            try:
                seen_ok: set[str] = set()
                for path, _line in vg._iter_path_citations(text):
                    key = path.lower()
                    if key in seen_ok or key in flagged_paths:
                        continue
                    seen_ok.add(key)
                    if vg._is_observed(path, observed):
                        _add(True, "file_reference", path)
            except Exception:
                pass
            try:
                for f in vg.scan_for_ungrounded_location_claims(
                        text, observed_files=observed, ledger_available=True,
                        max_flags=_VERIFY_MAX_PER_KIND):
                    if f.path and f.path.lower() in flagged_paths:
                        continue    # already reported as a file reference
                    _add(False, "location", f.claim,
                         "no read or grep of that location in this "
                         "sub-agent's trace",
                         "open the cited location yourself")
            except Exception:
                pass
            try:
                for f in vg.scan_for_unsourced_quantities(
                        text, observed_files=observed,
                        evidence_tools_used=frozenset(ev["tools"]),
                        max_flags=_VERIFY_MAX_PER_KIND):
                    _add(False, "quantity", f.quantity,
                         "no calculation output or lookup in this "
                         "sub-agent's trace",
                         "check where the number came from")
            except Exception:
                pass
            try:
                for f in vg.scan_for_unexercised_functional_claims(
                        text, exec_commands=list(ev["commands"]),
                        exec_ledger_available=True,
                        tools_used=frozenset(ev["tools"]),
                        max_flags=_VERIFY_MAX_PER_KIND):
                    why = ("this sub-agent never executed "
                           + (f"'{f.subject}'" if f.subject else "anything"))
                    if f.kind == "interactive":
                        why = ("interactive/browser behaviour was never "
                               "exercised in this sub-agent's run")
                    _add(False, "functional", f.claim, why,
                         "run it yourself before relying on this")
            except Exception:
                pass

        verdict["supported"] = supported
        verdict["unsupported"] = unsupported
        verdict["claims_checked"] = len(supported) + len(unsupported)
        if not verdict["claims_checked"]:
            # "checked" was assigned before any claim was evaluated. Leaving
            # it here would tell the parent the delegate's work was checked
            # and found sound, when in fact the report contained nothing this
            # check can rule on. That is the one thing a verifier must never
            # say: it would manufacture the confidence it exists to withhold.
            verdict["status"] = "no_claims"
        verdict["notice"] = _verification_notice(verdict)
        return verdict
    except Exception:
        # Degraded but well-formed: the parent still gets its delegate's work.
        #
        # The status has to be reset. By the time a scanner can raise, it is
        # already "checked", and returning that with zero claims and a blank
        # notice made a crash in the verifier indistinguishable from a clean
        # bill of health -- silently, because the notice is the only part the
        # parent reads and this path used to delete it.
        verdict["status"] = "unavailable"
        verdict["supported"] = []
        verdict["unsupported"] = []
        verdict["claims_checked"] = 0
        verdict["notice"] = (
            "[subagent-verify] the cross-check could not be completed, so "
            "none of this sub-agent's claims were checked — treat the report "
            "as unverified.")
        return verdict


def attach_verification(payload, *, tool_calls=None, repo_root=None):
    """Return the payload with the verdict attached — text never touched.

    The notice is placed FIRST so it heads the JSON the parent reads; the
    delegate's own fields (``result`` above all) are copied through
    byte-for-byte, and the machine-readable verdict lands under
    ``verification``. Idempotent: a payload that already carries a verdict is
    returned unchanged, so a second annotation pass cannot double-report.
    Non-dict input is returned as-is. Never raises.
    """
    if not isinstance(payload, dict):
        return payload
    if "verification" in payload:
        return payload
    try:
        verdict = verify_subagent_report(
            payload, tool_calls=tool_calls, repo_root=repo_root)
        notice = str(verdict.get("notice") or "")
        compact = {k: v for k, v in verdict.items() if k != "notice"}
        out: dict = {}
        if notice:
            out["verification_notice"] = notice
        out.update(payload)
        out["verification"] = compact
        return out
    except Exception:
        return payload


# Permission modes ordered from most to least restrictive. Used to clamp a
# child's mode so a subagent is never MORE permissive than its parent.
_MODE_RESTRICTIVENESS = {
    "plan": 0, "default": 1, "acceptEdits": 2, "bypassPermissions": 3,
}


def _clamp_child_mode(parent_mode: str, preset_mode: str) -> str:
    """Return the more restrictive of the parent's and the preset's mode.

    A subagent must never be able to do what its parent may not. Without this,
    a parent pinned to ``plan`` (the human-approval gate that ``exit_plan_mode``
    deliberately refuses to self-exit) could call ``subagent(...)`` whose preset
    carries ``mode="default"`` and obtain a child that writes/bashes on its
    behalf — an indirect plan-mode escape. Taking the lower (more restrictive)
    rank also preserves a preset that is intentionally stricter than the parent
    (e.g. read-only explore presets stay ``plan`` even from a ``default`` parent).
    """
    rank = _MODE_RESTRICTIVENESS
    if rank.get(parent_mode, 1) <= rank.get(preset_mode, 1):
        return parent_mode
    return preset_mode


def _derive_perms(parent_perms, mode: str, workspace=None):
    """Clone parent_perms with the sub-agent's mode + optional workspace.

    Both fields are optional — the caller passes ``workspace=None`` when
    no isolation is requested and the parent's workspace inherits. The child's
    mode is clamped so it is never more permissive than the parent's.
    """
    if parent_perms is None:
        return None
    mode = _clamp_child_mode(getattr(parent_perms, "mode", "") or "default", mode)
    # A locked session may not move its workspace. The child inherits the
    # lock through the role, but a workspace override (worktree isolation)
    # would relocate the boundary — the child would still be contained to
    # one folder, just not the folder the user confined the session to.
    if workspace is not None and getattr(parent_perms, "scope_locked", False):
        workspace = None
    # Bump nesting depth so the child's own subagent calls are refused once the
    # depth cap is hit (anti-recursion). Falls back gracefully if the perms
    # type predates the field.
    depth = int(getattr(parent_perms, "subagent_depth", 0)) + 1
    try:
        extra = {"subagent_depth": depth}
        # dataclasses.replace copies field VALUES, so a mutable one is
        # shared by reference. read_tracker is the stale-write guard's
        # memory of "this file was read at this mtime, so overwriting it is
        # safe" -- and that question is per agent. Shared, one child's write
        # bumped the entry and the next child's overwrite then saw an
        # unchanged mtime, passed, and clobbered it silently; the "no prior
        # read in this session" refusal was likewise satisfied by a
        # sibling's read. Give each child its own copy of what the parent
        # knew, so it inherits the history without writing into it.
        try:
            _tracker = getattr(parent_perms, "read_tracker", None)
            if isinstance(_tracker, dict):
                extra["read_tracker"] = dict(_tracker)
        except Exception:
            pass
        if workspace is not None:
            return dataclasses.replace(
                parent_perms, mode=mode, workspace=workspace, **extra)
        return dataclasses.replace(parent_perms, mode=mode, **extra)
    except Exception:
        try:
            if workspace is not None:
                return dataclasses.replace(
                    parent_perms, mode=mode, workspace=workspace)
            return dataclasses.replace(parent_perms, mode=mode)
        except Exception:
            return parent_perms


def _cheap_tier_enabled() -> bool:
    """Setting ``agent.subagents.cheap_tier`` — default ON.

    Gates the automatic cheap-tier routing of read-only subagents. Missing
    or unreadable settings mean the default (enabled) applies."""
    try:
        from delfin.user_settings import load_settings
        cfg = (((load_settings() or {}).get("agent") or {})
               .get("subagents") or {})
        return bool(cfg.get("cheap_tier", True))
    except Exception:
        return True


def _resolve_subagent_model(
    parent_client, subagent_type: str, model_override: str = "",
) -> tuple[str, str]:
    """Pick the model a subagent should run on: ``(model, tier)``.

    ``tier`` is ``"cheap"`` when the run routes to the parent provider's
    cheap tier, else ``"parent"``. Read-only presets (explore / plan /
    code-reviewer …) do their work fine on a cheap-tier model — running
    them on the parent's expensive model just burns budget. Rules:

      - writer presets ALWAYS keep the parent model (their edits are the
        part that is expensive to redo);
      - ``model_override="parent"`` pins the parent model for a hard
        read-only task; ``"cheap"`` requests cheap routing even when the
        ``agent.subagents.cheap_tier`` setting is off;
      - otherwise the setting (default ON) gates routing;
      - the cheap candidate must resolve for the parent's provider, differ
        from the parent model, not be marked broken, and support tool
        calling (subagents drive tools — a no-tools model returns prose
        instead of doing the work).

    Never raises: ANY failure keeps the parent model, so routing can
    never break a subagent run.
    """
    parent_model = str(getattr(parent_client, "model", "") or "")
    try:
        override = (model_override or "").strip().lower()
        if override == "parent":
            return parent_model, "parent"
        if is_writer_preset(subagent_type):
            return parent_model, "parent"
        if override != "cheap" and not _cheap_tier_enabled():
            return parent_model, "parent"
        provider = str(
            getattr(parent_client, "_provider", "")
            or getattr(parent_client, "provider", "") or ""
        ).strip().lower()
        if not provider:
            return parent_model, "parent"
        from .model_routing import is_known_broken, tier_model
        candidate = (tier_model(provider, "cheap", None) or "").strip()
        if not candidate or candidate == parent_model:
            return parent_model, "parent"
        if is_known_broken(candidate):
            return parent_model, "parent"
        from .model_capabilities import resolve as _resolve_caps
        base_url = str(getattr(parent_client, "_base_url", "") or "")
        caps = _resolve_caps(provider, candidate, base_url)
        if not getattr(caps, "supports_tools", False):
            return parent_model, "parent"
        return candidate, "cheap"
    except Exception:
        return parent_model, "parent"


def run_subagent(
    *,
    subagent_type: str,
    description: str,
    prompt: str,
    parent_client: "OpenAIClient",
    parent_perms: Optional["KitToolPermissions"],
    max_tool_calls: int | None = None,
    max_wall_s: float | None = None,
    max_output_tokens: int | None = None,
    isolation: str = "",
    resume_from: str = "",
    sa_id: str = "",
    model: str = "",
    inherit_context: bool = True,
    output_schema: Optional[dict] = None,
) -> SubagentResult:
    """Run a sub-agent loop and return its final assistant message.

    The sub-agent shares the parent's underlying OpenAI client (so the
    same KIT-Toolbox endpoint and API key) but with an isolated
    message list and (usually) tighter permissions.

    ``isolation="worktree"`` creates a fresh git worktree under
    ``$TMPDIR/delfin-wt-<hex>`` and points the sub-agent's workspace
    there, so any edits stay off the user's working tree until
    explicitly reviewed/merged. Falls back gracefully if the parent
    workspace isn't a git repo (just runs without isolation +
    reports a warning on the result).

    ``resume_from`` continues a FINISHED subagent: the stored
    conversation (user/assistant turns of all prior rounds) is replayed
    in front of the new prompt, and the stored ``subagent_type`` wins so
    permissions match the original preset. The session file accumulates
    across resumes under the same id.

    ``model`` is a per-call override for the cheap-tier routing of
    read-only presets: ``"parent"`` pins the parent's model (for a hard
    read-only task), ``"cheap"`` requests the cheap tier even when the
    ``agent.subagents.cheap_tier`` setting is off. Empty (the default)
    keeps the setting-driven behaviour. See ``_resolve_subagent_model``.

    ``inherit_context=True`` (default) appends the parent workspace's
    project rules (DELFIN.MD/AGENTS.md) and BM25-selected feedback/user
    typed memories to the child's system prompt — see
    ``_inherited_context``. Opt out with ``inherit_context=False``.
    Best-effort: memory failures never break the spawn.

    ``output_schema`` (a JSON-Schema dict, subset documented at
    ``validate_json_schema``) enforces a structured return: the child is
    instructed to end with exactly one matching JSON object; the final
    text is parsed and validated, with ONE correction round on mismatch.
    Success → ``SubagentResult.structured_output``; persistent mismatch →
    ``structured_output=None`` + ``schema_error`` (text stays available).

    Returns a SubagentResult; never raises.
    """
    resume_from = (resume_from or "").strip()
    prior: dict = {}
    if resume_from:
        prior = load_subagent_session(resume_from) or {}
        if not prior:
            return SubagentResult(
                subagent_type=subagent_type,
                description=description,
                final_text="",
                error=(
                    f"unknown resume id {resume_from!r}: no stored "
                    "subagent session found. Use a sa_id returned by a "
                    "previous subagent call."
                ),
            )
        subagent_type = prior.get("subagent_type") or subagent_type
        description = description or prior.get("description", "")

    preset = SUBAGENT_PRESETS.get(subagent_type)
    if preset is None:
        return SubagentResult(
            subagent_type=subagent_type,
            description=description,
            final_text="",
            error=f"unknown subagent_type: {subagent_type!r}. Pick one of {list(SUBAGENT_PRESETS)}.",
        )

    # Resolve limits (settings-configurable) unless the caller pinned them.
    _lim = _subagent_limits()
    if max_tool_calls is None:
        max_tool_calls = _lim["max_tool_calls"]
    if max_wall_s is None:
        max_wall_s = _lim["max_wall_s"]
    if max_output_tokens is None:
        max_output_tokens = _lim["max_output_tokens"]

    # Optionally spin up an isolated git worktree before deriving the
    # sub-agent's permissions. Defer the tear-down to the finally block.
    isolation = (isolation or "").strip().lower()
    worktree_info = None
    isolation_warning = ""
    if (isolation == "worktree" and parent_perms is not None
            and _worktree_would_be_voided(parent_perms)):
        # Do not create what the permission layer is about to discard, and
        # do not let the payload claim an isolation that did not happen.
        isolation = ""
        isolation_warning = _ISOLATION_VOIDED_BY_LOCK
    if isolation == "worktree" and parent_perms is not None:
        try:
            from .worktree import enter_worktree as _enter_wt
            parent_ws = getattr(parent_perms, "workspace", None)
            if parent_ws is not None:
                worktree_info = _enter_wt(parent_ws)
        except Exception as exc:
            isolation_warning = (
                f"worktree isolation requested but failed ({exc}); "
                "subagent ran in the parent workspace."
            )
            worktree_info = None

    sub_workspace = worktree_info.path if worktree_info else None
    sub_perms = _derive_perms(parent_perms, preset.mode, workspace=sub_workspace)

    # Run against an isolated SHALLOW COPY of the parent client: it shares the
    # same underlying OpenAI client (endpoint/key/model — thread-safe to share)
    # but carries its OWN _permissions. This makes concurrent fan-out subagents
    # safe: previously each run mutated the SHARED parent's _permissions (swap
    # + restore), a race that could run a subagent under a sibling's sandbox
    # (e.g. a read-only explorer transiently gaining a builder's write perms).
    # The parent client is never touched on this path.
    import copy as _copy
    saved_perms = None
    restore_parent = False
    try:
        sub_client = _copy.copy(parent_client)
        # A shallow copy shares the parent's MUTABLE evidence ledgers: the
        # child would union its own reads into the parent's session set in
        # place, so a file only the delegate opened grounded the parent's
        # citations -- even when the parent discarded the report (the
        # fan-out timeout path abandons the thread, it does not stop it).
        # The verified merge in _merge_delegate_evidence is the one route
        # a delegate's reads may take to the parent; give the child its own
        # ledger so that stays the only one.
        for _ledger in ("_observed_files", "_observed_files_session"):
            if hasattr(sub_client, _ledger):
                try:
                    setattr(sub_client, _ledger, set())
                except Exception:
                    pass
        sub_client.set_permissions(sub_perms)
    except Exception:
        # Fallback to the legacy swap-on-parent if copying fails.
        sub_client = parent_client
        saved_perms = getattr(parent_client, "_permissions", None)
        restore_parent = True
        if hasattr(parent_client, "set_permissions"):
            try:
                parent_client.set_permissions(sub_perms)
            except Exception:
                pass

    # Cheap-tier routing: read-only presets run fine on the provider's cheap
    # model. Switch the ISOLATED copy only — the legacy fallback above shares
    # the parent client, so it must never have its model touched. Any failure
    # keeps the parent model (routing must never break a subagent run).
    sub_model = str(getattr(parent_client, "model", "") or "")
    model_tier = "parent"
    if not restore_parent:
        try:
            _routed, _routed_tier = _resolve_subagent_model(
                parent_client, subagent_type, model)
            if _routed_tier == "cheap" and _routed:
                try:
                    if hasattr(sub_client, "switch_model"):
                        sub_client.switch_model(_routed)
                    else:
                        sub_client.model = _routed
                except Exception:
                    pass
                # Record the routed model only if the client actually took it
                # (the OpenAI path reads self.model per request).
                if str(getattr(sub_client, "model", "") or "") == _routed:
                    sub_model = _routed
                    model_tier = "cheap"
        except Exception:
            pass

    system_prompt = (
        preset.system_prompt
        + f"\n\nWorkspace: {sub_perms.workspace if sub_perms else '(none)'}"
        + f"\nTask label: {description}"
    )
    # Context inheritance: project rules + standing memories, best-effort
    # (a failure inside _inherited_context degrades to "" — spawning must
    # never depend on the memory subsystem being healthy).
    if inherit_context:
        try:
            system_prompt += _inherited_context(
                getattr(sub_perms, "workspace", None), prompt)
        except Exception:
            pass
    if output_schema and isinstance(output_schema, dict):
        system_prompt += _schema_instruction(output_schema)
    # On resume, inject a faithful recap (earlier conversation + the tool
    # outputs the subagent actually saw) as one context block, then the
    # new request. The system prompt is rebuilt fresh so workspace paths
    # stay current.
    messages = [{"role": "system", "content": system_prompt}]
    if prior:
        recap = _render_resume_recap(prior)
        if recap.strip():
            messages.append({"role": "user", "content": recap})
            messages.append({"role": "assistant", "content": (
                "Understood — I have my earlier findings and the tool "
                "outputs they produced. Continuing from there.")})
    messages.append({"role": "user", "content": prompt})

    final_text_parts: list[str] = []
    tool_calls_seen: list[dict] = []
    in_tokens = out_tokens = 0
    t0 = time.monotonic()
    truncated = False
    error = ""
    # Live-panel registry entry (removed in the finally below). Resumes
    # keep their original id so the session file accumulates.
    import uuid as _uuid
    # Honour a caller-supplied id (background runs reserve it up-front so the
    # parent can poll/retrieve the result); else resume keeps its id; else new.
    _sa_id = resume_from if prior else ((sa_id or "").strip() or _uuid.uuid4().hex[:8])
    _sa_started = time.time()
    _sa_actions: list[str] = []
    # Rich live transcript (text the subagent writes + tool calls with brief
    # in/out) so the dashboard can let you "go INTO" a running subagent and
    # watch its activity in the chat window.
    _sa_transcript: list = []
    _sa_text_buf: list[str] = []

    def _sa_push() -> None:
        _running_update(_sa_id, {
            "type": subagent_type,
            "description": (description or "")[:120],
            "started_at": _sa_started,
            "actions": _sa_actions[-12:],
            "last_action": _sa_actions[-1] if _sa_actions else "",
            "transcript": _sa_transcript[-40:],
        })

    def _sa_flush_text() -> None:
        if _sa_text_buf:
            _t = "".join(_sa_text_buf).strip()
            _sa_text_buf.clear()
            if _t:
                _sa_transcript.append({"k": "text", "v": _t[:1500]})

    _running_update(_sa_id, {
        "type": subagent_type,
        "description": (description or "")[:120],
        "started_at": _sa_started,
        "actions": [],
        "last_action": "",
        "transcript": [],
    })

    try:
        for event in sub_client.stream_message(
            messages=messages,
            system=system_prompt,
            max_tokens=max_output_tokens,
        ):
            if time.monotonic() - t0 > max_wall_s:
                truncated = True
                error = f"wall-clock budget exhausted ({max_wall_s:.0f}s)"
                break
            if len(tool_calls_seen) > max_tool_calls:
                truncated = True
                error = f"tool-call budget exhausted ({max_tool_calls})"
                break
            if event.type == "text_delta" and event.text:
                final_text_parts.append(event.text)
                _sa_text_buf.append(event.text)
            elif event.type == "tool_use":
                tool_calls_seen.append({
                    "name": event.tool_name,
                    "input": event.tool_input,
                })
                # Live drill-down: record this step so the dashboard can show
                # what THIS subagent is doing right now (read/write/bash …).
                _sa_flush_text()
                _sa_actions.append(
                    _format_action(event.tool_name, event.tool_input))
                _sa_transcript.append({
                    "k": "tool", "name": event.tool_name or "?",
                    "in": _format_action(event.tool_name, event.tool_input)})
                _sa_push()
            elif event.type == "tool_result":
                # Attach the output to the most recent tool call still
                # missing one (preserves what the subagent actually saw).
                _attached = False
                for tc in reversed(tool_calls_seen):
                    if "output" not in tc:
                        tc["output"] = event.tool_output
                        _attached = True
                        break
                if not _attached:
                    tool_calls_seen.append({
                        "name": event.tool_name or "?",
                        "input": "",
                        "output": event.tool_output,
                    })
                for _e in reversed(_sa_transcript):
                    if _e.get("k") == "tool" and "out" not in _e:
                        _e["out"] = str(event.tool_output or "")[:400]
                        break
                _sa_push()
            elif event.type == "message_delta":
                in_tokens = max(in_tokens, event.input_tokens)
                out_tokens = max(out_tokens, event.output_tokens)
    except Exception as exc:
        error = f"sub-agent stream raised: {exc}"
    finally:
        _running_update(_sa_id, None)
        # Restore parent permissions ONLY if we fell back to running on the
        # shared parent client (the normal isolated-copy path never touched it).
        if restore_parent and hasattr(parent_client, "set_permissions"):
            try:
                parent_client.set_permissions(saved_perms)
            except Exception:
                pass

    # If an isolated worktree was created, tear it down — keep the dir
    # when there are local changes so the user can review/merge.
    worktree_summary: dict = {}
    _held_jobs = (_jobs_holding_worktree(worktree_info.path)
                  if worktree_info is not None else [])
    if worktree_info is not None and _held_jobs:
        # The tree is still in use. Removing it would take the live
        # processes' working directory and their job registry with it —
        # see _kept_for_running_jobs.
        worktree_summary = _kept_for_running_jobs(worktree_info, _held_jobs)
    elif worktree_info is not None:
        try:
            from .worktree import exit_worktree as _exit_wt
            info = _exit_wt(worktree_info, keep_if_changed=True)
            worktree_summary = {
                "branch": info.branch,
                "had_changes": bool(info.had_changes),
                "final_path": str(info.final_path) if info.final_path else "",
                "cleaned_up": bool(info.cleaned_up),
            }
            # Surface WHAT the isolated subagent changed so the parent can
            # review parallel-writer work without hunting through the worktree.
            if info.had_changes:
                try:
                    from .worktree import diff_summary as _diff_sum
                    _ds = _diff_sum(info)
                    if _ds:
                        worktree_summary["diff_summary"] = _ds
                except Exception:
                    pass
        except Exception as exc:
            worktree_summary = {"error": f"worktree teardown failed: {exc}"}

    if isolation_warning and not error:
        # Surface as a soft warning when nothing else went wrong.
        worktree_summary["warning"] = isolation_warning

    final_text = "".join(final_text_parts).strip()
    if not final_text and not error:
        error = "sub-agent returned no text"

    # Schema-validated return: parse + validate the final message, with
    # exactly ONE correction round (a second text-only turn carrying the
    # validation errors) before giving up. A failed/truncated run skips
    # validation — there is nothing trustworthy to validate.
    structured_output: Optional[dict] = None
    schema_error = ""
    if (output_schema and isinstance(output_schema, dict)
            and final_text and not error and not truncated):
        obj = extract_json_object(final_text)
        errs = (validate_json_schema(obj, output_schema)
                if obj is not None
                else ["no JSON object found in the final message"])
        if not errs:
            structured_output = obj
        else:
            correction = (
                "Your final message did not satisfy the required output "
                "schema: " + "; ".join(errs)[:600]
                + "\nReply again with ONLY one JSON object matching the "
                "schema from your instructions — no prose, no explanation."
            )
            messages.append({"role": "assistant", "content": final_text})
            messages.append({"role": "user", "content": correction})
            fixed, c_in, c_out, c_err = _stream_text_only(
                sub_client, messages, system_prompt,
                max_output_tokens, deadline=t0 + max_wall_s)
            in_tokens += c_in
            out_tokens += c_out
            if fixed:
                final_text = fixed  # the corrected message IS the report
                obj2 = extract_json_object(fixed)
                errs2 = (validate_json_schema(obj2, output_schema)
                         if obj2 is not None
                         else ["no JSON object found in the corrected "
                               "message"])
                if not errs2:
                    structured_output = obj2
                else:
                    schema_error = "; ".join(errs2)
            else:
                schema_error = "; ".join(errs)
                if c_err:
                    schema_error += f" (correction round failed: {c_err})"

    elapsed_s = time.monotonic() - t0
    # Persist the conversation so the parent can resume this subagent
    # later via ``resume_id`` (resume-by-id). Store the
    # LOGICAL conversation (clean user/assistant turns) + the tool
    # interactions WITH outputs, accumulating across resumes — decoupled
    # from the recap-laden ``messages`` actually sent to the model.
    prior_messages = [
        {"role": m.get("role", ""), "content": m.get("content", "")}
        for m in (prior.get("messages") or [])
        if m.get("role") in ("user", "assistant")
    ]
    this_round = [{"role": "user", "content": prompt}]
    if final_text:
        this_round.append({"role": "assistant", "content": final_text})
    all_interactions = (prior.get("interactions") or []) + tool_calls_seen
    _save_subagent_session(
        _sa_id,
        subagent_type=subagent_type,
        description=description,
        messages=prior_messages + this_round,
        interactions=all_interactions,
        error=error,
        worktree=worktree_summary,
    )
    # Persist a telemetry record so the dashboard /agents stats command
    # and the subagent-pane can show real costs/durations across runs.
    _write_telemetry({
        "ts": time.time(),
        "subagent_type": subagent_type,
        "description": description,
        "isolation": isolation,
        "model": sub_model,
        "model_tier": model_tier,
        "elapsed_s": round(elapsed_s, 3),
        "input_tokens": in_tokens,
        "output_tokens": out_tokens,
        "tool_calls_count": len(tool_calls_seen),
        "truncated": truncated,
        "error": error,
    })
    return SubagentResult(
        subagent_type=subagent_type,
        description=description,
        final_text=final_text,
        tool_calls=tool_calls_seen,
        elapsed_s=elapsed_s,
        input_tokens=in_tokens,
        output_tokens=out_tokens,
        truncated=truncated,
        error=error,
        worktree=worktree_summary,
        sa_id=_sa_id,
        model=sub_model,
        model_tier=model_tier,
        workspace=str(getattr(sub_perms, "workspace", "") or ""),
        structured_output=structured_output,
        schema_error=schema_error,
    )


def _jobs_holding_worktree(path) -> list[dict]:
    """Background jobs still running inside a sub-agent's isolated tree.

    One implementation, in :mod:`delfin.agent.worktree`, so the sub-agent
    teardown, the ``exit_worktree`` tool and the post-merge cleanup cannot
    drift apart on which of them checks. Best-effort: a registry that
    cannot be read reports nothing, which keeps the teardown behaving
    exactly as it did before this check."""
    try:
        from .worktree import jobs_holding_worktree
        return jobs_holding_worktree(path) or []
    except Exception:
        return []


def _kept_for_running_jobs(info, jobs: list[dict]) -> dict:
    """Refuse the teardown and say why, with the job ids to act on.

    ``bash_background`` registers a job under the workspace it was started
    in, and for an isolated sub-agent that workspace IS the worktree. The
    tree was then removed with ``--force`` while the child process — in
    its own session, so it survives — kept running to its 24-hour cap:
    its working directory deleted underneath it, its completion event
    destroyed with the registry file, and ``bash_status`` answering
    ``unknown job_id``, which reads as a typo rather than as a live
    process nobody is counting.

    Refusing to remove a tree that is still in use is the honest half of
    that trade. Re-registering the record against the parent workspace
    would rescue the bookkeeping and still leave a running process whose
    cwd had been deleted — a record about a job whose tree no longer
    exists, filed under a workspace that never ran it. Keeping the tree
    costs a directory that the report names; the parent can end the jobs
    with ``bash_kill`` or let them finish.
    """
    ids = [str(j.get("job_id") or "?") for j in jobs][:8]
    return {
        "branch": getattr(info, "branch", ""),
        "had_changes": True,
        "final_path": str(getattr(info, "path", "")),
        "cleaned_up": False,
        "running_jobs": ids,
        "warning": (
            f"the isolated worktree was NOT removed: {len(jobs)} background "
            f"job(s) started by this sub-agent are still running in it "
            f"({', '.join(ids)}). Removing it would delete their working "
            f"directory and their completion records while the processes "
            f"keep running. End them with bash_kill(job_id) or wait for "
            f"them; the tree stays at {getattr(info, 'path', '')}."
        ),
    }


def auto_isolation_for(
    subagent_type: str, requested: str, *, background: bool,
) -> str:
    """The isolation a subagent should get when the call did not pin one.

    Auto-isolation for writers lived only in the parallel fan-out, which
    triggers at two or more subagent calls in one turn. A background writer
    is a single call and got none -- and it is the case that needs it most:
    fan-out siblings at least finish before the turn ends, while a
    background writer edits the tree while the parent is editing it, with
    no lock and no shared stale-write baseline between them.

    A single FOREGROUND writer is left alone on purpose: it has the tree to
    itself, so a worktree would only add a merge step nobody asked for.
    """
    requested = (requested or "").strip()
    if requested:
        return requested
    if background and is_writer_preset(subagent_type):
        return "worktree"
    return ""


# Stated, not implied: the caller asked for isolation and is not getting it.
_ISOLATION_VOIDED_BY_LOCK = (
    "worktree isolation was not applied: this session is confined to one "
    "folder, and moving the sub-agent into a worktree would move that "
    "boundary. The sub-agent ran in the parent workspace and is NOT "
    "isolated from it."
)


def _worktree_would_be_voided(parent_perms) -> bool:
    """Whether ``_derive_perms`` will drop a workspace override anyway.

    Checked BEFORE the worktree is created. It used to be created, left
    unused because the lock dropped the override, and torn down again --
    while the result still reported ``isolation: "worktree"``."""
    return bool(getattr(parent_perms, "scope_locked", False))


def is_writer_preset(subagent_type: str) -> bool:
    """True if a subagent type can WRITE (mutate the workspace) — i.e. its
    preset is NOT the read-only 'plan' permission profile. Used to auto-isolate
    parallel writers into separate worktrees so they can't clobber each other."""
    p = SUBAGENT_PRESETS.get((subagent_type or "").strip())
    return bool(p) and (p.mode or "").strip().lower() != "plan"


# ---------------------------------------------------------------------------
# Declarative orchestration — a deterministic multi-stage executor.
#
# The DRIVER is pure library code: the model only fills in the spec (via the
# 'orchestrate' tool) and does the work inside the spawned children. Stage
# order, barriers, template substitution, verification votes and majority
# arithmetic are all deterministic Python — no model in the driver seat.
# ---------------------------------------------------------------------------

_ORCH_MAX_STAGES = 3
_ORCH_MAX_CALLS_PER_STAGE = 6
_ORCH_MAX_VOTES = 3
# Mirrors the api_client fan-out ThreadPool cap (4 workers).
_ORCH_MAX_WORKERS = 4

# The skeptic ballot every verify vote must return (schema mechanics from
# run_subagent's output_schema path).
_ORCH_VERIFY_SCHEMA: dict = {
    "type": "object",
    "properties": {
        "refuted": {"type": "boolean"},
        "note": {"type": "string"},
    },
    "required": ["refuted"],
}


def _substitute_stage_refs(prompt: str, stage_results: dict) -> str:
    """Replace ``{{stage:NAME}}`` placeholders with that stage's results.

    Substitution value is the JSON-compact list of the named stage's
    payloads. Placeholders naming an unknown (or not-yet-run) stage are
    left untouched — stages only see PRIOR stages (barrier semantics).
    """
    out = str(prompt or "")
    if "{{stage:" not in out:
        return out
    for name, payloads in stage_results.items():
        token = "{{stage:" + name + "}}"
        if token in out:
            try:
                rendered = json.dumps(payloads, separators=(",", ":"),
                                      ensure_ascii=False)
            except (TypeError, ValueError):
                rendered = "[]"
            out = out.replace(token, rendered)
    return out


def _validate_orchestration_spec(spec, parent_perms) -> tuple[list, dict, str]:
    """Check spec shape + hard limits. Returns ``(stages, verify, error)``."""
    if not isinstance(spec, dict):
        return [], {}, "spec must be a JSON object with a 'stages' list"
    # No nested orchestration: children of an orchestration run at subagent
    # depth 1 (same cap as any spawned subagent) and an orchestration driven
    # FROM inside a subagent would multiply fan-out beyond the depth cap.
    if int(getattr(parent_perms, "subagent_depth", 0) or 0) >= 1:
        return [], {}, (
            "orchestration cannot be nested: it may only be driven by the "
            "top-level agent (orchestration children already run at "
            "subagent depth 1)")
    stages = spec.get("stages")
    if not isinstance(stages, list) or not stages:
        return [], {}, "spec.stages must be a non-empty list"
    if len(stages) > _ORCH_MAX_STAGES:
        return [], {}, (
            f"too many stages: {len(stages)} > max {_ORCH_MAX_STAGES}")
    seen_names: set[str] = set()
    for idx, stage in enumerate(stages):
        if not isinstance(stage, dict):
            return [], {}, f"stage #{idx + 1} must be an object"
        name = str(stage.get("name") or "").strip()
        if not name:
            return [], {}, f"stage #{idx + 1} needs a non-empty 'name'"
        if name in seen_names:
            return [], {}, (
                f"duplicate stage name {name!r} — template references "
                "would be ambiguous")
        seen_names.add(name)
        calls = stage.get("parallel")
        if not isinstance(calls, list) or not calls:
            return [], {}, (
                f"stage {name!r} needs a non-empty 'parallel' list")
        if len(calls) > _ORCH_MAX_CALLS_PER_STAGE:
            return [], {}, (
                f"stage {name!r} has {len(calls)} calls > max "
                f"{_ORCH_MAX_CALLS_PER_STAGE} per stage")
        for cidx, call in enumerate(calls):
            if not isinstance(call, dict):
                return [], {}, (
                    f"stage {name!r} call #{cidx + 1} must be an object")
            for req_key in ("subagent_type", "description", "prompt"):
                if not str(call.get(req_key) or "").strip():
                    return [], {}, (
                        f"stage {name!r} call #{cidx + 1} is missing "
                        f"{req_key!r}")
    verify = spec.get("verify") or {}
    if verify:
        if not isinstance(verify, dict):
            return [], {}, "spec.verify must be an object"
        if not str(verify.get("prompt_template") or "").strip():
            return [], {}, "spec.verify needs a 'prompt_template'"
    return stages, verify, ""


def _run_stage_calls(jobs: list[dict], parent_client, parent_perms) -> list[dict]:
    """Run one batch of subagent-call dicts, fanned out on a thread pool.

    ``jobs`` carry FINAL kwargs (prompts already substituted). Results come
    back in submission order — deterministic regardless of completion
    order. A single job runs inline (no pool). run_subagent never raises,
    but a raise anywhere still degrades to an error payload.
    """
    def _one(job: dict) -> dict:
        try:
            return run_subagent(
                subagent_type=str(job.get("subagent_type") or ""),
                description=str(job.get("description") or ""),
                prompt=str(job.get("prompt") or ""),
                parent_client=parent_client,
                parent_perms=parent_perms,
                isolation=str(job.get("isolation") or ""),
                output_schema=(job.get("output_schema")
                               if isinstance(job.get("output_schema"), dict)
                               else None),
            ).to_payload()
        except Exception as exc:
            return {
                "subagent_type": str(job.get("subagent_type") or ""),
                "description": str(job.get("description") or ""),
                "result": "",
                "error": f"orchestration call raised: {exc}",
            }

    return _collect_stage(jobs, _one)


def _collect_stage(jobs: list[dict], run_one) -> list[dict]:
    """Run one stage's delegates concurrently under ONE deadline.

    Named rather than inlined in the closure above so the policy has one
    home and a test can reach it: what matters here is a property of the
    stage, not of whichever caller happens to build the jobs.
    """
    if len(jobs) == 1:
        return [run_one(jobs[0])]
    import concurrent.futures as _cf
    # ONE deadline for the whole stage, in completion order.
    #
    # The tool-call fan-out learned this already: waiting on the futures in
    # submission order bought nothing (the pool runs them concurrently) and
    # a stage barrier with no deadline at all is worse -- one stalled
    # delegate holds every sibling's finished work, and the stage, and the
    # turn, for as long as it stalls. The child's own wall-clock guard
    # fires per streamed event, so a fully SILENT stream never trips it.
    #
    # A straggler is reported as one rather than dropped: its slot in the
    # returned list keeps the job it belonged to, so a caller can still
    # tell which of six calls did not come back.
    deadline = _orchestration_stage_timeout()
    results: list[dict] = [None] * len(jobs)          # type: ignore[list-item]
    # NOT a `with` block. Leaving one calls shutdown(wait=True), which
    # waits for every worker -- so the deadline would collect the results
    # early and the stage would sit there anyway until the straggler
    # returned. Measured at 30s against a 0.6s deadline before this was
    # written the other way; the timing assertion in the test is the only
    # thing that catches it, because every value in the result is already
    # correct by then.
    pool = _cf.ThreadPoolExecutor(
        max_workers=min(len(jobs), _ORCH_MAX_WORKERS),
        thread_name_prefix="orchestration",
    )
    try:
        futures = {pool.submit(run_one, job): index
                   for index, job in enumerate(jobs)}
        try:
            for fut in _cf.as_completed(list(futures), timeout=deadline):
                index = futures[fut]
                try:
                    results[index] = fut.result()
                except Exception as exc:
                    results[index] = {
                        "subagent_type": str(
                            (jobs[index] or {}).get("subagent_type") or ""),
                        "description": str(
                            (jobs[index] or {}).get("description") or ""),
                        "result": "",
                        "error": f"subagent failed: "
                                 f"{str(exc) or type(exc).__name__}",
                    }
        except Exception:
            pass          # the deadline; whatever landed is kept
    finally:
        # Queued work that never started is dropped; a delegate already
        # running keeps running and is reported below as unfinished --
        # the same contract the tool-call fan-out states for its own
        # stragglers.
        try:
            pool.shutdown(wait=False, cancel_futures=True)
        except TypeError:                      # pre-3.9 signature
            pool.shutdown(wait=False)
    for index, entry in enumerate(results):
        if entry is not None:
            continue
        job = jobs[index] or {}
        results[index] = {
            "subagent_type": str(job.get("subagent_type") or ""),
            "description": str(job.get("description") or ""),
            "result": "",
            "error": (f"did not finish within the stage deadline of "
                      f"{deadline:.0f}s and was left running; nothing it "
                      f"produced is in this result."),
        }
    return results



# How long a whole orchestration stage may wait for its delegates. Same
# shape as the tool-call fan-out's bound: a bit past the child's own
# wall-clock budget, because that guard fires per streamed event and a
# fully silent stream never trips it.
def _orchestration_stage_timeout() -> float:
    try:
        wall = float(_subagent_limits().get("max_wall_s", 300) or 300)
    except Exception:
        wall = 300.0
    return wall + 120.0


def run_orchestration(spec: dict, parent_client, parent_perms) -> dict:
    """Execute a declarative multi-stage subagent plan. Never raises.

    ``spec`` shape::

        {"stages": [{"name": str,
                     "parallel": [{"subagent_type": str, "prompt": str,
                                   "description": str, "output_schema"?: {},
                                   "isolation"?: "worktree"}, ...]},
                    ...],
         "verify": {"prompt_template": str,   # optional
                    "votes": int,             # clamped to 1..3
                    "subagent_type"?: str}}   # default "explore"

    Semantics:
      - stages run strictly in order with a barrier between them;
      - within a stage, calls fan out on a thread pool (cap 4 workers,
        matching the api_client fan-out); writers in a multi-call stage
        auto-isolate into worktrees unless the call pins isolation;
      - ``{{stage:NAME}}`` in a later stage's prompt is replaced with the
        JSON-compact results list of that finished stage;
      - the optional verify step runs ``votes`` parallel skeptic subagents
        per FINAL-stage result (prompt_template with ``{{result}}``
        substituted), each returning ``{refuted, note}`` via the
        output-schema mechanism; a strict majority of the valid ballots
        marking ``refuted=true`` moves the result to ``rejected``.

    Hard limits: max 3 stages, max 6 calls per stage, max 3 votes;
    children run at subagent depth 1 (so they cannot fan out further) and
    orchestration itself refuses to run from inside a subagent.

    Returns ``{"stages": {name: [payloads]}, "verified": [...],
    "rejected": [...]}`` — or ``{"error": ...}`` on a spec violation.
    Writes ONE summary telemetry record for the whole run.
    """
    t0 = time.monotonic()
    stages, verify, spec_err = _validate_orchestration_spec(spec, parent_perms)
    if spec_err:
        return {"error": spec_err}

    stage_results: dict[str, list[dict]] = {}
    stage_names: list[str] = []
    for stage in stages:
        name = str(stage.get("name")).strip()
        calls = stage.get("parallel")
        jobs: list[dict] = []
        for call in calls:
            job = {
                "subagent_type": str(call.get("subagent_type") or "").strip(),
                "description": str(call.get("description") or "").strip(),
                # Substitute against COMPLETED stages only (barrier: this
                # stage's own name is not in stage_results yet).
                "prompt": _substitute_stage_refs(
                    str(call.get("prompt") or ""), stage_results),
                "isolation": str(call.get("isolation") or "").strip(),
                "output_schema": call.get("output_schema"),
            }
            # Same auto-isolation rule as the api_client fan-out: parallel
            # WRITERS each get their own worktree so concurrent edits can't
            # clobber one another; read-only presets need none.
            try:
                if (len(calls) >= 2 and not job["isolation"]
                        and is_writer_preset(job["subagent_type"])):
                    job["isolation"] = "worktree"
            except Exception:
                pass
            jobs.append(job)
        stage_results[name] = _run_stage_calls(jobs, parent_client,
                                               parent_perms)
        stage_names.append(name)

    # ---- optional verification votes over the FINAL stage ---------------
    final_payloads = stage_results[stage_names[-1]]
    verified: list[dict] = []
    rejected: list[dict] = []
    vote_payloads_all: list[dict] = []
    if verify:
        try:
            votes = int(verify.get("votes") or 1)
        except (TypeError, ValueError):
            votes = 1
        votes = max(1, min(_ORCH_MAX_VOTES, votes))
        template = str(verify.get("prompt_template") or "")
        skeptic_type = (str(verify.get("subagent_type") or "").strip()
                        or "explore")
        vote_jobs: list[dict] = []
        for ridx, payload in enumerate(final_payloads):
            result_repr = json.dumps(
                {"description": payload.get("description", ""),
                 "result": payload.get("result", ""),
                 "structured_output": payload.get("structured_output")},
                separators=(",", ":"), ensure_ascii=False)
            for v in range(votes):
                vote_jobs.append({
                    "subagent_type": skeptic_type,
                    "description": f"verify result #{ridx + 1} "
                                   f"(vote {v + 1}/{votes})",
                    "prompt": template.replace("{{result}}", result_repr),
                    "isolation": "",
                    "output_schema": _ORCH_VERIFY_SCHEMA,
                    "_ridx": ridx,
                })
        ballots = _run_stage_calls(vote_jobs, parent_client, parent_perms)
        vote_payloads_all = ballots
        for ridx, payload in enumerate(final_payloads):
            mine = [b for j, b in zip(vote_jobs, ballots)
                    if j.get("_ridx") == ridx]
            valid = [b["structured_output"] for b in mine
                     if isinstance(b.get("structured_output"), dict)]
            refuted_n = sum(1 for s in valid if s.get("refuted") is True)
            notes = [str(s.get("note") or "") for s in valid
                     if s.get("note")]
            payload["verify"] = {
                "votes": votes,
                "valid_votes": len(valid),
                "refuted": refuted_n,
                "notes": notes[:_ORCH_MAX_VOTES],
            }
            # Strict majority of VALID ballots. No valid ballots → the
            # result stands (a broken skeptic must not veto real work),
            # but valid_votes=0 in the payload flags it for the caller.
            if valid and refuted_n * 2 > len(valid):
                rejected.append(payload)
            else:
                verified.append(payload)
    else:
        verified = list(final_payloads)

    # ---- one summary telemetry record for the whole orchestration --------
    all_payloads = [p for plist in stage_results.values() for p in plist]
    all_payloads += vote_payloads_all
    first_error = next((p.get("error") for p in all_payloads
                        if p.get("error")), "")
    _write_telemetry({
        "ts": time.time(),
        "subagent_type": "orchestration",
        "description": (f"{len(stage_names)} stage(s), "
                        f"{sum(len(stage_results[n]) for n in stage_names)} "
                        f"call(s), {len(vote_payloads_all)} vote(s)"),
        "isolation": "",
        "model": str(getattr(parent_client, "model", "") or ""),
        "model_tier": "parent",
        "elapsed_s": round(time.monotonic() - t0, 3),
        "input_tokens": sum(int(p.get("input_tokens") or 0)
                            for p in all_payloads),
        "output_tokens": sum(int(p.get("output_tokens") or 0)
                             for p in all_payloads),
        "tool_calls_count": sum(len(p.get("tool_calls") or [])
                                for p in all_payloads),
        "truncated": any(p.get("truncated") for p in all_payloads),
        "error": str(first_error or ""),
        "stages": stage_names,
        "verified": len(verified),
        "rejected": len(rejected),
    })
    return {
        "stages": stage_results,
        "verified": verified,
        "rejected": rejected,
    }


def get_subagent_result(sa_id: str) -> dict:
    """Fetch a background subagent's status/result by id.

    Returns ``{status: "running"|"finished"|"unknown", ...}``. For a finished
    run, ``final_text`` is the subagent's last report. Lets the parent agent
    (or the dashboard) collect a backgrounded subagent's output without
    resuming it.
    """
    sa_id = (sa_id or "").strip()
    if not sa_id:
        return {"error": "sa_id is required"}
    running = read_running()
    if sa_id in running:
        ent = running[sa_id] or {}
        return {"sa_id": sa_id, "status": "running",
                "subagent_type": ent.get("type", ""),
                "description": ent.get("description", ""),
                "started_at": ent.get("started_at", 0)}
    sess = load_subagent_session(sa_id)
    if not sess:
        # Distinguish "never existed" from "was running when its process
        # died". The caller is usually a parent waiting for a report: told
        # "running" it waits forever, told "unknown" it assumes it made the
        # id up. Only "died" lets it start the work again.
        dead = read_running(include_dead=True).get(sa_id)
        if dead:
            why = str(dead.get("error") or "").strip()
            # Collecting IS the delivery: whoever asked has now been told,
            # so the push channel must not announce it a second time.
            _claim_pending_report(sa_id)
            return {"sa_id": sa_id, "status": "died",
                    "subagent_type": dead.get("type", ""),
                    "description": dead.get("description", ""),
                    "started_at": dead.get("started_at", 0),
                    "error": ((why + " — the work has to be started again")
                              if why else
                              ("the process running this subagent is no "
                               "longer running, and it wrote no report — "
                               "the work has to be started again"))}
        return {"sa_id": sa_id, "status": "unknown",
                "error": "no running or finished subagent with this id"}
    # The untrimmed report when it was stored; the conversation is the
    # fallback for sessions written before the field existed.
    final = str(sess.get("final_report") or "")
    if not final:
        for m in reversed(sess.get("messages") or []):
            if m.get("role") == "assistant" and m.get("content"):
                final = m.get("content", "")
                break
    out = {"sa_id": sa_id, "status": "finished",
           "subagent_type": sess.get("subagent_type", ""),
           "description": sess.get("description", ""),
           "final_text": final,
           "error": sess.get("error", "")}
    # What the framework recorded about the run's isolated tree — kept
    # trees, and the reason one was kept, above all.
    wt = sess.get("worktree")
    if isinstance(wt, dict) and wt:
        out["worktree"] = wt
    _claim_pending_report(sa_id)
    # A collected background report gets the same cross-check as a foreground
    # one — the stored session keeps the tool calls WITH their outputs.
    return attach_verification(
        out, tool_calls=(sess.get("interactions") or []))


__all__ = [
    "SubagentPreset",
    "SUBAGENT_PRESETS",
    "list_subagents",
    "subagent_type_names",
    "reload_subagent_presets",
    "SubagentResult",
    "ReportClaim",
    "collect_report_evidence",
    "scan_report_test_claims",
    "scan_report_completion_claims",
    "verify_subagent_report",
    "attach_verification",
    "run_subagent",
    "run_orchestration",
    "validate_json_schema",
    "extract_json_object",
    "load_subagent_session",
    "list_finished",
    "get_subagent_result",
    "reserve_running",
    "mark_running_died",
    "drain_finished_subagents",
    "reap_pending_reports",
]
