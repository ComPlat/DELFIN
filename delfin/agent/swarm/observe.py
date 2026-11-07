"""Observe a swarm and one member session (package SWV) — read-only.

The ONE data API the dashboard (tab_swarm.py), the TUI (tui/src/swarm) and
an SW3 lead call to see two things the person asked for (06.11.): what a
swarm and its members are doing, and to step into one member and see what
it does. This module is pure: JSON-able dicts in, JSON-able dicts out; it
never starts, pauses, answers, prompts or writes anything.

Two functions:

* :func:`swarm_overview` — one row per member of a swarm (or of an
  operator wave given as a plain roster): package, role, model, state,
  commits on its branch beyond the base, last tool call, last commit
  subject, age of the last activity, pending approvals count.
* :func:`member_feed` — a member's recent steps (tool calls with result
  status, refusals, commits), newest last.

Data writers are EXISTING stores, never a second trace:

* ``audit_log.steps_since`` — the per-session step stream ``cmd_watch``
  already reads (cli.py:4690). ``member_feed`` reuses it; there is no new
  trace of its own.
* the ``session_health.collect()`` record shape (open_approvals, paused,
  idle_s, no_token_s, last_turn.stopped) — the same record the health
  classifier reads.

Both writers are injectable (``records=``, ``steps_fn=``, ``commits_fn=``)
so tests feed real-shaped fixtures; the production wrappers call
``session_health.collect()`` and ``audit_log.steps_since``.

Liveness is NEVER a process grep (agreed practice 6): it comes from state
files and the state record, exactly as session_health derives it. Nothing
secret is ever read beyond what the audit log already shows; the redaction
layer does the masking, this module does not lift it.

Ownership (a member's session is shown only to the user who owns it, the
same check the service uses) is enforced in the service layer; this module
renders only the data it is given. The QS ruling (relayed via the operator):
swarm.observe.* are pure reads only, ownership enforced server-side.
"""

from __future__ import annotations

import os
from typing import Any, Callable, Optional

#: The user-facing states of a member row (task SWV, target 5). They are a
#: faithful projection of the session record, in a fixed priority so the
#: loudest condition wins: an open approval beats everything, then a pause,
#: then an ended turn, then "waiting for model".
STATE_NEEDS_APPROVAL = "NEEDS APPROVAL"
STATE_PAUSED = "paused"
STATE_ENDED = "ended"
STATE_WAITING_MODEL = "waiting for model"
STATE_WORKING = "working"

#: How long (s) without a model token counts as "waiting for model". We do
#: not read settings here (pure module): this is a stable, documented floor.
WAIT_ON_MODEL_S = 10 * 60


def _state_of(record: Optional[dict]) -> str:
    """The SWV state of a member from its session-health-shaped record.

    ``None`` (no record at all) is treated as ``working``: we cannot claim
    a pause, an approval, a wait or an ended turn from absence, and a
    member that has simply not been seen is not invented into a liveness
    that the stores do not prove.
    """
    if not record:
        return STATE_WORKING
    if record.get("open_approvals"):
        return STATE_NEEDS_APPROVAL
    if record.get("paused"):
        return STATE_PAUSED
    turn = record.get("last_turn") or {}
    if turn.get("stopped") is True and not record.get("resume_planned"):
        return STATE_ENDED
    if (record.get("no_token_s") or 0) >= WAIT_ON_MODEL_S:
        return STATE_WAITING_MODEL
    return STATE_WORKING


def _members_of(source) -> list[dict]:
    """The member list a ``source`` describes: package + role per member.

    ``source`` is either an SW1 spec dict (``{"packages": [{"name", ...
    "roles"}]}``; the same shape ``status._status`` reads) or a plain
    roster file (path str/Path) mapping member name -> {package, role}
    (so an operator wave with no spec shows too). Unknown shapes raise a
    ValueError naming the problem; a swarm with no packages is empty.
    """
    if isinstance(source, dict):
        packages = source.get("packages")
        if isinstance(packages, list):
            out = []
            for _pkg in packages:
                if not isinstance(_pkg, dict) or not _pkg.get("name"):
                    continue
                roles = _pkg.get("roles") or {"builder": True}
                if not isinstance(roles, dict):
                    roles = {"builder": True}
                for role in roles:
                    out.append({
                        "name": f"{_pkg['name']}-{role}",
                        "package": _pkg["name"], "role": role,
                    })
            return out
        # A bare roster dict is also accepted: name -> {package, role}.
        out = []
        for name, meta in (source.items() or {}).items():
            if not isinstance(meta, dict):
                continue
            out.append({
                "name": str(name),
                "package": str(meta.get("package") or ""),
                "role": str(meta.get("role") or ""),
            })
        return out
    # A path (str or Path) to a roster file (JSON or YAML: name -> {..}).
    path = os.fspath(source) if not isinstance(source, str) else source
    if not path:
        raise ValueError("source must be a spec dict, a roster dict or a "
                         "roster file path")
    if not os.path.exists(path):
        raise ValueError(f"roster file does not exist: {path}")
    data = _read_roster(path)
    out = []
    for name, meta in (data or {}).items():
        if not isinstance(meta, dict):
            continue
        out.append({
            "name": str(name),
            "package": str(meta.get("package") or ""),
            "role": str(meta.get("role") or ""),
        })
    return out


def _read_roster(path: str) -> dict:
    """Load a roster file's member map: JSON, or YAML when pyyaml is here."""
    text = open(path, "r", encoding="utf-8").read()
    try:
        import json
        return json.loads(text) if text.strip().startswith("{") else {}
    except Exception:
        pass
    try:
        import yaml  # type: ignore
        data = yaml.safe_load(text)
        return data if isinstance(data, dict) else {}
    except Exception:
        raise ValueError(f"roster file is neither JSON nor YAML: {path}") \
            from None


# ---------------------------------------------------------------------------
# production writers (injectable for tests)
# ---------------------------------------------------------------------------

def _default_records():
    """The real health records for the open sessions (never raises)."""
    from .. import session_health
    try:
        return session_health.collect()
    except Exception:
        return []


def _default_steps(offset=0, *, session="", denied_only=False, log_path=None,
                   **kw):
    """The real step stream for one member, reusing audit_log (cmd_watch's
    source). Returns (rows, next_offset)."""
    from .. import audit_log as _al
    return _al.steps_since(offset, session=session, denied_only=denied_only,
                           log_path=log_path)


# ---------------------------------------------------------------------------
# swarm_overview
# ---------------------------------------------------------------------------

def swarm_overview(source,
                   *,
                   records: Optional[list] = None,
                   steps_fn: Optional[Callable] = None,
                   commits_fn: Optional[Callable] = None,
                   now: Optional[float] = None) -> list[dict]:
    """One row per member of a swarm (or a roster), read-only.

    ``records``: the session-health-shaped records for the member sessions
    (default: ``session_health.collect()``). ``steps_fn``: the step stream
    reader (default: ``audit_log.steps_since``), called per member with
    ``session=<member name>``. ``commits_fn(source)`` returns
    ``{member_name: [commit}]}`` (default: reads the member's branch
    behind its base via safe_git; never raises).

    Each row: ``{package, role, name, model, state, commits_on_branch,
    last_tool_call, last_commit_subject, last_activity_age_s,
    pending_approvals_count}`` — all JSON-able.
    """
    import time as _time
    from .. import audit_log as _al

    steps_fn = steps_fn or _default_steps
    commits_fn = commits_fn or _default_commits
    records = list(records or [])
    now = now if now is not None else _time.time()
    rec_by_name = {str(r.get("name") or r.get("session_id") or ""): r
                   for r in records if isinstance(r, dict)}

    rows = []
    for member in _members_of(source):
        name = member["name"]
        rec = rec_by_name.get(name)
        # The member's latest steps (reused from the audit log) drive the
        # last-tool-call field; the step rows are already newest-last when
        # read with a fresh cursor.
        steps, _offset = steps_fn(0, session=name)
        last = steps[-1] if steps else None

        commits = {}
        try:
            commits = commits_fn(name)
        except Exception:
            commits = {}
        member_commits = commits.get(name) or []

        rows.append({
            "package": member["package"],
            "role": member["role"],
            "name": name,
            "model": str(rec.get("title") or "") if rec else "",
            "state": _state_of(rec),
            "commits_on_branch": len(member_commits),
            "last_tool_call": {
                "tool": str(last.get("tool") or "") if last else "",
                "arg": " ".join(str(last.get("command") or "").split())[:96]
                       if last else "",
            } if last else None,
            "last_commit_subject": (str(member_commits[-1]) if member_commits
                                    else ""),
            "last_activity_age_s": float(rec.get("idle_s") or 0.0) if rec
                                   else 0.0,
            "pending_approvals_count": len(rec.get("open_approvals") or [])
                                       if rec else 0,
        })
    return rows


def _default_commits(member_name: str) -> dict:
    """Commit subjects of one member's branch behind its base, from safe_git.

    A swarm member's branch is ``swarm/<pkg>-<role>`` (launch.py:60) and its
    name is ``<pkg>-<role>``, so the branch is ``swarm/<name>``. The base is
    the branch's merge-base with ``main``; commits beyond it name what the
    member has added. Never raises, never writes: a missing branch, a closed
    repo or an error yields an empty list, not an exception.
    """
    from .. import safe_git
    if not member_name:
        return {member_name: []}
    branch = f"swarm/{member_name}"
    try:
        # The repo root: safe_git resolves the working tree internally; we
        # hand it the current directory.
        cwd = os.getcwd()
        if not safe_git.is_repository(cwd):
            return {member_name: []}
        base = safe_git.out(cwd, "merge-base", "main", branch).strip()
        if not base:
            return {member_name: []}
        log = safe_git.out(cwd, "log", "--format=%s", f"{base}..{branch}")
        subjects = [ln.strip() for ln in log.splitlines() if ln.strip()]
        return {member_name: subjects}
    except Exception:
        return {member_name: []}


# ---------------------------------------------------------------------------
# member_feed
# ---------------------------------------------------------------------------

def member_feed(session: str,
                since: int = 0,
                *,
                records: Optional[list] = None,
                steps_fn: Optional[Callable] = None,
                commits_fn: Optional[Callable] = None,
                log_path=None) -> list[dict]:
    """A member's recent steps, newest last.

    Steams tool calls (with result status) and REFUSALS from the audit-log
    step stream (reused from ``cmd_watch``) and adds the member's commits.
    Each entry is ``{at, kind, text}`` where ``kind`` is ``tool``,
    ``refusal`` or ``commit`` — JSON-able, never raising. ``records`` and
    ``steps_fn``/``commits_fn``/``log_path`` are the injectable writers
    (defaults: the real stores).
    """
    from .. import audit_log as _al

    steps_fn = steps_fn or _default_steps
    entries: list[dict] = []

    try:
        rows, _offset = steps_fn(int(since or 0), session=session,
                                 log_path=log_path)
    except Exception:
        rows = []
    for row in rows or []:
        if not isinstance(row, dict):
            continue
        denied = str(row.get("decision") or "") == "denied"
        kind = "refusal" if denied else "tool"
        text = str(row.get("command") or row.get("path")
                   or row.get("tool") or "")
        if denied:
            text = f"{text}\n        └─ {row.get('reason') or ''}"
        entries.append({
            "at": row.get("ts") or "",
            "kind": kind,
            "text": " ".join(text.split()),
        })

    try:
        commits = (commits_fn or _default_commits)(session)
    except Exception:
        commits = {}
    # Commits from git log are newest-first; the feed contract is newest
    # last, so invert ONLY the commits. The steps above were left in
    # audit_log.steps_since order (oldest-first, i.e. newest-last) and must
    # NOT be blanket-reversed -- a blanket entries.reverse() over both runs
    # the steps back to newest-first and jumbles the feed (F1).
    entries.extend(
        {"at": "", "kind": "commit", "text": str(subject)}
        for subject in (commits.get(session) or [])[::-1]
    )
    return entries
