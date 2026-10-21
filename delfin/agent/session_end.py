"""Shared session-end for the DELFIN agent.

One call site for everything that happens when a chat session ends —
today the dashboard's tab_agent and the CLI's cmd_chat/cmd_run each
grew their own. This module is the single place:

- skill learning (package 2): a finished session that qualified and
  carries evidence yields at most ONE skill proposal. The raw
  conversation alone is not enough: the CLI's saved ``chat_messages``
  keep no tool results, so the session's tool trace
  (~/.delfin/tool_traces/<session>.jsonl, written by engine for every
  tool call) is merged in as synthetic tool messages before
  qualification and evidence extraction.
- more end-of-session stages (session indexing, memory upkeep) attach
  here later; each runs in its own try/except so one stage's failure
  never breaks the others or the session save itself.

All text here is English; a learning failure never raises.
"""

from __future__ import annotations


def skill_learning_settings(settings: dict | None) -> dict:
    """The ``agent.skill_learning`` block, default enabled.

    Enabling only means proposals are written for a human to review —
    no code path ever activates a proposed skill automatically.
    """
    try:
        if settings is None:
            from delfin.user_settings import load_settings
            settings = load_settings()
        cfg = ((settings or {}).get("agent") or {}).get(
            "skill_learning") or {}
    except Exception:
        cfg = {}
    return {
        "enabled": bool(cfg.get("enabled", True)),
    }


def _synthetic_tool_messages(session_id: str,
                             max_entries: int = 400) -> list[dict]:
    """The session's tool trace as chat-message-shaped tool entries.

    The CLI's saved sessions keep user/assistant text only; the trace
    holds the tool names and outputs. Mapping trace entries to the
    same HTML-chip form the dashboard writes keeps ONE parser in
    skill_learning. Best-effort: any failure yields an empty list.
    """
    try:
        from delfin.agent import tool_trace as tt
        entries = tt.read(session_id or "", last_n=max_entries)
    except Exception:
        return []
    out: list[dict] = []
    for e in entries:
        tool = str(e.get("tool") or "")
        if not tool:
            continue
        # The trace stores the tool's OUTPUT as text already; embedding
        # it raw (not json.dumps'd) keeps summary JSON like
        # {"passed": 4} parseable by skill_learning's evidence
        # patterns — double-encoded quotes broke every match.
        output_text = str(e.get("output") or e.get("error") or "")
        out.append({
            "role": "tool",
            "content": (
                f'<span class="tool-name">{tool}</span>  '
                f'<span class="tool-param">{e.get("input") or ""}</span>'
                f'<details><summary> &rarr; {output_text}</details>'),
        })
    return out


def _session_test_runs(session_id: str) -> list[dict]:
    """The saved session's test ledger (``evidence["tests"]``), or []."""
    try:
        from delfin.agent import session_store
        data = session_store.load_session(session_id) if session_id else None
        tests = ((data or {}).get("evidence") or {}).get("tests") or []
        return [dict(t) for t in tests if isinstance(t, dict)]
    except Exception:
        return []


def _run_verdict(run) -> str:
    """"green", "red" or "" for one recorded test run of the session."""
    if not isinstance(run, dict):
        return ""
    def _n(key):
        try:
            return int(run.get(key) or 0)
        except (TypeError, ValueError):
            return 0
    status = str(run.get("status") or run.get("state")
                 or run.get("outcome") or "").strip().lower()
    exit_code = run.get("exit_code")
    if (_n("failed") or _n("errors") or _n("error")
            or status in ("failed", "fail", "error", "red")
            or (isinstance(exit_code, int) and exit_code != 0)):
        return "red"
    if _n("passed") or status in ("passed", "ok", "success", "green"):
        return "green"
    return ""


def resolve_pending_session_uses(session_id: str) -> int:
    """Give the skills THIS session invoked their outcome. Never raises.

    A skill is recorded as "started" when its text is handed to the agent
    (api_client._execute_skill); what came of it is known only now.
    success: at least one green test run and no red one. failure: a red
    run. Neither: "unverified" -- a session that ran no test proves
    nothing about the skill, and calling that a failure would be as much
    an invention as calling it a success. Only rows tagged with this
    session's id are touched, each "started" row gets one resolving row,
    and a second end of the same session resolves only what was started
    after the first. Returns how many uses were resolved.

    Limit: the outcome is the whole session's test record. A skill used
    for one thing in a session whose unrelated test went red is resolved
    "failure" all the same; nothing here links a test run to a skill.
    Cost: one session load plus one bounded read per ledger file.
    """
    try:
        sid = str(session_id or "").strip()
        if not sid:
            return 0
        from delfin.agent import skill_usage as _su
        names = _su.pending_for_session(sid)
        if not names:
            return 0
        verdicts = {_run_verdict(r) for r in _session_test_runs(sid)}
        outcome = ("failure" if "red" in verdicts
                   else "success" if "green" in verdicts else "unverified")
        return sum(_su.resolve_session(name, sid, outcome)
                   for name in names)
    except Exception:
        return 0


def learn_at_session_end(messages, *, session_id: str = "",
                         settings: dict | None = None, llm=None,
                         _propose=None):
    """Run skill learning for an ending session. Never raises.

    Merges the session's tool trace into the chat messages (the CLI
    path has no tool results in them), then delegates to
    skill_learning.learn_from_session. Returns the proposal or None.
    """
    # The usage record is not learning: it is kept whether or not proposals
    # are switched on, and here because this is the one hook the CLI run,
    # the CLI chat and the dashboard all call when a session ends.
    resolve_pending_session_uses(session_id)
    try:
        cfg = skill_learning_settings(settings)
        if not cfg["enabled"]:
            return None
        from delfin.agent.skill_learning import learn_from_session
        msgs = [m for m in (messages or []) if isinstance(m, dict)]
        msgs = msgs + _synthetic_tool_messages(session_id)
        return learn_from_session(msgs, settings=settings, llm=llm,
                                  _propose=_propose,
                                  runs=_session_test_runs(session_id))
    except Exception:
        return None


def index_at_session_end(session_id: str) -> bool:
    """Session indexing (package 5): the ended session becomes searchable.

    The same hook's second stage, in its own try/except like the skill
    learning above: an indexing failure never breaks the session save
    (or the other stages). Deliberately NOT gated on
    agent.skill_learning.enabled — that setting governs proposals, not
    the read-only search index.
    """
    try:
        from delfin.agent.session_end_index import reindex_finished_session
        return reindex_finished_session(session_id)
    except Exception:
        return False
