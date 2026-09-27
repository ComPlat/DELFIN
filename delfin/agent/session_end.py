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


def learn_at_session_end(messages, *, session_id: str = "",
                         settings: dict | None = None, llm=None,
                         _propose=None):
    """Run skill learning for an ending session. Never raises.

    Merges the session's tool trace into the chat messages (the CLI
    path has no tool results in them), then delegates to
    skill_learning.learn_from_session. Returns the proposal or None.
    """
    try:
        cfg = skill_learning_settings(settings)
        if not cfg["enabled"]:
            return None
        from delfin.agent.skill_learning import learn_from_session
        msgs = [m for m in (messages or []) if isinstance(m, dict)]
        msgs = msgs + _synthetic_tool_messages(session_id)
        return learn_from_session(msgs, settings=settings, llm=llm,
                                  _propose=_propose)
    except Exception:
        return None
