"""What a line starting with a punctuation mark means.

Four prefixes and one table. The table is data, not an if-chain, because
the help text is GENERATED from it — a hand-written help page and a
command list are one edit away from disagreeing, and `help_gen` exists
because that already happened once.

Nothing here talks to a terminal. A command handler takes a context and a
string and returns a CommandResult; the loop decides what to do with it.
That is what makes every routing rule testable without a PTY, and it is
the line that keeps this file from becoming the 4 500-line dispatcher it
replaces.
"""

from __future__ import annotations

import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

__all__ = [
    "CommandResult", "ReplCommand", "BUILTINS", "dispatch",
    "looks_like_command", "palette_rows", "expand_at_references",
    "complete_path", "SHELL_PREFIX", "MEMORY_PREFIX",
]

SHELL_PREFIX = "!"
MEMORY_PREFIX = "#"


@dataclass
class CommandResult:
    """What the loop should do next."""

    output: str = ""          # printed as chrome
    prompt: str = ""          # send THIS to the model instead of the raw line
    quit: bool = False
    handled: bool = True
    clear: bool = False       # start a fresh conversation


@dataclass(frozen=True)
class ReplCommand:
    name: str
    category: str
    summary: str
    handler: Callable[[Any, str], CommandResult]
    takes_args: bool = False


def looks_like_command(text: str) -> bool:
    """A leading slash is a command only when it cannot be a path.

    The dashboard's predicate, kept identical on purpose: `/home/u/x.py`
    is a path someone pasted, `/` alone is nothing, and `/help` is a
    command. Diverging here would mean the same line means two different
    things on the two surfaces.
    """
    text = (text or "").strip()
    if not text.startswith("/") or len(text) < 2:
        return False
    first = text.split()[0]
    if first.count("/") > 1:            # /home/user/... is a path
        return False
    if first.startswith("/."):
        return False
    return True


# ---------------------------------------------------------------------------
# @path and the file completer
# ---------------------------------------------------------------------------

def complete_path(fragment: str, workspace: Path) -> list[str]:
    """Paths under *workspace* matching *fragment*.

    Only under the workspace. Completing outside it would advertise files
    the agent may not read, which is a worse experience than no
    completion: it invites a request that is then refused.
    """
    fragment = (fragment or "").lstrip("@")
    try:
        base = (workspace / fragment).parent if fragment else workspace
        stem = Path(fragment).name if fragment else ""
        base = base.resolve()
        root = workspace.resolve()
        if base != root and root not in base.parents:
            return []
        out = []
        for entry in sorted(base.iterdir()):
            if entry.name.startswith(".") and not stem.startswith("."):
                continue
            if stem and not entry.name.startswith(stem):
                continue
            rel = entry.relative_to(root)
            out.append(f"{rel}/" if entry.is_dir() else str(rel))
        return out[:50]
    except Exception:
        return []


def expand_at_references(text: str, workspace: Path) -> str:
    """Turn `@path` into something the model can act on.

    The path is only ANNOTATED, never read here: reading it would put file
    contents into the conversation without going through the gate that
    decides whether they may be read at all.
    """
    out: list[str] = []
    found: list[str] = []
    for token in text.split():
        if token.startswith("@") and len(token) > 1:
            rel = token[1:]
            target = workspace / rel
            found.append(rel + ("" if target.exists() else "  (not found)"))
            out.append(rel)
        else:
            out.append(token)
    body = " ".join(out)
    if found:
        body += "\n\nFiles referenced: " + ", ".join(found)
    return body


# ---------------------------------------------------------------------------
# The builtin handlers
# ---------------------------------------------------------------------------

def _status(ctx, _args: str) -> CommandResult:
    engine = ctx.engine
    try:
        st = engine.get_status() or {}
    except Exception:
        st = {}
    lines = [
        f"model      {getattr(getattr(engine, 'client', None), 'model', '?')}",
        f"provider   {st.get('provider', '?')}",
        f"mode       {st.get('mode', '?')}   role {st.get('role', '?')}",
        f"approval   {_mode(engine) or '(no gate on this backend)'}",
        f"workspace  {ctx.workspace}",
        f"session    {st.get('session_id', '') or '(unsaved)'}",
        f"tokens     ↑{st.get('input_tokens', 0)} ↓{st.get('output_tokens', 0)}"
        f"  cached {st.get('cached_tokens', 0)}",
        f"cost       ${float(st.get('cost_usd', 0.0) or 0.0):.4f}",
    ]
    return CommandResult(output="\n".join(lines))


def _cost(ctx, _args: str) -> CommandResult:
    try:
        st = ctx.engine.get_status() or {}
    except Exception:
        return CommandResult(output="no cost recorded")
    return CommandResult(output=(
        f"${float(st.get('cost_usd', 0.0) or 0.0):.4f}  "
        f"↑{st.get('input_tokens', 0)} ↓{st.get('output_tokens', 0)} "
        f"cached {st.get('cached_tokens', 0)}"))


def _context(ctx, _args: str) -> CommandResult:
    try:
        return CommandResult(output=ctx.engine._build_context_status_block())
    except Exception:
        return CommandResult(output="context size is not available here")


def _compact(ctx, args: str) -> CommandResult:
    """The engine's own compactor, never a second cruder one."""
    try:
        ctx.engine._compact_history(force=True)
        return CommandResult(output=ctx.engine._compaction_status_line())
    except Exception as exc:
        return CommandResult(output=f"compaction failed: {exc}")


def _model(ctx, args: str) -> CommandResult:
    name = args.strip()
    if not name:
        current = getattr(getattr(ctx.engine, "client", None), "model", "?")
        return CommandResult(output=f"model is {current} — /model <name> to switch")
    client = getattr(ctx.engine, "client", None)
    switch = getattr(client, "switch_model", None)
    if not callable(switch):
        return CommandResult(output="this backend cannot switch model")
    try:
        switch(name)
    except Exception as exc:
        return CommandResult(output=f"could not switch: {exc}")
    return CommandResult(output=f"model → {name}")


def _mode_cmd(ctx, args: str) -> CommandResult:
    name = args.strip()
    if not name:
        return CommandResult(output=f"mode is {getattr(ctx.engine, 'mode', '?')}")
    try:
        ctx.engine._load_mode(name)
    except Exception as exc:
        return CommandResult(output=f"could not switch mode: {exc}")
    return CommandResult(output=f"mode → {name}")


def _permissions(ctx, args: str) -> CommandResult:
    name = args.strip()
    if not name:
        return CommandResult(output=(
            f"approval is {_mode(ctx.engine) or '(no gate)'} — "
            "/permissions plan|default|acceptEdits, or Shift+Tab"))
    if name == "bypassPermissions":
        return CommandResult(output=(
            "not from here. Unattended execution is reachable only by "
            "starting with --permission-mode bypassPermissions --unattended, "
            "so it stays something typed on purpose."))
    try:
        ctx.engine.set_kit_permission_mode(name)
    except Exception as exc:
        return CommandResult(output=f"could not switch: {exc}")
    return CommandResult(output=f"approval → {name}")


def _effort(ctx, args: str) -> CommandResult:
    level = args.strip()
    if not level:
        return CommandResult(output=f"effort is {getattr(ctx.engine, 'effort', '?')}")
    try:
        ctx.engine.effort = level
        client = getattr(ctx.engine, "client", None)
        if client is not None and hasattr(client, "effort"):
            client.effort = level
    except Exception as exc:
        return CommandResult(output=f"could not set effort: {exc}")
    return CommandResult(output=f"effort → {level}")


def _doctor(ctx, _args: str) -> CommandResult:
    try:
        from . import doctor
        return CommandResult(
            output=doctor.format_doctor(doctor.run_doctor(ctx.workspace)))
    except Exception as exc:
        return CommandResult(output=f"doctor failed: {exc}")


def _init(ctx, _args: str) -> CommandResult:
    try:
        from .project_init import init_project
        result = init_project(ctx.workspace)
        return CommandResult(output=str(result))
    except Exception as exc:
        return CommandResult(output=f"init failed: {exc}")


def _session(ctx, args: str) -> CommandResult:
    from . import session_store as ss
    parts = shlex.split(args) if args.strip() else ["ls"]
    action = parts[0]
    try:
        if action == "ls":
            rows = ss.list_sessions(limit=15, workspace=str(ctx.workspace))
            if not rows:
                return CommandResult(output="no sessions in this directory")
            return CommandResult(output="\n".join(
                f"  {str(r.get('updated_at', ''))[:16]}  "
                f"{str(r.get('session_id', ''))[:8]}  "
                f"{str(r.get('title', '') or '(untitled)')[:60]}"
                for r in rows))
        if action == "search" and len(parts) > 1:
            rows = ss.find_sessions(" ".join(parts[1:]), limit=15)
            if not rows:
                return CommandResult(output="nothing matched")
            return CommandResult(output="\n".join(
                f"  {str(r.get('session_id', ''))[:8]}  "
                f"{str(r.get('title', '') or '(untitled)')[:60]}" for r in rows))
    except Exception as exc:
        return CommandResult(output=f"session lookup failed: {exc}")
    return CommandResult(output="usage: /session ls | /session search <text>")


def _rewind(ctx, _args: str) -> CommandResult:
    """What the agent changed, and how to take it back."""
    try:
        from . import audit_log
        report = audit_log.build_changes_report(workspace=ctx.workspace)
        return CommandResult(output=audit_log.format_changes_report(report))
    except Exception as exc:
        return CommandResult(output=f"no change report available ({exc})")


def _mcp(ctx, _args: str) -> CommandResult:
    try:
        from . import mcp_client
        registry = mcp_client.get_registry(ctx.workspace)
        servers = mcp_client.effective_servers(registry)
        if not servers:
            return CommandResult(output="no MCP servers configured")
        return CommandResult(output="\n".join(f"  {s}" for s in servers))
    except Exception as exc:
        return CommandResult(output=f"MCP registry unavailable ({exc})")




def _trust(ctx, _args: str) -> CommandResult:
    """Read the offers before granting anything.

    Trusting a digest nobody has read is not consent, so this only ever
    SHOWS. Granting goes through hooks_editor / mcp_editor, the two call
    sites an existing test enumerates as the only ones allowed to.
    """
    try:
        from . import workspace_trust
        notices = workspace_trust.pending_notices(ctx.workspace)
        if not notices:
            return CommandResult(output="nothing is being withheld here")
        body = "\n".join(f"  {n}" for n in notices)
        return CommandResult(output=(
            f"{body}\n"
            "  Review them in the hooks / MCP editors before granting; "
            "only you can."))
    except Exception as exc:
        return CommandResult(output=f"trust state unavailable ({exc})")

# ---------------------------------------------------------------------------
# What the dashboard reaches, reached from here
#
# Each handler calls the module that already implements the thing. Where
# the dashboard's version is inline in a widget closure there is no
# function to call, and the terminal half is written out minimally and
# says so — importing the dashboard module would pull the notebook widget
# stack into a terminal that has none.
# ---------------------------------------------------------------------------

def _live_session(ctx) -> str:
    """The id the per-session stores are keyed by.

    The engine keeps one id and copies it onto the tool executor's
    permissions; reading both means a backend that only sets one still
    finds its own records instead of silently reporting an empty store.
    """
    engine = getattr(ctx, "engine", None)
    if engine is None:
        return ""
    sid = str(getattr(engine, "session_id", "") or "")
    if sid:
        return sid
    perms = getattr(engine, "kit_permissions", None)
    return str(getattr(perms, "task_session_id", "") or "")


def _tools(ctx, args: str) -> CommandResult:
    """The tool surface the model is offered, read from the catalogue.

    api_client owns both the catalogue and the filter that drops what
    this context could not execute, so both are called. The dashboard's
    /tools is a hand-kept table that has to be edited whenever a tool is
    added; a second copy here would drift the same way.
    """
    try:
        from . import api_client as ac
        role = ""
        try:
            role = str((ctx.engine.get_status() or {}).get("role", "") or "")
        except Exception:
            role = ""
        tools = ac.advertisable_tools(
            list(getattr(ac, "_DOC_TOOLS_OPENAI", None) or []),
            ac.ToolSurfaceContext(role=role))
    except Exception as exc:
        return CommandResult(output=f"tool catalogue unavailable ({exc})")
    query = args.strip().lower()
    rows = []
    for entry in tools:
        fn = entry.get("function", {}) if isinstance(entry, dict) else {}
        name = str(fn.get("name", "") or "")
        desc = " ".join(str(fn.get("description", "") or "").split())
        if query and query not in f"{name} {desc}".lower():
            continue
        rows.append(f"  {name:<24} {desc[:70]}")
    if not rows:
        return CommandResult(output=(
            f"no tool matches {query!r}" if query else "no tools advertised"))
    head = f"{len(rows)} tool(s)" + (f" matching {query!r}" if query else "")
    return CommandResult(output="\n".join([head, *rows]))


def _usage(ctx, _args: str) -> CommandResult:
    """Token and cost detail: /cost plus the rate that produced it.

    The rate comes from `pricing`, the one place model prices are written
    down, so the rate quoted here cannot contradict the total printed
    beside it. A model nobody has priced says so instead of showing a
    zero that reads like "free".
    """
    try:
        st = ctx.engine.get_status() or {}
    except Exception:
        return CommandResult(output="no usage recorded yet")
    inp = int(st.get("input_tokens", 0) or 0)
    out = int(st.get("output_tokens", 0) or 0)
    model = str(getattr(getattr(ctx.engine, "client", None), "model", "") or "")
    provider = str(st.get("provider", "") or "")
    lines = [
        f"model      {model or '?'}   provider {provider or '?'}",
        f"input      {inp:,} tokens",
        f"output     {out:,} tokens",
        f"cached     {int(st.get('cached_tokens', 0) or 0):,} tokens",
        f"messages   {len(getattr(ctx.engine, 'messages', []) or [])} in context",
    ]
    try:
        from . import pricing
        price = pricing.resolve(model, provider)
        if price.state == pricing.PRICED:
            lines.append(f"rate       ${price.input_per_mtok}/MTok in, "
                         f"${price.output_per_mtok}/MTok out")
            lines.append(
                f"cost       ${float(st.get('cost_usd', 0.0) or 0.0):.4f}")
        elif price.state == pricing.NON_BILLING:
            lines.append(f"rate       no per-token USD cost ({price.reason})")
            lines.append(f"spent      {inp + out:,} tokens")
        else:
            lines.append(f"rate       unmeasured — {price.reason}")
            lines.append(f"spent      {inp + out:,} tokens "
                         "(no rate to convert them)")
    except Exception as exc:
        lines.append(f"rate       unavailable ({exc})")
    return CommandResult(output="\n".join(lines))


def _export(ctx, _args: str) -> CommandResult:
    """The conversation as Markdown.

    The dashboard's export is inline in a closure over its own chat
    widget, so there is no function to call and this is the minimal
    terminal equivalent: the engine's own messages. It lands under the
    DELFIN state directory rather than in the workspace, so an export
    never turns up as a file the agent changed.
    """
    messages = getattr(getattr(ctx, "engine", None), "messages", None)
    if not isinstance(messages, list) or not messages:
        return CommandResult(output="nothing to export")
    lines = ["# DELFIN agent session", ""]
    for msg in messages:
        if not isinstance(msg, dict):
            continue
        role = str(msg.get("role", "") or "?")
        lines.append(f"### {role}")
        lines.append("")
        lines.append(_message_text(msg.get("content")))
        lines.append("")
    try:
        import time
        from . import state_paths
        target = (state_paths.ensure_dir(Path.home() / ".delfin" / "exports")
                  / f"session_{time.strftime('%Y%m%d_%H%M%S')}.md")
        state_paths.write_text(target, "\n".join(lines))
    except Exception as exc:
        return CommandResult(output=f"export failed: {exc}")
    return CommandResult(output=f"exported {len(messages)} message(s) → {target}")


def _message_text(content) -> str:
    """One message's text, whatever shape the backend stored it in.

    Content blocks arrive as a list on some backends and a plain string
    on others; rendering the list's repr would put JSON in the export
    where the user expects their conversation.
    """
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts = []
        for block in content:
            if isinstance(block, dict):
                parts.append(str(block.get("text", "") or ""))
            elif isinstance(block, str):
                parts.append(block)
        return "\n".join(p for p in parts if p.strip())
    return str(content or "")


def _memories(ctx, args: str) -> CommandResult:
    """The memory store, through memory_store's own reader."""
    scope = {"": "project", "global": "user", "all": "all"}.get(
        args.strip().lower())
    if scope is None:
        return CommandResult(output="usage: /memories [global|all]")
    try:
        from .memory_store import list_typed_memories
        mems = list_typed_memories(ctx.workspace, scope=scope)
    except Exception as exc:
        return CommandResult(output=f"memory store unavailable ({exc})")
    if not mems:
        return CommandResult(output=(
            "no memories stored — a line starting with # writes one"))
    lines: list[str] = []
    last_type = None
    for mem in mems:
        if mem.get("type") != last_type:
            last_type = mem.get("type")
            lines.append(f"  [{last_type}]")
        tag = " (global)" if mem.get("scope") == "user" else ""
        desc = str(mem.get("description") or (mem.get("body") or "")[:80] or "")
        lines.append(f"    {mem.get('name', '?')}{tag} — {desc.strip()}")
    lines.append("  /forget <name> deletes one")
    return CommandResult(output="\n".join(lines))


def _forget(ctx, args: str) -> CommandResult:
    """Delete one memory, by name.

    Both stores are searched and the deletion is reported with the path
    it removed: a "done" that names nothing cannot be told apart from a
    no-op, and this one deletes a file.
    """
    name = args.strip()
    if not name:
        return CommandResult(output="usage: /forget <name>  (see /memories)")
    try:
        from .memory_store import delete_typed_memory
        path = delete_typed_memory(ctx.workspace, name, scope="all")
    except Exception as exc:
        return CommandResult(output=f"could not delete: {exc}")
    if path is None:
        return CommandResult(output=(
            f"no memory named {name!r} — /memories lists them"))
    return CommandResult(output=f"deleted {name} → {path}")


def _undo(ctx, _args: str) -> CommandResult:
    """Drop the last turn from the context.

    The dashboard's /undo is inline in a closure that also rewinds its
    chat widget; this is the terminal half — the engine's message list
    only. It restores NO file, and says so, because the two undos are
    different acts: /undo-file is the one that touches the disk.
    """
    messages = getattr(getattr(ctx, "engine", None), "messages", None)
    if not isinstance(messages, list) or not messages:
        return CommandResult(output="nothing to undo")
    start = None
    for i in range(len(messages) - 1, -1, -1):
        msg = messages[i]
        if isinstance(msg, dict) and msg.get("role") == "user":
            start = i
            break
    if start is None:
        # No user turn to cut back to. Popping anyway would empty the
        # context, which is /clear's job and not what was asked for.
        return CommandResult(output="nothing to undo")
    removed = len(messages) - start
    del messages[start:]
    return CommandResult(output=(
        f"dropped {removed} message(s) from the context — no file was "
        "restored; /undo-file does that"))


def _undo_file(ctx, args: str) -> CommandResult:
    """Restore files from the undo journal — change_journal does the work.

    The rendering is minimal on purpose: the dashboard's formatters live
    in the widget module, and importing it here would pull the notebook
    stack into a terminal.
    """
    sid = _live_session(ctx)
    if not sid:
        return CommandResult(output=(
            "no active session — there is no change journal yet"))
    try:
        from . import change_journal as cj
    except Exception as exc:
        return CommandResult(output=f"undo journal unavailable ({exc})")
    scope = args.strip().lower() or "list"
    if scope in ("list", "ls"):
        return _undo_file_listing(cj, sid)
    if scope not in ("last", "turn", "session"):
        return CommandResult(output=(
            f"unknown scope {scope!r} — list | last | turn | session"))
    turn_seqs: list = []
    if scope == "turn":
        try:
            from .api_client import _doc_executor as executor
            turn_seqs = list(executor._turn_seqs_for(sid))
        except Exception:
            turn_seqs = []
        if not turn_seqs:
            return CommandResult(output=(
                "no file changes recorded in the current turn — try "
                "/undo-file last or /undo-file session"))
    try:
        res = cj.revert(sid, scope=scope, turn_seqs=turn_seqs,
                        workspace=Path(ctx.workspace))
    except Exception as exc:
        return CommandResult(output=f"undo failed: {exc}")
    _audit_undo(res, sid, scope)
    return CommandResult(output=_undo_file_result(scope, res))


def _undo_file_listing(cj, sid: str) -> CommandResult:
    """What is in the journal, and what of it can still be undone.

    An entry whose pre-image is gone is marked as such: listing it like
    any other offers an undo that cannot happen.
    """
    try:
        records = cj.list_changes(sid, last_n=20)
    except Exception as exc:
        return CommandResult(output=f"could not read the journal ({exc})")
    if not records:
        return CommandResult(output="no file changes recorded this session")
    lines = []
    for rec in records:
        if rec.get("dropped") or rec.get("truncated") or rec.get("lossy"):
            kind = "no restorable pre-image — cannot undo"
        elif rec.get("undone"):
            kind = "already undone"
        elif rec.get("deleted"):
            kind = "deleted"
        elif rec.get("created"):
            kind = "created"
        else:
            kind = "edited"
        lines.append(f"  [{rec.get('seq', '?')}] {rec.get('path', '?')} ({kind})")
    lines.append("  /undo-file last | turn | session")
    return CommandResult(output="\n".join(lines))


def _undo_file_result(scope: str, res: dict) -> str:
    """What the undo did, and for anything it did not do, why.

    A conflict and a missing pre-image are different refusals — one
    protects the user's own later edit, the other protects the file from
    being rewritten with content the journal never stored.
    """
    out = [f"undo ({scope})"]
    for key, label in (("reverted", "restored"),
                       ("conflicts", "conflict — NOT touched"),
                       ("skipped", "skipped")):
        for item in (res or {}).get(key) or []:
            if isinstance(item, dict):
                why = str(item.get("reason", "") or "").strip()
                out.append(f"  {label}: {item.get('path', '(no path)')}"
                           + (f" — {why}" if why else ""))
            else:
                out.append(f"  {label}: {item}")
    if len(out) == 1:
        out.append("  nothing to undo")
    return "\n".join(out)


def _audit_undo(res: dict, sid: str, scope: str) -> None:
    """An undo overwrites user files, so it belongs in /rewind's report.

    Without this the restore is invisible to the change report that is
    supposed to answer "what happened to my files".
    """
    try:
        from . import audit_log
        for path in (res or {}).get("reverted", []) or []:
            audit_log.append(audit_log.make_record(
                tool="undo_changes", decision="ok", path=str(path),
                session_id=sid, extra={"scope": scope, "source": "terminal"}))
    except Exception:
        pass


def _pending(ctx, _args: str) -> CommandResult:
    """The diff-approval queue, rendered by pending_changes itself."""
    sid = _live_session(ctx)
    if not sid:
        return CommandResult(output="no active session — nothing is pending")
    try:
        from . import pending_changes
        return CommandResult(output=pending_changes.render_pending(sid))
    except Exception as exc:
        return CommandResult(output=f"pending queue unavailable ({exc})")


def _approve(ctx, args: str) -> CommandResult:
    """Apply a staged diff. Named ids only, never "the newest one"."""
    arg = args.strip()
    sid = _live_session(ctx)
    if not sid or not arg:
        return CommandResult(output="usage: /approve <id|all>  (see /pending)")
    try:
        from . import pending_changes
        if arg.lower() == "all":
            res = pending_changes.approve_all(sid, workspace=ctx.workspace)
            lines = [f"applied {len(res.get('applied') or [])} change(s)"]
            for conflict in res.get("conflicts") or []:
                lines.append(f"  conflict #{conflict.get('id')} "
                             f"{conflict.get('path')}: {conflict.get('reason')}")
            for err in res.get("errors") or []:
                lines.append(f"  error: {err}")
            return CommandResult(output="\n".join(lines))
        res = pending_changes.approve(sid, arg, workspace=ctx.workspace)
    except Exception as exc:
        return CommandResult(output=f"approve failed: {exc}")
    if res.get("status") == "applied":
        return CommandResult(output=(
            f"applied #{res.get('id')} → {res.get('path', '')} "
            "(covered by the undo journal)"))
    return _change_verdict(arg, res)


def _reject(ctx, args: str) -> CommandResult:
    """Discard a staged diff without applying it."""
    arg = args.strip()
    sid = _live_session(ctx)
    if not sid or not arg:
        return CommandResult(output="usage: /reject <id|all>  (see /pending)")
    try:
        from . import pending_changes
        if arg.lower() == "all":
            res = pending_changes.reject_all(sid)
            return CommandResult(output=(
                f"rejected {len(res.get('rejected') or [])} change(s)"))
        res = pending_changes.reject(sid, arg)
    except Exception as exc:
        return CommandResult(output=f"reject failed: {exc}")
    return _change_verdict(arg, res)


def _change_verdict(change_id: str, res: dict) -> CommandResult:
    """The store's own word on an id it did or did not act on.

    An id that was not found reads as "not found", never as success:
    a staging queue that reports work it did not do is worse than one
    that refuses.
    """
    status = str((res or {}).get("status", "") or "unknown")
    reason = str((res or {}).get("reason", "") or "")
    return CommandResult(output=(
        f"#{change_id}: {status}" + (f" — {reason}" if reason else "")))


_GIT_READ_ONLY: dict[str, str] = {
    "status": "git status --short",
    "diff": "git diff --stat",
    "log": "git log --oneline -15",
    "branch": "git branch -a --no-color",
}


def _git(ctx, args: str) -> CommandResult:
    """Read-only git, through the same gate `!` uses.

    The four subcommands are a fixed table, not a pass-through: a /git
    that forwarded whatever followed it would be a second shell escape
    wearing a familiar name, and `!git ...` already exists for the rest.

    The gate is asked up front whether it would prompt. It parks its
    question for the MAIN thread, and a command handler runs ON the main
    thread — calling into a prompt from here would be the asker waiting
    for the answerer, i.e. itself. `!` does not have that problem because
    the loop runs it on a worker, so that is where this points.
    """
    sub = (args.strip().split() or [""])[0].lower()
    command = _GIT_READ_ONLY.get(sub)
    if command is None:
        return CommandResult(output="usage: /git status | diff | log | branch")
    run = getattr(getattr(ctx, "engine", None), "run_gated_bash", None)
    if not callable(run):
        return CommandResult(output=(
            "/git runs through the same gate the agent uses, and this "
            "backend has none"))
    perms = getattr(ctx.engine, "kit_permissions", None)
    try:
        would_ask = (perms is not None
                     and getattr(perms, "confirm_callback", None) is not None
                     and not perms.matches_bash_auto_allow(command))
    except Exception:
        would_ask = False
    if would_ask:
        return CommandResult(output=(
            f"the gate wants to approve this one — run `!{command}`, which "
            "the loop runs where the approval prompt can be answered"))
    try:
        out = run(command)
    except Exception as exc:
        return CommandResult(output=f"git {sub} failed: {exc}")
    from .repl import TerminalAgent
    body = "\n".join(f"  {line}" for line in TerminalAgent._shell_lines(out))
    return CommandResult(output=body or f"(no output from git {sub})")


def _agents(ctx, _args: str) -> CommandResult:
    """The subagent presets, from the registry that defines them."""
    try:
        from . import subagents
        presets = subagents.list_subagents()
    except Exception as exc:
        return CommandResult(output=f"subagent presets unavailable ({exc})")
    if not presets:
        return CommandResult(output="no subagent presets registered")
    lines = [
        f"  {str(p.get('subagent_type') or p.get('name') or '?'):<18} "
        f"{str(p.get('description') or '')[:60]}"
        for p in presets
    ]
    try:
        finished = subagents.list_finished(last_n=5)
    except Exception:
        finished = []
    if finished:
        lines.append("  recently finished (resume via resume_id):")
        for rec in finished:
            lines.append(
                f"    {str(rec.get('sa_id', '?')):<10} "
                f"{str(rec.get('subagent_type', '?')):<16} "
                f"{str(rec.get('description') or '')[:40]}")
    return CommandResult(output="\n".join(lines))


def _skills(ctx, args: str) -> CommandResult:
    """Discovered skills, and one skill's body on request."""
    try:
        from . import skills
    except Exception as exc:
        return CommandResult(output=f"skills unavailable ({exc})")
    name = args.strip()
    if name:
        try:
            skill = skills.get_skill(name, ctx.workspace)
        except Exception as exc:
            return CommandResult(output=f"could not read {name}: {exc}")
        if skill is None:
            return CommandResult(output=(
                f"no skill named {name!r} — /skills lists them"))
        return CommandResult(output=(
            f"{skill.name} — {skill.description}\n"
            f"source {skill.source}\n\n{skill.body}"))
    try:
        found = skills.discover_skills(ctx.workspace)
    except Exception as exc:
        return CommandResult(output=f"could not list skills ({exc})")
    if not found:
        return CommandResult(output=(
            "no skills discovered — drop a SKILL.md under "
            "~/.delfin/skills/<name>/ or <workspace>/.delfin/skills/"))
    lines = [f"  /{s.name:<22} {(s.description or '')[:60]}" for s in found]
    lines.append("  /skills <name> shows the body")
    return CommandResult(output="\n".join(lines))


def _hooks(ctx, _args: str) -> CommandResult:
    """Which hooks are registered, and what the load could not do.

    List only: adding, removing and dry-running a hook edit settings the
    user owns, and this surface never does. The warnings are printed
    because a hook withheld for lack of trust looks exactly like a hook
    that had nothing to say — see /trust grant hooks.
    """
    try:
        from . import hooks_editor
        rows = hooks_editor.list_hooks(ctx.workspace)
        warnings = hooks_editor.hook_warnings(ctx.workspace)
    except Exception as exc:
        return CommandResult(output=f"hooks unavailable ({exc})")
    if not rows and not warnings:
        return CommandResult(output="no hooks registered")
    lines: list[str] = []
    event = ""
    for row in rows:
        if row.get("event") != event:
            event = str(row.get("event", "") or "")
            lines.append(f"  {event}:")
        lines.append(f"    [{row.get('index')}] matcher="
                     f"{str(row.get('matcher', '') or ''):<20} "
                     f"{str(row.get('command', ''))[:60]}")
        lines.append(f"        source {row.get('source') or 'unknown'}")
    for note in warnings:
        lines.append(f"  ! {note}")
    return CommandResult(output="\n".join(lines))


def _bash(ctx, args: str) -> CommandResult:
    """Background shell jobs: list them, or kill one by id.

    Kill is the only mutation, and it needs the id: a /bash that killed
    "the last one" would reach a job the user is not looking at.
    """
    try:
        from . import bash_jobs
        registry = bash_jobs.get_registry()
    except Exception as exc:
        return CommandResult(output=f"job registry unavailable ({exc})")
    parts = args.strip().split()
    action = parts[0].lower() if parts else "ls"
    if action == "kill":
        if len(parts) < 2:
            return CommandResult(output="usage: /bash kill <job_id>")
        try:
            ok, note = registry.kill(parts[1])
        except Exception as exc:
            return CommandResult(output=f"kill failed: {exc}")
        return CommandResult(output=f"{parts[1]}: {note}" if ok
                             else f"not killed — {note}")
    if action not in ("ls", "jobs"):
        return CommandResult(output="usage: /bash [ls] | /bash kill <job_id>")
    try:
        jobs = sorted(registry.list_jobs(), key=lambda j: j.started_at,
                      reverse=True)
    except Exception as exc:
        return CommandResult(output=f"could not list jobs ({exc})")
    if not jobs:
        return CommandResult(output="no background bash jobs")
    lines = []
    for job in jobs[:20]:
        try:
            st = job.status_dict()
        except Exception:
            continue
        flag = "run " if st.get("running") else "done"
        lines.append(f"  {flag} {st.get('job_id')}  exit={st.get('exit_code')}  "
                     f"{float(st.get('elapsed_s', 0.0) or 0.0):>7.1f}s  "
                     f"{str(st.get('command', ''))[:50]}")
    lines.append("  /bash kill <job_id>")
    return CommandResult(output="\n".join(lines))


def _attention(ctx, args: str) -> CommandResult:
    """The attention inbox, rendered by the module that owns it.

    Read-only here. Answering and dismissing resolve an item the agent
    is waiting on, and those belong on the surface that parked it.
    """
    try:
        from . import attention
        kind = args.strip().lower()
        if kind and kind not in attention.ATTENTION_KINDS:
            return CommandResult(output=(
                "usage: /attention [" + "|".join(
                    sorted(attention.ATTENTION_KINDS)) + "]"))
        return CommandResult(output=attention.render_inbox(kind or None))
    except Exception as exc:
        return CommandResult(output=f"attention inbox unavailable ({exc})")


def _plans(ctx, args: str) -> CommandResult:
    """Saved plans, and one plan's body on request."""
    try:
        from . import memory_store
    except Exception as exc:
        return CommandResult(output=f"plan store unavailable ({exc})")
    name = args.strip()
    if name:
        try:
            rec = memory_store.get_plan(ctx.workspace, name)
        except Exception as exc:
            return CommandResult(output=f"could not read {name}: {exc}")
        if rec is None:
            return CommandResult(output=(
                f"no plan matching {name!r} — /plans lists them"))
        return CommandResult(output=(
            f"{rec.get('name')} — {rec.get('description', '')}\n"
            f"source {rec.get('path')}\n\n{rec.get('body', '')}"))
    try:
        plans = memory_store.list_plans(ctx.workspace)
    except Exception as exc:
        return CommandResult(output=f"could not list plans ({exc})")
    if not plans:
        return CommandResult(output=(
            "no saved plans — one is written when a plan-mode plan is "
            "approved"))
    import time
    lines = []
    for plan in plans[:25]:
        stamp = time.strftime("%Y-%m-%d %H:%M",
                              time.localtime(plan.get("created_at", 0) or 0))
        lines.append(f"  {stamp}  {str(plan.get('name', '?')):<28} "
                     f"{str(plan.get('description') or '')[:40]}")
    lines.append("  /plans <name> shows the body")
    return CommandResult(output="\n".join(lines))


def _commands(ctx, args: str) -> CommandResult:
    """User-defined slash commands, from the same discovery the router uses.

    Listing them from a second source would let /commands advertise a
    command the router does not have, or hide one it does.
    """
    try:
        from . import slash_commands
    except Exception as exc:
        return CommandResult(output=f"command store unavailable ({exc})")
    name = args.strip()
    if name:
        try:
            tpl = slash_commands.get_command(name, ctx.workspace)
        except Exception as exc:
            return CommandResult(output=f"could not read {name}: {exc}")
        if tpl is None:
            return CommandResult(output=(
                f"no command named {name!r} — /commands lists them"))
        return CommandResult(output=(
            f"/{tpl.name} — {tpl.description}\n"
            f"source {tpl.source}\n\n{tpl.body}"))
    try:
        found = slash_commands.discover_commands(ctx.workspace)
    except Exception as exc:
        return CommandResult(output=f"could not list commands ({exc})")
    if not found:
        return CommandResult(output=(
            "no custom commands — a markdown file in ~/.delfin/commands/ "
            "or <workspace>/.delfin/commands/ becomes one"))
    lines = []
    for tpl in found:
        hint = f" {tpl.argument_hint}" if tpl.argument_hint else ""
        lines.append(f"  /{tpl.name}{hint:<20} {(tpl.description or '')[:60]}")
    lines.append("  /commands <name> shows the body")
    return CommandResult(output="\n".join(lines))


def _trace(ctx, args: str) -> CommandResult:
    """The tool calls this session made, as tool_trace records them."""
    engine = getattr(ctx, "engine", None)
    session = getattr(engine, "trace_session", None)
    sid = ""
    try:
        sid = str(session() or "") if callable(session) else _live_session(ctx)
    except Exception:
        sid = _live_session(ctx)
    if not sid:
        return CommandResult(output="no session yet — nothing has been traced")
    try:
        from . import tool_trace
        entries = tool_trace.read(sid)
        if not entries:
            return CommandResult(output="no tool calls recorded this session")
        limit = int(args.strip()) if args.strip().isdigit() else 30
        return CommandResult(output=(
            tool_trace.format_summary(entries, limit=limit)
            + f"\nfull trace: {tool_trace.trace_path(sid)}"))
    except Exception as exc:
        return CommandResult(output=f"trace unavailable ({exc})")


def _clear(ctx, _args: str) -> CommandResult:
    return CommandResult(output="starting a fresh conversation", clear=True)


def _exit(ctx, _args: str) -> CommandResult:
    return CommandResult(quit=True)


def _help(ctx, args: str) -> CommandResult:
    from . import help_gen
    return CommandResult(output=help_gen.generate_help(
        palette_rows(), search=args.strip()))


def _mode(engine) -> str:
    from .repl import permission_mode
    return permission_mode(engine)


BUILTINS: dict[str, ReplCommand] = {
    c.name: c for c in [
        ReplCommand("/help", "session", "List these commands", _help, True),
        ReplCommand("/status", "session", "Model, mode, tokens, cost", _status),
        ReplCommand("/cost", "session", "What this session has cost", _cost),
        ReplCommand("/usage", "session", "Token and cost detail, with the rate",
                    _usage),
        ReplCommand("/context", "session", "How much window is left", _context),
        ReplCommand("/compact", "session", "Summarise history to free space",
                    _compact, True),
        ReplCommand("/export", "session", "Write this conversation as Markdown",
                    _export),
        ReplCommand("/undo", "session", "Drop the last turn from the context",
                    _undo),
        ReplCommand("/trace", "session", "Tool calls made this session",
                    _trace, True),
        ReplCommand("/clear", "session", "Start a fresh conversation", _clear),
        ReplCommand("/exit", "session", "Leave", _exit),
        ReplCommand("/quit", "session", "Leave", _exit),
        ReplCommand("/model", "setup", "Show or switch the model", _model, True),
        ReplCommand("/mode", "setup", "Show or switch the agent mode",
                    _mode_cmd, True),
        ReplCommand("/permissions", "setup", "Show or switch the approval "
                    "posture", _permissions, True),
        ReplCommand("/effort", "setup", "Show or set reasoning effort",
                    _effort, True),
        ReplCommand("/doctor", "setup", "Check the setup", _doctor),
        ReplCommand("/init", "setup", "Scaffold AGENTS.md and .delfin/", _init),
        ReplCommand("/mcp", "setup", "Configured MCP servers", _mcp),
        ReplCommand("/trust", "setup", "What this folder offers and withholds",
                    _trust),
        ReplCommand("/tools", "setup", "Tools the model is offered",
                    _tools, True),
        ReplCommand("/hooks", "setup", "Hooks registered for this workspace",
                    _hooks),
        ReplCommand("/agents", "setup", "Subagent presets", _agents),
        ReplCommand("/skills", "setup", "Discovered skills", _skills, True),
        ReplCommand("/commands", "setup", "User-defined slash commands",
                    _commands, True),
        ReplCommand("/session", "history", "List or search past sessions",
                    _session, True),
        ReplCommand("/rewind", "history", "What the agent changed here",
                    _rewind),
        ReplCommand("/memories", "history", "Stored memories", _memories, True),
        ReplCommand("/forget", "history", "Delete one memory by name",
                    _forget, True),
        ReplCommand("/plans", "history", "Saved plans", _plans, True),
        ReplCommand("/undo-file", "changes", "Restore files from the undo "
                    "journal", _undo_file, True),
        ReplCommand("/pending", "changes", "Diffs waiting for approval",
                    _pending),
        ReplCommand("/approve", "changes", "Apply a staged diff",
                    _approve, True),
        ReplCommand("/reject", "changes", "Discard a staged diff",
                    _reject, True),
        ReplCommand("/git", "workspace", "status | diff | log | branch",
                    _git, True),
        ReplCommand("/bash", "workspace", "Background jobs: list or kill one",
                    _bash, True),
        ReplCommand("/attention", "workspace", "The attention inbox",
                    _attention, True),
    ]
}


def palette_rows() -> list[tuple[str, str, str, bool]]:
    """The rows help_gen renders. Derived, so help cannot drift."""
    seen: set[str] = set()
    rows = []
    for cmd in BUILTINS.values():
        if cmd.summary in seen and cmd.name in ("/quit",):
            continue
        seen.add(cmd.summary)
        rows.append((cmd.category, cmd.name, cmd.summary, cmd.takes_args))
    return rows


# ---------------------------------------------------------------------------
# The route
# ---------------------------------------------------------------------------

@dataclass
class ReplContext:
    engine: Any = None
    workspace: Path = field(default_factory=Path.cwd)
    session_id: str = ""


def dispatch(text: str, ctx: ReplContext) -> CommandResult:
    """Builtin, then user command, then skill, then subagent, then the model.

    The dashboard's order, kept because a `/deploy` that means one thing
    in the browser and another in a terminal is worse than not having it.
    """
    raw = (text or "").strip()
    if not looks_like_command(raw):
        return CommandResult(handled=False)

    head, _, args = raw.partition(" ")
    name = head.lower()
    args = args.strip()

    cmd = BUILTINS.get(name)
    if cmd is not None:
        try:
            return cmd.handler(ctx, args)
        except Exception as exc:
            return CommandResult(output=f"{name} failed: {exc}")

    bare = name.lstrip("/")
    try:
        from . import slash_commands
        tpl = slash_commands.get_command(bare, ctx.workspace)
        if tpl is not None:
            return CommandResult(
                prompt=slash_commands.expand_template(tpl.body, args))
    except Exception:
        pass

    try:
        from . import skills
        skill = skills.get_skill(bare, ctx.workspace)
        if skill is not None:
            return CommandResult(
                prompt=skills.render_skill_invocation(skill, args))
    except Exception:
        pass

    try:
        from . import slash_commands
        expanded = slash_commands.builtin_subagent_command(bare, args)
        if expanded:
            return CommandResult(prompt=expanded)
    except Exception:
        pass

    # Unknown: the literal text goes to the model, which is what the
    # dashboard does. Refusing here would make a typo unanswerable.
    return CommandResult(handled=False)
