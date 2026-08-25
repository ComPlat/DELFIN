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
        ReplCommand("/context", "session", "How much window is left", _context),
        ReplCommand("/compact", "session", "Summarise history to free space",
                    _compact, True),
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
        ReplCommand("/session", "history", "List or search past sessions",
                    _session, True),
        ReplCommand("/rewind", "history", "What the agent changed here",
                    _rewind),
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
