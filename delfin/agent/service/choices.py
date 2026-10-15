"""What a window offers to pick from: commands, paths, models, efforts, modes.

The window completes ``/commands`` and ``@paths`` as the user types and
opens pickers for the model, the effort and the agent mode. It asks the
session for the lists, so they are the session's: its workspace's own
commands and skills, the paths under its workspace, the models its
provider serves. Picking one sends the command the REPL would get
(``/model <name>``), so a choice is applied exactly the way a typed one is.
"""

from __future__ import annotations

import re
import threading
import time
from pathlib import Path
from typing import Any

#: Reasoning effort as `delfin-agent chat --effort` takes it ("" is the
#: role's own default).
EFFORTS = ("", "low", "medium", "high", "xhigh")

#: How long a provider's model list is reused before it is asked again.
#: Listing the KIT models takes seconds; a picker opened twice in a minute
#: must not wait twice.
_MODELS_TTL_S = 600.0
_models_cache: dict[str, tuple[float, list[dict[str, str]]]] = {}
_models_lock = threading.Lock()


def command_list(workspace: Path) -> list[dict[str, Any]]:
    """Every ``/command`` a line may start with, in the order
    ``repl_commands.dispatch`` tries them: the builtins, the user's and the
    workspace's own commands, skills, then the subagent shortcuts. A name
    an earlier kind already takes is left out, because typing it runs the
    earlier one."""
    from .. import repl_commands as rc

    out: list[dict[str, Any]] = []
    seen: set[str] = set()

    def add(name: str, summary: str, category: str, takes_args: bool,
            source: str) -> None:
        key = name.lower()
        if not name.startswith("/") or key in seen:
            return
        seen.add(key)
        out.append({"name": name, "summary": " ".join(str(summary or "").split()),
                    "category": category, "takesArgs": bool(takes_args),
                    "source": source})

    for category, name, summary, takes_args in rc.palette_rows():
        add(name, summary, category, takes_args, "builtin")
    try:
        from .. import slash_commands
        for tpl in slash_commands.discover_commands(workspace):
            add(tpl.slash, tpl.description or tpl.argument_hint, "user", True, "user")
    except Exception:
        pass
    try:
        from .. import skills
        for skill in skills.discover_skills(workspace):
            add("/" + skill.name, skill.description, "skills", True, "skill")
    except Exception:
        pass
    try:
        from .. import slash_commands
        for key, (preset, gloss) in slash_commands._BUILTIN_SUBAGENT_COMMANDS.items():
            add("/" + key, f"{gloss} ({preset} subagent)", "subagents", True,
                "subagent")
    except Exception:
        pass
    return out


def complete_path(fragment: str, workspace: Path, limit: int = 50) -> list[str]:
    """Paths under the workspace that ``@fragment`` may complete to.

    The REPL completer's rules (``repl_commands.complete_path``): nothing
    outside the workspace, dot files only when asked for, folders with a
    trailing slash, at most ``limit``. One difference: a fragment that
    ends in a slash lists that folder -- the REPL's completer reads
    "delfin/" as the name "delfin" and offers the folder itself again, so
    completion could never go into one. What the user typed before the
    last slash is kept as typed.
    """
    fragment = str(fragment or "").lstrip("@")
    head, slash, stem = fragment.rpartition("/")
    try:
        root = workspace.resolve()
        base = (root / head).resolve() if slash else root
        if base != root and root not in base.parents:
            return []
        out: list[str] = []
        for entry in sorted(base.iterdir(), key=lambda e: e.name):
            if entry.name.startswith(".") and not stem.startswith("."):
                continue
            if stem and not entry.name.startswith(stem):
                continue
            name = f"{head}{slash}{entry.name}"
            out.append(f"{name}/" if entry.is_dir() else name)
            if len(out) >= limit:
                break
        return out
    except Exception:
        return []


def _first(seq: Any, default: str = "") -> str:
    return str(seq[0]) if seq else default


def _model_detail(description: str) -> str:
    """Where a model runs and what it costs, from a ``/v1/models``
    description: "KIT", "Azure (EU) · $1.38 in, $11 out /1M"."""
    text = str(description or "")
    host = re.search(r"Host:\s*([^•\n]+)", text)
    price = re.search(r"Price \(1M tokens\):\s*([^•\n]+)", text)
    parts = []
    if host:
        parts.append(host.group(1).strip())
    if price:
        rates = [r.strip() for r in re.sub(r"[⚠️]", "", price.group(1)).split(",")]
        rates = [r for r in rates if r.endswith((" in", " out"))] or rates[:2]
        parts.append(", ".join(rates) + " per 1M tokens")
    if parts:
        return " · ".join(parts)
    alias = re.search(r"Currently uses ([^.\n]+)", text)
    if alias:
        return f"alias of {alias.group(1).strip()}"
    return _first(text.split("\n")).strip()[:80]


def _served_models(client: Any) -> list[dict[str, str]]:
    """The chat models an OpenAI-compatible endpoint lists (KIT, OpenAI,
    a local vLLM or Ollama), or nothing when it cannot be asked."""
    base = str(getattr(client, "_base_url", "") or "")
    if not base:
        return []
    with _models_lock:
        hit = _models_cache.get(base)
        if hit and time.monotonic() - hit[0] < _MODELS_TTL_S:
            return list(hit[1])
    from .. import model_capabilities as mc
    entries = mc._fetch_openai_models(base, str(getattr(client, "_api_key", "") or "")) or []
    found: list[dict[str, str]] = []
    for entry in entries:
        info = entry.get("info") or {}
        model_id = str(info.get("id") or entry.get("id") or "")
        if not model_id or mc.nonchat_reason(model_id):
            continue
        meta = info.get("meta") or {}
        found.append({"value": model_id,
                      "detail": _model_detail(meta.get("description") or
                                              entry.get("description") or "")})
    if found:
        with _models_lock:
            _models_cache[base] = (time.monotonic(), list(found))
    return found


def _roster_fallback(engine: Any, client: Any) -> list[dict[str, str]]:
    """The curated KIT roster (model_capabilities._KIT_PREFERENCE, best
    first), for a KIT session whose ``/v1/models`` cannot be read: a
    picker that offers only the current model offers no switch at all.
    Nothing for another provider, whose names DELFIN does not curate."""
    provider = str(getattr(engine, "provider", "") or getattr(client, "_provider", "") or "")
    if provider.lower() != "kit":
        return []
    from .. import model_capabilities as mc
    return [{"value": m, "detail": "recommended · the model list could not be read"}
            for m in mc._KIT_PREFERENCE]


def choices(engine: Any, what: str) -> dict[str, Any]:
    """``{current, choices: [{value, detail?}]}`` for a picker."""
    client = getattr(engine, "client", None)
    if what == "model":
        current = str(getattr(client, "model", "") or "")
        listed = _served_models(client) if client is not None else []
        if not listed:
            listed = _roster_fallback(engine, client)
        if current and not any(c["value"] == current for c in listed):
            listed.insert(0, {"value": current, "detail": "current"})
        return {"current": current, "choices": listed}
    if what == "effort":
        current = str(getattr(client, "effort", "") or getattr(engine, "effort", "") or "")
        return {"current": current,
                "choices": [{"value": e, "detail": "the role's default" if not e else ""}
                            for e in EFFORTS]}
    if what == "mode":
        current = str(getattr(engine, "mode", "") or "")
        try:
            modes = [str(m) for m in (engine.available_modes() or [])]
        except Exception:
            modes = []
        if current and current not in modes:
            modes.insert(0, current)
        return {"current": current, "choices": [{"value": m} for m in modes]}
    raise ValueError("what must be model, effort, mode or args")


# -- what goes after a /command ------------------------------------------------

#: The fixed words a builtin takes (its handler's own grammar in
#: repl_commands), each with what it does.
_WORDS: dict[str, list[tuple[str, str]]] = {
    "/permissions": [("plan", "read and plan; every write and command asked"),
                     ("default", "ask before each write and each command"),
                     ("acceptEdits", "writes go through; commands still asked"),
                     ("rules", "the permission rules in force")],
    "/git": [("status", ""), ("diff", ""), ("log", ""), ("branch", ""),
             ("role", "maintainer or contributor here")],
    "/undo-file": [("list", "what the journal can restore"), ("last", "the last change"),
                   ("turn", "the last turn's changes"), ("session", "every change of this session")],
    "/memories": [("global", "the memories of every project"), ("all", "this project's and the global ones")],
    "/memory": [("private", "this machine only"), ("swarm", "this run's findings"),
                ("project", "the repository's .delfin/knowledge"),
                ("global", "every DELFIN user"), ("all", "every layer")],
    "/swarm": [("runs", "the runs on this machine"), ("harvest", "a run's findings with evidence"),
               ("propose", "keep one as your private lesson"), ("close", "end and archive a run")],
    "/jobs": [("ls", "the watched jobs"), ("watch on", "start the monitor"),
              ("watch off", "stop the monitor"), ("watch status", "is the monitor running")],
    "/trust": [("grant hooks", ""), ("grant mcp", ""), ("revoke hooks", ""), ("revoke mcp", "")],
    "/session": [("ls", "this folder's sessions")],
}

#: What a builtin takes when its values cannot be listed.
_USAGE: dict[str, str] = {
    "/editor": "/editor [text to start from]",
    "/compact": "/compact [what the summary must keep]",
    "/trace": "/trace [how many calls, default 30]",
    "/bug": "/bug <what went wrong>",
    "/session": "/session ls | /session search <words>",
    "/forget": "/forget <memory name>",
    "/approve": "/approve <id> | all  (/pending lists them)",
    "/reject": "/reject <id> | all  (/pending lists them)",
    "/bash": "/bash [ls] | <job id> | kill <job id>",
    "/attention": "/attention [clear [kind]]",
    "/tools": "/tools [search words] | groups | enable <group>",
    "/skills": "/skills [name] | proposals",
    "/commands": "/commands [name]",
    "/plans": "/plans [name]",
    "/help": "/help [search words]",
    "/share": "/share <lesson> [project|global] [scope=..] [topic=..] [as <handle>] | confirm <dialog> <code>",
    "/unshare": "/unshare <entry> [as <handle>]",
    "/lessons": "/lessons [show <id> | edit <id> field=text | delete <id>]",
    "/knowledge": "/knowledge [list | search | lint | dupes | stats | queue | accept | reject | merge | move | retire | rename | index | topic add]",
}


def _listed(engine: Any, workspace: Path, command: str) -> list[dict[str, str]]:
    """Values that depend on the session or its workspace."""
    if command == "/help":
        from .. import repl_commands as rc
        return [{"value": name.lstrip("/"), "detail": summary}
                for _cat, name, summary, _args in rc.palette_rows()]
    if command == "/skills":
        from .. import skills
        return [{"value": k.name, "detail": " ".join(str(k.description or "").split())[:80]}
                for k in skills.discover_skills(workspace)] + [{"value": "proposals", "detail": "skills proposed to you"}]
    if command == "/commands":
        from .. import slash_commands
        return [{"value": t.slash.lstrip("/"), "detail": str(t.description or "")[:80]}
                for t in slash_commands.discover_commands(workspace)]
    if command in ("/forget",):
        from ..memory_store import list_typed_memories
        return [{"value": str(m["name"]), "detail": str(m.get("description") or "")[:80]}
                for m in list_typed_memories(workspace, scope="all")]
    if command == "/plans":
        from .. import memory_store
        return [{"value": str(p["name"]), "detail": str(p.get("description") or "")[:80]}
                for p in memory_store.list_plans(workspace)[:25]]
    if command in ("/approve", "/reject"):
        from .. import pending_changes
        sid = str(getattr(engine, "session_id", "") or "")
        rows = [{"value": str(r.get("seq") or r.get("id")),
                 "detail": str(r.get("path") or r.get("file") or "")[:80]}
                for r in pending_changes.list_pending(sid)]
        return rows + ([{"value": "all", "detail": "every staged diff"}] if rows else [])
    if command == "/bash":
        from .. import bash_jobs
        jobs = bash_jobs.get_registry().list_jobs(include_finished=False)
        return [{"value": "ls", "detail": "the background jobs"}] + [
            {"value": f"kill {j.job_id}", "detail": " ".join(str(j.command).split())[:60]} for j in jobs]
    if command == "/attention":
        from .. import attention
        return [{"value": "clear", "detail": "empty the inbox"}] + [
            {"value": f"clear {k}", "detail": ""} for k in attention.ATTENTION_KINDS]
    if command == "/tools":
        from ..toolsets import GROUPS
        return [{"value": "groups", "detail": "the optional tool groups"}] + [
            {"value": f"enable {g}", "detail": f"offer the {g} tools from the next turn"} for g in GROUPS]
    return []


def argument_choices(engine: Any, workspace: Path, command: str) -> dict[str, Any]:
    """What may follow ``command`` (``/model``) on a typed line:
    ``{current, choices: [{value, detail?}], usage}``. A value is the text
    after the command and a space, as the user would type it; ``usage`` a
    line saying what the command takes when its values cannot all be
    listed (or "" when they can). A command the session does not run
    gets nothing. Values are read when asked: memories, plans, staged
    diffs and jobs change while the session runs."""
    command = str(command or "").strip().split(" ")[0].lower()
    if command in ("/model", "/effort", "/mode"):
        got = choices(engine, command[1:])
        listed = [c for c in got["choices"] if c.get("value")]
        return {"current": got["current"], "choices": listed, "usage": ""}
    known = {c["name"].lower(): c for c in command_list(workspace)}
    entry = known.get(command)
    if entry is None or not entry.get("takesArgs"):
        return {"current": "", "choices": [], "usage": ""}
    current = ""
    if command == "/permissions":
        try:
            current = str(getattr(engine.kit_permissions, "mode", "") or "")
        except Exception:
            current = ""
    values = [{"value": v, "detail": d} for v, d in _WORDS.get(command, [])]
    try:
        values += _listed(engine, workspace, command)
    except Exception:
        pass                                     # a list that cannot be read: the usage says it
    usage = _USAGE.get(command, "")
    if entry.get("source") == "user":
        from .. import slash_commands
        try:
            tpl = slash_commands.get_command(command.lstrip("/"), workspace)
            hint = str(getattr(tpl, "argument_hint", "") or "").strip()
        except Exception:
            hint = ""
        usage = f"{entry['name']} {hint or '[arguments]'}"
    elif entry.get("source") == "skill":
        usage = f"{entry['name']} [what to do with it]"
    elif entry.get("source") == "subagent":
        usage = f"{entry['name']} <the task>"
    elif not values and not usage:
        usage = f"{entry['name']} [arguments]"
    clean = [{"value": " ".join(str(c["value"]).split()), "detail": " ".join(str(c.get("detail") or "").split())}
             for c in values if str(c.get("value") or "").strip()]
    return {"current": current, "choices": clean, "usage": usage}
