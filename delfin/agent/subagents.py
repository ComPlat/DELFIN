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

  - max 30 tool calls per sub-agent
  - max 60 seconds wall-clock per sub-agent
  - max 8000 tokens output

Failures don't propagate: a crash inside the sub-agent returns an
error string, never an exception.
"""

from __future__ import annotations

import dataclasses
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .api_client import KitToolPermissions, OpenAIClient


# Per-run hard caps. These are FALLBACK defaults; settings
# ["agent"]["subagents"] overrides them (see _subagent_limits). Wall-clock
# raised 60→300s: 60s truncated real exploration/research runs (especially on
# slower KIT/Qwen models) before they could report back — subagents were too
# short-leashed to be useful for anything but trivial lookups.
_MAX_TOOL_CALLS = 40
_MAX_WALL_S = 300.0
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


def _running_update(sa_id: str, entry: dict | None) -> None:
    """Maintain the live per-subagent status file (entry=None removes).

    File-based so the dashboard can render a Claude-Code-style live drill-down
    (name · task · steps · status) without sharing memory with the worker
    thread."""
    try:
        _RUNNING_DIR.mkdir(parents=True, exist_ok=True)
        f = _RUNNING_DIR / f"{sa_id}.json"
        if entry is None:
            try:
                f.unlink()
            except FileNotFoundError:
                pass
        else:
            f.write_text(json.dumps(entry), encoding="utf-8")
    except Exception:
        pass


def read_running() -> dict:
    """Live registry: {id: {type, description, started_at, actions, last_action}}."""
    out: dict = {}
    try:
        for f in _RUNNING_DIR.glob("*.json"):
            try:
                out[f.stem] = json.loads(f.read_text(encoding="utf-8"))
            except Exception:
                continue
    except Exception:
        pass
    return out


# Finished-subagent sessions (Claude-Code SendMessage analog): each run
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


def _save_subagent_session(
    sa_id: str,
    *,
    subagent_type: str,
    description: str,
    messages: list[dict],
    interactions: list[dict],
    error: str = "",
) -> None:
    """Persist a finished subagent conversation for later resumption.

    Stores the logical user/assistant conversation AND the tool
    interactions WITH their (trimmed) outputs — so a resumed subagent
    sees what it actually read, not just its own conclusions. Best-effort,
    never raises.
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
        }
        _SESSIONS_DIR.mkdir(parents=True, exist_ok=True)
        path = _SESSIONS_DIR / f"{sa_id}.json"
        path.write_text(json.dumps(record, ensure_ascii=False),
                        encoding="utf-8")
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
        path = _TELEMETRY_PATH
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(record, ensure_ascii=False) + "\n")
        # Best-effort trim — only every 50 writes to avoid I/O thrash
        try:
            stat = path.stat()
            if stat.st_size > 1_000_000:  # ~5k records of 200 bytes
                lines = path.read_text(encoding="utf-8").splitlines()
                if len(lines) > _TELEMETRY_MAX_LINES:
                    tail = lines[-_TELEMETRY_MAX_LINES:]
                    path.write_text("\n".join(tail) + "\n", encoding="utf-8")
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

    def to_payload(self) -> dict:
        return {
            "subagent_type": self.subagent_type,
            "description": self.description,
            "sa_id": self.sa_id,
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


def _derive_perms(parent_perms, mode: str, workspace=None):
    """Clone parent_perms with the sub-agent's mode + optional workspace.

    Both fields are optional — the caller passes ``workspace=None`` when
    no isolation is requested and the parent's workspace inherits.
    """
    if parent_perms is None:
        return None
    try:
        if workspace is not None:
            return dataclasses.replace(parent_perms, mode=mode, workspace=workspace)
        return dataclasses.replace(parent_perms, mode=mode)
    except Exception:
        return parent_perms


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

    system_prompt = (
        preset.system_prompt
        + f"\n\nWorkspace: {sub_perms.workspace if sub_perms else '(none)'}"
        + f"\nTask label: {description}"
    )
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
    # watch its activity in the chat window, Claude-Code style.
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
    if worktree_info is not None:
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
    elapsed_s = time.monotonic() - t0
    # Persist the conversation so the parent can resume this subagent
    # later via ``resume_id`` (Claude-Code SendMessage analog). Store the
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
    )
    # Persist a telemetry record so the dashboard /agents stats command
    # and the subagent-pane can show real costs/durations across runs.
    _write_telemetry({
        "ts": time.time(),
        "subagent_type": subagent_type,
        "description": description,
        "isolation": isolation,
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
    )


def is_writer_preset(subagent_type: str) -> bool:
    """True if a subagent type can WRITE (mutate the workspace) — i.e. its
    preset is NOT the read-only 'plan' permission profile. Used to auto-isolate
    parallel writers into separate worktrees so they can't clobber each other."""
    p = SUBAGENT_PRESETS.get((subagent_type or "").strip())
    return bool(p) and (p.mode or "").strip().lower() != "plan"


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
        return {"sa_id": sa_id, "status": "unknown",
                "error": "no running or finished subagent with this id"}
    final = ""
    for m in reversed(sess.get("messages") or []):
        if m.get("role") == "assistant" and m.get("content"):
            final = m.get("content", "")
            break
    return {"sa_id": sa_id, "status": "finished",
            "subagent_type": sess.get("subagent_type", ""),
            "description": sess.get("description", ""),
            "final_text": final,
            "error": sess.get("error", "")}


__all__ = [
    "SubagentPreset",
    "SUBAGENT_PRESETS",
    "list_subagents",
    "subagent_type_names",
    "reload_subagent_presets",
    "SubagentResult",
    "run_subagent",
    "load_subagent_session",
    "list_finished",
    "get_subagent_result",
]
