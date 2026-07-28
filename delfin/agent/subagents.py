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
# Structured returns — a tiny dependency-free JSON-Schema subset validator.
#
# Supported subset (deliberately small; anything else in a schema is
# IGNORED, never enforced):
#   - "type": object / array / string / number / integer / boolean
#   - "required" (on objects)
#   - "properties" (on objects; validated only for keys that are present)
#   - "enum"
#   - "items" (single-schema form, applied to every array element)
# ---------------------------------------------------------------------------

_SCHEMA_TYPES: dict[str, Any] = {
    "object": dict,
    "array": list,
    "string": str,
    "boolean": bool,
    "integer": int,
    "number": (int, float),
}


def validate_json_schema(value, schema: dict, path: str = "$") -> list[str]:
    """Validate ``value`` against the supported JSON-Schema subset.

    Returns a list of human-readable error strings; empty list = valid.
    Unknown/unsupported schema keywords are ignored (subset documented
    above). Never raises on malformed schemas — they just validate less.
    """
    errors: list[str] = []
    if not isinstance(schema, dict):
        return errors
    t = schema.get("type")
    if isinstance(t, str) and t:
        expected = _SCHEMA_TYPES.get(t)
        if expected is not None:
            ok = isinstance(value, expected)
            # bool is an int subclass in Python — a boolean must not
            # satisfy integer/number.
            if t in ("integer", "number") and isinstance(value, bool):
                ok = False
            # ...and an integer/number must not satisfy boolean.
            if t == "boolean" and not isinstance(value, bool):
                ok = False
            if not ok:
                errors.append(
                    f"{path}: expected {t}, got {type(value).__name__}")
                return errors  # wrong type — descending would only cascade
    if "enum" in schema and isinstance(schema.get("enum"), list):
        if value not in schema["enum"]:
            errors.append(
                f"{path}: value {value!r} not in enum {schema['enum']!r}")
    if isinstance(value, dict):
        req = schema.get("required")
        if isinstance(req, list):
            for k in req:
                if k not in value:
                    errors.append(f"{path}: missing required property {k!r}")
        props = schema.get("properties")
        if isinstance(props, dict):
            for k, sub in props.items():
                if k in value and isinstance(sub, dict):
                    errors.extend(
                        validate_json_schema(value[k], sub, f"{path}.{k}"))
    if isinstance(value, list):
        items = schema.get("items")
        if isinstance(items, dict):
            for i, item in enumerate(value):
                errors.extend(
                    validate_json_schema(item, items, f"{path}[{i}]"))
    return errors


def extract_json_object(text: str) -> dict | None:
    """First JSON *object* embedded in ``text``, or None.

    Tolerates markdown fences and surrounding prose: scans for ``{``
    candidates and raw-decodes from each until one parses to a dict.
    """
    if not isinstance(text, str) or not text.strip():
        return None
    dec = json.JSONDecoder()
    i = text.find("{")
    while i != -1:
        try:
            obj, _end = dec.raw_decode(text, i)
        except json.JSONDecodeError:
            i = text.find("{", i + 1)
            continue
        if isinstance(obj, dict):
            return obj
        i = text.find("{", i + 1)
    return None


def _schema_instruction(schema: dict) -> str:
    """System-prompt clause enforcing a schema-shaped final message."""
    try:
        compact = json.dumps(schema, separators=(",", ":"),
                             ensure_ascii=False)
    except (TypeError, ValueError):
        compact = "{}"
    return (
        "\n\nStructured output contract: your FINAL message must be exactly "
        "one JSON object that validates against the JSON Schema below. No "
        "prose before or after it (a ```json fence around it is tolerated).\n"
        "Schema: " + compact
    )


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


def _running_update(sa_id: str, entry: dict | None) -> None:
    """Maintain the live per-subagent status file (entry=None removes).

    File-based so the dashboard can render a live drill-down
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
    model: str = ""          # model the run actually used ("" = unknown)
    model_tier: str = ""     # "cheap" when routed to the cheap tier, else "parent"
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
    # Bump nesting depth so the child's own subagent calls are refused once the
    # depth cap is hit (anti-recursion). Falls back gracefully if the perms
    # type predates the field.
    depth = int(getattr(parent_perms, "subagent_depth", 0)) + 1
    try:
        extra = {"subagent_depth": depth}
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
        structured_output=structured_output,
        schema_error=schema_error,
    )


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

    if len(jobs) == 1:
        return [_one(jobs[0])]
    import concurrent.futures as _cf
    with _cf.ThreadPoolExecutor(
        max_workers=min(len(jobs), _ORCH_MAX_WORKERS),
        thread_name_prefix="orchestration",
    ) as pool:
        futures = [pool.submit(_one, j) for j in jobs]
        return [f.result() for f in futures]


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
    "run_orchestration",
    "validate_json_schema",
    "extract_json_object",
    "load_subagent_session",
    "list_finished",
    "get_subagent_result",
]
