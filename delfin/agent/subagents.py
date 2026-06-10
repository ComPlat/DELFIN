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


_MAX_TOOL_CALLS = 30
_MAX_WALL_S = 60.0
_MAX_OUTPUT_TOKENS = 8000

_TELEMETRY_PATH = Path.home() / ".delfin" / "subagent_telemetry.jsonl"
_TELEMETRY_MAX_LINES = 5000


_RUNNING_PATH = Path.home() / ".delfin" / "subagent_running.json"


def _running_update(sa_id: str, entry: dict | None) -> None:
    """Maintain the live registry of running subagents (entry=None removes).

    File-based so the dashboard can render a Claude-Code-style live panel
    (name · task · status · elapsed) without sharing memory with the
    worker thread."""
    try:
        try:
            data = json.loads(_RUNNING_PATH.read_text(encoding="utf-8"))
        except Exception:
            data = {}
        if entry is None:
            data.pop(sa_id, None)
        else:
            data[sa_id] = entry
        _RUNNING_PATH.parent.mkdir(parents=True, exist_ok=True)
        _RUNNING_PATH.write_text(json.dumps(data), encoding="utf-8")
    except Exception:
        pass


def read_running() -> dict:
    """Live registry: {id: {type, description, started_at}}."""
    try:
        return json.loads(_RUNNING_PATH.read_text(encoding="utf-8"))
    except Exception:
        return {}


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

    def to_payload(self) -> dict:
        return {
            "subagent_type": self.subagent_type,
            "description": self.description,
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
    max_tool_calls: int = _MAX_TOOL_CALLS,
    max_wall_s: float = _MAX_WALL_S,
    max_output_tokens: int = _MAX_OUTPUT_TOKENS,
    isolation: str = "",
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

    Returns a SubagentResult; never raises.
    """
    preset = SUBAGENT_PRESETS.get(subagent_type)
    if preset is None:
        return SubagentResult(
            subagent_type=subagent_type,
            description=description,
            final_text="",
            error=f"unknown subagent_type: {subagent_type!r}. Pick one of {list(SUBAGENT_PRESETS)}.",
        )

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

    # Use the parent client's underlying OpenAI client + model, but
    # swap permissions for the duration of this call.
    saved_perms = getattr(parent_client, "_permissions", None)
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
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": prompt},
    ]

    final_text_parts: list[str] = []
    tool_calls_seen: list[dict] = []
    in_tokens = out_tokens = 0
    t0 = time.monotonic()
    truncated = False
    error = ""
    # Live-panel registry entry (removed in the finally below).
    import uuid as _uuid
    _sa_id = _uuid.uuid4().hex[:8]
    _running_update(_sa_id, {
        "type": subagent_type,
        "description": (description or "")[:120],
        "started_at": time.time(),
    })

    try:
        for event in parent_client.stream_message(
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
            elif event.type == "tool_use":
                tool_calls_seen.append({
                    "name": event.tool_name,
                    "input": event.tool_input,
                })
            elif event.type == "message_delta":
                in_tokens = max(in_tokens, event.input_tokens)
                out_tokens = max(out_tokens, event.output_tokens)
    except Exception as exc:
        error = f"sub-agent stream raised: {exc}"
    finally:
        _running_update(_sa_id, None)
        # Restore parent permissions.
        if hasattr(parent_client, "set_permissions"):
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
        except Exception as exc:
            worktree_summary = {"error": f"worktree teardown failed: {exc}"}

    if isolation_warning and not error:
        # Surface as a soft warning when nothing else went wrong.
        worktree_summary["warning"] = isolation_warning

    final_text = "".join(final_text_parts).strip()
    if not final_text and not error:
        error = "sub-agent returned no text"
    elapsed_s = time.monotonic() - t0
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
    )


__all__ = [
    "SubagentPreset",
    "SUBAGENT_PRESETS",
    "list_subagents",
    "subagent_type_names",
    "reload_subagent_presets",
    "SubagentResult",
    "run_subagent",
]
