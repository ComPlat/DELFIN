"""Settings-driven hooks (.delfin-native schema).

Loads hook definitions from JSON settings files and fires shell
commands at well-defined event points in the agent loop. The schema
mirrors the canonical CLI so existing hook configs port over with minor
tweaks::

    {
      "hooks": {
        "PreToolUse": [
          {
            "matcher": "edit_file|write_file",
            "hooks": [
              {"type": "command", "command": "ruff check ${file}",
               "timeout": 5000}
            ]
          }
        ],
        "PostToolUse": [...],
        "UserPromptSubmit": [...],
        "Stop": [...]
      }
    }

Settings are read from (later wins on overlap):

  1. ``~/.delfin/settings.json``                       — user-global
  2. ``<workspace>/.delfin/settings.json``             — project
  3. ``<workspace>/.delfin/settings.local.json``       — local override

A hook command may exit non-zero to *block* the upcoming tool call
(PreToolUse) or simply log to stderr (other events). It may also
emit a JSON object on stdout::

    {"decision": "block", "reason": "tests failed"}

to deliver a structured block reason back to the agent.

Environment variables passed to every hook:

    CLAUDE_HOOK_EVENT      e.g. "PreToolUse"
    CLAUDE_TOOL_NAME       tool name on tool-use events
    CLAUDE_TOOL_INPUT      JSON-serialised arguments
    CLAUDE_WORKSPACE       resolved workspace root
    CLAUDE_USER_PROMPT     user message on UserPromptSubmit

Inside the ``command`` template, ``${file}`` and ``${cwd}`` are
expanded from the tool arguments so simple linter hooks don't need
to parse JSON.

Failures are tolerated: a misconfigured hook never crashes the
agent, only logs to the audit log.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import audit_log as _audit


_DEFAULT_TIMEOUT_S = 30
_VALID_EVENTS = (
    "PreToolUse",
    "PostToolUse",
    "UserPromptSubmit",
    "Stop",
    "Notification",
)


@dataclass
class HookCommand:
    matcher: str = ""           # regex against tool name
    command: str = ""           # shell command to run
    timeout_s: float = _DEFAULT_TIMEOUT_S
    type: str = "command"

    @property
    def matcher_re(self) -> re.Pattern[str] | None:
        if not self.matcher:
            return None
        try:
            return re.compile(self.matcher)
        except re.error:
            return None


@dataclass
class HookResult:
    event: str
    matched: bool
    exit_code: int = 0
    stdout: str = ""
    stderr: str = ""
    decision: str | None = None     # "block" / "approve" / None
    reason: str = ""
    duration_s: float = 0.0
    command: str = ""

    @property
    def blocks(self) -> bool:
        return self.decision == "block" or self.exit_code not in (0, None)


@dataclass
class HooksConfig:
    by_event: dict[str, list[HookCommand]] = field(default_factory=dict)

    def for_event(self, event: str) -> list[HookCommand]:
        return self.by_event.get(event, [])

    def is_empty(self) -> bool:
        return not any(self.by_event.values())


def _user_settings_path() -> Path:
    return Path.home() / ".delfin" / "settings.json"


def _project_settings_paths(workspace: Path) -> list[Path]:
    base = Path(workspace) / ".delfin"
    return [base / "settings.json", base / "settings.local.json"]


def _read_json_safe(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {}


def _merge_hooks(into: HooksConfig, raw: dict[str, Any]) -> None:
    hooks_obj = raw.get("hooks") if isinstance(raw, dict) else None
    if not isinstance(hooks_obj, dict):
        return
    for event, entries in hooks_obj.items():
        if event not in _VALID_EVENTS:
            continue
        if not isinstance(entries, list):
            continue
        bucket = into.by_event.setdefault(event, [])
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            matcher = str(entry.get("matcher", ""))
            cmds = entry.get("hooks", [])
            if not isinstance(cmds, list):
                continue
            for c in cmds:
                if not isinstance(c, dict):
                    continue
                bucket.append(
                    HookCommand(
                        matcher=matcher,
                        command=str(c.get("command", "")),
                        timeout_s=float(c.get("timeout", _DEFAULT_TIMEOUT_S)) / 1000.0
                        if c.get("timeout", 0) and c.get("timeout", 0) > 100
                        else float(c.get("timeout", _DEFAULT_TIMEOUT_S)),
                        type=str(c.get("type", "command")),
                    )
                )


def load_hooks(
    workspace: Path | str | None = None,
    *,
    extra_paths: list[Path] | None = None,
) -> HooksConfig:
    """Read all settings files and return a merged HooksConfig."""
    cfg = HooksConfig()
    paths: list[Path] = [_user_settings_path()]
    if workspace is not None:
        paths.extend(_project_settings_paths(Path(workspace)))
    for p in extra_paths or ():
        paths.append(Path(p))
    for p in paths:
        _merge_hooks(cfg, _read_json_safe(p))
    return cfg


def _expand(template: str, arguments: dict[str, Any]) -> str:
    """Replace ${var} placeholders from arguments. Unknown vars left as-is."""
    def repl(m: re.Match[str]) -> str:
        key = m.group(1)
        val = arguments.get(key, m.group(0))
        if isinstance(val, (dict, list)):
            return json.dumps(val)
        return str(val)
    return re.sub(r"\$\{(\w+)\}", repl, template)


def _build_env(
    event: str,
    *,
    tool_name: str = "",
    arguments: dict[str, Any] | None = None,
    user_prompt: str = "",
    workspace: Path | str | None = None,
) -> dict[str, str]:
    env = dict(os.environ)
    env["CLAUDE_HOOK_EVENT"] = event
    if tool_name:
        env["CLAUDE_TOOL_NAME"] = tool_name
    if arguments is not None:
        try:
            env["CLAUDE_TOOL_INPUT"] = json.dumps(arguments, default=str)[:8000]
        except (TypeError, ValueError):
            env["CLAUDE_TOOL_INPUT"] = ""
    if user_prompt:
        env["CLAUDE_USER_PROMPT"] = user_prompt[:8000]
    if workspace:
        env["CLAUDE_WORKSPACE"] = str(workspace)
    return env


def _run_command(
    cmd: HookCommand,
    env: dict[str, str],
    cwd: Path | None,
    arguments: dict[str, Any] | None,
) -> tuple[int, str, str, float]:
    expanded = _expand(cmd.command, arguments or {})
    t0 = time.monotonic()
    try:
        proc = subprocess.run(
            expanded, shell=True, env=env,
            cwd=str(cwd) if cwd else None,
            capture_output=True, text=True,
            timeout=max(0.1, cmd.timeout_s),
        )
        return proc.returncode, proc.stdout or "", proc.stderr or "", time.monotonic() - t0
    except subprocess.TimeoutExpired as exc:
        return 124, exc.stdout or "", (exc.stderr or "") + "\n[hook timeout]", time.monotonic() - t0
    except Exception as exc:
        return 1, "", f"[hook error] {exc}", time.monotonic() - t0


def _parse_decision(stdout: str) -> tuple[str | None, str]:
    """If stdout is a JSON dict with a decision, parse it."""
    s = (stdout or "").strip()
    if not s.startswith("{"):
        return None, ""
    try:
        obj = json.loads(s)
    except (json.JSONDecodeError, ValueError):
        return None, ""
    if not isinstance(obj, dict):
        return None, ""
    decision = obj.get("decision")
    reason = str(obj.get("reason", ""))
    if decision in ("block", "approve"):
        return decision, reason
    return None, reason


def run_hooks(
    event: str,
    config: HooksConfig,
    *,
    tool_name: str = "",
    arguments: dict[str, Any] | None = None,
    user_prompt: str = "",
    workspace: Path | str | None = None,
) -> list[HookResult]:
    """Run all matching hooks for ``event``. Returns one HookResult per fire."""
    if event not in _VALID_EVENTS:
        return []
    hooks = config.for_event(event)
    if not hooks:
        return []
    env = _build_env(
        event, tool_name=tool_name, arguments=arguments,
        user_prompt=user_prompt, workspace=workspace,
    )
    results: list[HookResult] = []
    cwd = Path(workspace) if workspace else None
    for hk in hooks:
        if hk.type != "command" or not hk.command.strip():
            continue
        pat = hk.matcher_re
        # Empty matcher means "match all"; tool-less events ignore matcher.
        if event in ("PreToolUse", "PostToolUse"):
            if pat and tool_name and not pat.search(tool_name):
                continue
        rc, out, err, dur = _run_command(hk, env, cwd, arguments)
        decision, reason = _parse_decision(out)
        results.append(HookResult(
            event=event, matched=True,
            exit_code=rc, stdout=out, stderr=err,
            decision=decision,
            reason=reason or err.strip()[:240],
            duration_s=dur,
            command=hk.command,
        ))
        # Fire-and-forget audit
        try:
            _audit.append(_audit.make_record(
                tool="hook",
                decision="block" if (decision == "block" or rc != 0) else "ok",
                mode="hook",
                command=f"[{event}] {hk.command[:200]}",
                extra={
                    "matcher": hk.matcher,
                    "exit_code": rc,
                    "duration_s": round(dur, 3),
                },
            ))
        except Exception:
            pass
    return results


def first_block(results: list[HookResult]) -> HookResult | None:
    for r in results:
        if r.blocks:
            return r
    return None


__all__ = [
    "HookCommand",
    "HookResult",
    "HooksConfig",
    "load_hooks",
    "run_hooks",
    "first_block",
]
