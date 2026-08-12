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

The user's own file is always read: it is not a workspace and never
needed trusting. The two workspace files are read only when the USER
has trusted that directory for hooks — see ``workspace_trust``. A
repository the user merely checked out ships no commands this module
will run, and what it withheld is reported rather than dropped.

A hook command may exit non-zero to *block* the upcoming tool call
(PreToolUse) or simply log to stderr (other events). It may also
emit a JSON object on stdout::

    {"decision": "block", "reason": "tests failed"}

to deliver a structured block reason back to the agent.

Environment variables passed to every hook:

    DELFIN_HOOK_EVENT      e.g. "PreToolUse"
    DELFIN_TOOL_NAME       tool name on tool-use events
    DELFIN_TOOL_INPUT      JSON-serialised arguments
    DELFIN_WORKSPACE       resolved workspace root
    DELFIN_USER_PROMPT     user message on UserPromptSubmit

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
# Only events that something actually FIRES. "Notification" sat here and
# was raised nowhere: a user could configure a hook for it, the settings
# would validate, and nothing would ever run. A declared event that never
# arrives is a promise the configuration makes and the code does not keep,
# and it is worse than an absent one because it looks configured.
#
# The notification need itself is already served by push_notification and
# remote_trigger, which are tools rather than hooks; adding a third path
# to the same place would be the wrong way to fix a dead enum member.
_VALID_EVENTS = (
    "PreToolUse",
    "PostToolUse",
    "UserPromptSubmit",
    "Stop",
)


@dataclass
class HookCommand:
    matcher: str = ""           # regex against tool name
    command: str = ""           # shell command to run
    timeout_s: float = _DEFAULT_TIMEOUT_S
    type: str = "command"
    # Absolute path of the settings file this came from. A hook the user
    # wrote and a hook a repository shipped ran identically and printed
    # identically, so the listing could not tell them apart; carrying the
    # source makes /hooks, the audit record and the changes report all
    # able to say WHOSE hook ran.
    source: str = ""

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
    source: str = ""

    @property
    def blocks(self) -> bool:
        return self.decision == "block" or self.exit_code not in (0, None)


@dataclass
class HooksConfig:
    by_event: dict[str, list[HookCommand]] = field(default_factory=dict)
    # Everything the load could not do and did not crash over: a
    # malformed entry that was dropped, a workspace whose hooks were
    # withheld for lack of trust. Every call site swallowed the
    # exceptions these replace, so the failures were invisible.
    warnings: list[str] = field(default_factory=list)

    def for_event(self, event: str) -> list[HookCommand]:
        return self.by_event.get(event, [])

    def is_empty(self) -> bool:
        return not any(self.by_event.values())


def _user_settings_path() -> Path:
    return Path.home() / ".delfin" / "settings.json"


def _project_settings_paths(workspace: Path) -> list[Path]:
    """Every settings file this loader knows about inside *workspace*.

    Derived from the trust registry rather than spelled out here, so the
    set of files that RUN commands and the set of files that need a
    trust decision (and a protected glob) cannot drift apart. They had:
    ``settings.local.json`` was read here and protected nowhere.
    """
    from . import workspace_trust as _trust
    root = Path(workspace)
    return [root / rel
            for rel in _trust.get_kind(_trust.KIND_HOOKS).relative_paths]


def _read_json_safe(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {}


def _parse_timeout(raw: Any, warn: list[str], where: str) -> float | None:
    """Seconds for one hook entry, or None when the entry must be dropped.

    A ``timeout`` that is not a number raised TypeError out of the whole
    merge, and every call site swallows exceptions around ``load_hooks``
    -- so one typo in any settings layer silently disabled ALL hooks,
    including the user's own blocking ones. A configuration error must
    cost the entry that has it and nothing else, and must be said out
    loud: hooks that stop running look exactly like hooks that found
    nothing to complain about.

    Values above 100 are read as milliseconds (the canonical schema),
    at or below as seconds.
    """
    if raw is None or raw is False or raw == 0 or raw == "":
        return float(_DEFAULT_TIMEOUT_S)
    if isinstance(raw, bool):
        value = float(int(raw))
    else:
        try:
            value = float(raw)
        except (TypeError, ValueError):
            warn.append(
                f"{where}: timeout must be a number, got {raw!r} — this hook "
                f"entry was dropped; the others still load."
            )
            return None
    if value <= 0:
        warn.append(
            f"{where}: timeout must be positive, got {raw!r} — this hook "
            f"entry was dropped; the others still load."
        )
        return None
    return value / 1000.0 if value > 100 else value


def _merge_hooks(into: HooksConfig, raw: dict[str, Any],
                 source: Path | str = "") -> None:
    hooks_obj = raw.get("hooks") if isinstance(raw, dict) else None
    if not isinstance(hooks_obj, dict):
        return
    src = str(source)
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
                command = str(c.get("command", ""))
                timeout_s = _parse_timeout(
                    c.get("timeout"), into.warnings,
                    f"{src or '<settings>'} [{event}] {command[:60]}")
                if timeout_s is None:
                    continue
                bucket.append(
                    HookCommand(
                        matcher=matcher,
                        command=command,
                        timeout_s=timeout_s,
                        type=str(c.get("type", "command")),
                        source=src,
                    )
                )


def _project_hooks_allowed(workspace: Path | str | None) -> bool:
    """Whether hooks MAY be considered at all for *workspace*.

    A hook file is executable configuration: it runs a shell command
    before and after every tool call, outside the permission gate and
    outside filesystem isolation. That is right for a project the user
    chose to work on, and wrong for a folder that receives files from
    other people -- which is exactly what a locked scope means.

    The check lives HERE rather than at the call sites because it was
    written at one of four call sites and forgotten at the other three:
    the engine's own turn path and both dashboard paths loaded workspace
    hooks unguarded, so a locked office folder's settings file executed on
    every message. A guard a caller has to remember is a guard that gets
    forgotten.

    This answers only "is this folder categorically ineligible". The
    question that decides an ordinary directory -- "did the USER say they
    trust it?" -- is ``workspace_trust.gate`` and is asked below, for the
    same reason and in the same place. Being eligible was never consent:
    every directory that was not a registered office folder passed here,
    which includes one cloned a second ago.
    """
    if workspace is None:
        return False
    try:
        from .memory_store import is_office_workspace
        if is_office_workspace(workspace):
            return False
    except Exception:
        pass
    return True


def load_hooks(
    workspace: Path | str | None = None,
    *,
    extra_paths: list[Path] | None = None,
) -> HooksConfig:
    """Read all settings files and return a merged HooksConfig.

    The user's own ``~/.delfin/settings.json`` is always read: it is not
    a workspace and never needed trusting. A workspace's files are read
    only when the user has explicitly trusted that directory for hooks,
    with the exact content they trusted -- ``workspace_trust.gate``.
    What was withheld lands in ``cfg.warnings`` rather than being
    dropped, because hooks that stop running look exactly like hooks
    that had nothing to say.

    ``extra_paths`` is for a caller naming a file itself (the CLI, a
    test). It is deliberately NOT gated: naming a path is the user's own
    decision, the same one the trust store records. Nothing derives it
    from the workspace.
    """
    from . import workspace_trust as _trust

    cfg = HooksConfig()
    paths: list[Path] = [_user_settings_path()]
    decision = None
    if workspace is not None and _project_hooks_allowed(workspace):
        decision = _trust.gate(_trust.KIND_HOOKS, workspace)
        paths.extend(decision.paths)
    for p in extra_paths or ():
        paths.append(Path(p))
    for p in paths:
        _merge_hooks(cfg, _read_json_safe(p), source=p)
    if decision is not None and decision.notice:
        cfg.warnings.append(decision.notice)
    return cfg


def _expand(template: str, arguments: dict[str, Any]) -> str:
    """Replace ${var} placeholders from arguments. Unknown vars left as-is.

    Every substituted value is shell-quoted. The expanded string is run
    with ``shell=True``, so an unquoted value is command text: the
    documented example hook is ``ruff check ${file}``, and a file called
    ``x.py; curl … | sh`` — named by the model, or simply present in an
    untrusted repository — would then execute. Nothing the hook author
    wrote is wrong; the substitution was the sink. ``shlex`` has been
    imported for this since the module was written and was never called.

    An unknown placeholder is left as literal text and NOT quoted: it is
    not a value, and quoting it would change what the author sees.
    """
    def repl(m: re.Match[str]) -> str:
        key = m.group(1)
        if key not in arguments:
            return m.group(0)
        val = arguments[key]
        if isinstance(val, (dict, list)):
            val = json.dumps(val)
        return shlex.quote(str(val))
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
    env["DELFIN_HOOK_EVENT"] = event
    if tool_name:
        env["DELFIN_TOOL_NAME"] = tool_name
    if arguments is not None:
        try:
            env["DELFIN_TOOL_INPUT"] = json.dumps(arguments, default=str)[:8000]
        except (TypeError, ValueError):
            env["DELFIN_TOOL_INPUT"] = ""
    if user_prompt:
        env["DELFIN_USER_PROMPT"] = user_prompt[:8000]
    if workspace:
        env["DELFIN_WORKSPACE"] = str(workspace)
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
            source=hk.source,
        ))
        # Fire-and-forget audit. The SOURCE is part of the record: a hook
        # the user wrote and a hook a repository shipped executed
        # identically and logged identically, so "what ran a shell
        # command during this session" could not be answered.
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
                    "hook_event": event,
                    # Not ``path``: the changes report drops records whose
                    # absolute path lies outside the workspace, and the
                    # user's own settings file always does.
                    "source": hk.source,
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
