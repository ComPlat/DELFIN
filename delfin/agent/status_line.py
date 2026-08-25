"""Configurable status line (.delfin-native).

Reads ``statusLine`` from settings.json and renders the agent's
current state as a single line. The setting can be either:

  - a shell command (``"command": "..."``) that is run with the
    status payload available as JSON on stdin and the rendered
    line read from stdout;
  - a built-in template (``"template": "{branch} | {tokens} | {mode}"``)
    that we expand without spawning a subprocess.

Default template if no config exists::

    {tokens} tokens | mode={mode} | branch={branch}

Lookup order, later wins:

  1. ``~/.delfin/settings.json``
  2. ``<workspace>/.delfin/settings.json``
  3. ``<workspace>/.delfin/settings.local.json``

A misconfigured status line never crashes — render returns "" on
any failure.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


_DEFAULT_TEMPLATE = "{tokens} tokens | mode={mode} | branch={branch}"


def _read_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def _gather_status_lines(workspace: Path | None) -> list[dict]:
    """Return statusLine specs in user → project → local order."""
    paths = [Path.home() / ".delfin" / "settings.json"]
    if workspace is not None:
        paths.extend([
            workspace / ".delfin" / "settings.json",
            workspace / ".delfin" / "settings.local.json",
        ])
    out: list[dict] = []
    for idx, p in enumerate(paths):
        sl = _read_json(p).get("statusLine")
        if isinstance(sl, dict):
            spec = dict(sl)
            # A template is data; a command is code. The first path is the
            # user's own settings file, the rest are inside the workspace
            # -- and the winning spec's `command` is run with shell=True,
            # cwd set to that workspace, on every status refresh, before
            # the agent has taken a single action. Granting the agent a
            # colleague's directory, or opening a repository that ships
            # `.delfin/settings.local.json`, was enough to execute a
            # command of that folder's choosing, with its stdout becoming
            # the status line and its stderr discarded. No allow-list, no
            # confirmation, no security event, no audit record.
            #
            # Same reasoning as the hooks file, and the same rule: the
            # workspace may describe how the line LOOKS, and may not
            # decide what RUNS.
            if idx > 0:
                spec.pop("command", None)
            out.append(spec)
        elif isinstance(sl, str):
            out.append({"template": sl})
    return out


@dataclass
class StatusContext:
    workspace: Path | None = None
    mode: str = "default"
    branch: str = ""
    model: str = ""
    tokens: int = 0
    cost_usd: float = 0.0
    extras: dict[str, Any] = field(default_factory=dict)

    def to_payload(self) -> dict:
        return {
            "workspace": str(self.workspace) if self.workspace else "",
            "mode": self.mode,
            "branch": self.branch,
            "model": self.model,
            "tokens": self.tokens,
            "cost_usd": self.cost_usd,
            **self.extras,
        }


def _git_branch(workspace: Path | None) -> str:
    if workspace is None:
        return ""
    try:
        out = subprocess.run(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"],
            cwd=str(workspace), capture_output=True, text=True, timeout=2,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except (FileNotFoundError, subprocess.SubprocessError):
        pass
    return ""


def _expand_template(tpl: str, ctx: StatusContext) -> str:
    payload = ctx.to_payload()
    try:
        return tpl.format(**payload)
    except (KeyError, IndexError, ValueError):
        return tpl


def has_custom_status_line(workspace: Path | None) -> bool:
    """True when the user actually configured one.

    Callers that describe the line as the user's own need to be able to
    tell "configured" from "the built-in default fired": the terminal
    printed the default after every turn under a docstring saying it
    printed nothing unless configured, and the default repeats two fields
    the banner and the live turn line already carry.
    """
    return bool(_gather_status_lines(workspace))


def _render_default(ctx: StatusContext) -> str:
    """The built-in line, with unknown fields left out entirely.

    Formatting the default template against an empty branch produced
    ``0 tokens | mode=plan | branch=`` outside a git repository — a
    labelled field with nothing after it, which reads as a lookup that
    failed rather than as a directory that is not a repository. A user's
    own template still gets the empty string, because that is the truth
    and their template decides how to show it.
    """
    parts = [f"{ctx.tokens} tokens", f"mode={ctx.mode}"]
    if ctx.branch:
        parts.append(f"branch={ctx.branch}")
    return " | ".join(parts)


def render_status_line(ctx: StatusContext) -> str:
    """Render the active statusLine for the given context.

    Falls back to the default template if no config is found and
    fills in ``branch`` from git when missing.
    """
    if not ctx.branch:
        ctx.branch = _git_branch(ctx.workspace)
    specs = _gather_status_lines(ctx.workspace)
    if specs:
        # later (project / local) wins
        spec = specs[-1]
    else:
        return _render_default(ctx)[:240]
    if "command" in spec and isinstance(spec["command"], str):
        cmd = spec["command"]
        try:
            payload = json.dumps(ctx.to_payload(), default=str)
            proc = subprocess.run(
                cmd, shell=True, input=payload,
                capture_output=True, text=True,
                timeout=float(spec.get("timeout_s", 3.0)),
                cwd=str(ctx.workspace) if ctx.workspace else None,
            )
            return (proc.stdout or "").strip()[:240]
        except (subprocess.SubprocessError, OSError):
            # OSError too: cwd is passed unchecked, so a workspace that
            # has since been deleted or renamed raises FileNotFoundError
            # -- which is not a SubprocessError and escaped the render.
            # The caller wraps the whole status refresh in a bare except,
            # so the status line simply vanished with no explanation.
            return ""
    tpl = str(spec.get("template") or _DEFAULT_TEMPLATE)
    return _expand_template(tpl, ctx)[:240]


__all__ = ["StatusContext", "render_status_line", "has_custom_status_line"]
