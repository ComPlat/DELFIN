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
    """A single event from a streaming Claude response."""

    __slots__ = ("type", "text", "input_tokens", "output_tokens", "stop_reason",
                 "cost_usd", "tool_name", "tool_input", "tool_output")

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
                 extra_dirs: list[str] | None = None):
        self.model = model or self.DEFAULT_MODEL
        self.permission_mode = permission_mode
        self.cwd = cwd or None
        self.mcp_config = mcp_config  # path to MCP config JSON or empty
        self.allowed_tools = allowed_tools  # restrict CLI to these tools only
        self.extra_dirs = extra_dirs  # extra writable directories (--add-dir)
        self.claude_path = claude_path or shutil.which("claude") or "claude"
        if not shutil.which(self.claude_path):
            raise FileNotFoundError(
                f"Claude Code CLI not found at '{self.claude_path}'. "
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
                        "Claude Code CLI is not logged in. "
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
        resolved_key = api_key or os.environ.get("ANTHROPIC_API_KEY", "")
        if not resolved_key:
            raise ValueError(
                "No Anthropic API key found. Set the ANTHROPIC_API_KEY "
                "environment variable before launching the dashboard."
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
# KIT-Toolbox coding-agent permissions (Claude-Code-style safety layer)
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
    r"git\s+clean\s+-[a-zA-Z]*f[a-zA-Z]*d",
    r"git\s+update-ref\s+-d\b",
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
    r"^\s*find\b(?![^|;&]*-delete)(?![^|;&]*-exec[^|;&]*rm)",
    r"^\s*grep\b", r"^\s*rg\b", r"^\s*ag\b",
    r"^\s*sed\s+-n\b",                                       # read-only sed
    r"^\s*awk\s+",                                           # awk has no destructive default
    r"^\s*jq\b", r"^\s*yq\b",
    r"^\s*git\s+(?:status|diff|log|show|branch(?!\s+-D)|remote|config\s+--get|"
    r"rev-parse|describe|ls-files|ls-tree|blame|stash\s+list|tag\s*$|"
    r"shortlog|reflog|fetch|pull(?!\s+--rebase\s+--force)|switch|checkout|add|"
    r"restore(?!\s+--source)|commit\s+-m|commit\s+--message)\b",
    r"^\s*tar\s+-?t",                                        # tar list-only
    r"^\s*unzip\s+-l\b",
    # -- (b) coding workflow ---------------------------------------------
    r"^\s*python3?\s+-c\s+",
    r"^\s*python3?\s+--version\b",
    r"^\s*python3?\s+-m\s+(?:py_compile|pytest|unittest|doctest|timeit|venv|"
    r"pip\s+show|pip\s+list|pip\s+freeze|"
    r"compileall|json\.tool|http\.server|delfin(?:\.\w+)*)\b",
    r"^\s*python3?\s+\S+\.py\b",                             # run a script in repo
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
    r"^\s*xtb\s+\S+\.xyz\b",                                 # DELFIN-spezifisch: read-only xtb-Aufruf
    r"^\s*delfin(?:-\w+)?\b",                                # delfin-CLI-Wrapper
)

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
)


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
    bash_deny_patterns: tuple[str, ...] = _DEFAULT_BASH_DENY_PATTERNS
    bash_auto_allow_patterns: tuple[str, ...] = _DEFAULT_BASH_AUTO_ALLOW
    path_deny_globs: tuple[str, ...] = _DEFAULT_PATH_DENY_GLOBS
    path_protected_globs: tuple[str, ...] = _DEFAULT_PATH_PROTECTED_GLOBS
    bash_timeout_s: int = 120
    bash_max_timeout_s: int = 600
    max_output_chars: int = 12_000
    confirm_callback: Optional[Callable[[str, dict, str], bool]] = None
    pre_tool_hook: Optional[Callable[[str, dict], None]] = None
    post_tool_hook: Optional[Callable[[str, dict, str], None]] = None
    read_tracker: dict[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.workspace = Path(self.workspace).expanduser().resolve()
        valid = {"plan", "default", "acceptEdits", "bypassPermissions"}
        if self.mode not in valid:
            raise ValueError(f"mode must be one of {valid}, got {self.mode!r}")

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
        for pat in self.bash_auto_allow_patterns:
            if re.search(pat, cmd, re.IGNORECASE):
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
                "Read a file from the DELFIN repository. Returns the file content. "
                "Use for .py, .json, .md, .yaml, .txt, .xyz, .out files. "
                "For large files, use offset and limit to read a specific range."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "File path relative to repo root (e.g. 'delfin/agent/learned_profiles.json')",
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
                "Create a new file or fully overwrite an existing file in the "
                "workspace. Requires sandbox-relative path. For existing files, "
                "you MUST call read_file first (the tool tracks read mtimes). "
                "Returns a unified diff preview of the change. Use edit_file "
                "for partial changes — write_file replaces the entire file."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Workspace-relative file path.",
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
                "(strip the line-number prefix, keep leading whitespace)."
            ),
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Workspace-relative file path.",
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
                        "description": "Workspace-relative file path.",
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
                        "description": "Workspace-relative cwd (default = workspace root).",
                    },
                },
                "required": ["command", "description"],
            },
        },
    },
]


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
        return result

    def _dispatch(
        self,
        name: str,
        arguments: dict,
        permissions: Optional["KitToolPermissions"],
    ) -> str:
        # Coding-agent tools are only available with explicit permissions.
        if name in ("write_file", "edit_file", "multi_edit", "bash"):
            if permissions is None:
                return json.dumps({"error": (
                    f"Tool '{name}' requires permissions to be configured. "
                    "Pass a KitToolPermissions instance to OpenAIClient."
                )})
            gate_err = self._run_permission_gate(name, arguments, permissions)
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

        # Calc search tools
        if name in ("search_calcs", "get_calc_info", "calc_summary"):
            return self._execute_calc(name, arguments)

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

        # Repo file access tools (workspace-aware when permissions are set)
        if name == "read_file":
            return self._execute_read_file(arguments, permissions)
        elif name == "grep_file":
            return self._execute_grep_file(arguments, permissions)
        elif name == "list_files":
            return self._execute_list_files(arguments, permissions)

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
            return json.dumps(info, indent=2, ensure_ascii=False)

        elif name == "calc_summary":
            return json.dumps(
                self._calc_engine.summary(), indent=2, ensure_ascii=False
            )

        return json.dumps({"error": f"Unknown calc tool: {name}"})

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
        rel_path = arguments.get("path", "")
        if not rel_path:
            return json.dumps({"error": "path is required"})
        root = perms.workspace if perms is not None else self._repo_root()
        full = (root / rel_path) if not Path(rel_path).is_absolute() else Path(rel_path)
        if not full.exists():
            return json.dumps({"error": f"File not found: {rel_path}"})
        if full.is_dir():
            entries = sorted(p.name for p in full.iterdir())[:50]
            return json.dumps({"type": "directory", "entries": entries})
        try:
            lines = full.read_text(encoding="utf-8", errors="replace").splitlines()
        except Exception as exc:
            return json.dumps({"error": str(exc)})
        offset = arguments.get("offset", 0) or 0
        limit = arguments.get("limit", 200) or 200
        selected = lines[offset:offset + limit]
        result = "\n".join(f"{i + offset}  {line}" for i, line in enumerate(selected))
        if len(lines) > offset + limit:
            result += f"\n... ({len(lines)} lines total, showing {offset}-{offset + limit})"
        return result

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
        max_results = arguments.get("max_results", 30) or 30
        try:
            regex = _re.compile(pattern, _re.IGNORECASE)
        except _re.error as exc:
            return json.dumps({"error": f"Invalid regex: {exc}"})
        matches = []
        files = [search_path] if search_path.is_file() else sorted(search_path.rglob("*"))
        for fp in files:
            if not fp.is_file() or fp.suffix in (".pyc", ".so", ".gz", ".pdf", ".png", ".jpg"):
                continue
            if "/__pycache__/" in str(fp) or "/.git/" in str(fp):
                continue
            try:
                for i, line in enumerate(fp.read_text(encoding="utf-8", errors="replace").splitlines()):
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
        import fnmatch
        pattern = arguments.get("pattern", "*")
        root = perms.workspace if perms is not None else self._repo_root()
        matches = []
        for fp in sorted(root.rglob("*")):
            if not fp.is_file():
                continue
            rel = str(fp.relative_to(root))
            if "/__pycache__/" in rel or "/.git/" in rel:
                continue
            if fnmatch.fnmatch(rel, pattern):
                matches.append(rel)
                if len(matches) >= 50:
                    break
        return "\n".join(matches) if matches else "No files matching pattern."

    # -- Coding-agent helpers (sandbox + permission gating) -------------------

    def _workspace_root(self, perms: "KitToolPermissions") -> Path:
        return perms.workspace

    def _resolve_in_workspace(
        self, rel_path: str, perms: "KitToolPermissions"
    ) -> tuple[Optional[Path], Optional[str]]:
        """Resolve ``rel_path`` against the workspace and verify containment.

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
        try:
            resolved.relative_to(ws)
        except ValueError:
            return None, f"path escapes workspace sandbox ({ws}): {rel_path}"
        try:
            rel_str = str(resolved.relative_to(ws)).replace("\\", "/")
        except ValueError:
            rel_str = rel_path.replace("\\", "/")
        if perms.matches_path_deny(rel_str):
            return None, f"path is on the deny-list (.git, secrets, keys, .env, ...): {rel_str}"
        return resolved, None

    def _run_permission_gate(
        self, name: str, args: dict, perms: "KitToolPermissions"
    ) -> Optional[str]:
        """Run the policy + callback gate. Returns error string or None."""
        mode = perms.mode

        if mode == "plan":
            return f"plan mode is read-only — '{name}' rejected"

        if name in ("write_file", "edit_file", "multi_edit"):
            path_arg = args.get("path", "")
            resolved, err = self._resolve_in_workspace(path_arg, perms)
            if err:
                return err
            try:
                rel_str = str(resolved.relative_to(perms.workspace)).replace("\\", "/")
            except Exception:
                rel_str = path_arg.replace("\\", "/")
            is_protected = perms.matches_path_protected(rel_str)
            if is_protected:
                # Force explicit user confirmation regardless of mode.
                if perms.confirm_callback is None:
                    return (
                        f"'{name}' targets the agent's own safety layer "
                        f"('{rel_str}'). This always requires explicit user "
                        "confirmation but no confirm_callback is configured — "
                        "refusing to proceed."
                    )
                preview = (
                    "[SELF-MODIFICATION GUARD]\n"
                    "This file is part of the agent's own safety layer.\n"
                    "Approving this will let the agent rewrite the code that "
                    "controls its own permissions.\n\n"
                    + self._build_change_preview(name, args, resolved)
                )
                try:
                    ok = bool(perms.confirm_callback(name, args, preview))
                except Exception as exc:
                    return f"confirm_callback raised: {exc}"
                return None if ok else f"user denied '{name}' on protected path '{rel_str}'"
            if mode == "acceptEdits" or mode == "bypassPermissions":
                return None
            # mode == "default": ask the callback.
            if perms.confirm_callback is None:
                return (
                    f"'{name}' on '{path_arg}' requires user confirmation, but "
                    "no confirm_callback is configured. Set permissions.mode "
                    "to 'acceptEdits' or pass a confirm_callback."
                )
            preview = self._build_change_preview(name, args, resolved)
            try:
                ok = bool(perms.confirm_callback(name, args, preview))
            except Exception as exc:
                return f"confirm_callback raised: {exc}"
            return None if ok else f"user denied '{name}' on '{path_arg}'"

        if name == "bash":
            cmd = args.get("command", "") or ""
            if not cmd.strip():
                return "command is required"
            denied = perms.matches_bash_deny(cmd)
            if denied:
                return f"command rejected by deny-pattern {denied!r}: refusing to run."
            if mode == "bypassPermissions":
                return None
            if perms.matches_bash_auto_allow(cmd):
                return None
            if perms.confirm_callback is None:
                return (
                    f"bash command requires user confirmation, but no "
                    "confirm_callback is configured. Either configure one, "
                    "use mode='bypassPermissions' for trusted environments, "
                    f"or restrict to auto-allow patterns. Command: {cmd[:200]}"
                )
            preview = f"$ {cmd}\n(description: {args.get('description', '<none>')})"
            try:
                ok = bool(perms.confirm_callback(name, args, preview))
            except Exception as exc:
                return f"confirm_callback raised: {exc}"
            return None if ok else "user denied bash command"

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
        path_arg = arguments.get("path", "")
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

        diff = self._make_diff(old_text, content, str(resolved.relative_to(perms.workspace)))
        action = "created" if not existed else "overwritten"
        return f"File {action}: {resolved.relative_to(perms.workspace)}\n\n{diff}"

    def _execute_edit_file(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        path_arg = arguments.get("path", "")
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
        if count == 0:
            return json.dumps({"error": f"old_string not found in '{path_arg}'"})
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

        diff = self._make_diff(old_text, new_text, str(resolved.relative_to(perms.workspace)))
        replaced = count if replace_all else 1
        return f"Edited {resolved.relative_to(perms.workspace)} ({replaced} replacement(s)):\n\n{diff}"

    def _execute_multi_edit(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        path_arg = arguments.get("path", "")
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
                return json.dumps({"error": (
                    f"edit #{i+1}: old_string not found "
                    f"(after applying earlier edits in this batch)"
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

        diff = self._make_diff(old_text, text, str(resolved.relative_to(perms.workspace)))
        total = sum(per_edit_replacements)
        return (
            f"Multi-edited {resolved.relative_to(perms.workspace)} "
            f"({len(edits)} edit(s), {total} replacement(s) total):\n\n{diff}"
        )

    def _execute_bash(
        self, arguments: dict, perms: "KitToolPermissions"
    ) -> str:
        cmd = arguments.get("command", "") or ""
        description = arguments.get("description", "") or ""
        timeout = int(arguments.get("timeout_s", perms.bash_timeout_s) or perms.bash_timeout_s)
        timeout = max(1, min(timeout, perms.bash_max_timeout_s))
        cwd_arg = arguments.get("cwd", "") or ""

        if cwd_arg:
            cwd_resolved, err = self._resolve_in_workspace(cwd_arg, perms)
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
                ["/bin/bash", "-c", cmd],
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
        if len(out) > cap:
            out = out[:cap] + f"\n... (stdout truncated, {len(proc.stdout) - cap} chars omitted)"
        if len(err) > cap:
            err = err[:cap] + f"\n... (stderr truncated, {len(proc.stderr) - cap} chars omitted)"

        return json.dumps({
            "exit_code": proc.returncode,
            "elapsed_s": round(elapsed, 3),
            "stdout": out,
            "stderr": err,
            "command": cmd[:500],
            "description": description,
            "cwd": str(run_cwd.relative_to(perms.workspace)) or ".",
        }, ensure_ascii=False)


# Singleton — shared across all OpenAIClient instances.
_doc_executor = _DocToolExecutor()


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
                 permissions: Optional["KitToolPermissions"] = None):
        try:
            import openai  # noqa: F401
        except ImportError:
            _auto_install("openai", "openai>=1.0")
            import openai  # noqa: F401
        resolved_key = api_key or os.environ.get(key_env_var, "")
        if not resolved_key:
            raise ValueError(
                f"No API key found. Set the {key_env_var} "
                "environment variable before launching the dashboard."
            )
        import openai

        self.model = model or self.DEFAULT_MODEL
        kwargs: dict[str, Any] = {"api_key": resolved_key}
        if base_url:
            kwargs["base_url"] = base_url
        self.client = openai.OpenAI(**kwargs)
        # KIT-Toolbox coding-agent permissions (None disables write/edit/bash).
        self._permissions: Optional[KitToolPermissions] = permissions

    def set_permissions(self, permissions: Optional["KitToolPermissions"]) -> None:
        """Replace the KIT-Toolbox permissions policy at runtime."""
        self._permissions = permissions

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
        # Detect reasoning models (o3, o4-mini, azure.o3, azure.o4-mini)
        _base = self.model.split(".", 1)[-1] if self.model.startswith(("azure.", "kit.")) else self.model
        is_reasoning = _base.startswith("o")

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

        _CODING_TOOL_NAMES = {"write_file", "edit_file", "multi_edit", "bash"}
        if has_coding:
            advertised_tools = _DOC_TOOLS_OPENAI
        else:
            advertised_tools = [
                t for t in _DOC_TOOLS_OPENAI
                if t.get("function", {}).get("name") not in _CODING_TOOL_NAMES
            ]

        _total_in = 0
        _total_out = 0
        _MAX_TOOL_ROUNDS = 15

        for _round in range(_MAX_TOOL_ROUNDS + 1):
            kwargs: dict[str, Any] = {
                "model": self.model,
                "messages": api_messages,
                "stream": True,
                "stream_options": {"include_usage": True},
            }

            if is_reasoning:
                kwargs["max_completion_tokens"] = max_tokens
                if thinking_budget >= 64000:
                    kwargs["reasoning_effort"] = "high"
                elif thinking_budget >= 16000:
                    kwargs["reasoning_effort"] = "medium"
                else:
                    kwargs["reasoning_effort"] = "low"
            else:
                kwargs["max_tokens"] = max_tokens

            if (has_doc_tools or has_calc_tools or has_coding) and not is_reasoning:
                kwargs["tools"] = advertised_tools

            # Accumulate streamed tool calls (may arrive in chunks)
            _tool_calls: dict[int, dict] = {}  # index -> {id, name, arguments_parts}
            _text_chunks: list[str] = []

            stream = self.client.chat.completions.create(**kwargs)
            finish_reason = None
            try:
                for chunk in stream:
                    if chunk.usage:
                        _total_in += chunk.usage.prompt_tokens or 0
                        _total_out += chunk.usage.completion_tokens or 0

                    if not chunk.choices:
                        continue

                    choice = chunk.choices[0]
                    delta = choice.delta

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
                        "id": entry["id"],
                        "type": "function",
                        "function": {
                            "name": entry["name"],
                            "arguments": "".join(entry["arguments_parts"]),
                        },
                    })
                assistant_msg["tool_calls"] = tc_list
                api_messages.append(assistant_msg)

                # Execute each tool call and append results
                for tc in tc_list:
                    fn_name = tc["function"]["name"]
                    try:
                        fn_args = json.loads(tc["function"]["arguments"])
                    except json.JSONDecodeError:
                        fn_args = {}

                    # Coding-agent tools get a different namespace prefix so
                    # the UI/Whitelist layer can distinguish them from doc tools.
                    is_coding = fn_name in ("write_file", "edit_file", "multi_edit", "bash")
                    ns_prefix = "kit-coding" if is_coding else "delfin-docs"

                    # Emit tool_use event for UI display
                    yield StreamEvent(
                        type="tool_use",
                        tool_name=f"mcp__{ns_prefix}__{fn_name}",
                        tool_input=json.dumps(fn_args),
                    )

                    result = _doc_executor.execute(
                        fn_name, fn_args, permissions=self._permissions
                    )

                    # Emit tool_result event for UI display
                    yield StreamEvent(
                        type="tool_result",
                        tool_name=f"mcp__{ns_prefix}__{fn_name}",
                        tool_output=result[:2000],
                    )

                    api_messages.append({
                        "role": "tool",
                        "tool_call_id": tc["id"],
                        "content": result,
                    })

                # Loop back to get the model's next response
                continue

            # No tool calls — emit final message_delta and break
            cost = self._estimate_cost(_total_in, _total_out)
            yield StreamEvent(
                type="message_delta",
                input_tokens=_total_in,
                output_tokens=_total_out,
                cost_usd=cost,
                stop_reason=finish_reason or "end_turn",
            )
            break
        else:
            # Exhausted all tool rounds without a final text response.
            # Emit a message_delta so the engine knows streaming is done.
            cost = self._estimate_cost(_total_in, _total_out)
            yield StreamEvent(
                type="message_delta",
                input_tokens=_total_in,
                output_tokens=_total_out,
                cost_usd=cost,
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
        kit_key = api_key or os.environ.get("KIT_TOOLBOX_API_KEY", "")
        kit_workspace = Path(cwd).expanduser().resolve() if cwd else Path.cwd().resolve()
        kit_perms = KitToolPermissions(
            workspace=kit_workspace,
            mode=_map_kit_permission_mode(permission_mode),
            confirm_callback=kit_confirm_callback,
        )
        return OpenAIClient(
            api_key=kit_key, model=model,
            base_url="https://ki-toolbox.scc.kit.edu/api/v1",
            key_env_var="KIT_TOOLBOX_API_KEY",
            permissions=kit_perms,
        )
    if provider == "openai":
        if backend == "cli":
            return CodexCLIClient(model=model, cwd=cwd,
                                  permission_mode=permission_mode)
        openai_key = api_key or os.environ.get("OPENAI_API_KEY", "")
        return OpenAIClient(api_key=openai_key, model=model)
    if backend == "api":
        return APIClient(api_key=api_key, model=model)
    return CLIClient(model=model, claude_path=claude_path,
                     permission_mode=permission_mode, cwd=cwd,
                     mcp_config=mcp_config,
                     allowed_tools=allowed_tools,
                     extra_dirs=extra_dirs)
