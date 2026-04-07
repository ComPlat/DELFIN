"""Client backends for the DELFIN Agent: Claude Code CLI or Anthropic API."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from typing import Any, Generator


# ---------------------------------------------------------------------------
# Shared event type
# ---------------------------------------------------------------------------

class StreamEvent:
    """A single event from a streaming Claude response."""

    __slots__ = ("type", "text", "input_tokens", "output_tokens", "stop_reason",
                 "cost_usd", "tool_name", "tool_input")

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
    ):
        self.type = type
        self.text = text
        self.input_tokens = input_tokens
        self.output_tokens = output_tokens
        self.stop_reason = stop_reason
        self.cost_usd = cost_usd
        self.tool_name = tool_name
        self.tool_input = tool_input


# ---------------------------------------------------------------------------
# Backend base
# ---------------------------------------------------------------------------

class _BaseClient:
    """Common interface for both backends."""

    def stream_message(
        self,
        system: str,
        messages: list[dict[str, Any]],
        max_tokens: int = 4096,
        session_id: str = "",
        thinking_budget: int = 0,
    ) -> Generator[StreamEvent, None, None]:
        raise NotImplementedError


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
                 permission_mode: str = "", cwd: str = ""):
        self.model = model or self.DEFAULT_MODEL
        self.permission_mode = permission_mode
        self.cwd = cwd or None
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
            "--model", self.model,
            "--append-system-prompt", system,
        ]

        if self.permission_mode and self.permission_mode != "default":
            cmd.extend(["--permission-mode", self.permission_mode])

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
        max_tokens: int = 4096,
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
                    if block.get("type") == "tool_use":
                        current_tool_name = block.get("name", "")
                        current_tool_input_chunks = []
                    elif block.get("type") == "thinking":
                        in_thinking_block = True

                elif etype == "content_block_delta":
                    delta = evt.get("delta", {})
                    if delta.get("type") == "text_delta":
                        text = delta.get("text", "")
                        if text:
                            emitted_text = True
                            yield StreamEvent(type="text_delta", text=text)
                    elif delta.get("type") == "thinking_delta":
                        text = delta.get("thinking", "")
                        if text:
                            yield StreamEvent(type="thinking_delta", text=text)
                    elif delta.get("type") == "input_json_delta":
                        current_tool_input_chunks.append(
                            delta.get("partial_json", "")
                        )

                elif etype == "content_block_stop":
                    if in_thinking_block:
                        in_thinking_block = False
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
                    if block.get("type") == "tool_use":
                        tool_name = block.get("name", "")
                        tool_input = json.dumps(
                            block.get("input", {}), ensure_ascii=False
                        )
                        if not current_tool_name:
                            yield StreamEvent(
                                type="tool_use",
                                tool_name=tool_name,
                                tool_input=tool_input,
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

        # If we exit the loop, the process died
        rc = proc.poll()
        if rc is not None and rc != 0:
            stderr = proc.stderr.read() if proc.stderr else ""
            self._proc = None
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
            raise ImportError(
                "The 'anthropic' package is required for API mode. "
                "Install with: pip install 'anthropic>=0.40'"
            )
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
        max_tokens: int = 4096,
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
# Factory
# ---------------------------------------------------------------------------

def create_client(
    backend: str = "cli",
    api_key: str = "",
    model: str = "",
    claude_path: str = "",
    permission_mode: str = "",
    cwd: str = "",
) -> _BaseClient:
    """Create the appropriate client backend.

    Parameters
    ----------
    backend : str
        ``"cli"`` (default) or ``"api"``.
    api_key : str
        Only needed for ``"api"`` backend.
    model : str
        Model name/alias.
    claude_path : str
        Only for ``"cli"`` backend.
    permission_mode : str
        CLI permission mode (``"default"``, ``"acceptEdits"``, ``"plan"``,
        ``"auto"``, ``"bypassPermissions"``).
    cwd : str
        Working directory for the CLI process.
    """
    if backend == "api":
        return APIClient(api_key=api_key, model=model)
    return CLIClient(model=model, claude_path=claude_path,
                     permission_mode=permission_mode, cwd=cwd)
