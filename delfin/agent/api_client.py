"""Client backends for the DELFIN Agent: Claude CLI, Anthropic API, OpenAI API, or Codex CLI."""

from __future__ import annotations

import importlib
import json
import os
import shutil
import subprocess
import sys
from typing import Any, Generator


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
]


class _DocToolExecutor:
    """Lazy-loaded local executor for doc search tools."""

    def __init__(self) -> None:
        self._engine = None
        self._index: dict | None = None

    def _ensure_loaded(self) -> bool:
        """Load index and build search engine. Returns True if ready."""
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

    def execute(self, name: str, arguments: dict) -> str:
        """Execute a doc tool by name. Returns the result string."""
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
                 base_url: str = "", key_env_var: str = "OPENAI_API_KEY"):
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

        # Check if doc tools are available
        has_doc_tools = _doc_executor._ensure_loaded()

        _total_in = 0
        _total_out = 0
        _MAX_TOOL_ROUNDS = 5

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

            if has_doc_tools and not is_reasoning:
                kwargs["tools"] = _DOC_TOOLS_OPENAI

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

                    # Emit tool_use event for UI display
                    yield StreamEvent(
                        type="tool_use",
                        tool_name=f"mcp__delfin-docs__{fn_name}",
                        tool_input=json.dumps(fn_args),
                    )

                    result = _doc_executor.execute(fn_name, fn_args)

                    # Emit tool_result event for UI display
                    yield StreamEvent(
                        type="tool_result",
                        tool_name=f"mcp__delfin-docs__{fn_name}",
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
        "acceptEdits":         ["--full-auto"],  # = workspace-write + auto-approve
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
        return OpenAIClient(
            api_key=kit_key, model=model,
            base_url="https://ki-toolbox.scc.kit.edu/api/v1",
            key_env_var="KIT_TOOLBOX_API_KEY",
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
