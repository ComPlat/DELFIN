"""TURBOMOLE manual MCP -- a thin, gated, rate-limited client (package TM2, T2).

TURBOMOLE GmbH runs a public MCP server for its manual
(https://manual.turbomole.org/mcp, Streamable HTTP over DELFIN's own MCP
machinery in ``delfin.agent.mcp_client``). This module is a THIN wrapper on
top of that transport -- it does not re-implement an HTTP stack, it drives
``MCPServer.call_tool``. Its tools are ``input_choices`` and
``generate_input`` (a ``control`` file + run commands + success criteria for
one structure).

Licence and privacy rules this code enforces (the person, 06.11.):
- DELFIN never ships, downloads or copies TURBOMOLE binaries, basis sets or
  parameter files, and never copies manual text or an MCP-generated
  ``control`` into the repo. A generated ``control`` is used at runtime for
  that one calculation only (``write_generated_control`` writes into a NEW
  calculation folder and records that provenance).
- A structure goes to the external server ONLY when the user's settings allow
  it (``runtime.turbomole.mcp_allowed``, default OFF). When it is OFF, the
  client refuses to call ``generate_input`` -- nothing is sent, and DELFIN
  builds the ``control`` itself (see ``delfin/tools/adapters/turbomole.py``).
- The MCP describes <= 60 calls/minute; this client holds a >= 1 call / 1.1 s
  wire cap so it stays well under the limit, and surfaces a rate-limited
  answer rather than dropping it.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Optional, Sequence


class McpNotAllowed(RuntimeError):
    """Raised when a structure would be sent but ``turbomole.mcp_allowed`` is OFF."""


def mcp_allowed(settings_path: Optional[str] = None) -> bool:
    """Is the TURBOMOLE manual MCP allowed for this user?

    Reads ``runtime.turbomole.mcp_allowed`` from the user's settings through
    the established ``user_settings.load_runtime_settings`` path. Absent or
    unreadable -> False. The feature is OFF by default; nothing is sent while
    it is OFF.
    """
    from delfin import user_settings
    try:
        runtime = user_settings.load_runtime_settings(settings_path)
    except Exception:
        return False
    turbomole = runtime.get("turbomole", {}) if isinstance(runtime, dict) else {}
    return bool(turbomole.get("mcp_allowed", False))


class _RateLimiter:
    """Holds calls to >= ``min_interval_s`` apart (<= 60/min wire cap).

    ``now_fn`` and ``sleep_fn`` are injected so tests are deterministic and
    never touch real wall-clock -- no test may depend on machine speed.
    """

    def __init__(
        self,
        min_interval_s: float,
        *,
        now_fn: Optional[Callable[[], float]] = None,
        sleep_fn: Optional[Callable[[float], None]] = None,
    ) -> None:
        self.min_interval_s = float(min_interval_s)
        self._now = now_fn or time.monotonic
        self._sleep = sleep_fn or time.sleep
        self._last = 0.0

    def wait(self) -> float:
        """Block until the next call may go out; returns the call time."""
        now = self._now()
        elapsed = now - self._last
        if self._last and elapsed < self.min_interval_s:
            delay = self.min_interval_s - elapsed
            self._sleep(delay)
            now = self._now()
        self._last = now
        return self._last


@dataclass
class GeneratedInput:
    """The parsed result of an MCP ``generate_input`` call.

    ``control_text`` is the ``control`` DELFIN would write into a calculation
    folder; ``raw_text`` is whatever the wrapped client returned, kept so a
    rate-limit or error answer is never silently dropped. Only the parts
    DELFIN can map are surfaced -- an MCP answer proposing something outside
    DELFIN's vocabulary is caught later by the version/compatibility check.
    """

    control_text: str = ""
    raw_text: str = ""
    command_hints: tuple[str, ...] = field(default_factory=tuple)
    success_criteria: tuple[str, ...] = field(default_factory=tuple)


class TurbomoleManualMcp:
    """Thin client for the TURBOMOLE manual MCP (input_choices, generate_input).

    ``client`` is anything with ``call_tool(name, arguments) -> str`` -- in
    DELFIN that is ``delfin.agent.mcp_client.MCPServer`` (built from the
    user's ``mcp_servers.json`` entry for the manual server); in tests it is a
    fake that returns canned answers, so no live network is ever touched.
    """

    def __init__(
        self,
        client: Any,
        *,
        allowed: bool,
        min_interval_s: float = 1.1,
        now_fn: Optional[Callable[[], float]] = None,
        sleep_fn: Optional[Callable[[float], None]] = None,
    ) -> None:
        self.client = client
        self.allowed = bool(allowed)
        self.min_interval_s = float(min_interval_s)
        self._limiter = _RateLimiter(self.min_interval_s, now_fn=now_fn, sleep_fn=sleep_fn)

    def input_choices(self, arguments: Optional[dict] = None) -> str:
        """Ask the MCP for the valid choices for an input (read-only hint)."""
        return self._call("input_choices", dict(arguments or {}))

    def generate_input(self, arguments: Optional[dict] = None) -> GeneratedInput:
        """Ask the MCP to build a ``control`` for one structure.

        Refuses (McpNotAllowed) while the feature is OFF -- a structure is
        sent only when ``turbomole.mcp_allowed`` is on. The answer is parsed
        into a ``GeneratedInput``; a rate-limited or error answer is surfaced,
        not dropped.
        """
        if not self.allowed:
            raise McpNotAllowed(
                "turbomole.mcp_allowed is OFF: a structure would be sent to the "
                "external TURBOMOLE manual server. Enable it in settings first."
            )
        raw = self._call("generate_input", dict(arguments or {}))
        return _split_generated(raw)

    def _call(self, tool_name: str, arguments: dict) -> str:
        self._limiter.wait()
        return self.client.call_tool(tool_name, arguments)


def _split_generated(raw: str) -> GeneratedInput:
    """Parse the flattened MCP answer into a ``GeneratedInput``.

    The wrapped client returns the flattened ``content`` text. When that text
    is a JSON object (a tool result), the ``control`` / ``commands`` /
    ``success_criteria`` fields are pulled out; otherwise the whole text IS
    the control. Either way ``raw_text`` carries the full answer.
    """
    if raw is None:
        return GeneratedInput(raw_text="")
    text = str(raw)
    try:
        obj = json.loads(text)
    except (ValueError, TypeError):
        return GeneratedInput(control_text=text, raw_text=text)
    if isinstance(obj, dict):
        control = obj.get("control") or obj.get("control_file") or obj.get("control_text") or ""
        commands = obj.get("commands") or obj.get("command_lines") or ()
        success = obj.get("success_criteria") or obj.get("success") or ()
        return GeneratedInput(
            control_text=str(control) if control else "",
            raw_text=text,
            command_hints=tuple(commands),
            success_criteria=tuple(success),
        )
    return GeneratedInput(control_text=text, raw_text=text)


def write_generated_control(folder: Path, generated: GeneratedInput) -> Path:
    """Write ``generated.control_text`` into a NEW calculation ``folder``.

    Also records a provenance note naming the source, so a later reader knows
    the ``control`` came from the TURBOMOLE manual MCP and is for this one
    calculation only -- never something DELFIN ships. Returns the folder.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "control").write_text(generated.control_text, encoding="utf-8")
    provenance_lines = [
        "# This control was generated by the TURBOMOLE manual MCP (external "
        "service) for the calculation in this folder only.",
        "# It is not distributed with DELFIN.",
    ]
    (folder / "mcp_provenance.txt").write_text(
        "\n".join(provenance_lines) + "\n", encoding="utf-8"
    )
    return folder
