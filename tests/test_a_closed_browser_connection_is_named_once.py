"""A kernel message for a closed browser connection is named, once a burst.

Input: the server's stderr as the launcher forwards it. jupyter_server
logs ``self.log.warning(str(e))`` for every kernel message that reaches a
browser websocket that has already closed (a reload, a closed tab, a lost
connection) -- ``WebSocketClosedError`` carries no text, so each such
message printed an empty ``[W <time> ServerApp]`` header. Field report
2026-10-22 18:42 (KIT cluster): a run of them right before the reload,
read as a fault of DELFIN's.

After: the first header-only warning of a burst is printed with what it
means; the rest of that burst is not repeated. Every warning that carries
text still passes through unchanged.
"""
from __future__ import annotations

import io
import sys

from delfin import cli_voila


def _forward(lines: list[str], monkeypatch) -> str:
    raw = io.BytesIO()
    err = io.TextIOWrapper(raw, encoding="utf-8", write_through=True)
    monkeypatch.setattr(sys, "stderr", err)
    cli_voila._forward_filtered_stderr(
        io.BytesIO("".join(lines).encode("utf-8")))
    err.flush()
    return raw.getvalue().decode("utf-8")


def test_a_burst_of_empty_warnings_becomes_one_line_that_says_why(monkeypatch):
    burst = [f"[W 2026-10-22 20:40:50.{i:03d} ServerApp] \n" for i in range(12)]
    out = _forward(burst, monkeypatch)
    lines = [l for l in out.splitlines() if l.strip()]
    assert len(lines) == 1, out
    assert "browser connection" in lines[0] and "closed" in lines[0]
    assert lines[0].startswith("[W 2026-10-22 20:40:50.000 ServerApp]")


def test_a_warning_with_text_still_passes_unchanged(monkeypatch):
    real = "[W 2026-10-22 15:17:56.371 ServerApp] WebSocket ping timeout after 101040 ms.\n"
    out = _forward([real, "[W 2026-10-22 15:17:56.449 ServerApp] \n", real],
                   monkeypatch)
    assert out.count("WebSocket ping timeout after 101040 ms.") == 2
    assert out.count("browser connection") == 1


def test_an_error_header_without_text_is_not_reworded(monkeypatch):
    # Only the warning level is the closed-socket case; an empty error is
    # passed through as it came, not explained with a guess.
    line = "[E 2026-10-22 20:41:06.696 ServerApp] \n"
    assert _forward([line], monkeypatch) == line
