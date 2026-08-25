"""The terminal controller: one turn, rendered as it happens.

Division of labour, and the reason for it: ``repl_render`` decides what a
line says, this module decides when it is written and to which stream.
Keeping the two apart is what stops this file growing into the thing it
replaces — a single closure where formatting, threading and state are the
same 14 000 lines.

Stream discipline: **stdout is the answer, stderr is everything else.**
Tool lines, notices, thinking, the banner and the prompt all go to stderr,
so ``delfin-agent -p "..." > answer.txt`` produces exactly the answer and
nothing about how it was reached.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Any, Callable

from . import repl_render as rr

__all__ = ["RenderItem", "Transcript", "TurnResult", "run_turn", "TURN_KEYS"]


# The keys `cli._run_once` returns, and therefore what --output-format json
# promises. Named here so the two paths can be asserted equal instead of
# drifting into two shapes of the same answer.
TURN_KEYS: frozenset[str] = frozenset({
    "text", "tool_calls", "input_tokens", "output_tokens", "error",
})


@dataclass(frozen=True)
class RenderItem:
    """One thing that happened, on its way from the worker to the screen."""

    kind: str                      # text|thinking|tool_use|tool_result
                                   # |notice|denied|error|done
    text: str = ""
    name: str = ""
    data: dict | None = None


class Transcript:
    """Owns both streams and the cursor column.

    The column tracking is not decoration. The model streams
    ``"Let me check the tests"`` with no trailing newline and then calls a
    tool; without closing the line first, the tool line is appended to the
    middle of the sentence and the transcript stops being readable exactly
    when it matters most.
    """

    def __init__(self, out=None, err=None, *, theme: rr.Theme | None = None,
                 show_tools: bool = True, show_thinking: bool = False,
                 width: int | None = None) -> None:
        self.out = out if out is not None else sys.stdout
        self.err = err if err is not None else sys.stderr
        self.theme = theme if theme is not None else rr.theme_for(self.err)
        self.show_tools = show_tools
        self.show_thinking = show_thinking
        self.width = width or rr.terminal_width()
        # Two flags, because one had to answer two different questions and
        # got one of them wrong. `_out_open` is about the FILE: stdout ends
        # mid-line and needs closing before the process exits.
        # `_break_pending` is about the SCREEN: the next chrome line needs a
        # visual break first, and only the first one does. Sharing a flag
        # meant the first tool line cleared it and stdout was left without
        # its final newline.
        self._out_open = False
        self._break_pending = False

    # -- primitives ------------------------------------------------------
    def answer(self, delta: str) -> None:
        """Answer text. The only thing that ever reaches stdout."""
        if not delta:
            return
        self.out.write(delta)
        self._flush(self.out)
        self._out_open = not delta.endswith("\n")
        self._break_pending = self._out_open

    def chrome(self, line: str) -> None:
        """Everything that is not the answer, on stderr, on its own line.

        The break that closes a half-written answer line is written to
        STDERR, not stdout. Both streams share one cursor when they share a
        terminal, so the visual effect is identical — but a newline on
        stdout would land inside the answer, and `delfin-agent -p "..." >
        answer.txt` would capture a line break that the model never wrote.
        Redirecting stderr alone costs the visual break and keeps the
        answer exact, which is the right way round.
        """
        if not line:
            return
        if self._break_pending:
            self.err.write("\n")
            self._break_pending = False
        self.err.write(line + "\n")
        self._flush(self.err)

    def finish(self) -> None:
        """Close the answer stream, so a redirected stdout ends in a newline."""
        if self._out_open:
            self.out.write("\n")
            self._flush(self.out)
            self._out_open = False
            self._break_pending = False

    @staticmethod
    def _flush(stream) -> None:
        try:
            stream.flush()
        except Exception:
            pass

    # -- dispatch --------------------------------------------------------
    def render(self, item: RenderItem) -> None:
        kind = item.kind
        if kind == "text":
            self.answer(item.text)
        elif kind == "notice":
            self.chrome(rr.notice_line(item.text, theme=self.theme))
        elif kind == "thinking":
            if self.show_thinking:
                self.chrome(rr.thinking_line(
                    item.text, width=self.width, theme=self.theme))
        elif kind == "tool_use":
            if self.show_tools:
                self.chrome(rr.tool_headline(
                    item.name, item.text, width=self.width, theme=self.theme))
        elif kind == "tool_result":
            if self.show_tools:
                self.chrome(rr.tool_result_line(
                    item.name, item.text, meta=item.data,
                    width=self.width, theme=self.theme))
        elif kind == "denied":
            self.chrome(rr.denied_line(item.name, theme=self.theme))
        elif kind == "error":
            self.chrome(self.theme.red(f"[error] {item.text}"))


@dataclass
class TurnResult:
    text: str = ""
    tool_calls: list[dict] = field(default_factory=list)
    input_tokens: int = 0
    output_tokens: int = 0
    error: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "text": self.text,
            "tool_calls": self.tool_calls,
            "input_tokens": self.input_tokens,
            "output_tokens": self.output_tokens,
            "error": self.error,
        }


def _usage(engine) -> tuple[int, int]:
    usage = getattr(engine, "token_usage", {}) or {}
    return int(usage.get("input", 0) or 0), int(usage.get("output", 0) or 0)


def run_turn(engine, prompt: str, *, sink: Callable[[RenderItem], None],
             max_tokens: int = 0, memory_context: str = "") -> TurnResult:
    """One turn, with every callback turned into a RenderItem.

    The accounting mirrors ``cli._run_once`` deliberately — same token
    deltas, same result keys — so the interactive and the headless path
    report the same turn the same way.

    ``sink`` is called from whatever thread ``stream_response`` runs on.
    It must not touch the terminal; in the REPL it is a queue put, and the
    main thread does the writing.
    """
    chunks: list[str] = []
    tool_calls: list[dict] = []
    result = TurnResult()

    def _emit(item: RenderItem) -> None:
        try:
            sink(item)
        except Exception:
            pass

    def _on_token(text: str) -> None:
        if text:
            chunks.append(text)
            _emit(RenderItem("text", text=text))

    def _on_notice(text: str) -> None:
        if text:
            _emit(RenderItem("notice", text=text))

    def _on_thinking(text: str) -> None:
        if text:
            _emit(RenderItem("thinking", text=text))

    def _on_tool_use(name: str, input_json: str) -> None:
        tool_calls.append({"name": name, "input": rr._as_dict(input_json)})
        _emit(RenderItem("tool_use", name=name, text=input_json or ""))

    def _on_tool_result(name: str, output: str) -> None:
        _emit(RenderItem("tool_result", name=name, text=output or ""))

    def _on_tool_result_meta(name: str, meta: dict) -> None:
        _emit(RenderItem("tool_result", name=name, text="", data=dict(meta)))

    def _on_denied(name: str) -> None:
        _emit(RenderItem("denied", name=name))

    in_before, out_before = _usage(engine)
    try:
        full_text = engine.stream_response(
            user_message=prompt,
            memory_context=memory_context,
            on_token=_on_token,
            on_notice=_on_notice,
            on_thinking=_on_thinking,
            on_tool_use=_on_tool_use,
            on_tool_result=_on_tool_result,
            on_tool_result_meta=_on_tool_result_meta,
            on_permission_denied=_on_denied,
            max_tokens=max_tokens,
        ) or ""
    except Exception as exc:
        result.error = str(exc)
        full_text = ""
        _emit(RenderItem("error", text=result.error))

    in_after, out_after = _usage(engine)
    result.text = (full_text or "".join(chunks)).strip()
    result.tool_calls = tool_calls
    result.input_tokens = max(0, in_after - in_before)
    result.output_tokens = max(0, out_after - out_before)

    # The answer was streamed before the sanitiser ran on it, and a
    # terminal cannot un-print. Say so rather than let the screen and the
    # recorded answer quietly disagree.
    try:
        from .text_sanitize import sanitize_agent_text
        cleaned = sanitize_agent_text(result.text).text
        if cleaned != result.text:
            result.text = cleaned
            _emit(RenderItem("notice", text=(
                "The model leaked tool-channel text into the answer above; "
                "the recorded answer is the cleaned version.")))
    except Exception:
        pass

    _emit(RenderItem("done"))
    return result
