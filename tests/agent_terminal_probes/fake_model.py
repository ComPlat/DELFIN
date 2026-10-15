"""A scripted OpenAI-compatible endpoint, for driving the terminal offline.

Input: HTTP requests in the shape the ``ollama`` provider sends
(``POST /v1/chat/completions`` streamed or not, ``POST /api/show``,
``GET /api/tags``, ``GET /v1/models``). Output: replies chosen by a
marker in the LAST user message, so a test states the model's behaviour
in the prompt it types:

    do-bash <command>      one ``bash`` tool call, then a text answer
    do-bash-secret         a ``bash`` call carrying FAKE_SECRET in a header
    do-write <path> <text> one ``write_file`` tool call, then a text answer
    do-md                  a markdown answer (heading, list, code block)
    do-slow                ~30 s of streamed text, for interrupt tests
    do-delay <s>           headers at once, then <s> seconds before the
                           first token (a connected endpoint prefilling)
    do-hold <s>            <s> seconds before the response headers (a
                           queued endpoint), then a short answer
    do-fail500 / do-fail401   the endpoint answers with that HTTP status
    anything else          "received L lines, C chars; last line: ..."

Every chat request is appended as one JSON line to ``$FAKE_MODEL_LOG``
when that is set, so a test can assert what reached the model (a pasted
block arriving whole, a cancelled turn never arriving) rather than
reading it back off the screen.

Rejected: a monkeypatched client inside the child. It skips the HTTP
layer, the streaming parser and the wait indicator, which are the parts
of a turn a terminal user actually experiences.
"""

from __future__ import annotations

import json
import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

MODEL = "fake-model"
FAKE_SECRET = "sk-live-FAKE0123456789abcdef"


def _last(messages: list[dict], role: str) -> dict | None:
    for msg in reversed(messages):
        if msg.get("role") == role:
            return msg
    return None


def _text_of(msg: dict | None) -> str:
    if not msg:
        return ""
    content = msg.get("content")
    if isinstance(content, list):
        return "".join(str(p.get("text", "")) for p in content
                       if isinstance(p, dict))
    return str(content or "")


def plan(messages: list[dict]) -> dict:
    """What to answer: {"text": str} | {"tool": name, "args": dict}
    | {"status": int} | {"slow": True}."""
    if messages and messages[-1].get("role") == "tool":
        # Only a size: an answer quoting the result would put tool output
        # on screen as ANSWER text, and a test about how the terminal draws
        # the result could no longer tell the two apart.
        return {"text": f"Tool finished ({len(_text_of(messages[-1]))} "
                        "chars came back)."}
    user = _text_of(_last(messages, "user"))
    # The engine may append notes after the typed text; the marker is the
    # first one in the message.
    for line in user.splitlines():
        s = line.strip()
        if s == "do-bash-secret":
            # A credential the user never typed, so any copy of it on the
            # screen was put there by the terminal, not echoed from input.
            # At the END of a short command: a middle-truncated headline
            # keeps the tail, so this is where a narrow line still shows it.
            return {"tool": "bash", "args": {"command": (
                'curl -sH "Authorization: Bearer ' + FAKE_SECRET + '" x')}}
        if s.startswith("do-bash "):
            return {"tool": "bash", "args": {"command": s[8:]}}
        if s.startswith("do-write "):
            _, path, *rest = s.split(" ", 2)
            return {"tool": "write_file",
                    "args": {"path": path, "content": (rest[0] if rest else "")
                             + "\n"}}
        if s == "do-md":
            return {"text": "# Result\n\n- **first** point\n- second point\n\n"
                            "```python\nprint('hi')\n```\n"}
        if s == "do-slow":
            return {"slow": True}
        if s.startswith(("do-delay ", "do-hold ")):
            word, _, secs = s.partition(" ")
            try:
                wait = max(0.0, min(120.0, float(secs)))
            except ValueError:
                wait = 0.0
            key = "delay" if word == "do-delay" else "hold"
            return {"text": f"answered after {secs.strip()} s", key: wait}
        if s in ("do-fail500", "do-fail401"):
            return {"status": int(s[7:])}
    typed = typed_part(user)
    lines = typed.splitlines() or [""]
    return {"text": f"received {len(lines)} lines, {len(typed)} chars; "
                    f"last line: {lines[-1][:40]}"}


# Blocks the engine appends after what was typed (context status, file
# references, the CPU budget). Counting them would make "received N lines"
# describe the harness instead of the input.
_APPENDED = ("\n\nFiles referenced:", "\n\n# Context status",
             "\n\nThis session may use")


def typed_part(user: str) -> str:
    cut = len(user)
    for mark in _APPENDED:
        at = user.find(mark)
        if at != -1:
            cut = min(cut, at)
    return user[:cut]


class _Handler(BaseHTTPRequestHandler):
    def log_message(self, *_args) -> None:     # keep the pty clean
        pass

    def _json(self, code: int, body: dict) -> None:
        data = json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self) -> None:                   # noqa: N802
        if self.path.startswith("/api/tags"):
            self._json(200, {"models": [{"name": MODEL}]})
        elif self.path.startswith("/v1/models"):
            self._json(200, {"object": "list",
                             "data": [{"id": MODEL, "object": "model"}]})
        else:
            self._json(404, {"error": "not found"})

    def do_POST(self) -> None:                  # noqa: N802
        length = int(self.headers.get("Content-Length") or 0)
        try:
            body = json.loads(self.rfile.read(length) or b"{}")
        except ValueError:
            body = {}
        if self.path.startswith("/api/show"):
            self._json(200, {"model_info": {"fake.context_length": 32768},
                             "capabilities": ["completion", "tools"]})
            return
        if not self.path.endswith("/chat/completions"):
            self._json(404, {"error": "not found"})
            return
        messages = body.get("messages") or []
        step = plan(messages)
        log = self.server.log_path or os.environ.get("FAKE_MODEL_LOG")
        if log:
            with open(log, "a", encoding="utf-8") as fh:
                fh.write(json.dumps({"last_user": _text_of(
                    _last(messages, "user")), "plan": step,
                    "n_messages": len(messages),
                    "tool_result": (_text_of(messages[-1])
                                    if messages and messages[-1].get("role")
                                    == "tool" else None)}) + "\n")
        if "status" in step:
            self._json(step["status"], {"error": {
                "message": f"scripted HTTP {step['status']}",
                "type": "scripted"}})
            return
        if step.get("hold"):
            time.sleep(step["hold"])
        if body.get("stream"):
            self._stream(step)
        else:
            self._whole(step)

    # -- answers -----------------------------------------------------------

    def _whole(self, step: dict) -> None:
        msg: dict = {"role": "assistant", "content": step.get("text", "")}
        finish = "stop"
        if "tool" in step:
            msg = {"role": "assistant", "content": None, "tool_calls": [{
                "id": "call_1", "type": "function",
                "function": {"name": step["tool"],
                             "arguments": json.dumps(step["args"])}}]}
            finish = "tool_calls"
        if step.get("slow"):
            time.sleep(30)
            msg["content"] = "slow answer done"
        self._json(200, {"id": "x", "object": "chat.completion",
                         "created": 0, "model": MODEL,
                         "choices": [{"index": 0, "message": msg,
                                      "finish_reason": finish}],
                         "usage": {"prompt_tokens": 100,
                                   "completion_tokens": 10,
                                   "total_tokens": 110}})

    def _chunk(self, delta: dict, finish=None, usage=None) -> None:
        body = {"id": "x", "object": "chat.completion.chunk", "created": 0,
                "model": MODEL,
                "choices": [{"index": 0, "delta": delta,
                             "finish_reason": finish}]}
        if usage is not None:
            body["usage"] = usage
        self.wfile.write(b"data: " + json.dumps(body).encode() + b"\n\n")
        self.wfile.flush()

    def _stream(self, step: dict) -> None:
        self.send_response(200)
        self.send_header("Content-Type", "text/event-stream")
        self.end_headers()
        usage = {"prompt_tokens": 100, "completion_tokens": 10,
                 "total_tokens": 110}
        try:
            self._chunk({"role": "assistant"})
            if step.get("delay"):
                time.sleep(step["delay"])
            if "tool" in step:
                self._chunk({"tool_calls": [{
                    "index": 0, "id": "call_1", "type": "function",
                    "function": {"name": step["tool"],
                                 "arguments": json.dumps(step["args"])}}]})
                self._chunk({}, finish="tool_calls", usage=usage)
            elif step.get("slow"):
                for i in range(60):
                    self._chunk({"content": f"word{i} "})
                    time.sleep(0.5)
                self._chunk({}, finish="stop", usage=usage)
            else:
                text = step.get("text", "")
                for i in range(0, len(text), 16):
                    self._chunk({"content": text[i:i + 16]})
                self._chunk({}, finish="stop", usage=usage)
            self.wfile.write(b"data: [DONE]\n\n")
            self.wfile.flush()
        except (BrokenPipeError, ConnectionResetError):
            pass                                # the client walked away


def serve(port: int = 0, log_path: str = "") -> tuple[ThreadingHTTPServer, int]:
    """Start the endpoint on 127.0.0.1 in a daemon thread; (server, port)."""
    server = ThreadingHTTPServer(("127.0.0.1", port), _Handler)
    server.log_path = log_path
    server.daemon_threads = True
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, server.server_address[1]


if __name__ == "__main__":
    _srv, _port = serve(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    print(_port, flush=True)
    threading.Event().wait()
