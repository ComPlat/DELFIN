"""The TUI protocol's framing: one JSON object per line, nothing else.

The TypeScript client splits the stream on newlines, so a raw newline
inside a message would cut it in two, and one bad line must not end a
session. These are the rules both sides rely on.
"""

from __future__ import annotations

import io
import json
import threading

from delfin.agent.service import protocol


def test_a_newline_inside_a_value_never_ends_the_line():
    data = protocol.encode(protocol.event({"type": "notice",
                                           "text": "line one\nline two\r\n"}))
    assert data.endswith(b"\n")
    assert data.count(b"\n") == 1
    assert json.loads(data)["params"]["text"] == "line one\nline two\r\n"


def test_text_travels_as_text_not_as_escapes():
    data = protocol.encode({"text": "Energie ≈ −76.4 Eₕ 🚀"})
    assert "Energie ≈ −76.4 Eₕ 🚀".encode() in data


def test_a_bad_line_is_reported_and_the_next_one_still_read():
    stream = io.BytesIO(b'{"a": 1}\nnot json\n[1, 2]\n\n{"b": 2}\n')
    problems: list[str] = []
    got = list(protocol.read_messages(stream, on_bad_line=problems.append))
    assert got == [{"a": 1}, {"b": 2}]
    assert len(problems) == 2


def test_an_oversized_line_is_skipped_without_losing_the_next():
    big = b'{"x": "' + b"y" * 5000 + b'"}\n'
    stream = io.BytesIO(big + b'{"ok": true}\n')
    problems: list[str] = []
    got = list(protocol.read_messages(stream, on_bad_line=problems.append,
                                      max_bytes=1024))
    assert got == [{"ok": True}]
    assert problems and "longer than" in problems[0]


def test_params_are_checked_against_the_one_definition():
    assert protocol.invalid_params("session.prompt", {"text": "hi"}) == ""
    assert "missing" in protocol.invalid_params("session.prompt", {})
    assert "must be string" in protocol.invalid_params("session.prompt", {"text": 3})
    assert protocol.invalid_params("permission.reply",
                                   {"requestId": "p1", "answers": ["a", "b"]}) == ""
    assert "must be [string]" in protocol.invalid_params(
        "permission.reply", {"requestId": "p1", "answers": ["a", 2]})
    # A boolean is not an integer, whatever Python thinks.
    assert protocol.invalid_params("initialize", {"protocolVersion": True})
    assert "unknown method" in protocol.invalid_params("nope", {})


def test_the_schema_lists_every_message_and_every_durable_event_exists():
    schema = protocol.json_schema()
    assert set(schema["requests"]) == set(protocol.REQUESTS)
    assert set(schema["results"]) == set(protocol.REQUESTS)
    assert set(schema["events"]) == set(protocol.EVENTS)
    assert protocol.DURABLE_EVENTS <= set(protocol.EVENTS)
    for name, event_schema in schema["events"].items():
        assert event_schema["properties"]["type"] == {"const": name}
        for common in ("session", "ts"):
            assert common in event_schema["required"]
        # Only what is kept is numbered (eventlog.py).
        assert ("seq" in event_schema["required"]) == \
            (name in protocol.DURABLE_EVENTS)


def test_a_file_name_that_is_not_utf8_cannot_break_a_line():
    import os
    odd = os.fsdecode(b"caf\xe9.xyz")          # a lone surrogate inside
    data = protocol.encode(protocol.event({"type": "notice", "text": odd}))
    assert data.count(b"\n") == 1
    assert json.loads(data)["params"]["text"] == odd
    data.decode("utf-8")                       # the line itself is valid UTF-8


def test_a_value_json_has_no_type_for_travels_as_its_text():
    from pathlib import Path as _P
    data = protocol.encode({"path": _P("/tmp/x"), "raw": b"\x00"})
    assert json.loads(data) == {"path": "/tmp/x", "raw": "b'\\x00'"}


def test_a_line_written_in_pieces_arrives_whole():
    class Trickle(io.RawIOBase):
        """A pipe that takes at most 7 bytes a call, and now and then none."""

        def __init__(self):
            self.data = bytearray()
            self.calls = 0

        def writable(self):
            return True

        def write(self, b):
            self.calls += 1
            if self.calls % 3 == 0:
                return None                    # would block
            chunk = bytes(b[:7])
            self.data += chunk
            return len(chunk)

    sink = Trickle()
    writer = protocol.LineWriter(sink)
    assert writer.send(protocol.event({"type": "notice", "text": "x" * 100}))
    assert json.loads(bytes(sink.data))["params"]["text"] == "x" * 100


def test_threads_writing_at_once_never_interleave_a_line():
    sink = io.BytesIO()
    writer = protocol.LineWriter(sink)

    def burst(n: int) -> None:
        for i in range(200):
            writer.send(protocol.event({"type": "notice", "text": f"{n}:{i}" * 20}))

    threads = [threading.Thread(target=burst, args=(n,)) for n in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    lines = sink.getvalue().splitlines()
    assert len(lines) == 1600
    assert all(json.loads(line)["method"] == "event" for line in lines)


def test_a_writer_whose_reader_is_gone_stops_quietly():
    class Gone(io.RawIOBase):
        def writable(self):
            return True

        def write(self, _b):
            raise BrokenPipeError

    writer = protocol.LineWriter(Gone())
    assert writer.send({"a": 1}) is False
    assert writer.closed and writer.send({"a": 2}) is False
