"""R1 / phase 3 — wake on message for the readline (non-raw) terminal path.

Wave-12 finding 3: an incoming session_message did not start a turn in an
idle terminal session; it was only read after a key press. The raw-terminal
path already wakes (read_boxed polls every 0.5s -> _wake_text ->
_operator_messages), but the readline fallback (run() line 2135
read_block(self._read_line)) blocks on input() forever, so an idle session on
that path never sees a message until a key is pressed (verified repl.py:740
_input -> input(), no select/timer). These tests pin the new wakeable read for
that path.
"""

from delfin.agent.repl import read_block_wakeable


class Scripted:
    """A tiny scripted callable returning its list one item per call."""

    def __init__(self, values):
        self._values = list(values)
        self.calls = 0

    def __call__(self, *args, **kwargs):
        i = self.calls
        self.calls += 1
        if i < len(self._values):
            return self._values[i]
        return self._values[-1]


class TestWakeableRead:
    def test_wake_message_becomes_the_prompt(self):
        # Idle (ready always False); a session message arrives on the first
        # wake check -> it becomes the prompt, and read_line is never called.
        read_line = Scripted(["hello"])
        wake = Scripted(["✉ from s12: please look"])
        ready = Scripted([False, False])
        text, from_wake = read_block_wakeable(
            read_line, wake, ready=ready, tick=0)
        assert from_wake is True
        assert "s12" in text
        assert read_line.calls == 0      # never blocked on a key press

    def test_keypress_takes_a_normal_read(self):
        # A key is already waiting -> a real read_line, not a wake.
        read_line = Scripted(["hello"])
        wake = Scripted(["✉ should not win"])
        ready = Scripted([True])
        text, from_wake = read_block_wakeable(
            read_line, wake, ready=ready, tick=0)
        assert from_wake is False
        assert text == "hello"
        assert wake.calls == 0

    def test_polls_until_a_wake_arrives(self):
        wedges = Scripted([False, False, False])
        wake = Scripted(["", "", "✉ m"])     # nothing, nothing, then a message
        text, from_wake = read_block_wakeable(
            Scripted([]), wake, ready=wedges, tick=0)
        assert from_wake is True
        assert text == "✉ m"

    def test_wake_does_not_fire_when_a_message_string_is_blank(self):
        # Blank/empty wake results are not a prompt; the loop keeps waiting.
        wake = Scripted(["", "✉ go"])
        ready = Scripted([False, False, False])
        text, from_wake = read_block_wakeable(
            Scripted([]), wake, ready=ready, tick=0)
        assert from_wake is True and text == "✉ go"

    def test_backslash_continuation_after_a_real_first_line(self):
        # Only the FIRST line is wake-aware; the read_block tail still joins
        # backslash continuations through read_line.
        lines = Scripted(["abc\\", "def"])
        wake = Scripted([])
        ready = Scripted([True])
        text, from_wake = read_block_wakeable(
            lines, wake, ready=ready, tick=0)
        assert from_wake is False
        assert text == "abcdef"
        assert lines.calls == 2
