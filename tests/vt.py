"""A model terminal for judging what the box draws — with a right margin.

Built 2026-09-21 for the Tab/completion work, shared on purpose: an
 earlier model screen (``test_the_box_does_not_eat_the_transcript.Screen``)
 applied the box's escapes but not its WRAPPING, and on 2026-09-19 that
 gap let two faults through — the model saw a long line as one row, so
 it could not see that the real terminal wrapped it and a row-sized
 erase came up one short. This one wraps at *width* and knows ``CSI J``,
 the erase-from-cursor-to-end the drawing uses when a box shrinks.

It is a test instrument, not an emulator: it implements exactly the
operations the REPL's drawing emits (CR, LF, BS, TAB, SGR-ignore,
OSC-skip, CSI A/B/C/D/G/K/J), plus autowrap at the right margin.
``pyte`` is not installed here; when it is, prefer it.
"""

from __future__ import annotations


class Screen:
    """Rows of characters, a cursor, and the operations the box uses.

    Use it as the ``err`` stream of a ``TerminalAgent`` built with
    ``__new__`` (see ``test_the_box_does_not_eat_the_transcript`` for
    the wiring pattern): everything the agent writes lands on the grid,
    and ``text()`` reads it back row by row.
    """

    def __init__(self, width: int = 80, transcript=()):
        if width < 1:
            raise ValueError("a screen needs a positive width")
        self.width = width
        self.rows: list[list[str]] = [list(r) for r in transcript] or [[]]
        self.row = len(self.rows) - 1
        self.col = 0

    # -- bookkeeping ------------------------------------------------------

    def _fit(self, row: int) -> None:
        while len(self.rows) <= row:
            self.rows.append([])

    def _put(self, ch: str) -> None:
        """One character, wrapping at the right margin like a terminal."""
        if self.col >= self.width:      # the margin: wrap before writing
            self.row += 1
            self.col = 0
        self._fit(self.row)
        line = self.rows[self.row]
        if len(line) < self.col:
            line.extend(" " * (self.col - len(line)))
        if len(line) == self.col:
            line.append(ch)
        else:
            line[self.col] = ch         # overwrite ONE cell, keep the rest
        self.col += 1

    # -- the interpreter --------------------------------------------------

    def write(self, data: str) -> None:
        i = 0
        n = len(data)
        while i < n:
            ch = data[i]
            if ch == "\x1b":
                i = self._escape(data, i)
                continue
            if ch == "\r":
                self.col = 0
            elif ch == "\n":
                self.row += 1
                self._fit(self.row)
            elif ch == "\x7f" or ch == "\b":
                self.col = max(0, self.col - 1)
            elif ch == "\t":
                while self.col % 8:
                    self._put(" ")
                continue
            elif ch < " ":
                pass                    # other C0 controls print nothing
            else:
                self._put(ch)
            i += 1

    def _escape(self, data: str, i: int) -> int:
        """Apply one escape sequence starting at *i*; return the next index."""
        n = len(data)
        if i + 1 >= n:
            return n
        kind = data[i + 1]
        if kind == "]":                 # OSC: skip to BEL or ST
            j = i + 2
            while j < n and data[j] not in ("\x07",):
                if data[j] == "\x1b" and data[j + 1:j + 2] == "\\":
                    return j + 2
                j += 1
            return j + 1
        if kind != "[":
            return i + 2                # a two-character escape (Esc x)
        j = i + 2
        while j < n and not data[j].isalpha() and data[j] not in "?":
            j += 1
        if j >= n:
            return n                    # truncated sequence: drop it
        params = data[i + 2:j]
        verb = data[j]
        if "?" in params:               # private modes (?2004h etc.) — off
            return j + 1
        parts = [p for p in params.split(";")]
        nums = []
        for p in parts:
            nums.append(int(p) if p.isdigit() else 0)
        count = nums[0] if nums and nums[0] else 1
        if verb == "A":
            self.row = max(0, self.row - count)
        elif verb == "B":
            self.row += count
            self._fit(self.row)
        elif verb == "C":
            self.col = min(self.width, self.col + count)
        elif verb == "D":
            self.col = max(0, self.col - count)
        elif verb == "G":
            self.col = max(0, min(self.width - 1, (nums[0] or 1) - 1))
        elif verb == "K":
            self._fit(self.row)
            line = self.rows[self.row]
            mode = nums[0] if parts and parts[0].isdigit() else 0
            if mode == 0:               # cursor to end of line
                if len(line) > self.col:
                    del line[self.col:]
            elif mode == 1:             # start of line to cursor
                keep = line[self.col:]
                self.rows[self.row] = [" "] * self.col + keep
            else:                       # the whole line
                self.rows[self.row] = []
        elif verb == "J":
            mode = nums[0] if parts and parts[0].isdigit() else 0
            self._fit(self.row)
            if mode == 0:               # cursor to end of screen
                line = self.rows[self.row]
                if len(line) > self.col:
                    del line[self.col:]
                for r in range(self.row + 1, len(self.rows)):
                    self.rows[r] = []
            elif mode == 1:             # start of screen to cursor
                for r in range(0, self.row):
                    self.rows[r] = []
                line = self.rows[self.row]
                self.rows[self.row] = [" "] * self.col + line[self.col:]
            else:                       # everything
                for r in range(len(self.rows)):
                    self.rows[r] = []
        # SGR (m) and anything else: recognised and ignored
        return j + 1

    # -- the interface a stream needs -------------------------------------

    def flush(self):
        pass

    def isatty(self):
        return True

    @property
    def writable(self):
        return True

    # -- reading it back ----------------------------------------------------

    def text(self) -> list[str]:
        """The grid, trailing blanks stripped, wrapping already applied."""
        return ["".join(r).rstrip() for r in self.rows]

    def line(self, fragment: str) -> bool:
        """Whether any row *contains* *fragment* — the usual question."""
        return any(fragment in r for r in self.text())
