"""Pure rendering of the input area at the idle prompt.

No terminal, no I/O, no cursor movement: in goes (text, cursor offset,
width, hint), out comes a list of rows. Everything hard about an input
line is layout, and layout is decidable without a terminal — which is
the whole reason this module exists apart from ``repl.py``.

Two rules, and the text between them:

    ─────────────────────────────────────────────────────────
    > the message being typed
    ─────────────────────────────────────────────────────────
      esc interrupt · shift+tab approval mode · /help

It was a closed frame first. A frame has to defend its right edge on
every keystroke and at every width, and three of the four faults found
in the first week lived in exactly that: a border drawn from the new
height while the old one stood, a bottom border counted as content, a
cursor offset by a border that a narrow row did not draw. The rules say
the same thing to the eye and give the text four more columns.

Rules the renderer owns:

* Wrapping. Text longer than one row wraps onto further rows; the
  cursor travels with it, and ``BoxView.cursor`` says where on the
  screen it belongs, as (content-row index, column from the row's
  start), so the caller can move it without re-deriving the layout.
* Double-width characters. CJK and other East Asian wide characters
  occupy two columns; a wide char that would straddle the edge moves
  down to the next row whole rather than being split, because a
  terminal cannot put half of one in a column.
* Width below the minimum. Below ``MIN_WIDTH`` even the rules are not
  worth their rows: the area degenerates to a single row with the END
  of the text kept, so nothing is silently dropped and the cursor
  (where a person is typing) stays visible.
* A hint that does not fit. The hint is truncated from the left with an
  ellipsis, keeping its END — the key names, which is what a hint is
  for. It never wraps and never steals a row from the text.
* Nothing writes the last column. Every row is at most ``width - 1``
  wide, so no terminal wraps a row of its own accord.
"""

from __future__ import annotations

import unicodedata

__all__ = ["BoxView", "render_box", "MIN_WIDTH", "PROMPT",
           "char_width", "string_width"]

#: The prompt inside the box, drawn before the text on the first row.
PROMPT = "> "

#: Below this width the frame does not fit; a plain row is returned.
MIN_WIDTH = 20

_BOTTOM_LEFT, _BOTTOM_RIGHT = "╰", "╯"
_HORIZONTAL = "─"


def char_width(ch: str) -> int:
    """Columns *ch* occupies. Wide (CJK) and fullwidth characters are 2.

    Combining marks are 0 — they stack onto the character before them,
    and counting one for each would push the row past the rules
    exactly when the user types Vietnamese or Korean.
    """
    if unicodedata.combining(ch):
        return 0
    if ch == "\n":
        # A newline takes no column of its own: _wrap keeps it at the
        # end of a row as a marker, and the padding must not count it.
        return 0
    return 2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1


def string_width(text: str) -> int:
    return sum(char_width(ch) for ch in text)


def _visible_text(text: str) -> str:
    """Make pasted controls inert without changing character offsets.

    The decoder deliberately keeps every byte inside bracketed paste: the
    message sent to the model must remain what the user pasted.  The box is
    a different surface.  Writing ESC, CR or a tab into a terminal there can
    clear the screen, overwrite the prompt, or move the cursor somewhere the
    painter did not account for.  Each control therefore becomes exactly one
    printable character here, so cursor offsets into the original text still
    address the same position in the rendered copy.
    """
    shown: list[str] = []
    for ch in text or "":
        code = ord(ch)
        if ch == "\n":
            shown.append(ch)                 # a deliberate hard wrap
        elif ch == "\t":
            shown.append("⇥")                # visible horizontal tab
        elif code < 32:
            shown.append(chr(0x2400 + code))  # ESC -> ␛, CR -> ␍
        elif code == 127:
            shown.append("␡")
        elif 0x80 <= code < 0xA0:
            shown.append("�")                # C1 terminal controls
        else:
            shown.append(ch)
    return "".join(shown)


def _wrap(text: str, inner: int) -> list[tuple[str, int]]:
    """Wrap *text* to *inner* columns; one (row, width) pair per row.

    Wrapping is purely graphical — by columns, not by words — because
    this is an input line: the user's spacing is content, and reflowing
    it would move what the cursor points at.

    A word longer than the row is split across rows rather than
    overflowing, and a double-width character that does not fit in the
    last column moves down WHOLE: a terminal cannot render half of one,
    and the alternative — a one-column gap — desynchronises the right
    border from every other row.
    """
    if inner <= 0:
        return [("", 0)]
    rows: list[tuple[str, int]] = []
    current = ""
    used = 0
    for ch in text:
        if ch == "\n":
            # A pasted block keeps its newlines: they hard-wrap the row,
            # exactly as the terminal would render them. The newline
            # stays at the END of the row it ends — zero-width, so the
            # padding ignores it — because the cursor walk needs to tell
            # a hard break from a wrap edge.
            rows.append((current + "\n", used))
            current, used = "", 0
            continue
        w = char_width(ch)
        if w == 0:
            # A combining mark stacks onto the character before it; it
            # costs no column and can never force a wrap.
            current += ch
            continue
        if used + w > inner:
            rows.append((current, used))
            current, used = "", 0
        current += ch
        used += w
    rows.append((current, used))
    return rows


def _truncate_hint(hint: str, width: int) -> str:
    """Cut *hint* to *width* columns, keeping its END.

    The end of the hint is the key names — the part a person looks at —
    and the left side is prose, which is the cheaper half to lose.
    """
    if string_width(hint) <= width:
        return hint
    keep = "…" + hint
    while string_width(keep) > width and len(keep) > 1:
        keep = "…" + keep[2:]
    return keep


def _truncate_status(status: str, width: int) -> str:
    """Cut *status* to *width* columns, keeping its START.

    The status line is not prose-with-keys: its start carries the ⚙ and
    — while a walk is on — the ▶ that says which row Enter would take.
    Hint-style end-keeping moved that mark off the screen on exactly
    the busy line the walk is for. The handles later in the line are
    reached by walking to them; the mark is not.
    """
    if string_width(status) <= width:
        return status
    keep = status
    while string_width(keep) > width - 1 and len(keep) > 1:
        keep = keep[:-1]
    return keep + "…"


def _narrow_row(text: str, cursor: int, width: int) -> BoxView:
    """The below-MIN_WIDTH form: one safe, cursor-following row.

    It is a horizontal viewport, biased toward what precedes the cursor.
    Keeping the unconditional END made a cursor moved into the middle point
    at unrelated text.  The budget is ``width - 1`` for the same reason as
    the full box: writing the terminal's final column enables autowrap and
    silently turns this one-row fallback into two physical rows.
    """
    budget = max(0, int(width or 0) - 1)
    if budget <= 0:
        return BoxView([""], (0, 0), None, border=False)

    # A newline cannot exist in a one-row view.  Keep it visible and keep
    # its one-character offset, just as _visible_text does for other
    # controls.
    full = (PROMPT + (text or "")).replace("\n", "↵")
    stop = max(0, min(len(PROMPT) + cursor, len(full)))
    if string_width(full) <= budget:
        start, end = 0, len(full)
    else:
        start = stop
        used = 0
        while start > 0:
            step = char_width(full[start - 1])
            if used + step > budget:
                break
            start -= 1
            used += step
        end = stop
        while end < len(full):
            step = char_width(full[end])
            if used + step > budget:
                break
            end += 1
            used += step

    row = full[start:end]
    col = string_width(full[start:stop])
    return BoxView([row], (0, min(col, budget)), None, border=False)


class BoxView:
    """The rendered box: rows plus where the cursor belongs.

    ``rows`` is the full picture (the opening rule, the content rows,
    the closing rule, then the hint row when there is one). ``cursor``
    is a ``(row, column)`` into the CONTENT rows — row 0 is the first
    row under the opening rule — with the column in screen columns
    measured from the start of that row. The caller renders the rows and
    then places the cursor; it never re-derives the layout.

    ``border`` says whether a row carries a frame the caller must step
    over before placing the cursor. Both forms drawn today report False
    — the rules are their own rows, and the content rows start at column
    zero — and it is kept because the question is the caller's to ask,
    not to assume.
    """

    __slots__ = ("rows", "cursor", "hint_row", "border")

    def __init__(self, rows: list[str], cursor: tuple[int, int],
                 hint_row: int | None, border: bool = True):
        self.rows = rows
        self.cursor = cursor
        #: Index into ``rows`` of the hint, or None when there is none.
        self.hint_row = hint_row
        self.border = border

    def __eq__(self, other) -> bool:
        if not isinstance(other, BoxView):
            return NotImplemented
        return (self.rows == other.rows
                and self.cursor == other.cursor
                and self.hint_row == other.hint_row)

    def __repr__(self) -> str:  # pragma: no cover — debugging aid
        return (f"BoxView(rows={self.rows!r}, cursor={self.cursor!r},"
                f" hint_row={self.hint_row!r})")


def _cursor_position(wrapped: list[tuple[str, int]], stop: int, inner: int
                     ) -> tuple[int, int]:
    """Where character offset *stop* in the stream landed on screen.

    The cursor is located by WALKING the stream, not by summing row
    widths: a newline consumes a character but no column, and a wide
    character that moved down a row costs its columns where the layout
    put them. Summing widths desynchronises the cursor from its row in
    exactly those two cases.

    Returns (content-row index, column from the row's start). ``stop``
    past the end lands on the very end of the last row, clamped by the
    caller beforehand.
    """
    row = 0
    i = 0
    for row_text, _w in wrapped:
        col = 0
        ended_by_newline = row_text.endswith("\n")
        for ch in row_text:
            if i == stop:
                return (row, col)
            i += 1
            col += char_width(ch)
        # End of the row. A stop here sits after the row's last
        # character. A row ended by a NEWLINE has already consumed that
        # newline in the loop above, so the cursor belongs at the start
        # of the next row; a row ended by a WRAP EDGE consumed nothing
        # extra, so the cursor stays at this row's end. Advancing the
        # character counter here is the off-by-one that put every
        # end-of-row cursor one column and one row wrong.
        if i == stop:
            if ended_by_newline:
                return (row + 1, 0)
            return (row, col)
        row += 1
    return (row - 1 if row else 0, inner)


def render_box(text: str, cursor: int, width: int, hint: str = "",
               status: str = "") -> BoxView:
    """Render the framed input area. Pure; raises nothing.

    ``cursor`` is an offset into ``text`` (0..len). An offset outside
    the text is clamped, because a caller that moved the cursor and the
    buffer in the wrong order must not get a crash for it. A newline in
    ``text`` (a pasted block) hard-wraps the row, matching what the
    terminal would do with it anyway.
    """
    width = int(width or 0)
    text = text or ""
    cursor = max(0, min(int(cursor), len(text)))
    visible = _visible_text(text)
    if width < MIN_WIDTH:
        return _narrow_row(visible, cursor, width)

    # Two rules, not a frame. A closed box has to defend its right edge
    # on every keystroke and at every width, and it buys nothing the
    # rules do not: the eye reads "this is where you type" from the line
    # above and the line below. It also gives the text four more columns
    # and removes the whole class of faults that lived in the borders.
    inner = max(1, width - 1)
    rule = _HORIZONTAL * inner
    # The prompt rides at the start of the wrapped stream; it is part
    # of row 0 and the cursor offset is measured through it.
    stream = PROMPT + visible
    wrapped = _wrap(stream, inner)

    rows: list[str] = [rule]
    cursor_pos = _cursor_position(wrapped, len(PROMPT) + cursor, inner)
    for row_text, _row_w in wrapped:
        rows.append(row_text.removesuffix("\n"))
    rows.append(rule)

    hint_row = None
    # What is running, under the line you type on: a background suite, a
    # sub-agent still out. It was reachable only by asking (`/bash`), so
    # a run started twenty minutes ago was remembered or it was not.
    if status:
        rows.append("  " + _truncate_status(status, inner - 2))
    if hint:
        rows.append("  " + _truncate_hint(hint, inner - 2))
        hint_row = len(rows) - 1
    # border=False: there is no left border to step over, so the caller
    # places the cursor at the column this module reports.
    return BoxView(rows, cursor_pos, hint_row, border=False)


def viewport(view: BoxView, height: int) -> BoxView:
    """A fixed-height window over *view*, keeping the cursor visible.

    The renderer is honest about how tall the text is; the terminal is
    not; an unbounded paste into a 24-row window would push the frame
    off the top of the screen and the caller into cursor-position
    queries. So the wiring shows a window: the top and bottom borders
    and the hint always, and at most *height* content rows around the
    one the cursor is on — with ``…`` markers naming which side was
    cut, because a box that silently eats the top of a pasted block
    looks exactly like lost text.

    Pure, like everything else here. ``height`` below 1 shows one row.
    """
    height = max(1, int(height or 1))
    # Content is what sits BETWEEN the borders. Taking rows[1:] when
    # there is no hint row swept the bottom border in as content, and
    # the same border was appended again below — two "╰────╯" rows.
    content = (view.rows[1:view.hint_row - 1] if view.hint_row is not None
               else view.rows[1:-1])
    if len(content) <= height:
        return view
    cur_row = min(view.cursor[0], len(content) - 1)
    top = min(max(0, cur_row - height // 2),
              len(content) - height)
    shown = content[top:top + height]
    rows = [view.rows[0]]
    if top > 0:
        rows.append("…")
    rows.extend(shown)
    if top + height < len(content):
        rows.append("…")
    if view.hint_row is not None:
        rows.append(view.rows[view.hint_row - 1])   # bottom border
        rows.append(view.rows[view.hint_row])
    else:
        rows.append(view.rows[-1])
    # A marker row above `shown` is a row of the picture, so every
    # content row below it moves down by one. Without this the cursor
    # sat one row high the moment the window scrolled — on a full line
    # of text rather than at its end.
    cursor = (cur_row - top + (1 if top > 0 else 0), view.cursor[1])
    # The hint sits at the end of the NEW picture, which is shorter than
    # the one it was windowed from. Carrying the old index forward made
    # view.rows[view.hint_row] point past the end -- unread today, and
    # wrong the moment anything asks.
    new_hint = (len(rows) - 1) if view.hint_row is not None else None
    return BoxView(rows, cursor, new_hint, border=view.border)
