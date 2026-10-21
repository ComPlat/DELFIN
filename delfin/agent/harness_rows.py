"""User-role rows the harness wrote, marked where they are written.

A chat endpoint has two speaking roles, so everything the harness tells
the model between the user's messages -- turn-end follow-ups, verification
turns, job and run notes, wake-ups, other sessions' messages, compaction
summaries -- travels as role "user". The reader of the user's own words
(``api_client._classifier_user_texts``; through it the action classifier,
skill learning and the refusal lift) leaves a marked row out whole.

One mark, three carriers, chosen by where the row lives:

* ``HarnessRow`` -- a ``dict`` subclass for a row on its way to the
  endpoint (the tool loop's request list, the engine's wire view). No
  extra key: ``json.dumps``, the SDK's transform and the prompt bytes are
  those of a plain dict, so neither a strict endpoint nor the prefix cache
  sees it.
* ``HARNESS_KEY`` (``"_harness": True``) -- on a row kept in a history that
  is saved as JSON (``engine.messages``, the CLI's display projection). It
  survives ``json``, ``dict(row)`` and ``copy.deepcopy``. Private keys
  never reach a request: ``engine._wire_messages`` strips every ``_`` key
  (and renders a marked row as a ``HarnessRow``), the OpenAI-compatible
  client copies only role, content and the tool keys.
* ``HarnessText`` -- a ``str`` subclass for text that is not a row yet: a
  turn's message handed to ``stream_response`` and a queued steer.

A dashboard transcript row with ``origin == "event"`` (written by the
dashboard for a turn nobody typed) is the same mark and is read as one.
Rows restored from a session saved before the mark carry none; the
reader's text rules stay as their fallback.

Rejected: a side table of row ids or indices (indices shift under
compaction, the alternation sanitiser and trimming, and an id is a new key
on the row anyway); the key alone (it would reach the endpoint from the
tool loop's request list, which nothing strips); the subclass alone (a
JSON round trip and ``dict(row)`` drop the type, so a restored history
would have lost every mark).
"""

from __future__ import annotations

from typing import Any

HARNESS_KEY = "_harness"


class HarnessRow(dict):
    """A user-role row the harness wrote, as it is sent: a plain ``dict``
    to every consumer, identified by its type alone."""


class HarnessText(str):
    """Text the harness wrote that becomes a user-role row later.

    ``str`` methods return plain ``str`` (``.strip()``, concatenation,
    f-strings), so a caller that rebuilds the text re-wraps it."""


def harness_row(content: Any, role: str = "user") -> HarnessRow:
    """A row for the tool loop's request list; ``HarnessText`` content is
    stored as plain ``str``."""
    if isinstance(content, str):
        content = str.__str__(content)
    return HarnessRow(role=role, content=content)


def user_row(text: Any) -> dict:
    """The row a queued steer becomes: marked when the text is
    ``HarnessText``, a plain ``dict`` when somebody typed it."""
    if isinstance(text, HarnessText):
        return harness_row(text)
    return {"role": "user", "content": text}


def is_harness_text(value: Any) -> bool:
    return isinstance(value, HarnessText)


def is_harness_row(msg: Any) -> bool:
    """True for a row carrying any of the three marks, or whose content
    is ``HarnessText``."""
    if isinstance(msg, HarnessRow):
        return True
    if not isinstance(msg, dict):
        return False
    return (msg.get(HARNESS_KEY) is True or msg.get("origin") == "event"
            or isinstance(msg.get("content"), HarnessText))


def history_mark(row: dict, source: Any) -> dict:
    """``row`` with the history key set when ``source`` is a marked row;
    for a projection that copies the keys it knows."""
    if is_harness_row(source):
        row[HARNESS_KEY] = True
    return row
