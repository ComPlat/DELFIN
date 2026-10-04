"""Isolate text that enters the model context from outside.

Anything whose author is not the user of this session — another session's
message, a tool's output, a fetched page, a wake-up — is data, not an
instruction. The model has no marker telling it which text is data, so this
module gives it one.

``wrap(source, text)`` puts foreign text in a labelled block whose fence
carries a random per-call nonce. The opener and the closer share that nonce,
so the content cannot close the fence early or forge a second one: any
fence-looking line it emits has a nonce it could not know.

``unwrap(wrapped)`` recovers the inner text against the fence's own nonce —
used by tests to prove a content that contains fence-looking lines still
cannot escape its block.

``flags(text)`` classes instruction-like phrases the text carries. It matches
worded instruction shapes from an ALLOW-list per kind (approval claims,
"ignore ... instructions", git-push asks, destructive asks, permission asks),
never a deny-list of bare tokens, so "push onto the stack" or "the report
shows 3 deletions" are plain text, not instructions.
"""

from __future__ import annotations

import re
import secrets

_OPEN_RE = re.compile(
    r"<<<untrusted\s+source=\"[^\"]*\"\s+nonce=\"([^\"]*)\"\s*>>>")
_CLOSE_RE = re.compile(
    r"<<<end untrusted\s+source=\"[^\"]*\"\s+nonce=\"([^\"]*)\"\s*>>>")


def _nonce() -> str:
    return secrets.token_hex(8)


def wrap(source: str, text: str, *, nonce: str | None = None) -> str:
    """Fence ``text`` in a labelled block that names ``source``.

    The fence delimiter carries a fresh random nonce (unless one is given,
    for deterministic tests). Opener and closer share it, so the content
    cannot close the block early or emit text that the harness would read
    as its own.
    """
    nonce = nonce or _nonce()
    src = str(source or "")
    body = "" if text is None else str(text)
    return (f"<<<untrusted source=\"{src}\" nonce=\"{nonce}\">>>\n"
            f"{body}\n"
            f"<<<end untrusted source=\"{src}\" nonce=\"{nonce}\">>>")


def unwrap(wrapped: str) -> str:
    """Recover the text inside a well-formed wrapped block.

    The closer must carry the SAME nonce as the opener; the outermost (last)
    matching closer delimits the block, so a content that contains a
    fence-looking line cannot shorten the block and leak or hijack. Raises
    ValueError when no opener with a matching closer is found.
    """
    opener = _OPEN_RE.search(wrapped or "")
    if not opener:
        raise ValueError("no untrusted opener found")
    nonce = opener.group(1)
    start = opener.end()
    after = wrapped[start:]
    closers = list(_CLOSE_RE.finditer(after))
    end = None
    for m in closers:
        if m.group(1) == nonce:
            end = m.start()  # last matching closer wins (outermost)
    if end is None:
        raise ValueError("no matching nonced closer found")
    inner = after[:end]
    if inner.endswith("\n"):
        inner = inner[:-1]
    if inner.startswith("\n"):
        inner = inner[1:]
    return inner


# --------------------------------------------------------------------------
# flags(): instruction-like phrases, classified from allow-lists per kind.
# --------------------------------------------------------------------------

# Each entry is (kind, regex). The regex matches a WORDED instruction shape,
# not a bare token, so a word used descriptively ("push onto the stack",
# "3 deletions") is not an instruction. All patterns are case-insensitive.

_ALLOWED = (
    # approval claims: whoever wrote the text asserting that something is
    # approved, sanctioned, or already granted.
    ("approval", re.compile(r"\bi\s+approve\b|\b(we|the)\s+approve\b"
                            r"|\bapproved\b|\bapproval\s+granted\b"
                            r"|\bgrant(ed)?\s+approval\b|\bqs\s+approved\b"
                            r"|\boperator\s+approved\b", re.I)),
    # "ignore / disregard previous instructions (and ...)".
    ("ignore", re.compile(r"\bignore\s+(all|any|the|your)?\s*(previous|prior|"
                          r"earlier|given|current)\s+(instructions?|rules?)"
                          r"|\bignore\s+(all|any|the|your)"
                          r"|\bdisregard\s+(all|any|previous)"
                          r"|\bdo\s+not\s+follow\s+the\s+rules\b", re.I)),
    # git-push asks.
    ("push", re.compile(r"\bgit\s+push\b|\brun\s+git\s+push\b"
                        r"|\bpush\s+--force\b|\bpush\s+-f\b"
                        r"|\bforce\s*-?\s*push\b", re.I)),
    # destructive (deleting / overwriting / history-rewriting) asks.
    ("destructive", re.compile(r"\brm\s+-rf\b|\brm\s+-r\b|\brm\s+-f\b"
                               r"|\bgit\s+reset\s+--hard\b|\bgit\s+clean\s+-f"
                               r"|\bgit\s+rm\b|\bdrop\s+table\b|"
                               r"\bdrop\s+database\b", re.I)),
    # permission / authority requests aimed at the receiver.
    ("permission", re.compile(r"\bgive\s+me\s+permission\b|\bgrant\s+me\b"
                              r"|\ballow\s+me\s+to\b|\byou\s+(are|'re)\s+"
                              r"(permitted|allowed|authorized)\b"
                              r"|\bpermission\s+to\b|\bgive\s+me\s+access\b"
                              r"|\bi\s+need\s+permission\b", re.I)),
)


def flags(text: str | None) -> list[str]:
    """Classify instruction-like phrases in ``text``.

    Returns one label per kind found, e.g. ``["approval: approved",
    "push: git push"]`` — the kind first, so a caller can gate on the kind
    without parsing. Empty list for no instructions. Never raises; None and
    empty text return [].
    """
    if not text:
        return []
    found: list[str] = []
    for kind, pattern in _ALLOWED:
        m = pattern.search(text)
        if m:
            found.append(f"{kind}: {m.group(0).strip()}")
    return found
