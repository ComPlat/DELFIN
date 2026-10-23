"""Deterministic one-line digests ("wanted posters") for tool results.

When api_client._elide_old_tool_results drops an old tool result to free
context, the content is gone. Bug 172455: an agent that no longer knows
WHICH file it read or what a grep found re-runs the tool — 149 reads of
one file in a single run. This module builds a short, deterministic
replacement: tool name, the arguments that matter, and the few facts an
agent needs to decide whether re-running is worth it.

Rules (deliberate):
- No model calls. Same input -> same digest, byte for byte. Not cached:
  each result is elided once, and a cache keyed on whole tool outputs
  would hold megabytes of them.
- Never any secret: everything passes through output_guard.scrub_secrets
  before being embedded, and arguments go through it too.
- The digests encode the REAL output shapes of the api_client tools
  (numbered lines for read_file, "path:line: text" for grep_file, a JSON
  payload with exit_code/stdout/stderr for bash), not generic guesses.
- If anything surprises us, degrade to the generic digest rather than
  raising — an elision path must never crash the loop.
"""

from __future__ import annotations

import json
import re

from .output_guard import scrub_secrets

_MAX_CHARS = 300
_MAX_DEFS = 3
_MAX_HITS = 3
_MAX_FAILED = 3

# read_file output: "N  text" with two spaces after the number.
_READ_LINE_RE = re.compile(r"^(\d+)  ?(.*)$")
# grep_file output: "path:line: text" (see api_client._GREP_HIT_RE).
_HIT_RE = re.compile(r"^([^\s:][^:]*):(\d+):")
# read_file tail marker: "... (N lines total, showing A-B)"
_TOTALS_RE = re.compile(r"\((\d+) lines total, showing (\d+)-(\d+)\)")
# def/class at start of a line (after the read_file line number).
_DEF_RE = re.compile(r"^(?:def|class)\s+\w+")
# pytest verdict line, e.g. "3 failed, 12 passed in 4.5s".
_VERDICT_RE = re.compile(
    r"^\d+ (?:failed|passed|error)s?(?:, \d+ (?:failed|passed|error)s?)*"
    r"(?:, \d+ warning?s?)? in [\d.]+s?\b[^\n]*", re.MULTILINE)
_FAILED_RE = re.compile(r"^FAILED (\S+)", re.MULTILINE)

_TEST_TOOLS = {"run_tests", "gate"}
_GREP_TOOLS = {"grep_file", "grep", "search_files"}


def _clip(text: str, limit: int) -> str:
    text = text.strip()
    return text if len(text) <= limit else text[:limit - 1] + "…"


def _clean(value) -> str:
    """Scrub a value and make it safe to embed on one line."""
    return scrub_secrets(str(value)).replace("\n", " ").strip()


def _read_digest(arguments: dict, content: str) -> str:
    path = _clean(arguments.get("path", "?"))
    shown_start = shown_last = None
    total = None
    # The slice actually shown comes from the arguments when present, and
    # from the tail marker otherwise (offset/limit are optional; the
    # marker is what read_file always appends when it cut).
    offset = arguments.get("offset")
    limit = arguments.get("limit")
    m_total = _TOTALS_RE.search(content)
    lines = content.split("\n")
    if m_total:
        total = int(m_total.group(1))
        shown_start, shown_last = int(m_total.group(2)), int(m_total.group(3))
    else:
        n = len([ln for ln in lines if _READ_LINE_RE.match(ln)])
        if isinstance(offset, int):
            shown_start = max(1, offset)
            shown_last = shown_start + n - 1 if n else shown_start
        else:
            shown_start, shown_last = (1, n) if n else (None, None)
        total = shown_last if shown_last else n
    defs = []
    for ln in lines:
        m = _READ_LINE_RE.match(ln)
        if m and _DEF_RE.match(m.group(2) or ""):
            defs.append(_clip(m.group(2), 60))
            if len(defs) >= _MAX_DEFS:
                break
    parts = [f"read_file {path}"]
    if shown_start and shown_last:
        parts.append(f"lines {shown_start}-{shown_last}")
        if total and total > shown_last:
            parts.append(f"of {total} lines total")
    if defs:
        parts.append("defs: " + "; ".join(defs))
    return ", ".join(parts)


def _grep_digest(arguments: dict, content: str) -> str:
    pattern = _clip(_clean(arguments.get("pattern", "?")), 60)
    hits = []
    for ln in content.split("\n"):
        m = _HIT_RE.match(ln)
        if m:
            hits.append(f"{m.group(1)}:{m.group(2)}")
            if len(hits) >= _MAX_HITS:
                break
    if not hits and not content.strip():
        body = "0 hits"
    elif not hits:
        # "No matches found." or an unexpected shape: keep the first line
        # as the finding, scrubbed.
        body = _clip(scrub_secrets(content), 120)
    else:
        n_hits = sum(1 for ln in content.split("\n") if _HIT_RE.match(ln))
        body = f"{n_hits} hits at " + ", ".join(hits)
    return f"grep {pattern}: {body}"


def _payload_of(content: str):
    text = content.strip()
    if text.startswith("{") or text.startswith("["):
        try:
            return json.loads(text)
        except (json.JSONDecodeError, ValueError):
            return None
    return None


def _bash_digest(arguments: dict, content: str) -> str:
    cmd = _clip(_clean(arguments.get("command", "")), 80)
    payload = _payload_of(content)
    if isinstance(payload, dict):
        exit_code = payload.get("exit_code")
        out = payload.get("stdout") or payload.get("output") or ""
        err = payload.get("stderr") or ""
        tail = _clip(scrub_secrets(out.strip() or err.strip())
                     .split("\n")[-1] if (out or err) else "", 100)
        verdict = _VERDICT_RE.search(str(out) + "\n" + str(err))
        failed = [f.group(1) for f in _FAILED_RE.finditer(str(out))]
        parts = [f"bash {cmd}".rstrip()]
        if exit_code is not None:
            parts.append(f"exit {exit_code}")
        if verdict:
            parts.append(_clip(scrub_secrets(verdict.group(0)), 80))
        if failed:
            parts.append("FAILED " + ", ".join(
                _clip(f, 60) for f in failed[:_MAX_FAILED]))
        if tail:
            parts.append(f"last: {tail}")
        return ", ".join(parts)
    # Plain text (an older/foreign shape): keep the last line, which is
    # where a summary or traceback ends up.
    last = _clip(scrub_secrets(content).split("\n")[-1], 100)
    return ", ".join(p for p in (f"bash {cmd}".rstrip(),
                                 f"last: {last}") if p)


def _generic_digest(tool_name: str, arguments: dict, content: str) -> str:
    first = _clip(scrub_secrets(content).split("\n")[0] if content else "", 100)
    args = ""
    if isinstance(arguments, dict):
        for key in ("path", "query", "pattern", "target", "command", "doc_id"):
            if arguments.get(key):
                args = _clip(_clean(arguments[key]), 60)
                break
    parts = [tool_name]
    if args:
        parts.append(args)
    parts.append(f"{len(content)} chars")
    if first:
        parts.append(f"first: {first}")
    return ", ".join(parts)


def _build(tool_name: str, arguments, content: str) -> str:
    args = arguments if isinstance(arguments, dict) else {}
    if tool_name in _GREP_TOOLS:
        body = _grep_digest(args, content)
    elif tool_name in ("read_file", "read_text_file", "read_section"):
        body = _read_digest(args, content)
    elif tool_name == "bash" or tool_name in _TEST_TOOLS:
        body = _bash_digest(args, content)
    else:
        body = _generic_digest(tool_name, args, content)
    body = scrub_secrets(body)
    return body if len(body) <= _MAX_CHARS else body[:_MAX_CHARS - 1] + "…"


def digest(tool_name: str, arguments: dict, content: str) -> str:
    """A deterministic, secret-free, <=300-char wanted poster for one tool
    result. Same inputs -> same string. Never raises."""
    try:
        return _build(tool_name, arguments, content or "")
    except Exception:
        # The elision path must never crash on a surprising shape.
        try:
            return scrub_secrets(
                f"{tool_name}, {len(content or '')} chars (digest failed)")
        except Exception:
            return f"{tool_name}: digest failed"
