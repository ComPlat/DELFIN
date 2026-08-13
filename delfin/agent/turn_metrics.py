"""Per-turn timing metrics for the DELFIN agent.

Records one entry per agent turn — provider/model/role, total wall time,
**time-to-first-token (ttft)**, output size and tool-call count — to
``~/.delfin/turn_metrics/<session>.jsonl``.

The point is to make transient backend stalls *visible after the fact*. When a
trivial turn takes a minute and a half (e.g. the field report where Azure
GPT-5.4 answered "Hallo" in 92.7s), the breakdown tells whether the time went
into **waiting for the first token** (a backend/queue stall — high ttft, tiny
output, no tools) versus heavy generation or many tool rounds. That distinction
is exactly what a single "turn took 92.7s" number cannot give.

The severest form of that stall is the turn where the first token never
arrives at all, which records ``ttft_ms=None``. It was excluded from
:func:`is_stall` by a ``ttft is not None`` guard, so the summaries and the
eval report's "turn health" line called a backend that had stopped answering
healthy — and, since such turns count as turns while contributing no ttft,
the reported avg/p90 got *better* the more of them piled up. See
:func:`is_never_started` for how a silent turn is told apart from a record
that simply predates the field.

Three fields exist so that "the backend produced nothing" cannot absorb
every other way a turn can produce nothing:

``error``
    The exception the turn died on. The parameter existed from the start
    and the engine never passed one, so a turn that raised before any
    event — a connect/read timeout, an auth rejection — wrote the same
    record as a live backend that stayed silent, and the log is the only
    evidence there is afterwards. :func:`is_crashed` reads it and
    :func:`is_never_started` steps aside for it.
``first_event_ms``
    Time to the first stream EVENT of any kind, not to the first token.
    ``message_start`` and thinking deltas never stamp ttft, so a live
    connection to a silent model and a transport that delivered nothing
    at all both recorded ``ttft_ms=None``. ``first_event_ms`` separates
    them; :func:`silence_kind` names which one happened.
``retries``
    How many times the client re-issued the request. Three transient
    retries cost 1.5 + 3 + 6 s of sleep plus request latency; without the
    count, a turn that fought for that long looks exactly like a quiet
    one. ``None`` means the backend does not report retries — it is not
    a count of zero.

Provenance, for every field here: a MISSING key is a record written before
the field existed and must never be judged; present-and-``None`` is the
marker of an instrumented turn that has something to say.

Best-effort and dependency-free: never raises, caps file size, no-ops on IO
error. Mirrors :mod:`tool_trace`.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

_DIR = Path.home() / ".delfin" / "turn_metrics"
_MAX_BYTES = 2 * 1024 * 1024    # trim the file when it grows past this
_KEEP_TAIL = 2000               # lines kept when trimming

# A turn that took this long to its FIRST token, while emitting little and
# calling no tools, is flagged as a likely backend stall in summaries.
_SLOW_TTFT_MS = 20_000


def _clean_error(text: str) -> str:
    """An exception message, with credential material taken out.

    ``error`` is the one free-form field here, and an API failure quotes
    what it was given: a URL with a token in the query, an auth header,
    the key it rejected. These files are 0600 but a bug report bundles
    them, so the string is scrubbed before it reaches the disk rather
    than after. Best-effort: a guard that can break the recorder would
    cost more than it saves.
    """
    s = str(text or "")[:300]
    if not s:
        return ""
    try:
        from .output_guard import _redact_secrets
        return _redact_secrets(s, [])
    except Exception:
        return s


def _safe(session: str) -> str:
    s = "".join(c if (c.isalnum() or c in "-_.") else "_"
                for c in (session or "default"))
    return s[:80] or "default"


def metrics_path(session: str) -> Path:
    return _DIR / f"{_safe(session)}.jsonl"


def record(
    session: str,
    *,
    provider: str = "",
    model: str = "",
    role: str = "",
    total_ms: int = 0,
    ttft_ms: int | None = None,
    first_event_ms: int | None = None,
    output_chars: int = 0,
    tool_calls: int = 0,
    stopped: bool = False,
    error: str = "",
    retries: int | None = None,
    input_tokens: int = 0,
    output_tokens: int = 0,
    cached_tokens: int = 0,
) -> None:
    """Append one turn-timing entry to the session log. Never raises.

    The token counts are here because a claim about them could not be
    checked. The prompt is deliberately byte-stable and re-sent every turn,
    justified by prefix caching making that nearly free -- and a probe once
    found ``cached_tokens = 0`` on the KIT endpoint, which would make the
    justification wrong and real shortening the only lever that helps.

    Neither could be confirmed later, because the number reached no
    persistent store: the only consumer renders a cache line ONLY when the
    count is above zero, so on an endpoint that reports none it silently
    vanished. Recording it makes the question answerable from data instead
    of from one probe somebody remembers running.

    ``first_event_ms`` and ``retries`` are written on EVERY entry, ``None``
    included: an omitted key means "this recorder did not have the field",
    which is a different statement from "the field had nothing to report",
    and only the second may be judged.
    """
    try:
        from .state_paths import ensure_dir, open_append, secure_file
        from .state_paths import write_text as _write_secure
        ensure_dir(_DIR)
        p = metrics_path(session)
        try:
            if p.exists() and p.stat().st_size > _MAX_BYTES:
                lines = p.read_text(encoding="utf-8").splitlines()[-_KEEP_TAIL:]
                _write_secure(p, "\n".join(lines) + "\n")
        except Exception:
            pass
        entry = {
            "ts": time.time(),
            "provider": str(provider or ""),
            "model": str(model or ""),
            "role": str(role or ""),
            "total_ms": int(total_ms),
            "ttft_ms": (int(ttft_ms) if ttft_ms is not None else None),
            "first_event_ms": (int(first_event_ms)
                               if first_event_ms is not None else None),
            "output_chars": int(output_chars),
            "tool_calls": int(tool_calls),
            "stopped": bool(stopped),
            "error": _clean_error(error),
            "retries": (int(retries) if retries is not None else None),
            "input_tokens": int(input_tokens or 0),
            "output_tokens": int(output_tokens or 0),
            "cached_tokens": int(cached_tokens or 0),
        }
        # Created 0600 by ``open_append``: a chmod after the first write
        # leaves the file group-readable for the length of that write.
        with open_append(p) as f:
            f.write(json.dumps(entry, ensure_ascii=False) + "\n")
        secure_file(p)          # tighten a file that already existed
    except Exception:
        pass


def read(session: str, *, last_n: int | None = None) -> list[dict]:
    """Return the turn entries for ``session`` (optionally the last N)."""
    try:
        p = metrics_path(session)
        if not p.exists():
            return []
        out = []
        for line in p.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except Exception:
                continue
        return out[-last_n:] if last_n else out
    except Exception:
        return []


def _int(entry: dict, key: str) -> int:
    """Field as an int, treating absent/None/garbage as 0."""
    try:
        return int(entry.get(key) or 0)
    except (TypeError, ValueError):
        return 0


def is_crashed(entry: dict) -> bool:
    """The turn ended in an exception, recorded by the engine that caught it.

    The counterpart to :func:`is_never_started`: both describe a turn that
    delivered nothing, and only this one had somebody to blame other than
    a silent backend. A record with no ``error`` KEY predates the field
    and is not judged — absence of the key is not absence of a crash.
    """
    try:
        return bool(str(entry.get("error") or "").strip())
    except Exception:
        return False


def _delivered_nothing(entry: dict) -> bool:
    """The turn reached its end without giving the user anything.

    No text, no tool call, no counted output token, and no crash to
    explain it. Says nothing about whether the model was working — a turn
    that reasoned for three minutes and then stopped delivered nothing,
    and that is a fault worth naming even though the backend was busy.
    """
    if is_crashed(entry):
        return False
    return (_int(entry, "output_chars") == 0
            and _int(entry, "tool_calls") == 0
            and _int(entry, "output_tokens") == 0)


def _produced_no_token(entry: dict) -> bool:
    """The record says, and can prove, that the model never began.

    :func:`is_never_started` asks this one: an instrumented ``ttft_ms``
    that stayed ``None``, and the rest of the record agreeing.

    A reasoning delta STAMPS ttft — it is the model producing output —
    so a turn that thought and then said nothing fails this predicate and
    passes :func:`_delivered_nothing`. The two used to share one field
    and therefore one answer, which made a busy model that delivered
    nothing indistinguishable from an endpoint that never spoke.
    """
    if "ttft_ms" not in entry or entry.get("ttft_ms") is not None:
        return False
    return _delivered_nothing(entry)


def silence_kind(entry: dict) -> str:
    """Which kind of silence a token-less turn recorded.

    ``""``
        The turn produced tokens, or crashed — neither is silence.
    ``"transport"``
        ``first_event_ms is None``: not one event of any kind arrived, so
        the connection itself never delivered.
    ``"model"``
        A first event arrived and nothing was delivered — the stream was
        alive and the user got nothing. This is the discriminator that
        time-to-first-TOKEN cannot make, because ``message_start`` never
        stamps it and a reasoning delta stamps it without answering.
    ``"unclassified"``
        The record predates ``first_event_ms``. It is silent, and which
        kind cannot be recovered; guessing would invent evidence.

    Unlike :func:`is_never_started` this asks nothing about duration, and
    nothing about whether the model was working: a turn that reasoned at
    length and then stopped is silence to the person waiting for it.
    """
    try:
        if not _delivered_nothing(entry):
            return ""
        if "first_event_ms" not in entry:
            return "unclassified"
        return "transport" if entry.get("first_event_ms") is None else "model"
    except Exception:
        return ""


def is_never_started(entry: dict) -> bool:
    """A turn that spent the wall clock without a single token arriving.

    Telling this apart from "the field was not recorded" is the whole
    problem, and the record answers it in two independent ways.

    A recorded error settles it first: the turn did not wait out a silent
    backend, it died, and :func:`is_crashed` is the predicate for that.
    Both shapes write ``ttft_ms=None`` with an empty rest of the record,
    so without this the stall count absorbed every crash that took more
    than twenty seconds to fail — and the log, which is the only evidence
    left afterwards, could not tell the two apart at all.

    Provenance: :func:`record` always writes ``ttft_ms``, so the KEY is the
    marker. Present-and-``None`` comes from an instrumented turn whose
    first-token stamp was never taken — the backend stayed silent. A
    MISSING key is a record written before the field existed; it says
    nothing about what arrived, and judging it would rewrite yesterday's
    logs into stalls the moment this rule shipped.

    Corroboration: the rest of the record must agree that the turn really
    was silent — no text, no tool calls, no counted output tokens. If any
    of those is non-zero while ttft is ``None``, tokens demonstrably DID
    arrive and only the stamp was missed (thinking deltas, for one, never
    stamp it). That is a recording gap; calling it a backend stall would
    be a guess dressed as evidence.

    Wait: with no first token, the whole turn IS the wait, so ``total_ms``
    faces the same threshold a recorded ttft would. A turn that produced
    nothing within milliseconds never waited on the backend — that is an
    immediate failure, a different defect, and it must not inflate this
    count.
    """
    try:
        return (_produced_no_token(entry)
                and _int(entry, "total_ms") >= _SLOW_TTFT_MS)
    except Exception:
        return False


def is_stall(entry: dict) -> bool:
    """A turn dominated by waiting for the first token — little output, no
    tools — i.e. the backend, not the agent, ate the time.

    A turn where the first token never arrived at all is the extreme of
    that same failure and counts here too. The original ``ttft is not
    None`` guard excluded it, so the one shape this module was built to
    expose — and the one the dashboard watchdog kills on its first-token
    budget — was the single shape reported as healthy. Who ended such a
    turn does not matter: by the time anyone stopped it, the silence had
    already lasted longer than the threshold.

    A turn that raised is NOT counted here — it is a crash, reported by
    :func:`is_crashed` and counted separately in the roll-up, so neither
    failure hides inside the other's number.
    """
    try:
        if is_never_started(entry):
            return True
        ttft = entry.get("ttft_ms")
        return (ttft is not None and int(ttft) >= _SLOW_TTFT_MS
                and _int(entry, "tool_calls") == 0
                and _int(entry, "output_chars") <= 400)
    except Exception:
        return False


def _percentile(sorted_vals: list[int], q: float) -> int:
    """Nearest-rank percentile of a pre-sorted list (stdlib only)."""
    if not sorted_vals:
        return 0
    k = int(round((q / 100.0) * (len(sorted_vals) - 1)))
    k = max(0, min(len(sorted_vals) - 1, k))
    return sorted_vals[k]


def aggregate_turn_stats(
    window_days: float = 7,
    *,
    dir_path: Path | None = None,
    now: float | None = None,
) -> dict:
    """Windowed roll-up across ALL session files — the telemetry consumer
    the per-turn logs have been missing.

    Returns ``{turns, avg_ttft_ms, p90_ttft_ms, stalls, crashes,
    stopped_count}`` where ``stalls`` counts entries flagged by
    :func:`is_stall` and the ttft numbers cover only turns that recorded a
    first token.  ``crashes`` exists because the stall count stopped
    absorbing turns that raised: a run of crashing turns would otherwise
    have made this roll-up quieter than before, which is the exact failure
    mode the stall counter was fixed for.  Entries
    outside the last ``window_days`` are skipped (``window_days <= 0``
    disables the window); entries without a parseable ``ts`` stay
    visible.  Best-effort: never raises — corrupt lines/files skipped.

    A turn that never produced a token has no ttft to average, so it can
    only ever pull the reported percentiles DOWN by leaving their sample
    while still counting as a turn.  That is why it has to land in
    ``stalls``: the report renders the stall count on the same line as
    avg/p90, so a backend that answers less often now moves a number
    upward there instead of quietly flattering the latency figures.  The
    alternative — folding the never-started turns into the ttft sample at
    their ``total_ms`` — was rejected because it would report a
    first-token time for turns that never had one, turning a censored
    observation into an invented measurement; :func:`is_never_started`
    is public so a consumer can count the two kinds apart.
    """
    empty = {"turns": 0, "avg_ttft_ms": 0, "p90_ttft_ms": 0,
             "stalls": 0, "crashes": 0, "stopped_count": 0,
             "never_started": 0, "ttft_sample": 0}
    try:
        base = Path(dir_path) if dir_path else _DIR
        cutoff: float | None = None
        if window_days and float(window_days) > 0:
            cutoff = ((now if now is not None else time.time())
                      - float(window_days) * 86400.0)
        turns = stalls = crashes = stopped = never_started = 0
        ttfts: list[int] = []
        for fp in sorted(base.glob("*.jsonl")):
            try:
                lines = fp.read_text(encoding="utf-8").splitlines()
            except Exception:
                continue
            for ln in lines:
                try:
                    e = json.loads(ln)
                except Exception:
                    continue
                if not isinstance(e, dict):
                    continue
                if cutoff is not None:
                    try:
                        if float(e["ts"]) < cutoff:
                            continue
                    except (KeyError, TypeError, ValueError):
                        pass                    # no/bad ts → keep visible
                turns += 1
                ttft = e.get("ttft_ms")
                if ttft is not None:
                    try:
                        ttfts.append(int(ttft))
                    except (TypeError, ValueError):
                        pass
                if is_never_started(e):
                    never_started += 1
                if is_stall(e):
                    stalls += 1
                if is_crashed(e):
                    crashes += 1
                if e.get("stopped"):
                    stopped += 1
        ttfts.sort()
        return {
            "turns": turns,
            "avg_ttft_ms": (int(round(sum(ttfts) / len(ttfts)))
                            if ttfts else 0),
            "p90_ttft_ms": _percentile(ttfts, 90),
            "stalls": stalls,
            "crashes": crashes,
            "stopped_count": stopped,
            # A turn that waited out a silent backend. Counted apart from
            # a crash, which HAS an error, and apart from a stall, which
            # produced something late. It was only reachable through the
            # public predicate before, so the report could not say it --
            # and what the log cannot express, the report cannot say.
            "never_started": never_started,
            # How many turns the two ttft figures above are actually an
            # average OF. A mean over three samples out of ninety turns
            # reads exactly like a mean over ninety, and the difference
            # is whether the number means anything.
            "ttft_sample": len(ttfts),
        }
    except Exception:
        return empty


def format_summary(entries: list[dict], *, limit: int = 30) -> str:
    """One line per turn; flags stalls and crashes. Empty when no entries.

    This is what a bug report carries, so a line has to say which failure
    it was: a backend that answered nothing, a backend that was never
    reached, or a turn that died — with the retry count when the client
    reports one, because three retries and a quiet wait are the same
    number of seconds and not the same problem.
    """
    if not entries:
        return ""
    rows = []
    for e in entries[-limit:]:
        total = int(e.get("total_ms") or 0)
        ttft = e.get("ttft_ms")
        ttft_s = f"{int(ttft)/1000:.1f}s" if ttft is not None else "—"
        if is_crashed(e):
            flag = f"  ✖ crashed: {str(e.get('error') or '')[:80]}"
        elif is_stall(e):
            kind = {"transport": " (no transport)",
                    "model": " (model silent)"}.get(silence_kind(e), "")
            flag = f"  ⚠ backend-stall{kind}"
        else:
            flag = ""
        retries = e.get("retries")
        if retries:
            flag += f"  retries={int(retries)}"
        rows.append(
            f"{e.get('model','?')}  total={total/1000:.1f}s  "
            f"ttft={ttft_s}  out={int(e.get('output_chars') or 0)}ch  "
            f"tools={int(e.get('tool_calls') or 0)}{flag}"
        )
    return "\n".join(rows)


__all__ = ["metrics_path", "record", "read", "is_stall", "is_never_started",
           "is_crashed", "silence_kind", "format_summary",
           "aggregate_turn_stats"]
