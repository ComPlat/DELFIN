"""Per-skill usage ledger (wave 15, package L1).

A live skill is counted every time it is used, with the outcome, so a
skill that keeps failing is shown with its record instead of silently
staying in the list. This is the read side of "close the skill loop":
the proposal panel and the CLI skills view render the ledger, and the
skill-invocation path (api_client, protected) records into it on each
call.

Memory-poisoning guard: this ledger ONLY counts uses of skills. It
never reads session text, tool output or any content that could carry a
rule; `record_use` takes a skill name and an outcome, nothing else. A
name is validated as a plain identifier (no path separators, no '..',
no newline, no leading dot) -- the same guard skill_proposals uses for
proposal names, so an injected string can never escape the key slot or
carry free text. Stored under ``~/.delfin/skills/_usage``, invisible to
``discover_skills``: the usage ledger never activates anything.

Outcomes are strict: "success" or "failure". A use with any other
outcome raises ValueError (a contract bug must show in tests, never be
silently recorded). IO is best-effort: a write that cannot land must not
crash the invocation path it is wired into, so every read/write is
guarded and returns the data it has.
"""
from __future__ import annotations

import json
from pathlib import Path

from .memory_store import _atomic_write, _private_dir

_OK_OUTCOMES = frozenset({"success", "failure"})


def _usage_dir() -> Path:
    # Resolved per call, like every other path in memory_store: tests
    # (and any redirected HOME) must be able to move the store after
    # this module is loaded.
    return Path.home() / ".delfin" / "skills" / "_usage"


def __getattr__(name):
    # PEP 562: ``USAGE_DIR`` reads the CURRENT home, like the proposals
    # store. A module constant would pin the first HOME it ever saw.
    if name == "USAGE_DIR":
        return _usage_dir()
    raise AttributeError(name)


def _ledger_path(name: str) -> Path:
    return _usage_dir() / f"{name}.jsonl"


def _valid_name(name: str) -> str:
    """A pure identifier for a skill's usage key; raises ValueError."""
    n = (name or "").strip()
    if not n or "/" in n or n.startswith(".") or "\n" in n or "\r" in n:
        raise ValueError(f"invalid skill usage key: {name!r}")
    return n


# --------------------------------------------------------------------------
# record / read
# --------------------------------------------------------------------------


def record_use(name: str, outcome: str, *, detail: str = "",
               when: str = "") -> None:
    """Append one use of a skill to its ledger. Best-effort, never raises
    on IO.

    ``outcome`` must be "success" or "failure"; anything else raises
    ValueError so a wiring bug surfaces in tests instead of recording a
    meaningless row. ``when`` is an ISO timestamp (utcnow if omitted).
    The record carries no free advisory text beyond ``detail`` (a short
    context string for the human); it is display-only and can never turn
    into a rule.
    """
    import datetime as _dt
    n = _valid_name(name)
    if outcome not in _OK_OUTCOMES:
        raise ValueError(
            f"invalid skill outcome {outcome!r} for {name!r}: "
            f"expected one of {sorted(_OK_OUTCOMES)}")
    at = when or _dt.datetime.now(_dt.timezone.utc).isoformat(
        timespec="seconds")
    line = json.dumps({
        "skill": n,
        "outcome": outcome,
        "detail": str(detail or ""),
        "when": at,
    }, ensure_ascii=False)
    try:
        _private_dir(_usage_dir())
    except Exception:
        return
    # The append must be atomic against concurrent writers (a dashboard
    # session and a CLI session can share one repo): read the ledger,
    # append the line, and write it back with _atomic_write (temp file +
    # os.replace). Never a plain open(..., "a"), which truncates-then-
    # writes and lets a racing reader see a torn file.
    try:
        existing = ""
        try:
            existing = _ledger_path(n).read_text(encoding="utf-8")
        except OSError:
            pass
        _atomic_write(_ledger_path(n), existing + line + "\n")
    except Exception:
        return


def usage_records(name: str) -> list[dict]:
    """One dict per recorded use, newest first: skill, outcome, detail,
    when. Any unreadable line is skipped (never raises)."""
    n = _valid_name(name)
    path = _ledger_path(n)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    out = []
    for ln in lines:
        ln = ln.strip()
        if not ln:
            continue
        try:
            rec = json.loads(ln)
        except (ValueError, TypeError):
            continue
        if not isinstance(rec, dict) or rec.get("outcome") not in _OK_OUTCOMES:
            continue
        out.append({
            "skill": n,
            "outcome": rec.get("outcome"),
            "detail": str(rec.get("detail", "")),
            "when": str(rec.get("when", "")),
        })
    out.reverse()
    return out


def usage_for(name: str) -> dict:
    """One skill's roll-up: uses, successes, failures, last_used,
    last_outcome. Empty record when the skill was never used."""
    rows = usage_records(name)
    failures = sum(1 for r in rows if r["outcome"] == "failure")
    return {
        "skill": name,
        "uses": len(rows),
        "successes": len(rows) - failures,
        "failures": failures,
        "last_used": rows[0]["when"] if rows else "",
        "last_outcome": rows[0]["outcome"] if rows else "",
    }


def usage_summary() -> dict:
    """Roll-up over every skill that has a ledger entry, keyed by name."""
    base = _usage_dir()
    out = {}
    try:
        entries = sorted(base.glob("*.jsonl")) if base.is_dir() else []
    except OSError:
        return out
    for p in entries:
        name = p.stem
        try:
            out[name] = usage_for(name)
        except ValueError:
            continue
    return out
