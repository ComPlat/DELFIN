"""What the members of a wave ran into, counted, so the agents can be tuned
from counts rather than from impressions.

``delfin.agent.friction_report`` reads this machine's audit log and turn
metrics, which carry every refusal but not who in which wave met it, nor
the waits, retries, stalls and hints that only the session's own stream
shows. This reads the members' stream-json logs and their real runs:

  denial      a refused action, by tool and reason code
  stall       the progress watchdog's class (model-wait, gate-queue, ...)
  api_retry   a transient endpoint error, by its kind
  slot_wait   slot_timeout   no_first_byte   waiting for the endpoint
  hint        a correction the tools gave (shell, background, repeat, ...)
  budget      refusal     context     nudge     model     cache    notices
  stop        a run that did not end ``done`` (with its error)
  end         every run's end, with its tool calls (the denominator)
  crash       a traceback in the log; ``warning`` a SyntaxWarning and the like

``report`` groups them by kind and code, and per model as a rate per 100
tool calls. Details pass the bug-report redactor and lose home paths.

``report`` reads every log whole, which is fine for a person asking and not
for a pass every minute. The tick therefore keeps a running ledger:
``sweep`` reads each log from the byte the last sweep stopped at and appends
what is new to ``<root>/.friction.jsonl``, and ``write_report`` groups that
ledger into ``<root>/friction.json`` at most once an hour. Until 2026-11-25
nothing ran any of this unattended, so the counts existed only when someone
thought to ask.
"""
from __future__ import annotations

import collections
import json
import os
import re
import stat
from pathlib import Path

from . import wave

_NOTICES = [
    (r"waiting for a slot on", ("slot_wait", "endpoint_slot")),
    (r"^\s*stall: ", None),                          # repeats the stall event
    (r"auto-continue|Action watcher answers again", None),
    (r"Transient API error \((\w+)\)", ("api_retry", 1)),
    (r"Waited \d+ s for a local slot .*\((\w+)\)", ("slot_timeout", 1)),
    (r"No byte from the endpoint .*\((\w+)\)", ("no_first_byte", 1)),
    (r"\[Shell\]", ("hint", "shell")),
    (r"\[Background command\]", ("hint", "background")),
    (r"\[repeat\]", ("hint", "repeat")),
    (r"\[exists first\]", ("hint", "exists_first")),
    (r"^\s*\[git\] \d+ tracked files? changed and not committed", ("hint", "git_uncommitted")),
    (r"This session's work is in", ("hint", "work_dir")),
    (r"Plan mode is on", ("refusal", "plan_mode")),
    (r"Run budget (reached|exhausted)", ("budget", "run_budget")),
    (r"Action watcher down", ("watcher", "classifier_down")),
    (r"building its cache", ("cache", "first_turn_build")),
    (r"open task\(s\)", ("stop", "open_tasks")),
    (r"Context compacted", ("context", "compacted")),
    (r"announced an action", ("nudge", "announced_action")),
    (r"\[empty turn\]", ("model", "empty_turn")),
]
_STDERR = [
    (r"permission mode \w+ from --permission-mode", None),
    (r"SyntaxWarning", ("warning", "SyntaxWarning")),
    (r"request-head invariants broken", ("warning", "request_head_invariant")),
]


def _redact(text: str) -> str:
    try:
        from delfin.agent.bug_report import redact_for_report
        text = redact_for_report(text)[0]
    except Exception:
        text = re.sub(r"[A-Za-z0-9_\-]{32,}", "<long>", text)
    return re.sub(r"\s+", " ", re.sub(r"/home/[^/\s\"'`]*", "<home>", text)).strip()[:240]


def _classify(table, text):
    for pat, out in table:
        m = re.search(pat, text)
        if m:
            if out is None:
                return None
            kind, code = out
            return kind, (m.group(1) if code == 1 else code)
    return "", ""


def records(lines: list[str]):
    """``(kind, code, detail)`` for every friction in a stream-json log."""
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        if line.startswith("{"):
            try:
                e = json.loads(line)
            except ValueError:
                continue
            t = e.get("type")
            if t == "permission_denied":
                yield "denial", f"{e.get('tool')}:{e.get('reason_code')}", e.get("text") or ""
            elif t == "stall":
                yield "stall", str(e.get("class")), f"{float(e.get('duration_s') or 0):.0f} s"
            elif t == "notice":
                text = e.get("text") or ""
                c = _classify(_NOTICES, text)
                if c is None:
                    continue
                kind, code = c
                if not kind:
                    kind, code = "notice", re.sub(r"\d+", "N", text.strip())[:50]
                yield kind, code, text
            elif t == "result":
                stop = e.get("stop_reason") or "none"
                calls = e.get("tool_calls")
                n = len(calls) if isinstance(calls, list) else int(calls or 0)
                yield "end", stop, f"tool_calls={n}"
                if stop != "done" or e.get("error"):
                    yield "stop", stop, e.get("error") or ""
            continue
        if line.startswith("Traceback (most recent call last)"):
            exc = ""
            while i < len(lines):
                s = lines[i].rstrip()
                i += 1
                if s and not s[0].isspace():
                    exc = s.strip()
                    break
            yield "crash", exc.split(":")[0] or "unknown", exc
            continue
        c = _classify(_STDERR, line)
        if c and c[0]:
            yield c[0], c[1], line


def collect(run_names: list[str] | None = None) -> list[dict]:
    """Every friction record of the runs (default: all), from the members'
    logs and their real runs."""
    out = []
    for run in run_names if run_names is not None else wave.runs():
        for tid, a in wave.load_state(run)["agents"].items():
            model = str((a.get("task") or {}).get("model") or "")
            sources = [(Path(a["log"]), "member")]
            sources += [(p, "realrun:" + p.stem) for p in
                        sorted((Path(a["worktree"]) / "scratch" / "realrun").glob("*.jsonl"))]
            for path, src in sources:
                try:
                    lines = path.read_text(errors="replace").splitlines()
                except OSError:
                    continue
                for kind, code, detail in records(lines):
                    out.append({"run": run, "id": tid, "model": model, "src": src,
                                "kind": kind, "code": _redact(str(code))[:60],
                                "detail": _redact(str(detail))})
    return out


def group(recs: list[dict], top: int = 30) -> dict:
    """The records grouped: by reason (count, members, an example) and per
    model (rates per 100 tool calls)."""
    groups: dict = collections.defaultdict(list)
    for r in recs:
        if r["kind"] != "end":
            groups[(r["kind"], r["code"])].append(r)
    by_reason = [{"kind": k, "code": c, "count": len(rs),
                  "members": len({(r["run"], r["id"]) for r in rs}),
                  "example": next((r["detail"] for r in rs if r["detail"]), "")[:160]}
                 for (k, c), rs in sorted(groups.items(), key=lambda kv: -len(kv[1]))[:top]]
    per_model: dict = collections.defaultdict(collections.Counter)
    calls: collections.Counter = collections.Counter()
    for r in recs:
        per_model[r["model"]][r["kind"]] += 1
        if r["kind"] == "end":
            m = re.search(r"tool_calls=(\d+)", r["detail"])
            calls[r["model"]] += int(m.group(1)) if m else 0
    models = []
    for model, c in sorted(per_model.items(), key=lambda kv: -kv[1]["end"]):
        n = calls[model]
        rate = (lambda k: round(100 * c[k] / n, 1)) if n else (lambda k: None)
        models.append({"model": model or "?", "runs_ended": c["end"], "tool_calls": n,
                       **{f"{k}_per_100": rate(k) for k in ("denial", "api_retry", "stall", "hint")},
                       "crashes": c["crash"]})
    ends = collections.Counter(r["code"] for r in recs if r["kind"] == "end")
    return {"records": len(recs), "by_reason": by_reason, "per_model": models,
            "ended": dict(ends.most_common())}


def report(run_names: list[str] | None = None, top: int = 30) -> dict:
    """The records of the runs grouped: by reason (count, members, an
    example) and per model (rates per 100 tool calls).

    Reads every member log whole, so it is a person's or a tool's command,
    not something the tick runs every minute: the tick keeps the running
    ledger below (:func:`sweep`) and groups it with :func:`group`."""
    return group(collect(run_names), top)


# -- the running ledger the tick keeps -------------------------------------------------

#: How a ledger line's source is named. ``collect`` calls a member's own log
#: ``member``; the ledger adds the member, because a line that only says
#: "member" cannot be traced back to the log it came from.
_SRC = "member:"


def _ledger_path() -> Path:
    return wave.root() / ".friction.jsonl"


def _cursor_path() -> Path:
    return wave.root() / ".friction.cursor.json"


def _load_cursor() -> dict:
    try:
        cur = json.loads(_cursor_path().read_text())
        return cur if isinstance(cur, dict) else {}
    except (OSError, ValueError):
        return {}


def _replace_private(tmp: Path, p: Path) -> None:
    """Move *tmp* onto *p* so the result is readable by its owner only.

    The ledger, the cursor and the report sit beside the wave state, and the
    ledger carries member log excerpts: a world-readable copy in a shared
    waves root is other accounts' reading matter. DELFIN's other writers of
    this kind (attention.py, audit_log.py, approval_code.py) create 0o600.
    """
    os.replace(tmp, p)
    try:
        os.chmod(p, 0o600)
    except OSError:
        pass


def _save_cursor(cur: dict) -> None:
    p = _cursor_path()
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(f"{p.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(cur))
    _replace_private(tmp, p)


def _sources(run: str, tid: str, a: dict) -> list[tuple[Path, str]]:
    """The member's own log, then its real-run logs."""
    out = [(Path(a["log"]), f"{_SRC}{run}/{tid}")]
    try:
        wt = Path(a.get("worktree") or "")
    except TypeError:
        return out
    if wt.is_dir():
        out += [(p, "realrun:" + p.stem) for p in sorted((wt / "scratch" / "realrun")
                                                         .glob("*.jsonl"))]
    return out


def sweep(run_names: list[str] | None = None) -> dict:
    """Append what is new since the last sweep to ``<root>/.friction.jsonl``.

    Each member log is read from the byte the last sweep stopped at, not
    from the top: the tick calls this every pass, and re-reading every log
    of every run each time is what made an unattended count impossible.
    The cursor is written only after the lines reached the ledger, and a
    record that cannot be written leaves the cursor where it was, so the
    next sweep reads that stretch again rather than losing it.

    The ledger is the host's (one file over every run under the root), not
    per run: the per-model rates are across runs, and the root can hold a
    dozen runs whose friction is one picture.

    A stopped member's log is final; a running one's grows as it writes, and
    a line only half written is left for the next pass (the cursor stops at
    the last newline). A log that shrank -- a member resumed onto the same
    path -- is read again from the top, so its second life is counted too.

    Returns ``{"appended": n, "sources": n}``.
    """
    runs = run_names if run_names is not None else wave.runs()
    if not runs:
        return {"appended": 0, "sources": 0}
    with wave._flock(wave.root() / ".friction.lock"):
        cur = _load_cursor()
        new: dict = {}
        n_src = 0
        for run in runs:
            for tid, a in wave.load_state(run)["agents"].items():
                for path, src in _sources(run, tid, a):
                    key = str(path)
                    n_src += 1
                    try:
                        size = path.stat().st_size
                    except OSError:
                        continue
                    start = cur.get(key)
                    start = int(start) if isinstance(start, (int, float)) else 0
                    if start > size:                   # rewritten, not appended
                        start = 0
                    if start == size:
                        continue
                    try:
                        with open(path, "rb") as fh:
                            fh.seek(start)
                            chunk = fh.read()
                    except OSError:
                        continue
                    cut = chunk.rfind(b"\n")
                    if cut < 0:                        # a part line: next read
                        continue                       # takes it from byte 0
                    done = chunk[:cut + 1]
                    model = str((a.get("task") or {}).get("model") or "")
                    new[key] = {"next": start + len(done),
                                "recs": [{"run": run, "id": tid, "model": model, "src": src,
                                          "kind": k, "code": _redact(str(c))[:60],
                                          "detail": _redact(str(d))}
                                         for k, c, d in records(done.decode("utf-8", "replace")
                                                                .splitlines())]}
        appended = 0
        led = _ledger_path()
        led.parent.mkdir(parents=True, exist_ok=True)
        # Created owner-only, and an older world-readable ledger is corrected:
        # the records quote member logs (see _replace_private).
        fd = os.open(str(led), os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
        with os.fdopen(fd, "a") as fh:
            try:
                if stat.S_IMODE(os.fstat(fh.fileno()).st_mode) != 0o600:
                    os.fchmod(fh.fileno(), 0o600)
            except OSError:
                pass
            for key, e in new.items():
                for r in e["recs"]:
                    fh.write(json.dumps(r, default=str) + "\n")
                fh.flush()                             # the cursor must not run
                cur[key] = e["next"]                   # ahead of the ledger
                appended += len(e["recs"])
        _save_cursor(cur)
        return {"appended": appended, "sources": n_src}


def read_ledger() -> list[dict]:
    """The ledger's records ([] when there is none)."""
    try:
        lines = _ledger_path().read_text(errors="replace").splitlines()
    except OSError:
        return []
    out = []
    for line in lines:
        try:
            r = json.loads(line)
        except ValueError:
            continue
        if isinstance(r, dict) and r.get("kind"):
            out.append(r)
    return out


def write_report(top: int = 30) -> dict:
    """Rewrite ``<root>/friction.json`` from the ledger. The tick does this
    at most once an hour; the counts stay readable while the members run."""
    rep = group(read_ledger(), top)
    p = wave.root() / "friction.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(f"{p.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(rep, indent=1, default=str))
    _replace_private(tmp, p)
    return {"path": str(p), "records": rep["records"]}
