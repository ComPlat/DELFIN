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
"""
from __future__ import annotations

import collections
import json
import re
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


def report(run_names: list[str] | None = None, top: int = 30) -> dict:
    """The records grouped: by reason (count, members, an example) and per
    model (rates per 100 tool calls)."""
    recs = collect(run_names)
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
