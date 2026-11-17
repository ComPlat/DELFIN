"""Friction report: where real use lost turns and time, ranked, with tasks.

Input: two durable stores of this machine, read only.

* the audit log (``~/.delfin/audit.log`` and its weekly rotations
  ``audit-YYYY-Www.log``): every refused action (``decision: denied``)
  with its reason and command;
* the turn metrics (``agent_metrics.jsonl``): duration, model, silent
  exits and tool errors per turn.

Output: friction buckets ranked by count. Each bucket carries:
- the sessions it touched;
- the most frequent command heads;
- the models involved;
- one task proposal, a fixed text per bucket, for the person to approve
  or reject.

Nothing is changed or started. ``--save`` keeps the day's JSON under
``~/.delfin/friction/`` so the next run can show the change per bucket.

Semantics:
- A refusal is bucketed by the cause after "was not run:" in its reason.
  The first matching pattern in ``BUCKETS`` wins.
- Its command is reduced to a head: ``pytest``, ``git commit``,
  ``python -c``, ``grep``.
- Reason texts, arguments, paths and message texts never reach the
  output, so the report can be shown or sent anywhere.
- A turn of at least ``SLOW_TURN_S`` counts as slow.

Rejected alternatives:
- **Asking a model to summarise the log:** the report must be the same for
  the same input. It must also not depend on the endpoint whose stalls it
  counts.
- **Adding failure_log.jsonl:** it records the same refusals as the audit
  log and would count each of them twice.
"""

from __future__ import annotations

import json
import os
import re
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

#: A turn at or above this many seconds is a slow turn.
SLOW_TURN_S = 300.0
#: A bucket needs this many events before it gets a task proposal.
MIN_FOR_TASK = 3

# (key, pattern over the refusal cause, title, task proposal). Order matters.
BUCKETS: list[tuple[str, re.Pattern[str], str, str]] = [
    ("watcher_unavailable",
     re.compile(r"could not judge|confirming model could not answer", re.I),
     "The action watcher could not judge (model unreachable or slow)",
     "Find which of these commands the deterministic rules already know "
     "to stay in the workspace, and screen them without the model. "
     "Measure the confirming model's answer time under load."),
    ("watcher_refused",
     re.compile(r"action classifier refused|already refused by the action "
                r"classifier", re.I),
     "Allowed by the rules, then refused by the action watcher",
     "Replay the refused command heads against the watcher prompt; a "
     "read-only head refused more than once is a prompt or rule defect."),
    # Before "timeout": an approval window running out is a person away,
    # not a command that ran too long.
    ("approval_unanswered",
     # Upper case on purpose: the approval broker writes "TIMED OUT", the
     # command limit "timed out after".
     re.compile(r"TIMED OUT|(?i:the user is away|approval window)"),
     "Asked a person, nobody answered in time",
     "Sort these asks into those a rule could have decided (read-only, "
     "inside the workspace) and the rest; the rest should reach the "
     "person where they are (phone) instead of timing out."),
    ("guarded_location",
     re.compile(r"where approvals are answered|secret-deny|\.ssh|"
                r"credential", re.I),
     "Touched a guarded location (approvals, secrets)",
     "Refused by design: check that the refusal says what to do instead "
     "and that no session repeats it after the first refusal."),
    ("refusal_pause",
     re.compile(r"three actions were refused", re.I),
     "Session paused after three refusals",
     "Follows from the buckets above; check that the pause names the "
     "three refusals and that one answer lifts it."),
    ("read_outside_workspace",
     re.compile(r"outside the allowed workspace roots|escapes workspace",
                re.I),
     "Read or write outside the workspace roots",
     "List the outside paths by kind (temp files the agent wrote itself, "
     "the venv, the home directory) and decide per kind whether the agent "
     "should use a workspace path instead or the root should be granted."),
    ("unattended_needs_person",
     re.compile(r"no confirm_?callback is configured|nobody can be asked",
                re.I),
     "Needed a person, none was attached",
     "Route these to the approval queue (dashboard or phone) instead of "
     "refusing, or tell the agent up front which paths need a person."),
    ("deny_pattern",
     re.compile(r"deny-pattern", re.I),
     "A fixed deny pattern matched",
     "Check that the refusal text names the allowed way (for example "
     "naming files instead of rm -rf); count repeats in one session."),
    ("read_before_write",
     re.compile(r"read_file|read baseline", re.I),
     "Write without a prior read",
     "Check whether the agent read the file through bash first; if so, "
     "a bash read should count as the read baseline."),
    ("edit_no_match",
     re.compile(r"old_string|new_string must differ", re.I),
     "Edit text did not match",
     "Look at the edit tool's mismatch hint: does it show the nearest "
     "lines so the next attempt matches?"),
    ("timeout",
     re.compile(r"timed out|foreground ceiling", re.I),
     "Command ran into its time limit",
     "Move long runs to bash_background by default when the command is "
     "a test suite or a build."),
]
OTHER = ("other", None, "Other refusals",
         "Read a sample of these by hand and add a bucket if one repeats.")

_SLOW = ("slow_turns", None, "Turns of five minutes or more",
         "Split the time of these turns into model time and tool time, "
         "and compare with the same model's median turn.")
_SILENT = ("silent_exits", None, "Turns that ended without an answer",
           "Replay one such turn; a silent exit is an answer the person "
           "never got.")

_TITLES = {k: t for k, _, t, _ in BUCKETS + [OTHER, _SLOW, _SILENT]}
_TASKS = {k: p for k, _, _, p in BUCKETS + [OTHER, _SLOW, _SILENT]}

_SKIP_PREFIX = re.compile(
    r"^\s*(?:set\s+-[a-z]+\s+\w+\s*;\s*|cd\s+\S+\s*(?:&&|;)\s*|"
    r"env(?:\s+-u\s+\S+)+\s+|\w+=\S*\s+|time\s+|nice\s+(?:-n\s*\d+\s+)?|"
    r"timeout\s+(?:-\S+\s+)*\d+\w?\s+|(?:for|while)\s[^;]*;\s*do\s+)")
_SUBCOMMAND_TOOLS = {"git", "gh", "pip", "uv", "docker", "systemctl",
                     "npm", "cargo", "delfin-agent"}
_HEAD_CHARS = re.compile(r"[^\w.+-]")


def command_head(command: str) -> str:
    """The first program of ``command`` and, for git-like tools, its verb.

    ``set -o pipefail; ~/x/.venv/bin/python -m pytest -q tests/a.py``
    becomes ``pytest``; ``git commit -am "msg"`` becomes ``git commit``.
    No argument, path or quoted text is kept. An empty or unreadable
    command gives "".
    """
    s = str(command or "")
    for _ in range(6):
        t = _SKIP_PREFIX.sub("", s, count=1)
        if t == s:
            break
        s = t
    s = re.split(r"\s*(?:\|\||&&|[;|\n])\s*", s.strip(), maxsplit=1)[0]
    toks = s.split()
    if not toks:
        return ""
    prog = os.path.basename(toks[0].strip("\"'"))
    rest = toks[1:]
    if re.fullmatch(r"python[\d.]*", prog):
        if rest[:1] == ["-m"] and len(rest) > 1:
            head = rest[1]
        elif rest[:1] == ["-c"]:
            head = "python -c"
        elif rest[:1] == ["-"]:
            head = "python stdin"
        else:
            head = "python script" if rest else "python"
    elif prog in _SUBCOMMAND_TOOLS:
        verb, skip = "", False
        for r in rest:
            if skip or r.startswith("-"):
                # git -C DIR / -c KEY=VAL take a value; it is no verb
                skip = not skip and r in ("-C", "-c")
                continue
            verb = r
            break
        head = f"{prog} {verb}".strip()
    else:
        head = prog
    return " ".join(_HEAD_CHARS.sub("", p) for p in head.split())[:32]


def bucket_of(reason: str) -> str:
    """The bucket key for one refusal reason."""
    text = str(reason or "")
    cut = text.find(" was not run: ")
    cause = text[cut + len(" was not run: "):] if cut >= 0 else text
    for key, pat, _, _ in BUCKETS:
        if pat.search(cause):
            return key
    return OTHER[0]


def _ts(value: Any) -> float:
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return datetime.fromisoformat(str(value).rstrip("Z")).replace(
            tzinfo=timezone.utc).timestamp()
    except ValueError:
        return 0.0


def _audit_files(audit_path: Path, since_s: float) -> list[Path]:
    """The active log and every rotation written to since ``since_s``."""
    files = [p for p in audit_path.parent.glob("audit-*.log")
             if p.is_file() and p.stat().st_mtime >= since_s]
    if audit_path.is_file():
        files.append(audit_path)
    return files


def _jsonl(path: Path) -> Iterable[dict]:
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                try:
                    rec = json.loads(line)
                except ValueError:
                    continue
                if isinstance(rec, dict):
                    yield rec
    except OSError:
        return


def _new_bucket() -> dict:
    return {"count": 0, "sessions": set(), "heads": Counter(),
            "models": Counter()}


def collect(*, since_s: float, home: "Path | str | None" = None) -> dict:
    """Aggregate the friction since ``since_s`` (unix time) on this machine.

    ``home`` is the directory holding ``.delfin`` (default: the user's
    home). Missing or corrupt stores give empty counts, never an error.
    """
    base = Path(home or Path.home()) / ".delfin"
    buckets: dict[str, dict] = defaultdict(_new_bucket)
    actions = 0
    for path in _audit_files(base / "audit.log", since_s):
        for rec in _jsonl(path):
            if _ts(rec.get("ts")) < since_s or not rec.get("tool"):
                continue
            actions += 1
            if rec.get("decision") != "denied":
                continue
            b = buckets[bucket_of(rec.get("reason", ""))]
            b["count"] += 1
            if rec.get("session_id"):
                b["sessions"].add(str(rec["session_id"]))
            tool = str(rec.get("tool") or "")
            head = (command_head(rec.get("command", ""))
                    if tool.startswith("bash") else "")
            b["heads"][head or tool] += 1

    turns = 0
    durations: dict[str, list[float]] = defaultdict(list)
    for rec in _jsonl(base / "agent_metrics.jsonl"):
        if _ts(rec.get("ts")) < since_s:
            continue
        turns += 1
        model = str(rec.get("model") or "unknown")
        dur = float(rec.get("duration_s") or 0.0)
        durations[model].append(dur)
        sid = str(rec.get("session_id") or "")
        for key, hit in ((_SLOW[0], dur >= SLOW_TURN_S),
                         (_SILENT[0], bool(rec.get("silent_exit")))):
            if hit:
                b = buckets[key]
                b["count"] += 1
                b["models"][model] += 1
                if sid:
                    b["sessions"].add(sid)

    ranked = []
    for key, b in sorted(buckets.items(), key=lambda kv: (-kv[1]["count"],
                                                          kv[0])):
        ranked.append({
            "key": key, "title": _TITLES[key], "count": b["count"],
            "sessions": len(b["sessions"]),
            "heads": b["heads"].most_common(5),
            "models": b["models"].most_common(3),
            "task": _TASKS[key] if b["count"] >= MIN_FOR_TASK else "",
        })
    medians = {m: sorted(v)[len(v) // 2] for m, v in durations.items()}
    return {"since_s": since_s, "generated_s": time.time(),
            "actions": actions, "turns": turns,
            "median_turn_s": dict(sorted(medians.items())),
            "buckets": ranked}


def previous_report(save_dir: Path, before_day: str) -> "dict | None":
    """The newest saved report of a day earlier than ``before_day``."""
    days = sorted(p.stem for p in save_dir.glob("*.json")
                  if p.stem < before_day)
    if not days:
        return None
    try:
        return json.loads((save_dir / f"{days[-1]}.json").read_text())
    except (OSError, ValueError):
        return None


def save(data: dict, save_dir: Path) -> Path:
    """Write ``data`` as ``<save_dir>/YYYY-MM-DD.json``, owner-only."""
    save_dir.mkdir(parents=True, exist_ok=True)
    os.chmod(save_dir, 0o700)
    day = time.strftime("%Y-%m-%d", time.gmtime(data["generated_s"]))
    path = save_dir / f"{day}.json"
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1)
    return path


def render_text(data: dict, previous: "dict | None" = None) -> str:
    """Terminal rendering: the ranked buckets, then the proposed tasks."""
    before = {b["key"]: b["count"]
              for b in (previous or {}).get("buckets", [])}
    since = time.strftime("%Y-%m-%d %H:%M UTC",
                          time.gmtime(data.get("since_s") or 0))
    lines = [f"Friction since {since}: {data.get('actions', 0)} actions, "
             f"{data.get('turns', 0)} recorded turns", ""]
    buckets = data.get("buckets") or []
    if not buckets:
        lines.append("No refusals, slow turns or silent exits in the window.")
        return "\n".join(lines)
    for i, b in enumerate(buckets, 1):
        delta = ""
        if previous is not None:
            d = b["count"] - before.get(b["key"], 0)
            delta = f" ({'+' if d >= 0 else ''}{d} vs last report)"
        lines.append(f"{i}. {b['title']}: {b['count']}{delta}, "
                     f"{b['sessions']} session(s)")
        heads = ", ".join(f"{h}×{n}" for h, n in b.get("heads") or [])
        models = ", ".join(f"{m}×{n}" for m, n in b.get("models") or [])
        if heads:
            lines.append(f"   commands: {heads}")
        if models:
            lines.append(f"   models: {models}")
    if data.get("median_turn_s"):
        lines.append("")
        lines.append("Median turn: " + ", ".join(
            f"{m} {s:.0f} s" for m, s in data["median_turn_s"].items()))
    tasks = [b for b in buckets if b.get("task")]
    if tasks:
        lines += ["", "Proposed tasks (nothing starts without your yes):"]
        for i, b in enumerate(tasks, 1):
            lines.append(f"T{i} [{b['key']}, {b['count']}] {b['task']}")
    return "\n".join(lines)
