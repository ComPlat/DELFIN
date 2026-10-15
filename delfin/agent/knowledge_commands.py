"""The person's commands for the layered memory, one implementation for
every surface (REPL, terminal window, dashboard, CLI).

Commands (``handle(name, args, ...)``):

``/memory [private|swarm|project|global|all] [words]``
    every entry with its layer label -- ``private`` (typed memories and
    lessons, proposals marked with their status and ``public`` when they
    would go to a public repository), ``swarm <run>``, ``project``,
    ``global`` -- filterable by layer and words, with provenance for the
    shared ones and ``expired`` where ``review_by`` has passed.
``/share <id> [project|global] [scope=..] [topic=..] [as <handle>]
[release=<snippet>]``
    opens the share dialog (``knowledge_share``). With ``ask`` (REPL,
    window, CLI on /dev/tty) the code is read from the person's next line;
    without it (dashboard) the dialog stays open and
    ``/share confirm <dialog> <code>`` answers it.
``/unshare <entry>``
    withdraws one's own shared entry as a removal branch.
``/lessons [show <id> | edit <id> field=text ... | delete <id>]``
    the person's own lessons; an edit is taken as typed and re-checked
    by the next /share.
``/swarm [runs | harvest <run> | propose <run> <seq> | close <run>]``
    the run's findings, with evidence, for the supervising session.
``/knowledge [lint | dupes | stats | queue | search <words> | accept <q> |
reject <q> <reason> | edit <q> field=.. | merge <keep> <drop> |
merge-into <q> <entry> | move <entry> <topic> | retire <entry> <reason> |
rename <entry> <new> | index | topic add <slug> <text>]``
    the read-only views for everyone; every write is the admin's,
    through the admin's own confirmation.

Nothing here is reachable from a tool call: these are typed commands
(``repl_commands.BUILTINS``, the dashboard router, ``delfin-agent
knowledge``), and the agent's tool surface has no share or curation tool
(pinned by a test over the tool schemas and the executor).
"""
from __future__ import annotations

import re
import shlex
import time
from pathlib import Path
from typing import Callable, Optional

from . import knowledge as kn
from . import knowledge_private as kp
from . import knowledge_share as ks

#: Tests replace the clock the dialogs are timed with.
_CLOCK: Callable[[], float] = time.monotonic

COMMANDS = ("memory", "share", "unshare", "lessons", "swarm", "knowledge")
LAYERS = ("private", "swarm", "project", "global", "all")


def handle(name: str, args: str, *, workspace, ask=None, say=None,
           perms=None) -> str:
    name = name.lstrip("/").lower()
    try:
        words = shlex.split(args or "")
    except ValueError as exc:
        return f"could not read the arguments: {exc}"
    fn = {"memory": _memory, "share": _share, "unshare": _unshare,
          "lessons": _lessons, "swarm": _swarm,
          "knowledge": _knowledge}.get(name)
    if fn is None:
        return f"unknown command /{name}"
    try:
        return fn(words, workspace=Path(workspace or "."), ask=ask,
                  say=say, perms=perms)
    except ks.ShareRefused as exc:
        return f"refused: {exc}"
    except (ValueError, KeyError) as exc:
        return f"/{name}: {exc}"


# ---------------------------------------------------------------------------
# /memory
# ---------------------------------------------------------------------------

def memory_rows(workspace: Path, *, perms=None, today=None) -> list:
    """``(layer, label, text)`` for every entry this person can see."""
    rows = []
    try:
        from .memory_store import list_typed_memories
        for m in list_typed_memories(workspace, scope="all"):
            desc = str(m.get("description") or (m.get("body") or "")[:80])
            rows.append(("private", "private", f"memory {m.get('name')} -- "
                         f"{desc.strip()}"))
    except Exception:
        pass
    try:
        ks.sync_own_lessons()
    except Exception:
        pass
    for les in kp.list_lessons():
        if les.status == kp.WITHDRAWN:
            continue
        label = f"private · proposal {les.id}"
        if les.status == kp.PRIVATE:
            label += f" · would be {kp.target_label(les)}"
            tail = f"  (/share {les.id})"
        else:
            label += f" · {les.status}"
            tail = f"  [{les.shared_ref}]" if les.shared_ref else ""
        if les.note:
            tail += f"  reason: {les.note}"
        rows.append(("private", label, f"{les.title}{tail}"))
    try:
        from . import swarm_memory as sm
        rid = sm.current_run_id(perms)
        for e in (sm.read(rid) if rid else []):
            rows.append(("swarm", f"swarm {rid} #{e.seq}",
                         f"{e.author}: {e.text} [{e.evidence}]"))
    except Exception:
        pass
    from . import repo_visibility as rv
    in_delfin = kn.is_delfin_repo(workspace)
    if not in_delfin:
        root = kn.repo_root_of(workspace)
        vis = rv.label(rv.visibility(root))
        for e in kn.read_entries(kn.project_dir(root), layer=kn.PROJECT):
            rows.append(("project", _shared_label("project", vis, e, today),
                         _shared_text(e)))
    vis = rv.label(rv.visibility(ks.delfin_root()))
    for e in kn.read_entries(kn.global_dir()):
        if e.scope == "delfin-development" and not in_delfin:
            continue
        rows.append(("global", _shared_label("global", vis, e, today),
                     _shared_text(e)))
    return rows


def _shared_label(layer: str, vis: str, e: kn.Entry, today) -> str:
    label = f"{layer} · {vis}"
    if layer == "global":
        label += f" · {e.topic}"
    if e.expired(today):
        label += " · expired"
    return label


def _shared_text(e: kn.Entry) -> str:
    prov = "; ".join(f"{p.author} {p.date}" for p in e.provenance)
    return f"{e.slug}: {e.title}  ({prov})"


def _memory(words, *, workspace, perms=None, **_):
    layer = words[0].lower() if words and words[0].lower() in LAYERS \
        else "all"
    terms = [w.lower() for w in (words[1:] if words and words[0].lower()
                                 in LAYERS else words)]
    rows = [r for r in memory_rows(workspace, perms=perms)
            if layer == "all" or r[0] == layer]
    if terms:
        rows = [r for r in rows if all(t in (r[1] + " " + r[2]).lower()
                                       for t in terms)]
    lines = [f"  [{label}] {text}" for _l, label, text in rows] or [
        f"  no {layer + ' ' if layer != 'all' else ''}entries"]
    lines.append("  layers: private (this machine only) · swarm (this run) "
                 "· project (.delfin/knowledge, via git) · global (every "
                 "DELFIN user, via git)")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# dialogs
# ---------------------------------------------------------------------------

def _drive(d: "ks.Dialog", ask, say) -> str:
    """Show the dialog on the person's channel and read their answer."""
    say = say or (lambda _s: None)
    say(d.preview())
    if d.blocked:
        return "nothing was shared: the check blocked it (see above)"
    if ask is None:
        return (f"answer with: /share confirm {d.id} <code>   (or /share "
                f"cancel {d.id})")
    for _ in range(ks.MAX_WRONG + 2):
        try:
            typed = str(ask("code: ") or "")
        except (EOFError, KeyboardInterrupt):
            ks.cancel(d.id)
            return "cancelled; nothing was written"
        if not typed.strip():
            ks.cancel(d.id)
            return "cancelled; nothing was written"
        out = ks.confirm(d.id, typed, clock=_CLOCK)
        if out.status == "next":
            say(out.message)
            continue
        if out.status == "refused" and out.dialog is not None:
            say(out.message)
            continue
        return out.message
    ks.cancel(d.id)
    return "cancelled; nothing was written"


def _confirm_words(words) -> str:
    if len(words) < 3:
        return "usage: /share confirm <dialog> <code>"
    out = ks.confirm(words[1], words[2], clock=_CLOCK)
    if out.status == "next":
        return out.message + f"\nanswer with: /share confirm {words[1]} " \
                             "<code>"
    return out.message


def _share(words, *, workspace, ask=None, say=None, **_):
    if not words:
        return ("usage: /share <id> [project|global] [scope=always|"
                "delfin-development] [topic=<topic>] [as <handle>] "
                "[release=<snippet>]")
    if words[0] == "confirm":
        return _confirm_words(words)
    if words[0] == "cancel":
        return "cancelled" if len(words) > 1 and ks.cancel(words[1]) \
            else "no such dialog"
    les = kp.get(words[0])
    if les is None:
        raise KeyError(f"no lesson {words[0]} (/memory private)")
    layer, scope, topic, handle, release = "", "", "", "", []
    rest = words[1:]
    i = 0
    while i < len(rest):
        w = rest[i]
        if w in ("project", "global"):
            layer = w
        elif w.startswith("scope="):
            scope = w[6:]
        elif w.startswith("topic="):
            topic = w[6:]
        elif w.startswith("release="):
            release.append(w[8:])
        elif w == "as" and i + 1 < len(rest):
            handle = rest[i + 1]
            i += 1
        else:
            raise ValueError(f"unknown word {w!r}")
        i += 1
    if not layer:
        layer = les.prepared_layer or (
            "project" if (kn.repo_root_of(workspace) / ".git").exists()
            and not kn.is_delfin_repo(workspace) else "global")
    d = ks.start_share(les.id, layer, workspace=workspace, scope=scope,
                       topic=topic, handle=handle, releasable=release,
                       clock=_CLOCK)
    return _drive(d, ask, say)


def _unshare(words, *, workspace, ask=None, say=None, **_):
    if not words:
        return "usage: /unshare <entry> [as <handle>]"
    handle = words[2] if len(words) > 2 and words[1] == "as" else ""
    d = ks.start_unshare(words[0], workspace=workspace, handle=handle,
                         clock=_CLOCK)
    return _drive(d, ask, say)


# ---------------------------------------------------------------------------
# /lessons
# ---------------------------------------------------------------------------

def _fields(words) -> dict:
    out = {}
    key = None
    for w in words:
        m = re.match(r"^(title|situation|lesson|why|topic)=(.*)$", w)
        if m:
            key = m.group(1)
            out[key] = m.group(2)
        elif key:
            out[key] += " " + w
        else:
            raise ValueError(f"expected field=text, got {w!r}")
    return out


def _lessons(words, **_):
    if not words:
        rows = kp.list_lessons()
        if not rows:
            return "no lessons yet"
        return "\n".join(f"  {x.id:>3} [{x.status}] {x.title}"
                         for x in rows)
    verb = words[0]
    if verb == "show" and len(words) > 1:
        les = kp.get(words[1])
        if les is None:
            raise KeyError(f"no lesson {words[1]}")
        return (f"{les.text}\nevidence: {', '.join(les.evidence)}\n"
                f"topic: {les.topic or '-'}  status: {les.status}"
                + (f"  reason: {les.note}" if les.note else ""))
    if verb == "edit" and len(words) > 2:
        les = kp.edit(words[1], **_fields(words[2:]))
        return f"lesson {les.id} updated; /share {les.id} checks it again"
    if verb == "delete" and len(words) > 1:
        return "deleted" if kp.delete(words[1]) else "no such lesson"
    return "usage: /lessons [show <id> | edit <id> field=text | delete <id>]"


# ---------------------------------------------------------------------------
# /swarm
# ---------------------------------------------------------------------------

def _swarm(words, *, perms=None, **_):
    from . import swarm_memory as sm
    if not words or words[0] == "runs":
        rows = sm.runs()
        here = sm.current_run_id(perms)
        if not rows:
            return "no swarm runs on this machine"
        return "\n".join(f"  {r} ({state}){'  <- this session' if r == here else ''}"
                         for r, state in rows)
    if words[0] == "harvest" and len(words) > 1:
        rows = sm.harvest(words[1])
        if not rows:
            return f"no entries in run {words[1]}"
        lines = [f"  #{e.seq} {e.author} "
                 f"{time.strftime('%Y-%m-%d %H:%M', time.localtime(e.time))}"
                 f": {e.text}  [evidence {e.evidence}]" for e in rows]
        lines.append(f"  /swarm propose {words[1]} <seq> keeps one as your "
                     "private lesson; /share then offers it")
        return "\n".join(lines)
    if words[0] == "propose" and len(words) > 2:
        les = sm.to_lesson(words[1], int(words[2]))
        return f"private lesson {les.id}: {les.title} (/share {les.id})"
    if words[0] == "close" and len(words) > 1:
        return f"archived: {sm.close(words[1])}"
    return "usage: /swarm [runs | harvest <run> | propose <run> <seq> | " \
           "close <run>]"


# ---------------------------------------------------------------------------
# /knowledge
# ---------------------------------------------------------------------------

def _knowledge(words, *, workspace, ask=None, say=None, **_):
    from . import knowledge_curation as kc
    verb = words[0] if words else "list"
    args = words[1:]
    if verb == "list":
        by_topic: dict = {}
        for e in kn.read_entries(kn.global_dir()):
            by_topic.setdefault(e.topic, []).append(e)
        if not by_topic:
            return "the global knowledge has no entries yet"
        lines = []
        for topic in sorted(by_topic):
            lines.append(f"  {topic}")
            lines += [f"    {_shared_text(e)}" for e in by_topic[topic]]
        return "\n".join(lines)
    if verb == "search":
        terms = [a.lower() for a in args]
        hits = [e for e in kn.read_entries(kn.global_dir())
                if all(t in e.text.lower() for t in terms)]
        return "\n".join(f"  [{e.topic}] {_shared_text(e)}" for e in hits) \
            or "no match"
    if verb == "lint":
        problems = kc.lint()
        return "\n".join(f"  - {p}" for p in problems) or "lint: clean"
    if verb == "dupes":
        pairs = kc.dupes()
        return "\n".join(f"  {p.kind}: {p.a} / {p.b}" for p in pairs) \
            or "no duplicates or contradictions"
    if verb == "stats":
        st = kc.stats()
        lines = [f"  entries: {st['entries']}",
                 "  by topic: " + ", ".join(f"{k} {v}" for k, v in
                                            sorted(st["by_topic"].items())),
                 "  expired (re-review): " + (", ".join(st["expired"])
                                              or "none"),
                 "  never loaded here (retire candidates): "
                 + (", ".join(st["never_loaded"]) or "none")]
        lines += [f"  broken {b.kind} in {b.slug}: {b.ref}"
                  for b in st["broken"]]
        return "\n".join(lines)
    if verb == "queue":
        if not ks.is_admin():
            return "the review queue is the knowledge admin's"
        recs = [r for r in ks.read_queue() if r.get("status") == "pending"]
        if not recs:
            return "the review queue is empty"
        return "\n".join(
            f"  {r['id']} by {r.get('author')} [{r.get('topic')}, "
            f"{r.get('mode')}]" + (f" contradicts {', '.join(r['contradicts'])}"
                                   if r.get("contradicts") else "")
            + "\n" + "\n".join("      " + ln for ln in
                               str(r.get("entry", "")).splitlines())
            for r in recs)
    plan = _plan(kc, verb, args)
    if plan is None:
        return ("usage: /knowledge [list | search <words> | lint | dupes | "
                "stats | queue | accept <q> | reject <q> <reason> | edit <q> "
                "field=.. | merge <keep> <drop> | merge-into <q> <entry> | "
                "move <entry> <topic> | retire <entry> <reason> | rename "
                "<entry> <new> | index | topic add <slug> <text>]")
    d = kc.start(plan, clock=_CLOCK)
    return _drive(d, ask, say)


def _plan(kc, verb, args):
    if verb == "accept" and args:
        return kc.plan_accept(args[0])
    if verb == "reject" and len(args) > 1:
        return kc.plan_reject(args[0], " ".join(args[1:]))
    if verb == "edit" and len(args) > 1:
        return kc.plan_queue_edit(args[0], **_fields(args[1:]))
    if verb == "merge" and len(args) == 2:
        return kc.plan_merge(args[0], args[1])
    if verb == "merge-into" and len(args) == 2:
        return kc.plan_queue_merge(args[0], args[1])
    if verb == "move" and len(args) == 2:
        return kc.plan_move(args[0], args[1])
    if verb == "retire" and len(args) > 1:
        return kc.plan_retire(args[0], reason=" ".join(args[1:]))
    if verb == "rename" and len(args) == 2:
        return kc.plan_rename(args[0], args[1])
    if verb == "index":
        return kc.plan_index()
    if verb == "topic" and len(args) > 1 and args[0] == "add":
        return kc.plan_topic_add(args[1], " ".join(args[2:]))
    return None


def end_of_turn_hints() -> list:
    """The one-line share hints for lessons proposed this turn (each once).
    Never raises: a hint is never worth a failed turn."""
    try:
        return kp.pending_hints()
    except Exception:
        return []


def tty_io() -> tuple:
    """``(ask, say)`` bound to the controlling terminal, or ``(None, print)``
    when there is none -- then a write action stays unanswerable, which is
    the point: the CLI's stdout may be an agent's pipe."""
    try:
        fh = open("/dev/tty", "r+", encoding="utf-8")
    except OSError:
        return None, print

    def say(text: str) -> None:
        fh.write(str(text) + "\n")
        fh.flush()

    def ask(prompt: str) -> str:
        fh.write(prompt)
        fh.flush()
        line = fh.readline()
        if not line:
            raise EOFError
        return line.rstrip("\n")
    return ask, say


def run_cli(argv: list, *, workspace: Optional[Path] = None) -> int:
    """``delfin-agent knowledge <verb> ...``: the /knowledge command with
    the code dialog on /dev/tty."""
    ask, say = tty_io()
    out = handle("knowledge", " ".join(shlex.quote(a) for a in argv),
                 workspace=workspace or Path.cwd(), ask=ask, say=say)
    print(out)
    return 0 if not out.startswith(("refused", "/knowledge:")) else 1
