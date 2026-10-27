"""``/schedule`` in a session and ``delfin-agent schedule``: the words.

Input: the words after the command, the folder the schedule runs in (the
session's, or ``--cwd``), and who writes it. Output: the text to show.

    list                        every schedule, with its next run
    add <when> <task>           a new one, confirmed (the user wrote it)
    show N · ok N · pause N · resume N · delete N · run N · log [N]
    options for add: --usd X (a run's ceiling) --day-usd X --minutes N
                     --notify always|errors|never --cwd PATH

``ok N`` confirms a schedule the agent proposed, after showing what it
does. A session reads the user's answer through its own prompt (``ask``);
on a surface that cannot ask, ``ok N yes`` confirms what ``show N``
displayed. The command line refuses every change from inside an agent
session, and the store is masked there in any case (``schedules``).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Callable

from . import schedules
from .schedules import ScheduleError

USAGE = ("/schedule [list] · /schedule add <when> <task> · "
         "/schedule show|ok|pause|resume|delete|run N · /schedule log [N]")

_OPTIONS = {"--usd": "budget_run_usd", "--day-usd": "budget_day_usd",
            "--minutes": "max_minutes", "--notify": "notify", "--cwd": "workspace"}


def _service_line() -> str:
    """Whether the agent service that runs the schedules runs here."""
    try:
        from .service import client, rundir
        hello, _why = client.probe(timeout=1.0)
        host = rundir.hostname()
    except Exception:
        return ""
    if hello is not None:
        return f"runs in the agent service on {host} (pid {hello.get('pid')})"
    return (f"the agent service is not running on {host}: schedules run only "
            "while it does -- `delfin-agent schedule autostart` says how to "
            "keep it up")


def _status(e: dict, host: str) -> str:
    if e.get("host") != host:
        return f"on {e.get('host')}"
    if not e["confirmed"]:
        return "needs ok" if not e.get("confirmed_sha") else "changed: needs ok"
    if e.get("paused"):
        return "paused"
    if e.get("running"):
        return "running"
    return "active"


def listing(host: str | None = None) -> str:
    host = host or schedules._hostname()
    rows = schedules.entries()
    if not rows:
        return ("no schedules. /schedule add täglich 07:30 <task> makes one "
                f"(the time as you would say it; {schedules.EXAMPLES})")
    out = []
    for e in rows:
        last = f" · last {e['last_status']}" if e.get("last_status") else ""
        out.append(f"#{e['id']:<3} {e['when']} ({e['tz']}) · next "
                   f"{schedules.describe_next(e)} · {_status(e, host)}{last}\n"
                   f"     {e['task'][:70]}")
    return "\n".join(out)


def show(entry: dict) -> str:
    spent = schedules.spent_today(entry)
    return "\n".join([
        f"#{entry['id']} {entry['when']} ({entry['tz']}), next {schedules.describe_next(entry)}",
        f"  task:    {entry['task']}",
        f"  folder:  {entry['workspace']}",
        f"  host:    {entry['host']}",
        f"  limits:  ${entry['budget_run_usd']:.2f} a run, ${entry['budget_day_usd']:.2f} "
        f"a day (today ${spent:.2f}), {entry['max_minutes']} min a run",
        f"  report:  {entry['notify']}",
        f"  state:   {'confirmed' if entry['confirmed'] else 'NOT confirmed'}"
        f"{', paused' if entry.get('paused') else ''}"
        f" · by {entry.get('created_by') or '?'}"
        f" · check {schedules.fingerprint(entry)[:8]}",
    ])


def _options(words: list[str]) -> tuple[list[str], dict]:
    rest, opts = [], {}
    i = 0
    while i < len(words):
        key = words[i]
        if key in _OPTIONS and i + 1 < len(words):
            value = words[i + 1]
            name = _OPTIONS[key]
            try:
                opts[name] = (float(value) if name.startswith("budget")
                              else int(value) if name == "max_minutes" else value)
            except ValueError:
                raise ScheduleError(f"{key} needs a number", "budget") from None
            i += 2
            continue
        rest.append(key)
        i += 1
    return rest, opts


def handle(args: str, *, workspace: Path | str, created_by: str,
           ask: Callable[[str], str] | None = None) -> str:
    """One ``/schedule`` line. Never raises."""
    words = str(args or "").split()
    sub = words[0].lower() if words else "list"
    rest = words[1:]
    try:
        if sub in ("list", "ls"):
            return listing()
        if sub in ("help", "?"):
            return f"{USAGE}\n{schedules.EXAMPLES}"
        if sub in ("add", "new"):
            body, opts = _options(rest)
            ws = opts.pop("workspace", None) or str(workspace)
            entry = schedules.add(" ".join(body), workspace=str(Path(ws).expanduser()),
                                  created_by=created_by, confirm=True, **opts)
            line = _service_line()
            return (f"schedule #{entry['id']} made: {entry['when']} ({entry['tz']}), "
                    f"next {schedules.describe_next(entry)}\n  {entry['task'][:200]}"
                    + (f"\n  {line}" if line else ""))
        if sub == "log":
            rows = schedules.runs(limit=200)
            if rest:
                wanted = schedules._id(rest[0])
                rows = [r for r in rows if r["id"] == wanted]
            rows = rows[-10:]
            if not rows:
                return "no runs yet"
            return "\n".join(
                f"{time.strftime('%d.%m. %H:%M', time.localtime(r['ts']))} #{r['id']} "
                f"{r['status']}" + (f" ${r['costUsd']:.2f}" if r.get("costUsd") else "")
                + (f" · session {r['session'][:8]}" if r.get("session") else "")
                + (f"\n    {r['text'][:160]}" if r.get("text") else "")
                for r in rows)
        if not rest:
            return f"usage: /schedule {sub} N"
        entry_id = rest[0]
        if sub == "show":
            return show(schedules.get(entry_id))
        if sub in ("ok", "confirm"):
            entry = schedules.get(entry_id)
            if entry["confirmed"]:
                return f"#{entry['id']} is confirmed already"
            shown = schedules.fingerprint(entry)
            yes = len(rest) > 1 and rest[1].lower() in ("yes", "ja", "y")
            if not yes:
                if ask is None:
                    return (show(entry) + f"\nconfirm with: /schedule ok {entry['id']} yes")
                answer = str(ask(show(entry) + "\nLet this run as shown? [y/N] ") or "")
                if answer.strip().lower() not in ("y", "yes", "j", "ja"):
                    return f"#{entry['id']} left unconfirmed"
            done = schedules.confirm(entry_id, expect=shown)
            return f"#{done['id']} confirmed: next {schedules.describe_next(done)}"
        if sub == "pause":
            e = schedules.pause(entry_id)
            return f"#{e['id']} paused"
        if sub in ("resume", "continue"):
            e = schedules.resume(entry_id)
            return f"#{e['id']} resumed: next {schedules.describe_next(e)}"
        if sub in ("delete", "rm", "remove"):
            e = schedules.delete(entry_id)
            return f"#{e['id']} deleted ({e['task'][:60]})"
        if sub in ("run", "now"):
            e = schedules.get(entry_id)
            if e["host"] != schedules._hostname():
                return f"#{e['id']} runs on {e['host']}; start it there"
            if not e["confirmed"]:
                return f"#{e['id']} is not confirmed: /schedule ok {e['id']}"
            schedules.request_run(entry_id)
            line = _service_line()
            return f"#{e['id']} runs within a minute" + (f"; {line}" if line else "")
        return USAGE
    except ScheduleError as exc:
        return f"schedule: {exc}"
    except Exception as exc:                        # noqa: BLE001
        return f"schedule failed: {type(exc).__name__}: {exc}"


AUTOSTART = """The schedules run inside the agent service. To keep it up on this host:

  systemd (as your user, no root):
    mkdir -p ~/.config/systemd/user
    printf '[Service]\\nExecStart={python} -m delfin.agent.service.supervisor --idle-exit 315360000\\nRestart=on-failure\\n[Install]\\nWantedBy=default.target\\n' > ~/.config/systemd/user/delfin-agent-service.service
    systemctl --user enable --now delfin-agent-service
    (after a reboot it starts when you log in; `loginctl enable-linger` -- if the
     admin allows it -- starts it without a login)

  cron, where systemd --user is not available:
    crontab -e, and add:
    */10 * * * * {python} -m delfin.agent.service.launcher >/dev/null 2>&1

  On the Pi, scripts/pi/pi_messenger.sh installs it as a system unit.
Nothing was changed."""


def cli_main(argv: list[str] | None = None) -> int:
    """`delfin-agent schedule ...`."""
    # Not argparse: the task is free text, and a word in it that starts
    # with a dash is the task's, not an option.
    words = list(sys.argv[2:] if argv is None else argv)
    if words and words[0] in ("-h", "--help"):
        print(f"delfin-agent schedule {USAGE.replace('/schedule ', '')}\n"
              f"{schedules.EXAMPLES}\ndelfin-agent schedule autostart: keep the "
              "agent service up on this host")
        return 0
    yes = any(w in ("--yes", "-y") for w in words)
    words = [w for w in words if w not in ("--yes", "-y")]
    sub = (words[0].lower() if words else "list")
    if sub == "autostart":
        print(AUTOSTART.format(python=sys.executable))
        return 0
    from .approval_guard import inside_agent_session
    if sub not in ("list", "ls", "show", "log", "help") and inside_agent_session():
        print("ERROR: refused: this process was started by an agent session. "
              "Schedules are the user's to make and change; ask them to run it.",
              file=sys.stderr)
        return 1
    if sub in ("ok", "confirm") and yes:
        words.append("yes")
    cwd = ""
    if "--cwd" in words:
        at = words.index("--cwd")
        cwd = words[at + 1] if at + 1 < len(words) else ""
        del words[at:at + 2]
    workspace = Path(cwd or ".").expanduser().resolve()

    def ask(prompt: str) -> str:
        if not sys.stdin.isatty():
            return ""
        return input(prompt)

    text = handle(" ".join(words), workspace=workspace, created_by="cli", ask=ask)
    print(text)
    if sub in ("add", "new", "run", "now", "resume") and not text.startswith("schedule"):
        try:
            from .service import launcher
            launcher.ensure_service()
        except Exception as exc:                    # noqa: BLE001
            print(f"note: the agent service did not start: {exc}", file=sys.stderr)
    return 1 if text.startswith("schedule") else 0


__all__ = ["handle", "listing", "show", "cli_main", "USAGE"]
