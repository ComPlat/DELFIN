"""Watchers: the event a schedule waits for, and the cheap look at it.

Input: a watcher as written (``wenn job 123 fertig``, ``when file
calc.log changes``, ``wenn ausgabe "squeue -u max" sich ändert``) and,
when looking, the schedule's folder and the state of the last look.
Output: ``(spec, words used)`` when reading; ``(state, fired, note)`` when
looking -- the note names what happened, for the run's prompt and the
report.

Events:

    job N fertig|done|finishes      a Slurm job left the queue (squeue), with
                                    its end state (sacct) -- once
    pid|prozess N endet|exits       a process of this user ended -- once
    datei|file P erscheint|appears  a file inside the folder exists -- once
    datei|file P ändert sich|changes   its time stamp or size changed
    ausgabe|output "CMD" ändert sich|changes   a reading command's output
    ... alle|every N min            how often to look (default 60 s, >= 30 s)

The first look at a change only learns the state; a change fires on every
change after it.

Rules a watcher cannot pass, from whichever surface it was made:

* A command is one of ``READING_COMMANDS``, split without a shell, run
  with a clean environment (no credential variables), in the folder, under
  ``COMMAND_TIMEOUT_S``, its output cut to ``OUTPUT_LIMIT`` characters.
  Options that read a file of their own choosing (``grep -f``) are refused.
* A path -- the watched file, or a file argument of a command that reads
  files -- must lie inside the schedule's folder, links resolved; checked
  when the watcher is made and again at every look. Commands that take no
  paths (the queue, the disks) may not carry a ``/`` at all.

So a watcher can neither change anything nor read past its folder, and
nothing it runs needs an approval a session would have asked for.

Rejected: a shell command of the user's choosing (a phone message could
then run what a session would have to ask for); following symlinks out of
the folder; a check through the model (every look would cost a turn).
"""

from __future__ import annotations

import hashlib
import os
import re
import shlex
import subprocess
from pathlib import Path

from .schedules import ScheduleError

#: Default and shortest look interval.
WATCH_EVERY_S = 60
MIN_WATCH_S = 30
COMMAND_TIMEOUT_S = 20
OUTPUT_LIMIT = 4000

#: Reading commands that take no file paths: no argument may hold a "/".
NO_PATH_COMMANDS = frozenset({"squeue", "sacct", "sinfo", "sshare", "sprio", "qstat",
                              "bjobs", "df", "uptime", "free", "nvidia-smi", "date",
                              "whoami", "hostname"})
#: Reading commands whose arguments are files: each must lie in the folder.
PATH_COMMANDS = frozenset({"ls", "tail", "head", "wc", "grep", "du", "stat", "cat",
                           "md5sum", "sha256sum", "find"})
READING_COMMANDS = NO_PATH_COMMANDS | PATH_COMMANDS
#: Options refused because they read or write a file of their own choosing
#: or run something.
REFUSED_OPTIONS = {
    "grep": ("-f", "--file", "--exclude-from"),
    "find": ("-exec", "-execdir", "-ok", "-okdir", "-delete", "-fprint", "-fprintf",
             "-fls", "-fprint0"),
    "sacct": ("--output-file",),
    "du": ("--files0-from", "-X", "--exclude-from"),
    "wc": ("--files0-from",),
    "squeue": (),
}

LEAD_WORDS = frozenset({"wenn", "when", "sobald", "once", "nach", "after"})

_EVENT_DONE = r"(?:ist\s+)?(?:fertig|endet|beendet|vorbei|done|finishes|finished|ends|ended|exits|exited|is\s+done)"
_CHANGES = r"(?:sich\s+ändert|ändert\s+sich|sich\s+aendert|ändert|changes|changed)"
_APPEARS = r"(?:erscheint|existiert|da\s+ist|appears|exists|is\s+there)"

_LEAD = r"^\s*(?:wenn|when|sobald|once|nach|after)\s+"
_PATTERNS = (
    ("job", re.compile(_LEAD + r"(?:der\s+|the\s+)?job\s+(\S+?)(?:\s+" + _EVENT_DONE + r")?(?=[\s:;,]|$)", re.I)),
    ("pid", re.compile(_LEAD + r"(?:der\s+|the\s+)?(?:pid|prozess|process)\s+(\S+?)\s+" + _EVENT_DONE + r"(?=[\s:;,]|$)", re.I)),
    ("file", re.compile(_LEAD + r"(?:die\s+|the\s+)?(?:datei|file)\s+(\S+)\s+(" + _APPEARS + r"|" + _CHANGES + r")(?=[\s:;,]|$)", re.I)),
    ("output", re.compile(_LEAD + r"(?:die\s+|the\s+)?(?:ausgabe|output)(?:\s+(?:von|of))?\s+([\"'`])(.+?)\1\s+" + _CHANGES + r"(?=[\s:;,]|$)", re.I)),
)
_EVERY = re.compile(r"^\s*(?:,\s*)?(?:prüfe\s+|check\s+)?(?:alle|every)\s+(\d+)\s*([a-zäöü]+)(?=[\s:;,]|$)", re.I)


def split(text: str) -> tuple[dict, str]:
    """A watcher at the start of ``text`` and the task after it (may be
    empty: then the event is only reported)."""
    from .schedules import _UNIT
    raw = str(text or "").strip()
    for what, pattern in _PATTERNS:
        m = pattern.match(raw)
        if not m:
            continue
        if what == "output":
            target, event = m.group(2).strip(), "changes"
        elif what == "file":
            target = m.group(1)
            event = "changes" if re.fullmatch(_CHANGES, m.group(2), re.I) else "appears"
        else:
            target, event = m.group(1), "done"
        end = m.end()
        every = WATCH_EVERY_S
        e = _EVERY.match(raw[end:])
        if e and e.group(2).lower() in _UNIT:
            every = int(e.group(1)) * _UNIT[e.group(2).lower()]
            end += e.end()
            if every < MIN_WATCH_S:
                raise ScheduleError(f"a watcher looks at most every {MIN_WATCH_S} s",
                                    "interval")
        spec = {"kind": "watch", "what": what, "target": target, "event": event,
                "every": every, "text": raw[:end].strip().rstrip(",:;"), "tz": None}
        _validate(spec)
        rest = raw[end:].lstrip(" :;,-–—").strip()
        return spec, rest
    from .schedules import EXAMPLES
    raise ScheduleError(
        "no event found. Watchers: wenn job 123 fertig · when pid 4711 exits · "
        "wenn datei out.txt erscheint · when file calc.log changes · "
        f"wenn ausgabe \"squeue -u NAME\" sich ändert. {EXAMPLES}", "when")


def _validate(spec: dict) -> None:
    what, target = spec["what"], spec["target"]
    if what == "job" and not re.fullmatch(r"\d+(?:_\d+)?", target):
        raise ScheduleError(f"not a job number: {target}", "when")
    if what == "pid" and not target.isdigit():
        raise ScheduleError(f"not a process number: {target}", "when")
    if what == "output":
        argv = command_argv(target)
        name = argv[0]
        if name in NO_PATH_COMMANDS and any("/" in a for a in argv[1:]):
            raise ScheduleError(f"{name} takes no paths here", "path")


def command_argv(command: str) -> list[str]:
    """The command split without a shell, if it is one a watcher may run."""
    try:
        argv = shlex.split(command)
    except ValueError:
        raise ScheduleError("the command is not readable (quotes?)", "command") from None
    if not argv or argv[0] not in READING_COMMANDS:
        allowed = ", ".join(sorted(READING_COMMANDS))
        raise ScheduleError(f"a watcher runs only reading commands: {allowed}", "command")
    refused = REFUSED_OPTIONS.get(argv[0], ())
    for arg in argv[1:]:
        if any(arg == o or arg.startswith(o + "=") for o in refused):
            raise ScheduleError(f"{argv[0]} {arg} is not allowed in a watcher", "command")
    return argv


def _inside(workspace: str, path: str) -> Path:
    root = Path(os.path.realpath(workspace))
    target = Path(os.path.realpath(os.path.join(root, os.path.expanduser(path))))
    if target != root and root not in target.parents:
        raise ScheduleError(f"{path} is outside the schedule's folder", "path")
    return target


def _path_args(argv: list[str]) -> list[str]:
    args = [a for a in argv[1:] if not a.startswith("-")]
    if argv[0] == "grep" and args and not any(a in ("-e", "--regexp") or a.startswith("--regexp=")
                                              for a in argv[1:]):
        args = args[1:]                     # the pattern
    return args


def _command_paths(argv: list[str], workspace: str) -> None:
    """A file-reading command's paths, and the values of its options
    (``--opt=VALUE``), lie in ``workspace``."""
    for arg in _path_args(argv):
        _inside(workspace, arg)
    for arg in argv[1:]:
        if arg.startswith("-") and "=" in arg:
            _inside(workspace, arg.split("=", 1)[1])


def validate_paths(spec: dict, workspace: str) -> None:
    """Every path the watcher would touch lies in ``workspace``."""
    if spec.get("kind") != "watch":
        return
    if spec["what"] == "file":
        _inside(workspace, spec["target"])
    elif spec["what"] == "output":
        argv = command_argv(spec["target"])
        if argv[0] in PATH_COMMANDS:
            _command_paths(argv, workspace)


def _command_env() -> dict[str, str]:
    """What a watcher's command sees of the environment: nothing that could
    be a credential."""
    keep = ("PATH", "HOME", "USER", "LOGNAME", "LANG", "TZ", "SLURM_CONF")
    env = {k: os.environ[k] for k in keep if k in os.environ}
    env.setdefault("LANG", "C.UTF-8")
    return env


def _run(argv: list[str], cwd: str) -> str:
    try:
        done = subprocess.run(argv, cwd=cwd, env=_command_env(), stdin=subprocess.DEVNULL,
                              capture_output=True, text=True, errors="replace",
                              timeout=COMMAND_TIMEOUT_S)
    except FileNotFoundError:
        raise ScheduleError(f"{argv[0]} is not installed on this host", "command") from None
    except subprocess.TimeoutExpired:
        return f"(no answer within {COMMAND_TIMEOUT_S} s)"
    return ((done.stdout or "") + (done.stderr or ""))[:OUTPUT_LIMIT]


#: Asked through the throttled scheduler client (delfin.scheduler_client),
#: never directly: a cluster's operators count every scheduler query.
_SCHEDULER = frozenset({"squeue", "sacct", "sinfo", "sshare", "sprio", "scontrol"})
_ENDED = ("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
          "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE")


def _ask_scheduler(argv: list[str], cwd: str):
    """The scheduler's answer, through its throttle: the output, None when
    the command failed, or ``scheduler_client.THROTTLED`` (not asked now)."""
    import shutil
    from delfin import scheduler_client
    if shutil.which(argv[0]) is None:
        raise ScheduleError(f"{argv[0]} is not installed on this host", "command")

    def runner(cmd):
        try:
            done = subprocess.run(list(cmd), cwd=cwd, env=_command_env(),
                                  stdin=subprocess.DEVNULL, capture_output=True,
                                  text=True, errors="replace", timeout=COMMAND_TIMEOUT_S)
        except (OSError, subprocess.TimeoutExpired):
            return None
        return done.stdout if done.returncode == 0 else None
    return scheduler_client.query(scheduler_client.kind_of_argv(argv), argv, runner=runner)


def _job(target: str, workspace: str, previous: str | None) -> tuple[str, bool, str]:
    from delfin import scheduler_client
    queued = _ask_scheduler(["squeue", "-h", "-j", target, "-o", "%T"], workspace)
    if queued is scheduler_client.THROTTLED:
        return previous or "unknown", False, ""
    listed = [line.strip() for line in str(queued or "").splitlines() if line.strip()]
    if listed:
        return listed[0], False, ""
    # Not listed, or squeue failed (an id the queue has forgotten, or the
    # controller did not answer): the accounting says which.
    try:
        final = _ask_scheduler(["sacct", "-n", "-X", "-j", target, "-o", "State"], workspace)
    except ScheduleError:
        final = None
    states = str(final or "").split() if final is not scheduler_client.THROTTLED else []
    state = states[0] if states else ""
    if queued is None and not any(state.startswith(e) for e in _ENDED):
        return previous or "unknown", False, ""      # nobody could say it ended
    return "gone", True, f"job {target} has finished ({state or 'no longer in the queue'})"


def check(spec: dict, workspace: str, previous: str | None) -> tuple[str, bool, str]:
    """One look: ``(state, fired, note)``. Raises ScheduleError when the
    watcher cannot work here (a path that left the folder, a command that
    is missing)."""
    what, target = spec["what"], spec["target"]
    if what == "job":
        return _job(target, workspace, previous)
    if what == "pid":
        try:
            os.kill(int(target), 0)
        except ProcessLookupError:
            return "gone", True, f"process {target} has ended"
        except PermissionError:
            pass
        return "running", False, ""
    if what == "file":
        path = _inside(workspace, target)
        if spec["event"] == "appears":
            return ("there", True, f"{target} appeared") if path.exists() else ("missing", False, "")
        try:
            st = path.stat()
            state = f"{st.st_mtime_ns}:{st.st_size}"
        except OSError:
            state = "missing"
        if previous is None or previous == state:
            return state, False, ""
        size = state.split(":")[1] if ":" in state else "gone"
        return state, True, f"{target} changed (size {size})"
    argv = command_argv(target)
    if argv[0] in PATH_COMMANDS:
        _command_paths(argv, workspace)
    if argv[0] in _SCHEDULER:
        from delfin import scheduler_client
        answer = _ask_scheduler(argv, workspace)
        if answer is scheduler_client.THROTTLED:
            # Not asked this look. Return the previous state unchanged
            # (None stays None): storing "" here would make the next real
            # look fire "output changed" on output that did not change,
            # because "" is never equal to an output hash.
            return previous, False, ""
        out = str(answer or "")[:OUTPUT_LIMIT]
    else:
        out = _run(argv, workspace)
    state = hashlib.sha256(out.encode("utf-8", "replace")).hexdigest()
    if previous is None or previous == state:
        return state, False, ""
    return state, True, f"the output of `{target}` changed:\n{out[-1500:]}"


__all__ = ["split", "check", "validate_paths", "command_argv", "READING_COMMANDS",
           "WATCH_EVERY_S", "MIN_WATCH_S"]
