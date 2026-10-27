"""User schedules: a task that runs at set times, in the agent service.

Input: a schedule the user writes -- in a session (``/schedule``), on the
command line (``delfin-agent schedule``) or from the paired phone -- as the
time, then the task::

    täglich 07:30 Fasse die Slurm-Jobs zusammen
    mo-fr 18:00 tz=Europe/Berlin Prüfe die Rechnungen
    alle 2h Schau nach neuen Ergebnissen
    30 7 * * 1-5 Lies die neuen Logs

Output: an entry in ``~/.delfin/user_schedules/schedules.json`` and, for
every run, a line in ``runs.jsonl`` beside it (what ran, how it ended, its
cost and session), which the TUI's attention list and the phone's
schedules chat report.

Semantics.

* The time is kept as written plus the zone it is meant in (the host's
  named zone when none is given), and the next run is computed in that
  zone: 07:30 stays 07:30 across a change of summer and winter time. A
  time that falls into the skipped hour runs once, an hour later that day;
  a time inside the repeated hour runs once, at its first occurrence.
* Day-of-month and weekday in a cron line widen each other when both are
  restricted, as in cron.
* An interval (``alle 2h``) counts from the last run. The shortest is five
  minutes: every run is an agent turn and costs.
* A schedule runs only on the host it was made on (``host``): on a cluster
  the home directory is shared by every login node, and each runs its own
  agent service.
* A schedule runs only when the user confirmed exactly its content: the
  hash of what it does, where and when (``fingerprint``) is signed at
  confirmation with the account key (``trust_store``, HMAC; no process
  carrying the agent-session marker gets that key), and a change made
  since -- by anyone, through any route -- withdraws it. Where no account
  key can be had at all (a refused key file), the plain hash stands in,
  so the user's schedules still run. One the user writes is confirmed by writing it; one the
  agent proposes waits for ``/schedule ok N`` or the phone's ``ok N``.
* A missed run (the host was off) is made up once if it is less than
  ``CATCH_UP_S`` late, and skipped with a record otherwise; missed runs are
  never stacked. A run never starts while the schedule's previous one
  still runs.
* Each run has a cost ceiling and a time ceiling, and each schedule a daily
  cost ceiling; reached, the run is not started and the skip is recorded.

The store is the user's alone: a 0700 directory with 0600 files, masked in
every cage an agent's command runs in and refused to the file tools, like
the messenger's activation records (``sandbox._HOME_SECRET_DIRS``,
``api_client._DEFAULT_PATH_DENY_GLOBS``). The command line refuses to
change it from inside an agent session.

Rejected: keeping these in ``cron.json`` beside the agent's own wake-ups.
A build from before this one reads that file and, meeting a kind it does
not know, fires the entry and deletes it -- and on a cluster every login
node may run a different build.
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import re
import time
from contextlib import contextmanager
from datetime import date, datetime, timedelta, timezone, tzinfo
from pathlib import Path
from typing import Iterator

_DIR = Path.home() / ".delfin" / "user_schedules"

#: The shortest interval: every run is an agent turn.
MIN_INTERVAL_S = 300
#: A missed run later than this is skipped, not made up.
CATCH_UP_S = 6 * 3600.0
#: Defaults a schedule starts with; each is changeable per schedule.
DEFAULT_RUN_USD = 1.0
DEFAULT_DAY_USD = 2.0
DEFAULT_MAX_MINUTES = 30
NOTIFY = ("always", "errors", "never")
#: Runs kept in runs.jsonl.
KEEP_RUNS = 500
#: The longest task text.
MAX_TASK = 4000

EXAMPLES = ("Examples of when: täglich 07:30 · mo-fr 18:00 · alle 2h · "
            "30 7 * * 1-5 · daily 7:30 tz=Europe/Berlin")


class ScheduleError(ValueError):
    """A schedule that cannot be made or changed; ``code`` names why for a
    surface that words it itself (the phone): when, time, interval, zone,
    never, task, workspace, unknown, host, budget."""

    def __init__(self, message: str, code: str = "when"):
        super().__init__(message)
        self.code = code


# ---------------------------------------------------------------------------
# Time zones
# ---------------------------------------------------------------------------

def available_zone_names() -> set[str]:
    try:
        from zoneinfo import available_timezones
        return set(available_timezones())
    except Exception:
        return set()


def zone(name: str) -> tzinfo:
    """The tzinfo for a zone name; UTC needs no zone database."""
    if name in ("UTC", "Etc/UTC", "utc"):
        return timezone.utc
    from zoneinfo import ZoneInfo
    try:
        return ZoneInfo(name)
    except Exception:
        raise ScheduleError(f"unknown time zone: {name}", "zone") from None


def _known_zone(name: str) -> bool:
    try:
        zone(name)
        return True
    except ScheduleError:
        return False


def local_zone() -> str:
    """The host's zone by name (TZ, /etc/localtime, /etc/timezone); UTC
    when none can be named -- an offset alone does not know summer time."""
    candidates = [os.environ.get("TZ", "").lstrip(":")]
    try:
        link = os.readlink("/etc/localtime")
        if "zoneinfo/" in link:
            candidates.append(link.split("zoneinfo/", 1)[1])
    except OSError:
        pass
    try:
        candidates.append(Path("/etc/timezone").read_text().strip())
    except OSError:
        pass
    for name in candidates:
        if name in ("UTC", "Etc/UTC"):
            return "UTC"
        if name and "/" in name and _known_zone(name):
            return name
    return "UTC"


# ---------------------------------------------------------------------------
# Reading the time
# ---------------------------------------------------------------------------

_DAY_WORDS = {
    0: ("mo", "mon", "montag", "montags", "monday", "mondays"),
    1: ("di", "tue", "tues", "dienstag", "dienstags", "tuesday", "tuesdays"),
    2: ("mi", "wed", "mittwoch", "mittwochs", "wednesday", "wednesdays"),
    3: ("do", "thu", "thur", "thurs", "donnerstag", "donnerstags", "thursday",
        "thursdays"),
    4: ("fr", "fri", "freitag", "freitags", "friday", "fridays"),
    5: ("sa", "sat", "samstag", "samstags", "saturday", "saturdays"),
    6: ("so", "sun", "sonntag", "sonntags", "sunday", "sundays"),
}
_DAY = {word: day for day, words in _DAY_WORDS.items() for word in words}
_DAY_RE = "|".join(sorted(_DAY, key=len, reverse=True))
_DAYS_TOKEN = re.compile(rf"^(?:{_DAY_RE})\.?(?:-(?:{_DAY_RE})\.?)?"
                         rf"(?:,(?:{_DAY_RE})\.?(?:-(?:{_DAY_RE})\.?)?)*$")
_EVERY_DAY = {"täglich", "taeglich", "tägl", "tägl.", "daily"}
_WEEKDAYS = {"werktags", "wochentags", "weekdays", "weekday"}
_WEEKEND = {"wochenende", "wochenendes", "weekend", "weekends"}
_FILLER = {"jeden", "jede", "jedem", "every", "each", "am", "on", "um", "at",
           "und", "and", "uhr", "o'clock", "an", "den", "der", "die"}
_HOURLY = {"stündlich", "stuendlich", "hourly"}
_UNIT = {"s": 1, "sek": 1, "sec": 1,
         "m": 60, "min": 60, "mins": 60, "minute": 60, "minuten": 60,
         "minutes": 60,
         "h": 3600, "std": 3600, "stunde": 3600, "stunden": 3600,
         "hour": 3600, "hours": 3600,
         "d": 86400, "tag": 86400, "tage": 86400, "tagen": 86400,
         "day": 86400, "days": 86400}
_CRON_FIELD = re.compile(r"^[\d*/,\-]+$")
_TIME = re.compile(r"^(\d{1,2})(?::(\d{2}))?$")


def _days(token: str) -> list[int]:
    out: set[int] = set()
    for part in token.split(","):
        first, _, last = part.partition("-")
        a = _DAY[first.rstrip(".")]
        if not last:
            out.add(a)
            continue
        b = _DAY[last.rstrip(".")]
        d = a
        while True:
            out.add(d)
            if d == b:
                break
            d = (d + 1) % 7
    return sorted(out)


def _cron_field(text: str, low: int, high: int) -> list[int]:
    values: set[int] = set()
    for part in text.split(","):
        rng, _, step_s = part.partition("/")
        step = int(step_s) if step_s else 1
        if step < 1:
            raise ScheduleError(f"cron step out of range: {part}. {EXAMPLES}")
        if rng == "*":
            a, b = low, high
        elif "-" in rng:
            a_s, _, b_s = rng.partition("-")
            a, b = int(a_s), int(b_s)
        else:
            a = int(rng)
            b = high if step_s else a
        if not (low <= a <= high and low <= b <= high and a <= b):
            raise ScheduleError(f"cron field out of range: {part}. {EXAMPLES}")
        values.update(range(a, b + 1, step))
    return sorted(values)


def _cron(fields: list[str]) -> dict:
    try:
        minutes = _cron_field(fields[0], 0, 59)
        hours = _cron_field(fields[1], 0, 23)
        dom = _cron_field(fields[2], 1, 31)
        months = _cron_field(fields[3], 1, 12)
        dow_cron = _cron_field(fields[4], 0, 7)
    except ValueError as exc:
        if isinstance(exc, ScheduleError):
            raise
        raise ScheduleError(f"cannot read the cron line. {EXAMPLES}") from None
    # cron: 0 and 7 are Sunday; Python's weekday: Monday 0 .. Sunday 6.
    dow = sorted({(d - 1) % 7 for d in dow_cron})
    return {"kind": "cron", "minutes": minutes, "hours": hours, "dom": dom,
            "months": months, "dow": dow,
            "domAny": fields[2] == "*", "dowAny": fields[4] == "*"}


def split_when(text: str) -> tuple[dict, str]:
    """The schedule at the start of ``text`` and the task after it.

    Returns ``(spec, rest)``; ``spec["text"]`` is the time as written,
    ``spec["tz"]`` the zone when one was named (None otherwise). Raises
    ScheduleError when no time is found or it cannot be a schedule."""
    raw = str(text or "").strip()
    tokens = raw.split()
    tz: str | None = None
    # A cron line: five fields of digits, *, /, - and commas.
    if len(tokens) >= 5 and all(_CRON_FIELD.match(t) for t in tokens[:5]):
        spec = _cron(tokens[:5])
        i = 5
        i, tz = _zone_token(tokens, i)
        spec.update(text=" ".join(tokens[:5]), tz=tz)
        return spec, _rest(tokens, i)
    i = 0
    days: list[int] | None = None
    times: list[list[int]] = []
    every = 0
    said_every = False
    pending_at = False
    while i < len(tokens):
        word = tokens[i].lower()
        bare = word.rstrip(",:;")
        if every or times:
            # The time is read: from here on only more of it -- another
            # time, "Uhr", a zone, a day word -- and the task starts at the
            # first other word ("so", "an", "die" are task words here).
            if not (_times(bare) is not None or bare in ("uhr", "o'clock")
                    or bare.startswith("tz=") or bare.upper() == "UTC"
                    or ("/" in tokens[i] and _known_zone(tokens[i].rstrip(",:;")))
                    or (not every and (bare in _EVERY_DAY | _WEEKDAYS | _WEEKEND))
                    or (not every and bare in ("und", "and") and i + 1 < len(tokens)
                        and _times(tokens[i + 1].lower().rstrip(",:;")))):
                break
            if every and _times(bare) is not None:
                break
        if bare in ("alle", "every", "each") and i + 1 < len(tokens):
            n, used = _interval(tokens, i + 1)
            if n:
                every, i = n, i + 1 + used
                continue
        if bare in _HOURLY:
            every, i = 3600, i + 1
            continue
        if bare in _EVERY_DAY:
            days, i = list(range(7)), i + 1
            continue
        if bare in ("tag", "day") and said_every:
            days, i = list(range(7)), i + 1
            continue
        if bare in _WEEKDAYS:
            days, i = [0, 1, 2, 3, 4], i + 1
            continue
        if bare in _WEEKEND:
            days, i = [5, 6], i + 1
            continue
        if _DAYS_TOKEN.match(bare):
            days = sorted(set(days or []) | set(_days(bare)))
            i += 1
            continue
        found = _times(bare)
        if found is not None and (":" in bare or pending_at or _uhr_follows(tokens, i)):
            if not found:
                raise ScheduleError(f"{tokens[i].rstrip(',:;')} is not a time of "
                                    f"day. {EXAMPLES}", "time")
            times.extend(found)
            i += 1
            pending_at = False
            if word.endswith(":") and not word.endswith("::"):
                break                       # "18:00: task" -- the colon ends the time
            continue
        if bare in _FILLER:
            said_every = said_every or bare in ("jeden", "jede", "every", "each")
            pending_at = bare in ("um", "at")
            i += 1
            continue
        if "=" in bare and bare.startswith("tz="):
            name = tokens[i].split("=", 1)[1]
            if not _known_zone(name):
                raise ScheduleError(f"unknown time zone: {name}", "zone")
            tz, i = ("UTC" if name.upper() in ("UTC", "ETC/UTC") else name), i + 1
            continue
        if bare.upper() == "UTC" or ("/" in tokens[i] and _known_zone(tokens[i].rstrip(",:;"))):
            tz = "UTC" if bare.upper() == "UTC" else tokens[i].rstrip(",:;")
            i += 1
            continue
        break
    if every:
        if every < MIN_INTERVAL_S:
            raise ScheduleError("the shortest interval is 5 minutes: every run "
                                "is an agent turn", "interval")
        return {"kind": "every", "every": every, "text": _said(tokens, i),
                "tz": tz}, _rest(tokens, i)
    if times:
        spec = {"kind": "at", "days": days if days is not None else list(range(7)),
                "times": sorted(times), "text": _said(tokens, i), "tz": tz}
        return spec, _rest(tokens, i)
    if days is not None:
        raise ScheduleError(f"a time of day is missing (HH:MM). {EXAMPLES}", "time")
    raise ScheduleError(f"no schedule found at the start. {EXAMPLES}", "when")


def _uhr_follows(tokens: list[str], i: int) -> bool:
    return i + 1 < len(tokens) and tokens[i + 1].lower().rstrip(",:;") in ("uhr", "o'clock")


def _times(token: str) -> list[list[int]] | None:
    """[[h, m], ...] for "7:30" or "7:30,18:00"; [] when it is shaped like
    a time and is none (25:00); None when it is not shaped like one."""
    out: list[list[int]] = []
    for part in token.split(","):
        if not part:
            continue
        m = _TIME.match(part)
        if not m:
            return None
        h, mi = int(m.group(1)), int(m.group(2) or 0)
        if h > 23 or mi > 59:
            return []
        out.append([h, mi])
    return out or None


def _interval(tokens: list[str], i: int) -> tuple[int, int]:
    """(seconds, tokens used) for "2h", "2 h", "15 min", "3 Stunden"."""
    tok = tokens[i].lower().rstrip(",:;")
    m = re.match(r"^(\d+)\s*([a-zäöü]+)$", tok)
    if m and m.group(2) in _UNIT:
        return int(m.group(1)) * _UNIT[m.group(2)], 1
    if tok.isdigit() and i + 1 < len(tokens):
        unit = tokens[i + 1].lower().rstrip(",:;.")
        if unit in _UNIT:
            return int(tok) * _UNIT[unit], 2
    return 0, 0


def _zone_token(tokens: list[str], i: int) -> tuple[int, str | None]:
    if i < len(tokens):
        tok = tokens[i]
        if tok.lower().startswith("tz="):
            name = tok.split("=", 1)[1]
            if not _known_zone(name):
                raise ScheduleError(f"unknown time zone: {name}", "zone")
            return i + 1, "UTC" if name.upper() in ("UTC", "ETC/UTC") else name
        if tok.upper() == "UTC" or ("/" in tok and _known_zone(tok)):
            return i + 1, "UTC" if tok.upper() == "UTC" else tok
    return i, None


def _said(tokens: list[str], i: int) -> str:
    return " ".join(tokens[:i]).rstrip(",:;")


def _rest(tokens: list[str], i: int) -> str:
    return " ".join(tokens[i:]).lstrip(":-–— ").strip()


def parse_when(text: str) -> dict:
    """A time with nothing after it."""
    spec, rest = split_when(text)
    if rest:
        raise ScheduleError(f"not part of a time: {rest[:40]}. {EXAMPLES}", "when")
    return spec


# ---------------------------------------------------------------------------
# The next run
# ---------------------------------------------------------------------------

#: How far ahead a run is looked for: 29 February comes back within 8 years.
_HORIZON_DAYS = 8 * 366


def _day_matches(spec: dict, d: date) -> bool:
    if spec["kind"] == "at":
        return d.weekday() in spec["days"]
    if d.month not in spec["months"]:
        return False
    dom_ok = d.day in spec["dom"]
    dow_ok = d.weekday() in spec["dow"]
    if spec["domAny"] and spec["dowAny"]:
        return True
    if spec["domAny"]:
        return dow_ok
    if spec["dowAny"]:
        return dom_ok
    return dom_ok or dow_ok


def _clock_times(spec: dict) -> list[tuple[int, int]]:
    if spec["kind"] == "at":
        return [(h, m) for h, m in spec["times"]]
    return [(h, m) for h in spec["hours"] for m in spec["minutes"]]


def next_fire(spec: dict, after: float, tz_name: str) -> float:
    """The first run strictly after ``after`` (epoch seconds), in
    ``tz_name`` unless the spec names its own zone."""
    if spec["kind"] == "every":
        return float(after) + int(spec["every"])
    tzi = zone(spec.get("tz") or tz_name)
    start = datetime.fromtimestamp(after, tzi).date()
    clock = _clock_times(spec)
    for offset in range(_HORIZON_DAYS):
        d = start + timedelta(days=offset)
        if not _day_matches(spec, d):
            continue
        # Every wall-clock time of the day once, at its first occurrence; a
        # time inside a skipped hour maps to the instant after the jump. The
        # day's earliest instant wins: the mapping can reorder them.
        best = None
        for h, m in clock:
            ts = datetime(d.year, d.month, d.day, h, m, tzinfo=tzi, fold=0).timestamp()
            if ts > after and (best is None or ts < best):
                best = ts
        if best is not None:
            return best
    raise ScheduleError("this schedule never fires (no such date)", "never")


# ---------------------------------------------------------------------------
# The store
# ---------------------------------------------------------------------------

def _file() -> Path:
    return _DIR / "schedules.json"


def _runs_file() -> Path:
    return _DIR / "runs.jsonl"


def _ensure_dir() -> Path:
    _DIR.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        os.chmod(_DIR, 0o700)
    except OSError:
        pass
    return _DIR


@contextmanager
def _locked() -> Iterator[None]:
    _ensure_dir()
    fd = os.open(_DIR / ".lock", os.O_CREAT | os.O_RDWR, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        os.close(fd)


def _read() -> dict:
    try:
        data = json.loads(_file().read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {"next_id": 1, "entries": {}}
    if not isinstance(data, dict) or not isinstance(data.get("entries"), dict):
        return {"next_id": 1, "entries": {}}
    data.setdefault("next_id", 1)
    return data


def _write(data: dict) -> None:
    from .state_paths import write_text_atomic
    write_text_atomic(_file(), json.dumps(data, indent=1, ensure_ascii=False))


def fingerprint(entry: dict) -> str:
    """What the user confirms: what runs, where, when and at what cost."""
    keep = [entry.get(k) for k in ("task", "when", "tz", "workspace", "host",
                                   "budget_run_usd", "budget_day_usd",
                                   "max_minutes")]
    return hashlib.sha256(json.dumps(keep, ensure_ascii=False).encode()).hexdigest()


_MAC_PURPOSE = "user_schedules"


def _seal(entry: dict) -> str:
    """The confirmation of ``entry`` as it is now. Refused inside an agent
    session: only the user confirms."""
    from . import trust_store
    from .approval_guard import inside_agent_session
    if inside_agent_session():
        raise ScheduleError("refused inside an agent session: only the user "
                            "confirms a schedule", "unknown")
    tag = trust_store.mac(_MAC_PURPOSE, fingerprint(entry).encode())
    return tag or "sha:" + fingerprint(entry)


def confirmed(entry: dict) -> bool:
    """Whether the user confirmed exactly this content: the HMAC under the
    account key, or -- only where there is no account key at all -- the
    plain hash."""
    tag = str(entry.get("confirmed_sha") or "")
    if not tag:
        return False
    from . import trust_store
    if tag.startswith("sha:"):
        return (tag == "sha:" + fingerprint(entry)
                and trust_store.account_key(create=True) is None)
    return trust_store.check_mac(_MAC_PURPOSE, fingerprint(entry).encode(), tag) is True


def _spec(entry: dict) -> dict:
    spec = parse_when(entry["when"])
    spec["tz"] = entry.get("tz") or spec.get("tz")
    return spec


def _hostname() -> str:
    from .service import rundir
    return rundir.hostname()


def entries() -> list[dict]:
    """Every schedule, oldest first, each with ``id`` and ``confirmed``."""
    data = _read()
    out = []
    for key, entry in data["entries"].items():
        if isinstance(entry, dict):
            out.append({**entry, "id": int(key), "confirmed": confirmed(entry)})
    return sorted(out, key=lambda e: e["id"])


def get(entry_id: int | str) -> dict:
    data = _read()
    entry = data["entries"].get(str(_id(entry_id)))
    if not isinstance(entry, dict):
        raise ScheduleError(f"no schedule {entry_id}", "unknown")
    return {**entry, "id": int(_id(entry_id)), "confirmed": confirmed(entry)}


def _id(entry_id: int | str) -> int:
    try:
        return int(str(entry_id).lstrip("#"))
    except ValueError:
        raise ScheduleError(f"no schedule {entry_id}", "unknown") from None


def add(text: str, *, workspace: str, created_by: str, confirm: bool,
        host: str | None = None, tz: str | None = None,
        budget_run_usd: float = DEFAULT_RUN_USD,
        budget_day_usd: float = DEFAULT_DAY_USD,
        max_minutes: int = DEFAULT_MAX_MINUTES, notify: str = "always",
        now: float | None = None) -> dict:
    """A new schedule from "<when> <task>". ``confirm``: the user wrote
    it (a session, the command line, the phone after its ``ok``); False
    for a proposal, which waits for the user."""
    spec, task = split_when(text)
    task = task.strip()
    if not task:
        raise ScheduleError("the task is missing after the time", "task")
    if len(task) > MAX_TASK:
        raise ScheduleError(f"the task is longer than {MAX_TASK} characters", "task")
    ws = Path(str(workspace or "")).expanduser()
    if not ws.is_absolute() or not ws.is_dir():
        raise ScheduleError(f"no such folder: {workspace}", "workspace")
    if notify not in NOTIFY:
        raise ScheduleError(f"notify is one of {', '.join(NOTIFY)}", "budget")
    run_usd, day_usd = float(budget_run_usd), float(budget_day_usd)
    if not (0 < run_usd <= day_usd) or day_usd > 1000:
        raise ScheduleError("a run's ceiling must be above 0 and at most the "
                            "day's", "budget")
    if not 1 <= int(max_minutes) <= 24 * 60:
        raise ScheduleError("a run may take 1 to 1440 minutes", "budget")
    zone_name = spec.get("tz") or tz or local_zone()
    zone(zone_name)
    now = time.time() if now is None else now
    first = next_fire(spec, now, zone_name)
    entry = {
        "task": task, "when": spec["text"], "tz": zone_name,
        "workspace": str(ws.resolve()), "host": host or _hostname(),
        "created_at": now, "created_by": str(created_by)[:20],
        "budget_run_usd": run_usd, "budget_day_usd": day_usd,
        "max_minutes": int(max_minutes), "notify": notify,
        "paused": False, "next_fire_at": first, "last_fired_at": 0.0,
        "last_status": "", "last_session": "", "running": "",
        "run_requested": False, "spent": {}, "confirmed_sha": "",
    }
    if confirm:
        entry["confirmed_sha"] = _seal(entry)
    with _locked():
        data = _read()
        new_id = int(data.get("next_id") or 1)
        while str(new_id) in data["entries"]:
            new_id += 1
        data["entries"][str(new_id)] = entry
        data["next_id"] = new_id + 1
        _write(data)
    return {**entry, "id": new_id, "confirmed": confirmed(entry)}


def _change(entry_id: int | str, fn) -> dict:
    with _locked():
        data = _read()
        key = str(_id(entry_id))
        entry = data["entries"].get(key)
        if not isinstance(entry, dict):
            raise ScheduleError(f"no schedule {entry_id}", "unknown")
        fn(entry)
        _write(data)
        return {**entry, "id": int(key), "confirmed": confirmed(entry)}


def confirm(entry_id: int | str, *, expect: str = "",
            now: float | None = None) -> dict:
    """The user confirms this content. ``expect``: the fingerprint the user
    was shown; a different one (changed since) is refused. A first time
    that went by while it waited is not made up: it looks ahead."""
    def apply(entry: dict) -> None:
        sha = fingerprint(entry)
        if expect and expect != sha:
            raise ScheduleError("the schedule changed since it was shown; "
                                "look at it again", "unknown")
        entry["confirmed_sha"] = _seal(entry)
        moment = time.time() if now is None else now
        if float(entry.get("next_fire_at") or 0) < moment:
            entry["next_fire_at"] = next_fire(_spec(entry), moment, entry["tz"])
    return _change(entry_id, apply)


def pause(entry_id: int | str) -> dict:
    return _change(entry_id, lambda e: e.update(paused=True))


def resume(entry_id: int | str, *, now: float | None = None) -> dict:
    def apply(entry: dict) -> None:
        entry["paused"] = False
        moment = time.time() if now is None else now
        if float(entry.get("next_fire_at") or 0) < moment:
            entry["next_fire_at"] = next_fire(_spec(entry), moment, entry["tz"])
    return _change(entry_id, apply)


def request_run(entry_id: int | str) -> dict:
    return _change(entry_id, lambda e: e.update(run_requested=True))


def delete(entry_id: int | str) -> dict:
    with _locked():
        data = _read()
        key = str(_id(entry_id))
        entry = data["entries"].pop(key, None)
        if not isinstance(entry, dict):
            raise ScheduleError(f"no schedule {entry_id}", "unknown")
        _write(data)
        return {**entry, "id": int(key)}


def active_here(host: str | None = None) -> bool:
    """Whether this host has a schedule that will run: the agent service
    stays up for it."""
    host = host or _hostname()
    return any(e.get("host") == host and e["confirmed"] and not e.get("paused")
               for e in entries())


# ---------------------------------------------------------------------------
# What the agent service does with them
# ---------------------------------------------------------------------------

def _day_key(ts: float, tz_name: str) -> str:
    return datetime.fromtimestamp(ts, zone(tz_name)).strftime("%Y-%m-%d")


def spent_today(entry: dict, now: float | None = None) -> float:
    now = time.time() if now is None else now
    return float((entry.get("spent") or {}).get(_day_key(now, entry["tz"]), 0.0))


def take_due(host: str, now: float | None = None) -> list[dict]:
    """The schedules of ``host`` to start now, marked as starting; and the
    skips decided on the way (missed, over the day's ceiling, changed),
    each recorded as a run.

    Called by the agent service only. A schedule is due when it is
    confirmed, not paused, not running, and asked to run or past its time;
    its next time is set at once, so a slow start never fires it twice."""
    now = time.time() if now is None else now
    due: list[dict] = []
    skipped: list[tuple[dict, str, str]] = []
    with _locked():
        data = _read()
        changed = False
        for key, entry in data["entries"].items():
            if not isinstance(entry, dict) or entry.get("host") != host:
                continue
            if entry.get("paused") or entry.get("running"):
                continue
            asked = bool(entry.get("run_requested"))
            at = float(entry.get("next_fire_at") or 0)
            if not asked and now < at:
                continue
            if not confirmed(entry):
                if asked:
                    entry["run_requested"] = False
                    changed = True
                    skipped.append(({**entry, "id": int(key)}, "unconfirmed",
                                     "not confirmed: changed since confirmed, or "
                                     "proposed and never confirmed"))
                continue
            try:
                spec = _spec(entry)
                entry["next_fire_at"] = next_fire(spec, now, entry["tz"]) \
                    if spec["kind"] != "every" else now + spec["every"]
            except ScheduleError as exc:
                entry["paused"] = True
                changed = True
                skipped.append(({**entry, "id": int(key)}, "failed", str(exc)))
                continue
            entry["run_requested"] = False
            changed = True
            if not asked and now - at > CATCH_UP_S:
                skipped.append(({**entry, "id": int(key)}, "missed",
                                f"missed by {int((now - at) // 60)} min (host off "
                                "or the agent service not running)"))
                continue
            if spent_today(entry, now) >= float(entry.get("budget_day_usd") or 0):
                skipped.append(({**entry, "id": int(key)}, "budget",
                                f"today's ceiling of ${entry['budget_day_usd']:.2f} "
                                "is reached"))
                continue
            entry["running"] = "starting"
            entry["last_fired_at"] = now
            due.append({**entry, "id": int(key)})
        if changed:
            _write(data)
    for entry, status, why in skipped:
        record_run(entry, status=status, text=why, now=now)
    return due


def started(entry_id: int, session_id: str) -> None:
    try:
        _change(entry_id, lambda e: e.update(running=session_id or "starting",
                                             last_session=session_id or ""))
    except ScheduleError:
        pass


def finished(entry: dict, *, status: str, text: str = "", cost: float = 0.0,
             session_id: str = "", now: float | None = None) -> dict:
    """A run ended: its record, the schedule's last state and the day's
    spend. ``status``: done, failed, timeout, waiting (an approval held
    it until its time ran out)."""
    now = time.time() if now is None else now

    def apply(e: dict) -> None:
        e["running"] = ""
        e["last_status"] = status
        if session_id:
            e["last_session"] = session_id
        day = _day_key(now, e["tz"])
        spent = {k: v for k, v in (e.get("spent") or {}).items() if k >= _day_key(now - 8 * 86400, e["tz"])}
        spent[day] = round(float(spent.get(day, 0.0)) + max(0.0, float(cost)), 6)
        e["spent"] = spent
    try:
        _change(entry["id"], apply)
    except ScheduleError:
        pass                                  # deleted while it ran: the run is still recorded
    return record_run(entry, status=status, text=text, cost=cost,
                      session_id=session_id, now=now)


def release_all(host: str) -> None:
    """At the service's start: a run the previous service held is over."""
    with _locked():
        data = _read()
        stale = [e for e in data["entries"].values()
                 if isinstance(e, dict) and e.get("host") == host and e.get("running")]
        for e in stale:
            e["running"] = ""
            e["last_status"] = "interrupted"
        if stale:
            _write(data)


def record_run(entry: dict, *, status: str, text: str = "", cost: float = 0.0,
               session_id: str = "", now: float | None = None) -> dict:
    now = time.time() if now is None else now
    with _locked():
        lines = _read_runs()
        seq = (lines[-1]["seq"] + 1) if lines else 1
        record = {"seq": seq, "ts": round(now, 3), "id": int(entry["id"]),
                  "task": str(entry.get("task") or "")[:200],
                  "when": str(entry.get("when") or ""), "host": entry.get("host", ""),
                  "status": status, "text": str(text or "")[:4000],
                  "costUsd": round(float(cost or 0.0), 6), "session": session_id,
                  "notify": entry.get("notify", "always")}
        lines.append(record)
        from .state_paths import write_text_atomic
        write_text_atomic(_runs_file(), "".join(
            json.dumps(r, ensure_ascii=False) + "\n" for r in lines[-KEEP_RUNS:]))
    return record


def _read_runs() -> list[dict]:
    out = []
    try:
        text = _runs_file().read_text(encoding="utf-8")
    except OSError:
        return out
    for line in text.splitlines():
        try:
            row = json.loads(line)
        except ValueError:
            continue
        if isinstance(row, dict) and isinstance(row.get("seq"), int):
            out.append(row)
    return out


def runs(after_seq: int = 0, limit: int = 50) -> list[dict]:
    """Runs recorded after ``after_seq``, oldest first."""
    rows = [r for r in _read_runs() if r["seq"] > after_seq]
    return rows[-limit:]


def last_seq() -> int:
    rows = _read_runs()
    return rows[-1]["seq"] if rows else 0


def describe_next(entry: dict, now: float | None = None) -> str:
    """The next run as the schedule's own wall clock shows it."""
    ts = float(entry.get("next_fire_at") or 0)
    if not ts:
        return "-"
    return datetime.fromtimestamp(ts, zone(entry["tz"])).strftime("%a %d.%m. %H:%M")


__all__ = [
    "ScheduleError", "split_when", "parse_when", "next_fire", "local_zone",
    "available_zone_names", "zone", "entries", "get", "add", "confirm",
    "pause", "resume", "request_run", "delete", "active_here", "take_due",
    "started", "finished", "release_all", "record_run", "runs", "last_seq",
    "fingerprint", "confirmed", "spent_today", "describe_next",
    "MIN_INTERVAL_S", "CATCH_UP_S", "EXAMPLES",
]
