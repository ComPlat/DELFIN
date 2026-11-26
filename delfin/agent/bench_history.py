"""The benchmark history: every round kept, and the trend read back.

:mod:`bench_nightly` runs a graded round and hands back a ``RunRecord``,
and that was the end of it — nothing stored the record, so a week later
no one could say whether DELFIN got better. A round that is not kept
cannot show an improvement. This module is the store and the read side.

Two files, append-only, under one root (``agent.bench_history_root``,
default ``~/.delfin/bench_history`` — never inside a repository, which
:func:`root` refuses):

    rounds.jsonl      one line per round: round id, commit, model,
                      provider, seed mode, task-set version, budget
    instances.jsonl   one line per drawn instance: task id, capability
                      tag, pass/fail/score, unmeasured, duration, tokens

Append-only in the strong sense: a line is added with a single ``write``
to a descriptor opened ``O_APPEND``, and nothing here opens a file for
writing any other way, so an existing byte is never rewritten and a
concurrent round cannot interleave a partial line (a write smaller than
PIPE_BUF to an O_APPEND descriptor is atomic, and one JSON line of a
record is far smaller). ``round_id`` is therefore not a key — two rounds
may legitimately share one id, and both lines are kept and shown.

A record that cannot be written RAISES (:class:`HistoryWriteError`). It
is never swallowed: a round that ran and was not kept is the exact
silent failure this module exists to end, and a caller must find out
that the round it just paid for is gone.

Comparability
-------------
A pass rate means something only across rounds drawn the same way, of
the same tasks, of the same model. :func:`comparable` names the three
axes (seed mode, task-set version, model), and the read side drops every
round that differs from the reference on any of them, reporting the
dropped round with the axis that excluded it — so "the trend is flat"
can never quietly mean "three of my five rounds were incomparable".
Unmeasured rounds (every started instance unmeasured: the run never
happened) are excluded too, and say so: scoring one fabricates a
measurement.
"""

from __future__ import annotations

import json
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

from . import state_paths

#: Every field the read side can exclude a round on. A round that
#: differs from the reference on any of these is not the same
#: measurement, and a rate across two of them answers another question.
COMPARABLE_AXES = ("seed_mode", "task_set_version", "model")

ROUNDS_FILE = "rounds.jsonl"
INSTANCES_FILE = "instances.jsonl"

#: Longest setting value that still names one directory. A root is a
#: path, not a prose field; a newline in it would split a JSONL reader.
_MAX_ROOT_LEN = 512


class HistoryWriteError(RuntimeError):
    """A record could not be appended, so the round was NOT kept.

    The message names the file and the reason. Raised, not logged: the
    caller has to decide what to do about a round it can no longer show.
    """


# ---------------------------------------------------------------------------
# The root
# ---------------------------------------------------------------------------

#: The leaf the store lives in, under the DELFIN state root. ONE name for
#: both sides of the redirect: the resolver below and the
#: state_paths.USER_STATE_RESOLVERS entry that points the suite and
#: `RedirectedUserState` at a run's own copy. If the two ever differ, the
#: redirect sends a test's writes to the old leaf while the code writes to
#: the new one — and the fixture rounds land in the operator's history.
HISTORY_LEAF = "bench_history"


def default_root() -> Path:
    """``~/.delfin/bench_history``, with no setting and no argument.

    Registered in ``state_paths.USER_STATE_RESOLVERS`` under this name,
    which is how the suite and ``RedirectedUserState`` move the whole
    store into a run's own directory: the table replaces this function
    with a zero-argument callable that answers the redirected path, so
    it must take no argument and must NOT read the user's settings
    (settings are redirected separately, and a resolver that consulted
    them could point a test's writes back at the real home).
    """
    return state_paths.state_root() / HISTORY_LEAF


def _setting_root() -> str:
    """``agent.bench_history_root`` from the user's settings, or ""."""
    try:
        from delfin.user_settings import load_settings
        own = (load_settings() or {}).get("agent") or {}
        value = own.get("bench_history_root", "")
    except Exception:
        return ""
    return str(value or "").strip()


def _git_toplevel(path: Path) -> str:
    """Nearest directory holding a git ``.git`` entry at or above *path*.

    A ``.git`` counts when it is a file (a worktree's ``gitdir:`` link) or
    a directory with a ``HEAD`` in it, as git itself requires. An empty
    ``.git`` directory left in a shared parent such as ``/tmp`` is not a
    repository, and treating it as one refused every round below it.
    """
    cur = path
    for _ in range(24):
        try:
            dot = cur / ".git"
            if dot.is_file() or (dot / "HEAD").exists():
                return str(cur)
        except OSError:
            return ""
        parent = cur.parent
        if parent == cur:
            return ""
        cur = parent
    return ""


def root() -> Path:
    """Where the history lives (never inside a repository).

    Order: the ``agent.bench_history_root`` setting, else
    :func:`default_root`. An empty setting falls back to the default like
    every other DELFIN setting — "keep it nowhere" is not a thing this
    setting can say, because the round still has to be kept somewhere and
    silently dropping it is the defect.

    A root inside a git work tree is refused: benchmark records would
    land in the checkout, show up as untracked files in every
    ``git status`` a session runs, and get committed by someone.
    """
    raw = _setting_root()
    if raw:
        if len(raw) > _MAX_ROOT_LEN or "\n" in raw or "\x00" in raw:
            raise HistoryWriteError(
                "agent.bench_history_root is not a usable path (it is "
                f"{len(raw)} chars and must be a single directory path)")
        path = Path(raw).expanduser()
    else:
        path = default_root()
    top = _git_toplevel(path)
    if top:
        raise HistoryWriteError(
            f"the benchmark history root ({path}) lies inside the git "
            f"work tree {top}: benchmark records must never be written "
            "inside a repository. Point agent.bench_history_root at a "
            "directory outside any checkout.")
    return path


# ---------------------------------------------------------------------------
# Writing (append-only; a failure raises)
# ---------------------------------------------------------------------------

def _open_append(path: Path):
    """Open *path* append-only, creating it owner-only.

    ``O_APPEND`` is what makes "never rewritten" true rather than
    intended: the kernel seeks to the end as part of the write, so two
    processes cannot land on the same offset, and there is no open mode
    here that could truncate.
    """
    try:
        state_paths.ensure_dir(path.parent)
        fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_APPEND,
                     state_paths.FILE_MODE)
    except OSError as exc:
        raise HistoryWriteError(
            f"cannot open the benchmark history file {path}: "
            f"{type(exc).__name__}: {exc}") from exc
    try:
        state_paths.secure_file(path)
    except OSError:
        pass        # mode repair is best-effort; the handle is what matters
    return os.fdopen(fd, "a", encoding="utf-8")


def _append_line(path: Path, payload: dict, what: str) -> None:
    """Append exactly one JSON line. Any failure raises HistoryWriteError."""
    try:
        line = json.dumps(payload, sort_keys=True, default=str)
    except (TypeError, ValueError) as exc:
        raise HistoryWriteError(
            f"the {what} record cannot be serialised to JSON: "
            f"{type(exc).__name__}: {exc}") from exc
    try:
        with _open_append(path) as fh:
            fh.write(line + "\n")
            fh.flush()
            os.fsync(fh.fileno())
    except OSError as exc:
        raise HistoryWriteError(
            f"the {what} record could not be appended to {path}: "
            f"{type(exc).__name__}: {exc}. This round is NOT kept — its "
            "runs were paid for and would be unreportable.") from exc


def _root_base(history_root: "Path | str | None") -> Path:
    """The caller's root, or the settings-resolved default."""
    if history_root is None:
        return root()
    text = str(history_root)
    if not text.strip():
        raise HistoryWriteError(
            "the benchmark history root is set but empty: the round has "
            "nowhere to go. Name a directory, or pass nothing to use the "
            "default root.")
    path = Path(text).expanduser()
    # The no-repository rule applies to a caller-named root too: it is
    # the same write, arriving by a different argument, and a flag is
    # exactly how an operator would point it at a checkout by mistake.
    top = _git_toplevel(path)
    if top:
        raise HistoryWriteError(
            f"the benchmark history root ({path}) lies inside the git "
            f"work tree {top}: benchmark records must never be written "
            "inside a repository.")
    return path


def write_round(record_dict: dict, *,
                history_root: "Path | str | None" = None) -> Path:
    """Append the one round line. Raises HistoryWriteError if not written."""
    path = _root_base(history_root) / ROUNDS_FILE
    _append_line(path, record_dict, "round")
    return path


def write_instances(round_id: str, instance_dicts: Iterable[dict], *,
                    history_root: "Path | str | None" = None,
                    started_at: "float | None" = None) -> Path:
    """Append one line per drawn instance, in order.

    All-or-alarm on the file: a round whose instance lines are missing
    would be a round with a header and no measurements, which is exactly
    as unreportable as no record at all — so a write failure raises and
    the caller treats the whole round as not kept.

    ``started_at`` is stamped onto every line so the reader can tell two
    rounds that share a round id apart: the id is not a key here (a
    re-run asked for by name legitimately writes a second round with the
    same id), and without the stamp its instance lines would be merged
    into one of them and the pass rate of both would be wrong.
    """
    path = _root_base(history_root) / INSTANCES_FILE
    for inst in instance_dicts:
        payload = dict(inst)
        payload["round_id"] = round_id
        if started_at is not None:
            payload["round_started_at"] = float(started_at)
        _append_line(path, payload,
                     f"instance record ({payload.get('task_id', '?')})")
    return path


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

@dataclass
class Round:
    """One round line, plus the instance lines that name its round id."""

    round_id: str = ""
    commit: str = ""
    model: str = ""
    provider: str = ""
    seed_mode: str = ""
    task_set_version: str = ""
    started_at: float = 0.0
    ended_at: float = 0.0
    daily_budget_core_hours: float = 0.0
    spent_core_hours: float = 0.0
    refused_over_budget: int = 0
    n_instances: int = 0
    n_unmeasured: int = 0
    input_tokens: int = 0
    output_tokens: int = 0
    instances: list[dict] = field(default_factory=list)

    @property
    def axes(self) -> dict:
        return {"seed_mode": self.seed_mode,
                "task_set_version": self.task_set_version,
                "model": self.model}

    def graded(self) -> list[dict]:
        """The instances that ran AND were graded (the scored population).

        Mirrors ``bench_report._score``: a skipped instance never ran,
        and an unmeasured one ran nothing — counting either as a failure
        would fabricate a measurement.
        """
        return [i for i in self.instances
                if i.get("ran") and not i.get("unmeasured")
                and i.get("success") is not None]

    def pass_rate(self) -> Optional[float]:
        graded = self.graded()
        if not graded:
            return None
        return sum(1 for i in graded if i.get("success")) / len(graded)


def _read_lines(path: Path) -> list[dict]:
    """Every JSON object in *path*; a missing file is an empty history.

    A line that is not a JSON object is skipped and counted by
    :func:`unreadable_count` rather than dropping the whole history: one
    torn line (a machine lost mid-write) must not hide every round.
    """
    out: list[dict] = []
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return out
    except OSError as exc:
        raise HistoryWriteError(
            f"the benchmark history {path} cannot be read: "
            f"{type(exc).__name__}: {exc}") from exc
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            obj = json.loads(line)
        except ValueError:
            _UNREADABLE[path] = _UNREADABLE.get(path, 0) + 1
            continue
        if isinstance(obj, dict):
            out.append(obj)
        else:
            _UNREADABLE[path] = _UNREADABLE.get(path, 0) + 1
    return out


#: Torn/non-object lines seen per file since the process started, so the
#: read side can say "2 lines were unreadable" instead of losing them
#: silently.
_UNREADABLE: dict[Path, int] = {}


def unreadable_count(path: Path) -> int:
    return _UNREADABLE.get(Path(path), 0)


def load(history_root: "Path | str | None" = None) -> list[Round]:
    """Every round, oldest first, with its instances attached.

    ``history_root`` overrides the setting — for a test or a one-off
    directory. Rounds are keyed by (round_id, started_at): a re-run of
    the same round id is a second measurement, not a replacement.
    """
    base = Path(history_root).expanduser() if history_root else root()
    round_rows = _read_lines(base / ROUNDS_FILE)
    inst_rows = _read_lines(base / INSTANCES_FILE)

    by_key: dict[tuple, list[dict]] = {}
    unstamped: dict[str, list[dict]] = {}
    for inst in inst_rows:
        rid = str(inst.get("round_id", "") or "")
        by_key.setdefault((rid, _instance_stamp(inst)), []).append(inst)
        if "round_started_at" not in inst:
            unstamped.setdefault(rid, []).append(inst)

    rounds: list[Round] = []
    for row in round_rows:
        rnd = Round(
            round_id=str(row.get("round_id", "") or ""),
            commit=str(row.get("commit", "") or ""),
            model=str(row.get("model", "") or ""),
            provider=str(row.get("provider", "") or ""),
            seed_mode=str(row.get("seed_mode", "") or ""),
            task_set_version=str(row.get("task_set_version", "") or ""),
            started_at=_row_float(row.get("started_at")),
            ended_at=_row_float(row.get("ended_at")),
            daily_budget_core_hours=_row_float(
                row.get("daily_budget_core_hours")),
            spent_core_hours=_row_float(row.get("spent_core_hours")),
            refused_over_budget=int(row.get("refused_over_budget") or 0),
            n_instances=int(row.get("n_instances") or 0),
            n_unmeasured=int(row.get("n_unmeasured") or 0),
            input_tokens=int(row.get("input_tokens") or 0),
            output_tokens=int(row.get("output_tokens") or 0),
        )
        # Join on (round_id, round started_at). Two rounds that share an
        # id — a re-run the operator asked for by name — therefore keep
        # their own instances, and each pass rate is over its own run.
        key = (rnd.round_id, _round_stamp(row))
        rnd.instances = list(by_key.get(key, []))
        # A round whose lines were written before the stamp existed (or
        # by a caller that passed no started_at) has no key to match:
        # take its unstamped lines, which are the only ones that could
        # belong to it.
        if not rnd.instances and key[0]:
            rnd.instances = list(unstamped.get(key[0], []))
        rounds.append(rnd)
    rounds.sort(key=lambda r: (r.started_at, r.round_id))
    return rounds


def _row_float(value: Any) -> float:
    """A stored field as a float, never raising: 0.0 when it is not one.

    A fourth contract beside the deliberate ``_as_float`` collisions
    (delfin/agent/data_checks.py:47 rejects bool/NaN/inf because they
    poison a statistic; delfin/common/control_validator.py:1011 accepts
    them because CONTROL keys are a different question). Reading a
    history line is neither: one unreadable field in one line must not
    abort a history read, so a value that is not numeric is treated as
    the missing value it effectively is. Named apart so the difference
    is not mistaken for a duplicate.
    """
    try:
        return float(value or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _round_stamp(row: dict) -> float:
    """A round line's key stamp: its own started_at."""
    return _row_float(row.get("started_at"))


def _instance_stamp(row: dict) -> float:
    """An instance line's key stamp: the started_at of its ROUND.

    Deliberately not the instance's duration or any per-instance field:
    the two files are joined on (round_id, round started_at), and an
    instance line carries the round's stamp because that is the only
    thing that can tell two same-id rounds apart.
    """
    return _row_float(row.get("round_started_at"))


def comparable(target: Round, reference: Round) -> tuple[bool, str]:
    """Can *target* be put next to *reference*? Returns (yes, why not).

    The first differing axis names itself. An empty value on either side
    counts as a difference: a record that does not say which seed mode it
    used is not known to match one that does.
    """
    for axis in COMPARABLE_AXES:
        a, b = target.axes[axis], reference.axes[axis]
        if str(a or "") != str(b or ""):
            return False, f"{axis} differs ({b!r} vs {a!r})"
    return True, ""


def not_measurable_reason(rnd: Round) -> str:
    """Why this round cannot carry a pass rate, or "" when it can."""
    if rnd.graded():
        return ""
    if not rnd.instances:
        return ("no instance lines (only the round header was written: "
                "the round measured nothing)")
    started = [i for i in rnd.instances if i.get("ran")]
    if not started:
        return "no instance ran (every draw was refused before starting)"
    if all(i.get("unmeasured") for i in started):
        return ("every started instance was UNMEASURED (the runs never "
                "happened), so the round measured nothing")
    return "no started instance was graded (no pass/fail to count)"


# ---------------------------------------------------------------------------
# The read side's two questions
# ---------------------------------------------------------------------------

def inst_capability(inst: dict) -> str:
    """The capability tag of one stored instance line.

    The line IS an ``InstanceRecord`` plus its round key, so the tag is
    ``task_class`` — DELFIN's own capability vocabulary, the same field
    ``bench_report`` scores "per tool family" from. Named as one accessor
    so the read side and the ``--capability`` filter cannot drift from the
    writer, and so the brief's word ("capability tag") and the record's
    word (``task_class``) cannot be confused for two fields.
    """
    return str(inst.get("task_class") or "")


def inst_score(inst: dict):
    """The graded score of one stored instance: ``quality`` in the
    record's vocabulary, None when it was never graded."""
    return inst.get("quality")


def history_rows(rnd: Round, capability: str = "") -> list[dict]:
    """Per-task rows of one round: task, capability, pass/fail, score."""
    rows = []
    for inst in rnd.instances:
        cap = inst_capability(inst)
        if capability and cap != capability:
            continue
        rows.append({
            "task_id": inst.get("task_id", ""),
            "capability": cap,
            "ran": bool(inst.get("ran")),
            "unmeasured": bool(inst.get("unmeasured")),
            "success": inst.get("success"),
            "score": inst_score(inst),
            "skip_reason": inst.get("skip_reason", ""),
            "error": inst.get("error", ""),
        })
    return rows


def trend(history: list[Round], capability: str = "",
          reference: Optional[Round] = None) -> dict:
    """Pass rate per round (per task inside it), plus why rounds are out.

    ``reference`` defaults to the newest round that can carry a rate.
    Rounds are split three ways, each named:

      ``included``  comparable with the reference, and it measured something
      ``excluded``  dropped, with ``{round_id, reason}``
    """
    usable = [r for r in history
              if not not_measurable_reason(r)
              and (not capability or any(
                  inst_capability(i) == capability for i in r.instances))]
    excluded: list[dict] = []
    if not usable:
        for r in history:
            excluded.append({"round_id": r.round_id,
                             "commit": r.commit,
                             "reason": not_measurable_reason(r)
                             or "no instance of capability "
                                f"{capability!r} in this round"})
        return {"reference": None, "capability": capability,
                "included": [], "excluded": excluded}

    ref = reference or usable[-1]
    included: list[dict] = []
    for r in history:
        if not_measurable_reason(r):
            excluded.append({"round_id": r.round_id, "commit": r.commit,
                             "reason": not_measurable_reason(r)})
            continue
        if capability and not any(inst_capability(i) == capability
                                 for i in r.instances):
            excluded.append({"round_id": r.round_id, "commit": r.commit,
                             "reason": f"no instance of capability "
                                       f"{capability!r}"})
            continue
        ok, why = comparable(r, ref)
        if not ok:
            excluded.append({"round_id": r.round_id, "commit": r.commit,
                             "reason": why})
            continue
        rows = history_rows(r, capability)
        graded = [i for i in rows if i["ran"] and not i["unmeasured"]
                  and i["success"] is not None]
        included.append({
            "round_id": r.round_id,
            "commit": r.commit,
            "started_at": r.started_at,
            "n": len(graded),
            "pass_rate": (sum(1 for i in graded if i["success"]) / len(graded)
                          if graded else None),
            "tasks": rows,
        })
    return {"reference": ref.round_id, "capability": capability,
            "included": included, "excluded": excluded}


def diff_rounds(history: list[Round], a_id: str, b_id: str,
                capability: str = "") -> dict:
    """What changed between two named rounds, task by task.

    Refuses (``error``) when either round is absent or the two are not
    comparable — a delta across two seed modes or task sets is a number
    whose meaning nobody can state, and printing it is how a wrong
    conclusion gets a decimal point.
    """
    by_id: dict[str, Round] = {}
    for r in history:
        by_id.setdefault(r.round_id, r)
    a, b = by_id.get(a_id), by_id.get(b_id)
    if a is None or b is None:
        missing = a_id if a is None else b_id
        return {"error": f"no round with id {missing!r} in the history",
                "known": sorted(by_id)}
    if not_measurable_reason(a) or not_measurable_reason(b):
        return {"error": "one of the two rounds measured nothing: "
                + (not_measurable_reason(a) or not_measurable_reason(b)),
                "round_a": a_id, "round_b": b_id}
    ok, why = comparable(b, a)
    if not ok:
        return {"error": f"{a_id} and {b_id} are not comparable: {why}",
                "round_a": a_id, "round_b": b_id}

    def keyed(rnd: Round) -> dict:
        # Keyed on the BASE task id. bench_holdout_glue gives every draw
        # a unique id ("<instance>.<draw>") so learning state cannot
        # collide; two rounds of the same task then differ only in that
        # suffix, and a trend has to compare like tasks.
        out: dict[str, dict] = {}
        for row in history_rows(rnd, capability):
            out[str(row["task_id"]).rsplit(".", 1)[0]] = row
        return out

    ka, kb = keyed(a), keyed(b)
    tasks = []
    for tid in sorted(set(ka) | set(kb)):
        ra, rb = ka.get(tid), kb.get(tid)
        entry: dict[str, Any] = {"task_id": tid}
        if ra is None or rb is None:
            entry["change"] = ("only in " + (b_id if ra is None else a_id))
            tasks.append(entry)
            continue
        entry["capability"] = rb.get("capability", "")
        entry[a_id] = _cell(ra)
        entry[b_id] = _cell(rb)
        if ra["success"] is True and rb["success"] is False:
            entry["change"] = "regressed"
        elif ra["success"] is False and rb["success"] is True:
            entry["change"] = "fixed"
        else:
            entry["change"] = "same"
        tasks.append(entry)
    pa = _rate(ka)
    pb = _rate(kb)
    return {"round_a": a_id, "round_b": b_id, "commit_a": a.commit,
            "commit_b": b.commit, "capability": capability,
            "pass_rate_a": pa, "pass_rate_b": pb,
            "delta": (None if pa is None or pb is None else pb - pa),
            "comparable_axes": a.axes, "tasks": tasks}


def _cell(row: dict) -> str:
    if not row["ran"]:
        return "not run"
    if row["unmeasured"]:
        return "unmeasured"
    if row["success"] is None:
        return "ungraded"
    return ("pass" if row["success"] else "fail") + f"/{row['score']}"


def _rate(keyed: dict) -> Optional[float]:
    vals = [r for r in keyed.values()
            if r["ran"] and not r["unmeasured"] and r["success"] is not None]
    if not vals:
        return None
    return sum(1 for r in vals if r["success"]) / len(vals)


# ---------------------------------------------------------------------------
# Scheduling: the round as a systemd user service + timer (opt-in)
# ---------------------------------------------------------------------------
#
# The way delfin/agent/swarm/guard.py installs the wave tick: DELFIN's own
# installer writes the units and systemctl enables them. Off by default —
# neither importing this module nor running a round registers anything;
# `bench schedule --enable` is the one way a round lands in the schedule.
#
# A timer beside the service (not a looping service like the tick's): a
# benchmark round is a burst that should end, not a loop that must be
# kept alive, and `Persistent=true` makes a machine that was off at 03:00
# run the missed round on boot instead of silently losing a data point.
#
# The budget gate survives scheduling: `--budget` is part of ExecStart,
# so a fired round refuses instances over the daily core-hour limit
# exactly like a hand-run one. No scheduled path may drop it.

UNIT = "delfin-bench-nightly.service"
TIMER = "delfin-bench-nightly.timer"

#: Daily wall-clock time the round fires at (the operator's local zone;
#: systemd resolves OnCalendar against the machine's clock).
DEFAULT_ON_CALENDAR = "*-*-* 03:00:00"


def _unit_field(value, *, what: str) -> str:
    """One operator-supplied value that will sit on a unit's own line.

    A newline or a carriage return in a unit file ends the directive and starts
    another one: a ``--model`` value carrying a newline produced a SECOND
    ``ExecStart`` line in ``unit_text``'s output. A quote is refused for a
    different reason -- the file's own parser is then guessing where the
    argument ended, and a value that needs escaping is one whose escaping a
    later edit can break. Both are refused where the unit is written, which is
    where a person can still see the mistake; the alternative is a unit
    systemd rejects at 03:00 with nobody watching.

    ``timer_text`` carried its own copy of the newline/quote rule for
    ``OnCalendar``. It shares this one now: two validators for the same act
    drift, and the weaker one wins by being the one a new field is added to.
    """
    text = str(value if value is not None else "")
    if "\n" in text or "\r" in text or '"' in text or "\x00" in text:
        raise ValueError(f"not a usable {what}: {value!r}")
    return text


def unit_dir() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
    return Path(base) / "systemd" / "user"


def _q(text) -> str:
    """One systemd word: quoted, with % written as %%."""
    return '"' + str(text).replace("%", "%%").replace('"', '\\"') + '"'


def round_argv(*, model: str, provider: str = "", instances: int = 10,
               budget: float = 50.0, practice: bool = False,
               seed_path: str = "", python: str = "",
               every: str = "") -> list:
    """The command one round runs — the budget gate included, always.

    ``--budget`` is never omitted: a scheduled round that lost the gate
    would sample without limit, and the point of the gate (PL4: budget
    in code, not in prompts) is that nothing downstream can forget it.

    Every operator string is validated here, once, because this is the one
    place all of them pass through: a value that carries a newline would
    otherwise end the ``ExecStart`` directive that carries it.
    """
    model = _unit_field(model, what="--model value")
    provider = _unit_field(provider, what="--provider value")
    seed_path = _unit_field(seed_path, what="--seed-path value")
    argv = [python or sys.executable, "-P", "-W", "ignore",
            "-m", "delfin.agent.bench_nightly",
            "--model", str(model), "--instances", str(int(instances)),
            "--budget", repr(float(budget))]
    if provider:
        argv += ["--provider", str(provider)]
    if practice:
        argv.append("--practice")
    elif seed_path:
        argv += ["--seed-path", str(seed_path)]
    return argv


def unit_text(*, model: str, provider: str = "", instances: int = 10,
              budget: float = 50.0, practice: bool = False,
              seed_path: str = "", python: str = "") -> str:
    """The oneshot user service that runs exactly one round."""
    argv = round_argv(model=model, provider=provider, instances=instances,
                      budget=budget, practice=practice,
                      seed_path=seed_path, python=python)
    return "\n".join([
        "[Unit]",
        "Description=DELFIN benchmark: one graded nightly round "
        "(kept in the history)",
        "",
        "[Service]",
        "Type=oneshot",
        "Environment=" + _q("PYTHONUNBUFFERED=1"),
        "ExecStart=" + " ".join(_q(a) for a in argv),
        "Nice=10",
        "IOSchedulingClass=idle",
        "NoNewPrivileges=yes",
        "",
        "[Install]",
        "WantedBy=default.target",
        "",
    ])


def timer_text(on_calendar: str = DEFAULT_ON_CALENDAR) -> str:
    """The timer that fires the round service. A bad OnCalendar string is
    refused here, not by systemctl at 03:00 with nobody watching."""
    text = _unit_field(str(on_calendar or "").strip(),
                       what="OnCalendar value")
    if not text:
        raise ValueError(f"not a usable OnCalendar value: {on_calendar!r}")
    return "\n".join([
        "[Unit]",
        "Description=DELFIN benchmark: the nightly round timer",
        "",
        "[Timer]",
        "OnCalendar=" + text,
        # A machine that was off when the round was due runs it on boot:
        # a lost night is a gap in the trend, and a gap reads as "flat".
        "Persistent=true",
        "",
        "[Install]",
        "WantedBy=timers.target",
        "",
    ])


def _systemctl(*args: str):
    import subprocess

    return subprocess.run(["systemctl", "--user", *args],
                          capture_output=True, text=True, timeout=60)


def schedule_status() -> dict:
    """What is installed (nothing, by design, until `--enable` says so)."""
    d = unit_dir()
    return {"service": str(d / UNIT), "timer": str(d / TIMER),
            "installed": (d / UNIT).exists() or (d / TIMER).exists(),
            "service_file": (d / UNIT).exists(),
            "timer_file": (d / TIMER).exists(),
            "active": ""}


def install_schedule(*, model: str, provider: str = "", instances: int = 10,
                     budget: float = 50.0, practice: bool = False,
                     seed_path: str = "", on_calendar: str = "",
                     python: str = "") -> dict:
    """Write the units, reload, enable and start the TIMER (idempotent).

    Refuses when the host has no systemd user manager and says what to do
    instead — a round scheduled into a file no supervisor reads is a
    round that silently never runs, the same failure this package is
    about. Raises ValueError without a model: an ExecStart that cannot
    run is not worth enabling.
    """
    if not str(model or "").strip():
        raise ValueError("a scheduled benchmark round needs --model: the "
                         "round must say which model it measures")
    cal = on_calendar or DEFAULT_ON_CALENDAR
    s_text = unit_text(model=model, provider=provider, instances=instances,
                       budget=budget, practice=practice,
                       seed_path=seed_path, python=python)
    t_text = timer_text(cal)
    try:
        probe = _systemctl("show-environment")
    except Exception as exc:   # noqa: BLE001 — no systemctl to ask (OSError, timeout)
        probe = None
        why = f"{type(exc).__name__}"
    else:
        why = (probe.stderr or "").strip()[:200]
    if probe is None or probe.returncode != 0:
        raise HistoryWriteError(
            f"no systemd user manager here ({why}); run the round by hand "
            "or schedule it with your own supervisor — the exact command "
            "is printed by `bench schedule` with no arguments.")
    d = unit_dir()
    d.mkdir(parents=True, exist_ok=True)
    (d / UNIT).write_text(s_text)
    (d / TIMER).write_text(t_text)
    _systemctl("daemon-reload")
    on = _systemctl("enable", "--now", TIMER)
    active = ""
    try:
        active = _systemctl("is-active", TIMER).stdout.strip()
    except Exception:
        pass
    return {"service": str(d / UNIT), "timer": str(d / TIMER),
            "on_calendar": cal, "active": active,
            "error": (on.stderr or "").strip()[:200] if on.returncode else ""}


def remove_schedule() -> bool:
    """Disable, stop and delete both units. True when either was there."""
    d = unit_dir()
    existed = (d / UNIT).exists() or (d / TIMER).exists()
    for name in (TIMER, UNIT):
        if not (d / name).exists():
            continue
        try:
            _systemctl("disable", "--now", name)
        except Exception:
            pass
        try:
            (d / name).unlink()
        except OSError:
            pass
    try:
        _systemctl("daemon-reload")
    except Exception:
        pass
    return existed


# ---------------------------------------------------------------------------
# The CLI verbs (cli.py registers a subparser per handler and routes here,
# so the read side is reachable through the product's own entry point and
# is tested through it). The text they print comes from the format_*
# functions below, kept separate so both are testable on their own.
# ---------------------------------------------------------------------------

def cli_history(args) -> int:
    """``delfin-agent bench history`` — every kept round, per task.

    Exit 1 when nothing is kept yet: a script that pipes this into a
    comparison must not read "no output" as "flat".
    """
    rounds = load()
    if not rounds:
        print(format_history([]))
        return 1
    print(format_history(rounds, getattr(args, "capability", "") or ""))
    return 0


def cli_trend(args) -> int:
    """``delfin-agent bench trend --capability X`` — pass rate per round.

    --capability is required by the parser, so argparse refuses a bare
    `bench trend` with its usual exit 2 before this runs; no second
    refusal is kept here for a case the parser already turns away.
    """
    capability = str(args.capability or "").strip()
    if not capability:
        # Reachable even though --capability is required: argparse accepts
        # `--capability ""`, and trending over every capability mixed
        # together is the population error the axis exists to prevent.
        print("ERROR: bench trend needs a non-empty --capability (the task "
              "family whose trend you want; `bench history` lists what was "
              "ever run).", file=sys.stderr)
        return 2
    rounds = load()
    if not rounds:
        print(format_history([]))
        return 1
    print(format_trend(rounds, capability))
    return 0


def cli_diff(args) -> int:
    """``delfin-agent bench diff <a> <b>`` — what changed, task by task."""
    rounds = load()
    if not rounds:
        print(format_history([]))
        return 1
    text = format_diff(rounds, args.round_a, args.round_b,
                       getattr(args, "capability", "") or "")
    print(text)
    return 1 if text.startswith("cannot compare") else 0


def cli_schedule(args) -> int:
    """``delfin-agent bench schedule [--enable|--disable]`` — opt-in.

    Nothing is registered without --enable: a benchmark round costs real
    core-hours and model spend, so it must never land on a schedule
    because some other command happened to run.
    """
    if getattr(args, "disable", False):
        removed = remove_schedule()
        print("removed the benchmark round schedule."
              if removed else "no benchmark round schedule was installed.")
        return 0
    if not getattr(args, "enable", False):
        st = schedule_status()
        plan = round_argv(model=args.model or "<model>",
                          provider=args.provider or "",
                          instances=args.instances, budget=args.budget,
                          practice=args.practice,
                          seed_path=args.seed_path or "")
        return _print_schedule_state(st, plan, args)
    if not str(args.model or "").strip():
        print("ERROR: bench schedule --enable needs --model: the round must "
              "say which model it measures.", file=sys.stderr)
        return 2
    try:
        res = install_schedule(model=args.model, provider=args.provider or "",
                               instances=args.instances, budget=args.budget,
                               practice=args.practice,
                               seed_path=args.seed_path or "",
                               on_calendar=args.on_calendar or "")
    except (ValueError, HistoryWriteError) as exc:
        print(f"the benchmark round was NOT scheduled: {exc}",
              file=sys.stderr)
        return 1
    print(f"scheduled the benchmark round: {res['timer']}\n"
          f"  on:     {res['on_calendar']}\n"
          f"  active: {res['active'] or 'unknown'}")
    if res["error"]:
        print(f"  systemctl said: {res['error']}")
    return 0 if not res["error"] else 1


def _print_schedule_state(st: dict, plan: list, args) -> int:
    """The off-by-default state, plus the exact command it would run."""
    print(("benchmark round schedule: INSTALLED" if st["installed"]
           else "benchmark round schedule: off (nothing is scheduled)"))
    print(f"  service: {st['service']}\n  timer:   {st['timer']}")
    if st.get("active"):
        print(f"  active:  {st['active']}")
    print("  the round it would run:\n    "
          + " ".join(str(a) for a in plan))
    print(f"    (the --budget gate is part of that command: a fired round "
          f"still refuses instances over {args.budget:g} core-hours)")
    print("  install it with:  delfin-agent bench schedule --enable "
          "--model M\n  take it off with: delfin-agent bench schedule "
          "--disable")
    return 0


def format_history(rounds: list[Round], capability: str = "") -> str:
    """Per-round, per-task pass rates with the commit — and why the rest
    of the history is not shown."""
    if not rounds:
        return ("no benchmark rounds are kept yet. A round is stored "
                "automatically by `python -m delfin.agent.bench_nightly` "
                "(the module entry point that keeps its record). The "
                "scheduled bench watch cycle does not keep a round: it "
                "writes its run file and report under the runs directory "
                "instead, so this history stays empty afterwards.")
    tr = trend(rounds, capability)
    lines = []
    cap = f" (capability {capability})" if capability else ""
    if tr["reference"] and tr["included"]:
        lines.append(f"benchmark history{cap} — comparable with "
                     f"round {tr['reference']}")
    elif tr["reference"]:
        # A reference exists but nothing compares against it (every other
        # round was excluded): naming it would suggest a comparison that
        # is not there.
        lines.append(f"benchmark history{cap} — reference round "
                     f"{tr['reference']}, nothing comparable with it")
    else:
        lines.append(f"benchmark history{cap}")
    if not tr["included"]:
        lines.append("  no round is comparable; all of them are left out.")
    for row in tr["included"]:
        rate = ("n/a" if row["pass_rate"] is None
                else f"{row['pass_rate'] * 100:.0f}%")
        lines.append(f"  {row['round_id']:<22} {row['commit']:<12} "
                     f"pass {rate:>4}  n={row['n']}")
        for task in row["tasks"]:
            lines.append(f"      {task['task_id']:<30} "
                         f"{task['capability']:<20} {_cell(task)}")
    for row in tr["excluded"]:
        lines.append(f"  excluded {row['round_id']:<20} {row['reason']}")
    return "\n".join(lines)


def format_trend(rounds: list[Round], capability: str) -> str:
    """One line per comparable round for a capability + the excluded ones."""
    tr = trend(rounds, capability)
    lines = [f"trend for capability {capability!r} "
             f"(comparable with round {tr['reference'] or '?'})"]
    if not tr["included"]:
        lines.append("  nothing comparable in the history.")
    for row in tr["included"]:
        rate = ("n/a" if row["pass_rate"] is None
                else f"{row['pass_rate'] * 100:.0f}%")
        lines.append(f"  {row['round_id']:<22} {row['commit']:<12} {rate:>4}"
                     f"  n={row['n']}  "
                     + time.strftime("%Y-%m-%d %H:%M",
                                     time.localtime(row["started_at"])))
    for row in tr["excluded"]:
        lines.append(f"  excluded {row['round_id']:<20} {row['reason']}")
    return "\n".join(lines)


def format_diff(rounds: list[Round], a_id: str, b_id: str,
                capability: str = "") -> str:
    d = diff_rounds(rounds, a_id, b_id, capability)
    if d.get("error"):
        out = [f"cannot compare {a_id} with {b_id}: {d['error']}"]
        if d.get("known"):
            out.append("  known rounds: " + ", ".join(d["known"]))
        return "\n".join(out)
    def pct(x):
        return "n/a" if x is None else f"{x * 100:.0f}%"
    lines = [f"{a_id} ({d['commit_a']}) -> {b_id} ({d['commit_b']})"]
    lines.append(f"  pass rate {pct(d['pass_rate_a'])} -> "
                 f"{pct(d['pass_rate_b'])}"
                 + (f"  ({d['delta'] * 100:+.0f} points)"
                    if d["delta"] is not None else ""))
    lines.append("  same tasks, seed mode, task-set version and model: "
                 + ", ".join(f"{k}={v}"
                             for k, v in d["comparable_axes"].items()))
    for task in d["tasks"]:
        left = task.get(a_id, "-")
        right = task.get(b_id, "-")
        lines.append(f"  {task['task_id']:<30} {task.get('capability', ''):<20}"
                     f" {left:>12} -> {right:<12} {task['change']}")
    return "\n".join(lines)
