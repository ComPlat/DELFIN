"""Auto mode: the harness starts the next iteration when nobody types.

Input: one session's workspace and id, its engine at each call, and four
events from the surface (terminal, window worker, dashboard): the session
went idle (``next_prompt``), an auto turn ended (``finish_iteration``), any
turn ended (``note_turn_end``), a person sent a line
(``note_user_input``). Output: the iteration prompt (a ``HarnessText``,
which grants no push), the veto installed on the permission object for the
length of an iteration, notices, a "needs you" queue, a per-loop log and a
stop report.

Where the work comes from:

* a running loop (``start_loop``, a confirmed :mod:`loops` file): its steps,
  on its trigger -- ``manual`` back to back with ``loop_gap_seconds``
  between, ``idle`` after ``idle_seconds`` without input, ``every N`` at
  most once per interval, ``job`` on a finished-job wake;
* otherwise self-found work (:meth:`find_work`), ranked: a finished job
  (the wake), in-progress tasks of this session, an announced step not
  taken, the last test run failing, red CI on the current branch, pending
  tasks. Registers the user keeps (``TODO.md``, ``.delfin/todo.md``) and
  other sessions' tasks are out of scope: queued as questions, never
  started.

The frame, enforced here and in the gate rather than in the prompt:

* Posture: taken as it is. ``bypassPermissions`` refuses auto mode; a
  posture wider than the one auto was enabled under stops it
  (``auto_guard.mode_rank``). A backend without the native gate (the CLI
  backends, whose gate runs in another process) refuses it.
* Branch: never the default branch. A clean default branch is left for
  ``auto/<stamp>`` (or the loop's branch); a dirty one refuses.
* Per call (``veto``): ``auto_guard.classify_call``; outward and widening
  calls are queued as needs-you items and not run, measure and scope
  writes are refused.
* Per iteration: the changed paths (against the snapshot taken at its
  start; files dirty before it are left out) are committed with explicit
  paths; ``auto_guard.measure_findings`` judges them; the loop's check runs
  through the agent's own gate (``repl.shell_escape_result``); a hard
  finding or a worse measurement is undone with ``git revert`` (history
  kept). One JSON line per iteration in the loop log.
* Stops: iterations, wall time and tokens (the session's token counter
  since the start), a plateau of K iterations without progress (manual,
  idle and self-found work; event triggers are bounded by their events),
  the same failure ``repeat_failures`` times in a row, ``LOOP: blocked``,
  a goal met by the check; ``LOOP: needs-human`` pauses.
* The user first: a person's line refuses the next tool call of a running
  iteration and requests a stop between calls -- the call already running
  finishes -- and pauses auto mode until ``/auto on`` (or, for a loop with
  ``after_user: resume``, until the user's turn has ended).

Progress is measured, never taken from the answer: a metric that moved
the right way, a check that started passing, or (no check) a kept commit.
The ``LOOP:`` line can only stop or pause a loop, never count as progress.

Rejected: registering iterations in the scheduler (``~/.delfin/cron.json``)
-- a scheduled entry fires a fresh turn in whatever session holds the
listener, outside this session's budget and veto; a model-written
"progress" field (the measure would be the answer's own claim).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional

from . import auto_guard, loops
from .harness_rows import HarnessText

__all__ = ["AutoController", "WorkItem", "DEFAULTS"]

DEFAULTS = {
    "idle_seconds": 120,
    "loop_gap_seconds": 5,
    "max_iterations": 20,
    "wall": 4 * 3600,
    "tokens": 2_000_000,
    "plateau": 3,
    "repeat_failures": 2,
    "after_user": "pause",
    "registers": ["TODO.md", ".delfin/todo.md"],
    "check_timeout": 900,
}

_MARKER = re.compile(r"^\s*LOOP:\s*(progress|done|blocked|needs-human)\b:?\s*(.*)$",
                     re.IGNORECASE | re.MULTILINE)
_UNCHECKED = re.compile(r"^\s*[-*]\s+\[ \]\s+(\S.*)$", re.MULTILINE)


@dataclass
class WorkItem:
    rank: int
    kind: str
    text: str
    in_scope: bool
    ref: str = ""


@dataclass
class _Run:
    spec: loops.LoopSpec
    started_at: float
    tokens_at_start: int
    is_loop: bool
    iteration: int = 0
    kept: int = 0
    reverted: int = 0
    step_index: int = 0
    last_iteration_at: float = 0.0
    last_metric: Optional[float] = None
    last_exit: Optional[int] = None
    no_progress: int = 0
    fail_fp: str = ""
    fail_count: int = 0
    last_item: str = ""
    last_result: str = ""


@dataclass
class _Iteration:
    n: int
    started_at: float
    head: str
    before: dict
    item: str = ""
    preempted: bool = False
    stop_sent: bool = False
    engine: Any = None


# ---------------------------------------------------------------------------
# git, through subprocess, in the workspace only
# ---------------------------------------------------------------------------

class _Git:
    def __init__(self, workspace: Path):
        self.ws = Path(workspace)

    def run(self, *args: str, timeout: float = 30,
            raw: bool = False) -> tuple[int, str]:
        try:
            r = subprocess.run(["git", "-C", str(self.ws), *args],
                               capture_output=True, text=True, timeout=timeout,
                               stdin=subprocess.DEVNULL)
        except Exception as exc:
            return 1, str(exc)
        out = r.stdout if r.returncode == 0 else (r.stderr or r.stdout)
        return r.returncode, out if raw else out.strip()

    def is_repo(self) -> bool:
        code, out = self.run("rev-parse", "--is-inside-work-tree")
        return code == 0 and out == "true"

    def branch(self) -> str:
        code, out = self.run("branch", "--show-current")
        return out if code == 0 else ""

    def head(self) -> str:
        code, out = self.run("rev-parse", "HEAD")
        return out if code == 0 else ""

    def defaults(self) -> set[str]:
        names = {"main", "master"}
        code, out = self.run("symbolic-ref", "--short", "refs/remotes/origin/HEAD")
        if code == 0 and out:
            names.add(out.split("/", 1)[-1])
        return names

    def dirty(self) -> dict[str, str]:
        """{relpath: content hash or '' (deleted)} of every changed,
        added or untracked path."""
        code, out = self.run("status", "--porcelain=v1", "-z",
                             "--untracked-files=all", raw=True)
        if code != 0:
            return {}
        res: dict[str, str] = {}
        parts = out.split("\0")
        i = 0
        while i < len(parts):
            entry = parts[i]
            i += 1
            if len(entry) < 4:
                continue
            xy, path = entry[:2], entry[3:]
            if xy[0] in "RC":
                i += 1             # the source path follows
            res[path] = _hash(self.ws / path)
        return res

    def tracked(self, path: str) -> bool:
        code, _ = self.run("ls-files", "--error-unmatch", "--", path)
        return code == 0

    def show(self, rev: str, path: str) -> str:
        code, out = self.run("show", f"{rev}:{path}", raw=True)
        return out if code == 0 else ""

    def switch(self, branch: str) -> tuple[bool, str]:
        code, _ = self.run("rev-parse", "--verify", "--quiet",
                           f"refs/heads/{branch}")
        args = ("switch", branch) if code == 0 else ("switch", "-c", branch)
        code, out = self.run(*args)
        return code == 0, out

    def commit(self, paths: list[str], message: str) -> str:
        code, out = self.run("add", "-A", "--", *paths)
        if code != 0:
            return ""
        code, out = self.run("commit", "-q", "--no-verify", "-m", message,
                             "--", *paths)
        return self.head() if code == 0 else ""

    def revert(self, sha: str) -> bool:
        code, _ = self.run("revert", "--no-edit", sha)
        if code != 0:
            self.run("revert", "--abort")
        return code == 0


def _hash(path: Path) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return ""


def _read(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""


def _any_text(content: Any) -> str:
    """Text of a message's content: a string, or the text and content
    fields of its blocks (tool results included)."""
    if isinstance(content, str):
        return content
    out = []
    for b in content if isinstance(content, list) else []:
        if isinstance(b, dict):
            for k in ("text", "content", "output"):
                v = b.get(k)
                if isinstance(v, str):
                    out.append(v)
                elif isinstance(v, list):
                    out.append(_any_text(v))
    return "\n".join(out)


def _short(seconds: float) -> str:
    seconds = max(0, int(seconds))
    if seconds >= 3600:
        return f"{seconds // 3600}h{(seconds % 3600) // 60:02d}m"
    if seconds >= 60:
        return f"{seconds // 60}m"
    return f"{seconds}s"


def _tokens(engine) -> int:
    usage = getattr(engine, "token_usage", None) or {}
    try:
        return int(usage.get("input", 0) or 0) + int(usage.get("output", 0) or 0)
    except Exception:
        return 0


def _default_check(engine, command: str) -> tuple[int, str]:
    """The check through the agent's own gate, as ``!cmd`` runs."""
    from .repl import shell_escape_result
    run = shell_escape_result(engine, command)
    if run.get("refused") or run.get("error"):
        return 126, str(run.get("refused") or run.get("error"))
    out = "\n".join(s for s in (run.get("stdout"), run.get("stderr")) if s)
    code = run.get("exit_code")
    return (code if isinstance(code, int) else 1), out


def _default_notify(text: str, kind: str = "confirm_pending") -> None:
    """The attention inbox, which the window's attention feed, /attention
    and the dashboard read: a needs-you item is ``confirm_pending``, a stop
    report ``run_finished``."""
    try:
        from . import attention
        attention.emit_attention(kind, title="auto mode",
                                 detail=str(text)[:2000])
    except Exception:
        pass


def _settings() -> dict:
    try:
        from delfin.user_settings import load_settings
        raw = ((load_settings() or {}).get("agent") or {}).get("auto_mode")
        return raw if isinstance(raw, dict) else {}
    except Exception:
        return {}


# ---------------------------------------------------------------------------
# the controller
# ---------------------------------------------------------------------------

class AutoController:
    """One session's auto mode. See the module docstring."""

    def __init__(self, workspace: Path | str, *, session_id: str = "",
                 clock: Callable[[], float] = time.time,
                 run_check: Optional[Callable[[Any, str], tuple[int, str]]] = None,
                 notify: Optional[Callable[[str], None]] = None,
                 settings: Optional[dict] = None,
                 ci_probe: Optional[Callable[[str], str]] = None):
        self.workspace = Path(workspace)
        self.session_id = str(session_id or "")
        self.clock = clock
        self.run_check = run_check or _default_check
        self.notify = notify or _default_notify
        self.settings = {**DEFAULTS, **(_settings() if settings is None else {}),
                         **(settings or {})}
        self.ci_probe = ci_probe
        self.git = _Git(self.workspace)
        self.enabled = False
        self.paused = ""
        self.posture = ""
        self.loop: Optional[_Run] = None
        self._run: Optional[_Run] = None
        self.iteration: Optional[_Iteration] = None
        self.needs_you: list[dict] = []
        self.stop_report = ""
        self.last_turn_end = clock()
        self._notices: list[str] = []
        self._said = ""
        self._resume_after_user = False
        self._user_turn = False
        self.user_waiting: Callable[[], bool] = lambda: False
        self.log_path = self._log_path("auto")

    # -- small helpers -------------------------------------------------------
    def _log_path(self, name: str) -> Path:
        key = hashlib.sha256(str(self.workspace.resolve()).encode()).hexdigest()[:12]
        return Path.home() / ".delfin" / "loop_logs" / key / f"{name}.jsonl"

    def _log(self, row: dict) -> None:
        try:
            self.log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(self.log_path, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(row, sort_keys=True) + "\n")
            try:
                os.chmod(self.log_path, 0o600)
            except OSError:
                pass
        except OSError:
            pass

    def _tell(self, text: str, kind: str = "confirm_pending") -> None:
        try:
            try:
                self.notify(text, kind)
            except TypeError:
                self.notify(text)
        except Exception:
            pass

    def _say(self, text: str) -> None:
        self._notices.append(text)

    def pop_notices(self) -> list[str]:
        out, self._notices = self._notices, []
        return out

    def _queue(self, kind: str, text: str, why: str = "") -> dict:
        for item in self.needs_you:
            if item["text"] == text and item["kind"] == kind:
                return item
        item = {"id": len(self.needs_you) + 1, "kind": kind, "text": text,
                "why": why, "at": self.clock(),
                "loop": self._run.spec.name if self._run else "",
                "iteration": self.iteration.n if self.iteration else 0}
        self.needs_you.append(item)
        self._tell(f"needs you (#{item['id']}): {text}"
                   + (f" — {why}" if why else ""))
        return item

    def _perms(self, engine):
        try:
            return getattr(engine, "kit_permissions", None)
        except Exception:
            return None

    # -- enabling ------------------------------------------------------------
    def _frame_refusal(self, engine) -> str:
        perms = self._perms(engine)
        if perms is None:
            return ("auto mode needs DELFIN's own tool gate, and this backend "
                    "has none (a CLI backend's gate runs in another process "
                    "and cannot see an iteration). Not enabled.")
        mode = str(getattr(perms, "mode", "") or "")
        if mode == "bypassPermissions":
            return ("auto mode does not run under bypassPermissions: nothing "
                    "would ask before a write or a command while nobody "
                    "watches. Switch to acceptEdits or default first "
                    "(/mode), then /auto on.")
        return ""

    def _own_branch(self, wanted: str = "") -> tuple[bool, str]:
        if not self.git.is_repo():
            return True, ("not a git repository: iterations are not committed "
                          "and cannot be reverted; a worse result stops the "
                          "run instead.")
        current = self.git.branch()
        defaults = self.git.defaults()
        if wanted and current == wanted:
            return True, ""
        if not wanted and current and current not in defaults:
            return True, ""
        # Tracked changes only: untracked files (a new loop file among
        # them) stay where they are on a switch, and an iteration leaves
        # every path that was dirty before it out of its commit anyway.
        dirty = {p for p in self.git.dirty()
                 if self.git.tracked(p) and not p.startswith(".delfin/")}
        if dirty:
            return False, (f"the workspace is on {current or 'a detached HEAD'}"
                           f" with {len(dirty)} uncommitted change(s); auto "
                           "mode commits every iteration on its own branch and "
                           "will not carry them along. Commit or move them, "
                           "then try again.")
        target = wanted or time.strftime("auto/%Y%m%d-%H%M",
                                         time.localtime(self.clock()))
        if target in defaults:
            return False, f"{target} is the default branch; a loop never commits there."
        ok, out = self.git.switch(target)
        if not ok:
            return False, f"could not switch to {target}: {out[:200]}"
        return True, f"working on branch {target} (from {current or 'HEAD'})."

    def _new_run(self, spec: loops.LoopSpec, engine, is_loop: bool) -> _Run:
        return _Run(spec=spec, started_at=self.clock(),
                    tokens_at_start=_tokens(engine), is_loop=is_loop)

    def _auto_spec(self) -> loops.LoopSpec:
        s = self.settings
        return loops.LoopSpec(
            name="auto", goal="the next item of work in this session's scope",
            trigger="idle", scope="auto",
            max_iterations=int(s["max_iterations"]),
            wall_seconds=int(s["wall"]), max_tokens=int(s["tokens"]),
            plateau=int(s["plateau"]),
            repeat_failures=int(s["repeat_failures"]),
            after_user=str(s["after_user"]))

    def enable(self, engine) -> str:
        """``/auto on``. Returns the line to show."""
        if self.enabled and self.paused:
            return self.resume()
        if self.enabled:
            return "auto mode is already on. " + self.status_line()
        refusal = self._frame_refusal(engine)
        if refusal:
            return refusal
        ok, note = self._own_branch()
        if not ok:
            return "auto mode not enabled: " + note
        self.posture = str(getattr(self._perms(engine), "mode", "") or "")
        self.enabled = True
        self.paused = ""
        self.stop_report = ""
        self._said = ""
        self._run = self._new_run(self._auto_spec(), engine, False)
        self.log_path = self._log_path("auto")
        self.last_turn_end = self.clock()
        return (f"auto mode on (posture {self.posture}): after "
                f"{_short(self.settings['idle_seconds'])} without input the "
                "next item of work in scope starts by itself. Pushes, pull "
                "requests, merges and messages out are queued for you, never "
                "run. /auto off stops it; typing anything pauses it. "
                + note).strip()

    def start_loop(self, spec: Optional[loops.LoopSpec], engine) -> str:
        """``/loop run <name>``. Returns the line to show."""
        if spec is None:
            return "no such loop."
        if spec.errors:
            return f"loop {spec.name} is not valid: " + "; ".join(spec.errors)
        if not loops.is_confirmed(spec):
            return (f"loop {spec.name} is not confirmed for its current "
                    f"content. Read it with /loop show {spec.name}, then "
                    f"/loop confirm {spec.name}.")
        if self.loop is not None:
            return (f"loop {self.loop.spec.name} is running; /loop stop first.")
        refusal = self._frame_refusal(engine)
        if refusal:
            return refusal
        ok, note = self._own_branch(spec.branch)
        if not ok:
            return f"loop {spec.name} not started: " + note
        run = self._new_run(spec, engine, True)
        if spec.check:
            code, out = self._check(engine, spec.check)
            run.last_exit, run.last_metric = code, spec.metric_value(out)
            if code == 0:
                return (f"loop {spec.name}: the goal is already met (the check "
                        f"exits 0). Nothing started. {note}").strip()
        self.posture = str(getattr(self._perms(engine), "mode", "") or "")
        self.loop = self._run = run
        self.enabled = True
        self.paused = ""
        self.stop_report = ""
        self.log_path = self._log_path(spec.name)
        self.last_turn_end = self.clock() - float(self.settings["loop_gap_seconds"])
        base = ""
        if run.last_metric is not None:
            base = f" Baseline metric {run.last_metric:g}."
        return (f"loop {spec.name} started on branch "
                f"{self.git.branch() or '(no git)'}: {spec.summary()}.{base} "
                f"/loop stop ends it; typing anything pauses it. {note}").strip()

    def disable(self, reason: str = "turned off") -> str:
        if not self.enabled and self.loop is None:
            return "auto mode is off."
        return self.stop(reason)

    def resume(self) -> str:
        if not self.enabled:
            return "auto mode is off; /auto on or /loop run <name> starts it."
        self.paused = ""
        self._resume_after_user = False
        self.last_turn_end = self.clock() - float(self.settings["idle_seconds"])
        return "auto mode resumed. " + self.status_line()

    def stop(self, reason: str) -> str:
        run = self._run
        name = run.spec.name if run else "auto"
        lines = [f"auto stopped — {name}: {reason}."]
        if run:
            used = max(0, self._tokens_used_cache)
            lines.append(
                f"{run.iteration} iteration(s): {run.kept} kept, "
                f"{run.reverted} reverted; {_short(self.clock() - run.started_at)}"
                f" and ~{used} tokens used.")
            if run.last_metric is not None:
                lines.append(f"last metric {run.last_metric:g}.")
            lines.append(f"log: {self.log_path}")
        if self.needs_you:
            lines.append(f"needs you ({len(self.needs_you)}): "
                         + "; ".join(f"#{n['id']} {n['text'][:80]}"
                                     for n in self.needs_you[-5:]))
        report = " ".join(lines)
        if self.iteration is not None:
            self._disarm()
            self.iteration = None
        self.stop_report = report
        self.enabled = False
        self.loop = None
        self._run = None
        self.paused = ""
        self._say(report)
        self._tell(report, "run_finished")
        return report

    _tokens_used_cache = 0

    # -- what the surface reports ----------------------------------------------
    def note_turn_end(self, *, auto_turn: bool = False) -> None:
        self.last_turn_end = self.clock()
        if not auto_turn and self._user_turn:
            self._user_turn = False
            if self._resume_after_user and self.enabled:
                self.paused = ""
                self._resume_after_user = False

    def note_user_input(self, engine=None, *, steer: bool = False) -> bool:
        """A person sent a line. True when it preempts a running iteration.

        ``steer``: the line went INTO the running turn (the user directs
        it); the turn goes on, inside the frame, and auto mode pauses after
        it instead of refusing its next call."""
        if not self.enabled:
            return False
        after = self._run.spec.after_user if self._run else "pause"
        if after == "resume":
            self._resume_after_user = True
        self.paused = "the user sent a message"
        if self.iteration is not None and not steer:
            self.iteration.preempted = True
            self._user_turn = True
            return True
        self._user_turn = self.iteration is None
        return False

    # -- the veto ------------------------------------------------------------
    def veto(self, name: str, args: Optional[dict]) -> Optional[str]:
        it = self.iteration
        if it is None:
            return None
        waiting = False
        try:
            waiting = bool(self.user_waiting())
        except Exception:
            waiting = False
        if it.preempted or waiting:
            it.preempted = True
            self.paused = self.paused or "the user sent a message"
            if not it.stop_sent:
                it.stop_sent = True
                try:
                    if it.engine is not None:
                        it.engine.request_stop()
                except Exception:
                    pass
            return ("refused: the user sent a message, and the user comes "
                    "first. This auto iteration ends here: stop calling "
                    "tools and say in one line where the work stands.")
        spec = self._run.spec if self._run else None
        verdict = auto_guard.classify_call(
            name, args or {}, workspace=str(self.workspace),
            measure=spec.measure if spec else (),
            touch=spec.touch if spec else (),
            default_branches=sorted(self.git.defaults()))
        if verdict is None:
            return None
        if verdict.queue:
            what = str((args or {}).get("command") or name)
            item = self._queue(verdict.kind, f"{name}: {what[:200]}",
                               verdict.reason)
            return (f"blocked: auto mode does not do this ({verdict.reason}); "
                    f"it was not run and is queued for the user as needs-you "
                    f"item #{item['id']}. Go on with local work, or end the "
                    "iteration with 'LOOP: needs-human <question>' if the "
                    "work cannot continue without it.")
        return (f"blocked: {verdict.reason}. An auto iteration may not change "
                "the measure it is judged by, nor anything outside what the "
                "loop may touch. Not run.")

    def _arm(self, engine) -> None:
        perms = self._perms(engine)
        if perms is not None:
            try:
                setattr(perms, "auto_veto", self.veto)
            except Exception:
                pass

    def _disarm(self) -> None:
        eng = self.iteration.engine if self.iteration else None
        perms = self._perms(eng)
        if perms is not None and getattr(perms, "auto_veto", None) == self.veto:
            try:
                setattr(perms, "auto_veto", None)
            except Exception:
                pass

    # -- deciding the next iteration -------------------------------------------
    def _budget_stop(self, engine) -> str:
        run = self._run
        if run is None:
            return ""
        spec = run.spec
        used = _tokens(engine) - run.tokens_at_start
        self._tokens_used_cache = used
        if run.iteration >= spec.max_iterations:
            return (f"budget spent: {run.iteration} of {spec.max_iterations} "
                    "iterations")
        if self.clock() - run.started_at >= spec.wall_seconds:
            return f"budget spent: time ({_short(spec.wall_seconds)})"
        if used >= spec.max_tokens:
            return f"budget spent: tokens ({used} of {spec.max_tokens})"
        return ""

    def _posture_stop(self, engine) -> str:
        perms = self._perms(engine)
        if perms is None:
            return "the session no longer has DELFIN's tool gate"
        mode = str(getattr(perms, "mode", "") or "")
        if mode == "bypassPermissions" or auto_guard.widened(self.posture, mode):
            return (f"the posture changed from {self.posture} to {mode} since "
                    "auto was enabled; auto mode never runs wider than it "
                    "started")
        return ""

    def next_prompt(self, engine, *, typed: str = "", wake: str = "",
                    work: Optional[list[WorkItem]] = None,
                    pause_key: str = "") -> str:
        """The prompt of the next iteration, or "" (see the module docstring).
        Never raises."""
        try:
            return self._next(engine, typed, wake, work, pause_key)
        except Exception as exc:
            self._say(f"auto: the next iteration could not be prepared ({exc})")
            return ""

    def _next(self, engine, typed, wake, work, pause_key) -> str:
        if not self.enabled or self._run is None or self.iteration is not None:
            return ""
        if str(typed or "").strip() or self.paused:
            return ""
        if pause_key:
            try:
                from . import session_pause
                if session_pause.wake_blocked(pause_key):
                    return ""
            except Exception:
                pass
        why = self._posture_stop(engine) or self._budget_stop(engine)
        if why:
            self.stop(why)
            return ""
        run, spec, now = self._run, self._run.spec, self.clock()
        idle = now - self.last_turn_end
        ranked: list[WorkItem] = []
        item = ""
        if run.is_loop:
            if spec.trigger == "job":
                if not wake:
                    return ""
            elif not wake:
                if spec.trigger == "manual" and idle < self.settings["loop_gap_seconds"]:
                    return ""
                if spec.trigger == "idle" and idle < self.settings["idle_seconds"]:
                    return ""
                if spec.trigger == "every" and (
                        idle < self.settings["loop_gap_seconds"]
                        or (run.last_iteration_at
                            and now - run.last_iteration_at < spec.every_seconds)):
                    return ""
        else:
            if not wake and idle < self.settings["idle_seconds"]:
                return ""
            ranked = list(work) if work is not None else self.find_work(engine)
            if wake:
                ranked.insert(0, WorkItem(0, "job", "evaluate the finished job "
                                          "the wake below names", True))
            ranked.sort(key=lambda w: w.rank)
            for w in ranked:
                if not w.in_scope:
                    self._queue("out_of_scope", w.text,
                                "found by auto mode outside this session's "
                                "task; say if it should be started")
            inscope = [w for w in ranked if w.in_scope]
            if not inscope:
                state = f"none:{len(self.needs_you)}"
                if self._said != state:
                    self._said = state
                    self._say("auto: nothing in scope to do"
                              + (f"; {len(self.needs_you)} item(s) need you "
                                 "(/auto needs)" if self.needs_you else "."))
                return ""
            item = inscope[0].text
            self._said = ""
            self._say("auto: next → " + item + (
                "\n  also found: " + "; ".join(
                    ("" if w.in_scope else "[asks] ") + w.text[:70]
                    for w in ranked if w.text != item)[:600]
                if len(ranked) > 1 else ""))
        return self._begin(engine, wake, item)

    def _begin(self, engine, wake: str, item: str) -> str:
        run = self._run
        run.iteration += 1
        run.last_iteration_at = self.clock()
        self.iteration = _Iteration(
            n=run.iteration, started_at=self.clock(),
            head=self.git.head() if self.git.is_repo() else "",
            before=self.git.dirty() if self.git.is_repo() else {},
            item=item, engine=engine)
        self._arm(engine)
        return HarnessText(self._prompt(wake, item))

    def _budget_left(self) -> str:
        run, spec = self._run, self._run.spec
        return (f"{spec.max_iterations - run.iteration} iterations, "
                f"{_short(spec.wall_seconds - (self.clock() - run.started_at))}, "
                f"{max(0, spec.max_tokens - self._tokens_used_cache)} tokens")

    def _prompt(self, wake: str, item: str) -> str:
        run, spec = self._run, self._run.spec
        head = ("[auto — a system event, not the user] Nobody typed this. "
                + (f"Loop \"{spec.name}\", iteration {run.iteration} of "
                   f"{spec.max_iterations}" if run.is_loop else
                   f"Auto mode, iteration {run.iteration} of {spec.max_iterations}")
                + f" (left after this one: {self._budget_left()}).")
        parts = [head]
        if run.is_loop:
            parts.append(f"Goal: {spec.goal}")
            steps = spec.steps
            if spec.per_iteration == "one" and steps:
                k = min(run.step_index, len(steps) - 1)
                steps = [steps[k]]
                parts.append(f"This iteration: step {k + 1} of {len(spec.steps)}.")
            for i, st in enumerate(steps, 1):
                line = f"{i}. {st.text}"
                if st.run:
                    line += f"  (run: {st.run})"
                if st.check:
                    line += f"  (check: {st.check})"
                parts.append(line)
        else:
            parts.append("The harness picked this item from the session's own "
                         f"work: {item}. Do that item only.")
        if wake:
            parts.append(str(wake).strip())
        branch = self.git.branch() or "(no git)"
        frame = [f"Frame, enforced by the harness: work locally on branch "
                 f"{branch}; after this turn the harness commits what you "
                 "changed"]
        if spec.check:
            frame.append(f", runs `{spec.check}` and keeps or reverts the "
                         "change by its result")
        frame.append(". Pushes, pull requests, merges, branch deletions, "
                     "messages and uploads leaving this machine are refused "
                     "and queued for the user.")
        if spec.measure:
            frame.append(" Do not change: " + ", ".join(spec.measure[:8]) + ".")
        if spec.touch:
            frame.append(" Change only: " + ", ".join(spec.touch[:8]) + ".")
        frame.append(" More in the `loops` skill.")
        parts.append("".join(frame))
        parts.append("End with one line: `LOOP: progress <what changed>`, "
                     "`LOOP: done`, `LOOP: blocked <why>` or "
                     "`LOOP: needs-human <question>`. Progress is what the "
                     "harness measures, not what this line says.")
        return "\n".join(parts)

    def _check(self, engine, command: str) -> tuple[int, str]:
        try:
            code, out = self.run_check(engine, command)
            return int(code), str(out or "")
        except Exception as exc:
            return 125, f"the check could not run: {exc}"

    def cancel_iteration(self, why: str = "") -> None:
        """The surface could not start the turn the last prompt was for:
        take the iteration back (veto off, not counted). Never raises."""
        if self.iteration is None:
            return
        try:
            self._disarm()
        finally:
            self.iteration = None
            if self._run is not None and self._run.iteration > 0:
                self._run.iteration -= 1
            if why:
                self._say(f"auto: iteration not started ({why})")

    # -- the end of an iteration -----------------------------------------------
    def finish_iteration(self, engine, answer: str = "", error: str = "") -> str:
        """Judge, commit or revert, log; returns the digest ("" when no
        iteration was running). Never raises."""
        it = self.iteration
        if it is None:
            self.note_turn_end()
            return ""
        try:
            return self._finish(engine, it, str(answer or ""), str(error or ""))
        except Exception as exc:
            self._disarm()
            self.iteration = None
            self.note_turn_end(auto_turn=True)
            return self.stop(f"the iteration could not be judged ({exc})")

    def _changes(self, it: _Iteration) -> dict[str, tuple[str, str]]:
        after = self.git.dirty()
        out: dict[str, tuple[str, str]] = {}
        for path, digest in after.items():
            if path in it.before:
                continue                  # dirty before: not this iteration's
            if path.startswith(".delfin/") and not path.startswith(".delfin/loops/"):
                continue                  # the agent's own state files
            old = self.git.show(it.head, path) if it.head else ""
            new = _read(self.workspace / path) if digest else ""
            out[path] = (old, new)
        return out

    def _finish(self, engine, it: _Iteration, answer: str, error: str) -> str:
        run, spec = self._run, self._run.spec
        self._disarm()
        self.iteration = None
        self.note_turn_end(auto_turn=True)
        self._tokens_used_cache = _tokens(engine) - run.tokens_at_start
        marker = None
        for m in _MARKER.finditer(answer):
            marker = (m.group(1).lower(), m.group(2).strip())
        repo = bool(it.head)
        changes = self._changes(it) if repo else {}
        changed = sorted(changes)
        findings = auto_guard.measure_findings(changes, measure=spec.measure)
        hard = [f for f in findings if f.hard]
        soft = [f for f in findings if not f.hard]
        sha = ""
        if changed and repo:
            label = spec.name if run.is_loop else "auto"
            first = (marker[1] if marker and marker[1] else
                     (answer.strip().splitlines() or [""])[-1])[:60]
            sha = self.git.commit(changed, (
                f"{label} iteration {it.n}: {first or 'changes'}\n\n"
                f"Committed by DELFIN auto mode; files: {', '.join(changed)[:400]}"))
        kept, why = True, ""
        code: Optional[int] = None
        out = ""
        metric: Optional[float] = None
        progress = False
        if hard:
            kept, why = False, ("changed the measure: "
                                + ", ".join(f.path for f in hard))
        elif spec.check:
            code, out = self._check(engine, spec.check)
            metric = spec.metric_value(out)
            prev_m, prev_c = run.last_metric, run.last_exit
            worse = False
            if metric is not None and prev_m is not None:
                worse = (metric > prev_m if spec.better == "lower"
                         else metric < prev_m)
                progress = (metric < prev_m if spec.better == "lower"
                            else metric > prev_m)
            if prev_c == 0 and code != 0:
                worse = True
            if code == 0 and prev_c != 0:
                progress = True
            if worse and changed:
                kept = False
                why = (f"measured worse: metric {prev_m:g} → {metric:g}"
                       if metric is not None and prev_m is not None else
                       f"measured worse: the check went from exit {prev_c} to "
                       f"{code}")
        else:
            progress = bool(changed)
        if not kept and sha:
            if self.git.revert(sha):
                run.reverted += 1
            else:
                return self.stop(f"iteration {it.n} had to be reverted "
                                 f"({why}) and git revert failed; the commit "
                                 f"{sha[:10]} is still there")
        elif not kept and changed and not repo:
            return self.stop(f"iteration {it.n}: {why}; without git it cannot "
                             "be undone")
        elif not kept:
            pass
        if kept:
            if changed:
                run.kept += 1
            if spec.check and code is not None:
                run.last_exit, run.last_metric = code, (
                    metric if metric is not None else run.last_metric)
        else:
            progress = False
        for f in soft:
            self._queue("review", f"iteration {it.n} changed {f.path}: "
                        f"{f.detail}", "kept, but a person should look")
        # repeated identical failure
        fp = ""
        if error:
            fp = "error:" + error[:200]
        elif spec.check and code not in (None, 0) and metric is None:
            tail = [ln.strip() for ln in out.splitlines() if ln.strip()][-3:]
            fp = f"{code}:" + " | ".join(tail)[:300]
        elif not run.is_loop and not changed and it.item and it.item == run.last_item:
            fp = "item:" + it.item[:200]
        if fp and fp == run.fail_fp:
            run.fail_count += 1
        elif fp:
            run.fail_fp, run.fail_count = fp, 1
        else:
            run.fail_fp, run.fail_count = "", 0
        run.last_item = it.item
        run.no_progress = 0 if progress else run.no_progress + 1
        if spec.per_iteration == "one" and kept and spec.steps:
            st = spec.steps[min(run.step_index, len(spec.steps) - 1)]
            ok = True
            if st.check:
                ok = self._check(engine, st.check)[0] == 0
            if ok:
                run.step_index += 1
        result = ("kept" if kept else "reverted") + (
            f" ({sha[:10]})" if sha and kept else "")
        if not changed:
            result = "no change"
        if metric is not None:
            result += f" · metric {metric:g}"
        elif code is not None:
            result += f" · check exit {code}"
        run.last_result = result
        self._log({"at": round(self.clock(), 3), "loop": spec.name,
                   "iteration": it.n, "item": it.item,
                   "tried": answer.strip()[-400:], "changed": changed,
                   "commit": sha, "kept": kept, "why": why,
                   "check_exit": code, "metric": metric,
                   "flagged": [f"{f.path}: {f.detail}" for f in soft],
                   "marker": list(marker) if marker else None,
                   "tokens": self._tokens_used_cache,
                   "measure_source": auto_guard.measure_source(),
                   "preempted": it.preempted, "error": error[:300]})
        digest = (f"auto · {spec.name} · iteration {it.n}: {result}"
                  + (f" — {why}" if why else "")
                  + f" · left: {self._budget_left()}")
        self._say(digest)
        # stops, most specific first
        if spec.check and code == 0 and kept:
            self.stop(f"goal met (`{spec.check}` exits 0)")
            return digest
        if marker and marker[0] == "blocked":
            self._queue("blocked", marker[1] or "the agent reports it is blocked")
            self.stop(f"blocked: {marker[1] or 'no reason given'}")
            return digest
        if marker and marker[0] == "done" and not spec.check:
            self.stop("the agent reports the goal met (no check confirms it)")
            return digest
        if marker and marker[0] == "needs-human":
            self._queue("question", marker[1] or "the agent needs a decision")
            self.paused = "waiting for your answer"
        if run.fail_count >= spec.repeat_failures:
            self.stop(f"the same failure {run.fail_count} times in a row: "
                      f"{run.fail_fp[:160]}")
            return digest
        if (spec.trigger in ("manual", "idle") and
                run.no_progress >= spec.plateau):
            self.stop(f"plateau: {run.no_progress} iterations without "
                      "measured progress")
            return digest
        if run.is_loop and spec.per_iteration == "one" and \
                run.step_index >= len(spec.steps):
            self.stop("every step done")
            return digest
        why = self._budget_stop(engine)
        if why:
            self.stop(why)
        return digest

    # -- self-found work -----------------------------------------------------
    def find_work(self, engine) -> list[WorkItem]:
        """Ranked candidates; see the module docstring. Never raises."""
        out: list[WorkItem] = []
        sid = str(getattr(engine, "session_id", "") or self.session_id)
        try:
            from . import agent_tasks
            summ = agent_tasks.open_task_summary(self.workspace, sid or None)
            for t in summ.get("in_progress", []):
                out.append(WorkItem(2, "task", f"task {t.get('seq') or t.get('id')}"
                                    f" in progress: {t.get('subject')}", True,
                                    str(t.get("id"))))
            for t in summ.get("pending", []):
                out.append(WorkItem(5, "task", f"task {t.get('seq') or t.get('id')}"
                                    f" pending: {t.get('subject')}", True,
                                    str(t.get("id"))))
            if sid:
                foreign = agent_tasks.open_foreign_tasks(self.workspace, sid)
                for t in foreign.get("tasks", []):
                    out.append(WorkItem(8, "foreign", f"another session's task "
                                        f"{t.get('id')}: {t.get('subject')}",
                                        False, str(t.get("id"))))
        except Exception:
            pass
        try:
            note = str(getattr(engine, "pending_turn_continuation", "") or "")
            if note:
                out.append(WorkItem(3, "announced", "the step announced in the "
                                    "last answer and not taken", True))
        except Exception:
            pass
        last = self._last_test_verdict(engine)
        if last and re.search(r"\b\d+ (failed|errors?)\b", last, re.I):
            out.append(WorkItem(4, "tests", "the last test run failed: " + last,
                                True))
        red = self._ci_red()
        if red:
            out.append(WorkItem(4, "ci", red, True))
        for rel in self.settings.get("registers") or []:
            p = self.workspace / str(rel)
            for m in list(_UNCHECKED.finditer(_read(p)))[:10]:
                out.append(WorkItem(7, "register", f"{rel}: {m.group(1).strip()}",
                                    False))
        out.sort(key=lambda w: w.rank)
        return out

    @staticmethod
    def _last_test_verdict(engine) -> str:
        """The newest test verdict line ("2 failed, 5 passed") in the
        session's last 40 messages, any role; "" when there is none."""
        try:
            from .working_state import _TEST_VERDICT_RE
            msgs = list(getattr(engine, "messages", []) or [])[-40:]
            for msg in reversed(msgs):
                text = _any_text(msg.get("content") if isinstance(msg, dict)
                                 else "")
                for line in reversed(text.splitlines()):
                    if _TEST_VERDICT_RE.search(line):
                        return " ".join(line.split())[:160]
        except Exception:
            return ""
        return ""

    _ci_cache: tuple[float, str] = (0.0, "")

    def _ci_red(self) -> str:
        """Red CI on the current (non-default) branch, via gh, cached ten
        minutes; "" when there is none or it cannot be told."""
        try:
            branch = self.git.branch()
            if not branch or branch in self.git.defaults():
                return ""
            at, value = self._ci_cache
            if self.clock() - at < 600:
                return value
            if self.ci_probe is not None:
                value = str(self.ci_probe(branch) or "")
            elif shutil.which("gh"):
                r = subprocess.run(
                    ["gh", "run", "list", "--branch", branch, "--limit", "1",
                     "--json", "conclusion,displayTitle,url"],
                    cwd=str(self.workspace), capture_output=True, text=True,
                    timeout=10, stdin=subprocess.DEVNULL)
                rows = json.loads(r.stdout or "[]") if r.returncode == 0 else []
                value = ""
                if rows and rows[0].get("conclusion") == "failure":
                    value = (f"CI is red on {branch}: "
                             f"{rows[0].get('displayTitle', '')} {rows[0].get('url', '')}")
            else:
                value = ""
            self._ci_cache = (self.clock(), value)
            return value
        except Exception:
            return ""

    # -- visibility ----------------------------------------------------------
    def state(self) -> dict:
        run = self._run
        st = {"enabled": self.enabled, "paused": self.paused,
              "loop": run.spec.name if run and run.is_loop else "",
              "iteration": run.iteration if run else 0,
              "max_iterations": run.spec.max_iterations if run else 0,
              "last": run.last_result if run else "",
              "running": self.iteration is not None,
              "needs_you": [dict(n) for n in self.needs_you],
              "stop_report": self.stop_report,
              "branch": self.git.branch() if self.git.is_repo() else ""}
        if run:
            st["tokens_left"] = max(0, run.spec.max_tokens - self._tokens_used_cache)
            st["wall_left_s"] = max(0, int(run.spec.wall_seconds
                                           - (self.clock() - run.started_at)))
        else:
            st["tokens_left"] = 0
            st["wall_left_s"] = 0
        return st

    def status_line(self) -> str:
        """One line for a mode line; "" while auto mode is off."""
        run = self._run
        if not self.enabled or run is None:
            return ""
        name = f"loop {run.spec.name}" if run.is_loop else "auto"
        parts = [f"{name} {run.iteration}/{run.spec.max_iterations}"]
        if self.iteration is not None:
            parts.append("running")
        elif self.paused:
            parts.append("paused")
        if run.last_result:
            parts.append(run.last_result)
        if self.needs_you:
            parts.append(f"{len(self.needs_you)} need you")
        return " · ".join(parts)
