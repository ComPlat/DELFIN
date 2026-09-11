"""Owner-only creation, repair and retention for the agent's state tree.

Every recorder under ``~/.delfin`` used to create its directory with a bare
``mkdir`` and no mode, so the tree landed at whatever the process umask
allowed.  Observed on a shared machine (umask 002): ``~/.delfin`` itself at
0700 — inherited from an older umask, set by no code — and every
subdirectory it created at 0775 with most files at 0664.  Sub-agent session
records, which carry up to 60 tool outputs plus the delegate's full report,
had no ``chmod`` at all.

Two separate defects, and the second is why write-time ``chmod`` alone was
not enough:

1. CREATION.  A new directory/file must be owner-only from the start, on
   any umask.  :func:`ensure_dir` and :func:`secure_file` are the single
   pair every state writer calls.
2. LEGACY.  The per-file ``chmod`` calls added earlier only ever touch the
   file being written *now*.  A trace from April keeps its 0664 forever,
   because nothing appends to it again.  :func:`repair_tree` is the
   one-shot sweep that fixes what is already on disk.

Retention is the third piece, and it is what turns a one-day permission
slip into a four-month archive: without it the tree only grows.  Nothing
expired — ``cleanup_old_sessions`` existed with no production caller, and
sub-agent sessions, running markers, tool traces, turn metrics, transcript
archives, handoffs and bundles had no prune at all.

:func:`run_startup_maintenance` does the repair and the prune together,
once per process, and is deliberately total: it never raises, and a failure
on one directory does not stop the next.
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Any, Callable, Iterable, NamedTuple

# Owner-only, on any umask: neither value carries group/other bits, so a
# umask can only ever leave them as they are.
DIR_MODE = 0o700
FILE_MODE = 0o600

# Bits that must never be set on anything in the state tree.
_SHARED_BITS = 0o077

# A walk over a state tree is bounded so a pathological directory (a user
# pointing the archive at their home, a runaway trace dir) cannot turn
# start-up into a multi-minute stat storm.
_MAX_SWEEP_PATHS = 20_000

# Defaults, overridable from settings (see :func:`_retention_settings`).
#
# Two tiers on purpose. Derived state — traces, metrics, markers, archived
# transcripts, bundles, reports — is telemetry about work that is already
# finished, and a month of it is generous. Saved sessions are RESUMABLE
# user work: deleting one silently costs the user something they might
# still want, so it keeps a longer default and its own setting.
DEFAULT_RETENTION_DAYS = 30
DEFAULT_SESSION_RETENTION_DAYS = 90


class StateDir(NamedTuple):
    """One target in the agent's state tree.

    ``kind`` picks what happens to the target:

    ``"session"``  resumable user work — pruned on the longer tier
    ``"derived"``  telemetry — pruned on the short tier, modes repaired
    ``"keep"``     modes repaired, never pruned (a running record)
    ``"prune"``    pruned, modes left alone

    ``"prune"`` exists for the bug archive. A filed report is made
    group-readable ON PURPOSE, so that a maintainer can read it; a sweep
    that tightened it back would fight the writer and leave the mode
    depending on whether a report was written before or after the last
    start-up. Its growth still needs bounding.

    ``glob`` bounds the target to matching files directly inside ``path``
    instead of the whole subtree.  That is what lets the rotated audit logs
    be swept without touching the directory they share — ``~/.delfin``
    itself, which also holds the settings file, the credential store and
    the memory store.  Recursing there would be wrong for a chmod and
    catastrophic for a prune.
    """

    path: Path
    kind: str
    glob: str | None = None


# ---------------------------------------------------------------------------
# Creation
# ---------------------------------------------------------------------------

def _chmod(path: Path, mode: int) -> bool:
    """Best-effort ``chmod``. True when the mode was changed.

    Tightening permissions must never break the thing it protects: a state
    tree on a filesystem that rejects ``chmod`` (some NFS/CIFS mounts) has
    to keep working, just less privately.
    """
    try:
        current = os.stat(path).st_mode & 0o7777
    except OSError:
        return False
    if current == mode:
        return False
    try:
        os.chmod(path, mode)
    except OSError:
        return False
    return True


def ensure_dir(path: Any, *, mode: int = DIR_MODE) -> Path:
    """Create ``path`` (with parents) and make it owner-only.

    The ``chmod`` is unconditional rather than creation-only: every one of
    these directories already exists on installed machines, created by a
    bare ``mkdir`` under the process umask, so a mode applied only to fresh
    directories would fix nothing where it matters.  Never raises.
    """
    p = Path(path)
    try:
        p.mkdir(parents=True, exist_ok=True, mode=mode)
    except OSError:
        return p
    _chmod(p, mode)
    return p


def secure_file(path: Any, *, mode: int = FILE_MODE) -> Path:
    """Make an existing state file owner-only. Never raises."""
    p = Path(path)
    _chmod(p, mode)
    return p


def write_text(path: Any, text: str, *, encoding: str = "utf-8") -> Path:
    """Write a state file owner-only from the moment it exists.

    The file is created with 0600 by ``os.open`` before any content is
    written, so there is no window in which it sits on disk group-readable.
    """
    p = Path(path)
    ensure_dir(p.parent)
    fd = os.open(str(p), os.O_WRONLY | os.O_CREAT | os.O_TRUNC, FILE_MODE)
    try:
        with os.fdopen(fd, "w", encoding=encoding) as fh:
            fh.write(text)
    except Exception:
        try:
            os.close(fd)
        except OSError:
            pass
        raise
    secure_file(p)
    return p


def open_append(path: Any, *, encoding: str = "utf-8"):
    """Open a state file for appending, owner-only from creation.

    ``chmod`` after the first write leaves the file group-readable for the
    duration of that write; passing the mode to ``os.open`` closes that
    window.  On an existing file ``os.open`` ignores the mode, so the
    caller still tightens it afterwards via :func:`secure_file`.
    """
    p = Path(path)
    fd = os.open(str(p), os.O_WRONLY | os.O_CREAT | os.O_APPEND, FILE_MODE)
    return os.fdopen(fd, "a", encoding=encoding)


# ---------------------------------------------------------------------------
# Repair — what is already on disk
# ---------------------------------------------------------------------------

def _iter_files(root: Path, glob: str | None, max_paths: int):
    """Yield the files a sweep may touch under ``root``.

    With ``glob`` the walk is non-recursive and limited to matching names:
    that is the mode used for a directory the agent shares with state it
    does not own.  Symlinks are never yielded — following one would let a
    link planted inside the tree redirect a ``chmod`` or an ``unlink`` at
    an arbitrary file.
    """
    seen = 0
    if glob is not None:
        try:
            entries = sorted(root.glob(glob))
        except OSError:
            return
        for p in entries:
            seen += 1
            if seen > max_paths:
                return
            try:
                if p.is_symlink() or not p.is_file():
                    continue
            except OSError:
                continue
            yield p
        return
    try:
        for dirpath, _dirnames, filenames in os.walk(str(root),
                                                     followlinks=False):
            for name in filenames:
                seen += 1
                if seen > max_paths:
                    return
                p = Path(dirpath) / name
                try:
                    if p.is_symlink():
                        continue
                except OSError:
                    continue
                yield p
    except OSError:
        return


def repair_tree(
    root: Any,
    *,
    dir_mode: int = DIR_MODE,
    file_mode: int = FILE_MODE,
    glob: str | None = None,
    max_paths: int = _MAX_SWEEP_PATHS,
) -> int:
    """Tighten every path under ``root`` that is readable beyond its owner.

    Idempotent: a second run over a repaired tree changes nothing and
    returns 0.  Only paths that actually carry group/other bits are
    touched, so the sweep costs one ``stat`` per path and a ``chmod`` only
    where it is needed.

    ``glob`` restricts the sweep to matching files directly inside ``root``
    and leaves ``root`` itself alone — for a directory the agent shares
    with state it does not own.  Returns the number of paths whose mode was
    changed.  Never raises.
    """
    root_path = Path(root)
    changed = 0
    try:
        if not root_path.is_dir():
            return 0
    except OSError:
        return 0
    if glob is None and _chmod(root_path, dir_mode):
        changed += 1
    if glob is None:
        seen = 0
        try:
            for dirpath, dirnames, _files in os.walk(str(root_path),
                                                     followlinks=False):
                for name in dirnames:
                    seen += 1
                    if seen > max_paths:
                        break
                    p = Path(dirpath) / name
                    try:
                        if p.is_symlink():
                            continue
                    except OSError:
                        continue
                    if _needs_tightening(p) and _chmod(p, dir_mode):
                        changed += 1
        except OSError:
            pass
    for p in _iter_files(root_path, glob, max_paths):
        if _needs_tightening(p) and _chmod(p, file_mode):
            changed += 1
    return changed


def _needs_tightening(path: Path) -> bool:
    try:
        return bool(os.stat(path).st_mode & _SHARED_BITS)
    except OSError:
        return False


# ---------------------------------------------------------------------------
# Retention
# ---------------------------------------------------------------------------

def prune_old(
    root: Any,
    *,
    max_age_days: float,
    glob: str | None = None,
    now: float | None = None,
    on_remove: Callable[[Path], None] | None = None,
    max_paths: int = _MAX_SWEEP_PATHS,
) -> int:
    """Delete files under ``root`` last modified more than N days ago.

    ``max_age_days <= 0`` disables pruning entirely and returns 0 — that
    is the escape hatch for anyone who wants the archive to keep growing.
    ``glob`` bounds the prune to matching files directly inside ``root``.
    ``on_remove`` is called with each path just before it is unlinked, so a
    caller can drop companion state (a job's output tempfiles, say).

    ``root`` itself is never removed — an empty ``tool_traces/`` costs
    nothing and removing it would race the writer about to use it — but
    sub-directories left empty ARE, because the bug archive stores one
    directory per report and pruning only its files would leave an
    ever-growing tree of empty shells.  Returns the number of files
    removed.  Never raises.
    """
    try:
        days = float(max_age_days)
    except (TypeError, ValueError):
        return 0
    if days <= 0:
        return 0
    root_path = Path(root)
    try:
        if not root_path.is_dir():
            return 0
    except OSError:
        return 0
    cutoff = (now if now is not None else time.time()) - days * 86400.0
    removed = 0
    for p in _iter_files(root_path, glob, max_paths):
        try:
            if p.stat().st_mtime >= cutoff:
                continue
        except OSError:
            continue
        if on_remove is not None:
            try:
                on_remove(p)
            except Exception:
                pass
        try:
            p.unlink()
            removed += 1
        except OSError:
            pass
    if removed and glob is None:
        _drop_empty_subdirs(root_path)
    return removed


def _drop_empty_subdirs(root: Path) -> None:
    """Remove directories left empty under ``root``, deepest first.

    ``root`` itself always survives.  Best-effort: ``rmdir`` fails
    harmlessly on a directory a concurrent writer has just refilled.
    """
    try:
        walked = list(os.walk(str(root), topdown=False, followlinks=False))
    except OSError:
        return
    for dirpath, _dirnames, _filenames in walked:
        if Path(dirpath) == root:
            continue
        try:
            os.rmdir(dirpath)
        except OSError:
            continue


# ---------------------------------------------------------------------------
# The tree itself
# ---------------------------------------------------------------------------

# Each entry names the module and the attribute holding a state directory.
# Resolved LIVE on every call, never cached: the attributes are module-level
# constants, and the test suite redirects them per test. A snapshot taken at
# import time would point the sweep at the user's real home from inside a
# test run — which for a prune means deleting their sessions.
# ``file`` marks an attribute that resolves to a FILE, not a directory: the
# entry then targets that file's directory, bounded by ``glob``.
_STATE_DIR_SOURCES: tuple[tuple[str, str, str, str | None, bool], ...] = (
    ("delfin.agent.session_store", "_SESSIONS_DIR", "session", None, False),
    ("delfin.agent.session_store", "_transcript_archive_path", "derived",
     None, False),
    ("delfin.agent.session_store", "_handoffs_path", "derived", None, False),
    ("delfin.agent.session_store", "_bundles_path", "derived", None, False),
    ("delfin.agent.tool_trace", "_DIR", "derived", None, False),
    ("delfin.agent.turn_metrics", "_DIR", "derived", None, False),
    ("delfin.agent.subagents", "_SESSIONS_DIR", "derived", None, False),
    ("delfin.agent.subagents", "_RUNNING_DIR", "derived", None, False),
    ("delfin.agent.subagents", "_PENDING_DIR", "derived", None, False),
    ("delfin.agent.bug_report", "_FALLBACK_DIR", "prune", None, False),
    # Rotated audit logs: eight were observed at 0664. They sit directly in
    # ``~/.delfin`` next to the settings file and the credential store, so
    # the glob is what keeps the sweep off everything else in there.
    ("delfin.agent.audit_log", "_default_log_path", "derived",
     "audit-*.log", True),
    # The active audit log and the failure log are repaired, never pruned:
    # they are the running record, and both already self-trim.
    ("delfin.agent.audit_log", "_default_log_path", "keep",
     "audit.log", True),
    ("delfin.agent.failure_log", "_LOG_PATH", "keep",
     "failure_log.jsonl", True),
)


def state_dirs() -> list[StateDir]:
    """The agent's state targets, resolved from their owning modules.

    A module that cannot be imported simply drops out — this must work on a
    partial install.  Only existing directories are returned: the sweep has
    nothing to do for one that was never created.
    """
    import importlib

    out: list[StateDir] = []
    seen: set[tuple[str, str | None]] = set()
    for mod_name, attr, kind, glob, is_file in _STATE_DIR_SOURCES:
        try:
            mod = importlib.import_module(mod_name)
            raw = getattr(mod, attr, None)
        except Exception:
            continue
        if raw is None:
            continue
        if callable(raw):                 # a per-call resolver
            try:
                raw = raw()
            except Exception:
                continue
        try:
            p = Path(raw)
            if is_file:
                p = p.parent
            if not p.is_dir():
                continue
            key = (str(p.resolve()), glob)
        except (OSError, TypeError, ValueError):
            continue
        if key in seen:
            continue
        seen.add(key)
        out.append(StateDir(path=p, kind=kind, glob=glob))
    return out


def _retention_settings() -> tuple[float, float]:
    """``(derived_days, session_days)`` from settings. Never raises."""
    derived: Any = DEFAULT_RETENTION_DAYS
    sessions: Any = DEFAULT_SESSION_RETENTION_DAYS
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent") or {}
        derived = agent.get("state_retention_days", derived)
        sessions = agent.get("session_retention_days", sessions)
    except Exception:
        pass
    try:
        derived_f = float(derived)
    except (TypeError, ValueError):
        derived_f = float(DEFAULT_RETENTION_DAYS)
    try:
        sessions_f = float(sessions)
    except (TypeError, ValueError):
        sessions_f = float(DEFAULT_SESSION_RETENTION_DAYS)
    return derived_f, sessions_f


# One sweep per process. The work is idempotent, but repeating it on every
# session save would put a full tree walk on the hot path.
_MAINTENANCE_DONE = False


def run_startup_maintenance(
    *,
    dirs: Iterable[StateDir] | None = None,
    retention_days: float | None = None,
    session_retention_days: float | None = None,
    prune: bool = True,
    force: bool = False,
    now: float | None = None,
) -> dict:
    """Repair modes and apply retention across the state tree, once.

    Returns ``{"repaired": n, "pruned": n, "dirs": n, "ran": bool}``.
    ``ran`` is False when a previous call in this process already did the
    work — the caller can therefore be a hot path.

    Never raises: start-up maintenance that can abort start-up is worse
    than the permissions it fixes.
    """
    global _MAINTENANCE_DONE
    if _MAINTENANCE_DONE and not force:
        return {"repaired": 0, "pruned": 0, "dirs": 0, "ran": False}
    _MAINTENANCE_DONE = True

    try:
        targets = list(dirs) if dirs is not None else state_dirs()
    except Exception:
        return {"repaired": 0, "pruned": 0, "dirs": 0, "ran": True}

    if retention_days is None or session_retention_days is None:
        cfg_derived, cfg_sessions = _retention_settings()
        if retention_days is None:
            retention_days = cfg_derived
        if session_retention_days is None:
            session_retention_days = cfg_sessions

    repaired = 0
    pruned = 0
    for entry in targets:
        try:
            path, kind, glob = entry.path, entry.kind, entry.glob
        except AttributeError:            # a bare path is accepted too
            path, kind, glob = Path(entry), "derived", None  # type: ignore
        if prune and kind != "keep":
            days = (session_retention_days if kind == "session"
                    else retention_days)
            try:
                pruned += prune_old(path, max_age_days=days, glob=glob,
                                    now=now)
            except Exception:
                pass
        # Repair AFTER the prune: no point chmod'ing a file about to go.
        if kind != "prune":
            try:
                repaired += repair_tree(path, glob=glob)
            except Exception:
                pass
    return {"repaired": repaired, "pruned": pruned,
            "dirs": len(targets), "ran": True}


def reset_maintenance_flag() -> None:
    """Let the one-shot run again (tests; a long-lived daemon)."""
    global _MAINTENANCE_DONE
    _MAINTENANCE_DONE = False


# ---------------------------------------------------------------------------
# A run that is nobody's
# ---------------------------------------------------------------------------
#
# The agent reads its own past into every prompt: the session briefing
# quotes recent failures from ``~/.delfin/outcome_history.jsonl``, recall
# draws on the memory stores, the provider profile carries what earlier
# runs learned. That is right for a user. It is wrong for a benchmark
# attempt, a probe, or a suite: those must see nothing that an earlier
# run left, and leave nothing behind for a later one.
#
# Observed 2026-09-10. A benchmark task asking which of the user's
# calculations has the lowest energy answered about "the three runs in
# runs/" -- the wording of a DIFFERENT task, whose failed attempts sat in
# the user's history and reached the prompt as "failure lessons". Of 179
# records in that history, 59 were benchmark prompts verbatim and most of
# the rest were probes; the user's own work was a handful of lines.
#
# The table below is every module-level sink under ``~/.delfin`` that a
# run writes to. Each entry is resolved at import time in its module, so
# redirecting means setting the attribute -- which is also what the test
# suite does, per test, from this same table.
#
# Read-mostly stores are deliberately absent: the documentation index, the
# model-capability cache, credentials (isolated separately) and the
# settings files. Redirecting those would not stop a write, it would only
# hide the real content from a run that legitimately reads it.

#: ``(module, attribute, relative path under the redirected home)``.
USER_STATE_SINKS: tuple[tuple[str, str, str], ...] = (
    ("delfin.agent.bash_jobs", "_INDEX_PATH", "bash_jobs_index.json"),
    # A kept dashboard session announces itself so a later request can
    # find its kernel. A test that armed one wrote into the user's real
    # ~/.delfin and their next landing page then offered a session that
    # was a test fixture.
    ("delfin.dashboard.session", "RECORD_DIR", "kept_sessions"),
    # The benchmark's per-checkout run lock. It lived under tests/fixtures
    # first, where the checkout-leak guard would have caught it, and the
    # move to ~/.delfin brought it into this table's scope instead.
    ("delfin.agent.benchmark_runner", "_RUN_LOCK_DIR", "benchmark_locks"),
    ("delfin.agent.provider_profile", "_LOCAL_STATE_PATH",
     "provider_profile_state.json"),
    ("delfin.agent.job_fix", "_ATTEMPTS_PATH", "fix_attempts.json"),
    ("delfin.agent.session_store", "_SESSIONS_DIR", "agent_sessions"),
    ("delfin.agent.outcome_tracker", "_DEFAULT_PATH", "outcome_history.jsonl"),
    ("delfin.agent.agent_metrics", "_LOG_PATH", "agent_metrics.jsonl"),
    ("delfin.agent.context_tracker", "_DEFAULT_PATH", "context_usage.jsonl"),
    ("delfin.agent.failure_log", "_LOG_PATH", "failure_log.jsonl"),
    ("delfin.agent.eval_loop", "_REPORTS_DIR", "eval_reports"),
    ("delfin.agent.eval_loop", "_TASK_DRAFTS_DIR", "bug_tasks"),
    ("delfin.agent.bug_report", "_FALLBACK_DIR", "agent_bugs"),
    ("delfin.agent.bug_report", "_TASK_DRAFTS_DIR", "bug_tasks"),
    ("delfin.agent.benchmark", "_DEFAULT_RUNS_DIR", "benchmark_runs"),
    ("delfin.agent.scheduler", "_DEFAULT_PATH", "cron.json"),
    ("delfin.dashboard.schedules", "_DEFAULT_PATH", "schedules.json"),
    ("delfin.agent.memory_store", "_DEFAULT_PATH", "agent_memory.json"),
    ("delfin.agent.skill_registry", "_LOCAL_SKILLS_DIR", "skills"),
    ("delfin.agent.job_monitor", "_WATCHED_PATH", "watched_jobs.json"),
    ("delfin.agent.job_monitor", "_AGENT_WATCH_INDEX_PATH",
     "agent_watch_index.json"),
    ("delfin.agent.job_monitor", "_FINDINGS_PATH", "monitor_findings.jsonl"),
    ("delfin.agent.job_monitor", "_PID_PATH", "job_monitor.pid"),
    ("delfin.agent.bug_watcher", "_PID_PATH", "bug_watcher.pid"),
    ("delfin.agent.scheduler_daemon", "_PID_PATH", "scheduler_daemon.pid"),
    # The user's own settings file. A permission rule the agent persists on
    # approval is written here, so a test exercising that path edited the
    # real file -- and permission rules are exactly what must not be
    # granted by accident.
    ("delfin.agent.hooks_editor", "_USER_SETTINGS", "settings.json"),
    ("delfin.agent.kit_settings", "USER_SETTINGS_PATH", "settings.json"),
)

#: Sinks resolved per call rather than at import: ``(module, function,
#: relative path)``. The two project-store resolvers take the repo root
#: and land under ``projects/<slug>/<leaf>``; the rest take no argument.
USER_STATE_RESOLVERS: tuple[tuple[str, str, str], ...] = (
    ("delfin.agent.audit_log", "_default_log_path", "audit.log"),
    # Which directories the user has trusted to run commands. A run that
    # granted trust must never grant it in the real store: the entry would
    # outlive the run and let a later, real session honour a workspace's
    # hooks and MCP servers on the strength of a fixture directory.
    ("delfin.agent.workspace_trust", "_trust_store_path",
     "trusted_workspaces.json"),
    # The state-tree maintenance sweep walks these three, and its prune
    # DELETES: pointed at the real home from inside a test run it would
    # remove the user's archived transcripts, handoffs and bundles.
    ("delfin.agent.session_store", "_transcript_archive_path",
     "transcript_archive"),
    ("delfin.agent.session_store", "_handoffs_path", "handoffs"),
    ("delfin.agent.session_store", "_bundles_path", "bundles"),
    ("delfin.agent.attention", "_inbox_path", "attention_inbox.jsonl"),
    ("delfin.agent.change_journal", "_undo_root", "undo"),
    ("delfin.agent.memory_store", "_delfin_plans_dir", "projects"),
    ("delfin.agent.memory_store", "_delfin_memory_dir", "projects"),
    ("delfin.agent.memory_store", "_delfin_global_memory_dir", "memory"),
)

#: The resolvers that take a repo root and add a per-project leaf.
PROJECT_LEAVES: dict[str, str] = {
    "_delfin_memory_dir": "memory",
    "_delfin_plans_dir": "plans",
}

#: What a LIVE run keeps real even while everything else is redirected.
#: A benchmark attempt runs with the user's configuration -- their
#: settings, their provider -- it just does not run with their past. The
#: bench's own results directory and run lock are the bench's, not the
#: attempt's: both are written outside the guard.
#:
#: The audit log stays real too, for the opposite reason: it is where
#: the bench READS after the attempt. The gate denials an attempt cost
#: and the paths it wrote are both counted from the audit log once the
#: guard has closed -- and for one night they were counted from a
#: scratch log that no longer existed, so every block reported "denials
#: not observed" and the cost axis went dark (2026-09-11). The log does
#: not reach a prompt on its own, and its records carry the workspace
#: they were made in, so an attempt's lines never read as the user's.
KEPT_BY_A_LIVE_RUN: frozenset[tuple[str, str]] = frozenset({
    ("delfin.agent.hooks_editor", "_USER_SETTINGS"),
    ("delfin.agent.kit_settings", "USER_SETTINGS_PATH"),
    ("delfin.agent.benchmark", "_DEFAULT_RUNS_DIR"),
    ("delfin.agent.benchmark_runner", "_RUN_LOCK_DIR"),
    ("delfin.agent.audit_log", "_default_log_path"),
})

#: Set this to a directory and a ``delfin agent`` process keeps its
#: history, memory, profile learning, sessions and logs there instead of
#: under ``~/.delfin``. For probes and interviews driven through the CLI,
#: which otherwise write into the user's history exactly as a benchmark
#: attempt used to.
SCRATCH_STATE_ENV = "DELFIN_SCRATCH_STATE"


class RedirectedUserState:
    """Point every writable sink at ``home/.delfin`` for the duration.

    A context manager, re-entrant across separate instances: each one
    remembers what it replaced and puts exactly that back, so a guard
    inside a suite that already redirected restores the suite's redirect,
    not the real home.

    Nothing on disk is touched. The attributes are swapped in memory;
    what a run writes lands under ``home`` and what it reads is whatever
    ``home`` holds -- for a fresh directory, nothing. A process killed
    mid-run therefore leaves the user's files as they were, which the
    snapshot-and-restore guards this replaced could not promise.
    """

    def __init__(self, home: "Path | str", *,
                 keep: Iterable[tuple[str, str]] = ()) -> None:
        self.home = Path(home)
        self.root = self.home / ".delfin"
        self._keep = frozenset(keep)
        self._saved: list[tuple[Any, str, Any]] = []
        self.redirected: dict[tuple[str, str], Path] = {}

    def __enter__(self) -> "RedirectedUserState":
        import importlib

        for mod_name, attr, rel in USER_STATE_SINKS:
            if (mod_name, attr) in self._keep:
                continue
            try:
                mod = importlib.import_module(mod_name)
            except Exception:
                continue
            if not hasattr(mod, attr):
                continue
            self._saved.append((mod, attr, getattr(mod, attr)))
            target = self.root / rel
            setattr(mod, attr, target)
            self.redirected[(mod_name, attr)] = target

        for mod_name, attr, rel in USER_STATE_RESOLVERS:
            if (mod_name, attr) in self._keep:
                continue
            try:
                mod = importlib.import_module(mod_name)
            except Exception:
                continue
            original = getattr(mod, attr, None)
            if not callable(original):
                continue
            self._saved.append((mod, attr, original))
            leaf = PROJECT_LEAVES.get(attr)
            if leaf is not None:
                def _resolve(repo_root, _root=self.root, _leaf=leaf):
                    from delfin.agent.memory_store import _project_slug
                    return (_root / "projects" / _project_slug(repo_root)
                            / _leaf)
            else:
                def _resolve(_target=self.root / rel):
                    return _target
            setattr(mod, attr, _resolve)
            self.redirected[(mod_name, attr)] = self.root / rel
        return self

    def __exit__(self, *exc) -> None:
        # Reverse order, so an attribute redirected twice by mistake ends
        # up with its first original.
        for mod, attr, original in reversed(self._saved):
            try:
                setattr(mod, attr, original)
            except Exception:
                pass
        self._saved = []


def scratch_state_from_environment(
        environ: "dict | None" = None) -> "RedirectedUserState | None":
    """Honour :data:`SCRATCH_STATE_ENV` for the life of this process.

    Returns the redirect so a caller that wants to end it can; a CLI
    simply lets it stand. An unset or blank variable does nothing, and
    the directory is created so the first write does not have to.
    """
    env = os.environ if environ is None else environ
    raw = (env.get(SCRATCH_STATE_ENV) or "").strip()
    if not raw:
        return None
    home = Path(raw).expanduser()
    try:
        ensure_dir(home / ".delfin")
    except Exception:
        pass
    redirect = RedirectedUserState(home, keep=KEPT_BY_A_LIVE_RUN)
    redirect.__enter__()
    return redirect


__all__ = [
    "DIR_MODE", "FILE_MODE",
    "DEFAULT_RETENTION_DAYS", "DEFAULT_SESSION_RETENTION_DAYS",
    "StateDir", "ensure_dir", "secure_file", "write_text", "open_append",
    "repair_tree", "prune_old", "state_dirs", "run_startup_maintenance",
    "reset_maintenance_flag",
    "USER_STATE_SINKS", "USER_STATE_RESOLVERS", "PROJECT_LEAVES",
    "KEPT_BY_A_LIVE_RUN", "SCRATCH_STATE_ENV", "RedirectedUserState",
    "scratch_state_from_environment",
]
