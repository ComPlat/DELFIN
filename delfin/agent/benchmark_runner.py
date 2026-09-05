"""Live-engine runner for the canned-task benchmark suite.

The data plane (``benchmark.py``) defines Task / Trajectory / scoring; this
module bridges to the real ``AgentEngine`` and feeds it the prompts.  A
single function — ``run_suite`` — is the entire surface area:

    results = run_suite(load_tasks(), model="kit.qwen3.5-397b-A17b",
                        backend="api", provider="kit")
    path = write_run(results, model="kit.qwen3.5-397b-A17b")

An engine factory can be injected for unit-testing — the real factory
defers all heavy imports until invocation so a bench CLI parse fails
fast on a bad argument before pulling the engine stack.
"""

from __future__ import annotations

from pathlib import Path

import contextlib
import os
import re
import time
from typing import Any, Callable, Iterable, Optional

from .benchmark import (
    BenchmarkResult,
    Task,
    Trajectory,
    score_outcome,
)
from .benchmark_fixtures import ensure_office_fixtures


# ---------------------------------------------------------------------------
# Trajectory builder — convert _run_once output into a Trajectory.
# ---------------------------------------------------------------------------


# Same five ACTION variants tab_agent.py accepts in the dashboard.
_ACTION_PATTERNS: tuple[re.Pattern[str], ...] = (
    re.compile(r"^\s*ACTION\s*:\s*(/\S.*?)\s*$", re.MULTILINE),
    re.compile(r"^\s*ACTION\s+(/\S.*?)\s*$", re.MULTILINE),
    re.compile(r"^\s*Action\s*:\s*(/\S.*?)\s*$", re.MULTILINE),
    re.compile(r"`(/(?:tab|control|orca|effort|mode|provider|model)\s+\S[^`]*)`"),
)
_BARE_SLASH_PREFIXES = (
    "/tab ", "/control ", "/orca ", "/effort ", "/mode ", "/provider ",
    "/model ",
)


def extract_actions(text: str) -> list[str]:
    """Pull ACTION lines (canonical + tolerant variants) from text.

    Mirrors the dashboard's _action_cmd parser so the benchmark sees the
    same routing decisions an end user would.  Bare-slash lines (no
    ``ACTION:`` prefix) are accepted only when the slash command starts
    with a known dashboard prefix — that guards against false positives
    like ``/home/user/file.py``.
    """

    if not text:
        return []
    out: list[str] = []
    for rx in _ACTION_PATTERNS:
        for m in rx.finditer(text):
            cmd = (m.group(1) or "").strip()
            if cmd and cmd not in out:
                out.append(cmd)
    # Bare-slash lines (whitelisted prefixes only)
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith("/"):
            continue
        if any(s.startswith(p) for p in _BARE_SLASH_PREFIXES):
            if s not in out:
                out.append(s)
    return out


def trajectory_from_run(raw: dict, *, duration_s: float, cost_usd: float = 0.0,
                       checkout_root: str = "") -> Trajectory:
    """Convert ``_run_once``'s return dict into a Trajectory.

    ``checkout_root`` is where the repository under test lives. It is
    carried so the scorer can drop that prefix from tool inputs -- see
    ``benchmark._strip_checkout_prefix``. Optional, and an omitted root
    changes nothing, so older callers keep working."""
    text = str(raw.get("text") or "")
    return Trajectory(
        text=text,
        actions=extract_actions(text),
        tool_calls=list(raw.get("tool_calls") or []),
        duration_s=float(duration_s),
        cost_usd=float(cost_usd),
        input_tokens=int(raw.get("input_tokens") or 0),
        output_tokens=int(raw.get("output_tokens") or 0),
        error=str(raw.get("error") or ""),
        checkout_root=str(checkout_root or ""),
    )


# ---------------------------------------------------------------------------
# Engine factory + cost estimator
# ---------------------------------------------------------------------------


# (model, backend, provider, mode, task_class) -> engine. The task class
# is what decides the write boundary; the mode cannot, because behaviour
# and generic-project tasks both run as `solo`. Widened rather than
# guessed at the call site.
EngineFactory = Callable[..., Any]


def _iso(epoch: float) -> str:
    """Audit-log timestamps are ISO strings compared lexicographically."""
    import datetime as _dt
    try:
        return _dt.datetime.fromtimestamp(
            float(epoch), _dt.timezone.utc).isoformat()
    except (OSError, OverflowError, ValueError):
        return ""


def _denials_during(engine: Any, *, since_ts: str) -> Optional[int]:
    """How many tool calls the gates refused in this task's window.

    A cost-side number: guards are supposed to refuse the wrong thing, and
    a rising count over a suite that already passed means they started
    refusing the right thing too. That direction never shows up in a pass
    rate, and it is how a check teaches the model to phrase work so
    nothing keys on it.

    Returns ``None`` when there was no way to look -- no audit log, no
    workspace, a mocked engine in a unit test. Reporting 0 there would say
    "nothing was refused", which is a different sentence from "nobody
    looked", and this framework has already shipped that confusion once.
    """
    try:
        from delfin.agent import audit_log as _audit
        perms = getattr(engine, "kit_permissions", None)
        workspace = getattr(perms, "workspace", None) if perms else None
        if not workspace:
            return None
        report = _audit.build_changes_report(
            since_ts=since_ts or None, workspace=str(workspace))
        denied = report.get("denied")
        # An empty skeleton comes back on any failure, and its window
        # records count is 0 -- indistinguishable from a clean run that
        # simply refused nothing. Only a window that actually saw records
        # can say zero.
        window = report.get("window") or {}
        if not window.get("records") and not denied:
            return None
        return len(denied or [])
    except Exception:
        return None


def workspace_for(root: Path, *, mode: str = "", task_class: str = "") -> Optional[Path]:
    """The fixture workspace a task class writes into, or ``None``.

    Only the classes that WRITE get one. Office already had this, keyed
    on the mode; behaviour and generic-project tasks edit toy files just
    as concretely and had the whole checkout as their workspace, because
    the mode cannot tell them apart — both run as ``solo``. So the key is
    the task CLASS, which is what actually decides where a task belongs.

    Measured over the four packaged task files: office 11, behaviour 12,
    generic_project 8. Those 31 are every task that writes; the other 48
    read and are left alone deliberately — see the note in #111 on why
    moving their root is a separate problem.

    ``None`` when the class writes nothing or the directory is absent: a
    missing fixture must not silently redirect a run somewhere else.
    """
    cls = (task_class or "").strip().lower()
    rel: Optional[str] = None
    if (mode or "") == "office" or cls == "office":
        rel = "office_workspace"
    elif cls.startswith("behavior"):
        rel = "behavior_workspace"
    elif cls == "generic_project":
        rel = "user_project_workspace"
    if rel is None:
        return None
    candidate = Path(root) / "tests" / "fixtures" / rel
    return candidate if candidate.is_dir() else None


def _default_engine_factory(model: str, backend: str, provider: str,
                            mode: str, task_class: str = "") -> Any:
    """Build a real AgentEngine for the given config.

    AgentEngine creates its own client internally; we just hand it the
    resolved tuple and let it own the lifecycle.  Heavy imports happen
    here, not at module-import time.
    """
    import os as _os
    from pathlib import Path as _Path
    from .engine import AgentEngine

    # Office mode is DEFINED by working inside one folder, so a task in
    # that mode has to start there. Pointing it at the repository root
    # instead left the fixtures three levels down, and whether the model
    # found them was luck: one run cited the full path, the next answered
    # "the file is not in the working directory" — and that scored as a
    # failure of the rubric it never got to exercise. A benchmark that
    # intermittently measures file discovery instead of the behaviour
    # under test is worse than no benchmark.
    root = _Path(_os.getcwd()).resolve()
    fixtures = workspace_for(root, mode=mode, task_class=task_class)
    if fixtures is not None:
        root = fixtures

    return AgentEngine(
        repo_dir=root,
        backend=backend or "api",
        provider=provider,
        model=model,
        mode=mode or "solo",
    )


def _cost_delta(before: float, after: float) -> float:
    """Δ cost for a single turn — defends against engines that don't
    expose cost_usd."""
    try:
        return max(0.0, float(after) - float(before))
    except (TypeError, ValueError):
        return 0.0


# ---------------------------------------------------------------------------
# Public runner
# ---------------------------------------------------------------------------


# Behavior tasks edit toy files under this repo-relative dir. Attempts must
# start from a PRISTINE copy: an earlier attempt's edits would leak into the
# next one (observed live: a later repeats run scored against fixtures the
# first run had already modified, confounding the comparison). The runner
# snapshots the dir before each attempt and restores it afterwards — the
# task file's manual git-checkout instructions are thereby mechanised.
_BEHAVIOR_WS_RELS: tuple[Path, ...] = (
    Path("tests") / "fixtures" / "behavior_workspace",
    Path("tests") / "fixtures" / "user_project_workspace",
    # Office tasks read and edit the same small tables; without a restore
    # between attempts a later run scores against rows an earlier one
    # changed.
    Path("tests") / "fixtures" / "office_workspace",
)
_BEHAVIOR_WS_REL = _BEHAVIOR_WS_RELS[0]


def checkout_fingerprint(root: Path | str) -> dict[str, tuple[int, int]]:
    """(size, mtime) for every file under ``root``, excluding the guarded
    fixture workspaces and anything git already ignores by convention.

    Used to answer one question after a run: did the agent change the
    checkout outside the directories it is allowed to change? The guard
    below covers three fixture directories, while the agent's workspace is
    the whole tree -- so a run once rewrote a test file two directories
    away and nothing said so.
    """
    base = Path(root)
    guarded = [base / rel for rel in _BEHAVIOR_WS_RELS]
    # `.claude` holds the worktrees parallel agents check out into. A run
    # with four of them alive reported 4992 changed files and then "… and
    # 4992 more", which is not a warning anybody reads twice — and the
    # walk is O(every file in every worktree) before it says so. Observed
    # live, 2026-08-12.
    skip = {".git", "__pycache__", ".pytest_cache", "node_modules",
            ".claude", ".venv", ".mypy_cache", ".ruff_cache"}
    out: dict[str, tuple[int, int]] = {}
    for path in base.rglob("*"):
        try:
            if not path.is_file():
                continue
            if any(part in skip for part in path.parts):
                continue
            if any(g in path.parents for g in guarded):
                continue
            st = path.stat()
            out[str(path.relative_to(base))] = (st.st_size, int(st.st_mtime))
        except OSError:
            continue
    return out


def changed_outside_workspaces(
    root: Path | str, before: dict[str, tuple[int, int]],
) -> list[str]:
    """Paths that differ from ``before`` -- added, removed or rewritten.

    Reporting, not reverting: restoring the whole checkout would throw
    away whatever the developer had uncommitted, which is a worse failure
    than the one being fixed."""
    after = checkout_fingerprint(root)
    changed = {p for p in set(before) | set(after)
               if before.get(p) != after.get(p)}
    return sorted(changed)


def paths_this_run_wrote(root: Path | str, *,
                         since_ts: str) -> Optional[set[str]]:
    """Paths under ``root`` that this run's own tool trace recorded, or
    ``None`` when there was no way to look.

    Same source as :func:`_denials_during` -- the audit log -- and the
    same distinction: an empty window means nobody looked, and reporting
    an empty SET there would relabel every changed file as somebody
    else's writing. Returned relative to ``root``, because that is how
    the fingerprint is keyed; comparing absolute against relative would
    match nothing and quietly attribute the whole diff elsewhere.
    """
    try:
        from delfin.agent import audit_log as _audit
        base = Path(root).resolve()
        report = _audit.build_changes_report(
            since_ts=since_ts or None, workspace=str(base))
        window = report.get("window") or {}
        written = report.get("files_written") or []
        if not window.get("records") and not written:
            return None
        out: set[str] = set()
        for entry in written:
            raw = str((entry or {}).get("path", "") or "")
            if not raw:
                continue
            try:
                out.add(str(Path(raw).resolve().relative_to(base)))
            except (ValueError, OSError):
                continue        # written outside the checkout
        return out
    except Exception:
        return None


def attributed_changes(
    changed: list[str],
    written: Optional[set[str]],
) -> tuple[Optional[list[str]], Optional[list[str]]]:
    """Split a fingerprint diff into (this run wrote it, nobody here did).

    ``(None, None)`` when ``written`` is ``None``: with no trace to check
    against, neither column can be filled in honestly.
    """
    if written is None:
        return None, None
    mine = [p for p in changed if p in written]
    other = [p for p in changed if p not in written]
    return mine, other


_CHANGE_LIST_CAP = 20


def _capped(paths: list[str]) -> list[str]:
    lines = [f"     {p}" for p in paths[:_CHANGE_LIST_CAP]]
    if len(paths) > _CHANGE_LIST_CAP:
        lines.append(f"     … and {len(paths) - _CHANGE_LIST_CAP} more")
    return lines


def format_checkout_change_report(changed: list[str],
                                  written: Optional[set[str]]) -> str:
    """What the run may say about a checkout that is not as it was found.

    The diff is a fact about the TREE. Turning it into a sentence about
    the RUN was the defect: a suite writing into the same checkout put
    569 files into a warning that named this run as their cause. The two
    groups are printed apart, and the group the run cannot vouch for says
    so instead of being counted as its own.
    """
    if not changed:
        return ""
    mine, other = attributed_changes(changed, written)
    if mine is None:
        return "\n".join([
            f"\n⚠️  {len(changed)} file(s) OUTSIDE the fixture workspaces "
            f"changed while this run was going — the checkout is not as it "
            f"was found. This run's own write trace could not be read, so "
            f"which of them it caused is unknown:",
            *_capped(changed),
            "     review with `git status` before trusting the next run.",
        ])
    lines: list[str] = []
    if mine:
        lines.append(f"\n⚠️  {len(mine)} file(s) this run wrote OUTSIDE the "
                     f"fixture workspaces — the checkout is not as it was "
                     f"found:")
        lines.extend(_capped(mine))
    if other:
        lines.append(f"\nℹ️  {len(other)} file(s) changed during this run "
                     f"that nothing in this run's trace wrote. Another "
                     f"process writing into the same checkout is "
                     f"indistinguishable from this run in a size/mtime "
                     f"diff, so these are not counted against it:")
        lines.extend(_capped(other))
    lines.append("     review with `git status` before trusting the next run.")
    return "\n".join(lines)


_fixture_notice_shown = False


def _prepare_office_fixtures() -> str:
    """Build the workbook fixtures, and say once why if it cannot.

    The workbook tasks read files that are generated rather than
    committed. When they cannot be generated the tasks do not measure a
    model at all — they measure a missing dependency, and score it as a
    failing model. That is the one outcome worth being loud about, so the
    reason is printed rather than logged, once per process.
    """
    global _fixture_notice_shown
    try:
        _written, reason = ensure_office_fixtures()
    except Exception as exc:  # a broken builder must not take a run down
        _written, reason = [], f"the fixture builder raised: {exc}"
    if reason and not _fixture_notice_shown:
        _fixture_notice_shown = True
        print(f"\n⚠️  office workbook fixtures were NOT built: {reason}")
        print("     Tasks that read a .xlsx will fail for that reason and "
              "not because of the model — exclude them from any comparison.")
    return reason


class _PristineWorkspace:
    """Snapshot/restore guard for the fixture workspaces (a dir that does not
    exist under the current working directory is skipped).

    ``failed`` records a snapshot that did not happen. It used to swallow
    the error and set ``_pairs`` empty, so the restore then did nothing at
    all -- silently, at the one moment the caller most needed to hear it.
    """

    def __init__(self, root: Path | None = None) -> None:
        import os as _os
        base = root or Path(_os.getcwd())
        self._bases = [base / rel for rel in _BEHAVIOR_WS_RELS]
        self._pairs: list[tuple[Path, Path]] = []
        self._snap_root: Path | None = None
        self.failed = False

    def __enter__(self) -> "_PristineWorkspace":
        import shutil
        import tempfile
        live = [ws for ws in self._bases if ws.is_dir()]
        if not live:
            return self
        # One retry, then refusal. A file appearing or vanishing mid-copy
        # makes copytree raise, and that really happens: writing new
        # fixtures into the workspace while a run was in progress raced it
        # exactly this way. That is transient and deserves a second try.
        #
        # A snapshot that fails twice does not. The guard exists to keep
        # the user's files; a run that cannot keep them has no business
        # touching them. Carrying on unprotected and silent is what
        # happened before, and an attempt then deleted a fixture that
        # nothing put back.
        last: Exception | None = None
        for _attempt in (1, 2):
            self._pairs = []
            try:
                self._snap_root = Path(tempfile.mkdtemp(prefix="bench-ws-"))
                for i, ws in enumerate(live):
                    snap = self._snap_root / f"ws{i}"
                    shutil.copytree(ws, snap)
                    self._pairs.append((ws, snap))
                self.failed = False
                return self
            except Exception as exc:
                last = exc
                if self._snap_root is not None:
                    shutil.rmtree(self._snap_root, ignore_errors=True)
                    self._snap_root = None
        self._pairs = []
        self.failed = True
        names = ", ".join(ws.name for ws in live)
        raise RuntimeError(
            f"could not snapshot the fixture workspace(s) {names}: {last}. "
            f"Refusing to run: without the snapshot nothing can restore "
            f"them afterwards, and an attempt that deletes a fixture would "
            f"leave it deleted.")

    def __exit__(self, *exc) -> None:
        import shutil
        try:
            for ws, snap in self._pairs:
                shutil.rmtree(ws, ignore_errors=True)
                shutil.copytree(snap, ws)
        finally:
            if self._snap_root is not None:
                shutil.rmtree(self._snap_root, ignore_errors=True)


def _run_task_once(
    task: Task,
    *,
    model: str,
    backend: str,
    provider: str,
    profile_name: str,
    engine_factory: EngineFactory,
    max_tokens: int,
    run_once: Callable[..., dict],
    clock: Callable[[], float],
) -> BenchmarkResult:
    """Single attempt — kept private so retry-aggregation logic lives
    in one place at the public ``run_task`` entry."""
    try:
        engine = engine_factory(model, backend, provider, task.mode,
                                task.task_class)
    except Exception as exc:
        traj = Trajectory(error=f"engine init failed: {exc}")
        return score_outcome(
            task, traj, model=model, profile_name=profile_name, ts=clock(),
        )
    # Run budget from the task caps: a benchmark task must not burn an
    # unbounded multiple of its declared cost/duration budget (observed:
    # one task at 27x its cost cap). Generous factors keep legitimate
    # thoroughness possible; the engine's wind-down asks the agent to
    # wrap up at 80% and refuses new turns past 110%.
    try:
        if float(task.max_cost_usd or 0) > 0:
            engine.run_budget_usd = max(1.0, float(task.max_cost_usd) * 4.0)
        if float(task.max_duration_s or 0) > 0:
            engine.run_budget_s = max(120.0, float(task.max_duration_s) * 4.0)
    except Exception:
        pass
    cost_before = float(getattr(engine, "cost_usd", 0.0) or 0.0)
    t0 = clock()
    try:
        with _PristineWorkspace():
            raw = run_once(engine, task.prompt, max_tokens=max_tokens)
    except Exception as exc:
        raw = {"text": "", "tool_calls": [], "input_tokens": 0,
               "output_tokens": 0, "error": f"_run_once raised: {exc}"}
    t1 = clock()
    cost_after = float(getattr(engine, "cost_usd", 0.0) or 0.0)
    traj = trajectory_from_run(
        raw, duration_s=t1 - t0, cost_usd=_cost_delta(cost_before, cost_after),
        checkout_root=str(Path(os.getcwd()).resolve()),
    )
    traj.denials = _denials_during(engine, since_ts=_iso(t0))
    return score_outcome(
        task, traj, model=model, profile_name=profile_name, ts=clock(),
    )


def run_task(
    task: Task,
    *,
    model: str,
    backend: str = "api",
    provider: str = "",
    profile_name: str = "",
    engine_factory: EngineFactory | None = None,
    max_tokens: int = 1024,
    repeats: int = 1,
    run_once: Callable[..., dict] | None = None,
    clock: Callable[[], float] | None = None,
    on_replicate: Callable[[int, BenchmarkResult], None] | None = None,
) -> BenchmarkResult:
    """Execute one task and return a scored result.

    ``repeats=1`` (default): single sample, returned as-is.
    ``repeats>1``: run the task N times against a FRESH engine each time
    (so prior-turn state doesn't bias the result) and return a single
    aggregated BenchmarkResult with median-quality + majority-success
    + quality_stdev as a noise indicator.  ``on_replicate`` (if given)
    is called after each individual run with ``(idx, result)`` so a
    caller can stream progress.

    The defaults wire up the real ``AgentEngine`` + ``_run_once`` from
    ``cli.py``; tests can inject stubs.
    """
    now = clock or time.time
    factory = engine_factory or _default_engine_factory
    if run_once is None:
        from .cli import _run_once as _real_run_once
        run_once = _real_run_once

    # Before the snapshot, never after: the guard restores the workspace
    # to whatever it held when it was taken, so a workbook written later
    # is deleted again after the first attempt.
    _prepare_office_fixtures()

    n = max(1, int(repeats))
    replicates: list[BenchmarkResult] = []
    for i in range(n):
        r = _run_task_once(
            task,
            model=model, backend=backend, provider=provider,
            profile_name=profile_name, engine_factory=factory,
            max_tokens=max_tokens, run_once=run_once, clock=now,
        )
        replicates.append(r)
        if on_replicate is not None:
            try:
                on_replicate(i, r)
            except Exception:
                pass

    if n == 1:
        return replicates[0]
    from .benchmark import aggregate_replicates
    return aggregate_replicates(replicates)


# ---------------------------------------------------------------------------
# One run at a time in one checkout
# ---------------------------------------------------------------------------
#
# Demonstrated 2026-08-12, not argued: a full test-suite run failed two
# office fixture tests while a benchmark was running, and `git status`
# caught it mid-attempt -- `D tests/fixtures/office_workspace/buchungen.csv`,
# deleted by the agent under test, not yet restored by the fixture guard.
# Neither side was broken. The MEASUREMENT was, and it reported a product
# defect that did not exist.
#
# The same happened again on 2026-08-13 in the other direction: a suite
# and a live benchmark in one checkout, and the benchmark's own
# "this run changed N files" report blamed itself for the suite's writes.
#
# The lock does not stop anybody. It tells them, which is the whole
# difference between a run that is wrong and a run that says why.
# NOT inside the checkout. A sidecar under tests/fixtures/ is a new path
# appearing during a run, which is exactly what the suite's own leak guard
# exists to report -- and it would have been right to: the checkout must
# look the same after a run as before it. Runtime coordination state lives
# where the rest of this framework's runtime state lives, keyed by the
# checkout it protects so two checkouts never block each other.
_RUN_LOCK_DIR = Path.home() / ".delfin" / "benchmark_locks"


def _run_lock_path(root: Path) -> Path:
    import hashlib
    key = hashlib.sha1(str(Path(root).resolve()).encode("utf-8")).hexdigest()[:16]
    return _RUN_LOCK_DIR / f"{key}.lock"


class BenchmarkRunInProgress(RuntimeError):
    """Another run holds the fixture workspaces."""


@contextlib.contextmanager
def fixture_run_lock(root: Path | str | None = None, *, owner: str = ""):
    """Hold the fixture workspaces for the length of one run.

    An exclusive ``flock`` on a sidecar beside the fixtures, so a second
    benchmark -- or a test suite that takes the same lock -- is told what
    is happening instead of silently measuring a half-written state. The
    holder's pid and what it is doing are written into the file, because
    "resource busy" is not something anybody can act on.

    A filesystem without ``flock`` yields rather than refusing: this is a
    measurement aid, and refusing to measure at all would be worse than
    the interference it prevents.
    """
    base = Path(root) if root is not None else Path(os.getcwd())
    path = _run_lock_path(base)
    handle = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        handle = open(path, "a+")
        import fcntl
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            # "Somebody holds it" and "this filesystem has no flock" are
            # different answers and only one of them is a refusal. NFS
            # without lockd, some container overlays: conflating them
            # turns an unsupported filesystem into a permanent refusal to
            # measure at all, which is worse than the interference.
            import errno as _errno
            if exc.errno not in (_errno.EACCES, _errno.EAGAIN):
                raise
            handle.seek(0)
            held = (handle.read() or "").strip() or "another run"
            handle.close()
            raise BenchmarkRunInProgress(
                f"the fixture workspaces are held by {held}. Two runs in "
                f"one checkout measure each other's writes -- wait for it, "
                f"or use a separate checkout.")
        handle.seek(0)
        handle.truncate()
        handle.write(f"pid {os.getpid()}: {owner or 'a benchmark run'}")
        handle.flush()
    except BenchmarkRunInProgress:
        raise
    except Exception:
        if handle is not None:
            try:
                handle.close()
            except Exception:
                pass
        yield                       # no flock here: measure rather than refuse
        return
    try:
        yield
    finally:
        try:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        except Exception:
            pass
        try:
            handle.close()
        except Exception:
            pass


def run_suite(
    tasks: Iterable[Task],
    *,
    model: str,
    backend: str = "api",
    provider: str = "",
    profile_name: str = "",
    engine_factory: EngineFactory | None = None,
    max_tokens: int = 1024,
    repeats: int = 1,
    progress: Callable[[Task, BenchmarkResult], None] | None = None,
    on_replicate: Callable[[Task, int, BenchmarkResult], None] | None = None,
    run_once: Callable[..., dict] | None = None,
) -> list[BenchmarkResult]:
    """Run every task and return scored (optionally aggregated) results.

    ``repeats`` is forwarded to ``run_task`` so each task is sampled
    N times; the returned list still has one entry per task with the
    aggregate stats.  ``progress`` fires once per task with the final
    aggregated result; ``on_replicate`` fires for each individual
    sample as it lands.
    """

    import os as _os
    _root = Path(_os.getcwd())
    # Held for the whole suite, before the fingerprint is taken: the
    # baseline this run is compared against has to be a state nobody else
    # is writing into. Raises BenchmarkRunInProgress and names the holder
    # if somebody is.
    with fixture_run_lock(_root, owner=f"benchmark run, model {model}"):
        return _run_suite_locked(
            tasks, model=model, backend=backend, provider=provider,
            profile_name=profile_name, engine_factory=engine_factory,
            max_tokens=max_tokens, repeats=repeats, progress=progress,
            on_replicate=on_replicate, run_once=run_once, root=_root)


def _run_suite_locked(
    tasks: Iterable[Task],
    *,
    model: str,
    backend: str = "api",
    provider: str = "",
    profile_name: str = "",
    engine_factory: EngineFactory | None = None,
    max_tokens: int = 1024,
    repeats: int = 1,
    progress: Callable[[Task, BenchmarkResult], None] | None = None,
    on_replicate: Callable[[Task, int, BenchmarkResult], None] | None = None,
    run_once: Callable[..., dict] | None = None,
    root: Path | None = None,
) -> list[BenchmarkResult]:
    """The suite body, under the fixture lock."""
    _root = Path(root) if root is not None else Path(os.getcwd())
    _started_ts = _iso(time.time())
    try:
        _before = checkout_fingerprint(_root)
    except Exception:
        _before = {}

    out: list[BenchmarkResult] = []
    for task in tasks:
        per_task_replicate: Callable[..., None] | None = (
            (lambda t: (lambda i, r: on_replicate(t, i, r)))(task)
            if on_replicate is not None else None
        )
        result = run_task(
            task,
            model=model, backend=backend, provider=provider,
            profile_name=profile_name, engine_factory=engine_factory,
            max_tokens=max_tokens, repeats=repeats,
            run_once=run_once,
            on_replicate=per_task_replicate,
        )
        out.append(result)
        if progress is not None:
            try:
                progress(task, result)
            except Exception:
                pass  # best-effort progress reporting

    # Say what the run changed outside the directories it may change. The
    # fixture guard restores those three; everything else the agent touches
    # persists, and a run that quietly rewrites a test file leaves the next
    # measurement starting from a different tree than the last.
    try:
        _changed = changed_outside_workspaces(_root, _before) if _before else []
    except Exception:
        _changed = []
    if _changed:
        _written = paths_this_run_wrote(_root, since_ts=_started_ts)
        _report = format_checkout_change_report(_changed, _written)
        if _report:
            print(_report)
    return out


# ---------------------------------------------------------------------------
# Profile-name resolution
# ---------------------------------------------------------------------------


def resolve_profile_name(model: str) -> str:
    """Return the profile's ``notes`` field (or 'default') for stamping
    onto results — useful for filtering historical runs by profile."""
    try:
        from .model_profiles import get_profile
        p = get_profile(model)
        return (p.notes or "default")[:80]
    except Exception:
        return "default"


__all__ = [
    "extract_actions",
    "trajectory_from_run",
    "run_task",
    "run_suite",
    "resolve_profile_name",
    "EngineFactory",
]
