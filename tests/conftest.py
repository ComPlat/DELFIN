"""Shared test fixtures.

Redirect the subagent state files (live registry, telemetry, finished
sessions, outstanding background reports) away from the real
``~/.delfin`` so test runs never leave artifacts that show up in the
user's dashboard (live panel, ``/agents`` listing, stats) or get pushed
into a real session's next turn.
"""

from __future__ import annotations

import pytest


# ---------------------------------------------------------------------------
# The checkout itself is a sink, and it was invisible
# ---------------------------------------------------------------------------
# The user-home guard above covers ``~/.delfin``. It does not cover the
# other half of the same incident: a mock standing in for a workspace
# yields the RELATIVE path ``MagicMock/<attribute>/<id>``, so every store
# built from one wrote under the process CWD -- the repository checkout.
# An ordinary run left 3.1 MB in 560 files there, and
# ``git status --untracked-files=all`` reported nothing, because the
# ``.delfin/`` ignore rule structurally covers exactly those paths.
#
# ``TaskStore`` now refuses a non-absolute base directory, which stops
# that particular class at the source. This fixture is the other half:
# it does not know what the next leak will be, only that the checkout
# must look the same after a run as before it.
_CHECKOUT_ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]

# Directories the PRODUCT or the suite populates on purpose. Each is a
# declared output location, not test spill, and each is named here with
# the reason -- an entry may not be added on the grounds that git ignores
# it, because catching a gitignored write is exactly what this guard is
# for. The incident it was built for hid under the `.delfin/` rule.
#
# Enumerated from two CI failures rather than guessed a third time: this
# guard cannot fail on a machine where the toolchain is already installed
# and the fixtures already built, and a fresh clone -- which is what CI
# is -- creates every one of them. Roots, not leaf names, so the next
# file the installer adds does not need a third round.
_GENERATED_ROOTS = (
    # The chemistry toolchain, installed on demand by runtime_setup.
    "delfin/qm_tools",
    # The benchmark workbooks, materialised from a reviewable spec before
    # every run -- see delfin/agent/benchmark_fixtures.py.
    "tests/fixtures/office_workspace",
)

# Build noise: byte-code caches and tool caches created by running at all,
# plus the worktrees parallel agents check out under `.claude`, which are
# separate working copies with their own runs going on in them.
_CHECKOUT_NOISE = ("__pycache__", ".pytest_cache", ".git", ".mypy_cache",
                   ".ruff_cache", ".hypothesis",
                   ".claude", ".venv", "node_modules")


def _checkout_entries() -> frozenset:
    """Every path under the checkout, minus build noise and declared
    output locations."""
    import os
    roots = tuple(str(_CHECKOUT_ROOT / r) for r in _GENERATED_ROOTS)

    def _is_generated(path: str) -> bool:
        # Containment, not a string prefix. Comparing the raw strings would
        # exempt a sibling whose name merely begins with a declared root --
        # `qm_tools_backup` under `qm_tools` -- and an exemption this guard
        # never announces is the one place it must not be generous.
        return any(path == r or path.startswith(r + os.sep) for r in roots)

    out = set()
    for dirpath, dirnames, filenames in os.walk(_CHECKOUT_ROOT):
        dirnames[:] = [d for d in dirnames if d not in _CHECKOUT_NOISE]
        if _is_generated(dirpath):
            dirnames[:] = []
            continue
        for name in dirnames + filenames:
            if name.endswith(".pyc") or name in _CHECKOUT_NOISE:
                continue
            full = os.path.join(dirpath, name)
            if _is_generated(full):
                continue
            out.add(full)
    return frozenset(out)


def _unexpected_under_a_generated_root() -> frozenset:
    """What appeared inside a declared output location that its own ignore
    rules do not describe.

    Exempting the two roots wholesale would have cost real coverage: both
    hold tracked source -- the installer scripts and the CSV fixtures --
    and one of them is where the office benchmark writes. A leak there
    would have been the one place this guard stayed quiet.

    Git already carries the distinction, written by hand and reviewed:
    ``.gitignore`` names the seven directories the installer fills and the
    workbook patterns the generator writes, under the heading that says to
    keep the scripts and ignore the local install. So the rule inside a
    declared root is git's answer, not a second list here that would drift
    from the first. ``docs`` is created empty and never filled, and an
    empty directory does not exist as far as git is concerned.

    No git, or a checkout that is not a repository: the roots stay exempt,
    which is where this stood before. Reporting nothing beats failing a
    run over a missing tool.
    """
    import subprocess
    try:
        done = subprocess.run(
            ["git", "-C", str(_CHECKOUT_ROOT), "status", "--porcelain", "-z",
             "--", *_GENERATED_ROOTS],
            capture_output=True, text=True, timeout=60, check=False)
    except (OSError, subprocess.SubprocessError):
        return frozenset()
    if done.returncode != 0:
        return frozenset()
    # `XY <path>` per entry, NUL-separated; a rename adds its source as a
    # second entry, which is a path either way.
    return frozenset(
        entry[3:] for entry in done.stdout.split("\0") if entry[3:])


@pytest.fixture(scope="session", autouse=True)
def _the_suite_does_not_write_into_the_checkout():
    """Fail the run when new paths appeared in the checkout.

    Known limit, stated rather than implied: this compares the SET of
    paths, not their contents. A leak that already ran once in this
    checkout writes the same paths again and is not re-reported until the
    tree is cleaned. On a fresh clone -- which is what CI is -- every leak
    is a new path, so nothing escapes there. Comparing contents instead
    would mean hashing the tree twice per run, and would false-fail on
    every file the package legitimately rewrites.
    """
    before = _checkout_entries()
    before_generated = _unexpected_under_a_generated_root()
    yield
    new = sorted(p for p in _checkout_entries() - before)
    new += sorted(
        str(_CHECKOUT_ROOT / p)
        for p in _unexpected_under_a_generated_root() - before_generated)
    assert not new, (
        f"{len(new)} path(s) appeared in the checkout during the run; "
        "a git ignore rule may make them invisible to `git status`: "
        + ", ".join(str(p) for p in new[:20])
    )


_CWD_BEFORE_TEST: dict[str, str] = {}


@pytest.hookimpl(tryfirst=True)
def pytest_runtest_setup(item):
    import os
    _CWD_BEFORE_TEST[item.nodeid] = os.getcwd()


@pytest.hookimpl(trylast=True)
def pytest_runtest_teardown(item, nextitem):
    """A test may not leave the process somewhere else.

    The process CWD is shared by every test that follows, so one chdir
    that is not undone turns every later relative path into a different
    file -- and the failure then lands on whoever wrote the relative path
    rather than on whoever moved the ground under it. Observed on both
    sides of that: a file of mine that passes alone and fails eleven
    times inside the full run, and a test named for launch-directory
    independence failing only in the big block.

    A hook and not an autouse fixture, which is where the first attempt
    was wrong: fixture teardown runs BEFORE ``monkeypatch`` undoes its
    own chdir, so the guard accused a test that had done exactly the
    right thing. ``trylast`` puts this after every finalizer, so what it
    sees is what the next test will get. A guard that blames correct work
    is one somebody switches off.

    Restored as well as reported: a run that stops dead at the first
    polluter tells you less than one that names it and carries on.
    """
    import os
    before = _CWD_BEFORE_TEST.pop(item.nodeid, None)
    if before is None:
        return
    after = os.getcwd()
    if after != before:
        os.chdir(before)
        raise AssertionError(
            f"this test left the process in {after} instead of {before}; "
            "a chdir from a test must be undone (monkeypatch.chdir, or a "
            "try/finally), or every later relative path reads a "
            "different file")


@pytest.fixture(autouse=True)
def _reset_workspace_trust_caches():
    """Trust state is process-global: a parsed store and a record of which
    refusals have already been reported. Both would otherwise leak from one
    test's tmp workspace into the next one's assertions."""
    from delfin.agent import workspace_trust as wt
    wt.reset_announcements()
    yield
    wt.reset_announcements()


@pytest.fixture(autouse=True)
def _isolate_subagent_state(tmp_path, monkeypatch):
    from delfin.agent import subagents as sa
    monkeypatch.setattr(sa, "_RUNNING_DIR",
                        tmp_path / "subagent_running")
    monkeypatch.setattr(sa, "_TELEMETRY_PATH",
                        tmp_path / "subagent_telemetry.jsonl")
    monkeypatch.setattr(sa, "_SESSIONS_DIR",
                        tmp_path / "subagent_sessions")
    # Reserving a background id records a report the parent is owed; the
    # engine drains those into the next turn, so a test that spawned one
    # would otherwise announce itself in a real session.
    monkeypatch.setattr(sa, "_PENDING_DIR",
                        tmp_path / "subagent_pending")
    yield


@pytest.fixture(autouse=True)
def _isolate_user_wide_memory(tmp_path, monkeypatch):
    """Keep the user-scoped stores out of the real ``~/.delfin``.

    ``tmp_path`` bounds the PROJECT store, because that one is derived from
    the repo root a test passes in. The two USER-scoped stores are not: they
    resolve through ``Path.home()`` and ignore whatever root the test used.

    Both were observed writing into the real home directory during an
    ordinary run. The office registry had collected 178 entries, 176 of them
    reprs of mock objects a test had passed where a path was expected — and
    that registry is one of the carriers of the folder lock, so its contents
    decide which directory counts as an office workspace. The repository
    checkout itself was among them, which made a full-suite run treat DELFIN
    as an office folder and fail one memory test that passed in isolation.

    The redirect steps aside for a test that has pointed ``Path.home`` at a
    directory of its own: those tests were already isolated, by the more
    direct route, and silently overriding them would move the store out from
    under their assertions.
    """
    import pathlib

    from delfin.agent import memory_store as ms
    real_home = pathlib.Path.home()
    home = tmp_path / "user_home_delfin"

    def _redirected(name: str) -> pathlib.Path:
        current = pathlib.Path.home()
        if current != real_home:
            return current / ".delfin" / name
        return home / name

    monkeypatch.setattr(
        ms, "_delfin_global_memory_dir", lambda: _redirected("memory"))
    monkeypatch.setattr(
        ms, "_office_ws_file",
        lambda: _redirected("office_workspaces.json"))
    # The registry caches its parsed contents in a module global, so a stale
    # cache from an earlier test would outlive the redirect.
    monkeypatch.setattr(ms, "_office_ws_cache", None)
    yield
    monkeypatch.setattr(ms, "_office_ws_cache", None)


# Every module-level sink under the user's ~/.delfin that a test can write to.
# These are resolved at IMPORT time, so a test pointing Path.home elsewhere
# never redirected them -- which is why an ordinary run left 203 job records,
# 141 of them naming pytest tmp directories, in the user's real locator index.
#
# Read-mostly stores are deliberately absent: the documentation index, the
# model-capability cache, credentials (isolated separately) and the settings
# files. Redirecting those would not stop a write, it would only hide the
# real content from tests that legitimately read it.
_USER_STATE_SINKS: tuple[tuple[str, str, str], ...] = (
    ("delfin.agent.bash_jobs", "_INDEX_PATH", "bash_jobs_index.json"),
    ("delfin.agent.provider_profile", "_LOCAL_STATE_PATH",
     "provider_profile_state.json"),
    ("delfin.agent.job_fix", "_ATTEMPTS_PATH", "fix_attempts.json"),
    ("delfin.agent.session_store", "_SESSIONS_DIR", "agent_sessions"),
    ("delfin.agent.outcome_tracker", "_DEFAULT_PATH", "outcome_history.jsonl"),
    ("delfin.agent.agent_metrics", "_LOG_PATH", "agent_metrics.jsonl"),
    ("delfin.agent.context_tracker", "_DEFAULT_PATH", "context_usage.jsonl"),
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
    # real file -- and permission rules are exactly what must not be granted
    # by accident.
    ("delfin.agent.hooks_editor", "_USER_SETTINGS", "settings.json"),
    ("delfin.agent.kit_settings", "USER_SETTINGS_PATH", "settings.json"),
)


@pytest.fixture(scope="session")
def _user_state_targets():
    """Resolve the sink table once. A module that cannot be imported (an
    optional dependency is missing) simply drops out."""
    import importlib
    out = []
    for mod_name, attr, rel in _USER_STATE_SINKS:
        try:
            mod = importlib.import_module(mod_name)
        except Exception:
            continue
        if hasattr(mod, attr):
            out.append((mod, attr, rel))
    return tuple(out)


# The sinks that are resolved per call rather than at import. Redirecting
# these one by one, rather than swapping Path.home for the whole suite: that
# was measured, and it breaks 357 tests that legitimately read the real home.
_USER_STATE_RESOLVERS: tuple[tuple[str, str, str], ...] = (
    ("delfin.agent.audit_log", "_default_log_path", "audit.log"),
    # Which directories the user has trusted to run commands. A test that
    # granted trust must never grant it in the real store: the entry would
    # outlive the run and let a later, real session honour a workspace's
    # hooks and MCP servers on the strength of a pytest tmp directory.
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
)


@pytest.fixture(scope="session")
def _user_state_resolvers():
    import importlib
    out = []
    for mod_name, attr, rel in _USER_STATE_RESOLVERS:
        try:
            mod = importlib.import_module(mod_name)
        except Exception:
            continue
        if callable(getattr(mod, attr, None)):
            out.append((mod, attr, rel))
    return tuple(out)


@pytest.fixture(autouse=True)
def _isolate_user_state(tmp_path, monkeypatch, _user_state_targets,
                        _user_state_resolvers):
    """Point every writable user-state sink into the test's own directory.

    The per-project memory store belongs here too, and is the one that was
    least obvious: its path is ``~/.delfin/projects/<slug>/memory``, where
    the slug is derived from the repo root. Passing ``tmp_path`` as the repo
    root therefore did NOT bound it -- it only changed the slug, and every
    test that saved a memory left a directory named after itself in the
    user's real home.

    A test that patches one of these itself still wins: its ``setattr`` runs
    after this fixture's."""
    import pathlib

    real_home = pathlib.Path.home()
    fallback = tmp_path / "user_home"

    def _already_isolated() -> bool:
        """A test that has pointed ``Path.home`` somewhere of its own.

        Those were isolated by the more direct route, many of them assert on
        the resulting path, and some rely on side effects of the original
        resolver (the legacy-store migration). For them the fixture steps
        aside entirely and calls the original.
        """
        return pathlib.Path.home() != real_home

    for mod, attr, rel in _user_state_targets:
        monkeypatch.setattr(mod, attr, fallback / ".delfin" / rel)

    _PROJECT_LEAF = {"_delfin_memory_dir": "memory",
                     "_delfin_plans_dir": "plans"}
    for mod, attr, rel in _user_state_resolvers:
        original = getattr(mod, attr)
        leaf = _PROJECT_LEAF.get(attr)
        if leaf is not None:
            def _resolve(repo_root, _o=original, _leaf=leaf):
                if _already_isolated():
                    return _o(repo_root)
                from delfin.agent.memory_store import _project_slug
                return (fallback / ".delfin" / "projects"
                        / _project_slug(repo_root) / _leaf)
        else:
            def _resolve(_o=original, _rel=rel):
                if _already_isolated():
                    return _o()
                return fallback / ".delfin" / _rel
        monkeypatch.setattr(mod, attr, _resolve)


# ---------------------------------------------------------------------------
# The sandbox-escape gate
# ---------------------------------------------------------------------------
# The only tests that put a real command inside a real sandbox and check
# that it cannot get out were ``skipif``-gated on a functional bwrap. CI
# runs on a plain ubuntu runner with no bubblewrap installed and no apt
# step to add it, so the gate evaporates exactly where nobody is
# watching: the tests report as skipped, the run is green, and no
# assertion about confinement has been made in months. Everything else in
# those files inspects an argv LIST -- it asserts what the wrap would be
# asked to do, never that anything was confined.
#
# Skipping is still right on a developer laptop without bwrap. It is not
# right on a machine that is supposed to have it. ``DELFIN_EXPECT_ISOLATION=1``
# is the difference: with it set, an unusable sandbox FAILS instead of
# skipping, so an environment that loses bubblewrap says so.
_ISOLATION_EXPECTED_ENV = "DELFIN_EXPECT_ISOLATION"


def isolation_is_expected() -> bool:
    """Whether this environment promises a working sandbox."""
    import os
    return str(os.environ.get(_ISOLATION_EXPECTED_ENV, "")).strip().lower() in (
        "1", "true", "yes")


def sandbox_is_functional() -> bool:
    """Whether bwrap can actually confine anything here.

    Production's own probe, not a second one: a private copy in the test
    file probed the MINIMAL wrap (``bwrap --ro-bind / / true``) while the
    real wrap also asks for ``--dev /dev --proc /proc --tmpfs /tmp``. On a
    host where the minimal probe passes and the full one does not, the
    gate opened for a wrap that cannot be built.
    """
    try:
        from delfin.agent.api_client import _bwrap_functional
        return bool(_bwrap_functional())
    except Exception:
        return False


def requires_a_working_sandbox(fn):
    """Skip without a usable sandbox -- unless one was promised."""
    import functools

    if sandbox_is_functional():
        return fn
    if isolation_is_expected():
        @functools.wraps(fn)
        def _fail(*args, **kwargs):
            pytest.fail(
                f"{_ISOLATION_EXPECTED_ENV} is set, so this environment "
                "promises a working sandbox, but bwrap could not build the "
                "wrap the agent uses. Confinement is NOT being tested here."
            )
        return _fail
    return pytest.mark.skip(
        reason=("no functional bwrap (unprivileged container?); set "
                f"{_ISOLATION_EXPECTED_ENV}=1 to make this a failure")
    )(fn)


# Secret API keys delfin resolves, each from the env OR ~/.delfin/credentials.json.
_SECRET_KEY_VARS = ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY")


@pytest.fixture(autouse=True)
def _ci_parity_no_secrets(tmp_path, monkeypatch):
    """Make the local test environment match CI's: NO API key reachable from
    either the environment or the ``~/.delfin/credentials.json`` store.

    CI runs key-free and fully mocked (see ``.github/workflows/ci.yml`` —
    "no secrets are required"). A test that quietly relies on an ambient key
    therefore passes on a developer box that has one, but fails in CI — the
    exact "lokal grün, CI rot" trap. Stripping the keys here surfaces that
    whole class of failure on the developer's own machine, before the push.

    A test that genuinely needs a key still sets it itself (this autouse
    fixture runs first; the test's own ``monkeypatch.setenv`` then wins).
    """
    for var in _SECRET_KEY_VARS:
        monkeypatch.delenv(var, raising=False)
    try:
        from delfin.agent import credentials as _cred
        monkeypatch.setattr(_cred, "_DEFAULT_PATH",
                            tmp_path / "no_credentials.json")
    except Exception:
        pass
    yield


# ---------------------------------------------------------------------------
# CI tiering (2026-06-27). The fast gate (`pytest -m "not slow"`) is the merge
# check and must stay ~minutes. Two kinds of test are auto-marked `slow` and
# moved to the non-blocking slow-tests job:
#   (a) heavy MANTA build/integration tests (construct full complexes; the
#       slowest single test is ~130s — together they push the suite >60 min),
#   (b) ORDER-DEPENDENT tests that only pass in the full serial context (they
#       rely on module-level state another test sets up — sys.modules-driven
#       scans, runtime-mutated registries); in any subset they false-fail.
# Real fix for (b) = test isolation (tracked debt); until then they run in the
# (non-blocking) slow job so they neither block nor false-fail the gate.
# Add new build-heavy / order-dependent test files to this set.
# ---------------------------------------------------------------------------
_SLOW_TEST_FILES = frozenset({
    # (a) build/integration
    "test_all_features_contracts.py", "test_classify_and_router.py",
    "test_crest_censo.py", "test_decompose_kekulize_split.py",
    "test_embed_determinism.py", "test_equatorial_square_kappa4.py",
    "test_eta_ring_hydrogens.py", "test_fffree_backbone_reembed.py",
    "test_fffree_chelate_backbone.py", "test_fffree_conformer_complete.py",
    "test_fffree_conformer_coverage.py", "test_fffree_conformer_seating.py",
    "test_fffree_diatomic_orient.py", "test_fffree_hapto_recall.py",
    "test_fffree_interlig_vdw_gate.py", "test_fffree_joint_declash.py",
    "test_fffree_metalloid_donor.py", "test_fffree_nhc_carbene.py",
    "test_highcn_coverage.py", "test_isomer_benchmark.py",
    "test_isomer_enum_determinism.py", "test_ligand_rigid_seating.py",
    "test_metal_isomer_determinism.py", "test_multihapto_patches.py",
    "test_multi_sigma_path_v2.py", "test_platform.py",
    "test_ring_bounds_embed.py", "test_smiles_converter_regressions.py",
    "test_tool_contracts.py", "test_uff_soft_determinism.py",
    "test_user_smiles_suite.py",
    "test_welle5l_t3b_daqiwaz_pyridine_coverage.py",
    # (b) order-dependent (pass only in the full serial context)
    "test_no_regression_undefined_names.py", "test_post_optimizer.py",
    "test_theorem_d_asymmetric.py", "test_welle5p_c_hotfix.py",
})


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: heavy or order-dependent test, excluded from the fast CI gate",
    )


def pytest_collection_modifyitems(config, items):
    for item in items:
        fname = item.nodeid.split("::", 1)[0].rsplit("/", 1)[-1]
        if fname in _SLOW_TEST_FILES:
            item.add_marker(pytest.mark.slow)


@pytest.fixture(autouse=True)
def _mute_out_of_band_notifications(monkeypatch):
    """No test may reach the user's desktop or the network.

    The attention inbox itself is redirected to tmp above, but its
    delivery fan-out is not a file: it shells out to ``notify-send`` and
    POSTs to the configured webhook. Any test that emits an attention
    event — a confirm timing out, a scheduled run failing, and now a
    containment block — would pop a real notification on the developer's
    screen. Tests that assert on the transports patch these themselves,
    after this fixture.
    """
    from delfin.agent import notify as _n
    monkeypatch.setattr(_n, "send_notification", lambda *a, **k: False)
    monkeypatch.setattr(
        _n, "send_remote_trigger",
        lambda *a, **k: _n.TriggerResult(sent=False, error="muted in tests"))


@pytest.fixture(autouse=True)
def _isolate_agent_telemetry(monkeypatch, tmp_path):
    """Redirect the agent's per-user telemetry sinks into the test tmp dir.

    Engine fixtures with mocked clients still execute the real
    turn-metrics/tool-trace/failure-log writers at turn end; without this
    every test run leaves mock records in the user's real ~/.delfin and
    pollutes the eval loop's health statistics. Tests that monkeypatch
    these constants themselves simply override this baseline."""
    from delfin.agent import failure_log as _fl
    from delfin.agent import tool_trace as _tt
    from delfin.agent import turn_metrics as _tm
    monkeypatch.setattr(_tm, "_DIR", tmp_path / "turn_metrics")
    monkeypatch.setattr(_tt, "_DIR", tmp_path / "tool_trace")
    monkeypatch.setattr(_fl, "_LOG_PATH", tmp_path / "failure_log.jsonl")


@pytest.fixture(autouse=True)
def _state_maintenance_runs_once_per_test(monkeypatch):
    """Keep the one-shot state sweep from leaking across tests.

    ``run_startup_maintenance`` latches a module global so it costs nothing
    after the first session write. Left alone in a suite, whichever test
    ran first would consume the latch and every later test would silently
    skip the sweep — including the tests that assert it happens.
    """
    from delfin.agent import state_paths as _sp
    monkeypatch.setattr(_sp, "_MAINTENANCE_DONE", False)
    yield
