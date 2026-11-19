"""One function -- ``state_root()`` -- decides where DELFIN's state lives.

Task s09, first batch: ``~/.delfin`` used to be written directly by about
110 places in ``delfin/`` via ``Path.home() / ".delfin"``, so tests could not
redirect state without patching HOME, and a missed patch wrote the real user
store. This file pins the targets of the first conversion batch.

Part 1 (green on base): with ``DELFIN_STATE`` unset and HOME pointed at a
fixture root, each of the ten modules resolves a path under the state root
through its own function or constant, and that path is *exactly* the base
one: ``<home>/.delfin/<rel>``. This pins today's behaviour.

Part 2 (RED on base, green once converted): with ``DELFIN_STATE`` set to
another fixture root, the same modules resolve through it. On base they
write ``Path.home()`` and ignore the variable, so this fails there.

Part 3 (RED on base, green once converted): none of the ten modules
contains the literal ``Path.home() / ".delfin"`` any more.

Each probe reloads the module first so the value read is the module's own,
freshly computed under the current environment -- an autouse fixture that
monkeypatched the same attribute is overwritten by the reload and cannot
mask what the code actually produces.
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path

import pytest

# (module, how to read the path, relative path under the state root).
# how is ("attr", name) for a module constant, ("call", name) for a public
# path function read after reload.
MIN_MEMORY = "delfin.agent.memory_store"
MIN_SESSION = "delfin.agent.session_store"
MIN_SUBAGENTS = "delfin.agent.subagents"
MIN_SKILLS = "delfin.agent.skills"
MIN_LOOPS = "delfin.agent.loops"
MIN_EVAL = "delfin.agent.eval_loop"
MIN_BUG = "delfin.agent.bug_report"
MIN_STORE = "delfin.messenger.store"
MIN_PAIRING = "delfin.messenger.pairing"
MIN_SCHEDULES = "delfin.dashboard.schedules"

PROBES = (
    (MIN_MEMORY, ("attr", "_DEFAULT_PATH"), "agent_memory.json"),
    (MIN_MEMORY, ("attr", "_MEMORY_WRITE_LOCK"), "memory_writes"),
    (MIN_MEMORY, ("call", "_delfin_global_memory_dir"), "memory"),
    (MIN_MEMORY, ("call", "_office_ws_file"), "office_workspaces.json"),
    (MIN_SESSION, ("attr", "_SESSIONS_DIR"), "agent_sessions"),
    (MIN_SESSION, ("call", "_transcript_archive_path"), "transcript_archive"),
    (MIN_SESSION, ("call", "_handoffs_path"), "handoffs"),
    (MIN_SESSION, ("call", "_bundles_path"), "bundles"),
    (MIN_SUBAGENTS, ("attr", "_TELEMETRY_PATH"), "subagent_telemetry.jsonl"),
    (MIN_SUBAGENTS, ("attr", "_RUNNING_DIR"), "subagent_running"),
    (MIN_SUBAGENTS, ("attr", "_PENDING_DIR"), "subagent_pending"),
    (MIN_SUBAGENTS, ("attr", "_SESSIONS_DIR"), "subagent_sessions"),
    (MIN_SUBAGENTS, ("call", "_md_preset_search_dirs"), "subagents"),
    (MIN_SKILLS, ("call", "_skill_dirs"), "skills"),
    (MIN_LOOPS, ("call", "personal_dir"), "loops"),
    (MIN_LOOPS, ("call", "ledger_path"), "loops_confirmed.json"),
    (MIN_EVAL, ("attr", "_REPORTS_DIR"), "eval_reports"),
    (MIN_EVAL, ("attr", "_TASK_DRAFTS_DIR"), "bug_tasks"),
    (MIN_BUG, ("attr", "_FALLBACK_DIR"), "agent_bugs"),
    (MIN_BUG, ("attr", "_TASK_DRAFTS_DIR"), "bug_tasks"),
    (MIN_STORE, ("call", "state_dir"), "messenger"),
    (MIN_PAIRING, ("attr", "_DIR"), "messenger_pairing"),
    (MIN_SCHEDULES, ("attr", "_DEFAULT_PATH"), "schedules.json"),
)

# The ten modules under conversion, as source paths to grep for the missing
# ``Path.home() / ".delfin"`` literal.
TEN_MODULES = (
    "delfin/agent/memory_store.py",
    "delfin/agent/session_store.py",
    "delfin/agent/subagents.py",
    "delfin/agent/skills.py",
    "delfin/agent/loops.py",
    "delfin/agent/eval_loop.py",
    "delfin/agent/bug_report.py",
    "delfin/messenger/store.py",
    "delfin/messenger/pairing.py",
    "delfin/dashboard/schedules.py",
)


def _read_probe(module_name: str, how: tuple[str, str]) -> Path:
    mod = importlib.import_module(module_name)
    if how[0] == "attr":
        return getattr(mod, how[1])
    # "call" -- some path functions need an argument or return a list.
    fn = getattr(mod, how[1])
    if how[1] == "_md_preset_search_dirs":
        value = fn()[1]
    elif how[1] == "_skill_dirs":
        value = fn(None)[1]
    else:
        value = fn()
    return value


def _probe_fresh(module_name: str, how: tuple[str, str]) -> Path:
    """Re-run the module's code under the current environment, read the probe,
    then put the module's original namespace back.

    ``importlib.reload`` re-executes the module in the same module object, so
    without the restore every class and constant other modules already hold
    (``SubagentResult``, a path constant) is replaced for the rest of the test
    process, and a later ``isinstance`` check or write sees the probe's copy."""
    mod = importlib.import_module(module_name)
    saved = dict(vars(mod))
    try:
        importlib.reload(mod)
        return _read_probe(module_name, how)
    finally:
        vars(mod).clear()
        vars(mod).update(saved)


@pytest.fixture
def home_tmp(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("DELFIN_STATE", raising=False)
    monkeypatch.delenv("DELFIN_MESSENGER_DIR", raising=False)
    return home


def test_modules_resolve_under_home_when_state_unset(home_tmp):
    """With HOME pointed at a fixture root and DELFIN_STATE unset, each of
    the ten modules' paths is exactly the base one: <home>/.delfin/<rel>.

    Reloading first ensures the value read is the module's own code output,
    freshly computed under the exhausted environment, not an autouse
    fixture's monkeypatch of the same attribute."""
    for module_name, how, rel in PROBES:
        got = _probe_fresh(module_name, how)
        expected = Path(os.environ["HOME"]) / ".delfin" / rel
        assert got == expected, (
            f"{module_name} {how[1]} resolved to {got}; expected {expected}")


def test_modules_resolve_under_delfin_state_when_set(tmp_path, monkeypatch):
    """With DELFIN_STATE set to a fixture root, the converted modules use it.
    Red on base (they write ``Path.home()``); green once converted."""
    state = tmp_path / "state"
    monkeypatch.delenv("DELFIN_STATE", raising=False)
    monkeypatch.setenv("DELFIN_STATE", str(state))
    monkeypatch.setenv("HOME", str(tmp_path / "home_ignored"))
    monkeypatch.delenv("DELFIN_MESSENGER_DIR", raising=False)
    for module_name, how, rel in PROBES:
        got = _probe_fresh(module_name, how)
        expected = state / rel
        assert got == expected, (
            f"{module_name} {how[1]} resolved to {got}; expected {expected} ("
            "DELFIN_STATE is not honoured -- not yet converted to state_root())")


def test_no_module_writes_path_home_delfin_literal():
    """The conversion is source-visible: after it, none of the ten modules
    contains the literal ``Path.home() / ".delfin"`` any more."""
    import delfin.agent.state_paths as sp
    repo_root = Path(sp.__file__).resolve().parent.parent.parent
    missing = []
    for rel in TEN_MODULES:
        src = (repo_root / rel).read_text(encoding="utf-8")
        if 'Path.home() / ".delfin"' in src:
            missing.append(rel)
    assert not missing, (
        "these module(s) still write Path.home() / \".delfin\" directly; "
        "convert them to state_root(): " + ", ".join(missing))
