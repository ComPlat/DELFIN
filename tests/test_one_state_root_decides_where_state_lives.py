"""One function -- ``state_root()`` -- decides where DELFIN's state lives.

Task s09, first batch: ``~/.delfin`` used to be written directly by about
110 places in ``delfin/`` via ``Path.home() / ".delfin"``, so tests could not
redirect state without patching HOME, and a missed patch wrote the real user
store. This file pins the targets of the first conversion batch.

Part 1 (this commit, green on base): with ``DELFIN_STATE`` unset and HOME
pointed at a fixture root, each of the ten modules resolves a path under the
state root through its own function or constant, and that path is *exactly*
the base one: ``<home>/.delfin/<rel>``.

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


@pytest.fixture
def home_tmp(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("DELFIN_STATE", raising=False)
    os.environ.pop("DELFIN_MESSENGER_DIR", None)
    return home


def test_modules_resolve_under_home_when_state_unset(home_tmp):
    """With HOME pointed at a fixture root and DELFIN_STATE unset, each of
    the ten modules' paths is exactly the base one: <home>/.delfin/<rel>.

    Reloading first ensures the value read is the module's own code output,
    freshly computed under the exhausted environment, not an autouse
    fixture's monkeypatch of the same attribute."""
    for module_name, how, rel in PROBES:
        importlib.reload(importlib.import_module(module_name))
        got = _read_probe(module_name, how)
        expected = Path(os.environ["HOME"]) / ".delfin" / rel
        assert got == expected, (
            f"{module_name} {how[1]} resolved to {got}; expected {expected}")
