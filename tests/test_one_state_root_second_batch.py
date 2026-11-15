"""Second batch of step s09: ten more modules resolve their root through ``state_root()``.

Batch 1 (PR #735) made ``delfin.agent.state_paths.state_root()`` the single
function that decides where DELFIN's state lives (``DELFIN_STATE`` env var
if set, else ``~/.delfin``) and converted ten modules. These ten still spell
``Path.home() / ".delfin"`` themselves; this file pins their behaviour the
same way batch 1 pinned the first ten, so the conversion here is
characterized before it is made.

Part 1 (green on base): with ``DELFIN_STATE`` unset and HOME pointed at a
fixture root, each of the ten modules resolves its state path through its
own function or constant, and that path is *exactly* the base one:
``<home>/.delfin/<rel>``.

Part 2 (RED on base, green once converted): with ``DELFIN_STATE`` set to
another fixture root, the same modules resolve through it. On base they
write ``Path.home()`` and ignore the variable, so this fails there.

Part 3 (RED on base, green once converted): none of these ten modules
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
MIN_INDEX = "delfin.agent.session_index"
MIN_MESSAGES = "delfin.agent.session_messages"
MIN_PAUSE = "delfin.agent.session_pause"
MIN_PRESENCE = "delfin.agent.session_presence"
MIN_REPORT = "delfin.agent.session_report"
MIN_COMPACTION = "delfin.agent.compaction_log"
MIN_FAILURE = "delfin.agent.failure_log"
MIN_OUTCOME = "delfin.agent.outcome_tracker"
MIN_CHANGE = "delfin.agent.change_journal"
MIN_CONTEXT = "delfin.agent.context_tracker"

PROBES = (
    (MIN_INDEX, ("call", "_index_path"), "session_index.sqlite"),
    (MIN_MESSAGES, ("attr", "_DIR"), "session_inbox"),
    (MIN_PAUSE, ("attr", "_DIR"), "session_pause"),
    (MIN_PRESENCE, ("attr", "_DIR"), "session_presence"),
    (MIN_REPORT, ("call", "_report_dir"), "session_reports"),
    (MIN_COMPACTION, ("attr", "_LOG_DIR"), "compaction_log"),
    (MIN_FAILURE, ("attr", "_LOG_PATH"), "failure_log.jsonl"),
    (MIN_OUTCOME, ("attr", "_DEFAULT_PATH"), "outcome_history.jsonl"),
    (MIN_CHANGE, ("call", "_undo_root"), "undo"),
    (MIN_CONTEXT, ("attr", "_DEFAULT_PATH"), "context_usage.jsonl"),
)

# The ten modules under conversion, as source paths to grep for the missing
# ``Path.home() / ".delfin"`` literal.
TEN_MODULES = (
    "delfin/agent/session_index.py",
    "delfin/agent/session_messages.py",
    "delfin/agent/session_pause.py",
    "delfin/agent/session_presence.py",
    "delfin/agent/session_report.py",
    "delfin/agent/compaction_log.py",
    "delfin/agent/failure_log.py",
    "delfin/agent/outcome_tracker.py",
    "delfin/agent/change_journal.py",
    "delfin/agent/context_tracker.py",
)

# Every literal spelling of the HOME-anchored state root we must leave
# behind. Spacing variants are listed because the conversion is a source
# change and a reformat would otherwise slip one through.
_HOME_STATE_LITERALS = (
    'Path.home() / ".delfin"',
    'Path.home()/".delfin"',
    'expanduser("~/.delfin',
)


def _read_probe(module_name: str, how: tuple[str, str]) -> Path:
    mod = importlib.import_module(module_name)
    if how[0] == "attr":
        return getattr(mod, how[1])
    # "call" -- every path function in this batch takes no argument.
    fn = getattr(mod, how[1])
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
        if any(literal in src for literal in _HOME_STATE_LITERALS):
            missing.append(rel)
    assert not missing, (
        "these module(s) still write the HOME state path directly; "
        "convert them to state_root(): " + ", ".join(missing))
