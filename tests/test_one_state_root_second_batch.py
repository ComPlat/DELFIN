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


def _read_probe(module_name: str, how: tuple[str, str]) -> Path:
    mod = importlib.import_module(module_name)
    if how[0] == "attr":
        return getattr(mod, how[1])
    # "call" -- every path function in this batch takes no argument.
    fn = getattr(mod, how[1])
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
