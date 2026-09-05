"""One store, one session id, two surfaces that used to disagree.

THE INCIDENT. ``TaskStore.list`` treated ``""`` and ``None`` as
opposites: ``None`` meant "every session", ``""`` meant "nothing at
all". The model-facing reminder passed ``None`` when the id was empty
and therefore saw EVERY session's leftovers; the user's ticker passed
``""`` and therefore saw none of them and printed "No tasks yet". The
CLI backend mints no session id, so this was not a corner case — it was
the default state of that backend: the panel said the list was empty
while the prompt was listing it.

Both now go through ``resolve_session_scope``, and empty means unscoped
— the reading ``task_list``, ``task_adopt`` and the reminder already
documented.
"""

from __future__ import annotations

import re

import pytest

from delfin.agent.agent_tasks import TaskStore, get_store, resolve_session_scope
from delfin.agent import task_ticker as TT
from delfin.agent.engine import AgentEngine


@pytest.fixture
def store(tmp_path):
    from delfin.agent import agent_tasks
    agent_tasks._STORES.clear()
    return get_store(tmp_path)


def _block_for(monkeypatch, workspace, session_id: str) -> str:
    class _Perms:
        pass
    p = _Perms()
    p.workspace = workspace
    p.task_session_id = session_id
    p.mode = "default"
    eng = AgentEngine.__new__(AgentEngine)
    monkeypatch.setattr(AgentEngine, "kit_permissions",
                        property(lambda self: p))
    return eng._build_open_tasks_block()


def _panel_rows(html: str) -> int:
    return len(re.findall(r"font-family:monospace", html))


def _block_rows(block: str) -> int:
    return len([ln for ln in block.splitlines()
                if ln.startswith("- [")])


@pytest.mark.parametrize("session_id", ["", "sess-a"])
def test_the_prompt_and_the_panel_count_the_same_tasks(
        store, monkeypatch, session_id):
    store.create("first", session_id="sess-a")
    store.create("second", session_id="sess-b")
    block = _block_for(monkeypatch, store.base_dir, session_id)
    html = TT.render_html(store.base_dir, session_id=session_id)
    assert _block_rows(block) == _panel_rows(html) > 0


def test_an_empty_id_is_unscoped_on_both_surfaces(store, monkeypatch):
    """The CLI backend's normal state: no session id at all."""
    store.create("left behind", session_id="sess-a")
    assert "left behind" in _block_for(monkeypatch, store.base_dir, "")
    assert "left behind" in TT.render_html(store.base_dir, session_id="")


def test_the_resolver_is_the_single_definition():
    assert resolve_session_scope("") is None
    assert resolve_session_scope(None) is None
    assert resolve_session_scope("   ") is None
    assert resolve_session_scope("sess-a") == "sess-a"


def test_the_store_treats_empty_and_none_alike(tmp_path):
    s = TaskStore(tmp_path)
    s.create("a", session_id="x")
    s.create("b", session_id="y")
    assert len(s.list(session_id="")) == len(s.list(session_id=None)) == 2
    assert len(s.list(session_id="x")) == 1


def test_a_scoped_session_still_sees_only_its_own(store, monkeypatch):
    store.create("mine", session_id="sess-a")
    store.create("theirs", session_id="sess-b")
    block = _block_for(monkeypatch, store.base_dir, "sess-a")
    html = TT.render_html(store.base_dir, session_id="sess-a")
    assert "mine" in block and "theirs" not in block
    assert "mine" in html and "theirs" not in html
