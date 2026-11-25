"""`forget` must delete a memory that a later session still reads.

The recall path (prompt_loader) merges TWO stores into the External Memory
block: ``~/.delfin/memory`` (the user-wide store, injected first) and the
per-project store. The model therefore SEES global memories in its recalled
block and has their names to hand.

The ``forget`` tool deleted from the project store only -- ``scope`` was
never passed, so it defaulted to ``"project"``. For a global fact the call
answered "no memory named ..." and deleted nothing, while the fact stayed in
the prompt of every later session in every workspace. The tool that exists
to keep the store truthful was the one place that could not tell the truth
about it.

``/forget`` in the REPL passes ``scope="all"``; the dashboard's ``/forget``
searches the project store and then falls back to the global one. The tool
the AGENT calls was the only entry point that reached one store.

These tests drive the product's own tool dispatch -- ``_doc_executor.
execute(name, arguments, permissions)``, the path the model's call takes --
not the handler, so a future change to the dispatch is covered too.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions

PHRASE = "identity anchor phrase eta-six"


@pytest.fixture()
def home(tmp_path, monkeypatch):
    """A home of our own: the stores resolve off ``Path.home()``."""
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setattr(Path, "home", lambda: h)
    return h


@pytest.fixture()
def workspaces(tmp_path):
    a = tmp_path / "project_a"
    b = tmp_path / "project_b"
    a.mkdir()
    b.mkdir()
    return a, b


def _perms(ws: Path) -> KitToolPermissions:
    return KitToolPermissions(workspace=ws, mode="bypassPermissions")


def _call(tool: str, args: dict, ws: Path) -> dict:
    """One tool call through the dispatch the model's call goes through."""
    return json.loads(A._doc_executor.execute(tool, args, _perms(ws)))


def _recalled(workspace: Path) -> str:
    """What a session in `workspace` reads into its prompt."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader()
    loader.workspace_root = workspace
    return loader._load_external_memory_context(max_chars=6000, task_text="")


def test_forget_removes_a_global_memory_it_can_see(home, workspaces):
    """The defect: the name is in the recalled block, the delete misses."""
    from delfin.agent import memory_store as ms

    ws_a, ws_b = workspaces
    ms.save_typed_memory(f"global: user: {PHRASE}", repo_root=ws_a,
                         allow_scope_prefix=True)
    slug = ms.list_typed_memories(ws_a, scope="user")[0]["name"]

    # Precondition: the model really can see this name and this fact --
    # from a workspace that did not write it.
    assert PHRASE in _recalled(ws_b)

    out = _call("forget", {"name": slug}, ws_b)

    assert out.get("status") == "deleted", out
    assert ms.list_typed_memories(ws_a, scope="user") == []
    assert PHRASE not in _recalled(ws_b), (
        "the global fact is still recalled in another workspace after a "
        "forget that reported success")


def test_forget_still_removes_a_project_memory(home, workspaces):
    """The control: the behaviour that already worked must keep working.

    Green on origin/main and green after the fix -- a change that broke the
    project store would pass the test above and fail this one.
    """
    from delfin.agent import memory_store as ms

    ws_a, _ws_b = workspaces
    out = _call("remember", {"text": "project: the flag --fast exists"}, ws_a)
    slug = out["slug"]
    assert ms.list_typed_memories(ws_a)

    gone = _call("forget", {"name": slug}, ws_a)

    assert gone.get("status") == "deleted", gone
    assert ms.list_typed_memories(ws_a) == []


def test_forget_of_a_name_in_neither_store_still_says_so(home, workspaces):
    """A scope="all" search must not turn a real miss into a silent success.

    Green on origin/main too -- it pins that widening the search did not
    make the tool report a deletion that did not happen.
    """
    _ws_a, ws_b = workspaces
    out = _call("forget", {"name": "never-existed-anywhere"}, ws_b)
    assert "error" in out, out
    assert "never-existed-anywhere" in out["error"]
