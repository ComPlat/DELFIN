"""A greeting came back as fifteen tasks from other sessions.

Observed live, in a terminal session against the KIT backend: the user
typed a greeting and the turn ended with a list of open `validator_kit`
tasks that belonged to earlier sessions of the same workspace.

The scoping itself was sound. `open_task_summary` filters on
`perms.task_session_id`, and an EMPTY id means unscoped — the whole
workspace — which is documented and correct: `task_list`, `task_adopt`
and the model-facing reminder all rely on it.

What was wrong is how the id came to be empty. `_fresh_session_id` left
it empty for `backend == "cli"`, because the CLI subprocess backend
announces its own id on the stream. But `create_client` routes on the
PROVIDER: with `provider="kit"` it returns an OpenAIClient no matter what
`backend` says, and that client announces nothing. So the one backend
string that was allowed to defer the id was also reachable by a client
that never supplies one, and the deferral never ended.

The question is now asked of the object that was built.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from delfin.agent import agent_tasks as at
from delfin.agent.engine import AgentEngine


@pytest.fixture
def workspace_with_parked_work():
    """A workspace holding fifteen open tasks owned by an earlier session."""
    tmp = Path(tempfile.mkdtemp(prefix="task-scope-"))
    ws = tmp / "agent_workspace"
    ws.mkdir()
    store = at.get_store(ws)
    for i in range(15):
        store.create(subject=f"validator_kit step {i}", session_id="earlier")
    return ws


def _engine(ws: Path, backend: str) -> AgentEngine:
    # A key, not THE key: the client is built but never called here, and
    # the suite must not depend on a credential store being populated.
    return AgentEngine(repo_dir=ws, backend=backend, provider="kit",
                       api_key="test-key-not-used",
                       model="kit.qwen3.5-397b-A17b", mode="solo")


def test_a_kit_session_is_scoped_whatever_the_backend_string_says(
        workspace_with_parked_work):
    """The reproducer, at the level the defect lives on.

    `--backend cli --provider kit` is a real invocation: the backend
    default comes from a settings file, and nothing rejects the pair.
    """
    eng = _engine(workspace_with_parked_work, "cli")

    assert type(eng.client).__name__ == "OpenAIClient", (
        "the premise of this test is that the provider wins the routing")
    assert eng.session_id, (
        "a client that announces no session id must be given one, or every "
        "session in this workspace shares the unscoped bucket")

    perms = eng.kit_permissions
    assert getattr(perms, "task_session_id", "") == eng.session_id

    summary = at.open_task_summary(perms.workspace, perms.task_session_id)
    assert summary["state"] == "none", (
        f"fifteen tasks from another session were reported as this "
        f"session's unfinished work: {summary['counts']}")


def test_the_api_backend_was_never_affected(workspace_with_parked_work):
    """The control. It has always minted an id; this must not change."""
    eng = _engine(workspace_with_parked_work, "api")
    assert eng.session_id
    summary = at.open_task_summary(
        eng.kit_permissions.workspace, eng.kit_permissions.task_session_id)
    assert summary["state"] == "none"


def test_a_client_that_announces_its_id_still_defers(monkeypatch,
                                                     workspace_with_parked_work):
    """The deferral is the reason the empty id exists at all.

    Minting one for a backend that is about to be told a different id
    would split one conversation across two task buckets.
    """
    eng = _engine(workspace_with_parked_work, "api")
    monkeypatch.setattr(type(eng.client), "supplies_session_id", True,
                        raising=False)
    assert eng._fresh_session_id() == ""


def test_the_two_announcing_clients_declare_it():
    """Both backends that emit `session_init` say so on the class.

    A capability read off the class is checkable here; one inferred from a
    string at the call site is the thing that broke.
    """
    from delfin.agent import api_client as A

    assert A.CLIClient.supplies_session_id is True
    assert A.CodexCLIClient.supplies_session_id is True
    assert A.OpenAIClient.supplies_session_id is False
    assert A.APIClient.supplies_session_id is False


def test_an_empty_id_still_means_the_whole_workspace():
    """Unscoped is deliberate, and stays deliberate.

    This is the behaviour the fix must NOT change: `task_list` and
    `task_adopt` are documented to show everything when there is no
    session to scope to.
    """
    assert at.resolve_session_scope("") is None
    assert at.resolve_session_scope(None) is None
    assert at.resolve_session_scope("sid") == "sid"
