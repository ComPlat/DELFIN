"""The executor must ask the same question the prompt asks about attribution.

Measured while landing package h1: ``delfin/agent/commit_attribution.py`` says
one module answers "does this session attribute its commits?" for both the
prompt block and the mechanical splice, and the prompt half does ask it. The
executor half did not: ``api_client._commit_trailer_enabled`` read
``agent.commit_coauthor`` directly, so a wave member was told nothing about the
trailer while its ``git commit`` still got ``Co-Authored-By:`` and
``Agent: <run>/<id>`` spliced in. A reviewer of the merged history then reads
the product name, the model alias and the run id as public text -- which is the
defect landing rule E names.

This case is the contract for that second half. It is red today because the
executor is still wired to the bare key; the wiring lives in
``tests/fixtures/h1-trailer-off-for-members.patch`` (``api_client.py`` and
``user_settings.py`` both sit behind the self-modification guard, wave rule 13,
so the diff is carried in the tree for a session that may apply it).

It is written as ``xfail(strict=True)`` on purpose: the moment the patch is
applied the executor answers "no" and this case XPASSes, which fails the suite
and tells whoever applied the patch to delete the marker. A plain skip or a
non-strict xfail would go on passing after the fix and quietly measure nothing
-- the repository's convention for a contract that predates its wiring
(tests/test_elided_tool_result_digest.py).
"""
from __future__ import annotations

import os
import subprocess

import pytest

from delfin.agent import api_client as A

#: The marker ``wave._launch`` writes into a member's environment. Spelled here
#: rather than imported, so this file measures the executor's own behaviour on a
#: tree that has neither the decision module nor the wiring: on such a tree the
#: import would raise, and a collection error is not the failure this case names.
ENV_RUN = "DELFIN_SWARM_RUN"


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """A git repository and executor permissions, both session markers pinned.

    Mirrors tests/test_commit_trailer.py's fixture: ``DELFIN_AGENT_LABEL``
    unset so the ``Agent:`` line cannot confuse the count, and the git global
    config neutralised so a host's own trailer setting cannot leak in.
    """
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    monkeypatch.delenv("DELFIN_AGENT_LABEL", raising=False)
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", "/dev/null")
    monkeypatch.setenv("GIT_CONFIG_SYSTEM", "/dev/null")
    repo = tmp_path / "repo"
    repo.mkdir()
    env = {"GIT_CONFIG_GLOBAL": "/dev/null", "GIT_CONFIG_SYSTEM": "/dev/null",
           "PATH": os.environ["PATH"]}
    run = lambda *a: subprocess.run(["git", *a], cwd=repo, env=env,
                                    check=True, capture_output=True)
    run("init", "-q", "-b", "main")
    run("-c", "user.email=t@t", "-c", "user.name=t",
        "commit", "-q", "--allow-empty", "-m", "init")
    return repo


def _agent_commit(repo, monkeypatch, model="kit.glm-5.3-flash",
                  settings=None, name="f.txt"):
    """One model bash call through the real executor; the last commit message."""
    import delfin.user_settings as us
    base = {"agent": {"bash_isolation": "off"}}
    if settings:
        base["agent"].update(settings["agent"])
    monkeypatch.setattr(us, "load_settings", lambda *a, **k: base)
    A.set_commit_trailer_model(model)
    try:
        perms = A.KitToolPermissions(workspace=str(repo), mode="default",
                                     confirm_callback=lambda *a, **k: True)
        result = A._doc_executor.execute(
            "bash",
            {"command": f"echo hi > {name} && git add {name} && "
                        "git -c user.email=t@t -c user.name=t commit -q -m work"},
            perms)
        assert "error" not in result, result
    finally:
        A.set_commit_trailer_model("")
    msg = subprocess.run(["git", "log", "-1", "--format=%B"], cwd=repo,
                         capture_output=True, text=True).stdout
    return msg


@pytest.mark.xfail(
    strict=True,
    reason="The executor still reads agent.commit_coauthor directly instead of "
           "commit_attribution.commit_coauthor_for; the wiring is "
           "tests/fixtures/h1-trailer-off-for-members.patch. Delete this "
           "marker with that patch (wave rule 13).")
def test_the_executor_obeys_the_attribution_decision(repo, monkeypatch):
    """Both sides of the switch, through the splice a member's commit goes through.

    The knob has to stay a knob: the patch cannot simply delete the line, so
    the deployment-on case is pinned beside the default-off one. One case, so a
    failure names exactly one contract to satisfy.
    """
    monkeypatch.setenv(ENV_RUN, "w48")
    off = _agent_commit(repo, monkeypatch, name="off.txt")
    assert "Co-Authored-By:" not in off, off
    assert "Agent:" not in off, off

    on = _agent_commit(repo, monkeypatch, name="on.txt",
                       settings={"agent": {"swarm_commit_coauthor": True}})
    assert "Co-Authored-By:" in on, on


def test_an_ordinary_session_still_gets_it(repo, monkeypatch):
    """The never-worse case, green on both trees: the patch must not touch it.

    If this fails, the change reached a session it has no business changing.
    """
    monkeypatch.delenv(ENV_RUN, raising=False)
    msg = _agent_commit(repo, monkeypatch)
    assert "Co-Authored-By:" in msg, msg
    assert "kit.glm-5.3-flash" in msg, msg
