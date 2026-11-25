"""A wave member's commits carry no attribution trailer unless deployed to.

Measured on w48 while writing the landing rules: the product appends two
trailer lines to every ``git commit`` a model bash call runs --

    Co-Authored-By: DELFIN-Agent (<model>) <email>     (api_client._commit_trailer_value)
    Agent: <run>/<id>                                  (DELFIN_AGENT_LABEL, wave._launch)

-- and instructs the model to add the first itself (``prompt_loader.
_commit_attribution_block``). A reviewer of a merged branch then reads the
product name and the model alias as public text, and the landing rules say a
commit message is public text; a merged history cannot be reworded. The
decision therefore lives in ``delfin/agent/commit_attribution.py``: off for a
session inside a wave (``agent.swarm_commit_coauthor``, default false),
unchanged for every other session (``agent.commit_coauthor``, default true).

The executor splice and the prompt block are wired to it through
scratch/h1-trailer-off-for-members.patch (both files sit behind the
self-modification guard; wave rule 13). These tests pin the decision itself,
including the case that must NOT change: a plain session's sub-agent keeps its
attribution even though its permissions carry a swarm run id.
"""
from __future__ import annotations

import pytest

from delfin.agent import commit_attribution as ca


@pytest.fixture
def member(monkeypatch):
    """A process started as a member of a wave (what wave._launch writes)."""
    monkeypatch.setenv(ca.ENV_RUN, "w48")
    return None


@pytest.fixture
def plain(monkeypatch):
    monkeypatch.delenv(ca.ENV_RUN, raising=False)
    return None


def test_a_wave_member_is_off_by_default(member):
    assert ca.in_swarm_run() is True
    assert ca.commit_coauthor_for({}) is False
    assert ca.commit_coauthor_for(None) is False


def test_a_deployment_can_switch_the_member_back_on(member):
    """The knob is the point: configurable, not removed."""
    assert ca.commit_coauthor_for({"swarm_commit_coauthor": True}) is True


def test_an_ordinary_session_is_untouched(plain):
    assert ca.in_swarm_run() is False
    assert ca.commit_coauthor_for({}) is True
    assert ca.commit_coauthor_for({"commit_coauthor": False}) is False


def test_the_member_knob_does_not_leak_into_the_ordinary_answer(plain):
    """A member's key says nothing about a normal session, and the other way
    round: the two questions stay two keys. Fails if one key is read for both
    session kinds."""
    assert ca.commit_coauthor_for({"swarm_commit_coauthor": True}) is True
    assert ca.commit_coauthor_for({"commit_coauthor": False,
                                   "swarm_commit_coauthor": True}) is False


def test_the_member_answer_ignores_the_ordinary_key(member):
    """Turning the ordinary key off must not be the only way to silence a
    member (it is off anyway), and turning it on must not switch the member on:
    that is what would happen if the member branch read the wrong key."""
    assert ca.commit_coauthor_for({"commit_coauthor": True}) is False


def test_the_run_marker_is_the_one_the_swarm_writes():
    """The module names DELFIN_SWARM_RUN itself (importing swarm_memory here
    would drag the store's dependencies into every prompt build). This test is
    the reason that is safe: if the swarm ever changes its marker, this goes
    red instead of the member check silently matching nothing."""
    from delfin.agent import swarm_memory
    from delfin.agent.swarm import wave, job
    assert ca.ENV_RUN == swarm_memory.ENV_RUN == "DELFIN_SWARM_RUN"
    src = __import__("pathlib").Path(wave.__file__).read_text()
    assert f"{ca.ENV_RUN}=" in src, "wave._launch must still write the marker"
    assert ca.ENV_RUN in __import__("pathlib").Path(job.__file__).read_text()


def test_a_delegated_subagent_of_a_plain_session_is_not_a_member(
        plain, monkeypatch):
    """The trap this decision walks around: swarm_memory.current_run_id falls
    back to perms.swarm_run_id, which subagents._derive_perms sets for a plain
    session's child. Using it here would strip attribution from every
    delegation in the product. A child is not a member; the env marker is the
    line."""
    from delfin.agent import swarm_memory

    class _Perms:
        swarm_run_id = "s-derived"

    assert swarm_memory.current_run_id(_Perms()) == "s-derived"
    assert ca.in_swarm_run() is False
    assert ca.commit_coauthor_for({}) is True
