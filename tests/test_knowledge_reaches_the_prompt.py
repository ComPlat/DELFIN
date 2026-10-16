"""The read path in the composed prompt: repository lessons in the
volatile tail when they match the task, fenced as data; swarm findings at
each round boundary of every agent of the run; the stable head unchanged.
"""
from __future__ import annotations

from types import SimpleNamespace

import pytest

from delfin.agent import knowledge as kn
from delfin.agent import untrusted


@pytest.fixture(autouse=True)
def env(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("DELFIN_SWARM_RUN", raising=False)
    gdir = tmp_path / "pack" / "knowledge"
    kn.write_topics(gdir, {t: "" for t in kn.DEFAULT_TOPICS})
    kn.write_entry(gdir, kn.Entry(
        slug="sbatch-modules", title="Load modules inside the job script",
        topic="cluster-hpc", scope="always",
        provenance=[kn.Provenance("anna", "2026-09-01", "job:1")],
        review_by="2099-01-01",
        situation="A SLURM job submitted with sbatch from a login node.",
        lesson="Load the module environment inside the job script.",
        why="The job does not inherit the login shell's modules."))
    kn.write_index(gdir)
    monkeypatch.setattr(kn, "global_dir", lambda: gdir)
    ws = tmp_path / "ws"
    ws.mkdir()
    return {"ws": ws, "gdir": gdir}


def _loader(ws):
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    loader.workspace_root = ws
    loader.is_delfin_workspace = False
    return loader


def _sections(loader, task, role_id="solo_agent", mode_id="solo"):
    secs = loader.compose_sections(role_id=role_id, mode_id=mode_id,
                                   task_text=task, session_key="kn-1")
    return {s.name: s for s in secs}


def test_a_matching_lesson_reaches_the_volatile_tail(env):
    loader = _loader(env["ws"])
    secs = _sections(loader, "my sbatch job cannot find the module")
    assert "knowledge" in secs
    sec = secs["knowledge"]
    assert sec.layer == loader.LAYER_VOLATILE
    text = sec.content
    assert "Load modules inside the job script" in text
    assert "[UNTRUSTED EXTERNAL CONTENT" in text


def test_the_dashboard_role_reads_it_too(env):
    secs = _sections(_loader(env["ws"]), "sbatch job cannot find the module",
                     role_id="dashboard_agent", mode_id="dashboard")
    assert "knowledge" in secs


def test_an_unrelated_task_carries_no_lesson(env):
    secs = _sections(_loader(env["ws"]), "rename a variable in utils.py")
    assert "knowledge" not in secs or not secs["knowledge"].content.strip()


def test_bare_mode_reads_no_knowledge(env):
    loader = _loader(env["ws"])
    loader.skip_external_memory = True
    secs = _sections(loader, "my sbatch job cannot find the module")
    assert "knowledge" not in secs or not secs["knowledge"].content.strip()


def test_findings_arrive_at_each_round_for_the_other_agents(env):
    from delfin.agent import swarm_memory as sm
    perms_a = SimpleNamespace(swarm_run_id="wave-5", subagent_id="sa-a")
    perms_b = SimpleNamespace(swarm_run_id="wave-5", subagent_id="sa-b")
    seen_b: dict = {}
    assert sm.news_for(perms_b, seen_b) == ""
    sm.append("wave-5", author="sa-a", text="fixture tmp_home is shared",
              evidence="test:tests/conftest.py")
    block = sm.news_for(perms_b, seen_b)
    assert "fixture tmp_home is shared" in untrusted.unwrap(block)
    assert sm.news_for(perms_b, seen_b) == ""          # each entry once
    assert sm.news_for(perms_a, {}) == ""              # not its own


def test_a_delegate_inherits_the_matching_lessons(env):
    from delfin.agent.subagents import _inherited_context
    out = _inherited_context(env["ws"], "sbatch job cannot find the module")
    assert "Load modules inside the job script" in out
    assert "[UNTRUSTED EXTERNAL CONTENT" in out


def test_the_round_boundary_feeds_the_swarm_news():
    """The tool loop reads the run's new findings where it reads the
    delegation notes: once per round, for parent and delegates alike."""
    import inspect
    from delfin.agent import api_client as A
    src = inspect.getsource(A)
    at = src.index('_note_src = getattr(self, "_subagent_notes", None)')
    window = src[at - 4000:at + 4000]
    assert "news_for(" in window
    assert "not from the user" in window


def test_subagents_of_one_session_share_one_run(env):
    from delfin.agent import api_client as A
    from delfin.agent.subagents import _derive_perms
    parent = A.KitToolPermissions(workspace=env["ws"], mode="default")
    parent.task_session_id = "sess-42"
    child = _derive_perms(parent, "default")
    assert parent.swarm_run_id == "s-sess-42"
    assert child.swarm_run_id == parent.swarm_run_id
