"""Tests for delfin.agent.episodes — episodic session recall.

The session store used to be write-only: nothing recalled a finished
session. These tests cover the episode record writer, the deterministic
state-derivation, BM25 recall, pruning, the settings gate and the
headless-CLI save hook.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from delfin.agent import episodes as ep


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """Fake home + a repo directory, so the store lands under tmp."""
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    r = tmp_path / "repo"
    r.mkdir()
    return r


# ---------------------------------------------------------------------------
# save_episode
# ---------------------------------------------------------------------------

def test_save_episode_writes_frontmatter_and_sections(repo):
    p = ep.save_episode(
        "abc12345-6789",
        repo_root=repo,
        goal="Fix the ORCA convergence failure",
        outcome="Patched the input builder; tests pass",
        decisions=["Used tightscf keywords"],
        open_items=["#3 verify cluster run (pending)"],
        cost_usd=0.1234,
        verdict="PASS",
    )
    assert p.is_file()
    assert p.parent.name == "episodes"
    today = time.strftime("%Y-%m-%d")
    assert p.name == f"{today}_abc12345.md"
    text = p.read_text(encoding="utf-8")
    assert "session_id: abc12345-6789" in text
    assert f"date: {today}" in text
    assert "verdict: PASS" in text
    assert "cost: 0.1234" in text
    assert "## Goal\nFix the ORCA convergence failure" in text
    assert "## Outcome\nPatched the input builder; tests pass" in text
    assert "- Used tightscf keywords" in text
    assert "- #3 verify cluster run (pending)" in text


def test_save_episode_same_session_same_day_overwrites(repo):
    ep.save_episode(
        "samesid1", repo_root=repo, goal="first pass", outcome="",
        decisions=[], open_items=[],
    )
    p = ep.save_episode(
        "samesid1", repo_root=repo, goal="second pass", outcome="done",
        decisions=[], open_items=[],
    )
    files = list(ep._episodes_dir(repo).glob("*.md"))
    assert len(files) == 1
    assert "second pass" in p.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# build_episode_from_state
# ---------------------------------------------------------------------------

def test_build_episode_from_state_derives_fields():
    chat = [
        {"role": "user", "content": "Please fix the failing isomer test"},
        {"role": "assistant",
         "content": "I will look at the test first.\n\nReading tests now."},
        {"role": "assistant",
         "content": "Root cause found in isomer.py.\n\n"
                    "Patched the deterministic flag; the test passes."},
    ]
    state = {
        "cost_usd": 0.5,
        "todo_payload": [
            {"id": 1, "subject": "fix test", "status": "completed"},
            {"id": 2, "subject": "run suite", "status": "pending"},
            {"id": 3, "subject": "report", "status": "in_progress"},
        ],
    }
    fields = ep.build_episode_from_state(state, chat)
    assert fields["goal"] == "Please fix the failing isomer test"
    assert fields["outcome"] == (
        "Patched the deterministic flag; the test passes.")
    assert fields["decisions"] == [
        "I will look at the test first.",
        "Root cause found in isomer.py.",
    ]
    assert fields["open_items"] == [
        "#2 run suite (pending)",
        "#3 report (in_progress)",
    ]
    assert fields["cost_usd"] == 0.5


def test_build_episode_falls_back_to_engine_messages():
    """Headless CLI saves have no chat list — engine messages (with
    content-block lists) must still yield goal + outcome."""
    state = {
        "engine_messages": [
            {"role": "user",
             "content": [{"type": "text", "text": "summarise the repo"}]},
            {"role": "assistant",
             "content": [
                 {"type": "text", "text": "The repo builds ORCA inputs."},
             ]},
        ],
    }
    fields = ep.build_episode_from_state(state, [])
    assert fields["goal"] == "summarise the repo"
    assert fields["outcome"] == "The repo builds ORCA inputs."
    assert fields["open_items"] == []
    assert fields["cost_usd"] == 0.0


# ---------------------------------------------------------------------------
# recall_episodes
# ---------------------------------------------------------------------------

def _seed_two_episodes(repo):
    ep.save_episode(
        "aaaa1111", repo_root=repo,
        goal="Improve dashboard websocket latency",
        outcome="Tuned the polling interval",
        decisions=[], open_items=[],
    )
    ep.save_episode(
        "bbbb2222", repo_root=repo,
        goal="Debug orca geometry convergence failure",
        outcome="Added tightscf slowconv keywords",
        decisions=[], open_items=["rerun on cluster"], verdict="PASS",
    )


def test_recall_ranks_matching_episode_first(repo):
    _seed_two_episodes(repo)
    block = ep.recall_episodes(
        repo, "the orca geometry convergence job fails again")
    lines = block.splitlines()
    assert lines[0] == "# Similar past sessions"
    assert "orca geometry convergence" in lines[1]
    assert "(PASS)" in lines[1]
    assert "tightscf" in lines[1]
    assert "[open: rerun on cluster]" in lines[1]
    # The unrelated websocket episode has zero overlap — excluded.
    assert "websocket" not in block


def test_recall_respects_max_entries_and_max_chars(repo):
    # Two episodes that BOTH match the task text.
    ep.save_episode(
        "cccc3333", repo_root=repo,
        goal="orca convergence run one", outcome="tuned scf settings",
        decisions=[], open_items=[],
    )
    ep.save_episode(
        "dddd4444", repo_root=repo,
        goal="orca convergence run two", outcome="tuned grid settings",
        decisions=[], open_items=[],
    )
    block = ep.recall_episodes(
        repo, "another orca convergence problem", max_entries=1)
    assert len(block.splitlines()) == 2  # header + exactly one entry

    capped = ep.recall_episodes(
        repo, "another orca convergence problem", max_chars=80)
    assert capped
    assert len(capped) <= 80


def test_recall_empty_store_returns_empty(repo):
    assert ep.recall_episodes(repo, "anything at all") == ""


def test_recall_no_match_or_empty_task_returns_empty(repo):
    _seed_two_episodes(repo)
    assert ep.recall_episodes(
        repo, "quantum espresso phonon dispersion bands") == ""
    assert ep.recall_episodes(repo, "") == ""


# ---------------------------------------------------------------------------
# pruning
# ---------------------------------------------------------------------------

def test_prune_keeps_newest_100(repo):
    d = ep._episodes_dir(repo)
    d.mkdir(parents=True)
    for i in range(104):
        (d / f"old_{i:03d}.md").write_text(
            "---\n"
            f"session_id: s{i:03d}\n"
            f"date: 2019-{i:03d}\n"
            "verdict:\n"
            "cost: 0.0\n"
            "---\n\n"
            f"## Goal\nold goal {i}\n",
            encoding="utf-8",
        )
    p = ep.save_episode(
        "newsid99", repo_root=repo, goal="fresh work", outcome="",
        decisions=[], open_items=[],
    )
    files = {f.name for f in d.glob("*.md")}
    assert len(files) == 100
    assert p.name in files
    # The five oldest records were evicted; the newest old ones survive.
    for i in range(5):
        assert f"old_{i:03d}.md" not in files
    assert "old_103.md" in files


# ---------------------------------------------------------------------------
# settings gate (prompt loader unit)
# ---------------------------------------------------------------------------

def test_loader_episode_recall_respects_enabled_setting(repo, monkeypatch):
    from delfin.agent.prompt_loader import PromptLoader
    import delfin.user_settings as us

    ep.save_episode(
        "gate1234", repo_root=repo,
        goal="Debug orca geometry convergence failure",
        outcome="Added tightscf keywords", decisions=[], open_items=[],
    )
    loader = PromptLoader(repo)

    monkeypatch.setattr(
        us, "load_settings",
        lambda *a, **k: {"agent": {"episodes": {"enabled": True}}})
    out = loader._load_episode_recall_context(
        "the orca geometry convergence job fails")
    assert "Similar past sessions" in out

    monkeypatch.setattr(
        us, "load_settings",
        lambda *a, **k: {"agent": {"episodes": {"enabled": False}}})
    assert loader._load_episode_recall_context(
        "the orca geometry convergence job fails") == ""


# ---------------------------------------------------------------------------
# CLI hook
# ---------------------------------------------------------------------------

def test_cli_run_saves_episode(repo, monkeypatch, capsys):
    from delfin.agent import cli as agent_cli

    monkeypatch.chdir(repo)

    class FakeEngine:
        session_id = ""

        def export_state(self):
            return {
                "engine_messages": [
                    {"role": "user", "content": "check the failing tests"},
                    {"role": "assistant",
                     "content": "Two tests fail.\n\n"
                                "Both predate this branch."},
                ],
                "cost_usd": 0.02,
            }

        def record_cycle_outcome(self, *args, **kwargs):
            return None

    monkeypatch.setattr(
        agent_cli, "_build_engine", lambda args: FakeEngine())
    monkeypatch.setattr(
        agent_cli, "_run_once",
        lambda engine, prompt, max_tokens=4096: {
            "text": "done", "tool_calls": [],
            "input_tokens": 1, "output_tokens": 1, "error": "",
        })
    monkeypatch.setattr(
        agent_cli, "_save_session",
        lambda engine, repo_root: "epis-cli-1234")
    monkeypatch.setattr(
        "delfin.agent.eval_loop.maybe_run_scheduled",
        lambda *a, **k: None)

    rc = agent_cli.main(["run", "--cwd", str(repo), "check the failing tests"])
    assert rc == 0

    files = list(ep._episodes_dir(repo).glob("*.md"))
    assert len(files) == 1
    text = files[0].read_text(encoding="utf-8")
    assert "session_id: epis-cli-1234" in text
    assert "verdict: PASS" in text
    assert "check the failing tests" in text
    assert "Both predate this branch." in text
