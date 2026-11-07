"""The start form offers an own worktree in a plain folder, with consent.

"Own worktree" appeared only inside a git repository. A person's default
workspace is usually a plain folder, so the box is offered there too --
labelled with what it costs: the folder becomes a local git project. The
start asks once (the same consent the engine's NEEDS_LOCAL_PROJECT window
asks) and only with a yes makes the folder a repository and starts in a
worktree of it. A session worktree's changes come back by merge, shown --
never overwritten -- so a per-session "Merge back" action exists.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from delfin.dashboard import agent_sessions as AS

widgets = pytest.importorskip("ipywidgets")


def _git(repo, *args) -> str:
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout


def _form(tmp_path, monkeypatch, home):
    """A session tab whose new-session form is open, pointed at ``home``."""
    from delfin.agent import scheduler as S
    from delfin.dashboard import agent_sessions as AS
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    (tmp_path / "open.json").write_text("[]")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(home))
    # No session opens by itself (create_tab opens the default workspace
    # otherwise): the tests drive the start form only.
    monkeypatch.setattr(AS, "default_workspace", lambda ctx: str(home))

    def build(ctx):
        state = {"active_session_id": "", "chat_messages": [],
                 "streaming": False, "engine": None}
        return widgets.VBox(), {"state": state, "shutdown": lambda: None}

    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           agent_dir=str(home))
    _widget, refs = AS.create_tab(ctx, build=build)
    # create_tab opens the default workspace when nothing is open; the
    # tests only judge what the form adds, so the pre-opened one is closed.
    for rec in list(refs["sessions"]()):
        refs["close"](rec["key"])
    form = refs["form"]
    form["new"].click()
    return refs, form


@pytest.fixture
def plain(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    ws = home / "people" / "someone" / "agent_workspace"
    ws.mkdir(parents=True)
    (ws / "notes.md").write_text("# notes\n", encoding="utf-8")
    return ws


def test_the_box_is_offered_in_a_plain_folder_with_what_it_costs(
        plain, tmp_path, monkeypatch):
    refs, form = _form(tmp_path, monkeypatch, plain)
    form["workdir"].value = str(plain)
    box = form["own_worktree"]
    assert box.layout.display != "none"
    assert box.description != "Own worktree"          # says the cost
    assert "local project" in box.description
    # Unchanged: inside a repository the label stays the old one.
    assert (AS.OWN_WORKTREE_LABEL_IN_REPO == "Own worktree")


def _worktree_sessions(refs):
    """The sessions the form added: create_tab always keeps one default
    session open (closing the last re-opens it), so the tests pick the
    sessions whose workspace is a session worktree."""
    from delfin.dashboard import agent_sessions as AS
    return [rec for rec in refs["sessions"]()
            if AS.read_worktree_sidecar(rec["workspace"]) is not None
            or ".delfin/worktrees/" in str(rec["workspace"])]


def test_a_start_without_a_yes_changes_nothing(plain, tmp_path, monkeypatch):
    """No consent recorded, no repository, no worktree session."""
    refs, form = _form(tmp_path, monkeypatch, plain)
    form["workdir"].value = str(plain)
    form["own_worktree"].value = True
    form["start"].click()          # the confirm row appears, nothing done
    assert not (plain / ".git").exists()
    assert _worktree_sessions(refs) == []
    form["consent_no"].click()
    assert not (plain / ".git").exists()
    assert _worktree_sessions(refs) == []


def test_a_yes_makes_the_folder_a_project_and_starts_in_its_worktree(
        plain, tmp_path, monkeypatch):
    refs, form = _form(tmp_path, monkeypatch, plain)
    form["workdir"].value = str(plain)
    form["own_worktree"].value = True
    form["start"].click()          # asks: the confirm row
    form["consent_yes"].click()    # the yes
    sessions = _worktree_sessions(refs)
    assert len(sessions) == 1
    where = sessions[0]["workspace"]
    assert where.startswith(str(plain / ".delfin" / "worktrees"))
    assert _git(plain, "status", "--porcelain").strip() == ""
    assert _git(plain, "log", "--oneline").strip()
    for rec in sessions:
        refs["close"](rec["key"])


def test_merge_back_brings_the_changes_and_reports_a_conflict(
        plain, tmp_path, monkeypatch):
    """Two sessions on their own worktrees; the same file edited in both:
    the first merge lands in the main checkout, the second is reported as
    a conflict and the file is NOT overwritten."""
    refs, form = _form(tmp_path, monkeypatch, plain)
    form["workdir"].value = str(plain)
    form["own_worktree"].value = True
    form["start"].click()
    form["consent_yes"].click()
    form["own_worktree"].value = True
    form["start"].click()          # already a project: no consent needed
    wt_a, wt_b = [Path(rec["workspace"]) for rec in _worktree_sessions(refs)]
    (wt_a / "notes.md").write_text("from a\n", encoding="utf-8")
    (wt_b / "notes.md").write_text("from b\n", encoding="utf-8")
    recs = _worktree_sessions(refs)
    buttons = {rec["workspace"]: rec.get("merge_btn") for rec in recs}
    assert all(buttons.values()), "every session worktree row offers merge"
    buttons[str(wt_a)].click()
    assert "Merged back" in refs["notice"].value, refs["notice"].value
    # The button did the whole merge: the change is in the main checkout,
    # uncommitted, for review.
    assert (plain / "notes.md").read_text(encoding="utf-8") == "from a\n"
    assert "notes.md" in _git(plain, "status", "--porcelain")

    conflict = AS.merge_session_worktree(wt_b)
    assert conflict["status"] == "conflict"
    assert "notes.md" in conflict["message"]
    assert (plain / "notes.md").read_text(encoding="utf-8") == "from a\n"
