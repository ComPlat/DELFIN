"""A workspace folder without git becomes a local project on the user's yes.

Worktrees, sub-agents in their own worktrees and their merges need git; a
person's default workspace is a plain folder. The chain checked here: no
init without consent -> consent -> git init + .gitignore + first commit
(QM output, large files and secrets kept out) -> two sub-agents of one
orchestration in their own worktrees -> merged back, a conflict reported
and never overwritten. Refused: the home directory, a folder inside
another repository, system paths.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent import local_project as LP
from delfin.agent import subagents as SA
from delfin.agent import worktree as WT
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

ORCA_TEXT = ("                                 * O   R   C   A *\n"
             "                     ****ORCA TERMINATED NORMALLY****\n")


def _git(repo, *args) -> str:
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout


def _plain(tmp_path, monkeypatch) -> Path:
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    ws = home / "people" / "someone" / "agent_workspace"
    ws.mkdir(parents=True)
    (ws / "notes.md").write_text("# notes\n")
    (ws / "shared.py").write_text("A = 1\nB = 2\n")
    return ws


def _perms(ws, answer=None, mode="acceptEdits"):
    perms = KitToolPermissions(workspace=ws, mode=mode)
    perms.task_session_id = "local-project-test"
    asked: list = []
    if answer is not None:
        def _cb(name, args, preview):
            asked.append((name, preview))
            return answer
        perms.confirm_callback = _cb
    return perms, asked


def _commit_tree(ws) -> list[str]:
    return sorted(_git(ws, "ls-files").split())


# -- consent -----------------------------------------------------------------

def test_nothing_is_made_a_repository_without_someone_to_ask(tmp_path, monkeypatch):
    ws = _plain(tmp_path, monkeypatch)
    perms, _ = _perms(ws, answer=None)
    out = json.loads(_DocToolExecutor().execute("enter_worktree", {}, perms))
    assert "error" in out and "nobody can be asked" in out["error"]
    assert not (ws / ".git").exists() and not (ws / ".gitignore").exists()


@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_a_no_is_asked_once_and_leaves_the_folder_as_it_was(tmp_path, monkeypatch, mode):
    ws = _plain(tmp_path, monkeypatch)
    perms, asked = _perms(ws, answer=False, mode=mode)
    ex = _DocToolExecutor()
    first = json.loads(ex.execute("enter_worktree", {}, perms))
    again = json.loads(ex.execute("enter_worktree", {}, perms))
    assert "declined" in first["error"] and "not asked again" in again["error"]
    assert [a[0] for a in asked] == ["local_project_init"]
    assert "LOCAL git repository" in asked[0][1] and "No remote" in asked[0][1]
    assert not (ws / ".git").exists() and not (ws / ".gitignore").exists()


@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_a_yes_makes_one_first_commit_and_the_worktree(tmp_path, monkeypatch, mode):
    ws = _plain(tmp_path, monkeypatch)
    monkeypatch.setattr(LP, "SIZE_LIMIT_BYTES", 4096)
    (ws / "calc").mkdir()
    (ws / "calc" / "job.out").write_text(ORCA_TEXT)            # by name
    (ws / "calc" / "job.gbw").write_bytes(b"\0" * 64)
    (ws / "calc" / "renamed.log").write_text(ORCA_TEXT)        # by content
    (ws / "big.dat").write_bytes(b"0" * 5000)                  # by size
    (ws / ".env").write_text("TOKEN=x\n")
    (ws / "id_ed25519").write_text("key\n")
    (ws / "model.pt").write_bytes(b"w")
    (ws / "synthetic.out.txt").write_text(LP.QM_SYNTHETIC_MARKER + "\n" + ORCA_TEXT)
    perms, asked = _perms(ws, answer=True, mode=mode)
    out = json.loads(_DocToolExecutor().execute("enter_worktree", {}, perms))
    assert out.get("status") == "ok", out
    assert len(asked) == 1, [a[0] for a in asked]    # one question, not two
    assert "recognised as ORCA output" in asked[0][1]
    committed = _commit_tree(ws)
    assert committed == [".gitignore", "notes.md", "shared.py",
                         "synthetic.out.txt"], committed
    assert _git(ws, "branch", "--show-current").strip() == "main"
    assert _git(ws, "remote").strip() == ""                   # no remote
    assert _git(ws, "status", "--porcelain").strip() == ""    # kept out = ignored
    ignore = (ws / ".gitignore").read_text()
    assert LP.GITIGNORE_MARKER in ignore
    assert "/big.dat" in ignore and "/calc/renamed.log" in ignore
    wt_path = Path(out["path"])
    assert (wt_path / "shared.py").exists() and out["branch"].startswith("agent/")
    WT.exit_worktree(WT.WorktreeInfo(Path(out["repo_dir"]), wt_path,
                                     out["branch"], out["base_ref"], 0.0),
                     keep_if_changed=False)


def test_an_existing_gitignore_is_extended_not_replaced(tmp_path, monkeypatch):
    ws = _plain(tmp_path, monkeypatch)
    (ws / ".gitignore").write_text("scratch/\n")
    LP.init_local_project(ws)
    text = (ws / ".gitignore").read_text()
    assert text.startswith("scratch/\n") and LP.GITIGNORE_MARKER in text


def test_a_negation_that_re_includes_orca_output_rolls_back(tmp_path, monkeypatch):
    ws = _plain(tmp_path, monkeypatch)
    (ws / "x.out").write_text("text\n")
    # A negation after DELFIN's block (it is kept as it is, being marked)
    # re-includes *.out; the staged-list check catches it.
    (ws / ".gitignore").write_text(LP.DEFAULT_GITIGNORE + "!*.out\n")
    with pytest.raises(LP.LocalProjectError, match="must stay out"):
        LP.init_local_project(ws)
    assert not (ws / ".git").exists()
    assert (ws / ".gitignore").read_text() == LP.DEFAULT_GITIGNORE + "!*.out\n"


# -- refusals ----------------------------------------------------------------

def test_home_itself_is_refused(tmp_path, monkeypatch):
    ws = _plain(tmp_path, monkeypatch)
    home = Path.home()
    assert "home directory" in LP.refusal(home)
    with pytest.raises(LP.LocalProjectError, match="home directory"):
        LP.init_local_project(home)
    perms, asked = _perms(ws, answer=True)
    out = _DocToolExecutor()._offer_local_project(home, perms, "a worktree")
    assert out and asked == []
    assert not (home / ".git").exists() and not (ws / ".git").exists()


def test_a_folder_inside_another_repository_is_refused(tmp_path, monkeypatch):
    _plain(tmp_path, monkeypatch)
    outer = tmp_path / "outer"
    (outer / "inner").mkdir(parents=True)
    _git(outer, "init", "-q")
    assert "inside the repository" in LP.refusal(outer / "inner")
    with pytest.raises(LP.LocalProjectError):
        LP.init_local_project(outer / "inner")
    assert not (outer / "inner" / ".git").exists()


def test_system_paths_are_refused(tmp_path):
    root = Path(tmp_path.anchor)              # "/" here, a drive elsewhere
    assert LP.refusal(root)
    for name in ("etc", "usr"):
        if (root / name).is_dir():
            assert LP.refusal(root / name)


def test_plan_mode_is_refused_without_a_question(tmp_path, monkeypatch):
    ws = _plain(tmp_path, monkeypatch)
    perms, asked = _perms(ws, answer=True, mode="plan")
    ex = _DocToolExecutor()
    out = ex._offer_local_project(ws, perms, "a worktree")
    assert out and "Plan mode" in json.loads(out)["error"] and asked == []


def test_the_signatures_are_the_license_guards_own():
    """local_project ships in an installed DELFIN; scripts/ does not. The two
    lists must stay one list."""
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "scripts"))
    try:
        import license_guard
    finally:
        sys.path.pop(0)
    assert tuple(LP.QM_OUTPUT_SIGNATURES) == tuple(license_guard.QM_OUTPUT_SIGNATURES)
    assert LP.QM_SYNTHETIC_MARKER == license_guard.QM_SYNTHETIC_MARKER


# -- several agents ------------------------------------------------------------

def _fake_run_subagent(edits: dict):
    """run_subagent standing in for a model: in its worktree (when given
    one) it writes the files ``edits[description]`` names and commits."""
    def _run(*, subagent_type, description, prompt, parent_client,
             parent_perms, isolation="", output_schema=None, **_kw):
        ws = Path(parent_perms.workspace)
        target, info = ws, None
        if isolation == "worktree":
            info = WT.enter_worktree(ws)
            target = info.path
        for rel, text in edits.get(description, {}).items():
            (target / rel).write_text(text)
        summary = {}
        if info is not None:
            _git(target, "add", "-A")
            _git(target, "-c", "user.name=t", "-c", "user.email=t@t",
                 "commit", "-qm", description)
            info = WT.exit_worktree(info, keep_if_changed=True)
            summary = {"branch": info.branch, "had_changes": info.had_changes,
                       "final_path": str(info.final_path or ""),
                       "cleaned_up": info.cleaned_up}
        return SA.SubagentResult(subagent_type=subagent_type,
                                 description=description, final_text="done",
                                 worktree=summary)
    return _run


def _orchestrate(ws, perms, monkeypatch, edits, merge=True):
    monkeypatch.setattr(SA, "run_subagent", _fake_run_subagent(edits))
    perms.orchestration_runner = (
        lambda spec, **kw: SA.run_orchestration(spec, None, perms, **kw))
    spec = {"merge": merge, "stages": [{"name": "work", "parallel": [
        {"subagent_type": "general-purpose", "description": d,
         "prompt": f"edit the files for {d}, a long enough briefing"}
        for d in edits]}]}
    return json.loads(_DocToolExecutor().execute(
        "orchestrate", {"spec": spec}, perms))


def test_two_sub_agents_work_in_their_worktrees_and_are_merged(tmp_path, monkeypatch):
    ws = _plain(tmp_path, monkeypatch)
    perms, asked = _perms(ws, answer=True)
    edits = {"alpha": {"alpha.py": "ALPHA = 1\n"},
             "beta": {"shared.py": "A = 1\nB = 3\n"}}
    out = _orchestrate(ws, perms, monkeypatch, edits)
    assert [a[0] for a in asked] == ["local_project_init"]
    merges = out["merges"]
    assert [(m["description"], m["status"]) for m in merges] == [
        ("alpha", "ok"), ("beta", "ok")], merges
    # The main checkout got both, uncommitted, for review; no sub-agent
    # wrote it directly (its first commit is still the only one).
    assert (ws / "alpha.py").read_text() == "ALPHA = 1\n"
    assert (ws / "shared.py").read_text() == "A = 1\nB = 3\n"
    assert len(_git(ws, "log", "--format=%H").split()) == 1
    assert not [b for b in _git(ws, "branch").split() if b.startswith("agent/")]


def test_a_conflict_is_reported_and_nothing_is_overwritten(tmp_path, monkeypatch):
    ws = _plain(tmp_path, monkeypatch)
    perms, _ = _perms(ws, answer=True)
    edits = {"first": {"shared.py": "A = 10\nB = 2\n"},
             "second": {"shared.py": "A = 20\nB = 2\n"}}
    out = _orchestrate(ws, perms, monkeypatch, edits)
    status = {m["description"]: m["status"] for m in out["merges"]}
    assert status == {"first": "ok", "second": "conflict"}, out["merges"]
    assert (ws / "shared.py").read_text() == "A = 10\nB = 2\n"
    kept = next(m for m in out["merges"] if m["description"] == "second")
    assert (Path(kept["path"]) / "shared.py").read_text() == "A = 20\nB = 2\n"


def test_declined_writers_run_as_before_and_a_merge_is_refused(tmp_path, monkeypatch):
    ws = _plain(tmp_path, monkeypatch)
    perms, _ = _perms(ws, answer=False)
    edits = {"one": {"one.py": "1\n"}, "two": {"two.py": "2\n"}}
    refused = _orchestrate(ws, perms, monkeypatch, edits, merge=True)
    assert "error" in refused and not (ws / "one.py").exists()
    ran = _orchestrate(ws, perms, monkeypatch, edits, merge=False)
    assert "stages" in ran and (ws / "one.py").exists()
    assert not (ws / ".git").exists()


def test_a_writer_whose_worktree_cannot_be_made_does_not_run(tmp_path):
    """run_subagent itself: no fallback into the main checkout for a writer."""

    class _Perms:
        mode = "default"
        workspace = tmp_path
        session_allowed_tools = frozenset()

    res = SA.run_subagent(
        subagent_type="general-purpose", description="w",
        prompt="a self-contained briefing of more than twenty characters",
        parent_client=object(), parent_perms=_Perms(), isolation="worktree")
    assert res.error.startswith("refused: worktree isolation")
    assert "did not run" in res.error


def test_the_always_on_prompt_does_not_name_local_projects():
    """The skill note is loaded for such tasks only; the tool texts the
    model always sees carry the orchestrate merge key and nothing else."""
    from delfin.agent.api_client import _DOC_TOOLS_OPENAI
    text = json.dumps(_DOC_TOOLS_OPENAI)
    assert "local project" not in text.lower()
    assert re.search(r"'merge'\?: true", text)
