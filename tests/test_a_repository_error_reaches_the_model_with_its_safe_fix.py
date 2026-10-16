"""A git error text reaches the model together with its standard fix.

The pack skill ``git-errors`` holds, per situation, the text git prints,
the cause, the fix, what not to do and when to ask. ``repo_error_hints``
maps each text to its section. Pinned here:

* the skill loads through the existing loader, its index entry is one
  line under the listing's per-skill cap, and its sections and the
  trigger table name the same headings;
* each listed error text selects its section, and ordinary output
  (a green test run, a clean status, a log, a traceback, an ORCA line)
  selects none;
* the trigger case is produced by real git and reaches a real ``bash``
  tool result as its section, once per session;
* a user message naming such a text, or asking about checkouts or the
  dashboard's branch, gets a pointer in the turn, and other messages
  get nothing;
* every command the skill tells the agent to run passes the real gate
  where it is meant to run (``boundaries.explain`` and the push/merge
  gate with the git role set), and every command it calls refused is.
"""

from __future__ import annotations

import json
import subprocess
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from conftest import child_env
from delfin import user_settings as US
from delfin.agent import api_client as A
from delfin.agent import boundaries as B
from delfin.agent import job_monitor as jm
from delfin.agent import repo_error_hints as H
from delfin.agent import skills as S
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

pytestmark = pytest.mark.skipif(
    subprocess.run(["git", "--version"], capture_output=True).returncode != 0,
    reason="git is required to build the repositories")

SKILL_FILE = (Path(A.__file__).resolve().parent / "pack" / "skills"
              / "git-errors.md")
SHA = "0123456789abcdef0123456789abcdef01234567"
A_URL = "git@github.com:me/alpha.git"
A_KEY = "github.com/me/alpha"


def _skill():
    sk = S.get_skill(H.SKILL_NAME)
    assert sk is not None, "git-errors does not load"
    return sk


# ---------------------------------------------------------------------------
# the skill and its index line
# ---------------------------------------------------------------------------

def test_the_skill_loads_through_the_existing_loader():
    sk = _skill()
    assert Path(sk.source).resolve() == SKILL_FILE
    assert sk.domains == ("code",)
    assert sk.applies_to("code") and not sk.applies_to("office")
    assert 0 < len(sk.description) <= A._SKILL_LISTING_DESC_CHARS
    listing = A._skill_listing(S.discover_skills(domain="code"))
    entry = f"{sk.name} — {sk.description}"
    assert entry in listing and "\n" not in entry


def test_every_section_has_a_trigger_and_every_trigger_a_section():
    sk = _skill()
    assert H.headings(sk.body.splitlines()) == [h for h, _ in H.TRIGGERS]
    for heading, text in H.sections().items():
        assert text.startswith(f"## {heading}")
        assert "Never:" in text, f"{heading}: says nothing about what not to do"


def test_the_delfin_section_starts_from_ops_status():
    text = H.sections()["DELFIN's own checkouts and servers"]
    assert "`delfin-agent ops status`" in text
    first_fix = [ln for ln in text.splitlines() if ln.startswith("Fix:")][0]
    assert "delfin-agent ops status" in first_fix


# ---------------------------------------------------------------------------
# which text selects which section
# ---------------------------------------------------------------------------

ERRORS = [
    ("fatal: 'main' is already used by worktree at '/srv/x/dashboard'",
     "Branch held by another worktree"),
    ("fatal: 'main' is already checked out at '/srv/x/dashboard'",
     "Branch held by another worktree"),
    ("Note: switching to 'origin/main'.\n\nYou are in 'detached HEAD' state. "
     "You can look around, make experimental\n",
     "Detached HEAD"),
    ("HEAD detached at 1a2b3c4\nnothing to commit, working tree clean",
     "Detached HEAD"),
    ("fatal: Unable to create '/srv/x/.git/index.lock': File exists.\n\n"
     "Another git process seems to be running in this repository",
     "index.lock exists"),
    (" ! [rejected]        main -> main (fetch first)\n"
     "error: failed to push some refs to 'github.com:me/x.git'",
     "Push rejected or branches diverged"),
    (" ! [rejected]        feat -> feat (non-fast-forward)",
     "Push rejected or branches diverged"),
    ("hint: Updates were rejected because the tip of your current branch is "
     "behind", "Push rejected or branches diverged"),
    ("Your branch and 'origin/main' have diverged,\nand have 2 and 3 "
     "different commits each", "Push rejected or branches diverged"),
    ("fatal: Need to specify how to reconcile divergent branches.",
     "Push rejected or branches diverged"),
    ("error: The following untracked working tree files would be "
     "overwritten by checkout:\n\tdata.csv",
     "Untracked files would be overwritten"),
    ("error: Your local changes to the following files would be overwritten "
     "by merge:\n\tf.py\nPlease commit your changes or stash them before you "
     "merge.", "Local changes would be overwritten"),
    ("Auto-merging f.py\nCONFLICT (content): Merge conflict in f.py\n"
     "Automatic merge failed; fix conflicts and then commit the result.",
     "Merge conflict"),
    ("def f():\n<<<<<<< HEAD\n    return 1\n=======\n    return 2\n"
     ">>>>>>> feature\n", "Merge conflict"),
    ("fatal: not a git repository (or any of the parent directories): .git",
     "Not a git repository"),
    ("fatal: not a git repository: /srv/x/.git/worktrees/session-7",
     "Not a git repository"),
    ("! [remote rejected] main (shallow update not allowed)",
     "Shallow clone"),
    ("fatal: ambiguous argument 'v1.0..HEAD': unknown revision or path not "
     "in the working tree.", "Shallow clone"),
    ("1a2b3c4 (grafted, HEAD -> main) last commit", "Shallow clone"),
    ("version https://git-lfs.github.com/spec/v1\noid sha256:4d7a\n"
     "size 12345\n", "LFS pointer files"),
    ("git-lfs filter-process: git-lfs: command not found",
     "LFS pointer files"),
    ("Encountered 2 file(s) that should have been pointers, but weren't:",
     "LFS pointer files"),
    ("git@github.com: Permission denied (publickey).\nfatal: Could not read "
     "from remote repository.", "Permission denied (publickey)"),
    ("This conversation is already open in another session — open it there",
     "DELFIN's own checkouts and servers"),
    ("[Voila] ERROR | OSError: [Errno 98] Address already in use",
     "DELFIN's own checkouts and servers"),
    ("remote: error: GH006: Protected branch update failed for "
     "refs/heads/main.", "Contributing a change"),
    (" ! [remote rejected] main -> main (protected branch hook declined)",
     "Contributing a change"),
]


@pytest.mark.parametrize("text,heading", ERRORS)
def test_each_error_text_selects_its_section(text, heading):
    assert heading in H.match(text), H.match(text)
    hint = H.tool_hint(text)
    assert f"## {heading}" in hint and "skill(name='git-errors')" in hint


def test_every_section_is_reached_by_some_listed_text():
    reached = {h for text, _ in ERRORS for h in H.match(text)}
    assert reached == {h for h, _ in H.TRIGGERS}


UNRELATED = [
    "",
    "On branch main\nYour branch is up to date with 'origin/main'.\n\n"
    "nothing to commit, working tree clean",
    "a1b2c3d fix the energy table\n9f8e7d6 add the CSV export\n",
    "============ 412 passed, 3 skipped in 12.31s ============",
    "Traceback (most recent call last):\n  File \"x.py\", line 3\n"
    "KeyError: 'energy'",
    "FINAL SINGLE POINT ENERGY      -1234.567890123\n****ORCA TERMINATED "
    "NORMALLY****",
    "Switched to branch 'feature'\n",
    "To github.com:me/x.git\n   1a2b3c4..5d6e7f8  feature -> feature",
    "Already up to date.",
    "drwxr-xr-x 2 u u 4096 Oct 14 10:00 worktrees\n-rw-r--r-- 1 u u 3 x",
    "# Conflict of interest statement\nThe authors declare no conflict.",
]


@pytest.mark.parametrize("text", UNRELATED)
def test_ordinary_output_selects_nothing(text):
    assert H.match(text) == []
    assert H.tool_hint(text) == ""


def test_a_section_is_delivered_once_per_session_and_two_at_most():
    shown: set = set()
    many = "\n".join(t for t, _ in ERRORS)
    first = H.tool_hint(many, shown)
    assert first.count("\n## ") + first.startswith("## ") == 2
    assert len(shown) == 2
    again = H.tool_hint(ERRORS[0][0], shown)
    assert again == ""


# ---------------------------------------------------------------------------
# the trigger case, produced by real git, through the real bash tool
# ---------------------------------------------------------------------------

def _git(*args, cwd, env):
    return subprocess.run(["git", *args], cwd=str(cwd), capture_output=True,
                          text=True, env=env)


@pytest.fixture
def held_main(tmp_path):
    """``ws/main_checkout`` holds main; ``ws/dashboard`` is a second
    worktree on ``dash`` one commit behind."""
    env = child_env(tmp_path)
    env["GIT_CONFIG_GLOBAL"] = str(tmp_path / "gitconfig")
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    (tmp_path / "gitconfig").write_text("[user]\n\tname = t\n\temail = t@t\n")
    ws = tmp_path / "ws"
    main = ws / "main_checkout"
    main.mkdir(parents=True)
    assert _git("init", "-q", "-b", "main", str(main), cwd=ws, env=env
                ).returncode == 0
    (main / "f.txt").write_text("1\n")
    _git("add", "f.txt", cwd=main, env=env)
    _git("commit", "-qm", "one", cwd=main, env=env)
    _git("branch", "dash", cwd=main, env=env)
    (main / "f.txt").write_text("2\n")
    _git("commit", "-qam", "two", cwd=main, env=env)
    r = _git("worktree", "add", "-q", str(ws / "dashboard"), "dash",
             cwd=main, env=env)
    assert r.returncode == 0, r.stderr
    return {"ws": ws, "main": main, "dash": ws / "dashboard", "env": env}


def test_real_git_prints_the_text_the_first_section_names(held_main):
    r = _git("switch", "main", cwd=held_main["dash"], env=held_main["env"])
    assert r.returncode != 0
    assert H.match(r.stderr) == ["Branch held by another worktree"], r.stderr


def test_the_bash_tool_result_carries_the_section_once(held_main,
                                                       monkeypatch):
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", held_main["env"]["GIT_CONFIG_GLOBAL"])
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    perms = KitToolPermissions(workspace=str(held_main["ws"]),
                               mode="acceptEdits")
    ex = _DocToolExecutor()

    def run(cmd):
        return json.loads(ex.execute(
            "bash", {"command": cmd, "cwd": "dashboard"}, perms))

    # A held branch no longer fails: the executor checks out the same
    # commit detached (#238) and says so. git's own refusal text still
    # names the held-branch section, once.
    detached = run("git switch main")
    assert detached.get("exit_code") == 0, detached
    if "skill_hint" in detached:
        assert "## Branch held by another worktree" in detached["skill_hint"]
    assert "skill_hint" not in run("git switch main")
    # A real git error still carries its section once: a stale index.lock.
    gitdir = Path(_git("rev-parse", "--absolute-git-dir",
                       cwd=held_main["dash"], env=held_main["env"]
                       ).stdout.strip())
    (gitdir / "index.lock").write_text("")
    first = run("git switch dash")
    assert first.get("exit_code") not in (0, None), first
    assert "index.lock" in first["stderr"]
    hint = first.get("skill_hint", "")
    assert "## index.lock exists" in hint
    assert "Detached HEAD" not in hint
    assert "skill_hint" not in run("git switch dash")
    # The checkout is detached now; git status says so, and that section
    # too is carried once.
    status = run("git status")
    assert "## Detached HEAD" in status.get("skill_hint", ""), status
    assert "skill_hint" not in run("git status")


# ---------------------------------------------------------------------------
# a user message
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("text,heading", [
    ("Ich bekomme fatal: 'main' is already used by worktree at '/x/dash' "
     "- was tun?", "Branch held by another worktree"),
    ("Welcher Branch ist im Dashboard ausgecheckt?",
     "DELFIN's own checkouts and servers"),
    ("which checkout does the dashboard run from?",
     "DELFIN's own checkouts and servers"),
    ("Räum die alten worktrees auf", "DELFIN's own checkouts and servers"),
    ("Open a pull request for this", "Contributing a change"),
])
def test_a_user_message_calls_for_its_section(text, heading):
    assert heading in H.user_topics(text)
    assert "skill(name='git-errors')" in H.user_hint(text)


@pytest.mark.parametrize("text", [
    "Optimise the geometry of the Cr complex with r2SCAN-3c.",
    "Make the energy table a CSV, please.",
    "Wie viele Rechnungen sind noch nicht fertig?",
    "Report the transport properties.",
])
def test_other_user_messages_get_nothing(text):
    assert H.user_topics(text) == []
    assert H.user_hint(text) == ""


@pytest.fixture
def world(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("DELFIN_HOME", raising=False)
    ws = tmp_path / "ws"
    ws.mkdir()
    pack = tmp_path / "packs"
    shared = pack / "pack" / "shared"
    shared.mkdir(parents=True)
    agents = pack / "pack" / "agents"
    agents.mkdir()
    for name in ("delfin_context.md", "work_cycle_rules.md",
                 "goal_decomposition_rules.md",
                 "universal_input_template.md", "minimal_final_verdict.md"):
        (shared / name).write_text("# x")
    for name in ("session_manager", "builder_agent", "test_agent"):
        (agents / f"{name}.md").write_text(f"# {name}")
    modes = pack / "pack_lite" / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (pack / "pack_lite" / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
              - builder_agent
              - test_agent
    """))
    return {"ws": ws, "pack": pack}


def _engine(world):
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine
    seen: list = []

    def fake_stream(system, messages, max_tokens=4096, session_id="",
                    thinking_budget=0):
        seen.append(str(system) + "\n" + json.dumps(messages, default=str))
        yield StreamEvent(type="session_init", text="repo-err")
        yield StreamEvent(type="message_start", input_tokens=10)
        yield StreamEvent(type="text_delta", text="done")
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=fake_stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=world["pack"], backend="cli",
                             mode="quick", pack_dir=world["pack"])
    client._permissions = KitToolPermissions(workspace=str(world["ws"]),
                                             mode="acceptEdits")
    return engine, seen


def test_the_turn_carries_the_pointer_only_when_called_for(world):
    # The steering list of the interactive role, as the L2 correction test
    # reads it: the blocks the solo agent's turn tail is built from.
    engine, _seen = _engine(world)
    engine.stream_response(
        "git says: fatal: 'main' is already used by worktree at '/x/dash'")
    blocks = dict(engine._steering_blocks("solo_agent"))
    assert "[git-errors]" in blocks.get("repo_errors", "")
    assert "Branch held by another worktree" in blocks["repo_errors"]
    engine.stream_response("Optimise the geometry of the Cr complex.")
    assert "repo_errors" not in dict(engine._steering_blocks("solo_agent"))


def test_a_session_without_skills_gets_no_pointer(world):
    engine, _seen = _engine(world)
    engine.kit_permissions.skip_skill_discovery = True
    engine.stream_response(
        "git says: fatal: 'main' is already used by worktree at '/x/dash'")
    assert "repo_errors" not in dict(engine._steering_blocks("solo_agent"))


# ---------------------------------------------------------------------------
# the gate: what the skill says to run runs; what it calls refused is
# ---------------------------------------------------------------------------

#: (as written in the skill, as run in the fixture repository)
AGENT_RUNS = [
    ("`git worktree list`", "git worktree list"),
    ("`git fetch origin`", "git fetch origin"),
    ("`git switch --detach origin/<branch>`", "git switch --detach origin/main"),
    ("`git switch <branch>`", "git switch main"),
    ("`git branch <name>`", "git branch keep-work"),
    ("`git reflog`", "git reflog"),
    ("`pgrep -a git`", "pgrep -a git"),
    ("`git pull --rebase`", "git pull --rebase"),
    ("`git merge origin/<branch>`", "git merge origin/main"),
    ("`git diff --name-only --diff-filter=U`",
     "git diff --name-only --diff-filter=U"),
    ("`git diff --check`", "git diff --check"),
    ("`git add <file>`", "git add f.txt"),
    ("`git merge --abort`", "git merge --abort"),
    ("`git rebase --abort`", "git rebase --abort"),
    ("`pwd`", "pwd"),
    ("`cat .git`", "cat .git"),
    ("`git rev-parse --is-shallow-repository`",
     "git rev-parse --is-shallow-repository"),
    ("`git fetch --unshallow`", "git fetch --unshallow"),
    ("`git lfs pull`", "git lfs pull"),
    ("`git remote -v`", "git remote -v"),
    ("`delfin-agent ops status`", "delfin-agent ops status"),
    ("`git switch -c <user>/<topic>`", "git switch -c alice/topic"),
    ("`git add <paths>`", "git add f.txt"),
    ("`git commit -m \"<what and why>\"`", "git commit -m \"step one\""),
    ("`gh pr create --base main --fill`", "gh pr create --base main --fill"),
    ("`gh pr view <n> --json headRefOid,statusCheckRollup`",
     "gh pr view 12 --json headRefOid,statusCheckRollup"),
]

#: Commands the skill calls refused, or leaves to the user.
REFUSED = [
    ("(refused)", "git stash"),
    ("`git push --force` (refused)", "git push --force origin alice/topic"),
    ("`git reset --hard`", "git reset --hard origin/main"),
    ("`rm <path>/index.lock`", "rm .git/index.lock"),
]


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "idx.json")
    cfg = tmp_path / "gitconfig"
    cfg.write_text("[user]\n\tname = t\n\temail = t@t\n")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(cfg))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    monkeypatch.delenv("GH_REPO", raising=False)
    env = child_env(tmp_path)
    env["GIT_CONFIG_GLOBAL"] = str(cfg)
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    root = tmp_path / "alpha"
    root.mkdir()
    _git("init", "-q", "-b", "main", str(root), cwd=tmp_path, env=env)
    (root / "f.txt").write_text("x\n")
    _git("add", "f.txt", cwd=root, env=env)
    _git("commit", "-qm", "init", cwd=root, env=env)
    _git("remote", "add", "origin", A_URL, cwd=root, env=env)
    return root


def _allowed(verdict: dict) -> bool:
    return (verdict["verdict"] == "allow"
            or (verdict["verdict"] == "ask"
                and verdict["if_approved"] == "allow"))


@pytest.mark.parametrize("written,command", AGENT_RUNS,
                         ids=[c for _, c in AGENT_RUNS])
def test_what_the_skill_says_to_run_passes_the_gate(repo, written, command):
    assert written in SKILL_FILE.read_text(encoding="utf-8"), written
    v = B.explain("bash", {"command": command}, cwd=str(repo),
                  permission_mode="default", asker="terminal")
    assert _allowed(v), (command, v.get("verdict"), v.get("reason"))


@pytest.mark.parametrize("written,command", REFUSED,
                         ids=[c for _, c in REFUSED])
def test_what_the_skill_calls_refused_is_refused(repo, written, command):
    assert written in SKILL_FILE.read_text(encoding="utf-8"), written
    v = B.explain("bash", {"command": command}, cwd=str(repo),
                  permission_mode="bypassPermissions", asker="terminal")
    assert v["verdict"] == "refuse", (command, v)


def _gate_perms(ws, *, grant=True):
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    if grant:
        A._grant_push_from(perms, "commit and push it", new_request=True)
    return perms


def _gate(perms, cmd):
    return _DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd}, perms)


def test_a_contributor_pushes_the_branch_the_skill_names(repo):
    US.save_settings({"agent": {"git_role": "contributor"}})
    text = SKILL_FILE.read_text(encoding="utf-8")
    assert "`git push -u origin <branch>`" in text
    subprocess.run(["git", "switch", "-q", "-c", "alice/topic"],
                   cwd=str(repo), check=True, capture_output=True)
    assert _gate(_gate_perms(repo), "git push -u origin alice/topic") is None
    err = _gate(_gate_perms(repo), "git push origin main")
    assert err and err.startswith("blocked:"), err


def test_a_maintainer_merges_as_the_skill_says_and_only_so(repo, monkeypatch):
    US.save_settings({"agent": {"git_roles": {A_KEY: "maintainer"}}})
    text = SKILL_FILE.read_text(encoding="utf-8")
    assert ("`gh pr merge <n> --merge --match-head-commit <headRefOid>`"
            in text)
    monkeypatch.setattr(A, "_gh_pr_view", lambda where, r, pr: {
        "headRefOid": SHA, "statusCheckRollup": [
            {"__typename": "CheckRun", "name": "tests",
             "status": "COMPLETED", "conclusion": "SUCCESS"}]})
    cmd = f"gh pr merge 12 --merge --match-head-commit {SHA}"
    assert _gate(_gate_perms(repo), cmd) is None
    for flag in ("--admin", "--auto"):
        assert f"`{flag}`" in text
        err = _gate(_gate_perms(repo), f"{cmd} {flag}")
        assert err and flag in err, err
    err = _gate(_gate_perms(repo), "gh pr merge 12 --merge")
    assert err and "--match-head-commit" in err, err


def test_a_contributor_merge_is_left_to_the_user(repo, monkeypatch):
    US.save_settings({"agent": {"git_role": "contributor"}})
    monkeypatch.setattr(A, "_gh_pr_view", lambda where, r, pr: {
        "headRefOid": SHA, "statusCheckRollup": [
            {"__typename": "CheckRun", "name": "tests",
             "status": "COMPLETED", "conclusion": "SUCCESS"}]})
    err = _gate(_gate_perms(repo),
                f"gh pr merge 12 --merge --match-head-commit {SHA}")
    assert err and err.startswith("blocked:"), err
