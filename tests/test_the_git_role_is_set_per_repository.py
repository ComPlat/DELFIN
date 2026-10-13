"""The git role is resolved for the repository a command acts on.

``agent.git_role`` was one value for the installation: ``maintainer``
unlocked a push onto the default branch and ``gh pr merge`` in every
repository at once, ``contributor`` refused both everywhere. These tests
pin the per-repository table ``agent.git_roles`` through the real bash
gate (``_run_permission_gate``), with local repositories whose remotes are
fake URLs -- nothing is pushed and nothing reaches a network.

Also pinned, because the role decides only one question: the one-request
grant, the unattended refusal, the force refusal and the rejected-push
report are unchanged for a maintainer; a typed request for the
pull-request route holds even where the role would allow a push; and a
maintainer's merge is pinned to a head whose checks are green.
"""

from __future__ import annotations

import json
import os
import subprocess

import pytest

from conftest import child_env
from delfin import user_settings as US
from delfin.agent import api_client as A
from delfin.agent import harness_rows as HR
from delfin.agent import job_monitor as jm
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

pytestmark = pytest.mark.skipif(
    subprocess.run(["git", "--version"], capture_output=True).returncode != 0,
    reason="git is required to build the repositories")

A_URL = "git@github.com:me/alpha.git"
B_URL = "https://github.com/me/beta.git"
A_KEY = "github.com/me/alpha"
B_KEY = "github.com/me/beta"
SHA = "0123456789abcdef0123456789abcdef01234567"


def _git(*args, cwd, env):
    subprocess.run(["git", *args], cwd=str(cwd), check=True,
                   capture_output=True, env=env)


@pytest.fixture(autouse=True)
def _isolated(tmp_path, monkeypatch):
    monkeypatch.setattr(jm, "_AGENT_WATCH_INDEX_PATH", tmp_path / "index.json")
    cfg = tmp_path / "gitconfig"
    cfg.write_text("[user]\n\tname = t\n\temail = t@t\n")
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", str(cfg))
    monkeypatch.setenv("GIT_CONFIG_NOSYSTEM", "1")
    monkeypatch.delenv("GH_REPO", raising=False)


@pytest.fixture
def repos(tmp_path):
    """Two repositories side by side, each on ``main`` with one commit."""
    env = child_env(tmp_path)
    env["GIT_CONFIG_GLOBAL"] = str(tmp_path / "gitconfig")
    env["GIT_CONFIG_NOSYSTEM"] = "1"
    out = {}
    for name, url in (("alpha", A_URL), ("beta", B_URL)):
        root = tmp_path / name
        root.mkdir()
        _git("init", "-q", str(root), cwd=tmp_path, env=env)
        _git("checkout", "-q", "-b", "main", cwd=root, env=env)
        (root / "f.txt").write_text("x\n")
        _git("add", "f.txt", cwd=root, env=env)
        _git("commit", "-qm", "init", cwd=root, env=env)
        _git("remote", "add", "origin", url, cwd=root, env=env)
        out[name] = root
    out["env"] = env
    return out


def _settings(**agent):
    US.save_settings({"agent": agent})


def _perms(ws, *, mode="bypassPermissions", grant=True, answer=None):
    perms = KitToolPermissions(workspace=ws, mode=mode)
    if grant:
        A._grant_push_from(perms, "commit and push it", new_request=True)
    if answer is not None:
        perms.asked = []

        def _cb(*a, **k):
            perms.asked.append(a)
            return answer
        perms.confirm_callback = _cb
    return perms


def _gate(perms, cmd, **extra):
    return _DocToolExecutor()._run_permission_gate(
        "bash", {"command": cmd, **extra}, perms)


def _refused_as_contributor(err):
    return bool(err) and "pull request" in err and err.startswith("blocked:")


# ---------------------------------------------------------------------------
# Keys
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("url", [
    "git@github.com:me/alpha.git", "git@GitHub.com:me/alpha",
    "ssh://git@github.com/me/alpha.git", "ssh://git@github.com:22/me/alpha",
    "https://github.com/me/alpha.git", "https://user@github.com/me/alpha/",
    "HTTPS://GITHUB.COM/me/alpha",
])
def test_every_spelling_of_one_remote_is_one_key(url):
    from delfin.agent import git_role as GR
    assert GR.normalise_remote(url) == A_KEY


def test_a_key_keeps_the_path_case_and_rejects_nonsense():
    from delfin.agent import git_role as GR
    assert GR.normalise_remote("git@github.com:Me/Alpha.git") == "github.com/Me/Alpha"
    assert GR.normalise_remote("https://github.com/") is None
    assert GR.normalise_remote("https://github.com/me/../x") is None
    assert GR.normalise_remote("ftp://x") is None
    assert GR.normalise_remote("") is None


# ---------------------------------------------------------------------------
# One repository maintainer, the other not
# ---------------------------------------------------------------------------

def test_own_repository_maintainer_pushes_main_with_the_grant(repos):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    assert _gate(_perms(repos["alpha"]), "git push origin main") is None


def test_own_repository_maintainer_without_the_grant_is_refused(repos):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"], grant=False), "git push origin main")
    assert err and "has not asked for a push" in err, err


def test_the_other_repository_stays_contributor(repos):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["beta"]), "git push origin main")
    assert _refused_as_contributor(err), err
    assert B_KEY in err, err


def test_both_in_one_session(repos):
    """One grant each, one session: the role follows the directory."""
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    perms = _perms(repos["alpha"])
    assert _gate(perms, "git push origin main") is None
    err = _gate(perms, "git push origin main", cwd=str(repos["beta"]))
    assert _refused_as_contributor(err), err


@pytest.mark.parametrize("url", [B_URL, "git@github.com:me/beta",
                                 "ssh://git@github.com/me/beta.git"])
def test_a_push_by_url_resolves_by_that_url(repos, url):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]), f"git push {url} main")
    assert _refused_as_contributor(err), err


def test_a_push_by_url_to_the_maintainer_repository_is_allowed(repos):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    assert _gate(_perms(repos["beta"]),
                 "git push https://github.com/me/alpha main") is None


def test_git_dash_c_resolves_the_named_repository(repos):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]), "git -C ../beta push origin main")
    assert _refused_as_contributor(err), err
    assert _gate(_perms(repos["beta"]),
                 "git -C ../alpha push origin main") is None


def test_a_leading_cd_resolves_the_repository_it_enters(repos):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]), "cd ../beta && git push origin main")
    assert _refused_as_contributor(err), err


def test_a_global_maintainer_still_works(repos):
    _settings(git_role="maintainer")
    assert _gate(_perms(repos["alpha"]), "git push origin main") is None
    assert _gate(_perms(repos["beta"]), "git push origin main") is None


def test_a_repository_contributor_wins_over_a_global_maintainer(repos):
    _settings(git_role="maintainer", git_roles={B_KEY: "contributor"})
    err = _gate(_perms(repos["beta"]), "git push origin main")
    assert _refused_as_contributor(err), err
    assert _gate(_perms(repos["alpha"]), "git push origin main") is None


def test_an_entry_written_as_a_url_counts(repos):
    _settings(git_role="contributor", git_roles={A_URL: "maintainer"})
    assert _gate(_perms(repos["alpha"]), "git push origin main") is None


@pytest.mark.parametrize("table", ["maintainer", ["github.com/me/alpha"], 7])
def test_a_malformed_table_is_contributor(repos, table):
    _settings(git_role="maintainer", git_roles=table)
    err = _gate(_perms(repos["alpha"]), "git push origin main")
    assert _refused_as_contributor(err), err


def test_an_entry_that_is_no_role_is_contributor(repos):
    _settings(git_role="maintainer", git_roles={A_KEY: "owner"})
    err = _gate(_perms(repos["alpha"]), "git push origin main")
    assert _refused_as_contributor(err), err


def test_a_branch_push_is_untouched_by_the_role(repos):
    _settings(git_role="contributor")
    assert _gate(_perms(repos["beta"]), "git push -u origin me/topic") is None


def test_a_local_repository_is_keyed_by_its_path(repos, tmp_path):
    env = repos["env"]
    local = tmp_path / "local"
    local.mkdir()
    _git("init", "-q", str(local), cwd=tmp_path, env=env)
    _git("checkout", "-q", "-b", "main", cwd=local, env=env)
    (local / "f.txt").write_text("x\n")
    _git("add", "f.txt", cwd=local, env=env)
    _git("commit", "-qm", "init", cwd=local, env=env)
    _settings(git_role="contributor")
    assert _refused_as_contributor(_gate(_perms(local), "git push"))
    _settings(git_role="contributor",
              git_roles={str(local.resolve()): "maintainer"})
    from delfin.agent import git_role as GR
    assert GR.repo_key(local) == str(local.resolve())
    assert _gate(_perms(local), "git push") is None


# ---------------------------------------------------------------------------
# A command whose repository cannot be told is contributor
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("cmd", [
    "cd . && cd ../beta && git push origin main",
    "git -c remote.origin.pushurl=https://github.com/me/beta push origin main",
    "git --work-tree=../beta push origin main",
    "GIT_WORK_TREE=../beta git push origin main",
    "sh -c 'cd ../beta && git push origin main'",
    "git push origin main && git -C ../beta push origin main",
    "for d in . ../beta; do git -C $d push origin main; done",
    "pushd ../beta && git push origin main",
])
def test_an_unresolvable_shape_is_contributor(repos, cmd):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]), cmd)
    assert _refused_as_contributor(err), (cmd, err)


def test_every_push_url_of_the_remote_counts(repos):
    env = repos["env"]
    _git("remote", "set-url", "--add", "--push", "origin", A_URL,
         cwd=repos["alpha"], env=env)
    _git("remote", "set-url", "--add", "--push", "origin", B_URL,
         cwd=repos["alpha"], env=env)
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]), "git push origin main")
    assert _refused_as_contributor(err), err


def test_a_push_insteadof_rewrite_is_followed(repos):
    env = repos["env"]
    _git("config", "url.https://github.com/me/beta.pushInsteadOf",
         "https://github.com/me/alpha", cwd=repos["alpha"], env=env)
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]),
                "git push https://github.com/me/alpha main")
    assert _refused_as_contributor(err), err


def test_a_second_remote_counts_when_none_is_named(repos):
    _git("remote", "add", "upstream", B_URL, cwd=repos["alpha"],
         env=repos["env"])
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]), "git push")
    assert _refused_as_contributor(err), err
    assert _gate(_perms(repos["alpha"]), "git push origin main") is None


# ---------------------------------------------------------------------------
# What the role does not change
# ---------------------------------------------------------------------------

def test_a_maintainer_force_push_is_refused(repos):
    _settings(git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]), "git push --force origin main")
    assert err, "a force push went through for a maintainer"


def test_a_maintainer_attended_without_a_grant_is_asked(repos):
    _settings(git_roles={A_KEY: "maintainer"})
    perms = _perms(repos["alpha"], mode="default", grant=False, answer=False)
    err = _gate(perms, "git push origin main")
    assert perms.asked, "the user was not asked"
    assert err and "did not approve" in err, err


def test_a_rejected_maintainer_push_is_reported_and_spends_nothing(repos):
    _settings(git_roles={A_KEY: "maintainer"})
    perms = _perms(repos["alpha"])
    note = A._after_push(
        {"command": "git push origin main"},
        json.dumps({"exit_code": 0, "stdout": "", "stderr": (
            "To github.com:me/alpha.git\n"
            " ! [rejected]        main -> main (fetch first)\n"
            "error: failed to push some refs\n")}),
        perms)
    assert note == ""
    assert perms.push_grants.get("push") == 1


# ---------------------------------------------------------------------------
# gh pr merge
# ---------------------------------------------------------------------------

def _green(head=SHA, checks=None):
    rollup = checks if checks is not None else [
        {"__typename": "CheckRun", "name": "tests", "status": "COMPLETED",
         "conclusion": "SUCCESS"},
        {"__typename": "StatusContext", "context": "lint", "state": "SUCCESS"},
    ]
    return {"headRefOid": head, "statusCheckRollup": rollup}


@pytest.fixture
def pr_view(monkeypatch):
    seen = []

    def _set(payload):
        def _view(where, repo, pr):
            seen.append((str(where), repo, pr))
            return payload
        monkeypatch.setattr(A, "_gh_pr_view", _view)
    _set.seen = seen
    return _set


def _merge(perms, tail=""):
    return _gate(perms, f"gh pr merge 12 --squash --match-head-commit {SHA}{tail}")


def test_a_maintainer_merges_a_pinned_green_head(repos, pr_view):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green())
    assert _merge(_perms(repos["alpha"])) is None
    assert pr_view.seen and pr_view.seen[0][1] == A_KEY


def test_a_merge_in_the_contributor_repository_is_refused(repos, pr_view):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green())
    err = _merge(_perms(repos["beta"]))
    assert err and "maintainer's decision" in err, err


def test_a_merge_naming_another_repository_resolves_by_it(repos, pr_view):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green())
    for spec in ("-R me/beta", "--repo me/beta", "--repo=github.com/me/beta"):
        err = _merge(_perms(repos["alpha"]), " " + spec)
        assert err and "maintainer's decision" in err, (spec, err)
    err = _gate(_perms(repos["alpha"]),
                f"gh pr merge https://github.com/me/beta/pull/3 "
                f"--match-head-commit {SHA}")
    assert err and "maintainer's decision" in err, err


def test_a_merge_without_the_grant_is_refused(repos, pr_view):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green())
    err = _merge(_perms(repos["alpha"], grant=False))
    assert err and "has not asked" in err, err


def test_a_merge_without_the_head_pin_is_refused(repos, pr_view):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green())
    err = _gate(_perms(repos["alpha"]), "gh pr merge 12 --squash")
    assert err and "--match-head-commit" in err, err


@pytest.mark.parametrize("flag", ["--admin", "--auto"])
def test_a_merge_past_the_checks_is_refused(repos, pr_view, flag):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green())
    err = _merge(_perms(repos["alpha"]), " " + flag)
    assert err and flag in err, err


def test_a_merge_of_a_moved_head_is_refused(repos, pr_view):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green(head="f" * 40))
    err = _merge(_perms(repos["alpha"]))
    assert err and "head" in err and "f" * 12 in err, err


@pytest.mark.parametrize("check", [
    {"__typename": "CheckRun", "name": "tests", "status": "IN_PROGRESS",
     "conclusion": ""},
    {"__typename": "CheckRun", "name": "tests", "status": "COMPLETED",
     "conclusion": "FAILURE"},
    {"__typename": "StatusContext", "context": "tests", "state": "PENDING"},
])
def test_a_merge_with_a_check_not_green_is_refused(repos, pr_view, check):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green(checks=[check]))
    err = _merge(_perms(repos["alpha"]))
    assert err and "tests" in err and "green" in err, err


def test_a_merge_whose_checks_cannot_be_read_is_refused(repos, pr_view):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(None)
    err = _merge(_perms(repos["alpha"]))
    assert err and err.startswith("blocked:"), err


def test_the_checks_are_read_through_gh(repos, tmp_path, monkeypatch):
    """The real reader, against a gh that prints a fixed answer."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    log = tmp_path / "gh-args"
    gh = bindir / "gh"
    gh.write_text("#!/bin/sh\n"
                  f"echo \"$@\" >> '{log}'\n"
                  f"cat <<'EOF'\n{json.dumps(_green())}\nEOF\n")
    gh.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ['PATH']}")
    _settings(git_roles={A_KEY: "maintainer"})
    assert _merge(_perms(repos["alpha"])) is None
    args = log.read_text()
    assert "pr view 12" in args and f"-R {A_KEY}" in args, args
    assert "headRefOid" in args and "statusCheckRollup" in args, args


# ---------------------------------------------------------------------------
# The user asked for the pull-request route
# ---------------------------------------------------------------------------

def _typed(perms, text):
    A._grant_push_for_request(perms, [{"role": "user", "content": text}])


def _harness(perms, text):
    A._grant_push_for_request(perms, [HR.harness_row(HR.HarnessText(text))])


@pytest.mark.parametrize("phrase", [
    "Mach das bitte über PR und push es",
    "push das per Pull Request",
    "bitte nicht direkt auf main, push es",
    "kein direkter Push, aber push den Branch",
    "push it via a PR",
    "push it, but don't push to main",
    "no direct push to main please, push the branch",
])
def test_a_typed_request_for_the_pr_route_holds_for_a_maintainer(
        repos, pr_view, phrase):
    _settings(git_role="maintainer")
    pr_view(_green())
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, phrase)
    assert perms.push_grants.get("push"), phrase
    err = _gate(perms, "git push origin main")
    assert err and "pull request" in err and "asked" in err, (phrase, err)
    assert "git switch -c" in err and "link" in err, err
    err = _merge(perms)
    assert err and "pull request" in err, (phrase, err)
    assert _gate(perms, "git push -u origin me/topic") is None


def test_the_route_holds_for_the_rest_of_the_session(repos):
    _settings(git_role="maintainer")
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, "über PR bitte")
    _typed(perms, "jetzt push es")
    assert perms.push_grants.get("push")
    err = _gate(perms, "git push origin main")
    assert err and "pull request" in err, err


@pytest.mark.parametrize("lift", [
    "push it directly to main", "push direkt auf main",
    "jetzt direkt auf main pushen"])
def test_a_later_typed_request_for_the_direct_push_lifts_it(repos, lift):
    _settings(git_role="maintainer")
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, "über PR bitte")
    _typed(perms, lift)
    assert _gate(perms, "git push origin main") is None, lift


def test_harness_text_neither_sets_nor_lifts_the_route(repos):
    _settings(git_role="maintainer")
    perms = _perms(repos["alpha"])
    _harness(perms, "[watch] über PR, kein direkter Push")
    perms.push_grants["push"] = 1
    assert _gate(perms, "git push origin main") is None
    _typed(perms, "push it via a PR")
    _harness(perms, "push it directly to main")
    perms.push_grants["push"] = 1
    err = _gate(perms, "git push origin main")
    assert err and "pull request" in err, err


def test_the_route_is_held_for_the_repository_it_was_asked_in(repos):
    _settings(git_role="maintainer")
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, "push it via a PR")
    err = _gate(perms, "git push origin main")
    assert err and "pull request" in err, err
    assert _gate(perms, "git -C ../beta push origin main") is None


def test_a_sub_agent_prompt_does_not_lift_the_route(repos):
    import dataclasses
    _settings(git_role="maintainer")
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, "push it via a PR")
    child = dataclasses.replace(perms, subagent_depth=1)
    A._grant_push_from(child, "push it directly to main", new_request=True)
    err = _gate(perms, "git push origin main")
    assert err and "pull request" in err, err


# ---------------------------------------------------------------------------
# Only the user sets it
# ---------------------------------------------------------------------------

def test_the_settings_file_is_a_protected_path(repos):
    perms = KitToolPermissions(workspace=repos["alpha"])
    assert perms.matches_path_protected(".delfin_settings.json")
    assert perms.matches_path_protected("sub/.delfin_settings.json")


@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_the_agent_cannot_write_the_settings_file(repos, tmp_path,
                                                  monkeypatch, mode):
    """Through the real executor: the write tool and three shell forms."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    target = home / US.SETTINGS_FILE_NAME
    planted = json.dumps({"agent": {"git_roles": {B_KEY: "maintainer"}}})
    perms = KitToolPermissions(workspace=repos["alpha"], mode=mode)
    perms.confirm_callback = lambda *a, **k: True
    attempts = [
        ("write_file", {"path": str(target), "content": planted}),
        ("bash", {"command": f"echo '{planted}' > ~/{US.SETTINGS_FILE_NAME}"}),
        ("bash", {"command": f"cd ~ && echo '{planted}' > {US.SETTINGS_FILE_NAME}"}),
        ("bash", {"command": (
            "python3 -c \"import pathlib; pathlib.Path.home().joinpath("
            f"'{US.SETTINGS_FILE_NAME}').write_text('{{}}')\"")}),
    ]
    for tool, args in attempts:
        A._doc_executor.execute(tool, args, perms)
        assert not target.exists() or B_KEY not in target.read_text(), (
            tool, args)


def test_the_handover_still_names_the_role_as_the_users():
    from delfin.agent import handover
    assert handover.kinds().get("git_role") is True


# ---------------------------------------------------------------------------
# /git role
# ---------------------------------------------------------------------------

def test_git_role_shows_and_sets_the_role_of_this_repository(repos):
    from delfin.agent import repl_commands as rc
    _settings(git_role="contributor")
    ctx = rc.ReplContext(workspace=repos["alpha"])
    shown = rc.dispatch("/git role", ctx).output
    assert "contributor" in shown and A_KEY in shown, shown
    out = rc.dispatch("/git role maintainer", ctx).output
    assert "maintainer" in out and A_KEY in out, out
    agent = US.load_settings()["agent"]
    assert agent["git_roles"] == {A_KEY: "maintainer"}
    assert agent["git_role"] == "contributor"
    assert "maintainer" in rc.dispatch("/git role", ctx).output
    assert _gate(_perms(repos["alpha"]), "git push origin main") is None
    err = _gate(_perms(repos["beta"]), "git push origin main")
    assert _refused_as_contributor(err), err
    rc.dispatch("/git role contributor", ctx)
    assert US.load_settings()["agent"]["git_roles"] == {A_KEY: "contributor"}


def test_git_role_refuses_a_word_that_is_no_role(repos):
    from delfin.agent import repl_commands as rc
    ctx = rc.ReplContext(workspace=repos["alpha"])
    out = rc.dispatch("/git role owner", ctx).output
    assert "usage" in out, out
    assert "git_roles" not in (US.load_settings().get("agent") or {}) or \
        not US.load_settings()["agent"]["git_roles"]


def test_git_role_outside_a_repository_says_so(tmp_path):
    from delfin.agent import repl_commands as rc
    plain = tmp_path / "plain"
    plain.mkdir()
    out = rc.dispatch("/git role maintainer", rc.ReplContext(workspace=plain)).output
    assert "not inside a git repository" in out, out


# ---------------------------------------------------------------------------
# The role is read before the command runs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("before", [
    "git remote set-url origin https://github.com/me/beta",
    "git remote add up https://github.com/me/beta",
    "git config remote.origin.pushurl https://github.com/me/beta",
    "git config url.https://github.com/me/beta.pushInsteadOf "
    "git@github.com:me/alpha.git",
    "git config core.editor true",
    "git config user.name x",
    "git -c core.hooksPath=h commit -qm x",
    "git -C . commit -qm x",
    "python3 retarget.py",
    "git fetch --set-upstream https://github.com/me/beta main",
    "git commit -qm \"$(git remote set-url origin https://github.com/me/beta)\"",
    "(git status)",
])
@pytest.mark.parametrize("sep", [" && ", "; ", "\n"])
def test_an_earlier_step_that_may_move_the_push_is_contributor(
        repos, before, sep):
    """The push URL is read before the line runs; an earlier step that is
    not one of the few that keep it leaves the push unresolved."""
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    cmd = before + sep + "git push origin main"
    err = _gate(_perms(repos["alpha"]), cmd)
    assert _refused_as_contributor(err), (cmd, err)


@pytest.mark.parametrize("cmd", [
    "git add f.txt && git commit -qm 'x; y' && git push origin main",
    "cd . && git add -A && git status && git push origin main",
    "git fetch origin && git pull --ff-only origin main && "
    "git log -1 --oneline && git push origin main",
    "git add f.txt\ngit commit -qm x\ngit push origin main",
])
def test_the_usual_steps_before_a_push_keep_the_role(repos, cmd):
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    assert _gate(_perms(repos["alpha"]), cmd) is None, cmd


def test_an_earlier_step_before_a_merge_is_contributor(repos, pr_view):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green())
    err = _gate(_perms(repos["alpha"]),
                "git remote set-url origin https://github.com/me/beta && "
                f"gh pr merge 12 --squash --match-head-commit {SHA}")
    assert err and "maintainer's decision" in err, err


@pytest.mark.parametrize("key", [
    "branch.main.pushRemote", "remote.pushDefault", "branch.main.remote"])
def test_a_push_naming_no_remote_follows_the_configured_push_remote(
        repos, key):
    """``git push`` with no remote goes where branch.<b>.pushRemote,
    remote.pushDefault or branch.<b>.remote says, which may be a URL."""
    _git("config", key, B_URL, cwd=repos["alpha"], env=repos["env"])
    _settings(git_role="contributor", git_roles={A_KEY: "maintainer"})
    err = _gate(_perms(repos["alpha"]), "git push")
    assert _refused_as_contributor(err), (key, err)
    assert _gate(_perms(repos["alpha"]), "git push origin main") is None


@pytest.mark.parametrize("cmd", [
    "git config core.hooksPath h",
    "git config --local core.hookspath /tmp/h",
    "git config set core.hooksPath h",
    "git config --global core.fsmonitor ./m.sh",
])
@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_a_hook_path_write_is_asked_in_every_mode(repos, cmd, mode):
    """A hook directory or fsmonitor program set now runs from the next
    auto-allowed ``git commit``/``git status``: asked as a git
    configuration change, refused with nobody to ask."""
    perms = _perms(repos["alpha"], mode=mode, grant=False, answer=False)
    err = _gate(perms, cmd)
    assert perms.asked and "GIT CONFIGURATION" in str(perms.asked[0]), (
        cmd, mode, perms.asked)
    assert err, (cmd, mode)
    headless = _perms(repos["alpha"], mode=mode, grant=False)
    err = _gate(headless, cmd)
    assert err and "git configuration" in err.lower(), (cmd, mode, err)


@pytest.mark.parametrize("cmd", [
    "git config --get core.hooksPath", "git config -l", "git config user.name"])
def test_reading_the_configuration_is_not_asked_as_a_hook_write(repos, cmd):
    perms = _perms(repos["alpha"], grant=False, answer=True)
    assert _gate(perms, cmd) is None
    assert not perms.asked, (cmd, perms.asked)


# ---------------------------------------------------------------------------
# A host that is not ASCII
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("url", [
    "https://gitlab.Kit.edu/me/x",      # KELVIN SIGN lower-cases to k
    "git@gıthub.com:me/alpha.git",      # dotless i
    "https://GİTHUB.com/me/alpha",
    "https://ｇｉｔｈｕｂ.com/me/alpha",
])
def test_a_host_that_is_not_ascii_has_no_key(url):
    from delfin.agent import git_role as GR
    assert GR.normalise_remote(url) is None, url
    assert GR.entry_key(url) is None, url


def test_an_entry_with_a_non_ascii_host_counts_for_nothing(repos):
    from delfin.agent import git_role as GR
    assert GR.entry_key("gitlab.Kit.edu/me/x") is None


# ---------------------------------------------------------------------------
# A pull request with no checks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("checks", [[], ["not a check"]])
def test_a_merge_with_no_checks_asks_when_attended(repos, pr_view, checks):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green(checks=checks))
    perms = _perms(repos["alpha"], mode="default", answer=False)
    err = _merge(perms)
    assert perms.asked and "no checks" in str(perms.asked[-1]), perms.asked
    assert err and err.startswith("blocked:"), err
    perms = _perms(repos["alpha"], mode="default", answer=True)
    assert _merge(perms) is None
    assert perms.asked and "no checks" in str(perms.asked[-1])


@pytest.mark.parametrize("mode", ["default", "bypassPermissions"])
def test_a_merge_with_no_checks_is_refused_unattended(repos, pr_view, mode):
    _settings(git_roles={A_KEY: "maintainer"})
    pr_view(_green(checks=[]))
    perms = _perms(repos["alpha"], mode=mode)
    if mode == "bypassPermissions":
        perms.confirm_callback = lambda *a, **k: True
    err = _merge(perms)
    assert err and "no checks" in err and err.startswith("blocked:"), err


# ---------------------------------------------------------------------------
# A resumed session keeps the route
# ---------------------------------------------------------------------------

def test_a_resumed_session_keeps_the_pr_route(repos, tmp_path, monkeypatch):
    import types
    from delfin.agent import cli as agent_cli
    from delfin.agent import session_store as S
    from delfin.agent.engine import AgentEngine

    def _engine(perms):
        eng = AgentEngine.__new__(AgentEngine)
        eng.mode = "solo"
        eng.route = ["solo_agent"]
        eng.current_role_index = 0
        eng.role_outputs = {}
        eng.compaction_summaries = {}
        eng.messages = []
        eng.token_usage = {"input": 0, "output": 0}
        eng.cost_usd = 0.0
        eng.session_id = "pr-route-rt-1"
        for spec in AgentEngine._SESSION_FIELDS:
            setattr(eng, spec.attr, spec.reset())
        eng.client = types.SimpleNamespace(_permissions=perms)
        return eng

    monkeypatch.setattr(S, "_SESSIONS_DIR", tmp_path / "sessions")
    monkeypatch.setattr(AgentEngine, "_sync_task_session", lambda self: None)
    _settings(git_role="maintainer")
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, "über PR bitte")
    assert perms.pr_route_repos == {A_KEY}
    sid = agent_cli._save_session(_engine(perms), tmp_path)

    fresh = _perms(repos["alpha"], grant=False)
    resumed = _engine(fresh)
    agent_cli._resume_or_create(resumed, types.SimpleNamespace(
        session=sid, mode="solo", verbose=False, cwd=str(tmp_path)))
    assert fresh.pr_route_repos == {A_KEY}
    _typed(fresh, "jetzt push es")
    err = _gate(fresh, "git push origin main")
    assert err and "pull request" in err, err


# ---------------------------------------------------------------------------
# The route matcher
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("phrase", [
    "über einen PR", "mach einen PR", "lieber per PR",
    "bitte nicht direkt pushen", "open a pull request instead",
    "go through a PR", "erstell bitte einen PR dafür",
])
def test_the_route_is_read_in_its_common_phrasings(repos, phrase):
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, phrase)
    assert perms.pr_route_repos == {A_KEY}, phrase


@pytest.mark.parametrize("phrase", [
    "nicht über PR", "kein PR nötig, direkt auf main", "mach keinen PR",
    "no PR needed, push to main", "don't open a PR", "without a pull request",
])
def test_a_declined_pr_does_not_set_the_route(repos, phrase):
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, phrase)
    assert perms.pr_route_repos == set(), phrase


def test_kein_pr_noetig_lifts_the_route_only_as_a_typed_message(repos):
    _settings(git_role="maintainer")
    perms = _perms(repos["alpha"], grant=False)
    _typed(perms, "über PR bitte")
    _harness(perms, "kein PR nötig, direkt auf main")
    assert perms.pr_route_repos == {A_KEY}
    _typed(perms, "kein PR nötig, direkt auf main")
    assert perms.pr_route_repos == set()
