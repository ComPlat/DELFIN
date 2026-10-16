"""DELFIN's own git reads run no program the repository names.

DELFIN reads repository state from outside the sandbox (launch banner,
system prompt, evidence fingerprint, introspection, held-branch recovery,
worktrees). The repository is the session's workspace, so its config and
attributes are not DELFIN's to trust. ``safe_git`` overrides every
program-naming key the repository defines, plus a fixed set, and these
tests hold that: on the command line it builds, and on what git then does.
The behavioural checks use ``false`` as the configured program -- harmless,
and git fails loudly when it runs it, so "not run" is observable.
"""
from __future__ import annotations

import os
import subprocess

import pytest

from delfin.agent import safe_git as SG


@pytest.fixture
def repo(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    r = tmp_path / "repo"
    r.mkdir()
    env = {**os.environ, "HOME": str(home)}
    run = lambda *a: subprocess.run(["git", *a], cwd=r, env=env, check=True,
                                    capture_output=True)
    run("init", "-q")
    run("-c", "user.email=t@t", "-c", "user.name=t", "commit", "-q",
        "--allow-empty", "-m", "init")
    return r, run


def _overrides(argv):
    return [argv[i + 1] for i, a in enumerate(argv) if a == "-c"]


def test_the_fixed_overrides_come_before_the_subcommand(repo):
    r, _ = repo
    argv = SG.argv_for(["status", "--porcelain"], str(r))
    sub = argv.index("status")
    assert argv[:2] == ["git", "-C"]
    assert "--no-optional-locks" in argv[:sub]
    for kv in ("core.fsmonitor=false", "core.hooksPath=/dev/null",
               "core.pager=cat", "diff.external="):
        assert kv in _overrides(argv[:sub]), kv
    assert "--ignore-submodules=all" in argv[sub:]


@pytest.mark.parametrize("sub,flags", [
    ("diff", ["--no-textconv", "--no-ext-diff"]),
    ("log", ["--no-textconv", "--no-ext-diff"]),
    ("show", ["--no-textconv", "--no-ext-diff"]),
    ("blame", ["--no-textconv"]),
    ("grep", ["--no-textconv"]),
])
def test_a_reading_command_gets_no_textconv_and_no_external_diff(repo, sub, flags):
    r, _ = repo
    argv = SG.argv_for([sub], str(r))
    i = argv.index(sub)
    assert argv[i + 1:i + 1 + len(flags)] == flags


def test_a_repository_program_key_is_overridden_empty(repo):
    r, run = repo
    run("config", "filter.lfsish.clean", "false")
    run("config", "core.sshCommand", "false")
    run("config", "diff.pic.textconv", "false")
    ov = _overrides(SG.argv_for(["status"], str(r)))
    assert "filter.lfsish.clean=" in ov
    assert "filter.lfsish.required=false" in ov
    assert "core.sshcommand=" in ov    # git lists section and name lower-cased
    assert "diff.pic.textconv=" in ov


def test_an_included_file_is_read_as_the_repository_s_own(repo, tmp_path):
    r, run = repo
    inc = tmp_path / "more.cfg"
    inc.write_text("[diff \"x\"]\n\ttextconv = false\n")
    run("config", "include.path", str(inc))
    run("config", f"includeIf.gitdir:{r}/.git.path", str(inc))
    ov = _overrides(SG.argv_for(["status"], str(r)))
    assert "diff.x.textconv=" in ov
    # git refuses an include key on the command line; the included entry
    # is what is overridden, and git still runs.
    assert not any(o.startswith(("include.", "includeif.")) for o in ov)
    res = SG.run(["rev-parse", "--show-toplevel"], r)
    assert res.returncode == 0, res.stderr


def test_the_user_s_own_global_configuration_is_left_alone(repo, tmp_path):
    r, _ = repo
    (tmp_path / "home" / ".gitconfig").write_text(
        "[filter \"userlfs\"]\n\tclean = false\n")
    ov = _overrides(SG.argv_for(["status"], str(r)))
    assert not any(o.startswith("filter.userlfs.") for o in ov)


def test_a_key_the_command_line_cannot_name_refuses_the_call(repo):
    r, _ = repo
    cfg = r / ".git" / "config"
    cfg.write_text(cfg.read_text() + "[filter \"a=b\"]\n\tclean = false\n")
    with pytest.raises(SG.Refused):
        SG.argv_for(["status"], str(r))
    res = SG.run(["status", "--porcelain"], r)
    assert res.returncode == 128 and "cannot override" in res.stderr


def test_a_required_filter_the_repository_names_is_not_run(repo):
    """Control and check in one: a tracked file whose stat changed is passed
    through its clean filter by plain ``git status``; with ``false`` as a
    required filter that fails the command. Through safe_git the filter is
    blanked, so status compares the file as it is and reports it modified."""
    r, run = repo
    (r / "a.dat").write_text("one\n")
    run("add", "a.dat")
    run("-c", "user.email=t@t", "-c", "user.name=t", "commit", "-q", "-m", "a")
    run("config", "filter.f.clean", "false")
    run("config", "filter.f.required", "true")
    (r / ".gitattributes").write_text("*.dat filter=f\n")
    (r / "a.dat").write_text("two\n")
    plain = subprocess.run(["git", "status", "--porcelain"], cwd=r,
                           capture_output=True, text=True)
    assert plain.returncode != 0, plain.stdout
    res = SG.run(["status", "--porcelain"], r)
    assert res.returncode == 0, res.stderr
    assert " M a.dat" in res.stdout


def test_no_inherited_git_variable_reaches_git(repo, monkeypatch, tmp_path):
    r, _ = repo
    monkeypatch.setenv("GIT_DIR", str(tmp_path / "elsewhere"))
    monkeypatch.setenv("GIT_CONFIG_PARAMETERS", "'core.pager'='false'")
    res = SG.run(["rev-parse", "--git-dir"], r)
    assert res.returncode == 0, res.stderr
    assert res.stdout.strip() in (".git", str(r / ".git"))


def test_a_failure_is_a_result_not_an_exception(tmp_path):
    res = SG.run(["status"], tmp_path / "missing")
    assert res.returncode != 0
    assert SG.out(tmp_path / "missing", "status") == ""
    assert SG.is_repository(tmp_path / "missing") is False
