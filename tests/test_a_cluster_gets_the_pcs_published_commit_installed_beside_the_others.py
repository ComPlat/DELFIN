"""DELFIN installs itself on a cluster per commit: the PC's own, published, side by side.

Input: a throwaway DELFIN-shaped git checkout (pyproject.toml, and a stub
``delfin/installers/install_delfin.sh`` that records its arguments) with a
bare repository standing in for the install source; a cluster reached
through a fake ``ssh`` that runs every remote command line in a local home,
so ctun's rsync uploads and ``bash -s`` script feeds really run; a fake
``squeue`` that names the commits of queued jobs. Output: what
``delfin.remote.install`` returns or raises, and the files on the "cluster".

Semantics pinned here:

* a commit on a branch of the install source, or behind one, is uploaded as
  exactly that commit (``git archive``) into ``<install_root>/<sha>/`` with
  ``.delfin-commit``, and installed there with ``--repo`` that directory and
  ``--prefix`` the software dir; a second call finds it complete;
* a commit on no branch of the source is never uploaded;
* changes in the PC's working tree are noted and do not travel;
* a failed install is reported with its log and started again only on
  request;
* ``remove_unused`` keeps what it is told to keep and what a queued job names
  (``--comment=delfin-commit=<sha>``), removes the other installs of
  DELFIN's, and never touches a directory without ``.delfin-commit``;
* every install and removal is in the audit log.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.remote import connection, errors, install

_FAKE_SSH = """\
#!{python}
# Parses argv as ssh does (options, destination, then the remote words joined
# into one command line) and runs that line with sh in $HOME.
import os, subprocess, sys
WITH_VALUE = set("bcDEeFIiJLlmOoPpQRSWw")
args, i = sys.argv[1:], 0
while i < len(args) and args[i].startswith("-"):
    flag = args[i]
    if flag == "-O" and args[i + 1] == "check":
        sys.exit(0 if os.path.exists(os.environ["FAKE_MASTER"]) else 255)
    i += 2 if len(flag) == 2 and flag[1] in WITH_VALUE else 1
command = " ".join(args[i + 1:])
sys.exit(subprocess.run(["sh", "-c", command], cwd=os.environ["HOME"]).returncode)
"""

_FAKE_SQUEUE = """\
#!/bin/sh
cat "$FAKE_SQUEUE" 2>/dev/null
"""

_STUB_INSTALLER = """\
#!/usr/bin/env bash
# Stands in for install_delfin.sh: records its arguments and makes a venv
# whose python "imports DELFIN" (exits 0). FAIL: fails before any venv.
# TOOLFAIL: a tool fails, as std2 without gfortran. MPIFAIL: OpenMPI fails.
echo "stub install $*"
[ -f FAIL ] && { echo "stub failure" >&2; exit 3; }
mkdir -p "$2/.venv/bin" && printf '#!/bin/sh\\nexit 0\\n' > "$2/.venv/bin/python"
chmod +x "$2/.venv/bin/python"
mkdir -p "$4" && echo "$2" >> "$4/installed-by-stub"
[ -f TOOLFAIL ] && { echo "[delfin-install] FAILED qm: xtb std2 (missing: std2)"; \\
                     echo "[delfin-install] FAILED   install of tools (see the lines above)"; exit 1; }
[ -f MPIFAIL ] && { echo "[delfin-install] FAILED   openmpi"; exit 1; }
exit 0
"""


def _git(cwd, *args):
    return subprocess.run(["git", "-C", str(cwd), *args], check=True, capture_output=True,
                          text=True).stdout.strip()


class World:
    def __init__(self, tmp: Path):
        self.home = tmp / "home"
        self.source = tmp / "source.git"
        self.checkout = tmp / "checkout"
        self.squeue = tmp / "squeue.txt"
        self.audit = self.home / ".delfin" / "audit.log"

    def commit(self, message: str, files: dict[str, str]) -> install.Commit:
        for name, text in files.items():
            path = self.checkout / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")
        _git(self.checkout, "add", "-A")
        _git(self.checkout, "commit", "-q", "-m", message)
        return install.pc_commit(self.checkout)

    def publish(self, branch: str = "main") -> None:
        _git(self.checkout, "push", "-q", str(self.source), f"HEAD:refs/heads/{branch}")

    def installed(self, sha: str) -> Path:
        return self.home / "inst" / sha

    def conn(self):
        settings = {"remote": {"clusters": {"kcist": {"install_root": "~/inst",
                                                      "software_dir": "~/sw"}}}}
        return connection.connect("kcist", settings=settings)


@pytest.fixture
def world(tmp_path, monkeypatch):
    w = World(tmp_path)
    w.home.mkdir()
    ctun_dir = w.home / ".config" / "cluster-tunnel"
    ctun_dir.mkdir(parents=True)
    (ctun_dir / "config.yaml").write_text(
        "clusters:\n  kcist:\n    host: login.example\n    user: me\n"
        "    strict_host_key_checking: 'yes'\n"
        f"    known_hosts_file: {ctun_dir / 'known_hosts'}\n", encoding="utf-8")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, text in (("ssh", _FAKE_SSH.format(python=sys.executable)),
                       ("squeue", _FAKE_SQUEUE)):
        (bin_dir / name).write_text(text, encoding="utf-8")
        (bin_dir / name).chmod(0o755)
    master = tmp_path / "master"
    master.touch()
    for key, value in (("HOME", w.home), ("XDG_CONFIG_HOME", w.home / ".config"),
                       ("XDG_CACHE_HOME", w.home / ".cache"), ("FAKE_MASTER", master),
                       ("FAKE_SQUEUE", w.squeue),
                       ("PATH", f"{bin_dir}:/usr/bin:/bin"),
                       ("GIT_AUTHOR_NAME", "t"), ("GIT_AUTHOR_EMAIL", "t@example.org"),
                       ("GIT_COMMITTER_NAME", "t"), ("GIT_COMMITTER_EMAIL", "t@example.org")):
        monkeypatch.setenv(key, str(value))
    monkeypatch.delenv("CTUN_CONFIG", raising=False)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    subprocess.run(["git", "init", "-q", "--bare", str(w.source)], check=True)
    subprocess.run(["git", "init", "-q", str(w.checkout)], check=True)
    w.commit("first", {"pyproject.toml": 'name = "delfin-complat"\n',
                       "delfin/installers/install_delfin.sh": _STUB_INSTALLER,
                       "delfin/version.txt": "one\n"})
    return w


def _ensure(world, commit, **kwargs):
    return install.ensure_installed(world.conn(), commit=commit, sources=[str(world.source)],
                                    poll_seconds=0.05, wait_seconds=60, **kwargs)


def test_a_published_commit_is_installed_into_its_own_directory(world):
    commit = install.pc_commit(world.checkout)
    world.publish()
    outcome = _ensure(world, commit)
    assert (outcome.state, outcome.installed_now) == ("complete", True)
    assert outcome.path == f"~/inst/{commit.sha}" and outcome.source == str(world.source)
    tree = world.installed(commit.sha)
    marker = json.loads((tree / ".delfin-commit").read_text(encoding="utf-8"))
    assert marker["commit"] == commit.sha and marker["source"] == str(world.source)
    assert (tree / "delfin" / "version.txt").read_text(encoding="utf-8") == "one\n"
    assert not (tree / ".git").exists(), "the commit's tree, not a clone"
    log = (tree / ".delfin-install.log").read_text(encoding="utf-8")
    assert f"--repo {tree} --prefix {world.home / 'sw'} --profile standard --yes" in log
    assert (world.home / "sw" / "installed-by-stub").exists()
    again = _ensure(world, commit)
    assert (again.state, again.installed_now) == ("complete", False)
    assert "cluster_install" in world.audit.read_text(encoding="utf-8")


def test_a_second_commit_installs_beside_the_first(world):
    first = install.pc_commit(world.checkout)
    world.publish()
    _ensure(world, first)
    second = world.commit("second", {"delfin/version.txt": "two\n"})
    world.publish()
    _ensure(world, second)
    assert install.installs(world.conn()) == {first.sha: "complete", second.sha: "complete"}
    assert (world.installed(first.sha) / "delfin" / "version.txt").read_text() == "one\n"
    assert (world.installed(second.sha) / "delfin" / "version.txt").read_text() == "two\n"


def test_a_commit_behind_a_branch_head_counts_as_published(world):
    older = install.pc_commit(world.checkout)
    world.commit("newer", {"delfin/version.txt": "newer\n"})
    world.publish()
    assert install.published_source(older, [str(world.source)]) == str(world.source)


def test_a_commit_on_no_branch_of_the_source_is_never_uploaded(world):
    world.publish()
    local = world.commit("only on the PC", {"delfin/version.txt": "local\n"})
    with pytest.raises(errors.CommitNotPublished) as caught:
        _ensure(world, local)
    assert caught.value.code == "commit_not_published"
    assert not world.installed(local.sha).exists()


def test_a_commit_only_in_another_repository_is_not_published(world, tmp_path):
    other = tmp_path / "somebodys-fork.git"
    subprocess.run(["git", "init", "-q", "--bare", str(other)], check=True)
    world.publish()
    foreign = world.commit("pushed elsewhere", {"delfin/version.txt": "foreign\n"})
    _git(world.checkout, "push", "-q", str(other), "HEAD:refs/heads/main")
    assert install.published_source(foreign, [str(world.source)]) is None
    assert install.published_source(foreign, [str(world.source), str(other)]) == str(other)


def test_changes_in_the_pcs_working_tree_are_noted_and_do_not_travel(world):
    world.publish()
    (world.checkout / "delfin" / "version.txt").write_text("edited, not committed\n")
    commit = install.pc_commit(world.checkout)
    assert commit.dirty
    outcome = _ensure(world, commit)
    assert any("changes beyond commit" in n for n in outcome.notes)
    assert (world.installed(commit.sha) / "delfin" / "version.txt").read_text() == "one\n"


def test_a_failed_install_is_reported_and_started_again_only_on_request(world):
    commit = world.commit("breaks the install", {"FAIL": "1\n"})
    world.publish()
    with pytest.raises(errors.InstallFailed) as caught:
        _ensure(world, commit)
    assert caught.value.code == "install_failed" and "stub failure" in caught.value.log
    with pytest.raises(errors.InstallFailed, match="failed before"):
        _ensure(world, commit)
    log = world.installed(commit.sha) / ".delfin-install.log"
    log.write_text("", encoding="utf-8")
    with pytest.raises(errors.InstallFailed):
        _ensure(world, commit, retry=True)
    assert "stub install" in log.read_text(encoding="utf-8"), "retry ran the installer again"


def test_a_tool_that_fails_is_a_warning_of_a_complete_install(world):
    """std2 is built with gfortran, which kcist's login node lacks: DELFIN
    works without it, and the calculations that need it say so."""
    commit = world.commit("std2 cannot be built here", {"TOOLFAIL": "1\n"})
    world.publish()
    outcome = _ensure(world, commit)
    assert outcome.state == "complete"
    # The installer's own summary line ("see the lines above") is left out.
    assert outcome.warnings == ["FAILED qm: xtb std2 (missing: std2)"]
    assert _ensure(world, commit).warnings == outcome.warnings


@pytest.mark.parametrize("marker", ["MPIFAIL", "FAIL"])
def test_no_openmpi_or_no_working_venv_is_a_failed_install(world, marker):
    commit = world.commit(f"breaks with {marker}", {marker: "1\n"})
    world.publish()
    with pytest.raises(errors.InstallFailed):
        _ensure(world, commit)


def test_unused_installs_go_and_kept_or_queued_ones_stay(world):
    shas = []
    for n in range(3):
        if n:
            world.commit(f"c{n}", {"delfin/version.txt": f"{n}\n"})
        world.publish()
        commit = install.pc_commit(world.checkout)
        _ensure(world, commit)
        shas.append(commit.sha)
    foreign = world.home / "inst" / "not-delfins"
    foreign.mkdir()
    world.squeue.write_text(f"{install.job_comment(shas[0])}\n\nsomething else\n")
    removed = install.remove_unused(world.conn(), keep=[shas[2]])
    assert removed == [shas[1]]
    assert world.installed(shas[0]).exists() and world.installed(shas[2]).exists()
    assert not world.installed(shas[1]).exists()
    assert foreign.exists()
    assert '"auto_remove"' in world.audit.read_text(encoding="utf-8")


def test_the_script_removes_nothing_that_is_not_delfins(world, tmp_path):
    victim = world.home / "inst" / "precious"
    victim.mkdir(parents=True)
    (victim / "data").write_text("keep me")
    r = subprocess.run(["bash", str(install.SCRIPT), "remove", "~/inst/precious"],
                       capture_output=True, text=True, cwd=world.home)
    assert r.returncode == 2 and (victim / "data").exists()


def test_the_job_comment_names_the_commit():
    assert install.job_comment("ab" * 20) == "delfin-commit=" + "ab" * 20
