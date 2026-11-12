"""DELFIN's clusters: which ones the remote layer will drive, and how it says why not.

Input: a home under ``tmp_path`` with ctun's config at its default place
(``$XDG_CONFIG_HOME/cluster-tunnel/config.yaml``), the ``remote`` block of
DELFIN's settings as a dict, and a fake ``ssh`` first on PATH (it answers
``ssh -O check`` from a marker file and runs the remote command line with a
local ``sh``). Output: ``delfin.remote.clusters.load`` records and what
``delfin.remote.connection`` returns or raises.

Semantics pinned here:

* a cluster is DELFIN's when ``remote.clusters`` names it, and usable only
  when ctun has its host key pinned, a ctun budget guard reads every
  command (not only its first word), and every file ctun's ssh reads lies
  in ctun's directories or ~/.ssh, out of every agent's reach;
* the problems are listed, never raised, each with what to do;
* remote directories default to ~/delfin/jobs and ~/delfin/installs and are
  never the home itself, ``/`` or a path with ``..``;
* ``connect`` refuses an unknown name, an unusable cluster, and any call from
  inside an agent session; ctun's errors arrive as DELFIN's own.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from cluster_tunnel import api as ctun_api
from delfin.remote import clusters, connection, errors

_FAKE_SSH = """\
#!{python}
import json, os, subprocess, sys

state = os.environ["FAKE_SSH_STATE"]
args = sys.argv[1:]
with open(os.path.join(state, "argv.jsonl"), "a", encoding="utf-8") as log:
    log.write(json.dumps(args) + "\\n")
if "-O" in args:
    if args[args.index("-O") + 1] == "check" and os.path.exists(os.path.join(state, "master")):
        sys.exit(0)
    print("Control socket connect(sock): No such file or directory", file=sys.stderr)
    sys.exit(255)
sys.exit(subprocess.run(["sh", "-c", args[-1]], cwd=os.path.join(state, "home")).returncode)
"""


class Home:
    """A user's home with ctun configured, and DELFIN's settings for it."""

    def __init__(self, root: Path):
        self.root = root
        self.ctun_dir = root / ".config" / "cluster-tunnel"
        self.ctun_dir.mkdir(parents=True)
        self.config = self.ctun_dir / "config.yaml"
        self.known_hosts = self.ctun_dir / "known_hosts"
        self.state = root.parent / "ssh-state"
        (self.state / "home").mkdir(parents=True)

    def ctun(self, entry: str = "", *, pinned: bool = True, name: str = "kcist") -> None:
        pin = (f"    strict_host_key_checking: 'yes'\n"
               f"    known_hosts_file: {self.known_hosts}\n") if pinned else ""
        self.config.write_text(
            f"clusters:\n  {name}:\n    host: login.example\n    user: me\n{pin}{entry}",
            encoding="utf-8")

    @staticmethod
    def settings(**clusters_) -> dict:
        return {"remote": {"clusters": clusters_ or {"kcist": {}}}}

    def log_in(self) -> None:
        (self.state / "master").touch()


@pytest.fixture
def home(tmp_path, monkeypatch):
    root = tmp_path / "home"
    root.mkdir()
    monkeypatch.setenv("HOME", str(root))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(root / ".config"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(root / ".cache"))
    monkeypatch.delenv("CTUN_CONFIG", raising=False)
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    h = Home(root)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "ssh"
    fake.write_text(_FAKE_SSH.format(python=sys.executable), encoding="utf-8")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:/usr/bin:/bin")
    monkeypatch.setenv("FAKE_SSH_STATE", str(h.state))
    return h


def _problems(home: Home, **settings) -> tuple[str, ...]:
    return clusters.load(settings=home.settings(**settings))["kcist"].problems


# ---------------------------------------------------------------------------
# Which clusters are DELFIN's, and which are usable
# ---------------------------------------------------------------------------

def test_without_a_remote_block_delfin_has_no_clusters(home):
    home.ctun()
    assert clusters.load(settings={}) == {}
    assert clusters.load(settings={"remote": {"clusters": []}}) == {}


def test_a_pinned_cluster_is_usable_with_the_default_directories(home):
    home.ctun()
    cluster = clusters.load(settings=home.settings())["kcist"]
    assert cluster.usable, cluster.problems
    assert (cluster.target, cluster.pinned) == ("me@login.example", True)
    assert (cluster.work_dir, cluster.install_root) == ("~/delfin/jobs", "~/delfin/installs")
    assert cluster.allowed_partitions == ()


def test_a_cluster_ctun_does_not_know_is_named_with_the_fix(home):
    home.ctun(name="other")
    (problem,) = _problems(home)
    assert "has no cluster 'kcist'" in problem and "ctun config" in problem


def test_an_unpinned_host_key_keeps_the_cluster_out(home):
    home.ctun(pinned=False)
    (problem,) = _problems(home)
    assert "ctun -t kcist hostkey --pin" in problem


def test_a_budget_guard_that_reads_only_first_words_keeps_it_out(home):
    budget = home.ctun_dir / "budget" / "kcist.sh"
    budget.parent.mkdir()
    budget.write_text("echo 0\n", encoding="utf-8")
    home.ctun("    budget:\n      session_limit: 10\n")
    (problem,) = _problems(home)
    assert "guard_scope: every-command" in problem
    home.ctun("    budget:\n      session_limit: 10\n      guard_scope: every-command\n")
    assert _problems(home) == ()


@pytest.mark.parametrize("key", ["ssh_config", "known_hosts_file", "budget script"])
def test_a_file_ctun_reads_where_an_agent_could_write_it_keeps_it_out(home, tmp_path, key):
    """A ProxyCommand written into one of these would run in DELFIN's own
    process, outside every cage."""
    outside = tmp_path / "workspace" / "planted"
    outside.parent.mkdir()
    outside.write_text("", encoding="utf-8")
    if key == "ssh_config":
        home.ctun(f"    ssh_config: {outside}\n")
    elif key == "known_hosts_file":
        home.config.write_text(
            "clusters:\n  kcist:\n    host: login.example\n    user: me\n"
            f"    strict_host_key_checking: 'yes'\n    known_hosts_file: {outside}\n",
            encoding="utf-8")
    else:
        home.ctun(f"    budget:\n      script: {outside}\n      guard_scope: every-command\n")
    (problem,) = _problems(home)
    assert str(outside) in problem and "outside ctun's directories" in problem


def test_an_ssh_config_in_dot_ssh_is_out_of_reach_already(home):
    ssh_dir = home.root / ".ssh"
    ssh_dir.mkdir()
    (ssh_dir / "config_hpc").write_text("", encoding="utf-8")
    home.ctun(f"    ssh_config: {ssh_dir / 'config_hpc'}\n")
    assert _problems(home) == ()


def test_a_ctun_config_outside_its_directory_keeps_every_cluster_out(home, tmp_path):
    elsewhere = tmp_path / "workspace" / "ctun.yaml"
    elsewhere.parent.mkdir()
    home.ctun()
    elsewhere.write_text(home.config.read_text(encoding="utf-8"), encoding="utf-8")
    cluster = clusters.load(settings=home.settings(), ctun_config=elsewhere)["kcist"]
    assert any("ctun's config" in p and "outside" in p for p in cluster.problems)


def test_an_unreadable_ctun_config_is_a_problem_not_an_exception(home):
    home.config.write_text("clusters: [", encoding="utf-8")
    (problem,) = _problems(home)
    assert problem.startswith("ctun's config cannot be read")


@pytest.mark.parametrize("value, expected", [
    ("jobs", "~/jobs"),
    ("~/delfin/jobs/", "~/delfin/jobs"),
    ("/scratch/me/jobs", "/scratch/me/jobs"),
    ("~//a/./b", "~/a/b"),
])
def test_a_work_dir_is_kept_under_home_or_absolute(home, value, expected):
    home.ctun()
    cluster = clusters.load(settings=home.settings(kcist={"work_dir": value}))["kcist"]
    assert cluster.usable, cluster.problems
    assert cluster.work_dir == expected


@pytest.mark.parametrize("value", ["~", "~/", "/", "", "  ", "~/a/../..", "../x",
                                   "~bob/x", "a\nb", 7])
def test_a_work_dir_that_is_everything_or_escapes_is_refused(home, value):
    home.ctun()
    problems = _problems(home, kcist={"work_dir": value})
    assert len(problems) == 1 and "work_dir" in problems[0], problems


def test_unknown_settings_and_bad_partitions_are_named(home):
    home.ctun()
    problems = _problems(home, kcist={"workdir": "~/x", "allowed_partitions": ["small", "a b"]})
    assert any("unknown setting 'workdir'" in p for p in problems)
    assert any("allowed_partitions" in p and "a b" in p for p in problems)


# ---------------------------------------------------------------------------
# connect() and the connection
# ---------------------------------------------------------------------------

def test_connect_refuses_a_name_the_settings_do_not_list(home):
    home.ctun()
    with pytest.raises(errors.UnknownCluster, match="lists: kcist"):
        connection.connect("horeka", settings=home.settings())


def test_connect_refuses_an_unusable_cluster_with_its_problems(home):
    home.ctun(pinned=False)
    with pytest.raises(errors.ClusterNotUsable) as caught:
        connection.connect("kcist", settings=home.settings())
    assert caught.value.code == "cluster_not_usable"
    assert any("hostkey --pin" in p for p in caught.value.problems)


def test_a_command_an_agent_runs_cannot_connect(home, monkeypatch):
    home.ctun()
    home.log_in()
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "1234")
    with pytest.raises(errors.RemoteRefused) as caught:
        connection.connect("kcist", settings=home.settings())
    assert caught.value.code == "remote_refused"
    assert not (home.state / "argv.jsonl").exists(), "ssh was started anyway"


def test_without_a_login_the_connection_says_login_required(home):
    home.ctun()
    conn = connection.connect("kcist", settings=home.settings())
    assert conn.is_live() is False
    with pytest.raises(errors.LoginRequired) as caught:
        conn.run(["true"])
    assert caught.value.code == "login_required"


def test_a_logged_in_connection_runs_with_the_pinned_key_only(home):
    home.ctun()
    home.log_in()
    conn = connection.connect("kcist", settings=home.settings())
    result = conn.run(["sh", "-c", "echo hello; exit 4"])
    assert (result.exit_code, result.stdout) == (4, "hello\n")
    argvs = [json.loads(line) for line in
             (home.state / "argv.jsonl").read_text(encoding="utf-8").splitlines()]
    assert all("StrictHostKeyChecking=yes" in a for a in argvs), argvs


_SPENT = ctun_api.Decision(allowed=False, reason="used 12 of 10 jobh", used=12.0,
                           limit=10.0, unit="jobh")


@pytest.mark.parametrize("exc, mapped", [
    (ctun_api.LoginRequired("kcist"), errors.LoginRequired),
    (ctun_api.BudgetExhausted("kcist", ["sbatch", "job.sh"], _SPENT), errors.BudgetExhausted),
    (ctun_api.BudgetGuardError("kcist", ["sbatch", "job.sh"], _SPENT), errors.BudgetGuardError),
    (ctun_api.HostKeyError("kcist presents another key"), errors.HostKeyMismatch),
    (ctun_api.UnknownClusterError("kcist"), errors.ClusterNotUsable),
    (ctun_api.TransportError("rsync failed"), errors.TransportError),
    (OSError("ssh: not found"), errors.TransportError),
])
def test_ctuns_errors_arrive_as_delfins(home, monkeypatch, exc, mapped):
    home.ctun()
    conn = connection.connect("kcist", settings=home.settings())

    def boom(*args, **kwargs):
        raise exc

    monkeypatch.setattr(ctun_api, "run", boom)
    with pytest.raises(mapped):
        conn.run(["sbatch", "job.sh"])
