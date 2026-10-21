"""A git remote the workspace names widens the egress only for reads.

The workspace's repository configuration is the workspace's: a cloned or
unpacked repository can name any remote. Before this change every remote
host became "known" to the egress proxy (an HTTP remote included, so a
plain-HTTP POST to it passed in an attended session), and every SSH
remote became a direct TCP route that bypasses the proxy entirely, on
any port, attended or not.

Input: a repository under ``tmp_path`` with remotes added by
``git remote add`` (no network), name resolution replaced by a label
function, and the settings' ``allowed_domains``. Output:
``egress_proxy.configured_domains``, ``direct_tcp_endpoints``, and the
proxy's verdict (``EgressProxy.decide``) for a plain-HTTP write.

Semantics:

* an HTTPS remote's host is known for the proxy (reads pass as for any
  configured host; an upload is judged as for every known host); an
  ``http://`` remote's host is not;
* an SSH remote gets a direct TCP route only when its host is named by
  something the workspace does not control: the default allow-list
  (GitHub) or ``agent.sandbox_network.allowed_domains``; any other SSH
  remote gets none, and the reason is recorded once;
* the remotes are read with ``safe_git`` (no program the repository
  names runs).
"""
from __future__ import annotations

import subprocess

import pytest

from delfin.agent import egress_proxy as EP
from delfin.agent import security_events as SE


@pytest.fixture
def repo(tmp_path, monkeypatch):
    r = tmp_path / "ws"
    r.mkdir()
    env = {"GIT_CONFIG_NOSYSTEM": "1", "HOME": str(tmp_path),
           "PATH": "/usr/bin:/bin"}
    subprocess.run(["git", "init", "-q", str(r)], check=True, env=env)
    for name, url in (
            ("origin", "https://git.example.test/team/repo.git"),
            ("plain", "http://plain.example.test/team/repo.git"),
            ("sshport", "ssh://git@attacker.example.test:443/x.git"),
            ("scp", "git@scp.example.test:team/repo.git"),
            ("gh", "git@github.com:team/repo.git")):
        subprocess.run(["git", "-C", str(r), "remote", "add", name, url],
                       check=True, env=env)
    monkeypatch.setattr(EP, "_GIT_REMOTES", {})
    monkeypatch.setattr(EP, "_DIRECT_TCP", {})
    monkeypatch.setattr(EP, "slurm_controllers", lambda: [])
    monkeypatch.setattr(EP, "_resolve", lambda host, port: [f"{host}:{port}"])
    monkeypatch.setattr(EP, "_settings", lambda: {})
    for var in [v for v in __import__("os").environ
                if v.upper().endswith("LICENSE_FILE")]:
        monkeypatch.delenv(var)
    SE.clear()
    yield r
    SE.clear()


def test_an_https_remote_is_known_and_an_http_one_is_not(repo):
    hosts = EP.configured_domains(repo)
    assert "git.example.test" in hosts
    assert "plain.example.test" not in hosts


def test_an_ssh_remote_gets_a_route_only_where_the_user_named_it(repo,
                                                                 monkeypatch):
    direct = set(EP.direct_tcp_endpoints(repo))
    assert "github.com:22" in direct, direct
    assert not any(d.startswith(("attacker.example.test",
                                 "scp.example.test")) for d in direct), direct
    assert any(e.kind == "network" and "attacker.example.test" in e.detail
               for e in SE.recent(50)), [e.detail for e in SE.recent(50)]
    monkeypatch.setattr(EP, "_DIRECT_TCP", {})
    monkeypatch.setattr(EP, "_settings", lambda: {
        "allowed_domains": ["scp.example.test"]})
    direct = set(EP.direct_tcp_endpoints(repo))
    assert "scp.example.test:22" in direct, direct
    assert not any(d.startswith("attacker.example.test") for d in direct)


def test_a_plain_http_write_to_a_remote_host_is_not_passed_as_known(repo):
    class _Perms:
        workspace = repo
        mode = "default"
        confirm_callback = staticmethod(lambda *a: False)

    proxy = EP.EgressProxy.__new__(EP.EgressProxy)
    proxy.allowed_domains = tuple(EP.DEFAULT_ALLOWED_DOMAINS)
    proxy.allowed_ports = (80, 443)
    proxy.upload_ask_bytes = 1 << 20
    proxy.upload_domains = ()
    proxy._upload_from_settings = False
    verdict, _why, upload = proxy.decide("plain.example.test", 80,
                                         method="POST", body_bytes=10,
                                         perms=_Perms())
    assert verdict == "ask" and upload, (verdict, upload)
    verdict, _why, _ = proxy.decide("git.example.test", 443,
                                    perms=_Perms())
    assert verdict == "allow"


def test_the_remotes_are_read_without_the_repositorys_programs(repo,
                                                              monkeypatch):
    from delfin.agent import safe_git
    calls = []
    real = safe_git.run

    def _spy(args, cwd, **kw):
        calls.append(list(args))
        return real(args, cwd, **kw)
    monkeypatch.setattr(safe_git, "run", _spy)
    assert "https://git.example.test/team/repo.git" in EP.git_remotes(repo)
    assert calls, "git_remotes did not go through safe_git"
