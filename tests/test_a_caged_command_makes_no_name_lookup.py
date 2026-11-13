"""A caged command reaches no host-name resolver unless it needs one.

Input: ``socket_guard.name_lookup_sockets`` for a command's routes and
network mode and the host's nsswitch configuration, and the Unix-socket
allow-list ``landlock_exec.main`` hands the socket guard.

Semantics under test:

* systemd-resolved's sockets are not on the list every caged command
  gets: they answer host-name lookups only, and a lookup carries its
  name to the DNS (``<data>.attacker.example`` reaches the attacker's
  name server). A command behind the proxy needs no lookup, and its
  inet datagrams are refused already;
* they are offered to a command with a direct route (SLURM, a license
  server) only on a host whose ``hosts:`` line names ``resolve``, where
  that program could not resolve its route otherwise, and to an open
  network, which can send DNS anyway;
* the user-name services (sssd, nscd) and munge stay on the list.

Measured before the change (probe, bubblewrap with the socket guard,
bypass, unattended): a caged command connected to
/run/systemd/resolve/io.systemd.Resolve; UDP to any resolver was already
refused.
"""
from __future__ import annotations

from delfin.agent import landlock_exec as LE
from delfin.agent import socket_guard as SG


def _nss(tmp_path, hosts_line):
    p = tmp_path / "nsswitch.conf"
    p.write_text(f"passwd: files systemd\n{hosts_line}\n")
    return str(p)


def test_the_resolver_is_not_on_every_commands_list():
    assert not any(p.startswith("/run/systemd/resolve") for p in SG.DEFAULT_ALLOWED)
    for kept in ("/run/munge", "/var/lib/sss/pipes", "/run/nscd"):
        assert kept in SG.DEFAULT_ALLOWED


def test_a_command_behind_the_proxy_gets_no_resolver(tmp_path):
    nss = _nss(tmp_path, "hosts: files resolve [!UNAVAIL=return] dns")
    assert SG.name_lookup_sockets(False, "proxy", nss) == ()
    assert SG.name_lookup_sockets(False, "none", nss) == ()


def test_a_direct_route_gets_it_only_where_names_resolve_through_it(tmp_path):
    assert SG.name_lookup_sockets(
        True, "proxy", _nss(tmp_path, "hosts: files dns")) == ()
    assert SG.name_lookup_sockets(
        True, "proxy", _nss(tmp_path, "hosts:  files resolve [!UNAVAIL=return] dns")
    ) == SG.NAME_LOOKUP_SOCKETS
    assert SG.name_lookup_sockets(
        True, "proxy", _nss(tmp_path, "# hosts: resolve\nhosts: files dns")) == ()
    assert SG.name_lookup_sockets(True, "proxy", str(tmp_path / "missing")) == ()


def test_an_open_network_keeps_it(tmp_path):
    assert SG.name_lookup_sockets(
        False, "open", _nss(tmp_path, "hosts: files dns")) == SG.NAME_LOOKUP_SOCKETS


def test_the_helper_hands_the_guard_the_narrowed_list(tmp_path, monkeypatch):
    got = {}

    def fake_run(write, hide, allow_socket, command, strict, **kw):
        got["allow"] = list(allow_socket)
        return 0

    monkeypatch.setattr(LE, "_run_guarded", fake_run)
    ws = str(tmp_path)
    assert LE.main(["--write", ws, "--net", "proxy", "--", "true"]) == 0
    assert not any(p.startswith("/run/systemd/resolve") for p in got["allow"])
    assert "/run/munge" in got["allow"] and ws in got["allow"]
    assert LE.main(["--write", ws, "--net", "open", "--", "true"]) == 0
    assert "/run/systemd/resolve" in got["allow"]
