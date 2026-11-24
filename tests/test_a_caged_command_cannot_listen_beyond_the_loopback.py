"""A caged command cannot listen beyond the loopback.

Measured 2026-11-24 on a shared server without a host firewall: an
unattended agent ran a parallel ORCA (``%pal nprocs 2``) through its bash
tool. OpenMPI's TCP transport in the ORCA ranks listened on
0.0.0.0:1024-1025 and [::]:1024-1025 (``btl_tcp_port_min_v4`` is 1024),
and mpirun on four more wildcard ports, for three minutes, until the
swarm's own watch stopped every agent. The cage held the command's
OUTBOUND traffic to the egress proxy -- the socket guard judged every
connect() -- and nothing judged a bind() or a listen(), so the sockets
were reachable from the network. Environment variables do not close it:
with ``oob_tcp_if_include``, ``btl_tcp_if_include`` and PMIx's
``ptl_tcp_if_include`` all set to the loopback, mpirun still bound the
wildcard (measured in a network namespace of its own, where nothing is
reachable from outside).

Since then the supervisor performs every bind() and listen() itself: a
wildcard becomes the loopback of its family, any other address but the
loopback is refused, a listen() on an unbound socket binds it to the
loopback first, and a datagram from an unbound socket under a restricted
network leaves from a loopback-bound one. Only an open network inside a
SLURM allocation of several nodes keeps its binds (``--inbound any``).

The rules are checked on every machine; the tests that need the kernel's
seccomp user notification run the command under the real guard, and one
of them through the argv the bash tool builds. The live probe checks a
bind WITHOUT listen first and stops there unless it landed on the
loopback, so a broken guard never leaves a listening wildcard socket
behind, not even for the test's lifetime.
"""
from __future__ import annotations

import ipaddress
import json
import socket
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import child_env

import _egress_test_support as ES
from delfin.agent import api_client as A
from delfin.agent import landlock_exec as LL
from delfin.agent import socket_guard as SG

HELPER = ES.HELPER


def _held_by_a_namespace_instead(tmp_path, monkeypatch):
    """On a kernel without seccomp user notification the socket guard does
    not run; what holds the way in there is checked instead of skipping:
    the bash tool's cage puts the command in a network namespace of its
    own where bubblewrap works, and refuses an unattended command where
    nothing can hold it."""
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    argv = A._bash_isolation_argv("true", perms.workspace, perms)
    if A._bwrap_functional():
        assert "--unshare-net" in argv, argv
    elif not A._landlock_functional():
        assert argv[:2] == ["/bin/sh", "-c"] and "refused" in argv[-1], argv


# -- the rule, on every machine ------------------------------------------------

@pytest.mark.parametrize("asked, bound", [
    ("127.0.0.1", "127.0.0.1"),
    ("127.8.9.10", "127.8.9.10"),
    ("0.0.0.0", "127.0.0.1"),
    ("::1", "::1"),
    ("::", "::1"),
    ("::ffff:0.0.0.0", "::ffff:127.0.0.1"),
    ("::ffff:127.0.0.1", "::ffff:127.0.0.1"),
    ("10.0.0.1", None),
    ("192.168.1.5", None),
    ("198.51.100.7", None),
    ("fe80::1", None),
    ("2001:db8::1", None),
    ("::ffff:10.0.0.1", None),
])
def test_a_bind_lands_on_the_loopback_or_is_refused(asked, bound):
    target = SG.Supervisor(0, []).bind_target(ipaddress.ip_address(asked))
    assert (str(target) if target is not None else None) == bound


def test_a_multi_node_allocation_with_an_open_network_binds_as_asked():
    any_ = SG.Supervisor(0, [], inbound="any")
    for asked in ("0.0.0.0", "10.0.0.1", "::"):
        ip = ipaddress.ip_address(asked)
        assert any_.bind_target(ip) == ip
    # anything but "any" is the loopback rule
    assert SG.Supervisor(0, [], inbound="nonsense").inbound == "loopback"


@pytest.mark.parametrize("env, multi", [
    ({}, False),
    ({"SLURM_JOB_NUM_NODES": "4"}, False),
    ({"SLURM_JOB_ID": "7", "SLURM_JOB_NUM_NODES": "1"}, False),
    ({"SLURM_JOB_ID": "7", "SLURM_JOB_NUM_NODES": "2"}, True),
    ({"SLURM_JOB_ID": "7", "SLURM_NNODES": "3"}, True),
    ({"SLURM_JOB_ID": "7", "SLURM_JOB_NUM_NODES": "two"}, False),
])
def test_a_multi_node_allocation_is_recognised(env, multi):
    assert A._in_multi_node_allocation(env) is multi


@pytest.mark.parametrize("mode, multi, inbound", [
    ("open", True, True),
    ("open", False, False),
    ("none", True, False),
])
def test_only_an_open_network_in_a_multi_node_allocation_keeps_its_binds(
        monkeypatch, mode, multi, inbound):
    from delfin.agent import egress_proxy as EP
    monkeypatch.setattr(EP, "network_mode", lambda: mode)
    monkeypatch.setattr(A, "_in_multi_node_allocation", lambda env=None: multi)
    args, _proxy_dir = A._sandbox_network_args(None)
    assert (["--inbound", "any"] == args[2:4]) is inbound, args


def test_a_command_nobody_saw_never_keeps_its_binds(monkeypatch):
    """A caged_unasked command gets no open network (it becomes the proxy),
    so no --inbound any either, in a multi-node allocation too."""
    from delfin.agent import egress_proxy as EP
    monkeypatch.setattr(EP, "network_mode", lambda: "open")
    monkeypatch.setattr(A, "_in_multi_node_allocation", lambda env=None: True)
    monkeypatch.setattr(EP, "get_proxy", lambda: (_ for _ in ()).throw(OSError("none here")))
    args, _proxy_dir = A._sandbox_network_args(None, allow_slurm=False)
    assert "--inbound" not in args, args


def test_the_helper_hands_the_inbound_mode_to_the_supervisor(monkeypatch):
    seen = {}

    def fake(*a, **kw):
        seen.update(kw)
        return 0

    monkeypatch.setattr(LL, "_run_guarded", fake)
    assert LL.main(["--fs", "0", "--inbound", "any", "--", "true"]) == 0
    assert seen["inbound"] == "any"
    seen.clear()
    assert LL.main(["--fs", "0", "--", "true"]) == 0
    assert seen["inbound"] == "loopback"


def test_bind_and_listen_go_to_the_supervisor_in_every_network_mode():
    nr_bind, nr_listen = SG._arch_dgram()[0], SG._arch_listen()
    for restricted in (False, True):
        p = SG._program(block_inet_datagrams=restricted)
        for nr in (nr_bind, nr_listen):
            i = next(i for i, ins in enumerate(p) if ins[0] == SG._JEQ and ins[3] == nr)
            assert p[i + 1 + p[i][1]] == (SG._RET, 0, 0, SG._RET_USER_NOTIF)


# -- the guard, live ------------------------------------------------------------

_PROBE = r'''
import errno, json, os, socket, sys
args = json.loads(sys.argv[1])
out = {}
def err(e):
    return errno.errorcode.get(e.errno, str(e.errno))
# First, and alone: a wildcard bind WITHOUT listen. A broken guard leaves a
# bound socket nobody can connect to, and the probe stops before any listen.
s = socket.socket(); s.bind(("0.0.0.0", 0)); out["tcp 0.0.0.0"] = s.getsockname()[0]; s.close()
if out["tcp 0.0.0.0"] != "127.0.0.1":
    print(json.dumps(out)); sys.exit(0)
def name(fam, addr, listen):
    s = socket.socket(fam, socket.SOCK_STREAM)
    try:
        s.bind(addr)
        if listen:
            s.listen(1)
        return s.getsockname()[0]
    except OSError as e:
        return err(e)
    finally:
        s.close()
out["tcp 0.0.0.0 listen"] = name(socket.AF_INET, ("0.0.0.0", 0), True)
out["tcp 127.0.0.1 listen"] = name(socket.AF_INET, ("127.0.0.1", 0), True)
if args.get("external"):
    out["tcp external"] = name(socket.AF_INET, (args["external"], 0), False)
if socket.has_ipv6 and args.get("ipv6"):
    out["tcp6 ::"] = name(socket.AF_INET6, ("::", 0, 0, 0), True)
s = socket.socket(); s.listen(1); out["listen unbound"] = s.getsockname()[0]; s.close()
u = socket.socket(socket.AF_INET, socket.SOCK_DGRAM); u.bind(("", 0)); out["udp ''"] = u.getsockname()[0]; u.close()
# a server a tool binds to "" still serves its own client on the loopback
srv = socket.socket(); srv.bind(("", 0)); srv.listen(1)
c = socket.create_connection(("127.0.0.1", srv.getsockname()[1]), timeout=5)
a, _ = srv.accept(); c.sendall(b"ok"); out["own server"] = a.recv(2).decode()
for x in (a, c, srv):
    x.close()
if args.get("unix_ok"):
    def unix(path):
        s = socket.socket(socket.AF_UNIX)
        try:
            s.bind(path); s.listen(1); os.unlink(path); return "bound"
        except OSError as e:
            return err(e)
        finally:
            s.close()
    out["unix allowed"] = unix(os.path.join(args["unix_ok"], "a.sock"))
    if args.get("unix_other"):
        out["unix elsewhere"] = unix(os.path.join(args["unix_other"], "b.sock"))
    os.chdir(args["deep"])
    s = socket.socket(socket.AF_UNIX); s.bind("r.sock"); s.listen(1)
    c = socket.socket(socket.AF_UNIX); c.connect("r.sock")
    a, _ = s.accept(); c.send(b"d"); out["unix deep relative"] = a.recv(1).decode()
    for x in (a, c, s):
        x.close()
    os.unlink("r.sock")
print(json.dumps(out))
'''


def _external_address():
    """An IPv4 address of this host beyond the loopback, or "". Read from
    the kernel's routing trie (the addresses it marks LOCAL), so the test
    opens no socket to find it."""
    try:
        lines = Path("/proc/net/fib_trie").read_text().splitlines()
    except OSError:
        return ""
    for prev, line in zip(lines, lines[1:]):
        if line.strip() == "/32 host LOCAL":
            ip = prev.strip().lstrip("|-").strip()
            try:
                if not ipaddress.ip_address(ip).is_loopback:
                    return ip
            except ValueError:
                continue
    return ""


def _ipv6_loopback():
    try:
        with socket.socket(socket.AF_INET6) as s:
            s.bind(("::1", 0))
        return True
    except OSError:
        return False


def _run(argv, ws, tmp_path):
    out = subprocess.run(argv, cwd=str(ws), env=child_env(tmp_path),
                         capture_output=True, text=True, timeout=240)
    assert out.returncode == 0, out.stderr[-2000:]
    return json.loads(out.stdout.strip().splitlines()[-1])


def _expected(args):
    want = {"tcp 0.0.0.0": "127.0.0.1", "tcp 0.0.0.0 listen": "127.0.0.1",
            "tcp 127.0.0.1 listen": "127.0.0.1", "listen unbound": "127.0.0.1",
            "udp ''": "127.0.0.1", "own server": "ok"}
    if args.get("external"):
        want["tcp external"] = "EACCES"
    if args.get("ipv6"):
        want["tcp6 ::"] = "::1"
    if args.get("unix_ok"):
        want["unix allowed"] = "bound"
        want["unix deep relative"] = "d"
        if args.get("unix_other"):
            want["unix elsewhere"] = "EACCES"
    return want


def _args(ws, tmp_path, with_unix_other=True):
    deep = ws / ("d" * 60) / ("e" * 60)
    deep.mkdir(parents=True)
    other = tmp_path / "other"
    other.mkdir()
    return {"external": _external_address(), "ipv6": _ipv6_loopback(),
            "unix_ok": str(ws), "deep": str(deep),
            "unix_other": str(other) if with_unix_other else ""}


@pytest.fixture
def loopback_web():
    yield from ES.serve_web()


@pytest.fixture
def loopback_proxy(loopback_web, monkeypatch):
    yield from ES.session_proxy(loopback_web, monkeypatch)


def _probe_under_the_helper(tmp_path, net_args):
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    ws = Path(perms.workspace)
    (ws / "probe.py").write_text(_PROBE)
    args = _args(ws, tmp_path)
    argv = [sys.executable, "-I", HELPER, "--fs", "0", *net_args(perms),
            "--allow-socket", str(ws), "--",
            sys.executable, str(ws / "probe.py"), json.dumps(args)]
    assert _run(argv, ws, tmp_path) == _expected(args)


def test_the_guard_holds_every_bind_to_the_loopback_behind_the_proxy(
        tmp_path, loopback_proxy, monkeypatch):
    """The default network, the measured run's."""
    if not SG.available():
        return _held_by_a_namespace_instead(tmp_path, monkeypatch)
    _probe_under_the_helper(tmp_path, lambda perms: [
        "--net", "proxy", "--proxy-socket", loopback_proxy.socket_path,
        "--proxy-token-file", loopback_proxy.token_file_for(perms)])


def test_the_guard_holds_every_bind_to_the_loopback_with_an_open_network(
        tmp_path, monkeypatch):
    """An open network opens the way out, not the way in."""
    if not SG.available():
        return _held_by_a_namespace_instead(tmp_path, monkeypatch)
    _probe_under_the_helper(tmp_path, lambda perms: ["--net", "open"])


def test_with_inbound_any_a_bind_is_carried_out_as_asked(tmp_path, monkeypatch):
    """The multi-node allocation's mode, checked without any listen: a
    wildcard bind stays the wildcard."""
    if not SG.available():
        return _held_by_a_namespace_instead(tmp_path, monkeypatch)
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    ws = Path(perms.workspace)
    probe = ("import socket; s = socket.socket(); s.bind(('0.0.0.0', 0)); "
             "print(s.getsockname()[0]); s.close()")
    out = subprocess.run([sys.executable, "-I", HELPER, "--fs", "0", "--net", "open",
                          "--inbound", "any", "--", sys.executable, "-c", probe],
                         cwd=str(ws), env=child_env(tmp_path), capture_output=True,
                         text=True, timeout=120)
    assert out.stdout.split() == ["0.0.0.0"], out.stderr[-1000:]


def test_the_bash_tools_cage_holds_every_bind_to_the_loopback(tmp_path, monkeypatch):
    """The same through the argv the bash tool builds, bubblewrap with the
    helper inside where it works -- the cage the measured run was in."""
    if not SG.available() or not (A._bwrap_functional() or A._landlock_functional()):
        return _held_by_a_namespace_instead(tmp_path, monkeypatch)
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    ws = Path(perms.workspace)
    (ws / "probe.py").write_text(_PROBE)
    # A directory outside the roots is behind the cage's fresh /tmp under
    # bubblewrap, so "elsewhere" is not asked here.
    args = _args(ws, tmp_path, with_unix_other=False)
    argv = A._bash_isolation_argv(
        f"{sys.executable} {ws / 'probe.py'} '{json.dumps(args)}'", ws, perms)
    assert _run(argv, ws, tmp_path) == _expected(args)
