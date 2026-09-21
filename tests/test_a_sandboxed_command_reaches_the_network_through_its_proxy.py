"""An isolated command reaches the network only through DELFIN's proxy.

Modelled on the sandbox of a coding agent: filesystem limits alone let a
command in an unattended run post a workspace file to any host, reach an
internal service or ask the cloud metadata endpoint for credentials. The
command now gets a proxy that lets through the allowed domains; the socket
guard (or a network namespace) makes it the only way out, and refuses the
metadata address in every mode. Driven on 2026-09-16 through the agent's
bash tool under bubblewrap, under Landlock and uncaged: PyPI 200 through
the proxy, example.com refused, UDP refused, metadata refused.
"""
import http.server
import os
import socket
import subprocess
import sys
import threading
import urllib.error
import urllib.request

import pytest

from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP
from delfin.agent import socket_guard as SG

HELPER = str(__import__("pathlib").Path(A.__file__).with_name("landlock_exec.py"))
guard = pytest.mark.skipif(not SG.available(), reason="needs seccomp user notification")


@pytest.fixture
def web():
    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            body = f"hello {self.headers.get('Host')}".encode()
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *a):
            pass

    srv = http.server.HTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield srv.server_address[1]
    srv.shutdown()


@pytest.fixture
def proxy(web):
    refused = []
    p = EP.EgressProxy(allowed_domains=["allowed.test", "*.corp.test"],
                       allowed_ports=[web],
                       connector=lambda h, port, t: socket.create_connection(("127.0.0.1", port), timeout=t),
                       on_refusal=refused.append)
    p.refused = refused
    yield p
    p.close()


def _through(proxy, url, token=True):
    auth = f"delfin:{proxy.token}@" if token else ""
    fwd = socket.socket(socket.AF_INET)
    fwd.bind(("127.0.0.1", 0))
    fwd.listen(4)
    port = fwd.getsockname()[1]

    def pipe():
        conn, _ = fwd.accept()
        up = socket.socket(socket.AF_UNIX)
        up.connect(proxy.socket_path)
        import select
        while True:
            r, _, _ = select.select([conn, up], [], [], 10)
            if not r:
                return
            for s in r:
                d = s.recv(65536)
                if not d:
                    return
                (up if s is conn else conn).sendall(d)

    threading.Thread(target=pipe, daemon=True).start()
    opener = urllib.request.build_opener(urllib.request.ProxyHandler(
        {"http": f"http://{auth}127.0.0.1:{port}"}))
    try:
        return opener.open(url, timeout=10).read().decode()
    except urllib.error.HTTPError as exc:
        return f"{exc.code} {exc.read().decode()}"


def test_the_proxy_lets_through_allowed_domains_only(proxy, web):
    assert _through(proxy, f"http://allowed.test:{web}/").startswith("hello allowed.test")
    assert _through(proxy, f"http://api.corp.test:{web}/").startswith("hello api.corp.test")
    out = _through(proxy, f"http://evil.test:{web}/")
    assert out.startswith("403") and "allowed_domains" in out
    assert proxy.refused


def test_the_proxy_wants_its_token(proxy, web):
    assert _through(proxy, f"http://allowed.test:{web}/", token=False).startswith("407")


def test_the_proxy_itself_refuses_private_addresses(monkeypatch):
    p = EP.EgressProxy(allowed_domains=["internal.test"], allowed_ports=[80])
    try:
        monkeypatch.setattr(socket, "getaddrinfo", lambda h, port, *a, **k: [
            (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("10.0.0.5", port))])
        with pytest.raises(OSError, match="refused"):
            p._connector("internal.test", 80, 2.0)
    finally:
        p.close()


def test_domain_patterns():
    assert EP.domain_allowed("files.pythonhosted.org", EP.DEFAULT_ALLOWED_DOMAINS)
    assert EP.domain_allowed("a.b.corp.test", ["*.corp.test"])
    assert not EP.domain_allowed("corp.test", ["*.corp.test"])
    assert not EP.domain_allowed("evilpypi.org", ["pypi.org"])


_PROBE = r'''
import os, socket, sys, threading, urllib.request, urllib.error
port = int(sys.argv[1])
def get(url):
    try:
        return urllib.request.urlopen(url, timeout=10).read().decode()[:20]
    except urllib.error.HTTPError as e:
        return f"HTTP{e.code}"
    except Exception as e:
        return "error"
def tcp(host, p, timeout=3):
    s = socket.socket(); s.settimeout(timeout)
    try:
        s.connect((host, p)); return "connected"
    except OSError as e:
        return "refused" if e.errno == 13 else "error"
def udp():
    try:
        socket.socket(socket.AF_INET, socket.SOCK_DGRAM); return "created"
    except OSError:
        return "refused"
res = {}
t = threading.Thread(target=lambda: res.update(thread=tcp("127.0.0.1", port)))
t.start(); t.join()
print(get(f"http://allowed.test:{port}/"), get(f"http://evil.test:{port}/"),
      tcp("127.0.0.1", port), res["thread"], udp(), tcp("169.254.169.254", 80))
'''


def _helper(proxy, web, mode, tmp_path):
    probe = tmp_path / "probe.py"
    probe.write_text(_PROBE)
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    return subprocess.run(
        [sys.executable, "-I", HELPER, "--write", str(ws), "--net", mode,
         "--proxy-socket", proxy.socket_path, "--strict", "1", "--",
         sys.executable, str(probe), str(web)],
        capture_output=True, text=True, timeout=120)


@guard
def test_proxy_mode_is_the_only_way_out(proxy, web, tmp_path):
    if not A._landlock_functional():
        pytest.skip("needs Landlock")
    out = _helper(proxy, web, "proxy", tmp_path)
    fields = out.stdout.split()
    assert fields[0] == "hello" and fields[1].startswith("allowed.test"), out
    assert fields[2:] == ["HTTP403", "refused", "refused", "refused", "refused"], out


@guard
def test_open_mode_still_refuses_the_metadata_address(proxy, web, tmp_path):
    if not A._landlock_functional():
        pytest.skip("needs Landlock")
    out = _helper(proxy, web, "open", tmp_path)
    fields = out.stdout.split()
    assert fields[-4:] == ["connected", "connected", "created", "refused"], out


@guard
def test_a_connect_from_a_thread_is_answered(proxy, web, tmp_path):
    """The notification names the THREAD; the socket is the process's."""
    if not A._landlock_functional():
        pytest.skip("needs Landlock")
    out = _helper(proxy, web, "open", tmp_path)
    assert out.stdout.split()[-3] == "connected", out


def test_an_uncaged_attended_command_gets_the_guard(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_process_cage_functional", lambda: False)
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    monkeypatch.setattr(SG, "available", lambda: True)
    monkeypatch.delenv(A._PROCESS_CAGE_ENV, raising=False)
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="acceptEdits")
    argv = A._bash_isolation_argv("ls", tmp_path, perms, mode="auto")
    assert argv[2].endswith("landlock_exec.py") and "--deny-socket" in argv
    assert argv[argv.index("--net") + 1] == "open"
    monkeypatch.setenv(A._PROCESS_CAGE_ENV, "off")
    assert A._bash_isolation_argv("ls", tmp_path, perms, mode="auto") == ["/bin/bash", "-c", "ls"]


def test_an_isolated_command_under_bubblewrap_runs_the_guard_inside(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    monkeypatch.setattr(A.shutil, "which", lambda _x: "/usr/bin/bwrap")
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    monkeypatch.setattr(SG, "available", lambda: True)
    perms = A.KitToolPermissions(workspace=str(tmp_path), lock_workspace=True)
    argv = A._bash_isolation_argv("ls", tmp_path, perms)
    i = argv.index(HELPER)
    assert argv[i + 1:i + 5] == ["--fs", "0", "--net", "proxy"]
    proxy_dir = os.path.dirname(argv[argv.index("--proxy-socket") + 1])
    assert argv[argv.index(proxy_dir) - 1] == "--ro-bind"
    monkeypatch.setattr(SG, "available", lambda: False)
    argv = A._bash_isolation_argv("ls", tmp_path, perms)
    assert "--unshare-net" in argv and HELPER not in argv


def test_the_slurm_controller_is_reachable_from_the_sandbox(tmp_path, monkeypatch):
    conf = tmp_path / "slurm.conf"
    conf.write_text("# cluster\nSlurmctldHost=localhost(127.0.0.1)\nSlurmctldPort=6817\n")
    monkeypatch.setenv("SLURM_CONF", str(conf))
    assert "127.0.0.1:6817" in EP.slurm_controllers()


def test_the_network_mode_comes_from_the_settings(monkeypatch):
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "none"})
    assert A._sandbox_network_args() == (["--net", "none"], "")
