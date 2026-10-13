"""An isolated command reaches the network only through DELFIN's proxy.

Modelled on the sandbox of a coding agent: filesystem limits alone let a
command in an unattended run post a workspace file to any host, reach an
internal service or ask the cloud metadata endpoint for credentials. The
command now gets a proxy that lets through the allowed domains; the socket
guard (or a network namespace) makes it the only way out, and refuses the
metadata address in every mode. Driven on 2026-09-16 through the agent's
bash tool under bubblewrap, under Landlock and uncaged: PyPI 200 through
the proxy, example.com refused, UDP refused, metadata refused.

The tests at the end start a child under each real kernel layer with a
per-session proxy (the session's permissions decide); the decisions of
that proxy are tested host-independently in
test_a_child_process_reaches_the_network_only_through_the_proxy.py.
"""
import http.server
import os
import shutil
import socket
import subprocess
import sys
import threading
import urllib.error
import urllib.request
from pathlib import Path

import pytest

import _egress_test_support as ES
from conftest import child_env
from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP
from delfin.agent import landlock_exec as LL
from delfin.agent import socket_guard as SG

HELPER = ES.HELPER
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


def test_an_attended_command_is_caged_and_reaches_the_net_through_the_proxy(
        tmp_path, monkeypatch):
    """An attended session used to run its commands through a plain shell
    with only a socket guard: no filesystem containment, and `--net open`.
    Since "auto" means "wherever this host can hold it", the same command
    is now held to the workspace and its network goes through the proxy.
    Measured before that change: `echo x > "$HOME"/f` in acceptEdits wrote
    into the home directory, because the gate reads the command TEXT and
    the path was not in it."""
    monkeypatch.setattr(A, "_process_cage_functional", lambda: False)
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    monkeypatch.setattr(SG, "available", lambda: True)
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_landlock_functional", lambda: True)
    monkeypatch.delenv(A._PROCESS_CAGE_ENV, raising=False)
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="acceptEdits")
    argv = A._bash_isolation_argv("ls", tmp_path, perms, mode="auto")
    assert any(str(a).endswith("landlock_exec.py") for a in argv), argv
    assert argv[argv.index("--net") + 1] == "proxy", argv

    # With the cage switched off by name and no filesystem layer, the
    # filesystem is as it is -- and the network is still held by the
    # socket guard (no --fs policy, --net proxy).
    monkeypatch.setattr(A, "_landlock_functional", lambda: False)
    monkeypatch.setenv(A._PROCESS_CAGE_ENV, "off")
    argv = A._bash_isolation_argv("ls", tmp_path, perms, mode="auto")
    assert argv[-4:] == ["--", "/bin/bash", "-c", "ls"], argv
    assert argv[argv.index("--fs") + 1] == "0", argv
    assert argv[argv.index("--net") + 1] == "proxy", argv

    # And it lets go entirely when the user set the network open.
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "open"})
    assert A._bash_isolation_argv(
        "ls", tmp_path, perms, mode="auto") == ["/bin/bash", "-c", "ls"]


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
    # Without the guard: a network namespace, and the helper inside it
    # forwarding to the proxy (it used to be no network at all).
    monkeypatch.setattr(SG, "available", lambda: False)
    argv = A._bash_isolation_argv("ls", tmp_path, perms)
    assert "--unshare-net" in argv and HELPER in argv
    assert argv[argv.index("--netns") + 1] == "1"


def test_the_slurm_controller_is_reachable_from_the_sandbox(tmp_path, monkeypatch):
    conf = tmp_path / "slurm.conf"
    conf.write_text("# cluster\nSlurmctldHost=localhost(127.0.0.1)\nSlurmctldPort=6817\n")
    monkeypatch.setenv("SLURM_CONF", str(conf))
    assert "127.0.0.1:6817" in EP.slurm_controllers()


def test_the_network_mode_comes_from_the_settings(monkeypatch):
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "none"})
    assert A._sandbox_network_args() == (["--net", "none"], "")


# -- a child under a session's proxy cannot go around it ---------------------
# session_web serves a page, an upload sink and a package index;
# session_proxy is the process-wide proxy, deciding by the session's
# permissions (allowed.test allowed, everything else refused unasked).

@pytest.fixture
def session_web():
    yield from ES.serve_web()


@pytest.fixture
def session_proxy(session_web, monkeypatch):
    yield from ES.session_proxy(session_web, monkeypatch)


_ROUTE_PROBE = r'''
import socket, sys, urllib.request, urllib.error
port = int(sys.argv[1])
def get(url, data=None):
    try:
        return urllib.request.urlopen(urllib.request.Request(url, data=data),
                                      timeout=10).read().decode().split()[0]
    except urllib.error.HTTPError as e:
        return f"HTTP{e.code}"
    except Exception as e:
        return "error"
def tcp(host, p):
    s = socket.socket(); s.settimeout(3)
    try:
        s.connect((host, p)); return "connected"
    except OSError as e:
        return "refused"
print(tcp("127.0.0.1", port), get(f"http://allowed.test:{port}/"),
      get(f"http://evil.test:{port}/"))
'''


def _run(argv, ws, tmp_path, extra_env=None):
    env = child_env(tmp_path)
    env.update(extra_env or {})
    return subprocess.run(argv, cwd=str(ws), env=env, capture_output=True,
                          text=True, timeout=240)


def _helper_argv(proxy, perms, ws, extra):
    return [sys.executable, "-I", HELPER, "--fs", "0", "--net", "proxy",
            "--proxy-socket", proxy.socket_path,
            "--proxy-token-file", proxy.token_file_for(perms),
            "--allow-socket", "/", *extra, "--"]


@guard
@pytest.mark.parametrize("mode", ["default", "acceptEdits", "plan", "bypassPermissions"])
def test_a_child_python_cannot_connect_directly_in_any_mode(
        session_proxy, session_web, tmp_path, monkeypatch, mode):
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    if not (A._bwrap_functional() or A._landlock_functional()):
        pytest.skip("needs bubblewrap or Landlock for the default isolation")
    perms = ES.perms_in(tmp_path, mode=mode)
    probe = Path(perms.workspace) / "probe.py"
    probe.write_text(_ROUTE_PROBE)
    argv = A._bash_isolation_argv(
        f"{sys.executable} {probe} {session_web.server_address[1]}",
        perms.workspace, perms)
    out = _run(argv, perms.workspace, tmp_path)
    assert out.stdout.split() == ["refused", "hello", "HTTP403"], out


@pytest.mark.skipif(LL.abi_version() < 4, reason="needs Landlock ABI 4 (TCP connect)")
def test_landlock_tcp_holds_the_child_to_the_proxy_without_the_socket_guard(
        session_proxy, session_web, tmp_path):
    perms = ES.perms_in(tmp_path)
    ws = Path(perms.workspace)
    (ws / "probe.py").write_text(_ROUTE_PROBE)
    argv = _helper_argv(session_proxy, perms, ws, ["--guard", "0"]) + [
        sys.executable, str(ws / "probe.py"), str(session_web.server_address[1])]
    out = _run(argv, ws, tmp_path)
    assert out.stdout.split() == ["refused", "hello", "HTTP403"], out


@pytest.mark.skipif(shutil.which("bwrap") is None or not A._bwrap_functional(),
                    reason="needs a working bubblewrap")
def test_a_network_namespace_reaches_the_proxy_without_the_socket_guard(
        session_proxy, session_web, tmp_path):
    """bubblewrap without seccomp user notification used to mean no network
    at all (--unshare-net, proxy unreachable): pip and git broke. The
    forwarder now runs inside the namespace."""
    perms = ES.perms_in(tmp_path)
    ws = Path(perms.workspace)
    (ws / "probe.py").write_text(_ROUTE_PROBE)
    proxy_dir = os.path.dirname(session_proxy.socket_path)
    argv = (["bwrap", "--dev-bind", "/", "/", "--proc", "/proc", "--unshare-net",
             "--ro-bind", proxy_dir, proxy_dir, "--die-with-parent"]
            + _helper_argv(session_proxy, perms, ws, ["--guard", "0", "--netns", "1"])
            + [sys.executable, str(ws / "probe.py"), str(session_web.server_address[1])])
    out = _run(argv, ws, tmp_path)
    assert out.stdout.split() == ["refused", "hello", "HTTP403"], out


# -- a normal unattended run keeps working ----------------------------------

@guard
def test_an_unattended_run_installs_fetches_and_calls_its_model(
        session_proxy, session_web, tmp_path, monkeypatch):
    if not (A._bwrap_functional() or A._landlock_functional()):
        pytest.skip("needs bubblewrap or Landlock for the default isolation")
    monkeypatch.setattr(A, "_record_security_event", lambda *a, **k: None)
    port = session_web.server_address[1]
    index = f"http://pypi.test:{port}/simple"
    monkeypatch.setenv("PIP_INDEX_URL", index)
    EP.note_endpoint(f"http://models.test:{port}/v1")
    perms = ES.perms_in(tmp_path, mode="bypassPermissions")
    ws = Path(perms.workspace)

    git = shutil.which("git")
    remote = ws / "remote.git"
    clone = ws / "clone"
    if git:
        subprocess.run([git, "init", "-q", "--bare", str(remote)], check=True,
                       env=child_env(tmp_path))
        seed = ws / "seed"
        subprocess.run([git, "init", "-q", str(seed)], check=True, env=child_env(tmp_path))
        (seed / "f").write_text("x")
        for cmd in (["add", "f"], ["-c", "user.name=t", "-c", "user.email=t@t",
                                   "commit", "-qm", "m"],
                    ["push", "-q", str(remote), "HEAD:refs/heads/main"]):
            subprocess.run([git, "-C", str(seed), *cmd], check=True, env=child_env(tmp_path))
        subprocess.run([git, "clone", "-q", str(remote), str(clone)], check=True,
                       env=child_env(tmp_path))

    model = ws / "model.py"
    model.write_text(
        "import sys, urllib.request\n"
        "r = urllib.request.urlopen(urllib.request.Request(sys.argv[1], "
        "data=b'{\"messages\": []}'), timeout=10)\n"
        "print('model', r.read().decode())\n")
    cmd = (f"{sys.executable} -m pip install -q --no-cache-dir --target {ws}/site "
           f"tinypkg && {sys.executable} -c 'import sys; sys.path.insert(0, \"{ws}/site\"); "
           f"import tinypkg; print(\"pkg\", tinypkg.VALUE)' && "
           f"{sys.executable} {model} http://models.test:{port}/v1/chat")
    if git:
        cmd += f" && {git} -C {clone} fetch -q origin && echo fetched"
    argv = A._bash_isolation_argv(cmd, ws, perms)
    out = _run(argv, ws, tmp_path, {"PIP_INDEX_URL": index,
                                    "PIP_TRUSTED_HOST": "pypi.test",
                                    "PIP_DISABLE_PIP_VERSION_CHECK": "1"})
    assert out.returncode == 0, (out.stderr, session_web.seen, session_proxy.refused)
    assert "pkg 42" in out.stdout and "model got 16" in out.stdout, out
    if git:
        assert "fetched" in out.stdout, out
