"""The suite-wide guard of tests/external_agent_guard.py fires.

Each in-process case builds a stand-in "host" agent CLI in tmp_path and
registers it with the guard as a real binary, so the cases hold on a
machine with no agent CLI installed (CI) as on one with it. The case
drains what it caused from the guard's record before it ends, because
the conftest fixture would otherwise fail it -- the last test runs a
child pytest with the real conftest to show that this failure happens.
"""

from __future__ import annotations

import os
import shutil
import socket
import subprocess
import sys
import textwrap
import threading
import urllib.error
import urllib.request
from pathlib import Path

import pytest

import external_agent_guard as guard
from conftest import child_env

TESTS = Path(__file__).resolve().parent


@pytest.fixture()
def host_cli(tmp_path, monkeypatch):
    """A stand-in host ``claude`` on PATH, registered as a real binary."""
    bindir = tmp_path / "hostbin"
    bindir.mkdir()
    exe = bindir / "claude"
    exe.write_text("#!/bin/sh\necho REAL-CLI-RAN\n")
    exe.chmod(0o755)
    monkeypatch.setitem(guard._STATE, "real",
                        guard._STATE["real"] | {str(exe.resolve())})
    monkeypatch.setenv("PATH", os.environ["PATH"] + os.pathsep + str(bindir))
    yield exe
    guard.take_trips()


def test_the_guard_is_installed_for_the_session():
    assert guard._STATE["installed"]
    assert shutil.which is guard._which
    assert subprocess.Popen._execute_child is guard._execute_child
    assert socket.socket.connect is guard._connect


def test_discovery_by_name_finds_no_agent_cli(host_cli):
    assert shutil.which("claude") is None
    assert shutil.which("codex") is None
    # By explicit path: refused and recorded.
    assert shutil.which(str(host_cli)) is None
    assert any(str(host_cli) in t for t in guard.take_trips())


def test_no_session_variable_of_a_calling_agent_or_key_is_left():
    left = [k for k in os.environ
            if k.upper().startswith(guard.STRIPPED_PREFIXES)]
    assert left == []


def test_spawning_the_cli_by_path_is_refused(host_cli):
    with pytest.raises(guard.ExternalAgentTripwire):
        subprocess.run([str(host_cli), "-p", "hello"], capture_output=True)
    assert any("hello" in t for t in guard.take_trips())


def test_spawning_through_a_symlink_or_an_interpreter_is_refused(host_cli, tmp_path):
    link = tmp_path / "agent"
    link.symlink_to(host_cli)
    with pytest.raises(guard.ExternalAgentTripwire):
        subprocess.Popen([str(link)])
    with pytest.raises(guard.ExternalAgentTripwire):
        subprocess.Popen(["/bin/sh", str(host_cli)])
    assert len(guard.take_trips()) == 2


def test_spawning_the_cli_by_name_is_refused():
    with pytest.raises(guard.ExternalAgentTripwire):
        subprocess.run(["codex", "exec", "x"], capture_output=True)
    assert guard.take_trips()


def test_a_shell_that_resolves_the_name_meets_the_tripwire(host_cli):
    r = subprocess.run("claude -p from-a-shell", shell=True,
                       capture_output=True, text=True)
    assert r.returncode == 127
    assert "REAL-CLI-RAN" not in r.stdout
    assert "tripwire" in r.stderr
    trips = guard.take_trips()
    assert any("from-a-shell" in t for t in trips)


def test_os_exec_and_posix_spawn_are_refused(host_cli):
    with pytest.raises(guard.ExternalAgentTripwire):
        os.posix_spawn(str(host_cli), [str(host_cli)], dict(os.environ))
    assert len(guard.take_trips()) == 1
    # os.execv would replace the test process if the guard failed, so it is
    # checked by identity, not by a call.
    assert os.execv.__module__ == guard.__name__
    assert os.execve.__module__ == guard.__name__


def test_a_remote_endpoint_is_refused():
    # 192.0.2.0/24 is TEST-NET-1: refused before any packet leaves.
    with pytest.raises(guard.ExternalAgentTripwire):
        socket.create_connection(("192.0.2.1", 443), timeout=1)
    with pytest.raises(urllib.error.URLError):
        urllib.request.urlopen("http://192.0.2.1/v1/messages", timeout=1)
    s = socket.socket(socket.AF_INET6, socket.SOCK_STREAM)
    try:
        with pytest.raises(guard.ExternalAgentTripwire):
            s.connect_ex(("2001:db8::1", 443))
    finally:
        s.close()
    trips = guard.take_trips()
    assert len(trips) == 3 and all("192.0.2.1" in t or "2001:db8" in t
                                   for t in trips)


def test_a_model_client_meets_the_tripwire_not_the_endpoint():
    """The OpenAI-compatible client, given a key, still cannot reach a
    remote base URL."""
    from delfin.agent import api_client as A
    try:
        client = A.OpenAIClient(api_key="sk-test-not-a-key",
                                model="gpt-4o-mini",
                                base_url="http://192.0.2.1:9/v1")
    except TypeError:
        pytest.skip("OpenAIClient signature changed")
    with pytest.raises(Exception):
        list(client.stream_message("s", [{"role": "user", "content": "hi"}],
                                   max_tokens=8))
    assert any("192.0.2.1" in t for t in guard.take_trips())


def test_loopback_and_unix_sockets_are_untouched(tmp_path):
    srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    srv.bind(("127.0.0.1", 0))
    srv.listen(1)
    port = srv.getsockname()[1]
    t = threading.Thread(target=lambda: srv.accept()[0].close())
    t.start()
    socket.create_connection(("localhost", port), timeout=2).close()
    t.join(2)
    srv.close()
    path = str(tmp_path / "s")
    u = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    u.bind(path)
    u.listen(1)
    c = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    c.connect(path)
    c.close()
    u.close()
    assert guard.take_trips() == []


@pytest.mark.fake_agent_cli
def test_the_marker_gives_a_stub_and_still_refuses_the_real_one(host_cli):
    found = shutil.which("claude")
    assert found and Path(found).read_text().count("calls.log")
    subprocess.run(["claude", "-p", "stubbed"], check=True)
    assert "stubbed" in (Path(found).parent / "calls.log").read_text()
    with pytest.raises(guard.ExternalAgentTripwire):
        subprocess.run([str(host_cli)])
    with pytest.raises(guard.ExternalAgentTripwire):
        socket.create_connection(("192.0.2.1", 443), timeout=1)
    assert len(guard.take_trips()) == 2


def test_the_product_sees_no_cli_and_refuses_to_build_its_client():
    from delfin.agent import api_client as A
    with pytest.raises(FileNotFoundError):
        A.CLIClient(model="sonnet")
    with pytest.raises(FileNotFoundError):
        A.CodexCLIClient(model="gpt-5")
    assert guard.take_trips() == []


def test_a_test_that_swallows_the_refusal_still_fails(tmp_path):
    """The conftest fixture turns a recorded refusal into a failure even
    when the code under test caught the exception."""
    inner = tmp_path / "inner"
    inner.mkdir()
    (inner / "conftest.py").write_text(textwrap.dedent(f"""\
        import importlib.util, sys
        sys.path.insert(0, {str(TESTS)!r})
        _spec = importlib.util.spec_from_file_location(
            "delfin_suite_conftest", {str(TESTS / 'conftest.py')!r})
        _c = importlib.util.module_from_spec(_spec)
        sys.modules["delfin_suite_conftest"] = _c
        _spec.loader.exec_module(_c)
        globals().update({{k: getattr(_c, k) for k in dir(_c)
                          if not k.startswith("__")}})
        """))
    (inner / "test_inner.py").write_text(textwrap.dedent("""\
        import socket, subprocess
        import pytest

        def test_swallowed_spawn():
            try:
                subprocess.run(["claude", "-p", "x"])
            except Exception:
                pass

        def test_swallowed_connect():
            try:
                socket.create_connection(("192.0.2.1", 443), timeout=1)
            except Exception:
                pass

        def test_clean():
            assert True

        @pytest.mark.fake_agent_cli
        def test_stub_allowed():
            subprocess.run(["claude", "-p", "x"], check=True)
        """))
    env = child_env(tmp_path)
    # The child imports this checkout's delfin, not whatever an editable
    # install of the interpreter points at.
    env["PYTHONPATH"] = os.pathsep.join(
        [str(TESTS.parent)] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    r = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
         "-p", "no:randomly", "-o", "addopts=", "--rootdir", str(inner),
         str(inner / "test_inner.py")],
        cwd=str(inner), env=env, capture_output=True, text=True, timeout=300)
    out = r.stdout + r.stderr
    assert r.returncode != 0, out
    # the call phases pass; the two that reached out error at teardown
    assert "4 passed, 2 errors" in out, out
    assert "ERROR test_inner.py::test_swallowed_spawn" in out, out
    assert "ERROR test_inner.py::test_swallowed_connect" in out, out
    assert "external agent tripwire" in out, out
    assert "spawned an agent CLI by name: claude -p x" in out, out
    assert "network connection to 192.0.2.1:443" in out, out
