"""The guard refuses a remote destination reached through urllib.

A CHARACTERISATION test, not proof of the proxy fix: it passes on the
commit before the fix as well. It is here because the two cases it covers
had no test of their own, and because the shape of the case is what the
proxy blind spot needed (see below).

Measured 2026-11-26, in a swarm member's own log: 7,846 refusals in two
hours, every one a test process reaching the cage's egress proxy. The
member's shell carried ``https_proxy=http://127.0.0.1:<cage port>``; the
guard refuses a SOCKET to a non-loopback address, so urllib connected to
127.0.0.1 -- allowed -- and the request left the machine from the proxy,
where the in-process tripwire could not see it. That blind spot is closed
by stripping the proxy variables with the provider keys
(``external_agent_guard.STRIPPED_PREFIXES``), and the red-to-green
evidence for it is in the guard's own suite:

    tests/test_the_suite_reaches_no_external_agent.py::test_a_remote_endpoint_is_refused
    tests/test_the_suite_reaches_no_external_agent.py::test_a_model_client_meets_the_tripwire_not_the_endpoint

Both fail on the commit before the fix (the urllib leg records no trip) and
pass after it.

Runs the guard in a CHILD pytest process: ``install`` is idempotent per
process, so only a fresh process actually installs it. The destination is a
literal address from TEST-NET-1 (RFC 5737) -- no DNS, no route, the same
answer on every machine.
"""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap

TESTS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(TESTS)

#: 192.0.2.1 is TEST-NET-1: no route, no DNS, no host behind it anywhere.
#: The path is deliberately NOT ``/models``: conftest answers a model-list
#: URL offline before it is opened (_model_lists_are_offline), which would
#: hide the attempt and leave this test measuring that fixture.
REMOTE_IP = "192.0.2.1"

#: The inner run loads the SUITE's conftest, which is what installs the
#: guard (pytest_configure). Without it the inner process is unguarded and
#: every case here would pass for the wrong reason.
INNER_CONFTEST = """
import importlib.util, sys
sys.path.insert(0, %(tests)r)
_spec = importlib.util.spec_from_file_location(
    "delfin_suite_conftest", %(conftest)r)
_c = importlib.util.module_from_spec(_spec)
sys.modules["delfin_suite_conftest"] = _c
_spec.loader.exec_module(_c)
globals().update({k: getattr(_c, k) for k in dir(_c) if not k.startswith("__")})
"""

_CASE = """\
    import urllib.request

    def test_swallowed_open():
        try:
            urllib.request.urlopen("http://%(ip)s/v1/ping", timeout=2).close()
        except Exception:
            pass

    def test_clean():
        assert True
    """ % {"ip": REMOTE_IP}


def _run_child(tmp_path, env):
    inner = tmp_path / "inner"
    inner.mkdir()
    (inner / "conftest.py").write_text(INNER_CONFTEST % {
        "tests": TESTS, "conftest": os.path.join(TESTS, "conftest.py")})
    (inner / "test_inner.py").write_text(textwrap.dedent(_CASE))
    return subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
         "-o", "addopts=", "--rootdir", str(inner), str(inner / "test_inner.py")],
        cwd=str(inner), env=env, capture_output=True, text=True, timeout=300)


def test_urllib_to_a_remote_address_is_refused_and_named(tmp_path):
    """A product path that opens a URL on a remote host is stopped with the
    address in the refusal, even though the code under test swallowed the
    exception -- that is what makes an egress visible instead of silent."""
    (tmp_path / "home").mkdir()
    env = {
        "PATH": os.environ.get("PATH", "/usr/bin:/bin"),
        "HOME": str(tmp_path / "home"),
        "PYTHONPATH": os.pathsep.join([REPO, TESTS]),
    }
    r = _run_child(tmp_path, env)
    out = r.stdout + r.stderr
    assert r.returncode != 0, f"the guard refused nothing:\n{out}"
    assert "external agent tripwire" in out, out
    assert f"network connection to {REMOTE_IP}:80" in out, (
        f"the refusal did not name the destination:\n{out}")
