"""web_fetch: the documented refusal, not the network.

Rule 6: a test binds to DELFIN, never to this machine. A DNS failure is
what THIS cage shows; a machine with network answers differently. The
machine-independent documented behaviour is the URL deny-list: localhost,
RFC1918, cloud metadata and *.internal hosts are refused before any
network call, and the refusal names the host and the list. Live reachability
is a real-run question (scratch/realrun.json), not a unit test's.
"""

from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions

import json

import pytest


@pytest.fixture
def perms(tmp_path):
    perms = KitToolPermissions(workspace=str(tmp_path / "ws"))
    perms.mode = "acceptEdits"
    perms.task_session_id = "f3-web"
    return perms


@pytest.fixture
def ex():
    return _DocToolExecutor()


@pytest.mark.parametrize("url,host", [
    ("http://localhost/index.html", "localhost"),
    ("http://127.0.0.1/index.html", "127.0.0.1"),
    ("http://169.254.169.254/latest/meta-data", "169.254.169.254"),
    ("http://10.0.0.1/x", "10.0.0.1"),
    ("http://192.168.1.1/x", "192.168.1.1"),
    ("http://db.internal/x", "db.internal"),
])
def test_web_fetch_refuses_internal_hosts_and_names_the_cause(ex, perms, url, host):
    out = ex.execute("web_fetch", {"url": url}, perms)
    # Two layers refuse an internal host, and which answers first depends on
    # what ran earlier in the process (CI's shard 3 met the executor's
    # deny-list first, 2026-11-25). Both must refuse and name the host;
    # neither may hand back fetched content.
    assert host in out, out
    assert ("web_fetch refused" in out and "Nothing was sent" in out) or "deny-list" in out, out


@pytest.mark.parametrize("arguments", [
    {"query": ""},
    {"query": "   "},
    {},
])
def test_web_search_refuses_an_empty_query_before_any_lookup(ex, perms, arguments):
    """The half of ``web_search`` that is answerable without a network.

    ``_execute_web_search`` (api_client.py) refuses before it reaches
    ``web_tools``, so this outcome is the same on a machine with network and
    one without -- unlike the live search, which is a real-run question.
    Covered here because the package's own harness exercises ``web_fetch``
    and left ``web_search`` with no test at all.
    """
    out = json.loads(ex.execute("web_search", arguments, perms))
    assert "error" in out, out
    assert "query is required" in out["error"], out["error"]


def test_web_fetch_refusal_is_structured(ex, perms):
    out = ex.execute("web_fetch", {"url": "http://10.0.0.1/x"}, perms)
    body = out.split("\n", 1)[1] if out.startswith("[UNTRUSTED") else out
    parsed = json.loads(body)
    assert "web_fetch refused" in parsed["error"]
    assert "10.0.0.1" in parsed["error"]
    assert "Nothing was sent" in parsed["error"]
