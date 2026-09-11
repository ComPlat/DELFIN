"""web_search answered "network error: [Errno 104] Connection reset by
peer" twice in a field session and once from the office network, in the
same minute, and the agent concluded the thing did not exist "in the
web" (2026-09-11). A reset is not an answer: the search is tried again,
twice, before it is reported as a network error that names the count.
"""

import urllib.error

from delfin.agent import web_tools as w


def test_a_reset_is_retried_and_the_third_answer_is_used(monkeypatch):
    calls = {"n": 0}
    body = b'<div class="result__body"><a class="result__a" href="https://example.org/x">Example</a><a class="result__snippet">snip</a></div>'

    def fake_fetch(url, timeout_s, want_status=False):
        calls["n"] += 1
        if calls["n"] < 3:
            raise urllib.error.URLError(ConnectionResetError(104, "Connection reset by peer"))
        return (body, "text/html", 200)

    monkeypatch.setattr(w, "_fetch_bytes", fake_fetch)
    monkeypatch.setattr(w.time, "sleep", lambda s: None)
    out = w.web_search("tada molecule", max_results=3)
    assert not out.get("error"), out
    assert calls["n"] == 3


def test_a_reset_every_time_is_reported_with_the_count(monkeypatch):
    def always_reset(url, timeout_s, want_status=False):
        raise urllib.error.URLError(ConnectionResetError(104, "Connection reset by peer"))

    monkeypatch.setattr(w, "_fetch_bytes", always_reset)
    monkeypatch.setattr(w.time, "sleep", lambda s: None)
    out = w.web_search("tada molecule")
    assert "Connection reset" in out["error"] and "after 3 attempts" in out["error"]


def test_another_network_error_is_not_retried(monkeypatch):
    calls = {"n": 0}

    def refused(url, timeout_s, want_status=False):
        calls["n"] += 1
        raise urllib.error.URLError(OSError("Name or service not known"))

    monkeypatch.setattr(w, "_fetch_bytes", refused)
    out = w.web_search("tada molecule")
    assert "network error" in out["error"] and calls["n"] == 1


def test_the_reset_is_recognised_by_its_wording_too():
    assert w._is_connection_reset(urllib.error.URLError("[Errno 104] Connection reset by peer"))
    assert not w._is_connection_reset(urllib.error.URLError("timed out"))
