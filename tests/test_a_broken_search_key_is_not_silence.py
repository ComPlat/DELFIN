"""A configured search engine that failed must say so.

Input: ``web_tools.web_search`` with one keyed engine configured whose
function raises, and the unkeyed fallback failing too.

Output: the error the caller receives names the keyed failure.

Measured 2026-11-22: ``_keyed_search`` built an error string
(``"<engine>: <Error>: <message>"``) and the only mention of
``keyed_error`` in the module was its own assignment. A wrong, expired
or rate-limited key therefore looked exactly like a search that found
nothing, and the person was never told. The string is now appended to
whatever the fallback reports, and only when the fallback also fails: a
search that succeeds without the keyed engine needs no apology.
"""
from __future__ import annotations

import urllib.error

import pytest

from delfin.agent import web_tools


@pytest.fixture
def one_broken_engine(monkeypatch):
    """One configured engine, and calling it raises."""
    def _boom(query, max_results, timeout_s):
        raise RuntimeError("401 Unauthorized")

    monkeypatch.setattr(web_tools, "_KEYED_ENGINES", [("brave", _boom)])
    monkeypatch.setattr(web_tools, "keyed_engines_configured",
                        lambda: ["brave"])
    monkeypatch.setattr(web_tools, "_check_url", lambda url: "")


def _fallback_raises(exc):
    """Make the unkeyed fallback fail at the seam the product uses.

    ``_fetch_bytes`` is the one place this module opens a URL (the same
    seam tests/test_web_search_asks_which_engines_are_configured.py
    patches). Patching urllib directly misses it and reaches the
    suite's network tripwire instead, which is how this test first
    failed.
    """
    def _fetch(url, timeout_s, want_status=False, headers=None):
        raise exc
    return _fetch


def test_a_failing_key_is_named_when_the_fallback_also_fails(
        one_broken_engine, monkeypatch):
    monkeypatch.setattr(web_tools, "_fetch_bytes",
                        _fallback_raises(urllib.error.HTTPError(
                            "u", 503, "Service Unavailable", {}, None)))
    out = web_tools.web_search("beta HRS crystal", max_results=3)
    err = str(out.get("error") or "")
    assert "503" in err, err
    assert "brave" in err and "401 Unauthorized" in err, (
        "a configured engine that failed must be named, or a broken key "
        f"is indistinguishable from an empty result: {err!r}")


def test_a_network_error_also_carries_the_keyed_failure(
        one_broken_engine, monkeypatch):
    monkeypatch.setattr(web_tools, "_fetch_bytes",
                        _fallback_raises(
                            urllib.error.URLError("name resolution")))
    err = str(web_tools.web_search("x", max_results=1).get("error") or "")
    assert "network error" in err and "brave" in err, err


def test_without_a_configured_engine_no_apology_is_added(monkeypatch):
    """Control: nothing was tried, so there is nothing to report. The
    fallback's own error stands alone."""
    monkeypatch.setattr(web_tools, "_KEYED_ENGINES", [])
    monkeypatch.setattr(web_tools, "keyed_engines_configured", lambda: [])
    monkeypatch.setattr(web_tools, "_check_url", lambda url: "")
    monkeypatch.setattr(web_tools, "_fetch_bytes",
                        _fallback_raises(urllib.error.HTTPError(
                            "u", 503, "Service Unavailable", {}, None)))
    err = str(web_tools.web_search("x", max_results=1).get("error") or "")
    assert "503" in err
    assert "tried first and failed" not in err, err


def test_a_succeeding_fallback_says_nothing_about_the_key(
        one_broken_engine, monkeypatch):
    """Control: the search worked. A failed keyed attempt is then an
    implementation detail, not an error the caller must read."""
    page = b"<html><a class='result__a' href='http://x'>title</a></html>"

    def _ok(url, timeout_s, want_status=False, headers=None):
        # (body, content type, status): the shape the caller unpacks.
        return page, "text/html", 200

    monkeypatch.setattr(web_tools, "_fetch_bytes", _ok)
    out = web_tools.web_search("x", max_results=1)
    assert "error" not in out or "tried first and failed" not in str(
        out.get("error") or "")
