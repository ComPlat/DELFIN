"""A search that could not run must not look like a search that found nothing.

Measured on this host: DuckDuckGo answers every real query with HTTP 202
— an anti-bot challenge page with no result markup. The tool reported
``result_count: 0`` and no error, which is the shape of "nothing matched".

That is the worst possible failure mode for this model in this mode. The
office prompt forbids answering from memory; a tool that says "no results"
without saying it never ran is precisely what makes answering from memory
look like the only option left. The same field session shows the model
then reaching a WORKING search under another server, which is what the
error message now points it at.

Three outcomes have to stay distinguishable, and calling a genuinely
empty search "blocked" would be a different lie:

* 202                      -> refused by the backend
* 200 with a "no results"  -> genuinely empty, and that is a real answer
* 200 with neither         -> unparseable; whether anything matched is
                              unknown, and unknown is not zero
"""

from __future__ import annotations

import pytest

import delfin.agent.web_tools as W


@pytest.fixture
def offline(monkeypatch):
    """Serve a canned page. No network, no dependence on DuckDuckGo's mood."""
    def _serve(status: int, html: str):
        def _fetch(url, timeout_s, want_status=False):
            if want_status:
                return html.encode(), "text/html", status
            return html.encode(), "text/html"
        monkeypatch.setattr(W, "_fetch_bytes", _fetch)
        monkeypatch.setattr(W, "_ddg_instant_answer", lambda q, t: [])
        monkeypatch.setattr(W, "_check_url", lambda u: None)
    return _serve


_RESULT_MARKUP = ('<a class="result__a" href="https://example.com">Title</a>'
                  '<a class="result__snippet">Snippet</a>')


def test_a_challenge_is_reported_as_a_failure(offline):
    offline(202, "<html>please verify you are human</html>")
    out = W.web_search("Kostenstelle Dienstreise Formular")
    assert "error" in out
    assert out["result_count"] == 0
    assert "202" in out["error"]


def test_the_failure_says_none_should_be_inferred(offline):
    """The message has to close the door the model would otherwise walk
    through: no results retrieved is not licence to recall some."""
    offline(202, "<html>challenge</html>")
    assert "none should be inferred" in W.web_search("x")["error"]


def test_the_failure_names_the_alternative(offline):
    """A working search exists on this cluster under another server. A
    refusal that does not say so leaves the model with nothing to try."""
    out = W.web_search("x") if False else None
    offline(202, "<html>challenge</html>")
    err = W.web_search("x")["error"]
    assert "another server" in err
    assert "ask them for a URL" in err


def test_a_genuinely_empty_search_is_not_called_a_failure(offline):
    """The engine says so in words when it means it. Reporting that as
    'blocked' would be a different lie."""
    offline(200, "<html><p>No results found for asdkjhasd</p></html>")
    out = W.web_search("asdkjhasd")
    assert "error" not in out
    assert out["result_count"] == 0


@pytest.mark.parametrize("phrase", [
    "No results found", "keine Ergebnisse", "no results for",
])
def test_the_empty_notice_is_recognised_in_either_language(offline, phrase):
    offline(200, f"<html><p>{phrase}</p></html>")
    assert "error" not in W.web_search("q")


def test_an_unparseable_page_is_reported_as_unknown_not_zero(offline):
    offline(200, "<html>markup we cannot read</html>")
    out = W.web_search("q")
    assert "error" in out
    assert "unknown" in out["error"]


def test_real_results_still_come_back(offline):
    offline(200, _RESULT_MARKUP)
    out = W.web_search("q")
    assert "error" not in out
    assert out["result_count"] == 1
    assert out["results"][0]["url"] == "https://example.com"


def test_the_instant_answer_fallback_still_rescues_a_blocked_page(monkeypatch):
    """It was there for exactly this case and must keep priority over the
    new error: a real answer beats an honest failure."""
    def _fetch(url, timeout_s, want_status=False):
        html = b"<html>challenge</html>"
        return (html, "text/html", 202) if want_status else (html, "text/html")
    monkeypatch.setattr(W, "_fetch_bytes", _fetch)
    monkeypatch.setattr(W, "_check_url", lambda u: None)
    monkeypatch.setattr(
        W, "_ddg_instant_answer",
        lambda q, t: [{"title": "Caffeine", "url": "https://en.wikipedia.org/x",
                       "snippet": "..."}])
    out = W.web_search("caffeine")
    assert "error" not in out
    assert out["result_count"] == 1


def test_the_status_is_available_to_callers_that_ask(monkeypatch):
    """_fetch_bytes stays backward compatible: two values unless asked."""
    def _fetch_raw(url, timeout_s, want_status=False):
        return (b"x", "text/html", 200) if want_status else (b"x", "text/html")

    monkeypatch.setattr(W, "_fetch_bytes", _fetch_raw)
    assert len(W._fetch_bytes("u", 1)) == 2
    assert len(W._fetch_bytes("u", 1, want_status=True)) == 3


def test_an_unknown_status_is_not_treated_as_a_failure(monkeypatch):
    """A replacement predating the want_status parameter yields no status.
    Unknown is not evidence of a challenge, and claiming one would be the
    same overreach the fix exists to remove -- pointed the other way."""
    def _old_signature(url, timeout_s):
        return b"<html>nothing recognisable</html>", "text/html"

    monkeypatch.setattr(W, "_fetch_bytes", _old_signature)
    monkeypatch.setattr(W, "_ddg_instant_answer", lambda q, t: [])
    monkeypatch.setattr(W, "_check_url", lambda u: None)
    out = W.web_search("q")
    assert "error" not in out
    assert out["result_count"] == 0


def test_a_two_value_fetch_does_not_break_the_search(monkeypatch):
    """The parameter is new; adding it must not break existing callers."""
    def _old_signature(url, timeout_s):
        return _RESULT_MARKUP.encode(), "text/html"

    monkeypatch.setattr(W, "_fetch_bytes", _old_signature)
    monkeypatch.setattr(W, "_check_url", lambda u: None)
    out = W.web_search("q")
    assert "error" not in out
    assert out["result_count"] == 1
