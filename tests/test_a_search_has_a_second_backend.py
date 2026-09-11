"""DuckDuckGo challenged every query from here and from the cluster on
2026-09-11 (HTTP 202 on all three of its endpoints), and the only
answer web_search had was "the backend refused". OpenAlex and Wikipedia
answer keyless from both networks; a scientific question gets the
scholarly index first, and every result names the index it came from.
"""

import json
import urllib.error

from delfin.agent import web_tools as w

_CHALLENGE = (b"<html><body>anomaly detected, challenge</body></html>", "text/html", 202)


def _openalex_body():
    return json.dumps({"results": [
        {"title": "Triazatruxene: a rigid donor for TADF", "doi": "https://doi.org/10.1002/advs.201700989",
         "publication_year": 2018,
         "primary_location": {"landing_page_url": "https://x", "source": {"display_name": "Advanced Science"}},
         "authorships": [{"author": {"display_name": "A. Author"}}, {"author": {"display_name": "B. Author"}}]},
        {"title": "no link", "doi": "", "primary_location": {}},
    ]}).encode()


def test_a_challenged_query_falls_back_to_openalex(monkeypatch):
    calls = []

    def fetch(url, timeout_s, want_status=False):
        calls.append(url)
        if "duckduckgo.com/html" in url:
            return _CHALLENGE
        if "api.duckduckgo.com" in url:
            return (b"{}", "application/json", 200)
        if "api.openalex.org" in url:
            return (_openalex_body(), "application/json", 200)
        raise AssertionError(url)

    monkeypatch.setattr(w, "_fetch_bytes", fetch)
    out = w.web_search("triazatruxene TADF", max_results=5)
    assert not out.get("error"), out
    assert out["source"] == "openalex" and out["result_count"] == 1
    hit = out["results"][0]
    assert hit["url"] == "https://doi.org/10.1002/advs.201700989"
    assert "Advanced Science" in hit["snippet"] and "2018" in hit["snippet"]
    assert hit["source"] == "openalex"


def test_wikipedia_answers_when_openalex_has_nothing(monkeypatch):
    def fetch(url, timeout_s, want_status=False):
        if "duckduckgo.com/html" in url:
            return _CHALLENGE
        if "api.duckduckgo.com" in url:
            return (b"{}", "application/json", 200)
        if "api.openalex.org" in url:
            return (b'{"results": []}', "application/json", 200)
        if "wikipedia.org" in url:
            return (json.dumps(["tadf", ["Thermally activated delayed fluorescence"], ["A mechanism"],
                                ["https://en.wikipedia.org/wiki/TADF"]]).encode(), "application/json", 200)
        raise AssertionError(url)

    monkeypatch.setattr(w, "_fetch_bytes", fetch)
    out = w.web_search("thermally activated delayed fluorescence")
    assert out["source"] == "wikipedia" and out["results"][0]["url"].endswith("/TADF")


def test_every_backend_failing_is_still_named_as_a_refusal(monkeypatch):
    def fetch(url, timeout_s, want_status=False):
        if "duckduckgo.com/html" in url:
            return _CHALLENGE
        raise urllib.error.URLError("down")

    monkeypatch.setattr(w, "_fetch_bytes", fetch)
    out = w.web_search("anything")
    assert out.get("error") and "anti-bot challenge" in out["error"]
    assert out["result_count"] == 0


def test_a_working_duckduckgo_is_not_second_guessed(monkeypatch):
    body = b'<div class="result__body"><a class="result__a" href="https://example.org/x">Example</a><a class="result__snippet">snip</a></div>'

    def fetch(url, timeout_s, want_status=False):
        if "duckduckgo.com/html" in url:
            return (body, "text/html", 200)
        raise AssertionError("no fallback should be called: " + url)

    monkeypatch.setattr(w, "_fetch_bytes", fetch)
    out = w.web_search("example")
    assert out["source"] == "duckduckgo-html" and out["result_count"] == 1
