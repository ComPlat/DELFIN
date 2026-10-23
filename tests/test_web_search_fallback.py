"""web_search must degrade gracefully when the DuckDuckGo HTML endpoint is
blocked (a 202 anomaly page with no results, as happens from datacenter IPs):
it falls back to the Instant Answer JSON API. Offline — network is mocked.
"""

from __future__ import annotations

import json

from delfin.agent import web_tools


def _mock_fetch(html_body: bytes, ia_payload: dict):
    def fake(url, timeout_s):
        if "api.duckduckgo.com" in url:
            return (json.dumps(ia_payload).encode("utf-8"), "application/json")
        return (html_body, "text/html")
    return fake


def test_falls_back_to_instant_answer_when_html_empty(monkeypatch):
    monkeypatch.setattr(web_tools, "_fetch_bytes", _mock_fetch(
        b"<html>202 anomaly, no result__a here</html>",
        {"Heading": "Toluene",
         "AbstractText": "An aromatic hydrocarbon, C6H5CH3.",
         "AbstractURL": "https://en.wikipedia.org/wiki/Toluene",
         "RelatedTopics": [
             {"Text": "Benzene ring", "FirstURL": "https://example.com/benzene"},
         ]},
    ))
    r = web_tools.web_search("toluene")
    assert r["result_count"] >= 1
    assert r["source"] == "duckduckgo-instant-answer"
    assert "aromatic" in r["results"][0]["snippet"].lower()
    assert r["results"][0]["url"] == "https://en.wikipedia.org/wiki/Toluene"


def test_html_results_win_when_present(monkeypatch):
    html = (
        b'<a class="result__a" href="https://ex.com/x">Real Title</a>'
        b'<a class="result__snippet">real snippet</a>'
    )
    # IA payload present but must NOT be used when HTML has hits.
    monkeypatch.setattr(web_tools, "_fetch_bytes", _mock_fetch(
        html, {"AbstractText": "should not appear"}))
    r = web_tools.web_search("anything")
    assert r["source"] == "duckduckgo-html"
    assert r["results"][0]["title"] == "Real Title"


def test_no_results_anywhere_is_clean_empty(monkeypatch):
    monkeypatch.setattr(web_tools, "_fetch_bytes", _mock_fetch(
        b"<html>nothing</html>", {}))
    r = web_tools.web_search("zzznotathing")
    assert r["result_count"] == 0
    assert r["results"] == []
