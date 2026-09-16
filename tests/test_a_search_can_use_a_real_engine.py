"""With a key, web_search asks a real search engine first; without one it
searches as before.

Asked for on 2026-09-16 ("googeln wie du"): the scraped DuckDuckGo page is
challenged on this network, and OpenAlex/Wikipedia are indexes, not the
web. GOOGLE_CSE_API_KEY + GOOGLE_CSE_CX or BRAVE_SEARCH_API_KEY in the
credential store switch a real engine on; nothing changes without them.
"""
import json

from delfin.agent import web_tools as W


def _creds(monkeypatch, **values):
    monkeypatch.setattr(W, "_credential", lambda name: values.get(name, ""))


def test_without_a_key_no_engine_is_configured_and_the_chain_is_untouched(monkeypatch):
    _creds(monkeypatch)
    assert W.keyed_engines_configured() == []
    assert W._keyed_search("x", 5, 5) == ([], "", "")


def test_google_answers_first_when_its_key_is_there(monkeypatch):
    _creds(monkeypatch, GOOGLE_CSE_API_KEY="k", GOOGLE_CSE_CX="cx")
    seen = {}

    def fake_fetch(url, timeout_s, want_status=False, headers=None):
        seen["url"] = url
        body = json.dumps({"items": [{"title": "ORCA 6", "link": "https://orcaforum.kofo.mpg.de/x",
                                      "snippet": "released July 2024"}]}).encode()
        return body, "application/json"
    monkeypatch.setattr(W, "_fetch_bytes", fake_fetch)
    out = W.web_search("ORCA 6 release")
    assert out["source"] == "google"
    assert out["results"][0]["url"].startswith("https://orcaforum")
    assert "customsearch/v1" in seen["url"] and "cx=cx" in seen["url"]


def test_brave_sends_its_token_in_the_header(monkeypatch):
    _creds(monkeypatch, BRAVE_SEARCH_API_KEY="tok")
    seen = {}

    def fake_fetch(url, timeout_s, want_status=False, headers=None):
        seen["headers"] = headers or {}
        body = json.dumps({"web": {"results": [{"title": "t", "url": "https://a.b/c", "description": "d"}]}}).encode()
        return body, "application/json"
    monkeypatch.setattr(W, "_fetch_bytes", fake_fetch)
    out = W.web_search("anything")
    assert out["source"] == "brave" and out["results"][0]["url"] == "https://a.b/c"
    assert seen["headers"].get("X-Subscription-Token") == "tok"


def test_a_failing_keyed_engine_falls_through_to_the_keyless_chain(monkeypatch):
    _creds(monkeypatch, BRAVE_SEARCH_API_KEY="tok")

    def fake_fetch(url, timeout_s, want_status=False, headers=None):
        if "brave.com" in url:
            raise OSError("quota")
        return b"<html>no results found</html>", "text/html", 200
    monkeypatch.setattr(W, "_fetch_bytes", fake_fetch)
    monkeypatch.setattr(W, "_ddg_instant_answer", lambda q, t: [])
    monkeypatch.setattr(W, "_openalex_search", lambda q, m, t: [])
    monkeypatch.setattr(W, "_wikipedia_search", lambda q, m, t: [])
    out = W.web_search("anything")
    assert "results" in out or "error" in out
    assert out.get("source") != "brave"


def test_the_tool_description_names_the_engines():
    from delfin.agent import api_client as A
    import inspect
    src = inspect.getsource(A)
    assert "Uses Google or " in src and "Brave when a key is configured" in src
