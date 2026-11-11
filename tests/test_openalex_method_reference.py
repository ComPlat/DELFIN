"""find_method_reference: original-paper candidates for an ORCA keyword.

Recorded OpenAlex responses only; no network.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse

import pytest

from delfin.agent import scholarly_sources as ss
from delfin.agent import source_apis, web_tools

BECKE = {"id": "https://openalex.org/W1", "doi": "https://doi.org/10.1063/1.464913",
         "title": "Density-functional thermochemistry. III. The role of exact exchange",
         "publication_year": 1993,
         "authorships": [{"author": {"display_name": "Axel D. Becke"}}],
         "primary_location": {"source": {"display_name": "The Journal of Chemical Physics"}},
         "cited_by_count": 100000, "type": "article"}
LYP = {"id": "https://openalex.org/W2", "doi": "https://doi.org/10.1103/physrevb.37.785",
       "title": "Development of the Colle-Salvetti correlation-energy formula into a "
                "functional of the electron density",
       "publication_year": 1988,
       "authorships": [{"author": {"display_name": "Chengteh Lee"}}],
       "primary_location": {"source": {"display_name": "Physical Review B"}},
       "cited_by_count": 90000, "type": "article"}


class Net:
    def __init__(self, answer):
        self.answer = answer
        self.urls = []

    def __call__(self, url, *, timeout_s, headers=None):
        self.urls.append(url)
        return json.dumps({"results": self.answer(url)}).encode()

    def params(self, i):
        return {k: v[0] for k, v in
                urllib.parse.parse_qs(urllib.parse.urlsplit(self.urls[i]).query).items()}


@pytest.fixture(autouse=True)
def _isolated(monkeypatch):
    monkeypatch.setattr(ss, "contact_email", lambda: "")
    monkeypatch.setattr(ss, "_credential", lambda name: "")
    monkeypatch.setattr(web_tools, "outbound_secret_reason", lambda text: None)
    source_apis.clear_cache()
    yield
    source_apis.clear_cache()


def test_known_method_searches_the_original_titles(monkeypatch):
    def answer(url):
        return [BECKE] if "exact exchange" in urllib.parse.unquote_plus(url) else [LYP]
    net = Net(answer)
    monkeypatch.setattr(ss, "_get", net)
    out = ss.find_method_reference("B3LYP")
    assert out["status"] == "ok"
    dois = [c["doi"] for c in out["candidates"]]
    assert dois == ["10.1063/1.464913", "10.1103/physrevb.37.785"]
    assert all(c["match"] == "known-title" for c in out["candidates"])
    p = net.params(0)
    assert p["filter"].startswith("title.search:Density-functional thermochemistry")
    assert "primary_topic.field.id:fields/16" in p["filter"]
    assert p["sort"] == "cited_by_count:desc"
    assert "verify" in out["note"].lower()


def test_aliases_share_the_original_paper(monkeypatch):
    net = Net(lambda url: [BECKE])
    monkeypatch.setattr(ss, "_get", net)
    ss.find_method_reference("def2-TZVP")
    assert "Balanced basis sets" in net.params(0)["filter"]


def test_unknown_name_is_a_title_search_and_cached(monkeypatch):
    net = Net(lambda url: [BECKE])
    monkeypatch.setattr(ss, "_get", net)
    first = ss.find_method_reference("XYZ-99")
    assert first["candidates"][0]["match"] == "name-in-title"
    assert net.params(0)["filter"].startswith("title.search:XYZ-99,")
    again = ss.find_method_reference("XYZ-99")
    assert again == first and len(net.urls) == 1


def test_too_narrow_topic_filter_falls_back_to_title_only(monkeypatch):
    net = Net(lambda url: [] if "primary_topic" in urllib.parse.unquote(url) else [BECKE])
    monkeypatch.setattr(ss, "_get", net)
    out = ss.find_method_reference("PBE")
    assert out["status"] == "ok" and len(net.urls) == 2
    assert "primary_topic" not in net.params(1)["filter"]


def test_network_error_is_reported_not_raised(monkeypatch):
    def broken(url, *, timeout_s, headers=None):
        raise urllib.error.URLError("unreachable")
    monkeypatch.setattr(ss, "_get", broken)
    out = ss.find_method_reference("SMD")
    assert out["status"] == "error" and "network error" in out["note"]


def test_nothing_found_and_bad_input(monkeypatch):
    monkeypatch.setattr(ss, "_get", Net(lambda url: []))
    assert ss.find_method_reference("NoSuchMethod")["status"] == "not_found"
    assert ss.find_method_reference("")["status"] == "invalid_input"
    assert ss.find_method_reference("x" * 200)["status"] == "invalid_input"


def test_source_apis_re_exports_the_agent_functions(monkeypatch):
    monkeypatch.setattr(ss, "_get", Net(lambda url: [BECKE]))
    out = source_apis.find_method_reference("B3LYP", max_results=1)
    assert out["candidates"][0]["doi"] == "10.1063/1.464913"
    fed = source_apis.federated_search("exact exchange", sources=["openalex"],
                                       unpaywall=False)
    assert fed["status"] == {"openalex": "ok"}
