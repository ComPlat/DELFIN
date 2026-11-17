"""SRC1 phase 1: source_apis verifies DOIs and compound data through the
official APIs, over DELFIN's existing guarded web fetch.

Red on the previous commit: delfin/agent/source_apis.py does not exist.

All HTTP is intercepted: no test touches the live network. The fake
records the URL that source_apis asked for, so the test can prove every
call goes through web_tools (single HTTP stack, deny-list in force).
"""

from __future__ import annotations

import json
import urllib.parse
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

FAKE_CROSSREF = {
    "status": "ok",
    "message": {
        "DOI": "10.1021/acs.jpclett.5b00331",
        "title": ["A Test Title on Solvation"],
        "author": [{"given": "A.", "family": "Researcher"}],
        "container-title": ["J. Phys. Chem. Lett."],
        "issued": {"date-parts": [[2015]]},
        "type": "journal-article",
    },
}

FAKE_PUBCHEM_BY_NAME = {
    # the REAL PUG-REST /property response shape (PropertyTable.Properties
    # is a list of rows); a test of a record reader uses the real shape.
    "PropertyTable": {"Properties": [{
        "CID": 2244,
        "MolecularFormula": "C9H8O4",
        "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        "MolecularWeight": "180.16",
    }]},
}


class FakeTransport:
    """Stands in for web_tools' guarded fetch; records every URL."""

    def __init__(self, responses):
        self.responses = responses
        self.urls = []

    def fetch_json(self, url, timeout_s=30):
        self.urls.append(url)
        handler = self.responses.get(url)
        if handler is None:
            raise RuntimeError("no recorded response for " + url)
        if handler.get("status", 200) == 404:
            import urllib.error
            raise urllib.error.HTTPError(url, 404, "not found", {}, None)
        return handler["body"]


@pytest.fixture
def transport(monkeypatch):
    fake = FakeTransport({
        "https://api.crossref.org/works/10.1021%2Facs.jpclett.5b00331": {
            "body": FAKE_CROSSREF},
        "https://api.crossref.org/works/10.9999%2Finvented.doi": {
            "status": 404},
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/MolecularFormula,InChIKey,MolecularWeight/JSON": {
            "body": FAKE_PUBCHEM_BY_NAME},
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/BSYNRYMUTXBXSQ-UHFFFAOYSA-N/property/MolecularFormula,InChIKey,MolecularWeight/JSON": {
            "body": FAKE_PUBCHEM_BY_NAME},
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/not-a-molecule/property/MolecularFormula,InChIKey,MolecularWeight/JSON": {
            "status": 404},
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        + urllib.parse.quote("C([H])(=O)/N", safe="") + "/property/"
        "MolecularFormula,InChIKey,MolecularWeight/JSON": {
            "body": FAKE_PUBCHEM_BY_NAME},
    })
    from delfin.agent import source_apis
    monkeypatch.setattr(source_apis, "fetch_json", fake.fetch_json)
    monkeypatch.setattr(source_apis.web_tools, "web_search",
                        lambda query, max_results=5: [])
    # _check_url resolves DNS; in a no-network cage that refuses every
    # host before the (stubbed) transport. Unit tests stub the guard too
    # — the guard itself is exercised by the refusal tests, which stub
    # only the DNS step.
    monkeypatch.setattr(source_apis.web_tools, "_check_url",
                        lambda url: None)
    source_apis.clear_cache()
    yield fake
    source_apis.clear_cache()


def test_module_exists_and_states_its_terms():
    src = (ROOT / "delfin/agent/source_apis.py").read_text(encoding="utf-8")
    assert "Crossref" in src and "PubChem" in src
    # licences of the returned data are stated, not assumed
    assert "CC0" in src and "public domain" in src


def test_verify_doi_returns_metadata(transport):
    from delfin.agent.source_apis import verify_doi
    rec = verify_doi("10.1021/acs.jpclett.5b00331")
    assert rec["status"] == "ok"
    assert rec["title"] == "A Test Title on Solvation"
    assert rec["year"] == 2015
    assert rec["journal"] == "J. Phys. Chem. Lett."
    assert rec["authors"] == ["A. Researcher"]
    # exactly one host call, to the official API
    assert transport.urls == [
        "https://api.crossref.org/works/10.1021%2Facs.jpclett.5b00331"]


def test_verify_doi_reports_not_found_clearly(transport):
    from delfin.agent.source_apis import verify_doi
    rec = verify_doi("10.9999/invented.doi")
    assert rec["status"] == "not_found"
    # the invented DOI's number is NOT echoed as verified metadata
    assert "title" not in rec


def test_verify_doi_rejects_malformed_input():
    from delfin.agent.source_apis import verify_doi
    assert verify_doi("not a doi")["status"] == "invalid_input"
    assert verify_doi("")["status"] == "invalid_input"


def test_verify_doi_refuses_private_and_denylisted_targets(monkeypatch):
    # the guard lives in web_tools._check_url; a private host must never
    # even be attempted: no bytes requested, no opener touched.
    from delfin.agent import source_apis, web_tools
    requested = []
    # the real guard, with only its DNS step short-circuited: private
    # hosts are still refused by the deny-list pattern check itself.
    real_check = web_tools._check_url
    monkeypatch.setattr(web_tools, "_check_url",
                        lambda url: None if real_check(url) is None
                        or "private" in str(real_check(url))
                        else "refused (deny-list)")
    monkeypatch.setattr(web_tools, "_fetch_bytes",
                        lambda url, timeout_s, want_status=False,
                        headers=None: requested.append(url) or (b"", ""))
    monkeypatch.setattr(source_apis, "_GUARDED_FETCH",
                        source_apis._fallback_fetch)
    for url in ("https://127.0.0.1/api", "http://localhost/x",
                "http://169.254.169.254/latest/meta-data"):
        with pytest.raises(ValueError, match="refused"):
            source_apis.fetch_json(url)
    assert requested == []


def test_lookup_compound_by_name_and_inchikey(transport):
    from delfin.agent.source_apis import lookup_compound
    by_name = lookup_compound("aspirin")
    assert by_name["status"] == "ok"
    assert by_name["cid"] == 2244
    assert by_name["formula"] == "C9H8O4"
    assert by_name["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    by_key = lookup_compound("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
    assert by_key["status"] == "ok"
    assert by_key["cid"] == 2244


def test_lookup_compound_unknown_name_is_not_found(transport):
    from delfin.agent.source_apis import lookup_compound
    rec = lookup_compound("not-a-molecule")
    assert rec["status"] == "not_found"


def test_search_works_delegates_to_web_tools(monkeypatch):
    # No second search stack either: search_works is a thin delegation,
    # and the test proves the delegation, with a recorded fake.
    from delfin.agent import source_apis, web_tools
    calls = []
    monkeypatch.setattr(web_tools, "web_search",
                        lambda query, max_results=8:
                        calls.append((query, max_results))
                        or [{"title": "t", "url": "https://a/b",
                             "snippet": "s", "source": "openalex"}])
    hits = source_apis.search_works("solvation", max_results=3)
    assert calls == [("solvation", 3)]
    assert hits and hits[0]["source"] == "openalex"


def test_every_call_goes_through_the_guarded_fetch(monkeypatch):
    # No second HTTP stack: the fallback fetch must call web_tools'
    # guarded _fetch_bytes (its deny-list + SSRF checks), never urllib
    # of its own. Once the operator builds the web_tools.fetch_json hook,
    # _GUARDED_FETCH IS that hook.
    from delfin.agent import source_apis, web_tools
    seen = []
    monkeypatch.setattr(web_tools, "_check_url", lambda url: None)
    monkeypatch.setattr(web_tools, "_fetch_bytes",
                        lambda url, timeout_s, want_status=False,
                        headers=None: seen.append(url)
                        or (b'{"a": 1}', "application/json"))
    monkeypatch.setattr(source_apis, "_GUARDED_FETCH",
                        source_apis._fallback_fetch)
    assert source_apis.fetch_json("https://api.crossref.org/works/x") \
        == {"a": 1}
    assert seen == ["https://api.crossref.org/works/x"]
    monkeypatch.setattr(web_tools, "_check_url",
                        lambda url: "deny-listed")
    with pytest.raises(ValueError, match="refused"):
        source_apis.fetch_json("https://api.crossref.org/works/y")


def test_licence_and_terms_are_documented_in_the_docstring():
    import delfin.agent.source_apis as sa
    doc = sa.__doc__ or ""
    for term in ("terms", "rate", "contact"):
        assert term in doc.lower()


def test_non_json_body_is_an_error_not_an_empty_record(monkeypatch):
    from delfin.agent import source_apis, web_tools
    monkeypatch.setattr(web_tools, "_check_url", lambda url: None)
    monkeypatch.setattr(web_tools, "_fetch_bytes",
                        lambda url, timeout_s, want_status=False,
                        headers=None: (b"<html>challenge</html>",
                                       "text/html"))
    monkeypatch.setattr(source_apis, "_GUARDED_FETCH",
                        source_apis._fallback_fetch)
    with pytest.raises(ValueError, match="non-JSON"):
        source_apis.fetch_json("https://api.crossref.org/works/x")


def test_cache_returns_what_the_api_returned_nothing_else(transport):
    from delfin.agent.source_apis import verify_doi
    first = verify_doi("10.1021/acs.jpclett.5b00331")
    second = verify_doi("10.1021/acs.jpclett.5b00331")
    assert first == second
    # only ONE host call: the second came from the cache
    assert len(transport.urls) == 1


def test_pubchem_routes_are_chosen_by_identifier_shape(transport):
    # name / SMILES / InChIKey are three separate PUG-REST routes; the
    # SMILES value must be percent-quoted in the URL it produced.
    from delfin.agent.source_apis import lookup_compound
    rec = lookup_compound("C([H])(=O)/N")
    assert transport.urls and "/compound/smiles/" in transport.urls[-1]
    assert "%2F" in transport.urls[-1]  # the slash was quoted
    assert rec["status"] == "ok"
    assert not any("/compound/name/" in u and "=O" in u
                   for u in transport.urls)
