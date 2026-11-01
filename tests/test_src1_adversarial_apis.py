"""Adversarial tests for SRC1 P1 (reviewer nacht-s4).

Recorded responses only -- no live network. Attacks the shapes and
transports of delfin/agent/source_apis.py.
"""
from __future__ import annotations

import json
import urllib.error

import pytest

from delfin.agent import source_apis, web_tools


class FakeTransport:
    def __init__(self, routes):
        self.routes = routes
        self.calls = []

    def fetch_json(self, url, timeout_s=30):
        self.calls.append(url)
        for prefix, spec in self.routes.items():
            if url.startswith(prefix):
                if spec.get("code") == 404:
                    raise urllib.error.HTTPError(
                        url, 404, "not found", None, None)
                if spec.get("code"):
                    raise urllib.error.HTTPError(
                        url, spec["code"], "err", None, None)
                if spec.get("raw") is not None:
                    raise ValueError("non-JSON body: x")
                return spec["json"]
        raise AssertionError(f"unexpected url {url}")


@pytest.fixture()
def transport(monkeypatch):
    fake = FakeTransport({})
    monkeypatch.setattr(source_apis, "fetch_json", fake.fetch_json)
    # DNS resolution is refused in the cage; the guard is exercised
    # separately, the behaviour under test here is the module's own.
    monkeypatch.setattr(web_tools, "_check_url", lambda url: None)
    source_apis.clear_cache()
    return fake


def test_doi_wrappers_and_case_hit_one_key(transport):
    transport.routes["https://api.crossref.org/works/10.1039%2Fc8sc02158a"] = {
        "json": {"message": {"title": ["Real Paper"], "DOI": "10.1039/C8SC02158A"}}}
    for form in ("10.1039/C8SC02158A", "doi:10.1039/C8SC02158A",
                 "https://doi.org/10.1039/C8SC02158A",
                 "  10.1039/c8sc02158a  "):
        rec = source_apis.verify_doi(form)
        assert rec["status"] == "ok", form
        assert rec["doi"] == "10.1039/c8sc02158a"
    assert len(transport.calls) == 1, "wrappers must share one cache key"


def test_invented_wellformed_doi_is_not_found_not_echoed(transport):
    transport.routes["https://api.crossref.org/works/"] = {"code": 404}
    rec = source_apis.verify_doi("10.9999/INVENTED-PAPER-1234")
    assert rec["status"] == "not_found"
    for field in ("title", "authors", "journal", "year"):
        assert not rec.get(field), f"invented DOI echoed as {field}"


@pytest.mark.parametrize("bad", [None, "", "   ", "10.", "10.abc/x",
                                 "crossref.org/works/10.1/x", "10.1234567890/x"])
def test_malformed_doi_is_invalid_input(transport, bad):
    rec = source_apis.verify_doi(bad)
    assert rec["status"] == "invalid_input"
    assert transport.calls == []


def test_unicode_suffix_is_quoted_never_ok_from_nothing(transport):
    # Cyrillic look-alike suffix: must be percent-quoted (no control
    # characters in the request line) and come back as not_found, never
    # as a fabricated ok record.
    transport.routes["https://api.crossref.org/works/"] = {"code": 404}
    rec = source_apis.verify_doi("10.9999/официальный-дою")
    assert rec["status"] == "not_found"
    sent = transport.calls[0]
    assert "%" in sent, f"unicode suffix sent unquoted: {sent}"


def test_http_500_is_error_and_not_cached(transport):
    transport.routes["https://api.crossref.org/works/10.1039%2Fx"] = {"code": 503}
    first = source_apis.verify_doi("10.1039/x")
    assert first["status"] == "error" and "503" in first["error"]
    second = source_apis.verify_doi("10.1039/x")
    assert len(transport.calls) == 2, "an error must not poison the cache"


def test_pubchem_routes_by_shape(transport):
    transport.routes = {
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/":
            {"json": {"PropertyTable": {"Properties": [{"CID": 2519,
             "MolecularFormula": "C8H10N4O2", "InChIKey":
             "RYYVLZVUVIJVGH-UHFFFAOYSA-N", "MolecularWeight": 194.19}]}}},
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/":
            {"json": {"PropertyTable": {"Properties": [{"CID": 241}]}}},
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/":
            {"json": {"PropertyTable": {"Properties": [{"CID": 2244}]}}},
    }
    assert source_apis.lookup_compound("caffeine")["cid"] == 2519
    assert source_apis.lookup_compound("CN1C=NC2=C1C=O")["cid"] == 241
    assert source_apis.lookup_compound(
        "RYYVLZVUVIJVGH-UHFFFAOYSA-N")["cid"] == 2244
    smiles_call = [c for c in transport.calls if "/smiles/" in c]
    assert smiles_call, "SMILES shape must use the smiles route"
    assert "%" in smiles_call[0], "SMILES must be percent-quoted"


def test_pubchem_unknown_and_missing_fields(transport):
    transport.routes["https://pubchem.ncbi.nlm.nih.gov/"] = {"code": 404}
    rec = source_apis.lookup_compound("no-such-molecule-xyz")
    assert rec["status"] == "not_found"
    # A 200 with an empty Properties row must not pretend success.
    transport.routes["https://pubchem.ncbi.nlm.nih.gov/"] = {
        "json": {"PropertyTable": {"Properties": []}}}
    source_apis.clear_cache()
    rec = source_apis.lookup_compound("empty-table-xyz")
    assert rec["status"] == "ok" and rec["cid"] is None
    assert rec["formula"] is None


# --- host_for_action: the ONE mapping the security layer reads ---------------

def test_host_for_action_maps_known_and_never_guesses_unknown():
    # The rule layer (permission_rules._subjects) turns this value into
    # the host a web_fetch(domain:) rule may match. A wrong or guessed
    # host here would make the policy bind the WRONG host — so unknown
    # and hostile shapes must give None, never a fallback.
    from delfin.agent.source_apis import host_for_action
    assert host_for_action("doi") == "api.crossref.org"
    assert host_for_action("compound") == "pubchem.ncbi.nlm.nih.gov"
    assert host_for_action("works") == "api.openalex.org"
    # Case variants normalise to the SAME fixed host (closed mapping —
    # every result is one of the three documented API hosts).
    assert host_for_action("DOI") == "api.crossref.org"
    for hostile in (None, "", "  ", "unknown", "scrape",
                    "doi; rm -rf", {"action": "doi"}, ["doi"]):
        assert host_for_action(hostile) is None, f"{hostile!r} must not map"
