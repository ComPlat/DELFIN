"""Federated scholarly search: adapters, merge, cache, rate limit, errors.

No test touches the network: every source answers from a recorded,
synthetic response through a fake ``scholarly_sources._get``. The shapes
are the services' real response shapes (OpenAlex /works, Crossref
/works, Semantic Scholar /paper/search, arXiv Atom, Europe PMC /search,
Unpaywall /v2/<doi>); the contents are invented.
"""

from __future__ import annotations

import json
import threading
import time
import types
import urllib.error
import urllib.parse

import pytest

from delfin.agent import scholarly_sources as ss
from delfin.agent import source_apis, web_tools

OPENALEX = {"results": [
    {"id": "https://openalex.org/W1", "doi": "https://doi.org/10.1063/1.464913",
     "title": "Density-functional thermochemistry. III. The role of exact exchange",
     "publication_year": 1993,
     "authorships": [{"author": {"display_name": "Axel D. Becke"}}],
     "primary_location": {"source": {"display_name": "The Journal of Chemical Physics"},
                          "landing_page_url": "https://doi.org/10.1063/1.464913"},
     "best_oa_location": None, "open_access": {"is_oa": False},
     "cited_by_count": 100000, "type": "article"},
    {"id": "https://openalex.org/W2", "doi": "https://doi.org/10.1039/test.2",
     "title": "An open benchmark of spin-state energies",
     "publication_year": 2022,
     "authorships": [{"author": {"display_name": "Jane Doe"}},
                     {"author": {"display_name": "Rick Roe"}}],
     "primary_location": {"source": {"display_name": "Chemical Science"}},
     "best_oa_location": {"pdf_url": "https://pubs.example.org/oa/test2.pdf",
                          "license": "cc-by", "is_oa": True},
     "open_access": {"is_oa": True, "oa_status": "gold"},
     "cited_by_count": 12, "type": "article"},
]}

CROSSREF = {"status": "ok", "message": {"items": [
    {"DOI": "10.1063/1.464913",
     "title": ["Density-functional thermochemistry. III. The role of exact exchange"],
     "author": [{"given": "Axel D.", "family": "Becke"}],
     "container-title": ["The Journal of Chemical Physics"],
     "issued": {"date-parts": [[1993, 4, 1]]}, "type": "journal-article",
     "is-referenced-by-count": 90000},
    {"DOI": "10.26434/chemrxiv-2024-abcde",
     "title": ["A preprint on iron spin crossover"],
     "author": [{"given": "Ann", "family": "Smith"}],
     "institution": [{"name": "ChemRxiv"}],
     "issued": {"date-parts": [[2024]]}, "type": "posted-content",
     "license": [{"URL": "https://creativecommons.org/licenses/by-nc/4.0/"}]},
]}}

S2 = {"total": 1, "data": [
    {"paperId": "abc123",
     "title": "Density-functional thermochemistry. III. The role of exact exchange",
     "year": 1993, "authors": [{"name": "A. Becke"}], "venue": "J. Chem. Phys.",
     "externalIds": {"DOI": "10.1063/1.464913"}, "citationCount": 95000,
     "isOpenAccess": False, "openAccessPdf": None},
]}

ARXIV = b"""<?xml version="1.0" encoding="UTF-8"?>
<feed xmlns="http://www.w3.org/2005/Atom" xmlns:arxiv="http://arxiv.org/schemas/atom">
  <title>arXiv Query</title>
  <entry>
    <id>http://arxiv.org/abs/2301.00001v2</id>
    <published>2023-01-01T00:00:00Z</published>
    <title>Spin-state energetics of
      iron complexes</title>
    <author><name>John Doe</name></author>
    <author><name>Mary Major</name></author>
    <link href="http://arxiv.org/abs/2301.00001v2" rel="alternate" type="text/html"/>
    <link title="pdf" href="http://arxiv.org/pdf/2301.00001v2" rel="related" type="application/pdf"/>
  </entry>
</feed>
"""

EPMC = {"hitCount": 1, "resultList": {"result": [
    {"id": "123", "source": "MED", "pmid": "123", "pmcid": "PMC1",
     "title": "Spin-state energetics of iron complexes.",
     "authorString": "Doe J, Major M.",
     "journalInfo": {"journal": {"title": "J Comput Chem"}},
     "pubYear": "2024", "doi": "10.1000/epmc.1", "citedByCount": 5,
     "isOpenAccess": "Y", "license": "cc by",
     "fullTextUrlList": {"fullTextUrl": [
         {"availabilityCode": "OA", "documentStyle": "html", "url": "https://europepmc.org/articles/PMC1"},
         {"availabilityCode": "OA", "documentStyle": "pdf",
          "url": "https://europepmc.org/articles/PMC1/pdf"}]}},
]}}

UNPAYWALL = {"doi": "10.1063/1.464913", "is_oa": True,
             "best_oa_location": {"url_for_pdf": "https://repo.example.edu/becke.pdf",
                                  "url": "https://repo.example.edu/becke",
                                  "license": "cc-by", "version": "acceptedVersion"}}

HOST_BODY = {
    "api.openalex.org": json.dumps(OPENALEX).encode(),
    "api.crossref.org": json.dumps(CROSSREF).encode(),
    "api.semanticscholar.org": json.dumps(S2).encode(),
    "export.arxiv.org": ARXIV,
    "www.ebi.ac.uk": json.dumps(EPMC).encode(),
    "api.unpaywall.org": json.dumps(UNPAYWALL).encode(),
}


class FakeNet:
    def __init__(self):
        self.calls = []          # (url, headers)
        self.fail = {}           # host -> exception
        self.block = {}          # host -> threading.Event

    def __call__(self, url, *, timeout_s, headers=None):
        host = urllib.parse.urlsplit(url).hostname
        self.calls.append((url, dict(headers or {})))
        if host in self.block:
            self.block[host].wait(5)
        if host in self.fail:
            raise self.fail[host]
        return HOST_BODY[host]

    def hosts(self):
        return [urllib.parse.urlsplit(u).hostname for u, _ in self.calls]


@pytest.fixture
def net(monkeypatch):
    fake = FakeNet()
    monkeypatch.setattr(ss, "_get", fake)
    monkeypatch.setattr(ss, "contact_email", lambda: "")
    monkeypatch.setattr(ss, "_credential", lambda name: "")
    monkeypatch.setattr(web_tools, "outbound_secret_reason", lambda text: None)
    source_apis.clear_cache()
    ss._COOLDOWN.clear()
    yield fake
    for ev in fake.block.values():
        ev.set()
    source_apis.clear_cache()
    ss._COOLDOWN.clear()


# --- parsing, per adapter ----------------------------------------------------

def test_openalex_records_have_the_scholar_fields(net):
    recs = ss.search_openalex("exact exchange", 5)
    first, second = recs
    assert first["title"].startswith("Density-functional thermochemistry")
    assert first["authors"] == ["Axel D. Becke"]
    assert first["year"] == 1993 and first["venue"] == "The Journal of Chemical Physics"
    assert first["doi"] == "10.1063/1.464913"
    assert first["url"] == "https://doi.org/10.1063/1.464913"
    assert first["cited_by"] == 100000 and first["pdf_url"] == ""
    assert second["pdf_url"] == "https://pubs.example.org/oa/test2.pdf"
    assert second["pdf_license"] == "cc-by" and second["is_oa"] is True
    url = net.calls[0][0]
    q = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
    assert q["search"] == ["exact exchange"] and "mailto" not in q


def test_crossref_labels_chemrxiv_preprints(net):
    recs = ss.search_crossref("spin crossover", 5)
    assert recs[0]["year"] == 1993 and recs[0]["authors"] == ["Axel D. Becke"]
    pre = recs[1]
    assert pre["preprint"] == "ChemRxiv" and pre["venue"] == "ChemRxiv"
    assert pre["ids"]["crossref_license"] == "cc-by-nc"
    assert pre["pdf_url"] == ""      # Crossref never claims a PDF


def test_chemrxiv_source_filters_on_its_doi_prefix(net):
    recs = ss.search_chemrxiv("spin", 5)
    q = urllib.parse.parse_qs(urllib.parse.urlsplit(net.calls[0][0]).query)
    assert q["filter"] == ["prefix:10.26434"]
    assert all(r["sources"] == ["chemrxiv"] for r in recs)


def test_semantic_scholar_parses_and_sends_key_in_header_only(net, monkeypatch):
    monkeypatch.setattr(ss, "_credential",
                        lambda name: "s2-secret-key-123456" if name == "SEMANTIC_SCHOLAR_API_KEY" else "")
    recs = ss.search_semantic_scholar("exact exchange", 3)
    assert recs[0]["doi"] == "10.1063/1.464913" and recs[0]["cited_by"] == 95000
    url, headers = net.calls[0]
    assert headers.get("x-api-key") == "s2-secret-key-123456"
    assert "s2-secret-key" not in url


def test_openalex_key_travels_as_bearer_header(net, monkeypatch):
    monkeypatch.setattr(ss, "_credential",
                        lambda name: "oa-key-abcdefghijkl" if name == "OPENALEX_API_KEY" else "")
    ss.search_openalex("x", 1)
    url, headers = net.calls[0]
    assert headers["Authorization"] == "Bearer oa-key-abcdefghijkl"
    assert "oa-key" not in url and "api_key" not in url


def test_arxiv_atom_is_parsed_and_pdf_linked_without_licence(net):
    recs = ss.search_arxiv("spin state iron", 3)
    (rec,) = recs
    assert rec["title"] == "Spin-state energetics of iron complexes"
    assert rec["authors"] == ["John Doe", "Mary Major"] and rec["year"] == 2023
    assert rec["ids"]["arxiv"] == "2301.00001" and rec["preprint"] == "arXiv"
    assert rec["pdf_url"].startswith("https://arxiv.org/pdf/")
    assert rec["pdf_license"] == ""
    q = urllib.parse.parse_qs(urllib.parse.urlsplit(net.calls[0][0]).query)
    assert q["search_query"] == ["all:spin AND all:state AND all:iron"]


def test_arxiv_feed_with_entities_is_refused():
    with pytest.raises(ss.SourceError):
        ss.parse_arxiv_atom(b'<?xml version="1.0"?><!DOCTYPE x [<!ENTITY a "b">]><feed/>')


def test_europe_pmc_takes_the_oa_pdf_and_licence(net):
    (rec,) = ss.search_europe_pmc("spin", 3)
    assert rec["doi"] == "10.1000/epmc.1" and rec["venue"] == "J Comput Chem"
    assert rec["authors"] == ["Doe J", "Major M"]
    assert rec["pdf_url"] == "https://europepmc.org/articles/PMC1/pdf"
    assert rec["pdf_license"] == "cc-by" and rec["ids"]["pmcid"] == "PMC1"


@pytest.mark.parametrize("raw,code", [
    ("https://creativecommons.org/licenses/by/4.0/", "cc-by"),
    ("http://creativecommons.org/publicdomain/zero/1.0/", "cc0"),
    ("CCBY", "cc-by"), ("cc by-nc-nd", "cc-by-nc-nd"), ("cc-by-sa", "cc-by-sa"),
    ("pd", "public-domain"), ("publisher-specific-oa", "publisher-specific-oa"),
    ("", ""), (None, ""),
])
def test_licences_are_normalised(raw, code):
    assert ss.normalise_license(raw) == code


def test_only_open_licences_permit_storing_a_copy():
    assert ss.licence_permits_storage("cc-by")
    assert ss.licence_permits_storage("https://creativecommons.org/licenses/by-nc/4.0/")
    for lic in ("", "implied-oa", "publisher-specific-oa", "other-oa", None):
        assert not ss.licence_permits_storage(lic)


# --- federation, merge, de-dup -----------------------------------------------------

ALL = ["openalex", "crossref", "semantic_scholar", "arxiv", "europe_pmc"]


def test_federated_search_merges_and_deduplicates(net):
    out = ss.federated_search("spin", sources=ALL, unpaywall=False, timeout_s=5)
    assert out["status"] == {k: "ok" for k in ALL}
    recs = out["results"]
    becke = [r for r in recs if r["doi"] == "10.1063/1.464913"]
    assert len(becke) == 1
    assert becke[0]["sources"] == ["openalex", "crossref", "semantic_scholar"]
    assert becke[0]["cited_by"] == 100000
    spin = [r for r in recs if r["title"].startswith("Spin-state energetics")]
    assert len(spin) == 1, "arXiv preprint and Europe PMC paper are one work"
    assert set(spin[0]["sources"]) == {"europe_pmc", "arxiv"}
    assert spin[0]["pdf_license"] == "cc-by"        # the storable copy wins
    assert spin[0]["preprint"] == ""                # the published version wins
    # Becke was found by three sources: fused relevance ranks it first.
    assert recs[0]["doi"] == "10.1063/1.464913"


def test_sorting_by_year_and_citations(net):
    out = ss.federated_search("spin", sources=ALL, unpaywall=False)
    years = [r["year"] for r in ss.sort_records(out["results"], "year")]
    assert years == sorted(years, reverse=True)
    cites = [r["cited_by"] or 0 for r in ss.sort_records(out["results"], "citations")]
    assert cites == sorted(cites, reverse=True)


def test_title_dedup_needs_close_years():
    a = dict(ss._empty_record("arxiv"), title="Same Title", year=2001)
    b = dict(ss._empty_record("europe_pmc"), title="Same title.", year=2020)
    merged = ss.merge_records([("arxiv", [a]), ("europe_pmc", [b])])
    assert len(merged) == 2


def test_results_are_cached_per_source(net):
    ss.federated_search("spin", sources=["openalex", "crossref"], unpaywall=False)
    n = len(net.calls)
    again = ss.federated_search("SPIN", sources=["openalex", "crossref"], unpaywall=False)
    assert len(net.calls) == n == 2
    assert again["status"] == {"openalex": "ok", "crossref": "ok"}


def test_a_failing_source_does_not_hide_the_others(net):
    net.fail["export.arxiv.org"] = urllib.error.URLError("Name or service not known")
    out = ss.federated_search("spin", sources=["openalex", "arxiv"], unpaywall=False)
    assert out["status"]["openalex"] == "ok"
    assert out["status"]["arxiv"].startswith("error: network error")
    assert out["results"]
    line = ss.status_line(out["status"])
    assert line.startswith("arXiv: network error") and "\n" not in line


def test_a_slow_source_times_out_and_partial_results_show(net):
    net.block["www.ebi.ac.uk"] = threading.Event()
    t0 = time.monotonic()
    out = ss.federated_search("slow spin", sources=["openalex", "europe_pmc"],
                              unpaywall=False, timeout_s=0.5)
    assert time.monotonic() - t0 < 3
    assert out["status"] == {"openalex": "ok", "europe_pmc": "timed out"}
    assert len(out["results"]) == 2


def test_rate_limited_source_pauses_without_asking_again(net):
    net.fail["api.semanticscholar.org"] = urllib.error.HTTPError(
        "u", 429, "Too Many Requests", {}, None)
    out = ss.federated_search("spin", sources=["semantic_scholar"], unpaywall=False)
    assert "429" in out["status"]["semantic_scholar"]
    n = len(net.calls)
    out = ss.federated_search("other", sources=["semantic_scholar"], unpaywall=False)
    assert len(net.calls) == n
    assert "paused" in out["status"]["semantic_scholar"]


def test_errors_are_not_cached(net):
    net.fail["api.openalex.org"] = urllib.error.URLError("down")
    ss.federated_search("spin", sources=["openalex"], unpaywall=False)
    del net.fail["api.openalex.org"]
    out = ss.federated_search("spin", sources=["openalex"], unpaywall=False)
    assert out["status"] == {"openalex": "ok"} and out["results"]


def test_empty_and_secret_queries_send_nothing(net, monkeypatch):
    assert ss.federated_search("   ")["results"] == []
    monkeypatch.setattr(web_tools, "outbound_secret_reason",
                        lambda text: "it contains a credential DELFIN holds")
    out = ss.federated_search("sk-abc", sources=["openalex"])
    assert out["results"] == [] and "Not sent" in out["notes"][0]
    assert net.calls == []


# --- Unpaywall and the contact e-mail ------------------------------------------------

def test_unpaywall_is_skipped_without_a_contact_email(net):
    out = ss.federated_search("spin", sources=["openalex"], unpaywall=True)
    assert "unpaywall" not in out["status"]
    assert any(n.startswith("Unpaywall skipped") for n in out["notes"])
    assert "api.unpaywall.org" not in net.hosts()


def test_unpaywall_adds_the_oa_pdf_with_the_configured_address(net, monkeypatch):
    monkeypatch.setattr(ss, "contact_email", lambda: "someone@example.org")
    out = ss.federated_search("spin", sources=["openalex"], unpaywall=True)
    becke = next(r for r in out["results"] if r["doi"] == "10.1063/1.464913")
    assert becke["pdf_url"] == "https://repo.example.edu/becke.pdf"
    assert becke["pdf_license"] == "cc-by" and "unpaywall" in becke["sources"]
    up = [u for u, _ in net.calls if "unpaywall" in u]
    assert up and all("email=someone%40example.org" in u for u in up)
    oa_url = next(u for u, _ in net.calls if "openalex" in u)
    assert "mailto=someone%40example.org" in oa_url


def test_contact_email_comes_only_from_env_or_settings(monkeypatch, tmp_path):
    from delfin import user_settings
    settings = tmp_path / "settings.json"
    monkeypatch.setattr(user_settings, "get_settings_path", lambda *a, **k: settings)
    monkeypatch.delenv("DELFIN_CONTACT_EMAIL", raising=False)
    assert ss.contact_email() == ""
    settings.write_text(json.dumps({"literature": {"contact_email": "me@lab.example"}}))
    assert ss.contact_email() == "me@lab.example"
    monkeypatch.setenv("DELFIN_CONTACT_EMAIL", "env@lab.example")
    assert ss.contact_email() == "env@lab.example"
    monkeypatch.setenv("DELFIN_CONTACT_EMAIL", "not an address")
    assert ss.contact_email() == ""


def test_no_contact_address_is_built_into_the_code():
    import inspect
    for mod in (ss, web_tools, source_apis):
        src = inspect.getsource(mod)
        assert "complat.org" not in src
        assert "@gmail" not in src


# --- rate limiting -------------------------------------------------------------------

def test_rate_limit_spaces_requests_per_host_without_blocking_others(monkeypatch):
    clock = {"t": 1000.0}
    slept = []
    fake_time = types.SimpleNamespace(monotonic=lambda: clock["t"],
                                      sleep=lambda s: slept.append(round(s, 3)))
    monkeypatch.setattr(source_apis, "time", fake_time)
    monkeypatch.setattr(source_apis, "_LAST_REQUEST", {})
    source_apis._rate_limit("export.arxiv.org")
    source_apis._rate_limit("export.arxiv.org")
    assert slept == [3.0]                 # arXiv: one request every 3 s
    source_apis._rate_limit("api.openalex.org")
    assert slept == [3.0]                 # another host is not held up
    source_apis._rate_limit("export.arxiv.org")
    assert slept == [3.0, 6.0]            # the third waits for its own slot


def test_one_connection_at_a_time_where_the_terms_say_so():
    with source_apis._HostSlot("export.arxiv.org", 1.0):
        with pytest.raises(TimeoutError):
            with source_apis._HostSlot("export.arxiv.org", 0.1):
                pass
    with source_apis._HostSlot("export.arxiv.org", 0.1):
        pass                              # released again


def test_fetch_bytes_goes_through_the_guarded_stack(monkeypatch):
    seen = []
    monkeypatch.setattr(web_tools, "_check_url", lambda url: None)
    monkeypatch.setattr(web_tools, "_fetch_bytes",
                        lambda url, timeout_s, want_status=False, headers=None:
                        seen.append((url, headers)) or (b"<feed/>", "application/atom+xml"))
    monkeypatch.setattr(source_apis, "_LAST_REQUEST", {})
    assert source_apis.fetch_bytes("https://www.ebi.ac.uk/x",
                                   headers={"Accept": "application/json"}) == b"<feed/>"
    assert seen == [("https://www.ebi.ac.uk/x", {"Accept": "application/json"})]
    monkeypatch.setattr(web_tools, "_check_url", lambda url: "deny-listed")
    with pytest.raises(ValueError, match="refused"):
        source_apis.fetch_bytes("https://www.ebi.ac.uk/y")


# --- storing a copy --------------------------------------------------------------------

class _Resp:
    def __init__(self, body, ctype):
        self.body, self.headers = body, {"Content-Type": ctype}

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def read(self, n):
        return self.body[:n]


def _opener(body, ctype="application/pdf"):
    return types.SimpleNamespace(open=lambda req, timeout: _Resp(body, ctype))


def test_download_open_pdf_only_for_open_licences(monkeypatch):
    monkeypatch.setattr(web_tools, "_check_url", lambda url: None)
    monkeypatch.setattr(web_tools, "_guarded_opener", lambda: _opener(b"%PDF-1.7 x"))
    monkeypatch.setattr(source_apis, "_LAST_REQUEST", {})
    rec = {"pdf_url": "https://pubs.example.org/oa/a.pdf", "pdf_license": "cc-by"}
    assert ss.download_open_pdf(rec).startswith(b"%PDF-")
    with pytest.raises(ValueError, match="licence"):
        ss.download_open_pdf(dict(rec, pdf_license="publisher-specific-oa"))
    with pytest.raises(ValueError, match="arXiv"):
        ss.download_open_pdf({"pdf_url": "https://arxiv.org/pdf/1", "pdf_license": "cc-by"})
    with pytest.raises(ValueError, match="arXiv"):
        ss.download_open_pdf({"pdf_url": "https://export.arxiv.org/pdf/1", "pdf_license": "cc-by"})
    # A host that only ends in the same letters is not arXiv.
    assert ss.download_open_pdf({"pdf_url": "https://notarxiv.org/a.pdf",
                                 "pdf_license": "cc-by"}).startswith(b"%PDF-")


def test_download_refuses_a_landing_page(monkeypatch):
    monkeypatch.setattr(web_tools, "_check_url", lambda url: None)
    monkeypatch.setattr(web_tools, "_guarded_opener",
                        lambda: _opener(b"<html>challenge</html>", "text/html"))
    monkeypatch.setattr(source_apis, "_LAST_REQUEST", {})
    with pytest.raises(ValueError, match="rather than a PDF"):
        ss.download_open_pdf({"pdf_url": "https://p.example.org/a.pdf",
                              "pdf_license": "cc-by"})


def test_markdown_and_bibtex_records():
    rec = ss._openalex_record(OPENALEX["results"][1])
    md = ss.to_markdown(rec, retrieved="2026-11-11", pdf_file="Doe2022_An.pdf")
    assert "DOI: [10.1039/test.2](https://doi.org/10.1039/test.2)" in md
    assert "Retrieved: 2026-11-11" in md and "Found in: OpenAlex" in md
    bib = ss.to_bibtex(rec)
    assert bib.startswith("@article{Doe2022,") and "doi = {10.1039/test.2}" in bib
    assert ss.suggest_basename(rec) == "Doe2022_An_open_benchmark_of_spin_state"


def test_crossref_retries_without_select_on_http_400(net, monkeypatch):
    urls = []

    def get(url, *, timeout_s, headers=None):
        urls.append(url)
        if "select=" in url:
            raise urllib.error.HTTPError(url, 400, "Bad Request", {}, None)
        return HOST_BODY["api.crossref.org"]
    monkeypatch.setattr(ss, "_get", get)
    recs = ss.search_crossref("spin", 3)
    assert len(urls) == 2 and "select=" not in urls[1]
    assert recs[0]["doi"] == "10.1063/1.464913"


@pytest.mark.parametrize("author,stem", [
    ("Jane Doe", "Doe"), ("Doe J", "Doe"), ("Doe JA", "Doe"), ("Doe, Jane", "Doe"),
    ("Axel D. Becke", "Becke"), ("Ana Hernández", "Hernandez"),
])
def test_file_name_uses_the_surname(author, stem):
    rec = {"authors": [author], "year": 2020, "title": "A title"}
    assert ss.suggest_basename(rec) == f"{stem}2020_A_title"


def test_a_host_refusing_the_download_is_one_clear_line(monkeypatch):
    def refuse(req, timeout):
        raise urllib.error.HTTPError(req.full_url, 403, "Forbidden", {}, None)
    monkeypatch.setattr(web_tools, "_check_url", lambda url: None)
    monkeypatch.setattr(web_tools, "_guarded_opener",
                        lambda: types.SimpleNamespace(open=refuse))
    monkeypatch.setattr(source_apis, "_LAST_REQUEST", {})
    with pytest.raises(ValueError, match=r"refused the download \(HTTP 403\)"):
        ss.download_open_pdf({"pdf_url": "https://p.example.org/a.pdf",
                              "pdf_license": "cc-by"})
