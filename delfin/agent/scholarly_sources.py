"""Federated scholarly search over free, official APIs.

One query goes to several literature indexes at once, the way Zotero,
JabRef or Publish or Perish search: every source is a small adapter that
returns the same record shape, has its own rate limit (enforced per host
in ``source_apis``), timeout and cache, and fails on its own. The results
are merged and de-duplicated (by DOI, then by normalised title + year),
and each merged record says which sources found it.

All HTTP goes through ``source_apis.fetch_bytes`` -> ``web_tools``'
guarded stack (deny-list, SSRF checks, redirect re-validation). There is
no second opener; the PDF download uses the same guarded opener.

Terms, checked against each service's own documentation on 2026-11-11:

- **OpenAlex** -- https://help.openalex.org/how-to-use-the-api/rate-limits-and-authentication :
  key optional (``api_key`` or ``Authorization: Bearer``), 10 req/s
  without a key, 30 with a free key; a daily USD budget (keyless measured
  at $0.10/day, $0.001 per search). Data CC0. No polite pool any more in
  the docs; a ``mailto`` is sent only when configured.
- **Crossref** -- https://www.crossref.org/documentation/retrieve-metadata/rest-api/access-and-authentication/ :
  public pool 5 req/s, concurrency 1; polite pool (``mailto``) 10 req/s,
  concurrency 3. Metadata CC0. ChemRxiv preprints are Crossref DOIs with
  the prefix 10.26434 and are labelled, never scraped from chemrxiv.org.
- **Semantic Scholar Graph API** -- https://www.semanticscholar.org/product/api :
  without a key a pool of 1000 req/s shared by all keyless users (429s
  are common); with a key (``x-api-key`` header) 1 req/s.
- **arXiv API** -- https://info.arxiv.org/help/api/tou.html : one request
  every three seconds, one connection at a time; metadata CC0; e-prints
  may not be stored and served without the licence's permission, so an
  arXiv PDF is only ever linked here, never stored (the Atom feed carries
  no licence).
- **Europe PMC REST** (covers PubMed and PMC) --
  https://europepmc.org/RestfulWebService : free, no key; EBI states
  10 requests/s per client.
- **Unpaywall** -- https://unpaywall.org/products/api : ``email``
  parameter required, about 100 000 calls/day. Used only to find the
  open-access PDF of a DOI, and only when the person configured a
  contact e-mail (never a built-in address).
- **Google Scholar** has no API and its terms forbid automated querying,
  so it is not a source here.

Contact e-mail: ``DELFIN_CONTACT_EMAIL`` (environment) or the user
setting ``literature.contact_email``; optional. Optional API keys come
from DELFIN's credential store (env or credential file):
``OPENALEX_API_KEY``, ``SEMANTIC_SCHOLAR_API_KEY``; they travel in a
request header, never in a URL.
"""

from __future__ import annotations

import concurrent.futures as _cf
import copy
import datetime as _dt
import json
import os
import re
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Callable

from . import source_apis, web_tools

__all__ = [
    "SOURCES", "DEFAULT_SOURCES", "federated_search", "merge_records",
    "sort_records", "find_method_reference", "contact_email",
    "normalise_license", "licence_permits_storage", "unpaywall_lookup",
    "download_open_pdf", "to_bibtex", "to_markdown", "suggest_basename",
    "SourceError",
]

#: The record every adapter returns (and merge_records keeps).
RECORD_FIELDS = ("title", "authors", "year", "venue", "doi", "url",
                 "pdf_url", "pdf_license", "is_oa", "cited_by", "type",
                 "preprint", "ids", "sources")

_DEFAULT_TIMEOUT_S = 10.0
_CACHE_PREFIX = "fed:"
_COOLDOWN_S = 60.0          # after a 429, leave that source alone this long
_MAX_PDF_BYTES = 50 * 1024 * 1024


class SourceError(RuntimeError):
    """One source failed; the message is one line for the person."""


# --- configuration -------------------------------------------------------------

_EMAIL_RE = re.compile(r"^[^@\s<>\"']+@[^@\s<>\"']+\.[A-Za-z]{2,}$")


def contact_email() -> str:
    """The person's optional contact address for polite API use, or "".

    ``DELFIN_CONTACT_EMAIL`` wins; else ``literature.contact_email`` in the
    user settings file, read without writing it. Never a built-in value.
    """
    raw = os.environ.get("DELFIN_CONTACT_EMAIL", "").strip()
    if not raw:
        try:
            from delfin.user_settings import get_settings_path
            path = get_settings_path()
            if path.exists():
                data = json.loads(path.read_text(encoding="utf-8") or "{}")
                lit = data.get("literature") if isinstance(data, dict) else None
                if isinstance(lit, dict):
                    raw = str(lit.get("contact_email") or "").strip()
        except Exception:
            raw = ""
    return raw if _EMAIL_RE.match(raw) else ""


def _credential(name: str) -> str:
    return web_tools._credential(name)


# --- small helpers ---------------------------------------------------------------

def _get(url: str, *, timeout_s: float, headers: "dict | None" = None) -> bytes:
    """One guarded GET (tests replace this with recorded responses)."""
    return source_apis.fetch_bytes(url, timeout_s=timeout_s, headers=headers)


def _get_json(url: str, *, timeout_s: float, headers: "dict | None" = None):
    body = _get(url, timeout_s=timeout_s,
                headers={"Accept": "application/json", **(headers or {})})
    try:
        return json.loads(body.decode("utf-8", errors="replace"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise SourceError("answered with something that is not JSON") from exc


def _norm_doi(doi) -> str:
    text = str(doi or "").strip()
    text = re.sub(r"^(?:(?:https?://)?(?:dx\.)?doi\.org/|doi:)", "", text,
                  flags=re.IGNORECASE)
    text = text.strip().lower()
    return text if re.match(r"^10\.\d{4,9}/\S+$", text) else ""


def _http_url(value) -> str:
    text = str(value or "").strip()
    return text if text.startswith(("https://", "http://")) else ""


def _int_or_none(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _clean(text) -> str:
    return re.sub(r"\s+", " ", str(text or "")).strip()


def _empty_record(source: str) -> dict:
    return {"title": "", "authors": [], "year": None, "venue": "", "doi": "",
            "url": "", "pdf_url": "", "pdf_license": "", "is_oa": False,
            "cited_by": None, "type": "", "preprint": "", "ids": {},
            "sources": [source]}


def _finish(rec: dict) -> "dict | None":
    rec["title"] = _clean(rec.get("title"))
    if not rec["title"]:
        return None
    rec["doi"] = _norm_doi(rec.get("doi"))
    if rec["doi"] and not rec.get("url"):
        rec["url"] = "https://doi.org/" + rec["doi"]
    rec["url"] = _http_url(rec.get("url"))
    rec["pdf_url"] = _http_url(rec.get("pdf_url"))
    rec["pdf_license"] = normalise_license(rec.get("pdf_license"))
    if rec["doi"].startswith("10.26434/"):
        rec["preprint"] = "ChemRxiv"
        rec["venue"] = rec.get("venue") or "ChemRxiv"
    return rec


# --- licences --------------------------------------------------------------------

_LICENCE_COMPACT = {
    "ccby": "cc-by", "ccbysa": "cc-by-sa", "ccbync": "cc-by-nc",
    "ccbyncsa": "cc-by-nc-sa", "ccbynd": "cc-by-nd", "ccbyncnd": "cc-by-nc-nd",
    "cc0": "cc0", "pd": "public-domain", "publicdomain": "public-domain",
}

#: Licences under which keeping a verbatim copy in one's own literature
#: folder is allowed (every Creative Commons licence permits copying
#: verbatim; CC0 / public domain too). Anything else -- publisher-specific,
#: "implied OA", bronze, unknown -- is linked, never stored.
_STORABLE = frozenset(_LICENCE_COMPACT.values())


def normalise_license(raw) -> str:
    """A short licence code (``cc-by``, ``cc0`` ...) or the cleaned text."""
    text = str(raw or "").strip().lower()
    if not text:
        return ""
    m = re.search(r"creativecommons\.org/(licenses|publicdomain)/([a-z\-]+)", text)
    if m:
        if m.group(1) == "publicdomain":
            return "cc0" if m.group(2) == "zero" else "public-domain"
        return "cc-" + m.group(2)
    compact = re.sub(r"[^a-z0-9]", "", text)
    if compact in _LICENCE_COMPACT:
        return _LICENCE_COMPACT[compact]
    return re.sub(r"\s+", "-", text)[:60]


def licence_permits_storage(licence) -> bool:
    return normalise_license(licence) in _STORABLE


# --- adapters ----------------------------------------------------------------------

def _openalex_record(w: dict) -> "dict | None":
    rec = _empty_record("openalex")
    rec["title"] = w.get("title") or w.get("display_name") or ""
    rec["doi"] = w.get("doi") or ""
    rec["year"] = _int_or_none(w.get("publication_year"))
    rec["authors"] = [_clean((a.get("author") or {}).get("display_name"))
                      for a in (w.get("authorships") or [])
                      if isinstance(a, dict)]
    rec["authors"] = [a for a in rec["authors"] if a]
    prim = w.get("primary_location") or {}
    rec["venue"] = _clean((prim.get("source") or {}).get("display_name"))
    rec["url"] = "" if w.get("doi") else (prim.get("landing_page_url") or "")
    best = w.get("best_oa_location") or {}
    if best.get("pdf_url"):
        rec["pdf_url"] = best.get("pdf_url")
        rec["pdf_license"] = best.get("license") or ""
    oa = w.get("open_access") or {}
    rec["is_oa"] = bool(oa.get("is_oa"))
    rec["cited_by"] = _int_or_none(w.get("cited_by_count"))
    rec["type"] = str(w.get("type") or "")
    if rec["type"] == "preprint":
        rec["preprint"] = rec["venue"] or "preprint"
    if w.get("id"):
        rec["ids"]["openalex"] = str(w["id"])
    return _finish(rec)


_OPENALEX_SELECT = ("id,doi,title,display_name,publication_year,authorships,"
                    "primary_location,best_oa_location,open_access,"
                    "cited_by_count,type")


def _openalex_request(params: dict, timeout_s: float) -> list[dict]:
    params = dict(params)
    params.setdefault("select", _OPENALEX_SELECT)
    mail = contact_email()
    if mail:
        params["mailto"] = mail
    headers = {}
    key = _credential("OPENALEX_API_KEY")
    if key:
        headers["Authorization"] = f"Bearer {key}"
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
    data = _get_json(url, timeout_s=timeout_s, headers=headers)
    out = []
    for w in (data.get("results") or []) if isinstance(data, dict) else []:
        if isinstance(w, dict):
            rec = _openalex_record(w)
            if rec:
                out.append(rec)
    return out


def search_openalex(query: str, max_results: int = 10,
                    timeout_s: float = _DEFAULT_TIMEOUT_S) -> list[dict]:
    return _openalex_request({"search": query,
                              "per-page": max(1, min(int(max_results), 50))},
                             timeout_s)


def _crossref_record(item: dict) -> "dict | None":
    rec = _empty_record("crossref")
    rec["title"] = " ".join(str(t) for t in (item.get("title") or []))
    rec["doi"] = item.get("DOI") or ""
    rec["authors"] = [
        _clean(" ".join(x for x in (a.get("given"), a.get("family")) if x)
               or a.get("name"))
        for a in (item.get("author") or []) if isinstance(a, dict)]
    rec["authors"] = [a for a in rec["authors"] if a]
    rec["venue"] = _clean(" ".join(str(t) for t in
                                   (item.get("container-title") or [])))
    if not rec["venue"]:
        inst = item.get("institution") or []
        if inst and isinstance(inst[0], dict):
            rec["venue"] = _clean(inst[0].get("name"))
    parts = ((item.get("issued") or {}).get("date-parts") or [[]])[0] or []
    rec["year"] = _int_or_none(parts[0]) if parts else None
    rec["cited_by"] = _int_or_none(item.get("is-referenced-by-count"))
    rec["type"] = str(item.get("type") or "")
    if rec["type"] == "posted-content":
        rec["preprint"] = rec["venue"] or "preprint"
    # Crossref's licence applies to the work; Crossref gives no OA PDF
    # link we could trust, so it is recorded but no PDF is claimed.
    lic = [l.get("URL") for l in (item.get("license") or [])
           if isinstance(l, dict) and l.get("URL")]
    if lic:
        rec["ids"]["crossref_license"] = normalise_license(lic[0])
    return _finish(rec)


def search_crossref(query: str, max_results: int = 10,
                    timeout_s: float = _DEFAULT_TIMEOUT_S,
                    prefix: str = "") -> list[dict]:
    params = {"query": query, "rows": max(1, min(int(max_results), 50)),
              "select": ("DOI,title,author,container-title,issued,type,"
                         "is-referenced-by-count,license")}
    if prefix:
        params["filter"] = f"prefix:{prefix}"
    mail = contact_email()
    if mail:
        params["mailto"] = mail
    url = "https://api.crossref.org/works?" + urllib.parse.urlencode(params)
    try:
        data = _get_json(url, timeout_s=timeout_s)
    except urllib.error.HTTPError as exc:
        if exc.code != 400:
            raise
        # A field list the API no longer accepts must not cost the whole
        # source: once more without `select` (full records, same parser).
        params.pop("select")
        url = "https://api.crossref.org/works?" + urllib.parse.urlencode(params)
        data = _get_json(url, timeout_s=timeout_s)
    items = ((data.get("message") or {}).get("items") or []) \
        if isinstance(data, dict) else []
    return [r for r in (_crossref_record(i) for i in items
                        if isinstance(i, dict)) if r]


def search_chemrxiv(query: str, max_results: int = 10,
                    timeout_s: float = _DEFAULT_TIMEOUT_S) -> list[dict]:
    """ChemRxiv preprints: Crossref restricted to ChemRxiv's DOI prefix."""
    out = search_crossref(query, max_results, timeout_s, prefix="10.26434")
    for rec in out:
        rec["sources"] = ["chemrxiv"]
    return out


def _s2_record(p: dict) -> "dict | None":
    rec = _empty_record("semantic_scholar")
    rec["title"] = p.get("title") or ""
    ext = p.get("externalIds") or {}
    rec["doi"] = ext.get("DOI") or ""
    if ext.get("ArXiv"):
        rec["ids"]["arxiv"] = str(ext["ArXiv"])
    if p.get("paperId"):
        rec["ids"]["s2"] = str(p["paperId"])
    rec["year"] = _int_or_none(p.get("year"))
    rec["authors"] = [_clean(a.get("name")) for a in (p.get("authors") or [])
                      if isinstance(a, dict) and a.get("name")]
    rec["venue"] = _clean(p.get("venue"))
    rec["cited_by"] = _int_or_none(p.get("citationCount"))
    rec["is_oa"] = bool(p.get("isOpenAccess"))
    pdf = p.get("openAccessPdf") or {}
    if isinstance(pdf, dict) and pdf.get("url"):
        rec["pdf_url"] = pdf.get("url")
        rec["pdf_license"] = pdf.get("license") or ""
    rec["url"] = "" if rec["doi"] else (p.get("url") or "")
    return _finish(rec)


def search_semantic_scholar(query: str, max_results: int = 10,
                            timeout_s: float = _DEFAULT_TIMEOUT_S) -> list[dict]:
    params = {"query": query, "limit": max(1, min(int(max_results), 100)),
              "fields": ("title,year,authors,venue,externalIds,citationCount,"
                         "openAccessPdf,isOpenAccess,url")}
    headers = {}
    key = _credential("SEMANTIC_SCHOLAR_API_KEY")
    if key:
        headers["x-api-key"] = key
    url = ("https://api.semanticscholar.org/graph/v1/paper/search?"
           + urllib.parse.urlencode(params))
    data = _get_json(url, timeout_s=timeout_s, headers=headers)
    rows = (data.get("data") or []) if isinstance(data, dict) else []
    return [r for r in (_s2_record(p) for p in rows if isinstance(p, dict)) if r]


_ATOM = "{http://www.w3.org/2005/Atom}"
_ARXIV = "{http://arxiv.org/schemas/atom}"


def _arxiv_query(query: str) -> str:
    terms = re.findall(r"[\w\-]+", query, flags=re.UNICODE)[:8]
    return " AND ".join(f"all:{t}" for t in terms) or "all:" + query[:50]


def parse_arxiv_atom(body: bytes) -> list[dict]:
    if b"<!ENTITY" in body or b"<!DOCTYPE" in body:
        raise SourceError("answered with a feed that declares entities")
    try:
        root = ET.fromstring(body)
    except ET.ParseError as exc:
        raise SourceError("answered with something that is not Atom") from exc
    out = []
    for entry in root.findall(f"{_ATOM}entry"):
        ident = (entry.findtext(f"{_ATOM}id") or "").strip()
        if "/api/errors" in ident:
            raise SourceError("rejected the query: "
                              + _clean(entry.findtext(f"{_ATOM}summary"))[:120])
        rec = _empty_record("arxiv")
        rec["title"] = entry.findtext(f"{_ATOM}title") or ""
        rec["authors"] = [_clean(a.findtext(f"{_ATOM}name"))
                          for a in entry.findall(f"{_ATOM}author")]
        rec["authors"] = [a for a in rec["authors"] if a]
        published = entry.findtext(f"{_ATOM}published") or ""
        rec["year"] = _int_or_none(published[:4])
        rec["doi"] = entry.findtext(f"{_ARXIV}doi") or ""
        journal = _clean(entry.findtext(f"{_ARXIV}journal_ref"))
        rec["venue"] = journal or "arXiv"
        rec["preprint"] = "" if journal else "arXiv"
        rec["type"] = "article" if journal else "preprint"
        rec["url"] = ident.replace("http://", "https://", 1)
        m = re.search(r"arxiv\.org/abs/([^\s]+?)(?:v\d+)?$", ident)
        if m:
            rec["ids"]["arxiv"] = m.group(1)
        for link in entry.findall(f"{_ATOM}link"):
            if link.get("title") == "pdf" and link.get("href"):
                # Linked only: the feed states no licence (see module doc).
                rec["pdf_url"] = link.get("href").replace("http://", "https://", 1)
        rec["is_oa"] = True
        rec = _finish(rec)
        if rec:
            out.append(rec)
    return out


def search_arxiv(query: str, max_results: int = 10,
                 timeout_s: float = _DEFAULT_TIMEOUT_S) -> list[dict]:
    params = {"search_query": _arxiv_query(query), "start": 0,
              "max_results": max(1, min(int(max_results), 50))}
    url = "https://export.arxiv.org/api/query?" + urllib.parse.urlencode(params)
    body = _get(url, timeout_s=timeout_s,
                headers={"Accept": "application/atom+xml"})
    return parse_arxiv_atom(body)


def _epmc_record(r: dict) -> "dict | None":
    rec = _empty_record("europe_pmc")
    rec["title"] = r.get("title") or ""
    rec["doi"] = r.get("doi") or ""
    authors = str(r.get("authorString") or "").rstrip(".")
    rec["authors"] = [a.strip() for a in authors.split(",") if a.strip()]
    journal = ((r.get("journalInfo") or {}).get("journal") or {}).get("title")
    rec["venue"] = _clean(journal or r.get("journalTitle") or "")
    rec["year"] = _int_or_none(r.get("pubYear"))
    rec["cited_by"] = _int_or_none(r.get("citedByCount"))
    rec["is_oa"] = str(r.get("isOpenAccess") or "").upper() == "Y"
    for key in ("pmid", "pmcid"):
        if r.get(key):
            rec["ids"][key] = str(r[key])
    urls = ((r.get("fullTextUrlList") or {}).get("fullTextUrl") or [])
    for u in urls:
        if (isinstance(u, dict) and u.get("availabilityCode") == "OA"
                and u.get("documentStyle") == "pdf" and u.get("url")):
            rec["pdf_url"] = u["url"]
            rec["pdf_license"] = r.get("license") or ""
            break
    if not rec["doi"]:
        src, ident = r.get("source"), r.get("id")
        if src and ident:
            rec["url"] = (f"https://europepmc.org/article/"
                          f"{urllib.parse.quote(str(src))}/"
                          f"{urllib.parse.quote(str(ident))}")
    if str(r.get("source") or "") == "PPR":
        rec["preprint"] = "preprint"
    return _finish(rec)


def search_europe_pmc(query: str, max_results: int = 10,
                      timeout_s: float = _DEFAULT_TIMEOUT_S) -> list[dict]:
    params = {"query": query, "format": "json", "resultType": "core",
              "pageSize": max(1, min(int(max_results), 25))}
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
           + urllib.parse.urlencode(params))
    data = _get_json(url, timeout_s=timeout_s)
    rows = ((data.get("resultList") or {}).get("result") or []) \
        if isinstance(data, dict) else []
    return [r for r in (_epmc_record(x) for x in rows if isinstance(x, dict)) if r]


def unpaywall_lookup(doi: str, *, timeout_s: float = _DEFAULT_TIMEOUT_S) -> dict:
    """Open-access location of a DOI from Unpaywall.

    ``{"status": "ok", "pdf_url", "license", "is_oa", "url"}``,
    ``{"status": "skipped", "reason"}`` without a contact e-mail, or
    ``{"status": "not_found"|"error", ...}``. Cached like the searches.
    """
    norm = _norm_doi(doi)
    if not norm:
        return {"status": "invalid_input", "reason": "not a DOI"}
    mail = contact_email()
    if not mail:
        return {"status": "skipped", "reason": _UNPAYWALL_SKIP}
    key = _CACHE_PREFIX + "unpaywall:" + norm
    hit = source_apis._cached(key)
    if hit is not None:
        return dict(hit)
    url = (f"https://api.unpaywall.org/v2/{urllib.parse.quote(norm, safe='/')}"
           f"?{urllib.parse.urlencode({'email': mail})}")
    try:
        data = _get_json(url, timeout_s=timeout_s)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            rec = {"status": "not_found"}
            source_apis._store(key, rec)
            return dict(rec)
        return {"status": "error", "reason": f"HTTP {exc.code}"}
    except Exception as exc:
        return {"status": "error", "reason": _one_line(exc)}
    best = (data.get("best_oa_location") or {}) if isinstance(data, dict) else {}
    rec = {"status": "ok", "is_oa": bool(data.get("is_oa")),
           "pdf_url": _http_url(best.get("url_for_pdf")),
           "url": _http_url(best.get("url")),
           "license": normalise_license(best.get("license"))}
    source_apis._store(key, rec)
    return dict(rec)


_UNPAYWALL_SKIP = ("Unpaywall skipped: it needs a contact e-mail; set "
                   "DELFIN_CONTACT_EMAIL or the setting literature.contact_email")


@dataclass(frozen=True)
class Source:
    key: str
    label: str
    search: Callable[..., list]
    default: bool
    host: str


SOURCES: dict[str, Source] = {
    s.key: s for s in (
        Source("openalex", "OpenAlex", search_openalex, True, "api.openalex.org"),
        Source("crossref", "Crossref", search_crossref, True, "api.crossref.org"),
        Source("semantic_scholar", "Semantic Scholar", search_semantic_scholar,
               False, "api.semanticscholar.org"),
        Source("arxiv", "arXiv", search_arxiv, True, "export.arxiv.org"),
        Source("europe_pmc", "Europe PMC", search_europe_pmc, True,
               "www.ebi.ac.uk"),
        Source("chemrxiv", "ChemRxiv (via Crossref)", search_chemrxiv, False,
               "api.crossref.org"),
    )
}

#: Source priority when merging: the first record of a work fixes the
#: fields the later ones only fill in.
_PRIORITY = ("openalex", "crossref", "europe_pmc", "semantic_scholar",
             "chemrxiv", "arxiv")


def default_sources() -> list[str]:
    """The sensible default set: keyless sources that rarely refuse, plus
    Semantic Scholar when a key is configured."""
    out = [k for k, s in SOURCES.items() if s.default]
    if _credential("SEMANTIC_SCHOLAR_API_KEY"):
        out.append("semantic_scholar")
    return out


DEFAULT_SOURCES = tuple(k for k, s in SOURCES.items() if s.default)


# --- merge / de-duplicate -----------------------------------------------------------

def _norm_title(title: str) -> str:
    return re.sub(r"[^0-9a-z]+", "", str(title or "").casefold())


def _merge_into(target: dict, rec: dict) -> None:
    for src in rec.get("sources") or []:
        if src not in target["sources"]:
            target["sources"].append(src)
    for field in ("title", "venue", "doi", "url", "type"):
        if not target.get(field) and rec.get(field):
            target[field] = rec[field]
    if target.get("year") is None and rec.get("year") is not None:
        target["year"] = rec["year"]
    if len(rec.get("authors") or []) > len(target.get("authors") or []):
        target["authors"] = list(rec["authors"])
    if rec.get("cited_by") is not None:
        target["cited_by"] = max(target.get("cited_by") or 0, rec["cited_by"])
    target["is_oa"] = bool(target.get("is_oa") or rec.get("is_oa"))
    # A published version wins over a preprint label.
    if target.get("preprint") and not rec.get("preprint") and rec.get("venue"):
        target["preprint"] = ""
        target["venue"] = rec["venue"]
    have_storable = licence_permits_storage(target.get("pdf_license"))
    if rec.get("pdf_url") and (not target.get("pdf_url") or (
            not have_storable and licence_permits_storage(rec.get("pdf_license")))):
        target["pdf_url"] = rec["pdf_url"]
        target["pdf_license"] = rec.get("pdf_license") or ""
    for k, v in (rec.get("ids") or {}).items():
        target["ids"].setdefault(k, v)


def merge_records(per_source: "list[tuple[str, list[dict]]]") -> list[dict]:
    """Merge the sources' result lists into one de-duplicated list.

    Same work = same DOI; without a DOI on either side, the same
    normalised title and a year at most one apart (a preprint and its
    paper). Ranking is reciprocal-rank fusion over the sources' own
    relevance orders (``score``); ``sources`` lists who found it.
    """
    order = {k: i for i, k in enumerate(_PRIORITY)}
    per_source = sorted(per_source, key=lambda kv: order.get(kv[0], 99))
    merged: list[dict] = []
    by_doi: dict[str, dict] = {}
    by_title: dict[str, list[dict]] = {}
    for _name, records in per_source:
        for rank, rec in enumerate(records):
            rec = copy.deepcopy(rec)
            gain = 1.0 / (60 + rank)
            doi = rec.get("doi") or ""
            tkey = _norm_title(rec.get("title"))
            target = by_doi.get(doi) if doi else None
            if target is None and tkey:
                for cand in by_title.get(tkey, []):
                    same_doi_or_none = (not doi or not cand.get("doi")
                                        or cand.get("doi") == doi)
                    years = (cand.get("year"), rec.get("year"))
                    close = (None in years or abs(years[0] - years[1]) <= 1)
                    if same_doi_or_none and close:
                        target = cand
                        break
            if target is None:
                rec["score"] = gain
                rec["sources"] = list(rec.get("sources") or [])
                rec["ids"] = dict(rec.get("ids") or {})
                merged.append(rec)
                if doi:
                    by_doi[doi] = rec
                if tkey:
                    by_title.setdefault(tkey, []).append(rec)
                continue
            target["score"] = target.get("score", 0.0) + gain
            had_doi = target.get("doi")
            _merge_into(target, rec)
            if not had_doi and target.get("doi"):
                by_doi[target["doi"]] = target
    return sort_records(merged, "relevance")


def sort_records(records: list[dict], how: str = "relevance") -> list[dict]:
    if how == "year":
        return sorted(records, key=lambda r: (r.get("year") is None,
                                              -(r.get("year") or 0),
                                              -r.get("score", 0.0)))
    if how == "citations":
        return sorted(records, key=lambda r: (r.get("cited_by") is None,
                                              -(r.get("cited_by") or 0),
                                              -r.get("score", 0.0)))
    return sorted(records, key=lambda r: -r.get("score", 0.0))


# --- federation ------------------------------------------------------------------------

_EXECUTOR = _cf.ThreadPoolExecutor(max_workers=8,
                                   thread_name_prefix="delfin-litsearch")
_COOLDOWN: dict[str, float] = {}
_COOLDOWN_LOCK = threading.Lock()


def _one_line(exc: BaseException) -> str:
    if isinstance(exc, urllib.error.HTTPError):
        if exc.code == 429:
            return "rate limited (HTTP 429), try again later"
        return f"HTTP {exc.code}"
    if isinstance(exc, urllib.error.URLError):
        return f"network error: {_clean(exc.reason)[:80]}"
    if isinstance(exc, (TimeoutError, _cf.TimeoutError)):
        return "timed out"
    text = _clean(str(exc)) or type(exc).__name__
    return text[:120]


def _run_source(key: str, query: str, max_results: int,
                timeout_s: float) -> list[dict]:
    ckey = f"{_CACHE_PREFIX}{key}:{max_results}:{query.casefold()}"
    hit = source_apis._cached(ckey)
    if hit is not None:
        return copy.deepcopy(hit["records"])
    with _COOLDOWN_LOCK:
        until = _COOLDOWN.get(key, 0.0)
    if until > time.monotonic():
        raise SourceError("rate limited (HTTP 429), paused for a minute")
    try:
        records = SOURCES[key].search(query, max_results, timeout_s)
    except urllib.error.HTTPError as exc:
        if exc.code == 429:
            with _COOLDOWN_LOCK:
                _COOLDOWN[key] = time.monotonic() + _COOLDOWN_S
        raise SourceError(_one_line(exc)) from exc
    except SourceError:
        raise
    except Exception as exc:
        raise SourceError(_one_line(exc)) from exc
    source_apis._store(ckey, {"records": copy.deepcopy(records)})
    return records


def federated_search(query: str, *, sources=None, max_results: int = 10,
                     timeout_s: float = _DEFAULT_TIMEOUT_S, sort: str = "relevance",
                     unpaywall: "bool | None" = None) -> dict:
    """Search the chosen sources concurrently; merged, de-duplicated results.

    Every source runs in its own worker with its own timeout; whatever
    has answered after ``timeout_s`` is merged and shown, a source still
    running is reported as timed out (it keeps running and fills the cache
    for the next search). Returns ``{"query", "results", "status",
    "notes"}``; ``status`` maps each source to ``"ok"``, ``"error: ..."``
    or ``"timed out"``. Unpaywall (when on and a contact e-mail is
    configured) adds open-access PDF links for DOIs that have none.
    """
    text = _clean(query)
    if not text:
        return {"query": "", "results": [], "status": {},
                "notes": ["Type a query first."]}
    why = web_tools.outbound_secret_reason(text)
    if why:
        return {"query": text, "results": [], "status": {},
                "notes": [f"Not sent: {why}."]}
    keys = [k for k in (sources if sources is not None else default_sources())
            if k in SOURCES and k != "unpaywall"]
    want_unpaywall = (unpaywall if unpaywall is not None
                      else bool(sources is None or "unpaywall" in (sources or ())))
    deadline = time.monotonic() + float(timeout_s)
    futures = {_EXECUTOR.submit(_run_source, k, text, int(max_results),
                                float(timeout_s)): k for k in keys}
    done, pending = _cf.wait(futures, timeout=float(timeout_s))
    status: dict[str, str] = {}
    per_source = []
    for fut in done:
        key = futures[fut]
        try:
            per_source.append((key, fut.result()))
            status[key] = "ok"
        except Exception as exc:
            status[key] = "error: " + _one_line(exc)
    for fut in pending:
        status[futures[fut]] = "timed out"
    results = merge_records(per_source)
    notes: list[str] = []
    if want_unpaywall:
        if not contact_email():
            notes.append(_UNPAYWALL_SKIP)
        else:
            status["unpaywall"] = _enrich_with_unpaywall(results, deadline)
    return {"query": text, "results": sort_records(results, sort),
            "status": status, "notes": notes}


def _enrich_with_unpaywall(results: list[dict], deadline: float,
                           limit: int = 10) -> str:
    todo = [r for r in results
            if r.get("doi") and not licence_permits_storage(r.get("pdf_license"))][:limit]
    if not todo:
        return "ok"
    remaining = max(2.0, deadline - time.monotonic())
    futs = {_EXECUTOR.submit(unpaywall_lookup, r["doi"], timeout_s=remaining): r
            for r in todo}
    done, pending = _cf.wait(futs, timeout=remaining)
    errors = 0
    for fut in done:
        rec = futs[fut]
        try:
            got = fut.result()
        except Exception:
            errors += 1
            continue
        if got.get("status") == "error":
            errors += 1
        if got.get("status") != "ok" or not got.get("pdf_url"):
            continue
        if not rec.get("pdf_url") or licence_permits_storage(got.get("license")):
            rec["pdf_url"] = got["pdf_url"]
            rec["pdf_license"] = got.get("license") or ""
            rec["is_oa"] = True
            if "unpaywall" not in rec["sources"]:
                rec["sources"].append("unpaywall")
    if pending:
        return "timed out"
    return f"error: {errors} lookup(s) failed" if errors else "ok"


def status_line(status: dict, notes: "list[str] | None" = None) -> str:
    """One line naming the sources that failed (empty when all answered)."""
    bad = [f"{SOURCES[k].label if k in SOURCES else k.title()}: "
           f"{v[len('error: '):] if v.startswith('error: ') else v}"
           for k, v in status.items() if v != "ok"]
    return " · ".join(bad + list(notes or []))


# --- method references (agent side) -----------------------------------------------

#: Search phrases for the titles of the ORIGINAL papers of common ORCA
#: keywords. They are search terms only: the candidates and their DOIs
#: always come from OpenAlex, and a DOI is cited only after verify_doi.
_METHOD_TITLE_HINTS = {
    "b3lyp": ["Density-functional thermochemistry. III. The role of exact exchange",
              "Development of the Colle-Salvetti correlation-energy formula "
              "into a functional of the electron density"],
    "pbe": ["Generalized gradient approximation made simple"],
    "pbe0": ["Toward reliable density functional methods without adjustable "
             "parameters: The PBE0 model"],
    "bp86": ["Density-functional exchange-energy approximation with correct "
             "asymptotic behavior",
             "Density-functional approximation for the correlation energy of "
             "the inhomogeneous electron gas"],
    "tpss": ["Climbing the density functional ladder: nonempirical "
             "meta-generalized gradient approximation"],
    "tpssh": ["Comparative assessment of a new nonempirical density "
              "functional: Molecules and hydrogen-bonded complexes"],
    "r2scan": ["Accurate and numerically efficient r2SCAN meta-generalized "
               "gradient approximation"],
    "r2scan-3c": ["r2SCAN-3c: A Swiss army knife composite electronic-structure method"],
    "b97-3c": ["B97-3c: A revised low-cost variant of the B97-D density "
               "functional method"],
    "wb97x-d": ["Long-range corrected hybrid density functionals with damped "
                "atom-atom dispersion corrections"],
    "cam-b3lyp": ["A new hybrid exchange-correlation functional using the "
                  "Coulomb-attenuating method"],
    "m06": ["The M06 suite of density functionals for main group "
            "thermochemistry"],
    "def2": ["Balanced basis sets of split valence, triple zeta valence and "
             "quadruple zeta valence quality for H to Rn"],
    "d3": ["A consistent and accurate ab initio parametrization of density "
           "functional dispersion correction (DFT-D) for the 94 elements H-Pu"],
    "d3bj": ["Effect of the damping function in dispersion corrected density "
             "functional theory"],
    "d4": ["A generally applicable atomic-charge dependent London dispersion "
           "correction"],
    "smd": ["Universal solvation model based on solute electron density and "
            "on a continuum model of the solvent"],
    "cpcm": ["Quantum calculation of molecular energies and energy gradients "
             "in solution by a conductor solvent model"],
    "rijcosx": ["Efficient, approximate and parallel Hartree-Fock and hybrid "
                "DFT calculations. A chain-of-spheres algorithm"],
    "dlpno-ccsd(t)": ["Natural triple excitations in local coupled cluster "
                      "calculations with pair natural orbitals"],
    "zora": ["Relativistic regular two-component Hamiltonians"],
    "cc-pvtz": ["Gaussian basis sets for use in correlated molecular "
                "calculations. I. The atoms boron through neon and hydrogen"],
    "gfn2-xtb": ["GFN2-xTB An accurate and broadly parametrized "
                 "self-consistent tight-binding quantum chemical method"],
    "orca": ["The ORCA quantum chemistry program package",
             "Software update: The ORCA program system Version 5.0"],
}

_METHOD_ALIASES = {
    "def2-svp": "def2", "def2-tzvp": "def2", "def2-tzvpp": "def2",
    "def2-qzvp": "def2", "def2-qzvpp": "def2", "def2-svpd": "def2",
    "def2-tzvpd": "def2", "d3zero": "d3", "d3(bj)": "d3bj",
    "wb97x-d3": "wb97x-d", "ωb97x-d": "wb97x-d", "cc-pvdz": "cc-pvtz",
    "cc-pvqz": "cc-pvtz", "rij": "rijcosx", "cosx": "rijcosx",
    "cpcm(water)": "cpcm", "pbeh": "pbe0", "pbe1pbe": "pbe0",
    "x2c": "zora",
}

#: OpenAlex topic fields (Scopus ASJC areas) where method papers sit:
#: Chemistry (16), Physics and Astronomy (31), Materials Science (25),
#: Computer Science (17). Used as the topic/concept filter; OpenAlex's
#: legacy concepts are superseded by topics.
_METHOD_FIELD_FILTER = "primary_topic.field.id:fields/16|fields/31|fields/25|fields/17"


def _method_key(name: str) -> str:
    key = re.sub(r"\s+", "", str(name or "").strip().lower())
    key = key.replace("ω", "w")
    return _METHOD_ALIASES.get(key, key)


def find_method_reference(name: str, *, max_results: int = 5,
                          timeout_s: float = _DEFAULT_TIMEOUT_S) -> dict:
    """Original-paper candidates for an ORCA method, basis set or functional.

    Uses OpenAlex: a title search (the method name, or the known title of
    its original paper) restricted to the chemistry/physics topic fields,
    most-cited first -- the original paper of a method is nearly always
    its most-cited title match. Each candidate carries its DOI; cite one
    only after ``verify_doi`` confirmed it.

    Returns ``{"status": "ok"|"not_found"|"error"|"invalid_input",
    "method", "candidates": [record + "match"], "note"}``.
    """
    raw = _clean(name)
    if not raw or len(raw) > 80:
        return {"status": "invalid_input", "method": raw, "candidates": [],
                "note": "give one method, basis-set or functional name"}
    why = web_tools.outbound_secret_reason(raw)
    if why:
        return {"status": "invalid_input", "method": raw, "candidates": [],
                "note": f"not sent: {why}"}
    key = _method_key(raw)
    hints = _METHOD_TITLE_HINTS.get(key, [])
    searches = [(h, "known-title") for h in hints] or [(raw, "name-in-title")]
    ckey = f"{_CACHE_PREFIX}method:{key}:{max_results}"
    hit = source_apis._cached(ckey)
    if hit is not None:
        return copy.deepcopy(hit)
    found: list[dict] = []
    errors: list[str] = []
    per_search = max(1, min(int(max_results), 10))
    for phrase, how in searches:
        title = re.sub(r"[,:|]", " ", phrase)
        params = {"filter": f"title.search:{title},{_METHOD_FIELD_FILTER}",
                  "sort": "cited_by_count:desc", "per-page": per_search}
        try:
            recs = _openalex_request(params, timeout_s)
        except urllib.error.HTTPError as exc:
            if exc.code == 429:
                errors.append(_one_line(exc))
                continue
            recs = []            # e.g. 400 for a filter OpenAlex rejects
        except Exception as exc:
            errors.append(_one_line(exc))
            continue
        if not recs:
            # The topic filter can be too narrow (a method paper filed
            # under another field) or rejected: one retry on the title alone.
            params["filter"] = f"title.search:{title}"
            try:
                recs = _openalex_request(params, timeout_s)
            except Exception as exc:
                errors.append(_one_line(exc))
                continue
        for rec in recs:
            rec["match"] = how
            found.append(rec)
    merged = merge_records([("openalex", found)]) if found else []
    for rec in merged:
        rec.pop("score", None)
    merged = sort_records(merged, "citations")[: int(max_results)]
    if not merged and errors:
        return {"status": "error", "method": raw, "candidates": [],
                "note": "OpenAlex: " + errors[0]}
    out = {"status": "ok" if merged else "not_found", "method": raw,
           "candidates": merged,
           "note": ("Candidates from OpenAlex, most cited first. Verify the "
                    "DOI with verify_doi before citing; the original paper "
                    "of a composite method may be several papers.")}
    source_apis._store(ckey, out)
    return copy.deepcopy(out)


# --- storing a record -------------------------------------------------------------------

def suggest_basename(rec: dict) -> str:
    """``Author2021_First_title_words``: file-name safe, at most 70 chars."""
    first = (rec.get("authors") or ["Unknown"])[0]
    if "," in first:                       # "Doe, Jane"
        tokens = first.split(",")[0].split()
    else:
        tokens = first.split()
        # Europe PMC / PubMed style "Doe J" or "Doe JA": surname first.
        if len(tokens) > 1 and re.fullmatch(r"[A-Z]{1,3}\.?", tokens[-1]):
            tokens = tokens[:1]
    last = tokens[-1] if tokens else "Unknown"
    year = str(rec.get("year") or "nd")
    words = re.findall(r"[A-Za-z0-9]+", rec.get("title") or "")[:6]
    base = f"{last}{year}_" + "_".join(words)
    base = re.sub(r"[^A-Za-z0-9_\-]", "", base.encode("ascii", "ignore").decode())
    return (base or "reference")[:70].rstrip("_-")


def _bib_escape(text: str) -> str:
    return str(text or "").replace("\\", "").replace("{", "").replace("}", "")


def to_bibtex(rec: dict) -> str:
    key = re.sub(r"[^A-Za-z0-9]", "", suggest_basename(rec).split("_")[0]) or "ref"
    kind = "misc" if rec.get("preprint") else "article"
    fields = [("title", rec.get("title")),
              ("author", " and ".join(rec.get("authors") or [])),
              ("journal" if kind == "article" else "howpublished", rec.get("venue")),
              ("year", rec.get("year")), ("doi", rec.get("doi")),
              ("url", rec.get("url"))]
    body = ",\n".join(f"  {k} = {{{_bib_escape(v)}}}" for k, v in fields if v)
    return f"@{kind}{{{key},\n{body}\n}}\n"


def to_markdown(rec: dict, *, retrieved: "str | None" = None,
                pdf_file: str = "", pdf_note: str = "") -> str:
    retrieved = retrieved or _dt.date.today().isoformat()
    lines = [f"# {rec.get('title') or 'Untitled'}", ""]
    if rec.get("authors"):
        lines.append(f"- Authors: {', '.join(rec['authors'])}")
    if rec.get("venue") or rec.get("year"):
        lines.append(f"- Published: {rec.get('venue') or '?'}, {rec.get('year') or 'n.d.'}")
    if rec.get("preprint"):
        lines.append(f"- Preprint: {rec['preprint']}")
    if rec.get("doi"):
        lines.append(f"- DOI: [{rec['doi']}](https://doi.org/{rec['doi']})")
    elif rec.get("url"):
        lines.append(f"- Link: {rec['url']}")
    if rec.get("cited_by") is not None:
        lines.append(f"- Cited by: {rec['cited_by']} (at retrieval)")
    if rec.get("pdf_url"):
        lic = rec.get("pdf_license") or "licence unknown"
        lines.append(f"- Open-access PDF: {rec['pdf_url']} ({lic})")
    if pdf_file:
        lines.append(f"- Local copy: {pdf_file}")
    if pdf_note:
        lines.append(f"- PDF: {pdf_note}")
    names = [SOURCES[s].label if s in SOURCES else s.title()
             for s in rec.get("sources") or []]
    lines.append(f"- Found in: {', '.join(names)}")
    lines.append(f"- Retrieved: {retrieved}")
    lines += ["", "```bibtex", to_bibtex(rec).rstrip(), "```", ""]
    return "\n".join(lines)


def download_open_pdf(rec: dict, *, timeout_s: float = 30.0,
                      max_bytes: int = _MAX_PDF_BYTES) -> bytes:
    """The open-access PDF of ``rec``, when its licence permits keeping it.

    Refuses (ValueError, one line) a missing or non-permissive licence, a
    non-https link, an arXiv e-print (see module doc), anything that is
    not a PDF (a publisher landing or challenge page is never stored),
    and anything above ``max_bytes``. Same guarded opener as every other
    request: deny-list, SSRF checks, redirect re-validation.
    """
    url = _http_url(rec.get("pdf_url"))
    if not url:
        raise ValueError("no open-access PDF link")
    if not licence_permits_storage(rec.get("pdf_license")):
        raise ValueError("the licence ({}) does not allow keeping a copy; "
                         "link only".format(rec.get("pdf_license") or "unknown"))
    host = (urllib.parse.urlsplit(url).hostname or "").lower()
    if host.endswith("arxiv.org"):
        raise ValueError("arXiv e-prints are linked, not stored")
    why = web_tools._check_url(url)
    if why:
        raise ValueError(f"refused: {why}")
    req = urllib.request.Request(url, headers={
        "User-Agent": web_tools._USER_AGENT, "Accept": "application/pdf"})
    with source_apis._HostSlot(host, timeout_s):
        source_apis._rate_limit(host)
        with web_tools._guarded_opener().open(req, timeout=timeout_s) as resp:
            ctype = str(resp.headers.get("Content-Type", "")).lower()
            body = resp.read(max_bytes + 1)
    if len(body) > max_bytes:
        raise ValueError("PDF larger than the {} MB limit".format(max_bytes // 2**20))
    if not body.startswith(b"%PDF-"):
        raise ValueError("the link answered with {} rather than a PDF; open it "
                         "in the browser instead".format(ctype or "something else"))
    return body
