"""Verify literature sources and compound data through OFFICIAL APIs.

The agents verify DOIs and literature values all the time; publisher
pages block scrapers (403, JavaScript walls) and stealth browsing breaks
the sites' terms. This module therefore asks the official, free,
documented APIs instead, and every HTTP call goes through
``delfin.agent.web_tools``' guarded stack — the same single HTTP stack
with its deny-list, redirect re-validation and SSRF IP checks. There is
no second opener here.

Services used (terms checked against the service's own documentation,
not assumed):

- **Crossref** (``api.crossref.org/works/<doi>``): DOI -> title, authors,
  journal, year. Crossref metadata is distributed under CC0 1.0
  (https://www.crossref.org/documentation/metadata-stewardship/
  licensing-and-enforcement/). No key. DELFIN joins the "polite pool"
  with a contact address in the User-Agent (Crossref etiquette:
  https://github.com/CrossRef/rest-api-doc#etiquette).
- **PubChem PUG-REST** (``pubchem.ncbi.nlm.nih.gov/rest/pug``): compound
  lookup by name, SMILES or InChIKey -> CID, formula, InChIKey, molar
  mass. PubChem data is public domain (U.S. National Library of
  Medicine, https://pubchem.ncbi.nlm.nih.gov/docs/). Rate limit: at
  most 5 requests per second, at most 400 requests per minute
  (https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest#section=Access-Statistics).
- **OpenAlex** (already used by ``web_tools.web_search``): scholarly
  work search. No key; a mailto address joins the polite pool.

Behaviour rules enforced below, not just documented:

- Requests are rate-limited per host and cached (in-process TTL cache)
  so repeated verifications in one session do not re-download.
- A "not found" answer is cached like a hit; the cache never invents
  metadata — a record is only ever what an API returned.
- Identifier text never reaches a URL unescaped: path segments are
  percent-quoted, and the host is always the fixed API host — user
  input cannot steer the connection anywhere else.
"""

from __future__ import annotations

import json
import re
import threading
import time
import urllib.error
import urllib.parse

from . import web_tools

__all__ = ["verify_doi", "lookup_compound", "search_works", "clear_cache",
           "host_for_action"]


# --- the fixed host per action -----------------------------------------------

# ONE derivation place for the action -> API host mapping. The
# permission rule layer reads it (permission_rules._subjects) so the
# managed network policy and web_fetch(domain:) rules bind verify_source
# like they bind web_fetch (operator decision 2026-11-01, closing s4's
# P2 finding). An unknown action answers None — never a guessable
# fallback host.
_ACTION_HOSTS = {
    "doi": "api.crossref.org",
    "compound": "pubchem.ncbi.nlm.nih.gov",
    "works": "api.openalex.org",
}


def host_for_action(action) -> "str | None":
    """The fixed API host an action dials, or None for an unknown action."""
    return _ACTION_HOSTS.get(str(action or "").strip().lower())


# --- guarded transport -------------------------------------------------------

_FETCH_LOCK = threading.Lock()
_LAST_REQUEST: dict[str, float] = {}
_MIN_DELAY_S = 0.3  # per host: ~3 req/s, inside both APIs' stated limits


def _rate_limit(host: str) -> None:
    with _FETCH_LOCK:
        now = time.monotonic()
        last = _LAST_REQUEST.get(host)
        if last is not None and now - last < _MIN_DELAY_S:
            time.sleep(_MIN_DELAY_S - (now - last))
        _LAST_REQUEST[host] = time.monotonic()


def _fallback_fetch(url: str, *, timeout_s: int = 30) -> bytes:
    """Body bytes from web_tools' guarded stack (the pre-hook path).

    Used only while the operator has not yet built the public
    ``web_tools.fetch_json`` hook; same deny-list, same opener.
    """
    why = web_tools._check_url(url)
    if why:
        raise ValueError(f"refused to send this request: {why}")
    body, _ctype = web_tools._fetch_bytes(url, timeout_s)
    return body


#: The guarded JSON fetch source_apis uses for every call. Once the
#: web_tools hook exists it is that function (one stack, one place); a
#: monkeypatched or absent attribute falls back to the direct guarded
#: bytes path above. Never a urllib call of its own.
_GUARDED_FETCH = getattr(web_tools, "fetch_json", None) or _fallback_fetch


def fetch_json(url: str, *, timeout_s: int = 30) -> dict:
    """GET ``url`` through the guarded stack and parse the JSON.

    Raises ValueError with a clear message on a non-JSON body (a
    challenge page must never come back as an empty record).
    """
    parsed = urllib.parse.urlsplit(url)
    host = parsed.hostname or ""
    why = web_tools._check_url(url)
    if why:
        raise ValueError(f"refused to ask {host!r}: {why}")
    _rate_limit(host)
    try:
        data = _GUARDED_FETCH(url, timeout_s=timeout_s)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise
        raise ValueError(f"{host} answered HTTP {exc.code}") from exc
    if isinstance(data, dict):
        # the web_tools hook returns a parsed record already
        return data
    if not isinstance(data, (bytes, str)):
        raise ValueError(f"{host} answered with a non-JSON body")
    try:
        return json.loads(data if isinstance(data, str)
                          else data.decode("utf-8", errors="replace"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ValueError(f"{host} answered with non-JSON body: {exc}") from exc


# --- cache -------------------------------------------------------------------

_CACHE: dict[str, tuple[float, dict]] = {}
_CACHE_LOCK = threading.Lock()
_TTL_S = 3600.0


def clear_cache() -> None:
    with _CACHE_LOCK:
        _CACHE.clear()


def _cached(key: str):
    with _CACHE_LOCK:
        entry = _CACHE.get(key)
        if entry is None:
            return None
        stamp, value = entry
    if time.monotonic() - stamp > _TTL_S:
        return None
    return value


def _store(key: str, value: dict) -> None:
    with _CACHE_LOCK:
        _CACHE[key] = (time.monotonic(), value)


# --- Crossref: DOI verification ----------------------------------------------

_CROSSREF_HOST = "api.crossref.org"
# A DOI is a case-insensitive prefix 10.<registrar>/<suffix>; the suffix
# may contain most printable characters, so it is quoted, never filtered.
_DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")


def _normalise_doi(doi: str) -> str:
    text = str(doi or "").strip()
    # Strip the common "doi:" / "DOI:" / "https://doi.org/" wrappers so
    # the underlying identifier is what is compared.
    text = re.sub(r"^(?:(?:https?://)?doi\.org/|doi:)", "", text,
                  flags=re.IGNORECASE)
    return text.strip().lower()


def _crossref_record(payload: dict) -> dict:
    msg = payload.get("message") or {}
    title = " ".join(str(t) for t in (msg.get("title") or []))
    authors = [
        " ".join(x for x in (a.get("given"), a.get("family")) if x)
        for a in (msg.get("author") or [])
    ]
    issued = ((msg.get("issued") or {}).get("date-parts") or [[]])[0]
    return {
        "status": "ok",
        "doi": str(msg.get("DOI") or "").lower(),
        "title": title,
        "authors": [a for a in authors if a],
        "journal": " ".join(str(t) for t in (msg.get("container-title") or [])),
        "year": issued[0] if issued else None,
        "type": str(msg.get("type") or ""),
        "source": "crossref",
    }


def verify_doi(doi: str) -> dict:
    """Resolve a DOI through Crossref.

    Returns ``{"status": "ok", ...metadata}``, ``{"status": "not_found"}``
    (the DOI does not exist — an invented DOI is NOT echoed as verified
    metadata), or ``{"status": "invalid_input"}`` for anything not shaped
    like a DOI at all.
    """
    normalised = _normalise_doi(doi)
    if not normalised or not _DOI_RE.match(normalised):
        return {"status": "invalid_input",
                "error": "not shaped like a DOI (10.xxxx/suffix)"}
    key = "doi:" + normalised
    hit = _cached(key)
    if hit is not None:
        return dict(hit)
    path = urllib.parse.quote(normalised, safe="")
    url = f"https://{_CROSSREF_HOST}/works/{path}"
    try:
        payload = fetch_json(url)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            rec = {"status": "not_found", "source": "crossref"}
            _store(key, rec)
            return dict(rec)
        rec = {"status": "error", "error": f"HTTP {exc.code}",
               "source": "crossref"}
        return rec
    rec = _crossref_record(payload)
    rec["doi"] = normalised
    _store(key, rec)
    return dict(rec)


# --- PubChem: compound lookup -------------------------------------------------

_PUBCHEM_HOST = "pubchem.ncbi.nlm.nih.gov"
_PUBCHEM_PROPS = "MolecularFormula,InChIKey,MolecularWeight"
# An InChIKey is three hyphen-separated blocks of letters.
_INCHIKEY_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


def _pubchem_record(payload: dict) -> dict:
    # Real PUG-REST /property shape: {"PropertyTable": {"Properties": [
    #   {"CID": ..., "MolecularFormula": ..., ...}]}}
    rows = ((payload.get("PropertyTable") or {}).get("Properties")) or []
    row = rows[0] if rows else {}
    return {
        "status": "ok",
        "cid": row.get("CID"),
        "formula": row.get("MolecularFormula"),
        "inchikey": row.get("InChIKey"),
        "molar_mass_g_mol": row.get("MolecularWeight"),
        "source": "pubchem",
    }


def _pubchem_lookup(route: str, value: str) -> dict:
    """One PUG-REST route per identifier kind; the value is quoted."""
    path = urllib.parse.quote(value, safe="")
    url = (f"https://{_PUBCHEM_HOST}/rest/pug/compound/{route}/{path}"
           f"/property/{_PUBCHEM_PROPS}/JSON")
    key = "pubchem:" + route + ":" + value.lower()
    hit = _cached(key)
    if hit is not None:
        return dict(hit)
    try:
        payload = fetch_json(url)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            rec = {"status": "not_found", "source": "pubchem"}
            _store(key, rec)
            return dict(rec)
        rec = {"status": "error", "error": f"HTTP {exc.code}",
               "source": "pubchem"}
        return rec
    rec = _pubchem_record(payload)
    _store(key, rec)
    return dict(rec)


def lookup_compound(value: str) -> dict:
    """Look a compound up in PubChem by name, SMILES or InChIKey.

    The three identifier kinds are separate PUG-REST routes
    (``/compound/name``, ``/compound/smiles``, ``/compound/inchikey``);
    the route is chosen from the identifier's shape, and a SMILES is
    percent-quoted, never sent raw.
    """
    text = str(value or "").strip()
    if not text:
        return {"status": "invalid_input", "error": "empty identifier"}
    if _INCHIKEY_RE.match(text):
        return _pubchem_lookup("inchikey", text)
    # A SMILES carries bond/stereo/ring marks a bare name never has:
    # slashes, brackets, parentheses, ring closures, chirality.
    if re.search(r"[=\[\]()#/\\@]", text):
        return _pubchem_lookup("smiles", text)
    return _pubchem_lookup("name", text)


# --- OpenAlex: scholarly search ------------------------------------------------

def search_works(query: str, max_results: int = 5) -> list[dict]:
    """Scholarly works for a query, through web_tools' OpenAlex search.

    Same transport, same polite-pool mailto, same deny-list. Returns
    ``[{title, url, snippet, source}]``; an empty list means no (or no
    reachable) results.
    """
    return web_tools.web_search(query, max_results=max_results)
