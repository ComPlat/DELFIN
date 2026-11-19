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
  work search. Data CC0 (https://help.openalex.org/how-to-use-the-api/
  api-overview). Terms re-checked 2026-11-11
  (https://help.openalex.org/how-to-use-the-api/rate-limits-and-authentication):
  the API key is optional (``api_key`` parameter or an ``Authorization:
  Bearer`` header; DELFIN sends it as the header, never in a URL);
  keyless use is 10 requests/s with a small daily USD budget (measured:
  ``x-ratelimit-limit-usd: 0.1`` per day, a search costs $0.001), a free
  key gives 30/s and 10x the budget. The docs no longer mention a polite
  pool; a ``mailto`` contact is sent only when the person configured one.

The federated literature search (OpenAlex, Crossref, Semantic Scholar,
arXiv, Europe PMC, Unpaywall) lives in ``delfin.agent.scholarly_sources``
and is re-exported here (``federated_search``, ``find_method_reference``)
so the agent and the dashboard use one implementation.

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
           "host_for_action", "fetch_bytes", "federated_search",
           "find_method_reference"]


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

#: Per-host spacing where a service states a stricter limit than the
#: default (each value from that service's own terms, see
#: delfin.agent.scholarly_sources): arXiv one request every three
#: seconds; Semantic Scholar 1 request/s with a key, a shared pool
#: without; Crossref's public pool 5/s with concurrency 1.
_HOST_MIN_DELAY_S = {
    "export.arxiv.org": 3.0,
    "api.semanticscholar.org": 1.0,
    "api.crossref.org": 0.25,
    "api.openalex.org": 0.15,
    "www.ebi.ac.uk": 0.15,
    "api.unpaywall.org": 0.1,
}

#: Hosts whose terms allow only this many connections at once.
_HOST_CONCURRENCY = {"export.arxiv.org": 1, "api.crossref.org": 1}
_HOST_SLOTS: dict[str, threading.BoundedSemaphore] = {}


def _rate_limit(host: str) -> None:
    """Wait for this host's next request slot.

    The slot is reserved under the lock and the wait happens outside it,
    so a slow host (arXiv: 3 s) never holds up requests to another host
    that run concurrently in the federated search.
    """
    delay = _HOST_MIN_DELAY_S.get(host, _MIN_DELAY_S)
    with _FETCH_LOCK:
        now = time.monotonic()
        last = _LAST_REQUEST.get(host)
        slot = now if last is None else max(now, last + delay)
        _LAST_REQUEST[host] = slot
    if slot > now:
        time.sleep(slot - now)


class _HostSlot:
    """At most ``_HOST_CONCURRENCY[host]`` requests to ``host`` at once."""

    def __init__(self, host: str, timeout_s: float):
        self.sem = None
        limit = _HOST_CONCURRENCY.get(host)
        if limit:
            with _FETCH_LOCK:
                self.sem = _HOST_SLOTS.setdefault(
                    host, threading.BoundedSemaphore(limit))
        self.host = host
        self.timeout_s = timeout_s

    def __enter__(self):
        if self.sem is not None and not self.sem.acquire(
                timeout=max(0.1, float(self.timeout_s))):
            raise TimeoutError(f"{self.host} is busy with another request")
        return self

    def __exit__(self, *exc):
        if self.sem is not None:
            self.sem.release()
        return False


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
    try:
        with _HostSlot(host, timeout_s):
            _rate_limit(host)
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


def fetch_bytes(url: str, *, timeout_s: float = 15,
                headers: "dict | None" = None) -> bytes:
    """GET ``url`` through web_tools' guarded stack; the raw body.

    The same deny-list, SSRF checks, per-host rate limit and connection
    cap as :func:`fetch_json`, for the services that answer in another
    format (arXiv's Atom) or need a request header (an optional API key
    travels in a header, never in a URL). An HTTP error is re-raised as
    ``urllib.error.HTTPError`` so a caller can tell 404 from 429; a body
    larger than web_tools' cap is an error, never a silently cut record.
    """
    parsed = urllib.parse.urlsplit(url)
    host = parsed.hostname or ""
    why = web_tools._check_url(url)
    if why:
        raise ValueError(f"refused to ask {host!r}: {why}")
    with _HostSlot(host, timeout_s):
        _rate_limit(host)
        got = web_tools._fetch_bytes(url, timeout_s, headers=headers)
    body = got[0]
    if len(body) > web_tools._MAX_BYTES:
        raise ValueError(f"{host} answered with more than "
                         f"{web_tools._MAX_BYTES} bytes")
    return body


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


# --- federated literature search (delfin.agent.scholarly_sources) -------------

def federated_search(query: str, **kwargs) -> dict:
    """Search several scholarly APIs at once; see scholarly_sources."""
    from . import scholarly_sources
    return scholarly_sources.federated_search(query, **kwargs)


def find_method_reference(name: str, **kwargs) -> dict:
    """Original-paper candidates for a method/basis/functional name."""
    from . import scholarly_sources
    return scholarly_sources.find_method_reference(name, **kwargs)
