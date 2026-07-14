"""Lightweight web-search + web-fetch primitives for the KIT agent.

When Jerome's project pulls in an external library (BoTorch, Ax, SMAC,
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

import html as _html
import ipaddress
import json
import re
import socket
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Optional


_USER_AGENT = "Mozilla/5.0 (DELFIN-KIT-agent; +https://kit.edu)"
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


def _check_url(url: str) -> Optional[str]:
    """Return an error string if the URL is forbidden, else None."""
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
                    or ip.is_multicast or ip.is_reserved:
                return (
                    f"resolved IP {ip_str} for host {host!r} is "
                    "private/loopback/link-local — refused"
                )
    except socket.gaierror as exc:
        return f"DNS lookup failed: {exc}"
    return None


def _fetch_bytes(url: str, timeout_s: int) -> tuple[bytes, str]:
    """GET a URL with a strict cap. Returns (body, content_type)."""
    req = urllib.request.Request(
        url, headers={"User-Agent": _USER_AGENT, "Accept": "text/html,*/*"}
    )
    with urllib.request.urlopen(req, timeout=timeout_s) as resp:
        ctype = resp.headers.get("Content-Type", "")
        body = resp.read(_MAX_BYTES + 1)
    return body, ctype


def _html_to_text(body: bytes) -> str:
    """Minimal HTML → plain text. Drops <script>/<style> and tags."""
    try:
        text = body.decode("utf-8", errors="replace")
    except Exception:
        return ""
    # Strip script + style first (their contents are not user-readable).
    text = re.sub(r"<script[^>]*>.*?</script>", " ", text,
                  flags=re.DOTALL | re.IGNORECASE)
    text = re.sub(r"<style[^>]*>.*?</style>", " ", text,
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


def web_search(query: str, *, max_results: int = 8,
               timeout_s: int = _DEFAULT_TIMEOUT_S) -> dict:
    """DuckDuckGo search: HTML result list, with an Instant-Answer fallback
    when the HTML endpoint returns nothing (e.g. blocked by a 202 anomaly)."""
    if not query.strip():
        return {"error": "query must be non-empty"}
    encoded = urllib.parse.urlencode({"q": query})
    url = f"https://duckduckgo.com/html/?{encoded}"
    err = _check_url(url)
    if err:
        return {"error": err}
    try:
        body, ctype = _fetch_bytes(url, timeout_s)
    except urllib.error.HTTPError as exc:
        return {"error": f"HTTP {exc.code}: {exc.reason}"}
    except urllib.error.URLError as exc:
        return {"error": f"network error: {exc.reason}"}
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

    return {"query": query, "results": hits, "result_count": len(hits),
            "source": source}
