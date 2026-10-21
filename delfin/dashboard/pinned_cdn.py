"""The files the dashboard page loads from another origin: one place.

Input: none. Output: :data:`ASSETS`, each a URL pinned to an exact version
and the SRI hash of that file's bytes. The page loads each with
``integrity`` + ``crossorigin``, so the browser refuses a file whose bytes
differ, and ``web_security.PAGE_CSP`` names exactly these URLs (full
paths, not the host), so the page can load no other file from the same
CDN -- another version or another package -- even if markup that escaped
the output encoding asked for one.

Everything else the page runs comes from the dashboard's own origin: the
3Dmol bundle (``static/3Dmol-min.js``, inlined at start), Voila's and
Jupyter's bundles, and the drawing editor (``/static/ketcher``).

A new foreign script or stylesheet is added here, with its hash, or it is
refused by the browser; ``tests/test_the_dashboard_page_loads_code_only_
from_where_it_names.py`` fails for a foreign URL in the dashboard source
that is not pinned and named.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Asset:
    url: str
    sri: str
    kind: str  # "script" or "style"


XTERM_VERSION = "5.5.0"
XTERM_FIT_VERSION = "0.10.0"
MERMAID_VERSION = "11.17.2"

XTERM_JS = Asset(
    "https://cdn.jsdelivr.net/npm/@xterm/xterm@5.5.0/lib/xterm.js",
    "sha384-M169f14mRZOXm3hD/v2Ti0ThIT/RnAQagXA9nlE15yHAtrW19gdePJh/HaTzUOe/",
    "script")
XTERM_CSS = Asset(
    "https://cdn.jsdelivr.net/npm/@xterm/xterm@5.5.0/css/xterm.css",
    "sha384-8Xk9wy/gzEDUKrXtrmCFa2bBuK3BpjpDuL/p0SeKQX19Khl/M+lHOgD/CyYf7efP",
    "style")
XTERM_FIT_JS = Asset(
    "https://cdn.jsdelivr.net/npm/@xterm/addon-fit@0.10.0/lib/addon-fit.js",
    "sha384-iF+jqbuti4XlB64clWgFWYEscb+UnSRv3VgVikGYZu+otNFnSHr7y7NcKfBnGizn",
    "script")
MERMAID_JS = Asset(
    "https://cdn.jsdelivr.net/npm/mermaid@11.17.2/dist/mermaid.min.js",
    "sha384-EOXBFmc3gx5mb+vn0vPvvGqACToJD24hhacX5Yx+8NUUQrHIle/Qi5Bg9o3zKwW2",
    "script")

# The URLs are spelled out (not concatenated) so a search for a CDN finds
# every one; the versions above are the ones in them.
ASSETS = (XTERM_JS, XTERM_CSS, XTERM_FIT_JS, MERMAID_JS)


def urls(kind: str) -> tuple[str, ...]:
    """The pinned URLs of one kind, in a fixed order."""
    return tuple(a.url for a in ASSETS if a.kind == kind)
