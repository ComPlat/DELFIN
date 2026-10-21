"""The dashboard page loads script, style, images and connections only from
sources its Content-Security-Policy names.

Input: ``web_security.PAGE_CSP`` (what every dashboard response carries) and
the dashboard's own source files. Properties asserted:

* the page CSP has fetch directives for scripts, styles, images, fonts,
  workers, frames and connections, and none of them allows a whole host, a
  whole scheme or ``*``: a foreign script is named by its full, versioned
  path, so another file on the same CDN (another version, another package)
  is refused;
* ``connect-src`` is the dashboard's own origin only, so injected markup
  cannot send what it reads to another site with ``fetch``/``XMLHttpRequest``;
* every script or stylesheet URL the dashboard source loads from another
  origin is pinned to a version, carries an SRI hash, and is the one the CSP
  names; there is no unpinned fallback (3Dmol.org without a version used to
  be the fallback of ten hand-written viewers);
* the launcher's MathJax stand-in is served from the dashboard's own origin,
  so ``script-src`` needs no ``data:``.

The browser side of the same claim (every tab rendered under this policy
with no violation) is measured in headless Chromium; see the commit that
introduced the policy.
"""

from __future__ import annotations

import re
from pathlib import Path

from delfin.dashboard import web_security

REPO = Path(__file__).resolve().parents[1]
DASHBOARD = REPO / "delfin" / "dashboard"

FETCH_DIRECTIVES = ("script-src", "style-src", "img-src", "font-src",
                    "connect-src", "worker-src", "frame-src")


def _directives(csp: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for part in csp.split(";"):
        words = part.split()
        if words:
            out.setdefault(words[0].lower(), words[1:])
    return out


def _foreign(sources: list[str]) -> list[str]:
    return [s for s in sources if s.startswith(("http://", "https://"))]


def test_the_page_csp_names_a_source_list_for_every_kind_of_fetch():
    d = _directives(web_security.PAGE_CSP)
    for name in FETCH_DIRECTIVES:
        assert name in d, (name, web_security.PAGE_CSP)
        assert "'self'" in d[name], (name, d[name])


def test_no_directive_allows_a_whole_host_a_whole_scheme_or_everything():
    for name, sources in _directives(web_security.PAGE_CSP).items():
        for source in sources:
            assert source != "*", (name, source)
            assert source not in ("http:", "https:", "ws:", "wss:"), (name, source)
            if source.startswith(("http://", "https://")):
                path = source.split("://", 1)[1].partition("/")[2]
                # a full path to one file, never a host or a directory
                assert path and not path.endswith("/"), (name, source)
                assert "*" not in source, (name, source)


def test_connections_stay_on_the_dashboards_own_origin():
    assert _directives(web_security.PAGE_CSP)["connect-src"] == ["'self'"]


def test_scripts_accept_no_data_urls():
    sources = _directives(web_security.PAGE_CSP)["script-src"]
    assert "data:" not in sources and "blob:" not in sources, sources


_URL = re.compile(r"""["'](https://[^"'\s]+\.(?:js|css))["']""")
_CONCAT = re.compile(
    r"""["'](https://[^"'\s]+@)["']\s*\+\s*([A-Z_][A-Z0-9_]*)\s*\+\s*["']([^"']+)["']""")


def _external_assets() -> dict[str, str]:
    """``url -> file`` for every script/stylesheet URL in the dashboard source,
    with ``"...@" + VERSION + "/path"`` concatenations resolved."""
    found: dict[str, str] = {}
    for src in sorted(DASHBOARD.glob("*.py")):
        text = src.read_text(encoding="utf-8")
        for match in _URL.finditer(text):
            found[match.group(1)] = src.name
        consts = dict(re.findall(r"^\s*([A-Z_][A-Z0-9_]*)\s*=\s*['\"]([^'\"]+)['\"]",
                                 text, re.M))
        for head, const, tail in _CONCAT.findall(text):
            found[head + consts.get(const, "?") + tail] = src.name
    return found


def test_the_scan_sees_the_dashboards_foreign_scripts():
    urls = _external_assets()
    assert any("xterm" in u for u in urls), urls
    assert any("mermaid" in u for u in urls), urls


def test_every_foreign_script_or_stylesheet_is_pinned_hashed_and_named_by_the_csp():
    d = _directives(web_security.PAGE_CSP)
    allowed = set(_foreign(d["script-src"])) | set(_foreign(d["style-src"]))
    problems = []
    for url, name in sorted(_external_assets().items()):
        if not re.search(r"@\d+\.\d+\.\d+/", url):
            problems.append(f"{name}: {url} is not pinned to a version")
        if url not in allowed:
            problems.append(f"{name}: {url} is not named by the page CSP")
    assert not problems, "\n".join(problems)
    from delfin.dashboard import pinned_cdn
    for asset in pinned_cdn.ASSETS:
        assert asset.sri.startswith("sha384-") and len(asset.sri) > 60, asset
        assert asset.url in allowed, asset


def test_no_viewer_falls_back_to_an_unpinned_3dmol():
    hits = [p.name for p in DASHBOARD.glob("*.py")
            if "3Dmol.org/build" in p.read_text(encoding="utf-8")]
    assert hits == [], hits


def test_the_mathjax_stand_in_is_served_from_the_dashboards_origin():
    launcher = (REPO / "delfin" / "cli_voila.py").read_text(encoding="utf-8")
    uses_data_url = '"mathjax_url": "data:' in launcher
    names_the_route = "NO_TYPESETTING_ROUTE" in launcher
    assert not uses_data_url and names_the_route
