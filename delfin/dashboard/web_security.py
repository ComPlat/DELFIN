"""HTTP-level protections of the dashboard server.

Installed by the ``delfin.dashboard.server_guard`` server extension into
the running jupyter-server. Three parts:

``files route``
    Voila registers ``/voila/files/(.*)`` with ``AllowListFileHandler``, a
    plain tornado ``StaticFileHandler``: no authentication, the content
    type taken from the file name, links followed out of the root, hidden
    paths served. :func:`install_files_route` replaces the handler of that
    route (in place, in the router, so the order in which extensions load
    does not matter) with :class:`DelfinFileHandler`, built on
    jupyter-server's ``AuthenticatedFileHandler``. It keeps Voila's
    allow-list and deny-list, and adds:

    * authentication and Jupyter's XSRF/Referer check on every GET/HEAD;
    * the real path must stay inside the real root (no link out of it);
    * no path component may start with ``.``, except below an explicitly
      named directory of DELFIN's own (the drawing editor's launch copy);
    * outside those DELFIN directories only raster images and PDFs are
      served inline. Anything else, HTML/SVG/JS above all, is sent as
      ``text/plain`` with ``Content-Disposition: attachment`` and a
      ``sandbox`` CSP, so user data never runs as a page on the
      dashboard's origin, even if the allow-list is widened by config.

    Rejected: serving user HTML inline under ``sandbox allow-scripts`` (what
    Jupyter's ``/files/`` does). It confines the page to an opaque origin,
    but still renders attacker-chosen markup inside the dashboard's frame
    tree; nothing in the UI needs it.

``headers``
    :data:`SECURITY_HEADERS` (the page CSP among them, :func:`_page_csp`)
    are merged into ``settings["headers"]``, which
    every ``JupyterHandler`` applies; :class:`SecurityHeadersTransform` adds
    the same headers to every other response (Voila's own static handlers
    are plain tornado handlers and set none).

``session cookie``
    The login cookie gets ``SameSite=Lax`` and a bounded lifetime
    (:func:`login_days`, default 7 days, enforced on read through
    ``max_age_days``); the XSRF cookie gets ``SameSite=Lax`` as well.

Every part is best effort towards the caller (the server must start), but
each failure is logged at WARNING; the files route, if it cannot be
replaced by the authenticated handler, is replaced by one that refuses
every request (fail closed).
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterable
from typing import Any

from delfin.dashboard import pinned_cdn as _pinned_cdn


def _page_csp() -> str:
    """The Content-Security-Policy of every dashboard page.

    Measured in headless Chromium (every tab opened, every nested tab bar
    clicked, the drawing editor's frame included) against a report-only
    policy of ``'self'`` alone; each source below answers a violation that
    policy produced, or a load the page makes on demand:

    ``script-src 'unsafe-inline' 'unsafe-eval'``
        Inline ``<script>`` in widget HTML outputs and the Voila template,
        inline handlers (``<img onerror=...>`` is how an HTML widget runs
        code), and ``eval`` in Voila's own bundle, which is how every
        ``Javascript`` output -- the dashboard's start-up script among
        them -- is run. Without them the page does not start.
    ``script-src <pinned files>``
        xterm, its fit add-on and mermaid, each by its full versioned path
        and loaded with an SRI hash (:mod:`pinned_cdn`). Nothing else from
        cdn.jsdelivr.net: py3Dmol's own loader (``3dmol@<version>``) never
        runs because the bundled 3Dmol is defined first, and the CSP keeps
        it from running if it did.
    ``style-src 'unsafe-inline' <xterm.css>``
        Several hundred inline style attributes and ``<style>`` elements
        from ipywidgets, Voila and the dashboard.
    ``img-src data: blob:``, ``font-src data:``
        Inline images (logo, plots, the MANTA animation) and the fonts the
        bundles embed.
    ``worker-src blob:``
        The force-field worker is built from a ``Blob``.
    ``frame-src data: blob:``
        The literature viewer falls back to a ``data:`` PDF.
    ``connect-src 'self'``
        The kernel websocket, the contents API, the beacon. A browser that
        implements CSP Level 3 matches ``ws:``/``wss:`` to the page's own
        origin under ``'self'``.

    Not allowed, and measured to cost nothing: Voila's template links
    Font Awesome from ``cdn.jsdelivr.net/npm/@fortawesome/fontawesome-free@^5``
    (a version range, no SRI). The same font is bundled in Voila's and the
    widget manager's own scripts and renders from there.
    """
    scripts = " ".join(_pinned_cdn.urls("script"))
    styles = " ".join(_pinned_cdn.urls("style"))
    return "; ".join([
        f"script-src 'self' 'unsafe-inline' 'unsafe-eval' {scripts}",
        f"style-src 'self' 'unsafe-inline' {styles}",
        "img-src 'self' data: blob:",
        "font-src 'self' data:",
        "connect-src 'self'",
        "worker-src 'self' blob:",
        "frame-src 'self' data: blob:",
        "media-src 'self' data: blob:",
        "manifest-src 'self'",
        "object-src 'none'",
        "base-uri 'self'",
        "form-action 'self'",
        "frame-ancestors 'self'",
    ])


#: What every dashboard page may load and run (:func:`_page_csp`).
#:
#: No ``default-src``: Jupyter appends ``default-src 'none'`` to the policy
#: of its API responses, and the first occurrence of a directive is the one
#: a browser applies.
PAGE_CSP = _page_csp()

#: A same-origin, empty script the launcher hands Voila as ``mathjax_url``
#: (the dashboard typesets nothing; see ``cli_voila``), so ``script-src``
#: needs no ``data:``.
NO_TYPESETTING_ROUTE = "/delfin/no-typesetting.js"

#: For a file that is a download or a passive image: nothing in it runs.
FILE_CSP = "default-src 'none'; frame-ancestors 'self'; sandbox"

#: For a PDF shown in an iframe. ``sandbox`` is left out: browsers refuse
#: to start their PDF viewer in a sandboxed document.
PDF_CSP = "default-src 'none'; object-src 'self'; frame-ancestors 'self'"

#: ``Referrer-Policy: same-origin``, not ``no-referrer``: Jupyter's file
#: handlers accept a cookie-authenticated GET without an XSRF token only
#: when the Referer names the same host (``JupyterHandler.check_referer``),
#: so ``no-referrer`` turns every image and PDF the page loads into a 403.
#: ``same-origin`` sends nothing to any other site, which is what keeps a
#: token-bearing URL out of third-party logs.
SECURITY_HEADERS = {
    "Content-Security-Policy": PAGE_CSP,
    "X-Content-Type-Options": "nosniff",
    "X-Frame-Options": "SAMEORIGIN",
    "Referrer-Policy": "same-origin",
    "Permissions-Policy": "camera=(), microphone=(), geolocation=()",
}

#: Served inline from user data: they cannot carry script.
PASSIVE_IMAGE_SUFFIXES = frozenset(
    {".png", ".jpg", ".jpeg", ".gif", ".ico", ".webp", ".bmp"})
PDF_SUFFIXES = frozenset({".pdf"})

#: The allow-list ``delfin-voila`` passes to Voila: what the UI requests
#: from ``/voila/files/`` -- images (the MANTA animation), PDFs (the
#: literature viewer), and the drawing editor's launch copy, which is only
#: used when the published ``/static/ketcher`` route could not be set up.
EDITOR_LAUNCH_COPY = ".delfin/ketcher"
FILE_ALLOWLIST = (
    r"(\.delfin/ketcher/.*\.(html|js|css|svg|png|jpg|gif|ico|json|wasm"
    r"|woff|woff2|ttf))|(.*\.(png|jpg|jpeg|gif|ico|webp|pdf))"
)

LOGIN_DAYS_ENV = "DELFIN_DASHBOARD_LOGIN_DAYS"
DEFAULT_LOGIN_DAYS = 7.0


def login_days() -> float:
    """Lifetime of the dashboard's login cookie, in days (bounded 1 h..30 d)."""
    try:
        value = float(os.environ.get(LOGIN_DAYS_ENV, "") or DEFAULT_LOGIN_DAYS)
    except ValueError:
        value = DEFAULT_LOGIN_DAYS
    return min(30.0, max(1.0 / 24.0, value))


def classify(path: str, *, trusted: bool) -> str:
    """How a file is served: ``'page'``, ``'image'``, ``'pdf'`` or ``'download'``.

    *trusted* is True for a file below one of DELFIN's own directories.
    """
    if trusted:
        return "page"
    suffix = os.path.splitext(path)[1].lower()
    if suffix in PASSIVE_IMAGE_SUFFIXES:
        return "image"
    if suffix in PDF_SUFFIXES:
        return "pdf"
    return "download"


def _inside(path: str, root: str) -> bool:
    return path == root or path.startswith(root.rstrip(os.sep) + os.sep)


_HANDLER_CLASS: Any = None


def handler_class():
    """The authenticated replacement for Voila's ``AllowListFileHandler``."""
    global _HANDLER_CLASS
    if _HANDLER_CLASS is not None:
        return _HANDLER_CLASS

    from urllib.parse import quote

    from jupyter_server.base.handlers import AuthenticatedFileHandler
    from tornado import web

    class DelfinFileHandler(AuthenticatedFileHandler):
        """``/voila/files/``: authenticated, root-confined, no active content."""

        def initialize(self, path, allowlist=(), denylist=(),
                       trusted_dirs=(EDITOR_LAUNCH_COPY,),
                       default_filename=None, **_ignored):
            self.allowlist = list(allowlist or ())
            self.denylist = list(denylist or ())
            self.trusted_dirs = tuple(trusted_dirs or ())
            self._delfin_kind = "download"
            web.StaticFileHandler.initialize(
                self, path=path, default_filename=default_filename)

        def get_absolute_path(self, root, path):
            allowed = any(re.fullmatch(p, path, re.IGNORECASE)
                          for p in self.allowlist)
            denied = any(re.fullmatch(p, path, re.IGNORECASE)
                         for p in self.denylist)
            if not allowed:
                raise web.HTTPError(403, "File not allowlisted")
            if denied:
                raise web.HTTPError(403, "File denylisted")
            return web.StaticFileHandler.get_absolute_path(root, path)

        def validate_absolute_path(self, root, absolute_path):
            # tornado's own check: inside the root (lexically), exists,
            # directory -> default file. Jupyter's hidden-file rule is
            # replaced by the stricter one below.
            checked = web.StaticFileHandler.validate_absolute_path(
                self, root, absolute_path)
            if checked is None:
                return None
            real_root = os.path.realpath(root)
            real = os.path.realpath(checked)
            if not _inside(real, real_root):
                raise web.HTTPError(403, "outside the served directory")
            rel = os.path.relpath(real, real_root).replace(os.sep, "/")
            trusted = any(rel == d or rel.startswith(d.rstrip("/") + "/")
                          for d in self.trusted_dirs)
            if not trusted and any(part.startswith(".")
                                   for part in rel.split("/")):
                raise web.HTTPError(404)
            self._delfin_kind = classify(rel, trusted=trusted)
            return checked

        def get_content_type(self):
            if self._delfin_kind == "download":
                return "text/plain; charset=utf-8"
            return super().get_content_type()

        def set_headers(self):
            super().set_headers()
            kind = self._delfin_kind
            if kind == "page":
                csp = PAGE_CSP
            elif kind == "pdf":
                csp = PDF_CSP
            else:
                csp = FILE_CSP
            self.set_header("Content-Security-Policy", csp)
            self.set_header("X-Content-Type-Options", "nosniff")
            if kind == "download":
                name = os.path.basename(self.absolute_path or "") or "download"
                self.set_header(
                    "Content-Disposition",
                    "attachment; filename*=UTF-8''" + quote(name, safe=""))

    _HANDLER_CLASS = DelfinFileHandler
    return DelfinFileHandler


def refusing_handler_class():
    """What the files route becomes when the authenticated handler cannot be
    built: every request answered 403."""
    from tornado import web

    class RefuseFiles(web.RequestHandler):
        def prepare(self):
            raise web.HTTPError(403, "file route unavailable")

    return RefuseFiles


def _rules(router: Any) -> Iterable[Any]:
    """Every tornado ``Rule`` below *router*, depth first."""
    for rule in list(getattr(router, "rules", ()) or ()):
        yield rule
        target = getattr(rule, "target", None)
        if hasattr(target, "rules"):
            yield from _rules(target)


def install_files_route(web_app: Any, log: Any = None) -> int:
    """Replace every route served by Voila's ``AllowListFileHandler``.

    Returns how many routes were replaced. Idempotent: a route already
    served by :class:`DelfinFileHandler` is left as it is.
    """
    try:
        from voila.static_file_handler import AllowListFileHandler
    except Exception:
        return 0
    try:
        replacement = handler_class()
        failure = None
    except Exception as exc:          # pragma: no cover - import failure
        replacement = refusing_handler_class()
        failure = exc
    replaced = 0
    for router in (getattr(web_app, "default_router", None),
                   getattr(web_app, "wildcard_router", None)):
        for rule in _rules(router):
            target = getattr(rule, "target", None)
            if isinstance(target, type) and issubclass(target, AllowListFileHandler):
                rule.target = replacement
                replaced += 1
    if failure is not None and log is not None:
        log.warning("[delfin] /voila/files/ is refused: the authenticated "
                    "file handler could not be built (%s)", failure)
    return replaced


def install_headers(web_app: Any) -> None:
    """Default security headers on every response of *web_app*."""
    settings = web_app.settings
    headers = dict(settings.get("headers") or {})
    for name, value in SECURITY_HEADERS.items():
        headers.setdefault(name, value)
    settings["headers"] = headers
    if SecurityHeadersTransform not in getattr(web_app, "transforms", []):
        web_app.add_transform(SecurityHeadersTransform)


def install_cookie_policy(serverapp: Any) -> None:
    """``SameSite=Lax`` and a bounded lifetime for the login and XSRF cookies."""
    days = login_days()
    provider = getattr(serverapp, "identity_provider", None)
    if provider is not None:
        options = dict(provider.cookie_options or {})
        options.setdefault("samesite", "Lax")
        options.setdefault("expires_days", days)
        provider.cookie_options = options
        read = dict(provider.get_secure_cookie_kwargs or {})
        read.setdefault("max_age_days", days)
        provider.get_secure_cookie_kwargs = read
    settings = serverapp.web_app.settings
    xsrf = dict(settings.get("xsrf_cookie_kwargs") or {})
    xsrf.setdefault("samesite", "Lax")
    settings["xsrf_cookie_kwargs"] = xsrf


def _make_transform():
    from tornado.web import OutputTransform

    class _SecurityHeadersTransform(OutputTransform):
        """Adds :data:`SECURITY_HEADERS` a handler did not set itself."""

        def transform_first_chunk(self, status_code, headers, chunk, finishing):
            for name, value in SECURITY_HEADERS.items():
                if name not in headers:
                    headers[name] = value
            return status_code, headers, chunk

    return _SecurityHeadersTransform


SecurityHeadersTransform = _make_transform()


def install_no_typesetting_route(web_app: Any) -> None:
    """Serve :data:`NO_TYPESETTING_ROUTE`: a constant, empty script."""
    from jupyter_server.base.handlers import JupyterHandler
    from jupyter_server.utils import url_path_join
    from tornado import web

    class NoTypesetting(JupyterHandler):
        @web.authenticated
        def get(self):
            self.set_header("Content-Type", "text/javascript; charset=utf-8")
            self.set_header("Cache-Control", "max-age=86400")
            self.finish("// no typesetting\n")

    base = web_app.settings.get("base_url", "/")
    web_app.add_handlers(".*$", [(url_path_join(base, NO_TYPESETTING_ROUTE),
                                  NoTypesetting)])


def install(serverapp: Any) -> list[str]:
    """Install every part; return the names of the parts that failed."""
    failed: list[str] = []
    log = getattr(serverapp, "log", None)
    for name, step in (
        ("headers", lambda: install_headers(serverapp.web_app)),
        ("cookies", lambda: install_cookie_policy(serverapp)),
        ("typesetting stand-in",
         lambda: install_no_typesetting_route(serverapp.web_app)),
        ("files route", lambda: install_files_route(serverapp.web_app, log)),
    ):
        try:
            step()
        except Exception as exc:
            failed.append(name)
            if log is not None:
                log.warning("[delfin] web protection %r not installed: %s",
                            name, exc)
    return failed
