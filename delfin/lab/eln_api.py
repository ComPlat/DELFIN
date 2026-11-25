"""A client for the Chemotion ELN's REST API, and the answer it gives.

Input: a base URL, a collection of readings the caller asked for, and a
function that resolves the credential by NAME. Output: parsed JSON, or a
:class:`ClientError` that names what the server said and nothing about the
credential.

What this module is not. It holds no token, and it cannot mint one: the
credential is fetched through ``credential_fn(name)`` at the moment a
request needs it, so a URL, a log line, a refusal and a ``repr`` of any
object here carry no secret. ``display_url`` is the only way an address
leaves this module.

The server facts, read from the API's own source on 2026-11-24 (the
``@main`` of ComPlat/chemotion_ELN; the strings are cited in
:data:`SOURCES`):

* The API is Grape under ``/api`` + version ``v1`` (``app/api/api.rb``:
  ``prefix :api`` / ``version 'v1'``), and it mounts ``SampleAPI``,
  ``ReactionAPI``, ``ResearchPlanAPI``, ``CollectionAPI``,
  ``SearchAPI``, ``AttachmentAPI`` and ``AuthenticationAPI``.
* ``POST /api/v1/authentication/token`` with ``username`` and
  ``password`` answers ``{token: <chemtoken_...>}`` (unwrapped, no root
  key), and every other call authenticates with
  ``Authorization: Bearer <token>`` (``api.rb``,
  ``detect_current_user_from_api_token``).
* List endpoints are paginated by ``grape-kaminari``: the request carries
  ``page`` and ``per_page``, the answer is a bare JSON ARRAY, and the
  TOTAL and the page count arrive in the ``X-Total`` and
  ``X-Total-Pages`` response headers (``Grape::Kaminari``; the
  ``paginate per_page: 7, max_per_page: 100`` declarations in
  ``research_plan_api.rb`` and ``sample_api.rb``).
* Searches are POSTs: ``POST /api/v1/search/<plural>`` with the query in
  the body (``search_api.rb``). Their answer is
  ``{"<plural>": [...], "total_elements": N}``.
* A refusal is JSON ``{"error": "<text>"}`` with an HTTP status
  (``error!('403 Forbidden', 403)`` throughout; some endpoints answer
  ``{"error": ..., "code": ...}``), so the TEXT of the error is the
  server's own sentence and is quoted back -- never paraphrased into
  "request failed".

What was NOT verified, and is therefore not modelled: the exact shape of
a create body beyond the envelope (`name`/`body`/`collection_id` plus,
for reactions, the `materials` the source shows it reads), what a
materials entry must carry for the server to accept it, and the
`attachments` mapping request's body. :data:`UNVERIFIED` lists them and
:func:`draft_body` refuses instead of guessing.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from typing import Any, Callable, Optional
from urllib.parse import urlencode, urlsplit, urlunsplit

#: Where each fact above was read. Recorded rather than re-fetched: a
#: checkout of that commit is what a completeness pass re-derives them
#: from.
SOURCES = {
    "api.root": "ComPlat/chemotion_ELN app/api/api.rb @main 2026-11-24 "
                "(prefix :api, version 'v1', the mount list)",
    "authentication": "app/api/chemotion/authentication_api.rb "
                      "(POST /authentication/token -> {token: ...})",
    "auth.header": "app/api/api.rb detect_current_user_from_api_token "
                   "(Bearer; the token must start with 'chemtoken')",
    "pagination": "Grape::Kaminari; paginate per_page/max_per_page in "
                  "sample_api.rb and research_plan_api.rb (bare array + "
                  "X-Total / X-Total-Pages headers)",
    "search": "app/api/chemotion/search_api.rb (POST /api/v1/search/<plural>, "
              "body = query, answer {<plural>: [...], total_elements: N})",
    "samples": "app/api/chemotion/sample_api.rb and "
               "app/api/entities/sample_entity.rb",
    "reactions": "app/api/chemotion/reaction_api.rb and "
                 "app/api/entities/reaction_entity.rb",
    "research_plans": "app/api/chemotion/research_plan_api.rb and "
                      "app/api/entities/research_plan_entity.rb",
    "collections": "app/api/chemotion/collection_api.rb "
                   "(GET /collections -> {own: [...], shared_with_me: [...]})",
    "errors": "error!('...', status) throughout; body {error: text}",
}

#: The shapes NOT established from source. A caller that would depend on
#: one of these is refused with this text rather than served a guess.
UNVERIFIED = {
    "reaction.create": "the full set of attributes POST /reactions requires "
                       "besides the material lists; no draft is written "
                       "without them",
    "sample.create": "the container/analyses nesting a fresh sample needs; "
                     "a draft is written with name + collection only",
    "attachment.attach": "the body that maps uploaded keys onto an element "
                         "(modelled from the official client, not from the "
                         "route)",
}

#: The default per-request timeout, seconds. The server has no rate limit
#: (eln_routes.SERVER_FACTS), so a hung request is bounded by this and by
#: the page cap in :meth:`Client._paged`.
DEFAULT_TIMEOUT_S = 30.0

#: The largest page the server will serve (``max_per_page: 100``).
MAX_PER_PAGE = 100

#: The credential's NAME in DELFIN's credential store. The value is never
#: written here; this is the key a person stores under.
TOKEN_NAME = "CHEMOTION_API_TOKEN"

#: The name the username is stored under, for the token-minting call.
USER_NAME = "CHEMOTION_USER"

#: A SHA-256 hex digest, the shape an approval's digest has. Used to refuse
#: ``send_approved`` anything that is not one.
_DIGEST_RE = re.compile(r"^[0-9a-f]{64}$")


class ClientError(RuntimeError):
    """A call that did not produce an answer, or produced a refusal.

    ``status`` is the HTTP status when there was one, ``detail`` the
    server's own error text when it said anything, ``route`` the verb and
    path (never the URL with userinfo). The message never contains the
    credential: the transport strips it before raising, and
    ``display_url`` is the only address that reaches a message.
    """

    def __init__(self, message: str, *, status: Optional[int] = None,
                 detail: str = "", route: str = "") -> None:
        super().__init__(message)
        self.status = status
        self.detail = detail
        self.route = route


def display_url(url: str) -> str:
    """``url`` without userinfo or query, for a message a person reads.

    The same contract as ``delfin.agent.endpoint_gate.display_url``: a
    credential pasted into a URL must not reach a log, an error or a
    transcript. Reimplemented here rather than imported so this module
    keeps no dependency on the agent package (``delfin/lab`` is imported
    on machines with no agent).
    """
    try:
        parts = urlsplit(str(url or ""))
        host = parts.hostname or ""
        if parts.port:
            host = f"{host}:{parts.port}"
        return urlunsplit((parts.scheme, host, parts.path, "", ""))
    except Exception:                                   # noqa: BLE001
        return ""


@dataclass
class Page:
    """One page of a list endpoint, with what the server said about it."""

    rows: list
    total: Optional[int] = None
    total_pages: Optional[int] = None
    page: int = 1
    per_page: int = 0


@dataclass
class Answer:
    """One call's outcome: the parsed body plus the headers a caller needs.

    Only the headers :data:`KEPT_HEADERS` names are kept. Copying every
    response header into a return value is how a ``Set-Cookie`` or an
    ``Authorization`` echo ends up in a transcript.
    """

    route: str
    status: int
    body: Any = None
    headers: dict = field(default_factory=dict)

    @property
    def total(self) -> Optional[int]:
        return _as_int(self.headers.get("x-total"))

    @property
    def total_pages(self) -> Optional[int]:
        return _as_int(self.headers.get("x-total-pages"))


#: The response headers a caller may see. ``x-total``/``x-total-pages``
#: are Grape::Kaminari's; ``content-type`` decides JSON vs text.
KEPT_HEADERS = ("x-total", "x-total-pages", "x-page", "x-per-page",
                "content-type")


def _as_int(value: Any) -> Optional[int]:
    try:
        return int(str(value))
    except (TypeError, ValueError):
        return None


def _kept(headers: Any) -> dict:
    out: dict = {}
    try:
        items = headers.items()
    except AttributeError:
        return out
    for key, value in items:
        low = str(key).lower()
        if low in KEPT_HEADERS:
            out[low] = str(value)
    return out


# ---------------------------------------------------------------------------
# The transport: one seam, so a fake server can be the far end
# ---------------------------------------------------------------------------

def _urlopen_send(method: str, url: str, body: Optional[bytes],
                  headers: dict, timeout: float) -> tuple:
    """``(status, headers, body bytes)`` over urllib. The real transport.

    A loopback target (the fake server in the tests, an instance in
    Docker on the same host) is sent WITHOUT the environment's HTTP
    proxy, which would otherwise intercept it; a remote ELN keeps the
    default opener so an ``http_proxy`` is honoured. The same split
    ``delfin.lab.chemotion`` makes, for the same measured reason.
    """
    import urllib.error
    import urllib.request

    req = urllib.request.Request(url, data=body, headers=headers or {},
                                 method=method)
    opener = (urllib.request.build_opener(urllib.request.ProxyHandler({}))
              if is_loopback(url) else None)
    try:
        with (opener or urllib.request).open(req, timeout=timeout) as resp:
            return int(resp.status), _kept(resp.headers), resp.read()
    except urllib.error.HTTPError as exc:
        # The server answered; its status and body are the answer.
        return int(exc.code), _kept(getattr(exc, "headers", None)), \
            (exc.read() or b"")
    except OSError as exc:
        raise ClientError(
            f"{method} {display_url(url)}: {type(exc).__name__}",
            route=f"{method} {_path_of(url)}") from None


def is_loopback(url: str) -> bool:
    """Is ``url`` this machine? Decides whether a proxy may be used."""
    host = (urlsplit(str(url or "")).hostname or "").lower().rstrip(".")
    return host in ("localhost", "127.0.0.1", "::1") or host.startswith("127.")


def _path_of(url: str) -> str:
    return urlsplit(str(url or "")).path or "/"


# ---------------------------------------------------------------------------
# The client
# ---------------------------------------------------------------------------

class Client:
    """Calls the ELN's API. Sends nothing the route gate has not allowed.

    ``credential_fn(name) -> str`` is the ONLY way a secret enters. It is
    called lazily, once per client, and the result is held in
    ``self._token`` (``_secret`` in a ``repr``). ``transport`` is the
    seam the tests answer with a fake server; the default is
    :func:`_urlopen_send`.

    Every call is decided by :func:`delfin.lab.eln_routes.decide_route`
    BEFORE it is sent. The gate allows reads only, so a client built on it
    cannot write at all -- which is the point: :class:`DraftWriter` is the
    one writer, and it takes a person's approval per call.
    """

    def __init__(self, base_url: str, credential_fn: Callable[[str], str],
                 *, transport: Optional[Callable] = None,
                 timeout: float = DEFAULT_TIMEOUT_S) -> None:
        self.base_url = str(base_url or "").rstrip("/")
        self._credential_fn = credential_fn
        self._transport = transport or _urlopen_send
        self._timeout = float(timeout)
        self._token: Optional[str] = None

    def __repr__(self) -> str:                          # no secrets
        return (f"Client(base_url={display_url(self.base_url)!r}, "
                f"token={'set' if self._token else 'unset'})")

    # -- the credential ---------------------------------------------------

    def token(self) -> str:
        """The bearer token, read from the store on first use.

        Raised as a :class:`ClientError` naming the credential's NAME (not
        its value, and not the person) when the store has none: a caller
        that has to tell someone what to store is told the key, never the
        secret.
        """
        if self._token is None:
            value = ""
            try:
                value = str(self._credential_fn(TOKEN_NAME) or "")
            except Exception:                           # noqa: BLE001
                value = ""
            if not value:
                raise ClientError(
                    f"no {TOKEN_NAME} in the credential store; store it with "
                    "`delfin-agent key set` (the value never goes in a prompt)")
            self._token = value
        return self._token

    # -- one request ------------------------------------------------------

    def _send(self, method: str, path: str, *, query: dict | None = None,
              body: Any = None) -> Answer:
        """Gate, send, parse. The one READ path.

        The route gate allows reads only, so this raises for a write. That
        is the whole design: a caller holding this object cannot change
        anything on the server, and the one write DELFIN makes is made by
        :meth:`send_approved`, which requires a person's grant to exist
        first. Two methods rather than a flag, because a flag is a
        parameter a caller can set.
        """
        from . import eln_routes

        decision = eln_routes.decide_route(method, path, body if body else None)
        if not decision.allowed:
            raise ClientError(
                f"refused before sending: {decision.reason}",
                route=f"{method} {path}")
        return self._transfer(method, path, query=query, body=body)

    def send_approved(self, method: str, path: str, *, digest: str,
                      body: Any = None) -> Answer:
        """Send one WRITE that a person has already approved.

        ``digest`` is the approval's digest, and it is required: a caller
        that has not consumed a grant for this exact request has no digest
        to pass, so this method cannot be reached from a fresh
        conversation, from a tool argument, or from a permission mode.
        :func:`delfin.lab.eln_write.create_draft` is its only caller, and
        it passes the digest it just consumed.

        The route gate is not consulted again, and that is deliberate
        rather than a hole. The gate answers "may this be sent WITHOUT a
        human decision", and its answer for every write is no; this is the
        case its docstring names as a later layer's work -- "a person
        confirms, then the transport sends" -- with the confirmation made
        concrete in ``eln_write``. What this method DOES still require is
        that the table knows the route: a write to a path
        :data:`eln_routes.ELN_ROUTES` has never listed is refused here
        too, so the table stays the single list of what may be reached.
        """
        from . import eln_routes

        if not _DIGEST_RE.match(str(digest or "")):
            raise ClientError(
                "send_approved needs the digest of a consumed approval; an "
                "ELN write is sent only after a person approved it",
                route=f"{method} {path}")
        row = eln_routes.match_route(method, path)
        if row is None or not row.writes:
            raise ClientError(
                f"not a write the route table knows: {method} {path}",
                route=f"{method} {path}")
        return self._transfer(method, path, body=body)

    def _transfer(self, method: str, path: str, *, query: dict | None = None,
                  body: Any = None) -> Answer:
        """Put one already-decided request on the wire and parse the answer."""
        url = self.base_url + path
        if query:
            url = f"{url}?{urlencode(query)}"
        headers = {"Accept": "application/json"}
        payload = None
        if body is not None:
            payload = json.dumps(body).encode("utf-8")
            headers["Content-Type"] = "application/json"
        # Every route the gate allows is an authenticated read: the table
        # carries no unauthenticated entry (even the token mint it refuses
        # is a write). So the bearer header goes on, or the call does not.
        headers["Authorization"] = f"Bearer {self.token()}"
        status, raw_headers, raw = self._transport(
            method, url, payload, headers, self._timeout)
        text = raw.decode("utf-8", "replace") if isinstance(raw, bytes) else \
            str(raw or "")
        parsed: Any = None
        if text.strip():
            try:
                parsed = json.loads(text)
            except ValueError:
                parsed = None
        if status >= 400:
            # The server's own sentence, quoted. `error!('403 Forbidden')`
            # and the `{error:, code:}` form both land here.
            detail = ""
            if isinstance(parsed, dict):
                detail = str(parsed.get("error") or parsed.get("message") or "")
            if not detail:
                detail = text.strip()[:300]
            raise ClientError(
                f"{method} {path}: HTTP {status}"
                + (f" - {detail}" if detail else ""),
                status=status, detail=detail, route=f"{method} {path}")
        return Answer(route=f"{method} {path}", status=status, body=parsed,
                      headers=raw_headers)

    # -- reads ------------------------------------------------------------

    def collections(self) -> dict:
        """``GET /api/v1/collections`` -> ``{"own": [...], "shared_with_me": []}``.

        The server's own envelope, passed through. Resolving a collection
        NAME to its id is :func:`resolve_collection`'s job.
        """
        out = self._send("GET", "/api/v1/collections").body
        return out if isinstance(out, dict) else {
            "own": [], "shared_with_me": []}

    def search(self, element: str, query: dict, *,
               per_page: int = MAX_PER_PAGE, page: int = 1) -> Page:
        """``POST /api/v1/search/<element>`` -- the read path the gate allows.

        ``element`` is the plural the route spells (``samples``,
        ``reactions``, ``research_plans``, ``wellplates``, ``screens``).
        The answer is ``{"<element>": [...], "total_elements": N}``;
        ``total_elements`` becomes :attr:`Page.total`.
        """
        element = _element_name(element)
        body = dict(query or {})
        body.setdefault("page", page)
        body.setdefault("per_page", per_page)
        answer = self._send("POST", f"/api/v1/search/{element}", body=body)
        data = answer.body if isinstance(answer.body, dict) else {}
        rows = data.get(element)
        if rows is None:
            for value in data.values():
                if isinstance(value, list):
                    rows = value
                    break
        total = _as_int(data.get("total_elements")) or answer.total
        return Page(rows=list(rows or []), total=total,
                    total_pages=answer.total_pages, page=page,
                    per_page=per_page)

    def search_all(self, element: str, query: dict, *, max_pages: int = 20,
                   per_page: int = MAX_PER_PAGE) -> list:
        """Follow the pages of one search, up to ``max_pages``.

        A paginated read that stops at page one is a wrong answer waiting
        to be quoted ("none of these has an analysis" when page two had
        one), so the loop is here rather than at each caller. It stops at
        the server's ``total_pages`` when it says, and at ``max_pages``
        regardless -- a broken header cannot spin forever.
        """
        rows: list = []
        page = 1
        while page <= max(1, int(max_pages)):
            got = self.search(element, query, per_page=per_page, page=page)
            rows.extend(got.rows)
            if got.total_pages is not None and page >= got.total_pages:
                break
            if got.total is not None and len(rows) >= got.total:
                break
            if len(got.rows) < per_page:
                break
            page += 1
        return rows

    def element(self, element: str, element_id: Any) -> dict:
        """One element by id, from the search endpoint (a VERIFIED route).

        ``GET /api/v1/<element>/<id>`` is NOT used: no per-id sample read is
        established from the API source, and a path that 404s on the real
        server is worse than a slower one that works. The search body
        carries the ids (``search_api.rb``'s ``:search_params`` are the
        ``ui_state`` shape -- ``checkedIds``/``included_ids`` -- so a
        search scoped to one id is a read the table already allows), and
        the row is taken from that answer.

        Raises :class:`ClientError` when the id is not in the answer: an
        empty dict returned as if it were an element is how "no analyses"
        gets reported for a sample that was never read.
        """
        element = _element_name(element)
        wanted = _as_int(element_id)
        if wanted is None:
            raise ClientError(f"element id must be an integer, not "
                              f"{element_id!r}", route=f"search {element}")
        body = {"checkedAll": False, "checkedIds": [wanted],
                "uncheckedIds": [], "included_ids": [wanted],
                "excluded_ids": []}
        page = self.search(element, body, per_page=MAX_PER_PAGE)
        for row in page.rows:
            if _as_int(row.get("id")) == wanted:
                return row
        raise ClientError(
            f"no {_singular(element)} with id {wanted} in the answer; it may "
            "not exist, or this account may not read it",
            status=404, route=f"search {element}")


def _element_name(element: str) -> str:
    """Normalise ``sample``/``samples``/``Sample`` to the route's plural."""
    text = str(element or "").strip().lower().replace("-", "_")
    aliases = {"sample": "samples", "reaction": "reactions",
               "researchplan": "research_plans",
               "research_plan": "research_plans",
               "wellplate": "wellplates", "screen": "screens"}
    if text in aliases:
        return aliases[text]
    if not text.endswith("s"):
        return aliases.get(text, text + "s")
    return text


def _singular(element: str) -> str:
    if element.endswith("ies"):
        return element[:-3] + "y"
    return element[:-1] if element.endswith("s") else element


def resolve_collection(client: Client, wanted: str) -> dict:
    """The collection called ``wanted``, from the server's own list.

    Matched on ``label`` (the name a person types) then ``id``. Refuses
    with the list of names it DID see rather than "not found": the
    difference between "you spelled it wrong" and "you cannot see it" is
    the difference between a fixable answer and a wrong one.
    """
    listing = client.collections()
    rows = list(listing.get("own") or []) + \
        list(listing.get("shared_with_me") or [])
    key = str(wanted or "").strip()
    for row in rows:
        if not isinstance(row, dict):
            continue
        if key and key in (str(row.get("label") or ""),
                           str(row.get("name") or ""),
                           str(row.get("id") or "")):
            return row
    names = sorted({str(r.get("label") or r.get("name") or r.get("id") or "")
                    for r in rows if isinstance(r, dict)})
    raise ClientError(
        f"no collection named {key!r}; the server lists: "
        + (", ".join(n for n in names if n) or "(none)")
        + " (a collection you cannot read is not listed either)")
