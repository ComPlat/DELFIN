"""Allowlist of Chemotion ELN routes: which request may be made, and how that
is decided.

Purpose. Later work will talk to a Chemotion ELN instance -- the integration
target is a throwaway instance in Docker, NEVER a live one. This module is the
gate every request passes BEFORE it could be sent. It decides offline from a
table; nothing here sends, prepares or simulates a request (there is no HTTP
code in this file; a test asserts that as well).

Why not a verb gate. Facts about the server, verified in ComPlat/chemotion_ELN
at commit 3657b15 on 2026-11-19 (:data:`SERVER_FACTS`):

* Its tokens carry NO scopes. There is no read-only token, so the verb says
  nothing about trust.
* A verb-keyed gate is wrong there: searches are POSTs (``/api/v1/search/*``),
  and so are destructive operations (``/api/v1/versions/revert``,
  ``/api/v1/generic_elements/delete_revision``,
  ``/api/v1/molecules/delete_name``). "POST is a write" would refuse the
  searches and "GET is safe" would allow any GET. Both fail on this server.
* The blast radius can live in the BODY: ``DELETE /api/v1/ui_state`` without
  ``currentCollection.id`` in the body targets everything the user owns.
* Updates use ``declared(params, include_missing: false)``: a partial update
  merges -- but defaulted parameters and boiling/melting ranges are reset by
  the update.
* Every element model is ``acts_as_paranoid`` and there is NO undelete via the
  API. A delete here cannot be undone, which is why nothing in this gate
  reclassifies a delete toward "read" or "bounded" without evidence.
* There is no rate limiting at all.

So every entry carries its verdict as DATA: the verb (scoping here is per
verb+route, as the facts above show), whether it writes, and whether the
decision needs the request body. :func:`decide_route` consults nothing else.

The gate. :func:`decide_route` answers ONE question: may this request proceed
WITHOUT a human decision (``allowed``). Reads may; everything else may not --
the gate is read-only by construction, because this repository is foundations
for an integration and no approval flow exists yet. The alternative rejected:
auto-allowing a write whose scope is provably narrow (``refused_write_scoped``
below). With ``acts_as_paranoid`` making every delete final and no undo in the
API, a bounded write is still final; authorising one is a later layer's work
(a person confirms, then the transport sends). This layer only ever answers
"yes" to reads.

Refusals (``reason`` tokens; the tests name each one):

* ``read`` -- allowed read: nothing on the server changes.
* ``refused_write`` -- a write; refused here regardless of its target.
* ``refused_write_scoped`` -- a write whose narrowed scope WAS found in the
  body; still refused, with the scope named for whoever reviews it later.
* ``refused_write_unscoped`` -- a body-scoped write with its scope key
  missing, e.g. ``DELETE /api/v1/ui_state`` without ``currentCollection.id``:
  the blast radius is "everything the user owns".
* ``refused_needs_body`` -- the entry says the decision needs the body and
  the request had none.
* ``refused_not_allowlisted`` -- a (verb, path) no entry covers.
* ``refused_not_parsed`` -- unknown verb, a path this table cannot hold, or a
  body supplied that is not a JSON object.

FAIL CLOSED. A route not in the table is refused. An unparseable request is
refused. Nothing is allowed because it "looks like" an entry: matching is on
the exact path and verb of the entry. One wildcard exists -- a trailing ``*``
segment, copied verbatim from the verified family ``/api/v1/search/*`` --
and it matches EXACTLY ONE additional segment: ``/api/v1/search`` (no scope)
and ``/api/v1/search/a/b`` match nothing. Entries carry route strings as the
server spells them; placeholders like ``/:id`` or regexes were rejected
because they invent spellings the server never routes and widen matching past
what the source states.

``writes`` is per-ROUTE, not per-verb, deliberately: a POST search
(``/api/v1/search/*``, ``writes=False``) and a POST delete
(``.../delete_revision``, ``writes=True``) are both POSTs on this server and
the rows say so. The tests pin both directions.

Scope fields. :attr:`Route.scope_fields` names body keys (dotted for nesting)
that narrow the blast radius; today exactly one entry needs the body for its
verdict (``DELETE /api/v1/ui_state`` -- see the facts above): without
``currentCollection.id`` in the body the delete covers everything the user
owns, with it it is bounded. :func:`decide_route` reads those keys when the
entry says it needs the body and reports what it found in
:attr:`Decision.scope`. A missing scope key never mutes a write toward a
read; it names the unbounded blast radius instead.

Input, output. :func:`decide_route` (verb, path, body) -> :class:`Decision`.
It sends nothing and builds nothing.

Provenance. Route strings and verbs are recorded from ComPlat/chemotion_ELN at
commit 3657b15 (checked in source on 2026-11-19: ``SOURCES`` says for each
entry what its source was; the @3657b15 facts were transcribed from a source
check of that commit and this table is their record), or from the request
sequences DELFIN itself already writes (:mod:`delfin.lab.chemotion`, whose
``SOURCES`` names the endpoint files it read). The table is NOT complete: what
is not in it is refused by rule (fail closed), not documented here. A
completeness pass back over the server's route surface at 3657b15 is the named
next step before any request is made. The AGPL client library question is
settled in LICENSE_GUARD.md: nothing AGPL is imported or copied here, and no
new dependency is added at all (stdlib only).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Optional

#: The facts this table is shaped by: verified in ComPlat/chemotion_ELN at
#: commit 3657b15 on 2026-11-19. A re-derivation of the table from a checkout
#: of that commit is what the completeness pass (see the module docstring) will
#: do; until then these are the authority for WHY the entries are shaped the
#: way they are (per-route classification, body scope, refusal of writes).
SERVER_FACTS = (
    "tokens carry no scopes -- there is no read-only token "
    "(app/api/chemotion/authentication_api.rb)",
    "searches are POSTs under /api/v1/search/*; destructive operations are "
    "POST too: /api/v1/versions/revert, "
    "/api/v1/generic_elements/delete_revision, /api/v1/molecules/delete_name",
    "DELETE /api/v1/ui_state without currentCollection.id in the body "
    "targets everything the user owns",
    "updates use declared(params, include_missing: false): a partial update "
    "merges, but defaulted parameters and boiling/melting ranges are reset",
    "every element model is acts_as_paranoid and there is no undelete via "
    "the API",
    "there is no rate limiting at all",
)

#: What each entry's route string was taken from. ``@3657b15`` = a source check
#: of ComPlat/chemotion_ELN at that commit on 2026-11-19 (recorded, not
#: re-read while writing this file: the completeness pass names re-deriving
#: it). ``exporter`` = the request sequence DELFIN itself writes today
#: (:func:`delfin.lab.chemotion.request_sequence`), whose ``SOURCES`` names the
#: endpoint files behind it.
SOURCES = {
    "@3657b15": "ComPlat/chemotion_ELN @3657b15, checked in source "
                "2026-11-19",
    "exporter": "delfin.lab.chemotion request sequence (its SOURCES names "
                "app/api/api.rb and the endpoint files it read)",
}


@dataclass(frozen=True)
class Route:
    """One allowlisted (verb, path). One entry per (verb, path) pair.

    ``verb``       the HTTP verb, uppercase.
    ``path``       the route string exactly as the server spells it, minus
                   query string; a trailing ``*`` segment is the family
                   wildcard of ``/api/v1/search/*`` and matches exactly one
                   more segment.
    ``writes``     True when the call changes state on the server -- per
                   ROUTE, not per verb (see the module docstring).
    ``needs_body`` the decision needs the request body to classify this
                   request's blast radius (exactly the entries whose radius
                   lives in the body: ``DELETE /api/v1/ui_state``).
    ``scope_fields`` body keys (dotted for nesting) that narrow the radius;
                   empty when nothing in the body narrows it.
    ``reason``     what the entry is for, in one clause. A refusal echoes the
                   route and this clause together.
    ``source``     which SOURCES class the route string came from: "@3657b15"
                   or "exporter".
    """

    verb: str
    path: str
    writes: bool
    needs_body: bool = False
    scope_fields: tuple = ()
    reason: str = ""
    source: str = "@3657b15"


@dataclass
class Decision:
    """What :func:`decide_route` concluded about one request.

    ``allowed``       may this request proceed without a human decision. True
                      only for reads; the gate is read-only by construction.
    ``writes``        the entry's classification, copied through even when the
                      request is refused -- a caller's audit shows what it was
                      going to send. For a refusal from parsing or lookup (no
                      entry matched) it is True conservatively: an audit reads
                      "this could be a write", never "safe".
    ``body_consulted`` the entry needed the body and a dict body was present,
                      so the decision was made from it (``scope`` says what
                      was found).
    ``reason``        "<token>:<facts>" per the module docstring's table.
    ``route``         the matched entry, or None for a refusal before lookup.
    ``scope``         the scope key/value pairs actually present in the body.
    """

    allowed: bool
    verb: str
    path: str
    writes: bool
    body_consulted: bool = False
    reason: str = ""
    route: Optional[Route] = None
    scope: dict = field(default_factory=dict)


#: Verbs Grape routes here carry. An unknown verb is unparseable, not unknown.
VERBS = frozenset({"GET", "POST", "PUT", "PATCH", "DELETE"})

#: Refusal/read tokens, exactly as the tests name them.
REASON_READ = "read"
REASON_WRITE = "refused_write"
REASON_WRITE_SCOPED = "refused_write_scoped"
REASON_WRITE_UNSCOPED = "refused_write_unscoped"
REASON_NEEDS_BODY = "refused_needs_body"
REASON_NOT_ALLOWLISTED = "refused_not_allowlisted"
REASON_NOT_PARSED = "refused_not_parsed"


ELN_ROUTES: tuple = (
    # -- reads: a search is a POST whose body is the query: it writes nothing.
    Route("POST", "/api/v1/search/*", writes=False,
          reason="search: the POST body is the query, nothing is written",
          source="@3657b15"),
    # -- a read added for the ELN specialist (e1). It is the collection
    # listing the search path cannot replace: a collection's NAME has to be
    # resolved to an id before a search can be scoped to it, and the server
    # answers that only here (`GET /collections` -> own + shared_with_me).
    # A GET, so there is no body whose scope could be misread, and the
    # answer is scoped by the server to what this account may read
    # (collection_api.rb: `serialized_own_collections_for(current_user)`
    # plus the shared list).
    Route("GET", "/api/v1/collections", writes=False,
          reason="lists the collections this account may read; a read",
          source="@main 2026-11-24"),
    # -- writes: everything below changes state. All refused at this gate.
    Route("POST", "/api/v1/authentication/token", writes=True,
          reason="mints an API token for the account (and the token has no "
                 "scopes, see SERVER_FACTS)",
          source="exporter"),
    Route("POST", "/api/v1/samples", writes=True,
          reason="creates a sample under a collection (collection_id in the "
                 "body)",
          source="exporter"),
    Route("POST", "/api/v1/research_plans", writes=True,
          reason="creates a research plan under a collection (collection_id "
                 "in the body)",
          source="exporter"),
    # The reaction create the ELN specialist drafts through. It was missing
    # from this table, which is why a draft to it was refused as
    # `refused_not_allowlisted` rather than as a write -- and a refusal for
    # the wrong reason is how a route that is never reached looks the same
    # as a rule that holds. Established from reaction_api.rb (`resource
    # :reactions` + `post do`), read 2026-11-24.
    Route("POST", "/api/v1/reactions", writes=True,
          reason="creates a reaction under a collection (collection_id in the "
                 "body; materials are resolved by UpdateMaterials)",
          source="@main 2026-11-24"),
    Route("POST", "/api/v1/attachments/upload_chunk", writes=True,
          reason="writes bytes into a pending upload (its target is the key "
                 "in the body)",
          source="exporter"),
    Route("POST", "/api/v1/attachments/upload_chunk_complete", writes=True,
          reason="closes a pending upload (checksum must match the chunks)",
          source="exporter"),
    Route("POST", "/api/v1/attachments", writes=True,
          reason="maps uploaded keys onto an element (the element is named "
                 "in the body)",
          source="exporter"),
    Route("DELETE", "/api/v1/ui_state", writes=True, needs_body=True,
          scope_fields=("currentCollection.id",),
          reason="DELETE ui_state without currentCollection.id in the body "
                 "targets everything the user owns; with it, only that "
                 "collection's state",
          source="@3657b15"),
    Route("POST", "/api/v1/versions/revert", writes=True,
          reason="replaces an element with an earlier revision; elements are "
                 "acts_as_paranoid and there is no undelete via the API",
          source="@3657b15"),
    Route("POST", "/api/v1/generic_elements/delete_revision", writes=True,
          reason="deletes a revision; acts_as_paranoid means the delete "
                 "stands (no undelete via the API)",
          source="@3657b15"),
    Route("POST", "/api/v1/molecules/delete_name", writes=True,
          reason="deletes a name of a molecule; there is no undo in the API",
          source="@3657b15"),
)


def match_route(method: str, path: str) -> Optional[Route]:
    """Look one (verb, path) up. The matching entry or None.

    Exact match on verb and the entry's route string; the trailing ``*``
    segment of ``/api/v1/search/*`` is the only wildcard and it matches
    exactly one more non-empty segment. Unknown verbs and unparseable paths
    simply match nothing (the caller is :func:`decide_route`, which refuses
    those before this runs and says why).
    """
    if not isinstance(method, str) or not isinstance(path, str):
        return None
    verb = method.strip().upper()
    if verb not in VERBS:
        return None
    for row in ELN_ROUTES:
        if row.verb == verb and _path_matches(row.path, path):
            return row
    return None


def decide_route(method: str, path: str,
                 body: Optional[object] = None) -> Decision:
    """Decide whether one request may proceed without a human decision.

    Input: the request's verb (str), its path without the query string (str),
    and its body -- a dict (already parsed), a JSON object as str/bytes, or
    None when the request has no body. Output: a :class:`Decision`. This is
    a pure function: it sends nothing and builds nothing.

    Order of the decision (each step refuses on its own terms):

    1. Parse (REASON_NOT_PARSED): unknown/absent verb (anything outside
       :data:`VERBS`), a path that is not a str, is empty or carries spaces,
       has a query string in it, does not sit under ``/api/v1/``, has a
       doubled slash or a ``..`` segment, or a body that was supplied but is
       not a JSON object.
    2. Look up (REASON_NOT_ALLOWLISTED): a (verb, path) no entry covers is
       refused.
    3. Refuse writes (REASON_WRITE / REASON_WRITE_SCOPED /
       REASON_WRITE_UNSCOPED / REASON_NEEDS_BODY) per the module docstring;
       an entry with ``needs_body=True`` and no body is refused as
       REASON_NEEDS_BODY, and its scope keys are read when a body exists.
    4. Otherwise read (``allowed=True``, REASON_READ).

    ``writes``/``allowed`` in the result mirror the entry's classification; a
    prohibited request is never reported as a read.
    """
    reason, route = _parse(method, path, body)
    if route is None:
        return Decision(allowed=False, verb=str(method), path=str(path),
                        writes=True, reason=reason)
    return _decide_for(route, path, body)


def _parse(method: object, path: object, body: object) -> tuple:
    """Returns (reason, route): a non-empty reason is why lookup never ran."""
    if not isinstance(method, str) or not isinstance(path, str):
        return (f"{REASON_NOT_PARSED}: verb and path must be strings", None)
    verb = method.strip().upper()
    if verb not in VERBS:
        return (f"{REASON_NOT_PARSED}: unknown verb {method.strip()!r}", None)
    if not _path_ok(path):
        return (f"{REASON_NOT_PARSED}: path this table cannot hold", None)
    parsed, ok = _as_body(body)
    if not ok:
        return (f"{REASON_NOT_PARSED}: body is not a JSON object", None)
    row = match_route(verb, path)
    if row is None:
        return (f"{REASON_NOT_ALLOWLISTED}: {verb} {path}", None)
    return ("", row)


def _decide_for(route: Route, path: str, body: object) -> Decision:
    parsed, _ = _as_body(body)
    scope = _present_scope(parsed, route) if isinstance(parsed, dict) else {}
    missing = [k for k in route.scope_fields if k not in scope]
    context = {
        "verb": route.verb,
        "path": path,
        "writes": route.writes,
        "route": route,
        "scope": scope,
    }

    if not route.writes:
        return Decision(allowed=True,
                        reason=f"{REASON_READ}: {path}",
                        **context)

    if route.needs_body and parsed is None:
        return Decision(allowed=False, body_consulted=False,
                        reason=f"{REASON_NEEDS_BODY}: {path} needs "
                               f"{','.join(route.scope_fields)} in the body",
                        **context)

    if route.scope_fields:
        body_consulted = isinstance(parsed, dict)
        if missing:
            return Decision(allowed=False, body_consulted=body_consulted,
                            reason=f"{REASON_WRITE_UNSCOPED}: {path}: "
                                   f"no {','.join(missing)} in the body -- "
                                   f"{route.reason}",
                            **context)
        return Decision(allowed=False, body_consulted=body_consulted,
                        reason=f"{REASON_WRITE_SCOPED}: {path}: scope="
                               f"{','.join(route.scope_fields)}: "
                               f"{route.reason}",
                        **context)

    return Decision(allowed=False,
                    reason=f"{REASON_WRITE}: {path}: {route.reason}",
                    **context)


def _path_ok(path: str) -> bool:
    """A path this table can hold: plain, under /api/v1/, no smuggled bits."""
    if (not path or path != path.strip() or "?" in path or "//" in path
            or not path.startswith("/api/v1/")):
        return False
    return ".." not in path.split("/")


def _as_body(body: object) -> tuple:
    """(parsed_body_or_None, ok). A supplied body must be a JSON object."""
    if body is None:
        return None, True
    if isinstance(body, dict):
        return body, True
    if isinstance(body, (str, bytes)):
        try:
            loaded = json.loads(body)
        except ValueError:
            return None, False
        return (loaded, True) if isinstance(loaded, dict) else (None, False)
    return None, False


def _present_scope(body: Optional[dict], route: Route) -> dict:
    found: dict = {}
    if not isinstance(body, dict):
        return found
    for key in route.scope_fields:
        value = _dig(body, key)
        if value is not None and value != "":
            found[key] = value
    return found


def _dig(body: dict, key: str):
    """``a.b`` reads through nested dicts; absent at any level -> None."""
    node = body
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


def _path_matches(entry_path: str, request_path: str) -> bool:
    """Strict segment match, exact strings only.

    Segments split on "/" with the empties preserved, so a trailing or
    leading slash makes a different path that matches nothing. A trailing
    ``*`` (the search family) consumes exactly one non-empty segment.
    """
    entry_segs = entry_path.split("/")
    request_segs = request_path.split("/")
    if entry_segs[-1] == "*":
        if len(entry_segs) != len(request_segs) or request_segs[-1] == "":
            return False
        entry_segs, request_segs = entry_segs[:-1], request_segs[:-1]
    if len(entry_segs) != len(request_segs):
        return False
    return all(a == b for a, b in zip(entry_segs, request_segs))