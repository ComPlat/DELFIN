"""A fake Chemotion ELN server, in process, answering the real shapes.

Where the shapes come from. NOT from a running server: every body below is
transcribed from the API's own source at ``ComPlat/chemotion_ELN@main``,
read on 2026-11-24, and each is annotated with the file it came from. A
shape invented here would make the client pass a test and fail on the real
thing, which is the defect this file exists to prevent.

The routes it serves, and the source of each answer:

=======================================  ===================================
``POST /api/v1/search/<plural>``         ``search_api.rb``:
                                         ``{<plural>: [...],
                                         total_elements: N}``, honouring
                                         ``page``/``per_page``
``GET /api/v1/collections``              ``collection_api.rb``:
                                         ``{own: [...], shared_with_me:
                                         [...]}``
``GET /api/v1/samples/<id>``             ``sample_api.rb`` +
                                         ``entities/sample_entity.rb``:
                                         ``{sample: {...}}`` with
                                         ``container.children`` holding the
                                         analyses (``container_entity.rb``)
``GET /api/v1/reactions/<id>``           ``reaction_api.rb`` +
                                         ``entities/reaction_entity.rb``:
                                         ``{reaction: {...}}``
``POST /api/v1/reactions``               ``reaction_api.rb``: creates under
                                         ``collection_id`` and answers
                                         ``{reaction: {...}}``
``POST /api/v1/samples``                 ``sample_api.rb``: answers
                                         ``{sample: {...}}``
``POST /api/v1/authentication/token``    ``authentication_api.rb``:
                                         ``{token: "chemtoken_..."}``
=======================================  ===================================

Errors are the server's own forms: ``{"error": "403 Forbidden"}`` with
status 403 (the ``error!('403 Forbidden', 403)`` calls throughout
``sample_api.rb``), ``401 Unauthorized`` for a missing or wrong bearer
token (``api.rb``'s ``detect_current_user_from_api_token`` only accepts a
token that starts with ``chemtoken``).

What it is NOT. It is not a Chemotion instance and it does not claim to be:
``SHAPES_FROM`` states exactly which source file each answer comes from, and
:data:`NOT_MODELLED` lists what this fake does not reproduce (permission
detail levels, the paranoid-delete semantics, the analysis `dataset`
attachment wiring). A test that needs one of those is a test of the real
stack, not of this client.

The transport is a FUNCTION, not a socket: ``send(method, url, body,
headers, timeout) -> (status, headers, bytes)``, the contract
``delfin.lab.eln_api._urlopen_send`` implements. The tests inject
``FakeELN().send``; nothing here binds a port, so the suite cannot be
blocked by a firewall and cannot be aimed at anything real.
"""

from __future__ import annotations

import json
from urllib.parse import parse_qs, urlsplit

#: Which source file each served shape was transcribed from.
SHAPES_FROM = {
    "search": "app/api/chemotion/search_api.rb @main 2026-11-24",
    "collections": "app/api/chemotion/collection_api.rb @main 2026-11-24",
    "samples": "app/api/chemotion/sample_api.rb + app/api/entities/"
               "sample_entity.rb @main 2026-11-24",
    "reactions": "app/api/chemotion/reaction_api.rb + app/api/entities/"
                 "reaction_entity.rb @main 2026-11-24",
    "container": "app/api/entities/container_entity.rb @main 2026-11-24",
    "token": "app/api/chemotion/authentication_api.rb @main 2026-11-24",
    "auth": "app/api/api.rb detect_current_user_from_api_token @main "
            "2026-11-24",
}

#: What this fake deliberately does not reproduce. A client behaviour that
#: depends on one of these is verified against a real instance, not here.
NOT_MODELLED = (
    "collection share detail levels (an element below the detail level is "
    "served anonymised -- ElementDetailLevelCalculator)",
    "acts_as_paranoid: a delete is not undone by this fake because it serves "
    "no delete at all",
    "the analysis `dataset` container's attachment wiring",
    "reaction materials resolution (Usecases::Reactions::UpdateMaterials): "
    "the fake stores what was posted",
)

#: The bearer prefix ``api.rb`` requires (``current_token.start_with?
#: ('chemtoken')``).
TOKEN_PREFIX = "chemtoken"


def _analysis(kind: str, name: str) -> dict:
    """One analysis, as ``ContainerEntity`` renders a child of a container.

    ``container_type: 'analysis'`` and ``extended_metadata.kind`` are the
    two fields :func:`delfin.lab.eln_tools._analysis_names` reads; both are
    from ``container_entity.rb`` (``expose :children``, and the
    ``extended_metadata[:kind]`` mapping).
    """
    return {
        "id": abs(hash((kind, name))) % 10_000,
        "name": name,
        "container_type": "analysis",
        "description": "",
        "extended_metadata": {"kind": kind, "status": "done", "report": False},
        "children": [],
    }


def _dataset() -> dict:
    """The ``dataset`` child every analysis container carries."""
    return {"id": 0, "name": "dataset", "container_type": "dataset",
            "description": "", "extended_metadata": {},
            "attachments": [], "children": []}


class FakeELN:
    """An in-process ELN with two collections and five samples.

    ``calls`` records every request the client made, so a test can assert
    what was SENT (verb, path, body, and whether a bearer token rode
    along) rather than only what came back. That is what makes the
    approval and credential tests measurable.
    """

    def __init__(self, *, token: str = "chemtoken_test123") -> None:
        self.token = token
        self.calls: list[dict] = []
        self.refuse_with: dict | None = None
        self.own = [
            {"id": 3, "label": "Synthesis 2026", "permission_level": 10,
             "is_shared": False, "shared_by_id": None},
            {"id": 4, "label": "Analysis", "permission_level": 10,
             "is_shared": False, "shared_by_id": None},
        ]
        self.shared: list = []
        # Samples: two WITH an analysis, two WITHOUT, one in the other
        # collection. ``has_analysis`` is what the chemist asks about, so
        # both sides of it exist here.
        self.samples = [
            self._sample(11, "S-001", 3, [_analysis("nmr", "1H NMR")]),
            self._sample(12, "S-002", 3, [_analysis("ir", "IR")]),
            self._sample(13, "S-003", 3, []),
            self._sample(14, "S-004", 3, []),
            self._sample(15, "S-005", 4, [_analysis("ms", "MS")]),
        ]
        self.reactions = [
            self._reaction(21, "Oxidation of 4a", 3, "done"),
            self._reaction(22, "Boc deprotection", 3, "running"),
        ]
        self.next_id = 100
        # Attachment state: the chunked upload's own ledger, so a test can
        # assert the exporter's sequence arrived whole rather than that a
        # 200 came back.
        self.upload_chunks: dict = {}
        self.uploaded: list = []
        self.attached: list = []

    # -- fixtures ---------------------------------------------------------

    def _sample(self, sid: int, name: str, collection: int,
                analyses: list) -> dict:
        children = list(analyses)
        if children:
            children.append(_dataset())
        return {
            "id": sid, "name": name, "short_label": f"{name}-lab",
            "external_label": "", "type": "sample",
            "sum_formula": "C10H12O2", "is_restricted": False,
            "collection_id": collection,
            "container": {"id": 900 + sid, "name": "root",
                          "container_type": "root",
                          "description": "",
                          "extended_metadata": {},
                          "children": children},
        }

    def _reaction(self, rid: int, name: str, collection: int,
                  status: str) -> dict:
        return {
            "id": rid, "name": name, "short_label": f"R{rid}",
            "type": "reaction", "status": status,
            "collection_id": collection,
            "reactants": [{"id": 1, "name": "4a", "short_label": "4a"}],
            "products": [{"id": 2, "name": "5a", "short_label": "5a"}],
            "solvents": [],
            "starting_materials": [],
        }

    # -- the transport ----------------------------------------------------

    def send(self, method, url, body, headers, timeout) -> tuple:
        """``(status, headers, bytes)`` -- the contract ``eln_api`` calls."""
        parts = urlsplit(url)
        path = parts.path
        query = parse_qs(parts.query)
        payload = {}
        if body:
            try:
                payload = json.loads(body.decode("utf-8"))
            except Exception:                           # noqa: BLE001
                payload = {}
        auth = (headers or {}).get("Authorization", "")
        self.calls.append({"method": method, "path": path,
                           "query": {k: v[0] for k, v in query.items()},
                           "body": payload, "auth": auth,
                           "url": url, "headers": dict(headers or {})})
        if self.refuse_with is not None:
            return (self.refuse_with["status"],
                    {"content-type": "application/json"},
                    json.dumps({"error": self.refuse_with["error"]}).encode())

        if path == "/api/v1/authentication/token" and method == "POST":
            return self._json(200, {"token": self.token})

        # Everything below needs the bearer the server's own helper demands.
        if not auth.startswith(f"Bearer {TOKEN_PREFIX}"):
            return self._json(401, {"error": "401 Unauthorized"})

        if path == "/api/v1/collections" and method == "GET":
            return self._json(200, {"own": self.own,
                                    "shared_with_me": self.shared})

        if path.startswith("/api/v1/search/") and method == "POST":
            return self._search(path.rsplit("/", 1)[-1], payload)

        if path.startswith("/api/v1/samples/") and method == "GET":
            return self._one("samples", self.samples, path)

        if path.startswith("/api/v1/reactions/") and method == "GET":
            return self._one("reactions", self.reactions, path)

        if path == "/api/v1/reactions" and method == "POST":
            collection = payload.get("collection_id")
            if collection is not None and not self._writable(collection):
                return self._json(403, {"error": "403 Forbidden"})
            self.next_id += 1
            row = {"id": self.next_id, "type": "reaction",
                   "name": payload.get("name"),
                   "short_label": f"R{self.next_id}",
                   "status": "planned", "collection_id": collection,
                   "reactants": [], "products": [], "solvents": [],
                   "starting_materials": []}
            self.reactions.append(row)
            return self._json(200, {"reaction": row})

        if path == "/api/v1/samples" and method == "POST":
            collection = payload.get("collection_id")
            if collection is not None and not self._writable(collection):
                return self._json(403, {"error": "403 Forbidden"})
            self.next_id += 1
            row = self._sample(self.next_id, str(payload.get("name")), 0, [])
            row["collection_id"] = collection
            self.samples.append(row)
            return self._json(200, {"sample": row})

        if path == "/api/v1/attachments/upload_chunk" and method == "POST":
            # chemotion.build_import's chunk step. The server keeps the key
            # so the complete step can check one arrived.
            self.upload_chunks.setdefault(str(payload.get("key")), []).append(
                payload.get("counter"))
            return self._json(200, {"ok": True})

        if path == "/api/v1/attachments/upload_chunk_complete" \
                and method == "POST":
            key = str(payload.get("key"))
            if key not in self.upload_chunks:
                return self._json(422, {"error": "unknown upload key"})
            self.uploaded.append({"key": key,
                                  "filename": payload.get("filename"),
                                  "checksum": payload.get("checksum"),
                                  "chunks": len(self.upload_chunks[key])})
            return self._json(200, {"ok": True, "key": key})

        if path == "/api/v1/attachments" and method == "POST":
            self.attached.append({"element_type": payload.get("element_type"),
                                  "attachment_ids": payload.get(
                                      "attachment_ids")})
            return self._json(200, {"ok": True})

        return self._json(404, {"error": "404 Not Found"})

    # -- helpers ----------------------------------------------------------

    def _writable(self, collection_id) -> bool:
        """``writable_collection_for``: own collections with full access."""
        return any(str(c["id"]) == str(collection_id) and
                   int(c.get("permission_level", 0)) >= 10 for c in self.own)

    def _one(self, element: str, rows: list, path: str) -> tuple:
        try:
            wanted = int(path.rstrip("/").rsplit("/", 1)[-1])
        except ValueError:
            return self._json(404, {"error": "404 Not Found"})
        for row in rows:
            if row["id"] == wanted:
                return self._json(200, {element[:-1]: row})
        # sample_api.rb / reaction_api.rb: `Sample.find(id)` raises and
        # collection_api.rb's rescue_from answers the text below.
        return self._json(404, {"error": "Couldn't find Sample"})

    def _search(self, element: str, body: dict) -> tuple:
        """``search_api.rb``: paginate by ``page``/``per_page``, and say the
        total as ``total_elements``.

        The id selection is the ``ui_state`` shape ``search_api.rb``'s
        ``:search_params`` declare (``checkedIds``/``included_ids``, with
        ``uncheckedIds``/``excluded_ids`` removing), which is how a caller
        asks for ONE row without a per-id GET (none is established from the
        source).
        """
        source = {"samples": self.samples,
                  "reactions": self.reactions}.get(element, [])
        included = {str(i) for i in (body.get("included_ids")
                                     or body.get("checkedIds") or [])}
        excluded = {str(i) for i in (body.get("excluded_ids")
                                     or body.get("uncheckedIds") or [])}
        if included:
            source = [r for r in source if str(r.get("id")) in included]
        if excluded:
            source = [r for r in source if str(r.get("id")) not in excluded]
        collection = body.get("collection_id")
        if collection is not None:
            source = [r for r in source
                      if str(r.get("collection_id")) == str(collection)]
        name_filter = ""
        conditions = body.get("search_by_method") or ""
        if conditions == "substring" and body.get("value"):
            name_filter = str(body["value"]).lower()
        if name_filter:
            source = [r for r in source
                      if name_filter in str(r.get("name", "")).lower()]
        per_page = max(1, int(body.get("per_page") or 7))
        page = max(1, int(body.get("page") or 1))
        total = len(source)
        start = (page - 1) * per_page
        rows = source[start:start + per_page]
        headers = {"content-type": "application/json",
                   "x-total": str(total),
                   "x-total-pages": str(max(1, -(-total // per_page))),
                   "x-page": str(page), "x-per-page": str(per_page)}
        return (200, headers,
                json.dumps({element: rows,
                            "total_elements": total}).encode())

    @staticmethod
    def _json(status: int, obj) -> tuple:
        return (status, {"content-type": "application/json"},
                json.dumps(obj).encode("utf-8"))


#: The local URL a test points a client at. Loopback, so
#: ``eln_api.is_loopback`` answers True and the proxy is skipped -- the
#: same path a Docker instance on this host takes. Nothing binds it: the
#: transport is the in-process function above.
FAKE_URL = "http://127.0.0.1:8080"
