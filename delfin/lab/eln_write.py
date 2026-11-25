"""A write to the ELN, and the person's approval that has to precede it.

Input: what the caller wants written (a reaction or a sample, as a DRAFT)
and, on the second call, the approval id a person granted. Output: the
created element, or a refusal that says exactly what would have been
written and what is now missing.

WHY THIS IS NOT A PERMISSION-MODE QUESTION. DELFIN's confirmation dialog
follows the permission mode, and ``bypassPermissions`` answers it without
asking (``delfin/lab/approvals.py`` records the same reasoning for a
device write). A permission mode is a session posture a person chose once;
the ELN write is a different act every time, with a different body, and
an element the API cannot undelete (``acts_as_paranoid``, no undelete
route -- ``eln_routes.SERVER_FACTS``). So the mode is never consulted
here. ``approve`` reads no permission mode and cannot be satisfied by
one.

WHAT THE APPROVAL IS BOUND TO. The SHA-256 of the exact request: base
URL, verb, path, and the JSON body as it would be sent, with sorted keys
so an equivalent dict is the same request. A grant for one draft cannot
carry another -- a changed name, a changed collection, a changed material
list, a different server. The grant is signed with the account key
derived for the purpose :data:`PURPOSE` (``trust_store.mac``), which is
EMPTY inside an agent session, so a grant written by an agent's own shell
does not verify. It is single use and expires (:data:`TTL_S`).

WHAT A REFUSAL SAYS. The draft in full -- element, name, collection, the
number of materials -- and the exact command a person runs to approve it,
including the id and the digest's first 12 characters so a person can see
that the thing they approve is the thing that was shown. It never says
"permission denied" and stops.

WHAT IS NOT WRITTEN. A body :data:`delfin.lab.eln_api.UNVERIFIED` names
is refused rather than guessed at: a create whose required attributes
were not established from the API source would be a plausible-looking
element on a real server, which is exactly the failure this repository
refuses.
"""

from __future__ import annotations

import hashlib
import json
import secrets
from typing import Any, Callable, Optional

from . import eln_api

#: The key purpose the grant is signed under. One purpose, separate from
#: ``lab_approval``, so a grant for one cannot verify as the other.
PURPOSE = "eln_write"

#: How long a filed request stays approvable, seconds. Shorter than the
#: lab's 300 s: an ELN draft is composed in a conversation, and an
#: approval that outlives the conversation is an approval of something
#: the person no longer remembers seeing.
TTL_S = 240.0

#: The elements a DRAFT can be made of. A list, not a dict of shapes: the
#: shape of each is decided by :func:`draft_body`, and adding one there
#: without adding it here is what the test refuses.
DRAFT_ELEMENTS = ("reaction", "sample")

#: Where each draft element is created. Read from the API's mount and
#: resource (:data:`eln_api.SOURCES`), not guessed from the name.
_CREATE_PATH = {
    "reaction": "/api/v1/reactions",
    "sample": "/api/v1/samples",
}


class WriteRefused(RuntimeError):
    """A write that did not happen, with the draft that would have been."""

    def __init__(self, message: str, *, code: str = "",
                 draft: Optional[dict] = None,
                 request_id: str = "", digest: str = "") -> None:
        super().__init__(message)
        self.code = code
        self.draft = draft or {}
        self.request_id = request_id
        self.digest = digest


def draft_body(element: str, fields: dict) -> dict:
    """The exact JSON body a draft create would send.

    Only the fields established from source, plus the caller's own. The
    envelope is ``name`` + ``collection_id`` (both read from
    ``sample_api.rb``/``reaction_api.rb``), ``body`` for a research-plan
    style element, and ``materials`` for a reaction (the route reads
    ``materials`` and hands it to ``Usecases::Reactions::UpdateMaterials``).
    A key this function does not know is passed through UNCHANGED rather
    than dropped: a caller that learned a field from the API source must
    not have it silently removed here.
    """
    kind = str(element or "").strip().lower()
    if kind not in DRAFT_ELEMENTS:
        raise WriteRefused(
            f"cannot draft a {kind or 'nameless'} element: this writer "
            f"creates {', '.join(DRAFT_ELEMENTS)} as drafts",
            code="refused_unknown_element")
    data = dict(fields or {})
    if not str(data.get("name") or "").strip():
        raise WriteRefused(
            "a draft needs a name: it is what a person sees in the "
            "collection before anything else about it",
            code="refused_unnamed")
    out: dict = {"name": str(data.pop("name")).strip()}
    collection = data.pop("collection_id", None)
    if collection is not None:
        out["collection_id"] = collection
    if kind == "reaction":
        materials = data.pop("materials", None)
        if isinstance(materials, dict) and materials:
            out["materials"] = materials
        # A reaction with no materials is still a legal draft on the
        # server (its materials are an empty list); it is NOT treated as
        # an error, because a person may want the empty vessel first.
    out.update(data)
    return out


def digest_for(base_url: str, element: str, body: dict) -> str:
    """SHA-256 over the exact request a create would make.

    The body is serialised with sorted keys, so two callers that built
    the same request in a different key order agree; anything that
    changes the body changes the digest, which is what makes a grant
    single-purpose without a second table.
    """
    kind = str(element or "").strip().lower()
    path = _CREATE_PATH.get(kind, f"/api/v1/{kind}s")
    blob = json.dumps([eln_api.display_url(base_url), "POST", path, body],
                      sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def describe(element: str, body: dict) -> str:
    """What will be written, in one sentence a person can check."""
    kind = str(element or "").strip().lower()
    name = str(body.get("name") or "").strip() or "(unnamed)"
    bits = [f"{kind} draft {name!r}"]
    if body.get("collection_id") is not None:
        bits.append(f"in collection {body['collection_id']}")
    materials = body.get("materials")
    if isinstance(materials, dict) and materials:
        counts = ", ".join(f"{len(v) if isinstance(v, list) else 1} {k}"
                           for k, v in sorted(materials.items()))
        bits.append(f"with {counts}")
    return " ".join(bits)


# ---------------------------------------------------------------------------
# The approval ledger: file, grant, consume. The kernel's shape, ELN's keys
# ---------------------------------------------------------------------------

def _payload(rec: dict) -> bytes:
    return json.dumps([rec["id"], rec["digest"], rec["expires"],
                       rec["approver"]], separators=(",", ":")).encode()


def _now() -> float:
    from . import store
    return store.now()


def file_request(*, digest_: str, person: str, session: str,
                 element: str, body: dict, base_url: str) -> dict:
    """File (or return the still-pending) request for this exact write."""
    from . import store
    now = _now()
    with store.state() as doc:
        ledger = doc.setdefault("eln_writes", {})
        for rec in ledger.values():
            if (rec.get("digest") == digest_
                    and rec.get("status") in ("pending", "granted")
                    and rec.get("expires", 0) > now):
                return dict(rec)
        rid = "eln-" + secrets.token_hex(4)
        rec = {
            "id": rid, "digest": digest_, "person": person,
            "session": session, "element": str(element),
            "name": str(body.get("name") or ""),
            "collection_id": body.get("collection_id"),
            "base_url": eln_api.display_url(base_url),
            "created": now, "expires": now + TTL_S, "status": "pending",
            "approver": "", "mac": "",
        }
        ledger[rid] = rec
        return dict(rec)


def grant(request_id: str, *, approver: Optional[str] = None) -> dict:
    """The person's grant. Refused inside an agent session, and for a
    request filed by someone else.

    The same two gates as ``delfin.lab.approvals.grant``, for the same
    measured reason: an agent's own shell must not be able to sign its own
    write, and one person's grant must not release another person's.
    """
    from delfin.agent.approval_guard import inside_agent_session
    from delfin.agent import trust_store
    from . import store

    if inside_agent_session():
        raise WriteRefused(
            "only a person approves an ELN write; this process was started "
            "by an agent session", code="refused_agent_session")
    who = approver or store.person()
    now = _now()
    with store.state() as doc:
        ledger = doc.setdefault("eln_writes", {})
        rec = ledger.get(str(request_id or ""))
        if rec is None:
            raise WriteRefused(
                f"no ELN write request {request_id!r}", code="refused_no_such")
        if rec.get("person") != who:
            raise WriteRefused(
                f"{request_id} is {rec.get('person')}'s request; only that "
                "person approves it", code="refused_not_yours")
        if rec.get("status") != "pending":
            raise WriteRefused(
                f"{request_id} is {rec.get('status')}", code="refused_not_pending")
        if rec.get("expires", 0) <= now:
            rec["status"] = "expired"
            raise WriteRefused(f"{request_id} has expired",
                               code="refused_expired")
        rec["approver"] = who
        tag = trust_store.mac(PURPOSE, _payload(rec))
        if not tag:
            raise WriteRefused(
                "no signing key available here; nothing approved",
                code="refused_no_key")
        rec["mac"] = tag
        rec["status"] = "granted"
        rec["granted"] = now
        return dict(rec)


def consume(doc: dict, request_id: str, *, digest_: str, person: str,
            now: float) -> Optional[str]:
    """Check and consume one grant inside the caller's transaction.

    Returns None when the grant is valid (and marks it consumed), else the
    refusal reason. Mirrors ``delfin.lab.approvals.consume``: a replayed
    grant, a grant for another body, an expired one and an unsigned one
    each say which they are.
    """
    from delfin.agent import trust_store

    rec = doc.get("eln_writes", {}).get(str(request_id or ""))
    if rec is None:
        return (f"approval {request_id!r} does not exist; an ELN write is "
                "approved only by a person through `delfin-agent eln approve "
                "<id>`")
    if rec.get("status") == "consumed":
        return f"approval {request_id} was already used (replay refused)"
    if rec.get("status") != "granted":
        return f"approval {request_id} is {rec.get('status')}, not granted"
    if rec.get("expires", 0) <= now:
        rec["status"] = "expired"
        return f"approval {request_id} has expired"
    if rec.get("digest") != digest_:
        return (f"approval {request_id} was given for a different write "
                f"({rec.get('element')} {rec.get('name')!r} in collection "
                f"{rec.get('collection_id')}); the body has changed since")
    if rec.get("person") != person or rec.get("approver") != person:
        return f"approval {request_id} belongs to {rec.get('person')}"
    if trust_store.check_mac(PURPOSE, _payload(rec), rec.get("mac")) is not True:
        return (f"approval {request_id} does not carry a valid signature "
                "(not granted by a person)")
    rec["status"] = "consumed"
    rec["consumed"] = now
    return None


def pending() -> list[dict]:
    """The requests a person may approve now."""
    from . import store
    doc = store.snapshot()
    now = _now()
    return [dict(r) for r in doc.get("eln_writes", {}).values()
            if r.get("status") == "pending" and r.get("expires", 0) > now]


# ---------------------------------------------------------------------------
# The writer: the one route by which a draft reaches the server
# ---------------------------------------------------------------------------

def create_draft(client: "eln_api.Client", element: str, fields: dict,
                 *, approval: str = "", person: str = "",
                 session: str = "") -> dict:
    """Create one DRAFT, but only against a person's approval.

    Two calls, deliberately. Without ``approval`` the request is FILED and
    the refusal carries everything a person needs to decide (the draft in
    full, the id, the digest prefix, and the command). With a valid
    ``approval`` the same body is sent once, and the grant is consumed.

    The permission mode is not a parameter and is never read: a caller
    cannot pass one, and ``bypassPermissions`` reaches the same refusal as
    every other mode.
    """
    from . import store

    body = draft_body(element, fields)
    kind = str(element).strip().lower()
    digest_ = digest_for(client.base_url, kind, body)
    who = person or store.person()
    sid = session or ""

    if not approval:
        rec = file_request(digest_=digest_, person=who, session=sid,
                           element=kind, body=body, base_url=client.base_url)
        raise WriteRefused(
            f"needs a person's approval. It would write "
            f"{describe(kind, body)} to {eln_api.display_url(client.base_url)}"
            f" (request {rec['id']}, digest {digest_[:12]}). A person runs: "
            f"`delfin-agent eln approve {rec['id']}`, then this is repeated "
            f"with approval={rec['id']!r}. Nothing was sent.",
            code="needs_approval", draft={"element": kind, "body": body},
            request_id=rec["id"], digest=digest_)

    # Consumed in one short transaction and the send made OUTSIDE it: the
    # state lock is the lab's, held by every device operation on this
    # installation, and a create against a slow server inside it would
    # stall them all for the request's whole timeout. Consuming first keeps
    # the property that matters -- one grant, one send, even when the send
    # raises -- without holding a lock across the network.
    with store.state() as doc:
        why = consume(doc, approval, digest_=digest_, person=who, now=_now())
    if why:
        raise WriteRefused(
            f"{why}. Nothing was sent; the draft was "
            f"{describe(kind, body)}.",
            code="refused_approval", draft={"element": kind, "body": body},
            request_id=str(approval), digest=digest_)
    answer = client.send_approved("POST", _CREATE_PATH[kind],
                                  digest=digest_, body=body)
    return {"ok": True, "request_id": str(approval), "digest": digest_,
            "element": kind, "sent": body, "answer": answer.body,
            "status": answer.status}


def send_with_approval(client: "eln_api.Client", call: Callable[[], Any],
                       *, approval: str, digest_: str, person: str) -> Any:
    """Run ``call`` once, guarded by ``approval`` for ``digest_``.

    A narrower seam than :func:`create_draft`, for a write whose request is
    built elsewhere (an attachment upload of a run's export). The grant is
    consumed in its own transaction BEFORE ``call`` runs (never holding the
    state lock across the network), so a raising ``call`` does not leave a
    usable grant behind.
    """
    from . import store

    with store.state() as doc:
        why = consume(doc, approval, digest_=digest_, person=person, now=_now())
    if why:
        raise WriteRefused(why, code="refused_approval",
                           request_id=str(approval), digest=digest_)
    return call()
