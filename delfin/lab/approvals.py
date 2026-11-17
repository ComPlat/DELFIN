"""Human approval of one hazardous action.

Input: a request the kernel files when a hazardous write arrives without an
approval; a grant the person makes with ``delfin-agent lab approve <id>``.
Output: a single-use approval bound to exactly that action.

Binding. The request's digest is SHA-256 over the installation, the
manifest checksum, device, channel, the value in the channel's unit, the
person and the session. A grant is accepted only for a request with the
same digest, so an approval of 60 degC cannot carry a 140 degC write, and
one session's approval cannot carry another session's write.

Who may grant. Not a process started by an agent session (the
``DELFIN_AGENT_SESSION`` marker), and not another person: the grant names
the approver, and the approver must be the person the request names. The
grant's MAC is HMAC-SHA256 under ``trust_store.derived_key("lab_approval")``,
which returns no key inside an agent session, so a grant written into the
state file by an agent's shell does not verify. In people mode each person
has their own trust store, so another person's grant does not verify either.

Single use and expiry. A consumed grant is refused as a replay; a request
expires ``TTL_S`` after it was filed. An e-stop cancels pending requests.

Rejected: an approval that covers a device for a period (the model could
use it for any value in range); approvals answered through the agent's own
confirmation dialog (that dialog follows the permission mode, and bypass
answers it without asking).
"""

from __future__ import annotations

import hashlib
import json
import secrets
from typing import Any, Optional

from . import audit as _audit
from . import store as _store

TTL_S = 300.0
_PURPOSE = "lab_approval"


class ApprovalError(PermissionError):
    pass


def digest(*, installation: str, checksum: str, device: str, channel: str,
           value: Any, person: str, session: str) -> str:
    v = round(float(value), 9) if isinstance(value, (int, float)) \
        and not isinstance(value, bool) else str(value)
    blob = json.dumps([installation, checksum, device, channel, v, person,
                       session], separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _payload(rec: dict) -> bytes:
    return json.dumps([rec["id"], rec["digest"], rec["expires"],
                       rec["approver"]], separators=(",", ":")).encode()


def file_request(doc: dict, *, digest_: str, person: str, session: str,
                 device: str, channel: str, value: Any, unit: Optional[str],
                 now: float) -> dict:
    """File (or return the pending) request for this exact action."""
    for rec in doc["approvals"].values():
        if rec["digest"] == digest_ and rec["status"] in ("pending", "granted") \
                and rec["expires"] > now:
            return rec
    rid = "lab-" + secrets.token_hex(4)
    rec = {"id": rid, "digest": digest_, "person": person, "session": session,
           "device": device, "channel": channel, "value": value, "unit": unit,
           "created": now, "expires": now + TTL_S, "status": "pending",
           "approver": "", "mac": ""}
    doc["approvals"][rid] = rec
    return rec


def grant(request_id: str, *, approver: Optional[str] = None,
          now: Optional[float] = None) -> dict:
    """The person's grant. Refused inside an agent session."""
    from delfin.agent.approval_guard import inside_agent_session
    from delfin.agent import trust_store
    if inside_agent_session():
        _audit.alert("lab_approval_refused",
                     f"{request_id}: grant from inside an agent session")
        raise ApprovalError("only a person approves a lab action; this "
                            "process was started by an agent session")
    approver = approver or _store.person()
    now = _store.now() if now is None else now
    with _store.state() as doc:
        rec = doc["approvals"].get(request_id)
        if rec is None:
            raise ApprovalError(f"no approval request {request_id}")
        if rec["person"] != approver:
            _audit.alert("lab_approval_refused",
                         f"{request_id}: {approver} tried to approve an "
                         f"action of {rec['person']}")
            raise ApprovalError(f"{request_id} is {rec['person']}'s request; "
                                "only that person approves it")
        if rec["status"] != "pending":
            raise ApprovalError(f"{request_id} is {rec['status']}")
        if rec["expires"] <= now:
            rec["status"] = "expired"
            raise ApprovalError(f"{request_id} has expired")
        rec["approver"] = approver
        tag = trust_store.mac(_PURPOSE, _payload(rec))
        if not tag:
            raise ApprovalError("no signing key available; nothing approved")
        rec["mac"] = tag
        rec["status"] = "granted"
        rec["granted"] = now
        out = dict(rec)
    _audit.record("approval_granted", id=request_id, approver=approver,
                  device=out["device"], channel=out["channel"],
                  value=out["value"], unit=out["unit"],
                  session=out["session"])
    return out


def consume(doc: dict, request_id: str, *, digest_: str, person: str,
            now: float) -> Optional[str]:
    """Check and consume a grant inside the kernel's transaction.

    Returns None when the grant is valid (and marks it consumed), else the
    refusal reason."""
    from delfin.agent import trust_store
    rec = doc["approvals"].get(str(request_id or ""))
    if rec is None:
        return (f"approval {request_id!r} does not exist; a hazardous action "
                "is approved only by a person through `delfin-agent lab "
                "approve <id>`")
    if rec["status"] == "consumed":
        return f"approval {request_id} was already used (replay refused)"
    if rec["status"] != "granted":
        return f"approval {request_id} is {rec['status']}, not granted"
    if rec["expires"] <= now:
        rec["status"] = "expired"
        return f"approval {request_id} has expired"
    if rec["digest"] != digest_:
        return (f"approval {request_id} was given for a different action "
                f"({rec['device']}.{rec['channel']} = {rec['value']} "
                f"{rec['unit'] or ''}, session {rec['session']})")
    if rec["person"] != person or rec["approver"] != person:
        return f"approval {request_id} belongs to {rec['person']}"
    ok = trust_store.check_mac(_PURPOSE, _payload(rec), rec.get("mac"))
    if ok is not True:
        return (f"approval {request_id} does not carry a valid signature "
                "(not granted by a person)")
    rec["status"] = "consumed"
    rec["consumed"] = now
    return None


def pending() -> list[dict]:
    doc = _store.snapshot()
    now = _store.now()
    return [r for r in doc["approvals"].values()
            if r["status"] == "pending" and r["expires"] > now]
