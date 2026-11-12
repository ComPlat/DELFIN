"""Chemotion ELN exporter: turn a finished run record into ELN objects.

Input: a ``delfin.lab.run/1`` record (``delfin/lab/record.py``). Output:

* ``build_import()``  a structured model of what to write to the ELN: one
  research-plan style element carrying the run parameters as its body, one
  sample entry per run sample with its analyses, and the record / device
  log / measurements / parameters as attachments.
* ``request_sequence()`` the ordered HTTP requests an upload would send.
* ``dry_run()``  writes the exact request bodies to a folder, sending
  nothing. A person can review them before any upload.
* ``upload()``  POSTs the sequence to a Chemotion ELN through an injectable
  transport. It refuses unless a token (supplied by reference, never
  embedded) and an explicit confirmation are given.

The endpoint facts are the CURRENT official API, read from source this
session and cited in :data:`SOURCES`:

* Grape API under ``/api/v1`` (prefix ``api``, version ``v1``) mounting
  ``SampleAPI``, ``ResearchPlanAPI``, ``AttachmentAPI`` and
  ``AuthenticationAPI`` (chemotion_ELN ``app/api/api.rb``).
* ``POST /api/v1/authentication/token`` with ``username`` and ``password``
  returns ``{"token": "<chemtoken_...>"}``; the API-token auth branch
  requires the bearer token to start with ``chemtoken``
  (``authentication_api.rb``).
* Samples are created under ``/samples`` with a ``collection_id``, a nested
  ``container`` (their analyses live there) and ``segments``
  (``sample_api.rb``).
* Attachments upload in 5000-byte chunks: ``POST /attachments/upload_chunk``
  (multipart field ``file``, form ``key``, ``counter``) then
  ``POST /attachments/upload_chunk_complete`` (form ``key``, ``filename``,
  ``checksum`` = MD5 hex of the whole file) -- from the official Python
  client chemotionapi ``analyses_manager.py``.

What is NOT verified here (named in ``unverified`` and below so a deploy
confirms it against the ELN version): the exact create-request body of the
research-plan endpoint, how analyses nest into containers/datasets for a
fresh sample, token expiry, and the request that maps uploaded attachment
keys onto an element (modelled from the official client's
``Attachments.save``). Nothing here contacts a real server; that is only
``upload()``, and never from these tests.
"""

from __future__ import annotations

import hashlib
import json
import mimetypes
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

FORMAT = "delfin.lab.chemotion/1"

# Chunk size used by the official client (chemotionapi analyses_manager.py).
MAX_UPLOAD_CHUNK = 5000

# Endpoint facts read from source this session, with their primary URLs.
SOURCES = {
    "api.root": "ComPlat/chemotion_ELN app/api/api.rb "
                "(prefix api, version v1; mounts the *_API resources)",
    "authentication": "ComPlat/chemotion_ELN app/api/chemotion/authentication_api.rb",
    "samples": "ComPlat/chemotion_ELN app/api/chemotion/sample_api.rb",
    "research_plans": "ComPlat/chemotion_ELN app/api/chemotion/research_plan_api.rb",
    "attachments": "Chemotion/chemotionapi chemotion_api/elements/analyses_manager.py",
}

UNVERIFIED = [
    "exact request body of POST /research_plans (create params)",
    "how analyses nest into containers -> datasets for a fresh sample",
    "token expiry and the full /authentication mount path",
    "the request that maps uploaded attachment keys onto an element "
    "(modelled from the official client's Attachments.save)",
]


class ChemotionError(RuntimeError):
    """A refusal or failed upload; nothing leaves the process silently."""


@dataclass
class HttpRequest:
    """One exact HTTP request an upload would send (what dry_run writes).

    ``kind`` selects how the body is put on the wire:

    * ``"json"`` -- ``data`` is the JSON request body.
    * ``"multipart"`` -- ``files`` is a list of ``(field, filename, bytes)``
      parts; ``data`` holds the extra form fields.
    * ``"attach"`` -- JSON body that maps uploaded attachment keys onto an
      element (UNVERIFIED exact shape; see :data:`UNVERIFIED`).
    """

    method: str
    path: str
    data: Optional[dict] = None
    files: list = field(default_factory=list)
    kind: str = "json"
    auth: bool = True


def build_import(rec: dict) -> dict:
    """Map a ``delfin.lab.run/1`` record to the ELN import model.

    One research-plan element (body = the run parameters), one sample entry
    per run sample with its analyses, and the full record / device log /
    measurements / parameters as attachments.
    """
    samples = []
    by_sample: dict = {}
    for m in rec.get("measurements") or []:
        by_sample.setdefault(m.get("sample"), []).append(m)
    for sid, s in (rec.get("samples") or {}).items():
        samples.append({
            "name": sid,
            "external_label": s.get("container"),
            "location": f"{s.get('module')}:{s.get('position')}",
            "analyses": [
                {"kind": m.get("quantity"), "device": m.get("device"),
                 "value": m.get("value"), "unit": m.get("unit"),
                 "container": m.get("container")}
                for m in by_sample.get(sid, [])
            ],
        })
    parameters = {
        "run_id": rec.get("run_id"),
        "status": rec.get("status"),
        "plan": rec.get("plan"),
        "plan_digest": rec.get("plan_digest"),
        "manifest": rec.get("manifest"),
        "actor": rec.get("actor"),
        "steps": rec.get("steps"),
        "deviations": rec.get("deviations"),
        "fair": rec.get("fair"),
        "started": rec.get("started"),
        "ended": rec.get("ended"),
    }
    return {
        "format": FORMAT,
        "element": {
            "type": "research_plan",
            "name": f"Automated run {rec['run_id']}",
            "body": parameters,
        },
        "samples": samples,
        "attachments": [
            _attachment("record.json", json.dumps(rec, indent=1, sort_keys=True)),
            _attachment("device_log.json",
                        json.dumps(rec.get("device_log") or [], indent=1)),
            _attachment("measurements.json",
                        json.dumps(rec.get("measurements") or [], indent=1)),
            _attachment("parameters.json",
                        json.dumps(parameters, indent=1, sort_keys=True)),
        ],
        "unverified": list(UNVERIFIED),
    }


def _attachment(filename: str, text: str) -> dict:
    ctype = mimetypes.guess_type(filename)[0] or "application/octet-stream"
    return {"filename": filename, "content_type": ctype, "text": text}


def request_sequence(model: dict,
                     collection_id: int = 1) -> list[HttpRequest]:
    """Ordered HTTP requests that would write ``model`` to the ELN.

    Research-plan create, then one sample create per sample, then each
    attachment as a chunked upload (one request per 5000-byte chunk, plus a
    completion request), then the attach step mapping the upload keys onto
    the research-plan element.
    """
    seq = [HttpRequest("POST", "/api/v1/research_plans", kind="json", data={
        "name": model["element"]["name"],
        "body": model["element"]["body"],
        "collection_id": collection_id,
        "attachments": [a["filename"] for a in model["attachments"]],
    })]
    for s in model["samples"]:
        seq.append(HttpRequest("POST", "/api/v1/samples", kind="json", data={
            "name": s["name"],
            "external_label": s.get("external_label"),
            "collection_id": collection_id,
            "container": {"analyses": s.get("analyses") or []},
        }))
    keys = []
    for a in model["attachments"]:
        key = str(uuid.uuid4())
        keys.append(key)
        data = a["text"].encode("utf-8")
        digest = hashlib.md5(data).hexdigest()
        n = 0
        for i in range(0, len(data), MAX_UPLOAD_CHUNK):
            chunk = data[i:i + MAX_UPLOAD_CHUNK]
            seq.append(HttpRequest(
                "POST", "/api/v1/attachments/upload_chunk", kind="multipart",
                data={"key": key, "counter": str(n)},
                files=[("file", a["filename"], chunk)]))
            n += 1
        seq.append(HttpRequest(
            "POST", "/api/v1/attachments/upload_chunk_complete",
            kind="multipart", data={"key": key, "filename": a["filename"],
                                    "checksum": digest}))
    seq.append(HttpRequest(
        "POST", "/api/v1/attachments", kind="attach", data={
            "element_type": "ResearchPlan",
            "attachment_ids": keys,
        }))
    return seq


def dry_run(model: dict, out_dir: Path | str, collection_id: int = 1,
            ) -> list[Path]:
    """Write the exact request bodies to ``out_dir``; send nothing.

    Each request becomes a file ``<n>_<method>_<slug>.json`` holding the
    request (method, path, body or form fields and file parts) exactly as an
    upload would send it. Returns the paths written.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for i, req in enumerate(request_sequence(model, collection_id)):
        slug = req.path.rstrip("/").rsplit("/", 1)[-1] or "root"
        p = out / f"{i + 1:02d}_{req.method}_{slug}.json"
        p.write_text(json.dumps(_as_dict(req), indent=1, sort_keys=True),
                     encoding="utf-8")
        written.append(p)
    return sorted(written)


def _as_dict(req: HttpRequest) -> dict:
    d: dict = {"method": req.method, "path": req.path, "kind": req.kind}
    if req.data is not None:
        d["data"] = req.data
    d["files"] = [{"field": f, "filename": fn, "bytes(utf-8)": b.decode("utf-8")}
                  for f, fn, b in req.files]
    return d


def _multipart_body(req: HttpRequest, boundary: str) -> bytes:
    parts = []
    for f, fn, b in req.files:
        parts.append(
            f"--{boundary}\r\n"
            f'Content-Disposition: form-data; name="{f}"; filename="{fn}"\r\n'
            f"Content-Type: application/octet-stream\r\n\r\n".encode())
        parts.append(b)
        parts.append(b"\r\n")
    for k, v in (req.data or {}).items():
        parts.append(
            f"--{boundary}\r\n"
            f'Content-Disposition: form-data; name="{k}"\r\n\r\n'
            f"{v}\r\n".encode())
    parts.append(f"--{boundary}--\r\n".encode())
    return b"".join(parts)


def upload(model: dict, base_url: str, token: str, *,
           collection_id: int = 1, confirm: bool = False,
           transport: Optional[Callable] = None) -> list[dict]:
    """POST ``model`` to the ELN at ``base_url``.

    Refuses (``ChemotionError``) unless ``confirm`` is True and no request is
    sent. ``token`` is a value supplied by reference by the caller (the CLI
    reads it from an env var or a file, never embedded). ``transport`` is an
    injectable ``(req, base_url, token, boundary) -> response`` for tests; the
    default posts the requests over HTTP with ``urllib.request``.
    """
    if not isinstance(token, str) or not token:
        raise ChemotionError("upload needs a token, supplied by reference")
    if not confirm:
        raise ChemotionError("upload refused: explicit confirmation required")
    responses = []
    for req in request_sequence(model, collection_id):
        if transport is not None:
            responses.append(transport(req, base_url, token))
            continue
        boundary = f"----delfin{uuid.uuid4().hex}"
        if req.kind == "json":
            body = json.dumps(req.data).encode("utf-8")
            ctype = "application/json"
            headers = {"Content-Type": ctype}
        elif req.kind == "multipart":
            body = _multipart_body(req, boundary)
            headers = {"Content-Type": f"multipart/form-data; boundary={boundary}"}
        else:  # attach
            body = json.dumps(req.data).encode("utf-8")
            headers = {"Content-Type": "application/json"}
        if req.auth:
            headers["Authorization"] = f"Bearer {token}"
        responses.append(_urlopen(req.method, base_url + req.path, body,
                                  headers))
    return responses


def _urlopen(method: str, url: str, body: bytes, headers: dict) -> dict:
    import urllib.request
    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    try:
        # A loopback target (a test fake server) must not be routed through the
        # environment's HTTP proxy, which would intercept it. Remote ELN hosts
        # keep the default opener so an http_proxy is honoured.
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({})
        ) if _is_loopback(url) else None
        with (opener or urllib.request).open(req, timeout=30) as resp:
            status = resp.status
            payload = resp.read()
    except urllib.error.HTTPError as exc:            # noqa: N818
        raise ChemotionError(f"{method} {url}: HTTP {exc.code}") from exc
    except OSError as exc:
        raise ChemotionError(f"{method} {url}: {exc}") from exc
    return {"status": status, "body": payload.decode("utf-8", "replace")}


def _is_loopback(url: str) -> bool:
    from urllib.parse import urlsplit
    host = urlsplit(url).hostname or ""
    host = host.lower().rstrip(".").split(":")[0]
    return host in ("localhost", "127.0.0.1", "::1") or host.startswith("127.")
