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
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Optional


def _iso(ts) -> Optional[str]:
    """A timestamp as UTC ISO-8601.

    ``record.build`` stores ``started``/``ended`` as ISO strings while the
    steps and measurements carry epoch seconds, so both arrive here; a value
    that is already a string is passed through rather than converted twice.
    """
    if ts is None:
        return None
    if isinstance(ts, str):
        return ts
    try:
        return datetime.fromtimestamp(float(ts), timezone.utc).isoformat()
    except (TypeError, ValueError, OSError):
        return str(ts)


FORMAT = "delfin.lab.chemotion/1"

# Chunk size used by the official client (chemotionapi analyses_manager.py).
MAX_UPLOAD_CHUNK = 5000

# Endpoint and payload facts read from source, with their primary locations.
# Re-read on 2026-11-24 for the container tree and the analyses/dataset split.
SOURCES = {
    "api.root": "ComPlat/chemotion_ELN app/api/api.rb "
                "(prefix api, version v1; mounts the *_API resources; "
                "Detection of an API token: current_token starts with 'chemtoken')",
    "authentication": "ComPlat/chemotion_ELN app/api/chemotion/authentication_api.rb",
    "samples": "ComPlat/chemotion_ELN app/api/chemotion/sample_api.rb "
               "(sample.container = update_datamodel(params[:container]))",
    "reactions": "ComPlat/chemotion_ELN app/api/chemotion/reaction_api.rb "
                 "(reaction.container = update_datamodel(container_info))",
    "containers": "ComPlat/chemotion_ELN app/usecases/containers/update_datamodel.rb "
                  "(only description, children and attachments are read -- no "
                  "analyses key) and app/models/container.rb "
                  "(create_analysis_with_dataset!: analysis -> dataset -> files)",
    "attachments": "Chemotion/ChemotionAPI chemotion_api/elements/analyses_manager.py "
                   "and elements/empty_elements/{analyses,dataset}.json",
}

UNVERIFIED = [
    "exact create params of POST /research_plans (body, collection_id, "
    "attachments) -- its api.rb was not re-read for this export",
    "the ZIP LAYOUT the ELN's own file import expects: only the JSON import "
    "model and the server's container tree are verified here, so the bundle "
    "must be checked against one ELN version before it is relied on",
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


def _attachment(filename: str, text: str) -> dict:
    ctype = mimetypes.guess_type(filename)[0] or "application/octet-stream"
    return {"filename": filename, "content_type": ctype, "text": text}


#: Format of the file a person imports. The request shapes below were read
#: from the ELN's source and the official client on 2026-11-24
#: (:data:`SOURCES`); the FILE-import layout of a collection zip is not one
#: of them, which is why it stays in :data:`UNVERIFIED`.
IMPORT_FORMAT = "delfin.lab.chemotion-import/1"


def _measurement_file(m: dict) -> str:
    """The measurement as the data file its analysis carries.

    One JSON document per measurement rather than the run's whole
    measurements array: on the server a dataset's file belongs to that
    analysis, so a shared file would attach the same bytes to every dataset.
    """
    return json.dumps(m, indent=1, sort_keys=True, default=str)


def _conditions(rec: dict) -> list[dict]:
    """The ``set`` steps of the run: what it held its devices at.

    A ``measure`` step is not a condition (it reads), a ``wait`` step is not
    a condition (it is a guard on one), and a ``transfer`` is where the
    sample went -- that belongs to the sample, not to the reaction.

    Only a set that SUCCEEDED is a condition. The kernel records a refused
    attempt as a step of its own (``code: approval_required``, with no device
    and no value), and the same setpoint then appears again once it is
    approved: counted blindly, the refusal becomes a condition holding a
    device at nothing. The refusal is not lost -- it stays in ``steps`` and in
    the audit chain, which is where provenance reads it.
    """
    out = []
    for st in rec.get("steps") or []:
        if st.get("op") != "set":
            continue
        o = st.get("outcome") or {}
        if not o.get("ok"):
            continue
        out.append({
            "step": st.get("step"),
            "device": o.get("device"),
            "channel": o.get("channel"),
            "value": o.get("value"),
            "unit": o.get("unit"),
            "level": o.get("level"),
            "readback": o.get("readback"),
            "command": o.get("command"),
            "at": _iso(st.get("ts")),
        })
    return out


def analysis_container(measurement: dict, rec: dict) -> dict:
    """One analysis as the container tree the ELN persists.

    The shape is the client's own (``analyses.json`` / ``dataset.json``), not
    an invented one: an analysis is a container of type ``analysis`` whose
    child is a dataset, the files hang on the DATASET, and the values ride in
    the dataset's ``extended_metadata``, because the server writes container
    metadata from ``extended_metadata`` and a bare key would be dropped
    (``Usecases::Containers::UpdateDatamodel#update_description``, and
    ``Container#create_analysis_with_dataset!`` for the tree).
    """
    name = f"{measurement.get('quantity') or 'measurement'} on " \
           f"{measurement.get('device') or 'device'}"
    dataset = {
        "id": str(uuid.uuid4()),
        "is_new": True,
        "is_deleted": False,
        "name": name,
        "description": "",
        "container_type": "dataset",
        "extended_metadata": {
            "report": True,
            "delfin": {
                "quantity": measurement.get("quantity"),
                "value": measurement.get("value"),
                "unit": measurement.get("unit"),
                "device": measurement.get("device"),
                "container": measurement.get("container"),
                "measured_at": _iso(measurement.get("ts")),
                "run_id": rec.get("run_id"),
                "identifier": (rec.get("fair") or {}).get("identifier"),
            },
        },
        "attachments": [{
            "id": str(uuid.uuid4()),
            "is_new": True,
            # run AND sample AND quantity: two samples measured for the same
            # quantity both produced "<run>-signal.json" before, so the file
            # an analysis carries could not be told apart from its twin.
            "filename": f"{rec.get('run_id')}-"
                        f"{measurement.get('sample') or 'sample'}-"
                        f"{measurement.get('quantity') or 'measurement'}.json",
            "content_type": "application/json",
            "checksum": hashlib.md5(
                _measurement_file(measurement).encode("utf-8")).hexdigest(),
            "raw": _measurement_file(measurement),
        }],
        "children": [],
    }
    return {
        "id": str(uuid.uuid4()),
        "is_new": True,
        "is_deleted": False,
        "name": name,
        "description": "",
        "container_type": "analysis",
        "extended_metadata": {"report": True, "index": 1},
        # the client's own split: container attachments stay empty, the
        # files hang on the dataset
        "attachments": [],
        "children": [dataset],
    }


def _provenance(rec: dict) -> dict:
    """Manifest pin, the audit chain, who approved what, and the record's hash.

    The chain is REPORTED, not claimed. ``prev`` is the SHA-256 of the
    previous line's raw bytes as written (``delfin.agent.audit_log``), so it
    cannot be re-derived from the parsed lines a record carries: it has to be
    verified on the log's bytes. ``audit_chain`` is that verification of the
    installation's log -- through ``delfin.agent.audit_log.verify``, the one
    reader of the file -- taken when the export is built. A run whose lines
    are no longer in the log it came from therefore says so instead of
    presenting itself as chained.
    """
    audit = list(rec.get("device_log") or [])
    chain = _chain_report()

    approvals: dict = {}
    for line in audit:
        rid = line.get("id")
        if not rid or not str(line.get("event", "")).startswith("approval_"):
            continue
        ap = approvals.setdefault(rid, {
            "id": rid, "device": line.get("device"),
            "channel": line.get("channel"), "value": line.get("value"),
            "unit": line.get("unit"), "session": line.get("session"),
            "requested_by": None, "approved_by": None,
            "requested_at": None, "granted_at": None, "consumed_at": None,
        })
        ev = line.get("event")
        if ev == "approval_requested":
            ap["requested_by"] = line.get("person")
            ap["requested_at"] = line.get("ts")
        elif ev == "approval_granted":
            ap["approved_by"] = line.get("approver")
            ap["granted_at"] = line.get("ts")
        elif ev == "approval_consumed":
            ap["consumed_at"] = line.get("ts")

    return {
        "installation": rec.get("installation"),
        "manifest": rec.get("manifest"),
        "plan_digest": rec.get("plan_digest"),
        "actor": rec.get("actor"),
        "started": rec.get("started"),
        "ended": rec.get("ended"),
        "status": rec.get("status"),
        "deviations": rec.get("deviations") or [],
        "audit": audit,
        "audit_lines": len(audit),
        "audit_chain": chain,
        "audit_chain_intact": bool(chain and chain.get("ok")),
        "approvals": sorted(approvals.values(), key=lambda a: str(a["id"])),
        "record_sha256": hashlib.sha256(
            json.dumps(rec, indent=1, sort_keys=True).encode("utf-8")
        ).hexdigest(),
    }


def _chain_report() -> Optional[dict]:
    """``audit_log.verify`` of the installation's active log, or None.

    None means the question was not answered (no log, or verification could
    not run), which is different from ``ok: False`` -- a chain that was
    checked and is broken. The export distinguishes the two so a reader is
    never told a run is proven when nothing was proven.
    """
    try:
        from delfin.agent import audit_log
        from . import store
        path = store.lab_dir() / "audit.log"
        if not path.exists():
            return None
        report = audit_log.verify(path)
    except Exception:                             # noqa: BLE001 - reported as None
        return None
    if not isinstance(report, dict):
        return None
    return {k: report.get(k) for k in
            ("ok", "records", "chained", "legacy", "problems", "path")
            if k in report}


def _fair(rec: dict) -> dict:
    """The record's FAIR slots, plus the ones a repository asks for.

    ``license`` is reported as it is (the record leaves it unset) together
    with ``license_set``: an identifier without a licence is not reusable,
    and silently filling one in would claim a right nobody granted.
    """
    fair = dict(rec.get("fair") or {})
    fair.setdefault("identifier", None)
    fair.setdefault("creator", None)
    fair.setdefault("license", None)
    fair.setdefault("created", None)
    fair["license_set"] = fair.get("license") is not None
    fair["metadata_standard"] = IMPORT_FORMAT
    fair["schema"] = rec.get("schema")
    fair["searchable"] = {
        "run_id": rec.get("run_id"),
        "installation": rec.get("installation"),
        "plan_id": (rec.get("plan") or {}).get("id"),
        "samples": sorted((rec.get("samples") or {})),
        "quantities": sorted({m.get("quantity")
                              for m in rec.get("measurements") or []
                              if m.get("quantity")}),
    }
    return fair


def build_import(rec: dict) -> dict:
    """Map a ``delfin.lab.run/1`` record to the ELN import model.

    One research-plan element (body = the run parameters), one reaction
    carrying the run's conditions and the sample as its product, one sample
    entry per run sample with its analyses in the container tree the ELN
    persists, the full record / device log / measurements / parameters as
    attachments with their bytes, and the provenance and FAIR sections.
    """
    by_sample: dict = {}
    for m in rec.get("measurements") or []:
        by_sample.setdefault(m.get("sample"), []).append(m)

    samples = []
    for sid, s in (rec.get("samples") or {}).items():
        measurements = by_sample.get(sid, [])
        samples.append({
            "name": sid,
            "external_label": s.get("container"),
            "location": f"{s.get('module')}:{s.get('position')}",
            "container": {
                "id": str(uuid.uuid4()),
                "is_new": True,
                "is_deleted": False,
                "name": "root",
                "description": "",
                "container_type": "root",
                "extended_metadata": {"report": True},
                "attachments": [],
                "children": [analysis_container(m, rec)
                             for m in measurements],
            },
            "analyses": [
                {"kind": m.get("quantity"), "device": m.get("device"),
                 "value": m.get("value"), "unit": m.get("unit"),
                 "container": m.get("container")}
                for m in measurements
            ],
        })

    parameters = {
        "format": IMPORT_FORMAT,
        "run_id": rec.get("run_id"),
        "schema": rec.get("schema"),
        "status": rec.get("status"),
        "plan": rec.get("plan"),
        "plan_digest": rec.get("plan_digest"),
        "manifest": rec.get("manifest"),
        "installation": rec.get("installation"),
        "actor": rec.get("actor"),
        "started": rec.get("started"),
        "ended": rec.get("ended"),
        "measurements": rec.get("measurements") or [],
        "samples": sorted((rec.get("samples") or {})),
        "deviations": rec.get("deviations") or [],
        "fair": _fair(rec),
    }

    attachments = [
        _attachment("record.json", json.dumps(rec, indent=1, sort_keys=True)),
        _attachment("device_log.json",
                    json.dumps(rec.get("device_log") or [], indent=1)),
        _attachment("measurements.json",
                    json.dumps(rec.get("measurements") or [], indent=1)),
        _attachment("parameters.json",
                    json.dumps(parameters, indent=1, sort_keys=True)),
        _attachment("provenance.json",
                    json.dumps(_provenance(rec), indent=1, sort_keys=True)),
        _attachment("fair.json",
                    json.dumps(_fair(rec), indent=1, sort_keys=True)),
    ]
    # one data file per measurement, the bytes of the analysis that carries it
    for sid in sorted(by_sample):
        for i, m in enumerate(by_sample[sid], start=1):
            attachments.append(_attachment(
                f"{rec.get('run_id')}-{sid}-{m.get('quantity') or 'm'}-{i}.json",
                _measurement_file(m)))

    products = [{"sample": s["name"], "container": s.get("external_label"),
                 "location": s.get("location"), "module": (
                     (rec.get("samples") or {}).get(s["name"]) or {}
                 ).get("module")}
                for s in samples]
    reactions = [{
        "name": f"Automated run {rec.get('run_id')}",
        "type": "reaction",
        "status": rec.get("status"),
        "conditions": _conditions(rec),
        "products": products,
        "starting_materials": [],
        "reactants": [],
        "plan_digest": rec.get("plan_digest"),
        "duration_s": _duration_s(rec),
    }]

    return {
        "format": IMPORT_FORMAT,
        "element": {
            "type": "research_plan",
            "name": f"Automated run {rec.get('run_id')}",
            # the record itself, not a three-key summary of it: this is what
            # the plan element's body holds for a reader of the ELN.
            "body": dict(rec),
        },
        "reactions": reactions,
        "samples": samples,
        "attachments": attachments,
        "provenance": _provenance(rec),
        "fair": _fair(rec),
        "unverified": list(UNVERIFIED),
    }


def _duration_s(rec: dict) -> Optional[float]:
    """Wall time of the run, when both stamps are readable."""
    try:
        a = datetime.fromisoformat(str(rec["started"]))
        b = datetime.fromisoformat(str(rec["ended"]))
    except (KeyError, TypeError, ValueError):
        return None
    return round((b - a).total_seconds(), 3)


def write_import(model: dict, out: Path | str) -> tuple[Path, dict]:
    """Write ``model`` as one import file, plus its readable JSON twin.

    The zip is the artefact a person hands to the ELN. The twin
    (``chemotion_import.json``) is the same content uncompressed, because a
    reviewer cannot read a zip in a diff and because the export that already
    exists is JSON. Nothing here sends anything.
    """
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    members = {"chemotion_import.json": json.dumps(
        model, indent=1, sort_keys=True, default=str)}
    for a in model.get("attachments") or []:
        members[f"attachments/{a['filename']}"] = a.get("text") or ""
    # the manifest is a member too, so it must say so: a reader comparing its
    # list with the archive would otherwise find one name missing.
    manifest = {
        "format": IMPORT_FORMAT,
        "created": _iso(datetime.now(timezone.utc).timestamp()),
        "members": [],
        "record_sha256": (model.get("provenance") or {}).get("record_sha256"),
        "identifier": (model.get("fair") or {}).get("identifier"),
        "license_set": (model.get("fair") or {}).get("license_set"),
        "unverified": list(model.get("unverified") or []),
    }
    manifest["members"] = sorted(members) + ["manifest.json"]
    members["manifest.json"] = json.dumps(manifest, indent=1, sort_keys=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name in sorted(members):
            z.writestr(name, members[name])
    (path.parent / "chemotion_import.json").write_text(
        json.dumps(model, indent=1, sort_keys=True, default=str),
        encoding="utf-8")
    return path, manifest


def import_files(model: dict) -> dict:
    """The bundle's members as ``{name: text}``, without writing anything."""
    out = {"chemotion_import.json": json.dumps(
        model, indent=1, sort_keys=True, default=str)}
    for a in model.get("attachments") or []:
        out[f"attachments/{a['filename']}"] = a.get("text") or ""
    return out


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
