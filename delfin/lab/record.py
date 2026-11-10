"""The run record, and its export for Chemotion ELN.

Input: a finished run's checkpoint and the manifest. Output:
``runs/<run_id>/record.json`` (schema ``delfin.lab.run/1``) and
``runs/<run_id>/chemotion_export.json``.

The record is what a later reader -- a person, an ELN, the ``ml_train``
package closing the self-driving-lab loop -- needs without asking anyone:
the plan and its digest, the manifest version and checksum it ran against,
who ran it, every step's outcome, the measurements with units, every
deviation, where each sample went, and the lab audit lines of this run.
FAIR fields (Wilkinson et al. 2016): a persistent run identifier, the
creator, a licence slot, provenance (software and version, manifest
checksum) and timestamps in UTC ISO-8601.

Chemotion export. The payload follows the Chemotion data model as far as
it was verified (docs/design/lab-control.md section 8): a research-plan
style element with the run as body and the record attached as a JSON
file, and one entry per sample with its analyses. What is NOT verified
(the nesting of analyses/containers/datasets, the zip layout of a
collection import) is named in the payload's ``unverified`` list, and
nothing is uploaded: an upload needs the API token of a person.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from typing import Optional

from . import audit as _audit
from .workflow import run_dir

SCHEMA = "delfin.lab.run/1"


def _iso(ts: Optional[float]) -> Optional[str]:
    if ts is None:
        return None
    return datetime.fromtimestamp(float(ts), timezone.utc).isoformat()


def plan_digest(workflow: dict) -> str:
    blob = json.dumps(workflow, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


def build(cp: dict, manifest) -> dict:
    try:
        from importlib.metadata import version
        sw = version("delfin-complat")
    except Exception:
        sw = "unknown"
    run_id = cp["run_id"]
    audit_lines = [r for r in _audit.read_all()
                   if r.get("session") == run_id or r.get("run") == run_id]
    return {
        "schema": SCHEMA,
        "run_id": run_id,
        "status": cp.get("status"),
        "installation": manifest.installation,
        "manifest": {"version": manifest.version,
                     "checksum": manifest.checksum},
        "plan": cp.get("workflow"),
        "plan_digest": plan_digest(cp.get("workflow") or {}),
        "actor": {"person": cp.get("person"),
                  "coordinator": cp.get("coordinator")},
        "started": _iso(cp.get("started")),
        "ended": _iso(cp.get("ended")),
        "steps": cp.get("log") or [],
        "measurements": cp.get("measurements") or [],
        "deviations": cp.get("deviations") or [],
        "samples": cp.get("samples") or {},
        "device_log": audit_lines,
        "fair": {
            "identifier": f"delfin-lab-run:{manifest.installation}:{run_id}",
            "creator": cp.get("person"),
            "license": None,
            "provenance": {"software": "DELFIN", "version": sw,
                           "manifest_checksum": manifest.checksum},
            "created": _iso(cp.get("ended") or cp.get("started")),
        },
    }


def to_chemotion(rec: dict) -> dict:
    samples = []
    for sid, s in (rec.get("samples") or {}).items():
        samples.append({
            "name": sid,
            "external_label": s.get("container"),
            "location": f"{s.get('module')}:{s.get('position')}",
            "analyses": [{"kind": m.get("quantity"), "device": m.get("device"),
                          "value": m.get("value"), "unit": m.get("unit"),
                          "measured_at": _iso(m.get("ts"))}
                         for m in rec.get("measurements") or []
                         if m.get("sample") == sid],
        })
    return {
        "format": "delfin.lab.chemotion-export/0",
        "element": {"type": "research_plan",
                    "name": f"Automated run {rec['run_id']}",
                    "body": {"status": rec.get("status"),
                             "plan_digest": rec.get("plan_digest"),
                             "deviations": rec.get("deviations")}},
        "samples": samples,
        "attachments": [{"filename": "record.json",
                         "content_type": "application/json"}],
        "unverified": [
            "nesting of analyses -> containers -> datasets",
            "zip layout of a collection import",
            "LabIMotion generic-element field names for a run",
        ],
    }


def write(cp: dict, manifest) -> dict:
    rec = build(cp, manifest)
    d = run_dir(cp["run_id"])
    (d / "record.json").write_text(json.dumps(rec, indent=1, sort_keys=True),
                                   encoding="utf-8")
    (d / "chemotion_export.json").write_text(
        json.dumps(to_chemotion(rec), indent=1, sort_keys=True),
        encoding="utf-8")
    return rec
