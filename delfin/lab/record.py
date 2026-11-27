"""The run record, and its export for Chemotion ELN.

Input: a finished run's checkpoint and the manifest. Output:
``runs/<run_id>/record.json`` (schema ``delfin.lab.run/1``),
``runs/<run_id>/chemotion_import.json`` and
``runs/<run_id>/chemotion_import.zip``.

The record is what a later reader -- a person, an ELN, the ``ml_train``
package closing the self-driving-lab loop -- needs without asking anyone:
the plan and its digest, the manifest version and checksum it ran against,
who ran it, every step's outcome, the measurements with units, every
deviation, where each sample went, and the lab audit lines of this run.
FAIR fields (Wilkinson et al. 2016): a persistent run identifier, the
creator, a licence slot, provenance (software and version, manifest
checksum) and timestamps in UTC ISO-8601.

Chemotion import. ``chemotion.build_import()`` maps the record to the ELN and
``chemotion.write_import()`` writes it as one zip a person uploads: the
research plan, the reaction with its conditions, one sample per run sample
whose analyses carry their data files, the record and its provenance, and the
FAIR metadata. The shapes follow what was read in the ELN's own source and in
the official client (``chemotion.SOURCES``); what could NOT be verified stays
in the payload's ``unverified`` list. Nothing is uploaded: an upload needs the
API token of a person, and ``chemotion.upload()`` is the only code that sends.
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
    """The run as the ELN import model (``delfin.lab.chemotion-import/1``).

    Delegates to :func:`delfin.lab.chemotion.build_import`, so the export the
    run writes and the model the request sequence uploads are ONE shape. A
    second, thinner mapping here (format ``chemotion-export/0``) is what let a
    finished run leave a file that named one attachment and carried none.
    """
    from .chemotion import build_import
    return build_import(rec)


def write(cp: dict, manifest) -> dict:
    """Write the run's record and the file a person imports into the ELN.

    Three files: ``record.json`` (the record), ``chemotion_import.json`` (its
    readable twin) and ``chemotion_import.zip`` (the artefact). The import is
    derived and written on every call, so it can always be rebuilt from a
    record that already exists. A failure to write it is logged and does not
    lose the record: the record is the primary artefact and a missing import
    is reproducible from it, while an unwritten record is not.
    """
    rec = build(cp, manifest)
    d = run_dir(cp["run_id"])
    (d / "record.json").write_text(json.dumps(rec, indent=1, sort_keys=True),
                                   encoding="utf-8")
    _write_import(rec, d)
    return rec


def _write_import(rec: dict, d) -> None:
    from . import chemotion
    try:
        chemotion.write_import(chemotion.build_import(rec),
                               d / "chemotion_import.zip")
    except Exception as exc:                      # noqa: BLE001 - see docstring
        import logging
        logging.getLogger(__name__).error(
            "run %s: record written, Chemotion import not written: %s",
            rec.get("run_id"), exc)
