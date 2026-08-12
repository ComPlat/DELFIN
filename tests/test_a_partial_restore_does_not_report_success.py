"""An almost-empty restore looked exactly like a complete one.

The session file carried twenty-two keys and no version. ``restore_state``
returned None. ``_restore_evidence`` swallowed every per-field exception
into an empty default -- no counter, no log, no return value. So a
restore that recovered two fields out of twenty and a restore that
recovered all twenty were, to every caller, the same event: nothing was
raised, so it worked.

That is worse than a failure. The session goes on to answer with the
guards enabled over a record it did not actually get back, and the one
place a human might have noticed -- an unattended run's log -- said
nothing at all.

Bundles were the same defect with the evidence sitting right there: the
exporter wrote a manifest carrying ``kind`` and ``version`` from the
first release, and the importer never opened the file. Any zip
containing a ``session.json`` was ingested as ours, and a bundle from a
newer build was read with today's loaders.

So: a schema field on every file this build writes; a restore report
naming what came back, what was missing and what failed, surfaced by the
callers; a version NEWER than the code is a stated refusal rather than a
partial read; older is a NAMED migration rather than a silent default.
"""

from __future__ import annotations

import json
import zipfile

import pytest

from delfin.agent import session_store as ss
from delfin.agent.engine import AgentEngine


@pytest.fixture()
def sessions_dir(tmp_path, monkeypatch):
    d = tmp_path / "sessions"
    d.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", d)
    monkeypatch.setattr(ss, "_ensure_dir", lambda: d)
    monkeypatch.setattr(ss, "acquire_session_lock", lambda *a, **k: True)
    return d


def _engine() -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())
    eng.mode = "solo"
    eng.route = []
    eng.messages = []
    eng.session_id = ""
    eng.client = None
    return eng


# ---------------------------------------------------------------------------
# The file says which schema it is
# ---------------------------------------------------------------------------

def test_a_saved_session_carries_its_schema_version(sessions_dir):
    ss.save_session("s1", mode="solo")
    data = json.loads((sessions_dir / "s1.json").read_text())
    assert data["schema_version"] == ss.SESSION_SCHEMA_VERSION


# ---------------------------------------------------------------------------
# The restore says what it actually did
# ---------------------------------------------------------------------------

def test_a_complete_restore_reports_complete():
    src = _engine()
    src.session_id = "s1"
    src.mode = "solo"
    state = src.export_state()
    state["schema_version"] = ss.SESSION_SCHEMA_VERSION
    report = _engine().restore_state(state)
    assert report.complete, report.summary()
    assert not report.missing and not report.failed


def test_a_thin_restore_names_what_is_missing():
    """The case that used to be indistinguishable from the one above."""
    report = _engine().restore_state({
        "schema_version": ss.SESSION_SCHEMA_VERSION,
        "mode": "solo", "engine_messages": [], "session_id": "s1",
    })
    assert not report.complete
    assert "_exec_commands_session" in report.missing
    assert "role_test_evidence" in report.missing
    assert "_project_dir" in report.missing


def test_a_field_that_cannot_be_read_is_reported_as_failed(monkeypatch):
    """Swallowed into an empty default, with a counter this time."""
    boom = AgentEngine._SessionField(
        "_exec_commands_session", "exec_commands", "evidence",
        list, lambda v: (_ for _ in ()).throw(ValueError("bad")), list)
    others = tuple(f for f in AgentEngine._SESSION_FIELDS
                   if f.attr != "_exec_commands_session")
    monkeypatch.setattr(AgentEngine, "_SESSION_FIELDS", others + (boom,))
    report = _engine().restore_state({
        "schema_version": ss.SESSION_SCHEMA_VERSION,
        "mode": "solo", "engine_messages": [], "session_id": "s1",
        "evidence": {"exec_commands": ["pytest"]},
    })
    assert "_exec_commands_session" in report.failed
    assert not report.complete


def test_the_report_says_it_in_one_line():
    report = _engine().restore_state(
        {"mode": "solo", "engine_messages": [], "session_id": "s1"})
    text = report.summary()
    assert "schema v" in text and "restored" in text and "missing" in text


# ---------------------------------------------------------------------------
# A newer schema is refused, an older one is migrated by name
# ---------------------------------------------------------------------------

def test_a_session_from_a_newer_build_is_refused():
    with pytest.raises(ss.SessionSchemaError) as exc:
        _engine().restore_state({
            "schema_version": ss.SESSION_SCHEMA_VERSION + 1,
            "mode": "solo", "engine_messages": [], "session_id": "s1",
        })
    assert "refus" in str(exc.value).lower()


def test_an_older_session_is_migrated_and_says_so():
    report = _engine().restore_state({
        "schema_version": 1,
        "mode": "solo", "engine_messages": [], "session_id": "s1",
        "cost_usd": 4.0,
    })
    assert report.migrations
    assert any("v1_to_v2" in m for m in report.migrations)


def test_an_unversioned_session_is_read_as_the_first_schema():
    report = _engine().restore_state(
        {"mode": "solo", "engine_messages": [], "session_id": "s1"})
    assert report.schema_version == 1
    assert report.migrations


def test_the_migration_does_not_mutate_the_caller_dict():
    data = {"mode": "solo", "engine_messages": [], "session_id": "s1"}
    _engine().restore_state(data)
    assert "schema_version" not in data


def test_the_migration_gives_the_cost_baseline_the_current_total():
    """A v1 file has no baseline. Zero would book the whole session's
    spend as the first turn after the resume; the current total books
    nothing, which is the honest answer."""
    eng = _engine()
    eng.restore_state({
        "mode": "solo", "engine_messages": [], "session_id": "s1",
        "cost_usd": 4.0,
    })
    assert eng._last_outcome_cost == 4.0


# ---------------------------------------------------------------------------
# The bundle manifest is opened
# ---------------------------------------------------------------------------

def _bundle(tmp_path, *, manifest: dict | None, session: dict) -> "object":
    p = tmp_path / "b.delfin-session"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("session.json", json.dumps(session))
        if manifest is not None:
            z.writestr("manifest.json", json.dumps(manifest))
    return p


def test_a_bundle_from_a_newer_build_is_refused(tmp_path, sessions_dir):
    p = _bundle(tmp_path, manifest={
        "kind": ss.BUNDLE_KIND,
        "version": ss.BUNDLE_FORMAT_VERSION + 1,
    }, session={"session_id": "x"})
    with pytest.raises(ss.SessionSchemaError):
        ss.import_bundle(p)


def test_a_bundle_carrying_a_newer_session_is_refused(tmp_path, sessions_dir):
    p = _bundle(tmp_path, manifest={"kind": ss.BUNDLE_KIND, "version": 1},
                session={"session_id": "x",
                         "schema_version": ss.SESSION_SCHEMA_VERSION + 1})
    with pytest.raises(ss.SessionSchemaError):
        ss.import_bundle(p)


def test_a_zip_that_is_not_a_session_bundle_is_refused(tmp_path, sessions_dir):
    p = _bundle(tmp_path, manifest={"kind": "something-else", "version": 1},
                session={"session_id": "x"})
    with pytest.raises(ss.SessionSchemaError):
        ss.import_bundle(p)


def test_our_own_bundle_still_round_trips(tmp_path, sessions_dir, monkeypatch):
    monkeypatch.setattr(ss, "_bundles_dir", lambda: tmp_path)
    monkeypatch.setattr(ss, "_transcript_archive_dir", lambda: tmp_path)
    ss.save_session("s1", mode="solo", title="a session")
    path = ss.export_bundle("s1")
    assert path is not None
    new_id = ss.import_bundle(path)
    assert new_id
    assert ss.load_session(new_id)["mode"] == "solo"


def test_a_bundle_with_no_manifest_at_all_still_imports(tmp_path,
                                                       sessions_dir):
    """Bundles written before the manifest existed must stay readable."""
    p = _bundle(tmp_path, manifest=None, session={"session_id": "x"})
    assert ss.import_bundle(p)
