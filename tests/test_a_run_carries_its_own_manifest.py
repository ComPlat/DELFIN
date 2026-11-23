"""A result that cannot be traced to what produced it is an anecdote.

A session report answers "what did this session do". It cannot answer
"could someone else arrive here again": the commit that ran, the model
that answered, the prompt it read, the engine versions on the machine
and whether the claimed work is present in the workspace are not
activity, and none of them was recorded in one place.

What is asserted here is the property that makes the record worth
having: **a field the manifest could not resolve is named, not
defaulted.** A manifest whose gaps look like values invites a reader to
believe a run was reproducible when the record cannot say so, which is
worse than no manifest at all.

Input: a session id, a workspace, and the settings only the caller
knows. Output: the manifest, and the list of what it does not carry.

Nothing here reads or writes the person's own stores: every test gets a
temporary state root, and the one test that drives the session-end stage
substitutes the writer.
"""
from __future__ import annotations

import json
import os

import pytest

from delfin.agent import run_manifest as rm


@pytest.fixture(autouse=True)
def _own_state(tmp_path, monkeypatch):
    """A state root of this test's own, so no real report is touched."""
    monkeypatch.setenv("DELFIN_STATE", str(tmp_path / "state"))
    monkeypatch.delenv("DELFIN_AGENT_LABEL", raising=False)


@pytest.fixture()
def workspace(tmp_path):
    ws = tmp_path / "ws"
    (ws / "calc").mkdir(parents=True)
    (ws / "calc" / "job.inp").write_text(
        "! B3LYP def2-SVP\n* xyz 0 1\nH 0 0 0\n*\n", encoding="utf-8")
    (ws / "calc" / "job.out").write_text(
        "FINAL SINGLE POINT ENERGY -1.0\nTOTAL RUN TIME: 0 days 0 hours "
        "0 minutes 3 seconds\n", encoding="utf-8")
    return ws


# --- the property the record exists for -------------------------------

def test_a_field_it_could_not_resolve_is_named(workspace):
    """No caller passed a temperature or a seed, so neither is invented."""
    man = rm.build(session_id="s1", workspace=workspace)
    assert man["model"]["temperature"] is None
    assert man["model"]["seed"] is None
    assert "model.temperature" in man["unresolved"]
    assert "model.seed" in man["unresolved"]


def test_what_the_caller_knows_is_not_reported_missing(workspace):
    man = rm.build(session_id="s1", workspace=workspace,
                   settings={"provider": "kit", "temperature": 0.0,
                             "seed": 7},
                   termination="done")
    assert man["model"]["provider"] == "kit"
    assert man["model"]["temperature"] == 0.0
    assert man["model"]["seed"] == 7
    for name in ("model.provider", "model.temperature", "model.seed",
                 "run.termination"):
        assert name not in man["unresolved"]


def test_an_absent_termination_is_not_read_as_done(workspace):
    """The record must tell a run that finished from one that was
    stopped. Defaulting to "done" would make every killed run read as a
    completed one."""
    man = rm.build(session_id="s1", workspace=workspace)
    assert man["run"]["termination"] is None
    assert "run.termination" in man["unresolved"]
    assert "run.termination" in rm.missing(man)


def test_complete_is_about_the_manifest_not_the_run(workspace):
    """A complete manifest may record a failed contract; that is the
    distinction the name protects."""
    man = rm.build(session_id="s1", workspace=workspace,
                   settings={"provider": "kit", "model": "kit.glm"},
                   termination="done")
    man["contract"]["failed"] = True
    man["contract"]["flags"] = {"engine_invoked": False}
    assert rm.complete(man) is True, rm.missing(man)


def test_a_configured_model_is_not_passed_off_as_a_measured_one(workspace):
    """A short run has no per-turn metrics, so the configured id is the
    fallback -- and the record says which of the two it is, because a
    reader comparing two runs needs that distinction."""
    man = rm.build(session_id="s1", workspace=workspace,
                   settings={"model": "kit.glm"})
    assert man["model"]["id"] == "kit.glm"
    assert man["model"]["id_source"] == "configured"
    assert "model.id" not in man["unresolved"]

    none = rm.build(session_id="s1", workspace=workspace)
    assert none["model"]["id"] is None
    assert none["model"]["id_source"] is None
    assert "model.id" in none["unresolved"]


def test_missing_reports_a_required_field_left_none(workspace):
    man = rm.build(session_id="s1", workspace=workspace,
                   termination="done")
    man["code"]["commit"] = None
    man["unresolved"] = []          # a caller that cleared the list
    assert "code.commit" in rm.missing(man)
    assert rm.complete(man) is False


# --- what it carries ---------------------------------------------------

def test_the_contract_is_the_manifests_only_judgement(workspace):
    """The oracle reads the artifacts, not the model's account of them."""
    man = rm.build(session_id="s1", workspace=workspace,
                   expected_steps=1, termination="done")
    flags = man["contract"]["flags"]
    assert flags is not None, man["unresolved"]
    assert flags.get("engine_invoked") is True
    assert man["contract"]["steps_found"] == 1
    assert man["contract"]["expected_steps"] == 1


def test_an_unjudged_flag_is_carried_as_unjudged(workspace):
    """True, False and None are three answers and stay three."""
    man = rm.build(session_id="s1", workspace=workspace,
                   termination="done")
    assert set(man["contract"]["flags"]) >= {"engine_invoked",
                                            "grid_complete"}
    assert all(v is True or v is False or v is None
               for v in man["contract"]["flags"].values())


def test_every_evidence_file_carries_a_digest(workspace):
    man = rm.build(session_id="s1", workspace=workspace,
                   termination="done")
    paths = {f["path"]: f for f in man["workspace"]["files"]}
    assert "calc/job.inp" in paths and "calc/job.out" in paths
    for entry in paths.values():
        assert entry["sha256"] and len(entry["sha256"]) == 64
        assert entry["bytes"] > 0


def test_a_truncated_file_list_says_it_was_truncated(workspace,
                                                     monkeypatch):
    """A list cut at a cap and not marked reads as a complete one."""
    monkeypatch.setattr(rm, "_MAX_FILES", 1)
    man = rm.build(session_id="s1", workspace=workspace,
                   termination="done")
    assert len(man["workspace"]["files"]) == 1
    assert man["workspace"]["omitted"] == 1


def test_the_prompt_digest_changes_with_the_pack(workspace, monkeypatch):
    """A changed instruction changes the manifest, or the digest is
    decoration."""
    first = rm.build(session_id="s1", workspace=workspace)["prompt"]
    assert first["pack_digest"] and first["pack_files"] > 0

    from delfin.agent.prompt_loader import PromptLoader
    pack = PromptLoader().agent_dir
    probe = pack / "zz_manifest_probe.md"
    probe.write_text("probe\n", encoding="utf-8")
    try:
        second = rm.build(session_id="s1", workspace=workspace)["prompt"]
    finally:
        probe.unlink()
    assert second["pack_digest"] != first["pack_digest"]
    assert second["pack_files"] == first["pack_files"] + 1


def test_the_commit_and_whether_it_was_dirty_are_separate(workspace,
                                                          monkeypatch):
    monkeypatch.setenv("DELFIN_GIT_COMMIT", "abc1234-dirty")
    man = rm.build(session_id="s1", workspace=workspace)
    assert man["code"]["commit"] == "abc1234"
    assert man["code"]["dirty"] is True


@pytest.mark.parametrize("output,expected", [
    ("*****************\n   Program Version 6.1.1  -  RELEASE   -\n", "6.1.1"),
    ("   * xtb version 6.7.1 (edcfbbe) compiled by 'conda'\n", "6.7.1"),
    ("╔════\n   Version 3.0.2, Sun, 25 August\n",
     "3.0.2"),
    ("no version anywhere\n", None),
    ("*****************\n--------------\n", None),
])
def test_a_banner_line_is_not_a_version(workspace, monkeypatch, output,
                                        expected):
    """Measured on a real run: taking the first line of the output gave
    "*****************" for ORCA and a box rule for CREST. A version
    field holding decoration reads as a version to anyone comparing two
    manifests, so it is matched out of the text or left absent."""
    monkeypatch.setattr(rm.shutil, "which",
                        lambda name: "/bin/true" if name == "orca" else None)

    class _Res:
        stdout, stderr = output, ""

    monkeypatch.setattr(rm.subprocess, "run", lambda *a, **k: _Res())
    man = rm.build(session_id="s1", workspace=workspace)
    assert man["environment"]["engines"]["orca"]["version"] == expected
    named = "environment.engines.orca.version" in man["unresolved"]
    assert named is (expected is None)


def test_an_unknown_role_is_a_gap_not_a_null(workspace):
    """Which role prompt the model read decides what it was told it may
    do, so an unknown role is named."""
    assert "prompt.role" in rm.build(
        session_id="s1", workspace=workspace)["unresolved"]
    assert "prompt.role" not in rm.build(
        session_id="s1", workspace=workspace, role="solo_agent")["unresolved"]


def test_an_engine_that_is_not_installed_is_recorded_as_absent(workspace,
                                                               monkeypatch):
    """Absent is a fact about the machine and belongs in the record; it
    must not look like an engine whose version could not be read."""
    monkeypatch.setattr(rm.shutil, "which", lambda name: None)
    man = rm.build(session_id="s1", workspace=workspace)
    assert man["environment"]["engines"]["orca"] is None
    assert not [u for u in man["unresolved"]
                if u.startswith("environment.engines.orca")]


# --- storage -----------------------------------------------------------

def test_it_is_written_atomically_and_read_back(tmp_path, workspace):
    man = rm.build(session_id="s1", workspace=workspace,
                   termination="done")
    target = tmp_path / "out" / "s1.manifest.json"
    assert rm.write(man, target) == target
    assert not list(target.parent.glob("*.tmp"))
    assert rm.read(target) == man


def test_a_schema_it_does_not_know_is_refused(tmp_path, workspace):
    """A moved field would otherwise be read as an absent one."""
    target = tmp_path / "future.json"
    man = rm.build(session_id="s1", workspace=workspace)
    man["schema"] = rm.SCHEMA + 1
    target.write_text(json.dumps(man), encoding="utf-8")
    assert rm.read(target) is None


def test_the_file_is_not_world_readable(tmp_path, workspace):
    """It names every evidence file of the run."""
    target = tmp_path / "s1.manifest.json"
    rm.write(rm.build(session_id="s1", workspace=workspace), target)
    assert (os.stat(target).st_mode & 0o077) == 0


def test_a_write_that_cannot_happen_returns_none_and_does_not_raise(
        tmp_path, workspace):
    """The manifest records a run that already finished; failing to
    store it must not change that run's outcome."""
    blocked = tmp_path / "file"
    blocked.write_text("x", encoding="utf-8")
    assert rm.write(rm.build(session_id="s1"), blocked / "deeper.json") is None


def test_no_session_id_writes_nothing(workspace):
    assert rm.write_run_manifest("", workspace=workspace) is None


def test_the_session_end_call_writes_beside_the_report(workspace):
    written = rm.write_run_manifest("sess-42", workspace=workspace,
                                    termination="done")
    assert written is not None and written.is_file()
    assert written.name.endswith(".manifest.json")
    assert written.parent.name == "session_reports"
    assert rm.read(written)["run"]["session_id"] == "sess-42"


# --- it never raises ---------------------------------------------------

@pytest.mark.parametrize("field,blow_up", [
    ("code", "delfin.utils.get_git_commit_info"),
    ("ledger", "delfin.agent.change_journal.list_changes"),
])
def test_a_collector_that_fails_leaves_its_field_unresolved(
        workspace, monkeypatch, field, blow_up):
    """One broken source must not take the whole record down with it."""
    module, _, name = blow_up.rpartition(".")
    import importlib
    monkeypatch.setattr(importlib.import_module(module), name,
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError))
    man = rm.build(session_id="s1", workspace=workspace)
    assert man["schema"] == rm.SCHEMA
    assert any(u.startswith(field) for u in man["unresolved"]), man["unresolved"]


def test_a_workspace_that_is_not_there_is_said_not_assumed():
    man = rm.build(session_id="s1", workspace="/nonexistent/at/all")
    assert "workspace.root" in man["unresolved"]
    assert "contract" in rm.missing(man)
