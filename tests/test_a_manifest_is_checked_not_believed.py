"""A record nobody can test is a record nobody has to be right about.

The run manifest claims a result can be arrived at again. This is the
check of that claim. What is asserted here is the rule that makes the
check worth running: **a field that cannot be re-derived is never a
pass.** Three verdicts, never two — match, differs, unavailable — and
``reproduced()`` is true only when something was compared and nothing is
open.

The same rule as the execution contract's unjudged flag, for the same
reason: reading absence as agreement is what turns a provenance record
into decoration.

Also asserted: the comparison uses ONE collector. Replay builds a fresh
manifest with ``run_manifest.build`` and diffs it, so "replay disagrees"
can only mean the run changed, never that two readers drifted apart.

Input: a stored manifest and the tree as it is now. Output: the report.
Nothing here touches the person's own stores.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent import run_manifest as rm
from delfin.agent import run_replay as rr


@pytest.fixture(autouse=True)
def _own_state(tmp_path, monkeypatch):
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


@pytest.fixture()
def stored(workspace):
    return rm.build(session_id="s1", workspace=workspace, expected_steps=1,
                    termination="done", role="solo_agent",
                    settings={"model": "kit.glm", "provider": "kit",
                              "temperature": 0.0, "seed": 7})


def _verdict(report, name):
    for c in report.checks:
        if c.name == name:
            return c.verdict
    return None


# --- the rule the check exists for ------------------------------------

def test_an_unchanged_run_reproduces(stored, workspace):
    """Nothing moved, so every compared field agrees."""
    rep = rr.replay(stored, workspace=workspace)
    assert rep.differs() == [], [vars(c) for c in rep.differs()]
    assert rep.reproduced() is True, [vars(c) for c in rep.unavailable()]


def test_a_recorded_value_that_cannot_be_re_derived_is_not_a_pass(
        stored, workspace, monkeypatch):
    """The whole point. A claim the record made and this cannot check
    leaves the question open, and open is not agreement.

    The example matters: a field the record NEVER held is a gap in the
    record, not an unverifiable claim, and the two want different work.
    Here the record holds a prompt digest and the pack cannot be read,
    so the claim exists and cannot be checked.
    """
    def _no_pack(*a, **k):
        raise RuntimeError("pack unreadable")

    monkeypatch.setattr(
        "delfin.agent.prompt_loader.PromptLoader.__init__", _no_pack)
    rep = rr.replay(stored, workspace=workspace)
    assert _verdict(rep, "prompt.pack_digest") == rr.UNAVAILABLE
    assert "prompt.pack_digest" in {c.name for c in rep.open_questions()}
    assert rep.reproduced() is False
    assert rep.differs() == []          # not a difference: it is unknown


def test_a_gap_in_the_record_is_reported_as_a_gap_not_as_a_change(
        stored, workspace):
    """Whether the record is complete is the manifest's question. Replay
    carries the answer so a reader sees it, and does not turn it into a
    verification failure -- an incomplete record and a changed result
    call for different work."""
    stored["prompt"]["pack_digest"] = None
    rep = rr.replay(stored, workspace=workspace)
    assert "prompt.pack_digest" in rep.incomplete
    assert "prompt.pack_digest" in {c.name for c in rep.consistent_unknowns()}
    assert rep.reproduced() is True
    assert "incomplete" in rr.render(rep)


def test_a_changed_commit_is_a_difference(stored, workspace):
    stored["code"]["commit"] = "0000000"
    rep = rr.replay(stored, workspace=workspace)
    assert _verdict(rep, "code.commit") == rr.DIFFERS
    assert rep.reproduced() is False


def test_an_edited_evidence_file_is_caught_by_its_digest(stored, workspace):
    """The claim rests on these bytes; a digest is how that is checked."""
    (workspace / "calc" / "job.out").write_text(
        "FINAL SINGLE POINT ENERGY -99.0\n", encoding="utf-8")
    rep = rr.replay(stored, workspace=workspace)
    assert _verdict(rep, "workspace.calc/job.out") == rr.DIFFERS
    assert _verdict(rep, "workspace.calc/job.inp") == rr.MATCH


def test_a_deleted_evidence_file_is_a_difference_not_an_absence(
        stored, workspace):
    """The workspace is here and the file is not: the evidence a result
    was read from is gone, which is a difference in the result."""
    (workspace / "calc" / "job.out").unlink()
    rep = rr.replay(stored, workspace=workspace)
    assert _verdict(rep, "workspace.calc/job.out") == rr.DIFFERS


def test_a_missing_workspace_is_unavailable_not_a_difference(stored,
                                                             tmp_path):
    """On another machine the files are simply not here. Calling that a
    difference would accuse a run that may be perfectly reproducible."""
    rep = rr.replay(stored, workspace=tmp_path / "nowhere")
    assert _verdict(rep, "workspace.files") == rr.UNAVAILABLE
    assert rep.reproduced() is False


def test_an_engine_that_is_gone_is_a_difference(stored, workspace,
                                                monkeypatch):
    """A number that came out of a program no longer here cannot be
    trusted, so its absence is reported as a difference."""
    stored["environment"]["engines"]["orca"] = {"path": "/opt/orca/orca",
                                                "version": "6.1.1"}
    monkeypatch.setattr(rm.shutil, "which", lambda name: None)
    rep = rr.replay(stored, workspace=workspace)
    assert _verdict(rep, "environment.engines.orca") == rr.DIFFERS


def test_an_engine_that_appeared_since_is_not_a_difference(stored,
                                                           workspace):
    """It cannot have changed a result that was already recorded."""
    stored["environment"]["engines"]["orca"] = None
    rep = rr.replay(stored, workspace=workspace)
    assert _verdict(rep, "environment.engines.orca") == rr.MATCH


def test_each_contract_flag_is_its_own_check(stored, workspace):
    """"The contract differs" would hide which flag did."""
    names = {c.name for c in rr.replay(stored, workspace=workspace).checks}
    assert {"contract.engine_invoked", "contract.grid_complete"} <= names


def test_a_flag_unjudged_when_recorded_is_never_silently_a_match(
        stored, workspace):
    """None is not equal to None. The flag is judged NOW and was not
    then, so the record claimed nothing about it -- reported as unknown,
    never as agreement, and it does not fail the verification because
    there was no claim to fail."""
    stored["contract"]["flags"]["engine_invoked"] = None
    rep = rr.replay(stored, workspace=workspace)
    assert _verdict(rep, "contract.engine_invoked") == rr.UNAVAILABLE
    assert "contract.engine_invoked" not in {c.name for c in rep.matches()}
    assert "contract.engine_invoked" in {
        c.name for c in rep.consistent_unknowns()}


# --- it refuses what it cannot read -----------------------------------

def test_a_schema_it_does_not_know_compares_nothing(stored, workspace):
    """A moved field would be read as an absent one, so nothing is
    compared rather than compared wrongly -- and `reproduced` is false
    because an empty comparison is not a successful one."""
    stored["schema"] = rm.SCHEMA + 1
    rep = rr.replay(stored, workspace=workspace)
    assert rep.checks == []
    assert rep.reproduced() is False
    assert any("schema" in n for n in rep.notes)


def test_an_empty_comparison_never_reports_reproduced():
    rep = rr.Report()
    assert rep.checks == []
    assert rep.reproduced() is False


def test_an_unreadable_file_is_reported_and_does_not_raise(tmp_path):
    bad = tmp_path / "not.json"
    bad.write_text("{nope", encoding="utf-8")
    rep = rr.replay_file(bad)
    assert rep.checks == []
    assert rep.notes and rep.reproduced() is False


# --- one collector, not two -------------------------------------------

def test_replay_derives_through_the_manifests_own_collector(stored,
                                                            workspace,
                                                            monkeypatch):
    """Pinned so a second collector cannot be written for this module:
    if replay stopped going through run_manifest.build, "replay
    disagrees" would start meaning "the two readers disagree"."""
    called = {}
    real = rm.build

    def _spy(**kw):
        called.update(kw)
        return real(**kw)

    monkeypatch.setattr(rm, "build", _spy)
    rr.replay(stored, workspace=workspace)
    assert called, "replay did not build a fresh manifest"
    assert called.get("session_id") == "s1"
    assert str(called.get("workspace")) == str(workspace)


def test_the_fields_that_must_not_be_compared_are_named(stored, workspace):
    """written_at and the activity counters differ between two honest
    readings; the report says so rather than leaving a reader to wonder
    why they are missing."""
    rep = rr.replay(stored, workspace=workspace)
    names = {c.name for c in rep.checks}
    assert not any(n.startswith(("activity.", "written_at", "run.wall"))
                   for n in names)
    assert any("not compared" in n for n in rep.notes)


# --- what a person sees ------------------------------------------------

def test_the_rendering_puts_what_is_wrong_first(stored, workspace):
    """A reader scanning the top must see the differences and the open
    fields without reading the matches."""
    stored["code"]["commit"] = "0000000"
    stored["prompt"]["pack_digest"] = None
    text = rr.render(rr.replay(stored, workspace=workspace))
    first_differs = text.index(rr.DIFFERS)
    first_unavailable = text.index(rr.UNAVAILABLE)
    first_match = text.index("  " + rr.MATCH)
    assert first_differs < first_unavailable < first_match
    assert "NOT REPRODUCED" in text


def test_the_report_is_json_serialisable(stored, workspace):
    """A program reads this too, and the counts are the summary it needs."""
    data = rr.replay(stored, workspace=workspace).as_dict()
    json.dumps(data)
    assert data["counts"][rr.MATCH] > 0
    assert data["reproduced"] is True
    assert data["session_id"] == "s1"


# --- the distinction the first draft got wrong -------------------------

def test_unknown_then_and_now_is_consistent_not_a_failed_check(stored,
                                                                workspace):
    """`no_egress` is unjudged for a workspace that logged no egress,
    then and now. Calling that "not reproduced" made the predicate false
    for nearly every run, and a predicate that is always false says
    nothing. Nothing was claimed, so nothing failed to verify -- but it
    is still shown.
    """
    rep = rr.replay(stored, workspace=workspace)
    open_names = {c.name for c in rep.open_questions()}
    quiet_names = {c.name for c in rep.consistent_unknowns()}
    assert "contract.no_egress" in quiet_names
    assert "contract.no_egress" not in open_names
    assert rep.reproduced() is True


def test_a_recorded_judgement_that_cannot_be_judged_now_does_block(
        stored, tmp_path):
    """The other side of it: the record made a claim and there is no way
    to check it here, so the question stays open."""
    rep = rr.replay(stored, workspace=tmp_path / "gone")
    assert any(c.name == "contract.flags" for c in rep.open_questions())
    assert rep.reproduced() is False


def test_a_field_the_record_never_held_does_not_block(stored, workspace):
    """Nothing was claimed. It is reported so a reader sees the gap, and
    it does not turn into a verification failure."""
    stored["model"]["id"] = None
    rep = rr.replay(stored, workspace=workspace)
    assert "model.id" in {c.name for c in rep.consistent_unknowns()}
    assert rep.reproduced() is True
