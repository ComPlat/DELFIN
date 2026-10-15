"""The share path: target, visibility, the human-only dialog, the queue and
the branch. Every repository here is a temporary `git init`; nothing is
pushed and no network is used (the visibility probe is replaced).
"""
from __future__ import annotations

import json
import os
import stat
import subprocess
from pathlib import Path

import pytest

from delfin.agent import knowledge as kn
from delfin.agent import knowledge_private as kp
from delfin.agent import knowledge_share as ks
from delfin.agent import repo_visibility as rv

LESSON = ("Situation: a parallel pytest run of a large suite reports "
          "failures. Lesson: re-run the failing tests serially before "
          "believing them. Why: shared state between workers invents "
          "failures.")


class Clock:
    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t


def _git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


def _init_repo(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    _git(path, "init", "-q", "-b", "main")
    (path / "README.md").write_text("x\n")
    _git(path, "add", "README.md")
    _git(path, "commit", "-q", "-m", "init")
    return path


@pytest.fixture(autouse=True)
def env(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    for k in ("GIT_AUTHOR_NAME", "GIT_COMMITTER_NAME"):
        monkeypatch.setenv(k, "tester")
    for k in ("GIT_AUTHOR_EMAIL", "GIT_COMMITTER_EMAIL"):
        monkeypatch.setenv(k, "tester@example.invalid")
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    monkeypatch.setenv("DELFIN_KNOWLEDGE_QUEUE", str(tmp_path / "queue"))
    delfin = _init_repo(tmp_path / "delfin-src")
    gdir = delfin / "delfin" / "agent" / "pack" / "knowledge"
    gdir.mkdir(parents=True)
    (delfin / "delfin" / "agent" / "pack" / "manifest.yaml").write_text("x")
    kn.write_topics(gdir, {t: "" for t in kn.DEFAULT_TOPICS})
    kn.write_index(gdir)
    _git(delfin, "add", "-A")
    _git(delfin, "commit", "-q", "-m", "pack")
    monkeypatch.setattr(ks, "_DELFIN_ROOT", delfin)
    monkeypatch.setattr(kn, "global_dir", lambda: gdir)
    answers = {"value": "unknown"}
    monkeypatch.setattr(rv, "_RUNNER", lambda _root: answers["value"])
    project = _init_repo(tmp_path / "proj")
    return {"home": home, "delfin": delfin, "gdir": gdir, "vis": answers,
            "project": project, "queue": tmp_path / "queue"}


def _lesson(env, text=LESSON, topic="tools-and-sandbox"):
    return kp.propose(text, evidence=["test:tests/test_x.py::t"],
                      source="agent", topic=topic,
                      title="Serial re-run before believing failures",
                      workspace=str(env["project"]))


def _share(env, les, layer="global", **kw):
    kw.setdefault("handle", "maxh")
    kw.setdefault("workspace", env["project"])
    kw.setdefault("clock", Clock())
    return ks.start_share(les.id, layer, **kw)


# ---------------------------------------------------------------------------
# visibility
# ---------------------------------------------------------------------------

def test_unknown_visibility_counts_as_public(env):
    d = _share(env, _lesson(env))
    assert d.target.public
    assert ks.PUBLIC_BANNER in d.preview()
    assert "PUBLIC: anyone can read this, and it cannot be taken back " \
           "after merging" in d.preview()


def test_a_private_repository_gets_the_base_rules_and_one_step(env):
    env["vis"]["value"] = "private"
    les = _lesson(env, text=LESSON.replace("failures.", "failures at "
                                           "12.5 kcal/mol."))
    clock = Clock()
    d = _share(env, les, layer="project", clock=clock)
    assert not d.target.public and not d.blocked, d.preview()
    out = ks.confirm(d.id, d.code, clock=clock)
    assert out.status == "done", out.message


# ---------------------------------------------------------------------------
# the check blocks, names the span, and never rewrites
# ---------------------------------------------------------------------------

def test_a_leak_blocks_the_dialog_and_an_edit_rechecks(env):
    les = _lesson(env, text=LESSON.replace("believing them",
                                           "believing them (ask a@b.org)"))
    d = _share(env, les)
    assert d.blocked and d.code == ""
    assert "a@b.org" in d.preview() and "email" in d.preview()
    assert ks.confirm(d.id, "ANYTHING").status == "refused"
    # the stored text is unchanged: nothing was rewritten
    assert "a@b.org" in kp.get(les.id).lesson
    kp.edit(les.id, lesson="re-run the failing tests serially first.")
    d2 = _share(env, kp.get(les.id))
    assert not d2.blocked, d2.preview()


def test_a_topic_outside_the_taxonomy_cannot_be_created_by_sharing(env):
    d = _share(env, _lesson(env, topic="my-new-topic"))
    assert d.blocked and "taxonomy" in d.preview()


def test_an_entry_over_the_length_limit_is_blocked(env):
    long = LESSON.replace("believing them.", "believing them. " +
                          "Also check the logs carefully. " * 30)
    d = _share(env, _lesson(env, text=long))
    assert d.blocked and "words" in d.preview()


# ---------------------------------------------------------------------------
# the human-only, two-step confirmation
# ---------------------------------------------------------------------------

def test_public_needs_two_confirmations_with_a_pause(env):
    clock = Clock()
    d = _share(env, _lesson(env), clock=clock)
    first = d.code
    out = ks.confirm(d.id, first, clock=clock)
    assert out.status == "next"
    second = out.dialog
    assert second.stage == 2 and second.code and second.code != first
    assert kn.render_entry(second.entry).strip() in second.preview()
    # typed at once: not an answer
    assert ks.confirm(d.id, second.code, clock=clock).status == "refused"
    clock.t += ks.PAUSE_S + 0.1
    done = ks.confirm(d.id, second.code, clock=clock)
    assert done.status == "done", done.message


def test_a_wrong_code_does_nothing_and_three_expire(env):
    clock = Clock()
    d = _share(env, _lesson(env), clock=clock)
    for _ in range(ks.MAX_WRONG):
        assert ks.confirm(d.id, "WRONG1", clock=clock).status in (
            "refused", "expired")
    assert ks.confirm(d.id, d.code, clock=clock).status == "expired"
    assert list(env["queue"].glob("*.json")) == []


def test_the_code_is_never_written_anywhere(env, tmp_path):
    clock = Clock()
    d = _share(env, _lesson(env), clock=clock)
    code = d.code
    for p in tmp_path.rglob("*"):
        if p.is_file() and ".git" not in p.parts:
            assert code not in p.read_text(errors="replace"), p


def test_a_code_from_another_dialog_does_not_answer(env):
    clock = Clock()
    les = _lesson(env)
    a = _share(env, les, clock=clock)
    b = _share(env, les, clock=clock)
    assert ks.confirm(b.id, a.code, clock=clock).status == "refused"


def test_inside_an_agent_session_nothing_can_be_shared(env, monkeypatch):
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "123")
    with pytest.raises(ks.ShareRefused):
        _share(env, _lesson(env))


# ---------------------------------------------------------------------------
# where a confirmed share lands
# ---------------------------------------------------------------------------

def _confirm_all(d, clock):
    out = ks.confirm(d.id, d.code, clock=clock)
    if out.status == "next":
        clock.t += ks.PAUSE_S + 0.1
        out = ks.confirm(d.id, out.dialog.code, clock=clock)
    return out


def test_a_global_share_lands_in_the_admin_queue_not_in_a_branch(env):
    clock = Clock()
    les = _lesson(env)
    d = _share(env, les, clock=clock)
    out = _confirm_all(d, clock)
    assert out.status == "done"
    subs = list(env["queue"].glob("*.json"))
    assert len(subs) == 1
    sub = json.loads(subs[0].read_text())
    assert sub["status"] == "pending" and sub["author"] == "maxh"
    assert "Serial re-run" in sub["entry"]
    branches = _git(env["delfin"], "branch", "--list")
    assert branches.strip() == "* main"
    assert kp.get(les.id).status == kp.QUEUED


def test_a_project_share_prepares_a_branch_and_touches_no_working_tree(env):
    env["vis"]["value"] = "private"
    clock = Clock()
    proj = env["project"]
    (proj / "dirty.txt").write_text("the person's unstaged work\n")
    before_status = _git(proj, "status", "--porcelain")
    les = _lesson(env)
    d = _share(env, les, layer="project", clock=clock)
    out = _confirm_all(d, clock)
    assert out.status == "done", out.message
    branch = kp.get(les.id).shared_ref
    assert branch.startswith("delfin-knowledge/")
    files = _git(proj, "ls-tree", "-r", "--name-only", branch)
    assert ".delfin/knowledge/serial-re-run-before-believing-failures.md" \
        in files
    assert _git(proj, "status", "--porcelain") == before_status
    assert not (proj / ".delfin").exists()
    assert "git push" in out.message           # offered, not run
    assert _git(proj, "remote") == ""


def test_the_project_layer_of_the_delfin_tree_is_the_global_layer(env):
    d = _share(env, _lesson(env), layer="project", workspace=env["delfin"],
               scope="delfin-development")
    assert d.target.layer == kn.GLOBAL
    assert d.target.scope == "delfin-development"
    assert "scope" in d.preview()


def test_a_near_duplicate_extends_the_existing_entry(env):
    gdir = env["gdir"]
    existing = kn.Entry(
        slug="serial-first", title="Serial re-run before believing failures",
        topic="tools-and-sandbox", scope="always",
        provenance=[kn.Provenance("anna", "2026-09-01", "job:4711")],
        review_by="2027-03-01",
        situation="A parallel pytest run of a large suite reports failures.",
        lesson="Re-run the failing tests serially before believing them.",
        why="Shared state between workers invents failures.")
    kn.write_entry(gdir, existing)
    d = _share(env, _lesson(env))
    assert d.mode == "extend" and d.entry.slug == "serial-first"
    authors = {p.author for p in d.entry.provenance}
    assert authors == {"anna", "maxh"}
    assert "extend" in d.preview()


def test_a_contradiction_is_shown_and_goes_to_the_admin(env):
    gdir = env["gdir"]
    kn.write_entry(gdir, kn.Entry(
        slug="trust-parallel", title="Trust the parallel run",
        topic="tools-and-sandbox", scope="always",
        provenance=[kn.Provenance("anna", "2026-09-01", "job:4711")],
        review_by="2027-03-01",
        situation="A parallel pytest run of a large suite reports failures.",
        lesson="Never re-run the failing tests serially before believing "
               "them.",
        why="Re-runs waste time."))
    clock = Clock()
    d = _share(env, _lesson(env), clock=clock)
    assert "trust-parallel" in d.preview() and "contradict" in d.preview()
    _confirm_all(d, clock)
    sub = json.loads(next(env["queue"].glob("*.json")).read_text())
    assert sub["contradicts"] == ["trust-parallel"]


# ---------------------------------------------------------------------------
# unshare
# ---------------------------------------------------------------------------

def _commit_entry(repo, folder_rel, entry):
    folder = repo / folder_rel
    kn.write_entry(folder, entry)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "entry")


def test_unshare_withdraws_ones_own_entry_as_a_branch(env):
    env["vis"]["value"] = "private"
    proj = env["project"]
    _commit_entry(proj, ".delfin/knowledge", kn.Entry(
        slug="mine", title="My entry title here", topic="workflows",
        scope="project",
        provenance=[kn.Provenance("maxh", "2026-09-01", "job:1")],
        review_by="2027-01-01", situation="s.", lesson="l.", why="w."))
    clock = Clock()
    d = ks.start_unshare("mine", workspace=proj, handle="maxh", clock=clock)
    out = _confirm_all(d, clock)
    assert out.status == "done", out.message
    branch = out.ref
    files = _git(proj, "ls-tree", "-r", "--name-only", branch)
    assert ".delfin/knowledge/mine.md" not in files


def test_unshare_of_someone_elses_entry_is_refused(env):
    proj = env["project"]
    _commit_entry(proj, ".delfin/knowledge", kn.Entry(
        slug="theirs", title="Their entry title", topic="workflows",
        scope="project",
        provenance=[kn.Provenance("anna", "2026-09-01", "job:1")],
        review_by="2027-01-01", situation="s.", lesson="l.", why="w."))
    with pytest.raises(ks.ShareRefused):
        ks.start_unshare("theirs", workspace=proj, handle="maxh")


# ---------------------------------------------------------------------------
# the admin role
# ---------------------------------------------------------------------------

def test_admins_come_from_the_managed_policy_first(env, monkeypatch):
    monkeypatch.setattr(ks, "_settings_admins", lambda: ("alice",))
    monkeypatch.setattr(ks, "_policy_admins", lambda: ("bob",))
    assert ks.admins() == ("bob",)
    monkeypatch.setattr(ks, "_policy_admins", lambda: ())
    assert ks.admins() == ("alice",)


def test_the_admin_identity_is_the_account_not_the_chosen_handle(
        env, monkeypatch):
    monkeypatch.setattr(ks, "admins", lambda: ("bob",))
    monkeypatch.setattr(ks, "_identities", lambda: {"alice"})
    assert not ks.is_admin()
    monkeypatch.setattr(ks, "_identities", lambda: {"bob"})
    assert ks.is_admin()


def test_a_public_entry_needs_a_public_evidence_pointer(env):
    les = kp.propose(LESSON, evidence=["job:812345"], source="agent",
                     topic="tools-and-sandbox",
                     title="Serial re-run before believing failures")
    d = _share(env, les)
    assert d.blocked
    assert any("public pointer" in f for f in d.findings), d.findings
    kp.edit(les.id, evidence="commit:75caa295")
    assert kp.get(les.id).evidence == ["commit:75caa295"]
    d2 = _share(env, kp.get(les.id))
    assert not d2.blocked, d2.preview()


def test_an_edited_evidence_must_still_be_a_pointer(env):
    les = _lesson(env)
    with pytest.raises(ValueError):
        kp.edit(les.id, evidence="see the chat above")


# ---------------------------------------------------------------------------
# the queue record's file mode
# ---------------------------------------------------------------------------

def _queued_record(env):
    clock = Clock()
    out = _confirm_all(_share(env, _lesson(env), clock=clock), clock)
    assert out.status == "done", out.message
    (rec,) = ks.read_queue()
    return rec


@pytest.mark.skipif(not hasattr(os, "fchmod"), reason="POSIX file modes")
@pytest.mark.parametrize("umask", [0o022, 0o002, 0o077],
                         ids=["022", "002", "077"])
def test_a_rewritten_record_keeps_the_mode_it_was_created_with(env, umask):
    """An admin decision rewrites the record. The rewrite went through
    Path.write_text, whose mode is 0o666 minus the umask, so a 0o640
    record in a shared queue became 0o644 or 0o664: readable by every
    account on the host. The mode must not depend on the umask of
    whichever process rewrites."""
    old = os.umask(umask)
    try:
        rec = _queued_record(env)
        p = Path(rec["_path"])
        created = stat.S_IMODE(p.stat().st_mode)
        assert created == 0o640, oct(created)
        rec["status"] = "rejected"
        rec["reason"] = "already in the docs"
        ks.write_queue_record(rec)
    finally:
        os.umask(old)
    assert stat.S_IMODE(p.stat().st_mode) == 0o640
    assert json.loads(p.read_text())["status"] == "rejected"
    assert [x.name for x in p.parent.iterdir()] == [p.name]


@pytest.mark.skipif(not hasattr(os, "fchmod"), reason="POSIX file modes")
def test_group_read_goes_only_to_the_queue_folders_group(env, monkeypatch):
    """A shared queue is read by the admin through the queue folder's
    group. A record whose group is something else -- the author's
    primary group, on a folder without the setgid bit -- must not be
    group-readable: that group is whoever shares the author's primary
    group, not the admin. The record stays 0o600 and the share says why
    the admin cannot read it."""
    other = os.getegid() + 4242
    monkeypatch.setattr(ks, "_queue_gid", lambda _folder: other,
                        raising=False)
    clock = Clock()
    out = _confirm_all(_share(env, _lesson(env), clock=clock), clock)
    assert out.status == "done", out.message
    (rec,) = ks.read_queue()
    p = Path(rec["_path"])
    assert stat.S_IMODE(p.stat().st_mode) == 0o600
    assert "setgid" in out.message


@pytest.mark.skipif(not hasattr(os, "fchmod"), reason="POSIX file modes")
def test_the_per_account_queue_is_never_group_readable(env, monkeypatch):
    monkeypatch.delenv("DELFIN_KNOWLEDGE_QUEUE")
    clock = Clock()
    out = _confirm_all(_share(env, _lesson(env), clock=clock), clock)
    assert out.status == "done", out.message
    (rec,) = ks.read_queue()
    p = Path(rec["_path"])
    assert p.parent == env["home"] / ".delfin" / "knowledge" / "queue"
    assert stat.S_IMODE(p.stat().st_mode) == 0o600
    assert stat.S_IMODE(p.parent.stat().st_mode) == 0o700
