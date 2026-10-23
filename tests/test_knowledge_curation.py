"""Curation of the global knowledge: lint, index, duplicates, hygiene, and
the admin-only write actions. Temporary git repositories only.
"""
from __future__ import annotations

import datetime as dt
import json
import subprocess

import pytest

from delfin.agent import knowledge as kn
from delfin.agent import knowledge_curation as kc
from delfin.agent import knowledge_private as kp
from delfin.agent import knowledge_share as ks
from delfin.agent import repo_visibility as rv


class Clock:
    def __init__(self):
        self.t = 5000.0

    def __call__(self):
        return self.t


def _git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


def _e(slug, **kw):
    base = dict(slug=slug, title=f"Title of {slug}", topic="workflows",
                scope="always",
                provenance=[kn.Provenance("anna", "2026-09-01", "test:tests/test_x.py::t")],
                review_by="2027-06-01", situation=f"Situation for {slug}.",
                lesson=f"Lesson for {slug}.", why="Because.")
    base.update(kw)
    return kn.Entry(**base)


@pytest.fixture(autouse=True)
def env(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    monkeypatch.setenv("DELFIN_KNOWLEDGE_QUEUE", str(tmp_path / "queue"))
    for k in ("GIT_AUTHOR_NAME", "GIT_COMMITTER_NAME"):
        monkeypatch.setenv(k, "tester")
    for k in ("GIT_AUTHOR_EMAIL", "GIT_COMMITTER_EMAIL"):
        monkeypatch.setenv(k, "tester@example.invalid")
    repo = tmp_path / "delfin-src"
    repo.mkdir()
    _git(repo, "init", "-q", "-b", "main")
    gdir = repo / "delfin" / "agent" / "pack" / "knowledge"
    gdir.mkdir(parents=True)
    (repo / "delfin" / "agent" / "pack" / "manifest.yaml").write_text("x")
    (repo / "delfin" / "agent" / "tools.py").write_text(
        "def existing_function():\n    pass\n")
    kn.write_topics(gdir, {t: "" for t in kn.DEFAULT_TOPICS})
    a = _e("serial-first", topic="tools-and-sandbox",
           title="Re-run parallel failures serially",
           situation="A parallel pytest run of the suite reports failures.",
           lesson="Re-run the failing tests serially before believing them.",
           why="Workers share state.")
    b = _e("serial-again", topic="tools-and-sandbox",
           title="Serial re-run of parallel failures",
           situation="A parallel pytest run of the suite reports failures.",
           lesson="Re-run the failing tests serially before believing them "
                  "as real.",
           why="Workers share state between tests.",
           provenance=[kn.Provenance("bert", "2026-09-05", "commit:abc1234")])
    stale = _e("old-cluster-note", topic="cluster-hpc",
               review_by="2026-01-01",
               lesson="Use `delfin/agent/gone.py` and `vanished_helper()`.")
    fresh = _e("tools-note", topic="tools-and-sandbox",
               lesson="Call `existing_function()` in "
                      "`delfin/agent/tools.py`.")
    for e in (a, b, stale, fresh):
        kn.write_entry(gdir, e)
    kn.write_index(gdir)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "seed")
    monkeypatch.setattr(ks, "_DELFIN_ROOT", repo)
    monkeypatch.setattr(kn, "global_dir", lambda: gdir)
    monkeypatch.setattr(rv, "_RUNNER", lambda _r: "public")
    monkeypatch.setattr(ks, "admins", lambda: ("boss",))
    who = {"ids": {"someone"}}
    monkeypatch.setattr(ks, "_identities", lambda: who["ids"])
    return {"repo": repo, "gdir": gdir, "who": who, "tmp": tmp_path}


def _as_admin(env):
    env["who"]["ids"] = {"boss"}


def _run(dialog, clock):
    out = ks.confirm(dialog.id, dialog.code, clock=clock)
    if out.status == "next":
        clock.t += ks.PAUSE_S + 0.1
        out = ks.confirm(dialog.id, out.dialog.code, clock=clock)
    return out


# ---------------------------------------------------------------------------
# read-only tools: anyone
# ---------------------------------------------------------------------------

def test_lint_passes_a_clean_tree_and_names_each_problem(env):
    gdir = env["gdir"]
    assert kc.lint(gdir) == []
    kn.write_entry(gdir, _e("bad", topic="astrology", title=""))
    problems = kc.lint(gdir)
    assert any("bad" in p and "taxonomy" in p for p in problems)
    assert any("bad" in p and "title" in p for p in problems)
    assert any("INDEX.md" in p for p in problems)      # out of sync now


def test_lint_applies_the_public_content_rules(env):
    gdir = env["gdir"]
    kn.write_entry(gdir, _e("leak", lesson="Mail a@b.org for access."))
    kn.write_index(gdir)
    assert any("leak" in p and "email" in p for p in kc.lint(gdir))


def test_lint_flags_exact_duplicates(env):
    gdir = env["gdir"]
    twin = kn.parse_entry((gdir / "tools-note.md").read_text(), slug="twin")
    kn.write_entry(gdir, twin)
    kn.write_index(gdir)
    assert any("duplicate" in p and "twin" in p for p in kc.lint(gdir))


def test_dupes_lists_near_duplicates_and_contradictions(env):
    gdir = env["gdir"]
    kn.write_entry(gdir, _e(
        "never-serial", topic="tools-and-sandbox",
        title="Do not re-run serially",
        situation="A parallel pytest run of the suite reports failures.",
        lesson="Never re-run the failing tests serially before believing "
               "them.", why="Wasted time."))
    found = kc.dupes(gdir)
    pairs = {(f.kind, frozenset((f.a, f.b))) for f in found}
    assert ("near-duplicate",
            frozenset({"serial-first", "serial-again"})) in pairs
    assert any(k == "contradiction" and "never-serial" in p
               for k, p in pairs)


def test_stats_lists_expired_never_loaded_and_broken_references(env):
    usage = env["tmp"] / "usage.json"
    usage.write_text(json.dumps({"serial-first": {"loads": 3, "last": 1}}))
    st = kc.stats(env["gdir"], usage_path=usage, repo_root=env["repo"],
                  today=dt.date(2026, 10, 15))
    assert st["expired"] == ["old-cluster-note"]
    assert "serial-again" in st["never_loaded"]
    assert "serial-first" not in st["never_loaded"]
    broken = {(b.slug, b.ref) for b in st["broken"]}
    assert ("old-cluster-note", "delfin/agent/gone.py") in broken
    assert ("old-cluster-note", "vanished_helper()") in broken
    assert not [b for b in st["broken"] if b.slug == "tools-note"]
    assert st["by_topic"]["tools-and-sandbox"] == 3


# ---------------------------------------------------------------------------
# write actions: admin only, confirmed by the admin, a branch never a push
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("make", [
    lambda: kc.plan_merge("serial-first", "serial-again"),
    lambda: kc.plan_move("serial-first", "workflows"),
    lambda: kc.plan_retire("old-cluster-note", reason="site retired"),
    lambda: kc.plan_rename("serial-first", "serial-rerun"),
    lambda: kc.plan_index(),
    lambda: kc.plan_topic_add("data-formats", "file formats"),
])
def test_a_non_admin_cannot_start_a_curation_action(env, make):
    with pytest.raises(ks.ShareRefused):
        kc.start(make())


def test_an_agent_session_cannot_curate_even_as_admin(env, monkeypatch):
    _as_admin(env)
    monkeypatch.setenv("DELFIN_AGENT_SESSION", "1")
    with pytest.raises(ks.ShareRefused):
        kc.start(kc.plan_retire("old-cluster-note", reason="x"))


def test_merge_keeps_both_provenances_and_regenerates_the_index(env):
    _as_admin(env)
    clock = Clock()
    plan = kc.plan_merge("serial-first", "serial-again")
    assert "serial-again.md" in plan.diff and "bert" in plan.diff
    d = kc.start(plan, clock=clock)
    assert d.stage == 1 and d.public
    out = _run(d, clock)
    assert out.status == "done", out.message
    repo, branch = env["repo"], out.ref
    files = _git(repo, "ls-tree", "-r", "--name-only", branch)
    assert "serial-again.md" not in files
    merged = _git(repo, "show",
                  f"{branch}:delfin/agent/pack/knowledge/serial-first.md")
    assert "anna" in merged and "bert" in merged
    index = _git(repo, "show",
                 f"{branch}:delfin/agent/pack/knowledge/INDEX.md")
    assert "serial-again" not in index
    # nothing in the working tree moved
    assert (env["gdir"] / "serial-again.md").exists()
    assert _git(repo, "status", "--porcelain") == ""


def test_retire_and_move_and_rename(env):
    _as_admin(env)
    repo = env["repo"]
    for plan, check in [
        (kc.plan_retire("old-cluster-note", reason="site retired"),
         lambda b: "old-cluster-note.md" not in _git(
             repo, "ls-tree", "-r", "--name-only", b)),
        (kc.plan_move("tools-note", "workflows"),
         lambda b: "topic: workflows" in _git(
             repo, "show", f"{b}:delfin/agent/pack/knowledge/tools-note.md")),
        (kc.plan_rename("tools-note", "tools-call-note"),
         lambda b: "tools-call-note.md" in _git(
             repo, "ls-tree", "-r", "--name-only", b)),
    ]:
        clock = Clock()
        out = _run(kc.start(plan, clock=clock), clock)
        assert out.status == "done", out.message
        assert check(out.ref)


def test_a_move_to_an_unknown_topic_is_refused(env):
    with pytest.raises(ValueError):
        kc.plan_move("tools-note", "astrology")


def test_split_writes_two_entries_with_the_provenance(env):
    _as_admin(env)
    parts = [
        dict(slug="serial-a", title="Serial re-run part A",
             situation="Parallel failures.", lesson="Re-run serially.",
             why="Shared state."),
        dict(slug="serial-b", title="Serial re-run part B",
             situation="Flaky marker.", lesson="Mark it flaky.",
             why="Known flake."),
    ]
    plan = kc.plan_split("serial-first", parts)
    clock = Clock()
    out = _run(kc.start(plan, clock=clock), clock)
    files = _git(env["repo"], "ls-tree", "-r", "--name-only", out.ref)
    assert "serial-a.md" in files and "serial-b.md" in files
    assert "serial-first.md" not in files
    a = _git(env["repo"], "show",
             f"{out.ref}:delfin/agent/pack/knowledge/serial-a.md")
    assert "anna" in a


# ---------------------------------------------------------------------------
# the admin review queue
# ---------------------------------------------------------------------------

LESSON = ("Situation: a cluster job is submitted from a login node. "
          "Lesson: load the module environment inside the job script. "
          "Why: the job does not inherit the login shell's modules.")


def _queued(env):
    les = kp.propose(LESSON, evidence=["doc:docs/cluster.md"],
                     source="agent",
                     topic="cluster-hpc", title="Load modules inside the "
                     "job script")
    clock = Clock()
    d = ks.start_share(les.id, "global", workspace=env["repo"],
                       handle="carl", clock=clock)
    assert not d.blocked, d.preview()
    out = _run(d, clock)
    assert out.status == "done"
    return les, out.ref


def test_accept_prepares_the_branch_and_tells_the_author(env):
    les, qid = _queued(env)
    _as_admin(env)
    clock = Clock()
    out = _run(kc.start(kc.plan_accept(qid), clock=clock), clock)
    assert out.status == "done", out.message
    files = _git(env["repo"], "ls-tree", "-r", "--name-only", out.ref)
    assert "delfin/agent/pack/knowledge/load-modules-inside-the-job-script" \
        ".md" in files
    assert ks.sync_own_lessons()
    assert kp.get(les.id).status == kp.PREPARED


def test_reject_records_the_reason_the_author_sees(env):
    les, qid = _queued(env)
    _as_admin(env)
    clock = Clock()
    out = _run(kc.start(kc.plan_reject(qid, "already in the docs"),
                        clock=clock), clock)
    assert out.status == "done", out.message
    ks.sync_own_lessons()
    got = kp.get(les.id)
    assert got.status == kp.REJECTED and got.note == "already in the docs"


def test_queue_merge_into_an_existing_entry(env):
    les, qid = _queued(env)
    _as_admin(env)
    clock = Clock()
    out = _run(kc.start(kc.plan_queue_merge(qid, "old-cluster-note"),
                        clock=clock), clock)
    text = _git(env["repo"], "show", f"{out.ref}:delfin/agent/pack/"
                "knowledge/old-cluster-note.md")
    assert "carl" in text and "anna" in text


def test_queue_edit_rechecks_the_admins_wording(env):
    les, qid = _queued(env)
    _as_admin(env)
    with pytest.raises(ValueError):
        kc.plan_queue_edit(qid, lesson="Ask a@b.org to load the modules.")
    plan = kc.plan_queue_edit(qid, lesson="Load the modules in the script.")
    clock = Clock()
    out = _run(kc.start(plan, clock=clock), clock)
    assert out.status == "done"


def test_only_the_admin_adds_a_topic(env):
    _as_admin(env)
    clock = Clock()
    out = _run(kc.start(kc.plan_topic_add("data-formats", "file formats"),
                        clock=clock), clock)
    text = _git(env["repo"], "show",
                f"{out.ref}:delfin/agent/pack/knowledge/TOPICS.md")
    assert "data-formats" in text


def test_lint_reads_the_provenance_too(env):
    gdir = env["gdir"]
    kn.write_entry(gdir, _e("local-pointer", provenance=[
        kn.Provenance("anna", "2026-09-01", "job:812345")]))
    kn.write_index(gdir)
    problems = kc.lint(gdir)
    assert any("local-pointer" in p and "public pointer" in p
               for p in problems), problems


def test_an_admin_plan_whose_result_fails_the_check_is_blocked(env):
    _as_admin(env)
    parts = [dict(slug="part-a", title="Part A of the lesson",
                  situation="Parallel failures.",
                  lesson="Ask a@b.org before re-running.", why="Shared.")]
    d = kc.start(kc.plan_split("serial-first", parts), clock=Clock())
    assert d.blocked and "email" in d.preview()
    assert "check: passed" not in d.preview()
