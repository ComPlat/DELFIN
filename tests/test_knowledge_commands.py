"""The person's commands: /memory, /share, /unshare, /lessons, /swarm,
/knowledge -- driven as the REPL and the dashboard drive them.
"""
from __future__ import annotations

import re
import subprocess

import pytest

from delfin.agent import knowledge as kn
from delfin.agent import knowledge_commands as kcmd
from delfin.agent import knowledge_private as kp
from delfin.agent import knowledge_share as ks
from delfin.agent import repo_visibility as rv
from delfin.agent import swarm_memory as sm

LESSON = ("Situation: a parallel pytest run of a large suite reports "
          "failures. Lesson: re-run the failing tests serially before "
          "believing them. Why: shared state between workers.")


class Clock:
    def __init__(self):
        self.t = 100.0

    def __call__(self):
        return self.t


def _git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True).stdout.strip()


@pytest.fixture(autouse=True)
def env(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    monkeypatch.delenv("DELFIN_SWARM_RUN", raising=False)
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
    kn.write_topics(gdir, {t: "" for t in kn.DEFAULT_TOPICS})
    kn.write_entry(gdir, kn.Entry(
        slug="scratch-purge", title="Scratch is purged after two weeks",
        topic="cluster-hpc", scope="always",
        provenance=[kn.Provenance("anna", "2026-09-01",
                                  "doc:docs/cluster.md")],
        review_by="2027-06-01", situation="Results on scratch.",
        lesson="Copy them home within two weeks.", why="Site purge."))
    kn.write_index(gdir)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "seed")
    monkeypatch.setattr(ks, "_DELFIN_ROOT", repo)
    monkeypatch.setattr(kn, "global_dir", lambda: gdir)
    monkeypatch.setattr(rv, "_RUNNER", lambda _r: "unknown")
    clock = Clock()
    monkeypatch.setattr(kcmd, "_CLOCK", clock)
    ws = tmp_path / "ws"
    ws.mkdir()
    return {"repo": repo, "gdir": gdir, "ws": ws, "clock": clock,
            "queue": tmp_path / "queue"}


class Person:
    """A human at the terminal: reads the screen, types the code shown."""

    def __init__(self, clock, *, wait=True):
        self.screen: list = []
        self.clock = clock
        self.wait = wait

    def say(self, text):
        self.screen.append(str(text))

    def ask(self, _prompt):
        codes = re.findall(r"code ([A-Z0-9]{3}-[A-Z0-9]{3})",
                           "\n".join(self.screen))
        if self.wait:
            self.clock.t += ks.PAUSE_S + 0.5     # reads the text first
        return codes[-1] if codes else ""


def _cmd(env, name, args, person=None):
    return kcmd.handle(name, args, workspace=env["ws"],
                       ask=person.ask if person else None,
                       say=person.say if person else None)


def test_memory_lists_every_layer_with_its_label(env, monkeypatch):
    les = kp.propose(LESSON, evidence=["test:t"], source="agent",
                     topic="tools-and-sandbox", title="Serial re-run")
    monkeypatch.setenv("DELFIN_SWARM_RUN", "wave-1")
    sm.append("wave-1", author="agent-a", text="test X is flaky",
              evidence="test:x")
    out = _cmd(env, "memory", "")
    assert f"[private · proposal {les.id} · would be public]" in out
    assert f"/share {les.id}" in out
    assert "[swarm wave-1 #1] agent-a: test X is flaky" in out
    assert "[global · public · cluster-hpc] scratch-purge" in out
    assert "anna 2026-09-01" in out                   # provenance
    only = _cmd(env, "memory", "swarm")
    assert "swarm" in only and "scratch-purge" not in only
    assert "Serial re-run" not in only


def test_share_in_the_repl_reads_the_code_from_the_person(env):
    les = kp.propose(LESSON, evidence=["test:t"], source="agent",
                     topic="tools-and-sandbox", title="Serial re-run")
    p = Person(env["clock"])
    out = _cmd(env, "share", f"{les.id} global as maxh", p)
    assert "queued" in out, (out, p.screen)
    screen = "\n".join(p.screen)
    assert ks.PUBLIC_BANNER in screen
    assert "step 1 of 2" in screen and "step 2 of 2" in screen
    assert kp.get(les.id).status == kp.QUEUED


def test_typing_too_fast_on_the_second_step_is_not_an_answer(env):
    les = kp.propose(LESSON, evidence=["test:t"], source="agent",
                     topic="tools-and-sandbox", title="Serial re-run")
    p = Person(env["clock"], wait=False)
    out = _cmd(env, "share", f"{les.id} global as maxh", p)
    assert "queued" not in out
    assert "read the full text first" in "\n".join(p.screen)
    assert kp.get(les.id).status == kp.PRIVATE


def test_share_in_the_dashboard_is_two_typed_commands(env):
    les = kp.propose(LESSON, evidence=["test:t"], source="agent",
                     topic="tools-and-sandbox", title="Serial re-run")
    p = Person(env["clock"])
    out = _cmd(env, "share", f"{les.id} global as maxh")
    m = re.search(r"/share confirm (\w+) <code>", out)
    assert m, out
    dialog = ks.get_dialog(m.group(1))
    out2 = kcmd.handle("share", f"confirm {dialog.id} {dialog.code}",
                       workspace=env["ws"])
    assert "step 2 of 2" in out2
    env["clock"].t += ks.PAUSE_S + 1
    out3 = kcmd.handle("share", f"confirm {dialog.id} {dialog.code}",
                       workspace=env["ws"])
    assert "queued" in out3
    del p


def test_a_blocked_share_says_what_to_reword(env):
    les = kp.propose(LESSON.replace("believing them", "asking a@b.org"),
                     evidence=["test:t"], source="agent",
                     topic="tools-and-sandbox", title="Serial re-run")
    p = Person(env["clock"])
    out = _cmd(env, "share", f"{les.id} global as maxh", p)
    assert "blocked" in out
    assert "email" in "\n".join(p.screen)
    _cmd(env, "lessons", f"edit {les.id} lesson=re-run the failing tests "
         "serially first.")
    out = _cmd(env, "share", f"{les.id} global as maxh", Person(env["clock"]))
    assert "queued" in out


def test_knowledge_read_views_are_for_everyone_writes_for_the_admin(
        env, monkeypatch):
    assert "scratch-purge" in _cmd(env, "knowledge", "list")
    assert "lint: clean" in _cmd(env, "knowledge", "lint")
    assert "entries: 1" in _cmd(env, "knowledge", "stats")
    monkeypatch.setattr(ks, "admins", lambda: ("boss",))
    monkeypatch.setattr(ks, "_identities", lambda: {"someone"})
    out = _cmd(env, "knowledge", "retire scratch-purge site gone",
               Person(env["clock"]))
    assert out.startswith("refused") and "admin" in out
    monkeypatch.setattr(ks, "_identities", lambda: {"boss"})
    out = _cmd(env, "knowledge", "retire scratch-purge site gone",
               Person(env["clock"]))
    assert "prepared branch" in out, out


def test_swarm_harvest_and_propose(env):
    sm.append("wave-2", author="agent-a",
              text="module repo_map exists already",
              evidence="commit:75caa295")
    sm.close("wave-2")
    out = _cmd(env, "swarm", "harvest wave-2")
    assert "#1 agent-a" in out and "evidence commit:75caa295" in out
    out = _cmd(env, "swarm", "propose wave-2 1")
    assert "private lesson" in out
    assert kp.list_lessons()[0].source == "swarm"
