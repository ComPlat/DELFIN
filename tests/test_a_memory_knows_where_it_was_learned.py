"""Provenance on a memory: which branch and commit its body was measured on.

Input: the repository a memory is written in. Output: a ``learned_at:
<branch>@<commit12>`` line in the file's frontmatter.

Semantics. Stamped when a body is written and re-stamped when a body is
replaced, never on recall -- recall says a memory was useful, not that it
was learned again. Absent on every file written before the field existed,
and an absent stamp is never read as a defect.

What it buys: ``memory_tidy`` can name the notes measured on a commit that
never reached the default branch, on a branch that no longer exists.

Alternative considered and rejected: a separate memory store per branch.
It recreates the orphaning this store was just measured to suffer from --
542 stores holding 189 notes, 534 of them keyed to paths that no longer
exist -- and worse, ``git checkout`` changes the key with nothing to
announce it, so the agent would forget mid-session on a branch switch.
Provenance on the note answers the same question and costs one line.
"""

from __future__ import annotations

import subprocess

import pytest

from delfin.agent import memory_store as M
from delfin.agent import memory_tidy as T


def _repo(path, body="base\n"):
    path.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "init", "-q", "-b", "main"], cwd=str(path),
                   check=True)
    for k, v in (("user.email", "a@b.c"), ("user.name", "t")):
        subprocess.run(["git", "config", k, v], cwd=str(path), check=True)
    (path / "f.txt").write_text(body, encoding="utf-8")
    subprocess.run(["git", "add", "-A"], cwd=str(path), check=True)
    subprocess.run(["git", "commit", "-q", "-m", "base"], cwd=str(path),
                   check=True)
    return path


def _short(path, rev="HEAD"):
    return subprocess.run(["git", "-C", str(path), "rev-parse", "--short=12",
                           rev], capture_output=True, text=True).stdout.strip()


# -- reading the stamp -----------------------------------------------------

def test_the_stamp_names_the_branch_and_the_commit(tmp_path):
    repo = _repo(tmp_path / "r")
    assert M._head_ref(repo) == f"main@{_short(repo)}"


def test_a_detached_head_says_so_rather_than_inventing_a_branch(tmp_path):
    repo = _repo(tmp_path / "r")
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "--detach"],
                   check=True)
    assert M._head_ref(repo) == f"detached@{_short(repo)}"


def test_outside_a_repository_there_is_no_stamp(tmp_path):
    plain = tmp_path / "plain"
    plain.mkdir()
    assert M._head_ref(plain) == ""


def test_a_repository_with_no_commit_has_no_stamp(tmp_path):
    empty = tmp_path / "e"
    empty.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=str(empty), check=True)
    assert M._head_ref(empty) == ""


def test_git_failing_gives_no_stamp_rather_than_raising(tmp_path,
                                                        monkeypatch):
    def _boom(*a, **k):
        raise OSError("no git here")
    monkeypatch.setattr(M.subprocess, "run", _boom)
    assert M._head_ref(tmp_path) == ""


# -- the field survives every rewrite --------------------------------------
#
# `source` is a composer-owned field and was erased by a single recall
# until it was carried forward explicitly. `learned_at` is the same shape
# of field, so it gets the same test before it can repeat the defect.

def test_the_composer_writes_the_line():
    out = M._compose_frontmatter(
        name="n", description="d", created_at=1, updated_at=1, use_count=1,
        memory_type="project", body="b", learned_at="main@abc123abc123")
    assert "learned_at: main@abc123abc123" in out


def test_no_stamp_writes_no_line(tmp_path):
    out = M._compose_frontmatter(
        name="n", description="d", created_at=1, updated_at=1, use_count=1,
        memory_type="project", body="b")
    assert "learned_at" not in out


def test_the_field_is_composer_owned_not_carried_as_an_extra():
    """In `extras` it would be written twice on every rewrite."""
    assert "learned_at" in M._KNOWN_FRONT_FIELDS


def test_every_rewrite_carries_the_stamp_forward():
    """Each call site must pass it. A site that forgets erases provenance
    the next time that path runs -- which is exactly how `source` was
    lost once already."""
    import inspect
    import re

    src = inspect.getsource(M)
    # Call sites only -- the `def` line is not one of them, and counting it
    # made this assertion off by one the first time it ran.
    composes = len(re.findall(r"=\s*_compose_frontmatter\(", src))
    carries = src.count("learned_at=")
    assert composes >= 4, "the call sites moved; re-read them"
    assert carries >= composes, (
        f"{composes} call sites, {carries} pass learned_at")


# -- what tidy does with it ------------------------------------------------

@pytest.fixture()
def store(tmp_path):
    d = tmp_path / "memory"
    d.mkdir(parents=True)
    return d


def _note(store, name, stamp, *, body="a body", source="agent"):
    stamp_line = f"learned_at: {stamp}\n" if stamp else ""
    (store / f"project_{name}.md").write_text(
        "---\n"
        f"name: {name}\ndescription: d\ncreated_at: 1\nupdated_at: 9999999999\n"
        f"use_count: 1\nsource: {source}\n{stamp_line}"
        "metadata:\n  type: project\n---\n\n"
        f"{body}\n", encoding="utf-8")


def test_a_note_from_a_deleted_branch_that_never_landed_is_reported(
        tmp_path, store):
    repo = _repo(tmp_path / "r")
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "-b", "side"],
                   check=True)
    (repo / "g.txt").write_text("x", encoding="utf-8")
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "w"],
                   check=True)
    side = _short(repo)
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "main"],
                   check=True)
    subprocess.run(["git", "-C", str(repo), "branch", "-D", "-q", "side"],
                   check=True)
    _note(store, "orphan", f"side@{side}")
    p = T.propose(store, repo_root=repo)
    assert [n.name for n in p.unlanded] == ["orphan"]


def test_a_note_from_the_default_branch_is_not_reported(tmp_path, store):
    repo = _repo(tmp_path / "r")
    _note(store, "landed", f"main@{_short(repo)}")
    assert T.propose(store, repo_root=repo).unlanded == []


def test_a_note_whose_branch_still_exists_is_not_reported(tmp_path, store):
    """Unlanded work in progress is not abandoned work."""
    repo = _repo(tmp_path / "r")
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "-b", "wip"],
                   check=True)
    (repo / "g.txt").write_text("x", encoding="utf-8")
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "w"],
                   check=True)
    _note(store, "wip", f"wip@{_short(repo)}")
    assert T.propose(store, repo_root=repo).unlanded == []


def test_a_note_with_no_stamp_is_never_reported(tmp_path, store):
    repo = _repo(tmp_path / "r")
    _note(store, "legacy", "")
    assert T.propose(store, repo_root=repo).unlanded == []


def test_a_user_written_note_is_never_reported(tmp_path, store):
    """The user's own words are not judged by where the code was."""
    repo = _repo(tmp_path / "r")
    _note(store, "theirs", "gone@000000000000", source="user")
    assert T.propose(store, repo_root=repo).unlanded == []


def test_without_a_repository_the_proposal_is_what_it_always_was(store):
    _note(store, "x", "gone@000000000000")
    assert T.propose(store).unlanded == []


def test_a_repository_with_no_default_branch_judges_nothing(tmp_path, store):
    """A wrong default branch would rule every note unlanded at once, so
    an unresolvable one rules none."""
    repo = tmp_path / "nobranch"
    repo.mkdir()
    subprocess.run(["git", "init", "-q"], cwd=str(repo), check=True)
    _note(store, "x", "gone@000000000000")
    assert T.propose(store, repo_root=repo).unlanded == []


# -- the report is a question, not an action -------------------------------

def test_the_finding_is_never_proposed_for_retirement(tmp_path, store):
    """Measured: a branch merged with --squash and deleted leaves its
    commit unreachable, so work that DID land reads as work that did not.
    Whether a project squashes is not knowable here, so this list is
    reported and `apply` does not touch it."""
    repo = _repo(tmp_path / "r")
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "-b", "side"],
                   check=True)
    (repo / "g.txt").write_text("x", encoding="utf-8")
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "w"],
                   check=True)
    side = _short(repo)
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "main"],
                   check=True)
    subprocess.run(["git", "-C", str(repo), "merge", "--squash", "-q", "side"],
                   check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "squashed"],
                   check=True)
    subprocess.run(["git", "-C", str(repo), "branch", "-D", "-q", "side"],
                   check=True)
    _note(store, "squashed", f"side@{side}")
    p = T.propose(store, repo_root=repo)
    # It IS reported -- the check cannot tell this case apart ...
    assert [n.name for n in p.unlanded] == ["squashed"]
    # ... and precisely because it cannot, nothing acts on it.
    assert p.retire == []
    before = (store / "project_squashed.md").read_text(encoding="utf-8")
    T.apply(p)
    assert (store / "project_squashed.md").read_text(encoding="utf-8") == before


def test_the_report_says_why_it_is_only_a_report(tmp_path, store):
    repo = _repo(tmp_path / "r")
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "-b", "side"],
                   check=True)
    (repo / "g.txt").write_text("x", encoding="utf-8")
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "w"],
                   check=True)
    side = _short(repo)
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "main"],
                   check=True)
    subprocess.run(["git", "-C", str(repo), "branch", "-D", "-q", "side"],
                   check=True)
    _note(store, "orphan", f"side@{side}")
    text = T.propose(store, repo_root=repo).render()
    assert "squash" in text
    assert "not proposed for anything" in text


def test_the_count_line_does_not_count_the_report(tmp_path, store):
    """`after` is what applying would leave. Reporting leaves everything."""
    repo = _repo(tmp_path / "r")
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "-b", "side"],
                   check=True)
    (repo / "g.txt").write_text("x", encoding="utf-8")
    subprocess.run(["git", "-C", str(repo), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(repo), "commit", "-q", "-m", "w"],
                   check=True)
    side = _short(repo)
    subprocess.run(["git", "-C", str(repo), "checkout", "-q", "main"],
                   check=True)
    subprocess.run(["git", "-C", str(repo), "branch", "-D", "-q", "side"],
                   check=True)
    _note(store, "orphan", f"side@{side}")
    p = T.propose(store, repo_root=repo)
    assert p.after == p.before == 1


# -- the callers pass the repository ---------------------------------------

def test_the_cli_and_the_dashboard_both_pass_it():
    """A check nothing calls with a repository reports nothing, forever."""
    import pathlib

    base = pathlib.Path(M.__file__).resolve().parent
    cli = (base / "cli.py").read_text(encoding="utf-8")
    assert "propose(store, repo_root=root)" in cli
    assert "repo_root=Path(workspace)" in cli
    tab = (base.parent / "dashboard" / "tab_agent.py").read_text(
        encoding="utf-8")
    assert "propose(store, repo_root=root)" in tab
