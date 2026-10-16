"""Roll back 1 commit asks first, names both commits, and only a click answers.

Input: a checkout made here (main pushed to a bare origin, then commits on
top: one origin has, one only local, a merge), :mod:`delfin.dashboard.
rollback` with a hand-driven clock, the header's ``BranchPanel`` with the
confirmation under its More row, and the dashboard's source.

Output:

* the plan names the commit that leaves and the one HEAD moves to (short
  sha and subject), counts the commits a merge takes with it, and says how
  it is undone (``git reset --hard <old>``, the reflog, Pull when
  origin/main has it); it is refused with changes in the tree (untracked
  included) and on a root commit;
* the confirm button is disabled when the question opens and a click
  within the hold changes nothing; Cancel and a new question change
  nothing; a click after the hold resets to the planned commit, and the
  undo command restores it;
* the reset is refused when HEAD moved after the question was asked;
* Esc presses Cancel, and no key handler presses the confirm button;
* the confirmation lives under More and is withdrawn when More closes;
* the dashboard's handler asks, and resets only from the confirmation.
"""

from __future__ import annotations

import os
import subprocess

import pytest

pytest.importorskip("ipywidgets")

from delfin.dashboard import branch_list as BL  # noqa: E402
from delfin.dashboard import header_ui as HU  # noqa: E402
from delfin.dashboard import rollback as RB  # noqa: E402

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args) -> str:
    return subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=True, capture_output=True,
        text=True).stdout.strip()


def _commit(cwd, name, msg=None):
    (cwd / name).write_text(name + "\n")
    _git(cwd, "add", name)
    _git(cwd, "commit", "-q", "-m", msg or f"add {name}")


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", os.devnull)
    monkeypatch.setenv("GIT_CONFIG_SYSTEM", os.devnull)
    monkeypatch.setenv(BL.GH_ENV, "off")
    (tmp_path / "home").mkdir()
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    co = tmp_path / "checkout"
    _git(tmp_path, "clone", "-q", str(origin), str(co))
    _git(co, "checkout", "-q", "-b", "main")
    _commit(co, "a.txt", "first")
    _commit(co, "b.txt", "second, on origin")
    _git(co, "push", "-q", "-u", "origin", "main")
    return co


class _Clock:
    def __init__(self):
        self.t = 1000.0

    def __call__(self):
        return self.t


def _head(co, short=False):
    return _git(co, "rev-parse", *(["--short"] if short else []), "HEAD")


# -- the plan -------------------------------------------------------------------

def test_the_plan_names_both_commits_and_how_it_is_undone(repo):
    p, why = RB.plan(repo)
    assert why == "" and p is not None
    assert p.branch == "main" and p.leaving == 1
    assert p.head == _head(repo) and p.head_subject == "second, on origin"
    assert p.target == _git(repo, "rev-parse", "HEAD~1") and p.target_subject == "first"
    text = RB.describe(p)
    for s in (p.head_short, "second, on origin", p.target_short, "first",
              f"git reset --hard {p.target_short}", f"git reset --hard {p.head_short}",
              "git reflog", "Pull"):
        assert s in text, s


def test_a_local_only_commit_is_not_promised_back_by_pull(repo):
    _commit(repo, "c.txt", "third, only here")
    p, _ = RB.plan(repo)
    assert p.on_remote == () and "Pull" not in RB.describe(p)
    assert f"git reset --hard {p.head_short}" in RB.describe(p)


def test_a_merge_says_how_many_commits_leave_with_it(repo):
    _git(repo, "checkout", "-q", "-b", "side")
    _commit(repo, "s1.txt")
    _commit(repo, "s2.txt")
    _git(repo, "checkout", "-q", "main")
    _commit(repo, "m.txt")
    _git(repo, "merge", "-q", "--no-ff", "-m", "merge side", "side")
    p, _ = RB.plan(repo)
    assert p.head_subject == "merge side" and p.target_subject == "add m.txt"
    assert p.leaving == 3                      # the merge and the two side commits
    assert "3 commits leave main" in RB.describe(p)


@pytest.mark.parametrize("dirt", ["tracked", "untracked"])
def test_changes_in_the_tree_refuse_before_anything_is_asked(repo, dirt):
    if dirt == "tracked":
        (repo / "a.txt").write_text("changed\n")
    else:
        (repo / "new.txt").write_text("new\n")
    p, why = RB.plan(repo)
    assert p is None and "uncommitted changes" in why


def test_a_root_commit_has_nothing_to_roll_back_to(tmp_path, monkeypatch):
    monkeypatch.setenv("GIT_CONFIG_GLOBAL", os.devnull)
    d = tmp_path / "r"
    d.mkdir()
    _git(d, "init", "-q", "-b", "main")
    _commit(d, "a.txt")
    p, why = RB.plan(d)
    assert p is None and "no earlier commit" in why


# -- the confirmation ------------------------------------------------------------

def _ask(repo, clock):
    c = RB.Confirm(now=clock, timer=False)
    p, _ = RB.plan(repo)
    done = []

    def _yes():
        done.append(RB.execute(repo, p))
    c.ask(RB.describe(p), f"Roll back to {p.target_short}", _yes)
    return c, p, done


def test_the_confirm_is_held_and_a_click_inside_the_hold_is_not_an_answer(repo):
    clock = _Clock()
    before = _head(repo)
    c, p, done = _ask(repo, clock)
    assert c.is_open and c.yes.disabled
    assert c.yes.description == f"Roll back to {p.target_short}"
    assert p.head_short in c.text.value and p.target_short in c.text.value
    clock.t += RB.HOLD_S / 2
    c.yes.click()                              # a double click, a key repeat
    assert done == [] and _head(repo) == before and c.is_open
    c.release()
    assert not c.yes.disabled


def test_cancel_and_a_new_question_change_nothing(repo):
    clock = _Clock()
    before = _head(repo)
    c, p, done = _ask(repo, clock)
    clock.t += 5
    c.no.click()
    assert not c.is_open and done == [] and _head(repo) == before
    c.yes.click()                              # a stale click after Cancel
    assert done == [] and _head(repo) == before
    c, p, done = _ask(repo, clock)
    clock.t += 5
    c.ask(RB.describe(p), "Roll back", lambda: done.append("second"))
    c.yes.click()                              # the hold starts again
    assert done == [] and _head(repo) == before and c.yes.disabled


def test_a_click_after_the_hold_resets_and_the_undo_restores(repo):
    clock = _Clock()
    before = _head(repo)
    c, p, done = _ask(repo, clock)
    clock.t += RB.HOLD_S + 0.01
    c.yes.click()
    assert not c.is_open
    ok, text, details = done[0]
    assert ok and _head(repo) == p.target
    assert f"git reset --hard {p.head_short}" in text
    assert _git(repo, "status", "--porcelain") == ""
    _git(repo, *RB.undo_command(p).split()[1:])
    assert _head(repo) == before


def test_the_reset_is_refused_when_head_moved_after_the_question(repo):
    p, _ = RB.plan(repo)
    _commit(repo, "c.txt", "arrived meanwhile")
    moved = _head(repo)
    ok, text, _d = RB.execute(repo, p)
    assert not ok and "moved" in text and _head(repo) == moved


def test_esc_presses_cancel_and_no_key_presses_the_confirm_button():
    js = HU.keys_js()
    esc = js[js.index("ev.key === 'Escape'"):]
    assert ".delfin-confirm .delfin-confirm-cancel" in esc[:400]
    assert "delfin-confirm-yes" not in js
    c = RB.Confirm(timer=False)
    assert "delfin-confirm-cancel" in c.no._dom_classes
    assert "delfin-confirm-yes" in c.yes._dom_classes


# -- in the header ----------------------------------------------------------------

def test_the_confirmation_sits_under_more_and_closing_more_withdraws_it(repo):
    import ipywidgets as W
    from delfin.dashboard.branch_panel import BranchPanel
    c = RB.Confirm(timer=False)
    back = W.Button(description="Roll back")
    panel = BranchPanel(repo, gh=None, synchronous=True, pull=W.Button(),
                        more_extra=[back], more_below=[c.box])
    panel.on_more_closed(c.cancel)
    assert c.box in panel.more_below.children and panel.more_below in panel.box.children
    assert panel.more_below.layout.display == "none"
    panel.more_btn.value = True
    assert panel.more_below.layout.display != "none"
    c.ask("x", "Roll back", lambda: None)
    panel.more_btn.value = False
    assert panel.more_below.layout.display == "none" and not c.is_open


def test_the_dashboard_asks_and_resets_only_from_the_confirmation():
    import inspect
    from delfin import dashboard as d
    src = inspect.getsource(d.create_dashboard)
    handler = src[src.index("def handle_rollback_delfin"):src.index("def _do_rollback_delfin")]
    assert "rollback_confirm.ask(" in handler
    assert "'reset'" not in handler and "_run_git_write" not in handler
    assert "more_below=[rollback_confirm.box]" in src
    assert "on_more_closed(rollback_confirm.cancel)" in src
    assert "HEAD~1'" not in src                  # the target is the planned sha
