"""The header bar: few controls, the rest one click away, and nothing that stays.

Input: :class:`branch_panel.BranchPanel` on a checkout made here (main, a
merged branch, an unmerged one, a pinned one, a merged branch held by a
worktree, a branch only origin has, a tag), with preview start/close
replaced by stand-ins, a Pull and a roll-back button handed in as the
dashboard does, and a temporary HOME.

Output, by clicking the widgets:

* top row: the branch switcher, Test, Pull, the more toggle; no filter box,
  no Clean up, no roll-back -- those are in the more row, hidden until
  asked, and the open/closed state is kept in ``~/.delfin``;
* the switcher opens a menu with a search field and the views Mine, All,
  Hidden, Merged, Remote; typing narrows, Enter picks the first row; a
  typed tag becomes a candidate row;
* every row but main has a hide action that acts at once, leaves git
  as it was, and can be undone from the message; a delete action exists
  only on a local branch origin/main already contains and that is not
  main, current or held, and deletes after one confirmation with the undo
  command in the details;
* Test shows a chip with an Open link and a close action; a close that
  fails still ends in a final message;
* messages: the latest one only, dismissible, success ones expire, earlier
  ones on demand;
* Clean up: a search, one select-all toggle acting on the visible rows
  and never ticking pinned or unmerged ones, and an "n of m selected"
  counter.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import time
from pathlib import Path

import pytest

pytest.importorskip("ipywidgets")

from delfin.dashboard import branch_list as BL  # noqa: E402
from delfin.dashboard import branch_preview as BP  # noqa: E402
from delfin.dashboard import header_ui as HU  # noqa: E402
from delfin.dashboard.branch_panel import BranchPanel  # noqa: E402

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}


def _git(cwd, *args, check=True) -> str:
    out = subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", *args],
        cwd=str(cwd), env=_ENV, check=check, capture_output=True, text=True)
    return out.stdout.strip()


def _commit(cwd, name):
    (cwd / name).write_text(name + "\n")
    _git(cwd, "add", name)
    _git(cwd, "commit", "-q", "-m", name)


@pytest.fixture
def repo(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv(BL.GH_ENV, "off")
    (tmp_path / "home").mkdir()
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    co = tmp_path / "checkout"
    _git(tmp_path, "clone", "-q", str(origin), str(co))
    _git(co, "checkout", "-q", "-b", "main")
    _commit(co, "a.txt")
    _git(co, "push", "-q", "-u", "origin", "main")
    _git(co, "branch", "merged-one")
    _git(co, "branch", "merged-held")
    _git(co, "branch", "mine/pinned")
    _git(co, "worktree", "add", "-q", str(tmp_path / "wt-held"), "merged-held")
    _git(co, "checkout", "-q", "-b", "feat")
    _commit(co, "f.txt")
    _git(co, "checkout", "-q", "-b", "far/away")
    _commit(co, "r.txt")
    _git(co, "push", "-q", "origin", "far/away")
    _git(co, "checkout", "-q", "main")
    _git(co, "branch", "-q", "-D", "far/away")
    _git(co, "tag", "v1", "main")
    _git(co, "fetch", "-q", "origin")
    BL.pin(co, "mine/pinned")
    return co


class _Fake:
    def __init__(self):
        self.started, self.closed = [], []
        self.fail_close = None

    def start(self, repo, ref, **_kw):
        self.started.append(ref)
        from delfin.agent import proc_identity
        p = BP.Preview(
            slug=BP.slug(ref), ref=ref, sha=_git(repo, "rev-parse", ref),
            branch=ref, repo=str(repo), worktree=str(BP.worktree_path(repo, ref)),
            port=45678, token="tok-abc", url_path="/voila/render/x.ipynb",
            pid=os.getpid(), proc_start=proc_identity.process_start(os.getpid()),
            ticks=None, host=socket.gethostname(), started_at=time.time(),
            delfin_file="", log="")
        BP._save(p)
        return True, f"Preview of {ref}", p

    def close(self, s, **_kw):
        self.closed.append(s)
        if self.fail_close:
            raise self.fail_close
        BP._record_path(s).unlink()
        return BP.CloseResult(True, "closed", "Closed.")


def _panel(repo, **kw):
    import ipywidgets as W
    fake = _Fake()
    pull = W.Button(description="Pull")
    back = W.Button(description="Roll back")
    panel = BranchPanel(repo, gh=None, start=fake.start, close=fake.close,
                        synchronous=True, pull=pull, more_extra=[back], **kw)
    panel.refresh()
    return panel, fake, pull, back


def _shown(panel):
    return [r for r in panel.shown_refs()]


def _git_state(co):
    return (_git(co, "for-each-ref"), _git(co, "status", "--porcelain"),
            _git(co, "branch", "--show-current"), _git(co, "stash", "list"))


# -- structure ----------------------------------------------------------------

def test_the_top_row_has_the_switcher_test_pull_and_more_only(repo):
    import ipywidgets as W
    panel, _f, pull, back = _panel(repo)
    top = list(panel.row.children)
    assert top[:2] == [panel.branch_wrap, panel.test_btn]
    assert pull in top and panel.more_btn in top
    assert back not in top and panel.cleanup_btn not in top
    assert not any(isinstance(w, (W.Dropdown, W.Text)) for w in top)
    assert panel.branch_btn.description.startswith("⎇ main")
    assert panel.test_btn.description == "▶ Test"
    assert list(panel.more_row.children) == [panel.cleanup_btn, back]
    assert panel.more_row.layout.display == "none"
    assert panel.search in panel.pop.children[0].children
    assert panel.pop.layout.display == "none"


def test_more_opens_and_its_state_is_remembered(repo):
    panel, *_ = _panel(repo)
    panel.more_btn.value = True
    assert panel.more_row.layout.display != "none"
    assert HU.load_prefs().get("more_open") is True
    again, *_ = _panel(repo)
    assert again.more_btn.value is True
    assert again.more_row.layout.display != "none"


# -- the switcher -------------------------------------------------------------

def test_search_inside_the_menu_narrows_and_enter_picks(repo):
    panel, *_ = _panel(repo)
    panel.open_menu()
    assert panel.pop.layout.display != "none"
    assert _shown(panel) == ["main", "mine/pinned", "feat"]
    panel.search.value = "fe"
    assert _shown(panel) == ["main", "feat"]
    panel.search_submit()
    assert panel.value == "main" or panel.value == "feat"
    panel.open_menu()
    panel.search.value = "feat"
    panel.view.value = "All"
    assert _shown(panel) == ["feat"]
    panel.search_submit()
    assert panel.value == "feat"
    assert panel.pop.layout.display == "none"
    assert panel.branch_btn.description.startswith("⎇ feat")


def test_the_views_show_what_the_default_list_leaves_out(repo):
    panel, *_ = _panel(repo)
    panel.open_menu()
    panel.view.value = "Remote"
    assert _shown(panel) == ["origin/far/away"]
    panel.view.value = "Merged"
    assert set(_shown(panel)) == {"merged-one", "merged-held", "mine/pinned"}
    panel.view.value = "All"
    panel.search.value = "v1"
    assert _shown(panel) == ["v1"]                  # a typed tag
    assert panel.row_widgets("v1").pin is not None


def test_hide_acts_at_once_leaves_git_alone_and_can_be_undone(repo):
    before = _git_state(repo)
    panel, *_ = _panel(repo)
    panel.open_menu()
    assert panel.row_widgets("main").hide is None
    panel.row_widgets("feat").hide.click()
    assert "feat" not in _shown(panel)
    assert _git_state(repo) == before
    msg = panel.toasts.latest()
    assert "feat" in msg["text"] and [a for a, _ in msg["actions"]] == ["Undo"]
    panel.toasts.run_action(0)
    assert "feat" in _shown(panel)
    panel.row_widgets("mine/pinned").hide.click()
    panel.view.value = "Hidden"
    assert _shown(panel) == ["mine/pinned"]
    panel.row_widgets("mine/pinned").show.click()
    panel.view.value = "Mine"
    assert "mine/pinned" in _shown(panel)
    assert _git_state(repo) == before


def test_delete_only_where_allowed_and_after_one_confirmation(repo):
    panel, *_ = _panel(repo)
    panel.open_menu()
    for ref in ("main", "feat"):
        assert panel.row_widgets(ref).delete is None, ref
    panel.view.value = "All"
    assert panel.row_widgets("feat").delete is None           # unmerged
    assert panel.row_widgets("merged-held").delete is None    # held
    assert panel.row_widgets("origin/far/away").delete is None
    sha = _git(repo, "rev-parse", "merged-one")
    panel.row_widgets("merged-one").delete.click()
    assert _git(repo, "rev-parse", "merged-one") == sha       # asked first
    yes = panel.confirm_box.children[1]
    yes.click()
    assert _git(repo, "rev-parse", "--verify", "--quiet", "merged-one",
                check=False) == ""
    msg = panel.toasts.latest()
    assert "merged-one" in msg["text"] and sha in msg["details"]


def test_a_remote_only_branch_is_added_to_the_list(repo):
    panel, *_ = _panel(repo)
    panel.open_menu()
    panel.view.value = "Remote"
    panel.row_widgets("origin/far/away").pin.click()
    panel.view.value = "Mine"
    assert "origin/far/away" in _shown(panel)


def test_new_branch_is_the_last_row_and_does_not_switch(repo):
    before = _git(repo, "branch", "--show-current")
    panel, *_ = _panel(repo)
    panel.open_menu()
    assert panel.rows_box.children[-1] is panel.new_row
    panel.new_text.value = "max/idea"
    panel.new_btn.click()
    assert panel.value == "max/idea"
    assert _git(repo, "branch", "--show-current") == before
    assert "max/idea" in BL.load_prefs(repo)["pinned"]


# -- test and the chip --------------------------------------------------------

def test_test_shows_a_chip_with_open_and_close(repo):
    panel, fake, *_ = _panel(repo)
    panel.select("feat")
    panel.test_btn.click()
    assert fake.started == ["feat"]
    chip = panel.chip_widgets("feat")
    assert "Open" in chip.link.value and 'target="_blank"' in chip.link.value
    assert "45678" not in chip.text.value and "tok-abc" not in chip.text.value
    chip.close.click()
    assert fake.closed == ["feat"]
    assert panel.chip_widgets("feat") is None
    assert panel.toasts.latest()["text"].startswith("Closed")


def test_a_failing_close_still_ends_in_a_final_message(repo):
    panel, fake, *_ = _panel(repo)
    panel.select("feat")
    panel.test_btn.click()
    fake.fail_close = OSError("stale file handle")
    panel.chip_widgets("feat").close.click()
    msg = panel.toasts.latest()
    assert msg["level"] == "error" and "Couldn't close" in msg["text"]
    assert "stale file handle" in msg["text"] + msg["details"]
    chip = panel.chip_widgets("feat")
    assert chip is None or "Closing" not in chip.text.value


# -- messages -----------------------------------------------------------------

def test_messages_show_the_latest_dismiss_and_expire():
    t = HU.Toasts(timer=False)
    t.show("first", level="ok")
    t.show("second", level="error", details="why")
    assert t.latest()["text"] == "second"
    assert len(t.history) == 2
    t.dismiss()
    assert t.latest() is None and t.box.children == ()
    t.show("done", level="ok")
    t.expire()                          # what the timer does after a few s
    assert t.latest() is None
    t.show("broken", level="error")
    t.expire()
    assert t.latest()["text"] == "broken"   # errors wait for the person


# -- clean up -----------------------------------------------------------------

def test_clean_up_search_select_all_visible_and_counter(repo):
    _git(repo, "branch", "merged-two")
    panel, *_ = _panel(repo)
    panel.more_btn.value = True
    panel.cleanup_btn.value = True
    names = [c.name for c in panel._cands]
    assert {"merged-one", "merged-two", "mine/pinned"} <= set(names)
    panel.include_unmerged.value = True
    assert "feat" in [c.name for c in panel._cands]
    ticked = {c.name for c, cb in zip(panel._cands, panel._checks) if cb.value}
    assert "feat" not in ticked and "mine/pinned" not in ticked
    panel.select_toggle.click()                    # select none
    assert not any(cb.value for cb in panel._checks)
    assert panel.cleanup_count.value.startswith("0 of ")
    panel.cleanup_search.value = "merged-t"
    vis = [c.name for c in panel.visible_cands()]
    assert vis == ["merged-two"]
    panel.select_toggle.click()                    # select all (visible)
    ticked = {c.name for c, cb in zip(panel._cands, panel._checks) if cb.value}
    assert ticked == {"merged-two"}
    assert panel.cleanup_count.value.startswith("1 of ")
    panel.cleanup_search.value = ""
    panel.select_toggle.click()                    # all visible eligible
    ticked = {c.name for c, cb in zip(panel._cands, panel._checks) if cb.value}
    assert "feat" not in ticked and "mine/pinned" not in ticked
    assert {"merged-one", "merged-two"} <= ticked


# -- the dashboard's assembly -------------------------------------------------

def _create_dashboard_source() -> str:
    import inspect
    from delfin import dashboard as d
    return inspect.getsource(d.create_dashboard)


def test_the_dashboard_puts_pull_in_the_bar_and_roll_back_behind_more():
    src = _create_dashboard_source()
    assert "pull=pull_delfin_btn" in src
    assert "more_extra=[rollback_delfin_btn]" in src
    assert "_preview_bar(header_toasts)" in src     # a preview: banner only
    assert "header_toasts.box" in src               # one message area
    assert "pull_delfin_output" not in src          # nothing piles up under the bar
    assert "_hu.keys_js()" in src


def test_the_status_hint_reads_up_to_date_behind_and_ahead(repo, tmp_path):
    from delfin.dashboard import _git_status_hint
    h = _git_status_hint(repo)
    assert h["text"] == "main ✓ up to date" and h["state"] == "ok"
    other = tmp_path / "other"
    _git(tmp_path, "clone", "-q", str(tmp_path / "origin.git"), str(other))
    _commit(other, "n1.txt")
    _commit(other, "n2.txt")
    _git(other, "push", "-q", "origin", "HEAD:main")
    _git(repo, "fetch", "-q", "origin")
    h = _git_status_hint(repo)
    assert h["text"] == "main ↓ 2 new commits" and h["state"] == "behind"
    assert h["icon"] == "↓2" and "fast-forward" in h["title"]
    _git(repo, "merge", "-q", "--ff-only", "origin/main")
    _commit(repo, "local.txt")
    h = _git_status_hint(repo)
    assert h["text"] == "main ↑ 1 local commit" and h["state"] == "warn"


def test_keep_session_is_drawn_as_a_switch_with_its_label():
    from delfin.dashboard import session as S
    strip = S.build_status_strip()
    toggle = next(c for c in strip.children if hasattr(c, "value")
                  and isinstance(getattr(c, "value"), bool))
    assert "delfin-switch" in toggle._dom_classes
    assert toggle.description == "Keep session"
    assert "browser" in toggle.tooltip
    assert ".jupyter-button.delfin-switch.mod-active::after" in HU.HEADER_CSS


def test_the_header_style_never_sets_visibility_with_important():
    """Widgets show and hide through layout.display (an inline style); a
    rule with !important would win over it and pin a panel open."""
    import re
    css = HU.HEADER_CSS
    assert not re.search(r"display\s*:[^;{}]*!important", css)
    assert "var(--jp-layout-color1" in css and "var(--jp-ui-font-color1" in css
