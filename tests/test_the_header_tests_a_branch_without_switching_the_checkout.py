"""The dashboard header tests a branch in a preview and never switches the checkout.

Input: the header's branch controls (``branch_panel.BranchPanel``) on a
checkout made here, with the preview start/stop replaced by stand-ins that
record what they were asked. Output, by clicking the widgets: Test branch
asks for a preview of the selected ref and the running previews are listed
as links that open a new tab; Close preview stops it; Manage list adds,
creates from main, hides, shows again and deletes (a second confirmation
for unmerged commits); Clean up lists, asks, deletes and shows the undo
commands. The checkout's branch, head and stash stay as they were. In a
preview the banner names the branch, the commit and where the code comes
from, and warns when that is not the preview worktree.

The dashboard source is read for the removed SWITCH path: no stash, no
switch of the shared checkout, PULL is a fast-forward of main.
"""

from __future__ import annotations

import os
import re
import socket
import subprocess
import time
from pathlib import Path

import pytest

pytest.importorskip("ipywidgets")

import delfin  # noqa: E402
from delfin.dashboard import branch_cleanup as BC  # noqa: E402
from delfin.dashboard import branch_list as BL  # noqa: E402
from delfin.dashboard import branch_preview as BP  # noqa: E402
from delfin.dashboard.branch_panel import BranchPanel  # noqa: E402

_ENV = {**os.environ, "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull}
_SRC = Path(delfin.__file__).resolve().parent / "dashboard" / "__init__.py"


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
    (tmp_path / "home").mkdir()
    origin = tmp_path / "origin.git"
    _git(tmp_path, "init", "-q", "--bare", "-b", "main", str(origin))
    co = tmp_path / "checkout"
    _git(tmp_path, "clone", "-q", str(origin), str(co))
    _git(co, "checkout", "-q", "-b", "main")
    _commit(co, "a.txt")
    _git(co, "push", "-q", "-u", "origin", "main")
    _git(co, "branch", "merged-one")
    _git(co, "checkout", "-q", "-b", "feat")
    _commit(co, "f.txt")
    _git(co, "checkout", "-q", "main")
    (co / "wip.txt").write_text("uncommitted\n")
    return co


def _state(co):
    return (_git(co, "branch", "--show-current"), _git(co, "rev-parse", "HEAD"),
            _git(co, "stash", "list"), (co / "wip.txt").read_text())


class _Fake:
    def __init__(self, repo):
        self.repo = repo
        self.started, self.stopped = [], []

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

    def stop(self, s, **_kw):
        self.stopped.append(s)
        BP._record_path(s).unlink()
        return True, f"Stopped {s}", False


def _panel(repo):
    fake = _Fake(repo)
    panel = BranchPanel(repo, gh=None, start=fake.start, stop=fake.stop,
                        synchronous=True)
    panel.refresh()
    return panel, fake


def _values(panel):
    return [v for _l, v in panel.dropdown.options]


def test_the_list_hides_merged_and_filters(repo):
    panel, _ = _panel(repo)
    assert _values(panel) == ["main", "feat"]
    panel.filter.value = "fe"
    assert _values(panel) == ["main", "feat"]
    panel.filter.value = "zzz"
    assert _values(panel) == ["main"]


def test_test_branch_opens_a_preview_link_and_close_stops_it(repo):
    before = _state(repo)
    panel, fake = _panel(repo)
    panel.dropdown.value = "feat"
    panel.test_btn.click()
    assert fake.started == ["feat"]
    assert len(panel.previews_box.children) == 1
    link = panel.previews_box.children[0].children[0].value
    assert 'target="_blank"' in link and "45678" in link and "tok-abc" in link
    assert "location.hostname" in link             # the page's own host
    close = panel.previews_box.children[0].children[1]
    close.click()
    assert fake.stopped == ["feat"]
    assert panel.previews_box.children == ()
    assert _state(repo) == before


def test_manage_add_new_hide_show_again(repo):
    before = _state(repo)
    panel, _ = _panel(repo)
    panel.add_text.value = "merged-one"
    panel.add_btn.click()
    assert "merged-one" in _values(panel)
    panel.new_text.value = "max/idea"
    panel.new_btn.click()
    assert panel.dropdown.value == "max/idea"
    assert "Test branch" in panel.status.value
    assert _state(repo) == before                  # no switch, files kept
    panel.dropdown.value = "feat"
    panel.hide_btn.click()
    assert "feat" not in _values(panel)
    assert [v for _l, v in panel.hidden_dd.options] == ["feat"]
    panel.hidden_dd.value = "feat"
    panel.unhide_btn.click()
    assert "feat" in _values(panel)


def test_delete_from_the_list_asks_twice_for_unmerged_commits(repo):
    panel, _ = _panel(repo)
    panel.dropdown.value = "feat"
    panel.delete_btn.click()
    yes = panel.confirm_box.children[1]
    yes.click()                                    # first confirmation
    assert _git(repo, "rev-parse", "--verify", "feat")
    assert "not in origin/main" in panel.confirm_box.children[0].value
    panel.confirm_box.children[1].click()          # second confirmation
    assert _git(repo, "rev-parse", "--verify", "--quiet", "feat",
                check=False) == ""
    assert "Undo:" in panel.status.value
    # main cannot be deleted from here
    panel.dropdown.value = "main"
    panel.delete_btn.click()
    assert panel.confirm_box.children == ()
    assert "cannot be deleted" in panel.status.value


def test_clean_up_lists_asks_deletes_and_shows_the_undo(repo):
    panel, _ = _panel(repo)
    panel.cleanup_btn.value = True
    labels = [cb.description for cb in panel._checks]
    assert any("merged-one" in x for x in labels)
    assert not any(" feat " in x for x in labels)
    panel.cleanup_delete_btn.click()
    assert "Delete 1 local branch" in panel.cleanup_confirm.children[0].value
    panel.cleanup_confirm.children[1].click()
    assert "Deleted 1" in panel.cleanup_result.value
    assert "git -C" in panel.cleanup_result.value
    assert str(BC._log_path()) in panel.cleanup_result.value
    assert _git(repo, "rev-parse", "--verify", "--quiet", "merged-one",
                check=False) == ""


def test_the_preview_banner_names_branch_commit_and_code(tmp_path):
    root = Path(delfin.__file__).resolve().parent.parent
    env = {BP.PREVIEW_ENV_BRANCH: "feat", BP.PREVIEW_ENV_WORKTREE: str(root),
           BP.PREVIEW_ENV_SHA: "0123456789"}
    html = BP.banner_html(env)
    assert "PREVIEW: feat @" in html and f"code from {root}" in html
    assert "WARNING" not in html
    env[BP.PREVIEW_ENV_WORKTREE] = str(tmp_path)
    assert "WARNING" in BP.banner_html(env)
    assert BP.in_preview(env) and not BP.in_preview({})


def test_the_dashboard_no_longer_switches_or_stashes_the_shared_checkout():
    src = _SRC.read_text()
    assert "handle_switch_branch" not in src
    assert "description='SWITCH'" not in src
    assert not re.search(r"_run_git_\w+\([^)]*'stash'", src)
    assert not re.search(r"_run_git_write\([^)]*'switch'", src)
    assert "'pull', '--ff-only', 'origin', 'main'" in src
    assert "BranchPanel(" in src and "banner_html()" in src


def test_compute_branch_options_follows_the_new_rule(repo):
    from delfin.dashboard import _compute_branch_options
    values = [v for _l, v in _compute_branch_options(repo)]
    assert values == ["main", "feat"]
    assert BL.branch_entries(None) [0].ref == "main"


def test_clean_up_does_not_preselect_a_branch_the_user_pinned(repo):
    """A branch just made from main has no commits of its own, so git calls
    it merged; it is the user's own entry, and is offered unselected."""
    panel, _ = _panel(repo)
    panel.new_text.value = "max/fresh"
    panel.new_btn.click()
    panel.cleanup_btn.value = True
    boxes = {cb.description.split(" ")[1]: cb for cb in panel._checks
             if cb.description.startswith("branch ")}
    assert boxes["max/fresh"].value is False
    assert "pinned" in boxes["max/fresh"].description
    assert boxes["merged-one"].value is True
    panel.cleanup_delete_btn.click()
    panel.cleanup_confirm.children[1].click()
    assert "Deleted 1." in panel.cleanup_result.value
    assert _git(repo, "rev-parse", "--verify", "max/fresh")
