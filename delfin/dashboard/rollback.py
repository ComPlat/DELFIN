"""Roll back 1 commit: the plan, the confirmation, the reset.

Input: the shared checkout (``repo``); for :func:`execute`, the git
runners of the dashboard (read: no optional locks; write: index-lock
recovery), so the reset goes through the same paths as Pull. Output:

* :func:`plan` -- ``(Plan, "")`` or ``(None, reason)``. A plan names the
  branch, the commit that leaves it (short sha, subject), the commit it
  moves to (the first parent), how many commits leave the branch with it
  (more than one when the commit is a merge), the remote branches that
  still contain it, and whether this checkout writes a reflog. Refused:
  no repository, uncommitted or untracked changes, no earlier commit.
* :func:`describe` -- the confirmation text (HTML): what moves where,
  that the action is ``git reset --hard <target>`` on the shared checkout
  (every session on it runs the older code), and how it is undone:
  ``git reset --hard <old sha>``; ``git reflog`` keeps the old commit
  (an unreachable entry expires after ``gc.reflogExpireUnreachable``,
  30 days by default); ↓ Pull brings it back when ``origin/main``
  contains it and the branch is main.
* :func:`execute` -- the reset to the planned target, refused when HEAD
  is no longer the planned commit (a Pull, a commit or another rollback
  happened in between) or the tree has changes by then.
* :class:`Confirm` -- the question in the More row: the text, a danger
  button that names the action, and Cancel. Only a deliberate click
  answers: the confirm button is disabled for :data:`HOLD_S` after the
  question opens, and a click that arrives earlier (a double click, a key
  repeat on a focused button) is ignored by the handler as well. Esc
  presses Cancel (``header_ui.keys_js``); no key presses the confirm
  button, and a new question replaces an open one without answering it.

Rejected: ``git revert`` (adds a commit on main that nobody wrote on a
branch, and Pull's fast-forward against origin/main fails afterwards);
asking in a message (the toast area expires success messages and is
shared with Pull's output); a second click on the same button as the
confirmation (agent_sessions' Stop all pattern): a double click would
roll back without the text having been shown.

Cost: one more click and :data:`HOLD_S` seconds per rollback; seven git
reads when the question opens (70 ms on the DELFIN checkout, 2026-10-15).
"""

from __future__ import annotations

import html
import subprocess
import threading
import time
from dataclasses import dataclass, field
from typing import Callable, Optional

#: The confirm button answers this long after the question opened, the
#: hold of the terminal window's approvals (tui footer/approval.tsx
#: HOLD_MS 500): a click or key meant for what was there before never
#: lands on it.
HOLD_S = 0.5


@dataclass(frozen=True)
class Plan:
    branch: str              # "" when detached
    head: str                # the commit that leaves (full sha)
    head_short: str
    head_subject: str
    target: str              # the first parent HEAD moves to (full sha)
    target_short: str
    target_subject: str
    leaving: int = 1         # commits no longer on the branch afterwards
    on_remote: tuple = field(default=())   # remote branches containing head
    reflog: bool = True      # core.logAllRefUpdates is not false


def _read(repo, *args) -> subprocess.CompletedProcess:
    from delfin.agent import safe_git
    return safe_git.run(args, repo, timeout=60, merge_stderr=True)


def _write(repo, *args) -> subprocess.CompletedProcess:
    from delfin.agent import safe_git
    return safe_git.run(args, repo, timeout=120, merge_stderr=True)


def plan(repo, run: Optional[Callable] = None) -> tuple[Optional[Plan], str]:
    """What Roll back 1 commit would do in *repo* now (see the module doc)."""
    run = run or _read
    if not repo:
        return None, "No DELFIN checkout found for this dashboard."
    st = run(repo, "status", "--porcelain")
    if st.returncode != 0:
        return None, (st.stdout or "git status failed").strip()
    if st.stdout.strip():
        return None, ("Not rolled back: the shared checkout has uncommitted changes. "
                      "Commit, stash or discard them first.")
    fmt = "%H%x00%h%x00%s"
    head = run(repo, "log", "-1", f"--format={fmt}", "HEAD")
    target = run(repo, "log", "-1", f"--format={fmt}", "HEAD~1")
    if head.returncode != 0 or target.returncode != 0 or not target.stdout.strip():
        return None, "Not rolled back: there is no earlier commit."
    h, hs, hsub = (head.stdout.strip().split("\x00") + ["", "", ""])[:3]
    t, ts, tsub = (target.stdout.strip().split("\x00") + ["", "", ""])[:3]
    branch = run(repo, "branch", "--show-current").stdout.strip()
    try:
        leaving = int(run(repo, "rev-list", "--count", f"{t}..{h}").stdout.strip() or 1)
    except ValueError:
        leaving = 1
    remote = run(repo, "branch", "-r", "--contains", h, "--format=%(refname:short)")
    on_remote = tuple(sorted(x.strip() for x in remote.stdout.splitlines()
                             if x.strip() and not x.strip().endswith("/HEAD")
                             )) if remote.returncode == 0 else ()
    cfg = run(repo, "config", "--bool", "core.logAllRefUpdates").stdout.strip()
    return Plan(branch=branch, head=h, head_short=hs, head_subject=hsub, target=t,
                target_short=ts, target_subject=tsub, leaving=leaving,
                on_remote=on_remote, reflog=cfg != "false"), ""


def undo_command(p: Plan) -> str:
    return f"git reset --hard {p.head_short}"


def describe(p: Plan) -> str:
    """The confirmation text (HTML) for *p*."""
    e = html.escape
    where = f"<b>{e(p.branch)}</b>" if p.branch else "the detached checkout"
    lines = [
        f"Roll back {where} by one commit?",
        f"Leaves: <code>{e(p.head_short)}</code> {e(p.head_subject)}",
        f"Back to: <code>{e(p.target_short)}</code> {e(p.target_subject)}",
    ]
    if p.leaving > 1:
        lines.append(f"It is a merge: {p.leaving} commits leave {e(p.branch or 'HEAD')} "
                     "with it.")
    lines.append(f"Runs <code>git reset --hard {e(p.target_short)}</code> in the shared "
                 "checkout; every session on it runs the older code. Nothing "
                 "uncommitted is lost: the checkout is clean.")
    undo = f"Undo: <code>{e(undo_command(p))}</code>"
    if p.reflog:
        undo += ("; <code>git reflog</code> lists it until gc expires the entry "
                 "(30 days by default)")
    else:
        undo += "; this checkout keeps no reflog, so note the sha"
    if p.branch == "main" and "origin/main" in p.on_remote:
        undo += ". origin/main has it, so ↓ Pull also brings it back"
    lines.append(undo + ".")
    return "<br>".join(lines)


def execute(repo, p: Plan, *, read: Optional[Callable] = None,
            write: Optional[Callable] = None) -> tuple[bool, str, list]:
    """Reset *repo* to ``p.target`` when it is still at ``p.head`` and clean.

    ``read(repo, *args)`` / ``write(repo, *args)`` return a
    CompletedProcess with the output in ``stdout``. Returns ``(ok, text,
    details)``."""
    read = read or _read
    write = write or _write
    now = read(repo, "rev-parse", "HEAD").stdout.strip()
    if now != p.head:
        return False, "Not rolled back: the checkout moved since you were asked.", [
            f"asked about {p.head_short}, now at {now[:9] or '?'}; open Roll back again"]
    st = read(repo, "status", "--porcelain")
    if st.returncode != 0 or st.stdout.strip():
        return False, "Not rolled back: the shared checkout has uncommitted changes.", [
            (st.stdout or "").strip()]
    result = write(repo, "reset", "--hard", p.target)
    details = [f"{p.head_short} -> {p.target_short}", (result.stdout or "").strip()
               or "(no output)"]
    if result.returncode != 0:
        return False, "Couldn't roll back.", details
    details.append(f"Undo: {undo_command(p)}")
    return True, (f"Rolled back {p.branch or 'HEAD'} to {p.target_short}. "
                  f"Undo: {undo_command(p)}"), details


class Confirm:
    """The confirmation box (see the module docstring).

    ``now`` and ``timer`` are for tests: with ``timer=False`` no thread
    enables the button (:meth:`release` does); the handler's own check
    of the hold is the same either way."""

    def __init__(self, *, hold_s: float = HOLD_S, now: Callable = time.monotonic,
                 timer: bool = True):
        import ipywidgets as W
        from delfin.dashboard import header_ui as HU
        self.W, self.HU = W, HU
        self.hold_s = float(hold_s)
        self.now = now
        self.timer = timer
        self.shown_at = 0.0
        self._seq = 0
        self._on_yes: Optional[Callable] = None
        self.text = W.HTML("")
        self.yes = HU.icon_button("Roll back", "", "delfin-btn", "is-danger", "is-filled",
                                  "delfin-confirm-yes")
        self.no = HU.icon_button("Cancel", "Esc", "delfin-btn", "delfin-confirm-cancel")
        self.box = HU.classed(W.HBox([self.text, self.yes, self.no]), "delfin-confirm",
                              "delfin-confirm-wide")
        self.box.layout.display = "none"
        self.yes.on_click(lambda _b: self._yes())
        self.no.on_click(lambda _b: self.cancel())

    @property
    def is_open(self) -> bool:
        return self.box.layout.display != "none"

    def ask(self, text_html: str, label: str, on_yes: Callable) -> None:
        """Show the question; *on_yes* runs on a confirm after the hold."""
        self._seq += 1
        seq = self._seq
        self._on_yes = on_yes
        self.text.value = text_html
        self.yes.description = label
        self.yes.disabled = True
        self.shown_at = self.now()
        self.box.layout.display = ""
        if self.timer and self.hold_s > 0:
            t = threading.Timer(self.hold_s, self._release_if, args=(seq,))
            t.daemon = True
            t.start()
        elif self.hold_s <= 0:
            self.yes.disabled = False

    def _release_if(self, seq: int) -> None:
        if seq == self._seq and self.is_open:
            self.yes.disabled = False

    def release(self) -> None:
        """Enable the confirm button now (the hold's timer, by hand)."""
        self._release_if(self._seq)

    def cancel(self) -> None:
        self._on_yes = None
        self.box.layout.display = "none"
        self.text.value = ""
        self.yes.disabled = True

    def _yes(self) -> None:
        if not self.is_open or self._on_yes is None:
            return
        if self.now() - self.shown_at < self.hold_s:
            return              # within the hold: not an answer
        fn = self._on_yes
        self.cancel()
        fn()
