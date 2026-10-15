"""The header's branch controls: switcher, Test, preview chips, more, clean up.

Input: the shared checkout (``repo``); for tests, the functions that start
and close previews and the ``gh`` to ask; from the dashboard, its Pull and
roll-back buttons and the status hint. Output: :class:`BranchPanel` with

* ``row`` -- the controls in the header line: the branch switcher (a pill
  that opens a menu: search, the views Mine / All / Hidden / Merged /
  Remote, a reload, one row per branch with its actions, and "New branch
  from main" last), ``▶ Test``, the running previews as chips, Pull, the
  more toggle, the status hint;
* ``box`` -- the more row (Clean up, roll back) and the Clean up panel,
  under the bar;
* ``toasts`` -- the message area (:class:`header_ui.Toasts`).

Rules kept from the list this replaces: nothing here switches the shared
checkout (Test opens a preview in a worktree, ``branch_preview.start``);
hiding changes only the per-user list in ``~/.delfin`` and is undone from
its message; a branch is deleted only through ``branch_cleanup`` with its
rules, re-checked at deletion. The menu offers delete only on a local
branch ``origin/main`` already contains that is not main, current or held
(one confirmation); a branch with commits of its own is deleted from Clean
up with "Also branches with their own commits" and the second
confirmation.

Long work (gh, fetch, starting or closing a preview, the clean-up listing)
runs in a thread and ends in a message, so the kernel keeps answering.
"""

from __future__ import annotations

import html
import os
import subprocess
import sys
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

from delfin.dashboard import branch_cleanup as BC
from delfin.dashboard import branch_list as BL
from delfin.dashboard import branch_preview as BP
from delfin.dashboard import header_ui as HU

SWITCH_NOTE = ("Test opens the selected branch in its own preview; the shared "
               "checkout stays on main.")

VIEWS = [("Mine", "mine"), ("All", "all"), ("Hidden", "hidden"),
         ("Merged", "merged"), ("Remote", "remote")]
_VIEW_TIPS = ["main, your pinned branches, open pull requests, branches with new commits",
              "every branch, here and on origin", "branches you hid",
              "branches already in main", "branches only origin has"]
MAX_ROWS = 40


def _esc(text) -> str:
    return html.escape(str(text or ""))


def preview_href(p) -> str:
    """Where the chip's Open link points: through the dashboard's own server
    when it proxies previews, else the preview's port on the page's host."""
    rel = BP.proxy_href(p)
    if rel:
        return rel
    return f"http://localhost:{int(p.port)}{p.url_path}?token={p.token}"


def preview_link_html(p, label: str = "Open ↗") -> str:
    """An anchor that opens the preview in a new tab. A direct-port link is
    rewritten on click to the host the page was loaded from."""
    href = preview_href(p)
    js = ""
    title = "Open the preview in a new tab"
    if href.startswith("http://localhost:"):
        path = href.split(f":{int(p.port)}", 1)[1]
        js = (f' onclick="this.href=location.protocol+\'//\'+location.hostname+'
              f'\':{int(p.port)}\'+{_esc(repr(path))};"')
        # Only without the proxy: then the preview has a port of its own.
        title += f" (through an SSH tunnel, forward port {int(p.port)} too)"
    return (f'<a class="delfin-preview-link" href="{_esc(href)}" target="_blank" '
            f'rel="noopener noreferrer"{js} title="{_esc(title)}">'
            f'{_esc(label)}</a>')


@dataclass
class RowW:
    box: object
    main: object
    test: object
    hide: object = None
    show: object = None
    pin: object = None
    delete: object = None


@dataclass
class ChipW:
    box: object
    text: object
    link: object = None
    close: object = None


def _on_enter(text_widget, fn) -> None:
    """Enter in a text box. ``Text.on_submit`` is deprecated in ipywidgets 8
    but is the only event that tells Enter from a value change."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        text_widget.on_submit(lambda _w: fn())


def _row_icon(r: BL.BranchRow) -> str:
    if r.kind == "main":
        return "⎇"
    if r.pr is not None:
        return "⇄"
    if r.pinned:
        return "📌"
    return {"remote": "☁", "tag": "🏷", "commit": "◆", "lookup": "🔍"}.get(r.kind, "•")


def _ci_html(ci: str) -> str:
    return {"passing": '<span class="ok">✓</span>', "failing": '<span class="bad">✗</span>',
            "pending": "…"}.get(ci, "")


class BranchPanel:
    def __init__(self, repo, *, gh: Optional[str] = "auto",
                 start: Callable = None, close: Callable = None,
                 finish_close: Callable = None, python: Optional[str] = None,
                 base_env: Optional[dict] = None, synchronous: bool = False,
                 on_change: Optional[Callable] = None, pull=None,
                 more_extra=(), hint=None, toasts: Optional[HU.Toasts] = None):
        import ipywidgets as W
        self.W = W
        self.repo = str(repo) if repo else ""
        self.gh = gh
        self._start = start or BP.start
        self._close = close or BP.close
        self._finish_close = finish_close or BP.finish_close_result
        self.python = python or sys.executable
        self.base_env = base_env
        self.synchronous = synchronous
        self.on_change = on_change
        self.toasts = toasts or HU.Toasts(timer=not synchronous)
        self.value = "main"
        self._all: list = [BL.BranchRow(ref="main", kind="main", mine=True)]
        self._rows: dict = {}
        self._shown: list = []
        self._refusals: dict = {}
        self._starting: dict = {}
        self._closing: dict = {}
        self._chips: dict = {}
        self._previews: list = []
        self._pending: list = []
        prefs = HU.load_prefs()

        # -- the switcher
        self.branch_btn = HU.icon_button("⎇ main ▾", "Branch to test (/ or Ctrl+K)",
                                         "delfin-btn", "delfin-branch-btn")
        self.search = W.Text(placeholder="Search branches, or type a tag or commit",
                             continuous_update=False)
        HU.classed(self.search, "delfin-branch-search", "delfin-debounced")
        self.reload_btn = HU.icon_button("⟳", "Reload: fetch from origin and the pull "
                                         "request states", "delfin-icon-btn")
        self.view = W.ToggleButtons(options=[v for v, _k in VIEWS], value="Mine",
                                    tooltips=_VIEW_TIPS)
        HU.classed(self.view, "delfin-views")
        head = HU.classed(W.HBox([self.search, self.reload_btn]), "delfin-pop-head")
        self.confirm_box = HU.classed(W.HBox([]), "delfin-confirm")
        self.confirm_box.layout.display = "none"
        self.new_open_btn = HU.icon_button("＋ New branch from main …",
                                           "Create a branch at main without switching",
                                           "delfin-btn", "delfin-brow-main")
        self.new_text = W.Text(placeholder="name, e.g. max/idea", continuous_update=False)
        HU.classed(self.new_text, "delfin-branch-search")
        self.new_btn = HU.icon_button("Create", "Create the branch and add it to your list",
                                      "delfin-btn", "is-primary")
        self.new_text.layout.display = "none"
        self.new_btn.layout.display = "none"
        self.new_row = HU.classed(W.HBox([self.new_open_btn, self.new_text, self.new_btn]),
                                  "delfin-brow")
        self.rows_box = HU.classed(W.VBox([self.new_row]), "delfin-branch-rows")
        self.pop_note = HU.classed(W.HTML(""), "delfin-pop-note")
        self.pop = HU.classed(W.VBox([head, self.view, self.confirm_box, self.rows_box,
                                      self.pop_note]), "delfin-branch-pop")
        self.pop.layout.display = "none"
        self.branch_wrap = HU.classed(W.VBox([self.branch_btn, self.pop]),
                                      "delfin-branch-wrap")

        # -- test, chips, more
        self.test_btn = HU.icon_button("▶ Test", SWITCH_NOTE, "delfin-btn", "is-primary")
        self.chips = HU.classed(W.HBox([]), "delfin-group")
        self.more_btn = W.ToggleButton(value=bool(prefs.get("more_open")),
                                       description="⋯ More",
                                       tooltip="Clean up and rarely used actions")
        HU.classed(self.more_btn, "delfin-btn")
        self.cleanup_btn = W.ToggleButton(value=False, description="🧹 Clean up …",
                                          tooltip="Delete merged branches and unused "
                                          "preview folders; every deletion can be undone")
        HU.classed(self.cleanup_btn, "delfin-btn")
        self.more_row = HU.classed(W.HBox([self.cleanup_btn, *more_extra]),
                                   "delfin-more-row")
        self.more_row.layout.display = "" if self.more_btn.value else "none"
        for b in more_extra:
            HU.classed(b, "delfin-btn")
        top = [self.branch_wrap, self.test_btn, self.chips]
        if pull is not None:
            HU.classed(pull, "delfin-btn")
            top.append(pull)
        top.append(self.more_btn)
        if hint is not None:
            top.append(hint)
        self.row = HU.classed(W.HBox(top), "delfin-group")

        self._build_cleanup()
        self.box = W.VBox([self.more_row, self.cleanup_box])

        self.branch_btn.on_click(lambda _b: self.toggle_menu())
        self.search.observe(lambda _c: self.render_rows(), names="value")
        _on_enter(self.search, self.search_submit)
        self.view.observe(lambda _c: self.render_rows(), names="value")
        self.reload_btn.on_click(lambda _b: self.reload())
        self.new_open_btn.on_click(lambda _b: self._open_new())
        self.new_btn.on_click(lambda _b: self.new_branch())
        _on_enter(self.new_text, self.new_branch)
        self.test_btn.on_click(lambda _b: self.test_selected())
        self.more_btn.observe(self._toggle_more, names="value")
        self.cleanup_btn.observe(self._toggle_cleanup, names="value")

    # -- plumbing ------------------------------------------------------------

    def _run(self, fn) -> None:
        if self.synchronous:
            fn()
        else:
            threading.Thread(target=fn, daemon=True).start()

    def say(self, text: str, *, level: str = "info", details: str = "",
            actions=None, key: str = "") -> None:
        self.toasts.show(text, level=level, details=details, actions=actions, key=key)

    def selected(self) -> str:
        return str(self.value or "")

    def _changed(self) -> None:
        if self.on_change:
            try:
                self.on_change()
            except Exception:
                pass

    # -- the list ------------------------------------------------------------

    def refresh(self, *, gh_refresh: bool = False, select: str = "") -> None:
        """Git first (fast), then again with gh in the background."""
        if select:
            self.value = select
        if not self.repo:
            self.render_rows()
            return
        try:
            self._all = BL.all_rows(self.repo, gh=None)
        except Exception as exc:
            self.say("The branch list could not be read.", level="error", details=str(exc))
        self._update_refusals()
        self.render_rows()
        self.refresh_previews()
        if self.gh:
            def _with_gh():
                try:
                    rows = BL.all_rows(self.repo, gh=self.gh, refresh=gh_refresh)
                except Exception:
                    return
                self._all = rows
                self.render_rows()
            self._run(_with_gh)

    def _update_refusals(self) -> None:
        names = [r.ref for r in self._all if r.kind == "local" and r.merged]
        try:
            self._refusals = BC.delete_refusals(self.repo, names) if names else {}
        except Exception:
            self._refusals = {}

    def _view_key(self) -> str:
        return dict(VIEWS).get(self.view.value, "mine")

    def render_rows(self) -> None:
        view = self._view_key()
        text = (self.search.value or "").strip()
        rows = BL.view_rows(self._all, view, text)
        if text and view in ("all", "remote") and self.repo:
            known = {r.ref for r in self._all}
            exact = [r for r in rows if r.ref == text or r.ref == f"origin/{text}"]
            if not exact:
                try:
                    cand = BL.lookup_row(self.repo, text, known)
                except Exception:
                    cand = None
                if cand is not None and (cand.kind != "lookup" or not rows):
                    rows = [cand] + rows
        more = max(0, len(rows) - MAX_ROWS)
        rows = rows[:MAX_ROWS]
        if self.repo:
            try:
                BL.fill_ahead(self.repo, [r for r in rows if r.ahead is None])
            except Exception:
                pass
        self._shown = rows
        widgets = [self._row_for(r).box for r in rows]
        if not rows:
            widgets.append(self.W.HTML('<span class="delfin-pop-note">No branch matches.</span>'))
        widgets.append(self.new_row)
        self.rows_box.children = tuple(widgets)
        self.pop_note.value = (f"{more} more — type to narrow the list." if more else
                               "Your list is personal; the shared checkout stays on main.")
        self.branch_btn.description = f"⎇ {self.value} ▾"
        self.test_btn.tooltip = (f"Open {self.value} in its own preview "
                                 "(the shared checkout stays on main)")

    def shown_refs(self) -> list:
        return [r.ref for r in self._shown]

    def row_widgets(self, ref: str) -> Optional[RowW]:
        for r in self._shown:
            if r.ref == ref:
                return self._row_for(r)
        return None

    def _row_for(self, r: BL.BranchRow) -> RowW:
        refusal = self._refusals.get(r.ref, "x") if (r.kind == "local" and r.merged) else "x"
        key = (r.ref, r.kind, r.pinned, r.hidden, r.mine, r.merged, r.ahead, r.pr,
               r.ci, refusal, r.ref == self.value, r.date, r.author)
        hit = self._rows.get(r.ref)
        if hit is not None and hit[0] == key:
            return hit[1]
        W = self.W
        main = HU.icon_button(f"{_row_icon(r)} {r.ref}", self._row_tip(r),
                              "delfin-btn", "delfin-brow-main")
        main.on_click(lambda _b, ref=r.ref: self.select(ref))
        meta = []
        if r.ahead:
            meta.append(f"+{int(r.ahead)}")
        if r.pr is not None:
            meta.append(f"#{int(r.pr)} {_ci_html(r.ci)}")
        if r.hidden:
            meta.append("hidden")
        if r.kind == "local" and r.merged:
            meta.append("in main")
        if r.kind == "lookup":
            meta.append("look up on origin")
        if r.date:
            meta.append(_esc(r.date) + (f" · {_esc(r.author)}" if r.author else ""))
        meta_w = HU.classed(W.HTML(" · ".join(meta)), "delfin-brow-meta")
        test = HU.icon_button("▶", f"Test {r.ref} in a preview", "delfin-icon-btn",
                              "delfin-brow-act")
        test.on_click(lambda _b, ref=r.ref: self.test_selected(ref))
        row = RowW(box=None, main=main, test=test)
        acts = [test]
        if r.kind not in ("main", "lookup"):
            if r.hidden:
                row.show = HU.icon_button("↺", "Show again in your list",
                                          "delfin-icon-btn")
                row.show.on_click(lambda _b, ref=r.ref: self.show_again(ref))
                acts.append(row.show)
            elif r.mine:
                row.hide = HU.icon_button("✕", "Hide from your list (git is not touched; "
                                          "undo from the message)", "delfin-icon-btn",
                                          "delfin-brow-act")
                row.hide.on_click(lambda _b, ref=r.ref: self.hide(ref))
                acts.append(row.hide)
        if r.kind != "main" and not r.mine and not r.hidden:
            row.pin = HU.icon_button("📌", "Add to your list" + (
                " (fetches it from origin)" if r.kind in ("remote", "lookup") else ""),
                "delfin-icon-btn")
            row.pin.on_click(lambda _b, rr=r: self.add(rr))
            acts.append(row.pin)
        if refusal == "":
            row.delete = HU.icon_button("🗑", "Delete this local branch (it is already "
                                        "in main; can be restored)", "delfin-icon-btn",
                                        "delfin-brow-act")
            row.delete.on_click(lambda _b, ref=r.ref: self.ask_delete(ref))
            acts.append(row.delete)
        box = HU.classed(W.HBox([main, meta_w, *acts]), "delfin-brow")
        if r.ref == self.value:
            box.add_class("is-selected")
        row.box = box
        self._rows[r.ref] = (key, row)
        return row

    @staticmethod
    def _row_tip(r) -> str:
        what = {"main": "main", "remote": "only on origin", "tag": "tag", "commit": "commit",
                "lookup": "not here yet: adding it fetches it from origin"}.get(r.kind, "branch")
        bits = [f"{r.ref} — {what}"]
        if r.ahead:
            bits.append(f"{r.ahead} commit(s) not in main")
        if r.pr is not None:
            bits.append(f"pull request #{r.pr}" + (f", checks {r.ci}" if r.ci else ""))
        return "; ".join(bits)

    # -- the menu ------------------------------------------------------------

    def open_menu(self) -> None:
        self.confirm_box.children = ()
        self.confirm_box.layout.display = "none"
        self.pop.layout.display = ""
        self.refresh()

    def close_menu(self) -> None:
        self.pop.layout.display = "none"
        self.confirm_box.children = ()
        self.confirm_box.layout.display = "none"
        self.new_text.layout.display = "none"
        self.new_btn.layout.display = "none"
        self.new_open_btn.layout.display = ""

    def toggle_menu(self) -> None:
        if self.pop.layout.display == "none":
            self.open_menu()
        else:
            self.close_menu()

    def select(self, ref: str) -> None:
        self.value = ref
        self.close_menu()
        self.render_rows()

    def search_submit(self) -> None:
        text = (self.search.value or "").strip().lower()
        rows = self._shown
        if not rows:
            return
        pick = next((r for r in rows if text and text in r.ref.lower()), rows[0])
        if pick.kind == "lookup":
            self.add(pick)
            return
        self.select(pick.ref)

    def reload(self) -> None:
        self.reload_btn.disabled = True
        self.reload_btn.description = "⟳ …"

        def _go():
            try:
                if self.repo:
                    BL._git(self.repo, "fetch", "--quiet", "--no-tags", "origin",
                            timeout=90.0)
                self.refresh(gh_refresh=True)
                self.say("Branch list updated.", level="ok", key="reload")
            except Exception as exc:
                self.say("Couldn't reload the branch list.", level="error", details=str(exc))
            finally:
                self.reload_btn.disabled = False
                self.reload_btn.description = "⟳"
        self._run(_go)

    # -- row actions -----------------------------------------------------------

    def hide(self, ref: str) -> None:
        was_pinned = ref in BL.load_prefs(self.repo)["pinned"]
        ok, msg = BL.hide(self.repo, ref)
        if not ok:
            self.say(msg, level="error")
            return

        def _undo():
            BL.unhide(self.repo, ref)
            if was_pinned:
                BL.pin(self.repo, ref)
            self.refresh()
        if self.value == ref:
            self.value = "main"
        self.refresh()
        self.say(f"Hid {ref} from your list.", level="ok", actions=[("Undo", _undo)])

    def show_again(self, ref: str) -> None:
        BL.unhide(self.repo, ref)
        BL.pin(self.repo, ref)
        self.refresh()

        def _undo():
            BL.hide(self.repo, ref)
            self.refresh()
        self.say(f"{ref} is back in your list.", level="ok", actions=[("Undo", _undo)])

    def add(self, r) -> None:
        def _go():
            if r.kind == "lookup":
                ok, msg, ref = BL.add_to_list(self.repo, r.ref)
                if not ok:
                    self.say(f"Couldn't find {r.ref} here or on origin.", level="error",
                             details=msg)
                    return
            else:
                ref = r.ref
                if r.kind == "remote" and ref.startswith("origin/"):
                    BL.fetch_remote_branch(self.repo, ref[len("origin/"):])
                BL.pin(self.repo, ref)
            self.refresh()

            def _undo():
                BL.unpin(self.repo, ref)
                self.refresh()
            self.say(f"Added {ref} to your list.", level="ok", actions=[("Undo", _undo)])
        self._run(_go)

    def ask_delete(self, ref: str) -> None:
        W = self.W
        yes = HU.icon_button("Delete", f"Delete {ref}", "delfin-btn", "is-danger", "is-filled")
        no = HU.icon_button("Cancel", "", "delfin-btn")

        def _cancel(_b=None):
            self.confirm_box.children = ()
            self.confirm_box.layout.display = "none"
        no.on_click(_cancel)
        yes.on_click(lambda _b: self._do_delete(ref))
        self.confirm_box.children = (W.HTML(
            f"Delete the local branch <b>{_esc(ref)}</b>? It is already in main; "
            "origin is not touched, and it can be restored."), yes, no)
        self.confirm_box.layout.display = ""

    def _do_delete(self, ref: str) -> None:
        self.confirm_box.children = ()
        self.confirm_box.layout.display = "none"
        ok, msg, report = BC.delete_branch(self.repo, ref, confirm_unmerged=False)
        if ok:
            self.say(f"Deleted {ref}.", level="ok",
                     details="Restore with:\n" + "\n".join(report.undo))
            if self.value == ref:
                self.value = "main"
        else:
            self.say(f"{ref} was not deleted.", level="error", details=msg)
        self.refresh()

    def _open_new(self) -> None:
        self.new_open_btn.layout.display = "none"
        self.new_text.layout.display = ""
        self.new_btn.layout.display = ""

    def new_branch(self) -> None:
        ok, msg, ref = BL.create_branch_from_main(self.repo, self.new_text.value)
        if not ok:
            self.say(msg, level="error")
            return
        self.new_text.value = ""
        self.close_menu()
        self.refresh(select=ref)
        self.say(f"Created {ref} from main and added it to your list.", level="ok",
                 actions=[("▶ Test", lambda: self.test_selected(ref))])

    # -- previews --------------------------------------------------------------

    def test_selected(self, ref: str = "") -> None:
        ref = ref or self.selected()
        if not ref:
            return
        self.value = ref
        self.close_menu()
        self.render_rows()
        self._starting[ref] = True
        self.render_chips()

        def _go():
            p = None
            try:
                ok, msg, p = self._start(self.repo, ref, python=self.python,
                                         base_env=self.base_env,
                                         progress=lambda _m: None)[:3]
            except Exception as exc:
                ok, msg = False, f"{type(exc).__name__}: {exc}"
            self._starting.pop(ref, None)
            self.refresh_previews()
            if ok and p is not None:
                self.say(f"Preview of {ref} is ready — open it from its chip.",
                         level="ok", details=msg)
            else:
                self.say(f"Couldn't start the preview of {ref}.", level="error",
                         details=msg, actions=[("Try again", lambda: self.test_selected(ref))])
            self._changed()
        self._run(_go)

    def refresh_previews(self) -> None:
        try:
            self._previews = BP.list_previews(self.repo)
        except Exception:
            self._previews = []
        try:
            self._pending = BP.pending_closes(self.repo)
        except Exception:
            self._pending = []
        self.render_chips()

    def _chip_state(self) -> tuple:
        return (tuple(sorted((p.slug, p.port) for p in self._previews)),
                tuple(sorted(p.slug for p in self._pending)))

    def poll_previews(self) -> bool:
        """Re-read the previews; redraw only when they changed (a preview
        closed from its own tab, or started from another window). True
        when the chips were redrawn."""
        if self._starting or self._closing:
            return False
        before = self._chip_state()
        try:
            self._previews = BP.list_previews(self.repo)
            self._pending = BP.pending_closes(self.repo)
        except Exception:
            return False
        if self._chip_state() == before:
            return False
        self.render_chips()
        return True

    def start_polling(self, every_s: float = 10.0) -> None:
        """Keep the chips current without a reload (a daemon thread)."""
        import time as _time

        def _loop():
            while True:
                _time.sleep(every_s)
                try:
                    self.poll_previews()
                except Exception:
                    pass
        threading.Thread(target=_loop, daemon=True, name="delfin-preview-chips").start()

    def chip_widgets(self, ref: str) -> Optional[ChipW]:
        return self._chips.get(ref)

    def render_chips(self) -> None:
        W = self.W
        chips: dict = {}
        for p in self._previews:
            if p.slug in self._closing:
                continue
            text = W.HTML(f'<span class="dot"></span>Preview: <b>{_esc(p.ref)}</b>')
            link = W.HTML(preview_link_html(p))
            close = HU.icon_button("✕", "Close this preview", "delfin-icon-btn")
            close.on_click(lambda _b, s=p.slug, ref=p.ref: self.close_preview(s, ref))
            box = HU.classed(W.HBox([text, link, close]), "delfin-chip")
            chips[p.ref] = ChipW(box, text, link, close)
        for ref in list(self._starting):
            if ref in chips:
                continue
            text = W.HTML(f'<span class="dot"></span>Starting <b>{_esc(ref)}</b> …')
            chips[ref] = ChipW(HU.classed(W.HBox([text]), "delfin-chip", "is-busy"), text)
        for _s, ref in list(self._closing.items()):
            text = W.HTML(f'<span class="dot"></span>Closing <b>{_esc(ref)}</b> …')
            chips[ref] = ChipW(HU.classed(W.HBox([text]), "delfin-chip", "is-busy"), text)
        for p in self._pending:
            if p.ref in chips:
                continue
            text = W.HTML(f'<span class="dot"></span><b>{_esc(p.ref)}</b> closed; its '
                          'folder has changes')
            rm = HU.icon_button("Remove", "Delete the preview folder with its changes",
                                "delfin-link-btn")
            keep = HU.icon_button("Keep", "Keep the preview folder", "delfin-link-btn")
            rm.on_click(lambda _b, s=p.slug: self.answer_close(s, True))
            keep.on_click(lambda _b, s=p.slug: self.answer_close(s, False))
            chips[p.ref] = ChipW(HU.classed(W.HBox([text, rm, keep]), "delfin-chip",
                                            "is-error"), text)
        self._chips = chips
        self.chips.children = tuple(c.box for c in chips.values())

    def close_preview(self, s: str, ref: str = "") -> None:
        ref = ref or s
        self._closing[s] = ref
        self.render_chips()

        def _go():
            try:
                r = self._close(s, progress=lambda _m: None)
            except Exception as exc:
                r = BP.CloseResult(False, "failed", f"{type(exc).__name__}: {exc}",
                                   "See what runs: python -m "
                                   "delfin.dashboard.branch_preview list")
            finally:
                self._closing.pop(s, None)
            try:
                self.refresh_previews()
            except Exception:
                self.render_chips()
            if r.state == "ask":
                self.say(f"Closed {ref}. Its folder has changes — remove it anyway?",
                         level="ask", details=r.details,
                         actions=[("Remove", lambda: self.answer_close(s, True)),
                                  ("Keep", lambda: self.answer_close(s, False))])
            elif r.ok:
                self.say(r.text, level="ok" if r.state == "closed" else "warn",
                         details=r.details)
            else:
                text = r.text if r.text.startswith("Couldn't") else (
                    f"Couldn't close the preview of {ref}: {r.text}")
                self.say(text, level="error", details=r.details,
                         actions=[("Try again", lambda: self.close_preview(s, ref))])
            self._changed()
        self._run(_go)

    def answer_close(self, s: str, discard: bool) -> None:
        def _go():
            try:
                r = self._finish_close(s, discard=discard)
            except Exception as exc:
                r = BP.CloseResult(False, "failed", f"{type(exc).__name__}: {exc}")
            self.refresh_previews()
            self.say(r.text if r.ok else f"Couldn't finish: {r.text}",
                     level="ok" if r.ok else "error", details=r.details)
        self._run(_go)

    # -- more ------------------------------------------------------------------

    def _toggle_more(self, change) -> None:
        on = bool(change.get("new"))
        self.more_row.layout.display = "" if on else "none"
        if not on and self.cleanup_btn.value:
            self.cleanup_btn.value = False
        HU.save_prefs(more_open=on)

    # -- clean up --------------------------------------------------------------

    def _build_cleanup(self) -> None:
        W = self.W
        self.cleanup_search = W.Text(placeholder="Search", continuous_update=False)
        HU.classed(self.cleanup_search, "delfin-branch-search", "delfin-debounced")
        self.select_toggle = HU.icon_button("Select none", "Tick or untick the rows shown "
                                            "(never ticks pinned branches or ones with "
                                            "their own commits)", "delfin-btn")
        self.cleanup_count = HU.classed(W.HTML(""), "delfin-count")
        self.include_unmerged = W.Checkbox(value=False, indent=False,
                                           description="Also branches with their own commits")
        self.include_unmerged.layout.width = "auto"
        self.stale_days = W.BoundedIntText(value=BC.DEFAULT_STALE_DAYS, min=1, max=3650,
                                           description="Unused folders after (days)",
                                           style={"description_width": "initial"})
        self.stale_days.layout.width = "220px"
        self.scan_btn = HU.icon_button("⟳", "List again", "delfin-icon-btn")
        self.cleanup_close = HU.icon_button("✕", "Close Clean up", "delfin-icon-btn")
        self.cleanup_list = HU.classed(W.VBox([]), "delfin-cleanup-rows")
        self.cleanup_kept = W.HTML("")
        self.cleanup_delete_btn = HU.icon_button("Delete selected", "", "delfin-btn",
                                                 "is-danger")
        self.cleanup_confirm = HU.classed(W.HBox([]), "delfin-confirm")
        self.cleanup_confirm.layout.display = "none"
        title = W.HTML('<span class="delfin-panel-title">Clean up</span> '
                       '<span class="delfin-pop-note">merged branches and unused preview '
                       'folders; every deletion can be restored</span>')
        self.cleanup_box = HU.classed(W.VBox([
            HU.classed(W.HBox([title, self.cleanup_close]), "delfin-pop-head"),
            HU.classed(W.HBox([self.cleanup_search, self.select_toggle, self.cleanup_count,
                               self.scan_btn]), "delfin-pop-head"),
            HU.classed(W.HBox([self.include_unmerged, self.stale_days]), "delfin-pop-head"),
            self.cleanup_list, self.cleanup_kept,
            HU.classed(W.HBox([self.cleanup_delete_btn]), "delfin-pop-head"),
            self.cleanup_confirm,
        ]), "delfin-panel")
        self.cleanup_box.layout.display = "none"
        self._checks: list = []
        self._cands: list = []
        self._pinned: set = set()
        self.cleanup_close.on_click(lambda _b: setattr(self.cleanup_btn, "value", False))
        self.cleanup_search.observe(lambda _c: self._filter_cleanup(), names="value")
        self.select_toggle.on_click(lambda _b: self.toggle_select())
        self.include_unmerged.observe(lambda _c: self.scan_cleanup(), names="value")
        self.scan_btn.on_click(lambda _b: self.scan_cleanup())
        self.cleanup_delete_btn.on_click(lambda _b: self.ask_cleanup())

    def _toggle_cleanup(self, change) -> None:
        on = bool(change.get("new"))
        self.cleanup_box.layout.display = "" if on else "none"
        if on:
            self.scan_cleanup()

    def scan_cleanup(self) -> None:
        self.cleanup_confirm.children = ()
        self.cleanup_confirm.layout.display = "none"
        self.cleanup_list.children = (self.W.HTML('<span class="delfin-pop-note">'
                                                  'Looking …</span>'),)

        def _go():
            try:
                cands, kept = BC.candidates(
                    self.repo, gh=self.gh, stale_days=int(self.stale_days.value),
                    protect_paths=BP.running_paths(self.repo),
                    include_unmerged=bool(self.include_unmerged.value))
            except Exception as exc:
                self.cleanup_list.children = ()
                self.say("Couldn't list what can be cleaned up.", level="error",
                         details=str(exc))
                return
            self._show_candidates(cands, kept)
        self._run(_go)

    @staticmethod
    def _label(c) -> str:
        if c.kind == "worktree":
            return f"Folder {Path(c.path).name} — {c.reason}"
        reason = (c.reason.replace(f"merged into {BL.BASE}", "already in main")
                  .replace("merged into main", "already in main")
                  .replace(f"not in {BL.BASE}", "not in main"))
        if c.unmerged and "not in" not in reason:
            reason += f"; {c.unmerged} commit(s) not in main"
        return f"Branch {c.name} — {reason}"

    def _show_candidates(self, cands, kept) -> None:
        W = self.W
        self._cands = cands
        self._checks = []
        try:
            self._pinned = set(BL.load_prefs(self.repo)["pinned"])
        except Exception:
            self._pinned = set()
        rows = []
        for c in cands:
            label = self._label(c) + (" — pinned" if self._is_pinned(c) else "")
            cb = W.Checkbox(value=self._eligible(c), description=label, indent=False,
                            tooltip=c.path or c.name)
            cb.layout.width = "auto"
            cb.observe(lambda _c: self._update_count(), names="value")
            self._checks.append(cb)
            rows.append(cb)
        if not rows:
            rows = [W.HTML('<span class="delfin-pop-note">Nothing to clean up.</span>')]
        self.cleanup_list.children = tuple(rows)
        if kept:
            items = "".join(
                f"<li>{_esc(Path(k.name).name if k.kind == 'worktree' else k.name)}"
                f" — {_esc(k.why)}</li>" for k in kept)
            self.cleanup_kept.value = (f'<details><summary class="delfin-pop-note">Kept '
                                       f'({len(kept)}): in use or with changes</summary>'
                                       f'<ul>{items}</ul></details>')
        else:
            self.cleanup_kept.value = ""
        self._filter_cleanup()

    def _is_pinned(self, c) -> bool:
        return c.kind == "branch" and c.name in self._pinned

    def _eligible(self, c) -> bool:
        """Ticked by default and by Select all: not pinned, no own commits."""
        return not c.unmerged and not self._is_pinned(c)

    def _filter_cleanup(self) -> None:
        needle = (self.cleanup_search.value or "").strip().lower()
        for c, cb in zip(self._cands, self._checks):
            hay = f"{c.name} {c.branch} {c.path}".lower()
            cb.layout.display = "" if (not needle or needle in hay) else "none"
        self._update_count()

    def visible_cands(self) -> list:
        return [c for c, cb in zip(self._cands, self._checks) if cb.layout.display != "none"]

    def _visible_state(self):
        vis = [(c, cb) for c, cb in zip(self._cands, self._checks)
               if cb.layout.display != "none"]
        elig = [cb for c, cb in vis if self._eligible(c)]
        all_on = bool(elig and all(cb.value for cb in elig)) or (
            not elig and any(cb.value for _c, cb in vis))
        return vis, elig, all_on

    def _update_count(self) -> None:
        n = sum(1 for cb in self._checks if cb.value)
        self.cleanup_count.value = f"{n} of {len(self._checks)} selected"
        _vis, _elig, all_on = self._visible_state()
        self.select_toggle.description = "Select none" if all_on else "Select all"
        self.cleanup_delete_btn.description = f"Delete selected ({n})"
        self.cleanup_delete_btn.disabled = n == 0

    def toggle_select(self) -> None:
        vis, elig, all_on = self._visible_state()
        if all_on:
            for _c, cb in vis:
                cb.value = False
        else:
            for cb in elig:
                cb.value = True
        self._update_count()

    def _selected_cands(self) -> list:
        return [c for c, cb in zip(self._cands, self._checks) if cb.value]

    def ask_cleanup(self) -> None:
        W = self.W
        sel = self._selected_cands()
        if not sel:
            return
        nb = sum(1 for c in sel if c.kind == "branch")
        nw = len(sel) - nb
        risky = [c for c in sel if c.unmerged]
        yes = HU.icon_button("Delete", "", "delfin-btn", "is-danger", "is-filled")
        no = HU.icon_button("Cancel", "", "delfin-btn")

        def _cancel(_b=None):
            self.cleanup_confirm.children = ()
            self.cleanup_confirm.layout.display = "none"
        no.on_click(_cancel)

        def _second(_b):
            all_btn = HU.icon_button("Delete those too", "", "delfin-btn", "is-danger",
                                     "is-filled")
            skip_btn = HU.icon_button("Skip those", "", "delfin-btn")
            all_btn.on_click(lambda _b2: self._do_cleanup(sel, True))
            skip_btn.on_click(lambda _b2: self._do_cleanup(
                [c for c in sel if not c.unmerged], False))
            names = ", ".join(c.name for c in risky[:6])
            self.cleanup_confirm.children = (W.HTML(
                f'<b>{len(risky)} selected branch(es) have commits that are not in main: '
                f'{_esc(names)}</b>'), all_btn, skip_btn, no)

        yes.on_click(_second if risky else (lambda _b: self._do_cleanup(sel, False)))
        parts = []
        if nb:
            parts.append(f"{nb} branch(es)")
        if nw:
            parts.append(f"{nw} folder(s)")
        self.cleanup_confirm.children = (W.HTML(
            f"Delete {' and '.join(parts)}? Nothing on origin changes; each one can be "
            "restored."), yes, no)
        self.cleanup_confirm.layout.display = ""

    def _do_cleanup(self, sel, confirm_unmerged: bool) -> None:
        self.cleanup_confirm.children = ()
        self.cleanup_confirm.layout.display = "none"
        report = BC.delete(self.repo, sel, confirm_unmerged=confirm_unmerged)
        details = []
        if report.undo:
            details.append("Restore with (also in "
                           + (report.log_path or str(BC._log_path())) + "):")
            details.extend(report.undo)
        if report.skipped:
            details.append("Not deleted:")
            details.extend(f"  {n} — {w}" for n, w in report.skipped)
        level = "ok" if not report.skipped else ("warn" if report.deleted else "error")
        text = f"Deleted {len(report.deleted)}." + (
            f" {len(report.skipped)} not deleted." if report.skipped else "")
        self.say(text, level=level, details="\n".join(details))
        self.refresh()
        self.scan_cleanup()


# -- inside a preview -----------------------------------------------------------

def preview_bar(toasts: Optional[HU.Toasts] = None, env=None):
    """The header of a preview dashboard: what runs here, and Close preview."""
    import ipywidgets as W
    env = os.environ if env is None else env
    try:
        info = BP.banner_info(env)
    except Exception as exc:
        info = {"branch": env.get(BP.PREVIEW_ENV_BRANCH, "?"), "sha": "", "code": "",
                "delfin_file": str(exc), "worktree": "", "inside": False}
    wrong = not info.get("inside")
    detail = (f'code from {_esc(info.get("code"))}<br>delfin.__file__ = '
              f'{_esc(info.get("delfin_file"))}')
    warn = " — this code is NOT from the preview's folder" if wrong else ""
    text = W.HTML(
        f'<div class="delfin-preview-banner" data-delfin-file="{_esc(info.get("delfin_file"))}">'
        f'<b>PREVIEW</b> · {_esc(info.get("branch"))} @ {_esc(info.get("sha"))}{warn}'
        f'<details><summary>Details</summary>{detail}</details></div>')
    close = HU.icon_button("Close preview", "Stop this preview and remove its folder "
                           "(asks in the main dashboard if it has changes)", "delfin-btn")
    bar = HU.classed(W.HBox([text, close]), "delfin-preview-bar")
    if wrong:
        bar.add_class("is-wrong")

    def _close(_b):
        ok, msg = close_this_preview(env)
        close.disabled = ok
        if ok:
            text.value = (f'<div class="delfin-preview-banner"><b>PREVIEW</b> · '
                          f'{_esc(info.get("branch"))} — closing … you can close this '
                          'tab.</div>')
        if toasts is not None:
            toasts.show(msg, level="info" if ok else "error")
    close.on_click(_close)
    return bar


def close_this_preview(env=None) -> tuple[bool, str]:
    """Ask a detached helper to close the preview this kernel runs in: the
    kernel belongs to the preview server's process group, so the stop must
    not run in it. The helper imports the shared checkout's code."""
    env = dict(os.environ if env is None else env)
    ref = env.get(BP.PREVIEW_ENV_BRANCH) or ""
    if not ref:
        return False, "This is not a preview."
    s = BP.slug(ref)
    code = env.get(BP.SHARED_CHECKOUT_ENV) or env.get(BP.PREVIEW_ENV_WORKTREE) or ""
    # Not tied to the preview's lifeline (lifeline.child_env): the helper
    # has to outlive the server it ends.
    child = {k: v for k, v in env.items()
             if not k.startswith(("DELFIN_PREVIEW_", "DELFIN_LIFELINE_"))}
    child["PYTHONPATH"] = code
    try:
        subprocess.Popen([sys.executable, "-m", "delfin.dashboard.branch_preview", "stop", s],
                         cwd=str(Path.home()), env=child, stdin=subprocess.DEVNULL,
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                         start_new_session=True)
    except OSError as exc:
        return False, f"Couldn't close this preview: {exc}"
    return True, "Closing this preview … the main dashboard shows the result."
