"""The header's branch controls: list, Test branch, previews, manage, clean up.

Input: the shared checkout (``repo``) and, for tests, the functions that
start and stop previews and the ``gh`` to ask. Output: :class:`BranchPanel`
with ``row`` (the controls that sit in the header line) and ``box`` (the
running previews, the manage and clean-up panels, the messages) for the
line below it.

Nothing here switches the shared checkout. Testing a branch opens a
preview (``branch_preview.start``) in a worktree and shows a link to it;
the link opens a tab only when the person clicks it. Deleting goes
through ``branch_cleanup`` with its rules, one confirmation, and a second
one for anything with commits not in ``origin/main``.

Long work (gh, starting a server, clean-up listing) runs in a thread and
reports into an HTML widget, so the kernel keeps answering the page.
"""

from __future__ import annotations

import html
import sys
import threading
from typing import Callable, Optional

from delfin.dashboard import branch_cleanup as BC
from delfin.dashboard import branch_list as BL
from delfin.dashboard import branch_preview as BP

SWITCH_NOTE = ("The shared checkout stays on main. To try a branch, pick it "
               "and press Test branch: it opens in its own worktree and server.")


def _esc(text) -> str:
    return html.escape(str(text or ""))


def preview_link_html(p) -> str:
    """A link that opens the preview in a new tab, on the host the page was
    loaded from (``localhost`` until the click says otherwise)."""
    path = f"{p.url_path}?token={p.token}"
    href = f"http://localhost:{p.port}{path}"
    js = (f"this.href=location.protocol+'//'+location.hostname+':{int(p.port)}'"
          f"+{_esc(repr(path))};")
    sha = (p.sha or "")[:8]
    return (f'<a class="delfin-preview-link" data-port="{int(p.port)}" '
            f'href="{_esc(href)}" target="_blank" rel="noopener noreferrer" '
            f'onclick="{js}" style="font-family:monospace; font-weight:600;">'
            f'Open preview: {_esc(p.ref)} @ {_esc(sha)} (port {int(p.port)})</a>'
            f'<span style="color:#666; font-size:11px; margin-left:8px;">'
            f'{_esc(p.worktree)} &middot; through an SSH tunnel, forward port '
            f'{int(p.port)} too</span>')


class BranchPanel:
    def __init__(self, repo, *, gh: Optional[str] = "auto",
                 start: Callable = None, stop: Callable = None,
                 finish_close: Callable = None, python: Optional[str] = None,
                 base_env: Optional[dict] = None, synchronous: bool = False,
                 on_change: Optional[Callable] = None):
        import ipywidgets as W
        self.W = W
        self.repo = str(repo) if repo else ""
        self.gh = gh
        self._start = start or BP.start
        self._stop = stop or BP.stop
        self._finish_close = finish_close or BP.finish_close
        self.python = python or sys.executable
        self.base_env = base_env
        self.synchronous = synchronous
        self.on_change = on_change
        self.entries: list = [BL.BranchEntry(ref="main", kind="main")]

        lay = W.Layout
        self.dropdown = W.Dropdown(options=[("main", "main")], value="main",
                                   description="Branch:", layout=lay(width="340px"),
                                   style={"description_width": "initial"})
        self.filter = W.Text(placeholder="filter", layout=lay(width="110px"),
                             continuous_update=True)
        self.refresh_btn = W.Button(description="↻", tooltip="Refresh the list",
                                    layout=lay(width="36px"))
        self.test_btn = W.Button(description="Test branch", button_style="success",
                                 tooltip=SWITCH_NOTE, layout=lay(width="110px"))
        self.manage_btn = W.ToggleButton(description="Manage list", value=False,
                                         layout=lay(width="105px"))
        self.cleanup_btn = W.ToggleButton(description="Clean up", value=False,
                                          layout=lay(width="90px"))
        self.row = W.HBox([self.dropdown, self.filter, self.refresh_btn,
                           self.test_btn, self.manage_btn, self.cleanup_btn],
                          layout=lay(align_items="center", gap="4px"))

        self.status = W.HTML("")
        self.previews_box = W.VBox([])
        # -- manage
        self.add_text = W.Text(placeholder="branch, origin/<branch>, tag or commit",
                               layout=lay(width="300px"))
        self.add_btn = W.Button(description="Add to list", layout=lay(width="110px"))
        self.new_text = W.Text(placeholder="new branch name, e.g. max/idea",
                               layout=lay(width="300px"))
        self.new_btn = W.Button(description="New branch from main",
                                layout=lay(width="170px"))
        self.hide_btn = W.Button(description="Hide from list",
                                 layout=lay(width="120px"))
        self.delete_btn = W.Button(description="Delete branch…", button_style="danger",
                                   layout=lay(width="125px"))
        self.hidden_dd = W.Dropdown(options=[], description="Hidden:",
                                    layout=lay(width="300px"),
                                    style={"description_width": "initial"})
        self.unhide_btn = W.Button(description="Show again", layout=lay(width="100px"))
        self.confirm_box = W.HBox([])
        self.manage_box = W.VBox([
            W.HTML(f'<span style="color:#555;">{_esc(SWITCH_NOTE)} Your list is '
                   'yours: kept in ~/.delfin, not in the repository.</span>'),
            W.HBox([self.add_text, self.add_btn]),
            W.HBox([self.new_text, self.new_btn]),
            W.HBox([W.HTML("Selected entry:"), self.hide_btn, self.delete_btn]),
            W.HBox([self.hidden_dd, self.unhide_btn]),
            self.confirm_box,
        ], layout=lay(display="none", border="1px solid #cfd8dc", padding="6px"))
        # -- clean up
        self.stale_days = W.BoundedIntText(value=BC.DEFAULT_STALE_DAYS, min=1,
                                           max=3650, description="stale after days:",
                                           style={"description_width": "initial"},
                                           layout=lay(width="190px"))
        self.scan_btn = W.Button(description="List again", layout=lay(width="100px"))
        self.cleanup_list = W.VBox([])
        self.cleanup_kept = W.HTML("")
        self.cleanup_delete_btn = W.Button(description="Delete selected",
                                           button_style="danger",
                                           layout=lay(width="140px"))
        self.cleanup_confirm = W.HBox([])
        self.cleanup_result = W.HTML("")
        self.cleanup_box = W.VBox([
            W.HBox([self.stale_days, self.scan_btn]),
            self.cleanup_list, self.cleanup_kept,
            W.HBox([self.cleanup_delete_btn]), self.cleanup_confirm,
            self.cleanup_result,
        ], layout=lay(display="none", border="1px solid #cfd8dc", padding="6px"))
        self._checks: list = []
        self._cands: list = []

        self.box = W.VBox([self.status, self.previews_box, self.manage_box,
                           self.cleanup_box])

        self.filter.observe(lambda _c: self._apply_filter(), names="value")
        self.refresh_btn.on_click(lambda _b: self.refresh(gh_refresh=True))
        self.test_btn.on_click(lambda _b: self.test_selected())
        self.manage_btn.observe(self._toggle_manage, names="value")
        self.cleanup_btn.observe(self._toggle_cleanup, names="value")
        self.add_btn.on_click(lambda _b: self.add())
        self.new_btn.on_click(lambda _b: self.new_branch())
        self.hide_btn.on_click(lambda _b: self.hide_selected())
        self.unhide_btn.on_click(lambda _b: self.unhide_selected())
        self.delete_btn.on_click(lambda _b: self.ask_delete_selected())
        self.scan_btn.on_click(lambda _b: self.scan_cleanup())
        self.cleanup_delete_btn.on_click(lambda _b: self.ask_cleanup())

    # -- plumbing ------------------------------------------------------------

    def _run(self, fn) -> None:
        if self.synchronous:
            fn()
        else:
            threading.Thread(target=fn, daemon=True).start()

    def say(self, text: str, *, color: str = "#37474f", pre: bool = False) -> None:
        body = (f'<pre style="margin:0; white-space:pre-wrap;">{_esc(text)}</pre>'
                if pre else _esc(text))
        self.status.value = (f'<div style="color:{color}; font-size:12px;">'
                             f'{body}</div>' if text else "")

    def selected(self) -> str:
        return str(self.dropdown.value or "")

    # -- the list ------------------------------------------------------------

    def refresh(self, *, gh_refresh: bool = False, select: str = "") -> None:
        """Git-only first (instant), then again with gh in the background."""
        if not self.repo:
            return
        keep = select or self.selected()
        try:
            self.entries = BL.branch_entries(self.repo, gh=None)
        except Exception as exc:
            self.say(f"Branch list unavailable: {exc}", color="#b71c1c")
        self._apply_filter(keep)
        self._refresh_hidden()
        self.refresh_previews()
        if self.gh:
            def _with_gh():
                try:
                    entries = BL.branch_entries(self.repo, gh=self.gh,
                                                refresh=gh_refresh)
                except Exception:
                    return
                self.entries = entries
                self._apply_filter(keep if keep else self.selected())
            self._run(_with_gh)

    def _apply_filter(self, keep: str = "") -> None:
        keep = keep or self.selected()
        shown = BL.filter_entries(self.entries, self.filter.value)
        opts = [(e.label, e.ref) for e in shown]
        if keep and keep not in [v for _l, v in opts]:
            match = next((e for e in self.entries if e.ref == keep), None)
            if match is not None:
                opts.append((match.label, match.ref))
        if list(self.dropdown.options) != opts:
            self.dropdown.options = opts
        values = [v for _l, v in opts]
        target = keep if keep in values else ("main" if "main" in values else
                                              (values[0] if values else None))
        if target is not None and self.dropdown.value != target:
            self.dropdown.value = target

    def _refresh_hidden(self) -> None:
        hidden = BL.hidden_entries(self.repo) if self.repo else []
        self.hidden_dd.options = [(h, h) for h in hidden]
        self.unhide_btn.disabled = not hidden

    # -- previews ------------------------------------------------------------

    def test_selected(self) -> None:
        ref = self.selected()
        if not ref:
            self.say("Pick a branch first.")
            return
        self.test_btn.disabled = True
        self.say(f"Opening a preview of {ref} ... (the shared checkout is not touched)")

        def _go():
            try:
                ok, msg = self._start(self.repo, ref, python=self.python,
                                      base_env=self.base_env,
                                      progress=lambda m: self.say(m))[:2]
            except Exception as exc:
                ok, msg = False, f"Error: {exc}"
            self.test_btn.disabled = False
            self.say(msg, color="#1b5e20" if ok else "#b71c1c", pre=not ok)
            self.refresh_previews()
            self._changed()
        self._run(_go)

    def refresh_previews(self) -> None:
        W = self.W
        rows = []
        try:
            previews = BP.list_previews(self.repo)
        except Exception:
            previews = []
        for p in previews:
            close = W.Button(description="Close preview", button_style="warning",
                             layout=W.Layout(width="120px"))
            close.on_click(lambda _b, s=p.slug: self.close_preview(s))
            rows.append(W.HBox([W.HTML(preview_link_html(p)), close],
                               layout=W.Layout(align_items="center", gap="8px")))
        self.previews_box.children = tuple(rows)

    def close_preview(self, s: str) -> None:
        self.say(f"Stopping the preview {s} ...")

        def _go():
            ok, msg, ask = self._stop(s)
            self.refresh_previews()
            self._changed()
            if not ask:
                self.say(msg, color="#1b5e20" if ok else "#b71c1c")
                return
            self.say(msg, color="#e65100", pre=True)
            W = self.W
            rm = W.Button(description="Remove it anyway", button_style="danger")
            keep = W.Button(description="Keep the worktree")

            def _answer(discard):
                self.confirm_box.children = ()
                self.previews_box.children = tuple(
                    c for c in self.previews_box.children if c is not row)
                ok2, msg2 = self._finish_close(s, discard=discard)
                self.say(msg2, color="#1b5e20" if ok2 else "#b71c1c")
            rm.on_click(lambda _b: _answer(True))
            keep.on_click(lambda _b: _answer(False))
            row = W.HBox([W.HTML(f"Worktree of {_esc(s)} has changes:"), rm, keep])
            self.previews_box.children = tuple(self.previews_box.children) + (row,)
        self._run(_go)

    def _changed(self) -> None:
        self.refresh()
        if self.on_change:
            try:
                self.on_change()
            except Exception:
                pass

    # -- manage --------------------------------------------------------------

    def _toggle_manage(self, change) -> None:
        self.manage_box.layout.display = "" if change.get("new") else "none"
        if change.get("new"):
            self._refresh_hidden()

    def add(self) -> None:
        ok, msg, ref = BL.add_to_list(self.repo, self.add_text.value)
        self.say(msg, color="#1b5e20" if ok else "#b71c1c")
        if ok:
            self.add_text.value = ""
            self.refresh(select=ref)

    def new_branch(self) -> None:
        ok, msg, ref = BL.create_branch_from_main(self.repo, self.new_text.value)
        if ok:
            self.new_text.value = ""
            self.refresh(select=ref)
            msg += " Press Test branch to open it."
        self.say(msg, color="#1b5e20" if ok else "#b71c1c")

    def hide_selected(self) -> None:
        ref = self.selected()
        ok, msg = BL.hide(self.repo, ref)
        self.say(msg, color="#1b5e20" if ok else "#b71c1c")
        if ok:
            self.refresh(select="main")

    def unhide_selected(self) -> None:
        ref = str(self.hidden_dd.value or "")
        if not ref:
            return
        BL.unhide(self.repo, ref)
        BL.pin(self.repo, ref)
        self.say(f"{ref} is back in your list.", color="#1b5e20")
        self.refresh(select=ref)

    def ask_delete_selected(self) -> None:
        W = self.W
        ref = self.selected()
        why = BC._branch_refusal(self.repo, ref)
        if why:
            self.say(f"{ref} cannot be deleted here: {why}.", color="#b71c1c")
            return
        n = BC.unmerged_count(self.repo, ref)
        sha = BL._out(self.repo, "rev-parse", "--short", f"refs/heads/{ref}")
        yes = W.Button(description="Yes, delete", button_style="danger")
        no = W.Button(description="Cancel")
        no.on_click(lambda _b: setattr(self.confirm_box, "children", ()))
        text = (f"Delete the local branch {_esc(ref)} ({_esc(sha)})? "
                "The remote is not touched; the undo command is logged.")

        def _second(_b):
            yes2 = W.Button(description=f"Delete with {n} unmerged commit(s)",
                            button_style="danger")
            yes2.on_click(lambda _b2: self._do_delete(ref, True))
            self.confirm_box.children = (W.HTML(
                f'<b style="color:#b71c1c;">{_esc(ref)} has {n} commit(s) not '
                'in origin/main.</b>'), yes2, no)

        yes.on_click(_second if n else (lambda _b: self._do_delete(ref, False)))
        self.confirm_box.children = (W.HTML(text), yes, no)

    def _do_delete(self, ref: str, confirm_unmerged: bool) -> None:
        self.confirm_box.children = ()
        ok, msg, _r = BC.delete_branch(self.repo, ref,
                                       confirm_unmerged=confirm_unmerged)
        self.say(msg, color="#1b5e20" if ok else "#b71c1c")
        if ok:
            self.refresh(select="main")

    # -- clean up ------------------------------------------------------------

    def _toggle_cleanup(self, change) -> None:
        self.cleanup_box.layout.display = "" if change.get("new") else "none"
        if change.get("new"):
            self.scan_cleanup()

    def scan_cleanup(self) -> None:
        self.cleanup_result.value = ""
        self.cleanup_confirm.children = ()
        self.cleanup_list.children = (self.W.HTML("Listing ..."),)

        def _go():
            try:
                cands, kept = BC.candidates(
                    self.repo, gh=self.gh, stale_days=int(self.stale_days.value),
                    protect_paths=BP.running_paths(self.repo))
            except Exception as exc:
                self.cleanup_list.children = (self.W.HTML(f"Error: {_esc(exc)}"),)
                return
            self._show_candidates(cands, kept)
        self._run(_go)

    def _show_candidates(self, cands, kept) -> None:
        W = self.W
        self._cands = cands
        self._checks = []
        rows = []
        try:
            pinned = set(BL.load_prefs(self.repo)["pinned"])
        except Exception:
            pinned = set()
        for c in cands:
            mine = c.kind == "branch" and c.name in pinned
            label = (f"worktree {c.path} — {c.reason}" if c.kind == "worktree"
                     else f"branch {c.name} ({c.sha[:8]}) — {c.reason}"
                     + (f" — {c.unmerged} commit(s) NOT in origin/main"
                        if c.unmerged else "")
                     + (f" — after its worktree {c.after}" if c.after else "")
                     + (" — pinned in your list" if mine else ""))
            cb = W.Checkbox(value=not c.unmerged and not mine, description=label,
                            indent=False,
                            layout=W.Layout(width="auto"))
            self._checks.append(cb)
            rows.append(cb)
        if not rows:
            rows = [W.HTML("Nothing to clean up.")]
        self.cleanup_list.children = tuple(rows)
        self.cleanup_delete_btn.disabled = not cands
        if kept:
            items = "".join(f"<li>{_esc(k.kind)} {_esc(k.name)} — {_esc(k.why)}</li>"
                            for k in kept)
            self.cleanup_kept.value = (
                f"<details><summary>Kept ({len(kept)})</summary><ul>{items}</ul>"
                "</details>")
        else:
            self.cleanup_kept.value = ""

    def _selected_cands(self) -> list:
        return [c for c, cb in zip(self._cands, self._checks) if cb.value]

    def ask_cleanup(self) -> None:
        W = self.W
        sel = self._selected_cands()
        if not sel:
            self.cleanup_result.value = "Nothing selected."
            return
        nb = sum(1 for c in sel if c.kind == "branch")
        nw = len(sel) - nb
        risky = [c for c in sel if c.unmerged]
        yes = W.Button(description="Yes, delete", button_style="danger")
        no = W.Button(description="Cancel")
        no.on_click(lambda _b: setattr(self.cleanup_confirm, "children", ()))

        def _second(_b):
            all_btn = W.Button(description="Delete those too", button_style="danger")
            skip_btn = W.Button(description="Skip those")
            all_btn.on_click(lambda _b2: self._do_cleanup(sel, True))
            skip_btn.on_click(lambda _b2: self._do_cleanup(
                [c for c in sel if not c.unmerged], False))
            names = ", ".join(c.name for c in risky[:6])
            self.cleanup_confirm.children = (W.HTML(
                f'<b style="color:#b71c1c;">{len(risky)} selected branch(es) have '
                f'commits not in origin/main: {_esc(names)}</b>'), all_btn,
                skip_btn, no)

        yes.on_click(_second if risky else (lambda _b: self._do_cleanup(sel, False)))
        self.cleanup_confirm.children = (W.HTML(
            f"Delete {nb} local branch(es) and remove {nw} worktree(s)? Nothing "
            "on the remote is touched; every deletion is logged with its undo "
            "command."), yes, no)

    def _do_cleanup(self, sel, confirm_unmerged: bool) -> None:
        self.cleanup_confirm.children = ()
        report = BC.delete(self.repo, sel, confirm_unmerged=confirm_unmerged)
        parts = [f"<b>Deleted {len(report.deleted)}.</b> "]
        if report.undo:
            parts.append("Undo (also in " + _esc(report.log_path or str(BC._log_path()))
                         + "):<pre style='white-space:pre-wrap;'>"
                         + _esc("\n".join(report.undo)) + "</pre>")
        if report.skipped:
            parts.append("Not deleted:<ul>" + "".join(
                f"<li>{_esc(n)} — {_esc(w)}</li>" for n, w in report.skipped)
                + "</ul>")
        self.cleanup_result.value = "".join(parts)
        self.refresh()
        self.scan_cleanup_keep_result()

    def scan_cleanup_keep_result(self) -> None:
        result = self.cleanup_result.value
        self.scan_cleanup()
        self.cleanup_result.value = result
