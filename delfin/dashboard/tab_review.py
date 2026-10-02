"""The Review tab: look at every frame of a multi-frame XYZ and rate it.

The ratings themselves -- loading, order, saving, resuming -- live in
:mod:`delfin.review`; this is only the screen on top of it.  The 3D view is
the dashboard's own (:func:`molecule_viewer.render_xyz_in_output`), so the
viewer settings apply here exactly as they do in the Calculations tab.

What is on the screen
---------------------

* Top: what to open (a multi-frame XYZ or a folder of them), the review
  name, the blinded switch with its seed, an optional findings file.
* Left: the files (normal mode) and the frames of the current file, rated
  ones marked; in the blinded mode only anonymous numbers.
* Middle: the structure, its position ("frame 3/12"), Prev / Next / Next
  unrated, Pass / Block, the categories and the comment on this frame.
* Right: progress, the SMILES and findings for this frame (normal mode
  only), a comment on the whole file and one on the whole review.

Keys (while the tab is visible and no text field has the focus): Left/Right
move, N jumps to the next unrated frame, P passes, B blocks, 1-9 toggle the
categories in their listed order.  The browser hands each key to a hidden
text field, the same bridge the Ketcher panel uses, so every key also works
headless through :meth:`ReviewPanel.handle_key`.
"""

from __future__ import annotations

import html
from pathlib import Path
from typing import Any

from delfin import review as rv

try:
    import ipywidgets as widgets
    HAS_WIDGETS = True
except ImportError:                                  # pragma: no cover
    widgets = None
    HAS_WIDGETS = False

__all__ = ['ReviewPanel', 'create_tab']

PANEL_CLASS = 'delfin-review-panel'
KEYS_CLASS = 'delfin-review-keys'
FRAME_WINDOW = 100
_MARK = {'pass': '✓', 'block': '✗', None: '·'}


def _keys_js() -> str:
    """Hand arrow keys, P, B, N and digits to the hidden field while the panel is shown."""
    return (
        "(function(){\n"
        "  if(window.__delfinReviewKeys) return;\n"
        "  window.__delfinReviewKeys=true;\n"
        "  var WANT={ArrowLeft:'left',ArrowRight:'right',p:'p',P:'p',b:'b',B:'b',n:'n',N:'n'};\n"
        "  document.addEventListener('keydown',function(ev){\n"
        "    if(ev.ctrlKey||ev.metaKey||ev.altKey) return;\n"
        "    var el=ev.target, tag=(el&&el.tagName||'').toUpperCase();\n"
        "    if(tag==='INPUT'||tag==='TEXTAREA'||tag==='SELECT'||(el&&el.isContentEditable)) return;\n"
        "    var key=WANT[ev.key]||(/^[1-9]$/.test(ev.key)?ev.key:null);\n"
        "    if(!key) return;\n"
        "    var panel=document.querySelector('." + PANEL_CLASS + "');\n"
        "    if(!panel||panel.offsetParent===null) return;\n"
        "    var box=panel.querySelector('." + KEYS_CLASS + "');\n"
        "    var input=box&&box.querySelector('input, textarea');\n"
        "    if(!input) return;\n"
        "    ev.preventDefault();\n"
        "    var setter=Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype,'value');\n"
        "    var line=Date.now()+' '+key;\n"
        "    if(setter&&setter.set) setter.set.call(input,line); else input.value=line;\n"
        "    input.dispatchEvent(new Event('input',{bubbles:true}));\n"
        "    input.dispatchEvent(new Event('change',{bubbles:true}));\n"
        "  },true);\n"
        "})();\n"
    )


class ReviewPanel:
    """Widgets and behaviour of the Review tab."""

    def __init__(self, start_dir: Any = None, ctx: Any = None):
        self.ctx = ctx
        self.session: rv.ReviewSession | None = None
        self._busy = False
        start = str(start_dir or Path.home())
        W, L = widgets, widgets.Layout

        # -- what to open --------------------------------------------------
        self.path_input = W.Text(value=start, description='Open:',
                                 placeholder='multi-frame .xyz or a folder of .xyz files',
                                 layout=L(width='420px'))
        self.name_input = W.Text(value='', description='Name:',
                                 placeholder='default: folder or file name',
                                 layout=L(width='260px'))
        self.open_btn = W.Button(description='Open', button_style='info',
                                 layout=L(width='80px'))
        self.blind_cb = W.Checkbox(value=False, description='Blinded', indent=False,
                                   layout=L(width='90px'))
        self.seed_input = W.IntText(value=0, description='Seed:', layout=L(width='140px'))
        self.findings_input = W.Text(value='', description='Findings:',
                                     placeholder='optional JSON; default findings.json',
                                     layout=L(width='340px'))

        # -- left: files and frames ----------------------------------------
        self.file_list = W.Select(options=[], rows=14, layout=L(width='100%'))
        self.frame_list = W.Select(options=[], rows=14, layout=L(width='100%'))
        self.left_title = W.HTML('<b>Files</b>')

        # -- middle: viewer and rating -------------------------------------
        self.viewer = W.Output(layout=L(width='100%', min_height='360px',
                                        border='1px solid #ddd'))
        self.position_html = W.HTML('')
        self.prev_btn = W.Button(description='◀ Prev', layout=L(width='90px'))
        self.next_btn = W.Button(description='Next ▶', layout=L(width='90px'))
        self.unrated_btn = W.Button(description='Next unrated', layout=L(width='120px'))
        self.pass_btn = W.Button(description='Pass (P)', button_style='success',
                                 layout=L(width='110px'))
        self.block_btn = W.Button(description='Block (B)', button_style='danger',
                                  layout=L(width='110px'))
        self.multi_cb = W.Checkbox(value=True, description='Several categories',
                                   indent=False, layout=L(width='170px'))
        self.category_box = W.HBox([], layout=L(flex_flow='row wrap', gap='4px'))
        self.category_buttons: list[Any] = []
        self.note_input = W.Text(value='', description='Comment:',
                                 placeholder='your comment on this frame',
                                 continuous_update=False, layout=L(width='100%'))

        # -- right: progress, context, comments ----------------------------
        self.progress_html = W.HTML('')
        self.context_html = W.HTML('')
        self.file_comment = W.Textarea(value='', placeholder='comment on this file',
                                       continuous_update=False,
                                       layout=L(width='100%', height='70px'))
        self.session_comment = W.Textarea(value='', placeholder='comment on the whole review',
                                          continuous_update=False,
                                          layout=L(width='100%', height='70px'))
        self.categories_input = W.Textarea(value='\n'.join(rv.DEFAULT_CATEGORIES),
                                           layout=L(width='100%', height='130px'))
        self.categories_btn = W.Button(description='Apply categories',
                                       layout=L(width='150px'))
        self.export_btn = W.Button(description='Export CSV', layout=L(width='110px'))
        self.status_html = W.HTML('')
        self.keys_input = W.Text(value='', layout=L(display='none'))
        self.keys_input.add_class(KEYS_CLASS)

        self._build_categories(rv.DEFAULT_CATEGORIES)
        self._wire()
        self.widget = self._layout()
        self.widget.add_class(PANEL_CLASS)
        self._set_enabled(False)
        self._install_keys()

    # -- layout ------------------------------------------------------------
    def _layout(self):
        W, L = widgets, widgets.Layout
        top = W.VBox([
            W.HBox([self.path_input, self.open_btn, self.name_input],
                   layout=L(flex_flow='row wrap', gap='6px')),
            W.HBox([self.blind_cb, self.seed_input, self.findings_input],
                   layout=L(flex_flow='row wrap', gap='6px')),
        ])
        self.files_box = W.VBox([self.left_title, self.file_list])
        left = W.VBox([self.files_box, W.HTML('<b>Frames</b>'), self.frame_list],
                      layout=L(width='22%', min_width='200px'))
        middle = W.VBox([
            self.position_html,
            self.viewer,
            W.HBox([self.prev_btn, self.next_btn, self.unrated_btn,
                    self.pass_btn, self.block_btn],
                   layout=L(flex_flow='row wrap', gap='6px')),
            W.HBox([W.HTML('<b>Categories</b>'), self.multi_cb], layout=L(gap='10px')),
            self.category_box,
            self.note_input,
        ], layout=L(width='50%', min_width='340px'))
        self.file_comment_box = W.VBox([W.HTML('<b>Comment on this file</b>'),
                                        self.file_comment])
        right = W.VBox([
            self.progress_html,
            self.context_html,
            self.file_comment_box,
            W.HTML('<b>Comment on this review</b>'), self.session_comment,
            W.HTML('<b>Categories</b> (one per line)'), self.categories_input,
            W.HBox([self.categories_btn, self.export_btn], layout=L(gap='6px')),
        ], layout=L(width='28%', min_width='220px'))
        return W.VBox([
            top,
            W.HBox([left, middle, right], layout=L(gap='12px', width='100%',
                                                   flex_flow='row wrap')),
            self.status_html,
            self.keys_input,
        ], layout=L(width='100%'))

    def _wire(self):
        self.open_btn.on_click(lambda _b: self.open())
        self.prev_btn.on_click(lambda _b: self._move('prev'))
        self.next_btn.on_click(lambda _b: self._move('next'))
        self.unrated_btn.on_click(lambda _b: self._move('unrated'))
        self.pass_btn.on_click(lambda _b: self.rate('pass'))
        self.block_btn.on_click(lambda _b: self.rate('block'))
        self.export_btn.on_click(lambda _b: self.export())
        self.categories_btn.on_click(lambda _b: self.apply_categories())
        self.blind_cb.observe(self._on_mode, names='value')
        self.seed_input.observe(self._on_mode, names='value')
        self.file_list.observe(self._on_file_pick, names='value')
        self.frame_list.observe(self._on_frame_pick, names='value')
        self.note_input.observe(self._on_note, names='value')
        self.file_comment.observe(self._on_file_comment, names='value')
        self.session_comment.observe(self._on_session_comment, names='value')
        self.keys_input.observe(self._on_key, names='value')

    def _install_keys(self):
        if self.ctx is None or getattr(self.ctx, 'js_output', None) is None:
            return
        try:
            from delfin.dashboard.helpers import _append_js
            _append_js(self.ctx, _keys_js())
        except Exception:                            # pragma: no cover
            pass

    def _build_categories(self, categories):
        self.category_buttons = []
        for i, cat in enumerate(categories):
            hint = f' ({i + 1})' if i < 9 else ''
            btn = widgets.ToggleButton(value=False, description=f'{cat}{hint}',
                                       tooltip=cat, layout=widgets.Layout(width='auto'))
            btn._review_category = cat
            btn.observe(self._on_category, names='value')
            self.category_buttons.append(btn)
        self.category_box.children = tuple(self.category_buttons)

    def _set_enabled(self, on: bool):
        for w in (self.prev_btn, self.next_btn, self.unrated_btn, self.pass_btn,
                  self.block_btn, self.export_btn, self.note_input,
                  self.file_comment, self.session_comment, self.categories_btn,
                  *self.category_buttons):
            w.disabled = not on

    def _say(self, text: str, error: bool = False):
        color = '#b00020' if error else '#555'
        self.status_html.value = f'<span style="color:{color}">{html.escape(text)}</span>'

    # -- opening -------------------------------------------------------------
    def open(self, path: str | None = None) -> rv.ReviewSession | None:
        if path is not None:
            self.path_input.value = str(path)
        raw = self.path_input.value.strip()
        if not raw:
            self._say('Give a multi-frame .xyz file or a folder.', error=True)
            return None
        findings = self.findings_input.value.strip() or None
        try:
            session = rv.ReviewSession.open(
                raw, name=self.name_input.value, blind=self.blind_cb.value,
                seed=int(self.seed_input.value or 0), findings_path=findings)
        except (OSError, ValueError) as exc:
            self._say(f'Cannot open {raw}: {exc}', error=True)
            return None
        if not len(session):
            self._say(f'No XYZ frames in {raw}.', error=True)
            return None
        self.session = session
        self._busy = True
        try:
            self.categories_input.value = '\n'.join(session.categories)
            self._build_categories(session.categories)
            self.session_comment.value = session.data.get('comment', '')
        finally:
            self._busy = False
        session.save()
        # Resume at the first frame still to be rated.
        if session.record(*session.order[0]) and session.next_unrated() is None:
            session.go(0)
        self._set_enabled(True)
        self.refresh(full=True)
        msg = f'Review file: {session.review_path}'
        if session.warnings:
            msg += '  |  ' + '; '.join(session.warnings)
        self._say(msg)
        return session

    # -- moving and rating ------------------------------------------------------
    def _move(self, where: str):
        s = self.session
        if s is None:
            return
        if where == 'prev':
            s.prev()
        elif where == 'next':
            s.next()
        elif where == 'unrated' and s.next_unrated() is None:
            self._say('Every frame is rated.')
        self.refresh()

    def selected_categories(self) -> list[str]:
        return [b._review_category for b in self.category_buttons if b.value]

    def rate(self, verdict: str):
        s = self.session
        if s is None:
            return None
        rec = s.rate(verdict, self.selected_categories(), self.note_input.value)
        self.refresh()
        p = s.progress()
        if p['rated'] == p['total']:
            self._say(f'All {p["total"]} frames rated. Review file: {s.review_path}')
        return rec

    def export(self) -> Path | None:
        if self.session is None:
            return None
        out = self.session.export_csv()
        self._say(f'CSV written: {out}')
        return out

    def apply_categories(self):
        cats = [c.strip() for c in self.categories_input.value.splitlines() if c.strip()]
        if not cats:
            self._say('The category list is empty.', error=True)
            return
        self._build_categories(cats)
        if self.session is not None:
            self.session.set_categories(cats)
            self.refresh()

    def handle_key(self, key: str):
        key = (key or '').strip().lower()
        if self.session is None:
            return
        if key == 'left':
            self._move('prev')
        elif key == 'right':
            self._move('next')
        elif key == 'n':
            self._move('unrated')
        elif key == 'p':
            self.rate('pass')
        elif key == 'b':
            self.rate('block')
        elif key.isdigit() and 1 <= int(key) <= len(self.category_buttons):
            btn = self.category_buttons[int(key) - 1]
            btn.value = not btn.value

    # -- observers ---------------------------------------------------------------
    def _on_key(self, change):
        value = str(change.get('new') or '')
        if not value:
            return
        self.handle_key(value.split(' ', 1)[-1])

    def _on_mode(self, _change):
        if self.session is None or self._busy:
            return
        self.session.set_mode(self.blind_cb.value, int(self.seed_input.value or 0))
        self.session.save()
        self.refresh(full=True)

    def _on_file_pick(self, change):
        if self._busy or self.session is None or not change.get('new'):
            return
        self.session.go_to_file(change['new'])
        self.refresh()

    def _on_frame_pick(self, change):
        if self._busy or self.session is None or change.get('new') is None:
            return
        self.session.go(int(change['new']))
        self.refresh()

    def _on_category(self, change):
        if self._busy or not change.get('new') or self.multi_cb.value:
            return
        self._busy = True
        try:
            for b in self.category_buttons:
                if b is not change['owner']:
                    b.value = False
        finally:
            self._busy = False

    def _on_note(self, change):
        """A changed comment on a frame that is already rated is stored at once."""
        s = self.session
        if self._busy or s is None:
            return
        frame = s.current
        rec = s.record(frame.file, frame.index)
        if rec is not None and rec.get('note', '') != change['new'].strip():
            rec['note'] = change['new'].strip()
            rec['note_updated'] = rv._now()
            s.save()

    def _on_file_comment(self, change):
        if self._busy or self.session is None or self.session.blind:
            return
        self.session.set_file_comment(self.session.current.file, change['new'])
        self.refresh_lists()

    def _on_session_comment(self, change):
        if self._busy or self.session is None:
            return
        self.session.set_session_comment(change['new'])

    # -- drawing -------------------------------------------------------------------
    def refresh(self, full: bool = False):
        s = self.session
        if s is None:
            return
        frame = s.current
        rec = s.record(frame.file, frame.index)
        self._busy = True
        try:
            chosen = set(rec.get('categories', [])) if rec else set()
            for b in self.category_buttons:
                b.value = b._review_category in chosen
            self.note_input.value = rec.get('note', '') if rec else ''
            meta = s.data['files'].get(frame.file, {})
            self.file_comment.value = '' if s.blind else meta.get('comment', '')
        finally:
            self._busy = False
        self.refresh_lists(full=full)
        self._draw_position(frame, rec)
        self._draw_progress()
        self._draw_context(frame)
        self._draw_viewer(frame)

    def refresh_lists(self, full: bool = False):
        s = self.session
        if s is None:
            return
        frame = s.current
        self._busy = True
        try:
            self.files_box.layout.display = 'none' if s.blind else ''
            self.file_comment_box.layout.display = 'none' if s.blind else ''
            if not s.blind:
                files = list(s.data['files'])
                present = {f.file for f in s.frames}
                opts = []
                for fname in files:
                    if fname not in present:
                        continue
                    n = s.data['files'][fname].get('n_frames', 0)
                    done = sum(1 for i in range(n) if s.record(fname, i))
                    mark = '✓' if n and done == n else ' '
                    opts.append((f'{mark} {fname}  {done}/{n}', fname))
                self.file_list.options = opts
                self.file_list.value = frame.file
                positions = [p for p, (f, _) in enumerate(s.order) if f == frame.file]
            else:
                positions = list(range(len(s)))
            # Thousands of options are rebuilt after every rating, so the list
            # shows a window around the current frame.
            if len(positions) > 2 * FRAME_WINDOW:
                at = positions.index(s.position)
                lo = max(0, min(at - FRAME_WINDOW, len(positions) - 2 * FRAME_WINDOW))
                positions = positions[lo:lo + 2 * FRAME_WINDOW]
            opts = []
            for pos in positions:
                fr = s.frame_at(pos)
                r = s.record(fr.file, fr.index)
                mark = _MARK[r['verdict'] if r else None]
                text = (f'{mark} #{pos + 1}' if s.blind
                        else f'{mark} frame {fr.index + 1}  {fr.label}')
                opts.append((text, pos))
            self.frame_list.options = opts
            self.frame_list.value = s.position
        finally:
            self._busy = False

    def _draw_position(self, frame, rec):
        s = self.session
        if s.blind:
            where = f'structure #{s.position + 1} of {len(s)}'
        else:
            n = s.data['files'].get(frame.file, {}).get('n_frames', 0)
            where = (f'<b>{html.escape(frame.file)}</b> &nbsp; frame {frame.index + 1}/{n}'
                     f' &nbsp; <i>{html.escape(frame.label)}</i>')
        if rec:
            color = '#1b7f3b' if rec['verdict'] == 'pass' else '#b00020'
            cats = ', '.join(rec.get('categories', []))
            state = (f'<span style="color:{color}"><b>{rec["verdict"]}</b></span>'
                     + (f' ({html.escape(cats)})' if cats else ''))
        else:
            state = '<span style="color:#888">not rated</span>'
        self.position_html.value = f'{where} &nbsp;&mdash;&nbsp; {state}'

    def _draw_progress(self):
        p = self.session.progress()
        pct = 100.0 * p['rated'] / p['total'] if p['total'] else 0.0
        self.progress_html.value = (
            f'<b>Progress</b> {p["rated"]} / {p["total"]} ({pct:.0f} %)<br>'
            f'<span style="color:#1b7f3b">pass {p["pass"]}</span> &nbsp; '
            f'<span style="color:#b00020">block {p["block"]}</span>'
            f'<div style="background:#eee;height:6px;margin-top:4px">'
            f'<div style="background:#4a90d9;height:6px;width:{pct:.1f}%"></div></div>'
        )

    def _draw_context(self, frame):
        s = self.session
        if s.blind:
            self.context_html.value = ''
            return
        parts = []
        meta = s.data['files'].get(frame.file, {})
        smiles = (meta.get('index') or {}).get('smiles')
        if smiles:
            parts.append(f'<b>SMILES</b><br><code style="word-break:break-all">'
                         f'{html.escape(smiles)}</code>')
        findings = s.findings_for(frame)
        if findings:
            items = ''.join(f'<li>{html.escape(t)}</li>' for t in findings)
            parts.append(f'<b>Findings</b><ul style="margin:2px 0 2px 18px">{items}</ul>')
        self.context_html.value = '<br>'.join(parts)

    def _draw_viewer(self, frame):
        comment = '' if self.session.blind else frame.comment
        try:
            from delfin.dashboard.molecule_viewer import render_xyz_in_output
            render_xyz_in_output(self.viewer, frame.xyz(comment))
        except Exception as exc:                     # pragma: no cover - no browser
            with self.viewer:
                print(f'Viewer unavailable: {exc}')


def create_tab(ctx) -> Any:
    """Build the Review tab.  Returns ``(tab_widget, refs_dict)``."""
    start = getattr(ctx, 'calc_dir', None) or Path.home()
    panel = ReviewPanel(start, ctx=ctx)
    box = widgets.VBox([panel.widget], layout=widgets.Layout(width='100%', padding='10px'))
    return box, {'review_panel': panel}


try:                                                 # pragma: no cover
    from delfin.dashboard.tab_registry import register_tab
    register_tab('review', 'Review', create_tab, order=8150)
except Exception:                                    # pragma: no cover
    pass
