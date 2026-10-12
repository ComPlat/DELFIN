"""The Review tab: look at every frame of a multi-frame XYZ and rate it.

The ratings themselves -- loading, order, saving, resuming -- live in
:mod:`delfin.review`; this is only the screen on top of it.  The 3D view is
the dashboard's own (:func:`molecule_viewer.render_xyz_in_output`), so the
viewer settings apply here exactly as they do in the Calculations tab.

What is on the screen
---------------------

* **Open a review** (a section that closes once a review is open): what to
  open (a multi-frame XYZ or a folder of them), the review name, an
  optional findings file, the blinded switch with its seed and the list of
  categories.
* A bar across the top: show/hide the list, Prev / Next / Next unrated,
  which frame is shown, and the progress (rated / total, a bar, pass and
  block counts).
* Left, narrow and foldable: the files (normal mode) and the frames, rated
  ones marked; in the blinded mode only anonymous numbers.
* Middle, large: the structure.
* Right: the rating -- a large Pass and a large Block button, the
  categories as toggles, the comment on this frame, a line that says when
  the rating was saved and a short key legend; below that the SMILES and
  findings (normal mode only).
* **More**: a comment on the whole file and on the whole review, Export CSV.

Keys (while the tab is visible and no text field has the focus; the tab's
lists and checkboxes do not hold them): Left/Right move, N jumps to the next
unrated frame, P passes, B blocks, 1-9 toggle the categories in their listed
order, Esc leaves a text field of the tab.  A rating or category key carries
the position that was on screen and is dropped when the kernel has moved on;
auto-repeat is dropped except for the arrows.  The
browser hands each key to a hidden text field, the same bridge the Ketcher
panel uses, so every key also works headless through
:meth:`ReviewPanel.handle_key`.  The script that listens for the keys is
registered with the dashboard's start-up scripts (``ctx.add_init_js``):
written into the shared script output on its own, it was cleared again by
the dashboard's one start-up send before the browser ran it, and no key
reached the tab.
"""

from __future__ import annotations

import html
from datetime import datetime
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
AT_CLASS = 'delfin-review-at-'
PASS_CLASS = 'delfin-review-pass'
BLOCK_CLASS = 'delfin-review-block'
FRAME_WINDOW = 100
VIEWER_WIDTH = '100%'   # of the middle column; a fixed 720 px ran under the rating column
VIEWER_HEIGHT = 520
_MARK = {'pass': '✓', 'block': '✗', None: '·'}
_HELP = "<span style='color:#757575;font-size:11.5px'>{}</span>"
_PASS_COLOR = '#1b7f3b'
_BLOCK_COLOR = '#b00020'
KEY_LEGEND = ('Keys: <b>P</b> pass &middot; <b>B</b> block &middot; <b>1</b>&ndash;<b>9</b> '
              'categories &middot; <b>N</b> next unrated &middot; <b>&larr;</b> / <b>&rarr;</b> '
              'previous / next &middot; <b>Esc</b> leaves a text field')
_CSS = """<style>
.delfin-review-panel .delfin-review-big {font-size:16px;font-weight:600;}
.delfin-review-panel .delfin-review-cat {text-align:left;justify-content:flex-start;}
.delfin-review-panel .widget-hbox {gap:8px;}
.delfin-review-panel .widget-vbox {gap:4px;}
.delfin-review-panel .delfin-review-bar {border-bottom:1px solid #e0e0e0;padding:2px 0 6px 0;}
</style>"""


def _keys_js() -> str:
    """Hand arrow keys, P, B, N and digits to the hidden field while the panel is shown.

    The line sent is ``<ms> <key> <position>``: the position is read from the
    ``delfin-review-at-<n>`` class the kernel sets on the hidden field after it
    has drawn frame *n*, so a rating key names the frame that was on screen.
    Auto-repeat is dropped for every key but the arrows (a held P would rate
    each frame it advanced to).  A text field, a text input and a select
    outside the panel keep the keys; the panel's own lists and checkboxes do
    not, so P after a click into the frame list still rates.

    A click on Pass or Block is taken in the capture phase and sent the same
    way (as ``p`` / ``b`` with the position), and does not reach the
    button's own handler: a button click carries no data to the kernel, so a
    second click before the next frame was drawn would rate that frame.
    Without the script (a plain notebook) the buttons' ``on_click`` rates.
    """
    return (
        "(function(){\n"
        "  if(window.__delfinReviewKeys) return;\n"
        "  window.__delfinReviewKeys=true;\n"
        "  var WANT={ArrowLeft:'left',ArrowRight:'right',p:'p',P:'p',b:'b',B:'b',n:'n',N:'n'};\n"
        "  function send(panel,key){\n"
        "    var box=panel.querySelector('." + KEYS_CLASS + "');\n"
        "    var input=box&&box.querySelector('input, textarea');\n"
        "    if(!input) return false;\n"
        "    var at=/(?:^|\\s)" + AT_CLASS + "(\\d+)(?:\\s|$)/.exec(box.className);\n"
        "    var setter=Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype,'value');\n"
        "    var line=Date.now()+' '+key+(at?' '+at[1]:'');\n"
        "    if(setter&&setter.set) setter.set.call(input,line); else input.value=line;\n"
        "    input.dispatchEvent(new Event('input',{bubbles:true}));\n"
        "    input.dispatchEvent(new Event('change',{bubbles:true}));\n"
        "    return true;\n"
        "  }\n"
        "  document.addEventListener('keydown',function(ev){\n"
        "    if(ev.ctrlKey||ev.metaKey||ev.altKey) return;\n"
        "    var panel=document.querySelector('." + PANEL_CLASS + "');\n"
        "    if(!panel||panel.offsetParent===null) return;\n"
        "    var el=ev.target, tag=(el&&el.tagName||'').toUpperCase();\n"
        "    var kind=(el&&el.type||'').toLowerCase();\n"
        "    var typing=tag==='TEXTAREA'||(el&&el.isContentEditable)||\n"
        "      (tag==='INPUT'&&['checkbox','radio','button','submit','reset','range'].indexOf(kind)<0)||\n"
        "      (tag==='SELECT'&&!panel.contains(el));\n"
        "    if(typing){\n"
        "      if(ev.key==='Escape'&&panel.contains(el)){el.blur();ev.preventDefault();}\n"
        "      return;\n"
        "    }\n"
        "    var key=WANT[ev.key]||(/^[1-9]$/.test(ev.key)?ev.key:null);\n"
        "    if(!key) return;\n"
        "    if(ev.repeat&&key!=='left'&&key!=='right'){ev.preventDefault();return;}\n"
        "    if(!send(panel,key)) return;\n"
        "    ev.preventDefault();\n"
        "  },true);\n"
        "  document.addEventListener('click',function(ev){\n"
        "    var panel=document.querySelector('." + PANEL_CLASS + "');\n"
        "    var btn=ev.target&&ev.target.closest&&ev.target.closest('."
        + PASS_CLASS + ", ." + BLOCK_CLASS + "');\n"
        "    if(!panel||!btn||!panel.contains(btn)) return;\n"
        "    var el=btn.tagName==='BUTTON'?btn:btn.querySelector('button');\n"
        "    if((el&&el.disabled)||btn.disabled) return;\n"
        "    var key=btn.classList.contains('" + PASS_CLASS + "')?'p':'b';\n"
        "    if(!send(panel,key)) return;\n"
        "    ev.preventDefault();ev.stopPropagation();\n"
        "  },true);\n"
        "})();\n"
    )


def _paint_category(change):
    """A chosen category stands out in colour, not only by a slightly darker grey."""
    change['owner'].button_style = 'warning' if change.get('new') else ''


class ReviewPanel:
    """Widgets and behaviour of the Review tab."""

    def __init__(self, start_dir: Any = None, ctx: Any = None):
        self.ctx = ctx
        self.session: rv.ReviewSession | None = None
        self._busy = False
        start = str(start_dir or Path.home())
        W, L = widgets, widgets.Layout
        label = {'description_width': '70px'}

        # -- open a review ---------------------------------------------------
        self.path_input = W.Text(value=start, description='Open:', style=label,
                                 placeholder='multi-frame .xyz or a folder of .xyz files',
                                 layout=L(width='560px', max_width='100%'))
        self.open_btn = W.Button(description='Open', button_style='info', icon='folder-open',
                                 layout=L(width='100px'))
        self.name_input = W.Text(value='', description='Name:', style=label,
                                 placeholder='default: folder or file name',
                                 layout=L(width='360px', max_width='100%'))
        self.findings_input = W.Text(value='', description='Findings:', style=label,
                                     placeholder='optional JSON; default findings.json',
                                     layout=L(width='560px', max_width='100%'))
        self.blind_cb = W.Checkbox(value=False, description='Blinded', indent=False,
                                   layout=L(width='90px'))
        self.seed_input = W.IntText(value=0, description='Seed:',
                                    style={'description_width': '40px'},
                                    layout=L(width='140px'))
        self.categories_input = W.Textarea(value='\n'.join(rv.DEFAULT_CATEGORIES),
                                           layout=L(width='360px', max_width='100%',
                                                    height='170px'))
        self.categories_btn = W.Button(description='Apply categories',
                                       layout=L(width='150px'))

        # -- the bar: list, navigation, position, progress -------------------
        self.list_toggle = W.ToggleButton(value=True, description='List', icon='list',
                                          tooltip='show or hide the files and frames',
                                          layout=L(width='80px'))
        self.prev_btn = W.Button(description='Prev', icon='arrow-left',
                                 tooltip='previous frame (←)', layout=L(width='90px'))
        self.next_btn = W.Button(description='Next', icon='arrow-right',
                                 tooltip='next frame (→)', layout=L(width='90px'))
        self.unrated_btn = W.Button(description='Next unrated', tooltip='next unrated (N)',
                                    layout=L(width='120px'))
        self.position_html = W.HTML('')
        self.progress_html = W.HTML('')

        # -- left: files and frames ------------------------------------------
        self.left_title = W.HTML('<b>Files</b>')
        self.file_list = W.Select(options=[], rows=10, layout=L(width='100%'))
        self.frame_list = W.Select(options=[], rows=16, layout=L(width='100%'))

        # -- middle: the structure ---------------------------------------------
        self.viewer = W.Output(layout=L(width='100%', min_height=f'{VIEWER_HEIGHT}px',
                                        border='1px solid #ddd'))

        # -- right: the rating ---------------------------------------------------
        self.pass_btn = W.Button(description='Pass (P)', button_style='success', icon='check',
                                 layout=L(width='48%', height='52px'))
        self.block_btn = W.Button(description='Block (B)', button_style='danger', icon='times',
                                  layout=L(width='48%', height='52px'))
        for b in (self.pass_btn, self.block_btn):
            b.add_class('delfin-review-big')
        self.pass_btn.add_class(PASS_CLASS)
        self.block_btn.add_class(BLOCK_CLASS)
        self.verdict_html = W.HTML('')
        self.multi_cb = W.Checkbox(value=True, description='Several categories',
                                   indent=False, layout=L(width='170px'))
        self.category_box = W.VBox([], layout=L(width='100%'))
        self.category_buttons: list[Any] = []
        self.note_input = W.Textarea(value='', placeholder='your comment on this frame',
                                     continuous_update=False,
                                     layout=L(width='100%', height='60px'))
        self.saved_html = W.HTML('')
        self.help_html = W.HTML(_HELP.format(KEY_LEGEND))
        self.context_html = W.HTML('')

        # -- more: comments on the file and the review, export ------------------
        self.file_comment = W.Textarea(value='', placeholder='comment on this file',
                                       continuous_update=False,
                                       layout=L(width='100%', height='60px'))
        self.session_comment = W.Textarea(value='', placeholder='comment on the whole review',
                                          continuous_update=False,
                                          layout=L(width='100%', height='60px'))
        self.export_btn = W.Button(description='Export CSV', icon='download',
                                   layout=L(width='130px'))
        self.status_html = W.HTML('')
        self.keys_input = W.Text(value='', layout=L(display='none'))
        self.keys_input.add_class(KEYS_CLASS)
        self._at_class = ''

        self._build_categories(rv.DEFAULT_CATEGORIES)
        self._wire()
        self.widget = self._layout()
        self.widget.add_class(PANEL_CLASS)
        self._set_enabled(False)
        self._install_keys()

    # -- layout ------------------------------------------------------------
    def _layout(self):
        W, L = widgets, widgets.Layout
        row = dict(flex_flow='row wrap', align_items='center')
        setup_body = W.VBox([
            W.HTML(_HELP.format(
                'Open a multi-frame .xyz or a folder of .xyz files. Every rating is saved at '
                'once to review_&lt;name&gt;.json in that folder; opening the same folder and '
                'name again resumes. Set the categories before you open a new review (one per '
                'line; keys 1&ndash;9 in this order); a resumed review keeps its own.')),
            W.HBox([self.path_input, self.open_btn], layout=L(**row)),
            W.HBox([self.name_input], layout=L(**row)),
            W.HBox([self.findings_input], layout=L(**row)),
            W.HBox([self.blind_cb, self.seed_input,
                    W.HTML(_HELP.format('Blinded: all frames shuffled by the seed; file '
                                        'names, labels and findings hidden'))],
                   layout=L(**row)),
            W.HBox([W.HTML('<b>Categories</b>'), self.categories_input, self.categories_btn],
                   layout=L(flex_flow='row wrap', align_items='flex-start')),
        ], layout=L())
        self.setup = W.Accordion(children=[setup_body])
        self.setup.set_title(0, 'Open a review')
        self.setup.selected_index = 0

        self.bar = W.HBox([
            self.list_toggle, self.prev_btn, self.next_btn, self.unrated_btn,
            W.Box([self.position_html], layout=L(flex='1 1 260px', padding='0 8px')),
            W.Box([self.progress_html], layout=L(width='300px')),
        ], layout=L(width='100%', **row))
        self.bar.add_class('delfin-review-bar')

        self.files_box = W.VBox([self.left_title, self.file_list])
        self.sidebar = W.VBox([self.files_box, W.HTML('<b>Frames</b>'), self.frame_list],
                              layout=L(width='200px', min_width='160px'))
        middle = W.VBox([self.viewer], layout=L(flex='1 1 480px', min_width='240px'))
        self.rating_box = W.VBox([
            W.HBox([self.pass_btn, self.block_btn],
                   layout=L(width='100%', justify_content='space-between')),
            self.verdict_html,
            W.HBox([W.HTML('<b>Categories</b>'), self.multi_cb],
                   layout=L(align_items='center')),
            self.category_box,
            W.HTML('<b>Comment</b>'),
            self.note_input,
            self.saved_html,
            self.help_html,
            self.context_html,
        ], layout=L(width='300px', min_width='260px'))

        self.file_comment_box = W.VBox([W.HTML('<b>Comment on this file</b>'),
                                        self.file_comment])
        more_body = W.VBox([
            self.file_comment_box,
            W.HTML('<b>Comment on this review</b>'), self.session_comment,
            self.export_btn,
        ], layout=L(width='100%'))
        self.more = W.Accordion(children=[more_body])
        self.more.set_title(0, 'More: comments on the file and the review, export')
        self.more.selected_index = None

        self.work = W.VBox([
            self.bar,
            W.HBox([self.sidebar, middle, self.rating_box],
                   layout=L(width='100%', flex_flow='row wrap',
                            align_items='flex-start')),
            self.more,
        ], layout=L(width='100%', display='none'))
        return W.VBox([
            W.HTML(_CSS),
            W.HTML('<h3 style="margin:0 0 4px 0">Review structures</h3>'),
            self.setup,
            self.work,
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
        self.list_toggle.observe(self._on_list_toggle, names='value')
        self.blind_cb.observe(self._on_mode, names='value')
        self.seed_input.observe(self._on_mode, names='value')
        self.file_list.observe(self._on_file_pick, names='value')
        self.frame_list.observe(self._on_frame_pick, names='value')
        self.note_input.observe(self._on_note, names='value')
        self.file_comment.observe(self._on_file_comment, names='value')
        self.session_comment.observe(self._on_session_comment, names='value')
        self.keys_input.observe(self._on_key, names='value')

    def _install_keys(self):
        """Register the key script with the dashboard's start-up scripts.

        The dashboard sends every tab's start-up script in one go and clears
        the shared script output before it writes; a script appended there
        on its own while the tab is built is gone before the browser runs it.
        """
        if self.ctx is None:
            return
        add = getattr(self.ctx, 'add_init_js', None)
        if callable(add):
            add(_keys_js())
            return
        if getattr(self.ctx, 'js_output', None) is None:
            return
        try:
            from delfin.dashboard.helpers import _append_js
            _append_js(self.ctx, _keys_js())
        except Exception:                            # pragma: no cover
            pass

    def _build_categories(self, categories):
        self.category_buttons = []
        for i, cat in enumerate(categories):
            hint = f'{i + 1}  ' if i < 9 else ''
            btn = widgets.ToggleButton(value=False, description=f'{hint}{cat}', tooltip=cat,
                                       disabled=self.session is None,
                                       layout=widgets.Layout(width='100%', height='30px'))
            btn.add_class('delfin-review-cat')
            btn._review_category = cat
            btn.observe(self._on_category, names='value')
            btn.observe(_paint_category, names='value')
            self.category_buttons.append(btn)
        self.category_box.children = tuple(self.category_buttons)

    def _set_enabled(self, on: bool):
        for w in (self.prev_btn, self.next_btn, self.unrated_btn, self.pass_btn,
                  self.block_btn, self.export_btn, self.note_input,
                  self.file_comment, self.session_comment, *self.category_buttons):
            w.disabled = not on

    def _say(self, text: str, error: bool = False):
        color = _BLOCK_COLOR if error else '#555'
        self.status_html.value = (f'<span style="color:{color};font-size:12px">'
                                  f'{html.escape(text)}</span>')

    # -- opening -------------------------------------------------------------
    def open(self, path: str | None = None) -> rv.ReviewSession | None:
        if path is not None:
            self.path_input.value = str(path)
        raw = self.path_input.value.strip()
        if not raw:
            self._say('Give a multi-frame .xyz file or a folder.', error=True)
            return None
        findings = self.findings_input.value.strip() or None
        wanted = [c.strip() for c in self.categories_input.value.splitlines() if c.strip()]
        try:
            session = rv.ReviewSession.open(
                raw, name=self.name_input.value, blind=self.blind_cb.value,
                seed=int(self.seed_input.value or 0), findings_path=findings,
                default_categories=wanted or None)
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
        self.work.layout.display = ''
        self.setup.selected_index = None
        self.saved_html.value = ''
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
        where = s.display_name()
        rec = s.rate(verdict, self.selected_categories(), self.note_input.value)
        color = _PASS_COLOR if rec['verdict'] == 'pass' else _BLOCK_COLOR
        self.saved_html.value = (
            f'<span style="color:{color};font-size:12px">&#10003; saved '
            f'{datetime.now().strftime("%H:%M:%S")}: {html.escape(where)} &rarr; '
            f'<b>{rec["verdict"]}</b></span>')
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
        if self.session is None:
            self._say(f'{len(cats)} categories set; a new review opened next uses them.')
            return
        self.session.set_categories(cats)
        self.refresh()
        self._say(f'{len(cats)} categories applied and saved in {self.session.review_path.name}.')

    def handle_key(self, key: str, at: int | None = None):
        """Act on one key; *at* is the position the browser showed when it was pressed.

        A rating or category key whose *at* is not the current position is
        dropped: the kernel moved on (after a rating) before the screen did.
        Navigation keys are relative and always act.
        """
        key = (key or '').strip().lower()
        if self.session is None:
            return
        if (at is not None and at != self.session.position
                and (key in ('p', 'b') or key.isdigit())):
            self._say(f'Key {key.upper()} ignored: it was pressed for a frame no longer shown.')
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
        parts = str(change.get('new') or '').split()
        if len(parts) < 2:
            return
        at = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else None
        self.handle_key(parts[1], at)

    def _on_list_toggle(self, change):
        self.sidebar.layout.display = '' if change.get('new') else 'none'

    def _on_mode(self, _change):
        if self.session is None or self._busy:
            return
        self.session.set_mode(self.blind_cb.value, int(self.seed_input.value or 0))
        self.session.save()
        self.saved_html.value = ''
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
        """A changed comment is stored at once: in the rating, or for an unrated frame in ``notes``."""
        s = self.session
        if self._busy or s is None:
            return
        frame = s.current
        rec = s.record(frame.file, frame.index)
        text = change['new'].strip()
        if rec is None:
            if text == s.pending_note(frame.file, frame.index):
                return
            s.set_pending_note(text)
        elif rec.get('note', '') != text:
            rec['note'] = text
            rec['note_updated'] = rv._now()
            s.save()
        else:
            return
        self.saved_html.value = (
            f'<span style="color:#555;font-size:12px">&#10003; comment saved '
            f'{datetime.now().strftime("%H:%M:%S")}</span>')

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
            self.note_input.value = (rec.get('note', '') if rec
                                     else s.pending_note(frame.file, frame.index))
            meta = s.data['files'].get(frame.file, {})
            self.file_comment.value = '' if s.blind else meta.get('comment', '')
        finally:
            self._busy = False
        self.refresh_lists(full=full)
        self._draw_position(frame, rec)
        self._draw_progress()
        self._draw_context(frame)
        self._draw_viewer(frame)
        self._mark_shown(s.position)

    def _mark_shown(self, position: int):
        """Tell the key script which frame is on screen (sent after the frame itself)."""
        if self._at_class:
            self.keys_input.remove_class(self._at_class)
        self._at_class = f'{AT_CLASS}{position}'
        self.keys_input.add_class(self._at_class)

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
                self.file_list.options = []
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
            where = f'<b>structure #{s.position + 1}</b> of {len(s)}'
        else:
            n = s.data['files'].get(frame.file, {}).get('n_frames', 0)
            where = (f'<b>{html.escape(frame.file)}</b> &nbsp; frame {frame.index + 1}/{n}'
                     f' &nbsp; <i>{html.escape(frame.label)}</i>')
        if rec:
            color = _PASS_COLOR if rec['verdict'] == 'pass' else _BLOCK_COLOR
            cats = ', '.join(rec.get('categories', []))
            state = (f'<span style="color:{color}"><b>{rec["verdict"]}</b></span>'
                     + (f' ({html.escape(cats)})' if cats else ''))
        else:
            state = '<span style="color:#888">not rated</span>'
        self.position_html.value = f'<span style="font-size:14px">{where} &nbsp;&mdash;&nbsp; {state}</span>'
        self.verdict_html.value = f'<span style="font-size:12px">This frame: {state}</span>'

    def _draw_progress(self):
        p = self.session.progress()
        pct = 100.0 * p['rated'] / p['total'] if p['total'] else 0.0
        self.progress_html.value = (
            f'<div style="font-size:13px"><b>{p["rated"]} / {p["total"]}</b> rated '
            f'({pct:.0f} %) &nbsp; <span style="color:{_PASS_COLOR}">pass {p["pass"]}</span>'
            f' &nbsp; <span style="color:{_BLOCK_COLOR}">block {p["block"]}</span></div>'
            f'<div style="background:#eee;height:6px;margin-top:3px;border-radius:3px">'
            f'<div style="background:#4a90d9;height:6px;width:{pct:.1f}%;'
            f'border-radius:3px"></div></div>'
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
            parts.append(f'<b>SMILES</b><br><code style="word-break:break-all;font-size:12px">'
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
            render_xyz_in_output(self.viewer, frame.xyz(comment),
                                 width=VIEWER_WIDTH, height=VIEWER_HEIGHT)
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
