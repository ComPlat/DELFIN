"""The Review tab: look at every frame of a multi-frame XYZ and rate it.

The ratings themselves -- loading, order, saving, resuming -- live in
:mod:`delfin.review`; this is only the screen on top of it.  The 3D view is
the dashboard's own (:func:`molecule_viewer.render_xyz_in_output`), so the
viewer settings apply here exactly as they do in the Calculations tab.

What is on the screen
---------------------

* **Open a review** (a compact section that closes once a review is open):
  what to open (a multi-frame XYZ or a folder of them), the review name, an
  optional findings file, the blinded switch with its seed and the list of
  categories.  A ``kategorien.txt`` beside the structures (or one folder up)
  fills the list for a new review while it still holds the built-in one.
* A bar across the top: show/hide the list, Prev / Next / Next unrated,
  which frame is shown, show/hide the H atoms and reset the view.
* Left, narrow and foldable: the files (normal mode) and the frames, rated
  ones marked; in the blinded mode only anonymous numbers (folded away by
  default there).
* Middle, as large as the window allows: the structure.  Its height is a
  share of the window height (``calc(100vh - ...)``, as the other tab
  bodies), and the dashboard's viewer factory keeps the 3D canvas sized to
  its box (``__delfinKeepSized``), so the picture follows the window.
* Right, compact: the progress, a large Pass and a large Block button, the
  categories as a grid of toggles that show their key, the comment on this
  frame with an Unsure toggle, a line that says when the rating was saved, a
  short key legend, and **More** (folded): SMILES and findings (normal mode),
  a comment on the whole file and on the whole review, adding categories,
  Export CSV.

Keys (while the tab is visible and no text field has the focus; the tab's
lists and checkboxes do not hold them): Left/Right move, N jumps to the next
unrated frame, P passes, B blocks, 1-9 and 0 toggle the first ten
categories, Shift+1..Shift+9 the categories 11-19 (read from the key's
place, ``event.code``, so the keyboard layout does not matter), U adds or
removes the word "unsicher" in the comment, H shows/hides the H atoms, R
resets the view, Esc leaves a text field of the tab.  A rating, category or
U key carries the position that was on screen and is dropped when the
kernel has moved on; auto-repeat is dropped except for the arrows.  The
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
import re
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

__all__ = ['ReviewPanel', 'create_tab', 'category_key', 'category_index']

PANEL_CLASS = 'delfin-review-panel'
KEYS_CLASS = 'delfin-review-keys'
AT_CLASS = 'delfin-review-at-'
PASS_CLASS = 'delfin-review-pass'
BLOCK_CLASS = 'delfin-review-block'
SIDE_CLASS = 'delfin-review-side'
VIEWER_CLASS = 'delfin-review-viewer'
FRAME_WINDOW = 100
VIEWER_WIDTH = '100%'   # of the middle column; a fixed 720 px ran under the rating column
# The dashboard's tab bodies are calc(100vh - 145px) tall (helpers.create_page_css);
# the folded setup, the status line and the bar above the viewer take about 130 px
# more.  A CSS length, not pixels: the viewer box follows the window, and the
# viewer factory's ResizeObserver (__delfinKeepSized) resizes the canvas with it.
VIEWER_HEIGHT = 'max(360px, calc(100vh - 275px))'
MAX_KEYED_CATEGORIES = 19
_MARK = {'pass': '✓', 'block': '✗', None: '·'}
_HELP = "<span style='color:#757575;font-size:11.5px'>{}</span>"
_PASS_COLOR = '#1b7f3b'
_BLOCK_COLOR = '#b00020'
KEY_LEGEND = ('<b>P</b> pass &middot; <b>B</b> block &middot; <b>1</b>&ndash;<b>9</b>, '
              '<b>0</b>, <b>&#8679;1</b>&ndash;<b>&#8679;9</b> categories &middot; '
              '<b>U</b> unsure &middot; <b>N</b> next unrated &middot; <b>&larr;</b> / '
              '<b>&rarr;</b> previous / next &middot; <b>H</b> H atoms &middot; '
              '<b>R</b> reset view &middot; <b>Esc</b> leaves a text field')
_CSS = """<style>
.delfin-review-panel .delfin-review-big {font-size:16px;font-weight:600;}
.delfin-review-panel .delfin-review-cat {text-align:left;justify-content:flex-start;
  overflow:hidden;text-overflow:ellipsis;white-space:nowrap;padding:0 6px;font-size:12.5px;}
.delfin-review-panel .widget-hbox {gap:6px;}
.delfin-review-panel .widget-vbox {gap:3px;}
.delfin-review-panel .delfin-review-bar {border-bottom:1px solid #e0e0e0;padding:0 0 4px 0;}
.delfin-review-panel .delfin-review-side .widget-select {flex:1 1 0 !important;min-height:0 !important;}
.delfin-review-panel .delfin-review-side .widget-select select {height:100% !important;}
.delfin-review-panel .delfin-review-viewer {overflow:hidden !important;}
.delfin-review-panel .delfin-review-rating {overflow-y:auto !important;overflow-x:hidden !important;}
</style>"""


# ---------------------------------------------------------------------------
# Keys of the categories
# ---------------------------------------------------------------------------

def category_key(i: int) -> str:
    """The key sent for the category at 0-based place *i*: 1-9, 0, s1-s9 (Shift), else ''."""
    if 0 <= i < 9:
        return str(i + 1)
    if i == 9:
        return '0'
    if 10 <= i < MAX_KEYED_CATEGORIES:
        return f's{i - 9}'
    return ''


def category_key_label(i: int) -> str:
    """The key as shown on the button (Shift as an arrow)."""
    key = category_key(i)
    return f'⇧{key[1:]}' if key.startswith('s') else key


def category_index(key: str) -> int | None:
    """The 0-based place of the category a key names, or None."""
    key = (key or '').strip().lower()
    if re.fullmatch(r'[1-9]', key):
        return int(key) - 1
    if key == '0':
        return 9
    if re.fullmatch(r's[1-9]', key):
        return 9 + int(key[1])
    return None


def _keys_js() -> str:
    """Hand arrow keys, P, B, N, U, H, R and digits to the hidden field while the panel is shown.

    The line sent is ``<ms> <key> <position>``: the position is read from the
    ``delfin-review-at-<n>`` class the kernel sets on the hidden field after it
    has drawn frame *n*, so a rating key names the frame that was on screen.
    Digits are read from the key's place (``event.code`` ``Digit0``-``Digit9``,
    ``Numpad0``-``Numpad9``): with Shift they are sent as ``s1``-``s9`` (what
    Shift+1 types differs from layout to layout).  Auto-repeat is dropped for
    every key but the arrows (a held P would rate each frame it advanced to).
    A text field, a text input and a select outside the panel keep the keys;
    the panel's own lists and checkboxes do not, so P after a click into the
    frame list still rates.

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
        "  var WANT={ArrowLeft:'left',ArrowRight:'right',p:'p',P:'p',b:'b',B:'b',n:'n',N:'n',\n"
        "    u:'u',U:'u',h:'h',H:'h',r:'r',R:'r'};\n"
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
        "    var key=null, dm=/^(Digit|Numpad)([0-9])$/.exec(ev.code||'');\n"
        "    if(dm){\n"
        "      if(!ev.shiftKey) key=dm[2];\n"
        "      else if(dm[1]==='Digit'&&dm[2]!=='0') key='s'+dm[2];\n"
        "    } else {\n"
        "      key=WANT[ev.key]||(/^[0-9]$/.test(ev.key)?ev.key:null);\n"
        "    }\n"
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


def _without_h(xyz: str) -> str:
    """An XYZ block without its hydrogen atoms (for the view only)."""
    lines = xyz.splitlines()
    if len(lines) < 2:
        return xyz
    atoms = [ln for ln in lines[2:] if ln.split() and ln.split()[0].upper() not in ('H', 'D')]
    return f'{len(atoms)}\n{lines[1]}\n' + '\n'.join(atoms) + '\n'


def _split_categories(text: str) -> list[str]:
    return rv.clean_categories(text.splitlines())


class ReviewPanel:
    """Widgets and behaviour of the Review tab."""

    def __init__(self, start_dir: Any = None, ctx: Any = None):
        self.ctx = ctx
        self.session: rv.ReviewSession | None = None
        self._busy = False
        self._prefilled: list[str] | None = None
        start = str(start_dir or Path.home())
        W, L = widgets, widgets.Layout
        label = {'description_width': '62px'}

        # -- open a review ---------------------------------------------------
        self.path_input = W.Text(value=start, description='Open:', style=label,
                                 placeholder='multi-frame .xyz or a folder of .xyz files',
                                 layout=L(width='560px', max_width='100%'))
        self.open_btn = W.Button(description='Open', button_style='info', icon='folder-open',
                                 layout=L(width='90px'))
        self.name_input = W.Text(value='', description='Name:', style=label,
                                 placeholder='default: folder or file name',
                                 layout=L(width='300px', max_width='100%'))
        self.findings_input = W.Text(value='', description='Findings:', style=label,
                                     placeholder='optional JSON; default findings.json',
                                     layout=L(width='360px', max_width='100%'))
        self.blind_cb = W.Checkbox(value=False, description='Blinded', indent=False,
                                   layout=L(width='80px'))
        self.seed_input = W.IntText(value=0, description='Seed:',
                                    style={'description_width': '36px'},
                                    layout=L(width='120px'))
        self.categories_input = W.Textarea(value='\n'.join(rv.DEFAULT_CATEGORIES),
                                           layout=L(width='300px', max_width='100%',
                                                    height='130px'))
        self.categories_btn = W.Button(description='Apply categories',
                                       layout=L(width='140px'))

        # -- the bar: list, navigation, position, view -----------------------
        self.list_toggle = W.ToggleButton(value=True, description='List', icon='list',
                                          tooltip='show or hide the files and frames',
                                          layout=L(width='72px'))
        self.prev_btn = W.Button(description='Prev', icon='arrow-left',
                                 tooltip='previous frame (←)', layout=L(width='80px'))
        self.next_btn = W.Button(description='Next', icon='arrow-right',
                                 tooltip='next frame (→)', layout=L(width='80px'))
        self.unrated_btn = W.Button(description='Next unrated', tooltip='next unrated (N)',
                                    layout=L(width='112px'))
        self.position_html = W.HTML('')
        self.h_toggle = W.ToggleButton(value=True, description='H atoms', icon='eye',
                                       tooltip='show or hide the hydrogen atoms (H)',
                                       layout=L(width='96px'))
        self.reset_btn = W.Button(description='Reset view', icon='refresh',
                                  tooltip='orientation and zoom as at the start (R)',
                                  layout=L(width='108px'))

        # -- left: files and frames ------------------------------------------
        self.left_title = W.HTML('<b>Files</b>')
        self.file_list = W.Select(options=[], rows=8, layout=L(width='100%'))
        self.frame_list = W.Select(options=[], rows=16, layout=L(width='100%'))

        # -- middle: the structure ---------------------------------------------
        self.viewer = W.Output(layout=L(width='100%', min_height='360px',
                                        border='1px solid #ddd'))
        self.viewer.add_class(VIEWER_CLASS)

        # -- right: the rating ---------------------------------------------------
        self.progress_html = W.HTML('')
        self.pass_btn = W.Button(description='Pass (P)', button_style='success', icon='check',
                                 layout=L(width='50%', height='54px'))
        self.block_btn = W.Button(description='Block (B)', button_style='danger', icon='times',
                                  layout=L(width='50%', height='54px'))
        for b in (self.pass_btn, self.block_btn):
            b.add_class('delfin-review-big')
        self.pass_btn.add_class(PASS_CLASS)
        self.block_btn.add_class(BLOCK_CLASS)
        self.verdict_html = W.HTML('')
        self.multi_cb = W.Checkbox(value=True, description='several', indent=False,
                                   tooltip='several categories per frame',
                                   layout=L(width='90px'))
        self.category_box = W.GridBox([], layout=L(
            width='100%', grid_template_columns='repeat(2, minmax(0, 1fr))',
            grid_gap='3px 4px'))
        self.category_buttons: list[Any] = []
        self.unsure_btn = W.ToggleButton(value=False, description='Unsure (U)',
                                         icon='question',
                                         tooltip='adds or removes the word "unsicher" '
                                                 'in the comment (U)',
                                         layout=L(width='110px', height='26px'))
        self.note_input = W.Textarea(value='', placeholder='your comment on this frame',
                                     continuous_update=False,
                                     layout=L(width='100%', height='58px'))
        self.saved_html = W.HTML('')
        self.help_html = W.HTML(_HELP.format(KEY_LEGEND))
        self.context_html = W.HTML('')

        # -- more: findings, comments on the file and the review, categories, export
        self.file_comment = W.Textarea(value='', placeholder='comment on this file',
                                       continuous_update=False,
                                       layout=L(width='100%', height='56px'))
        self.session_comment = W.Textarea(value='', placeholder='comment on the whole review',
                                          continuous_update=False,
                                          layout=L(width='100%', height='56px'))
        self.add_cat_input = W.Text(value='', placeholder='new categories, comma separated',
                                    layout=L(flex='1 1 auto', width='auto', min_width='0'))
        self.add_cat_btn = W.Button(description='Add', icon='plus', layout=L(width='70px'),
                                    tooltip='append categories; the others keep their keys')
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
        self._prefill_categories()

    # -- layout ------------------------------------------------------------
    def _layout(self):
        W, L = widgets, widgets.Layout
        row = dict(flex_flow='row wrap', align_items='center')
        setup_body = W.VBox([
            W.HBox([self.path_input, self.open_btn], layout=L(**row)),
            W.HBox([self.name_input, self.findings_input, self.blind_cb, self.seed_input],
                   layout=L(**row)),
            W.HBox([
                W.VBox([W.HTML('<b>Categories</b> ' + _HELP.format(
                            'one per line; keys 1&ndash;9, 0, &#8679;1&ndash;&#8679;9 in this '
                            'order; a kategorien.txt beside the structures fills a new '
                            'review')),
                        self.categories_input]),
                W.VBox([self.categories_btn, W.HTML(_HELP.format(
                    'Every rating is saved at once to review_&lt;name&gt;.json beside the '
                    'structures; the same folder and name again resumes, with the '
                    'review&rsquo;s own categories. Blinded: all frames shuffled by the '
                    'seed; file names, labels and findings hidden.'))],
                    layout=L(width='300px', max_width='100%')),
            ], layout=L(flex_flow='row wrap', align_items='flex-start')),
        ], layout=L())
        self.setup = W.Accordion(children=[setup_body])
        self.setup.set_title(0, 'Open a review')
        self.setup.selected_index = 0

        self.bar = W.HBox([
            self.list_toggle, self.prev_btn, self.next_btn, self.unrated_btn,
            W.Box([self.position_html], layout=L(flex='1 1 240px', padding='0 6px',
                                                 min_width='0', overflow='hidden')),
            self.h_toggle, self.reset_btn,
        ], layout=L(width='100%', **row))
        self.bar.add_class('delfin-review-bar')

        self.files_box = W.VBox([self.left_title, self.file_list])
        self.sidebar = W.VBox([self.files_box, W.HTML('<b>Frames</b>'), self.frame_list],
                              layout=L(width='170px', min_width='140px',
                                       height=VIEWER_HEIGHT))
        self.sidebar.add_class(SIDE_CLASS)
        middle = W.VBox([self.viewer], layout=L(flex='1 1 420px', min_width='240px'))

        self.file_comment_box = W.VBox([W.HTML('<b>Comment on this file</b>'),
                                        self.file_comment])
        more_body = W.VBox([
            self.context_html,
            self.file_comment_box,
            W.HTML('<b>Comment on this review</b>'), self.session_comment,
            W.HTML('<b>Add categories</b> ' + _HELP.format(
                'appended at the end; ratings and keys stay')),
            W.HBox([self.add_cat_input, self.add_cat_btn], layout=L(width='100%')),
            self.export_btn,
        ], layout=L(width='100%'))
        self.more = W.Accordion(children=[more_body])
        self.more.set_title(0, 'More: findings, comments, categories, export')
        self.more.selected_index = None

        self.rating_box = W.VBox([
            self.progress_html,
            W.HBox([self.pass_btn, self.block_btn],
                   layout=L(width='100%', justify_content='space-between')),
            self.verdict_html,
            W.HBox([W.HTML('<b>Categories</b>'), self.multi_cb],
                   layout=L(align_items='center', justify_content='space-between')),
            self.category_box,
            W.HBox([W.HTML('<b>Comment</b>'), self.unsure_btn],
                   layout=L(align_items='center', justify_content='space-between')),
            self.note_input,
            self.saved_html,
            self.help_html,
            self.more,
        ], layout=L(width='340px', min_width='260px', max_height=VIEWER_HEIGHT))
        self.rating_box.add_class('delfin-review-rating')

        self.work = W.VBox([
            self.bar,
            W.HBox([self.sidebar, middle, self.rating_box],
                   layout=L(width='100%', flex_flow='row wrap',
                            align_items='flex-start')),
        ], layout=L(width='100%', display='none'))
        return W.VBox([
            W.HTML(_CSS),
            self.setup,
            self.status_html,
            self.work,
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
        self.add_cat_btn.on_click(lambda _b: self.add_categories(self.add_cat_input.value))
        self.reset_btn.on_click(lambda _b: self.reset_view())
        self.h_toggle.observe(self._on_h_toggle, names='value')
        self.unsure_btn.observe(self._on_unsure_btn, names='value')
        self.list_toggle.observe(self._on_list_toggle, names='value')
        self.path_input.observe(lambda _c: self._prefill_categories(), names='value')
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
        """One toggle per category, its key in front; a choice on the frame is kept by name."""
        chosen = set(self.selected_categories()) if self.category_buttons else set()
        self.category_buttons = []
        for i, cat in enumerate(categories):
            key = category_key_label(i)
            hint = f'{key}  ' if key else ''
            tip = f'{cat}  (key {key})' if key else f'{cat}  (no key)'
            btn = widgets.ToggleButton(value=cat in chosen, description=f'{hint}{cat}',
                                       tooltip=tip, disabled=self.session is None,
                                       layout=widgets.Layout(width='100%', height='28px'))
            btn.add_class('delfin-review-cat')
            btn._review_category = cat
            btn.observe(self._on_category, names='value')
            btn.observe(_paint_category, names='value')
            if btn.value:
                btn.button_style = 'warning'
            self.category_buttons.append(btn)
        self.category_box.children = tuple(self.category_buttons)

    def _set_enabled(self, on: bool):
        for w in (self.prev_btn, self.next_btn, self.unrated_btn, self.pass_btn,
                  self.block_btn, self.export_btn, self.note_input, self.unsure_btn,
                  self.h_toggle, self.reset_btn, self.add_cat_input, self.add_cat_btn,
                  self.file_comment, self.session_comment, *self.category_buttons):
            w.disabled = not on

    def _say(self, text: str, error: bool = False):
        color = _BLOCK_COLOR if error else '#555'
        self.status_html.value = (f'<span style="color:{color};font-size:12px">'
                                  f'{html.escape(text)}</span>')

    # -- the category list of a new review -------------------------------------
    def _field_untouched(self) -> bool:
        """The categories field still holds what the tab put there, not what was typed."""
        cats = _split_categories(self.categories_input.value)
        return cats == list(rv.DEFAULT_CATEGORIES) or (
            self._prefilled is not None and cats == self._prefilled)

    def _prefill_categories(self, raw: str | None = None) -> list[str] | None:
        """Fill the field from a ``kategorien.txt`` beside what is to be opened.

        Only while the field holds the built-in list (or an earlier fill of
        this kind); a list the reader typed is never replaced.  A resumed
        review keeps its own categories regardless (see :meth:`open`).
        """
        raw = (self.path_input.value if raw is None else raw).strip()
        if not raw or not self._field_untouched():
            return None
        try:
            path = Path(raw).expanduser()
            cats = rv.read_category_file(path) if path.exists() else None
        except OSError:
            cats = None
        target = cats or list(rv.DEFAULT_CATEGORIES)
        if _split_categories(self.categories_input.value) != target:
            self.categories_input.value = '\n'.join(target)
            if self.session is None:
                self._build_categories(target)
        self._prefilled = cats
        return cats

    # -- opening -------------------------------------------------------------
    def open(self, path: str | None = None) -> rv.ReviewSession | None:
        if path is not None:
            self.path_input.value = str(path)
        raw = self.path_input.value.strip()
        if not raw:
            self._say('Give a multi-frame .xyz file or a folder.', error=True)
            return None
        self._prefill_categories(raw)
        findings = self.findings_input.value.strip() or None
        wanted = _split_categories(self.categories_input.value)
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
            self._prefilled = list(session.categories)
            self.category_buttons = []
            self._build_categories(session.categories)
            self.session_comment.value = session.data.get('comment', '')
            self.list_toggle.value = not session.blind
        finally:
            self._busy = False
        self.sidebar.layout.display = '' if self.list_toggle.value else 'none'
        session.save()
        # Resume at the first frame still to be rated.
        if session.record(*session.order[0]) and session.next_unrated() is None:
            session.go(0)
        self._set_enabled(True)
        self.work.layout.display = ''
        self.setup.selected_index = None
        self.setup.set_title(0, f'Open a review  —  now: {session.source.name} '
                                f'→ {session.review_path.name}')
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
        """Use the field's categories.

        Before a review is open, and for an open review nobody has rated
        anything in yet, the field's list is taken as it is.  Once ratings
        exist, categories are only appended: the listed ones keep their place
        and so their key, and a category left out of the field stays.
        """
        cats = _split_categories(self.categories_input.value)
        if not cats:
            self._say('The category list is empty.', error=True)
            return
        if self.session is None:
            self._build_categories(cats)
            self._say(f'{len(cats)} categories set; a new review opened next uses them.')
            return
        s = self.session
        if not s.has_ratings():
            s.set_categories(cats)
            self._show_categories()
            self._say(f'{len(cats)} categories applied and saved in {s.review_path.name}.')
            return
        before = s.categories
        added = s.append_categories(cats)
        self._show_categories()
        msg = (f'{len(added)} categories appended' if added else 'No new category')
        msg += f'; saved in {s.review_path.name}.'
        if [c for c in cats if c in before] != before:
            msg += (' A review with ratings keeps its categories in their order '
                    '(their keys stay); new ones go to the end.')
        self._say(msg)

    def add_categories(self, names) -> list[str]:
        """Append categories (a list, or text separated by commas or lines) to the review."""
        if isinstance(names, str):
            names = re.split(r'[,\n]', names)
        names = rv.clean_categories(names)
        if not names:
            return []
        if self.session is None:
            have = _split_categories(self.categories_input.value)
            new = [c for c in names if c not in have]
            self.categories_input.value = '\n'.join(have + new)
            self._build_categories(have + new)
            return new
        added = self.session.append_categories(names)
        self._show_categories()
        self.add_cat_input.value = ''
        if added:
            self._say(f'{len(added)} categories appended; saved in '
                      f'{self.session.review_path.name}.')
        else:
            self._say('No new category.')
        return added

    def _show_categories(self):
        s = self.session
        self._busy = True
        try:
            self.categories_input.value = '\n'.join(s.categories)
            self._prefilled = list(s.categories)
        finally:
            self._busy = False
        self._build_categories(s.categories)
        n = len(s.categories)
        if n > MAX_KEYED_CATEGORIES:
            self._say(f'{n} categories: those after the {MAX_KEYED_CATEGORIES}th '
                      'have no key and are chosen with the mouse.')

    # -- the unsure mark and the view --------------------------------------------
    def set_unsure(self, on: bool):
        """Add or remove the word "unsicher" in this frame's comment (saved at once)."""
        if self.session is None:
            return
        text = self.note_input.value
        new = rv.with_unsure(text, on)
        if new != text:
            self.note_input.value = new                  # _on_note stores it
        self._sync_unsure()

    def toggle_unsure(self):
        self.set_unsure(not rv.is_unsure(self.note_input.value))

    def _sync_unsure(self):
        on = rv.is_unsure(self.note_input.value)
        was, self._busy = self._busy, True
        try:
            self.unsure_btn.value = on
            self.unsure_btn.button_style = 'warning' if on else ''
        finally:
            self._busy = was

    def reset_view(self):
        """Draw the frame again: orientation and zoom as when it was first shown."""
        if self.session is not None:
            self._draw_viewer(self.session.current)

    def handle_key(self, key: str, at: int | None = None):
        """Act on one key; *at* is the position the browser showed when it was pressed.

        A rating, category or unsure key whose *at* is not the current
        position is dropped: the kernel moved on (after a rating) before the
        screen did.  Navigation and view keys are relative and always act.
        """
        key = (key or '').strip().lower()
        if self.session is None:
            return
        cat = category_index(key)
        if (at is not None and at != self.session.position
                and (key in ('p', 'b', 'u') or cat is not None)):
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
        elif key == 'u':
            self.toggle_unsure()
        elif key == 'r':
            self.reset_view()
        elif key == 'h':
            self.h_toggle.value = not self.h_toggle.value
        elif cat is not None and cat < len(self.category_buttons):
            btn = self.category_buttons[cat]
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

    def _on_h_toggle(self, change):
        self.h_toggle.description = 'H atoms' if change.get('new') else 'H hidden'
        self.h_toggle.button_style = '' if change.get('new') else 'info'
        if self.session is not None:
            self._draw_viewer(self.session.current)

    def _on_unsure_btn(self, change):
        if self._busy or self.session is None:
            return
        self.set_unsure(bool(change.get('new')))

    def _on_mode(self, _change):
        if self.session is None or self._busy:
            return
        self.session.set_mode(self.blind_cb.value, int(self.seed_input.value or 0))
        self.session.save()
        self.list_toggle.value = not self.session.blind
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
        self._sync_unsure()
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
        self._sync_unsure()
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
        self.position_html.value = (f'<span style="font-size:14px;white-space:nowrap">{where}'
                                    f' &nbsp;&mdash;&nbsp; {state}</span>')
        self.verdict_html.value = f'<span style="font-size:12px">This frame: {state}</span>'

    def _draw_progress(self):
        p = self.session.progress()
        pct = 100.0 * p['rated'] / p['total'] if p['total'] else 0.0
        self.progress_html.value = (
            f'<div style="font-size:13.5px"><b>{p["rated"]} / {p["total"]} rated</b> '
            f'({pct:.0f} %), <span style="color:{_PASS_COLOR}"><b>{p["pass"]}</b> pass</span>, '
            f'<span style="color:{_BLOCK_COLOR}"><b>{p["block"]}</b> block</span></div>'
            f'<div style="background:#eee;height:5px;margin-top:2px;border-radius:3px">'
            f'<div style="background:#4a90d9;height:5px;width:{pct:.1f}%;'
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
        xyz = frame.xyz(comment)
        if not self.h_toggle.value:
            xyz = _without_h(xyz)
        try:
            from delfin.dashboard.molecule_viewer import render_xyz_in_output
            render_xyz_in_output(self.viewer, xyz, width=VIEWER_WIDTH, height=VIEWER_HEIGHT)
        except Exception as exc:                     # pragma: no cover - no browser
            with self.viewer:
                print(f'Viewer unavailable: {exc}')


def create_tab(ctx) -> Any:
    """Build the Review tab.  Returns ``(tab_widget, refs_dict)``."""
    start = getattr(ctx, 'calc_dir', None) or Path.home()
    panel = ReviewPanel(start, ctx=ctx)
    box = widgets.VBox([panel.widget], layout=widgets.Layout(width='100%', padding='6px 10px'))
    return box, {'review_panel': panel}


try:                                                 # pragma: no cover
    from delfin.dashboard.tab_registry import register_tab
    register_tab('review', 'Review', create_tab, order=8150)
except Exception:                                    # pragma: no cover
    pass
