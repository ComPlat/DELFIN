"""The whole of Ketcher as a panel: draw it, keep it, open it again.

The structure editor has had a small Ketcher beside its input box for a while:
a frame that opens under a DRAW toggle and one button that hands a molecule
back as a SMILES.  That is the right size for "I would rather draw this than
type it" and the wrong size for everything else the editor can do.

This is the editor as its own thing.  The frame fills the panel, what is drawn
is kept in ``<calc>/Ketcher`` under a name and opened again from there, and a
drawing with an arrow in it comes back as a reaction SMILES rather than as an
error -- Ketcher refuses to write a molfile for a reaction, which is why the
small version could not read one at all.

Three things travel between here and the browser, and each of them goes the way
the structure editor already established:

* out through ``ctx.run_js``, which is fine for a question but useless for an
  answer -- it clears the previous script before sending the next, so a reply
  written that way can be thrown away before the page has run it;
* back through a hidden textarea whose value the script sets, which is ordered,
  survives a background thread and arrives as an ordinary widget change;
* stamped with a serial, because a widget only reports a value that *changed*
  and drawing the same thing twice would otherwise look like an answer that
  never came.
"""

from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import Any, Dict

import ipywidgets as widgets
from IPython import get_ipython

from . import ketcher as _ketcher

__all__ = ['Panel', 'build', 'read_js', 'copy_js']

#: The classes the page finds the two halves of a panel by.  A dashboard can
#: hold more than one -- the Ketcher tab has one and each structure editor has
#: another -- so every selector is qualified by the panel's own scope, the way
#: the structure editor qualifies its frame.  Asking for "the" frame is how a
#: press in one tab was answered into another.
FRAME_CLASS = 'delfin-ketcher-frame'
SYNC_CLASS = 'delfin-ketcher-sync'

#: How long an answer from the page is waited for before the panel stops
#: claiming it is reading.  A frame that has been reloaded never answers at
#: all, and a status line that says "Reading the drawing..." for good is a lie
#: about what is happening.
_ANSWER_LEASH = 25.0

#: What to ask the editor for, per saved format.  ``auto`` is what the SMILES
#: button uses: a molfile normally, an RXN file when there is an arrow on the
#: canvas, because Ketcher throws "The structure cannot be saved as *.MOL due
#: to reaction" rather than writing one.
_ASK_FOR = {
    '.ket': 'ket', '.mol': 'mol', '.rxn': 'rxn', '.smi': 'smi',
    '.cdxml': 'cdxml',
}


def read_js(scope: str, kind: str, want: str, *,
            frame_class: str = FRAME_CLASS,
            sync_class: str = SYNC_CLASS) -> str:
    """Ask the editor for what has been drawn, in the form *want*.

    *kind* travels back with the answer so the one channel can carry both
    questions -- what to put in the SMILES box, and what to write to disk --
    without a second hidden widget for each.

    The two classes are arguments because the structure editor has a Ketcher
    of its own, under names of its own that predate this module and that its
    tests read.  The question being asked is the same one, so it is asked from
    one place; only which elements on the page it is asked through differs.
    """
    return (
        "(function(){\n"
        "  var scope=" + json.dumps(str(scope)) + ";\n"
        "  var kind=" + json.dumps(str(kind)) + ";\n"
        "  var want=" + json.dumps(str(want)) + ";\n"
        "  var MARK=" + json.dumps(_ketcher.KET_MARK) + ";\n"
        "  var box=document.querySelector('." + sync_class + ".'+scope);\n"
        "  var input=box&&box.querySelector('textarea, input');\n"
        "  function hand(text){\n"
        "    if(!input) return;\n"
        "    var proto=(input.tagName==='TEXTAREA')\n"
        "      ? window.HTMLTextAreaElement.prototype\n"
        "      : window.HTMLInputElement.prototype;\n"
        "    var setter=Object.getOwnPropertyDescriptor(proto,'value');\n"
        "    var line=(Date.now())+'\\n'+kind+'\\n'+text;\n"
        "    if(setter&&setter.set) setter.set.call(input,line);\n"
        "    else input.value=line;\n"
        "    input.dispatchEvent(new Event('input',{bubbles:true}));\n"
        "    input.dispatchEvent(new Event('change',{bubbles:true}));\n"
        "  }\n"
        "  var host=document.querySelector('." + frame_class + ".'+scope);\n"
        "  var frame=host&&host.querySelector('iframe');\n"
        "  var api=null;\n"
        "  try{ api=frame&&frame.contentWindow&&frame.contentWindow.ketcher; }\n"
        "  catch(e){ api=null; }\n"
        "  if(!api){ hand('!no-editor'); return; }\n"
        "  function ask(){\n"
        "    if(want==='ket') return api.getKet();\n"
        "    if(want==='smi') return api.getSmiles();\n"
        "    if(want==='cdxml') return api.getCDXml();\n"
        "    if(want==='rxn') return api.getRxn();\n"
        "    if(want==='mol') return api.getMolfile();\n"
        "    var arrow=false;\n"
        "    try{ arrow=!!(api.containsReaction&&api.containsReaction()); }\n"
        "    catch(e){ arrow=false; }\n"
        "    /* Always both. For a reaction because an RXN file holds one\n"
        "       arrow, so Indigo flattens a scheme drawn in three steps into\n"
        "       \"the first thing, into everything else\" and the arrows\n"
        "       survive only in Ketcher's own KET; and for a structure so\n"
        "       that the drawing itself can be kept beside the job it was\n"
        "       drawn for. */\n"
        "    return Promise.all([arrow ? api.getRxn() : api.getMolfile(),\n"
        "                        api.getKet()])\n"
        "      .then(function(both){ return both[0]+MARK+both[1]; });\n"
        "  }\n"
        "  try{\n"
        "    Promise.resolve(ask()).then(function(text){ hand(text||''); },\n"
        "      function(err){ hand('!'+err); });\n"
        "  }catch(e){ hand('!'+e); }\n"
        "})();"
    )


def copy_js(text: str) -> str:
    """Put *text* on the clipboard, however this browser will let us.

    The same ladder the structure editor's XYZ copy uses: the clipboard API
    where there is one, a hidden textarea and ``execCommand`` where there is
    not -- a dashboard reached over plain HTTP at a machine name has no
    clipboard API at all, because it is not a secure context.
    """
    payload = json.dumps(str(text or ''))
    return (
        "(function(){\n"
        "  var text=" + payload + ";\n"
        "  function old(){\n"
        "    try{\n"
        "      var ta=document.createElement('textarea');\n"
        "      ta.value=text; ta.style.position='fixed'; ta.style.opacity='0';\n"
        "      document.body.appendChild(ta); ta.focus(); ta.select();\n"
        "      document.execCommand('copy'); document.body.removeChild(ta);\n"
        "    }catch(e){ window.prompt('Copy:', text); }\n"
        "  }\n"
        "  try{\n"
        "    if(navigator.clipboard&&navigator.clipboard.writeText){\n"
        "      navigator.clipboard.writeText(text).catch(old);\n"
        "    } else { old(); }\n"
        "  }catch(e){ old(); }\n"
        "})();"
    )


class Panel:
    """One Ketcher, with somewhere to keep what is drawn in it."""

    def __init__(self, parts: Dict[str, Any]):
        self.__dict__.update(parts)

    @property
    def widget(self):
        return self.box


def build(ctx, *, height: str = '72vh', scope: str = 'delfin-ketcher-tab',
          title: str = 'Ketcher', folder=None, compact: bool = False,
          fill: bool = False) -> Panel:
    """One editor panel, ready to be placed.

    *folder* is where Ketcher's own Save writes and where its own Open reads,
    given as a callable because it moves: the Ketcher tab keeps its drawings
    in one place, and a drawing opened from the Calculations tab is saved back
    into the folder it came out of, the way a document is.

    *compact* leaves out what belongs to a drawing board rather than to a file
    being looked at -- the SMILES row and the handover to Submit.

    *fill* makes the frame take the height it is given rather than a fixed
    one, so that it reaches the bottom of the pane and follows it, which is
    what the text view, the grid and the document beside it do.
    """
    main_io_loop = getattr(getattr(get_ipython(), 'kernel', None),
                           'io_loop', None)

    def schedule(func, *args, **kwargs):
        """Touch a widget from the thread that owns it."""
        if main_io_loop is not None:
            main_io_loop.add_callback(lambda: func(*args, **kwargs))
            return
        func(*args, **kwargs)

    state: Dict[str, Any] = {'asked': None, 'serial': 0}
    frame_selector = f'.{FRAME_CLASS}.{scope}'
    sync_selector = f'.{SYNC_CLASS}.{scope}'

    # -- the parts ------------------------------------------------------
    # It says whatever it has to say, and gives way rather than pushing the
    # row wider than the pane.
    status = widgets.HTML(value='', layout=widgets.Layout(
        flex='1 1 0', min_width='0'))

    frame = widgets.HTML(value='', layout=widgets.Layout(
        width='100%',
        **({'flex': '1 1 0', 'min_height': '0', 'height': '100%'}
           if fill else {})))
    frame.add_class(FRAME_CLASS)
    frame.add_class(scope)

    # What the editor hands back.  A widget value rather than a script,
    # because a script sent through run_js can be replaced before the page
    # has run it.
    sync = widgets.Textarea(value='', layout=widgets.Layout(display='none'))
    sync.add_class(SYNC_CLASS)
    sync.add_class(scope)

    install_btn = widgets.Button(
        description='FETCH KETCHER', icon='download', button_style='warning',
        tooltip='Fetch the newest published Ketcher (about 32 MB).',
        layout=widgets.Layout(width='180px'),
    )

    smiles_btn = widgets.Button(
        description='TO SMILES', icon='arrow-down', button_style='success',
        tooltip=('Read what is drawn.  A structure comes back as a SMILES, a '
                 'drawing with an arrow in it as a reaction SMILES.'),
        layout=widgets.Layout(width='140px'),
    )
    # One box, whatever was drawn.  A structure and a reaction are read out
    # by the same press and go to the same places; two boxes meant one of them
    # was always holding something stale from an earlier drawing.
    # A box rather than a line: a scheme drawn in several steps comes back as
    # one reaction per step, one to a line.
    # A box rather than a line, and one that can be made narrower: at a flat
    # 100% beside two buttons of a fixed width the row came to more than the
    # width it had, and the pane grew a sideways scrollbar for the difference.
    smiles_out = widgets.Textarea(
        value='', placeholder='what was drawn, as a SMILES',
        description='SMILES:', rows=2,
        layout=widgets.Layout(width='auto', flex='1 1 0', min_width='0'),
        style={'description_width': '80px'},
    )
    smiles_copy_btn = widgets.Button(
        description='COPY', icon='copy', layout=widgets.Layout(width='90px'),
        tooltip='Put the SMILES on the clipboard.',
    )
    to_submit_btn = widgets.Button(
        description='TO SUBMIT', icon='arrow-right', button_style='info',
        tooltip='Put it in the Submit tab\'s input box and go there.',
        layout=widgets.Layout(width='130px'),
    )

    # -- saying things --------------------------------------------------
    def _say(text: str, colour: str = '#333') -> None:
        import html as _html
        status.value = (f"<span style='color:{colour}'>"
                        f"{_html.escape(str(text))}</span>")

    def _say_later(text: str, colour: str = '#333') -> None:
        schedule(_say, text, colour)

    # -- the folder -----------------------------------------------------
    def _folder() -> Path:
        """Where this panel keeps things, which is not always one place."""
        if folder is not None:
            try:
                chosen = folder() if callable(folder) else folder
            except Exception:                           # noqa: BLE001
                chosen = None
            if chosen:
                return Path(chosen)
        return _ketcher.drawings_directory(
            getattr(ctx, 'calc_dir', None) or (Path.home() / 'calc'))

    def _scan_files() -> list:
        """What is kept, remembered by name so an answer can be resolved."""
        kept = _ketcher.list_in(_folder())
        state['files'] = {item.name: item for item in kept}
        return [item.name for item in kept]

    def _refresh_files() -> None:
        """Tell the editor what is kept, so its own Open list can show it."""
        _send(_ketcher.files_js(frame_selector, _scan_files()))

    # -- talking to the page --------------------------------------------
    def _wire() -> str:
        """The wiring, restated. It marks the window it bound to, so arriving
        again costs nothing -- and the frame may have been built since the
        page's startup script ran, which is the case the moment Ketcher is
        fetched for the first time."""
        return '\n'.join([_ketcher.focus_js(frame_selector),
                          _ketcher.wire_js(frame_selector, sync_selector)])

    def _send(*scripts) -> None:
        """Everything this panel says to the page, said once.

        ``ctx.run_js`` clears its output before it displays the next script,
        so two calls in a row can mean the first is thrown away before the
        browser has run it.  Whatever has to happen together is therefore
        joined and sent as one -- which is also why the focus script travels
        with every question rather than being sent beside it.  It marks the
        element it bound to, so arriving again costs nothing.
        """
        joined = '\n'.join(str(part) for part in scripts if part)
        if joined:
            ctx.run_js(joined)

    # -- the frame ------------------------------------------------------
    def _show_frame() -> bool:
        """Put the editor on the page, if there is one to put there.

        Widgets only.  Nothing is sent to the page from here: this is called
        on the way into asking a question, and a script sent immediately
        before another one is a script that may never run.
        """
        url = _ketcher.app_url()
        if not url:
            frame.value = ''
            install_btn.description = 'FETCH KETCHER'
            install_btn.button_style = 'warning'
            for button in (smiles_btn, to_submit_btn):
                button.disabled = True
            return False
        # Kept rather than hidden: an editor that is there is an editor that
        # can be brought up to date, and the tab is the place to do it from.
        install_btn.description = 'UPDATE KETCHER'
        install_btn.button_style = ''
        for button in (smiles_btn, to_submit_btn):
            button.disabled = False
        # Set once.  Setting the same frame again reloads it, and a reloaded
        # editor is an empty one -- which is a drawing thrown away every time
        # anything on this panel is refreshed.
        if url not in (frame.value or ''):
            frame.value = _ketcher.frame_html(url, height=height)
        return True

    # -- asking the page ------------------------------------------------
    def _ask(kind: str, want: str, saying: str) -> None:
        state['asked'] = kind
        state['serial'] += 1
        mine = state['serial']
        _say(saying)
        _send(_wire(), read_js(scope, kind, want))

        def _leash():
            if state['serial'] == mine and state['asked'] == kind:
                state['asked'] = None
                _say_later('The editor did not answer.  If it was reloaded, '
                           'draw again and ask once more.', '#ef6c00')

        timer = threading.Timer(_ANSWER_LEASH, _leash)
        timer.daemon = True
        timer.start()

    def _on_smiles(_button=None) -> None:
        if not _show_frame():
            return
        _ask('smiles', 'auto', 'Reading the drawing...')

    # -- what came back -------------------------------------------------
    def _on_sync(change) -> None:
        if change.get('name') != 'value':
            return
        raw = sync.value or ''
        parts = raw.split('\n', 2)
        if len(parts) < 3:
            return
        _serial, kind, payload = parts
        state['asked'] = None
        if payload.startswith('!'):
            trouble = payload[1:]
            _say('The editor is not open yet, so there is nothing to read.'
                 if trouble == 'no-editor'
                 else f'The drawing could not be read: {trouble}', '#d32f2f')
            return
        if kind == 'save':
            # The name and the format are the ones typed into Ketcher's own
            # Save dialog: "my aspirin.mol" comes back whole.
            filename, _, body = payload.partition('\n')
            suffix = Path(filename).suffix.lower() or '.ket'
            outcome = _ketcher.save_into(
                _folder(), Path(filename).stem, body, suffix)
            _say(outcome['status'], '#2e7d32' if outcome['ok'] else '#d32f2f')
            if outcome['ok']:
                _refresh_files()
            return
        if kind == 'save-failed':
            _say(f'The drawing could not be read out of the editor: {payload}',
                 '#d32f2f')
            return
        if kind == 'open-list':
            _refresh_files()
            return
        if kind == 'open':
            where = (state.get('files') or {}).get(payload)
            if where is None:
                _refresh_files()
                _say(f'{payload} is not there any more.', '#d32f2f')
                return
            got = _ketcher.read_drawing(where)
            if not got['ok']:
                _refresh_files()
                _say(got['status'], '#d32f2f')
                return
            open_text(got['text'], got['name'])
            return
        outcome = _ketcher.smiles_from_drawing(payload)
        if not outcome['ok']:
            _say(outcome['status'], '#d32f2f')
            return
        # Two boxes, because a reaction and a structure are read differently
        # by whoever picks them up -- and because one of them would otherwise
        # be overwritten by the other the next time TO SMILES is pressed.
        smiles_out.value = outcome['smiles']
        state['reaction'] = bool(outcome.get('reaction'))
        _say(outcome['status'], '#2e7d32')

    # -- opening one that was kept --------------------------------------
    def open_text(text: str, name: str = '') -> bool:
        """Put *text* into the editor, whatever format it is written in.

        The entry point the Calculations browser calls, so that a drawing
        double-clicked over there arrives here as the drawing it was.
        """
        body = str(text or '')
        if not body.strip():
            _say(f'{name or "That file"} is empty.', '#d32f2f')
            return False
        if not _show_frame():
            _say('Ketcher is not installed yet, so there is nowhere to open '
                 'it.  Press FETCH KETCHER first.', '#d32f2f')
            return False
        _send(_wire(), _ketcher.load_js(frame_selector, body))
        _say(f'{name or "The drawing"} is in the editor.', '#2e7d32')
        return True

    def _on_to_submit(_button=None) -> None:
        """The structure into the Submit tab, where the job is set up.

        Carried across rather than copied here: every field a calculation
        needs is over there, and a second, staler copy of them would be worse
        than a tab switch.  The same handover the reaction graph makes.
        """
        reaction = bool(state.get('reaction'))
        drawn = str(smiles_out.value or '').strip()
        if not drawn:
            _say('Press TO SMILES first -- there is nothing to send.', '#d32f2f')
            return
        box = (getattr(ctx, 'submit_refs', None) or {}).get('coords_widget')
        if box is None:
            _say('There is no Submit tab in this dashboard to send it to.',
                 '#d32f2f')
            return
        box.value = drawn
        _say(f'{drawn} is in the Submit tab' + (
            ' -- a reaction, so Convert, which builds one structure, has '
            'nothing to do with it.' if reaction
            else ' -- Convert turns it into coordinates.'), '#2e7d32')
        try:
            ctx.select_tab('Submit Job')
        except Exception:                               # noqa: BLE001
            pass

    def _on_copy_smiles(_button=None) -> None:
        _send(copy_js(smiles_out.value or ''))
        _say('The SMILES is on the clipboard.')

    # -- fetching it ----------------------------------------------------
    def _on_install(_button=None) -> None:
        install_btn.disabled = True
        _say('Fetching Ketcher...')

        def _work():
            try:
                outcome = _ketcher.install(
                    on_line=lambda line: _say_later(f'Ketcher: {line}'))
            except Exception as exc:                    # noqa: BLE001
                outcome = {'ok': False,
                           'status': f'Ketcher could not be installed: {exc}'}

            def _done():
                install_btn.disabled = False
                _say(outcome['status'],
                     '#2e7d32' if outcome.get('ok') else '#d32f2f')
                if _show_frame():
                    _send(_wire(),
                          _ketcher.files_js(frame_selector, _scan_files()))

            schedule(_done)

        thread = threading.Thread(target=_work, daemon=True)
        thread.start()

    # -- wiring ---------------------------------------------------------
    smiles_btn.on_click(_on_smiles)
    smiles_copy_btn.on_click(_on_copy_smiles)
    to_submit_btn.on_click(_on_to_submit)
    install_btn.on_click(_on_install)
    sync.observe(_on_sync, names='value')

    def _row(members):
        return widgets.HBox(members, layout=widgets.Layout(
            gap='8px', align_items='center', flex_wrap='wrap',
            width='100%', max_width='100%', overflow='hidden'))

    box = widgets.VBox(
        [
            widgets.HTML(f'<b>{title}</b>') if title else widgets.HTML(''),
            _row([smiles_btn, install_btn, status]),
            frame, sync,
        ] + ([] if compact else
             [_row([smiles_out, smiles_copy_btn, to_submit_btn])]),
        layout=widgets.Layout(
            width='100%', gap='6px',
            **({'flex': '1 1 0', 'min_height': '0', 'height': '100%',
                'flex_flow': 'column'} if fill else {})),
    )

    ready = _show_frame()
    # Registered on the context rather than sent: create_dashboard collects
    # every tab's startup script and sends them as one, after all of the tabs
    # are built.  Sending it here would be a second run_js against a page that
    # has not been drawn yet -- and it would clear whatever the tab before it
    # had sent.
    #
    # The wiring goes out with it: Ketcher's Save and Open have to be answered
    # from the first press, not from the first time somebody happens to touch
    # a control beside the frame.
    try:
        ctx.add_init_js('\n'.join([
            _ketcher.focus_js(frame_selector),
            _ketcher.wire_js(frame_selector, sync_selector),
            _ketcher.files_js(frame_selector, _scan_files()),
        ]))
    except Exception:                                   # noqa: BLE001
        pass
    if ready:
        version = _ketcher.installed_version()
        _say(f'Ketcher {version}: draw it, keep it, open it again.')
    else:
        _say('Ketcher is not here yet.  It is about 32 MB, fetched once and '
             'then it works without a network -- press FETCH KETCHER.',
             '#ef6c00')

    return Panel(locals())
