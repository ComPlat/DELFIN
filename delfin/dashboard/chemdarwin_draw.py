"""Drawing the rules ChemDarwin runs on, instead of typing them.

Four boxes in that tab take chemistry as text -- a seed, the reaction SMARTS,
and the two pattern boxes that say what a product may not gain and may not
lose.  The SMARTS ones are the reason this module exists: writing one by hand
means stating the fragment that changes *and* the connectivity it has to sit
in, as a query expression, while counting the atom maps yourself.

So there is one Ketcher for all four, and a DRAW button beside each box that
points it at that box.  One rather than four: the editor is a 30 MB React and
WASM application, and four frames would be four of it.  It is built on the
first press rather than with the tab, for the same reason -- the Calculations
browser builds its panel that way too, and sends the startup script the
dashboard has already flushed.

The part with teeth is the two pattern boxes.  ``apply_custom_reaction_iter``
zips them against the reaction SMARTS *by line number*: line 3 of Forbidden
constrains rule 3, and a rule appended without a line appended beside it
quietly inherits the filter of a rule that has nothing to do with it.  So
appending a rule here pads both pattern boxes with the ``-`` that box already
understands to mean "nothing", and a pattern is never added on its own -- it
is added to the rule it belongs to, picked from a list.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import ipywidgets as widgets

from rdkit import Chem

from . import ketcher as _ketcher
from . import ketcher_panel as _panel
from . import ketcher_smarts as _smarts

__all__ = ['build_section', 'rule_lines', 'aligned', 'add_to_line',
           'set_line', 'SCOPE', 'LEGEND']

#: This tab's own Ketcher.  Every selector the panel uses is qualified by it,
#: so the editor in the Ketcher tab and this one never answer each other's
#: questions.
SCOPE = 'delfin-ketcher-chemdarwin'

#: What ``_parse_forbidden_line`` in the ChemDarwin engine reads as "no
#: constraint here".  Repeated rather than imported: importing the tab module
#: from a module the tab module imports is a cycle.
_NOTHING = ('-', 'none', 'n/a', '.', '*')

#: Per box: what it is called, what to ask the editor for, and what to make of
#: the answer.  ``structure`` is the road that was already there -- a molfile
#: read back by RDKit -- because a seed is a SMILES and nothing more.
_FIELDS: Dict[str, Dict[str, str]] = {
    'seed':      {'label': 'Seed SMILES',        'want': 'auto',   'as': 'structure'},
    'rxn':       {'label': 'Reaction SMARTS',    'want': 'smarts', 'as': 'reaction'},
    'forbidden': {'label': 'Forbidden patterns', 'want': 'smarts', 'as': 'pattern'},
    'protected': {'label': 'Protected patterns', 'want': 'smarts', 'as': 'pattern'},
}

#: The second road, when Ketcher will not write a SMARTS for what is drawn.
_FALLBACK_WANT = {'reaction': 'rxn', 'pattern': 'mol'}

_COLOURS = {'ok': '#2e7d32', 'note': '#ef6c00', 'bad': '#d32f2f', '': '#666'}

LEGEND = '''
<div style="font-size:12px; line-height:1.5;">
<p style="margin:0 0 6px 0;"><b>The bonds you draw define the pattern; Any
bonds and query properties define how it is allowed to be attached.</b></p>
<table style="border-collapse:collapse;">
<tr><th style="text-align:left; padding:2px 12px 2px 0;">In Ketcher</th>
    <th style="text-align:left; padding:2px 12px 2px 0;">In the SMARTS</th>
    <th style="text-align:left; padding:2px 0;">What it is for</th></tr>
<tr><td style="padding:2px 12px 2px 0;">Reaction Arrow Tool</td>
    <td style="padding:2px 12px 2px 0;"><code>&gt;&gt;</code></td>
    <td>separates before from after</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Reaction Mapping Tool</td>
    <td style="padding:2px 12px 2px 0;"><code>:1</code> <code>:2</code> …</td>
    <td><b>what changes</b> — the same number is the same atom</td></tr>
<tr><td style="padding:2px 12px 2px 0;">mapped on the left, not on the right</td>
    <td style="padding:2px 12px 2px 0;">—</td>
    <td>the atom is deleted (ring contraction)</td></tr>
<tr><td style="padding:2px 12px 2px 0;">unmapped on the right</td>
    <td style="padding:2px 12px 2px 0;">—</td>
    <td>the atom is built (annulation)</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Bond: Single / Aromatic</td>
    <td style="padding:2px 12px 2px 0;"><code>-</code> <code>:</code></td>
    <td>connectivity as drawn</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Bond: <b>Any</b></td>
    <td style="padding:2px 12px 2px 0;"><code>~</code></td>
    <td>attached, no matter how</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Bond: Single/Aromatic</td>
    <td style="padding:2px 12px 2px 0;"><code>-,:</code></td>
    <td>tolerates either Kekulé form</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Query property “Aromaticity”</td>
    <td style="padding:2px 12px 2px 0;"><code>a</code> / <code>A</code></td>
    <td>aromatic or aliphatic</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Query property “H count”</td>
    <td style="padding:2px 12px 2px 0;"><code>H1</code></td>
    <td>a free position — something may be attached here</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Query property “Ring bond count”</td>
    <td style="padding:2px 12px 2px 0;"><code>x&lt;n&gt;</code></td>
    <td>ring context</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Query property “Connectivity”</td>
    <td style="padding:2px 12px 2px 0;"><code>X&lt;n&gt;</code></td>
    <td>how substituted</td></tr>
<tr><td style="padding:2px 12px 2px 0;">Generic atom A / Q / X / M</td>
    <td style="padding:2px 12px 2px 0;">atom list</td>
    <td>leave the surroundings open</td></tr>
</table>
<p style="margin:6px 0 0 0; color:#666;">What is drawn arrives below as a
SMARTS and can be sharpened by hand there — the line is checked again on
every change.</p>
</div>
'''


# -- the text, on its own ------------------------------------------------
# Kept free of widgets so the line arithmetic can be tested without a kernel.

def rule_lines(text: Any) -> List[str]:
    """The lines the engine will actually treat as rules.

    Written to match ``apply_custom_reaction_iter``: blanks are dropped, a
    ``name:`` in front of a line that has no arrow is the name of nothing and
    the line is dropped with it, and only what still holds ``>>`` counts.  If
    this and the engine ever disagree, every pattern line below is off by the
    difference, so it is deliberately the same reading.
    """
    kept: List[str] = []
    for raw in str(text or '').splitlines():
        line = raw.strip()
        if not line:
            continue
        if '>>' not in line and ':' in line:
            line = line.split(':', 1)[1].strip()
        if '>>' in line:
            kept.append(line)
    return kept


def aligned(filter_text: Any, rules: int) -> str:
    """Pad a pattern box so line *n* still belongs to rule *n*.

    Only ever pads.  A box with more lines than there are rules is harmless --
    the engine reads past them -- and cutting it would throw away a filter
    somebody typed for a rule they are about to write.
    """
    lines = str(filter_text or '').splitlines()
    while len(lines) < rules:
        lines.append('-')
    return '\n'.join(lines)


def set_line(text: Any, index: int, value: str) -> str:
    """Put *value* on line *index*, making the line exist if it does not."""
    lines = str(text or '').splitlines()
    while len(lines) <= index:
        lines.append('-')
    lines[index] = value
    return '\n'.join(lines)


def add_to_line(text: Any, index: int, pattern: str) -> str:
    """Add *pattern* to the filter line for rule *index*.

    ``;`` is what the engine splits a line on, and ``-`` is what it reads as
    an empty one -- so a first pattern replaces the placeholder rather than
    being appended to it.
    """
    lines = str(text or '').splitlines()
    while len(lines) <= index:
        lines.append('-')
    standing = lines[index].strip()
    if not standing or standing.lower() in _NOTHING:
        lines[index] = pattern
        return '\n'.join(lines)
    parts = [p.strip() for p in standing.split(';') if p.strip()]
    if pattern not in parts:
        parts.append(pattern)
    lines[index] = ';'.join(parts)
    return '\n'.join(lines)


def _first_pattern(text: Any, index: int) -> str:
    """The first real pattern on a filter line, for putting back in the editor."""
    lines = str(text or '').splitlines()
    if index >= len(lines):
        return ''
    for part in lines[index].split(';'):
        part = part.strip()
        if part and part.lower() not in _NOTHING:
            return part
    return ''


# -- the panel -----------------------------------------------------------

def build_section(ctx, targets: Dict[str, Any]) -> Dict[str, Any]:
    """One Ketcher for the four boxes, and the buttons that point it at them.

    *targets* maps ``seed`` / ``rxn`` / ``forbidden`` / ``protected`` onto the
    Textareas in the tab.  Returns the section to place, the four DRAW buttons
    to put beside those boxes, and the things the tests reach for.
    """
    state: Dict[str, Any] = {'target': 'rxn', 'panel': None, 'quiet': False}

    head = widgets.HTML()
    holder = widgets.Box(layout=widgets.Layout(
        width='100%', min_width='0', overflow='hidden'))

    legend = widgets.Accordion(
        children=[widgets.HTML(LEGEND)],
        selected_index=None,
        layout=widgets.Layout(width='100%'))
    legend.set_title(0, 'How a drawing becomes a SMARTS')

    # Not continuously updated: every keystroke would re-run the validator,
    # and it runs a reaction over a dozen probe molecules.
    preview = widgets.Textarea(
        value='', placeholder='what was drawn, as a SMARTS',
        rows=3, continuous_update=False,
        layout=widgets.Layout(width='auto', flex='1 1 0', min_width='0'))
    verdict = widgets.HTML(value='')

    # On by default: ChemDarwin's seeds are aromatic scaffolds, and Ketcher
    # draws a benzene fragment with single and double bonds, which match none
    # of them.  Off gives exactly what was drawn.
    aromatic_box = widgets.Checkbox(
        value=True, description='bonds also match aromatic', indent=False,
        layout=widgets.Layout(width='auto', min_width='0'))
    aromatic_box.add_class('chemdarwin-draw-switch')

    line_pick = widgets.Dropdown(
        options=[], description='Rule:',
        layout=widgets.Layout(width='auto', flex='1 1 0', min_width='0'),
        style={'description_width': '50px'})

    read_btn = widgets.Button(
        description='READ DRAWING', icon='arrow-down',
        button_style='success', layout=widgets.Layout(width='190px'),
        tooltip='Fetch what is in the editor, as a SMARTS.')
    apply_btn = widgets.Button(
        description='APPEND', icon='plus', button_style='primary',
        layout=widgets.Layout(width='220px'))
    replace_btn = widgets.Button(
        description='REPLACE RULE', icon='pencil',
        layout=widgets.Layout(width='170px'))
    load_btn = widgets.Button(
        description='OPEN IN EDITOR', icon='arrow-up',
        layout=widgets.Layout(width='160px'),
        tooltip='Put the chosen line back in the editor as a drawing.')
    close_btn = widgets.Button(
        description='CLOSE', icon='times',
        layout=widgets.Layout(width='140px'))

    # -- saying things ---------------------------------------------------
    def _verdict(level: str, text: str) -> None:
        import html as _html
        if not text:
            verdict.value = ''
            return
        verdict.value = (f"<span style='color:{_COLOURS.get(level, '#666')}'>"
                         f"{_html.escape(str(text))}</span>")

    def _show(text: str) -> None:
        """Put *text* in the preview without the validator answering itself."""
        state['quiet'] = True
        try:
            preview.value = text
        finally:
            state['quiet'] = False

    # -- what came back --------------------------------------------------
    def _answer(kind: str, payload: str) -> bool:
        if not str(kind).startswith('cd-'):
            return False
        name = kind[3:]
        retried = name.endswith('!again')
        if retried:
            name = name[:-6]
        field = _FIELDS.get(name)
        if field is None:
            return False

        if payload.startswith('!'):
            trouble = payload[1:]
            if trouble == 'no-editor':
                _verdict('bad', 'The editor is not open yet.')
                return True
            # Ketcher's SMARTS export is not infallible and the version is not
            # pinned, so a refusal is a reason to ask a different way, once.
            if not retried and field['want'] == 'smarts':
                _verdict('note', f'Ketcher wrote no SMARTS ({trouble}) -- '
                                 f'asking again through the drawing file.')
                _panel_now().ask(f'cd-{name}!again',
                                 _FALLBACK_WANT[field['as']],
                                 'Reading the drawing again ...')
                return True
            _verdict('bad', f'The drawing could not be read: {trouble}')
            return True

        outcome = _read(field['as'], payload, retried)
        _show(outcome.get('smarts') or '')
        _verdict(outcome.get('level') or ('ok' if outcome.get('ok') else 'bad'),
                 outcome.get('status') or '')
        return True

    def _with_trial(outcome: Dict[str, Any]) -> Dict[str, Any]:
        """Say now what Run would only say afterwards, and without a reason."""
        if not outcome.get('ok'):
            return outcome
        tried = _smarts.trial_on_seed(outcome['smarts'],
                                      targets['seed'].value)
        if tried.get('status'):
            outcome['status'] = f"{outcome['status']} · {tried['status']}"
            if tried['level'] == 'note':
                outcome['level'] = 'note'
        return outcome

    def _read(shape: str, payload: str, retried: bool) -> Dict[str, Any]:
        """The drawing, as the text the box it is bound for holds."""
        if shape == 'structure':
            drawn = _ketcher.smiles_from_drawing(payload)
            if drawn.get('reaction'):
                # The seed box is read as one SMILES; a reaction SMILES put in
                # it fails much later, in the engine, as an unreadable seed.
                return {'ok': False, 'level': 'bad',
                        'smarts': drawn.get('smiles') or '',
                        'status': ('That is a reaction -- the seed is a '
                                   'single structure. Take the arrow out of '
                                   'the drawing.')}
            return {'ok': bool(drawn.get('ok')),
                    'smarts': drawn.get('smiles') or '',
                    'level': 'ok' if drawn.get('ok') else 'bad',
                    'status': drawn.get('status') or ''}
        if shape == 'reaction':
            raw = (_smarts.reaction_smarts_from_rxn_block(payload) if retried
                   else payload)
            if not raw:
                return {'ok': False, 'level': 'bad', 'smarts': '',
                        'status': 'No reaction could be read from the drawing file.'}
            outcome = _smarts.normalize_reaction_smarts(
                raw, aromatic=aromatic_box.value)
            return _with_trial(outcome)
        raw = (_smarts.query_smarts_from_molblock(payload) if retried
               else payload)
        if not raw:
            return {'ok': False, 'level': 'bad', 'smarts': '',
                    'status': 'No pattern could be read from the drawing file.'}
        return _smarts.normalize_query_smarts(
            raw, aromatic=aromatic_box.value)

    # -- the editor, once somebody wants it ------------------------------
    def _panel_now():
        """Build the editor on first use, and give it its startup script.

        The dashboard collects every tab's startup JavaScript and sends it as
        one, before this runs -- so what this panel registered has to be sent
        after the fact.  The Calculations browser builds its panel the same
        way and sends the same last entry.
        """
        if state['panel'] is None:
            panel = _panel.build(ctx, height='76vh', scope=SCOPE, title='',
                                 compact=True, on_answer=_answer)
            state['panel'] = panel
            # The panel's own TO SMILES reads the drawing into a box this
            # section does not use.  READ DRAWING is the same press with the
            # answer in the right place, so only one of them is shown.
            panel.smiles_btn.layout.display = 'none'
            holder.children = [panel.widget]
            fresh = '\n'.join(getattr(ctx, 'init_js_parts', [])[-1:])
            if fresh.strip():
                try:
                    ctx.run_js(fresh)
                except Exception:                           # noqa: BLE001
                    pass
        return state['panel']

    # -- the shape of the row, per box -----------------------------------
    def _relist() -> None:
        """The rules to pick from, as they stand in the box right now."""
        rules = rule_lines(targets['rxn'].value)
        options = [(f'{i + 1}: {line[:70]}', i) for i, line in enumerate(rules)]
        standing = line_pick.value
        line_pick.options = options
        if options:
            line_pick.value = (standing if standing in range(len(options))
                               else len(options) - 1)
        line_pick.disabled = not options

    def _retarget(name: str) -> None:
        state['target'] = name
        field = _FIELDS[name]
        head.value = (f"<b>Ketcher →</b> {field['label']}")
        pattern = field['as'] == 'pattern'
        seed = field['as'] == 'structure'
        line_pick.layout.display = 'none' if seed else ''
        replace_btn.layout.display = 'none' if seed else ''
        if seed:
            apply_btn.description = 'USE AS SEED'
        elif pattern:
            apply_btn.description = 'ADD TO THIS RULE'
            replace_btn.description = 'REPLACE LINE'
        else:
            apply_btn.description = 'APPEND AS NEW RULE'
            replace_btn.description = 'REPLACE RULE'
        _relist()
        _show('')
        _verdict('', '')

    def _picked() -> Optional[int]:
        chosen = line_pick.value
        return chosen if isinstance(chosen, int) else None

    # -- the buttons -----------------------------------------------------
    def _on_read(_button=None) -> None:
        panel = _panel_now()
        if not panel._show_frame():
            _verdict('bad', 'Ketcher is not here yet -- press FETCH '
                            'KETCHER, which fetches it once.')
            return
        name = state['target']
        _verdict('note', 'Reading the drawing ...')
        panel.ask(f'cd-{name}', _FIELDS[name]['want'], 'Reading the drawing ...')

    def _align() -> None:
        rules = len(rule_lines(targets['rxn'].value))
        for name in ('forbidden', 'protected'):
            box = targets.get(name)
            if box is None:
                continue
            padded = aligned(box.value, rules)
            if padded != (box.value or ''):
                box.value = padded

    def _checked(text: str) -> Optional[Dict[str, Any]]:
        name = state['target']
        if name == 'seed':
            mol = Chem.MolFromSmiles(text)
            if mol is None:
                mol = Chem.MolFromSmiles(text, sanitize=False)
            if mol is None:
                _verdict('bad', 'RDKit cannot read this SMILES.')
                return None
            return {'ok': True, 'status': f'{mol.GetNumAtoms()} atoms'}
        found = _smarts.inspect(text, reaction=(name == 'rxn'))
        if not found['ok']:
            _verdict('bad', found['status'])
            return None
        return found

    def _on_apply(_button=None) -> None:
        text = (preview.value or '').strip()
        if not text:
            _verdict('bad', 'There is nothing in the box below to take.')
            return
        found = _checked(text)
        if found is None:
            return
        name = state['target']
        if name == 'seed':
            targets['seed'].value = text
            _verdict('ok', f"Taken as the seed · {found['status']}")
            return
        if name == 'rxn':
            rules = rule_lines(targets['rxn'].value)
            rules.append(text)
            targets['rxn'].value = '\n'.join(rules)
            _align()
            _relist()
            _verdict('ok', f"Appended as rule {len(rules)} · {found['status']}")
            return
        index = _picked()
        if index is None:
            _verdict('bad', 'A pattern always belongs to a rule -- add a '
                            'reaction SMARTS first.')
            return
        box = targets[name]
        box.value = add_to_line(box.value, index, text)
        _align()
        _verdict('ok', f'Added to rule {index + 1}.')

    def _on_replace(_button=None) -> None:
        text = (preview.value or '').strip()
        if not text:
            _verdict('bad', 'There is nothing in the box below to take.')
            return
        index = _picked()
        if index is None:
            _verdict('bad', 'No rule chosen.')
            return
        found = _checked(text)
        if found is None:
            return
        name = state['target']
        if name == 'rxn':
            rules = rule_lines(targets['rxn'].value)
            if index >= len(rules):
                _verdict('bad', f'Rule {index + 1} is no longer there.')
                return
            rules[index] = text
            targets['rxn'].value = '\n'.join(rules)
        else:
            targets[name].value = set_line(targets[name].value, index, text)
        _align()
        _relist()
        _verdict('ok', f"Rule {index + 1} replaced · {found['status']}")

    def _on_load(_button=None) -> None:
        name = state['target']
        if name == 'seed':
            text = (targets['seed'].value or '').strip()
        else:
            index = _picked()
            if index is None:
                _verdict('bad', 'No rule chosen.')
                return
            if name == 'rxn':
                rules = rule_lines(targets['rxn'].value)
                text = rules[index] if index < len(rules) else ''
            else:
                text = _first_pattern(targets[name].value, index)
        if not text:
            _verdict('bad', 'There is nothing there that could be drawn.')
            return
        # The line stays in the preview, whatever the editor makes of it: the
        # trip in is lossy and this is the only copy of what was written.
        _show(text)
        panel = _panel_now()

        # A rule goes over as an RXN, never as the SMARTS it is written as.
        # ``setMolecule`` reads a SMARTS and drops every atom map on the way
        # in, and a rule whose maps are gone is a different rule -- the maps
        # would then be guessed back by the MCS on the way out.
        sending, lossy = text, False
        if name == 'rxn':
            block = _smarts.rxn_block_from_smarts(text)
            if block:
                sending = block
                lossy = not _smarts.survives_the_editor(text)
        if not panel.open_text(sending, name='The line'):
            _verdict('bad', 'The editor could not open the line -- the status '
                            'line above the editor says why.')
            return
        if name == 'rxn' and lossy:
            _verdict('note', 'In the editor -- but Ketcher cannot hold this '
                             'rule\'s query terms, the hydrogen count among '
                             'them. Read back it would be a different rule. '
                             'The original is untouched in the box below.')
            return
        _verdict('ok', 'In the editor. Change it, then READ DRAWING.')

    def _on_close(_button=None) -> None:
        box.layout.display = 'none'

    def _revalidate(change=None) -> None:
        if state['quiet']:
            return
        text = (preview.value or '').strip()
        if not text:
            _verdict('', '')
            return
        name = state['target']
        if name == 'seed':
            mol = Chem.MolFromSmiles(text) or Chem.MolFromSmiles(text, sanitize=False)
            if mol is None:
                _verdict('bad', 'RDKit cannot read this SMILES.')
            else:
                _verdict('ok', f'{mol.GetNumAtoms()} atoms')
            return
        found = _smarts.inspect(text, reaction=(name == 'rxn'))
        if found['ok'] and name == 'rxn':
            found = _with_trial(dict(found, smarts=text))
        _verdict(found['level'], found['status'])

    read_btn.on_click(_on_read)
    apply_btn.on_click(_on_apply)
    replace_btn.on_click(_on_replace)
    load_btn.on_click(_on_load)
    close_btn.on_click(_on_close)
    preview.observe(_revalidate, names='value')
    aromatic_box.observe(lambda _c: _revalidate(), names='value')
    targets['rxn'].observe(lambda _c: _relist(), names='value')

    def _row(members):
        return widgets.HBox(members, layout=widgets.Layout(
            gap='8px', align_items='center', flex_wrap='wrap',
            width='100%', max_width='100%', overflow='hidden'))

    box = widgets.VBox(
        [
            _row([head, close_btn]),
            legend,
            holder,
            _row([read_btn, aromatic_box, line_pick]),
            preview,
            verdict,
            _row([apply_btn, replace_btn, load_btn]),
        ],
        layout=widgets.Layout(width='100%', gap='6px', display='none',
                              border='1px solid #d0d0d0', border_radius='6px',
                              padding='8px', overflow_x='hidden'),
    )

    box.add_class('chemdarwin-draw')

    def _opener(name: str):
        def go(_button=None) -> None:
            # Under the box it belongs to, not at the foot of the tab.  The
            # tab hands in how to do that, because only it knows what the two
            # containers hold; without one the section stays where it is.
            place = state.get('place')
            if place is not None:
                try:
                    place(name)
                except Exception:                           # noqa: BLE001
                    pass
            box.layout.display = ''
            _retarget(name)
            _panel_now()._show_frame()
        return go

    def _set_place(where) -> None:
        """Told, once the tab has built the containers, how to move."""
        state['place'] = where

    buttons = {
        name: widgets.Button(
            description='DRAW', icon='pencil', button_style='info',
            layout=widgets.Layout(width='100px'),
            tooltip=f"Draw {field['label']} instead of typing it.")
        for name, field in _FIELDS.items()
    }
    for name, button in buttons.items():
        button.on_click(_opener(name))

    _retarget('rxn')

    return {'widget': box, 'buttons': buttons, 'state': state,
            'preview': preview, 'verdict': verdict, 'line_pick': line_pick,
            'read': _on_read, 'apply': _on_apply, 'replace': _on_replace,
            'load': _on_load, 'retarget': _retarget, 'answer': _answer,
            'align': _align, 'relist': _relist, 'set_place': _set_place,
            'aromatic': aromatic_box}
