"""A viewer must begin at its own frame, not a prompt-width inside it.

An ipywidgets Output renders as a table: ``.jp-OutputArea-child`` is
``display: table`` with ``table-layout: fixed``, its first cell is the
``Out[7]:`` prompt and its second holds the content. Nothing in this dashboard
ever writes a prompt, but the theme still gives that empty cell
``--jp-cell-prompt-width`` -- 64px under Voila's light theme -- so every Output
panel started 64px in from its own border. On the structure viewers that is a
white band down the left inside the blue frame, with the picture pushed the
same 64px off the right edge and clipped there. In fullscreen, where there is
nothing beside it, it is all you see.

Measured in headless Chromium against Voila's own stylesheet and a real 3Dmol
canvas: 64px before, 0px after.

The rule lives in ``create_page_css`` and is deliberately not scoped to any
viewer class. Two of the Output viewers carry no usable class at all -- the
Turbomole one has no class whatsoever -- so a list of selectors cannot reach
them, and the next viewer someone adds would start out broken again.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

_HELPERS_PATH = (
    Path(__file__).resolve().parents[1] / 'delfin' / 'dashboard' / 'helpers.py'
)


def _load_create_page_css():
    """The real helpers module, whatever else is in sys.modules.

    Several tests in this suite install a three-function stub over
    delfin.dashboard.helpers and leave it there, so a plain import gives
    the stub or the real thing depending on what ran first. The file is
    executed under its proper name so that its `from .constants import ...`
    still resolves, and sys.modules is put back afterwards.
    """
    try:
        from delfin.dashboard.helpers import create_page_css
        return create_page_css
    except ImportError:
        pass

    name = 'delfin.dashboard.helpers'
    saved = sys.modules.get(name)
    spec = importlib.util.spec_from_file_location(name, _HELPERS_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
        return module.create_page_css
    finally:
        if saved is not None:
            sys.modules[name] = saved
        else:
            sys.modules.pop(name, None)


create_page_css = _load_create_page_css()

PROMPT_SELECTORS = ('.jp-OutputPrompt', '.jp-OutputArea-prompt')


def _rules(css):
    """(selector, body) for every rule in the stylesheet."""
    css = re.sub(r'</?style>', '', css)
    return [(m.group(1).strip(), m.group(2))
            for m in re.finditer(r'([^{}]+)\{([^{}]*)\}', css)]


def _prompt_rules():
    return [(sel, body) for sel, body in _rules(create_page_css().value)
            if any(p in sel for p in PROMPT_SELECTORS)]


def test_the_prompt_gutter_is_taken_away():
    rules = _prompt_rules()
    assert rules, (
        'nothing hides the Out[7] prompt cell, so every Output panel starts '
        'a prompt-width in from its own frame'
    )
    assert any('display' in body and 'none' in body for _sel, body in rules), (
        'the prompt cell is still laid out, so it still reserves its column'
    )


def test_the_rule_reaches_a_viewer_that_carries_no_class():
    """Two Output viewers have no class to hang a selector on."""
    for sel, _body in _prompt_rules():
        for part in sel.split(','):
            part = part.strip()
            if not part:
                continue
            # A bare .jp-OutputPrompt reaches every Output. A descendant
            # selector reaches only what someone remembered to mark.
            assert part in PROMPT_SELECTORS or part.lstrip().startswith('.jp-Output'), (
                f'{part!r} is scoped to an ancestor, so it misses the viewers '
                'that carry no class of their own'
            )


def test_the_content_cell_takes_the_whole_width():
    """With the prompt gone the picture must claim the space it freed."""
    rules = _rules(create_page_css().value)
    child = [body for sel, body in rules if '.jp-OutputArea-child' in sel]
    assert child, '.jp-OutputArea-child is never widened'
    assert any('width' in body and '100%' in body for body in child)


def test_the_separator_border_above_the_output_is_gone():
    """A notebook separates outputs with a transparent border; a viewer
    inside a frame does not want a line of its frame spent on that."""
    rules = _rules(create_page_css().value)
    output = [body for sel, body in rules if '.jp-OutputArea-output' in sel]
    assert output, '.jp-OutputArea-output is never adjusted'
    assert any('border-top' in body and '0' in body for body in output)
