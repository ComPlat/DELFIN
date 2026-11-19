"""Answer-time hint when the model hands the user hand-written viewer HTML.

The observed defect (2026-11-12): "show me benzene in 3D" was answered with
a raw ``<iframe srcdoc>`` embedding an NGL/3Dmol viewer loaded from a CDN,
written into the answer text. The agent tab renders its own interactive 3D
viewer card behind the ``show_molecule`` tool -- markup in the answer shows
as text on the card path and loads remote code on the HTML path, and the
model was never told either thing.

Three pieces live here:

* :func:`scan_raw_viewer_html` -- recognise a hand-written 3D viewer block
  (iframe with a srcdoc carrying a viewer script, or a script tag pointing
  at a known viewer library). Deliberately narrow: a ``<b>emphasis</b>`` or
  a code snippet must never fire it.
* :func:`raw_viewer_hint` -- the one visible hint naming the right tool.
* :func:`caveat_for` -- the combination the engine appends: only when the
  answer embeds a viewer AND ``show_molecule`` was not called this turn;
  the answer text is preserved, never rewritten.

Rejected: a forced correction turn (the card is reachable in the same turn
by a follow-up tool call, and a retry cannot un-print what the user already
saw); deleting the HTML from the answer (a rewrite of model text the user
is reading).
"""

from __future__ import annotations

import re
from typing import Optional

#: An <iframe> (any attribute order) -- the observed failure carried its
#: whole viewer inside ``srcdoc``, so the tag alone is the trigger.
_IFRAME = re.compile(r"<iframe\b", re.IGNORECASE)

#: A <script src> pointing at a known 3D-viewer library, with or without
#: an iframe around it.
_VIEWER_SCRIPT = re.compile(
    r"<script[^>]*src\s*=\s*[\"'][^\"']*(3dmol|ngl|glmol|jmol|molstar)"
    r"[^\"']*[\"']", re.IGNORECASE)

_HINT = (
    "[Show in 3D] This answer embeds a hand-written 3D viewer "
    "(<iframe>/<script>). The chat renders HTML as plain text and a viewer "
    "loaded from a CDN does not run. To show the molecule, call the "
    "`show_molecule` tool with the structure file (e.g. a written .xyz) -- "
    "the agent tab then renders its own interactive 3D viewer card."
)


def raw_viewer_hint() -> str:
    """The one hint, for tests and for the engine."""
    return _HINT


def scan_raw_viewer_html(text: str) -> Optional[str]:
    """Return the trigger found in ``text`` (\"iframe\" or \"script\"), or
    None. Narrow on purpose; see the module docstring."""
    if not text:
        return None
    if _IFRAME.search(text):
        return "iframe"
    if _VIEWER_SCRIPT.search(text):
        return "script"
    return None


def caveat_for(text: str, *, show_molecule_called: bool) -> str:
    """The hint for this answer, or \"\".

    Fires only when the answer embeds a raw viewer AND the turn never
    called ``show_molecule``: the tool result carries the card, and a hint
    on top of a delivered card would teach the model the working path is
    wrong.
    """
    if show_molecule_called:
        return ""
    if scan_raw_viewer_html(text) is None:
        return ""
    return "\n\n" + _HINT
