"""The card tool results (``show_molecule``, ``make_plot``): what they
carry, and how they are read.

Input: a molecule file the executor's read gate has already admitted.
Output: the result string ``DELFIN_CARD:`` + JSON::

    {"card": "molecule", "path": <absolute resolved path>, "title": <str>,
     "isovalue": <float or null>, "text": <plain summary line>}

The result carries no markup. The dashboard draws the card itself from
``path`` with the chat's existing 3D viewer (``tab_agent._mol3d_card_html``:
bundled 3Dmol, so it works on a node without outbound network; trajectories
animated; cubes drawn with atoms and both isosurface signs). The model and
the terminal are given ``text`` only (:func:`card_model_text`), never the
file content.

``make_plot`` (``delfin.agent.chat_plots.card_result``) writes a figure and
returns the same marker with ``{"card": "plot", "path", "title", "text"}``;
the dashboard shows that file with the chat's artifact preview. A card may
be followed by a newline and further lines for the model (an energy
profile's source notes): the dashboard parses the JSON object alone, so the
object stays within the 2000-character preview of a result however long
those lines are, and :func:`card_model_text` appends them to ``text``.

Rejected: an ``<iframe sandbox srcdoc>`` card built here. The sandboxed
document cannot see the dashboard's bundled 3Dmol and fetched an unpinned
build from 3Dmol.org, so on a node without network the card was an empty
460 px box; it also copied the whole file (escaped twice) into the tool
result, and duplicated the chat's existing viewer with a weaker one (first
frame only, cube surface without atoms or the negative lobe).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from delfin.agent import chat_media

#: Marker prefix of a card result. Only a result of the ``show_molecule``
#: or ``make_plot`` tool is read as a card, each as its own kind (see
#: tab_agent._show_molecule_card_html, _make_plot_card_html); the same
#: prefix in any other tool's output is ordinary text.
CARD_MARKER = "DELFIN_CARD:"

#: The suffixes the chat's 3D card draws as XYZ or cube
#: (tab_agent._mol3d_payload; .pdb is drawn there but not parsed here).
DRAWN_SUFFIXES = (".xyz", ".trj", ".traj", ".cube")


def _isovalue(value) -> "float | None":
    """A usable isosurface level, or None for the default.

    The model may send a string, a bool, zero, a negative or a non-finite
    number. Both signs are drawn, so the magnitude is what matters.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    if not math.isfinite(value) or value == 0:
        return None
    return abs(float(value))


def render_molecule_tool_result(
    path,
    *,
    isovalue=None,
    workspace_root: "str | None" = None,
    title: "str | None" = None,
) -> str:
    """The ``DELFIN_CARD:`` result string for one molecule file.

    Raises:
        ValueError: path outside ``workspace_root``, a file over the size
            cap, or content that is neither XYZ nor cube.
    """
    if "\n" in str(path):
        raise ValueError("takes a file path, not file content")
    suffix = Path(str(path)).suffix.lower()
    if suffix not in DRAWN_SUFFIXES:
        raise ValueError(
            f"draws {', '.join(DRAWN_SUFFIXES)} files, "
            f"not {suffix or 'a file without a suffix'}")
    media = chat_media.molecule(
        path, kind="cube" if suffix == ".cube" else None,
        workspace_root=workspace_root)
    resolved = str(Path(media["path"]).expanduser().resolve())
    payload = {
        "card": "molecule",
        "path": resolved,
        "title": str(title or media["path"]),
        "isovalue": _isovalue(isovalue),
        "text": media["text"],
    }
    return CARD_MARKER + json.dumps(payload)


def card_payload(result, card: str = "molecule") -> "dict | None":
    """The validated payload of a ``card`` result, or None.

    ``card`` is the kind the caller expects from the tool that produced
    ``result``: ``"molecule"`` (show_molecule) or ``"plot"`` (make_plot,
    ``delfin.agent.chat_plots.card_result``). None for anything that does
    not start with the marker, is of another kind, or whose JSON is not
    exactly that kind's shape (wrong types, missing keys).
    """
    split = _split_card(result)
    if split is None:
        return None
    payload = split[0]
    if not isinstance(payload, dict) or payload.get("card") != card:
        return None
    for key in ("path", "title", "text"):
        if not isinstance(payload.get(key), str):
            return None
    if card == "molecule":
        iso = payload.get("isovalue")
        if iso is not None and _isovalue(iso) is None:
            return None
    return payload


def _split_card(result) -> "tuple[object, str] | None":
    """``(decoded JSON, lines after it)`` of a card result, or None.

    The JSON object follows the marker and ends the result or a line; what
    follows that newline is plain text for the model. Anything else after
    the object makes the result no card.
    """
    if not isinstance(result, str) or not result.startswith(CARD_MARKER):
        return None
    try:
        payload, end = json.JSONDecoder().raw_decode(
            result, len(CARD_MARKER))
    except (ValueError, TypeError):
        return None
    tail = result[end:]
    if tail and not tail.startswith("\n"):
        return None
    return payload, tail[1:]


def card_model_text(result, card: str = "molecule") -> "str | None":
    """The plain summary the model receives for a card result, or None:
    ``text``, then the lines that follow the JSON object."""
    payload = card_payload(result, card)
    if payload is None:
        return None
    tail = _split_card(result)[1]
    return payload["text"] + ("\n" + tail if tail else "")
