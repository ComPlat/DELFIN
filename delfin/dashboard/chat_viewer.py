"""Package V1 — dashboard rendering for molecule chat media.

The shared chat hook (agreed with package V2) is the ``role="tool"`` raw-HTML
render path: a tool result that is HTML is inlined in the chat verbatim. This
module builds the string that hook shows for a molecule, and produces the
plain-text fallback for a terminal.

``render_molecule`` is the one entry point the protected ``show_molecule``
tool handler calls: it either returns the self-contained, rotatable 3Dmol
HTML card (dashboard, ``as_html=True``) or the plain-text summary a terminal
can show (``as_html=False``). The heavy lifting (read gate, XYZ/cube parsing,
HTML block) lives in :mod:`delfin.agent.chat_media`.
"""
from __future__ import annotations

from delfin.agent import chat_media

_DRAG_HINT = "drag to rotate · scroll to zoom · right-drag to translate"


def _card_html(media: dict, title: str | None = None) -> str:
    """Wrap a chat_media HTML block in a titled, hint-carrying chat card."""
    display = title or (media.get("path") is not None and media.get("path")) \
        or media["formula"]
    safe_title = (display or "").replace("<", "&lt;").replace(">", "&gt;")
    kind_label = {
        "xyz": "XYZ",
        "multixyz": "multi-frame XYZ",
        "cube": "cube isosurface",
    }.get(media["kind"], media["kind"])
    return (
        '<div class="delfin-molecule" style="margin:6px 0;padding:6px;'
        'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
        f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
        f'&#128231; <code>{safe_title}</code> &middot; {kind_label}</div>'
        f'{media["html"]}'
        f'<div style="font-size:10px;color:#9ca3af;margin-top:4px;">'
        f'{_DRAG_HINT}</div>'
        '</div>'
    )


def render_molecule(
    path_or_xyz: str,
    *,
    kind: str | None = None,
    isovalue: float | None = None,
    workspace_root: str | None = None,
    as_html: bool = True,
    title: str | None = None,
) -> str:
    """The string the chat (or terminal) shows for one molecule.

    Args:
        path_or_xyz: a path to an XYZ/cube file, or the raw XYZ/cube text.
        kind: ``"xyz"`` | ``"multixyz"`` | ``"cube"`` (auto-detected when None).
        isovalue: cube isosurface level (default 0.02).
        workspace_root: root file paths must stay inside (the read gate).
        as_html: True -> a self-contained, rotatable 3Dmol HTML card for the
            dashboard's ``role="tool"`` hook; False -> the plain-text fallback
            (formula, atom count, path) a terminal prints.
        title: optional title shown in the card header (defaults to the file
            path or formula).

    Raises:
        ValueError: unknown kind, path outside ``workspace_root``, or a file
            over the size cap.
    """
    media = chat_media.molecule(
        path_or_xyz,
        kind=kind,
        isovalue=isovalue,
        workspace_root=workspace_root,
    )
    if not as_html:
        return media["text"]
    return _card_html(media, title)
