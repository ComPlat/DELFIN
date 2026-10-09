"""Package V1 — dashboard rendering for molecule chat media.

The shared chat hook (agreed with package V2) is the ``role="tool"`` hook in
``delfin/dashboard/tab_agent.py``: a tool result that begins with the
``DELFIN_CARD:`` marker is inlined as a self-contained HTML card; everything
else is escaped like normal tool output (so an error that echoes a hostile
argument can never execute). This module builds that card and the plain-text
terminal fallback, and produces the ``DELFIN_CARD:`` result string the
protected ``show_molecule`` handler returns.

The card is an isolated ``<iframe sandbox="allow-scripts" srcdoc="...">``
WITHOUT ``allow-same-origin``: the 3Dmol viewer runs inside the iframe but
cannot reach the parent dashboard, and the ``srcdoc`` value is attribute-
escaped so file content (quotes, ``</script>``, ``</iframe>``) can never
break out of the iframe. The heavy lifting (read gate, XYZ/cube parsing) lives
in :mod:`delfin.agent.chat_media`.
"""
from __future__ import annotations

import html as _html
import json
import re

from delfin.agent import chat_media

#: Fixed marker prefix of the shared V1+V2 tool-card contract. The tab_agent
#: hook inlines ONLY output that starts with this and parses to the object
#: below; anything else is escaped like any tool result.
CARD_MARKER = "DELFIN_CARD:"

#: The shared card shape both show_molecule (V1) and make_plot (V2) emit.
#: ``escape`` is always "html"; ``html`` is the handler-built, already-escaped
#: sandboxed iframe card; ``text`` is the one-line plain terminal fallback.
CARD_KEYS = ("escape", "html", "text")

_DRAG_HINT = "drag to rotate · scroll to zoom · right-drag to translate"


def tool_result(html: str, text: str) -> str:
    """The tool result string for the shared contract: ``CARD_MARKER`` + JSON.

    The tab_agent hook recognises the marker, inlines ``html`` (the sandboxed
    iframe card) and ignores the rest for display. ``text`` is what a
    terminal / text-only reader shows.
    """
    payload = json.dumps({"escape": "html", "html": html, "text": text})
    return CARD_MARKER + payload


def _card_html(media: dict, title: str | None = None) -> str:
    """A titled, hint-carrying chat card holding a sandboxed 3Dmol viewer.

    The viewer media (a div + <script> that loads 3Dmol from the CDN) runs
    inside ``<iframe sandbox="allow-scripts" srcdoc="...">``. ``srcdoc`` is an
    attribute, so the whole nested document is HTML-escaped (quotes included):
    file content can never terminate the attribute or the iframe, and the
    iframe has no ``allow-same-origin`` so it cannot touch the parent page.
    """
    display = title or (media.get("path") is not None and media.get("path")) \
        or media["formula"]
    safe_title = _html.escape(str(display or ""))
    kind_label = {
        "xyz": "XYZ",
        "multixyz": "multi-frame XYZ",
        "cube": "cube isosurface",
    }.get(media["kind"], media["kind"])
    viewer_doc = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'></head>"
        f"<body style='margin:0'>{media['html']}</body></html>"
    )
    srcdoc = _html.escape(viewer_doc, quote=True)
    return (
        '<div class="delfin-molecule" style="margin:6px 0;padding:6px;'
        'border:1px solid #e5e7eb;border-radius:6px;background:#fff;">'
        f'<div style="font-size:11px;color:#6b7280;margin-bottom:4px;">'
        f'&#128231; <code>{safe_title}</code> &middot; '
        f'<span style="color:#9ca3af">{kind_label}</span></div>'
        '<iframe sandbox="allow-scripts" srcdoc="' + srcdoc + '" '
        'style="width:100%;height:460px;border:0;border-radius:4px;'
        'display:block;" title="molecule viewer"></iframe>'
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
        as_html: True -> the sandboxed-iframe HTML card for the dashboard's
            ``role="tool"`` hook; False -> the plain-text fallback (formula,
            atom count, path) a terminal prints.
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


def render_molecule_tool_result(
    path_or_xyz: str,
    *,
    kind: str | None = None,
    isovalue: float | None = None,
    workspace_root: str | None = None,
    title: str | None = None,
) -> str:
    """The full ``DELFIN_CARD:`` result string for the ``show_molecule`` tool.

    Wraps the sandboxed iframe card (dashboard) and the plain-text fallback
    (terminal) in the shared V1+V2 contract. The protected handler returns
    this; tab_agent inlines ``html`` when the marker is present and escapes
    everything else.
    """
    media = chat_media.molecule(
        path_or_xyz,
        kind=kind,
        isovalue=isovalue,
        workspace_root=workspace_root,
    )
    return tool_result(_card_html(media, title), media["text"])


def _is_safe_inline_card(html: str) -> bool:
    """Defence-in-depth gate on the inline card html (operator finding #2).

    Accepts ONLY the handler-built card shape: exactly one
    ``<iframe sandbox="allow-scripts" srcdoc="...">`` with no ``src=``, no
    ``on*`` event handlers and no ``allow-same-origin`` on that iframe, and no
    injection-capable tag (script/img/svg/object/embed/...) anywhere OUTSIDE
    the fully attribute-escaped ``srcdoc``. The file content lives only inside
    ``srcdoc``, escaped, so nothing here can break out of the iframe. Every
    other shape is rejected so tab_agent never inlines attacker-shaped html.
    """
    if html.count("<iframe") != 1:
        return False
    _m = re.search(r"<iframe\b([^>]*)>", html)
    if not _m:
        return False
    if 'sandbox="allow-scripts"' not in _m.group(1):
        return False
    # Strip the srcdoc attribute before inspecting the iframe's own
    # attributes: its escaped value legitimately contains "src=" for the
    # embedded 3Dmol script and must not trip the bans below. Only a real
    # attribute ON the iframe itself (allow-same-origin, src=, on*=) is
    # grounds for rejection.
    _attrs_plain = re.sub(r'srcdoc="[^"]*"', "", _m.group(1))
    if "allow-same-origin" in _attrs_plain:
        return False
    # The iframe must carry ONLY srcdoc + style/title — never src or a JS
    # event handler attribute.
    if re.search(r"\b(?:src|on[a-z]+)\s*=", _attrs_plain):
        return False
    # Drop the escaped srcdoc (quotes are &quot;, so it contains no literal
    # "</"), then reject any remaining injection-capable tag in the wrapper.
    # The single validated iframe (count==1 + per-attribute checks above) is
    # the only <iframe> left and is already proven safe.
    _stripped = re.sub(r"srcdoc=\"[^\"]*\"", "", html)
    for _tag in ("<script", "<img", "<svg", "<object", "<embed",
                 "<link", "<meta", "<form", "<style"):
        if _tag in _stripped:
            return False
    return True


def _card_payload(result: str) -> "dict | None":
    """The parsed, validated contract object for a ``DELFIN_CARD:`` result.

    ``None`` when the result is not a card (error output that echoes a hostile
    path, a non-card shape, or malformed content).
    """
    if not result.startswith(CARD_MARKER):
        return None
    try:
        payload = json.loads(result[len(CARD_MARKER):])
    except (ValueError, TypeError):
        return None
    if not isinstance(payload, dict):
        return None
    if payload.get("escape") != "html" or not isinstance(payload.get("html"), str):
        return None
    if not _is_safe_inline_card(payload["html"]):
        return None
    return payload


def parse_card_result(result: str) -> "str | None":
    """Return the inline html card from a ``DELFIN_CARD:`` result, else None.

    This is the reader half of the shared V1+V2 contract. Only a result that
    starts exactly with :data:`CARD_MARKER`, parses to
    ``{"escape":"html","html":<str>, ...}`` AND whose ``html`` is the exact
    single-sandboxed-iframe card shape (see :func:`_is_safe_inline_card`) is
    accepted; an error output (JSON that echoes a hostile path), a non-card
    shape, or malformed content returns None so the caller escapes the output
    rather than inlining it.
    """
    payload = _card_payload(result)
    return payload["html"] if payload is not None else None


def card_model_text(result: str) -> "str | None":
    """The short model-facing summary of a card result, or None.

    The model's tool-result context must NEVER receive the embedding html or
    raw file content (operator finding #1): only the one-line plain text
    summary (formula, atom count, format, path) built by chat_media from
    parsed data. For a result that is not a valid card this returns None so
    the caller keeps its existing (escaped) fallback behaviour.
    """
    payload = _card_payload(result)
    if payload is None:
        return None
    text = payload.get("text")
    if not isinstance(text, str):
        return None
    return text
