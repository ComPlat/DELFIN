"""Three findings of the 2026-10-22 17:42 field report (KIT cluster, GLM).

1. An image in the chat has its size before it is decoded. Python
   replaces the whole transcript on every update; an inline <img> without
   width and height is 0 px tall until its data URI is decoded, so every
   update collapsed the plots, pictures and PDF previews above the end
   for a moment (the report carried a 19 MB PDF's page preview). With
   width and height from the image's own header the browser reserves the
   box at once and scales it with max-width and height: auto.
2. The working line does not lie on the transcript's scrollbar. It is
   laid over the transcript (#362) across the full cell, scrollbar
   included; the transcript's scrollbar is now a thin 10 px bar and the
   line keeps 18 px from the right edge.
3. Offering to read a file, or asking about it, is not a claim about it.
   "Im Workspace liegt u. a. `AUFGABE.md` ... -- soll ich da
   reinschauen?" drew "Note: cited but not read this turn: AUFGABE.md".
"""
from __future__ import annotations

import base64
import re
import struct
import zlib

from delfin.agent import verify_guard as vg


def _png(w: int, h: int) -> bytes:
    raw = b"".join(b"\x00" + b"\x00\x00\x00" * w for _ in range(h))

    def chunk(kind, data):
        return (struct.pack(">I", len(data)) + kind + data
                + struct.pack(">I", zlib.crc32(kind + data) & 0xffffffff))
    return (b"\x89PNG\r\n\x1a\n"
            + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


def _img_tags(html: str) -> list[str]:
    return re.findall(r"<img\b[^>]*>", html)


def test_a_picture_has_its_size_before_it_is_decoded(tmp_path):
    from delfin.dashboard import tab_agent
    p = tmp_path / "picture.png"
    p.write_bytes(_png(40, 25))
    [tag] = _img_tags(tab_agent._render_artifact_body(p))
    assert 'width="40"' in tag and 'height="25"' in tag


def test_a_vector_plot_has_its_size_before_it_is_decoded(tmp_path):
    from delfin.dashboard import tab_agent
    p = tmp_path / "plot.svg"
    p.write_text('<svg xmlns="http://www.w3.org/2000/svg" '
                 'width="576pt" height="345.6pt" viewBox="0 0 576 345.6">'
                 '<rect width="10" height="10"/></svg>')
    [tag] = _img_tags(tab_agent._render_artifact_body(p))
    assert 'width="576"' in tag and 'height="346"' in tag


def test_the_box_is_read_from_the_header_alone():
    from delfin.dashboard import tab_agent
    assert tab_agent._img_box(_png(7, 3)) == ' width="7" height="3"'
    assert tab_agent._img_box(b"GIF89a" + struct.pack("<HH", 5, 9)) == \
        ' width="5" height="9"'
    assert tab_agent._img_box(b"not an image") == ""
    svg = b'<svg viewBox="0 0 120 80" xmlns="http://www.w3.org/2000/svg"/>'
    assert tab_agent._img_box(svg) == ' width="120" height="80"'


def _rule(selector: str) -> str:
    from delfin.dashboard import tab_agent
    css = " ".join(re.sub(r"/\*.*?\*/", "", tab_agent._AGENT_CSS,
                          flags=re.S).split())
    out = ""
    for head, body in re.findall(r"([^{}]+)\{([^{}]*)\}", css):
        if selector in [s.strip() for s in head.split(",")]:
            out += body + ";"
    return out


def test_the_working_line_keeps_off_the_scrollbar():
    assert "width: 10px" in _rule(
        ".delfin-agent-chat-host::-webkit-scrollbar")
    assert "margin-right: 18px" in _rule(
        ".delfin-agent-chat-frame > .delfin-agent-working-slot")


def test_offering_to_read_a_file_is_not_a_claim_about_it(tmp_path):
    (tmp_path / "AUFGABE.md").write_text("task\n")
    (tmp_path / "manuscript.pdf").write_bytes(b"%PDF-1.4\n")
    offer = ("Was hast du vor? Im Workspace liegt u. a. `AUFGABE.md` und "
             "`manuscript.pdf` — soll ich da reinschauen?")
    assert vg.scan_for_ungrounded_code_claims(
        offer, repo_root=tmp_path, observed_files=set()) == []
    claim = "`AUFGABE.md` beschreibt eine Auswertung der Spektren."
    flags = vg.scan_for_ungrounded_code_claims(
        claim, repo_root=tmp_path, observed_files=set())
    assert [f.kind for f in flags] == ["unread"]


def test_the_png_helper_is_a_png():
    assert base64.b64encode(_png(1, 1))[:8] == b"iVBORw0K"
