"""A view_image result is drawn in the dashboard chat, and its note says
who sees the image.

The model is sent the image as an image part of its next message; the
chat showed only the tool line, so an answer saying "the image is in the
chat" was false. The card is the chat's artifact preview of the resolved
path the executor read; anything else (an error, another shape, a path
that is not the resolved one, a type view_image does not accept) is no
card and is escaped like any tool output.
"""

from __future__ import annotations

import base64
import json
import struct
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions
from delfin.dashboard import tab_agent as T

# A 3x2 PNG, written byte for byte so the test needs no imaging library.
_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAMAAAACCAIAAAASFvFNAAAAEklEQVR4nGP4z8AARBAS"
    "TAAAzfAL9W1sCL0AAAAASUVORK5CYII=")


@pytest.fixture(autouse=True)
def _clear_pending():
    A._doc_executor.__dict__.pop("_pending_view_images", None)
    yield
    A._doc_executor.__dict__.pop("_pending_view_images", None)


def _view(ws: Path, rel: str) -> str:
    perms = KitToolPermissions(workspace=ws, mode="bypassPermissions")
    return A._doc_executor._execute_view_image({"path": rel}, perms)


def test_a_viewed_image_is_a_card_with_its_box_reserved(tmp_path):
    (tmp_path / "fig.png").write_bytes(_PNG)
    out = _view(tmp_path, "fig.png")
    card = T._view_image_card_html(out)
    assert card is not None
    assert "data:image/png;base64," in card
    # The width/height from the header reserve the box before the image
    # decodes, so the chat height does not change when it loads.
    assert 'width="3" height="2"' in card
    assert "fig.png" in card


def test_the_result_names_the_resolved_path(tmp_path):
    (tmp_path / "real").mkdir()
    (tmp_path / "real" / "a.png").write_bytes(_PNG)
    (tmp_path / "link").symlink_to(tmp_path / "real")
    out = json.loads(_view(tmp_path, "link/../real/a.png"))
    assert out["path"] == str((tmp_path / "real" / "a.png").resolve())


def test_an_error_result_is_no_card(tmp_path):
    (tmp_path / "a.txt").write_text("hello")
    assert T._view_image_card_html(_view(tmp_path, "a.txt")) is None
    assert T._view_image_card_html(_view(tmp_path, "nope.png")) is None
    assert T._view_image_card_html('{"error": "path is required"}') is None
    assert T._view_image_card_html("not json") is None
    assert T._view_image_card_html(None) is None


def test_a_forged_result_reaches_no_other_file(tmp_path):
    (tmp_path / "a.png").write_bytes(_PNG)
    (tmp_path / "notes.txt").write_text("secret")
    real = str((tmp_path / "a.png").resolve())

    def res(**kw):
        base = {"status": "ok", "path": real, "mime": "image/png"}
        base.update(kw)
        return json.dumps(base)

    assert T._view_image_card_html(res()) is not None
    # Another type than the tool accepts, or one that disagrees with the
    # file's suffix.
    assert T._view_image_card_html(
        res(path=str(tmp_path / "notes.txt"), mime="text/plain")) is None
    assert T._view_image_card_html(
        res(path=str(tmp_path / "notes.txt"))) is None
    assert T._view_image_card_html(res(mime="image/svg+xml")) is None
    # A relative path, or one that is not already resolved.
    assert T._view_image_card_html(res(path="a.png")) is None
    (tmp_path / "ln.png").symlink_to(tmp_path / "a.png")
    assert T._view_image_card_html(res(path=str(tmp_path / "ln.png"))) is None
    assert T._view_image_card_html(
        res(path=str(tmp_path / "x" / ".." / "a.png"))) is None
    assert T._view_image_card_html(res(status="error")) is None


def test_a_large_image_is_a_one_line_card(tmp_path):
    big = tmp_path / "big.png"
    big.write_bytes(_PNG + b"\0" * 2_100_000)
    card = T._view_image_card_html(_view(tmp_path, "big.png"))
    assert card is not None and "too large to embed" in card
    # No inline image; the existing download link is still offered.
    assert "<img" not in card


def test_a_webp_image_gets_its_box():
    # Lossless (VP8L) header of a 5x7 image: width-1 and height-1 in 14 bits.
    bits = (5 - 1) | ((7 - 1) << 14)
    data = (b"RIFF" + struct.pack("<I", 30) + b"WEBP" + b"VP8L"
            + struct.pack("<I", 10) + b"\x2f" + bits.to_bytes(4, "little")
            + b"\0" * 8)
    assert T._img_box(data) == ' width="5" height="7"'


def test_the_note_does_not_claim_the_user_sees_it(tmp_path):
    (tmp_path / "fig.png").write_bytes(_PNG)
    note = json.loads(_view(tmp_path, "fig.png"))["note"]
    assert note == A.VIEW_IMAGE_NOTE
    assert "shown to you" in note
    assert "dashboard chat" in note
    assert "do not see it unless" in note


def test_the_chat_draws_the_card_only_for_view_image():
    src = Path(T.__file__).read_text(encoding="utf-8")
    assert ('if tool_name == "view_image" and _raw_name in (\n'
            '                            "view_image", '
            '"mcp__delfin-docs__view_image"):') in src
    assert "_view_image_card_html(tool_output)" in src
