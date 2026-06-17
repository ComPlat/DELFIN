"""Agent view_image tool — let a vision model look at an image from disk."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions


def _png(path: Path) -> None:
    Image = pytest.importorskip("PIL.Image")
    Image.new("RGB", (8, 8), "red").save(path)


def _perms(ws: Path) -> KitToolPermissions:
    return KitToolPermissions(workspace=ws, mode="bypassPermissions")


@pytest.fixture(autouse=True)
def _clear_pending():
    A._doc_executor.__dict__.pop("_pending_view_images", None)
    yield
    A._doc_executor.__dict__.pop("_pending_view_images", None)


def test_view_image_loads_validates_and_stashes(tmp_path):
    _png(tmp_path / "x.png")
    ex = A._doc_executor
    out = json.loads(ex._execute_view_image({"path": "x.png"}, _perms(tmp_path)))
    assert out["status"] == "ok" and out["mime"] == "image/png"
    # The image is stashed so the stream loop can inject it as visual content.
    assert ex._pending_view_images
    assert ex._pending_view_images[-1].source_path.name == "x.png"


def test_view_image_rejects_non_image(tmp_path):
    (tmp_path / "a.txt").write_text("hello")
    out = json.loads(
        A._doc_executor._execute_view_image({"path": "a.txt"}, _perms(tmp_path)))
    assert "error" in out


def test_view_image_missing_file(tmp_path):
    out = json.loads(
        A._doc_executor._execute_view_image({"path": "nope.png"}, _perms(tmp_path)))
    assert "error" in out


def test_read_file_redirects_images_to_view_image(tmp_path):
    _png(tmp_path / "y.png")
    out = json.loads(
        A._doc_executor._execute_read_file({"path": "y.png"}, _perms(tmp_path)))
    assert "error" in out and "view_image" in out["error"]


def test_view_image_wired_end_to_end():
    src = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
           / "api_client.py").read_text(encoding="utf-8")
    assert '"name": "view_image"' in src          # advertised tool
    assert 'name == "view_image"' in src          # dispatch route
    assert "_pending_view_images" in src          # loop injection
    assert "image_url" in src                       # injected as visual content
