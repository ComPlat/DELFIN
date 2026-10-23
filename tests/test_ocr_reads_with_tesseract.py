"""OCR through pytesseract, the engine ocr_availability names first.

ocr_availability reported pytesseract as the engine wherever tesseract was
on PATH, and _ocr_image only knew easyocr: it refused, and every OCR read
of a scan came back empty -- also where easyocr was installed beside it.

Driven with a stand-in pytesseract, so it runs without tesseract.
"""

from __future__ import annotations

import io
import sys
import types

import pytest

from delfin.agent import office


def _png() -> bytes:
    from PIL import Image
    buf = io.BytesIO()
    Image.new("L", (8, 8), 255).save(buf, format="PNG")
    return buf.getvalue()


def _fake_tesseract(monkeypatch, *, languages=("deu", "eng", "osd")):
    calls: dict = {}
    mod = types.ModuleType("pytesseract")
    mod.Output = types.SimpleNamespace(DICT="dict")
    mod.get_languages = lambda config="": list(languages)

    def image_to_data(image, lang="", output_type=None):
        calls["lang"] = lang
        return {"text": ["", "Beleg", "4711", "", "zweite", "Zeile"],
                "conf": ["-1", "96", "90", "-1", "40", "30"],
                "block_num": [1, 1, 1, 2, 2, 2],
                "par_num": [1, 1, 1, 1, 1, 1],
                "line_num": [0, 1, 1, 0, 1, 1]}

    mod.image_to_data = image_to_data
    monkeypatch.setitem(sys.modules, "pytesseract", mod)
    monkeypatch.setattr(office, "ocr_availability", lambda: {
        "available": True, "engine": "pytesseract", "detail": [], "next_step": ""})
    return calls


def test_tesseract_reads_lines_with_their_confidence(monkeypatch):
    calls = _fake_tesseract(monkeypatch)

    text, mean, weak = office._ocr_image(_png())

    assert text == "Beleg 4711\nzweite Zeile"
    assert calls["lang"] == "deu+eng"
    assert mean == pytest.approx((0.93 + 0.35) / 2)
    assert weak == 1


def test_a_missing_language_is_left_out(monkeypatch):
    calls = _fake_tesseract(monkeypatch, languages=("eng",))
    office._ocr_image(_png())
    assert calls["lang"] == "eng"


def test_no_installed_language_names_the_package(monkeypatch):
    _fake_tesseract(monkeypatch, languages=("osd",))
    with pytest.raises(office.OfficeError, match="tesseract-ocr-deu tesseract-ocr-eng"):
        office._ocr_image(_png())
