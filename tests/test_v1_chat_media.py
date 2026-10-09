"""Package V1 — chat_media.molecule(): molecules/XYZ/cube shown in the chat.

Phase 2 tests for the pure HTML/text builder in ``delfin/agent/chat_media.py``.
Every case here is exercised against the REAL module through its public
``molecule(path_or_xyz, *, kind, workspace_root=None)`` entry point. Fixtures
are embedded so nothing on disk outside the workspace is read.

These are the red control for phase 2: on the baseline commit (no
``delfin/agent/chat_media.py``) every case fails on ``import``.
"""
import json

import pytest

import delfin.agent.chat_media as chat_media

H2O_XYZ = """3
water
O  0.000000  0.000000  0.000000
H  0.000000  0.000000  0.957200
H  0.957200  0.000000  0.000000
"""

MULTI_XYZ = (
    "3\nframe 1\n"
    "O  0.000000  0.000000  0.000000\n"
    "H  0.000000  0.000000  0.957200\n"
    "H  0.957200  0.000000  0.000000\n"
    "\n"
    "3\nframe 2\n"
    "O  0.050000  0.010000  0.000000\n"
    "H  0.010000  0.120000  0.957200\n"
    "H  0.957200  0.020000  0.100000\n"
)

# Small but genuinely parseable Gaussian cube: 1 hydrogen atom on a 2x2x2 grid.
H2_CUBE = (
    "test cube\n"
    "generated for v1 tests\n"
    "1   0.0   0.0   0.0\n"
    "2   1.0   0.0   0.0\n"
    "2   0.0   1.0   0.0\n"
    "2   0.0   0.0   1.0\n"
    "1  0.0   0.0   0.0   0.0\n"
    "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8\n"
)

_CDN = "3Dmol.org/build/3Dmol-min.js"


class TestMoleculeXYZ:
    def test_xyz_returns_html_block(self):
        res = chat_media.molecule(H2O_XYZ, kind="xyz")
        html = res["html"]
        # Self-contained, CDN 3Dmol, loads the model as xyz.
        assert "3Dmol-min.js" in html or _CDN in html
        assert "addModel" in html
        assert "width:100%" in html
        # A unique viewer id (several molecules in one chat must not clash).
        assert res["viewer_id"]
        assert res["kind"] == "xyz"

    def test_xyz_text_fallback_formula_and_count(self):
        res = chat_media.molecule(H2O_XYZ, kind="xyz")
        text = res["text"]
        assert "H2" in text and "O" in text  # formula
        assert "3" in text  # atom count
        # The terminal must get text, not HTML.
        assert "<" not in text and "script" not in text


class TestMoleculeMultiXyz:
    def test_multifxyz_html_and_frame_count(self):
        res = chat_media.molecule(MULTI_XYZ, kind="multixyz")
        assert res["kind"] == "multixyz"
        assert res["frames"] == 2
        assert "addModel" in res["html"]

    def test_multifxyz_text_fallback_mentions_frames(self):
        res = chat_media.molecule(MULTI_XYZ, kind="multixyz")
        assert "2" in res["text"]


class TestMoleculeCube:
    def test_cube_isosurface_block(self):
        res = chat_media.molecule(H2_CUBE, kind="cube", isovalue=0.02)
        assert res["kind"] == "cube"
        # The isosurface is drawn via 3Dmol's volumetric cube API.
        assert "addVolumetricData" in res["html"]
        assert "cube" in res["html"]
        assert "0.02" in res["html"] or "0.02" in json.dumps(res["html"])

    def test_cube_text_fallback(self):
        res = chat_media.molecule(H2_CUBE, kind="cube")
        text = res["text"]
        assert "H" in text and "1" in text  # one atom, hydrogen


class TestMoleculeGuardrails:
    def test_path_outside_workspace_refused(self, tmp_path):
        # The workspace_root does NOT contain the file -> the read gate must
        # refuse it even though it exists on disk.
        root = tmp_path / "ws"
        root.mkdir()
        outside = tmp_path / "elsewhere"
        outside.mkdir()
        file = outside / "molec.xyz"
        file.write_text(H2O_XYZ, encoding="utf-8")
        with pytest.raises(ValueError):
            chat_media.molecule(str(file), kind="xyz", workspace_root=str(root))

    def test_file_size_capped(self, tmp_path, monkeypatch):
        # Cap the module limit so a tiny file is already "too large" -> the
        # oversize guard trips without a megabytes-hard disk write.
        monkeypatch.setattr(chat_media, "_MAX_FILE_BYTES", 10)
        big = tmp_path / "big.xyz"
        big.write_text("1\nbig\n" + "H 0 0 0\n" * 100, encoding="utf-8")
        with pytest.raises(ValueError):
            chat_media.molecule(str(big), kind="xyz", workspace_root=str(tmp_path))

    def test_unknown_kind_refused(self):
        with pytest.raises(ValueError):
            chat_media.molecule(H2O_XYZ, kind="nope")

    def test_xyz_path_read_from_workspace(self, tmp_path):
        good = tmp_path / "water.xyz"
        good.write_text(H2O_XYZ, encoding="utf-8")
        res = chat_media.molecule(str(good), kind="xyz", workspace_root=tmp_path)
        assert res["kind"] == "xyz"
        assert str(good) in res["text"]  # path surfaces in the fallback


# ---------------------------------------------------------------------------
# Adversarial + PathLike regression tests (reviewer nacht-s32 findings).
# ---------------------------------------------------------------------------


class TestInjectionGuards:
    """A file whose content tries to break out of the surrounding <script>."""

    MALICIOUS = '3\n</script><script>window.__v1pwned=1</script>\nO 0 0 0\nH 0 0 1\nH 0 1 0\n'

    def test_xyz_comment_cannot_break_out_of_script(self):
        res = chat_media.molecule(self.MALICIOUS, kind="xyz")
        html = res["html"]
        # Only the one closing tag _viewer_html itself emits may appear; the
        # injected "</script>" must be escaped (e.g. as <\\/script>), so the
        # malicious payload cannot terminate the script element and become
        # executable markup.
        assert html.count("</script>") == 1
        # After JS-escape the payload is inert data inside a string literal,
        # never a live element.
        assert "<\\/script>" in html

    def test_cube_comment_cannot_break_out_of_script(self):
        cube = (
            "cube </script><script>window.__v1pwned=1</script>\n"
            "generated\n"
            "1 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n"
            "1 0 0 0 0\n0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8\n"
        )
        res = chat_media.molecule(cube, kind="cube")
        assert res["html"].count("</script>") == 1


class TestPathLike:
    def test_accepts_pathlib_path(self, tmp_path):
        import pathlib
        good = tmp_path / "water.xyz"
        good.write_text(H2O_XYZ, encoding="utf-8")
        # A caller (the dashboard) hands a Path, not a str — must not raise.
        res = chat_media.molecule(
            pathlib.Path(good), kind="xyz", workspace_root=tmp_path)
        assert res["kind"] == "xyz"
        assert "H2" in res["text"]

    def test_accepts_os_pathlike(self, tmp_path):
        import os
        good = tmp_path / "water.xyz"
        good.write_text(H2O_XYZ, encoding="utf-8")
        res = chat_media.molecule(
            os.fspath(good), kind="xyz", workspace_root=tmp_path)
        assert res["kind"] == "xyz"
