"""Adversarial tests for package V1: chat_media.molecule.

Against the builder's phase-2 API (delfin/agent/chat_media.py, hash af249946):
`molecule(path_or_xyz: str, *, kind=None, isovalue=None, workspace_root=None)`
returns a dict with keys html/text/kind/viewer_id/frames/atom_count/formula/path.

These tests probe the invariants the existing dashboard pipeline (tab_agent.py
role=tool verbatim HTML) requires: file content must not break out of the HTML
block, paths outside the workspace must be refused, size is capped, and the
terminal fallback is plain text. A red case here is a finding.
"""
from pathlib import Path

import pytest


WATER_XYZ = """3
water
O  0.0 0.0 0.1173
H  0.0 0.7572 -0.4692
H  0.0 -0.7572 -0.4692
"""

MULTI_XYZ = """2
frame1
H  0 0 0
H  0 0 1
2
frame2
H  0 0 0
H  1 0 0
"""

# A valid 3-atom XYZ whose comment line is hostile markup. role=tool renders
# tool results as verbatim HTML (tab_agent.py:10505), and data is embedded into
# a <script> block via json.dumps, which does NOT escape '</script>'.
MALICIOUS_XYZ_SCRIPT = """3
</script><script>window.__v1pwned=1</script>
C 0 0 0
O 1 0 0
H 0 1 0
"""

# Attribute-breakout attempt in a comment.
MALICIOUS_XYZ_ATTR = """3
" onmouseover="window.__v1pwned=1
C 0 0 0
O 1 0 0
H 0 1 0
"""

# Minimal but valid Gaussian cube (1 atom, 1 box per axis), single voxel.
CUBE_ONE = (
    "cube comment\n"
    "source comment\n"
    " 1  0.0 0.0 0.0\n"
    " 1  1.0 0.0 0.0\n"
    " 1  0.0 1.0 0.0\n"
    " 1  0.0 0.0 1.0\n"
    " 1  0.0 0.0 0.0\n"
    "1\n"
)

# A cube whose comment carries an HTML/JS breakout.
MALICIOUS_CUBE = (
    'x" onload="alert(1)"></script><script>v=1</script>\n'
    "src\n"
    " 1  0 0 0\n"
    " 1  1 0 0\n"
    " 1  0 1 0\n"
    " 1  0 0 1\n"
    " 1  0 0 0\n"
    "1\n"
)


def _write(tmp: Path, name: str, text: str) -> str:
    p = tmp / name
    p.write_text(text, encoding="utf-8")
    return str(p)


def _script_breakout(html: str) -> bool:
    """True if hostile content added a second '</script>' — a real breakout.

    The viewer block itself always ends with exactly ONE literal '</script>'
    (the close of its own <script> element, chat_media.py _viewer_html). A
    hostile file whose '../script>' survived unescaped would add a second one
    (or more), terminating the element early. So >1 means a genuine breakout.
    """
    return html.count("</script>") > 1


class TestInjection:
    """Data is embedded into a <script> block via json.dumps (chat_media.py
    _xyz_body/_cube_body). json.dumps does NOT escape '</script>', so hostile
    content must be neutralised or the file can execute JS when the block is
    inserted as verbatim chat HTML."""

    def test_malicious_xyz_script_breakout(self, tmp_path):
        from delfin.agent import chat_media
        p = _write(tmp_path, "bad.xyz", MALICIOUS_XYZ_SCRIPT)
        html = chat_media.molecule(p, kind="xyz")["html"]
        assert not _script_breakout(html), "hostile </script> reached the HTML"

    def test_malicious_cube_script_breakout(self, tmp_path):
        from delfin.agent import chat_media
        p = _write(tmp_path, "bad.cube", MALICIOUS_CUBE)
        html = chat_media.molecule(p, kind="cube")["html"]
        assert not _script_breakout(html), "hostile </script> reached the HTML"

    def test_attribute_breakout_from_comment(self, tmp_path):
        from delfin.agent import chat_media
        p = _write(tmp_path, "attr.xyz", MALICIOUS_XYZ_ATTR)
        html = chat_media.molecule(p, kind="xyz")["html"]
        # an unescaped double-quote + onmouseover would escape the data-*
        # attribute / JS string context.
        assert 'onmouseover="window' not in html


class TestReadGate:
    """A path is refused when workspace_root is given and the path is outside
    it (mirror read_file's _check_read_access at api_client.py:14769)."""

    def test_outside_workspace_root_refused(self, tmp_path):
        from delfin.agent import chat_media
        p = _write(tmp_path, "w.xyz", WATER_XYZ)
        other = tmp_path / "other"
        other.mkdir()
        with pytest.raises(ValueError):
            chat_media.molecule(p, kind="xyz", workspace_root=str(other))

    def test_path_inside_workspace_root_ok(self, tmp_path):
        from delfin.agent import chat_media
        p = _write(tmp_path, "w.xyz", WATER_XYZ)
        r = chat_media.molecule(p, kind="xyz", workspace_root=str(tmp_path))
        assert r["path"] == p
        assert r["atom_count"] == 3


class TestSizeCap:
    """A file over _MAX_FILE_BYTES is refused, not embedded (a huge payload
    would travel inside the transcript)."""

    def test_huge_file_refused(self, tmp_path):
        from delfin.agent import chat_media
        p = _write(tmp_path, "big.xyz", "C 0 0 0\n" * 600_000)
        with pytest.raises(ValueError):
            chat_media.molecule(p, kind="xyz", workspace_root=str(tmp_path))


class TestTerminalFallback:
    """The terminal gets plain TEXT (formula, atom count, path), not HTML."""

    def test_text_field_is_plain(self, tmp_path):
        from delfin.agent import chat_media
        p = _write(tmp_path, "w.xyz", WATER_XYZ)
        text = chat_media.molecule(p, kind="xyz", workspace_root=str(tmp_path))["text"]
        assert "Molecule:" in text
        assert "Atoms: 3" in text
        assert "<" not in text  # no HTML in the terminal text


class TestKindValidation:
    def test_bad_kind_rejected(self, tmp_path):
        from delfin.agent import chat_media
        p = _write(tmp_path, "w.xyz", WATER_XYZ)
        with pytest.raises(ValueError):
            chat_media.molecule(p, kind="nope", workspace_root=str(tmp_path))

    def test_none_path(self):
        from delfin.agent import chat_media
        with pytest.raises((ValueError, TypeError)):
            chat_media.molecule(None, kind="xyz")


class TestPathLike:
    """The function is documented as taking str, but the dashboard hands Path
    objects around (files_named_in_answer). Failing on a Path is a real
    friction the caller must catch."""

    def test_accepts_pathlike(self, tmp_path):
        from delfin.agent import chat_media
        p = Path(_write(tmp_path, "w.xyz", WATER_XYZ))
        r = chat_media.molecule(p, kind="xyz", workspace_root=str(tmp_path))
        assert r["atom_count"] == 3


class TestChatViewer:
    """phase-3a dashboard wrapper (delfin/dashboard/chat_viewer.py) reuses
    chat_media. It inherits whatever chat_media does — so to approve it the
    phase-2 invariants must hold; these pin that inheritance."""

    def test_terminal_text_has_no_html(self, tmp_path):
        from delfin.dashboard import chat_viewer
        p = _write(tmp_path, "w.xyz", WATER_XYZ)
        out = chat_viewer.render_molecule(
            p, kind="xyz", workspace_root=str(tmp_path), as_html=False)
        assert "<" not in out
        assert "Molecule:" in out

    def test_viewer_block_inherits_xss_breakout(self, tmp_path):
        # If chat_media still leaks `</script>`, the dashboard card that embeds
        # media["html"] verbatim (chat_viewer.py:36) carries the same XSS.
        from delfin.dashboard import chat_viewer
        p = _write(tmp_path, "bad.xyz", MALICIOUS_XYZ_SCRIPT)
        html = chat_viewer.render_molecule(
            p, kind="xyz", workspace_root=str(tmp_path), as_html=True)
        assert not _script_breakout(html), \
            "chat_viewer card inherits the chat_media `</script>` breakout"

    def test_viewer_accepts_pathlike(self, tmp_path):
        from delfin.dashboard import chat_viewer
        p = Path(_write(tmp_path, "w.xyz", WATER_XYZ))
        out = chat_viewer.render_molecule(
            p, kind="xyz", workspace_root=str(tmp_path), as_html=False)
        assert "Atoms: 3" in out
