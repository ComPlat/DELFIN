"""Package V1 phase 3 — chat_viewer.render_molecule(): dashboard HTML vs terminal text.

Tests the dashboard-side HTML card builder (the string the shared
``role="tool"`` chat hook inlines) and its plain-text terminal fallback. The
red control for this phase is that ``delfin.dashboard.chat_viewer`` did not
exist on the phase-2 commit: every case fails on ``import`` there.
"""
import json

import pytest

from delfin.dashboard import chat_viewer

H2O_XYZ = """3
water
O  0.000000  0.000000  0.000000
H  0.000000  0.000000  0.957200
H  0.957200  0.000000  0.000000
"""

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


class TestDashboardHTML:
    def test_render_molecule_dashboard_gets_html(self):
        html = chat_viewer.render_molecule(H2O_XYZ, kind="xyz")
        # The dashboard inline card: a 3Dmol viewer with the drag hint, and the
        # self-contained block from chat_media.
        assert "drag to rotate" in html
        assert "addModel" in html
        assert "3Dmol-min.js" in html
        assert "delfin-molecule" in html

    def test_cube_gives_isosurface_block(self):
        html = chat_viewer.render_molecule(H2_CUBE, kind="cube", isovalue=0.02)
        assert "addVolumetricData" in html
        assert "cube" in html
        assert "cube isosurface" in html
        assert "0.02" in html

    def test_title_surfaces_in_card(self):
        html = chat_viewer.render_molecule(H2O_XYZ, kind="xyz", title="my water")
        assert "my water" in html


class TestTerminalFallback:
    def test_terminal_gets_text_not_html(self):
        text = chat_viewer.render_molecule(H2O_XYZ, kind="xyz", as_html=False)
        assert "<" not in text
        assert "script" not in text
        assert "H2" in text and "O" in text  # formula
        assert "Atoms: 3" in text  # atom count

    def test_cube_terminal_fallback_is_text(self):
        text = chat_viewer.render_molecule(H2_CUBE, kind="cube", as_html=False)
        assert "<" not in text
        assert "H" in text and "1" in text


class TestGuardrails:
    def test_path_outside_workspace_refused(self, tmp_path):
        root = tmp_path / "ws"
        root.mkdir()
        outside = tmp_path / "elsewhere"
        outside.mkdir()
        file = outside / "molec.xyz"
        file.write_text(H2O_XYZ, encoding="utf-8")
        with pytest.raises(ValueError):
            chat_viewer.render_molecule(
                str(file), kind="xyz", workspace_root=str(root))

    def test_unknown_kind_refused(self):
        with pytest.raises(ValueError):
            chat_viewer.render_molecule(H2O_XYZ, kind="nope")


class TestSandboxedCard:
    """The card must be an isolated <iframe sandbox> whose srcdoc carries the
    3D viewer (operator security findings on the phase-3b patch)."""

    def test_card_is_a_sandboxed_iframe(self):
        html = chat_viewer.render_molecule(H2_CUBE, kind="cube")
        assert "<iframe" in html
        assert 'sandbox="allow-scripts"' in html
        # No allow-same-origin: the viewer must not reach the parent dashboard.
        assert "allow-same-origin" not in html

    def test_card_renders_via_srcdoc(self):
        html = chat_viewer.render_molecule(H2O_XYZ, kind="xyz")
        assert "srcdoc=" in html

    def test_content_cannot_break_out_of_srcdoc(self):
        # Quotes + </script>/</iframe> in the file must be neutralised inside
        # the srcdoc so they can't break the attribute or the iframe element.
        cube = (
            'cube "quoted" </script><script>x=1</script></iframe><iframe>\n'
            "generated\n"
            "1 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n"
            "1 0 0 0 0\n0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8\n"
        )
        html = chat_viewer.render_molecule(cube, kind="cube")
        # Only the card's own closing </iframe> may appear raw; the injected
        # one is escaped (entity-encoded) inside srcdoc.
        assert html.count("</iframe>") == 1
        # No raw executable <script> outside the (escaped) srcdoc.
        assert html.count("</script>") == 0
        # The data still reaches the viewer (escaped, so inert).
        assert "addVolumetricData" in html


class TestCardContract:
    """The shared V1+V2 tool-result contract (DELFIN_CARD: marker + JSON)."""

    def test_tool_result_marker_and_keys(self):
        out = chat_viewer.tool_result("<iframe>card</iframe>", "some text")
        assert out.startswith(chat_viewer.CARD_MARKER)
        payload = json.loads(out[len(chat_viewer.CARD_MARKER):])
        assert payload["escape"] == "html"
        assert payload["html"] == "<iframe>card</iframe>"
        assert payload["text"] == "some text"

    def test_render_molecule_tool_result_contract(self):
        out = chat_viewer.render_molecule_tool_result(H2_CUBE, kind="cube")
        assert out.startswith(chat_viewer.CARD_MARKER)
        payload = json.loads(out[len(chat_viewer.CARD_MARKER):])
        # The card is the sandboxed iframe; the terminal line is plain text.
        assert 'sandbox="allow-scripts"' in payload["html"]
        assert "allow-same-origin" not in payload["html"]
        assert "<" not in payload["text"]
        assert "addVolumetricData" in payload["html"]
