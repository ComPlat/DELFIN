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


class TestParseCardResult:
    """The reader half of the shared contract: only a real card is inlined."""

    def test_roundtrip_extracts_the_card(self):
        out = chat_viewer.render_molecule_tool_result(H2_CUBE, kind="cube")
        card = chat_viewer.parse_card_result(out)
        assert card is not None
        assert "<iframe" in card
        assert 'sandbox="allow-scripts"' in card

    def test_non_marker_returns_none(self):
        # An error (which may echo a hostile path) is never a card.
        assert chat_viewer.parse_card_result('{"error": "boom"}') is None
        evil = '<img src=x onerror="window.__pwned=1">'
        err = json.dumps({"error": evil})
        assert chat_viewer.parse_card_result(err) is None

    def test_malformed_or_wrong_shape_returns_none(self):
        assert chat_viewer.parse_card_result("DELFIN_CARD:not-json") is None
        assert chat_viewer.parse_card_result("DELFIN_CARD:{}") is None
        assert chat_viewer.parse_card_result(
            'DELFIN_CARD:{"escape":"text","html":"<img>"}') is None
        assert chat_viewer.parse_card_result(
            'DELFIN_CARD:{"escape":"html","text":"x"}') is None
    def test_nonce_fence_guard(self):
        # A wrapped result is not a card: the nonce fence means the input no
        # longer starts with the bare DELFIN_CARD: marker, so it must NOT be
        # inlined — parse_card_result returns None.
        wrapped = "BEGIN external>\nDELFIN_CARD:{}"
        assert chat_viewer.parse_card_result(wrapped) is None


class TestStrictCardSecurity:
    """Operator Finding-2 attack cases: parse_card_result must REJECT any html
    that is not the exact single-sandboxed-srcdoc-iframe card, even when the
    marker and {escape:"html", html:<str>} shape are present. Full pads of the
    operator's attack set — these were red on the loose parser and are green
    here. Builds a syntactically-valid card JSON whose html is attacker-shaped
    but wrapped in the DELFIN_CARD: marker."""

    def _card(self, html: str) -> str:
        return chat_viewer.CARD_MARKER + json.dumps(
            {"escape": "html", "html": html, "text": "x"})

    def test_selective_html_with_onload_and_src(self):
        # onload/src iframe: the escaped viewer may carry src= inside srcdoc,
        # but a REAL src= or on*= attribute on the iframe itself must reject.
        evil = ('<iframe sandbox="allow-scripts" src="//evil" '
                'srcdoc="<script>alert(1)</script>" onload="pwn()"></iframe>')
        assert chat_viewer.parse_card_result(self._card(evil)) is None

    def test_allow_same_origin_card_refused(self):
        evil = '<iframe sandbox="allow-scripts allow-same-origin" srcdoc="a"></iframe>'
        assert chat_viewer.parse_card_result(self._card(evil)) is None

    def test_non_sandboxed_card_refused(self):
        evil = '<iframe srcdoc="<script>alert(1)</script>"></iframe>'
        assert chat_viewer.parse_card_result(self._card(evil)) is None

    def test_extra_content_besides_iframe_refused(self):
        # A sibling <script> outside the escaped srcdoc is an injection vector.
        evil = '<iframe sandbox="allow-scripts" srcdoc="a"></iframe><script>alert(1)</script>'
        assert chat_viewer.parse_card_result(self._card(evil)) is None

    def test_two_iframes_refused(self):
        evil = ('<iframe sandbox="allow-scripts" srcdoc="a"></iframe>'
                '<iframe sandbox="allow-scripts" srcdoc="b"></iframe>')
        assert chat_viewer.parse_card_result(self._card(evil)) is None
