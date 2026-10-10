"""Package V1 — red test for the protected `show_molecule` tool wiring.

Asserts the shared V1+V2 contract: a valid result is `DELFIN_CARD:` + JSON
{"escape":"html","html":"<sandboxed iframe card>","text":"<plain line>"},
the card is an isolated iframe (sandbox=allow-scripts, no allow-same-origin,
srcdoc escaped against quotes/</script>/</iframe>), and an ERROR output
carries no raw HTML (so a hostile path echoed in an error can never be
inlined or executed). tab_agent inlines ONLY the marker card and escapes
everything else.

The dispatch cases are RED until the operator builds
`.gate/v1_show_molecule.patch` (protected-file wiring); the terminal-text
case (committed chat_viewer) is green already. Run only through the gate.
"""
from __future__ import annotations

import json

import pytest

from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions
from delfin.dashboard import chat_viewer

# A tiny but parseable Gaussian cube: 1 hydrogen on a 2x2x2 grid.
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

_MARKER = chat_viewer.CARD_MARKER


def _perms(ws):
    return KitToolPermissions(workspace=str(ws), mode="default")


def _payload(out):
    """Parse a DELFIN_CARD: result; fail if the marker is absent."""
    assert out.startswith(_MARKER), f"no card marker in {out[:120]!r}"
    return json.loads(out[len(_MARKER):])


class TestShowMoleculeWiring:
    def test_tool_is_registered_and_dispatchable(self, tmp_path):
        cube = tmp_path / "mole.cube"
        cube.write_text(H2_CUBE, encoding="utf-8")
        ex = _DocToolExecutor()
        out = ex._dispatch("show_molecule", {"path": str(cube)}, _perms(tmp_path))
        # The handler must be reached (not "unknown tool") and return a card.
        assert "unknown tool" not in out.lower()
        assert out.startswith(_MARKER)

    def test_cube_card_is_a_sandboxed_iframe(self, tmp_path):
        cube = tmp_path / "mole.cube"
        cube.write_text(H2_CUBE, encoding="utf-8")
        ex = _DocToolExecutor()
        out = ex._dispatch(
            "show_molecule", {"path": str(cube), "isovalue": 0.02},
            _perms(tmp_path))
        html = _payload(out)["html"]
        assert "<iframe" in html
        assert 'sandbox="allow-scripts"' in html
        assert "allow-same-origin" not in html  # isolated from the dashboard

    def test_cube_card_has_isosurface_block(self, tmp_path):
        cube = tmp_path / "mole.cube"
        cube.write_text(H2_CUBE, encoding="utf-8")
        ex = _DocToolExecutor()
        out = ex._dispatch(
            "show_molecule", {"path": str(cube), "isovalue": 0.02},
            _perms(tmp_path))
        html = _payload(out)["html"]
        assert "addVolumetricData" in html
        assert "0.02" in html

    def test_content_cannot_break_out_of_srcdoc(self, tmp_path):
        cube = tmp_path / "mole.cube"
        cube.write_text(
            'cube "quoted" </script><script>x=1</script></iframe><iframe>\n'
            "generated\n1 0 0 0\n2 1 0 0\n2 0 1 0\n2 0 0 1\n"
            "1 0 0 0 0\n0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8\n",
            encoding="utf-8",
        )
        ex = _DocToolExecutor()
        out = ex._dispatch("show_molecule", {"path": str(cube)}, _perms(tmp_path))
        html = _payload(out)["html"]
        # Only the card's own closing </iframe> may be raw; the injected one
        # is escaped inside sircdoc.
        assert html.count("</iframe>") == 1
        assert html.count("</script>") == 0


class TestErrorsAreEscaped:
    def test_path_outside_workspace_is_refused(self, tmp_path):
        root = tmp_path / "ws"
        root.mkdir()
        outside = tmp_path / "elsewhere"
        outside.mkdir()
        cube = outside / "mole.cube"
        cube.write_text(H2_CUBE, encoding="utf-8")
        ex = _DocToolExecutor()
        out = ex._dispatch("show_molecule", {"path": str(cube)}, _perms(root))
        assert not out.startswith(_MARKER)  # an error, not a card
        assert "error" in out

    def test_html_path_in_error_is_never_inlined(self, tmp_path):
        # A hostile path echoed by an error must NOT come back as a card that
        # tab_agent inlines: the error output is escaped JSON, never HTML.
        evil = '<img src=x onerror="window.__pwned=1">'
        ex = _DocToolExecutor()
        out = ex._dispatch("show_molecule", {"path": evil}, _perms(tmp_path))
        assert not out.startswith(_MARKER)
        assert "<img" not in out  # escaped: no raw HTML executable


class TestTerminalFallback:
    def test_terminal_gets_text_not_html(self, tmp_path):
        cube = tmp_path / "mole.cube"
        cube.write_text(H2_CUBE, encoding="utf-8")
        ex = _DocToolExecutor()
        out = ex._dispatch("show_molecule", {"path": str(cube)}, _perms(tmp_path))
        text = _payload(out)["text"]
        assert "<" not in text
        assert "script" not in text
        assert "H" in text  # formula


@pytest.mark.parametrize(
    "kind, expected",
    [("xyz", "XYZ"), ("cube", "cube"), ("pdb", "pdb")],
)
def test_kind_is_validated(kind, expected):
    ex = _DocToolExecutor()
    out = ex._dispatch("show_molecule", {"path": "x." + kind}, _perms("ws"))
    assert not out.startswith(_MARKER)
    assert expected in out.lower()
