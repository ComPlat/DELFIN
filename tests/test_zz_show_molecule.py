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

    def test_hostile_path_never_inlined(self, tmp_path):
        # A hostile path echoed by an error must NOT become a card tab_agent
        # inlines. The dispatch returns escaped JSON; parse_card_result returns
        # None for it, so the chat hook escapes it like any other tool output
        # (never inlined, never executes).
        evil = '<img src=x onerror="window.__pwned=1">'
        ex = _DocToolExecutor()
        out = ex._dispatch("show_molecule", {"path": evil}, _perms(tmp_path))
        assert not out.startswith(_MARKER)  # never a card
        assert chat_viewer.parse_card_result(out) is None  # never inlined

    def test_molecule_in_non_repo_workspace_is_accepted(self, tmp_path):
        # Operator-fix regression: the containment root is the readable root the
        # file lies in (the session workspace), not the DELFIN checkout -- a
        # molecule in a user's own folder must be accepted, not refused.
        cube = tmp_path / "mole.cube"
        cube.write_text(H2_CUBE, encoding="utf-8")
        ex = _DocToolExecutor()
        out = ex._dispatch("show_molecule", {"path": str(cube)}, _perms(tmp_path))
        assert out.startswith(_MARKER)  # accepted and rendered as a card

    def test_relative_path_resolves_against_the_workspace(self, tmp_path):
        # Operator-fix regression: a relative path resolves against
        # perms.workspace, like every other file tool.
        cube = tmp_path / "mole.cube"
        cube.write_text(H2_CUBE, encoding="utf-8")
        ex = _DocToolExecutor()
        out = ex._dispatch("show_molecule", {"path": "mole.cube"}, _perms(tmp_path))
        assert out.startswith(_MARKER)


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
    "kind, content",
    [
        (
            "xyz",
            "2\nwater\nO  0.000000  0.000000  0.000000\n"
            "H  0.000000  0.000000  0.957200\n",
        ),
        ("cube", H2_CUBE),
    ],
)
def test_kind_is_validated(tmp_path, kind, content):
    # A real existing file of each SUPPORTED kind passes the read gate first,
    # so the dispatch validates the format rather than refusing on a path
    # error. (The committed version used workspace "ws" that does not exist:
    # the read gate rightly refused before kind-validation was reached.)
    f = tmp_path / ("x." + kind)
    f.write_text(content, encoding="utf-8")
    ex = _DocToolExecutor()
    out = ex._dispatch("show_molecule", {"path": str(f)}, _perms(tmp_path))
    assert out.startswith(_MARKER)  # each supported kind is recognised


def test_unsupported_kind_is_refused(tmp_path):
    # PDB is NOT a supported kind for show_molecule (the handler parses XYZ and
    # cube only); a .pdb file must be refused as a card, never inlined.
    pdb = tmp_path / "x.pdb"
    pdb.write_text(
        "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00"
        "  0.00           O\nEND\n",
        encoding="utf-8",
    )
    ex = _DocToolExecutor()
    out = ex._dispatch("show_molecule", {"path": str(pdb)}, _perms(tmp_path))
    assert not out.startswith(_MARKER)  # not a card
    assert "error" in out  # refused loudly
