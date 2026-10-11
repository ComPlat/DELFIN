"""show_molecule: the read gate, the result, and the chat card.

Input: a path from the model. Output: ``DELFIN_CARD:`` + JSON with the
resolved path, a title and a plain summary, or ``{"error": ...}``.

Probed at the layers that decide: the executor's dispatch (the read gate
the tool really calls), the dashboard function that turns a result into
chat HTML, and the two places that act on the marker.

Measured against the first version: the card was an
``<iframe sandbox srcdoc>`` that loaded an unpinned 3Dmol build from
3Dmol.org -- the sandboxed document cannot see the dashboard's bundled copy,
so on a node without network the card was an empty 460 px box -- showed
only the first frame of a trajectory, drew a cube's positive lobe without
atoms, and copied the whole file into the tool result. The card now is the
chat's existing 3D card (bundled 3Dmol), drawn from the resolved path.
"""
from __future__ import annotations

import ast
import inspect
import json
import pathlib

import pytest

from delfin.agent import api_client
from delfin.agent.api_client import _DocToolExecutor, KitToolPermissions
from delfin.dashboard import chat_viewer
from delfin.dashboard import tab_agent as T

H_CUBE = (
    "test cube\n"
    "generated\n"
    "1   0.0   0.0   0.0\n"
    "2   1.0   0.0   0.0\n"
    "2   0.0   1.0   0.0\n"
    "2   0.0   0.0   1.0\n"
    "1  0.0   0.0   0.0   0.0\n"
    "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8\n"
)
WATER = "3\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n"
_MARKER = chat_viewer.CARD_MARKER


def _perms(ws):
    return KitToolPermissions(workspace=str(ws), mode="default")


def _show(ws, path, **extra):
    return _DocToolExecutor()._dispatch(
        "show_molecule", {"path": str(path), **extra}, _perms(ws))


def _payload(out):
    assert out.startswith(_MARKER), out[:200]
    return json.loads(out[len(_MARKER):])


# -- the read gate -------------------------------------------------------


@pytest.fixture
def ws(tmp_path):
    root = tmp_path / "ws"
    root.mkdir()
    out = tmp_path / "out"
    out.mkdir()
    (out / "secret.xyz").write_text(WATER, encoding="utf-8")
    (root / "link.xyz").symlink_to(out / "secret.xyz")
    (root / "ok.xyz").write_text(WATER, encoding="utf-8")
    return root


@pytest.mark.parametrize("path", [
    "link.xyz",                 # symlink to a file outside
    "../out/secret.xyz",        # parent traversal
    "{out}/secret.xyz",         # absolute path outside
    "/etc/hostname",
])
def test_a_file_outside_the_roots_is_refused(ws, path):
    path = path.format(out=ws.parent / "out")
    out = _show(ws, path)
    assert not out.startswith(_MARKER)
    assert "read denied" in json.loads(out)["error"]


def test_a_secret_glob_is_refused(ws):
    (ws / ".env.xyz").write_text(WATER, encoding="utf-8")
    assert "deny-glob" in json.loads(_show(ws, ".env.xyz"))["error"]


def test_a_relative_path_is_the_workspaces(ws):
    p = _payload(_show(ws, "ok.xyz"))
    assert p["path"] == str((ws / "ok.xyz").resolve())
    assert p["title"] == "ok.xyz"


def test_a_relative_path_without_explicit_perms_uses_the_sessions(ws):
    ex = _DocToolExecutor()
    ex._permissions = _perms(ws)
    out = ex._execute_show_molecule({"path": "ok.xyz"}, None)
    assert out.startswith(_MARKER), out


def test_a_suffix_the_card_cannot_draw_is_refused(ws):
    (ws / "m.pdb").write_text("ATOM      1  O   HOH A   1\nEND\n",
                              encoding="utf-8")
    assert "error" in json.loads(_show(ws, "m.pdb"))


# -- the result ------------------------------------------------------------


def test_the_result_carries_no_markup_and_no_file_content(ws):
    (ws / "evil.xyz").write_text(
        "1\n</script><img src=x onerror=alert(1)> SENTINEL\nH 0 0 0\n",
        encoding="utf-8")
    out = _show(ws, "evil.xyz")
    p = _payload(out)
    assert set(p) == {"card", "path", "title", "isovalue", "text"}
    assert "SENTINEL" not in out and "<" not in p["text"]


@pytest.mark.parametrize("given, kept", [
    (0.05, 0.05), (-0.05, 0.05), (0, None), ("0.1", None), (True, None),
    (float("nan"), None),
])
def test_the_isovalue_is_a_positive_finite_number_or_the_default(
        ws, given, kept):
    (ws / "h.cube").write_text(H_CUBE, encoding="utf-8")
    assert _payload(_show(ws, "h.cube", isovalue=given))["isovalue"] == kept


def test_a_payload_of_another_shape_is_no_card():
    for bad in ("DELFIN_CARD:{", 'DELFIN_CARD:{"card":"molecule"}',
                'DELFIN_CARD:{"escape":"html","html":"<iframe>","text":"t"}',
                '{"error": "x"}'):
        assert chat_viewer.card_payload(bad) is None


# -- the chat card -----------------------------------------------------------


def test_the_card_is_the_chats_3d_card(ws):
    (ws / "traj.xyz").write_text(WATER * 3, encoding="utf-8")
    html = T._show_molecule_card_html(_show(ws, "traj.xyz"))
    assert "delfin-chat-mol3d" in html
    assert 'data-mol3d-frames="3"' in html
    assert "<iframe" not in html and "3Dmol.org" not in html


def test_the_isovalue_reaches_the_card(ws):
    (ws / "h.cube").write_text(H_CUBE, encoding="utf-8")
    html = T._show_molecule_card_html(_show(ws, "h.cube", isovalue=0.05))
    assert 'data-mol3d-iso="0.05"' in html and "±0.05" in html


def test_the_page_script_draws_the_cards_isovalue():
    src = pathlib.Path(inspect.getfile(T)).read_text(encoding="utf-8")
    assert "getAttribute('data-mol3d-iso')" in src
    assert "isoval: iso," in src and "isoval: -iso," in src


def test_a_hostile_title_is_escaped(ws):
    name = '"><img src=x onerror=alert(1)>.xyz'
    (ws / name).write_text(WATER, encoding="utf-8")
    html = T._show_molecule_card_html(_show(ws, name))
    assert "<img src=x" not in html
    assert html.count("<img") == 1          # the card's own starter


def test_with_viewers_off_the_card_says_so(ws, monkeypatch):
    from delfin.dashboard import molecule_viewer
    monkeypatch.setattr(molecule_viewer, "get_viewer_profile",
                        lambda: {"enabled": False})
    html = T._show_molecule_card_html(_show(ws, "ok.xyz"))
    assert "delfin-chat-mol3d" not in html
    assert "no 3D view" in html and "H2 O" in html


def test_an_error_is_no_card(ws):
    assert T._show_molecule_card_html(_show(ws, "link.xyz")) is None


# -- only show_molecule's result is read as a card ---------------------------


def _fn_src(module, name):
    tree = ast.parse(pathlib.Path(inspect.getfile(module)).read_text(
        encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.unparse(node)
    raise AssertionError(name)


def test_the_chat_draws_a_card_only_for_show_molecule():
    src = _fn_src(T, "_on_tool_result")
    assert "_show_molecule_card_html(tool_output)" in src
    assert ("if tool_name == 'show_molecule' and _raw_name in "
            "('show_molecule', 'mcp__delfin-docs__show_molecule'):"
            in src), "an MCP tool named show_molecule, or any output with " \
                     "the marker, must not be drawn as a card"
    assert "parse_card_result" not in src


def test_the_loop_reports_show_molecule_under_the_name_the_chat_accepts():
    """The OpenAI loop emits its own tools as mcp__<ns>__<name>; the chat
    must accept the namespace show_molecule actually lands in."""
    import re
    src = inspect.getsource(api_client)
    m = re.search(r'ns_prefix = "kit-coding" if is_coding else "delfin-docs"',
                  src)
    assert m, "the namespace rule moved; re-check the chat's accepted names"
    coding = src[:m.start()].rsplit("is_coding = fn_name in (", 1)[1]
    assert '"show_molecule"' not in coding.split(")")[0]


def test_the_ui_preview_of_a_large_file_is_still_a_card(ws):
    """The dashboard sees the first 2000 characters of a result. The first
    version's result carried the whole file twice-escaped: a 12-atom
    benzene was 2070 characters, cut mid-JSON, and never drawn."""
    atoms = "".join(f"C {i:.6f} 1.234567 -0.987654\n" for i in range(500))
    (ws / "big.xyz").write_text(f"500\nbig\n{atoms}", encoding="utf-8")
    out = _show(ws, "big.xyz")
    assert len(out) < 2000
    assert T._show_molecule_card_html(out[:2000]) is not None


def test_the_model_split_is_keyed_on_the_tool_name():
    import re
    src = " ".join(inspect.getsource(api_client).split())
    assert re.search(
        r'if \(fn_name == "show_molecule" and str\(result\)\.startswith\('
        r'"DELFIN_CARD:"\)\):', src)
