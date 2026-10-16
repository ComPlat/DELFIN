"""Text from outside DELFIN reaches the dashboard's HTML as text, not markup.

Two kinds of sink:

* HTML built with f-strings (``widgets.HTML(value=f"...")``). A file
  name, a job directory, an exception message or a cluster command's
  stderr is interpolated into markup; without ``html.escape`` a ``<`` in
  it starts an element. Pinned by an AST lint over ``delfin/dashboard``
  for the values that always come from outside -- exception variables,
  ``.stdout``/``.stderr``, a path's ``.name``/``.stem`` -- and strictly
  (every interpolation) for the three tabs the gap analysis named.
* Data embedded in an inline ``<script>`` (a 3D viewer's XYZ/cube text).
  ``json.dumps`` leaves ``</script>`` intact, which ends the script element
  early. ``helpers.js_string_literal`` writes ``<``, ``>`` and ``&`` as
  JavaScript unicode escapes; py3Dmol's own embedding is routed through
  the same rule.

The behavioural checks use a harmless marker (``</script><b>marker</b>``)
and look only at the generated string.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

DASHBOARD = Path(__file__).resolve().parents[1] / "delfin" / "dashboard"

#: Calls whose result is safe to interpolate into HTML.
ESCAPERS = {
    "escape", "_esc", "esc", "_escape", "_html_escape", "html_escape",
    "quote", "_e", "_h", "e", "_attr", "_calc_display_path",
    "js_string_literal", "_cell",
}

#: Files whose HTML f-strings are checked interpolation by interpolation.
STRICT = ("tab_turbomole_builder.py", "tab_application.py",
          "tab_pipeline_status.py")


def _escaped(expr: ast.AST) -> bool:
    for node in ast.walk(expr):
        if isinstance(node, ast.Call):
            func = node.func
            name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
            if name in ESCAPERS:
                return True
    return False


def _html_fstrings(tree: ast.AST):
    """(function node, f-string node) for f-strings that contain markup."""
    for fn in ast.walk(tree):
        if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Module)):
            continue
        for node in ast.walk(fn):
            if not isinstance(node, ast.JoinedStr):
                continue
            text = "".join(v.value for v in node.values
                           if isinstance(v, ast.Constant) and isinstance(v.value, str))
            if "<" in text and ">" in text:
                yield fn, node


def _external(expr: ast.AST, exc_names: set) -> str:
    if isinstance(expr, ast.Attribute) and expr.attr == "__name__":
        return ""                     # a class name is code, not text
    for node in ast.walk(expr):
        if isinstance(node, ast.Name) and node.id in exc_names:
            return "exception text"
        if isinstance(node, ast.Attribute) and node.attr in ("stdout", "stderr"):
            return "process output"
        if isinstance(node, ast.Attribute) and node.attr in ("name", "stem"):
            return "a path's name"
    return ""


def _sites():
    out = []
    for path in sorted(DASHBOARD.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for fn, node in _html_fstrings(tree):
            exc_names = {h.name for h in ast.walk(fn)
                         if isinstance(h, ast.ExceptHandler) and h.name}
            for value in node.values:
                if not isinstance(value, ast.FormattedValue):
                    continue
                if value.format_spec is not None or _escaped(value.value):
                    continue
                why = _external(value.value, exc_names)
                if why:
                    out.append(f"{path.name}:{value.lineno} {why}: "
                               f"{ast.unparse(value.value)}")
    return sorted(set(out))


def test_no_html_fstring_in_the_dashboard_interpolates_outside_text_raw():
    sites = _sites()
    assert sites == [], "\n".join(sites)


#: Interpolations in the strict files that are DELFIN's own values: CSS
#: colours, counters, the dropdown's fixed backend names, run ids, and
#: fragments whose name ends in ``_html`` (built from escaped parts).
_OWN = {"color", "status_badge", "elapsed", "data_info",
        "len", "i", "int", "cores", "backend", "run_id"}


@pytest.mark.parametrize("name", STRICT)
def test_the_named_tabs_escape_every_interpolation(name):
    tree = ast.parse((DASHBOARD / name).read_text(encoding="utf-8"))
    raw = []
    for _fn, node in _html_fstrings(tree):
        for value in node.values:
            if not isinstance(value, ast.FormattedValue):
                continue
            expr = value.value
            if value.format_spec is not None or _escaped(expr):
                continue
            if isinstance(expr, ast.Constant):
                continue
            root = expr
            while isinstance(root, (ast.Call, ast.Attribute, ast.BinOp)):
                root = (root.func if isinstance(root, ast.Call)
                        else root.value if isinstance(root, ast.Attribute)
                        else root.left)
            if isinstance(root, ast.Name) and (
                    root.id in _OWN or root.id.endswith("_html")):
                continue
            raw.append(f"{name}:{value.lineno}: {ast.unparse(expr)}")
    raw = sorted(set(raw))
    assert raw == [], "\n".join(raw)


# -- data inside an inline <script> --------------------------------------------

MARKER = "</script><b>marker</b>"


def test_the_fukui_viewer_keeps_file_text_inside_its_script(monkeypatch):
    from delfin.dashboard import molecule_viewer as MV

    monkeypatch.setattr(MV, "get_viewer_profile", lambda: {
        "enabled": True, "width": 400, "height": 400, "zoom": 1.0,
        "style_js": "{}", "viewer_config_js": "{}"})
    xyz = f"1\n{MARKER}\nH 0 0 0\n"
    cube = f"{MARKER}\n"
    html = MV.build_fukui_viewer_html(xyz, labels=[MARKER], cube_text=cube)
    assert "<b>marker</b>" not in html
    assert html.count("</script>") == 1


def test_py3dmol_embeds_model_text_without_a_script_end(monkeypatch):
    import py3Dmol

    from delfin.dashboard import molecule_viewer  # noqa: F401 - installs the rule

    view = py3Dmol.view(width=100, height=100)
    view.addModel(f"1\n{MARKER}\nH 0 0 0\n", "xyz")
    html = view._make_html()
    assert "<b>marker</b>" not in html
