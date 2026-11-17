"""Bring-your-own-licence tools stay outside the package (LICENSE_GUARD.md)."""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BYOL_TOOLS = ("export_csd_crystals",)


def test_no_module_under_delfin_imports_a_byol_tool():
    hits = []
    for p in (ROOT / "delfin").rglob("*.py"):
        try:
            tree = ast.parse(p.read_text(encoding="utf-8", errors="replace"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""] + [a.name for a in node.names]
            if any(t in n for n in names for t in BYOL_TOOLS):
                hits.append(str(p.relative_to(ROOT)))
    assert not hits, hits


def test_the_byol_tools_live_in_tools_only():
    for t in BYOL_TOOLS:
        assert (ROOT / "tools" / f"{t}.py").exists()
        assert not list((ROOT / "delfin").rglob(f"{t}.py"))
