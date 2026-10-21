"""The repo map names the session's workspace, in any language.

Input: a PromptLoader whose ``workspace_root`` is a foreign repository
(TypeScript, Go, Rust) and whose ``repo_root`` is the prompt pack's tree.
Output: the ``--- Repo Map ---`` text.

Before the fix the map was built from ``repo_root`` -- the DELFIN source
tree -- so every foreign session was shown ``delfin/agent/...`` files,
and the index read only ``delfin/`` and ``tests/`` with a Python-only
symbol reader, so a TS/Go/Rust workspace produced nothing structural.
"""
from __future__ import annotations

from pathlib import Path

from delfin.agent import repo_map
from delfin.agent.prompt_loader import PromptLoader


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _foreign_repo(root: Path) -> Path:
    _write(root / "package.json", '{"name": "x", "scripts": {"test": "vitest"}}')
    _write(root / "src" / "math.ts",
           "export function addNumbers(a: number, b: number) {\n"
           "  return a + b;\n}\n"
           "export class Calculator {}\n")
    _write(root / "src" / "math.test.ts", "import { addNumbers } from './math';\n")
    _write(root / "app" / "main.go",
           "package main\n\nfunc HandleRequest() {}\n\ntype Server struct{}\n")
    _write(root / "crates" / "core" / "src" / "lib.rs",
           "pub fn parse_header(x: &str) -> u8 { 0 }\npub struct Frame;\n")
    _write(root / "node_modules" / "dep" / "index.js",
           "function addNumbers() {}\n")
    _write(root / ".env", "TOKEN=x\n")
    repo_map._build_repo_index.cache_clear()
    return root


def _loader_for(ws: Path) -> PromptLoader:
    loader = PromptLoader()
    loader.workspace_root = ws
    loader.is_delfin_workspace = False
    return loader


def test_the_map_is_built_from_the_workspace_not_the_pack(tmp_path):
    ws = _foreign_repo(tmp_path / "proj")
    ctx = _loader_for(ws)._load_repo_map_context(
        "Fix the failing test; the handler returns the wrong value")
    assert ctx, "a code task in a foreign workspace produced no map"
    assert "delfin/" not in ctx, ctx
    assert "app/" in ctx and "src/" in ctx, ctx


def test_the_layout_lists_top_level_dirs_and_sources_by_extension(tmp_path):
    ws = _foreign_repo(tmp_path / "proj")
    ctx = _loader_for(ws)._load_repo_map_context("fix the bug")
    assert ".ts" in ctx and ".go" in ctx and ".rs" in ctx, ctx
    assert "package.json" in ctx, ctx


def test_a_typescript_symbol_ranks_its_file(tmp_path):
    ws = _foreign_repo(tmp_path / "proj")
    ctx = _loader_for(ws)._load_repo_map_context(
        "fix addNumbers, it overflows")
    assert "src/math.ts" in ctx, ctx
    assert "addNumbers" in ctx, ctx
    # The test beside it is named as its test.
    assert "src/math.test.ts" in ctx, ctx


def test_go_and_rust_symbols_are_read(tmp_path):
    ws = _foreign_repo(tmp_path / "proj")
    ctx = _loader_for(ws)._load_repo_map_context(
        "fix HandleRequest and parse_header")
    assert "app/main.go" in ctx and "HandleRequest" in ctx, ctx
    assert "crates/core/src/lib.rs" in ctx and "parse_header" in ctx, ctx


def test_vendored_and_hidden_files_are_not_mapped(tmp_path):
    ws = _foreign_repo(tmp_path / "proj")
    ctx = _loader_for(ws)._load_repo_map_context("fix addNumbers")
    assert "node_modules" not in ctx, ctx
    assert ".env" not in ctx, ctx


def test_a_non_code_task_still_gets_no_map(tmp_path):
    ws = _foreign_repo(tmp_path / "proj")
    assert _loader_for(ws)._load_repo_map_context("hello there") == ""
