"""The role prompt carries DELFIN/chemistry blocks only where they apply,
and its verification advice matches the workspace's stack.

Measured on a Python, a TypeScript, a Rust and a Go repository (solo
role, GLM, every lazy module kept): about 1,544 tokens of the role prompt
were DELFIN/chemistry blocks (ORCA decision tree, data-search and ops-MCP
sections, the two-contexts preamble, dashboard access, the DELFIN tail of
"After every code edit"), 379 more in the provider profile, and about 770
were Python-only advice that misdirects on TS/Rust/Go. Built through the
real engine, with no network.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import engine as E

GLM = "kit.glm-5.3"           # every lazy module kept (all_prompt_modules)

DELFIN_ONLY = (
    "### Quick decision tree",
    "## DELFIN ops MCP tools",
    "## Data search tools",
    "## ORCA / chemistry questions",
    "## Two work contexts",
    "## Dashboard access",
    "backend_slurm.py",
    "--- Chemistry Protocol ---",
)
# Safety text that must survive in every workspace.
SAFETY = (
    "## Permission boundary",
    "## Sandbox security boundary",
    "## Directory permissions",
    "## Confirm before mutating",
    "## Never fabricate tool results",
)


def _delfin(base: Path) -> Path:
    (base / "delfin" / "agent").mkdir(parents=True)
    (base / "delfin" / "__init__.py").write_text("")
    (base / "delfin" / "agent" / "__init__.py").write_text("")
    (base / "pyproject.toml").write_text('[project]\nname = "delfin"\n')
    return base


def _ts(base: Path) -> Path:
    base.mkdir(parents=True)
    (base / "package.json").write_text(
        '{"name": "app", "scripts": {"test": "vitest", "build": "tsc"}}')
    (base / "tsconfig.json").write_text("{}")
    return base


def _rust(base: Path) -> Path:
    base.mkdir(parents=True)
    (base / "Cargo.toml").write_text('[package]\nname = "crate"\n')
    return base


def _python(base: Path) -> Path:
    base.mkdir(parents=True)
    (base / "pyproject.toml").write_text('[project]\nname = "pkg"\n')
    return base


def _engine(ws: Path, model: str = GLM) -> E.AgentEngine:
    eng = E.AgentEngine(repo_dir=ws, backend="api", provider="kit",
                        api_key="stub-key", model=model, mode="solo",
                        permission_mode="default")
    eng._load_mode("solo")
    return eng


def _prompt(eng: E.AgentEngine, task: str) -> str:
    eng.messages = [{"role": "user", "content": task}]
    return eng._build_current_system_prompt(task_text=task)


def _sections(eng: E.AgentEngine, task: str) -> dict[str, str]:
    eng.messages = [{"role": "user", "content": task}]
    captured = {}
    orig = eng.loader.compose_sections

    def _spy(**kw):
        secs = orig(**kw)
        captured.update({s.name: s.content for s in secs})
        return secs

    eng.loader.compose_sections = _spy
    try:
        eng._build_current_system_prompt(task_text=task)
    finally:
        eng.loader.compose_sections = orig
    return captured


TASK = "The test suite fails; find the cause and fix it, then run the tests."


@pytest.mark.parametrize("make", [_ts, _rust, _python])
def test_a_foreign_repo_gets_no_delfin_block(tmp_path, make):
    text = _prompt(_engine(make(tmp_path / "repo")),
                   "Refactor this method into a pure function and update "
                   "the functional component; compare the energy method.")
    leaked = [m for m in DELFIN_ONLY if m in text]
    assert leaked == [], leaked
    missing = [m for m in SAFETY if m not in text]
    assert missing == [], missing


def test_the_delfin_tree_keeps_every_block(tmp_path):
    text = _prompt(_engine(_delfin(tmp_path / "d")), TASK)
    for marker in DELFIN_ONLY[:-1]:
        assert marker in text, marker
    for marker in SAFETY:
        assert marker in text, marker


def test_the_profile_drops_delfin_rules_on_a_foreign_repo(tmp_path):
    foreign = _sections(_engine(_ts(tmp_path / "app")), TASK).get(
        "profile", "")
    home = _sections(_engine(_delfin(tmp_path / "d")), TASK).get(
        "profile", "")
    assert "DELFIN" in home or "search_docs" in home, (
        "the shipped profile carries no DELFIN rule; the test sees nothing")
    for word in ("DELFIN", "search_docs", "backend_slurm", "ORCA"):
        assert word not in foreign, word


def test_verification_advice_follows_the_stack(tmp_path):
    ts = _prompt(_engine(_ts(tmp_path / "app")), TASK)
    rust = _prompt(_engine(_rust(tmp_path / "crate")), TASK)
    py = _prompt(_engine(_python(tmp_path / "pkg")), TASK)
    home = _prompt(_engine(_delfin(tmp_path / "d")), TASK)
    for text in (ts, rust):
        assert "ast.parse" not in text
        assert ".venv-<projekt>" not in text
    assert "tsc --noEmit" in ts
    assert "cargo test" in rust and "tsc --noEmit" not in rust
    assert "ast.parse" in py and "tsc --noEmit" not in py
    assert "ast.parse" in home
