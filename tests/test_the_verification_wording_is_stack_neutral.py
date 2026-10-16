"""The prompt's verification text names no Python runner outside Python, and
asks for the output to be checked against the request.

Built through the real engine for a Python, a TypeScript, a Rust and a Go
repository and for the DELFIN tree, on both roster models (GLM keeps every
lazy module; DeepSeek runs the slim selection), with no network:

* "Test and lint runners (`pytest`, `ruff`, `mypy`) are auto-allowed"
  reached the Rust, Go and TS prompts; the sentence is now stack-neutral;
* the profile rule "pytest in background after every code change"
  contradicted the synchronous verification rule and reached a Rust
  build; it is gone;
* Verify-after-modify now asks for the run's output to be checked against
  each requirement, not only the exit code (science_summarises: 7 of 7
  failing attempts printed the wrong row and called the summary correct).

Paid for by deletion: solo_agent.md 10,778 -> 10,733 tokens, cacheable
head 11,553 -> 11,508 (tests/test_prompt_token_budget.py).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import engine as E

MODELS = ("kit.glm-5.3", "kit.deepseek-v4-flash")
CHECK = "against\neach requirement of the request, not only the exit code"


def _repo(base: Path, kind: str) -> Path:
    base.mkdir(parents=True)
    if kind == "python":
        (base / "pyproject.toml").write_text('[project]\nname = "pkg"\n')
    elif kind == "ts":
        (base / "package.json").write_text(
            '{"name": "app", "scripts": {"test": "vitest"}}')
        (base / "tsconfig.json").write_text("{}")
    elif kind == "rust":
        (base / "Cargo.toml").write_text('[package]\nname = "c"\n')
    elif kind == "go":
        (base / "go.mod").write_text("module m\n")
    elif kind == "delfin":
        (base / "delfin" / "agent").mkdir(parents=True)
        (base / "delfin" / "__init__.py").write_text("")
        (base / "delfin" / "agent" / "__init__.py").write_text("")
        (base / "pyproject.toml").write_text('[project]\nname = "delfin"\n')
    return base


def _prompt(ws: Path, model: str) -> str:
    eng = E.AgentEngine(repo_dir=ws, backend="api", provider="kit",
                        api_key="stub-key", model=model, mode="solo",
                        permission_mode="default")
    eng._load_mode("solo")
    task = "The test suite fails; find the cause and fix it, then run the tests."
    eng.messages = [{"role": "user", "content": task}]
    return eng._build_current_system_prompt(task_text=task)


def _flat(text: str) -> str:
    return " ".join(text.split())


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("kind", ["python", "ts", "rust", "go", "delfin"])
def test_the_built_prompt_asks_for_the_output_against_the_request(
        tmp_path, kind, model):
    text = _flat(_prompt(_repo(tmp_path / kind, kind), model))
    assert _flat(CHECK) in text
    assert "are auto-allowed: run them" not in text
    assert "pytest in background after every code change" not in text
    assert "Test, lint and type-check runs need no asking" in text


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("kind,runner", [
    ("ts", "npm test"), ("rust", "cargo test"), ("go", "go test")])
def test_a_foreign_stack_gets_its_runner_and_no_python_runner(
        tmp_path, kind, runner, model):
    text = _prompt(_repo(tmp_path / kind, kind), model)
    assert runner in text
    assert "`mypy`) are auto-allowed" not in text
    assert "running pytest on the affected module" not in text


def test_the_profile_rule_is_gone_from_every_profile():
    import json
    profiles = json.loads((Path(E.__file__).parent / "learned_profiles.json")
                          .read_text())
    assert "pytest in background after every code change" not in \
        json.dumps(profiles)
