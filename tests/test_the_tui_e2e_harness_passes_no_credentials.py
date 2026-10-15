"""The window's live checks and the service tests start DELFIN without credentials.

Input: the environment builders the harnesses start DELFIN from --
tui/spike/check_tmux.py (scrubbed_env, and the tmux server every
tui/test/e2e check types its commands into) and tests/conftest.py
(child_env, the service and worker tests' children).
Output: none of them hands a credential variable on.

Why: DELFIN copies the provider keys it finds in its environment into
$HOME/.delfin/credentials.json. Under a check's temporary HOME that was a
copy of the developer's real key in /tmp after every run. A check that
needs a key sets a fake one itself (scrubbed_env(**extra)).
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
E2E = ROOT / "tui" / "test" / "e2e"
sys.path.insert(0, str(ROOT / "tui" / "spike"))

import check_tmux  # noqa: E402
from conftest import child_env, is_credential  # noqa: E402

FAKE = "fake-value-that-must-not-pass"
CREDENTIALS = ("KIT_TOOLBOX_API_KEY", "OPENAI_API_KEY", "ANTHROPIC_API_KEY", "OPENAI_BASE_URL",
               "KIT_TOOLBOX_URL", "GH_TOKEN", "GITHUB_TOKEN", "HF_TOKEN", "MY_SERVICE_TOKEN",
               "DEEPSEEK_API_KEY", "SOME_API_KEY", "AWS_SECRET_ACCESS_KEY", "PGPASSWORD",
               "GOOGLE_APPLICATION_CREDENTIALS", "anthropic_api_key")
ORDINARY = ("PATH", "HOME", "LANG", "TERM", "COLORTERM", "PYTHONPATH", "DELFIN_RUNTIME_DIR",
            "DELFIN_SCRATCH_STATE", "DELFIN_TUI_LOGO", "SSH_AUTH_SOCK", "USER", "TOKENIZERS_PARALLELISM")


@pytest.fixture()
def keys(monkeypatch):
    for name in CREDENTIALS:
        monkeypatch.setenv(name, FAKE)
    monkeypatch.setenv("LANG", "C.UTF-8")


def test_the_two_rules_agree_and_name_every_credential():
    for name in CREDENTIALS:
        assert check_tmux.is_credential(name) and is_credential(name), name
    for name in ORDINARY:
        assert not check_tmux.is_credential(name) and not is_credential(name), name
    for attr in ("CREDENTIAL_PARTS", "CREDENTIAL_PREFIXES", "CREDENTIAL_NAMES"):
        import conftest
        assert getattr(check_tmux, attr) == getattr(conftest, attr), attr


def test_scrubbed_env_drops_them_and_keeps_the_rest(keys):
    env = check_tmux.scrubbed_env()
    assert not [k for k in env if k in CREDENTIALS]
    assert FAKE not in env.values()
    assert env["LANG"] == "C.UTF-8" and "PATH" in env
    # A check that needs a key names a fake one.
    assert check_tmux.scrubbed_env(OPENAI_API_KEY="sk-test")["OPENAI_API_KEY"] == "sk-test"


def test_the_tmux_server_starts_without_them(keys, monkeypatch):
    seen: list[dict] = []

    class Done:
        stdout = ""

    def run(argv, **kw):
        seen.append(kw.get("env"))
        return Done()

    monkeypatch.setattr(check_tmux.subprocess, "run", run)
    tmux = check_tmux.Tmux(80, 24)
    tmux.screen()
    tmux.close()
    assert len(seen) == 3
    for env in seen:
        assert env is not None, "a tmux call inherits the caller's environment"
        assert FAKE not in env.values()


def _calls(path: Path):
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "subprocess" \
                and node.func.attr in ("run", "Popen", "check_call", "check_output", "call"):
            yield node


def test_every_e2e_check_starts_its_children_from_a_scrubbed_environment():
    checks = sorted(E2E.glob("check_*.py"))
    assert len(checks) >= 10
    for path in checks:
        text = path.read_text(encoding="utf-8")
        for spelling in ("**os.environ", "dict(os.environ", "os.environ.copy(", "os.environ.items("):
            assert spelling not in text, f"{path.name}: {spelling} copies the caller's credentials"
        for call in _calls(path):
            assert any(k.arg == "env" for k in call.keywords), \
                f"{path.name}:{call.lineno}: a child without env= inherits the caller's credentials"


def test_child_env_drops_them(keys, tmp_path):
    env = child_env(tmp_path)
    assert not [k for k in env if k in CREDENTIALS]
    assert FAKE not in env.values()
    assert env["HOME"] == str(tmp_path / "child_home")
