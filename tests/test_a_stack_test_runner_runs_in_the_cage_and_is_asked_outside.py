"""A stack's test or check runner runs unasked inside the process cage, and
is asked everywhere else; its install, publish and network verbs are never
admitted.

Diagnosis 2026-10-15: 0 of 12 cargo test, 0 of 8 go test and 0 of 10 npm
test runs executed in 1,018 benchmark attempts -- each was refused with
"... Then STOP and wait", while pytest, which runs repository code the same
way, is on the auto-allow list. The rule (caged_test_runners, consulted by
``_caged_command_may_run_unasked``) admits the runner of the stack whose
marker file is in the run's directory, and only where bubblewrap, the
process cage, Landlock and a non-open egress policy hold. Outside the cage
the refusal no longer says STOP: it asks for the exact command once and
names the change UNTESTED until it ran.

The cage is pinned by monkeypatch, as in
test_a_caged_command_runs_without_asking.py, so the tables hold on every
host; the ones that need a real runner are not here.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent import caged_test_runners as ctr
from delfin.agent.api_client import KitToolPermissions, _doc_executor


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "home").mkdir()
    d = tmp_path / "ws"
    (d / "rs" / "src").mkdir(parents=True)
    (d / "rs" / "Cargo.toml").write_text('[package]\nname = "rs"\n')
    (d / "gomod").mkdir()
    (d / "gomod" / "go.mod").write_text("module gomod\n")
    (d / "ts" / "node_modules" / ".bin").mkdir(parents=True)
    (d / "ts" / "package.json").write_text(json.dumps(
        {"name": "ts", "scripts": {"test": "node --test test/"}}))
    for tool in ("vitest", "tsc"):
        (d / "ts" / "node_modules" / ".bin" / tool).write_text("")
    (d / "bare").mkdir()
    (d / "bare" / "package.json").write_text(json.dumps({"name": "bare"}))
    return d


def _cage(monkeypatch, kind: str) -> None:
    monkeypatch.setattr(A, "_bash_isolation_decision",
                        lambda perms, mode=None: (kind, "auto"))
    monkeypatch.setattr(A, "_process_cage_enabled", lambda: kind == "bwrap")
    monkeypatch.setattr(A, "_landlock_functional", lambda: kind == "bwrap")


def _gate(ws: Path, cmd: str, cwd: str, mode: str = "default",
          attended: bool = False):
    asked: list[str] = []

    def _cb(tool, args, preview):
        asked.append(preview)
        return False

    perms = KitToolPermissions(mode=mode, workspace=str(ws),
                               confirm_callback=_cb if attended else None)
    admitted: list[bool] = []
    args = {"command": cmd}
    if cwd:
        args["cwd"] = str(ws / cwd)
    verdict = _doc_executor._run_permission_gate(
        "bash", args, perms, caged_admission=admitted)
    return verdict, asked, admitted


# (command, cwd relative to the workspace)
RUNNERS = [
    ("cargo test", "rs"),
    ("cargo test --offline", "rs"),
    ("cargo check", "rs"),
    ("cargo check --offline --all-targets", "rs"),
    ("cargo clippy --offline -- -D warnings", "rs"),
    ("cargo test parse_roman -- --nocapture", "rs"),
    ("cargo test --manifest-path rs/Cargo.toml", ""),
    ("cd rs && cargo test --offline", ""),
    ("cargo test 2>&1 | tail -40", "rs"),
    ("go test ./...", "gomod"),
    ("go vet ./...", "gomod"),
    ("go test -run TestParsePort -count=1 -v ./...", "gomod"),
    ("go test ./... 2>&1 | head -50", "gomod"),
    ("cd gomod && go test ./...", ""),
    ("npm test", "ts"),
    ("npm run test", "ts"),
    ("npm test -- --test-reporter=spec", "ts"),
    ("npx vitest run", "ts"),
    ("npx tsc --noEmit", "ts"),
    ("npx tsc --noEmit -p tsconfig.json", "ts"),
    ("node --test test/", "ts"),
]

# Never admitted by the cage rule, wherever the cage holds.
NEVER = [
    ("cargo install ripgrep", "rs"),
    ("cargo publish", "rs"),
    ("cargo add serde", "rs"),
    ("cargo update", "rs"),
    ("cargo run", "rs"),
    ("cargo +nightly test", "rs"),
    ("cargo test --config 'target.x.runner=\"sh\"'", "rs"),
    ("cargo test -Z unstable-options", "rs"),
    ("cargo test && curl https://example.com", "rs"),
    ("cargo test > out.txt", "rs"),
    ("RUSTC_WRAPPER=sh cargo test", "rs"),
    ("cargo test", "bare"),
    ("go get example.com/x", "gomod"),
    ("go install example.com/x@latest", "gomod"),
    ("go mod download", "gomod"),
    ("go run .", "gomod"),
    ("go test -exec sh ./...", "gomod"),
    ("go test -toolexec sh ./...", "gomod"),
    ("go test -o bin ./...", "gomod"),
    ("go test $(echo ./...)", "gomod"),
    ("GOFLAGS=-mod=mod go test ./...", "gomod"),
    ("go test ./...", "bare"),
    ("npm install", "ts"),
    ("npm ci", "ts"),
    ("npm publish", "ts"),
    ("npm run build", "ts"),
    ("npm exec vitest", "ts"),
    ("npm test", "bare"),
    ("npm test; rm -rf src", "ts"),
    ("npm test &", "ts"),
    ("npx cowsay hi", "ts"),
    ("npx -y vitest run", "ts"),
    ("npx --package vitest vitest run", "ts"),
    ("npx vitest run", "bare"),
    ("npx tsc", "ts"),
    ("cd .. && cargo test", "rs"),
]


@pytest.mark.parametrize("attended", [True, False])
@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
@pytest.mark.parametrize("cmd,cwd", RUNNERS)
def test_a_stack_runner_runs_unasked_in_the_cage(cmd, cwd, mode, attended,
                                                 ws, monkeypatch):
    _cage(monkeypatch, "bwrap")
    verdict, asked, admitted = _gate(ws, cmd, cwd, mode, attended)
    assert verdict is None, verdict
    assert asked == []
    assert admitted == [True]


@pytest.mark.parametrize("cmd,cwd", [r for r in RUNNERS
                                     if not r[0].startswith("node ")])
def test_without_a_cage_the_runner_is_asked(cmd, cwd, ws, monkeypatch):
    _cage(monkeypatch, "none")
    verdict, asked, admitted = _gate(ws, cmd, cwd, attended=True)
    assert admitted == []
    assert len(asked) == 1 and verdict is not None


@pytest.mark.parametrize("cmd,cwd", [r for r in RUNNERS
                                     if not r[0].startswith("node ")])
def test_without_a_cage_the_headless_refusal_asks_once_and_names_untested(
        cmd, cwd, ws, monkeypatch):
    _cage(monkeypatch, "none")
    verdict, asked, admitted = _gate(ws, cmd, cwd, attended=False)
    assert admitted == [] and asked == []
    assert "not on the auto-allow list" in verdict
    assert "UNTESTED" in verdict
    assert "STOP and wait" not in verdict


@pytest.mark.parametrize("attended", [True, False])
@pytest.mark.parametrize("cmd,cwd", NEVER)
def test_an_install_publish_or_network_variant_is_never_admitted(
        cmd, cwd, attended, ws, monkeypatch):
    _cage(monkeypatch, "none")
    before = _gate(ws, cmd, cwd, attended=attended)
    _cage(monkeypatch, "bwrap")
    after = _gate(ws, cmd, cwd, attended=attended)
    assert after[2] == []
    assert (after[0] is None) == (before[0] is None)
    assert len(after[1]) == len(before[1])


@pytest.mark.parametrize("cmd,cwd", NEVER)
def test_the_grammar_admits_none_of_the_never_table(cmd, cwd, ws):
    # A runner chained to something else is admitted as a segment only;
    # the rest is what the cage rule's exclusions then judge.
    rest = ctr.without_admitted_runners(cmd, str(ws / cwd), str(ws))
    assert rest is None or A._caged_command_names_excluded(rest), rest


def test_a_runner_admission_is_a_channel(ws, monkeypatch):
    from delfin.agent import channels
    _cage(monkeypatch, "bwrap")
    channels.reset()
    _gate(ws, "cargo test --offline", "rs")
    assert channels.snapshot().get("caged_test_runner") == 1
    channels.reset()


def test_the_cage_setting_switched_off_keeps_asking(ws, monkeypatch):
    _cage(monkeypatch, "bwrap")
    from delfin import user_settings
    monkeypatch.setattr(user_settings, "load_settings", lambda *a, **k: {
        "agent": {"caged_commands_run_unasked": False}})
    verdict, asked, admitted = _gate(ws, "cargo test", "rs", attended=True)
    assert admitted == [] and len(asked) == 1


def test_an_open_network_keeps_asking(ws, monkeypatch):
    _cage(monkeypatch, "bwrap")
    from delfin.agent import egress_proxy
    monkeypatch.setattr(egress_proxy, "network_mode", lambda: "open")
    verdict, asked, admitted = _gate(ws, "go test ./...", "gomod",
                                     attended=True)
    assert admitted == [] and len(asked) == 1


def test_another_refusal_still_says_stop(ws, monkeypatch):
    _cage(monkeypatch, "none")
    verdict, _, _ = _gate(ws, "cargo install ripgrep", "rs")
    assert "STOP and wait" in verdict
    assert "UNTESTED" not in verdict
