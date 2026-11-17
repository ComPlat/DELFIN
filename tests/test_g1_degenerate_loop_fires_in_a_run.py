"""The degenerate coding loop becomes visible to the loop guards.

The failure mode: a model edits a line, reruns the same test, fails the same
way, twenty times. All four existing guards are blind to it -- a failing
``run_tests`` is a SUCCESSFUL tool call returning ``{"summary":
{"failed": 1}, "failures": [...]}`` (not an ``{"error"...}``), and the model's
edit differs every round (so the identical ``(name, args)`` check resets
every round). This file drives the real headless engine through the scripted
fake_model endpoint and asserts the fifth guard (test_loop_guard wiring in
api_client) stops the loop: steer once at 3 identical failures, abort at 5
with the ``degenerate_test_loop`` stop reason (exit code 10).

The wiring lives in api_client.py, which is a protected file; it is applied
by scratch/apply_api_client.py in an integration worktree. These tests pass
ONLY after that script has run. On a tree without the wiring the guard does
not fire and the run ends some other way -- list them in PR.md as failing on
the base tree and green after the apply.
"""

from __future__ import annotations

from tests.agent_terminal_probes import fake_model


# The fixture source: a function that returns the wrong constant. The model
# "edits" it each round (a comment line carrying the attempt number), but the
# edit does not change what the failing test observes, so the run_tests
# failure message is byte-identical every round while the edit arguments
# differ every round.
def _source(i: int) -> str:
    return ("def answer():\n"
            "    # attempt %d\n"
            "    return 3\n" % i)


def _test_source() -> str:
    return ("from source import answer\n"
            "\n"
            "\n"
            "def test_answer_is_forty_two():\n"
            "    assert answer() == 42\n")


class _Launch:
    """Scripted model: holds per-run state and records the system harness
    rows the engine sends, so the steer can be asserted."""

    def __init__(self, plan_fn):
        plan_fn.state = {"n": 0}
        self._plan_fn = plan_fn
        self.seen_system: list[str] = []

    def plan(self, messages):
        for m in messages or []:
            if m.get("role") == "system":
                text = str(m.get("content") or "")
                if text and text not in self.seen_system:
                    self.seen_system.append(text)
        return self._plan_fn(messages)


def _plan_degenerate(_messages) -> dict:
    """The degenerate loop: alternate an edit (attempt 0,1,2,...) with a
    run_tests on the always-failing test. The edit content differs every
    round, so the no-progress (identical tool args) guard cannot trip; the
    test result is identical, so the same-failure guard must."""
    st = _plan_degenerate.state
    last = _messages[-1] if _messages else {}
    if last.get("role") != "tool":
        # first model turn after the prompt: edit, then run it.
        return {"tool": "edit_file", "args": {
            "path": "source.py", "content": _source(st["n"])}}
    st["n"] += 1
    prev = None
    for m in reversed(_messages[:-1]):
        if m.get("role") == "assistant" and m.get("tool_calls"):
            prev = (m["tool_calls"][0].get("function") or {}).get("name")
            break
    if prev == "run_tests":
        return {"tool": "edit_file", "args": {
            "path": "source.py", "content": _source(st["n"])}}
    return {"tool": "run_tests", "args": {"target": "tests/test_source.py"}}


def _plan_changing(_messages) -> dict:
    """Control: the failing node set CHANGES every round. The model alternates
    between two always-failing tests, so the signature differs each round and
    the same-failure counter resets -- the loop must NOT be aborted."""
    st = _plan_changing.state
    st["n"] += 1
    last = _messages[-1] if _messages else {}
    if last.get("role") != "tool":
        return {"tool": "edit_file", "args": {
            "path": "source.py", "content": _source(999)}}
    return {"tool": "run_tests", "args": {
        "target": ("tests/test_a.py" if st["n"] % 2 else "tests/test_b.py")}}


def _cli_rc_events(capsys) -> list:
    """The JSON result object the headless run wrote to stdout."""
    out = capsys.readouterr().out
    result = None
    for line in out.splitlines():
        line = line.strip()
        if not line or not line.startswith("{"):
            continue
        try:
            result = __import__("json").loads(line)
        except ValueError:
            continue
    return result


def _run_engine(monkeypatch, tmp_path, capsys, plan_fn):
    """Launch one real headless turn against fake_model's scripted endpoint.
    The workspace holds an always-failing test the model edits + re-runs."""
    import io
    from pathlib import Path

    from delfin.agent import cli as agent_cli
    from delfin.agent import api_client as _A
    from delfin.agent import mcp_client as M
    from delfin.agent import model_capabilities as mc
    from delfin.agent import provider_profile as PP

    home = tmp_path / "home"
    ws = home / "ws"
    ws.mkdir(parents=True)
    (ws / "source.py").write_text(_source(0))
    (ws / "tests").mkdir()
    (ws / "tests" / "test_source.py").write_text(_test_source())
    (ws / "tests" / "test_a.py").write_text("def test_a():\n    assert 1 == 2\n")
    (ws / "tests" / "test_b.py").write_text("def test_b():\n    assert 3 == 4\n")
    (ws / "conftest.py").write_text("")       # workspace root on sys.path

    monkeypatch.chdir(ws)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))

    class _EmptyReg:
        def discover_all(self):
            return []

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    cap = mc.ModelCapabilities(model=fake_model.MODEL, provider="ollama",
                               context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: cap)
    monkeypatch.setattr(M, "get_registry", lambda *a, **k: _EmptyReg())
    baseline = tmp_path / "learned_profiles.json"
    baseline.write_bytes(PP._DEFAULT_PATH.read_bytes())
    monkeypatch.setattr(PP, "_DEFAULT_PATH", baseline)
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda engine, root, **kw: "sid-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create",
                        lambda engine, args: "")

    stdin = io.StringIO("")
    stdin.isatty = lambda: False              # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stdin)

    launch = _Launch(plan_fn)
    monkeypatch.setattr(fake_model, "plan", launch.plan)
    srv, port = fake_model.serve(log_path=str(tmp_path / "model.log"))
    monkeypatch.setenv("OLLAMA_HOST", f"http://127.0.0.1:{port}")
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    # Round budget above the abort point (5 test runs ~ round 10) so the
    # wired guard aborts first; small enough that an unwired run ends quickly.
    monkeypatch.setattr(_A, "_resolve_max_tool_rounds", lambda *a, **k: 12)

    rc = agent_cli.main(["-p", "make the tests pass", "--cwd", str(ws),
                         "--provider", "ollama",
                         "--model", fake_model.MODEL,
                         "--permission-mode", "default",
                         "--output-format", "json"])
    srv.shutdown()
    result = _cli_rc_events(capsys)
    return rc, result, launch


def test_the_degenerate_loop_steers_then_aborts(monkeypatch, tmp_path, capsys):
    """Edit+rerun the same failing test: steer once at 3, abort at 5 with the
    degenerate_test_loop stop reason and its own exit code 10."""
    rc, result, launch = _run_engine(monkeypatch, tmp_path, capsys,
                                     _plan_degenerate)
    assert result is not None, "no result object from the headless run"
    assert result["stop_reason"] == "degenerate_test_loop"
    assert rc == 10
    # The steer (injected once at 3 consecutive identical failures) reached
    # the model as a system harness row.
    assert any("SAME failures" in t for t in launch.seen_system)


def test_a_changing_failure_set_runs_on(monkeypatch, tmp_path, capsys):
    """Control: the failing node set CHANGES each round -> it is NOT a same
    failure -> the loop runs on and is never aborted by the new guard."""
    rc, result, launch = _run_engine(monkeypatch, tmp_path, capsys,
                                     _plan_changing)
    assert result is not None, "no result object from the headless run"
    assert result["stop_reason"] != "degenerate_test_loop"
    assert not any("SAME failures" in t for t in launch.seen_system)
