"""Red control for PR #321 CodeQL py/clear-text-storage-sensitive-data.

The facts-file write used _m.get("steer_text", _m.get("text", "")): when
steer_text was absent the run's OWN output text (the model's answer,
which can carry anything) flowed into .longrun_facts.json — clear text
in the workspace. The patch (operator builds .gate/kg4_codeql_steer_
text.patch) writes ONLY the scripted steer text from the task YAML.

Both directions (operator checklist): fires when the run output would
leak (pre-patch), does not fire when only the scripted steer text is
present (post-patch).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


class _SpyEngine:
    """Duck-typed engine whose OUTPUT text carries a secret marker."""

    def __init__(self):
        self.session_id = "sess-codeql"
        self.token_usage = {"input": 0, "output": 0, "cached": 0}
        self.cost_usd = 0.0
        self.mode = "solo"
        self._model_name = "fake-model"

    @property
    def model(self):
        return self._model_name

    def stream_response(self, **kwargs):
        on_token = kwargs.get("on_token")
        if on_token:
            on_token("SECRET-RUN-OUTPUT-9f3c")
        return "SECRET-RUN-OUTPUT-9f3c"

    def estimate_context_tokens(self):
        return 0

    def delegate_spend(self):
        from delfin.agent import agent_metrics
        return agent_metrics.DelegateLedger()


def test_run_output_text_never_reaches_the_facts_file(tmp_path, monkeypatch):
    """The steer is DECLARED (so the facts file is written) but the run's
    output text must not appear in it — only the scripted steer text.

    The facts file exists only INSIDE the pristine-workspace guard (it
    is written before verify runs and removed with the guard), so the
    observation point is a spy on run_acceptance — the same call the
    real _run_task_once makes there (proven shape in
    test_kg4_operator_findings.py)."""
    # Resolved per call, not at import: the greenproof loads patched
    # copies into delfin.agent.benchmark_runner and this test must run
    # against whatever module sits there (pre-patch = red, patched =
    # green).
    br = sys.modules["delfin.agent.benchmark_runner"]
    task = br.Task(id="lr_steered_goal", task_class="office",
                mode="solo", prompt="old goal", verify="never_called.py",
                mid_run_messages=(("after",
                                   "PERSON (mid-task steering, new goal): "
                                   "follow the new goal"),))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "tests" / "fixtures" / "office_workspace").mkdir(
        parents=True, exist_ok=True)
    captured: dict = {}

    def spy_acceptance(name, workspace, **kw):
        p = Path(workspace) / ".longrun_facts.json"
        captured["exists"] = p.is_file()
        if p.is_file():
            captured["facts"] = json.loads(p.read_text(encoding="utf-8"))
        return None, "spy: not a real accept"

    monkeypatch.setattr(br, "run_acceptance", spy_acceptance)

    def _run_once(engine, prompt, **kw):
        captured.setdefault("prompts", []).append(prompt)
        return {"text": "SECRET-RUN-OUTPUT-9f3c", "error": "",
                "tool_calls": [{"name": "write_file"}],
                "input_tokens": 1, "output_tokens": 1}

    br._run_task_once(
        task, model="fake-model", backend="test", provider="test",
        profile_name="test", engine_factory=lambda *a, **kw: _SpyEngine(),
        max_tokens=256, run_once=_run_once, clock=lambda: 0.0)
    facts = captured.get("facts")
    assert facts is not None, captured
    texts = json.dumps(facts)
    assert "SECRET-RUN-OUTPUT-9f3c" not in texts, texts
    # The scripted steer text IS there (that is what the accept judges).
    assert any("new goal" in m.get("text", "")
               for m in facts["steer_delivered"]), facts
