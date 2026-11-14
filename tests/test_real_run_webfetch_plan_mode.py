"""Real run: one row of the admission matrix through the real entry point.

Input: ``python -P -m delfin.agent.cli chat --provider ollama --model
fake-model --new --permission-mode plan -p "do-webfetch <url>"`` in a
child process, against the scripted OpenAI-compatible endpoint
(tests/agent_terminal_probes/fake_model.py, its ``do-webfetch`` marker)
on 127.0.0.1. The endpoint answers the first round with one
``web_fetch`` call of an outside URL and every later round with text.

Output pinned: what the session sent back to the model for that call is
the plan-mode refusal -- the same decision the matrix records for
``native`` / ``plan_*`` / ``web_fetch_outside``
(tests/test_every_admission_path_decides_a_call_the_same_way.py) -- and
the run ends with the exit code of a refused turn
(turn_outcome.EXIT_CODES["refused"]).

The child gets a private HOME and no credentials (conftest.child_env),
and its interpreter runs with -P, so the delfin it imports is the one on
PYTHONPATH (this checkout), not one found from the working directory.
Nothing leaves the machine: the refusal comes before the URL is resolved.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from conftest import child_env  # noqa: E402
from delfin.agent.turn_outcome import EXIT_CODES  # noqa: E402
from tests.agent_terminal_probes import fake_model  # noqa: E402

URL = "https://sink.example/collect"


def test_real_session_refuses_webfetch_in_plan_mode(tmp_path, monkeypatch):
    seen: list = []
    original = fake_model.plan

    def recording_plan(messages):
        seen.append(json.loads(json.dumps(messages)))
        return original(messages)

    monkeypatch.setattr(fake_model, "plan", recording_plan)
    server, port = fake_model.serve(log_path=str(tmp_path / "model.log"))
    try:
        ws = tmp_path / "ws"
        ws.mkdir()
        env = child_env(tmp_path)
        env.update(PYTHONPATH=str(REPO),
                   OLLAMA_HOST=f"http://127.0.0.1:{port}")
        for var in ("NO_COLOR", "OLLAMA_BASE_URL", "DELFIN_LIFELINE_PID"):
            env.pop(var, None)
        done = subprocess.run(
            [sys.executable, "-P", "-m", "delfin.agent.cli", "chat",
             "--provider", "ollama", "--model", fake_model.MODEL, "--new",
             "--permission-mode", "plan", "-p", f"do-webfetch {URL}"],
            cwd=str(ws), env=env, capture_output=True, text=True,
            timeout=300, stdin=subprocess.DEVNULL)
    finally:
        server.shutdown()
        server.server_close()
    out = done.stdout + done.stderr

    assert seen, "the endpoint was never asked:\n" + out[-2000:]
    first_user = fake_model._text_of(fake_model._last(seen[0], "user"))
    assert first_user.startswith(f"do-webfetch {URL}"), first_user[:200]
    results = [m for msgs in seen[1:] for m in msgs
               if m.get("role") == "tool"]
    assert results, "no tool result reached the model:\n" + out[-2000:]
    body = fake_model._text_of(results[0])
    assert json.loads(body)["error"].startswith("plan mode (read-only)"), body
    assert "'web_fetch' is not on the plan-mode read-only tool list" in body
    assert done.returncode == EXIT_CODES["refused"], out[-2000:]
    # The child's interpreter, environment and working directory import
    # this checkout's delfin.
    where = subprocess.run(
        [sys.executable, "-P", "-c", "import delfin; print(delfin.__file__)"],
        cwd=str(ws), env=env, capture_output=True, text=True, timeout=60)
    assert Path(where.stdout.strip()).resolve() == (
        REPO / "delfin" / "__init__.py").resolve(), where.stdout + where.stderr
