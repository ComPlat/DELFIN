"""The container image, the MCP example and the answering model are named exactly.

Input: ``delfin.agent.container_exec`` (``DEFAULT_OCI_IMAGE``,
``_SIF_SOURCE``), .github/ci/agent-container/Dockerfile, .github/dependabot.yml,
the container job of .github/workflows/ci.yml, docs/USER_MANUAL.md, the MCP
examples in ``delfin.agent.mcp_client`` and the dashboard's ``/mcp`` help,
and ``turn_metrics.record_request`` as the OpenAI-compatible client calls it.

Output: failures for
- a default image named by tag: ``python:3.11-slim`` is a pointer the
  registry can move, so the isolation backend would run whatever it points
  at that day (OWASP LLM03, CIS SSC 3.1.7). The reference must be
  ``repository@sha256:<64 hex>``, the same digest as the Dockerfile that
  Dependabot's ``docker`` ecosystem updates, the image CI pulls, and the
  commands the manual gives;
- an MCP example that runs ``npx -y`` on a package without a version:
  whoever copies it runs the newest upload of that package on every start;
- a request line in the metrics that does not carry the model id the
  provider says answered (``served_model``), when the response names one:
  a routing alias such as ``kit.deepseek-v4-flash`` can be pointed at
  another model without any change here, and the record is the only place
  that difference would show.
"""

from __future__ import annotations

import re
import threading
import types
from pathlib import Path

import pytest
import yaml

from delfin.agent import container_exec as ce
from delfin.agent import turn_metrics as tm

ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / ".github" / "ci" / "agent-container" / "Dockerfile"
_DIGEST = re.compile(r"@(sha256:[0-9a-f]{64})$")


def _digest(reference: str) -> str:
    match = _DIGEST.search(reference)
    assert match, f"{reference} is not pinned by digest"
    return match.group(1)


def test_the_default_image_is_a_digest():
    assert re.fullmatch(r"docker\.io/library/python@sha256:[0-9a-f]{64}",
                        ce.DEFAULT_OCI_IMAGE), ce.DEFAULT_OCI_IMAGE
    assert ce._SIF_SOURCE == "docker://" + ce.DEFAULT_OCI_IMAGE


def test_dependabot_keeps_the_digest_current_through_the_dockerfile():
    lines = [ln for ln in DOCKERFILE.read_text(encoding="utf-8").splitlines()
             if ln.startswith("FROM ")]
    assert len(lines) == 1, lines
    reference = lines[0].split()[1]
    assert reference.startswith("python:3.11-slim@"), reference
    assert _digest(reference) == _digest(ce.DEFAULT_OCI_IMAGE), (
        "Dependabot moved the Dockerfile; set container_exec.DEFAULT_OCI_IMAGE "
        "to the same digest")
    config = yaml.safe_load((ROOT / ".github" / "dependabot.yml").read_text(encoding="utf-8"))
    docker = [u for u in config["updates"] if u["package-ecosystem"] == "docker"]
    assert docker and docker[0]["directory"] == "/.github/ci/agent-container"
    assert docker[0].get("cooldown", {}).get("default-days", 0) >= 7


def test_ci_and_the_manual_pull_the_pinned_image():
    wf = yaml.safe_load((ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8"))
    pull = next(s["run"] for s in wf["jobs"]["container"]["steps"]
                if "Pull the default image" in s.get("name", ""))
    assert f"docker pull {ce.DEFAULT_OCI_IMAGE}" in pull
    assert f"podman pull {ce.DEFAULT_OCI_IMAGE}" in pull
    manual = (ROOT / "docs" / "USER_MANUAL.md").read_text(encoding="utf-8")
    assert "python:3.11-slim`" not in manual and "pull python:3.11-slim" not in manual
    assert ce.DEFAULT_OCI_IMAGE in manual
    assert ce._SIF_SOURCE in manual


@pytest.mark.parametrize("path", ["delfin/agent/mcp_client.py",
                                  "delfin/dashboard/tab_agent.py"])
def test_the_mcp_example_names_a_version(path):
    text = (ROOT / path).read_text(encoding="utf-8")
    found = re.findall(r"@modelcontextprotocol/server-filesystem(@[\w.-]+)?", text)
    assert found, f"{path}: the example is gone"
    assert all(found), f"{path}: npx -y without a version runs the newest upload"


def _client(session: str, served: str):
    from delfin.agent import api_client as ac

    text = types.SimpleNamespace(
        model=served, usage=None, choices=[types.SimpleNamespace(
            delta=types.SimpleNamespace(content="ok", tool_calls=None,
                                        reasoning_content=None),
            finish_reason="stop")])
    usage = types.SimpleNamespace(
        model=served, usage=types.SimpleNamespace(
            prompt_tokens=10, completion_tokens=1,
            prompt_tokens_details=types.SimpleNamespace(cached_tokens=0)),
        choices=[])

    class _Stream:
        def __iter__(self):
            return iter((text, usage))

        def close(self):
            pass

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    return _Stream()

    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = _Stub()
    c.model = "kit.deepseek-v4-flash"
    c._provider = "kit"
    c._base_url = "https://ki-toolbox.scc.kit.edu/api/v1"
    c._api_key = "x"
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._stop_flag = False
    c.metrics_session = session
    return c


def _one_request(client) -> None:
    try:
        for ev in client.stream_message(
                messages=[{"role": "user", "content": "hi"}],
                system="t", max_tokens=8):
            if getattr(ev, "type", "") == "message_start":
                break
    except Exception:
        pass


def test_the_request_line_names_the_model_that_answered(tmp_path, monkeypatch):
    monkeypatch.setattr(tm, "_DIR", tmp_path / "turn_metrics")
    _one_request(_client("served", "deepseek-v4-flash-2026-09-01"))
    (entry,) = tm.read_requests("served")
    assert entry["model"] == "kit.deepseek-v4-flash"
    assert entry["served_model"] == "deepseek-v4-flash-2026-09-01"


def test_a_response_without_a_model_adds_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(tm, "_DIR", tmp_path / "turn_metrics")
    _one_request(_client("anon", ""))
    (entry,) = tm.read_requests("anon")
    assert "served_model" not in entry


def test_the_served_model_is_bounded_text(tmp_path, monkeypatch):
    """The value comes from the endpoint: it is cut and stripped of control
    characters before it is written into a file a report quotes."""
    monkeypatch.setattr(tm, "_DIR", tmp_path / "turn_metrics")
    tm.record_request("bounded", model="m", served_model="x\n\x1b[2Jy" + "z" * 500)
    (entry,) = tm.read_requests("bounded")
    assert "\n" not in entry["served_model"] and "\x1b" not in entry["served_model"]
    assert len(entry["served_model"]) <= 120


def test_a_report_shows_which_model_answered(tmp_path, monkeypatch):
    monkeypatch.setattr(tm, "_DIR", tmp_path / "turn_metrics")
    tm.record_request("shown", model="kit.deepseek-v4-flash", ttft_ms=900,
                      served_model="deepseek-v4-flash-2026-09-01")
    tm.record_request("shown", model="same", ttft_ms=900, served_model="same")
    _, first, second = tm.format_requests(tm.read_requests("shown")).splitlines()
    assert "served=deepseek-v4-flash-2026-09-01" in first
    assert "served=" not in second
