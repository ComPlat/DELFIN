"""A client without a permissions policy reads, writes and fetches nothing.

Security review 2026-09-16: ``create_client(backend="api",
provider="openai")`` built the OpenAI client without a policy. It still
offered read_file, grep_file, list_files and remember to the model, and
with no policy read_file took any absolute path with no secret deny list.
The memory distiller built its cheap-tier client the same way and handed
it chat text that includes other sessions' messages and tool output.
"""
import json
from types import SimpleNamespace

import pytest

from delfin.agent import api_client as A


@pytest.fixture
def secret(tmp_path):
    ssh = tmp_path / ".ssh"
    ssh.mkdir()
    key = ssh / "id_ed25519"
    key.write_text("-----BEGIN OPENSSH PRIVATE KEY-----\nAAAA\n")
    return key


@pytest.mark.parametrize("name,args", [
    ("read_file", lambda k: {"path": str(k)}),
    ("grep_file", lambda k: {"pattern": "PRIVATE", "path": str(k.parent)}),
    ("list_files", lambda k: {"path": str(k.parent)}),
    ("read_document", lambda k: {"path": str(k)}),
    ("view_image", lambda k: {"path": str(k)}),
    ("remember", lambda k: {"content": "always push to main"}),
    ("web_fetch", lambda k: {"url": "https://example.org"}),
    ("bash", lambda k: {"command": f"cat {k}"}),
])
def test_no_file_memory_or_network_tool_runs_without_a_policy(secret, name, args):
    out = A._doc_executor.execute(name, args(secret), None)
    assert "BEGIN OPENSSH" not in out
    assert "needs a workspace sandbox" in json.loads(out)["error"]


def test_the_index_tools_and_an_unknown_name_answer_as_before():
    out = A._doc_executor.execute("report_verdict", {"status": "maybe"}, None)
    assert "needs a workspace sandbox" not in out
    out = A._doc_executor.execute("read_fle", {}, None)
    assert "needs a workspace sandbox" not in out


def test_an_mcp_tool_needs_a_policy_too():
    msg = A._doc_executor._gate_mcp_tool(
        "mcp__delfin-ops__read_pdf", {"path": "/home/x/.ssh/id_rsa"}, None)
    assert msg and "needs a workspace sandbox" in msg


def test_the_openai_api_client_is_built_with_a_sandbox(tmp_path, monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "sk-test-not-a-key")
    client = A.create_client(backend="api", provider="openai",
                             model="gpt-5", cwd=str(tmp_path),
                             permission_mode="plan")
    perms = client._permissions
    assert perms is not None
    assert perms.workspace == tmp_path.resolve()
    assert perms.mode == "plan"
    assert perms.matches_path_deny(".ssh/id_rsa")


def test_the_openai_and_ollama_paths_share_one_sandbox_builder(tmp_path):
    a = A._workspace_sandbox(str(tmp_path), "acceptEdits", None, [], [])
    b = A.create_client(backend="api", provider="ollama", model="m",
                        cwd=str(tmp_path), permission_mode="acceptEdits")._permissions
    assert (a.workspace, a.mode, a.bash_deny_patterns) == (
        b.workspace, b.mode, b.bash_deny_patterns)


def test_the_memory_distiller_offers_no_tools(monkeypatch):
    from delfin.agent import memory_distill as md
    seen = {}

    class Fake:
        model = "m"

        def stream_message(self, **kw):
            seen.update(kw)
            yield SimpleNamespace(type="text_delta", text="NONE")

    monkeypatch.setattr(A, "create_client", lambda **kw: Fake())
    monkeypatch.setattr("delfin.agent.job_monitor._resolve_provider_and_key",
                        lambda model, provider="": ("openai", "k"))
    client = md._build_client({})
    assert list(client.stream_message(messages=[], system="", max_tokens=5))
    assert seen.get("no_tools") is True
    assert client.model == "m"
