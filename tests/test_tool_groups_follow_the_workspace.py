"""The chemistry and office tool groups are offered where they apply.

Measured on four foreign repositories (Python, TypeScript, Rust, Go):
114 tools / 16,456 schema tokens per request, of which 43 tools
(mcp__delfin-ops__*, mcp__delfin-tools__*) and 6,661 tokens are DELFIN
chemistry and 12 tools / 2,028 tokens are office documents. Each test
reads the tool list off the request the client really builds.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import api_client as AC
from delfin.agent.api_client import KitToolPermissions, _OFFICE_TOOL_NAMES


def _make_delfin_workspace(base: Path) -> Path:
    (base / "delfin" / "agent").mkdir(parents=True)
    (base / "delfin" / "__init__.py").write_text("")
    (base / "delfin" / "agent" / "__init__.py").write_text("")
    return base


def _make_ts_repo(base: Path) -> Path:
    base.mkdir(parents=True)
    (base / "package.json").write_text('{"name": "app", "scripts": {"test": "vitest"}}')
    (base / "src").mkdir()
    (base / "src" / "index.ts").write_text("export const x = 1;\n")
    return base


_MCP_NAMES = ("mcp__delfin-ops__parse_orca_output",
              "mcp__delfin-tools__catalog",
              "mcp__userserver__lookup")


def _client(workspace, monkeypatch, **perm_kwargs):
    from delfin.agent import mcp_client as _mcp
    from delfin.agent import model_capabilities as _mc

    class _Registry:
        def discover_all(self):
            out = []
            for full in _MCP_NAMES:
                _, server, name = full.split("__")
                out.append(_mcp.MCPTool(server=server, name=name,
                                        description=f"{name} tool",
                                        schema={"type": "object"}))
            return out

        def discover_resources(self):
            return []

        def discover_prompts(self):
            return []

    caps = _mc.ModelCapabilities(model="m", provider="ollama",
                                 context_window=200_000, supports_tools=True)
    monkeypatch.setattr(AC._doc_executor, "_ensure_loaded", lambda: True)
    monkeypatch.setattr(AC._doc_executor, "_ensure_calc_loaded", lambda: True)
    monkeypatch.setattr(AC, "_office_backends_available", lambda: True)
    monkeypatch.setattr(AC, "_office_backend_set", lambda: None)
    monkeypatch.setattr(_mc, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(_mcp, "get_registry", lambda *a, **k: _Registry())
    monkeypatch.setattr(AC.time, "sleep", lambda *a, **k: None)
    client = AC.create_client(backend="api", provider="ollama",
                              model="qwen2.5-coder:32b", cwd=str(workspace))
    perm_kwargs.setdefault("mode", "acceptEdits")
    client.set_permissions(KitToolPermissions(workspace=workspace,
                                              **perm_kwargs))
    return client


def _offered(client, messages=None) -> list[str]:
    seen: dict = {}

    class _Delta:
        content, tool_calls, reasoning_content = "ok", None, None

    class _Choice:
        index, delta, finish_reason = 0, _Delta(), "stop"

    class _Chunk:
        choices, usage = [_Choice()], None

    class _Stream:
        def __iter__(self):
            return iter([_Chunk()])

        def close(self):
            pass

    def _create(**kwargs):
        seen["tools"] = kwargs.get("tools") or []
        return _Stream()

    client.client.chat.completions.create = _create
    msgs = messages or [{"role": "user", "content": "fix the failing test"}]
    list(client.stream_message("sys", msgs, max_tokens=32))
    assert "tools" in seen, "the client sent no request"
    return [t["function"]["name"] for t in seen["tools"]]


def _chemistry(names):
    return {n for n in names
            if n.startswith(("mcp__delfin-ops__", "mcp__delfin-tools__"))
            or n == "show_molecule"}


def test_a_foreign_repo_is_offered_neither_group(tmp_path, monkeypatch):
    names = _offered(_client(_make_ts_repo(tmp_path / "app"), monkeypatch))
    assert _chemistry(names) == set()
    assert set(names) & _OFFICE_TOOL_NAMES == set()
    # A third-party MCP server is not a DELFIN group and stays.
    assert "mcp__userserver__lookup" in names
    assert "bash" in names and "edit_file" in names


def test_a_delfin_workspace_keeps_both_groups(tmp_path, monkeypatch):
    names = _offered(_client(_make_delfin_workspace(tmp_path / "d"),
                             monkeypatch))
    assert {"mcp__delfin-ops__parse_orca_output",
            "mcp__delfin-tools__catalog", "show_molecule"} <= set(names)
    assert _OFFICE_TOOL_NAMES <= set(names)


def test_chemistry_files_in_the_workspace_bring_the_chemistry_group(
        tmp_path, monkeypatch):
    ws = tmp_path / "study"
    (ws / "chem" / "benzene").mkdir(parents=True)
    (ws / "chem" / "benzene" / "start.xyz").write_text("1\n\nC 0 0 0\n")
    names = _offered(_client(ws, monkeypatch))
    assert "mcp__delfin-ops__parse_orca_output" in names
    assert set(names) & _OFFICE_TOOL_NAMES == set()


def test_an_orca_output_is_recognised_by_its_content(tmp_path, monkeypatch):
    ws = tmp_path / "runs"
    ws.mkdir()
    (ws / "job.out").write_text(
        "                                 *****************\n"
        "                                 * O   R   C   A *\n")
    (ws / "pytest.out").write_text("collected 3 items\n")
    assert "mcp__delfin-ops__parse_orca_output" in _offered(
        _client(ws, monkeypatch))


def test_a_test_log_named_out_is_not_chemistry(tmp_path, monkeypatch):
    ws = _make_ts_repo(tmp_path / "app")
    (ws / "vitest.out").write_text("Test Files  3 passed (3)\n")
    assert _chemistry(_offered(_client(ws, monkeypatch))) == set()


def test_a_calc_directory_of_the_session_brings_the_chemistry_group(
        tmp_path, monkeypatch):
    ws = _make_ts_repo(tmp_path / "app")
    calc = tmp_path / "calc"
    calc.mkdir()
    names = _offered(_client(ws, monkeypatch, confirm_write_dirs=(calc,)))
    assert "mcp__delfin-ops__parse_orca_output" in names


def test_office_files_bring_the_office_group(tmp_path, monkeypatch):
    ws = _make_ts_repo(tmp_path / "app")
    (ws / "docs").mkdir()
    (ws / "docs" / "budget.xlsx").write_bytes(b"PK\x03\x04")
    names = _offered(_client(ws, monkeypatch))
    assert _OFFICE_TOOL_NAMES <= set(names)
    assert _chemistry(names) == set()


def test_the_office_role_always_has_the_office_group(tmp_path, monkeypatch):
    ws = tmp_path / "office"
    ws.mkdir()
    (ws / "notes.txt").write_text("x")
    client = _client(ws, monkeypatch)
    client._permissions.agent_role = "office_agent"
    names = _offered(client)
    assert "read_document" in names


def test_a_setting_turns_a_group_on_and_off(tmp_path, monkeypatch):
    from delfin import user_settings as US
    real = US.load_settings

    def _with(value):
        def _load(*a, **k):
            data = dict(real(*a, **k) or {})
            data["agent"] = {**(data.get("agent") or {}),
                             "toolsets": value}
            return data
        return _load

    ws = _make_ts_repo(tmp_path / "app")
    monkeypatch.setattr(US, "load_settings", _with({"chemistry": "on"}))
    assert "mcp__delfin-tools__catalog" in _offered(_client(ws, monkeypatch))

    d = _make_delfin_workspace(tmp_path / "d")
    monkeypatch.setattr(US, "load_settings", _with({"office": "off"}))
    names = _offered(_client(d, monkeypatch))
    assert set(names) & _OFFICE_TOOL_NAMES == set()
    assert "show_molecule" in names


def test_a_tool_the_session_already_used_is_never_withdrawn(
        tmp_path, monkeypatch):
    """A resumed history that called a chemistry tool keeps the group."""
    ws = _make_ts_repo(tmp_path / "app")
    history = [
        {"role": "user", "content": "parse that output"},
        {"role": "assistant", "content": "",
         "tool_calls": [{"id": "c1", "type": "function", "function": {
             "name": "mcp__delfin-ops__parse_orca_output",
             "arguments": "{}"}}]},
        {"role": "tool", "tool_call_id": "c1", "content": "{}"},
        {"role": "user", "content": "and now the tests"},
    ]
    names = _offered(_client(ws, monkeypatch), history)
    assert "mcp__delfin-ops__parse_orca_output" in names


def test_tools_enable_offers_the_group_from_the_next_turn(
        tmp_path, monkeypatch):
    from types import SimpleNamespace

    from delfin.agent import repl_commands as RC

    ws = _make_ts_repo(tmp_path / "app")
    client = _client(ws, monkeypatch)
    first = _offered(client)
    assert _chemistry(first) == set()

    engine = SimpleNamespace(kit_permissions=client._permissions,
                             get_status=lambda: {"role": ""})
    out = RC._tools(SimpleNamespace(engine=engine), "enable chemistry").output
    assert "chemistry" in out and "next turn" in out

    second = _offered(client)
    assert "mcp__delfin-ops__parse_orca_output" in second
    # Appending a group must not reorder what was already offered.
    assert [n for n in second if n in set(first)] == first


def test_files_written_mid_session_change_nothing(tmp_path, monkeypatch):
    """The decision is taken once; the head stays byte-stable."""
    ws = _make_ts_repo(tmp_path / "app")
    client = _client(ws, monkeypatch)
    first = _offered(client)
    (ws / "mol.xyz").write_text("1\n\nH 0 0 0\n")
    (ws / "report.docx").write_bytes(b"PK")
    assert _offered(client) == first


def test_a_hidden_group_is_named_in_one_line_with_its_switch(tmp_path):
    from delfin.agent.toolsets import WorkspaceToolsets
    ts = WorkspaceToolsets(_make_ts_repo(tmp_path / "app"),
                           settings={"chemistry": "auto", "office": "auto"})
    note = ts.hidden_note()
    assert note.count("\n") == 0
    assert "/tools enable chemistry" in note
    assert "/tools enable office" in note
