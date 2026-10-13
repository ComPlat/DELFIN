"""What the TUI's window completes and picks from comes from the session.

The composer completes /commands and @paths as the user types, and the
pickers offer models, efforts and agent modes. The lists are the
session's own: the commands a typed line would run (in the order the REPL
tries them), the paths under its workspace and nowhere else, the models
its provider serves. The resume picker lists the saved sessions.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent.service import choices, protocol
from tui_service import ServiceProcess
from tui_worker_client import WorkerClient, WorkerError


def _write(path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_the_commands_are_the_ones_a_typed_line_would_run(tmp_path):
    ws = tmp_path / "workspace"
    _write(ws / ".delfin" / "commands" / "deploy.md",
           "---\ndescription: Ship the build\n---\nDeploy $ARGUMENTS\n")
    # A workspace command named like a builtin never runs: /help is the REPL's.
    _write(ws / ".delfin" / "commands" / "help.md", "---\ndescription: mine\n---\nx\n")
    _write(ws / ".delfin" / "skills" / "crest-scan" / "SKILL.md",
           "---\nname: crest-scan\ndescription: Conformer search with CREST\n---\n# CREST\n")
    with WorkerClient(tmp_path) as w:
        commands = w.request("session.commands")["commands"]
    names = [c["name"] for c in commands]
    from delfin.agent import repl_commands as rc
    builtins = [row[1] for row in rc.palette_rows()]
    assert names[:len(builtins)] == builtins
    by_name = {c["name"]: c for c in commands}
    assert by_name["/help"]["source"] == "builtin"
    assert by_name["/help"]["summary"] == "List these commands"
    assert by_name["/deploy"] == {"name": "/deploy", "summary": "Ship the build",
                                  "category": "user", "takesArgs": True,
                                  "source": "user"}
    assert by_name["/crest-scan"]["source"] == "skill"
    assert by_name["/explore"]["source"] == "subagent"
    assert names.count("/help") == 1
    assert names.index("/deploy") < names.index("/crest-scan") < names.index("/explore")


def test_paths_complete_inside_the_workspace_only(tmp_path):
    ws = tmp_path / "workspace"
    _write(ws / "delfin" / "xtb.py", "")
    _write(ws / "delfin" / "orca.py", "")
    _write(ws / "docs" / "guide.md", "")
    _write(ws / ".hidden", "")
    _write(tmp_path / "outside.txt", "")
    with WorkerClient(tmp_path) as w:
        def complete(fragment: str) -> list[str]:
            return w.request("session.completePath", {"fragment": fragment})["paths"]

        assert complete("") == ["delfin/", "docs/"]          # no dot files unasked
        assert complete(".h") == [".hidden"]                  # asked for: shown
        assert complete("delfin/x") == ["delfin/xtb.py"]
        # A folder lists what it holds (the REPL's completer offered the
        # folder itself again, so completion never went into one).
        assert complete("delfin/") == ["delfin/orca.py", "delfin/xtb.py"]
        assert complete("del") == ["delfin/"]
        assert complete("@do") == ["docs/"]                   # the @ may come along
        assert complete("../") == []                          # never outside
        assert complete("../outside") == []
        assert complete("/etc/pass") == []
        with pytest.raises(WorkerError) as bad:
            w.request("session.completePath", {})
        assert bad.value.code == protocol.INVALID_PARAMS


def test_the_pickers_offer_the_sessions_efforts_modes_and_model(tmp_path):
    with WorkerClient(tmp_path) as w:
        effort = w.request("session.choices", {"what": "effort"})
        assert [c["value"] for c in effort["choices"]] == list(choices.EFFORTS)
        assert effort["current"] == ""
        mode = w.request("session.choices", {"what": "mode"})
        assert mode["current"] == "solo"
        assert "solo" in [c["value"] for c in mode["choices"]]
        # The fake engine has no endpoint to ask: its own model, and only that.
        model = w.request("session.choices", {"what": "model"})
        assert model == {"current": "fake-model",
                         "choices": [{"value": "fake-model", "detail": "current"}]}
        with pytest.raises(WorkerError) as bad:
            w.request("session.choices", {"what": "colour"})
        assert bad.value.code == protocol.INVALID_PARAMS
        # A picker's request never holds up the next one (it has a thread
        # of its own): a turn still runs and answers come back.
        w.request("session.choices", {"what": "model"})
        start, end = w.turn("hello")
        assert w.types(start, end + 1)[-1] == "turn.ended"


def test_a_providers_model_list_shows_where_each_runs_and_what_it_costs(monkeypatch):
    from types import SimpleNamespace
    from delfin.agent import model_capabilities as mc

    calls = []
    entries = [
        {"id": "kit.qwen3-embedding-8b", "info": {"meta": {"description": ""}}},
        {"id": "kit.deepseek-v4-flash", "info": {"id": "kit.deepseek-v4-flash", "meta": {
            "description": "Large coding model\nHost: KIT • Intelligence: 💡💡💡💡 • Input: Text"}}},
        {"id": "azure.gpt-5.6-sol", "info": {"meta": {"description":
            "Flagship\nHost: Azure (EU) • Speed: ⚡️ • Price (1M tokens): ⚠️ $5.50 in, "
            "$33.00 out, $0.55 cache read • Input: Text"}}},
        {"id": "alias.medium-local", "info": {"meta": {"description":
            "This is an alias for a medium model, hosted locally at KIT. Currently uses "
            "DeepSeek V4 Flash.\nPlease use this only for custom models."}}},
    ]
    monkeypatch.setattr(mc, "_fetch_openai_models",
                        lambda base, key: calls.append(base) or entries)
    monkeypatch.setattr(choices, "_models_cache", {})
    engine = SimpleNamespace(client=SimpleNamespace(model="kit.glm-5.3-flash",
                                                    _base_url="https://kit.example/v1",
                                                    _api_key="k"))
    got = choices.choices(engine, "model")
    assert got["current"] == "kit.glm-5.3-flash"
    assert got["choices"] == [
        {"value": "kit.glm-5.3-flash", "detail": "current"},
        {"value": "kit.deepseek-v4-flash", "detail": "KIT"},
        {"value": "azure.gpt-5.6-sol",
         "detail": "Azure (EU) · $5.50 in, $33.00 out per 1M tokens"},
        {"value": "alias.medium-local", "detail": "alias of DeepSeek V4 Flash"},
    ]
    choices.choices(engine, "model")
    assert calls == ["https://kit.example/v1"]              # asked once, then reused


def _store_dir(service: ServiceProcess) -> Path:
    """Where the service's session store is: asked the way the service
    finds it, scratch state and all."""
    out = subprocess.run([sys.executable, "-c",
                          "from delfin.agent import state_paths; "
                          "state_paths.scratch_state_from_environment(); "
                          "from delfin.agent import session_store as s; "
                          "print(s._SESSIONS_DIR)"],
                         env=service.env, capture_output=True, text=True, check=True)
    return Path(out.stdout.strip())


def test_the_resume_picker_lists_saved_sessions_newest_first(tmp_path):
    with ServiceProcess(tmp_path) as service:
        ws = service.workspace()
        other = service.workspace("other")
        store = _store_dir(service)
        store.mkdir(parents=True, exist_ok=True)
        now = time.time()
        for sid, title, where, age in (("aaa1", "older one", ws, 300),
                                       ("bbb2", "newer one", ws, 60),
                                       ("ccc3", "elsewhere", other, 10)):
            (store / f"{sid}.json").write_text(json.dumps({
                "session_id": sid, "title": title, "workspace": str(where),
                "model": "kit.deepseek-v4-flash", "provider": "kit",
                "updated_at": now - age, "chat_messages": [
                    {"role": "user", "content": f"about {title}"}]}))
        c = service.client()
        rows = c.request("sessions.past", {"workspace": str(ws)})["sessions"]
        assert [r["sessionId"] for r in rows] == ["bbb2", "aaa1"]
        assert rows[0]["title"] == "newer one" and rows[0]["running"] is False
        assert rows[0]["messages"] == 1
        everywhere = c.request("sessions.past", {})["sessions"]
        assert [r["sessionId"] for r in everywhere] == ["ccc3", "bbb2", "aaa1"]
        found = c.request("sessions.past", {"query": "about older"})["sessions"]
        assert [r["sessionId"] for r in found] == ["aaa1"]
        assert c.request("sessions.past", {"limit": 1})["sessions"][0]["sessionId"] == "ccc3"


def test_ctrl_t_shows_or_hides_the_task_list_as_in_the_repl(tmp_path):
    with WorkerClient(tmp_path) as w:
        mark = len(w.events)
        assert w.request("session.toggleTasks") == {"showTasks": True}
        on, _ = w.wait_for("notice", after=mark)
        assert on["text"] == "task list on"
        mark = len(w.events)
        assert w.request("session.toggleTasks") == {"showTasks": False}
        off, _ = w.wait_for("notice", after=mark)
        assert off["text"] == "task list off"


def test_the_scripted_tools_turn_has_every_kind_of_result(tmp_path):
    with WorkerClient(tmp_path) as w:
        start, end = w.turn("show me the tools")
        ended = [e for e in w.events[start:end] if e["type"] == "tool.ended"]
        names = [e["name"] for e in ended]
        assert names == ["grep_file", "web_fetch", "bash", "write_file", "read_file",
                         "mcp__delfin-docs__search_docs", "subagent"]
        blocked = next(e for e in ended if e["name"] == "read_file")
        assert blocked["meta"]["ok"] is False and "outside the workspace" in blocked["meta"]["error"]
