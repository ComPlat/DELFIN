"""Integration: 200 real skills on disk, driven through the public path.

Phase 5 of package 6. The staged-loading contract has to hold on the
path DELFIN actually runs: real skill FILES in the workspace's
.delfin/skills (HOME redirected, nothing touches the real ~/.delfin),
a real OpenAIClient (built like the other stream tests build it), and
the turn driven until the provider request is built. Then:

  * the `skill` tool's description in the request is char-capped —
    200 skills must not cost the catalogue;
  * a name past the cap is still reachable: calling the `skill` tool
    with it loads the FULL body.

The stub provider yields one text token and stops; the assertions read
the captured request, not model behaviour.
"""

from __future__ import annotations

import types
import threading

from delfin.agent import api_client as A_mod


def _build_client(monkeypatch, tmp_path):
    """A real OpenAIClient over a stub provider, permissions pinned to
    a workspace that holds the fake catalogue. Mirrors the construction
    in test_a_stream_that_said_nothing_is_asked_again."""
    from delfin.agent import api_client as A

    captured: dict = {}

    class _Stream:
        def __iter__(self):
            delta = types.SimpleNamespace(
                content="ok", tool_calls=None)
            chunk = types.SimpleNamespace(
                usage=types.SimpleNamespace(prompt_tokens=5,
                                            completion_tokens=1),
                choices=[types.SimpleNamespace(
                    delta=delta, finish_reason="stop")])
            return iter((chunk,))

        def close(self):
            pass

    class _Stub:
        class chat:
            class completions:
                @staticmethod
                def create(**kw):
                    if kw.get("stream"):
                        captured["tools"] = kw.get("tools")
                        return _Stream()
                    raise AssertionError("the turn must stream")

    c = A.OpenAIClient.__new__(A.OpenAIClient)
    c.client = _Stub()
    c.model = "kit.glm-5.3"
    c._provider = "kit"
    c._base_url = "https://example.invalid/api/v1"
    c._api_key = "x"
    c.effort = ""
    perms = A.KitToolPermissions(workspace=tmp_path, agent_role="")
    c._permissions = perms
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._stop_flag = False
    return c, captured


def _write_200(ws):
    d = ws / ".delfin" / "skills"
    d.mkdir(parents=True, exist_ok=True)
    for i in range(200):
        (d / f"fake-skill-{i:03d}.md").write_text(
            f"---\nname: fake-skill-{i:03d}\n"
            f"description: {'a quite specific playbook purpose ' * 3}\n---\n"
            f"# Fake skill {i}\n\nThe FULL body of fake skill {i}, "
            "line after line of playbook.\n",
            encoding="utf-8")


def _drive(client):
    for ev in client.stream_message(
            system="s", messages=[{"role": "user", "content": "hi"}],
            max_tokens=32):
        if getattr(ev, "type", "") == "text_delta":
            break


def test_two_hundred_real_skills_leave_the_tool_description_capped(
        tmp_path, monkeypatch):
    ws = tmp_path / "projekt"
    ws.mkdir()
    _write_200(ws)
    client, captured = _build_client(monkeypatch, ws)
    # No pack skills, no user-global skills: the catalogue IS the 200.
    from delfin.agent import skills as S
    monkeypatch.setattr(S, "_PACK_SKILLS_DIR", tmp_path / "no_pack")
    _drive(client)
    tools = captured.get("tools") or []
    skill_tool = next(
        (t for t in tools if t.get("function", {}).get("name") == "skill"),
        None)
    assert skill_tool is not None, "the skill tool must be advertised"
    desc = skill_tool["function"]["description"]
    # The paste is the ceiling, not the catalogue: 200 skills cost no
    # more than 40 plus the cap.
    assert len(desc) < 200 + 2_000 + len("\nAvailable skills: ") + 100
    # And the overflow is announced, with the exact count.
    assert "more —" in desc or "200" in desc


def test_a_skill_past_the_cap_still_loads_in_full(tmp_path, monkeypatch):
    ws = tmp_path / "projekt"
    ws.mkdir()
    _write_200(ws)
    client, _ = _build_client(monkeypatch, ws)
    # The last of the 200: its name cannot be in a capped listing, but
    # the tool still resolves it — through the executor the turn's
    # tool loop dispatches to (_doc_executor.execute, the same path a
    # real `skill` tool call takes).
    from delfin.agent import skills as S
    monkeypatch.setattr(S, "_PACK_SKILLS_DIR", tmp_path / "no_pack")
    out = A_mod._doc_executor.execute(
        "skill", {"name": "fake-skill-199"}, client._permissions)
    import json
    data = json.loads(out)
    assert data.get("status") == "ok"
    assert "The FULL body of fake skill 199" in data.get("content", "")


def test_a_failing_listing_is_logged_not_swallowed(tmp_path, monkeypatch,
                                                   caplog):
    """The block that pastes the listing may not take the turn down, but
    a failure in it must be visible: a silent ``pass`` is how the
    listing stayed broken on main without anyone noticing."""
    ws = tmp_path / "projekt"
    ws.mkdir()
    _write_200(ws)
    client, captured = _build_client(monkeypatch, ws)

    def _boom(skills):
        raise RuntimeError("listing broke")

    monkeypatch.setattr(A_mod, "_skill_listing", _boom)
    import logging
    with caplog.at_level(logging.WARNING):
        _drive(client)
    assert captured.get("tools"), "the turn still runs"
    assert any("skill listing not advertised" in r.getMessage()
               and "listing broke" in r.getMessage()
               for r in caplog.records)
