"""Changing the permissions mid-turn told the gate and not the model.

``set_kit_permission_mode`` writes ``perms.mode`` from whichever thread
the dashboard runs on. The tool gate reads it on its very next check, so
the new rules are in force immediately -- and the model is still working
to the ones it was told at the top of the turn, because the system prompt
is built once and frozen for the whole turn.

Both directions of that gap are bad in their own way. Switched to
``plan``, the model keeps trying to write and every attempt comes back
refused for a reason it was never given. Switched to ``acceptEdits``, it
keeps proposing and asking, and the user watches it request permission it
has already been granted.

The fix is the rail this framework already uses for the two other things
that happen mid-turn: what the user said, and what a background job
finished. A permission change is the third, and it now rides the same
one -- injected between tool rounds, and again if the turn ends without
a tool call.
"""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent.api_client import StreamEvent


@pytest.fixture
def agent_tree(tmp_path):
    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# solo mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def _engine(agent_tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick",
                           pack_dir=agent_tree)


# ---------------------------------------------------------------------------
# The queue itself
# ---------------------------------------------------------------------------

def _real_client(tmp_path, monkeypatch):
    """A real client, so the queue under test is the shipped one."""
    from delfin.agent import model_capabilities as mc
    import delfin.agent.api_client as A
    caps = mc.ModelCapabilities(model="m", provider="ollama",
                                context_window=200_000, supports_tools=True)
    monkeypatch.setattr(mc, "resolve", lambda *a, **k: caps)
    return A.create_client(backend="api", provider="ollama",
                           model="qwen2.5-coder:7b", cwd=str(tmp_path))


def test_a_run_note_is_queued_and_drained_once(tmp_path, monkeypatch):
    client = _real_client(tmp_path, monkeypatch)
    client.push_run_note("the rules changed")
    assert client._drain_run_notes() == ["the rules changed"]
    # Exactly once: a note delivered twice would have the model announce
    # the same change on two rounds running.
    assert client._drain_run_notes() == []


def test_an_empty_note_is_not_queued(tmp_path, monkeypatch):
    client = _real_client(tmp_path, monkeypatch)
    client.push_run_note("")
    client.push_run_note("   ")
    assert client._drain_run_notes() == []


def test_steering_and_run_notes_do_not_consume_each_other(tmp_path,
                                                          monkeypatch):
    """They share a lock and nothing else. A mid-run message from the user
    and a permission change are different facts and both have to arrive."""
    client = _real_client(tmp_path, monkeypatch)
    client.push_steer("mach weiter")
    client.push_run_note("[permissions] switched to plan")
    assert client._drain_steer() == ["mach weiter"]
    assert client._drain_run_notes() == ["[permissions] switched to plan"]


# ---------------------------------------------------------------------------
# ... and the change that fills it
# ---------------------------------------------------------------------------

def test_switching_the_mode_mid_turn_queues_a_note(agent_tree):
    fake = MagicMock()
    fake._observed_files_session = set()
    fake._permissions = MagicMock(mode="default")
    engine = _engine(agent_tree, fake)

    assert engine.set_kit_permission_mode("plan") is True
    fake.push_run_note.assert_called_once()
    said = fake.push_run_note.call_args[0][0]
    assert "'default'" in said and "'plan'" in said
    assert "mid-turn" in said


def test_setting_the_mode_it_already_has_says_nothing(agent_tree):
    """A no-op is not news. Announcing it would train the reader to skip
    the announcements that matter."""
    fake = MagicMock()
    fake._observed_files_session = set()
    fake._permissions = MagicMock(mode="acceptEdits")
    engine = _engine(agent_tree, fake)

    assert engine.set_kit_permission_mode("acceptEdits") is True
    fake.push_run_note.assert_not_called()


def test_an_unknown_mode_changes_nothing_and_says_nothing(agent_tree):
    fake = MagicMock()
    fake._observed_files_session = set()
    fake._permissions = MagicMock(mode="default")
    engine = _engine(agent_tree, fake)

    assert engine.set_kit_permission_mode("whatever") is False
    fake.push_run_note.assert_not_called()


def test_a_client_without_the_rail_still_switches(agent_tree):
    """The mode change is the point; the notice is the courtesy. An older
    client must not turn a permission change into an exception."""
    fake = MagicMock()
    fake._observed_files_session = set()
    fake._permissions = MagicMock(mode="default")
    fake.push_run_note = MagicMock(side_effect=AttributeError("no rail"))
    engine = _engine(agent_tree, fake)

    assert engine.set_kit_permission_mode("plan") is True
    assert fake._permissions.mode == "plan"


# ---------------------------------------------------------------------------
# ... reaches the person too
# ---------------------------------------------------------------------------

def test_the_note_is_shown_and_is_not_the_model_answer(agent_tree):
    """It rides the notice type, like every other sentence the framework
    speaks: visible to the user, and never part of what gets scored as
    the model's answer."""
    seen: list[str] = []
    fake = MagicMock()
    fake._observed_files_session = set()
    fake.stream_message = MagicMock(side_effect=lambda *a, **k: iter([
        StreamEvent(type="notice",
                    text="\n\n[permissions] switched to 'plan' mid-turn\n"),
        StreamEvent(type="text_delta", text="Verstanden, ich plane nur."),
        StreamEvent(type="message_delta", output_tokens=6, cost_usd=0.0),
    ]))
    engine = _engine(agent_tree, fake)
    out = engine.stream_response("mach was", on_token=seen.append)

    assert any("[permissions]" in s for s in seen)
    assert "[permissions]" not in (out or "")
    assert "ich plane nur" in out


# ---------------------------------------------------------------------------
# One call, one rulebook
# ---------------------------------------------------------------------------

def test_a_flip_between_two_gate_checks_cannot_decide_one_call(tmp_path):
    """The gate read ``perms.mode`` two or three times per call -- once at
    the top of _run_permission_gate, again for the network check, and a
    third time inside apply_patch and undo_changes AFTER the gate had
    already decided. A dashboard flip landing between two of those reads
    decided one call by two different policies, and the direction that
    matters is the one that lets a write through.

    The mode is pinned for the length of the call now, so a change made
    from another thread applies to the NEXT call and never to half of
    this one.
    """
    from delfin.agent import api_client as ac

    perms = ac.KitToolPermissions(workspace=tmp_path, mode="acceptEdits")
    seen: list[str] = []
    with ac._pin_permission_mode(perms):
        seen.append(ac.effective_mode(perms))
        perms.mode = "plan"                 # the dashboard, mid-call
        seen.append(ac.effective_mode(perms))
    # ... and the next call sees the change.
    seen.append(ac.effective_mode(perms))
    assert seen == ["acceptEdits", "acceptEdits", "plan"], seen


def test_the_agent_changing_it_itself_takes_effect_at_once(tmp_path):
    """exit_plan_mode is the agent deciding, not somebody else deciding
    around it -- so the rest of the call runs under the new mode."""
    from delfin.agent import api_client as ac

    perms = ac.KitToolPermissions(workspace=tmp_path, mode="plan")
    with ac._pin_permission_mode(perms):
        assert ac.effective_mode(perms) == "plan"
        perms.mode = "acceptEdits"
        ac.repin_permission_mode("acceptEdits")
        assert ac.effective_mode(perms) == "acceptEdits"


def test_without_a_pin_the_live_value_is_the_answer(tmp_path):
    """Outside a tool call there is nothing to be consistent WITH, and a
    stale pin would be worse than none."""
    from delfin.agent import api_client as ac

    perms = ac.KitToolPermissions(workspace=tmp_path, mode="default")
    assert ac.effective_mode(perms) == "default"
    perms.mode = "plan"
    assert ac.effective_mode(perms) == "plan"


# ---------------------------------------------------------------------------
# The other two things the dashboard changes mid-turn
# ---------------------------------------------------------------------------

def test_a_directory_granted_mid_turn_reaches_the_model(agent_tree, tmp_path):
    """An agent refused a write two rounds ago has already told the user it
    cannot do the thing. Granting the directory silently leaves that
    standing -- the gate opens and nobody tries the door again."""
    from delfin.agent import api_client as ac

    target = tmp_path / "daten"
    target.mkdir()
    fake = MagicMock()
    fake._observed_files_session = set()
    fake._permissions = ac.KitToolPermissions(workspace=agent_tree)
    engine = _engine(agent_tree, fake)

    ok, _msg = engine.add_kit_workspace_dir(str(target), persist=False)
    assert ok
    fake.push_run_note.assert_called_once()
    said = fake.push_run_note.call_args[0][0]
    assert "[permissions]" in said and str(target) in said
    # ... and it says not to redo what already worked.
    assert "redo" in said


def test_a_command_allowed_for_good_reaches_the_model(agent_tree,
                                                      monkeypatch):
    from delfin.agent import api_client as ac
    from delfin.agent import kit_settings as ks

    monkeypatch.setattr(ks, "persist_pattern", lambda *a, **k: None)
    fake = MagicMock()
    fake._observed_files_session = set()
    fake._permissions = ac.KitToolPermissions(workspace=agent_tree)
    engine = _engine(agent_tree, fake)

    ok, _msg = engine.persist_kit_pattern("pytest -q", kind="allow")
    assert ok
    said = fake.push_run_note.call_args[0][0]
    assert "allowed" in said and "pytest -q" in said


def test_a_command_denied_for_good_reaches_the_model(agent_tree, monkeypatch):
    """The direction that matters more: the agent must stop trying."""
    from delfin.agent import api_client as ac
    from delfin.agent import kit_settings as ks

    monkeypatch.setattr(ks, "persist_pattern", lambda *a, **k: None)
    fake = MagicMock()
    fake._observed_files_session = set()
    fake._permissions = ac.KitToolPermissions(workspace=agent_tree)
    engine = _engine(agent_tree, fake)

    ok, _msg = engine.persist_kit_pattern("rm -rf", kind="deny")
    assert ok
    said = fake.push_run_note.call_args[0][0]
    assert "denied" in said and "not attempt it again" in said


def test_the_rail_holds_under_several_writers_at_once(tmp_path, monkeypatch):
    """The dashboard, a confirm dialog and a settings change can all push
    while the turn is reading. A lost note is a rule the model never hears
    about; a doubled one has it announce the same change twice."""
    import threading

    client = _real_client(tmp_path, monkeypatch)
    drained: list[str] = []

    def _write(worker: int) -> None:
        for i in range(200):
            client.push_run_note(f"note {worker}-{i}")

    def _read() -> None:
        for _ in range(400):
            drained.extend(client._drain_run_notes())

    threads = [threading.Thread(target=_write, args=(k,)) for k in range(4)]
    threads.append(threading.Thread(target=_read))
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    drained.extend(client._drain_run_notes())

    assert len(drained) == 800
    assert len(set(drained)) == 800


# ---------------------------------------------------------------------------
# A model switch carries its own consequences
# ---------------------------------------------------------------------------

def test_switching_the_model_resizes_the_window_and_says_so(agent_tree,
                                                            monkeypatch):
    """Two things depended on every caller remembering them, and a
    sentence in a docstring is not a mechanism: re-size the compaction
    budget to the new model's real context window, and tell a turn in
    flight that it is now talking to something else."""
    from delfin.agent import api_client as ac

    client = ac.OpenAIClient.__new__(ac.OpenAIClient)
    client.model = "small-model"
    fake_engine_calls = {"refresh": 0, "notes": []}

    def _on_switch(previous, current):
        fake_engine_calls["refresh"] += 1
        fake_engine_calls["notes"].append((previous, current))

    client.on_model_switched = _on_switch
    client.switch_model("big-model")

    assert client.model == "big-model"
    assert fake_engine_calls["refresh"] == 1
    assert fake_engine_calls["notes"] == [("small-model", "big-model")]


def test_switching_to_the_same_model_does_nothing(agent_tree):
    from delfin.agent import api_client as ac

    client = ac.OpenAIClient.__new__(ac.OpenAIClient)
    client.model = "same"
    called: list[tuple] = []
    client.on_model_switched = lambda p, c: called.append((p, c))
    client.switch_model("same")
    client.switch_model("")
    assert called == []


def test_a_client_nobody_installed_a_probe_on_still_switches():
    """The switch is the point; the consequences are the courtesy."""
    from delfin.agent import api_client as ac

    client = ac.OpenAIClient.__new__(ac.OpenAIClient)
    client.model = "a"
    client.switch_model("b")
    assert client.model == "b"


def test_a_probe_that_raises_does_not_break_the_switch():
    from delfin.agent import api_client as ac

    client = ac.OpenAIClient.__new__(ac.OpenAIClient)
    client.model = "a"
    client.on_model_switched = MagicMock(side_effect=RuntimeError("boom"))
    client.switch_model("b")
    assert client.model == "b"


def test_the_engine_installs_the_probe_and_it_reaches_the_turn(agent_tree,
                                                               monkeypatch):
    """End to end: a switch on the client the engine owns re-sizes the
    window and puts a note in front of the running turn."""
    from delfin.agent import api_client as ac
    from delfin.agent.engine import AgentEngine

    real = ac.OpenAIClient.__new__(ac.OpenAIClient)
    real.model = "small-model"
    real.push_run_note = MagicMock()
    real._observed_files_session = set()

    with patch("delfin.agent.engine.create_client", return_value=real):
        engine = AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick",
                             pack_dir=agent_tree)
    refreshed: list[int] = []
    monkeypatch.setattr(engine, "_refresh_context_window",
                        lambda *a, **k: refreshed.append(1))

    real.switch_model("big-model")

    assert refreshed == [1], "the compaction budget was not re-sized"
    said = real.push_run_note.call_args[0][0]
    assert "[model]" in said
    assert "small-model" in said and "big-model" in said


# ---------------------------------------------------------------------------
# Stop, with something still in the queue
# ---------------------------------------------------------------------------

def _stop_notice(client) -> str:
    """Drive one round of the tool loop with a stop pending."""
    client.should_stop = lambda: True
    events = list(client.stream_message(
        "sys", [{"role": "user", "content": "go"}], max_tokens=64))
    return "".join(e.text or "" for e in events if e.type == "notice")


def test_a_message_typed_just_before_a_stop_is_handed_back(tmp_path,
                                                            monkeypatch):
    """Measured before the fix: it stayed in the queue and was injected
    into a later, unrelated turn, out of order with whatever the user
    typed after pressing Stop. Dropping it instead would be a message the
    person wrote that nobody ever saw."""
    client = _real_client(tmp_path, monkeypatch)
    client.push_steer("mach stattdessen X")

    notice = _stop_notice(client)
    assert "Stopped" in notice
    assert "mach stattdessen X" in notice
    assert "send it again" in notice
    # ... and it is gone from the queue, so no later turn replays it.
    assert client._drain_steer() == []


def test_a_run_note_pending_at_a_stop_is_not_destroyed_by_it(tmp_path,
                                                             monkeypatch):
    """A run note was treated as the same kind of thing as a typed
    message, and it is not.

    Both queues were drained by one line and quoted in one sentence, so a
    stop produced: Not delivered … “[permissions] The user switched this
    session from 'plan' to 'acceptEdits' just now, mid-turn.” — send it
    again if it still applies. The user is invited to re-send a line they
    never wrote; and the note, which is a fact about the session and is
    still true after the stop, is thrown away -- taking with it the only
    notice the model was ever going to get that its permissions changed.

    The note outlives the stop and is delivered by the next turn."""
    client = _real_client(tmp_path, monkeypatch)
    client.push_run_note("[permissions] switched to plan")

    notice = _stop_notice(client)
    assert client._drain_run_notes() == ["[permissions] switched to plan"]
    # Not quoted at the user as something they wrote.
    assert "[permissions] switched to plan" not in notice
    assert "send it again" not in notice


def test_a_stop_says_that_something_changed_underneath_it(tmp_path,
                                                          monkeypatch):
    """Keeping the note is not the same as keeping the user in the dark:
    the turn they stopped was running under rules that changed, and that
    is worth one clause."""
    client = _real_client(tmp_path, monkeypatch)
    client.push_run_note("[permissions] switched to plan")

    notice = _stop_notice(client)
    assert "1 change" in notice and "next message" in notice


def test_the_two_queues_are_reported_apart(tmp_path, monkeypatch):
    """Both pending at once: the typed message comes back to the user to
    re-send, the run note does not."""
    client = _real_client(tmp_path, monkeypatch)
    client.push_steer("mach stattdessen X")
    client.push_run_note("[permissions] switched to plan")

    notice = _stop_notice(client)
    assert "mach stattdessen X" in notice and "send it again" in notice
    assert "[permissions] switched to plan" not in notice
    assert client._drain_steer() == []
    assert client._drain_run_notes() == ["[permissions] switched to plan"]


def test_a_stop_with_an_empty_queue_says_only_that(tmp_path, monkeypatch):
    """The half that keeps the stop message readable."""
    client = _real_client(tmp_path, monkeypatch)
    notice = _stop_notice(client)
    assert "Stopped" in notice
    assert "Not delivered" not in notice
