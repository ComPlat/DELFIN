"""A turn that measures, compares or lands a change points at the
``experiment-discipline`` skill; other turns carry nothing.

Pinned here:

* the skill loads through the existing loader, its listing entry is one
  line under the per-skill cap, and it names the commands that exist
  (``delfin-agent provenance`` / ``compare-arms`` with flags the parser
  accepts);
* the rules the skill exists for are in it: pre-registration, one
  instrument state, null -> value, never-worse cases, the user's own
  message as the only approval, the measure versus the object, the root
  cause;
* messages that ask for an A/B run, a benchmark, "is it better now", a
  landing or a green test suite get the pointer, in English and German;
  a chemistry optimisation and ordinary requests get nothing;
* the engine puts the pointer into the turn's steering blocks only for a
  message that called for it, and never in a session without skills.
"""

from __future__ import annotations

import re
import shlex
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import api_client as A
from delfin.agent import cli as C
from delfin.agent import experiment_hints as H
from delfin.agent import skills as S
from delfin.agent.api_client import KitToolPermissions

SKILL_FILE = (Path(A.__file__).resolve().parent / "pack" / "skills"
              / "experiment-discipline.md")


def _skill():
    sk = S.get_skill(H.SKILL_NAME)
    assert sk is not None
    return sk


def test_the_skill_loads_through_the_existing_loader():
    sk = _skill()
    assert Path(sk.source).resolve() == SKILL_FILE
    assert sk.domains == ("code",)
    assert 0 < len(sk.description) <= A._SKILL_LISTING_DESC_CHARS
    listing = A._skill_listing(S.discover_skills(domain="code"))
    assert f"{sk.name} — {sk.description}" in listing


@pytest.mark.parametrize("rule", [
    "Pre-register", "One instrument state", "null -> value",
    "never-worse", "Only the user's own message approves",
    "Changing the check", "Repair, never discard",
    "symptom -> cause -> mechanism", "states the root cause",
    "Done by hand three times", "Search first",
])
def test_the_rules_are_in_the_skill(rule):
    assert rule in _skill().body


def test_every_command_the_skill_names_parses():
    parser = C.build_parser()
    commands = re.findall(r"`(delfin-agent (?:provenance|compare-arms)[^`]*)`",
                          _skill().body)
    assert len(commands) >= 3
    for cmd in commands:
        words = shlex.split(re.sub(r"<[^>]+>", "x", cmd))[1:]
        parser.parse_args(words)


@pytest.mark.parametrize("text", [
    "Run an A/B of the new cache against main.",
    "Benchmark both branches and tell me which is faster.",
    "Is it better now?",
    "Ist es jetzt besser?",
    "Did the retry change help?",
    "Compare the two runs before and after the fix.",
    "Vergleich die beiden Läufe mit der Baseline.",
    "Make the failing tests pass.",
    "Mach die Tests wieder grün machen bitte",
    "Can we land this change? The numbers look good.",
    "Start the optimisation loop on the scoring function.",
    "Please pre-register the criterion before we build it.",
    "Run a control run first.",
])
def test_an_experiment_message_gets_the_pointer(text):
    assert H.matches(text), text
    assert f"skill(name='{H.SKILL_NAME}')" in H.user_hint(text)


@pytest.mark.parametrize("text", [
    "Optimise the geometry of the Cr complex with r2SCAN-3c.",
    "Optimiere die Struktur mit xtb.",
    "Make the energy table a CSV, please.",
    "Wie viele Rechnungen sind noch nicht fertig?",
    "Fix the typo in the README.",
    "What does the function parse_block do?",
    "Write a test for the new parser.",
    "Hello!",
])
def test_other_messages_get_nothing(text):
    assert not H.matches(text), text
    assert H.user_hint(text) == ""


@pytest.fixture
def world(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.delenv("DELFIN_HOME", raising=False)
    ws = tmp_path / "ws"
    ws.mkdir()
    pack = tmp_path / "packs"
    shared = pack / "pack" / "shared"
    shared.mkdir(parents=True)
    agents = pack / "pack" / "agents"
    agents.mkdir()
    for name in ("delfin_context.md", "work_cycle_rules.md",
                 "goal_decomposition_rules.md",
                 "universal_input_template.md", "minimal_final_verdict.md"):
        (shared / name).write_text("# x")
    for name in ("session_manager", "builder_agent", "test_agent"):
        (agents / f"{name}.md").write_text(f"# {name}")
    modes = pack / "pack_lite" / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (pack / "pack_lite" / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
              - builder_agent
              - test_agent
    """))
    return {"ws": ws, "pack": pack}


def _engine(world):
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    def fake_stream(system, messages, max_tokens=4096, session_id="",
                    thinking_budget=0):
        yield StreamEvent(type="session_init", text="exp")
        yield StreamEvent(type="message_start", input_tokens=10)
        yield StreamEvent(type="text_delta", text="done")
        yield StreamEvent(type="message_delta", output_tokens=5, cost_usd=0.0)

    client = MagicMock()
    client.stream_message = MagicMock(side_effect=fake_stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=world["pack"], backend="cli",
                             mode="quick", pack_dir=world["pack"])
    client._permissions = KitToolPermissions(workspace=str(world["ws"]),
                                             mode="acceptEdits")
    return engine


def test_the_turn_carries_the_pointer_only_when_called_for(world):
    engine = _engine(world)
    engine.stream_response("Run an A/B of the cache switch and compare the runs.")
    blocks = dict(engine._steering_blocks("solo_agent"))
    assert f"[{H.SKILL_NAME}]" in blocks.get("experiment", "")
    engine.stream_response("Optimise the geometry of the Cr complex.")
    assert "experiment" not in dict(engine._steering_blocks("solo_agent"))


def test_a_session_without_skills_gets_no_pointer(world):
    engine = _engine(world)
    engine.kit_permissions.skip_skill_discovery = True
    engine.stream_response("Is it better now?")
    assert "experiment" not in dict(engine._steering_blocks("solo_agent"))


def test_the_pointer_costs_what_its_docstring_says():
    assert len(H.user_hint("Is it better now?")) // 4 == 71
