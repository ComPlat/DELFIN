"""KG3 guardrail 4: stable-prefix measurement on a RECORDED request.

The rows are NOT re-composed fresh: a real AgentEngine (throwaway pack,
fake client) runs two consecutive turns through stream_message; a
per-instance wrapper around engine.loader.compose_sections records the
section lists the engine actually built and sent per turn. Those two
recorded lists feed prompt_head.stable_prefix_len / check_ordering.
"""

import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import prompt_head as ph


@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def _engine(agent_tree):
    from delfin.agent.engine import AgentEngine
    fake = MagicMock()
    fake.stream_message = MagicMock(side_effect=lambda *a, **k: iter(()))
    with patch("delfin.agent.engine.create_client", return_value=fake):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                           mode="quick", pack_dir=agent_tree)


def test_stable_prefix_measured_on_recorded_engine_turns(agent_tree):
    engine = _engine(agent_tree)
    loader = engine.loader
    recorded = []
    original = loader.compose_sections

    def recorder(**kwargs):
        sections = original(**kwargs)
        # Record what the engine actually built THIS turn: the
        # PromptSection objects themselves, not a reshaped copy.
        recorded.append(list(sections))
        return sections

    loader.compose_sections = recorder

    # Two consecutive real turns (the client streams nothing, the turn
    # machinery still runs the prompt build).
    engine.stream_response("first turn")
    engine.stream_response("second turn")

    assert len(recorded) == 2, "two engine turns must each record a head"
    turn_a, turn_b = recorded[0], recorded[1]

    report = ph.stable_prefix_len(turn_a, turn_b)
    # Byte counts of the recorded rows — the before/after the reviewer
    # wants stated (turn A prefix bytes, turn B full head bytes).
    report_b = ph.stable_prefix_len(turn_b, turn_a)
    head_bytes = sum(len(r.content) for r in turn_b)
    print(f"recorded: {len(turn_a)} vs {len(turn_b)} rows; "
          f"shared prefix {report['rows']} rows / {report['bytes']} bytes; "
          f"turn-B head total {head_bytes} bytes; "
          f"rows_all_stable={report['rows_all_stable']}")

    # Safety ordering on the RECORDED rows.
    assert ph.check_ordering(turn_b) == []

    # The stable rows the engine composed sit inside the shared prefix —
    # the held-prompt mechanism keeps the system prompt byte-stable, so
    # every stable row is cached.
    stable_names_b = {
        r.name for r in turn_b if getattr(r, "layer", None) == 0}
    shared_names = {turn_b[i].name for i in range(report["rows"])}
    assert stable_names_b <= shared_names, (
        stable_names_b - shared_names)
    assert report["rows_all_stable"] is True, report
