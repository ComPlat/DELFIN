"""A count the answer states is a quote when a tool returned it.

Field case, 2026-09-11: read_file on a dataset description was cut
short; the answer said "800 Stichproben (800 Samples)"; the note under
it called the figure an estimate. The file's own line said 800. The
truncation guard had asked only whether SOMETHING was cut, never
whether the number itself was in what came back.
"""
import textwrap
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import verify_guard as vg
from delfin.agent.api_client import StreamEvent


@pytest.fixture(autouse=True)
def _fresh_pool():
    vg.reset_observed_numbers()
    yield
    vg.reset_observed_numbers()


TOOLS = ["mcp__delfin-docs__read_file"]


def test_a_number_the_cut_result_carried_is_not_flagged():
    vg.record_tool_numbers("title: shifts\nn_samples = 800\n"
                           "... [truncated, 40000 chars total]",
                           truncated=True)
    answer = "Der Datensatz enthält 800 Stichproben (800 Samples)."
    assert vg.scan_for_counts_over_truncated_output(answer, TOOLS) == []


def test_a_number_no_result_carried_is_still_flagged():
    vg.record_tool_numbers("a.pdf\nb.pdf\n... [truncated, 9000 chars total]",
                           truncated=True)
    answer = "Ich habe 31 PDF-Dateien verifiziert."
    assert vg.scan_for_counts_over_truncated_output(answer, TOOLS) == [
        "31 PDF-Dateien"]


def test_only_the_quoted_count_is_dropped_from_a_mixed_answer():
    vg.record_tool_numbers("n_samples = 800\n... [truncated, 40000 chars]",
                           truncated=True)
    answer = "800 Samples, verteilt auf 12 Ordner."
    assert vg.scan_for_counts_over_truncated_output(answer, TOOLS) == [
        "12 Ordner"]


def test_the_pool_answers_even_with_a_hole_in_it():
    """observed_numbers() says None after a cut; the seen-set does not."""
    vg.record_tool_numbers("x = 7\n[truncated, 5000 chars]", truncated=True)
    assert vg.observed_numbers() is None
    assert 7.0 in vg.numbers_the_tools_returned()


def test_a_fresh_turn_has_seen_nothing():
    assert vg.numbers_the_tools_returned() == frozenset()


# --- end to end through the engine -------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# quick mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


def _client(reply, tool_events):
    fake = MagicMock()
    fake._observed_files_session = set()
    calls = {"n": 0}

    def _stream(*_a, **_k):
        calls["n"] += 1
        if calls["n"] == 1:
            for ev in tool_events:
                yield ev
        yield StreamEvent(type="text_delta", text=reply)
        yield StreamEvent(type="message_delta", output_tokens=5,
                          cost_usd=0.0)

    fake.stream_message = MagicMock(side_effect=_stream)
    return fake


def _engine(tree, client):
    from delfin.agent.engine import AgentEngine
    with patch("delfin.agent.engine.create_client", return_value=client):
        return AgentEngine(repo_dir=tree, backend="cli", mode="quick",
                           pack_dir=tree)


_CUT_DATASET = StreamEvent(
    type="tool_result", tool_name="mcp__delfin-docs__read_file",
    tool_output="# NMR set\nn_samples = 800\nnuclei: 13C\n",
    output_truncated=True, output_chars=40000)


def test_the_note_stays_off_a_quoted_figure(agent_tree):
    engine = _engine(agent_tree, _client(
        "Der Datensatz hat 800 Samples.", (_CUT_DATASET,)))
    out = engine.stream_response("wie groß ist der datensatz?")
    assert "800 Samples" in out
    assert "an estimate, not a count" not in out


def test_the_note_still_reaches_a_figure_nothing_returned(agent_tree):
    engine = _engine(agent_tree, _client(
        "Ich habe 31 Rechnungen geprüft.", (_CUT_DATASET,)))
    out = engine.stream_response("wie viele rechnungen?")
    assert "an estimate, not a count" in out
