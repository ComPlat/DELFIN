"""The sickest turn produced the healthiest-looking record.

Live, 2026-08-12: two office tasks came back at quality 35, rate 0%, no
tool call, ~1.4 s. Their entire recorded model output was

    ⏳ Transient API error (BadRequestError); retrying 1/3 in 2s…
    ⏳ Transient API error (BadRequestError); retrying 2/3 in 3s…
    ⏳ Transient API error (BadRequestError); retrying 3/3 in 6s…

with "No deployments available for selected model" in the error field.
The endpoint had no capacity, nothing about the model was observed, and
the suite booked it as a model that answered badly -- at a quality
number, into the file baselines are compared against.

The cause is one line of plumbing: the harness's own sentences were
emitted as ``text_delta``, the same event the model's words arrive on.
So they became the assistant's message, were scored as output, and
stamped time-to-first-token. A turn that produced nothing looked like a
turn that produced something, and every consumer downstream drew the
wrong conclusion from the same wrong premise.

They now carry their own event type. Shown to the user, and to nothing
else.
"""

from __future__ import annotations

import pytest

from delfin.agent.api_client import StreamEvent


# ---------------------------------------------------------------------------
# Which sentences are the harness's own
# ---------------------------------------------------------------------------

# The retry banner -- the sentence from the live incident -- is asserted
# on the real client's event stream in test_transient_retry.py, where a
# request is actually made to fail. What follows is the other half: what
# the ENGINE does with a notice once one arrives, which is where the
# scoring, the answer and the timing were taken from.


# ---------------------------------------------------------------------------
# ... and what the engine does with them
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    import textwrap
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


def _engine_over(events, agent_tree):
    from unittest.mock import MagicMock, patch
    from delfin.agent.engine import AgentEngine

    fake = MagicMock()
    fake._observed_files_session = set()
    fake.stream_message = MagicMock(side_effect=lambda *a, **k: iter(events))
    with patch("delfin.agent.engine.create_client", return_value=fake):
        return AgentEngine(repo_dir=agent_tree, backend="cli", mode="quick",
                           pack_dir=agent_tree)


_RETRIES = [
    StreamEvent(type="notice", text="\n⏳ Transient API error "
                "(BadRequestError); retrying 1/3 in 2s…\n"),
    StreamEvent(type="notice", text="\n⏳ Transient API error "
                "(BadRequestError); retrying 2/3 in 3s…\n"),
    StreamEvent(type="notice", text="\n⏳ Transient API error "
                "(BadRequestError); retrying 3/3 in 6s…\n"),
    StreamEvent(type="message_delta", output_tokens=0, cost_usd=0.0),
]


def test_a_turn_of_nothing_but_banners_answers_nothing(agent_tree):
    """This is the whole finding. Three banners and no model word must not
    add up to an answer -- because an answer is what gets scored."""
    engine = _engine_over(_RETRIES, agent_tree)
    out = engine.stream_response("was ist die S1-Energie?")
    assert "Transient API error" not in (out or "")


def test_the_user_still_sees_every_one_of_them(agent_tree):
    """Excluded from the answer, never from the person waiting. A silent
    retry would be the same defect wearing the other coat."""
    engine = _engine_over(_RETRIES, agent_tree)
    seen: list[str] = []
    engine.stream_response("was ist die S1-Energie?", on_token=seen.append)
    assert sum("Transient API error" in s for s in seen) == 3


def test_a_banner_does_not_stamp_time_to_first_token(agent_tree, monkeypatch):
    """The record said the model started answering promptly. It never
    answered at all.

    Read back from the turn LOG. The first version of this asked the
    engine for an attribute that does not exist, so getattr returned None
    and the assertion passed against nothing -- caught while writing the
    log tests, and worth the note: a test that cannot fail is the thing
    it was written to prevent.
    """
    import json
    from delfin.agent import turn_metrics as tm

    monkeypatch.setattr(tm, "_DIR", agent_tree / "turns")
    engine = _engine_over(_RETRIES, agent_tree)
    engine.stream_response("was ist die S1-Energie?")

    entries = []
    for fp in sorted((agent_tree / "turns").glob("*.jsonl")):
        for line in fp.read_text(encoding="utf-8").splitlines():
            try:
                entries.append(json.loads(line))
            except Exception:
                continue
    assert entries, "the turn was not logged at all"
    assert entries[-1].get("ttft_ms") is None, entries[-1]
    assert entries[-1].get("output_chars") == 0, entries[-1]


def test_the_model_answer_still_arrives_normally(agent_tree):
    """The half that keeps this from being a way to lose answers."""
    engine = _engine_over([
        StreamEvent(type="notice", text="\n⏳ Transient API error; retrying…\n"),
        StreamEvent(type="text_delta", text="Die S1-Energie ist 2.31 eV."),
        StreamEvent(type="message_delta", output_tokens=9, cost_usd=0.0),
    ], agent_tree)
    out = engine.stream_response("was ist die S1-Energie?")
    assert "2.31 eV" in out
    assert "Transient API error" not in out
