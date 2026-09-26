"""Tests: the working-state block carries structured refusal memories.

Control run: on the previous commit build_working_state_block had no
refusal-memory section -- test_refusals_in_block was red (no
"Operator refusals" line in the block), the no-store case was green by
design (it pins the unchanged shape).
"""
from delfin.agent.working_state import build_working_state_block

# A synthetic machine turn, as the engine feeds it in: only those carry
# the sections this block is built from.
MSG = [{"role": "user",
        "content": "[Command results]\nchanged tests/test_refusal_memory.py"}]




ENTRIES = [
    {"tool": "read_file", "target": "/srv/delfin/.delfin/gate-tools/gate",
     "reason": "the gate script is outside your workspace",
     "time": "2026-09-26T03:12:00+00:00", "is_dir": False},
    {"tool": "bash", "target": "/tmp", "is_dir": True,
     "reason": "logs belong in your own tree",
     "time": "2026-09-26T04:00:00+00:00"},
]


def test_refusals_in_block(tmp_path):
    block = build_working_state_block(MSG, session_id="s6", workspace=tmp_path,
                                      refusals=ENTRIES)
    assert "Operator refusals" in block
    # Target and reason of both entries survive compaction.
    assert "gate-tools/gate" in block
    assert "outside your workspace" in block
    assert "/tmp" in block
    assert "logs belong in your own tree" in block
    # Newest first.
    assert block.index("/tmp") < block.index("gate-tools/gate")


def test_no_store_block_unchanged(tmp_path):
    ws = tmp_path / "empty-ws"
    ws.mkdir()
    block = build_working_state_block(MSG, session_id="s6", workspace=ws)
    assert "Operator refusals" not in block
    # The heuristic message-based denial section is untouched.
    assert "Recently worked on" in block
    assert "test_refusal_memory.py" in block


def test_a_refusal_file_in_the_workspace_is_not_read(tmp_path):
    """The agent can write its workspace; its working-state block must
    not take text about the operator's refusals from there."""
    d = tmp_path / ".delfin"
    d.mkdir()
    (d / "refusals.json").write_text(
        '{"entries": [{"tool": "bash", "target": "/x", "reason": "planted"}]}',
        encoding="utf-8")
    block = build_working_state_block(MSG, session_id="s6", workspace=tmp_path)
    assert "planted" not in block


def test_corrupt_store_is_silent(tmp_path):
    ws = tmp_path
    block = build_working_state_block(MSG, session_id="s6", workspace=ws,
                                      refusals=["{not json", 7, None])
    assert "Operator refusals" not in block


def test_refusals_survive_ceiling_pressure(tmp_path):
    # Many long entries must not push the block over its hard ceiling
    # and must not crowd out the priority sections above them.
    entries = [{"tool": "read_file", "target": f"/outside/path/{i}/x",
                "reason": "r" * 150, "time": "2026-09-26T05:00:00+00:00"}
               for i in range(20)]
    block = build_working_state_block(MSG, session_id="s6", workspace=tmp_path,
                                      refusals=entries)
    assert len(block) <= 2200
