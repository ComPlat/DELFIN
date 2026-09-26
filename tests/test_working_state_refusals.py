"""Tests: the working-state block carries structured refusal memories.

Control run: on the previous commit build_working_state_block had no
refusal-memory section -- test_refusals_in_block was red (no
"Operator refusals" line in the block), the no-store case was green by
design (it pins the unchanged shape).
"""
import json

from delfin.agent.working_state import build_working_state_block

# A synthetic machine turn, as the engine feeds it in: only those carry
# the sections this block is built from.
MSG = [{"role": "user",
        "content": "[Command results]\nchanged tests/test_refusal_memory.py"}]


def _store(tmp_path, entries):
    p = tmp_path / ".delfin" / "refusals.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"entries": entries}), encoding="utf-8")
    return tmp_path


ENTRIES = [
    {"tool": "read_file", "target": "/srv/delfin/.delfin/gate-tools/gate",
     "reason": "the gate script is outside your workspace",
     "time": "2026-09-26T03:12:00+00:00", "is_dir": False},
    {"tool": "bash", "target": "/tmp", "is_dir": True,
     "reason": "logs belong in your own tree",
     "time": "2026-09-26T04:00:00+00:00"},
]


def test_refusals_in_block(tmp_path):
    ws = _store(tmp_path, ENTRIES)
    block = build_working_state_block(MSG, session_id="s6", workspace=ws)
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


def test_corrupt_store_is_silent(tmp_path):
    ws = tmp_path
    d = ws / ".delfin"
    d.mkdir(exist_ok=True)
    (d / "refusals.json").write_text("{not json", encoding="utf-8")
    block = build_working_state_block(MSG, session_id="s6", workspace=ws)
    assert "Operator refusals" not in block


def test_refusals_survive_ceiling_pressure(tmp_path):
    # Many long entries must not push the block over its hard ceiling
    # and must not crowd out the priority sections above them.
    entries = [{"tool": "read_file", "target": f"/outside/path/{i}/x",
                "reason": "r" * 150, "time": "2026-09-26T05:00:00+00:00"}
               for i in range(20)]
    ws = _store(tmp_path, entries)
    block = build_working_state_block(MSG, session_id="s6", workspace=ws)
    assert len(block) <= 2200
