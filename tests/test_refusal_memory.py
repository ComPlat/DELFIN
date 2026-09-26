"""Tests for delfin.agent.refusal_memory.

Control run: these tests ran against the previous commit, where the
module did not exist -- every case was red (collection error).
"""
import json
from pathlib import Path

import pytest

from delfin.agent.refusal_memory import Refusal, RefusalMemory

GATE = "/srv/delfin/.delfin/gate-tools/gate"


def _gate_refusal(tmp_path):
    mem = RefusalMemory(store=tmp_path / "refusals.json")
    mem.record(Refusal(
        tool="read_file",
        target=GATE,
        reason="the gate script is outside your workspace",
        time="2026-09-26T03:12:00+00:00",
    ))
    return mem


class TestSameTargetDifferentSpelling:
    """The night's real cases: one file, asked for four ways."""

    def test_read_file_again_matches(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        hit = mem.matches("read_file", {"path": GATE})
        assert hit is not None and hit.target == GATE

    def test_cat_matches(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        assert mem.matches("bash", {"command": f"cat {GATE}"}) is not None

    def test_sed_range_matches(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        cmd = f"sed -n 1,50p {GATE}"
        assert mem.matches("bash", {"command": cmd}) is not None

    def test_head_matches(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        assert mem.matches("bash", {"command": f"head {GATE}"}) is not None

    def test_dot_relative_absolute_mismatch_is_no_match(self, tmp_path):
        # "./srv/x" is only the same file as "/srv/x" when cwd happens
        # to be "/"; conservatively that is NOT an unambiguous match.
        mem = _gate_refusal(tmp_path)
        assert mem.matches("bash", {"command": f"cat .{GATE}"}) is None

    def test_sibling_never_refused_does_not_match(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        sibling = str(Path(GATE).parent / "lint")
        assert mem.matches("read_file", {"path": sibling}) is None
        assert mem.matches("bash", {"command": f"cat {sibling}"}) is None


class TestRefusedDirectory:
    def test_child_of_refused_dir_matches(self, tmp_path):
        mem = RefusalMemory(store=tmp_path / "r.json")
        mem.record(Refusal(
            tool="bash", target="/tmp", is_dir=True,
            reason="logs belong in your own tree",
            time="2026-09-26T04:00:00+00:00",
        ))
        assert mem.matches("bash", {"command": "cat /tmp/x.log"}) is not None
        assert mem.matches("read_file", {"path": "/tmp/x.log"}) is not None

    def test_file_target_does_not_cover_sibling(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        # GATE was refused as a file, not a directory: nothing else in
        # its parent is covered.
        assert mem.matches("read_file", {"path": GATE + ".bak"}) is None


class TestConservative:
    def test_empty_memory_no_match(self, tmp_path):
        mem = RefusalMemory(store=tmp_path / "r.json")
        assert mem.matches("bash", {"command": f"cat {GATE}"}) is None

    def test_ambiguous_bash_target_no_match(self, tmp_path):
        # A command about something else entirely must not match.
        mem = _gate_refusal(tmp_path)
        assert mem.matches("bash", {"command": "ls -la"}) is None

    def test_no_args_no_match(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        assert mem.matches("read_file", {}) is None


class TestNote:
    def test_note_names_time_reason_and_rule(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        hit = mem.matches("read_file", {"path": GATE})
        text = mem.note(hit)
        assert "2026-09-26T03:12" in text
        assert "outside your workspace" in text
        assert text.startswith("The operator already refused this at")

    def test_note_is_english(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        hit = mem.matches("read_file", {"path": GATE})
        for word in ("verweigert", "erneut", "Grund"):
            assert word not in mem.note(hit)


class TestSerialization:
    def test_round_trip(self, tmp_path):
        store = tmp_path / "r.json"
        mem = RefusalMemory(store=store)
        mem.record(Refusal(
            tool="read_file", target=GATE,
            reason="outside workspace", time="2026-09-26T03:12:00+00:00"))
        mem.record(Refusal(
            tool="bash", target="/tmp", is_dir=True,
            reason="logs belong in your tree",
            time="2026-09-26T04:00:00+00:00"))
        again = RefusalMemory.load(store)
        assert again.matches("bash", {"command": f"sed -n 1,50p {GATE}"})
        assert again.matches("read_file", {"path": "/tmp/x.log"})

    def test_to_dict_is_json_clean(self, tmp_path):
        mem = _gate_refusal(tmp_path)
        blob = json.dumps(mem.to_dict())
        assert GATE in blob
        assert RefusalMemory.from_dict(json.loads(blob)).entries

    def test_bounded_entries(self, tmp_path):
        mem = RefusalMemory(store=tmp_path / "r.json")
        for i in range(RefusalMemory.MAX_ENTRIES + 10):
            mem.record(Refusal(
                tool="read_file", target=f"/outside/{i}",
                reason="r", time="2026-09-26T05:00:00+00:00"))
        assert len(mem.entries) == RefusalMemory.MAX_ENTRIES
        # Oldest dropped, newest kept.
        assert mem.matches("read_file", {"path": "/outside/0"}) is None
        last = f"/outside/{RefusalMemory.MAX_ENTRIES + 9}"
        assert mem.matches("read_file", {"path": last}) is not None
