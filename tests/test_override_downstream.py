"""Tests for occupier-override downstream chain detection and override marker."""
import re
import tempfile
from pathlib import Path

import pytest

from delfin.cli import (
    _downstream_stages,
    _invalidate_stage_files,
    _OVERRIDE_MARKER_RE,
    _PREFERRED_INDEX_RE,
    _parse_override_step_list,
)


class TestParseOverrideStepList:
    def test_empty(self):
        assert _parse_override_step_list("") == []
        assert _parse_override_step_list(None) == []

    def test_comma_separated(self):
        assert _parse_override_step_list("1,2") == [1, 2]

    def test_semicolon_separated(self):
        assert _parse_override_step_list("1;2;3") == [1, 2, 3]

    def test_int_input(self):
        assert _parse_override_step_list(2) == [2]

    def test_deduplication(self):
        assert _parse_override_step_list("1,2,1") == [1, 2]


class TestDownstreamStages:
    CONFIG = {"oxidation_steps": "1,2", "reduction_steps": "1,2,3"}

    def test_initial_override(self):
        ds = _downstream_stages("initial_OCCUPIER", self.CONFIG)
        assert "ox_step_1" in ds
        assert "ox_step_2" in ds
        assert "red_step_1" in ds
        assert "red_step_2" in ds
        assert "red_step_3" in ds

    def test_ox_step_1_override(self):
        ds = _downstream_stages("ox_step_1_OCCUPIER", self.CONFIG)
        assert ds == ["ox_step_2"]

    def test_red_step_1_override(self):
        ds = _downstream_stages("red_step_1_OCCUPIER", self.CONFIG)
        assert ds == ["red_step_2", "red_step_3"]

    def test_last_ox_step(self):
        ds = _downstream_stages("ox_step_2_OCCUPIER", self.CONFIG)
        assert ds == []

    def test_last_red_step(self):
        ds = _downstream_stages("red_step_3_OCCUPIER", self.CONFIG)
        assert ds == []

    def test_no_steps_configured(self):
        ds = _downstream_stages("initial_OCCUPIER", {})
        assert ds == []


class TestInvalidateStageFiles:
    def test_removes_out_inp_fprint(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "red_step_1.out").write_text("data")
            (root / "red_step_1.inp").write_text("data")
            (root / "red_step_1.inp.fprint").write_text("hash")

            _invalidate_stage_files("red_step_1", root)

            assert not (root / "red_step_1.out").exists()
            assert not (root / "red_step_1.inp").exists()
            assert not (root / "red_step_1.inp.fprint").exists()

    def test_nonexistent_files_no_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            _invalidate_stage_files("red_step_99", root)


def _apply_override_marker(content: str, new_index: int, date_str: str) -> str:
    """Simulate the override marker logic from cli.py."""
    new_content = _PREFERRED_INDEX_RE.sub(rf"\g<1>{new_index}", content, count=1)
    new_content = _OVERRIDE_MARKER_RE.sub("", new_content)
    marker = f"(Manual Override Applied: Index {new_index}, Date {date_str})"
    lines = new_content.splitlines(True)
    out_lines = []
    for line in lines:
        out_lines.append(line)
        if _PREFERRED_INDEX_RE.search(line):
            out_lines.append(marker + "\n")
    return "".join(out_lines)


class TestOverrideMarker:
    SAMPLE_OCCUPIER = (
        "Some header\n"
        "(Preferred Index: 1)\n"
        "(Electron number: even)\n"
    )

    def test_marker_inserted_after_preferred_index(self):
        result = _apply_override_marker(self.SAMPLE_OCCUPIER, 3, "2026-04-14")
        assert "(Preferred Index: 3)" in result
        marker = "(Manual Override Applied: Index 3, Date 2026-04-14)"
        assert marker in result
        lines = result.splitlines()
        pi_line = next(i for i, l in enumerate(lines) if "Preferred Index" in l)
        assert lines[pi_line + 1] == marker

    def test_marker_replaced_on_second_override(self):
        content = (
            "Some header\n"
            "(Preferred Index: 2)\n"
            "(Manual Override Applied: Index 2, Date 2026-04-10)\n"
            "(Electron number: even)\n"
        )
        result = _apply_override_marker(content, 5, "2026-04-14")
        assert result.count("Manual Override") == 1
        assert "Index 5" in result
        assert "(Preferred Index: 5)" in result

    def test_electron_number_preserved(self):
        result = _apply_override_marker(self.SAMPLE_OCCUPIER, 2, "2026-04-14")
        assert "(Electron number: even)" in result
