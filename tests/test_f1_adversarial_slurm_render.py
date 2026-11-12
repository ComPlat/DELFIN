"""Adversarial tests on slurm_tests.build_job — the render model Target 2
(durable host) must copy (F1 reviewer s2).

These pin the guarantees an inherited renderer must keep and probe the
allow-lists for bypasses (newline smuggle, leading dash, .., empty
segments). GREEN here = the model is sound; a durable-host renderer that
drops one of these is a red finding against that hash.
"""

from __future__ import annotations

import pytest

from delfin.agent import slurm_tests


def _render(**kw) -> str:
    defaults = dict(
        repo="/repo", ref="main", test_paths=["tests/test_x.py"],
        partition="cpu", minutes=25,
    )
    defaults.update(kw)
    return slurm_tests.build_job(**defaults)


class TestExportRule:
    def test_export_is_exactly_ALL(self):
        text = _render()
        assert "#SBATCH --export=ALL\n" in text
        assert "--export=ALL," not in text, (
            "an explicit --export list must never be rendered")

    def test_partition_newline_cannot_smuggle_an_export_line(self):
        with pytest.raises(ValueError):
            _render(partition="cpu\n#SBATCH --export=ALL,EVIL=1")

    def test_ref_newline_cannot_smuggle_a_command_line(self):
        with pytest.raises(ValueError):
            _render(ref="main\ntouch /tmp/x")

    def test_leading_dash_refused_everywhere(self):
        for field in ("ref", "partition", "mem", "job_name"):
            with pytest.raises(ValueError):
                _render(**{field: "-evil"})


class TestValidation:
    def test_dotdot_in_test_path_refused(self):
        with pytest.raises(ValueError):
            _render(test_paths=["tests/../etc/passwd"])

    def test_ref_dotdot_refused(self):
        with pytest.raises(ValueError):
            _render(ref="a/../b")

    def test_empty_segment_partition_refused(self):
        with pytest.raises(ValueError):
            _render(partition="cpu,,gpu")

    def test_empty_test_path_list_refused(self):
        with pytest.raises(ValueError):
            _render(test_paths=[])

    def test_minutes_zero_or_negative_refused(self):
        for m in (0, -5):
            with pytest.raises(ValueError):
                _render(minutes=m)

    def test_mem_suffix_only(self):
        with pytest.raises(ValueError):
            _render(mem="8GB; rm -rf /")

    def test_job_name_with_space_refused(self):
        with pytest.raises(ValueError):
            _render(job_name="my job")


class TestRenderedShape:
    def test_rendered_script_has_set_guard_and_signals(self):
        text = _render()
        assert text.startswith("#!/bin/bash")
        assert "set -euo pipefail" in text

    def test_same_args_same_script(self):
        assert _render() == _render()

    def test_job_name_reaches_the_header(self):
        assert "#SBATCH --job-name=delfin-agent-tests" in _render()
