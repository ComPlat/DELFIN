"""Adversarial tests for TM1 phase 2: the turbomole contrib workflow.

Targets the physical-realism reporting duty of
delfin/workflows/contrib/w21_turbomole_workflow.py::_count_imaginary:
a non-minimum (imaginary frequencies) must be counted and surfaced in
the result, never silently dropped. The builder's own tests never
exercise _count_imaginary (their fake aoforce writes no aoforce.out).
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.tools._registry import get as get_adapter


@pytest.fixture
def tm_fakes(tmp_path, monkeypatch):
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    for name in ("ridft", "dscf", "jobex"):
        f = bin_dir / name
        f.write_text("#!/bin/sh\nexit 0\n")
        f.chmod(f.stat().st_mode | stat.S_IEXEC)
    aof = bin_dir / "aoforce"
    aof.write_text(f"#!/bin/sh\necho '$hessian' > {work_dir}/hessian\nexit 0\n")
    aof.chmod(aof.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep + os.environ.get("PATH", ""))
    return work_dir


def _xyz(work_dir: Path) -> Path:
    p = work_dir / "mol.xyz"
    p.write_text("3\nwater\nO 0.0 0.0 0.0\nH 0.9 0.0 0.0\nH -0.3 0.8 0.0\n")
    return p


class TestImaginaryReporting:
    """A non-minimum must be reported, not hidden."""

    def test_negative_frequency_counted_and_warned(self, tm_fakes):
        """aoforce.out with one negative frequency -> count 1 + warning."""
        from delfin.workflows.contrib.w21_turbomole_workflow import (
            TurbomoleWorkflow,
        )

        (tm_fakes / "aoforce.out").write_text(
            " n(t-i)     frequency (cm-1)\n"
            "   1        -123.45\n"
            "   2         1650.00\n"
        )
        count = TurbomoleWorkflow._count_imaginary(tm_fakes)
        assert count == 1, (
            "one negative frequency in aoforce.out must be counted as 1 "
            f"imaginary mode, got {count!r}"
        )

    def test_positive_only_counts_zero(self, tm_fakes):
        from delfin.workflows.contrib.w21_turbomole_workflow import (
            TurbomoleWorkflow,
        )

        (tm_fakes / "aoforce.out").write_text(
            "   1        1650.00\n   2        3400.00\n"
        )
        assert TurbomoleWorkflow._count_imaginary(tm_fakes) == 0

    def test_run_surfaces_warning_on_non_minimum(self, tm_fakes):
        """End to end: a non-minimum run must carry the warning."""
        from delfin.workflows.registry import get as get_workflow

        (tm_fakes / "aoforce.out").write_text("   1        -99.00\n")
        wf = get_workflow("turbomole")
        assert wf is not None, "workflow not registered on this branch"
        result = wf.run(config={}, geometry=_xyz(tm_fakes), charge=0,
                        work_dir=tm_fakes)
        assert result.get("imaginary_frequencies") == 1, result
        assert any("imaginary" in w.lower() for w in result.get("warnings", [])), (
            "a non-minimum must be surfaced as a warning, not hidden"
        )

    def test_no_log_reported_as_unknown_not_zero(self, tm_fakes):
        """No aoforce.out -> the count is unknown (None), never claimed 0."""
        from delfin.workflows.contrib.w21_turbomole_workflow import (
            TurbomoleWorkflow,
        )

        assert TurbomoleWorkflow._count_imaginary(tm_fakes) is None
