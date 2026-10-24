"""Red tests for TM1 phase 2: turbomole contrib workflow through the WB2 path.

The user workflow: request -> spec -> resolve -> validate -> realism lint
-> dry run -> cost, all via delfin.agent.workflow_builder.run_checks /
delfin.tools.user_workflows.check_spec, saved as a contributed workflow
under delfin/workflows/contrib/ with a contract test.
"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from delfin.tools._registry import get as get_adapter


@pytest.fixture
def tm_fakes(tmp_path, monkeypatch):
    """Fakes for every Turbomole binary the workflow touches."""
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    bin_dir = tmp_path / "fakebin"
    bin_dir.mkdir()
    for name in ("ridft", "dscf", "jobex", "aoforce"):
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


class TestWorkflowRegistered:
    def test_w21_turbomole_workflow_is_discovered(self):
        """The contributed workflow must be registered and discoverable.

        Depends on the operator-built hook in delfin/workflows/registry.py
        importing delfin.workflows.contrib.w21_turbomole_workflow (the
        module's registration at import time is what makes this pass).
        """
        from delfin.workflows.registry import get as get_workflow

        wf = get_workflow("turbomole")
        assert wf is not None, (
            "contrib workflow 'turbomole' not discovered by the registry — "
            "the w21_turbomole_workflow module never registered")


class TestWorkflowSpecPassesWB2:
    """The workflow's spec passes the full WB2 check chain."""

    def _spec(self, work_dir: Path):
        return {
            "name": "turbomole_opt_freq",
            "kind": "template",
            "defaults": {"charge": "{charge}", "mult": "{mult}",
                         "method": "{method}", "basis": "{basis}"},
            "steps": [
                {"step": "turbomole_opt", "label": "tm_opt",
                 "kwargs": {"charge": "{charge}", "mult": "{mult}",
                            "method": "{method}", "basis": "{basis}"}},
                {"step": "turbomole_freq", "label": "tm_freq",
                 "kwargs": {"charge": "{charge}", "mult": "{mult}",
                            "method": "{method}", "basis": "{basis}"}},
            ],
            "template": True,
        }

    def test_spec_validates_through_wb2_chain(self):
        from delfin.tools import user_workflows as uw

        result = uw.check_spec(self._spec(Path(".")), geometry=True)
        assert result["ok"], result["errors"]

    def test_hostile_spec_fails_the_chain(self):
        from delfin.tools import user_workflows as uw

        spec = self._spec(Path("."))
        spec["steps"][0]["kwargs"]["mult"] = 0
        result = uw.check_spec(spec, geometry=True)
        assert not result["ok"], "a physically impossible spec must not pass"


class TestWorkflowRuns:
    def test_run_executes_opt_then_freq_with_fakes(self, tm_fakes):
        from delfin.workflows.registry import get as get_workflow

        wf = get_workflow("turbomole")
        assert wf is not None
        result = wf.run(config={}, geometry=_xyz(tm_fakes), charge=0,
                        mult=1, work_dir=tm_fakes)
        assert result is not None
        assert result["ok"], result.get("errors")

    def test_run_rejects_missing_geometry(self, tm_fakes):
        from delfin.workflows.registry import get as get_workflow

        wf = get_workflow("turbomole")
        with pytest.raises(Exception):
            wf.run(config={}, geometry=None, charge=0, work_dir=tm_fakes)
