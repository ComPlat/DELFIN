"""Tests for the delfin.tools adapter layer.

Covers: StepResult, StepStatus, StepError, StepAdapter, registry,
run_step() dispatch, Pipeline building/execution, PipelineTemplate
resolution, and defaults propagation.

All tests use a fake adapter — no chemistry tools needed.
"""

from __future__ import annotations

import textwrap
import time
from pathlib import Path
from typing import Any, Optional

import pytest

from delfin.tools._types import StepError, StepResult, StepStatus
from delfin.tools._base import StepAdapter


# ======================================================================
#  Fake adapter for testing (no external dependencies)
# ======================================================================

class FakeAdapter(StepAdapter):
    """Adapter that writes a dummy XYZ and returns configurable results."""

    name = "fake_step"
    description = "Fake adapter for testing"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if kwargs.get("fail_validation"):
            raise ValueError("validation failed on purpose")

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()

        if kwargs.get("fail"):
            return self._make_result(
                self.name, StepStatus.FAILED, work_dir, start,
                error="forced failure",
            )

        # Write a dummy output geometry
        out_xyz = work_dir / "output.xyz"
        if geometry and geometry.is_file():
            content = geometry.read_text()
        else:
            content = "1\nfake\nH  0.0  0.0  0.0\n"
        out_xyz.write_text(content)

        # Collect kwargs into data for inspection
        data = {"cores": cores, **{k: v for k, v in kwargs.items() if not k.startswith("_")}}

        # Fake artifacts
        artifacts = {}
        if kwargs.get("produce_gbw"):
            gbw = work_dir / "calc.gbw"
            gbw.write_bytes(b"fake-gbw")
            artifacts["gbw"] = gbw

        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            geometry=out_xyz,
            data=data,
            artifacts=artifacts,
        )


class FakeAnalysisAdapter(StepAdapter):
    """Adapter that does analysis only (no geometry output)."""

    name = "fake_analysis"
    description = "Fake analysis adapter"
    produces_geometry = False

    def execute(
        self,
        work_dir: Path,
        *,
        geometry: Optional[Path] = None,
        cores: int = 1,
        **kwargs: Any,
    ) -> StepResult:
        start = time.monotonic()
        return self._make_result(
            self.name, StepStatus.SUCCESS, work_dir, start,
            data={"analysis_result": 42.0},
        )


# ======================================================================
#  Fixtures
# ======================================================================

@pytest.fixture(autouse=True)
def _register_fake_adapters(monkeypatch):
    """Register fake adapters for every test, clean up after."""
    from delfin.tools import _registry

    # Save and restore registry state
    old_reg = dict(_registry._REGISTRY)
    old_disc = _registry._DISCOVERED

    _registry._DISCOVERED = True  # prevent real discovery
    _registry._REGISTRY.clear()
    _registry.register(FakeAdapter())
    _registry.register(FakeAnalysisAdapter())

    yield

    _registry._REGISTRY.clear()
    _registry._REGISTRY.update(old_reg)
    _registry._DISCOVERED = old_disc


@pytest.fixture
def xyz_file(tmp_path):
    """Create a minimal XYZ file."""
    p = tmp_path / "input.xyz"
    p.write_text("2\ntest\nH  0.0  0.0  0.0\nH  0.0  0.0  0.74\n")
    return p


# ======================================================================
#  StepResult / StepStatus / StepError
# ======================================================================

class TestStepResult:
    def test_ok_property(self):
        r = StepResult(step_name="x", status=StepStatus.SUCCESS)
        assert r.ok is True

        r2 = StepResult(step_name="x", status=StepStatus.FAILED)
        assert r2.ok is False

        r3 = StepResult(step_name="x", status=StepStatus.SKIPPED)
        assert r3.ok is False

    def test_require_success(self):
        r = StepResult(step_name="x", status=StepStatus.SUCCESS)
        assert r.require() is r

    def test_require_failure_raises(self):
        r = StepResult(step_name="x", status=StepStatus.FAILED, error="boom")
        with pytest.raises(StepError, match="boom"):
            r.require()

    def test_default_fields(self):
        r = StepResult(step_name="test", status=StepStatus.SUCCESS)
        assert r.geometry is None
        assert r.output_file is None
        assert r.work_dir is None
        assert r.data == {}
        assert r.artifacts == {}
        assert r.error is None
        assert r.elapsed_seconds == 0.0


class TestStepError:
    def test_message_format(self):
        e = StepError("orca_sp", "SCF did not converge")
        assert str(e) == "[orca_sp] SCF did not converge"
        assert e.step_name == "orca_sp"


# ======================================================================
#  StepAdapter base class
# ======================================================================

class TestStepAdapter:
    def test_copy_geometry_to_workdir(self, tmp_path, xyz_file):
        dest = StepAdapter._copy_geometry_to_workdir(xyz_file, tmp_path, "mol.xyz")
        assert dest == tmp_path / "mol.xyz"
        assert dest.read_text() == xyz_file.read_text()

    def test_make_result_timing(self, tmp_path):
        start = time.monotonic()
        time.sleep(0.01)
        r = StepAdapter._make_result("test", StepStatus.SUCCESS, tmp_path, start)
        assert r.elapsed_seconds >= 0.01
        assert r.step_name == "test"
        assert r.work_dir == tmp_path


# ======================================================================
#  Registry
# ======================================================================

class TestRegistry:
    def test_list_steps(self):
        from delfin.tools._registry import list_steps
        steps = list_steps()
        assert "fake_step" in steps
        assert "fake_analysis" in steps

    def test_get_existing(self):
        from delfin.tools._registry import get
        adapter = get("fake_step")
        assert adapter is not None
        assert adapter.name == "fake_step"

    def test_get_nonexistent(self):
        from delfin.tools._registry import get
        assert get("nonexistent_step") is None

    def test_register_no_name_raises(self):
        from delfin.tools._registry import register

        class BadAdapter(StepAdapter):
            name = ""
            description = "no name"
            def execute(self, work_dir, **kw):
                pass

        with pytest.raises(ValueError, match="no 'name'"):
            register(BadAdapter())


# ======================================================================
#  run_step() dispatch
# ======================================================================

class TestRunStep:
    def test_basic_execution(self, tmp_path, xyz_file):
        from delfin.tools._runner import run_step
        r = run_step("fake_step", geometry=xyz_file, work_dir=tmp_path / "work")
        assert r.ok
        assert r.geometry is not None
        assert r.geometry.is_file()
        assert r.work_dir == tmp_path / "work"

    def test_creates_workdir_if_none(self, tmp_path, xyz_file, monkeypatch):
        import os
        monkeypatch.chdir(tmp_path)
        from delfin.tools._runner import run_step
        r = run_step("fake_step", geometry=xyz_file)
        assert r.ok
        assert r.work_dir is not None
        assert r.work_dir.is_dir()
        assert r.work_dir.parent == tmp_path

    def test_unknown_step_raises(self):
        from delfin.tools._runner import run_step
        with pytest.raises(ValueError, match="Unknown step 'nope'"):
            run_step("nope")

    def test_missing_geometry_raises(self, tmp_path):
        from delfin.tools._runner import run_step
        with pytest.raises(FileNotFoundError):
            run_step("fake_step", geometry="/no/such/file.xyz", work_dir=tmp_path / "w")

    def test_validation_failure(self, tmp_path, xyz_file):
        from delfin.tools._runner import run_step
        with pytest.raises(ValueError, match="validation failed"):
            run_step("fake_step", geometry=xyz_file, work_dir=tmp_path / "w",
                     fail_validation=True)

    def test_adapter_failure_returns_failed_result(self, tmp_path, xyz_file):
        from delfin.tools._runner import run_step
        r = run_step("fake_step", geometry=xyz_file, work_dir=tmp_path / "w", fail=True)
        assert not r.ok
        assert r.status == StepStatus.FAILED
        assert r.error == "forced failure"

    def test_cores_passed_through(self, tmp_path, xyz_file):
        from delfin.tools._runner import run_step
        r = run_step("fake_step", geometry=xyz_file, cores=8, work_dir=tmp_path / "w")
        assert r.data["cores"] == 8

    def test_kwargs_passed_through(self, tmp_path, xyz_file):
        from delfin.tools._runner import run_step
        r = run_step("fake_step", geometry=xyz_file, work_dir=tmp_path / "w",
                     charge=2, method="B3LYP")
        assert r.data["charge"] == 2
        assert r.data["method"] == "B3LYP"

    def test_prev_artifacts_stripped(self, tmp_path, xyz_file):
        """_prev_artifacts should not leak to the adapter."""
        from delfin.tools._runner import run_step
        r = run_step("fake_step", geometry=xyz_file, work_dir=tmp_path / "w",
                     _prev_artifacts={"gbw": Path("/fake.gbw")})
        assert r.ok
        assert "_prev_artifacts" not in r.data


# ======================================================================
#  Pipeline building
# ======================================================================

class TestPipelineBuilding:
    def test_add_returns_self(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test")
        ret = p.add("fake_step", smiles="C")
        assert ret is p

    def test_trunk_order(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test")
        p.add("fake_step", label="step1")
        p.add("fake_analysis", label="step2")
        assert len(p._trunk) == 2
        assert p._trunk[0].step_name == "fake_step"
        assert p._trunk[1].step_name == "fake_analysis"

    def test_branch_creates_child(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test")
        b = p.branch("ox")
        assert "ox" in p._branches
        assert isinstance(b, Pipeline)
        assert b.name == "test/ox"

    def test_branch_same_name_returns_existing(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test")
        b1 = p.branch("ox")
        b2 = p.branch("ox")
        assert b1 is b2

    def test_repr(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test")
        p.add("fake_step")
        p.branch("b")
        assert "test" in repr(p)
        assert "1 steps" in repr(p)
        assert "1 branches" in repr(p)


# ======================================================================
#  Pipeline defaults
# ======================================================================

class TestPipelineDefaults:
    def test_defaults_merge_into_steps(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", defaults={"charge": 2, "method": "B3LYP"})
        p.add("fake_step", basis="def2-SVP")
        assert p._trunk[0].kwargs == {"charge": 2, "method": "B3LYP", "basis": "def2-SVP"}

    def test_step_kwargs_override_defaults(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", defaults={"charge": 2, "method": "B3LYP"})
        p.add("fake_step", charge=0)
        assert p._trunk[0].kwargs["charge"] == 0
        assert p._trunk[0].kwargs["method"] == "B3LYP"

    def test_defaults_propagate_to_branches(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", defaults={"charge": 2, "method": "B3LYP"})
        b = p.branch("ox")
        assert b._defaults == {"charge": 2, "method": "B3LYP"}

    def test_branch_steps_get_defaults(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", defaults={"charge": 2})
        b = p.branch("ox")
        b.add("fake_step", mult=2)
        assert b._trunk[0].kwargs == {"charge": 2, "mult": 2}

    def test_empty_defaults(self):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test")
        assert p._defaults == {}
        p.add("fake_step", charge=0)
        assert p._trunk[0].kwargs == {"charge": 0}


# ======================================================================
#  Pipeline execution
# ======================================================================

class TestPipelineExecution:
    def test_sequential_run(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add("fake_analysis")
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert len(result.results) == 2
        assert result.results[0].ok
        assert result.results[1].ok

    def test_geometry_propagation(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")  # produces geometry
        p.add("fake_step")  # should receive geometry from step 0
        result = p.run(geometry=xyz_file)
        assert result.ok
        # Step 1 should have used step 0's output geometry
        assert result.results[1].geometry is not None
        assert result.results[1].geometry.is_file()

    def test_stop_on_failure(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add("fake_step", fail=True)  # will fail
        p.add("fake_step")  # should not run
        result = p.run(geometry=xyz_file, stop_on_failure=True)
        assert not result.ok
        assert len(result.results) == 2  # stopped after failure

    def test_continue_on_failure(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add("fake_step", fail=True)
        p.add("fake_step")
        result = p.run(geometry=xyz_file, stop_on_failure=False)
        assert not result.ok
        assert len(result.results) == 3  # all ran

    def test_branch_execution(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")

        b1 = p.branch("alpha")
        b1.add("fake_step")
        b1.add("fake_analysis")

        b2 = p.branch("beta")
        b2.add("fake_step")

        result = p.run(geometry=xyz_file, cores=4)
        assert result.ok
        assert len(result.branch_results) == 2
        assert result.branch("alpha").ok
        assert result.branch("beta").ok
        assert len(result.branch("alpha").results) == 2

    def test_artifact_propagation(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step", produce_gbw=True)  # produces gbw artifact
        p.add("fake_step")  # should receive _prev_artifacts
        result = p.run(geometry=xyz_file)
        assert result.ok
        # The second step should have had _prev_artifacts injected (stripped before data)

    def test_on_step_callback(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        collected = []
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add("fake_analysis")
        p.on_step(lambda r: collected.append(r.step_name))
        p.run(geometry=xyz_file)
        assert collected == ["fake_step", "fake_analysis"]

    def test_defaults_applied_at_runtime(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path, defaults={"charge": 2})
        p.add("fake_step")
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert result.results[0].data["charge"] == 2


# ======================================================================
#  PipelineResult
# ======================================================================

class TestPipelineResult:
    def test_ok_all_success(self):
        from delfin.tools.pipeline import PipelineResult
        r = PipelineResult(
            name="test",
            results=[
                StepResult(step_name="a", status=StepStatus.SUCCESS),
                StepResult(step_name="b", status=StepStatus.SUCCESS),
            ],
            branch_results={},
        )
        assert r.ok

    def test_ok_false_on_failure(self):
        from delfin.tools.pipeline import PipelineResult
        r = PipelineResult(
            name="test",
            results=[
                StepResult(step_name="a", status=StepStatus.SUCCESS),
                StepResult(step_name="b", status=StepStatus.FAILED),
            ],
            branch_results={},
        )
        assert not r.ok

    def test_ok_false_on_branch_failure(self):
        from delfin.tools.pipeline import PipelineResult
        r = PipelineResult(
            name="test",
            results=[StepResult(step_name="a", status=StepStatus.SUCCESS)],
            branch_results={
                "b1": PipelineResult(
                    name="b1",
                    results=[StepResult(step_name="x", status=StepStatus.FAILED)],
                    branch_results={},
                )
            },
        )
        assert not r.ok

    def test_last(self):
        from delfin.tools.pipeline import PipelineResult
        r = PipelineResult(
            name="test",
            results=[
                StepResult(step_name="a", status=StepStatus.SUCCESS),
                StepResult(step_name="b", status=StepStatus.SUCCESS),
            ],
            branch_results={},
        )
        assert r.last.step_name == "b"

    def test_last_empty(self):
        from delfin.tools.pipeline import PipelineResult
        r = PipelineResult(name="test", results=[], branch_results={})
        assert r.last is None

    def test_all_results_flat(self):
        from delfin.tools.pipeline import PipelineResult
        r = PipelineResult(
            name="test",
            results=[StepResult(step_name="trunk", status=StepStatus.SUCCESS)],
            branch_results={
                "b1": PipelineResult(
                    name="b1",
                    results=[StepResult(step_name="b1_step", status=StepStatus.SUCCESS)],
                    branch_results={},
                ),
            },
        )
        flat = r.all_results
        assert len(flat) == 2
        assert flat[0].step_name == "trunk"
        assert flat[1].step_name == "b1_step"

    def test_summary_format(self):
        from delfin.tools.pipeline import PipelineResult
        r = PipelineResult(
            name="test",
            results=[
                StepResult(step_name="a", status=StepStatus.SUCCESS, elapsed_seconds=1.5),
                StepResult(step_name="b", status=StepStatus.FAILED, error="oops", elapsed_seconds=0.3),
            ],
            branch_results={},
        )
        s = r.summary()
        assert "FAILED" in s
        assert "oops" in s
        assert "1.5s" in s


# ======================================================================
#  PipelineTemplate — value resolution
# ======================================================================

class TestResolveValue:
    def test_simple_placeholder(self):
        from delfin.tools.pipeline import _resolve_value
        assert _resolve_value("{charge}", {"charge": 2}) == 2

    def test_arithmetic_plus(self):
        from delfin.tools.pipeline import _resolve_value
        assert _resolve_value("{charge}+1", {"charge": 2}) == 3

    def test_arithmetic_minus(self):
        from delfin.tools.pipeline import _resolve_value
        assert _resolve_value("{charge}-1", {"charge": 2}) == 1

    def test_string_interpolation(self):
        from delfin.tools.pipeline import _resolve_value
        assert _resolve_value("CPCM({solvent})", {"solvent": "water"}) == "CPCM(water)"

    def test_non_string_passthrough(self):
        from delfin.tools.pipeline import _resolve_value
        assert _resolve_value(42, {"x": 1}) == 42
        assert _resolve_value(None, {"x": 1}) is None

    def test_no_placeholder_passthrough(self):
        from delfin.tools.pipeline import _resolve_value
        assert _resolve_value("B3LYP", {"method": "PBE0"}) == "B3LYP"

    def test_unresolved_placeholder(self):
        from delfin.tools.pipeline import _resolve_value
        assert _resolve_value("{unknown}", {"charge": 2}) == "{unknown}"


class TestResolveKwargs:
    def test_flat_dict(self):
        from delfin.tools.pipeline import _resolve_kwargs
        r = _resolve_kwargs(
            {"charge": "{charge}", "method": "B3LYP", "mult": "{mult}"},
            {"charge": 2, "mult": 3},
        )
        assert r == {"charge": 2, "method": "B3LYP", "mult": 3}

    def test_nested_dict(self):
        from delfin.tools.pipeline import _resolve_kwargs
        r = _resolve_kwargs(
            {"metal_basis": {"Fe": "{fe_basis}", "Cu": "{cu_basis}"}},
            {"fe_basis": "def2-TZVP", "cu_basis": "def2-SVP"},
        )
        assert r == {"metal_basis": {"Fe": "def2-TZVP", "Cu": "def2-SVP"}}

    def test_list_values(self):
        from delfin.tools.pipeline import _resolve_kwargs
        r = _resolve_kwargs(
            {"items": ["{a}", "{b}", "literal"]},
            {"a": 1, "b": 2},
        )
        assert r == {"items": [1, 2, "literal"]}


# ======================================================================
#  PipelineTemplate building and execution
# ======================================================================

class TestPipelineTemplate:
    def test_build_resolves_placeholders(self):
        from delfin.tools.pipeline import PipelineTemplate
        t = PipelineTemplate("tpl")
        t.add("fake_step", charge="{charge}", method="{method}")
        pipe = t.build(charge=2, method="B3LYP")
        assert pipe._trunk[0].kwargs["charge"] == 2
        assert pipe._trunk[0].kwargs["method"] == "B3LYP"

    def test_build_with_branches(self):
        from delfin.tools.pipeline import PipelineTemplate
        t = PipelineTemplate("tpl")
        t.add("fake_step", charge="{charge}")
        b = t.branch("ox")
        b.add("fake_step", charge="{charge}+1")
        pipe = t.build(charge=2)
        assert pipe._trunk[0].kwargs["charge"] == 2
        assert pipe._branches["ox"]._trunk[0].kwargs["charge"] == 3

    def test_template_defaults(self):
        from delfin.tools.pipeline import PipelineTemplate
        t = PipelineTemplate("tpl", defaults={"method": "{method}", "basis": "{basis}"})
        t.add("fake_step", charge="{charge}")
        # defaults should be merged in
        assert t._trunk[0].kwargs == {"method": "{method}", "basis": "{basis}", "charge": "{charge}"}
        pipe = t.build(method="B3LYP", basis="def2-SVP", charge=0)
        assert pipe._trunk[0].kwargs == {"method": "B3LYP", "basis": "def2-SVP", "charge": 0}

    def test_template_defaults_propagate_to_branches(self):
        from delfin.tools.pipeline import PipelineTemplate
        t = PipelineTemplate("tpl", defaults={"method": "{method}"})
        b = t.branch("b1")
        b.add("fake_step", charge="{charge}")
        # branch inherits parent defaults
        assert b._trunk[0].kwargs == {"method": "{method}", "charge": "{charge}"}

    def test_run_builds_and_executes(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import PipelineTemplate
        t = PipelineTemplate("tpl")
        t.add("fake_step", charge="{charge}")
        t.add("fake_analysis")
        result = t.run(charge=0, cores=2, geometry=xyz_file, work_dir=tmp_path)
        assert result.ok
        assert len(result.results) == 2
        assert result.results[0].data["charge"] == 0

    def test_repr(self):
        from delfin.tools.pipeline import PipelineTemplate
        t = PipelineTemplate("tpl")
        t.add("fake_step")
        t.branch("b")
        r = repr(t)
        assert "tpl" in r
        assert "1 steps" in r
        assert "1 branches" in r


# ======================================================================
#  Public API imports
# ======================================================================

class TestPublicAPI:
    def test_all_exports_importable(self):
        from delfin.tools import (
            run_step,
            step_as_workflow_job,
            Pipeline,
            PipelineTemplate,
            PipelineResult,
            StepResult,
            StepStatus,
            StepError,
            list_steps,
            register,
        )
        # Just verify they're all callable/classes
        assert callable(run_step)
        assert callable(step_as_workflow_job)
        assert callable(list_steps)

    def test_list_steps_returns_dict(self):
        from delfin.tools import list_steps
        steps = list_steps()
        assert isinstance(steps, dict)


# ======================================================================
#  Pipeline advanced: add_if (conditional)
# ======================================================================

class TestPipelineConditional:
    def test_add_if_true(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add_if(lambda results, last: True, "fake_analysis")
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert len(result.results) == 2
        assert result.results[1].step_name == "fake_analysis"

    def test_add_if_false_skips(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add_if(lambda results, last: False, "fake_analysis")
        result = p.run(geometry=xyz_file)
        # SKIPPED counts as not-ok, but the pipeline ran to completion
        assert len(result.results) == 2
        assert result.results[0].ok
        assert result.results[1].status == StepStatus.SKIPPED

    def test_add_if_uses_last_result(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")  # produces data with cores=1
        p.add_if(
            lambda results, last: last.data.get("cores") == 1,
            "fake_analysis",
        )
        result = p.run(geometry=xyz_file)
        assert len(result.results) == 2
        assert result.results[1].step_name == "fake_analysis"
        assert result.results[1].status == StepStatus.SUCCESS


# ======================================================================
#  Pipeline advanced: add_loop
# ======================================================================

class TestPipelineLoop:
    def test_loop_stops_on_condition(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        call_count = [0]

        def count_and_stop(result, iteration):
            call_count[0] += 1
            return iteration >= 2  # stop after 3 iterations (0, 1, 2)

        p = Pipeline("test", base_dir=tmp_path)
        p.add_loop("fake_step", until=count_and_stop, max_iter=10)
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert call_count[0] == 3  # iterations 0, 1, 2

    def test_loop_respects_max_iter(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add_loop("fake_step", until=lambda r, i: False, max_iter=3)
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert len(result.results) == 1  # loop produces one final result

    def test_loop_propagates_geometry(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add_loop("fake_step", until=lambda r, i: i >= 1, max_iter=5)
        p.add("fake_analysis")  # should get geometry from loop
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert len(result.results) == 2


# ======================================================================
#  Pipeline advanced: add_transform
# ======================================================================

class TestPipelineTransform:
    def test_transform_modifies_kwargs(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add_transform(
            "fake_step",
            lambda kw, results, last: {**kw, "charge": 42},
        )
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert result.results[1].data["charge"] == 42

    def test_transform_uses_previous_data(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")  # data has cores=1
        p.add_transform(
            "fake_step",
            lambda kw, results, last: {**kw, "charge": last.data["cores"] * 10},
        )
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert result.results[1].data["charge"] == 10

    def test_transform_error_fails_step(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add_transform(
            "fake_step",
            lambda kw, results, last: 1/0,  # will raise
        )
        result = p.run(geometry=xyz_file, stop_on_failure=True)
        assert not result.ok
        assert result.results[1].status == StepStatus.FAILED


# ======================================================================
#  Pipeline advanced: add_compute
# ======================================================================

class TestPipelineCompute:
    def test_compute_stores_data(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add_compute(
            lambda results, last, wd: {"answer": 42, "pi": 3.14},
            label="my_calc",
        )
        result = p.run(geometry=xyz_file)
        assert result.ok
        assert result.results[1].data["answer"] == 42
        assert result.results[1].data["pi"] == 3.14

    def test_compute_has_access_to_results(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add("fake_step")
        p.add_compute(
            lambda results, last, wd: {"n_prev": len(results)},
        )
        result = p.run(geometry=xyz_file)
        assert result.results[1].data["n_prev"] == 1

    def test_compute_error_fails(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add_compute(lambda r, l, w: 1/0, label="bad_calc")
        result = p.run(geometry=xyz_file, stop_on_failure=True)
        assert not result.ok
        assert "division by zero" in result.results[0].error


# ======================================================================
#  Pipeline advanced: add_fan_out
# ======================================================================

class TestPipelineFanOut:
    def test_fan_out_runs_on_multiple_geometries(self, tmp_path):
        from delfin.tools.pipeline import Pipeline
        # Create multiple xyz files
        geoms = []
        for i in range(3):
            g = tmp_path / f"mol_{i}.xyz"
            g.write_text(f"1\nmol {i}\nH  0.0  0.0  {float(i)}\n")
            geoms.append(g)

        p = Pipeline("test", base_dir=tmp_path / "work")
        p.add_fan_out(
            "fake_step",
            geometries_from=lambda results, last: geoms,
            label="fan",
        )
        result = p.run()
        assert result.ok
        assert result.results[0].data["fan_out_count"] == 3
        assert result.results[0].data["fan_out_ok"] == 3

    def test_fan_out_empty_skips(self, tmp_path, xyz_file):
        from delfin.tools.pipeline import Pipeline
        p = Pipeline("test", base_dir=tmp_path)
        p.add_fan_out(
            "fake_step",
            geometries_from=lambda results, last: [],
            label="empty_fan",
        )
        result = p.run(geometry=xyz_file)
        assert len(result.results) == 1
        assert result.results[0].status == StepStatus.SKIPPED
