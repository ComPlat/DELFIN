"""A step that declares it produces a capability actually leaves that
capability in its artifacts.

The framework's wiring trusts the ``produces`` declaration: when a step
declares ``produces=("X",)`` and its result carries no artifact under key X,
the next step that consumes X silently gets nothing and fails only at its
consumer, far from the mis-keyed producer.  ``run_step()`` therefore reports
such a step: its result carries ``_missing_declared_products`` naming the
step and the keys, and the log says the same.  Declared products that live
outside the artifacts dict by design are not reported: a ``geometry`` port is
satisfied by ``StepResult.geometry`` and ``qc_output`` by
``StepResult.output_file``.
"""

import pytest

from delfin.tools._registry import _REGISTRY, register
from delfin.tools._runner import run_step
from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus
from pathlib import Path


class _MisKeyed(StepAdapter):
    """Declares a product it never writes — the defect this check reports."""

    name = "test_miskeyed_producer"
    description = "declares special_key, produces nothing"
    category = "meta"
    produces_geometry = False
    produces = ("special_key",)

    def validate_params(self, **kwargs):
        return None

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        return StepResult(step_name=self.name, status=StepStatus.SUCCESS,
                          work_dir=Path(work_dir))


class _WellKeyed(StepAdapter):
    """Declares a product and writes it — must stay unreported."""

    name = "test_wellkeyed_producer"
    description = "declares special_key and writes it"
    category = "meta"
    produces_geometry = False
    produces = ("special_key",)

    def validate_params(self, **kwargs):
        return None

    def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
        path = Path(work_dir) / "special.out"
        path.write_text("ok\n", encoding="utf-8")
        return StepResult(step_name=self.name, status=StepStatus.SUCCESS,
                          work_dir=Path(work_dir),
                          artifacts={"special_key": path})


@pytest.fixture()
def _registered():
    register(_MisKeyed())
    register(_WellKeyed())
    yield
    _REGISTRY.pop("test_miskeyed_producer", None)
    _REGISTRY.pop("test_wellkeyed_producer", None)


def test_a_mis_keyed_producer_is_reported(tmp_path, _registered):
    result = run_step("test_miskeyed_producer", work_dir=tmp_path)
    assert result.ok
    assert result.data["_missing_declared_products"] == ("special_key",)


def test_a_well_keyed_producer_is_not_reported(tmp_path, _registered):
    result = run_step("test_wellkeyed_producer", work_dir=tmp_path)
    assert result.ok
    assert "_missing_declared_products" not in result.data


def test_a_geometry_producer_with_a_geometry_is_not_reported(tmp_path,
                                                             _registered):
    """``geometry`` lives in StepResult.geometry, not the artifacts dict."""
    from delfin.tools._types import StepStatus

    class GeometryOK(StepAdapter):
        name = "test_geometry_producer"
        produces_geometry = True

        def validate_params(self, **kwargs):
            return None

        def execute(self, work_dir, *, geometry=None, cores=1, **kwargs):
            path = Path(work_dir) / "out.xyz"
            path.write_text("2\nwater\nO 0 0 0\nH 0 0 1\nH 0 1 0\n",
                            encoding="utf-8")
            return StepResult(step_name=self.name, status=StepStatus.SUCCESS,
                              work_dir=Path(work_dir), geometry=path)

    register(GeometryOK())
    try:
        result = run_step("test_geometry_producer", work_dir=tmp_path)
        assert result.ok
        assert "_missing_declared_products" not in result.data
    finally:
        _REGISTRY.pop("test_geometry_producer", None)


def test_every_registered_adapter_declaring_a_product_satisfies_it(tmp_path):
    """The real registry, contract by contract: every adapter that declares a
    product and can run on this machine (no external binaries, no extra
    Python modules) is executed and its declared keys checked against its
    artifacts.  Adapters that need an external program are not run here —
    their check would fail on the missing program, not on a mis-keyed
    result — but each is named in the output so the coverage is visible."""
    from delfin.tools._registry import list_steps

    failures = []
    not_run = []
    # A tiny valid PDB so structure adapters with a required file parameter
    # actually run; same for the report adapters' project_dir.
    structure = tmp_path / "tiny.pdb"
    structure.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00"
        "           N\nEND\n",
        encoding="utf-8",
    )
    extra_kwargs: dict = {
        "biostructure_summary": {"structure_file": str(structure)},
        "report_json": {"project_dir": str(tmp_path)},
        "report_docx": {"project_dir": str(tmp_path)},
    }
    for name, adapter in sorted(list_steps().items()):
        # Only the built-in adapters; a test that registered a fixture
        # adapter earlier in the same process is not a contract of DELFIN's.
        if not adapter.__module__.startswith("delfin.tools.adapters"):
            continue
        contract = adapter.contract()
        declared = sorted(contract.produces)
        if not declared:
            continue
        if contract.requires_binaries or contract.requires_python:
            not_run.append(f"{name} (needs {sorted(contract.requires_binaries)}"
                           f"{sorted(contract.requires_python)})")
            continue
        try:
            result = run_step(name, work_dir=tmp_path / name,
                              **extra_kwargs.get(name, {}))
        except Exception as exc:  # a required param is missing — not the defect
            not_run.append(f"{name} (needs params: {exc})")
            continue
        if not result.ok:
            not_run.append(f"{name} (failed: {result.error})")
            continue
        missing = result.data.get("_missing_declared_products", ())
        if missing:
            failures.append(f"{name}: declares {declared}, "
                            f"artifacts carry {sorted(result.artifacts)} — "
                            f"missing {missing}")
        import shutil
        shutil.rmtree(tmp_path / name, ignore_errors=True)
    assert not failures, failures
    print(f"adapters declaring products, not run here: {len(not_run)}")
    for line in not_run:
        print(f"  {line}")
