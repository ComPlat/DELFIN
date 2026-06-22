"""Core dispatch logic for ``run_step()``."""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any, Optional

from delfin.common.logging import get_logger
from delfin.tools._types import StepResult, StepStatus

logger = get_logger(__name__)


def run_step(
    step_name: str,
    *,
    geometry: Optional[str | Path] = None,
    cores: int = 1,
    work_dir: Optional[str | Path] = None,
    **kwargs: Any,
) -> StepResult:
    """Run a single chemistry tool step and return a :class:`StepResult`.

    Parameters
    ----------
    step_name:
        Registered adapter name (e.g. ``"xtb_opt"``, ``"orca_sp"``).
    geometry:
        Path to input XYZ file.  Passed through to the adapter.
    cores:
        Number of CPU cores the tool may use.
    work_dir:
        Explicit working directory.  If *None* a temporary directory
        ``<step_name>_XXXX/`` is created under the current directory.
    **kwargs:
        Adapter-specific parameters (charge, mult, method, …).
    """
    from delfin.tools._registry import get as get_adapter, list_steps

    adapter = get_adapter(step_name)
    if adapter is None:
        available = ", ".join(sorted(list_steps().keys()))
        raise ValueError(
            f"Unknown step {step_name!r}. Available: {available}"
        )

    # Resolve geometry to absolute path
    geom_path: Optional[Path] = None
    if geometry is not None:
        geom_path = Path(geometry).resolve()
        if not geom_path.is_file():
            raise FileNotFoundError(f"Geometry file not found: {geom_path}")

    # Strip internal pipeline keys before validation/execution
    kwargs.pop("_prev_artifacts", None)

    # Validate before creating work dir
    adapter.validate_params(**kwargs)

    # Create isolated working directory
    if work_dir is not None:
        wdir = Path(work_dir).resolve()
        wdir.mkdir(parents=True, exist_ok=True)
    else:
        wdir = Path(
            tempfile.mkdtemp(prefix=f"{step_name}_", dir=Path.cwd())
        )

    logger.info("Running step '%s' in %s (cores=%d)", step_name, wdir, cores)

    try:
        result = adapter.execute(wdir, geometry=geom_path, cores=cores, **kwargs)
    except Exception as exc:
        logger.error("Step '%s' raised: %s", step_name, exc)
        result = StepResult(
            step_name=step_name,
            status=StepStatus.FAILED,
            work_dir=wdir,
            error=str(exc),
        )

    return result


def step_as_workflow_job(
    job_id: str,
    step_name: str,
    *,
    geometry: Optional[str | Path] = None,
    cores_min: int = 1,
    cores_optimal: int = 2,
    cores_max: int = 4,
    dependencies: Optional[set[str]] = None,
    result_callback: Optional[Any] = None,
    **kwargs: Any,
) -> Any:
    """Wrap a ``run_step()`` call as a :class:`WorkflowJob` for the scheduler.

    The returned job's ``work`` callable receives the allocated core count
    from the dynamic pool and forwards it to ``run_step(cores=...)``.
    """
    from delfin.tools._types import StepError
    from delfin.workflows.engine.classic import WorkflowJob

    def work(allocated_cores: int) -> None:
        result = run_step(step_name, geometry=geometry, cores=allocated_cores, **kwargs)
        if result_callback is not None:
            result_callback(result)
        if not result.ok:
            raise StepError(step_name, result.error or "failed")

    return WorkflowJob(
        job_id=job_id,
        work=work,
        description=f"step:{step_name}",
        dependencies=dependencies or set(),
        cores_min=cores_min,
        cores_optimal=cores_optimal,
        cores_max=cores_max,
    )
