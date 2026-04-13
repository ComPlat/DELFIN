"""Pipeline builder for chaining tool steps into workflows.

Provides a declarative way to build sequential and branching workflows
from ``run_step()`` calls, with automatic geometry propagation.

Sequential example::

    from delfin.tools import Pipeline

    pipe = Pipeline("basic_opt")
    pipe.add("smiles_to_xyz", smiles="CCO")
    pipe.add("xtb_opt", charge=0)
    pipe.add("orca_sp", charge=0, method="B3LYP", basis="def2-SVP")
    results = pipe.run(cores=4)
    print(results.last.data["energy_Eh"])

Branching example (parallel oxidation + reduction)::

    pipe = Pipeline("ox_red")
    pipe.add("smiles_to_xyz", smiles="CCO")
    pipe.add("xtb_opt", charge=0)

    ox = pipe.branch("oxidation")
    ox.add("orca_opt", charge=1, mult=2, method="B3LYP", basis="def2-SVP")
    ox.add("orca_freq", charge=1, mult=2, method="B3LYP", basis="def2-SVP")

    red = pipe.branch("reduction")
    red.add("orca_opt", charge=-1, method="B3LYP", basis="def2-SVP")
    red.add("orca_freq", charge=-1, method="B3LYP", basis="def2-SVP")

    results = pipe.run(cores=16)  # branches run in parallel

Parametric template — define once, run on any system::

    from delfin.tools import PipelineTemplate

    # Define a reusable redox workflow template
    redox = PipelineTemplate("redox_potential")
    redox.add("smiles_to_xyz", smiles="{smiles}")
    redox.add("xtb_opt", charge="{charge}")
    redox.add("orca_opt", charge="{charge}", mult="{mult}",
              method="{method}", basis="{basis}",
              ri="RIJCOSX", aux_basis="def2/J", dispersion="D4",
              solvent="{solvent}", solvent_model="CPCM")

    ox = redox.branch("oxidation")
    ox.add("orca_opt", charge="{charge}+1", mult="{mult_ox}",
           method="{method}", basis="{basis}")
    ox.add("orca_freq", charge="{charge}+1", mult="{mult_ox}",
           method="{method}", basis="{basis}")

    red = redox.branch("reduction")
    red.add("orca_opt", charge="{charge}-1", mult="{mult_red}",
            method="{method}", basis="{basis}")
    red.add("orca_freq", charge="{charge}-1", mult="{mult_red}",
            method="{method}", basis="{basis}")

    # Apply to any system
    results = redox.run(
        smiles="[Fe+2](N)(N)(N)(N)(N)N",
        charge=2, mult=5, mult_ox=4, mult_red=6,
        method="B3LYP", basis="def2-TZVP", solvent="water",
        cores=16,
    )

Pipeline-level defaults (like CONTROL.txt for a pipeline)::

    pipe = Pipeline("metal_opt", defaults={
        "charge": 2, "mult": 3,
        "method": "B3LYP", "basis": "def2-SVP",
        "ri": "RIJCOSX", "aux_basis": "def2/J",
        "solvent": "water", "solvent_model": "CPCM",
    })
    pipe.add("xtb_opt")              # inherits charge, mult
    pipe.add("orca_opt")             # inherits all defaults
    pipe.add("orca_freq")            # inherits all defaults
    pipe.add("orca_sp", basis="def2-TZVP")  # override basis only
    results = pipe.run(cores=8, geometry="input.xyz")

Conditional steps — only run if condition is met::

    pipe = Pipeline("smart_opt")
    pipe.add("xtb_opt", charge=0)
    pipe.add("orca_freq", charge=0)
    pipe.add_if(
        lambda results, last: last.data.get("has_imaginary", False),
        "imag_fix", charge=0,
    )

Loop until convergence::

    pipe = Pipeline("imag_loop")
    pipe.add("orca_opt", charge=0)
    pipe.add_loop(
        "imag_fix", charge=0,
        until=lambda result, i: result.data.get("n_imaginary", 0) == 0,
        max_iter=5,
    )

Dynamic kwargs from previous results::

    pipe = Pipeline("adaptive")
    pipe.add("smiles_to_xyz", smiles="CCO")
    pipe.add_transform(
        "orca_sp",
        lambda kw, results, last: {**kw, "mult": 1 if last.data["n_atoms"] % 2 == 0 else 2},
        charge=0, method="B3LYP",
    )

Post-processing — compute derived values::

    def calc_redox(results, last, work_dir):
        E0 = results[1].data["energy_Eh"]
        E_ox = results[2].data["energy_Eh"]
        return {"redox_V": -((E0 - E_ox) * 2625.5) / 96485 - 4.28}

    pipe.add_compute(calc_redox, label="redox_potential")

Fan-out — run on multiple geometries in parallel::

    pipe.add("crest_conformers", charge=0)
    pipe.add_fan_out(
        "orca_sp",
        geometries_from=lambda results, last: list(last.work_dir.glob("*.xyz")),
        charge=0, method="B3LYP",
        label="screen_conformers",
    )

Scheduler integration::

    pipe = Pipeline("scheduled")
    pipe.add("smiles_to_xyz", smiles="C")
    pipe.add("xtb_opt", charge=0)
    pipe.run_scheduled(cores=8)  # uses DynamicCorePool
"""

from __future__ import annotations

import threading
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from delfin.common.logging import get_logger
from delfin.tools._types import StepResult, StepStatus, StepError

logger = get_logger(__name__)


def _resolve_cores(cores: int | str) -> int:
    """Resolve ``cores`` argument — supports ``"auto"`` for cluster detection."""
    if isinstance(cores, str) and cores.lower() == "auto":
        from delfin.cluster_utils import detect_cluster_environment
        info = detect_cluster_environment()
        detected = info.get("cpus_available") or 1
        logger.info(
            "Auto-detected %d cores (%s scheduler)",
            detected, info.get("scheduler", "unknown"),
        )
        return detected
    return int(cores)


@dataclass
class _StepSpec:
    """Specification for a single step in a pipeline."""
    step_name: str
    kwargs: Dict[str, Any]
    geometry_override: Optional[Path] = None  # explicit input geometry
    label: str = ""  # human-readable label (defaults to step_name)
    # --- Advanced flow control ---
    condition: Optional[Callable] = None       # (results, last) -> bool
    transform: Optional[Callable] = None       # (kwargs, results, last) -> kwargs
    loop_until: Optional[Callable] = None      # (result, iteration) -> bool (True=stop)
    loop_max: int = 1                          # max iterations (1 = no loop)
    fan_out_fn: Optional[Callable] = None      # (results, last) -> list[Path]
    compute_fn: Optional[Callable] = None      # (results, last, work_dir) -> dict
    retry_max: int = 1                         # max attempts (1 = no retry)
    retry_delay: float = 0.0                   # seconds between retries
    map_items: Optional[Callable] = None       # (results, last) -> list[dict]
    checkpoint: bool = False                   # save checkpoint after this step
    sub_pipeline: Optional[Any] = None         # Pipeline to run as a step
    sub_pipeline_kwargs: Optional[Dict] = None # kwargs for sub-pipeline run()
    reactive_fn: Optional[Callable] = None     # (results, last, pipeline) -> None


class PipelineInserter:
    """Helper passed to ``add_reactive()`` callbacks to insert steps dynamically.

    Steps added via the inserter are collected and executed immediately
    after the reactive callback returns.
    """

    def __init__(self, defaults: Dict[str, Any]):
        self._steps: List[_StepSpec] = []
        self._defaults = defaults

    def add(self, step_name: str, *, label: str = "", geometry: Optional[str | Path] = None,
            **kwargs: Any) -> "PipelineInserter":
        """Insert a normal step."""
        merged = {**self._defaults, **kwargs}
        self._steps.append(_StepSpec(
            step_name=step_name,
            kwargs=merged,
            label=label or step_name,
            geometry_override=Path(geometry) if geometry else None,
        ))
        return self

    def add_sub_pipeline(self, sub: "Pipeline", *, label: str = "",
                         geometry: Optional[str | Path] = None,
                         **run_kwargs: Any) -> "PipelineInserter":
        """Insert a sub-pipeline."""
        self._steps.append(_StepSpec(
            step_name="_sub_pipeline",
            kwargs={},
            label=label or f"sub:{sub.name}",
            geometry_override=Path(geometry) if geometry else None,
            sub_pipeline=sub,
            sub_pipeline_kwargs=run_kwargs,
        ))
        return self


class Pipeline:
    """Declarative builder for chaining tool steps.

    Steps added via :meth:`add` are executed sequentially; the output
    geometry of each step automatically becomes the input for the next.

    Use :meth:`branch` to fork the pipeline into parallel sub-pipelines
    that share the geometry produced by the trunk.
    """

    def __init__(
        self,
        name: str,
        *,
        base_dir: Optional[Path] = None,
        defaults: Optional[Dict[str, Any]] = None,
    ):
        self.name = name
        self.base_dir = base_dir
        self._defaults: Dict[str, Any] = dict(defaults or {})
        self._trunk: List[_StepSpec] = []
        self._branches: Dict[str, Pipeline] = {}
        self._on_step: Optional[Callable[[StepResult], None]] = None
        self._inherited_artifacts: Dict[str, Path] = {}  # from parent trunk

    # ------------------------------------------------------------------
    # Building
    # ------------------------------------------------------------------

    def add(
        self,
        step_name: str,
        *,
        geometry: Optional[str | Path] = None,
        label: str = "",
        **kwargs: Any,
    ) -> "Pipeline":
        """Append a step to the trunk of this pipeline.

        Parameters
        ----------
        step_name : registered adapter name (e.g. ``"xtb_opt"``)
        geometry : explicit input geometry (overrides auto-propagation)
        label : human-readable label for logging/UI
        **kwargs : adapter-specific parameters (merged with pipeline defaults;
                   explicit kwargs override defaults)
        """
        # Merge: defaults ← step kwargs (step wins)
        merged = {**self._defaults, **kwargs}
        spec = _StepSpec(
            step_name=step_name,
            kwargs=merged,
            geometry_override=Path(geometry) if geometry else None,
            label=label or step_name,
        )
        self._trunk.append(spec)
        return self

    def add_if(
        self,
        condition: Callable,
        step_name: str,
        *,
        label: str = "",
        **kwargs: Any,
    ) -> "Pipeline":
        """Add a step that only runs if ``condition(results, last_result)`` is True.

        Example::

            # Only run IMAG fix if imaginary frequencies were found
            pipe.add("orca_freq", charge=0)
            pipe.add_if(
                lambda results, last: last.data.get("has_imaginary", False),
                "imag_fix", charge=0,
            )

            # Only run expensive DFT if xTB energy is below threshold
            pipe.add_if(
                lambda results, last: last.data.get("energy_Eh", 0) < -100,
                "orca_opt", charge=0,
            )
        """
        merged = {**self._defaults, **kwargs}
        spec = _StepSpec(
            step_name=step_name,
            kwargs=merged,
            label=label or f"{step_name}?",
            condition=condition,
        )
        self._trunk.append(spec)
        return self

    def add_loop(
        self,
        step_name: str,
        *,
        until: Callable,
        max_iter: int = 10,
        label: str = "",
        **kwargs: Any,
    ) -> "Pipeline":
        """Add a step that repeats until ``until(result, iteration)`` returns True.

        The output geometry of each iteration becomes the input for the next.

        Example::

            # IMAG elimination loop — repeat until no imaginary frequencies
            pipe.add_loop(
                "imag_fix", charge=0,
                until=lambda result, i: result.data.get("n_imaginary", 0) == 0,
                max_iter=5,
            )

            # Optimize until energy change < threshold
            pipe.add_loop(
                "orca_opt", charge=0, method="B3LYP",
                until=lambda result, i: abs(result.data.get("energy_change", 1)) < 1e-6,
                max_iter=20,
            )
        """
        merged = {**self._defaults, **kwargs}
        spec = _StepSpec(
            step_name=step_name,
            kwargs=merged,
            label=label or f"{step_name}×{max_iter}",
            loop_until=until,
            loop_max=max_iter,
        )
        self._trunk.append(spec)
        return self

    def add_transform(
        self,
        step_name: str,
        transform: Callable,
        *,
        label: str = "",
        **kwargs: Any,
    ) -> "Pipeline":
        """Add a step whose kwargs are modified at runtime based on previous results.

        ``transform(kwargs, results, last_result)`` should return the modified kwargs dict.

        Example::

            # Set charge from previous step's computed optimal charge
            pipe.add_transform(
                "orca_opt",
                lambda kw, results, last: {**kw, "charge": last.data["optimal_charge"]},
            )

            # Set multiplicity based on electron count
            pipe.add_transform(
                "orca_sp",
                lambda kw, results, last: {
                    **kw,
                    "mult": 1 if last.data["n_electrons"] % 2 == 0 else 2,
                },
                charge=0, method="B3LYP",
            )
        """
        merged = {**self._defaults, **kwargs}
        spec = _StepSpec(
            step_name=step_name,
            kwargs=merged,
            label=label or f"{step_name}~",
            transform=transform,
        )
        self._trunk.append(spec)
        return self

    def add_compute(
        self,
        compute: Callable,
        *,
        label: str = "compute",
    ) -> "Pipeline":
        """Add a pure-Python computation step (no tool, just a function).

        ``compute(results, last_result, work_dir)`` should return a dict
        that is stored in ``result.data``.

        Example::

            # Compute redox potential from branch energies
            def calc_redox(results, last, work_dir):
                E_neutral = results[2].data["energy_Eh"]
                E_ox = results[3].data["energy_Eh"]
                F = 96485.3329
                E_redox = -((E_neutral - E_ox) * 2625.5) / F - 4.28
                return {"E_redox_V": E_redox}

            pipe.add_compute(calc_redox, label="redox_potential")

            # Save CSV of all energies
            def save_csv(results, last, work_dir):
                import csv
                path = work_dir / "energies.csv"
                with open(path, "w", newline="") as f:
                    w = csv.writer(f)
                    w.writerow(["step", "energy_Eh"])
                    for r in results:
                        e = r.data.get("energy_Eh")
                        if e is not None:
                            w.writerow([r.step_name, e])
                return {"csv_path": str(path)}

            pipe.add_compute(save_csv, label="export_csv")
        """
        spec = _StepSpec(
            step_name="_compute",
            kwargs={},
            label=label,
            compute_fn=compute,
        )
        self._trunk.append(spec)
        return self

    def add_fan_out(
        self,
        step_name: str,
        geometries_from: Callable,
        *,
        label: str = "",
        **kwargs: Any,
    ) -> "Pipeline":
        """Run a step in parallel on multiple geometries (one-to-many).

        ``geometries_from(results, last_result)`` should return a list of
        Path objects.  The step runs once per geometry in parallel threads.
        Results are collected as a list in ``result.data["fan_out_results"]``.
        The best (lowest energy) geometry becomes the output geometry.

        Example::

            # CREST produces ensemble → run ORCA SP on each conformer
            pipe.add("crest_conformers", charge=0)
            pipe.add_fan_out(
                "orca_sp",
                geometries_from=lambda results, last: list(
                    last.work_dir.glob("crest_conformers.xyz")  # or parse ensemble
                ),
                charge=0, method="B3LYP", basis="def2-SVP",
                label="screen_conformers",
            )
        """
        merged = {**self._defaults, **kwargs}
        spec = _StepSpec(
            step_name=step_name,
            kwargs=merged,
            label=label or f"{step_name}[*]",
            fan_out_fn=geometries_from,
        )
        self._trunk.append(spec)
        return self

    def add_retry(
        self,
        step_name: str,
        *,
        max_attempts: int = 3,
        delay: float = 0.0,
        label: str = "",
        **kwargs: Any,
    ) -> "Pipeline":
        """Add a step that retries on failure up to *max_attempts* times.

        Example::

            # Retry ORCA SP up to 3 times (transient SCF failures)
            pipe.add_retry("orca_sp", max_attempts=3, charge=0, method="B3LYP")

            # Retry with delay between attempts
            pipe.add_retry("xtb_opt", max_attempts=5, delay=2.0, charge=0)
        """
        merged = {**self._defaults, **kwargs}
        spec = _StepSpec(
            step_name=step_name,
            kwargs=merged,
            label=label or f"{step_name}↻{max_attempts}",
            retry_max=max_attempts,
            retry_delay=delay,
        )
        self._trunk.append(spec)
        return self

    def add_map(
        self,
        step_name: str,
        items_from: Callable,
        *,
        label: str = "",
        **kwargs: Any,
    ) -> "Pipeline":
        """Run a step once per item from a list (high-throughput screening).

        ``items_from(results, last_result)`` returns a list of dicts.
        Each dict is merged into the step kwargs for that run.
        If the dict contains ``"geometry"``, it overrides the input geometry.

        Example::

            # Screen 1000 SMILES through the same workflow
            pipe.add_map(
                "smiles_to_xyz",
                items_from=lambda r, l: [{"smiles": s} for s in smiles_list],
                label="screen_smiles",
            )

            # Run ORCA on multiple geometries with different charges
            pipe.add_map(
                "orca_sp",
                items_from=lambda r, l: [
                    {"geometry": p, "charge": 0} for p in Path("inputs").glob("*.xyz")
                ],
                method="B3LYP", basis="def2-SVP",
                label="batch_sp",
            )
        """
        merged = {**self._defaults, **kwargs}
        spec = _StepSpec(
            step_name=step_name,
            kwargs=merged,
            label=label or f"{step_name}[map]",
            map_items=items_from,
        )
        self._trunk.append(spec)
        return self

    def add_checkpoint(self, label: str = "checkpoint") -> "Pipeline":
        """Save pipeline state after this point for crash recovery.

        Creates a ``checkpoint.json`` file in the working directory.
        Use ``Pipeline.resume(checkpoint_path)`` to restart from last
        successful checkpoint.

        Example::

            pipe.add("smiles_to_xyz", smiles="CCO")
            pipe.add("xtb_opt", charge=0)
            pipe.add_checkpoint()  # save state here
            pipe.add("orca_opt", charge=0, method="B3LYP")  # expensive
            pipe.add_checkpoint()  # save again after ORCA
        """
        spec = _StepSpec(
            step_name="_checkpoint",
            kwargs={},
            label=label,
            checkpoint=True,
        )
        self._trunk.append(spec)
        return self

    def add_sub_pipeline(
        self,
        sub: "Pipeline",
        *,
        label: str = "",
        pass_geometry: bool = True,
        **run_kwargs: Any,
    ) -> "Pipeline":
        """Run an entire Pipeline as a single step (composition).

        The sub-pipeline receives the current geometry and runs to
        completion.  Its final geometry becomes this pipeline's current
        geometry, and its full results are available in
        ``result.data["sub_results"]``.

        Example::

            # Build a reusable IMAG-fix sub-pipeline
            imag_fix = Pipeline("imag_fix")
            imag_fix.add("orca_freq", charge=0)
            imag_fix.add_loop("imag_fix", charge=0,
                              until=lambda r, i: r.data.get("n_imaginary", 0) == 0,
                              max_iter=5)

            # Use it inside a larger workflow
            main = Pipeline("main")
            main.add("smiles_to_xyz", smiles="CCO")
            main.add("xtb_opt", charge=0)
            main.add("orca_opt", charge=0)
            main.add_sub_pipeline(imag_fix, label="fix_imag")
            main.add("orca_sp", charge=0, method="B3LYP", basis="def2-TZVP")

            # OCCUPIER pattern: sub-pipeline per FoB stage
            for fob in [0.0, 0.25, 0.5, 0.75, 1.0]:
                stage = Pipeline(f"fob_{fob}")
                stage.add("orca_opt", charge=0, fob=fob)
                stage.add("orca_freq", charge=0, fob=fob)
                main.add_sub_pipeline(stage, label=f"fob_{fob}")
        """
        spec = _StepSpec(
            step_name="_sub_pipeline",
            kwargs={},
            label=label or f"sub:{sub.name}",
            sub_pipeline=sub,
            sub_pipeline_kwargs=run_kwargs,
        )
        self._trunk.append(spec)
        return self

    def add_reactive(
        self,
        callback: Callable,
        *,
        label: str = "reactive",
    ) -> "Pipeline":
        """Dynamically insert steps based on previous results.

        ``callback(results, last_result, pipeline_inserter)`` receives a
        :class:`PipelineInserter` that can add steps to be executed
        immediately.  This enables patterns like OCCUPIER's
        contamination-triggered branching or ESD's conditional state
        population.

        Example::

            # OCCUPIER: add broken-symmetry job if contamination detected
            def check_contamination(results, last, inserter):
                s2 = last.data.get("s2_deviation", 0)
                if s2 > 0.1:
                    inserter.add("orca_sp", charge=0, broken_symmetry=True,
                                 label="bs_correction")

            pipe.add("orca_sp", charge=0)
            pipe.add_reactive(check_contamination, label="contamination_check")

            # ESD: populate state jobs based on config
            def populate_states(results, last, inserter):
                for state in ["S1", "S2", "T1"]:
                    inserter.add("orca_opt", charge=0, root=state,
                                 label=f"opt_{state}")

            pipe.add_reactive(populate_states, label="state_population")

            # Dynamic fan-out: create sub-pipelines on the fly
            def dynamic_stages(results, last, inserter):
                conformers = list(last.work_dir.glob("*.xyz"))
                for conf in conformers:
                    stage = Pipeline(f"refine_{conf.stem}")
                    stage.add("orca_opt", charge=0)
                    stage.add("orca_freq", charge=0)
                    inserter.add_sub_pipeline(stage, geometry=conf)

            pipe.add_reactive(dynamic_stages, label="dynamic_refine")
        """
        spec = _StepSpec(
            step_name="_reactive",
            kwargs={},
            label=label,
            reactive_fn=callback,
        )
        self._trunk.append(spec)
        return self

    # ------------------------------------------------------------------
    # Internal: fan-out execution
    # ------------------------------------------------------------------

    @staticmethod
    def _execute_fan_out(
        spec: _StepSpec,
        step_idx: int,
        base: Path,
        cores: int,
        results: list,
        last_result: Optional[StepResult],
        current_artifacts: dict,
        run_step_fn: Callable,
    ) -> StepResult:
        """Run a single step on multiple geometries in parallel."""
        import time as _time
        start = _time.monotonic()

        try:
            geom_list = spec.fan_out_fn(results, last_result)
        except Exception as exc:
            return StepResult(
                step_name=spec.step_name, status=StepStatus.FAILED,
                error=f"fan_out geometries_from raised: {exc}",
                elapsed_seconds=_time.monotonic() - start,
            )

        if not geom_list:
            return StepResult(
                step_name=spec.step_name, status=StepStatus.SKIPPED,
                error="fan_out returned empty geometry list",
                elapsed_seconds=_time.monotonic() - start,
            )

        fan_results: List[StepResult] = []
        cores_per = max(1, cores // len(geom_list))

        def _run_one(idx: int, geom: Path) -> None:
            r = run_step_fn(
                spec.step_name,
                geometry=geom,
                cores=cores_per,
                work_dir=base / f"{step_idx:02d}_{spec.label}" / f"geom_{idx:03d}",
                **spec.kwargs,
            )
            fan_results.append(r)

        threads = []
        for idx, geom in enumerate(geom_list):
            t = threading.Thread(target=_run_one, args=(idx, Path(geom)))
            threads.append(t)
            t.start()
        for t in threads:
            t.join()

        # Find best result (lowest energy)
        best = None
        best_energy = float("inf")
        for r in fan_results:
            if r.ok:
                e = r.data.get("energy_Eh") or r.data.get("energy_eV") or float("inf")
                if e < best_energy:
                    best_energy = e
                    best = r

        n_ok = sum(1 for r in fan_results if r.ok)
        data = {
            "fan_out_count": len(geom_list),
            "fan_out_ok": n_ok,
            "fan_out_failed": len(geom_list) - n_ok,
        }
        if best and best_energy < float("inf"):
            data["best_energy"] = best_energy

        return StepResult(
            step_name=spec.step_name,
            status=StepStatus.SUCCESS if n_ok > 0 else StepStatus.FAILED,
            geometry=best.geometry if best else None,
            work_dir=base / f"{step_idx:02d}_{spec.label}",
            data=data,
            artifacts=best.artifacts if best else {},
            error=None if n_ok > 0 else "all fan-out jobs failed",
            elapsed_seconds=_time.monotonic() - start,
        )

    # ------------------------------------------------------------------
    # Internal: map execution
    # ------------------------------------------------------------------

    @staticmethod
    def _execute_map(
        spec: _StepSpec,
        step_idx: int,
        base: Path,
        cores: int,
        results: list,
        last_result: Optional[StepResult],
        current_geom: Optional[Path],
        current_artifacts: dict,
        run_step_fn: Callable,
    ) -> StepResult:
        """Run a single step once per item in parallel."""
        import time as _time
        start = _time.monotonic()

        try:
            items = spec.map_items(results, last_result)
        except Exception as exc:
            return StepResult(
                step_name=spec.step_name, status=StepStatus.FAILED,
                error=f"map items_from raised: {exc}",
                elapsed_seconds=_time.monotonic() - start,
            )

        if not items:
            return StepResult(
                step_name=spec.step_name, status=StepStatus.SKIPPED,
                error="map returned empty item list",
                elapsed_seconds=_time.monotonic() - start,
            )

        map_results: List[StepResult] = []
        cores_per = max(1, cores // len(items))

        def _run_one(idx: int, item: dict) -> None:
            merged = {**spec.kwargs, **item}
            geom = merged.pop("geometry", current_geom)
            if geom is not None:
                geom = Path(geom)
            r = run_step_fn(
                spec.step_name,
                geometry=geom,
                cores=cores_per,
                work_dir=base / f"{step_idx:02d}_{spec.label}" / f"item_{idx:04d}",
                **merged,
            )
            map_results.append(r)

        threads = []
        for idx, item in enumerate(items):
            t = threading.Thread(target=_run_one, args=(idx, dict(item)))
            threads.append(t)
            t.start()
        for t in threads:
            t.join()

        # Find best result (lowest energy) for geometry propagation
        best = None
        best_energy = float("inf")
        for r in map_results:
            if r.ok:
                e = r.data.get("energy_Eh") or r.data.get("energy_eV") or float("inf")
                if e < best_energy:
                    best_energy = e
                    best = r

        n_ok = sum(1 for r in map_results if r.ok)
        data = {
            "map_count": len(items),
            "map_ok": n_ok,
            "map_failed": len(items) - n_ok,
            "map_results": [
                {"ok": r.ok, "data": r.data, "geometry": str(r.geometry) if r.geometry else None}
                for r in map_results
            ],
        }
        if best and best_energy < float("inf"):
            data["best_energy"] = best_energy

        return StepResult(
            step_name=spec.step_name,
            status=StepStatus.SUCCESS if n_ok > 0 else StepStatus.FAILED,
            geometry=best.geometry if best else None,
            work_dir=base / f"{step_idx:02d}_{spec.label}",
            data=data,
            artifacts=best.artifacts if best else {},
            error=None if n_ok > 0 else "all map jobs failed",
            elapsed_seconds=_time.monotonic() - start,
        )

    # ------------------------------------------------------------------
    # Internal: checkpoint save/restore
    # ------------------------------------------------------------------

    @staticmethod
    def _save_checkpoint(
        checkpoint_path: Path,
        pipeline_name: str,
        step_idx: int,
        results: List[StepResult],
        current_geom: Optional[Path],
        current_artifacts: Dict[str, Path],
    ) -> None:
        """Save pipeline state to JSON for crash recovery."""
        import json
        state = {
            "pipeline_name": pipeline_name,
            "step_idx": step_idx,
            "current_geometry": str(current_geom) if current_geom else None,
            "current_artifacts": {k: str(v) for k, v in current_artifacts.items()},
            "results": [
                {
                    "step_name": r.step_name,
                    "status": r.status.value,
                    "geometry": str(r.geometry) if r.geometry else None,
                    "work_dir": str(r.work_dir) if r.work_dir else None,
                    "data": r.data,
                    "artifacts": {k: str(v) for k, v in r.artifacts.items()},
                    "error": r.error,
                    "elapsed_seconds": r.elapsed_seconds,
                }
                for r in results
            ],
        }
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        checkpoint_path.write_text(json.dumps(state, indent=2, default=str))
        logger.info("Checkpoint saved: %s (after step %d)", checkpoint_path, step_idx)

    @staticmethod
    def _load_checkpoint(checkpoint_path: Path) -> dict:
        """Load checkpoint state from JSON."""
        import json
        state = json.loads(checkpoint_path.read_text())
        # Reconstruct StepResult objects
        results = []
        for rd in state["results"]:
            results.append(StepResult(
                step_name=rd["step_name"],
                status=StepStatus(rd["status"]),
                geometry=Path(rd["geometry"]) if rd["geometry"] else None,
                work_dir=Path(rd["work_dir"]) if rd["work_dir"] else None,
                data=rd.get("data", {}),
                artifacts={k: Path(v) for k, v in rd.get("artifacts", {}).items()},
                error=rd.get("error"),
                elapsed_seconds=rd.get("elapsed_seconds", 0.0),
            ))
        state["results"] = results
        if state["current_geometry"]:
            state["current_geometry"] = Path(state["current_geometry"])
        state["current_artifacts"] = {
            k: Path(v) for k, v in state.get("current_artifacts", {}).items()
        }
        return state

    @classmethod
    def resume(cls, checkpoint_path: Path | str) -> dict:
        """Load checkpoint data for resuming a pipeline.

        Returns a dict with keys: ``pipeline_name``, ``step_idx``,
        ``current_geometry``, ``current_artifacts``, ``results``.
        Pass ``resume_from`` to :meth:`run` to continue from the checkpoint.

        Example::

            state = Pipeline.resume("work_dir/checkpoint.json")
            pipe.run(
                resume_from=state,
                cores=8,
            )
        """
        return cls._load_checkpoint(Path(checkpoint_path))

    def branch(self, name: str) -> "Pipeline":
        """Create a named branch that starts from the trunk's final geometry.

        Branches run in parallel when :meth:`run` is called.
        Returns the new branch Pipeline for chaining.
        """
        if name in self._branches:
            return self._branches[name]
        child = Pipeline(
            f"{self.name}/{name}",
            base_dir=self.base_dir,
            defaults=self._defaults,  # inherit parent defaults
        )
        self._branches[name] = child
        return child

    def on_step(self, callback: Callable[[StepResult], None]) -> "Pipeline":
        """Set a callback invoked after each step completes."""
        self._on_step = callback
        return self

    # ------------------------------------------------------------------
    # Execution — simple sequential
    # ------------------------------------------------------------------

    def run(
        self,
        *,
        cores: int | str = 1,
        geometry: Optional[str | Path] = None,
        stop_on_failure: bool = True,
        work_dir: Optional[Path] = None,
        resume_from: Optional[dict] = None,
        provenance: bool = False,
    ) -> PipelineResult:
        """Execute the pipeline sequentially, then branches in parallel.

        Parameters
        ----------
        cores : total CPU cores available.  Pass ``"auto"`` to auto-detect
                from the cluster environment (SLURM, PBS, or local system).
        geometry : initial geometry for the first step
        stop_on_failure : abort remaining steps on first failure
        work_dir : root directory for all step working directories
        resume_from : checkpoint state dict from ``Pipeline.resume()``
                      — skips steps already completed in the checkpoint
        provenance : if True, write a ``provenance.json`` log to *work_dir*
                     after execution (all steps, timing, parameters)

        Returns
        -------
        PipelineResult with all step results
        """
        from delfin.tools._runner import run_step

        cores = _resolve_cores(cores)

        import time as _pipeline_time

        pipeline_start = _pipeline_time.monotonic()

        base = Path(work_dir) if work_dir else (self.base_dir or Path.cwd())
        base.mkdir(parents=True, exist_ok=True)

        results: List[StepResult] = []
        current_geom: Optional[Path] = Path(geometry) if geometry else None
        current_artifacts: Dict[str, Path] = dict(self._inherited_artifacts)  # start with inherited
        skip_until = -1  # for resume: skip steps already done

        # --- Resume from checkpoint ---
        if resume_from is not None:
            skip_until = resume_from.get("step_idx", -1)
            results = list(resume_from.get("results", []))
            if resume_from.get("current_geometry"):
                current_geom = resume_from["current_geometry"]
            if resume_from.get("current_artifacts"):
                current_artifacts.update(resume_from["current_artifacts"])
            logger.info("[%s] resuming from checkpoint (step %d)", self.name, skip_until)

        # --- Execute trunk sequentially ---
        for i, spec in enumerate(self._trunk):
            # Skip already-completed steps when resuming
            if i <= skip_until:
                continue

            last_result = results[-1] if results else None

            # --- Checkpoint: save state ---
            if spec.checkpoint:
                import time as _time
                _start = _time.monotonic()
                cp_path = base / f"{i:02d}_{spec.label}" / "checkpoint.json"
                try:
                    self._save_checkpoint(cp_path, self.name, i, results,
                                          current_geom, current_artifacts)
                    result = StepResult(
                        step_name="_checkpoint", status=StepStatus.SUCCESS,
                        work_dir=cp_path.parent, data={"checkpoint_path": str(cp_path)},
                        elapsed_seconds=_time.monotonic() - _start,
                    )
                except Exception as exc:
                    result = StepResult(
                        step_name="_checkpoint", status=StepStatus.FAILED,
                        error=str(exc), elapsed_seconds=_time.monotonic() - _start,
                    )
                results.append(result)
                if self._on_step:
                    self._on_step(result)
                continue

            # --- Sub-pipeline: run an entire pipeline as one step ---
            if spec.sub_pipeline is not None:
                import time as _time
                _start = _time.monotonic()
                sub_dir = base / f"{i:02d}_{spec.label}"
                sub_geom = spec.geometry_override or current_geom
                sub_kwargs = dict(spec.sub_pipeline_kwargs or {})
                sub_kwargs.setdefault("stop_on_failure", stop_on_failure)

                try:
                    sub_result = spec.sub_pipeline.run(
                        cores=cores,
                        geometry=sub_geom,
                        work_dir=sub_dir,
                        **sub_kwargs,
                    )
                    # Extract final geometry and data from sub-pipeline
                    sub_last = sub_result.last
                    sub_data = {
                        "sub_pipeline": spec.sub_pipeline.name,
                        "sub_ok": sub_result.ok,
                        "sub_steps": len(sub_result.results),
                        "sub_results": [
                            {"step": r.step_name, "ok": r.ok,
                             "data": r.data, "elapsed": r.elapsed_seconds}
                            for r in sub_result.results
                        ],
                    }
                    if sub_last and sub_last.data:
                        sub_data.update(sub_last.data)

                    result = StepResult(
                        step_name=f"_sub:{spec.sub_pipeline.name}",
                        status=StepStatus.SUCCESS if sub_result.ok else StepStatus.FAILED,
                        geometry=sub_last.geometry if sub_last else None,
                        work_dir=sub_dir,
                        data=sub_data,
                        artifacts=sub_last.artifacts if sub_last else {},
                        error=None if sub_result.ok else f"sub-pipeline '{spec.sub_pipeline.name}' failed",
                        elapsed_seconds=_time.monotonic() - _start,
                    )
                except Exception as exc:
                    result = StepResult(
                        step_name=f"_sub:{spec.sub_pipeline.name}",
                        status=StepStatus.FAILED,
                        work_dir=sub_dir,
                        error=f"sub-pipeline raised: {exc}",
                        elapsed_seconds=_time.monotonic() - _start,
                    )

                results.append(result)
                if self._on_step:
                    self._on_step(result)
                if result.ok:
                    if result.geometry:
                        current_geom = result.geometry
                    if result.artifacts:
                        current_artifacts.update(result.artifacts)
                elif stop_on_failure:
                    return PipelineResult(name=self.name, results=results, branch_results={})
                continue

            # --- Reactive: dynamically insert steps ---
            if spec.reactive_fn is not None:
                import time as _time
                _start = _time.monotonic()
                inserter = PipelineInserter(self._defaults)
                try:
                    spec.reactive_fn(results, last_result, inserter)
                except Exception as exc:
                    logger.error("[%s] reactive '%s' raised: %s", self.name, spec.label, exc)
                    results.append(StepResult(
                        step_name="_reactive", status=StepStatus.FAILED,
                        error=f"reactive callback raised: {exc}",
                        elapsed_seconds=_time.monotonic() - _start,
                    ))
                    if stop_on_failure:
                        return PipelineResult(name=self.name, results=results, branch_results={})
                    continue

                if not inserter._steps:
                    # Reactive produced no steps — record as success with no-op
                    results.append(StepResult(
                        step_name="_reactive", status=StepStatus.SUCCESS,
                        data={"reactive_label": spec.label, "steps_inserted": 0},
                        elapsed_seconds=_time.monotonic() - _start,
                    ))
                    if self._on_step:
                        self._on_step(results[-1])
                    continue

                # Execute inserted steps immediately
                reactive_results = []
                for j, rspec in enumerate(inserter._steps):
                    r_geom = rspec.geometry_override or current_geom

                    if rspec.sub_pipeline is not None:
                        # Inserted sub-pipeline
                        sub_dir = base / f"{i:02d}_{spec.label}" / f"sub_{j:02d}_{rspec.label}"
                        sub_kwargs = dict(rspec.sub_pipeline_kwargs or {})
                        sub_kwargs.setdefault("stop_on_failure", stop_on_failure)
                        sub_result = rspec.sub_pipeline.run(
                            cores=cores, geometry=r_geom, work_dir=sub_dir, **sub_kwargs,
                        )
                        sub_last = sub_result.last
                        r = StepResult(
                            step_name=f"_sub:{rspec.sub_pipeline.name}",
                            status=StepStatus.SUCCESS if sub_result.ok else StepStatus.FAILED,
                            geometry=sub_last.geometry if sub_last else None,
                            work_dir=sub_dir,
                            data={"sub_pipeline": rspec.sub_pipeline.name, "sub_ok": sub_result.ok},
                            artifacts=sub_last.artifacts if sub_last else {},
                            elapsed_seconds=sum(sr.elapsed_seconds for sr in sub_result.results),
                        )
                    else:
                        # Inserted normal step
                        r_dir = base / f"{i:02d}_{spec.label}" / f"step_{j:02d}_{rspec.label}"
                        r = run_step(
                            rspec.step_name, geometry=r_geom, cores=cores,
                            work_dir=r_dir, **rspec.kwargs,
                        )

                    reactive_results.append(r)
                    if self._on_step:
                        self._on_step(r)
                    if r.ok:
                        if r.geometry:
                            current_geom = r.geometry
                        if r.artifacts:
                            current_artifacts.update(r.artifacts)
                    elif stop_on_failure:
                        results.extend(reactive_results)
                        return PipelineResult(name=self.name, results=results, branch_results={})

                # Record all reactive results
                results.extend(reactive_results)
                continue

            # --- Conditional: skip if condition returns False ---
            if spec.condition is not None:
                try:
                    should_run = spec.condition(results, last_result)
                except Exception as exc:
                    logger.warning("[%s] condition for '%s' raised: %s", self.name, spec.label, exc)
                    should_run = False
                if not should_run:
                    logger.info("[%s] skipping '%s' (condition not met)", self.name, spec.label)
                    results.append(StepResult(
                        step_name=spec.step_name, status=StepStatus.SKIPPED,
                        error="condition not met",
                    ))
                    continue

            # --- Compute: pure Python function, no tool ---
            if spec.compute_fn is not None:
                import time as _time
                _start = _time.monotonic()
                compute_dir = base / f"{i:02d}_{spec.label}"
                compute_dir.mkdir(parents=True, exist_ok=True)
                try:
                    data = spec.compute_fn(results, last_result, compute_dir)
                    if not isinstance(data, dict):
                        data = {"return_value": data}
                    result = StepResult(
                        step_name="_compute", status=StepStatus.SUCCESS,
                        work_dir=compute_dir, data=data,
                        elapsed_seconds=_time.monotonic() - _start,
                    )
                except Exception as exc:
                    result = StepResult(
                        step_name="_compute", status=StepStatus.FAILED,
                        work_dir=compute_dir, error=str(exc),
                        elapsed_seconds=_time.monotonic() - _start,
                    )
                results.append(result)
                if self._on_step:
                    self._on_step(result)
                if not result.ok and stop_on_failure:
                    return PipelineResult(name=self.name, results=results, branch_results={})
                continue

            # --- Fan-out: run step on multiple geometries in parallel ---
            if spec.fan_out_fn is not None:
                result = self._execute_fan_out(
                    spec, i, base, cores, results, last_result, current_artifacts, run_step,
                )
                results.append(result)
                if self._on_step:
                    self._on_step(result)
                if result.ok:
                    if result.geometry:
                        current_geom = result.geometry
                    if result.artifacts:
                        current_artifacts.update(result.artifacts)
                elif stop_on_failure:
                    return PipelineResult(name=self.name, results=results, branch_results={})
                continue

            # --- Map: run step once per item in parallel ---
            if spec.map_items is not None:
                result = self._execute_map(
                    spec, i, base, cores, results, last_result,
                    current_geom, current_artifacts, run_step,
                )
                results.append(result)
                if self._on_step:
                    self._on_step(result)
                if result.ok:
                    if result.geometry:
                        current_geom = result.geometry
                    if result.artifacts:
                        current_artifacts.update(result.artifacts)
                elif stop_on_failure:
                    return PipelineResult(name=self.name, results=results, branch_results={})
                continue

            geom = spec.geometry_override or current_geom

            # --- Transform: modify kwargs based on previous results ---
            merged_kwargs = dict(spec.kwargs)
            if spec.transform is not None:
                try:
                    merged_kwargs = spec.transform(merged_kwargs, results, last_result)
                except Exception as exc:
                    logger.error("[%s] transform for '%s' raised: %s", self.name, spec.label, exc)
                    results.append(StepResult(
                        step_name=spec.step_name, status=StepStatus.FAILED,
                        error=f"transform raised: {exc}",
                    ))
                    if stop_on_failure:
                        return PipelineResult(name=self.name, results=results, branch_results={})
                    continue

            # Auto-inject artifacts from previous steps (MOREAD chaining)
            if current_artifacts:
                if "moread" not in merged_kwargs and "gbw" in current_artifacts:
                    if spec.step_name.startswith("orca_"):
                        merged_kwargs["moread"] = str(current_artifacts["gbw"])
                merged_kwargs.setdefault("_prev_artifacts", dict(current_artifacts))

            # --- Loop: repeat step until condition met ---
            if spec.loop_until is not None and spec.loop_max > 1:
                loop_geom = geom
                loop_result = None
                for iteration in range(spec.loop_max):
                    step_dir = base / f"{i:02d}_{spec.label}_iter{iteration}"
                    logger.info("[%s] step %d/%d: %s (iter %d/%d)",
                                self.name, i + 1, len(self._trunk), spec.label,
                                iteration + 1, spec.loop_max)
                    loop_result = run_step(
                        spec.step_name, geometry=loop_geom, cores=cores,
                        work_dir=step_dir, **merged_kwargs,
                    )
                    if self._on_step:
                        self._on_step(loop_result)
                    if not loop_result.ok:
                        break
                    if loop_result.geometry:
                        loop_geom = loop_result.geometry
                    if loop_result.artifacts:
                        current_artifacts.update(loop_result.artifacts)
                    try:
                        if spec.loop_until(loop_result, iteration):
                            break
                    except Exception as exc:
                        logger.warning("[%s] loop_until raised: %s", self.name, exc)
                        break

                result = loop_result
                if result is not None:
                    results.append(result)
                    if result.ok and result.geometry:
                        current_geom = result.geometry
                    elif not result.ok:
                        logger.warning("[%s] loop '%s' failed at iter %d", self.name, spec.label, iteration)
                        if stop_on_failure:
                            return PipelineResult(name=self.name, results=results, branch_results={})
                continue

            # --- Normal step execution (with retry support) ---
            step_dir = base / f"{i:02d}_{spec.label}"
            logger.info("[%s] step %d/%d: %s", self.name, i + 1, len(self._trunk), spec.label)

            result = None
            for attempt in range(spec.retry_max):
                attempt_dir = step_dir if spec.retry_max == 1 else step_dir / f"attempt_{attempt}"
                if attempt > 0:
                    logger.info("[%s] retrying '%s' (attempt %d/%d)",
                                self.name, spec.label, attempt + 1, spec.retry_max)
                    if spec.retry_delay > 0:
                        import time as _time
                        _time.sleep(spec.retry_delay)

                result = run_step(
                    spec.step_name,
                    geometry=geom,
                    cores=cores,
                    work_dir=attempt_dir,
                    **merged_kwargs,
                )
                if result.ok:
                    break

            results.append(result)

            if self._on_step:
                self._on_step(result)

            if result.ok:
                if result.geometry:
                    current_geom = result.geometry
                if result.artifacts:
                    current_artifacts.update(result.artifacts)
            elif not result.ok:
                logger.warning("[%s] step '%s' failed: %s", self.name, spec.label, result.error)
                if stop_on_failure:
                    return PipelineResult(
                        name=self.name,
                        results=results,
                        branch_results={},
                    )

        # --- Execute branches in parallel ---
        branch_results: Dict[str, PipelineResult] = {}
        if self._branches:
            cores_per_branch = max(1, cores // len(self._branches))
            threads: Dict[str, threading.Thread] = {}
            branch_out: Dict[str, PipelineResult] = {}

            def _run_branch(bname: str, bpipe: Pipeline) -> None:
                # Branches inherit trunk artifacts
                bpipe._inherited_artifacts = dict(current_artifacts)
                branch_out[bname] = bpipe.run(
                    cores=cores_per_branch,
                    geometry=current_geom,
                    stop_on_failure=stop_on_failure,
                    work_dir=base / bname,
                )

            for bname, bpipe in self._branches.items():
                # Propagate callback
                if self._on_step and bpipe._on_step is None:
                    bpipe._on_step = self._on_step
                t = threading.Thread(target=_run_branch, args=(bname, bpipe), daemon=True)
                threads[bname] = t
                t.start()

            for bname, t in threads.items():
                t.join()

            branch_results = branch_out

        pipeline_result = PipelineResult(
            name=self.name,
            results=results,
            branch_results=branch_results,
        )

        # --- Provenance logging ---
        if provenance:
            _write_provenance(base, self.name, pipeline_result,
                              cores, _pipeline_time.monotonic() - pipeline_start)

        return pipeline_result

    # ------------------------------------------------------------------
    # Execution — scheduler-integrated
    # ------------------------------------------------------------------

    def run_scheduled(
        self,
        *,
        cores: int | str = 1,
        geometry: Optional[str | Path] = None,
        work_dir: Optional[Path] = None,
        config: Optional[Dict[str, Any]] = None,
    ) -> PipelineResult:
        """Execute using the WorkflowManager/DynamicCorePool scheduler.

        This integrates with the existing DELFIN scheduling infrastructure
        for proper resource management alongside other running jobs.

        Parameters
        ----------
        cores : CPU cores (used if config not provided)
        geometry : initial geometry for the first step
        work_dir : root directory for step working directories
        config : DELFIN config dict (from CONTROL.txt); if provided,
                 the pipeline shares the global DynamicCorePool with
                 all other running DELFIN jobs. If ``None``, a minimal
                 config with ``PAL=cores`` is created.
        """
        from delfin.tools._runner import step_as_workflow_job
        from delfin.workflows.engine.classic import _WorkflowManager, WorkflowRunResult
        from delfin.workflows.scheduling.manager import get_global_manager

        cores = _resolve_cores(cores)

        base = Path(work_dir) if work_dir else (self.base_dir or Path.cwd())
        base.mkdir(parents=True, exist_ok=True)

        # Use provided config or build a minimal one.
        # When config comes from CONTROL.txt the GlobalJobManager is
        # already initialised and the _WorkflowManager will share its
        # DynamicCorePool — no resource conflicts.
        if config is None:
            config = {"PAL": cores}
        # Ensure the global manager is initialised for this config
        global_mgr = get_global_manager()
        if not global_mgr.is_initialized():
            global_mgr.ensure_initialized(config)

        manager = _WorkflowManager(config, label=self.name)

        results_map: Dict[str, StepResult] = {}
        # NOTE: results are funneled through _make_callback into results_map;
        # the previous `trunk_geom` dict (step index → output geometry) was
        # declared here but never populated, so it has been removed.

        def _make_callback(step_idx: str):
            def cb(result: StepResult):
                results_map[step_idx] = result
                if self._on_step:
                    self._on_step(result)
            return cb

        # --- Trunk: chain dependencies ---
        prev_job_id = None
        for i, spec in enumerate(self._trunk):
            job_id = f"trunk_{i:02d}_{spec.label}"
            deps = {prev_job_id} if prev_job_id else set()

            # For scheduled execution, geometry propagation is trickier
            # because steps run async. We use a closure to capture the
            # geometry from the previous step's result.
            def _make_work(idx, sp, prev_idx):
                def work(allocated_cores: int) -> None:
                    from delfin.tools._runner import run_step as _run
                    geom = sp.geometry_override
                    if geom is None and prev_idx is not None:
                        prev_result = results_map.get(f"trunk_{prev_idx:02d}_{self._trunk[prev_idx].label}")
                        if prev_result and prev_result.geometry:
                            geom = prev_result.geometry
                    result = _run(
                        sp.step_name,
                        geometry=geom,
                        cores=allocated_cores,
                        work_dir=base / f"{idx:02d}_{sp.label}",
                        **sp.kwargs,
                    )
                    results_map[f"trunk_{idx:02d}_{sp.label}"] = result
                    if self._on_step:
                        self._on_step(result)
                    if not result.ok:
                        raise StepError(sp.step_name, result.error or "failed")
                return work

            from delfin.workflows.engine.classic import WorkflowJob
            job = WorkflowJob(
                job_id=job_id,
                work=_make_work(i, spec, i - 1 if i > 0 else None),
                description=f"step:{spec.step_name}",
                dependencies=deps,
                cores_min=1,
                cores_optimal=min(cores, 2),
                cores_max=cores,
            )
            manager.add_job(job)
            prev_job_id = job_id

        last_trunk_id = prev_job_id

        # --- Branches: depend on last trunk step ---
        for bname, bpipe in self._branches.items():
            prev_branch_id = last_trunk_id
            for j, spec in enumerate(bpipe._trunk):
                job_id = f"{bname}_{j:02d}_{spec.label}"
                deps = {prev_branch_id} if prev_branch_id else set()

                def _make_branch_work(bn, idx, sp, prev_id):
                    def work(allocated_cores: int) -> None:
                        from delfin.tools._runner import run_step as _run
                        geom = sp.geometry_override
                        if geom is None and prev_id:
                            prev_result = results_map.get(prev_id)
                            if prev_result and prev_result.geometry:
                                geom = prev_result.geometry
                        result = _run(
                            sp.step_name,
                            geometry=geom,
                            cores=allocated_cores,
                            work_dir=base / bn / f"{idx:02d}_{sp.label}",
                            **sp.kwargs,
                        )
                        results_map[job_id] = result
                        if self._on_step:
                            self._on_step(result)
                        if not result.ok:
                            raise StepError(sp.step_name, result.error or "failed")
                    return work

                job = WorkflowJob(
                    job_id=job_id,
                    work=_make_branch_work(bname, j, spec, prev_branch_id),
                    description=f"step:{spec.step_name}",
                    dependencies=deps,
                    cores_min=1,
                    cores_optimal=min(cores // max(1, len(self._branches)), 2),
                    cores_max=cores // max(1, len(self._branches)),
                )
                manager.add_job(job)
                prev_branch_id = job_id

        # Execute all jobs — results are gathered via _make_callback into
        # results_map, so the manager.execute() return value is intentionally
        # discarded here.
        manager.execute()

        # Collect results in order
        trunk_results = [
            results_map.get(f"trunk_{i:02d}_{spec.label}",
                           StepResult(step_name=spec.step_name, status=StepStatus.SKIPPED, error="not executed"))
            for i, spec in enumerate(self._trunk)
        ]
        branch_res = {}
        for bname, bpipe in self._branches.items():
            branch_res[bname] = PipelineResult(
                name=bname,
                results=[
                    results_map.get(f"{bname}_{j:02d}_{spec.label}",
                                   StepResult(step_name=spec.step_name, status=StepStatus.SKIPPED, error="not executed"))
                    for j, spec in enumerate(bpipe._trunk)
                ],
                branch_results={},
            )

        return PipelineResult(
            name=self.name,
            results=trunk_results,
            branch_results=branch_res,
        )

    # ------------------------------------------------------------------
    # Conversion to registered workflow
    # ------------------------------------------------------------------

    def as_workflow(self, description: str = "") -> object:
        """Return a Workflow-protocol-compatible object for the registry.

        Usage::

            from delfin.workflows import register_workflow
            pipe = Pipeline("my_workflow")
            pipe.add("smiles_to_xyz", smiles="C")
            pipe.add("xtb_opt", charge=0)
            register_workflow(pipe.as_workflow("My custom workflow"))
        """
        return _PipelineWorkflow(self, description or f"Pipeline: {self.name}")

    def __repr__(self) -> str:
        n_steps = len(self._trunk)
        n_branches = len(self._branches)
        parts = [f"Pipeline({self.name!r}, {n_steps} steps"]
        if n_branches:
            parts.append(f", {n_branches} branches")
        parts.append(")")
        return "".join(parts)


@dataclass
class PipelineResult:
    """Result of executing a Pipeline."""

    name: str
    results: List[StepResult]
    branch_results: Dict[str, "PipelineResult"]

    @property
    def ok(self) -> bool:
        """True if all trunk steps and all branches succeeded."""
        if any(not r.ok for r in self.results):
            return False
        return all(br.ok for br in self.branch_results.values())

    @property
    def last(self) -> Optional[StepResult]:
        """Last trunk result (convenience for sequential pipelines)."""
        return self.results[-1] if self.results else None

    @property
    def all_results(self) -> List[StepResult]:
        """Flat list of all results (trunk + all branches)."""
        out = list(self.results)
        for br in self.branch_results.values():
            out.extend(br.all_results)
        return out

    def branch(self, name: str) -> "PipelineResult":
        """Get results for a named branch."""
        return self.branch_results[name]

    def collect(self, key: str = "energy_Eh") -> Dict[str, Any]:
        """Gather a value from all branches for comparison.

        Returns a dict mapping branch name to the value of ``key``
        from the last result of each branch.

        Example::

            result = pipe.run(cores=8)
            energies = result.collect("energy_Eh")
            # {"oxidation": -1234.56, "reduction": -1234.78}
            best = min(energies, key=energies.get)
        """
        collected = {}
        for bname, br in self.branch_results.items():
            last = br.last
            if last and last.ok:
                collected[bname] = last.data.get(key)
        return collected

    def collect_all(self, *keys: str) -> Dict[str, Dict[str, Any]]:
        """Gather multiple values from all branches.

        Returns ``{branch_name: {key1: val1, key2: val2, ...}}``.

        Example::

            data = result.collect_all("energy_Eh", "dipole_Debye")
            # {"ox": {"energy_Eh": -1234.56, "dipole_Debye": 3.2}, ...}
        """
        collected = {}
        for bname, br in self.branch_results.items():
            last = br.last
            if last and last.ok:
                collected[bname] = {k: last.data.get(k) for k in keys}
        return collected

    def summary(self) -> str:
        """Human-readable summary of the pipeline execution."""
        lines = [f"Pipeline '{self.name}': {'OK' if self.ok else 'FAILED'}"]
        for r in self.results:
            status = "OK" if r.ok else "FAIL"
            elapsed = f"{r.elapsed_seconds:.1f}s" if r.elapsed_seconds else ""
            lines.append(f"  [{status}] {r.step_name} {elapsed}")
            if not r.ok and r.error:
                lines.append(f"         error: {r.error}")
        for bname, br in self.branch_results.items():
            lines.append(f"  branch '{bname}': {'OK' if br.ok else 'FAILED'}")
            for r in br.results:
                status = "OK" if r.ok else "FAIL"
                elapsed = f"{r.elapsed_seconds:.1f}s" if r.elapsed_seconds else ""
                lines.append(f"    [{status}] {r.step_name} {elapsed}")
                if not r.ok and r.error:
                    lines.append(f"           error: {r.error}")
        return "\n".join(lines)


# ======================================================================
#  Provenance logging
# ======================================================================

def _write_provenance(
    base: Path,
    pipeline_name: str,
    result: "PipelineResult",
    cores: int,
    total_seconds: float,
) -> None:
    """Write a JSON provenance log after pipeline execution."""
    import json
    import datetime
    import platform

    def _result_to_record(r: StepResult) -> dict:
        return {
            "step_name": r.step_name,
            "status": r.status.value,
            "geometry_in": None,  # not tracked yet
            "geometry_out": str(r.geometry) if r.geometry else None,
            "work_dir": str(r.work_dir) if r.work_dir else None,
            "data": r.data,
            "artifacts": {k: str(v) for k, v in r.artifacts.items()},
            "error": r.error,
            "elapsed_seconds": round(r.elapsed_seconds, 3),
        }

    prov = {
        "pipeline_name": pipeline_name,
        "ok": result.ok,
        "timestamp": datetime.datetime.now().isoformat(),
        "hostname": platform.node(),
        "python_version": platform.python_version(),
        "cores": cores,
        "total_seconds": round(total_seconds, 3),
        "steps": [_result_to_record(r) for r in result.results],
        "branches": {},
    }
    for bname, br in result.branch_results.items():
        prov["branches"][bname] = {
            "ok": br.ok,
            "steps": [_result_to_record(r) for r in br.results],
        }

    prov_path = base / "provenance.json"
    prov_path.write_text(json.dumps(prov, indent=2, default=str))
    logger.info("Provenance log written: %s", prov_path)


class _PipelineWorkflow:
    """Adapts a Pipeline to the Workflow protocol for registry integration."""

    def __init__(self, pipeline: Pipeline, description: str):
        self.name = pipeline.name
        self.description = description
        self._pipeline = pipeline

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> PipelineResult:
        cores = config.get("PAL", kwargs.get("cores", 1))
        return self._pipeline.run(cores=cores, **kwargs)

    def run_cli(self, argv: List[str]) -> int:
        import argparse
        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("--cores", type=int, default=1, help="CPU cores")
        parser.add_argument("--geometry", type=str, default=None, help="Initial XYZ file")
        parser.add_argument("--work-dir", type=str, default=None, help="Working directory")
        args = parser.parse_args(argv)

        result = self._pipeline.run(
            cores=args.cores,
            geometry=args.geometry,
            work_dir=Path(args.work_dir) if args.work_dir else None,
        )
        print(result.summary())
        return 0 if result.ok else 1


# ======================================================================
#  PipelineTemplate — parametric workflow templates
# ======================================================================

def _resolve_value(value: Any, params: Dict[str, Any]) -> Any:
    """Resolve a template value against provided parameters.

    - ``"{charge}"`` → ``params["charge"]``
    - ``"{charge}+1"`` → ``params["charge"] + 1``
    - ``"{charge}-1"`` → ``params["charge"] - 1``
    - Non-string values pass through unchanged.
    - Strings without ``{`` pass through unchanged.
    """
    if not isinstance(value, str) or "{" not in value:
        return value

    # Simple replacement: "{key}"
    import re
    match = re.fullmatch(r"\{(\w+)\}", value)
    if match:
        key = match.group(1)
        if key in params:
            return params[key]
        return value  # unresolved placeholder

    # Arithmetic: "{key}+N" or "{key}-N"
    match = re.fullmatch(r"\{(\w+)\}\s*([+-])\s*(\d+)", value)
    if match:
        key, op, num = match.group(1), match.group(2), int(match.group(3))
        if key in params:
            base = params[key]
            return base + num if op == "+" else base - num
        return value

    # String interpolation for complex cases: "CPCM({solvent})"
    def _replacer(m):
        k = m.group(1)
        return str(params[k]) if k in params else m.group(0)

    return re.sub(r"\{(\w+)\}", _replacer, value)


def _resolve_kwargs(kwargs: Dict[str, Any], params: Dict[str, Any]) -> Dict[str, Any]:
    """Resolve all template placeholders in a kwargs dict."""
    resolved = {}
    for k, v in kwargs.items():
        if isinstance(v, dict):
            resolved[k] = {dk: _resolve_value(dv, params) for dk, dv in v.items()}
        elif isinstance(v, list):
            resolved[k] = [_resolve_value(item, params) for item in v]
        else:
            resolved[k] = _resolve_value(v, params)
    return resolved


class PipelineTemplate:
    """Parametric pipeline template — define once, apply to any system.

    Template parameters are specified as ``"{param_name}"`` strings in
    step kwargs.  When :meth:`run` is called, all placeholders are
    resolved against the provided keyword arguments.

    Supports arithmetic on integer parameters::

        charge="{charge}+1"  →  charge + 1
        charge="{charge}-1"  →  charge - 1

    Example::

        tpl = PipelineTemplate("opt_freq")
        tpl.add("xtb_opt", charge="{charge}")
        tpl.add("orca_opt", charge="{charge}", method="{method}",
                basis="{basis}", solvent="{solvent}")
        tpl.add("orca_freq", charge="{charge}", method="{method}",
                basis="{basis}", solvent="{solvent}")

        # Run on system A
        results_A = tpl.run(smiles="CCO", charge=0, method="B3LYP",
                            basis="def2-SVP", solvent="water", cores=8)

        # Run on system B — same workflow, different molecule
        results_B = tpl.run(geometry="complex.xyz", charge=2, mult=3,
                            method="PBE0", basis="def2-TZVP",
                            solvent="dmf", cores=16)
    """

    def __init__(self, name: str, *, defaults: Optional[Dict[str, Any]] = None):
        self.name = name
        self._defaults: Dict[str, Any] = dict(defaults or {})
        self._trunk: List[_StepSpec] = []
        self._branches: Dict[str, PipelineTemplate] = {}

    def add(self, step_name: str, *, label: str = "", **kwargs: Any) -> "PipelineTemplate":
        """Add a step to the template trunk (kwargs may contain ``{placeholders}``)."""
        merged = {**self._defaults, **kwargs}
        self._trunk.append(_StepSpec(
            step_name=step_name,
            kwargs=merged,
            label=label or step_name,
        ))
        return self

    def branch(self, name: str) -> "PipelineTemplate":
        """Create a named branch template."""
        if name not in self._branches:
            self._branches[name] = PipelineTemplate(
                f"{self.name}/{name}",
                defaults=self._defaults,  # inherit parent defaults
            )
        return self._branches[name]

    def build(self, **params: Any) -> Pipeline:
        """Resolve all placeholders and return a concrete :class:`Pipeline`."""
        pipe = Pipeline(self.name)

        for spec in self._trunk:
            resolved = _resolve_kwargs(spec.kwargs, params)
            pipe.add(spec.step_name, label=spec.label, **resolved)

        for bname, btpl in self._branches.items():
            child = pipe.branch(bname)
            for spec in btpl._trunk:
                resolved = _resolve_kwargs(spec.kwargs, params)
                child.add(spec.step_name, label=spec.label, **resolved)

        return pipe

    def run(
        self,
        *,
        cores: int = 1,
        geometry: Optional[str | Path] = None,
        work_dir: Optional[Path] = None,
        stop_on_failure: bool = True,
        **params: Any,
    ) -> PipelineResult:
        """Resolve template and execute immediately.

        All keyword arguments not consumed by ``run()`` itself are treated
        as template parameters.
        """
        pipe = self.build(**params)
        return pipe.run(
            cores=cores,
            geometry=geometry,
            work_dir=work_dir,
            stop_on_failure=stop_on_failure,
        )

    def as_workflow(self, description: str = "") -> object:
        """Return a Workflow-protocol-compatible object."""
        return _TemplateWorkflow(self, description or f"Template: {self.name}")

    def __repr__(self) -> str:
        n = len(self._trunk)
        b = len(self._branches)
        parts = [f"PipelineTemplate({self.name!r}, {n} steps"]
        if b:
            parts.append(f", {b} branches")
        parts.append(")")
        return "".join(parts)


class _TemplateWorkflow:
    """Adapts a PipelineTemplate to the Workflow protocol."""

    def __init__(self, template: PipelineTemplate, description: str):
        self.name = template.name
        self.description = description
        self._template = template

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> PipelineResult:
        cores = config.get("PAL", kwargs.pop("cores", 1))
        return self._template.run(cores=cores, **kwargs)

    def run_cli(self, argv: List[str]) -> int:
        import argparse
        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("--cores", type=int, default=1)
        parser.add_argument("--geometry", type=str, default=None)
        parser.add_argument("--work-dir", type=str, default=None)
        # Template params are passed as --key=value
        args, remaining = parser.parse_known_args(argv)
        params: Dict[str, Any] = {}
        for item in remaining:
            if item.startswith("--") and "=" in item:
                k, v = item[2:].split("=", 1)
                # Try int, then float, then string
                try:
                    params[k] = int(v)
                except ValueError:
                    try:
                        params[k] = float(v)
                    except ValueError:
                        params[k] = v

        result = self._template.run(
            cores=args.cores,
            geometry=args.geometry,
            work_dir=Path(args.work_dir) if args.work_dir else None,
            **params,
        )
        print(result.summary())
        return 0 if result.ok else 1
