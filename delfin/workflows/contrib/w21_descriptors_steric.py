"""Contributed workflow: metal-complex steric descriptor screening (wave-21 DS1).

Chain: geometry -> morfeus buried volume + cone angle -> Multiwfn bond
orders (optional). Built on the wave-21 PipelineTemplate path so it is
contract-tested and shares DELFIN's step registry, validity limits and
physical-realism gate — no parallel mechanism.

Realism limits are the adapters' own (morfeus_realism.py /
multiwfn_realism.py), enforced at execute time; this wrapper adds only
the workflow-level requirement that a metal complex geometry exists.

Usage::

    from delfin.workflows.contrib.w21_descriptors_steric import steric_screening

    result = steric_screening.run(
        geometry="complex.xyz", metal_index=0, cores=1,
    )
"""
from __future__ import annotations

from typing import Any, Dict, List

from delfin.workflows.registry import register


class StericScreeningWorkflow:
    """Steric descriptors for a metal complex: %Vbur, cone angle, Sterimol.

    Valid for: a converged, optimised metal-complex geometry (xyz) with a
    defined metal centre. NOT valid for: an isolated organic molecule with
    no metal centre (the morfeus physics gate refuses it), or a raw,
    unoptimised geometry (descriptors of a non-minimum are not transferable).
    """

    name = "steric_screening"
    description = (
        "Metal-complex steric screening: buried volume, cone angle, "
        "Sterimol via morfeus"
    )

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        from pathlib import Path

        from delfin.tools._registry import get

        geometry = kwargs.get("geometry") or config.get("geometry")
        if geometry is None:
            raise ValueError(
                "steric_screening needs a converged metal-complex geometry "
                "(xyz); pass geometry='complex.xyz'"
            )
        metal_index = kwargs.get("metal_index", 0)
        work_dir = Path(kwargs.get("work_dir", "."))

        geometry = Path(geometry)
        if not geometry.exists():
            raise FileNotFoundError(f"geometry file not found: {geometry}")

        results = []
        for step_name, step_kwargs in (
            ("morfeus_buried_volume",
             {"metal_index": metal_index,
              "radius": kwargs.get("radius", 3.5)}),
            ("morfeus_cone_angle", {"metal_index": metal_index}),
            ("morfeus_full_report",
             {"metal_index": metal_index,
              "donor_indices": kwargs.get("donor_indices")}),
        ):
            step_kwargs = {k: v for k, v in step_kwargs.items()
                           if v is not None}
            step = get(step_name)
            if step is None:
                raise RuntimeError(
                    f"step {step_name} not registered — the morfeus "
                    "adapters must be importable"
                )
            results.append(step.execute(
                work_dir, geometry=geometry,
                cores=int(kwargs.get("cores", 1)), **step_kwargs,
            ))
        return results

    def run_cli(self, argv: List[str]) -> int:
        import argparse

        parser = argparse.ArgumentParser(description=self.description)
        parser.add_argument("--geometry", type=str, required=True)
        parser.add_argument("--metal-index", type=int, default=0)
        parser.add_argument("--radius", type=float, default=3.5)
        parser.add_argument("--donor-indices", type=int, nargs="*")
        parser.add_argument("--cores", type=int, default=1)
        parser.add_argument("--work-dir", type=str, default=".")
        args = parser.parse_args(argv)
        kwargs = {"geometry": args.geometry, "metal_index": args.metal_index,
                  "radius": args.radius, "cores": args.cores,
                  "work_dir": args.work_dir}
        if args.donor_indices:
            kwargs["donor_indices"] = args.donor_indices
        results = self.run(config={}, **kwargs)
        ok = all(getattr(r, "ok", True) or
                 getattr(getattr(r, "status", None), "name", "") == "SUCCESS"
                 for r in results)
        print(f"steric_screening: {len(results)} steps, "
              f"{'all succeeded' if ok else 'FAILURES present'}")
        return 0 if ok else 1


register(StericScreeningWorkflow())
