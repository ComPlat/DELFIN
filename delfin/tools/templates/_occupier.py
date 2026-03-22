"""OCCUPIER-style staged optimization using sub-pipelines and reactive steps.

The OCCUPIER workflow gradually relaxes a structure through stages (FoB =
fraction of bonds), where each stage runs a full opt+freq sub-pipeline.
This template recreates that pattern using ``add_sub_pipeline()`` and
``add_reactive()`` for contamination-triggered corrections.

Usage::

    from delfin.tools.templates import occupier_stages

    # Simple staged optimization
    pipe = occupier_stages(
        charge=0, mult=1,
        method="B3LYP", basis="def2-SVP",
        fob_values=[0.0, 0.25, 0.5, 0.75, 1.0],
        cores=8,
    )
    result = pipe.run(geometry="start.xyz", cores=8)

    # With contamination check (broken-symmetry correction)
    pipe = occupier_stages(
        charge=2, mult=5,
        method="B3LYP", basis="def2-SVP",
        fob_values=[0.0, 0.5, 1.0],
        check_contamination=True,
        contamination_threshold=0.1,
    )
    result = pipe.run(geometry="metal_complex.xyz", cores=16)
"""

from __future__ import annotations

from typing import List, Optional

from delfin.tools.pipeline import Pipeline


def _make_fob_stage(
    fob: float,
    *,
    charge: int,
    mult: int = 1,
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    **orca_kwargs,
) -> Pipeline:
    """Build a sub-pipeline for a single FoB stage."""
    stage = Pipeline(f"fob_{fob:.2f}", defaults={
        "charge": charge, "mult": mult,
        "method": method, "basis": basis,
        **orca_kwargs,
    })
    stage.add("orca_opt", label=f"opt_fob{fob:.2f}")
    stage.add("orca_freq", label=f"freq_fob{fob:.2f}")
    return stage


def occupier_stages(
    *,
    charge: int,
    mult: int = 1,
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    fob_values: Optional[List[float]] = None,
    check_contamination: bool = False,
    contamination_threshold: float = 0.1,
    **orca_kwargs,
) -> Pipeline:
    """Build an OCCUPIER-style staged optimization pipeline.

    Parameters
    ----------
    charge, mult : molecular charge and spin multiplicity
    method, basis : DFT method and basis set
    fob_values : list of FoB (fraction of bonds) values to iterate
                 (default: [0.0, 0.25, 0.5, 0.75, 1.0])
    check_contamination : if True, add reactive contamination checks
                          after each stage
    contamination_threshold : S² deviation threshold for triggering
                              broken-symmetry correction
    **orca_kwargs : additional ORCA parameters (ri, dispersion, solvent, ...)

    Returns
    -------
    Pipeline ready to run with ``pipe.run(geometry=..., cores=...)``
    """
    if fob_values is None:
        fob_values = [0.0, 0.25, 0.5, 0.75, 1.0]

    pipe = Pipeline("occupier", defaults={
        "charge": charge, "mult": mult,
        "method": method, "basis": basis,
        **orca_kwargs,
    })

    for fob in fob_values:
        # Each FoB stage is a sub-pipeline
        stage = _make_fob_stage(
            fob, charge=charge, mult=mult,
            method=method, basis=basis, **orca_kwargs,
        )
        pipe.add_sub_pipeline(stage, label=f"fob_{fob:.2f}")

        # Optional: check for spin contamination and add BS correction
        if check_contamination:
            threshold = contamination_threshold

            def _check_contam(results, last, inserter, _fob=fob, _thr=threshold):
                if last is None or not last.ok:
                    return
                s2_dev = last.data.get("s2_deviation", 0)
                if abs(s2_dev) > _thr:
                    inserter.add(
                        "orca_sp",
                        charge=charge, mult=mult,
                        method=method, basis=basis,
                        broken_sym=True,
                        label=f"bs_fob{_fob:.2f}",
                        **orca_kwargs,
                    )

            pipe.add_reactive(_check_contam, label=f"contam_check_{fob:.2f}")

    return pipe
