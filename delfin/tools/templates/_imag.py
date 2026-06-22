"""IMAG elimination template: opt → freq → fix imaginary loop.

Mirrors the behavior of DELFIN's IMAG workflow: detect imaginary
frequencies and iteratively displace + re-optimize until clean.

Usage::

    from delfin.tools.templates import imag_elimination

    result = imag_elimination.run(
        geometry="optimized.xyz",
        charge=0, method="B3LYP", basis="def2-SVP",
        max_imag_iter=5,
        cores=8,
    )

As a composable sub-pipeline::

    from delfin.tools.templates import imag_sub_pipeline
    from delfin.tools import Pipeline

    main = Pipeline("full_workflow")
    main.add("smiles_to_xyz", smiles="CCO")
    main.add("xtb_opt", charge=0)
    main.add("orca_opt", charge=0, method="B3LYP", basis="def2-SVP")
    main.add_sub_pipeline(
        imag_sub_pipeline(charge=0, method="B3LYP", basis="def2-SVP", max_iter=5),
        label="imag_cleanup",
    )
    main.add("orca_sp", charge=0, method="B3LYP", basis="def2-TZVP")
"""

from delfin.tools.pipeline import Pipeline, PipelineTemplate


def imag_sub_pipeline(
    *,
    charge: int = 0,
    mult: int = 1,
    method: str = "B3LYP",
    basis: str = "def2-SVP",
    max_iter: int = 5,
    **orca_kwargs,
) -> Pipeline:
    """Build a reusable IMAG-fix sub-pipeline with concrete parameters.

    Returns a Pipeline that can be used with ``add_sub_pipeline()``.
    """
    pipe = Pipeline("imag_fix", defaults={
        "charge": charge, "mult": mult,
        "method": method, "basis": basis,
        **orca_kwargs,
    })

    # Frequency check
    pipe.add("orca_freq", label="freq_check")

    # IMAG fix loop: repeat until no imaginary frequencies
    pipe.add_loop(
        "imag_fix",
        until=lambda result, iteration: result.data.get("n_imaginary", 0) == 0,
        max_iter=max_iter,
        label="imag_loop",
    )

    return pipe


# ── Template version (parametric) ──────────────────────────────────────

imag_elimination = PipelineTemplate("imag_elimination", defaults={
    "charge": "{charge}",
    "method": "{method}",
    "basis": "{basis}",
})

# Step 1: ORCA frequency calculation
imag_elimination.add("orca_freq", label="freq_check")

# Note: The template version cannot use add_loop (needs concrete callables).
# For the loop version, use imag_sub_pipeline() which returns a concrete
# Pipeline. The template is for simple freq-then-fix workflows.
imag_elimination.add("imag_fix", label="imag_fix")
