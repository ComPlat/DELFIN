"""OPI-driven multi-step workflow template: r2SCAN-3c opt+freq → DLPNO SP.

Composes the OPI-backed ORCA building blocks that feed each other through
the Pipeline engine's geometry / artifact wiring::

    orca_opt_opi   -> optimized geometry + .gbw
    orca_freq_opi  -> verifies the minimum, gives thermochemistry
    orca_dlpno_opi -> DLPNO-CCSD(T) single point on the optimized geometry,
                      started from the chain's orbitals (MOREAD)

Usage::

    from delfin.tools.templates import opi_opt_freq_dlpno

    result = opi_opt_freq_dlpno.run(
        geometry="input.xyz", charge=0,
        dft_method="r2SCAN-3c", dft_basis="def2-SVP",
        wft_basis="def2-TZVPP", cores=8,
    )
"""

from delfin.tools.pipeline import PipelineTemplate

# ── OPI chain: DFT opt+freq → DLPNO-CCSD(T) single point ──────────────

opi_opt_freq_dlpno = PipelineTemplate("opi_opt_freq_dlpno", defaults={
    "charge": "{charge}",
    "dft_method": "{dft_method}",
    "dft_basis": "{dft_basis}",
})

# Step 1: 3D geometry from SMILES (skipped if a geometry is provided).
opi_opt_freq_dlpno.add("smiles_to_xyz", smiles="{smiles}", label="smiles")

# Step 2: DFT optimization with the OPI-backed block. Its StepResult carries
# the optimized geometry and the .gbw, which the engine wires forward.
opi_opt_freq_dlpno.add("orca_opt_opi", label="dft_opt")

# Step 3: frequencies verify the minimum and give the thermal correction.
opi_opt_freq_dlpno.add(
    "orca_freq_opi", label="dft_freq")

# Step 4: DLPNO-CCSD(T) single point on the optimized geometry. The engine
# injects the previous step's .gbw as ``moread`` start orbitals; the DFT
# ``method`` broadcast must not leak into the WFT step.
opi_opt_freq_dlpno.add("orca_dlpno_opi", label="dlpno_sp", method="")
