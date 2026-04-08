"""Classic workflow template: SMILES → xTB pre-opt → ORCA opt → freq.

Mirrors the behavior of DELFIN's built-in Classic workflow engine,
but as a composable PipelineTemplate.

Usage::

    from delfin.tools.templates import classic_opt_freq

    result = classic_opt_freq.run(
        smiles="CCO", charge=0, method="B3LYP", basis="def2-SVP",
        cores=8,
    )

    # With all options
    result = classic_opt_freq.run(
        smiles="[Fe+2](N)(N)(N)(N)(N)N",
        charge=2, mult=5,
        method="B3LYP", basis="def2-SVP",
        ri="RIJCOSX", aux_basis="def2/J",
        dispersion="D4",
        solvent="water", solvent_model="CPCM",
        metal_basis={"Fe": "def2-TZVP"},
        cores=16,
    )

    # From existing geometry (skip SMILES step)
    result = classic_opt_freq.run(
        geometry="input.xyz",
        charge=0, method="B3LYP", basis="def2-SVP",
        cores=4,
    )
"""

from delfin.tools.pipeline import PipelineTemplate

# ── Classic: SMILES → xTB → ORCA opt → ORCA freq ──────────────────────

classic_opt_freq = PipelineTemplate("classic_opt_freq", defaults={
    "charge": "{charge}",
    "method": "{method}",
    "basis": "{basis}",
})

# Step 1: SMILES → 3D geometry (skipped if geometry provided directly)
classic_opt_freq.add("smiles_to_xyz", smiles="{smiles}", label="smiles")

# Step 2: xTB pre-optimization (fast, gets close to minimum)
classic_opt_freq.add("xtb_opt", label="xtb_preopt")

# Step 3: ORCA geometry optimization (DFT)
classic_opt_freq.add("orca_opt", label="orca_opt")

# Step 4: ORCA frequency calculation (verify minimum)
classic_opt_freq.add("orca_freq", label="orca_freq")
