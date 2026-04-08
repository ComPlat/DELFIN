"""Screening templates: conformer ranking and multi-level optimization.

Conformer screening::

    from delfin.tools.templates import conformer_screening

    result = conformer_screening.run(
        smiles="CCO", charge=0,
        method="B3LYP", basis="def2-SVP",
        cores=16,
    )

Multi-level optimization (basis set ladder)::

    from delfin.tools.templates import multi_level_opt

    result = multi_level_opt.run(
        smiles="CCO", charge=0,
        method="B3LYP",
        small_basis="def2-SVP",
        large_basis="def2-TZVP",
        cores=8,
    )
"""

from delfin.tools.pipeline import PipelineTemplate

# ── Conformer screening: CREST → ORCA SP ranking ──────────────────────

conformer_screening = PipelineTemplate("conformer_screening", defaults={
    "charge": "{charge}",
    "method": "{method}",
    "basis": "{basis}",
})

# Step 1: Generate 3D from SMILES
conformer_screening.add("smiles_to_xyz", smiles="{smiles}", label="smiles")

# Step 2: CREST conformer search
conformer_screening.add("crest_conformers", label="crest")

# Step 3: Rank conformers with ORCA SP (best energy wins)
# Uses fan_out in the template — users run via build() + add_fan_out


# ── Multi-level optimization: xTB → small DFT → large DFT ─────────────

multi_level_opt = PipelineTemplate("multi_level_opt", defaults={
    "charge": "{charge}",
    "method": "{method}",
})

multi_level_opt.add("smiles_to_xyz", smiles="{smiles}", label="smiles")
multi_level_opt.add("xtb_opt", label="xtb_preopt")
multi_level_opt.add("orca_opt", basis="{small_basis}", label="small_basis_opt")
multi_level_opt.add("orca_opt", basis="{large_basis}", label="large_basis_opt")
multi_level_opt.add("orca_freq", basis="{large_basis}", label="freq")
