"""Flagship composite templates — the richest shipped workflows.

These compose the derived building blocks into complete, named workflows that
the platform registers as Applications (typed inputs/outputs, shown in the
Pipelines tab):

* :data:`xtb_thermochemistry` — a fully open-source path (SMILES → native xTB
  opt → xTB Hessian/thermo). Runnable end-to-end without any licensed engine.
* :data:`conformer_energy` — SMILES → xTB pre-opt → CREST conformer search →
  xTB refinement of the best conformer.
* :data:`full_workflow` — a CONTROL.txt-class reference: SMILES → xTB pre-opt →
  DFT opt → DFT freq → imaginary-frequency cleanup → large-basis single point.
  Demonstrates that the foundation can express a workflow of the same complexity
  as DELFIN's classic engine, purely by composing derived blocks (the real
  CONTROL.txt-driven workflow is never touched).
"""

from delfin.tools.pipeline import PipelineTemplate

# ── xTB thermochemistry: SMILES → xTB opt → xTB Hessian (open-source) ──

xtb_thermochemistry = PipelineTemplate("xtb_thermochemistry", defaults={
    "charge": "{charge}",
    "mult": "{mult}",
    "solvent": "{solvent}",
})
xtb_thermochemistry.add("smiles_to_xyz", smiles="{smiles}", label="smiles")
xtb_thermochemistry.add("xtb_optimize", method="{method}", label="xtb_opt")
xtb_thermochemistry.add("xtb_hessian", method="{method}", label="xtb_thermo")


# ── Conformer energy: SMILES → xTB pre-opt → CREST → xTB refine ───────

conformer_energy = PipelineTemplate("conformer_energy", defaults={
    "charge": "{charge}",
    "mult": "{mult}",
    "solvent": "{solvent}",
})
conformer_energy.add("smiles_to_xyz", smiles="{smiles}", label="smiles")
conformer_energy.add("xtb_optimize", method="{method}", label="xtb_preopt")
conformer_energy.add("crest_conformers", label="crest")
conformer_energy.add("xtb_optimize", method="{method}", label="refine_best")


# ── full_workflow: CONTROL.txt-class reference (classic engine shape) ──

full_workflow = PipelineTemplate("full_workflow", defaults={
    "charge": "{charge}",
    "mult": "{mult}",
    "method": "{method}",
    "basis": "{basis}",
    "solvent": "{solvent}",
    # legacy IMAG param vocabulary, broadcast for the imag_fix stage (ignored by
    # the ORCA/xTB steps, which take only the kwargs they declare):
    "metals": "{metals}",
    "main_basisset": "{main_basisset}",
    "metal_basisset": "{metal_basisset}",
})
# 1. SMILES → 3D geometry
full_workflow.add("smiles_to_xyz", smiles="{smiles}", label="smiles")
# 2. fast xTB pre-optimization (pin the xTB method, do not inherit the DFT one)
full_workflow.add("xtb_opt", method="XTB2", label="xtb_preopt")
# 3. DFT geometry optimization
full_workflow.add("orca_opt", label="dft_opt")
# 4. DFT frequencies (produces the hessian that step 5 auto-wires)
full_workflow.add("orca_freq", label="dft_freq")
# 5. eliminate imaginary frequencies (hess_file auto-wired from dft_freq)
full_workflow.add("imag_fix", label="imag_cleanup")
# 6. large-basis single-point refinement
full_workflow.add("orca_sp", basis="{refine_basis}", label="refine")
