"""Solvation-analysis workflow templates."""

from delfin.tools.pipeline import PipelineTemplate


solution_entropy_workflow = PipelineTemplate("solution_entropy_workflow", defaults={
    "charge": "{charge}",
    "mult": "{mult}",
    "method": "{method}",
    "basis": "{basis}",
    "solvent": "{solvent}",
    "temperature": "{temperature}",
    "concentration_standard": "{concentration_standard}",
})
solution_entropy_workflow.add("smiles_to_xyz", smiles="{smiles}", label="smiles")
solution_entropy_workflow.add("xtb_opt", method="XTB2", label="xtb_preopt")
solution_entropy_workflow.add("orca_opt", label="orca_opt")
solution_entropy_workflow.add("orca_freq", label="orca_freq")
solution_entropy_workflow.add("solution_entropy", label="solution_entropy")


reaction_solution_entropy = PipelineTemplate("reaction_solution_entropy", defaults={
    "species": "{species}",
    "stoichiometry": "{stoichiometry}",
    "temperature": "{temperature}",
})
reaction_solution_entropy.add("reaction_solution_entropy", label="reaction_solution_entropy")
