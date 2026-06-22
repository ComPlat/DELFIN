"""Redox potential template: neutral + oxidized + reduced branches.

Computes oxidation and reduction potentials by running geometry
optimizations and frequency calculations for the neutral, oxidized,
and reduced species in parallel branches.

Usage::

    from delfin.tools.templates import redox_potential

    result = redox_potential.run(
        smiles="[Fe+2](N)(N)(N)(N)(N)N",
        charge=2, mult=5,
        mult_ox=4, mult_red=6,
        method="B3LYP", basis="def2-SVP",
        cores=16,
    )

    # Collect energies across branches
    energies = result.collect("electronic_energy_Eh")
    print(energies)
    # {"oxidation": -1234.56, "reduction": -1234.78}
"""

from delfin.tools.pipeline import PipelineTemplate

# ── Redox potential: trunk + oxidation/reduction branches ──────────────

redox_potential = PipelineTemplate("redox_potential", defaults={
    "charge": "{charge}",
    "method": "{method}",
    "basis": "{basis}",
})

# Trunk: SMILES → xTB pre-opt → ORCA opt + freq (neutral species)
redox_potential.add("smiles_to_xyz", smiles="{smiles}", label="smiles")
redox_potential.add("xtb_opt", label="xtb_preopt")
redox_potential.add("orca_opt", label="neutral_opt")
redox_potential.add("orca_freq", label="neutral_freq")

# Branch: oxidation (+1 charge)
ox = redox_potential.branch("oxidation")
ox.add("orca_opt", charge="{charge}+1", mult="{mult_ox}", label="ox_opt")
ox.add("orca_freq", charge="{charge}+1", mult="{mult_ox}", label="ox_freq")

# Branch: reduction (-1 charge)
red = redox_potential.branch("reduction")
red.add("orca_opt", charge="{charge}-1", mult="{mult_red}", label="red_opt")
red.add("orca_freq", charge="{charge}-1", mult="{mult_red}", label="red_freq")
