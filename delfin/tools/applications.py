"""Built-in applications — workflows exposed as callable functions.

These derive Applications from the shipped pipeline templates, giving the
platform's application registry real, demonstrable content.  Outputs are mapped
to truthful values the underlying steps actually produce (energies), not derived
properties the template does not compute.  Define your own with
:meth:`Application.from_pipeline` and :func:`register_application`.
"""

from __future__ import annotations

from delfin.tools._application import Application, OutputSpec, register_application
from delfin.tools._spec import ParamSpec
from delfin.tools._keys import key
from delfin.tools.templates import classic_opt_freq, multi_level_opt, redox_potential

# SMILES → xTB pre-opt → ORCA opt → ORCA freq, returning the final energies.
opt_freq_energy = Application.from_pipeline(
    classic_opt_freq,
    name="opt_freq_energy",
    description="Optimise a molecule (xTB pre-opt → ORCA opt) and run a frequency "
                "job; returns the final electronic and Gibbs free energies.",
    category="dft",
    inputs=(
        ParamSpec("smiles", "str", required=True, description="Input SMILES string"),
        ParamSpec("charge", "int", required=True, description="Molecular charge"),
        ParamSpec("method", "str", default="B3LYP", description="DFT functional"),
        ParamSpec("basis", "str", default="def2-SVP", description="Basis set"),
    ),
    outputs=(
        OutputSpec("energy_Eh", step="orca_freq", key="energy_Eh",
                   type="float", unit="Eh", description="Final electronic energy"),
        OutputSpec("gibbs_Eh", step="orca_freq", key="gibbs_Eh",
                   type="float", unit="Eh", description="Gibbs free energy"),
    ),
)
register_application(opt_freq_energy)


# Redox: optimise + freq for the neutral, oxidized (+1) and reduced (−1) species
# in parallel branches; returns the three Gibbs free energies (the raw redox
# workflow output — the potential in V is derived from these via DELFIN's
# convention downstream).
redox = Application.from_pipeline(
    redox_potential,
    name="redox_potential",
    description="Redox workflow: neutral + oxidized + reduced opt/freq; returns the "
                "Gibbs free energies of all three species.",
    category="dft",
    inputs=(
        ParamSpec("smiles", "str", required=True, description="Input SMILES string"),
        key("charge", required=True),
        key("method"),
        key("basis"),
        ParamSpec("mult_ox", "int", default=2, description="Multiplicity of the oxidized species"),
        ParamSpec("mult_red", "int", default=2, description="Multiplicity of the reduced species"),
    ),
    outputs=(
        OutputSpec("e_neutral_Eh", step="neutral_freq", key="gibbs_Eh",
                   unit="Eh", description="Neutral Gibbs free energy"),
        OutputSpec("e_oxidation_Eh", step="ox_freq", branch="oxidation", key="gibbs_Eh",
                   unit="Eh", description="Oxidized Gibbs free energy"),
        OutputSpec("e_reduction_Eh", step="red_freq", branch="reduction", key="gibbs_Eh",
                   unit="Eh", description="Reduced Gibbs free energy"),
    ),
)
register_application(redox)


# Multi-level: SMILES → xTB pre-opt → small-basis DFT opt → large-basis DFT opt →
# freq; returns the final energy/Gibbs at the large basis.
multi_level_energy = Application.from_pipeline(
    multi_level_opt,
    name="multi_level_energy",
    description="Basis-set ladder: xTB → small-basis DFT opt → large-basis DFT opt → "
                "freq; returns the final energies at the large basis.",
    category="dft",
    inputs=(
        ParamSpec("smiles", "str", required=True, description="Input SMILES string"),
        key("charge", required=True),
        key("method"),
        ParamSpec("small_basis", "str", default="def2-SVP", description="Small (opt) basis"),
        ParamSpec("large_basis", "str", default="def2-TZVP", description="Large (refine) basis"),
    ),
    outputs=(
        OutputSpec("energy_Eh", step="freq", key="energy_Eh", unit="Eh",
                   description="Electronic energy at the large basis"),
        OutputSpec("gibbs_Eh", step="freq", key="gibbs_Eh", unit="Eh",
                   description="Gibbs free energy at the large basis"),
    ),
)
register_application(multi_level_energy)


__all__ = ["opt_freq_energy", "redox", "multi_level_energy"]
