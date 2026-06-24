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
from delfin.tools.templates import (
    classic_opt_freq,
    conformer_dG,
    conformer_energy,
    full_workflow,
    multi_level_opt,
    redox_potential,
    xtb_thermochemistry,
)

# A method input for the native-xTB workflows (GFN level, not a DFT functional).
_XTB_METHOD = ParamSpec("method", "str", default="gfn2",
                        enum=("gfn0", "gfn1", "gfn2", "gfnff"),
                        description="xTB level (GFN0/1/2-xTB or GFN-FF)")
_SOLVENT_OPT = ParamSpec("solvent", "str", default="",
                         description="Implicit solvent name (empty = gas phase)")

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


# Native xTB thermochemistry — a fully open-source path (no licensed engine):
# SMILES → xTB optimization → xTB Hessian/thermo. Runnable end-to-end.
xtb_thermo_app = Application.from_pipeline(
    xtb_thermochemistry,
    name="xtb_thermochemistry",
    description="Open-source thermochemistry: SMILES → native xTB optimization → "
                "xTB Hessian; returns energy, free energy, ZPVE and the imaginary "
                "frequency count. Needs only the free xtb binary.",
    category="semiempirical",
    inputs=(
        ParamSpec("smiles", "str", required=True, description="Input SMILES string"),
        key("charge", required=True),
        key("mult"),
        _XTB_METHOD,
        _SOLVENT_OPT,
    ),
    outputs=(
        OutputSpec("energy_Eh", step="xtb_thermo", key="energy_Eh",
                   unit="Eh", description="Total xTB energy"),
        OutputSpec("free_energy_Eh", step="xtb_thermo", key="free_energy_Eh",
                   unit="Eh", description="Total free energy (RRHO)"),
        OutputSpec("zpve_Eh", step="xtb_thermo", key="zpve_Eh",
                   unit="Eh", description="Zero-point vibrational energy"),
        OutputSpec("n_imaginary", step="xtb_thermo", key="n_imaginary",
                   type="int", description="Number of imaginary frequencies"),
    ),
)
register_application(xtb_thermo_app)


# Conformer energy — SMILES → xTB pre-opt → CREST conformer search → xTB
# refinement of the best conformer; returns its energy.
conformer_energy_app = Application.from_pipeline(
    conformer_energy,
    name="conformer_energy",
    description="Conformer search: SMILES → xTB pre-opt → CREST ensemble → xTB "
                "refinement of the best conformer; returns its energy.",
    category="semiempirical",
    inputs=(
        ParamSpec("smiles", "str", required=True, description="Input SMILES string"),
        key("charge", required=True),
        key("mult"),
        _XTB_METHOD,
        _SOLVENT_OPT,
    ),
    outputs=(
        OutputSpec("energy_Eh", step="refine_best", key="energy_Eh",
                   unit="Eh", description="Best-conformer xTB energy"),
    ),
)
register_application(conformer_energy_app)


# Conformer ΔG ranking — SMILES → xTB pre-opt → CREST ensemble → CENSO
# free-energy sorting (the ensemble auto-wires from CREST). Returns the conformer
# count and the directory holding the CENSO-ranked ensemble.
conformer_dG_app = Application.from_pipeline(
    conformer_dG,
    name="conformer_dG",
    description="Conformer free-energy ranking: SMILES → xTB pre-opt → CREST "
                "ensemble → CENSO ΔG sorting; returns the conformer count and the "
                "ranked-ensemble directory.",
    category="semiempirical",
    inputs=(
        ParamSpec("smiles", "str", required=True, description="Input SMILES string"),
        key("charge", required=True),
        key("mult"),
        _SOLVENT_OPT,
        ParamSpec("crest_method", "str", default="gfn2",
                  enum=("gfn2", "gfn1", "gfnff", "gfn2//gfnff"),
                  description="GFN level for the CREST search"),
        ParamSpec("functional", "str", default="r2scan-3c",
                  description="CENSO functional for the ΔG ranking"),
    ),
    outputs=(
        OutputSpec("n_conformers", step="crest", key="n_conformers",
                   type="int", description="Number of conformers found"),
        OutputSpec("ranked_dir", step="censo", key="output_dir",
                   type="str", description="Directory with the CENSO-ranked ensemble"),
    ),
)
register_application(conformer_dG_app)


# full_workflow — CONTROL.txt-class reference: SMILES → xTB pre-opt → DFT opt →
# DFT freq → imaginary-frequency cleanup → large-basis single point. Composes the
# derived blocks into one named workflow of classic-engine complexity (the real
# CONTROL.txt-driven workflow is untouched). hess_file for the imag_cleanup stage
# auto-wires from the dft_freq hessian.
full_workflow_app = Application.from_pipeline(
    full_workflow,
    name="full_workflow",
    description="CONTROL.txt-class reference workflow: SMILES → xTB pre-opt → DFT "
                "opt → DFT freq → imaginary-frequency cleanup → large-basis single "
                "point. Returns the refined electronic energy and the Gibbs free "
                "energy.",
    category="dft",
    inputs=(
        ParamSpec("smiles", "str", required=True, description="Input SMILES string"),
        key("charge", required=True),
        key("mult"),
        key("method"),
        ParamSpec("basis", "str", default="def2-SVP", description="Opt/freq basis set"),
        ParamSpec("refine_basis", "str", default="def2-TZVP",
                  description="Large basis for the final single point"),
        key("solvent", default=""),
        ParamSpec("metals", "list", default=[],
                  description="Metal centers for the imag-cleanup stage (usually [])"),
        ParamSpec("main_basisset", "str", default="def2-SVP",
                  description="imag_fix main basis set"),
        ParamSpec("metal_basisset", "str", default="def2-TZVP",
                  description="imag_fix metal basis set"),
    ),
    outputs=(
        OutputSpec("energy_Eh", step="refine", key="energy_Eh",
                   unit="Eh", description="Refined electronic energy (large basis)"),
        OutputSpec("gibbs_Eh", step="dft_freq", key="gibbs_Eh",
                   unit="Eh", description="Gibbs free energy (opt-level frequencies)"),
    ),
)
register_application(full_workflow_app)


__all__ = [
    "opt_freq_energy",
    "redox",
    "multi_level_energy",
    "xtb_thermo_app",
    "conformer_energy_app",
    "conformer_dG_app",
    "full_workflow_app",
]
