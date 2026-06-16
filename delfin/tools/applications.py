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
from delfin.tools.templates import classic_opt_freq

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


__all__ = ["opt_freq_energy"]
