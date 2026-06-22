"""Pre-built workflow templates for common computational chemistry tasks.

These templates demonstrate how to use the Pipeline API to build
workflows equivalent to DELFIN's built-in engines (Classic, IMAG,
OCCUPIER-style) while being fully customizable.

Quick start::

    from delfin.tools.templates import classic_opt_freq, imag_elimination

    # Run a classic optimization + frequency workflow
    result = classic_opt_freq.run(
        smiles="CCO", charge=0, method="B3LYP", basis="def2-SVP",
        cores=8,
    )
    print(result.summary())

    # Customize: add dispersion, solvent, RI
    result = classic_opt_freq.run(
        smiles="[Fe+2](N)(N)(N)(N)(N)N",
        charge=2, mult=5,
        method="B3LYP", basis="def2-SVP",
        ri="RIJCOSX", aux_basis="def2/J",
        dispersion="D4", solvent="water", solvent_model="CPCM",
        cores=16,
    )

Available templates:

- :data:`classic_opt_freq` — SMILES → xTB pre-opt → ORCA opt → freq (Classic workflow)
- :data:`imag_elimination` — ORCA opt → freq → IMAG fix loop (IMAG workflow)
- :data:`redox_potential` — Neutral + oxidized + reduced → redox potentials
- :data:`conformer_screening` — CREST conformers → ORCA SP ranking
- :data:`multi_level_opt` — xTB → small basis → large basis (refinement)
- :func:`occupier_stages` — Build OCCUPIER-style staged optimization
- :func:`esd_states` — Build ESD-style excited state pipeline
"""

from delfin.tools.templates._classic import classic_opt_freq
from delfin.tools.templates._imag import imag_elimination, imag_sub_pipeline
from delfin.tools.templates._redox import redox_potential
from delfin.tools.templates._screening import conformer_screening, multi_level_opt
from delfin.tools.templates._occupier import occupier_stages
from delfin.tools.templates._esd import esd_states

__all__ = [
    "classic_opt_freq",
    "imag_elimination",
    "imag_sub_pipeline",
    "redox_potential",
    "conformer_screening",
    "multi_level_opt",
    "occupier_stages",
    "esd_states",
]
