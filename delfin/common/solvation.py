"""Canonical construction of the ORCA implicit-solvation keyword.

Kept in one place on purpose. The same three-line "model, optionally with
solvent" snippet used to be duplicated across the ESD input generator,
OCCUPIER, the IMAG workflow and xyz_io, and every copy carried the same
defect: a model without a solvent was emitted as a bare ``CPCM``. ORCA does
not read that as "no solvation" - it reads it as a conductor with
epsilon = infinity, i.e. an infinitely polarizable medium. That silently
changes the physics of every excited-state energy in the run without
showing up anywhere in the log.
"""

from delfin.common.logging import get_logger

logger = get_logger(__name__)


def build_solvation_keyword(implicit_solvation_model, solvent) -> str:
    """Build the ORCA solvation keyword, or an empty string for gas phase.

    Args:
        implicit_solvation_model: Solvation model (e.g. 'CPCM', 'SMD', 'ALPB')
        solvent: Solvent name (e.g. 'toluene', 'water')

    Returns:
        ``'CPCM(toluene)'`` when both are set, otherwise ``''``.
    """
    model = str(implicit_solvation_model or "").strip()
    solvent_name = str(solvent or "").strip()

    if not model:
        return ""

    if not solvent_name:
        logger.warning(
            "implicit_solvation_model=%s is set but no solvent is given - "
            "dropping the solvation keyword and running in gas phase. "
            "A bare '%s' does NOT mean 'no solvation': ORCA reads it as a "
            "conductor with epsilon = infinity. Set 'solvent=' in CONTROL.txt "
            "if you want an actual solvent.",
            model,
            model,
        )
        return ""

    return f"{model}({solvent_name})"
