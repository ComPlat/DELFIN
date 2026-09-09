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


def solvation_keyword_for_method(implicit_solvation_model, solvent, method) -> str:
    """The solvation keyword a given method will actually accept.

    ORCA's xTB does not implement CPCM or SMD.  Handed one it does not warn and
    carry on -- it aborts the whole run:

        WARNING: Found SMD or SMDSolvent or CPCM keyword with XTB calculation.
                 This is not implemented.
        Error (ORCA_MAIN): ... aborting the run

    So an xTB step takes ALPB with the same solvent, which is what the rest of
    DELFIN already writes by hand (``! XTB2 ALPB(DMF)``, ``! XTB2 GOAT
    ALPB(DMF)``).  Everything else keeps the model CONTROL asked for.

    The failure this prevents is total, not subtle: every frame optimisation in
    a solvated MANTA run died at exit code 25 and the run produced nothing.
    """
    name = str(solvent or "").strip()
    if not name:
        return ""
    token = str(method or "").strip().upper()
    if token.startswith("XTB") or token in ("GFN2-XTB", "GFN1-XTB", "GFNFF"):
        return f"ALPB({name})"
    return build_solvation_keyword(implicit_solvation_model, name)
