"""Contributed workflow: CASSCF active-space setup (package CAS1).

A thin shell over delfin.casscf_active_space. The sequence is the
casscf-setup skill's: propose CAS(n,m) with a chemical justification ->
SHOW the proposal to the user -> the CASSCF step starts only after the
user's explicit confirmation (CasscfConfirmationGate; a fresh gate per
workflow run, confirmation is never inherited).
"""

from __future__ import annotations

from typing import Any, Dict

from delfin.workflows.registry import register


class CasscfWorkflow:
    """CASSCF: proposal, validation and occupation checks for the
    active space; the calculation itself starts only on confirmation."""

    name = "casscf"
    description = (
        "CASSCF active-space setup: propose CAS(n,m), validate the "
        "%casscf block, check natural occupations after convergence")

    def propose(self, **kwargs: Any):
        from delfin.casscf_active_space import ActiveSpaceProposal

        proposal = ActiveSpaceProposal(**kwargs)
        proposal.validate()
        return proposal

    @staticmethod
    def show(proposal: Any) -> str:
        """The proposal as shown to the user before any CASSCF run."""
        lines = [
            f"Proposed active space: CAS({proposal.n_electrons},"
            f"{proposal.n_orbitals})",
            f"Justification: {proposal.justification}",
        ]
        if proposal.named_orbital_classes:
            lines.append(
                "Active orbitals from: "
                + ", ".join(proposal.named_orbital_classes))
        if proposal.orbital_rotations:
            lines.append(f"Orbital rotations: {proposal.orbital_rotations}")
        lines.append(
            "Confirm before the CASSCF step starts "
            "(confirmation_gate().confirm(user_confirmed=True)).")
        return "\n".join(lines)

    @staticmethod
    def confirmation_gate():
        from delfin.casscf_active_space import CasscfConfirmationGate

        return CasscfConfirmationGate()

    @staticmethod
    def validate_block(block: Dict[str, Any]) -> bool:
        from delfin.casscf_active_space import validate_casscf_block

        return validate_casscf_block(block)

    @staticmethod
    def check_occupations(occupations):
        from delfin.casscf_active_space import check_active_occupations

        return check_active_occupations(occupations)

    def run(self, *, config: Dict[str, Any], **kwargs: Any) -> Any:
        """Run entry for the registry: refuses to start without the
        confirmation flag in config (the user's explicit confirmation)."""
        if not config.get("CASSCF_user_confirmed"):
            raise PermissionError(
                "the CASSCF step needs the user's explicit confirmation "
                "of the proposal before it starts")
        raise NotImplementedError(
            "actual CASSCF/NEVPT2 submission goes through DELFIN's "
            "ORCA submit path; this workflow provides proposal, "
            "validation and occupation checks")


register(CasscfWorkflow())
