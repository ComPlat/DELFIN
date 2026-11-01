# PROVENANCE — cas1_n2_2p5_cas (synthetic, licence-safe)

The committed fixture `cas1_n2_2p5_cas_synth.out` is SYNTHETIC:
written by the CAS1 builder for DELFIN's tests, not ORCA output.
Values invented to carry the real physics pattern: all six active
NOs clearly partial (1.273..0.727, healthy CAS(6,6)), two
iteration blocks so LAST-block (converged) semantics is tested.

The real CASSCF output this pattern is modeled on was produced by
the practice run of the CAS1 workflow (SLURM job 8178878, CAS(6,6),
E(CAS)=-108.666200 Eh; UNO reference job 8178849). On the real
output all six NOs are partial — no resize needed. ORCA EULA
(June 2025, sections 3c-3e) restricts redistribution of
ORCA-generated data, so the real files are NOT committed; they
live under the private practice area (.gate/) of the CAS1
worktree, and the operator holds a private copy.

Tests: the synthetic files are always run; the real-output checks
skip unless DELFIN_PRIVATE_FIXTURES points at the folder holding
the real files.
