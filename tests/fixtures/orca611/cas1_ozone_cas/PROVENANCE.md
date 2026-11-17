# PROVENANCE — cas1_ozone_cas (synthetic, licence-safe)

The committed fixture `cas1_ozone_cas_synth.out` is SYNTHETIC:
written by the CAS1 builder for DELFIN's tests, not ORCA output.
Values invented to carry the real physics pattern: five
occupations > 1.98 (O 2s/sigma-like), the 1.78/0.26 diradical
pair, two near-empty; two iteration blocks so the LAST-block
(converged) semantics of `parse_casscf_occupations` is tested.

The real CASSCF output this pattern is modeled on was produced by
the practice run of the CAS1 workflow (SLURM job 8178877, CAS(18,12),
E(CAS)=-224.298002 Eh; UNO reference job 8178848). On that real
output the check flags exactly the five near-2.0 orbitals and
proposes shrinking (CAS(12,9), the usual full-pi + lone-pair
choice). ORCA EULA (June 2025, sections 3c-3e) restricts
redistribution of ORCA-generated data, so the real files are NOT
committed; they live under the private practice area (.gate/) of
the CAS1 worktree, and the operator holds a private copy.

Tests: the synthetic files are always run; the real-output checks
skip unless DELFIN_PRIVATE_FIXTURES points at the folder holding
the real files.
