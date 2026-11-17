# PROVENANCE — cas1_n2_2p5_uno (synthetic, licence-safe)

The committed fixture `cas1_n2_2p5_ref_synth.out` is SYNTHETIC:
written by the CAS1 builder for DELFIN's tests, not ORCA output.
Values invented to carry the physics pattern: exactly six
fractional UNO occupations (indices 4-9, 1.194..0.806), the rest
2.0/0.0 — the shape that justifies CAS(6,6) (manual 3.16.1's own
N2 dissociation example).

The real ORCA output this pattern is modeled on comes from the
CAS1 practice runs (SLURM jobs 8178849 reference, 8178878 CASSCF
step, both completed normally; N2 at r(NN)=2.5 A). ORCA EULA
(June 2025, sections 3c-3e) restricts redistribution of
ORCA-generated data, so the real files are NOT committed; they
live under the private practice area (.gate/) of the CAS1
worktree, and the operator holds a private copy.

Tests: the synthetic files are always run; the real-output checks
skip unless DELFIN_PRIVATE_FIXTURES points at the folder holding
the real files.
