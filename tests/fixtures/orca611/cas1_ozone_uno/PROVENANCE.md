# PROVENANCE — cas1_ozone_uno (synthetic, licence-safe)

The committed fixture `cas1_ozone_ref_synth.out` is SYNTHETIC: written
by the CAS1 builder for DELFIN's tests, not ORCA output. First line
says so. The parser needs only the `UHF NATURAL ORBITALS` /
`Occupation Numbers:` header and `N[i]= value` entries; those lines
follow ORCA's documented print format (manual 3.1.2, 5.1.11), the
values are invented to carry the physics pattern (pi pair ~1.8/0.2,
cores 2.0).

The real ORCA reference output this pattern is modeled on was
produced by practice runs of the CAS1 workflow (SLURM jobs 8178848
reference, 8178877 CASSCF step, both completed normally; ozone
geometry NIST CCCBDB experimental, Herzberg 1966: r(OO)=1.278 A,
angle 116.8 deg). ORCA EULA (June 2025, sections 3c-3e) restricts
redistribution of ORCA-generated data, so the real files are NOT
committed; they live under the private practice area (.gate/) of
the CAS1 worktree, and the operator holds a private copy.

Tests: the synthetic files are always run; the real-output checks
skip unless the environment variable DELFIN_PRIVATE_FIXTURES points
at the folder holding the real files.
