---
name: descriptors-steric
description: >
  Compute steric descriptors of a metal complex (buried volume %Vbur,
  cone angle, Sterimol, full report) through DELFIN's step registry with
  physically enforced validity limits.
domains: [descriptors, workflows]
---

# Steric descriptors of a metal complex (morfeus)

## What this skill does

Runs the contributed workflow `steric_screening` (registered in DELFIN's
workflow registry like every contrib workflow): geometry →
`morfeus_buried_volume` → `morfeus_cone_angle` → `morfeus_full_report`,
each through the adapters' execute-time physics gate.

## When to use it

- You have a converged, optimised geometry (xyz) of a metal complex with
  a defined metal centre.
- You need %Vbur, Tolman cone angle, bite angle (with 2 donor indices)
  or Sterimol parameters.

## When NOT to use it

- The molecule has no metal centre: the physics gate refuses buried
  volume and cone angle (no ligand sphere around a metal).
- The geometry is unoptimised: descriptors of a non-minimum are not
  transferable — optimise first (classic_opt_freq template).
- Radius outside 1–10 Å: refused (physically sensible ligand-shell range).

## How to run it (workflow path, preferred)

    from delfin.workflows.registry import get
    wf = get("steric_screening")
    results = wf.run(config={}, geometry="complex.xyz", metal_index=0,
                     donor_indices=[1, 2])  # donors enable the bite angle

CLI:

    delfin-agent workflow run steric_screening --geometry complex.xyz \
        --metal-index 0 --donor-indices 1 2

## Rules the user cannot override

- Atom indices are 0-based and must be >= 0; a negative index is refused
  before any compute is spent.
- Sphere radius must be within [1, 10] Å.
- Bite angle needs exactly 2 donor indices.
- If morfeus itself rejects an index (e.g. it points past the last atom),
  the step fails with a named physics error, not a bare crash.

## Sources

- %Vbur / cone angle / Sterimol definitions: morfeus-ml documentation
  (https://morfeus-ml.readthedocs.io), Sterimol: Verloop 1987.
- DELFIN validity limits: delfin/tools/adapters/morfeus_realism.py
  (machine-readable `validity_limits` on each adapter).
