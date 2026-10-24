# Descriptors in DELFIN (wave-21 DS1)

Physical-realism status of every descriptor step, and how to run them.

## morfeus steric descriptors

Steps: `morfeus_buried_volume`, `morfeus_cone_angle`,
`morfeus_sterimol`, `morfeus_full_report` (registry names). Validity
limits are machine-readable on each adapter (`validity_limits`) and
**enforced at execute time** by the pre-execute physics gate in
`delfin/tools/adapters/morfeus_realism.py` — not just documented.

| Limit | Value | Why |
|---|---|---|
| `metal_index` / `atom*_index` / `donor_indices` | int, >= 0 | 0-based; wrapper converts to 1-based for morfeus. A negative index would silently become 0 after the +1 conversion and die in morfeus with a misleading error. |
| `radius` | 1–10 Å | Physically sensible ligand-shell range for the Poisson-Boltzmann sphere (%Vbur convention). |
| bite angle | exactly 2 `donor_indices` | Definition of the bite angle. |
| Sterimol axis | `atom1_index` != `atom2_index` | Two identical indices define a degenerate axis. |
| system | metal complex with a metal centre | Buried volume / cone angle are meaningless without a ligand sphere around a metal. |

If an out-of-range index still reaches morfeus (e.g. past the last
atom), the step fails with a ValueError naming the physics, not a bare
IndexError.

## Workflow: `steric_screening`

Contributed workflow (wave-21) chaining the three analysis steps on one
geometry:

    from delfin.workflows.registry import get
    get("steric_screening").run(
        config={}, geometry="complex.xyz", metal_index=0,
        donor_indices=[1, 2],  # optional, enables the bite angle
    )

Requires a converged, optimised metal-complex geometry; refuses a
negative index and a missing geometry file.

## Multiwfn bond orders and extra analyses

`multiwfn_bond_orders` (Wiberg/Mayer/etc. on a converged
single-determinant input, e.g. a .molden file), plus the wave-21
additions `multiwfn_population` (Hirshfeld/Mulliken/Löwdin/Becke) and
`multiwfn_orbital_composition`. Mulliken charges are flagged as
basis-partition dependent in the adapters' validity limits.

## cclib parsing

`cclib_parse`, `cclib_energies`, `cclib_vibrations` read ORCA output
files; they add no physics beyond what the output contains.
