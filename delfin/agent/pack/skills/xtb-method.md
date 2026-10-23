---
name: xtb-method
description: >
  Choose the right xTB level (GFN0/1/2, GFN-FF), pass charge and unpaired
  electrons correctly, pick ALPB vs GBSA solvation, and know when semiempirical
  tight-binding is a pre-screen only and DFT must follow. Read this before
  running any of DELFIN's xtb_* adapters.
---

# Choosing an xTB method in DELFIN

DELFIN drives xtb two ways — never invent a raw `xtb` command line:

- **Native adapters** (`delfin/tools/adapters/xtb_native.py`, standalone `xtb`
  binary): `xtb_sp`, `xtb_optimize`, `xtb_hessian`, `xtb_md`, `xtb_fukui`.
  Parameters: `method` (gfn0 | gfn1 | gfn2 | gfnff), `charge`, `mult`
  (multiplicity → `--uhf mult-1`), `solvent` + `solvent_model` (alpb | gbsa),
  `accuracy`, `etemp`, `opt_level`.
- **ORCA-driven** (`delfin/tools/adapters/xtb.py`): `xtb_opt`/`xtb_goat` run
  ORCA's built-in `XTB1`/`XTB2` keywords; `xtb_solvator.py` builds an explicit
  solvent shell via SOLVATOR. There the method enum is strictly XTB1|XTB2 —
  GFN0/GFN-FF exist only in the native adapters.

## Which GFN level when

| Level | Flag | Self-consistent | Fits | Does NOT fit |
|---|---|---|---|---|
| GFN2-xTB | `--gfn 2` | yes | structures, vibrational frequencies, non-covalent (NCI) energies, HOMO/LUMO gaps, opt+Hess | high-accuracy energies, multi-reference cases |
| GFN1-xTB | `--gfn 1` | yes | older pipeline compatibility, very fast pre-screens | electrostatics of polar systems (no multipoles in SCF) |
| GFN0-xTB | `--gfn 0` | **no** (non-SCF) | huge systems (thousands of atoms), rough geometries/frequencies in seconds | anything needing accurate energies or gaps |
| GFN-FF | `--gfnff` | force field (no electrons) | very large systems, MD, relaxation before xTB | electronic properties, charges beyond EEQ model |

Ground rules:

1. **Default is GFN2** (`gfn2`) for anything ≤ a few hundred atoms where a
   structure, spectrum or relative energy is the answer. GFN2 with ALPB
   solvation is DELFIN's chemistry-task reference protocol
   (`tasks_chem.yaml:11-17`, xtb 6.7.1, opt+hess).
2. **GFN1/GFN0/GFN-FF** are pre-screen levels: geometry or sampling whose
   output is re-refined at GFN2/DFT. Their absolute energies, gaps and
   frequencies are not comparable to GFN2's.
3. **Levels are not comparable to each other.** Never mix GFN levels across a
   series, an energy table or before/after a retry — same level, same
   solvation, same accuracy, or the comparison is void.

## Charge and unpaired electrons

- `charge` maps to `--chrg`; `mult` (2S+1) maps to `--uhf mult-1` — so singlet
  = `--uhf 0`, doublet = `--uhf 1`, triplet = `--uhf 2`.
- Count electrons first: Z − charge; odd electron count ⇒ mult ≥ 2, and for
  transition-metal/open-shell systems state the multiplicity explicitly.
- Getting `--uhf` wrong silently gives a wrong (or impossible) spin state —
  the benchmark has a trap task for exactly this
  (`tasks_w20_xt1.yaml`, wrong-`--uhf` trap).
- The adapters validate both: a non-integer/negative/unphysical `--uhf` is
  refused, and the validated values cannot be shadowed via `extra_args`
  (neither bare nor `--flag=value`).

## Implicit solvation: ALPB vs GBSA

- `solvent_model=alpb` (default) — analytic linearized Poisson–Boltzmann with
  fast Born radii; generally the better-behaved of the two for typical
  solutes, and DELFIN's reference protocol uses it (`tasks_chem.yaml:12`).
- `solvent_model=gbsa` — generalized Born surface area; kept for
  compatibility with older xTB workflows.
- CPCM/COSMO are **ORCA-only** and are refused for native xtb —
  `build_xtb_cmd` raises on anything but alpb/gbsa.
- The `solvent` value is a name from xtb's list (water, acetone, ch2cl2,
  …). It is passed through verbatim — xtb validates it at runtime — but a
  flag-shaped value is refused by the adapter.

## Run types and what to report

- `xtb_sp` — single point; report energy (Eh) and HOMO–LUMO gap (eV) with
  method, level and solvation stated.
- `xtb_optimize` (`--opt`, levels crude…extreme, default normal) — follow with
  `xtb_hessian` (`--hess`): a structure is a minimum only with **zero**
  imaginary frequencies (report the count).
- `xtb_hessian` — thermochemistry in the log; quote it with the electronic
  temperature used.
- `xtb_md` (`--omd`) — sampling; the trajectory is a pre-screen, not
  equilibria.
- `xtb_fukui` (`--vfukui`) — Fukui indices; interpret as electrophilic/
  nucleophilic sites, not as charges.

## Known failure modes — when xTB is not enough

- **Carbocations can rearrange.** GFN2 optimises the isopropyl cation to
  protonated propene — a real, recorded artefact (`tasks_chem.yaml:19-21`).
  Grade/grade-off carbocation tasks on species where the ordering is
  unambiguous; check any "rearrangement" you see against chemistry before
  trusting it.
- **Multi-reference and transition-metal spectra**: gaps and relative
  energies near degeneracies are unreliable; use DFT (ORCA) for the numbers,
  xTB only for structures.
- **High-accuracy energetics** (barriers < ~1 kcal/mol, binding energies for
  publication): xTB is a pre-screen; DFT must follow. State "pre-screen only"
  in the answer.
- **Non-self-consistent levels** (GFN0) and **force-field levels** (GFN-FF)
  give no meaningful electronic gaps at all — never report a gap from them.

## References

- GFN1-xTB: Grimme, Hansen, Brandenburg, Bannwarth, *Chem. Rev.* **2016**,
  116, 5105–5141, DOI 10.1021/acs.chemrev.5b00547
- GFN2-xTB: Bannwarth, Ehlert, Grimme, *J. Chem. Theory Comput.* **2019**, 15,
  1652–1671, DOI 10.1021/acs.jctc.8b01176
- GFN0-xTB: Brandenburg et al., preprint DOI 10.26434/chemrxiv.8326202
  (published-version DOI to be pinned by QS)
- GFN-FF: Spicher, Grimme, *Angew. Chem. Int. Ed.* **2020**, 59, 10465–10473,
  DOI 10.1002/anie.202004239
- sTDA-xTB (screening spectra): Grimme, Hansen, *Angew. Chem. Int. Ed.*
  **2017**, 56, 10441, DOI 10.1002/anie.201704075

Every quantitative statement in this skill is a method-characterisation claim
backed by the reference next to it, not a measured DELFIN number; the measured
benchmark bands live in `tasks_chem.yaml` and `tasks_w20_xt1.yaml`.

