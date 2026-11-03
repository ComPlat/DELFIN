---
domains: code
---
# CASSCF setup
> Set up a CASSCF/NEVPT2 calculation with a justified active space — every keyword grounded in the manual.

The user wants a CASSCF (optionally NEVPT2) calculation. This is the
workflow where hallucinated keywords historically did the most damage
(`Nactel`/`Nactorb` do NOT exist) — ground everything.

## 1. Clarify the chemistry first
1. System, charge, multiplicity, and WHY multireference is needed
   (near-degeneracy? TM d-shell? bond breaking? excited states?).
2. Propose the active space from the chemistry, not from habit:
   - TM complexes: metal d-orbitals + directly interacting ligand
     orbitals (e.g. Fe(II): start from CAS(6,5)).
   - π-systems: full π/π* space when feasible.
   State electrons/orbitals explicitly: CAS(n,m) = n electrons, m orbitals.

## 2. Ground the input in the manual
Look up section 6.15 ("The CASSCF and NEVPT2 Modules") via `search_docs`
before writing the block. The canonical %casscf keywords are
`nel`, `norb`, `mult`, `nroots` (+ `ptmethod` for NEVPT2) — verify
against the manual extract, cite the section in your answer.

## 3. Build the input
- Start from a converged RHF/UHF or DFT orbital set; for TM systems
  consider `!UNO` natural orbitals as the starting guess.
- State-averaging (`nroots` > 1 + weights) when excited states matter.
- NEVPT2 on top: `ptmethod sc_nevpt2` (verify spelling in the manual).
- Sensible resources: %maxcore per core × cores ≤ SLURM memory request
  (leave the configured headroom).

## 4. Sanity checks before handing over
- Active-space occupations after convergence should be clearly partial
  (≈0.02–1.98); occupations pinned at 0/2 mean wasted orbitals — say so.
- Warn when CAS size explodes (>14 orbitals → CASSCF infeasible; suggest
  RAS/ICE or a smaller, chemically argued space).

## Safety rules
- NEVER invent keyword names — `nel`/`norb`/`mult`/`nroots` are real;
  anything else must come from a `search_docs` hit you can cite.
- Do not submit; present the input and let the user submit (or ask for
  explicit confirmation if they want you to).

## 5. Use the tooling, not your judgement alone (package CAS1)

DELFIN now ships a tool layer that enforces the checks of sections 1-4:
- `delfin.casscf_active_space` — `validate_casscf_block` (real keys
  only, `nel`/`norb` mandatory, `0 < nel <= 2*norb`, `norb <= 16`),
  `check_active_occupations` (pin check, 0.02-1.98 inclusive bounds,
  resize recommendation), `ActiveSpaceProposal` (chemical justification
  and named orbital classes REQUIRED; `nel <= 2*norb` is the physics
  limit — `nel > norb` is a valid open-shell space, e.g. CAS(6,5)),
  `CasscfConfirmationGate` (`confirm(user_confirmed=True)` /
  `revoke()`; a fresh gate per run).
- Workflow `casscf` (registry): `propose()` -> `show()` (the proposal
  is SHOWN TO THE USER) -> `confirmation_gate().confirm(...)` before
  the CASSCF step starts. Oversized spaces (>16 orbitals) are refused
  with RAS / ICE / DMRG named as the alternatives.

Run the proposal and the block through these before writing any
ORCA input; the workflow's `run()` itself never starts a calculation
without `CASSCF_user_confirmed` — the calculation goes through
DELFIN's normal ORCA submit path after the user confirms.

## Output format
**Chemical justification** → **Active space (CAS(n,m), why)** →
**Input block** (manual section cited) → **Expected checks**.
