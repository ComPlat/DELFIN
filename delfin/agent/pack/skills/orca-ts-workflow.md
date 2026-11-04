---
domains: code
---
# Transition states, reaction paths and scans (ORCA)
> Full TS workflow: guess structures, OptTS, verify the single imaginary mode, IRC to both minima, relaxed scans, barrier with ZPE at one level of theory.

The user wants a barrier height, a mechanism, or a scan — not just an
optimized geometry.

## 1. Level of theory is chosen ONCE and stated everywhere
- Pick ONE level (functional + basis + dispersion + solvent) for the
  whole barrier: reactant, TS, product AND ZPE. A barrier assembled
  from single points at different levels is not a barrier.
- Cheap-first where valid: pre-optimise reactant/product/TS-guess with
  `xtb_opt` (GFN2-xTB via ORCA, delfin/tools/adapters/xtb.py), then run
  the TS work at the final level. The final level is what goes in the
  answer, with the xTB step named as pre-optimisation only.

## 2. The workflow, in DELFIN steps
1. **Guess structures**: endpoints from `xtb_opt`/`orca_opt`; a TS
   guess from a linear-interpolated path between them (NEB-TS is
   ORCA's keyword for the automated variant), or a constrained
   pre-optimisation near the suspected saddle.
2. **OptTS**: `run_step("orca_ts_opt", geometry=<guess>.xyz, charge=…,
   method=…, basis=…, hessian_mode=…)`. The Hessian strategy must be
   named explicitly: `AnFreq` (analytic, cheaper when available) or
   `NumFreq` (numeric, needed when analytics are unavailable). Never
   let it default silently — the adapter refuses an unknown mode.
3. **Verify the TS**: `run_step("orca_ts_freq_check", …)`. SUCCESS
   means exactly one imaginary mode; the step reports `n_imag` and
   `imag_mode` (negative, as ORCA prints it).
   - 0 imaginary modes: that is a minimum, not a TS — go back to 2
     with a different guess.
   - 2 or more: a higher-order saddle — displace along the spurious
     mode and re-optimise; never delete the mode from the analysis.
   - The mode must MOVE ALONG THE REACTION COORDINATE — describe the
     motion from the normal mode (orca_pltvib) before calling it "the"
     TS. An HCN↔HNC TS must show the H transfer.
4. **IRC both directions**: `run_step("orca_ts_irc", geometry=<ts>.xyz,
   moread=<ts>.gbw, …)`. It refuses to run without the TS orbitals —
   that refusal is the safety rail, do not route around it.
   Optimize both endpoint geometries and run `orca_ts_freq_check`
   thinking mode on them: the IRC ends must be REAL minima, and they
   must be the stated reactant and product.
5. **Relaxed scans**: ORCA `%geom` Constraint/Scan blocks via the
   input's extra blocks (e.g. a methyl-rotation dihedral scan);
   constraints named in the input, not edited post hoc.

## 3. Barrier, honest about method and ZPE
- ΔE‡ = E(TS) − E(reactant), ΔG‡ = G(TS) − G(reactant) — all four
  numbers from the SAME level; the frequencies for ZPE/thermal
  corrections from the same level too (freq-thermochemistry skill).
- State unit, level, solvent, and whether it is ΔE‡ or ΔG‡. Compare
  with literature only within a band, cited; a single number without
  its method line is not a result.

## Safety rules
- Exactly one imaginary mode or the structure is not called a TS.
- Imaginary-mode counts, spin contamination, unconverged SCF are
  reported, never hidden.
- IRC endpoints are checked as minima before the mechanism is claimed.
- No barrier compared across different levels of theory.
- Keywords cited from the ORCA manual via search_docs, not from memory.
- Before the TS verdict, run the `result_checks` tool (delfin-ops) on
  the frequency run with `job_type="ts"`: it confirms exactly one
  imaginary mode with the mode index/frequency named, converged SCF,
  and that every number you report traces to the output file it was
  read from (a value with no source is flagged `not_computed`).

## Output format
**Reaction** → **Guess strategy** → **TS verification** (n_imag, mode
motion) → **IRC endpoints** (minima check) → **Barrier** (ΔE‡/ΔG‡,
one level, ZPE same level) → **Caveats**.
