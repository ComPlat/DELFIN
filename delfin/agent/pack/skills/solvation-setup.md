# Solvation setup
> Choose and configure implicit solvation (CPCM/SMD) correctly — and say when implicit is not enough.

The user wants solvent effects in a calculation.

## 1. Choose the model for the QUESTION, not by default
- **CPCM**: geometries, relative energies in one solvent, spectra trends —
  cheap and robust.
- **SMD**: solvation free energies / partition equilibria — parameterised
  for ΔG_solv; say explicitly when SMD is the right tool.
- **Implicit is NOT enough** when the solvent coordinates (H-bonds to the
  active site, coordinating solvents like MeCN at open metal sites) —
  then recommend explicit solvent molecules + implicit bulk, and say why.

## 2. Ground the input in the manual
Look up the CPCM/SMD sections via `search_docs` before writing the block.
Verify the solvent NAME against the manual's solvent list — solvent names
are a classic silent failure (wrong/unknown name → calculation runs in
gas phase or aborts). Cite the section.

## 3. Build the input
- Solvent via the grounded name; SMD via the documented switch inside
  the CPCM block (verify exact syntax via search_docs — do not recite).
- Consistency rule: if geometries were optimised with solvation, single
  points for the same comparison should use the SAME solvation setup —
  flag mixed setups.

## 4. Sanity checks
- ΔG_solv magnitudes: neutral organics ≈ −1…−15 kcal/mol; ions tens of
  kcal/mol — values far outside → suspect setup.
- Charged species in low-ε solvents: warn about reliability.

## Safety rules
- Solvent names and switches only with a manual citation.
- No submits without explicit confirmation.

## Output format
**Modell-Wahl (CPCM/SMD/explizit, warum)** → **Input-Block** (Sektion
zitiert) → **Konsistenz-/Plausibilitäts-Checks**.
