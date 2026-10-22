---
domains: code
---
# Conformer search (GOAT / CREST)
> Run a requested conformer or global-minimum search (ORCA GOAT, CREST) with the right charge, spin, solvent and topology, and report the ensemble with its degeneracies.

The user asked for a conformer search, a global minimum, or an ensemble.
Do what was asked, with the tool that was named; if none was named, ORCA's
GOAT where ORCA is installed, CREST otherwise. Refinement (DFT
re-optimisation, frequencies, higher-level single points) is a separate
step: offer it in the answer, run it only when it was asked for.

## 1. Charge and spin, the same at every step
- Determine the total charge and the number of unpaired electrons from
  the structure before anything runs; say both in the answer. A radical,
  an ion, a metal centre: never assume neutral singlet.
- Each program spells it differently -- check it in each input:
  - ORCA: `* xyzfile <charge> <multiplicity>`, multiplicity = 2S+1;
  - xtb and CREST: `--chrg <charge> --uhf <n>`, n = unpaired electrons
    (2S), NOT the multiplicity;
  - xtb also reads `.CHRG` / `.UHF` files in the working directory, which
    silently override a missing flag: look for them first.
- Open shell: the search runs on one spin state. GFN2-xTB does not rank
  spin states reliably (transition metals above all): search on the
  expected ground state, then compare spin states only at DFT on the
  low-lying conformers, and say so.

## 2. Environment and topology
- An ion, a zwitterion or anything polar is searched in implicit solvent
  (ALPB/CPCM matching the experiment): in the gas phase a zwitterion
  collapses or transfers a proton, and the ranking of conformers changes.
- A conformer search keeps the bonding. A proton that moved, a broken or
  formed bond, or a changed metal coordination is a different species
  (tautomer, isomer), not a conformer: compare the connectivity of every
  kept structure with the input and set those apart. GOAT and CREST check
  topology by default; for metal complexes check it yourself.
- Searching for new bonding on purpose is another request (GOAT-EXPLORE,
  GOAT-REACT): run it only when it was asked for.

## 3. Running it
- Keywords and block options from the manual (`search_docs`: "GOAT",
  "%goat") or `crest --help`; cite them, never recite from memory.
- GOAT: `! XTB GOAT` with `%pal` for the workers. Results, next to the
  input:
  - `<basename>.globalminimum.xyz`: the lowest structure found;
  - `<basename>.finalensemble.xyz`: every kept conformer, one frame
    each, its energy in the comment line;
  - in the output, the table after `# Final ensemble info #`: relative
    energy (kcal/mol), degeneracy, % population, then the lowest energy
    (Eh), Sconf and Gconf.
- CREST: `crest_conformers.xyz` (unique conformers, energies in the
  comment lines), `crest_best.xyz` (the lowest), `crest_rotamers.xyz`
  (with rotamers), and the conformer table in its output.
- Read the numbers from these files, never from the log of one worker.
  `show_molecule` on the ensemble file plays the conformers frame by
  frame; on the global minimum it shows that one in 3D.
- Cores: threads set to the cores asked for (`%pal`, `-T`,
  `OMP_NUM_THREADS`); a search on a login node stays small, anything long
  goes to the queue.

## 4. Reading the ensemble
- Energy window and duplicates: the filters (energy, RMSD, rotational
  constants) decide what counts as one conformer; state them.
- Degeneracy: GOAT gives every conformer degeneracy 1 (unless CONFDEGEN
  is set) and no factor for mirror images; CREST counts rotamers into the
  degeneracy. Mirror images (gauche+/gauche- of butane) have the same
  energy: say whether the ensemble lists both, and count the degeneracy g
  explicitly when populations are reported: p_i ~ g_i exp(-E_i/RT).
- Name the quantity: xtb energies are ΔE on the xtb surface. Populations
  from ΔE are not populations from ΔG; an experiment measures ΔH or ΔG.
- A global minimum from one search is the lowest found, not proven: say
  how complete the search was (cycles, conformers in the window).

## 5. With a rotation scan, when one is asked for too
The ensemble gives the minima; a relaxed scan along a dihedral gives the
path and the barriers between them. Run the scan as its own calculation
(constrained optimisations at each angle, same method as the search),
show its frames as one multi-frame .xyz, and check the two against each
other: every minimum of the scan must be a conformer of the ensemble at
the same energy, and a conformer the scan never reaches lies on another
coordinate -- say which.

## Output format
**Request and system** (charge, unpaired electrons, solvent) ->
**Input** (manual section cited) -> **Ensemble** (count, window, table of
relative energies, degeneracies) -> **What it means and what it does not**
-> **Possible next steps** (only as an offer).
