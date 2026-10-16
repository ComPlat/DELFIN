---
domains: code
---
# Energy profile (PES / HPES)
> Build a reaction energy diagram (PES, HPES, free-energy profile) in the group's fixed layout with make_plot kind energy_profile -- from given values or from outputs of any program it finds itself.

The user wants a reaction energy diagram: states as flat bars, dashed
connectors, one colour per pathway. The layout is fixed and drawn by the
`make_plot` tool -- do NOT write a matplotlib script for it, and do not
restyle it. Your job is the chemistry and the numbers.

## 1. Where the energies come from
- **The user gives them** (in the chat, or their `Index;deltaG` CSV):
  draw them as given -- `y` values in kcal/mol relative to the first state
  plus `ylabel` (`ΔG`, `ΔE`, `ΔH` -- what they ARE) and `"from_user": true`,
  or `"data": {"file": "<name>.csv"}`. Values in Hartree are refused, and so
  are typed-in values without `from_user` -- numbers YOU read from outputs
  always go through `files`.
- **You collect them from calculations**: find the outputs (`search_calcs`,
  `list_files`, `read_file`) and pass the PATHS as `files` -- the tool reads
  the energies and converts them itself. Do not read the numbers out and do
  the Hartree arithmetic yourself; one wrong digit gives a plausible wrong
  profile.
- Any program: ORCA through DELFIN's own parser, Gaussian, xtb, Psi4,
  GAMESS, Molpro, TURBOMOLE, Q-Chem, NWChem ... through cclib. The caption
  names the program per state. An output the tool cannot read is refused
  by name -- then ask the user for the values, do not guess them.
- Every state from the SAME program and method, same solvent and
  temperature. Outputs of two programs in one profile are refused; if the
  functional, basis or solvent model differ within one program, say so
  before drawing.
- `"energy": "gibbs"` (default, needs a frequency run) or `"electronic"`
  (ΔE). Some outputs carry no free energy (Molpro, NWChem): then ΔE, and
  say so.

## 2. Mass balance
A profile follows ONE system through the reaction: every state holds the
same atoms and the same total charge.
- CO2 binding to a metal starts at **cat + free CO2**, then cat-CO2. A start
  at the bare catalyst is not on the surface -- that HPES is unrealistic.
- A molecule that leaves (H2O, CO, a ligand) stays in every later state as
  a separate fragment.
- **From calculations (`files`)** the tool counts the atoms and the charge
  of every state from the outputs (last geometry, total charge) and REFUSES
  an unbalanced profile, naming what each state has too much or too little.
  Then explain it to the user in chemical terms (e.g. "the start lacks the
  free CO2: it has to be catalyst + CO2, computed separately") and build the
  state from fragments, a list whose energies add up:
  `[["cat.out", "co2.out"], "cat-co2.out"]`. If a fragment was never
  computed, say which calculation is missing -- do not drop the state's
  partner or the check.
- **Values the user gives** cannot be checked by the tool. Draw them; if
  the state names already show a gap (start "cat", later "cat-CO2"), point
  it out.
- When the number of separate molecules changes (association,
  dissociation) in solution, the 1 atm -> 1 M standard-state correction
  (RT ln 24.46 = 1.89 kcal/mol per molecule at 298.15 K) matters. Point it
  out and ask; the tool does not apply it (with `files` its caption says
  when the count changes).

## 3. Transition states
- **Computed TS**: a state like any other, named `TS1`, `TS2-3` ... and
  joined straight. It must have exactly ONE imaginary mode; with `files` the
  tool checks the last frequency block and warns in its caption -- pass the
  warning on, never hide it. A minimum with an imaginary mode is flagged too.
- **TS not computed** but known to be there ("unknown TS between INT1 and
  INT2"): `null` in `y` or in `files`. It is drawn as a dashed hump between
  its neighbours, without a value. Never make up a number for it.

## 4. Call make_plot
One object per state, in order. A state of several molecules lists all
their outputs in its `files` -- that is ONE state, the energies add up:

```json
{"kind": "energy_profile", "title": "CO2 binding at Co-salen",
 "states": [
   {"name": "cat + CO2", "files": ["cat/cat.out", "co2/co2.out"]},
   {"name": "TS1",       "files": "ts1/ts1.out"},
   {"name": "cat-CO2",   "files": "catco2/catco2.out"},
   {"name": "TS?",       "files": null},
   {"name": "P + CO",    "files": ["p/p.out", "co/co.out"]}]}
```

Values the user gave: `{"name": "TS1", "value": 18.2}` per state, plus
`"ylabel": "ΔG"` and `"from_user": true`. Several pathways in one figure:
`"series": [{"label": "S=0", "states": [...]}, {"label": "S=1", "states":
[...]}]` -- never one series per state.

- Energies from `files` are relative to the first output of the first
  pathway (or `"reference": ...`), so branches share one zero.
- Several pathways in ONE figure: several series. One pathway per figure:
  one call each, each with its own `filename`.
- A state named in two pathways is one column: a branch lines up by its
  names, no counting; a shorter route ends under the shared product.
- Any number of steps, 2 or 80: the figure widens and the names wrap or
  turn by themselves.
- `title`: the user's heading; otherwise the reaction in a few words.
- `label` only with names the user gave; without them leave it out.
  Never invent "ground state"/"excited" labels.
- Axis in quantity calculus, `ΔG / kcal mol⁻¹`, built for you: `ylabel` is
  the quantity only (`ΔG`, `ΔE`, `ΔH`, `ΔG‡`), `units` (`kcal/mol`,
  `kJ/mol`, `eV`) only what the numbers really are.
- `format: "svg"` for a vector figure, `"png"` (300 dpi) for slides.

## 5. Answer
The file path; per pathway the highest point (state and value); every
warning from the caption (imaginary modes, unchecked states); every state
left out and why; and where the numbers came from (the outputs, or "as
the user gave them").
