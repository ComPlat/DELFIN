---
domains: code
---
# Energy profile (PES / HPES)
> Draw a reaction energy diagram (PES, HPES, free-energy profile) in the group's fixed layout with make_plot kind energy_profile.

The user wants a reaction energy diagram: states as flat bars, dashed
connectors, one colour per pathway. The layout is fixed and drawn by the
`make_plot` tool -- do NOT write a matplotlib script for it, and do not
restyle it. Your job is the numbers.

## 1. Get the energies right first
- Use the energy the user asks for (ΔG by default, ΔE / ΔH if they say so)
  and the SAME level of theory, solvent and temperature for every state.
- Relative to the first state of the route (that state is 0.00). From
  Hartree: (E_i - E_ref) * 627.509 kcal/mol.
- Every number must come from a calculation you read (`/analyze energy`,
  the output file, a table the user gave). A state with no result is left
  OUT and named in your answer -- never estimate or interpolate one.
- A transition state counts only with exactly one imaginary frequency.

## 2. Call make_plot
One series per pathway; `x` names every state, `y` gives its energy:

```json
{"kind": "energy_profile", "format": "png", "filename": "HPES",
 "title": "CO insertion at Co-S-salen",
 "series": [
   {"label": "singlet", "x": ["R", "TS1", "INT1", "TS2", "P"],
    "y": [0, 18.2, -3.1, 12.4, -15.0]},
   {"label": "triplet", "x": ["R", "TS1", "INT1", "TS?", "P"],
    "y": [0, 22.0, 5.5, null, -15.0]},
   {"label": "via TS3", "x": ["INT1", "TS3", "P'"],
    "y": [-3.1, 25.0, -8.0]}]}
```

- **Computed TS** (exactly one imaginary mode): a state like any other,
  with its name (`TS1`, `TS2-3`) and energy -- joined straight.
- **TS not computed** but known to be there (the user says "unknown TS
  between INT1 and P"): `null` as its energy. It is drawn as a dashed hump
  between its neighbours, without a value. Never make up a number for it.
- Any number of steps, 2 or 80: the figure widens and the names turn by
  themselves. Name every state the way the user names it.
- A state named in two pathways is one column: a branch lines up by its
  names, no counting. Begin the branch with the state it leaves from
  (`INT1` above), and a shorter route ends under the shared product.
- A state two pathways reach at the same energy is drawn once, split in
  their colours -- give both pathways the identical value.
- `title`: the plot heading. Use the user's; otherwise the reaction or
  system in a few words.
- `label` only with the names the user gave the pathways; without names
  leave it out (no legend). Never invent "ground state"/"excited" labels.
- Axis: quantity over unit, `ΔG / kcal mol⁻¹`, built for you. Set `ylabel`
  to the quantity only (`ΔG`, `ΔE`, `ΔH`, `ΔG‡`) and `units` (`kcal/mol`,
  `kJ/mol`, `eV`, `Eh`) ONLY to what the numbers really are -- relabelling
  is not converting.
- The user's own `Index;deltaG` CSV (semicolons, comma decimals, a `;` row
  between pathways, optional third column with the state name, `?` for a
  TS that was not computed) is drawn as is: `"data": {"file": "<name>.csv"}`.
- `format: "svg"` for a vector figure, `"png"` (300 dpi) for slides.

## 3. Answer
The file path, the highest point of each pathway (state and value), and
every state you left out and why. No restating of the whole table.
