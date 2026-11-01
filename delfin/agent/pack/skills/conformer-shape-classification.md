---
domains: code
---
# Conformer shape classification

> Classify and rank conformers of an existing ensemble by geometric shape
> (stretched, compact, planar, isotropic) and geometric chemical classes
> (E/Z-like, s-cis/s-trans, rotamers) — using the deterministic engine, never
> the model's own estimates.

The user asks WHICH conformers are most stretched / compact / planar /
linear / isotropic, asks for E/Z or s-cis/s-trans states, wants a conformer
described by its shape measures, or wants selected conformers exported.
Trigger phrases (German and English): "am stärksten gestreckt", "kompaktest(e)",
"am planarsten", "linearste", "isotrop", "Z-Konformere", "s-cis/s-trans",
"rotamer", "most stretched/compact/planar conformers", "exportiere",
"untersuche mein (final)ensemble", "klassifizierung(en) vorschlagen",
"langgestreckt", "speichere die … konformere", "zeige mir die … konformere".

## 0. Proactive profiling, NL selection, durable export

For "Untersuche mein Finalensemble und schlage mir sinnvolle
Klassifizierungen vor." use `delfin.analysis_tools.ensemble_profiling`
(same engine values, cached, status discipline included):

```python
from delfin.analysis_tools.ensemble_profiling import (
    profile_ensemble, suggest_classifications, resolve_selection,
    _select_and_rank, export_selection)
p = profile_ensemble(path_or_None)          # n, composition, energies,
                                            # descriptor spreads, fragments,
                                            # ez_cip (varies flag), diene
for s in suggest_classifications(p):        # 3-5 suggestions, each naming
    ...                                     # the measured spread it rests on
r = resolve_selection("die zehn am stärksten langgestreckten Konformere")
rows = _select_and_rank(path, r["descriptor"], r["selection"],
                        r["descending"], r["n"])     # [(id, value, relE)]
out = export_selection(path, [c for c, _, _ in rows], "top10_extended_dmax")
# -> xyz/csv/json in <calculations dir>/geometry_classification/;
#    overwrite needs overwrite=True; originals never modified.
```

Show an export with the `show_molecule` tool and state the
frame-to-conformer mapping (frame i = `conformer_id_to_export_frame`).
E/Z in every output is labelled a coordinate-based reconstruction
(RDKit bond orders + CIP), NOT independently validated stereochemistry.

### Context gating (dialog entry) — OBLIGATORY

Activation is CONTEXT-DEPENDENT. Reach for this skill only when the
user's message clearly concerns geometric conformer classification:
the trigger phrases above, or general interest in classifying/
comparing molecular shapes of conformers ("ich möchte meine Konformere
geometrisch klassifizieren", "welche geometrischen Unterschiede gibt
es in meinem Finalensemble", "was kann ich mit der geometry
classification machen", "ich möchte die molekülformen meiner
GOAT-Konformere vergleichen").

NEVER activate for: a plain greeting or chat start, unrelated DELFIN
tasks (quantum-chemical runs, job management, file operations without
a classification intent, general questions), or unprompted. Do NOT run
a full ensemble analysis automatically at chat start and do NOT inject
classification example prompts into unrelated answers.

Two-step entry when the interest is GENERAL (no concrete analysis
ordered): answer briefly inside the classification context with the
capability list (profile + suggestions; top-N rankings by dmax / npr1
/ plane_rms_A; E/Z-like, s-cis/s-trans, rotamer classes; export as
multi-XYZ with CSV+JSON metadata; 3D view), then ask whether the
ensemble should be profiled. If the user has ALREADY ordered a
concrete analysis ("die zehn am stärksten gestreckten …"), skip the
introduction and execute directly.

Follow-ups inside the classification context use the conversation
state: the last ensemble path, the numbered suggestions just given,
and the current selection. "Die zweite Möglichkeit" refers to the
second suggestion EXACTLY as numbered in your last answer. "Speichere
diese Auswahl" / "Zeige sie im 3D-Viewer" act on the most recent
selection. When the user changes topic, stop volunteering
classification suggestions; if they return with a new classification
question about the same ensemble, reuse the known path.

## 1. The engine does the science

All values come from `delfin.analysis_tools.geometry_classification_tools`
(over `conformer_shape_classification`). Never invent or estimate a number:
every measure you report must be a value that call returned. If a call fails,
say so instead of approximating.

```python
from delfin.analysis_tools import geometry_classification_tools as g

# 1) find the ensemble (explicit path or bounded, deterministic search)
ens = g.find_ensemble(path_or_None, search_dirs=None)   # ambiguity -> ValueError with candidates

# 2) descriptor catalogue / definitions
g.list_descriptors()                       # name -> definition, unit, range
g.describe_descriptor("dmax")              # descriptor definition and metadata

# 3) analyse once, then reuse the SHA256 cache automatically
g.analyze(path=ens)                        # cached; recompute only on change

# 4) ranking and filtering (values, units, direction are returned explicitly)
g.rank(path=ens, descriptor="dmax", reverse=True, top=5)
# heuristics (E-like/Z-like, s-cis/s-trans) -- geometric categories only:
g.filter_chemical(path=ens, classes=[("double_bond_geometry", "Z-like")])
# CHEMICALLY DEFINED E/Z (bond orders + CIP via RDKit; reconstruction from
# coordinates, NOT independently validated stereochemistry -- keep the
# caveat in every answer). Optional: restrict to one double bond:
g.filter_chemical(path=ens, classes=[("ez_cip", "E")], bond=(27, 0))
g.combined(path=ens, filters=[("double_bond_geometry", "Z-like")],
           ranking={"descriptor": "compacity", "top": 3})
g.describe_conformer(4, path=ens)          # one conformer, all its measures
g.export_multi_xyz([0, 4], out_path, path=ens)  # coordinates byte-identical
```

`find_ensemble` searches configured calc dirs and the current workspace only
(bounded depth). With several `finalensemble.xyz` files it raises a
ValueError listing the candidates — show the list and ask, do not guess.
The full preferred file is `GOAT_<name>.finalensemble.xyz`; the
`goat.N.N.finalensemble.xyz` siblings hold only subsets, avoid them when the
full one exists.

## 2. Which descriptor for which question

Shape terms are ambiguous. State the definition you used and, when the word
is genuinely ambiguous, offer the alternatives:

| German term | Primary descriptor | Alternatives to mention |
|---|---|---|
| am stärksten gestreckt / linearste | `npr1` ASCENDING (npr1 = I1/I3 in [0,1]; smallest = most rod-like; direction: `reverse=False`) | `dmax` (largest = longest), `npr2`, `eccentricity`; for "linear" also state npr2 (low npr1 + high npr2 = true rod) |
| am kompaktesten | `compacity` | `rg` (smallest), `kappa2` |
| am planarsten | `planarity` = (i2-i1)/(i3-i1), 1 = planar lamella | NOT `asphericity` (that measures rod-likeness) |
| isotrop / sphärisch | `kappa2` ASCENDING (relative shape anisotropy, 1.5*sum (tk-s/3)^2/s^2; 0 = isotropic, 1 = rod; `reverse=False`) | `npr1` AND `npr2` close to 1 (ideal isotropy) — npr1 ≈ 0.5 is NOT evidence of sphericity (a flat disc also has npr1 ≈ 0.5). CAUTION: `kappa2` and `asphericity` are MATHEMATICALLY IDENTICAL (both = 1.5*sum(t^2)/s^2 - 0.5) — never present them as independent cross-checks. |
| gestreckteste fragment | fragment-level: `shape_descriptors(..., selections={"fragment": [i, j, ...]})` | ask for an explicit atom selection |

Fragment questions: never invent atom indices. Since Phase 4 the
automatic recognition is available
(`from delfin.analysis_tools.conformer_fragment_recognition import
recognize_fragments, fragment_geometry, combined_ranking`):

- `recognize_fragments(path)` -> chemical fragments with status:
  "recognized" (bond perception succeeded and is identical on ALL
  frames), "plausible" (perception differs between frames), "unknown"
  (no consistent graph) / "ambiguous" (several backbones equally
  aromatic-adjacent). NEVER present a "plausible"/"unknown"/"ambiguous"
  fragment as a unique chemical definition. Bond perception derives
  the graph from COORDINATES via rdDetermineBonds(charge=0): the
  fragments are GRAPH-BASED RECONSTRUCTIONS, not independently
  validated chemistry — if the user has a SMILES/SDF reference, prefer
  it and say so.
- `fragment_geometry(path, atoms, kind="selection")` -> per-frame
  descriptors for EXACTLY those atoms, plus `plane_rms_A` (RMS
  distance to the best-fit plane, in A — the true out-of-plane
  measure, NOT the same property as the engine's PMI `planarity`).
- `combined_ranking(path, atoms, descriptor_a="plane_rms_A",
  ascending_a=True, descriptor_b="dmax", descending_b=True, top=5,
  tie="a_first")` -> transparent lexicographic two-criteria ranking
  over that atom selection (no silent whole-molecule fallback);
  state both criteria, both directions and the tie rule; no hidden
  weights. `atoms` is REQUIRED — a ranking without an explicit
  selection is refused.
- `pareto_front(path, atoms, minimize=("plane_rms_A",),
  maximize=("dmax",), top=None)` -> non-dominated SET, not a ranking;
  never present it as "the best".
- German terms: "das konjugierte Grundgerüst" -> conjugated_backbone
  fragment; "der aromatische Ring" -> aromatic_ring;
  "die konjugierte Dieneinheit" -> the conjugated_pi component. If the
  status is not "recognized", disclose it and ask.

Terminology discipline (keep these apart; the German terms are what the
user types, the engine identifier is what you use):
- "maximale Ausdehnung" = `dmax` (A, longest end-to-end distance);
- "Stabförmigkeit/Linearität" = `npr1` (I1/I3, small = rod) with `npr2`
  as complement; a long molecule is not automatically rod-like;
- "echte Fragmentplanarität" = `plane_rms_A` (A, best-plane RMS);
  the engine's `planarity` is a PMI ratio — different property,
  never call them equivalent.

## 3. Chemical classes — status discipline

E/Z and s-cis/s-trans are GEOMETRIC classifications from coordinates without
bond orders or CIP priorities. Report them as "Z-like"/"E-like" and
"s-cis/s-trans (geometric)", never as confirmed stereochemistry. Conformers
whose dihedral falls in the unclassified gap (45-135 deg for E/Z) have
status `not_assigned` — they are NOT matches for a Z-like or E-like filter
and count neither way; say how many were excluded as unassigned.

## 4. Presenting results

A complete answer includes: input file, number of conformers analysed, for
each hit: conformer id + original frame index, the descriptor used, its
value and unit, the ranking direction, relative energy in kcal/mol when
available (`rel_energy_kcal_mol`), and the caveats above. Keep the actual
numbers exactly as returned — no rounding beyond what the engine reports,
no invented precision.

## 5. Export and showing in the 3D viewer

`export_multi_xyz(indices, out_path, path=ens)` writes the selected frames
with unchanged coordinates into a multi-XYZ file. Name the output path in
the answer. The original ensemble file is never modified.

To SHOW a selection, call the `show_molecule` tool with the exported
file's path: the dashboard renders an interactive 3D card in the agent
chat (frame navigation included; limits: 4 MB and 500 frames per card,
a 38-atom 10-frame export is ~30 KB). Frame order in the card equals
the `indices` order — state the mapping "card frame i = conformer
<indices[i]>" in the answer so the original IDs stay traceable.

The Calculations-browser viewer (xyz_frames / xyz_current_frame) is a
SEPARATE viewer inside the dashboard state; there is no agent-side hook
that sets it. If the user asks for the browser viewer, point them to
the exported file to open there — do not claim the browser viewer was
driven remotely.
