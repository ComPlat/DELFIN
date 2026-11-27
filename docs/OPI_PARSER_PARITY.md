# Which ORCA quantities DELFIN reads through OPI, and which stay on its own parsers

DELFIN reads ORCA output itself (`delfin.energies`, `delfin.parser`,
`delfin.imag`, `delfin.ir_spectrum`, `delfin.esd_results`, and the readers
behind `delfin.api`) and it can also read the same quantities through ORCA's
own Python interface, OPI (`pip install "delfin[orca]"`). `delfin.orca_read`
is the one place that decides, per quantity, which of the two answers.

This file records that decision and its evidence, so a later reader does not
have to re-measure it. Every row was measured on 2026-11-26; where a value
could not be compared the row says so instead of implying agreement.

## How the decision is made

`orca_read.read(work_dir, base_name, quantity)` returns a `Reading` whose
`source` is `"opi"`, `"text"` or `"none"`. OPI answers only when all of these
hold, and the text parser answers otherwise:

1. OPI is importable (it is an extra; CI's main job has no OPI).
2. The `.out` header names ORCA >= 6.1.1, OPI's own contract.
3. `<base>.property.json` is present — the file OPI actually reads.
4. The run terminated normally and the SCF converged, checked on the text.

When both sides answer, the text value is the record: it is returned even
when OPI agrees (`source` is then `"opi"` and the note says "confirmed by
text shadow check"), and on a mismatch the text value is returned with a
warning naming the quantity, both values and the base name.

OPI is asked for the property half only. Its `Output.parse()` defaults to
also building gbw JSON, which runs `orca_2json` whenever a `<base>.gbw` sits
beside the output — and every real calculation folder keeps one. See the
first two entries under "Measured defects" below.

## The table

Population: the 18 committed `*.out` files under `tests/fixtures/`, of which 3
carry a `property.json` and only one (`orca611/synth_h2o`) carries a
thermochemistry block that both readers can answer. That is the honest size of
the evidence: the value comparisons below rest on 1 fixture × 5 quantities,
not on a corpus. Consumers are counted as modules under `delfin/` that name
the reader.

| Quantity | DELFIN parser | OPI getter (source) | Consumers | Decision | Why |
|---|---|---|---|---|---|
| Final single-point energy | `energies.find_electronic_energy` | `get_final_energy` (property JSON) | 11 modules | **OPI first** | Identical on the one fixture where both answer: `-76.32126902150691` Eh, exactly. Text shadows it. |
| Gibbs free energy | `energies.find_gibbs_energy` | `get_free_energy` (property JSON) | 11 modules | **OPI first** | Identical: text `-76.3177312`, OPI `-76.3177312027805` — the text prints 8 decimals, the JSON carries the double. Same number. |
| ZPE | `energies.find_ZPE` | `get_zpe` (property JSON) | 6 modules | **OPI first** | Identical: `0.02120238` vs `0.021202376295716842`; same sign (a positive correction). |
| Total enthalpy H | *(tail scan)* | `get_enthalpy` (property JSON) | — | **OPI first** | Identical: `-76.29628713` vs `-76.29628713260958`. |
| Inner energy U | *(tail scan)* | `get_inner_energy` (property JSON) | — | **OPI first** | Identical: `-76.29723134` vs `-76.29723134165184`. Note H = U + kT holds to 5 decimals, so the two are genuinely different quantities and neither is the other under a second name. |
| `<S**2>` deviation | `parser.extract_last_uhf_deviation` | `get_s2` (**`.out` grep**) | `occupier.py`, `reporting/occupier_reports.py`, `tools/result_checks.py` | **OPI first, text shadows** | Both grep the same `UHF SPIN CONTAMINATION` block; OPI returns `(expected, ideal)` and `read` returns their difference, same sign convention. That relation was measured in phase 1 on the real `o2_triplet` run (UKS triplet): text `Deviation` = 0.007206, OPI (2.007206, 2.0) — see `tests/test_opi_equivalence_q04.py`. **Not re-measured here**: no committed fixture carries the block (grep: none of the 18), so on this population the S2 row rests on that earlier measurement, not on a fresh one. On a closed-shell run neither side has a block and `read` answers `0.0`. |
| Multiplicity | `.inp` header | `get_mult` (property JSON) | `tools/result_checks.py` | **Text** | Text reads the *input's* intended multiplicity, which survives a failed run; OPI reads what ORCA *did*. Different questions, so neither is a replacement. On `synth_h2o` both answer `1`; on a run that crashed before SCF only the text side has a value, and `result_checks` wants the intended one. |
| Total charge | `.inp` header | `get_charge` (property JSON) | (via `result_checks`) | **Text** | Same reason as multiplicity. |
| Number of electrons | derived from `.inp` geometry + charge | `get_nelectrons` (property JSON) | — | **Text** | Measured: OPI answered `(None, None)` on the fixtures that do carry a property JSON, and `(28, None)` on the DLPNO one. DELFIN's derivation answers on every fixture. OPI lacks the quantity as read here. |
| Terminated normally | `smart_recalc.has_ok_marker` | `terminated_normally()` (`.out` grep) | `tools/result_checks.py` | **Text** | Same grep, but DELFIN's reader is the one `orca_recovery` and `smart_recalc` already treat as the record, and it is one call with no OPI install. No gain from routing; OPI is not faster (see Timing). |
| Failure reason | *(outcome phrase of the folder)* | none | — | **Text** | OPI says *that* a run failed, never *why*. It has no quantity to compare. |
| Geometry-optimisation converged | text markers, `None` when the run did not optimise | `geometry_optimization_converged()` (`.out` grep, "HURRAY") | `tools/result_checks.py` | **Text** | OPI answers `False` for a run that never optimised (the "HURRAY" marker is simply absent), which would turn "not applicable" into "failed". That is a semantic difference, not a missing value, and `text` is the answer DELFIN's callers already branch on. |
| Final geometry | `<base>.xyz` written by ORCA | `get_structure` (reads `results_properties…geometry.coordinates`; needs the gbw JSON in practice) | `imag.py`, `orca_recovery.py` | **Text** | DELFIN's reader is the chaining hand-off and must return exactly what the sibling `.xyz` carries. Measured on `synth_h2o`: `get_structure()` answered `None` in **all four** combinations of a stand-in `.gbw` and `.json` present/absent, so OPI gave nothing to compare either way — and my `read_gbw_json=False` change therefore removed no geometry that was reachable (that matters: it was the one quantity the gbw half might have served). |
| IR modes | `ir_spectrum.parse_ir_spectrum` | `get_ir` (`.out` grep) | `ir_spectrum.py` | **Text** | No committed fixture carries an `IR SPECTRUM` block, so no comparison was possible. Routing it would be a change with zero evidence behind it. |
| Imaginary modes | `imag.collect_imaginary_modes` | none | `imag.py`, `api.py` | **Text** | OPI has no getter for it. |
| Vibrational frequencies (`.hess`) | `imag.hessian_frequencies` | none | `imag.py` | **Text** | OPI reads JSON, never `.hess`; the run folder's `.hess` is not in its model. |
| Dipole moment | `api.extract_dipole` | `get_dipole` (property JSON) | 5 modules | **Text for now** | OPI has the quantity and it is in the property JSON, but **no committed fixture carries a dipole block in either form** — OPI answered `None` and DELFIN's reader found nothing to read. A decision without a comparison is a guess; left on text until a fixture with both exists. |
| Polarizability / quadrupole | `parser.parse_polarizability` | `get_polarizability`, `get_quadrupole` (property JSON) | 1 module | **Text for now** | Same: OPI has the getters, no fixture exercises them. |
| Mulliken / Loewdin charges | `api.extract_mulliken_charges`, `api.extract_loewdin_charges` | `get_mulliken`, `get_loewdin` (property JSON) | 5 modules each | **Text for now** | Same: no population block in any committed fixture, on either side. |
| Mayer / MBIS / Hirshfeld / CHELPG | via Multiwfn (`analysis_tools/multiwfn_wrapper.py`), not from ORCA output | OPI has all four getters (property JSON) | 3 modules (Mayer, Hirshfeld) | **No change** | DELFIN's route for these goes through Multiwfn, not through this reader; `chelpg` appears only as a `%chelpg` block-name constant in `common/orca_input.py`, so there is no ORCA-reading parser to compare. Worth knowing OPI could serve them directly. |
| HOMO / LUMO / gap | `api.extract_orbital_energies` | `get_homo`, `get_lumo`, `get_hl_gap` (gbw JSON) | 6 modules | **Text** | OPI's route needs the gbw half, i.e. `orca_2json` on the run folder. That is a heavier and host-dependent path than the text block DELFIN reads today, for a quantity with no parity evidence on this population. |
| TDDFT / absorption spectrum | `api.extract_excited_states`, `energies.find_state1_ohne_SOC` | **no getter** | 7 modules | **Text** | OPI 2.0.0 exposes no absorption-spectrum getter, even though its property model has an `absorption_spectrum` field. The quantity DELFIN's ESD path needs is not reachable through OPI's public API. |
| NMR chemical shifts | `nmr_spectrum.parse_nmr_orca` (`CHEMICAL SHIELDING SUMMARY`) | **no getter** | `reporting/nmr_report.py`, `ensemble_nmr.py` | **Text** | OPI's property model has a `chem_shift` record but `Output` exposes no getter for it (checked: no name in `dir(Output)` mentions nmr/shift/spectrum). |
| Entropy term | tail scan, `Final entropy term` (T·S, Eh) | `get_entropy` (property JSON) | — | **Text** | *Naming* is unresolved: ORCA prints T·S in Eh, not a molar entropy, and OPI's docstring says "entropy (in a.u.)". Both carry the same number (`0.0214440701709188` vs `0.02144407`), but a value whose unit is argued about is not one to route on. `read_thermochemistry` already handles ORCA's other label (`Total entropy correction`, −T·S) by negating it. |
| Failure/unterminated runs | text | partial | — | **Text** | Encoding a failed run through OPI loses the reason. |

## Measured defects found while doing this

Both were fixed in the commit that added `tests/test_opi_gate_reads_the_files_opi_reads.py`,
and both had silently sent every read to the text parser:

1. **The gate required `<base>.property.txt`, a file OPI never opens.**
   Instrumenting `open` across a full `Output.parse()` counted **0** opens of
   that file. It appears in OPI only as the *input* of `orca_2json` in
   `opi.execution.core.Runner.create_property_json`, which DELFIN does not
   call. On `synth_h2o` the reader answered `source="text"` with the note
   "property txt missing" although the property JSON was correct.
2. **`parse()` ran `orca_2json` for a `<base>.gbw`.** `parse()` defaults to
   `read_gbw_json=True`, so with a `.gbw` in the folder (every real run has
   one) it called `create_missing_gbw_json`. Without `orca_2json` in OPI's
   ORCA bin folder that raised `FileNotFoundError`, and the read fell back to
   text. Measured on a copy of `synth_h2o` with a stand-in `.gbw`: default
   `parse()` raised, `parse(read_gbw_json=False)` returned
   `G = -76.3177312027805`. The reader now parses the property half only.

**Still open, deliberately:** DELFIN's SCF-convergence check looks for
`SCF CONVERGED`; OPI's looks for `SUCCESS`. On `synth_h2o` only the second is
present, so even after fix 1 DELFIN's gate refuses that run (`"SCF not
converged"`) while OPI answers. Which marker is the contract is a separate
decision — DELFIN's `api.parse_orca_output` uses the same `SCF CONVERGED`
string — so it was left alone and is recorded here rather than changed in a
parser-parity change.

## Timing

Measured on `synth_h2o` with a `.gbw` present and the gate passing, median of
7 runs, Python 3.12 / OPI 2.0.0:

| What | Time |
|---|---|
| One quantity through `read` (OPI-backed) | 0.58 ms |
| Six quantities through `read` (OPI-backed) | 3.16 ms |
| One `_opi_values` parse | 0.18 ms |
| One text tail scan (`_text_values`, one quantity) | 0.054 ms |

OPI costs about 10× the text scan per quantity, and `read` re-parses the run
for each quantity — six calls cost six parses. That is acceptable where the
answer feeds a report and not where it is called in a loop over a folder; it
was not a reason to move any quantity by itself, since the energies above are
already routed.

## Not measured here

- **Real ORCA outputs in an archive.** The package calls for a local-only
  parity run over them (counts into the PR, never the content). This sandbox
  denies the shared archive, so it was **not run**; the numbers above
  rest on committed fixtures only. The probe is documented in `scratch/PR.md`.
- **A real finished calculation read through both readers for a report.** The
  one fixture with a report-relevant thermochemistry block is `synth_h2o`,
  and it is synthetic (its `PROVENANCE.md` says so). "Real run" evidence
  requires a real output, which is in the archive above.
