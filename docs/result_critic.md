# `/check` — the scientific-correctness critic (`result_critic.py`)

`/check` scans a finished ORCA output for results that are numerically
"done" yet scientifically wrong. It is **report-only**: it never edits
an output, never re-runs a calculation and never loosens a criterion.

## What it checks

| Check | Code | Level | Source |
|---|---|---|---|
| Run finished cleanly | `termination` / `no-termination` | ok / error | ORCA output contract |
| SCF converged | `scf-not-converged` / `scf-not-fully-converged` | error / warn | DELFIN's own wording (api.py:516) |
| Geometry converged; opt has frequencies | `geom-converged` / `geom-not-converged` / `opt-no-freq` | ok / warn | minimum needs Hessian confirmation |
| Imaginary modes on a minimum / TS count | `no-imag` / `min-has-imag` / `ts-one-imag` / `ts-no-imag` / `ts-many-imag` | ok/warn/error | saddle point = exactly 1 imaginary mode |
| Spin contamination (ground state) | `spin-contamination` / `spin-clean` | warn ≥ 0.10, error ≥ 1.0 | \|⟨S²⟩ − S(S+1)\| |
| Energy sign | `energy-nonneg` / `energy-ok` | warn | electronic energies < 0 |
| **Mode count** | `freq-count-mismatch` | error | 3N−6 (non-linear) / 3N−5 (linear rotor, zero first rotational constant) — ORCA manual 6.1.1 §4.6/§4.7 |
| **IR intensities** | `ir-neg-intensity` / `ir-zero-intensity` | error / ok | absorption cross-sections are non-negative (Beer–Lambert) |
| **IR imaginary peaks** | `ir-imaginary-as-peak` | error | imaginary frequencies are not physical lines |
| **TDDFT roots vs request** | `uv-few-roots` | warn | `%tddft nroots` — requested window truncated |
| **Excited-state ⟨S²⟩** | `uv-spin-contamination` | warn ≥ 0.10, error ≥ 1.0 | per-state column, same scale as ground state |
| **Triplet in singlet spectrum** | `uv-triplet-in-singlet` | error | spin-flip mixture (only for a mult-1 reference) |
| **Broadening stated** | `uv-no-broadening` | ok | a spectrum compared with experiment needs its FWHM |

Bold entries are the MS1 additions (wave 22).

## Wiring: what gets parsed, and by whom

The critic reuses DELFIN's own parsers — it never re-parses ORCA
output:

- `delfin.api.extract_vibrational_modes` (api.py:2063) for the mode
  list and IR intensities;
- `delfin.api.extract_excited_states` (api.py:1517) for the TDDFT
  transition table;
- `delfin.imag.collect_imaginary_modes`, `delfin.parser.
  extract_last_uhf_deviation`, `delfin.energies.find_electronic_energy`
  as before.

These parsers take a **folder** and read its **largest** `.out`,
while `/check` scans **one** file. To keep "file flagged" and "file
parsed" from diverging silently, the spectra checks run only when the
critiqued file **is** the largest `.out` of its folder; otherwise the
report carries a visible `spectra-skipped` warning naming the file
that would have been parsed instead.

## Known limits (stated, not hidden)

- **Raman**: DELFIN has no Raman parser; `/check` reports Raman as
  unparsed rather than inventing a check.
- **Excited-state ⟨S²⟩ without a header**: the per-state ⟨S²⟩ is read
  only when the absorption table's header names the column (S^2,
  S**2, <S**2>). An unnamed last numeric column is never attributed —
  fosc/D²/XYZ values resemble ⟨S²⟩ and guessing flagged clean states
  in testing. The ground-state ⟨S²⟩ check remains the guard there.
- **Linear molecules**: linearity comes from the thermochemistry
  block (zero first rotational constant). Without that block the
  output is assumed non-linear (3N−6).

## For the agent using `/check`

Treat ❌ as "do not report this result to the user without saying
why". Treat ⚠️ as "report with the caveat, verbatim from the check
message". The messages are written to be quoted — they carry the
number, the physical reason and the unit.
