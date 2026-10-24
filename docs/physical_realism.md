# Physical realism rules (WB1)

Every pipeline/application spec that DELFIN runs through the workflow
engine is checked by the **scientific lint** — one entry point,
`delfin.tools._science.scientific_lint`, with rule modules under
`delfin/tools/_science_rules/`. A finding is a `warning` unless the
input is physically impossible (then `error`). Each finding names the
offending step and a concrete fix, so an agent or a person can correct
the spec before spending compute.

Severity convention: **error only for physically impossible input**;
everything that makes a number less trustworthy but not impossible is a
warning.

## Pre-existing rules (delfin/tools/_science.py)

| Rule | Level | Fires when |
|---|---|---|
| `composite-with-basis` | warning | a composite method (…-3c) plus an explicit basis it would ignore |
| `composite-with-dispersion` | warning | a composite method plus an explicit dispersion it already includes |
| `relativity-without-rel-basis` | warning | a relativistic Hamiltonian (ZORA/DKH/X2C) with a non-recontracted basis |
| `solvent-model-without-solvent` | warning | a solvation model with no solvent named |
| `invalid-multiplicity` | error | multiplicity < 1 |

## Electronic rules (delfin/tools/_science_rules/electronic.py)

**Rule 3 — `diffuse-basis-required` (warning).** An anion (`charge < 0`)
or a step marked as a Rydberg / dipole-bound / diffuse excited state
needs diffuse basis functions. Silent on `ma-def2-*`, `maug-cc-*`,
`aug-cc-*`, `pcH-*` bases (matched case-insensitively — `AUG-cc-pVDZ` is
legal) and on composite methods, whose basis the method defines itself.
Fix: use a minimally-augmented or augmented basis.
*Sources:* Rienstra-Kiracofe et al., *J. Phys. Chem. A* 2013 — diffuse
functions for anions, DOI 10.1021/jp312755z; Zheng et al., *J. Chem.
Theory Comput.* 2009 — minimally augmented basis sets, DOI
10.1021/ct900566x; Kunitsa & Bravaya, *Theor. Chem. Acc.* 2015 —
diffuse functions for dipole-bound/Rydberg states, DOI
10.1007/s00214-015-1647-1.

**Rule 4 — `relativistic-required` (warning).** A step whose `metals`
kwarg contains a 4d/5d element, a lanthanide/actinide, or a mixed set
needs a relativistic Hamiltonian and a relativistically recontracted
basis. The decision is DELFIN's own policy
`delfin.utils._should_use_rel` (`delfin/utils.py:203`), reused by
import, never copied; a pin test ties the rule's verdict to it.
Fix: set `relativity` (ZORA/X2C/DKH) and a ZORA-/DKH-/X2C-/SARC-
recontracted basis.

**Rule 10 — `electron-multiplicity-parity` (error).** Spin multiplicity
and electron count of the *same* parity (even mult with even electrons,
odd with odd) are physically impossible: ORCA stops such a run with
`multiplicity (…) is even and number of electrons (…) is even ->
impossible`. This is DELFIN's own rule (`delfin/config.py:709-764`,
146 archived run stops; parity-derived multiplicity
`delfin/utils.py:414-428`). Fires only when an explicit `electrons`
kwarg is present — the static lint counts no electrons it was not
given. *Note:* parity alone is a necessary, not sufficient condition;
a full multiplicity-vs-electrons check is a possible follow-up.

## Chain-consistency rules (delfin/tools/_science_rules/chain.py)

Spec-level: run over each walked step sequence (trunk and each branch
independently).

**Rule 1 — `level-consistency` (warning).** One energy difference uses
one level of theory: two level-setting steps in one chain with
concretely stated, differing method/basis/dispersion are flagged. ZPE
and thermal corrections belong to the level the geometry was optimised
at. Steps that omit a value inherit it at runtime — absence is not a
different level. *Sources:* DELFIN's "compare methods, not programs"
discipline (wave-20 rule 1); geometry/level consistency of ZPE,
*J. Chem. Theory Comput.* 2021, 17 — DOI 10.1021/acs.jctc.1c00474.

**Rule 2 — `freq-after-opt-missing` / `ts-freq-missing` /
`ts-irc-missing` (warnings).** Every geometry optimisation is followed
by a frequency calculation (`orca_freq` or `orca_opt_freq`): a minimum
has 0 imaginary modes. A TS optimisation needs a frequency check
(exactly 1 imaginary mode, `orca_ts_freq_check`) and an IRC to both
minima (`orca_ts_irc`). *Source:* DELFIN's own TS adapters
(`delfin/tools/adapters/orca_ts_adapters.py`; the IRC adapter fails
closed with `irc_endpoints=None`, lines 294–298).

**Rule 8 — `solvent-inconsistent` (warning).** Solvent and solvation
model must be consistent through the whole chain; CPCM and SMD free
energies are not comparable. Silent on inherited values.

**Rule 9 — `thermochemistry-conditions-missing` (warning).** A
frequency calculation without explicit `freq_temp` computes
thermochemistry at ORCA's default 298.15 K / 1 atm — the conditions of
the reported Gibbs energy must be stated. *Source:* ORCA 6.1.1 manual
§4.7 Thermochemistry (ideal-gas statistical mechanics; defaults
298.15 K, 1 atm); the adapter knob is `freq_temp`
(`delfin/tools/adapters/orca.py:196-197,347`).

## System / pre-DFT rules (delfin/tools/_science_rules/system.py)

These rules read an optional `system` hint kwarg
(`noncovalent` | `flexible` | `3d_complex`) and fire only when it — or,
for rule 7, the `metals` kwarg — is present. The static lint invents no
chemistry.

**Rule 5 — `dispersion-recommended` (warning).** A non-covalent or
flexible system without a dispersion correction (or a composite method
that includes one): London dispersion dominates non-covalent binding
energies and conformational energy differences. *Source:* Caldeweyher,
Ehlert, Hansen, Neugebauer, Spicher, Bannwarth, Grimme, *J. Chem.
Phys.* 2019, 150, 154122 — DFT-D4, DOI 10.1063/1.5090222 (verified
against the dftd4 reference list).

**Rule 6 — `conformer-search-missing` (warning).** A flexible molecule
has many low-energy conformers; without a conformer search before DFT
the DFT geometry may be the first one drawn, not the global minimum.
Conformer steps recognised: `crest_conformers`, `crest`,
`conformer_search`, `genarris`, `mlp_screen`. *Source:* Pracht, Bohle,
Grimme, *Phys. Chem. Chem. Phys.* 2020, 22, 7169–7192 — CREST,
DOI 10.1039/C9CP06869D.

**Rule 7 — `spin-state-scan-missing` (warning).** A 3d metal complex
(`system="3d_complex"` **and** a `metals` kwarg) with no spin-state
scan and no broken-symmetry statement: one multiplicity is not evidence
of the ground state; spin-state energetics are strongly
method-sensitive. Scan steps recognised: `occupier`, `spin_scan`,
`guppy` (DELFIN's own OCCUPIER module, `delfin/occupier.py`).
*Sources:* Wappett & Goerigk, *J. Chem. Theory Comput.* 2023, 19 —
MME55 benchmark, DOI 10.1021/acs.jctc.3c00558; DELFIN's
broken-symmetry sequence rule (`delfin/config.py:709-764`).

## Testing discipline

Every rule is tested in both directions (`tests/test_wb1_*.py`): it
fires on the wrong spec and stays silent on the right one, through the
real entry point `scientific_lint` (never an inner function). All rules
fire only on kwargs actually present; unresolved template placeholders
are skipped. The existing rule set and the platform's dict shape
(`{step, level, rule, message}`) are pinned by `tests/test_science.py`
and `tests/test_wb1_baseline.py`.
