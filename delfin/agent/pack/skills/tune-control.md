---
domains: code
---
# Tune CONTROL parameters
> Read the active CONTROL.txt, research what it should say for this system,
> propose it, and apply only what the user confirms -- through DELFIN's own
> `/control` path.

CONTROL.txt decides what DELFIN's classic workflows (OCCUPIER, ESD, TADF,
IMAG, redox, thermodynamics) compute. Every key must be understood, every
choice must have a source, and the user confirms before anything is applied.

1. Run `/control show` to read the current CONTROL content.
2. Understand the system. Identify the system class (organic,
   3d-metal complex, 4d/5d or 4f/5f complex, ESD/UV-Vis, NMR, redox,
   thermodynamics). Look at SMILES if visible, otherwise ask the user one
   short question about the target system.
3. Look each key up, in one place: the per-key notes in
   `delfin/common/control_key_notes.json` (meaning, valid values, default,
   which code reads it, and where a choice matters scientifically the
   best-practice value per system class with its source). The key LIST and
   the defaults come from the validator schema at runtime
   (`delfin.common.control_validator.CONTROL_FIELD_SPECS`) -- never from
   memory, never from a second list.
4. For choices the notes do not settle, use `search_docs` (ORCA manual,
   methodology documents) or DELFIN's web tools to verify a source; cite it
   in the proposal (manual section or DOI). Do not rely on memory.
5. Propose only the keys that should change. Each proposal carries:
   old value, new value, one-line reason, and its source.
6. Check each proposal through the validator BEFORE showing it to the user
   (validate_control_config / the coercers): a proposal the validator
   refuses never reaches the user. If the validator's cross-key
   consistency notes (`cross_key_consistency_notes`) fire on the current
   file, name them in the proposal -- they are findings, not silent fixes.
7. Stop and wait for the user's confirmation. Present the proposals as
   `/control key <key> <value>` commands. Apply a key only after the user
   confirms it, and only through DELFIN's `/control` path -- never by
   editing CONTROL.txt by hand.
8. After applying, run `/control show` once and read back every changed key,
   so the user sees what the file now says.

Edge cases:
- Relativity is DECIDED AUTOMATICALLY: delfin/cli.py forces relativity=none
  for light-element systems and the auto policy picks the relativistic
  branch for 4d/5d/4f/5f metals (delfin/utils._should_use_rel). Do not
  propose a relativity value unless the policy would not produce one and
  there is a specific reason; state that reason.
- If `freq_type` must become numFREQ, say that it is far more expensive and
  name the method that has no analytical Hessian.
- If `solvent` is set but `implicit_solvation_model` is not, the run goes
  gas phase silently; propose the model (or clearing the solvent) and say
  which.

Output format: a short bulleted list in the user's language, exact numerical
values, one line of reason and source per proposal.
