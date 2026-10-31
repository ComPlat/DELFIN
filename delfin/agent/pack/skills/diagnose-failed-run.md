---
domains: code
---
# Diagnose failed run
> Systematically diagnose a failed/stalled calculation, propose a concrete fix — repair only with confirmation.

The user has a failed or stalled calculation (SLURM job, ORCA run, or a
DELFIN CONTROL.txt workflow) and wants to know WHY it failed and HOW to
fix it.  Follow this protocol strictly.

## 1. Collect evidence (read-only)

1. Read the tails of `*.err`, `*.out`, `*.log` and `slurm-*.out` in the
   calc folder (newest first).  Quote the exact failing line with file
   + line number — never paraphrase an error from memory.
2. If it is a DELFIN workflow, also read `CONTROL.txt` and the DELFIN
   log to see which step died.
3. `sacct -j <id> -o JobID,State,Elapsed,MaxRSS,ReqMem,Timelimit` for the
   SLURM view (OOM and TIMEOUT are often only visible here).

## 2. Classify the failure — two fundamentally different classes

**A. Mechanical / infrastructure** (deterministic, mechanically fixable):
- `bin/activate: No such file or directory` → venv/tarball unpack failed
- `oom-kill` / `OUT_OF_MEMORY` / MaxRSS ≈ ReqMem → memory
- `DUE TO TIME LIMIT` → walltime
- `Disk quota exceeded`, `Permission denied`, `No such file or directory`
- `NODE_FAIL`, segfault on one node → infrastructure flake

**B. Chemical / scientific** (requires human judgment):
- SCF not converged, geometry optimisation oscillates/diverges,
  imaginary frequencies, active-space problems, basis-set issues.

## 3. Ground the fix in the manual

For every ORCA/program keyword you propose, look it up first
(`search_docs` against the indexed manual from the Literature tab) and
cite the section.  Never propose a keyword from memory.

## 4b. DELFIN's own healing routes — use them before proposing an edit

DELFIN heals failures itself.  A failed run is repaired through DELFIN's
paths, never by rewriting the ORCA input by hand:

- **Recovery classes the engine heals itself**: SCF/DIIS/TRAH/MPI/memory
  and plain-geometry failures are handled by ORCA's recovery with
  modified parameters on retry — run the same job again through the
  recalc path (`/recalc` or the dashboard's Recalc) and DELFIN applies
  them.  Never edit `%scf`, `%geom` or keywords in the input yourself.
- **`stop_for_user` classes** (`RecoveryStrategy.reason_for_user()` in
  `delfin/orca_recovery.py`): SCAN_CONSTRAINT_ERROR (ORCA: "GSTEP: could
  not impose initial constraints") and SCAN_STEP_NOT_CONVERGED (a
  constrained scan step ran out of cycles, often with a SHARK/PROPERTIES
  abort).  No retry or parameter change helps — present the reason and
  the user's options (different start geometry/constraints; more or
  smaller scan steps; start the scan from the other side).
- **Input errors** (INPUT_ERROR, e.g. a keyword ORCA refuses): classify,
  quote ORCA's line, and give the user the corrected keyword — but the
  change goes through DELFIN's own input writer/CONTROL, never a
  hand-edited .inp.
- **Imaginary modes → IMAG**: a converged optimisation with an
  imaginary mode is IMAG's case, not a repair job.  Route it through
  IMAG (`delfin --imag` for the pipeline, or the `imag` workflow with
  the calculation's .inp/.out).  IMAG's default noise floor is
  `allow_imaginary_freq` = -50 cm-1 (the shipped default — imag.py
  NOISE_FLOOR, mirrored in control_validator.py; a measured value
  recorded in its comment).  A mode shallower than the floor is skipped
  BY DESIGN and silently.  If the user asks to heal a named mode (e.g.
  -21.86 cm-1), propose an explicit `allow_imaginary_freq` SHALLOWER
  than the measured mode (20 for -21.86, because the filter is
  `f < -abs(value)`), state the value and its reason, show it to the
  user and apply it through the `/control` path after confirmation —
  never silently, and never by changing the -50 default.

## 4. Propose — and only repair with explicit confirmation

- **Class A**: prepare the COMPLETE fix — the exact edit (diff of the
  submit script / CONTROL.txt / input) AND the resubmit command. Present
  both, then STOP and ask:
  *"Apply the fix and resubmit the job? (y/n)"*
  Only execute after explicit yes.
- **Class B**: present the analysis, 2-3 ranked options with trade-offs
  (e.g. damping vs. level shift vs. better guess), and what each implies
  scientifically. Do NOT apply chemistry changes autonomously — loosening
  convergence or shrinking an active space to "make it run" silently
  produces wrong science. The user decides.

## Safety rules
- Never `sbatch`/resubmit without the explicit confirmation above.
- Never weaken convergence criteria, never change method/basis/active
  space on your own.
- Max ONE repair attempt per confirmation — if it fails again, re-diagnose
  instead of looping.

## Output format
English prose: **Findings** (quote + file:line) → **Class** (A/B) →
**Root cause** → **Fix proposal** (concrete, with diff/command) → question.
