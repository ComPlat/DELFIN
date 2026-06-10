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
