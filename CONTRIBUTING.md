# Contributing to DELFIN

DELFIN automates quantum-chemical calculations that people publish results
from. A change that is merely plausible is not good enough here: a bug can
produce a number that looks entirely reasonable and is wrong. Most of what
follows comes from that.

## Getting set up

```bash
git clone git@github.com:ComPlat/DELFIN.git
cd DELFIN
python3.11 -m venv .venv && . .venv/bin/activate
pip install -e ".[agent,docs,dev]"
pytest -q -m "not slow"
```

Python 3.11 is what production runs and what CI gates on. `pyproject.toml`
also permits 3.10, but nothing verifies it — do not rely on it.

## Before you open a pull request

- **Run the fast suite.** `pytest -q -m "not slow"` is the merge gate.
- **Prove it against the real thing.** If your change affects how ORCA, xTB or
  CREST are driven, run it. Unit tests do not catch that ORCA looks for a
  restart file in its working directory rather than next to the input, or that
  it words an out-of-memory abort differently than the pattern expected. Both
  were found by running, not by reading.
- **Add the test that would have caught it.** Name it after the behaviour it
  protects, not after the function it calls.
- **Say what you measured.** "Fixes the retry" is not reviewable. "Benzene
  Hessian: 161 s cold, 8 s resumed, 0 displacements recomputed" is.

## Things that are easy to get wrong here

- **Do not change what is computed while fixing how it is scheduled.** Cores,
  restarts, logging and recovery may be reorganised freely. Functionals, basis
  sets, solvation and convergence criteria change results and belong in their
  own, separately justified change.
- **CONTROL keys are a user interface.** Renaming or repurposing one breaks
  input files people have kept for years.
- **A killed job is not a failed job.** Anything that decides a job is finished
  has to distinguish "ORCA said it converged", "ORCA gave up" and "the walltime
  ran out". They look alike and mean entirely different things.
- **Outputs get large.** Hundreds of MB is normal. Read a bounded tail, never
  the whole file.

## Style

`ruff` and `black` are configured in `pyproject.toml`. The lint job in CI is
informational for now while the existing findings are worked down — please do
not add new ones.

Comments should explain why something is the way it is, especially where the
obvious approach fails. The repository has a lot of hard-won detail in its
comments; keep that habit.

## Reporting a security problem

Not through an issue — see [SECURITY.md](SECURITY.md).
