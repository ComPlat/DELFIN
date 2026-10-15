# generalist_workspace

Workspace root for the `generalist` benchmark category
(`delfin/agent/pack/benchmark/tasks_generalist.yaml`).

The directory is empty on purpose. Each task's setup script
(`delfin/agent/pack/benchmark/setup/generalist_*.py`) builds one git
repository in a subdirectory (`textkit/`, `invoice-ts/`, `datewin/`, ...)
before the attempt. The acceptance script
(`delfin/agent/pack/benchmark/accept/generalist_*_verify.py`) then runs that
repository's tests and inspects its files and history. The benchmark runner
snapshots this directory before an attempt and restores it afterwards, so
nothing a run builds survives into the next one.

Nothing here is DELFIN or chemistry: the repositories are ordinary Python and
Node/TypeScript projects. Their purpose is to measure the agent the way a
general coding agent is judged.
