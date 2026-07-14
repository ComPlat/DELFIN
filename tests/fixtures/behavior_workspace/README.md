# Behaviour-benchmark throwaway workspace

Self-contained toy files that the **behavioural parity** benchmark
(`delfin/agent/pack/benchmark/tasks_auto_behavior.yaml`) operates on to
measure whether the agent plans / scouts / verifies / asks the way an expert
engineer does.

These files are **meant to be edited and run by the agent** during a
benchmark. They are committed so the live runner can reset them to a clean
state **between replicates**:

```sh
git checkout -- tests/fixtures/behavior_workspace/
git clean -fd  tests/fixtures/behavior_workspace/
```

Nothing here is imported by the real DELFIN package — it is fixture data only.
The `thermo.py` and `pipeline.py` toys carry a deliberate bug / gap so a
scout/plan task has something concrete to read-before-editing.
