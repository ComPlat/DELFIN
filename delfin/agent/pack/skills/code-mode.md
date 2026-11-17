---
name: code-mode
description: >
  Write one short Python program against DELFIN's tools (as a typed API) instead
  of one tool call per model turn. Best for fan-out work: checking many job
  folders, tabulating states and energies, summarising a run tree.
domains: code
---

# Code Mode

Instead of one model turn per tool call, write ONE short program and run it:

```
python -m delfin.agent.code_mode run my_program.py
```

## How to get the API

```
python -m delfin.agent.code_mode stub > tools_api.py
```

The stub is a typed module: one function per DELFIN read-only tool
(`job_status`, `read_calc_output`, `list_calcs`, `energy_table`, ...), each
returning parsed data. Import it in the program, call it, print the result.

## What the program may call

Read-only tools only (job status, output reading/parsing, listing, tables).
A call to a write/submit/delete tool raises an error naming the direct route —
submits and writes stay individually approved direct tool calls.

## Rules

- The program runs through the normal bash path, inside DELFIN's existing
  cage: no network beyond the existing proxy, no writes outside the work
  tree, wall-clock capped (120 s default).
- Every tool call goes through the same central gate a direct call takes —
  Code Mode is another way to reach the SAME checked tools, never a way
  around them. A gate refusal surfaces as an error in your program.
- Results are size-capped; a result over the cap comes back as an error
  naming the cap — filter or paginate in the program instead.
- Print a compact table or a short summary; only what you print reaches
  the model.

## Example

```python
import tools_api as api
failed, rows = [], []
for i in range(1, 21):
    info = api.read_calc_output(f"run_tree/job_{i:02d}")
    if info.get("outcome") != "succeeded":
        failed.append(info["job"])
    else:
        rows.append((info["job"], info["method"], info["energy_eh"]))
print("failed:", failed)
for job, method, energy in sorted(rows, key=lambda r: r[2]):
    print(f"{job}  {method}  {energy:.6f}")
```

Save it as `census.py`, then `python -m delfin.agent.code_mode run census.py`.

## When to use it

Many independent reads over many folders: 20 job folders' states + energies,
every conformer of an ensemble, all .out files in a run tree. For one or two
calls, call the tool directly.

## Status

Code Mode is OFF by default until the adoption measurement (the benchmark
family in `delfin/agent/pack/benchmark/tasks_w26_cm1.yaml`, preregistered
band) shows the gain; until then it ships disabled and this skill is the
reference for what it will do.
