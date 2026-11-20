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
python -m delfin.agent.code_mode stub > tools_api.py
python -m delfin.agent.code_mode run my_program.py
```

## How to get the API

```
python -m delfin.agent.code_mode stub > tools_api.py
```

The stub is a typed module: one function per DELFIN read-only tool, each
returning parsed data. Import it in the program, call it, print the result.

The functions it carries are `search_calcs`, `get_calc_info`, `list_docs`,
`read_section`, `bash_status` and `bash_output` — the read-only allow-list.
Read the names off the module rather than guessing at tool names.

## What the program may call

Read-only tools only (job status, output reading/parsing, listing, tables).
A call to a write/submit/delete tool raises an error naming the direct route —
submits and writes stay individually approved direct tool calls.

## Rules

- The program runs through the normal bash path, inside DELFIN's existing
  cage: no network beyond the existing proxy, no writes outside the work
  tree.
- Every tool call goes through the same central gate a direct call takes —
  Code Mode is another way to reach the SAME checked tools, never a way
  around them. A gate refusal surfaces as an error in your program.
- Results are size-capped; a result over the cap comes back as an error
  naming the cap — filter or paginate in the program instead.
- The wall-clock ceiling is the bash tool's, not Code Mode's: a program is
  one bash command, so `bash_timeout_s` (120 s by default, at most
  `bash_max_timeout_s`, 600 s) bounds it like every other command.
- Print a compact table or a short summary; only what you print reaches
  the model.

## Example

```python
import tools_api as api
failed, rows = [], []
for i in range(1, 21):
    info = api.get_calc_info(calc_id=f"job_{i:02d}")
    if info.get("outcome") != "succeeded":
        failed.append(info["calc_id"])
    else:
        rows.append((info["calc_id"], info.get("functional"),
                     info.get("final_energy_eh")))
print("failed:", failed)
for name, method, energy in sorted(rows, key=lambda r: r[2]):
    print(f"{name}  {method}  {energy:.6f}")
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
