# PACKAGE p3-a1 — review as a flow

The measurement for the review package: does a review reached through the
product produce findings that carry a `file:line`, a concrete failure scenario
and a verification verdict, and does it refuse to report a finding it could not
verify?

## Run it

```sh
PYTHONPATH=. python experiments/parity/p3-a1/measure.py \
    --out experiments/parity/p3-a1/before.json
```

`PYTHONPATH=.` is not optional. Run as a plain script, `sys.path[0]` is this
directory, so `import delfin` resolves to whatever copy is installed rather than
to the tree under test — and the run then reports a number about code that is
not this branch. The script refuses to measure when the loaded `delfin` is not
the tree it lives in, and `measured_tree` in the JSON records which one it was.

## What is driven

`delfin/agent/cli.py:main` — the front door `delfin-agent chat -p … --output-format
stream-json` uses — with the model's HTTP stream scripted. Everything between
the entry point and the finding is real: the engine's turn, the tool dispatch,
the `subagent` executor, `run_subagent`, the preset registry, the schema
validation. Only the model's words are scripted.

## The endpoints

| endpoint | meaning | how it is proven |
|---|---|---|
| `REACHED` | the findings-contract code executed | a `sys.settrace` call hook that fires only for frames whose file is the contract module |
| `USED` | the product path carried the contract's answer to the parent | the `subagent` tool payload the parent received holds the normalised review |
| `CORRECT` | the outcome is right | every verifiable finding survives with its line and scenario; every unverifiable one is dropped **and counted**; survivors are severity-ranked |

## The 12 variants

A small module is copied into the scratch workspace so each finding cites a
location that really exists. Each variant scripts a review that mixes findings
which can be checked against the workspace with findings that cannot: a missing
file, a line past EOF, a line of 0, no scenario, `verified: false`, a path
outside the workspace, an unknown severity, and an empty review.

## What this does not measure

Recall and precision of a review against real defects. Both are properties of
the model reading the code; scripting the model would make the number a
measurement of this script. They are named as probes in `scratch/realrun.json`
for whoever can run against a model.
