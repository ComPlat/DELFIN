# Validator A (protocol author), 20 drawn tasks, 2026-11-20

Rubric: experiments/2026-11-20-specialist-vs-generalist.md Addendum 4,
committed at 07e7ab94f BEFORE any row was read. Judged against origin/main
at 5f1975075. Validator A has a prediction on record (Addendum 3) and is
therefore not disinterested; a second validator judges the same 20 blind.

| # | name | verdict | reason |
|---|---|---|---|
| 1 | `_max_heavy_bond_order` | usable | clear docstring; the only product hit is a comment saying it is "retained", no own tests, so deletion is checkable |
| 2 | `_orbit_burnside_selftest` | usable | a self-test of `_enumerate_orbits`; a self-test does not belong in the product, so both routes are determinable |
| 3 | `_orca_ok` | usable (weakest) | no docstring; a `-> bool` success check on an ORCA output path in an xtb/CREST module, intent legible from the signature alone |
| 4 | `describe_mlp_backends` | usable | docstring names its intended consumer ("for agent/workflow use"); see the instrument finding below |
| 5 | `log_orca_job_plan` | usable (weak) | no docstring; `(label, jobs) -> None` in the occupier, so the call site is the point where jobs are planned |
| 6 | `check_goat_e2e` | usable | clear docstring, 7 test files reference it, only product hit is a comment |
| 7 | `est_reference_comparison` | usable | clear docstring, 2 own test files |
| 8 | `load_fukui_result` | usable | clear docstring, 1 own test file |
| 9 | `validate_relspin` | usable | no docstring but 1 own test; an adapter validator, the same shape a tool-contract gap was just fixed in |
| 10 | `write_co2_scan_opi_input` | usable | no docstring, 2 own test files |
| 11 | `_project_config_path` | usable | clear docstring, no own tests |
| 12 | `_step_finished` | usable (weak) | no docstring; the durable step machinery is where it belongs |
| 13 | `controller_hint` | usable | docstring names the consumer ("for a surface's mode line") |
| 14 | `expected_digest` | **rubric gap** | the work is already done in flight: PR #834 deletes it and is in the landing queue |
| 15 | `is_harness_text` | usable | clear docstring |
| 16 | `bwrap_substitution` | **blocked** | cage machinery: its docstring has it replacing `--tmpfs /tmp` with `--bind cage_mount(...)` "so the cage build can assert the swap" |
| 17 | `distill_for_shipping` | usable | no docstring, 3 own test files |
| 18 | `pdf_xfa_data` | usable | clear docstring, 1 own test file |
| 19 | `status_validate` | usable | docstring names the contract it checks, 1 own test |
| 20 | `step_resume_note` | usable | clear docstring, 2 own test files |

**Counts: 18 usable, 1 blocked, 1 rubric gap.** Two to replace, under the
5-of-20 threshold at which the corpus would have been the finding.

## The rubric needed a category, and that is reported rather than hidden

Row 14 fits none of the four verdicts. It is not underspecified, not
mis-scoped, and not blocked: it is **already landed** — PR #834 deletes
`principles_guard.expected_digest` and sits in the landing queue, so by the
time an arm runs, the task is done and the arm would be graded on work it
did not do.

The needed category is **already-landed**: the capability is removed or
connected by work in flight. Adding a category after reading a row is
exactly what fixing the rubric first was meant to prevent, so the fact that
it had to be extended is recorded here instead of the row being filed under
the nearest existing label.

It also implies a mechanical pre-flight that does not exist: **check each
drawn row against every open PR**, not only against `main`. A sample drawn
from `main` goes stale at the speed of the landing queue, and 17 PRs landed
on the day this one was drawn.

## Two corrections made during the validation, both from checking a match

- `expected_digest` first looked like it had **4 product callers**. A
  word-boundary search finds **none**: the hits were `expected_digests`, the
  plural parameter of `check()`. The scan was right and the first grep was
  sloppy.
- `_max_heavy_bond_order` and `check_goat_e2e` also looked like they had a
  product caller each. Both hits are **comments**. A comment is not a caller,
  and a substring is not a reference; getting this wrong in either direction
  is the main way this validation could fail.

## A finding about the instrument, not about the tasks

The apparatus reported that **exactly one** of the 20 has no consuming
suite (`describe_mlp_backends`, `delfin/mlp_tools/__init__.py`). A
word-boundary search finds **four** test files importing from
`delfin.mlp_tools`, which load that `__init__.py`. So the consuming-suite
resolver **under-reports**, and check 2 would have been treated as vacuous
for a row where it is not. The resolver mirrors the product's
`_related_test_files` and caps at 12 files; whether the cap or the mapping
is the cause is not established here.

## Is this set of 20 a fair instrument?

With the two replacements, yes for the work type it measures, and with one
reservation that is not about fairness but about power: **five of the twenty
have no docstring**, and for three of those the intended call site has to be
inferred from a signature and a module name. Those are the rows where two
competent agents would most plausibly disagree about what "connected"
means, and where the scorer's judgement does the most work.
