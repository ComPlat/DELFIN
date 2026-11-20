# The sub-agent specialist foundation

Measured at `93c498537..f02099cf9` on branch
`guard/one-protected-list-for-the-preset-check`, on the group's workstation
(Linux, 384 cores, bubblewrap available, ORCA installed, 102 ops-server
tools importable).

Four of these numbers CI cannot produce: the registry size depends on a
settings file, the catalogue size on an importable ops server, and the two
control runs compare two trees rather than assert a green one.

## What was driven

```
~/micromamba/envs/delfin-ci/bin/python -m pytest -q -p no:randomly <files>
```

from inside a worktree, and for the registry counts
`subagents._build_preset_registry()` called directly in that worktree with
`HOME` pointed at an empty directory.

## Numbers

| what | before | after |
|---|---|---|
| live preset registry, default installation | 12 | **17** |
| tool-name catalogue (`_known_preset_tool_names`) | 141 | **300** |
| protected prefixes a write scope may not enter | 17 | **43** |
| gate globs the preset check did not cover | 33 | **0** |
| shipped presets with a validation problem | 13 of 13 (unseen) | **0** |
| ops-server tools the catalogue derives | 19 of 102 | **102 of 102** |
| presets validated at load | 0 | all |
| per-request tool schema | 11,678 tokens | 11,721 (budget paid, not raised) |

The registry count was taken with an **empty HOME**. This host's own
`~/.delfin_settings.json` already sets both department switches true, so an
in-place measurement reads 17 before the change as well and proves nothing.
12 was reproduced at `93c498537`, `8b8cc0e68`, `da82bd09a` and `d554c53a7`.

## Control runs

Each new test file was run against the tree *without* its change, in that
tree, so the failure is the measurement and not a claim about one:

| file | at the baseline | with the change |
|---|---|---|
| `test_the_preset_check_uses_the_gates_protected_list.py` | 11 failed, 7 passed | 18 passed |
| `test_the_subagent_schema_follows_the_registry.py` | 5 failed, 1 passed | 6 passed |
| `test_every_shipped_preset_is_valid.py` | 13 failed, 2 passed | 15 passed |
| `test_the_chemistry_departments_are_on_by_default.py` | 5 failed | 5 passed |

The seven that pass at the baseline in the first row are the ordinary write
scopes and the frame files the old copy already covered; the one in the
second row is the non-empty enum, which was already true. Stating which
assertions were already true is the point of running the control at all.

## Across every suite the four streams touch

319 passed: the token budget, the md loader, `sa0_preset_fields`, the
schema-follows-registry pair, the validation pair, the shipped-preset CI
check, the label helper, the departments default, the protected-list pair,
the preset tools suite, `sa0_validate_preset`, `sa0_invariants`,
`sa0_adversarial_invariants`, `test_subagents`. Ruff clean on
`F821,F811,F601,F402,F823`.

## What was NOT verifiable from here

- **Whether a specialist is better than a generalist.** Nothing here
  measures that. 0 of 61 sub-agent calls across three agent waves went to a
  domain specialist, and for five of them the reason was that they were not
  loaded at all. The A/B measurement — same tasks, same model, same brief,
  one arm with the specialist named, scored on work landed and tests green
  rather than on the number of specialist calls — has not been run.
- **The cluster side.** The other wave operator's agents make no sub-agent
  calls by rule, not by preference, so that half of the record says nothing
  either way.
- **A real delegation chain.** `may_delegate_to` is empty in every shipped
  preset, so the graph pass added here has nothing to walk in production;
  it was driven against stand-in presets that declare a ring.
- **`display_name` end to end.** The helper and the seven display sites are
  tested, and the longest shipped label still renders cut to 16 characters
  by a pre-existing control width, which was reported rather than widened.
