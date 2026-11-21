# The verdict handoff gets a schema the check already reads

Measured at the baseline `5f1975075` and the branch head `50153975a`
(`swarm/night5/n5o`), on the group's workstation (Linux, 384 cores,
`bubblewrap` present but the sub-agent "cage supervisor" cannot take over
a child here — the child reports no listener number and the command comes
back `exit_code 126, "the command was not run"`).

One of these numbers CI cannot produce: the three *driven* verifier tests
run `run_subagent` through the real executor, and one of the scripted
tool rounds is a live `run_tests` call. That call goes to the test's own
cage supervisor, which this host does not hand back to the child, so the
round fails with `exit_code 126` whether or not the change is present.
CI runs those three against a working supervisor and gets a different
number; on this host the only honest record is that the two arms agree.

## What was driven

```
python -m pytest tests/test_the_verifier_tries_to_refute_a_claim.py -k driven -q
```

from inside a worktree checked out at each tree, the loaded module's
`__file__` confirmed to be inside that worktree before each run. The other
numbers in this change (the red/green of the new tests, the
`validate_preset` result, the preset load) are produced by CI on every
commit and are not repeated here.

## Numbers

| what | at the baseline 5f1975075 | at the head 50153975a |
|---|---|---|
| driven verifier tests, `-k driven` | **3 failed**, 1 passed, 11 deselected | **3 failed** (same three), 1 passed |
| the three, named | `test_driven_it_runs_the_check_and_cannot_edit_a_file`, `test_driven_without_a_callback_a_shell_write_does_not_land`, `test_driven_without_a_callback_a_writing_preset_keeps_the_contract` | the same three, same failure text |
| failure signature (all three) | cage supervisor no listener number; `exit_code 126, "the command was not run"`; or `planted*.txt` not created | identical |
| every other consumer test (18 files, 19 tests' worth of the rest) | 455 passed of 458 | 458 passed, **0 new failures** |

The two arms fail the *same three tests with the same text*. A change
that broke a driven test would fail it at the head and pass it at the
baseline; these fail at both, so the failure is the host's supervisor,
not this change.

## What was NOT verifiable from here

- **Whether the three driven tests pass on CI.** This host's cage
  supervisor cannot take over a child, so I could not drive the `run_tests`
  round to a real `"passed": 1` here. CI, with a working supervisor, is
  the only place that number will exist; it is not recorded as green.
- **A real delegation to a verdict preset.** The schema's effect on the
  return is driven with a scripted fake client (the mechanism CI itself
  uses); a live model call that fills the `live_run` block by running a
  suite was not made, and the four wave failures that motivated the field
  are other operators' record, not reproduced here.
- **The head_sha the caller re-checks.** `head_sha` is a string the
  delegate must fill; nothing here runs the caller's `git rev-parse`
  against it, because there is no live run to take a head from.
