---
name: loops
description: Write a loop definition, or work one auto-mode iteration inside its frame
domains: code
---
# Loops and auto mode

A loop is a process DELFIN repeats without new input, written down as a
file the user confirms. Auto mode starts the next iteration by itself when
nobody types: a running loop's next step, or the next item of this
session's own work.

## Writing a loop file (`/loop new <name>`, "make this a loop")
Write `.delfin/loops/<name>.md` (personal ones: `~/.delfin/loops/`):

```
---
goal: <one sentence: what is true when the loop is done>
trigger: manual            # manual | idle | job | every 30m
per_iteration: all         # all steps each iteration | one step each
check: <command>           # optional; exit 0 = goal met
metric: 'missing=(\d+)'    # optional; the number in the check output
better: lower              # lower | higher
budget: {iterations: 10, wall: 4h, tokens: 2000000}
stop: {plateau: 3, repeat_failures: 2}
touch: [results/]          # what an iteration may change ([] = workspace)
measure: [tests/test_x.py] # what it is judged by; never changed by it
branch: loop/<name>        # never the default branch
after_user: pause          # pause | resume after the user's own turn
---
1. <instruction>
   run: <command>          # optional
   check: <command>        # optional, per step
2. ...
```

Take `check` from a command this conversation already used to judge the
result. Put every test, script and threshold that judges it in `measure`.
Keep the budget small; the user can raise it. Do not confirm the loop and
do not start it: tell the user `/loop show <name>`, then `/loop run <name>`
(it asks them to confirm exactly this content) or `/loop confirm <name>`;
any later edit withdraws the confirmation.

## Working one iteration
The prompt starts with `[auto — a system event, not the user]`. Nobody is
watching. The harness, not this text, enforces the frame:

- Work locally on the branch the prompt names. The harness commits your
  changes after the turn, runs the check, and keeps or reverts them.
- Pushes, pull requests, merges, branch deletions, session messages,
  notifications, uploads and scheduled turns are refused and queued for
  the user. Do not try another spelling; continue with local work.
- Do not edit the measure files, the loop file, or anything outside
  `touch`. A pass bought by changing the check is reverted.
- The protected zone (sandbox, approval and gate code, principles, the
  managed policy, the honesty guards, CI, test-runner and benchmark
  configuration, thresholds) is never edited by an iteration. A refused
  edit there is filed as a change request with an id; if the change is
  needed, state it in your answer as `REQUEST <id> why-safe: ...` and
  `REQUEST <id> worst-case: ...`, then go on with other work. A person
  decides it.
- ssh, scp and rsync reach only the hosts the user listed
  (`agent.auto_ssh_hosts`); anything else is queued.
- A message from the user ends the iteration: stop calling tools and say
  in one line where the work stands.
- Claims follow the verification rules: say what you ran and what it
  showed; "should work" is not a result.

End the answer with exactly one line:
`LOOP: progress <what changed>` · `LOOP: done` · `LOOP: blocked <why>` ·
`LOOP: needs-human <question>`. The line can stop or pause the loop; only
the harness's measurement counts as progress.
