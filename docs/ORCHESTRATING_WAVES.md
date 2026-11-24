# Orchestrating waves of DELFIN agents

This is for whoever runs many DELFIN agents at once: a person at the
terminal, a DELFIN lead session, or another coding agent connected over
MCP. It describes how waves have been run on this project since
2026-11-05, and what each step guards against. Every rule below was
written after a measured defect.

## When a wave pays off

Run a wave when the work splits into packages that touch different files,
each small enough for one session, each with a measurable result. Do it
yourself when the work is one change in shared files, smaller than the
coordination it needs, or a question one measurement answers.

Give well-specified work to the wave and keep the expensive reasoning for
yourself: finding the root cause where members failed, and checking what
they produced.

## Connecting

`delfin-agent swarm wave mcp` serves every verb over MCP on stdio (server
name `delfin-swarm`). Register that command with your client once; the
same verbs exist on the command line as `delfin-agent swarm wave VERB`.

Settings (`agent.swarm` in the DELFIN settings):

| key | meaning | default |
|---|---|---|
| `waves_root` | where runs, worktrees and logs live (never under `~/.delfin`) | `~/delfin-waves` |
| `harness` | a separate checkout members run DELFIN from, kept on `origin/main` | the installed DELFIN |
| `max_parallel` | members at once per run | 4 |
| `provider_caps` | members at once per provider, across all runs, e.g. `{"kit": 3}` | none |
| `member_cpu`, `member_memory`, `member_tasks` | each member's scope ceilings | 200 %, 6G, 1024 |
| `git_identity` | `{"name", "email"}` every landed commit must carry | the swarm job identity |
| `test_python` | interpreter a landed branch's tests run on | this one |
| `github_repo`, `land_base` | where PRs are opened | the repository's own, `main` |

## The shape of a wave

`wave_compose` (Python: `delfin.agent.swarm.compose.compose`) turns a
list of packages into the whole task list. It starts nothing:

1. **Arms**: every package is worked on by members of *different* model
   families, independently. Where they agree, the result is robust; where
   they differ, the difference is the finding.
2. **Integrator**, one per package: it starts when its arms have ended, or
   when a quorum has ended and waited. It re-runs every arm's tests on the
   base (they must fail there) and on the branch (they must pass), keeps
   the better-measured work item by item, and asks for real runs.
3. **Lander**: bundles the landing branches into one branch for one full
   CI run. A serial lander starves every other merger.
4. **Reviewers**: four, of another family than the lander. They are
   read-only, each with one focus: behaviour, tests, safety, claims.
5. **Fixer**: works through every confirmed finding and writes the final
   pull-request text.

A package is `{"id", "title", "body"}`. The body says:

- what to do, and which files to touch;
- which suites to run before and after;
- what to measure;
- which real-run probes would show the change.

The delivery rules (`compose.DELIVERY_RULES`) are appended to every
member's brief.

## Running and watching

1. `wave_start(run, tasks, repo)` starts what is ready.
2. `wave_guard_install(systemd=true)` runs the tick every minute as a
   systemd user service, which restarts when it fails and comes back after a
   reboot. Without systemd, `wave_guard_install()` makes it an entry of the
   scheduler daemon, and `delfin-agent swarm wave tick --every 60` runs it
   under any supervisor (tmux on a login node, a batch job). The tick:
   - starts tasks as their dependencies end;
   - runs the stage: real runs outside the member's cage, the resume that
     reads them, the final mark;
   - resumes members an endpoint outage ended;
   - records events;
   - guards the host.
3. `wave_watch(runs, timeout_s, cursors)` blocks until something happens:
   a member ended (and how), went silent, left a traceback, or became
   final. Loop on it and pass the cursors back.
4. `wave_status`, `wave_log` and `wave_review` show what a member did and
   wrote.
5. `wave_resume(run, member, note)` steers a member with your review.
   `wave_add_tasks` appends a later stage.

Capacity is shared. Every orchestrator on an account sees the same runs,
and the provider caps count all of them. Agree on a split with the others
before you raise a limit, and never touch another orchestrator's runs or
worktrees.

## Your review: the members' reports are not evidence

Before anything lands, check it yourself:

- Re-run each new test without its fix: it must fail for the reason it
  names. Then run it with the fix: it must pass.
- Re-measure every number that is claimed.
- Read the code. A test that calls a helper instead of the product entry
  point proves that the helper exists.
- Compare the real-run table with what the text says it shows.
- Search the diff and the PR text for account names, home paths, host
  names and keys.

## Landing

`wave_land()` is a dry run:

1. the branch checks: commits, identity, home paths, forbidden paths;
2. the tests the branch adds, on `test_python`;
3. an independent review by a calibrated pair of reviewers from another
   model family.

`wave_land(push=True)` then pushes the ready branches and opens their
PRs. The reviewers are calibrated first: `delfin-agent swarm wave land
--calibrate`, then `--calibrate-pair`. An uncalibrated reviewer refuses
every branch. That is the design: a review that did not happen is never
reported as one.

Merging stays with a person.

## Safety

- Members run caged: no keys in their shell, no push URL, outbound traffic
  only through the egress proxy, and every socket they bind on the
  loopback.
- The tick checks this every minute. If a member's process listens beyond
  the loopback, or the host leaves its load and memory envelope on two
  ticks, every member is stopped. The swarm then holds until a person runs
  `wave_clear_stop`.
- `wave_emergency_stop(why)` does the same by hand. Use it on any doubt;
  worktrees, commits and logs stay.

## Lessons the tools now carry

- A member whose launcher has not reaped it is a zombie: it has ended.
  Count it as ended.
- Take the state lock for each change, not for a starter's whole life. A
  lifetime lock deadlocked the stage.
- Never edit a tree while a suite runs over it.
- Never use `git stash` in a shared checkout.
- Never run two writing agents in one tree.
- Commit before every validation run, so its code state is known.
- A real run imports DELFIN from its code copy (`python -P`). Without
  that, a repository-copy workspace imports itself.
