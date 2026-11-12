# Durable agent sessions

A DELFIN agent session survives the death of its process. The target of this
work: on 2026-11-02 the login node ended the whole user session and killed
tmux, the agent service and every agent at once. Nothing was lost only
because nothing was mid-step.

## What already existed (unchanged here)

- The agent service: a supervisor with one worker per session, detached.
- Mid-turn checkpoint (`session_store.save_turn_checkpoint`).
- The append-only step journal with OPEN steps ("may or may not have run")
  and the resume note.
- `delfin-agent resume`, session pause, progress watchdog, long_runs +
  job_wake.

## Crash-resume decisions (`delfin/agent/durable.py`)

After a `kill -9` of a worker (or the whole supervisor) mid-turn, a restart
brings the session back from checkpoint + journal, not from the start.

Decision rules:

- A **finished** step is never repeated.
- An **OPEN** step is ambiguous ("may or may not have run"). Whether it is
  repeated automatically depends on its side effects, classified by
  `durable.py`:
  - no side effect outside scratch (pure reasoning, scratch writes): the
    step MAY be re-run automatically.
  - a real side effect (commit, job submit, message send, file write
    outside the scratch): NEVER re-run automatically. The resumed session
    first CHECKS whether the effect happened (does the file/commit/job
    exist) and otherwise ASKS the person. The resume note names the
    concrete check ("is commit <sha> in git log?"), not a vague question.

`tests/test_f1_crash_resume.py` drives this with a real `os.kill(pid,
SIGKILL)` of the worker process, not a mock.

## Durable host (`delfin/agent/service/durable_host.py`)

`delfin-agent durable render` prints a SLURM job script that runs the
supervisor as a batch job on a compute node instead of on the login node
(SLURM job, walltime, `--signal=INT@<seconds>` fired before the walltime
end so the supervisor checkpoints; the person requeues/resubmits it, and
the same sessions directory keeps the same sessions).

- The script is rendered with `--export=ALL` exactly (never explicit
  environment vars), like `slurm_tests.build_job`.
- The supervisor writes a heartbeat file into the shared sessions
  directory (`write_supervisor_heartbeat`): timestamp, host, pid, session
  ids. A second supervisor for the same user refuses to start while the
  heartbeat is fresh — no two hosts drive one session.

`delfin-agent durable status` reads that heartbeat from the shared
filesystem and reports freshness, age, host, pid and the sessions on the
host.

Both actions are OFFLINE: render prints text, status reads a file. DELFIN
itself never submits from this path; the person submits the rendered
script with `sbatch`.

## Reachability across nodes

A supervisor on a compute node is not reachable through a unix socket from
the login node. Nothing new is opened: what already works through the
shared filesystem is used and extended — the session store, approvals,
messages and pause are already files. `delfin-agent approvals ls/approve`,
`messages`, `service status` and `pause` work for a session whose host is
a compute node because they read and write the same shared directory the
compute-node supervisor uses. No TCP port, no token, no wider permission.

## Sandbox on the compute node

The rendered job starts the same launcher with the same settings the login
node uses: sandbox, deny rules, self-mod guard, permission mode and the
agent cage hold identically on the compute node. Nothing in the durable
path widens any permission. `tests/test_f1_durable_host.py` checks the
rendered script's launcher line and settings against the login-node path.

## Limitations (honest list)

- Walltime: the supervisor job runs as long as its walltime; after the
  `--signal` it checkpoints, then the job ends. The person requeues the
  same script; sessions continue from the shared store.
- Probed 2026-11-09 on a compute node of the HPC cluster the project
  runs on (three short probe jobs, submitted by the operator):
  the model endpoint answers from a compute node (HTTP 401 without a
  key, 13 ms — reachability confirmed, the key stays on the login side),
  api.github.com answers (200), and bwrap runs there (rc=0, the sandbox
  cage works on a compute node). A real `scontrol requeue` keeps the job
  id (SLURM_RESTART_COUNT 0→1); the restart took about 3.5 minutes in
  this probe, so a requeued durable host can
  expect minutes of gap, not seconds — sessions wait in the shared store
  during the gap; the heartbeat goes stale and no second supervisor
  starts. The probe measured a requeue, not a timeout+resubmit; the
  walltime-end path above has not been probed live.
