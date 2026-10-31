# Schedules

A schedule runs a task at set times: in the agent service, on the host where it
was made, as a session of its own. It works the same on a Pi, a server and a
cluster login node, and it needs no root.

## Making one

| Where | How |
|---|---|
| A session (TUI, `delfin-agent chat`) | `/schedule add täglich 07:30 Fasse die Slurm-Jobs zusammen` |
| A shell | `delfin-agent schedule add mo-fr 18:00 --cwd ~/proj Prüfe die Rechnungen` |
| The phone, schedules chat | `neu alle 2h Schau nach neuen Ergebnissen`, then `ok N` |
| The phone, a session's chat | `plan täglich 07:30 Lies die Logs` (runs in that session's folder), then `ok N` in the schedules chat |
| Asking the agent | "jeden Tag um 7:30 ..." -- it proposes one (`cron_create` with `at`); you confirm with `/schedule ok N` |

Times, German or English, then the task:

    täglich 07:30 · daily 7:30 · jeden Tag um 7 Uhr · 07:30, 18:00
    mo-fr 18:00 · werktags · wochenende · mo,mi,fr 07:00 · montags 8 uhr
    alle 2h · every 30m · alle 15 min · stündlich       (at least 5 minutes)
    30 7 * * 1-5                                         (cron: min hour dom month dow)
    ... tz=Europe/Berlin                                 (default: the host's zone)

The time is kept in its zone: 07:30 stays 07:30 across summer and winter time.
A time inside the skipped hour runs once that day, an hour later; one inside the
repeated hour runs once.

## Watchers: when something happens

Instead of a time, an event; the task after it runs when it happens, and
without a task the event is only reported:

    wenn job 123456 fertig Fasse die Ergebnisse zusammen   (Slurm: squeue, then sacct)
    when pid 4711 exits
    wenn datei out/result.xyz erscheint Werte sie aus
    when file calc.log changes every 2 min Prüfe auf SCF-Fehler
    wenn ausgabe "squeue -u max" sich ändert

* The agent service looks every 60 s (`alle|every N min`, at least 30 s)
  without a model; only the event starts a run, its prompt naming what was
  seen. "fertig", "endet", "erscheint" fire once (then the watcher is listed
  as done); "ändert sich" fires on every change, the first look only learns
  the state, and with a task at most once per 5 minutes.
* A command is one of a fixed set of reading commands (squeue, sacct, sinfo,
  qstat, bjobs, df, free, uptime, nvidia-smi, ls, tail, head, wc, grep, du,
  stat, cat, find, checksums), split without a shell, run with a clean
  environment (no credential variables), in the folder, for at most 20 s;
  options that read or write a file of their own choosing (`grep -f`,
  `find -exec`, ...) are refused.
* Every path -- the watched file, a command's file argument -- must lie in
  the schedule's folder, links resolved, when it is made and at every look;
  a watcher that cannot work is paused and that is reported. Commands without
  paths may not carry a `/`. So a watcher, from whichever surface, can neither
  change anything nor read past its folder.

## Managing

`/schedule` (or `delfin-agent schedule`) lists them with the next run, and:

    show N · ok N · pause N · resume N · delete N · run N (now) · log [N]
    add options: --usd 1 (a run's ceiling) --day-usd 2 --minutes 30
                 --notify always|errors|never --cwd PATH

On the phone, in the chat "<machine> · Zeitpläne" (`zeitpläne` opens it):
`pläne`, `neu`, `ok N`, `pause N`, `weiter N`, `jetzt N`, `löschen N`,
`zeige N`, `hilfe` -- in English `list`, `new`, `ok`, `pause`, `resume`, `now`,
`delete`, `show`, `help`.

## What a run is

* One session in the agent service, in the schedule's folder, the task as its
  first prompt. The agent view lists it (`schedule #N · ...`); a window can
  open it, read it and go on with it. Finished, it sleeps; the run before it
  of the same schedule leaves the list (it stays in the session store).
* How it ended goes to `~/.delfin/user_schedules/runs.jsonl`, the TUI's
  attention list and the phone's schedules chat (`notify`: always, errors,
  never). `/schedule log` shows the last runs.
* A missed run (host or service off) is made up once if it is less than six
  hours late; later, it is skipped and that is reported. Missed runs are never
  stacked, and a run never starts while the previous one still runs.

## Security

* **Only you make a schedule run.** What you write in a session, a shell or on
  the paired phone is yours. The agent can only propose one; it runs after your
  `ok`, which records a hash of what it does, where, when and at what cost. Any
  change after that -- through any route, an edited file included -- withdraws
  it until you confirm again. A phone `ok` confirms only the content the phone
  was shown. At most five proposals wait at once.
* **The store is yours alone.** `~/.delfin/user_schedules` is 0700 with 0600
  files, masked in every cage an agent's command runs in (like the credential
  store), refused to the agent's file tools, and `delfin-agent schedule`
  refuses every change from inside an agent session.
* **The same cage as any session.** A run is a session with the scheduled
  posture (`agent.scheduler.permission_mode`, default `default`), never bypass;
  it asks before what it would ask in a window. An approval nobody answers is
  reported ("waits for an approval") and the run is stopped at its time limit.
* **Ceilings.** Each run: $1 and 30 minutes by default (`--usd`, `--minutes`),
  enforced by the session and by the service. Each schedule: $2 a day
  (`--day-usd`); reached, no further run that day, and that is reported. A
  model with no published rate has no measured spend: the time ceiling bounds
  it.
* **One host.** A schedule runs only on the host it was made on, so a shared
  home on a cluster does not run it on every login node.
* **The phone** can make schedules only in the host's phone folder
  (`DELFIN_MESSENGER_WORKSPACES`/zeitplaene) or in an enabled session's own
  folder -- it never names a folder -- at most five in ten minutes.

## Keeping the service up

The service runs the schedules and stays up while this host has one to run.

* Pi: `scripts/pi/pi_messenger.sh` installs it as a system unit.
* Server or cluster, as your user: `delfin-agent schedule autostart` prints a
  systemd `--user` unit and, where that is not available, a crontab line that
  starts the service when it is not running. It changes nothing itself.

Rejected: keeping these in `~/.delfin/cron.json` beside the agent's own
wake-ups. A build from before this one reads that file and, meeting a kind it
does not know, fires the entry and deletes it -- and on a cluster every login
node may run a different build.
