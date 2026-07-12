# RESURRECT — DELFIN self-optimizing loop + agent memory on a NEW machine

> Disaster recovery. Reconstitute the FULL setup (DELFIN construction + MANTA2 harness + the agent's
> memory) on a fresh machine, so the self-optimizing loop and the accumulated project knowledge continue.
> **PRIVATE / infra — never ship in a public DELFIN release.**

## The backed-up state (all PRIVATE)
1. **`hmaximilian/delfin-backup`** (GitHub, PRIVATE) — this DELFIN construction repo, full git history.
   (DELFIN's public remote is `ComPlat/DELFIN`; NOTHING is pushed there without explicit user approval.)
2. **`hmaximilian/agent_workspace`** (GitHub, PRIVATE) — `MANTA2/` harness, `quality_framework/`,
   `CLAUDE_MEMORY/` (the agent's memory), `Batch.txt`, SMILES pools, `BACKUP/backup.sh`.
3. **`/storage/MANTA_SICHERUNGEN/`** (machine-local) — 2-hourly git bundles + session transcripts.
   Transcripts stay local-only (may contain secrets); memory + code go to GitHub.

## Steps on a new machine
1. **Clone the two private repos** (adjust the home layout — see paths below):
   - `git clone git@github.com:hmaximilian/agent_workspace.git ~/agent_workspace`
   - `git clone https://github.com/hmaximilian/delfin-backup.git ~/ComPlat/DELFIN`
   - add DELFIN's remotes back: `origin` = `git@github.com:ComPlat/DELFIN.git` (PUBLIC, never push),
     `private-backup` = the delfin-backup URL.
2. **Restore the agent's memory** (THIS is the continuity — a fresh Claude session reads it and continues):
   - `cp -r ~/agent_workspace/CLAUDE_MEMORY/* ~/.claude/projects/<project-dir>/memory/`
   - `MEMORY.md` is the index; the `*.md` topic files hold the detail + all the rules/findings/plan.
3. **Recreate the Python env** (from this dir):
   - `micromamba create -n delfin -f RESTORE/environment_delfin.yml`
     (or `environment_delfin_minimal.yml` = explicit installs only; `pip_freeze_delfin.txt` = pip layer).
   - loop.py expects `~/micromamba/envs/delfin/bin/python`.
4. **Paths:** ~31 code files hardcode `/home/qmchem_max/...`. Easiest = create user `qmchem_max` with the
   same home layout (`~/agent_workspace`, `~/ComPlat/DELFIN`, `~/CCDC/...`, `~/micromamba/...`, `/storage/...`).
   Otherwise sed-remap those roots. (Long-term: make the code path-free.)
5. **CCDC (ORACLE ONLY — not needed to CONSTRUCT):** install the CSD Python API + a valid CCDC licence.
   `ccdc_tools/ccdc_env.py` resolves it via `CCDC_WRAPPER`/`CSD_RUN_PYTHON_API` env → PATH → `CSD_HOME`.
   Only the 5th-signal (mogul realism) + CCDC ground-truth need it; shipped DELFIN construction is CCDC-free.
6. **Ground-truth data (validation only):** CCDC `clean_v2` (~307k) is re-derivable from the CSD;
   `Batch.txt` + SMILES pools are in `agent_workspace/MANTA2`.
7. **Restart the 2-hourly backup:** re-add the cron `17 */2 * * * ~/agent_workspace/BACKUP/backup.sh`.
8. **Smoke test:** `cd ~/agent_workspace/MANTA2 && bash EYE/capped.sh -- python harness/loop.py --pool full:5 --label smoke --on "" --off "" --det-n 0`.

## What "resurrect the agent" really means
Claude starts fresh each session — there is no model checkpoint to restore. Continuity lives ENTIRELY in
`CLAUDE_MEMORY/` (MEMORY.md + topic files). Restore those and a new session picks up the project state,
the user's rules, the findings, and the plan. The memory IS the thread; keep it backed up (it is, 2-hourly).
