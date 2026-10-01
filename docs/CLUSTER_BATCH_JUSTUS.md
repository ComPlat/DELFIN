# Batch-Läufe auf JUSTUS 2 mit `delfin cluster`

`delfin cluster` baut eine Liste von `ID;SMILES`-Zeilen auf einem Slurm-Cluster: mit MANTA
(der ausgelieferten Konstruktion) und, zum Vergleich, mit Architector, molSimplify und epic-MACE.
Jede Array-Task baut einen Shard auf einem vollen Knoten (48 Kerne, bis 72 h), setzt nach einem
Abbruch fort und schreibt ins Archivformat von DELFIN (mehrteilige xyz, Kommentar
`<ID> frame<k> <label>`).

> **Daten:** Das Repository enthält keine Eingaben. Listen, Auswahlen, Specs und Ergebnisse liegen
> nur im Laufverzeichnis (Workspace). Lizenzierte Eingaben (z. B. aus einer Strukturdatenbank)
> gehören nie in ein öffentliches Repository.

## 1. Was die Befehle tun

| Befehl | Zweck |
|---|---|
| `delfin cluster prepare` | Liste (+ optional Auswahl) → Shards, je Shard die Specs (externe Werkzeuge), `manifest.json` mit sha256, Code-Commit, Env-Versionen, Einstellungen |
| `delfin cluster slurm` | schreibt das sbatch-Array-Skript (48 Kerne, `--time=72:00:00`, `%N`), mit `--submit` auch einreichen |
| `delfin cluster run-shard` | die Array-Task: baut einen Shard, ein Kindprozess je System, Zeitlimit je System, fortsetzbar |
| `delfin cluster status` | Fortschritt je Shard und Klasse, unvollständige Shards als `--array`-Angabe |
| `delfin cluster collect` | führt die fertigen Shards zu `archive_<label>/` zusammen, mit `build_/buildtime_/buildmem_/summary_<label>.json` |
| `delfin cluster repeat-stats` | Byte-Gleichheit einer zweimal gebauten Teilmenge (`prepare --repeat N`) |

Gleichwertig: `delfin-cluster …` und `python -m delfin.cluster_bench …` (das Slurm-Skript nutzt
die letzte Form).

**Zeitlimit je System** = ⌈`--timeout` (Vorgabe 21 600 s) × `SPEED_FACTOR`⌉, exakt dezimal
gerechnet. Es ist nicht die Wandzeit des Jobs (72 h). Ein System am Limit wird mit seiner ganzen
Prozessgruppe beendet und als `timeout` geführt; das Limit entscheidet nur, *welche* Systeme
fertig werden, nie die Bytes eines fertigen Systems.

**Klassen** (`build_<label>.json`, `_meta/<ID>.json`):

| Klasse | Bedeutung |
|---|---|
| `ok` | mindestens ein Frame |
| `timeout` | Limit je System erreicht |
| `empty` | Werkzeug lief, kein Frame (auch Absturz) |
| `fail` | das Werkzeug lehnt die Eingabe ab (CN/Zähnigkeit ohne Geometrie …) – Werkzeug-Fehlschlag |
| `not_expressible` | nicht in der Eingabesprache des Werkzeugs ausdrückbar (Schnitt Metall + Liganden + OZ 0..8 unmöglich; MACE: Stellenzahl ohne Geometrie) – **kein** Werkzeug-Fehlschlag, nicht im Nenner |

MANTA liest SMILES direkt und kennt nur `ok/empty/timeout/fail`.

**Determinismus und Provenienz:** `prepare` hält Commit (plus Hash ungesicherter Änderungen unter
`delfin/`), Interpreter und Paketversionen des DELFIN-Envs und des Werkzeug-Envs, sha256 der
Worker, die MANTA-Bauschalter und alle Einstellungen im Manifest fest. `run-shard` prüft vor dem
Bauen sha256 des Shards und der Specs sowie Code und Envs gegen das Manifest und **baut nicht**,
wenn etwas abweicht (Exit 2, Begründung in `refused_<zeit>.json`). `collect` meldet Chunks mit
verschiedenem Code, Env oder Limit als Problem (`n_problems` muss 0 sein).

MANTA baut je System genau wie die Kampagnen-Harness: Bauschalter aus
`delfin.cli_manta.construction_env("champion")`, Nachbearbeitungspässe aus,
`smiles_to_xyz_isomers(..., deterministic=True, max_isomers=100000, quality_mode="extreme")`,
7 Threads je Bau (`--threads`, für Byte-Gleichheit nicht ändern), BLAS/OpenMP 1 Thread,
`PYTHONHASHSEED=0`, `DELFIN_DETERMINISTIC=1`. Gleiche SMILES werden einmal gebaut; die übrigen
IDs bekommen die Frames mit eigener ID im Kopf (`cached` im Status).

## 2. Einmalig auf JUSTUS 2: Workspace, DELFIN, Envs

```bash
ssh <user>@justus2.uni-ulm.de
ws_allocate delfin_batch 60          # Workspace, 60 Tage; ws_extend delfin_batch 60 verlängert
WS=$(ws_find delfin_batch); cd $WS

# DELFIN vom privaten Zweig (nicht von ComPlat/DELFIN)
git clone -b work/2026-09-30-cluster-batch https://github.com/hmaximilian/delfin-backup.git DELFIN
git -C DELFIN log --oneline -1       # Commit notieren, er steht später in jedem Manifest
```

**DELFIN-Env** (Python 3.12.11, exakt die Versionen der MANTA-Referenz):

```bash
# micromamba einmal holen, falls nicht vorhanden:
#   curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba create -y -p $WS/delfin_env -c conda-forge python=3.12.11
$WS/delfin_env/bin/pip install -e "$WS/DELFIN" -c "$WS/DELFIN/env/reference-environment.lock"
DPY=$WS/delfin_env/bin/python
PYTHONNOUSERSITE=1 $DPY -m delfin cluster --help
```

Ohne Netz auf dem Login-Knoten: das Offline-Env des privaten MANTA-Pakets
(`justus_batch/env/install_env.sh $WS/delfin_env`, 143 Pins, Wheelhouse) installieren und danach
`$WS/delfin_env/bin/pip install --no-deps -e $WS/DELFIN`.

**Werkzeug-Envs** (eigene Interpreter, DELFIN wird dort nicht installiert): aus den privaten
Offline-Paketen nach `$WS` kopieren (`rsync`, siehe 3.) und installieren:

```bash
module load compiler/gnu                               # gcc für pynauty (Architector-Env)
bash competitors/env/install_env.sh $WS/bench_env       # Architector 0.0.10 + molSimplify 2.0.0, endet mit "ok"
bash competitors/env_mace/install_mace_env.sh $WS/mace_env   # Python 3.7 + RDKit 2020.09 + epic-MACE efb5778e
```

## 3. Eingaben kopieren (vom Arbeitsrechner)

```bash
R=<user>@justus2.uni-ulm.de:<WS>
rsync -av --progress ~/agent_workspace/SMILES2XYZ/Batch_v2.txt $R/inputs/
rsync -av ~/agent_workspace/justus_batch/competitors/{heldout_20k.txt,repro500.txt} $R/inputs/
rsync -av ~/agent_workspace/justus_batch/competitors/specs/specs_heldout_20k.jsonl $R/inputs/   # optional, s. u.
rsync -av ~/agent_workspace/justus_batch/competitors/{env,env_mace} $R/competitors/          # Offline-Installer
```

Specs: `prepare` schneidet jedes System mit `delfin.common.external_builders.split_complex_smiles`
im DELFIN-Env. Die kanonische Schreibweise der Ligand-SMILES hängt von der RDKit-Version ab
(lokal gesehen: 2 von 10 Systemen mit anderer, chemisch gleicher Schreibweise zwischen RDKit
2025.9 und 2026.03). Wer exakt die Specs des lokalen Zensus will, übergibt sie mit `--specs`.

## 4. Vorbereiten

```bash
cd $WS; export PYTHONNOUSERSITE=1; DPY=$WS/delfin_env/bin/python
B=$WS/inputs/Batch_v2.txt

# MANTA, ganze Liste (Shards à 500, gleiche SMILES zusammen), 300 Systeme doppelt für die Byte-Statistik
$DPY -m delfin cluster prepare --tool manta --input $B --run-dir $WS/runs/manta_batch --repeat 300

# Konkurrenten auf der Held-out-Auswahl (Reihenfolge der Auswahl bleibt, Shards à 250)
for T in architector molsimplify; do
  $DPY -m delfin cluster prepare --tool $T --input $B --select $WS/inputs/heldout_20k.txt \
       --specs $WS/inputs/specs_heldout_20k.jsonl \
       --tool-python $WS/bench_env/bin/python --run-dir $WS/runs/${T}_h20k --repeat 500
done
$DPY -m delfin cluster prepare --tool mace --input $B --select $WS/inputs/heldout_20k.txt \
     --specs $WS/inputs/specs_heldout_20k.jsonl \
     --tool-python $WS/mace_env/bin/python --run-dir $WS/runs/mace_h20k --repeat 500
```

Optionen: `--mode` (Architector `full`/`default`, MACE `paper`/`extended`), `--timeout`,
`--shard-size`, `--workers` (MANTA 36 × 7 Threads, sonst 48), `--label`. Ein Laufverzeichnis wird
nie überschrieben.

## 5. Einreichen: erst Pilot, dann Hauptlauf

```bash
# Pilot: 2 Shards je Werkzeug, eigener Ausgabeordner (--run pilot), ungeskaliertes Limit
for R in manta_batch architector_h20k molsimplify_h20k mace_h20k; do
  $DPY -m delfin cluster slurm $WS/runs/$R --run pilot --array 0-1 --submit
done
# Wandzeit/Status ansehen, Faktor F gegen lokal bestimmen, dann Hauptlauf mit demselben F für alle:
F=1.3
for R in manta_batch architector_h20k molsimplify_h20k mace_h20k; do
  $DPY -m delfin cluster slurm $WS/runs/$R --speed-factor $F --throttle 20 --submit
  $DPY -m delfin cluster slurm $WS/runs/$R --set repeat --speed-factor $F --submit   # zweiter Bau
done
```

- `--throttle N` begrenzt gleichzeitige Tasks (`%N`; 1920 Kerne je Nutzer = 40 Knoten insgesamt,
  auf die Werkzeuge aufteilen).
- `--partition`, `--account`, `--mem`, `--time` (Vorgabe `72:00:00`), `--setup 'module load …'`
  nach Bedarf. Ohne `--partition` fragt `--submit` mit `sbatch --test-only`, welche konfigurierten
  Partitionen den Job nehmen (`delfin.slurm_submit`).
- Ohne `--submit` wird nur `RUN/slurm/<tool>_<set>_<run>.sbatch` geschrieben; einreichen mit
  `sbatch <datei>`.

## 6. Überwachen, fortsetzen

```bash
$DPY -m delfin cluster status $WS/runs/architector_h20k          # -v: je Shard
squeue -u $USER -r | head; sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,ExitCode
```

Eine Task, die an der Job-Wandzeit oder an einem Knotenfehler stirbt, behält alle fertigen
Systeme. `status` zeigt die unvollständigen Shards als `--array`-Angabe; einfach erneut
einreichen (mit **demselben** `--speed-factor`, sonst verweigert `run-shard`):

```bash
$DPY -m delfin cluster slurm $WS/runs/architector_h20k --speed-factor $F --array 3,17-19 --submit
```

## 7. Zusammenführen und zurückholen

```bash
for R in manta_batch architector_h20k molsimplify_h20k mace_h20k; do
  $DPY -m delfin cluster collect $WS/runs/$R                       # -> RUN/collected/<label>/
  $DPY -m delfin cluster collect $WS/runs/$R --set repeat          # zweiter Bau der Teilmenge
  $DPY -m delfin cluster repeat-stats $WS/runs/$R                  # identisch / verschieden
done
```

`summary_<label>.json`: Klassen, `coverage_of_expressible`, fehlende Shards (`resubmit_array`),
Limit und Provenienz; `n_problems` muss 0 sein.

Zurück auf den Arbeitsrechner (nie in ein bestehendes Ziel):

```bash
rsync -av <user>@justus2.uni-ulm.de:<WS>/runs/ ~/agent_workspace/justus_returned/runs/
```

## 8. Hinweise

- **Workspace-Ablauf:** Workspaces auf JUSTUS 2 verfallen (`ws_list` zeigt die Restzeit);
  rechtzeitig `ws_extend delfin_batch 60` oder die Ergebnisse zurückholen. Nach dem Ablauf ist
  der Inhalt weg.
- **Kontingent:** MANTA auf ganz Batch_v2 rund 130 000 Kernstunden (Spanne 85–150 k);
  Architector auf heldout_20k bis ≈ 55 000, molSimplify ≈ 3 000, MACE ≤ 24 000 (Pilot zeigt,
  wo es liegt).
- **Nicht deterministisch:** Architector und molSimplify liefern schon lokal von Lauf zu Lauf
  andere Isomere bzw. Koordinaten; `repeat-stats` misst das. epic-MACE und MANTA waren lokal
  byte-gleich reproduzierbar.
- **Gleiches Limit, ungleicher Aufwand:** MANTA baut mit 7 Threads je System, die anderen
  einfädig.
- Arbeitsverzeichnisse der Werkzeuge liegen unter `chunk_NNNN/work/` (nach Erfolg gelöscht);
  `DELFIN_CLUSTER_WORK_ROOT` legt sie z. B. auf `$TMPDIR`.
- Eine Dashboard-Ansicht gibt es noch nicht; `delfin cluster status --json` liefert dieselben
  Zahlen maschinenlesbar.
