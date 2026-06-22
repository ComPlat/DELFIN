# DELFIN auf BwUniCluster 3.0 - Komplette Anleitung

Diese Anleitung beschreibt, wie du DELFIN auf dem **BwUniCluster 3.0** (KIT) optimal einrichtest und verwendest.

## Cluster-Informationen

- **Name**: bwUniCluster 3.0 + KIT-GFA-HPC 3
- **Standort**: Karlsruher Institut für Technologie (KIT)
- **Login**: `ssh username@uc3.scc.kit.edu`
- **Dokumentation**: https://wiki.bwhpc.de/e/BwUniCluster3.0

### Hardware-Spezifikationen

- **CPU Nodes**: 40 Cores pro Node (einige mit 28 Cores)
- **Memory**: Standard 96-128 GB, High-Memory bis 1 TB
- **Storage**: BeeOND für schnelles lokales I/O
- **Hyperthreading**: Verfügbar (sollte für HPC mit `--threads-per-core=1` deaktiviert werden)

## Schritt 1: Erste Installation

### 1.1 Auf den Cluster einloggen

```bash
ssh username@uc3.scc.kit.edu
```

### 1.2 Verfügbare Module prüfen

```bash
# ORCA suchen
module avail orca

# Python suchen
module avail python

# Beispielausgabe könnte sein:
# chem/orca/6.0.0
# chem/orca/6.1.1
# devel/python/3.11
```

### 1.3 DELFIN installieren

```bash
# Python-Modul laden
module load devel/python/3.11  # Version anpassen!

# Virtual Environment erstellen
cd ~
python -m venv delfin_env

# Aktivieren
source delfin_env/bin/activate

# DELFIN installieren
pip install --upgrade pip
pip install delfin-complat

# Prüfen
delfin --version
```

### 1.4 Installation dauerhaft machen

Füge zu deiner `~/.bashrc` hinzu:

```bash
# DELFIN Setup
alias delfin-activate='source ~/delfin_env/bin/activate'
```

Dann:
```bash
source ~/.bashrc
delfin-activate  # Aktiviert DELFIN Environment
```

## Schritt 2: Projekt einrichten

### 2.1 Arbeitsverzeichnis erstellen

```bash
# In deinem Home oder Workspace
cd $HOME  # oder cd /path/to/your/workspace

# Projektordner erstellen
mkdir -p delfin_projects/mein_molekuel
cd delfin_projects/mein_molekuel
```

### 2.2 CONTROL.txt und input erstellen

```bash
# Wenn du eine XYZ-Datei hast
delfin --define=meine_struktur.xyz

# Oder manuell
delfin --define
```

### 2.3 CONTROL.txt für BwUniCluster anpassen

Wichtige Einstellungen für den Cluster:

```ini
# ========================================
# CONTROL.txt für BwUniCluster 3.0
# ========================================

# Ressourcen (wird automatisch von SLURM erkannt)
pal_jobs=40                     # Wird von SLURM_CPUS_PER_TASK überschrieben
orca_parallel_strategy=auto     # MPI + OpenMP

# Workflow
method=OCCUPIER
OCCUPIER_method=auto
OCCUPIER_tree=deep
calc_initial=yes
oxidation_steps=1,2,3
reduction_steps=1,2,3
parallel_workflows=yes

# Error Recovery (WICHTIG für Cluster!)
enable_auto_recovery=yes
max_recovery_attempts=1
enable_job_timeouts=no          # Für lange Rechnungen

# Optional: XTB/CREST
XTB_OPT=yes
CREST=no
```

## Schritt 3: Job-Submission

### 3.1 Submit-Skript kopieren

```bash
# Standard Job (40 Cores, 120GB, 48h)
cp ~/ComPlat/DELFIN/examples/example_Job_Submission_Scripts/bwunicluster3_delfin.sh .

# Oder High-Memory Job (40 Cores, 500GB, 72h)
cp ~/ComPlat/DELFIN/examples/example_Job_Submission_Scripts/bwunicluster3_delfin_highmem.sh .
```

### 3.2 Skript anpassen (falls nötig)

Öffne das Skript und prüfe:

```bash
nano bwunicluster3_delfin.sh
```

Wichtige Parameter:

```bash
#SBATCH --partition=cpu         # Optionen: dev_cpu, cpu, highmem
#SBATCH --cpus-per-task=40      # 40 oder 28 je nach Node
#SBATCH --mem=120G              # Speicherbedarf
#SBATCH --time=48:00:00         # Max. Laufzeit

# Module Namen prüfen/anpassen
module load chem/orca/6.1.1     # Deine ORCA Version
module load devel/python/3.11   # Deine Python Version
```

### 3.3 Job einreichen

```bash
# Job starten
sbatch bwunicluster3_delfin.sh

# Output sollte sein:
# Submitted batch job 12345678
```

### 3.4 Job überwachen

```bash
# Job-Status prüfen
squeue -u $USER

# Output live anschauen
tail -f delfin_*.out

# Job abbrechen falls nötig
scancel 12345678  # Job-ID einsetzen
```

## Schritt 4: Partitions-Übersicht

### Entwicklung (kurze Tests)

```bash
#SBATCH --partition=dev_cpu
#SBATCH --time=00:30:00
```

- Schneller Queue-Durchsatz
- Max. 30 Minuten
- Für Tests und kurze Rechnungen

### Standard (normale Rechnungen)

```bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=40
#SBATCH --mem=120G
#SBATCH --time=48:00:00
```

- Für die meisten DELFIN-Rechnungen
- 40 Cores, 120 GB Speicher
- Bis 48 Stunden

### High-Memory (große Systeme)

```bash
#SBATCH --partition=highmem
#SBATCH --cpus-per-task=40
#SBATCH --mem=500G
#SBATCH --time=72:00:00
```

- Für sehr große Moleküle
- Korrelierte Methoden
- Bis 1 TB Speicher möglich

## Schritt 5: Best Practices

### 5.1 BeeOND für schnelles I/O nutzen

BeeOND erstellt ein schnelles lokales Dateisystem aus den SSDs der Compute-Nodes:

```bash
#SBATCH --constraint=BEEOND
#SBATCH --exclusive
```

DELFIN nutzt automatisch `$BEEOND_MOUNTPOINT` wenn verfügbar.

### 5.2 Ressourcen optimal nutzen

```bash
# In CONTROL.txt:
orca_parallel_strategy=auto     # DELFIN verteilt Cores optimal
parallel_workflows=yes          # Ox/Red parallel (bei 40 Cores)
enable_auto_recovery=yes        # Automatische Fehlerbehandlung
```

### 5.3 Lange Rechnungen

Für schwierige Systeme, die >48h brauchen:

```bash
# In CONTROL.txt:
enable_job_timeouts=no

# Im Submit-Skript:
#SBATCH --time=120:00:00  # 5 Tage
```

### 5.4 Mehrere Systeme parallel (Job Array)

Für High-Throughput Screening:

```bash
# Struktur:
# system_1/CONTROL.txt + input.txt
# system_2/CONTROL.txt + input.txt
# ...

# Job Array Submit-Skript anpassen und:
sbatch --array=1-10%3 array_script.sh
```

## Schritt 6: Troubleshooting

### Module nicht gefunden

```bash
# Alle verfügbaren Module anzeigen
module avail

# Modul suchen
module avail orca
module avail python
```

### ORCA nicht im PATH

```bash
# Nach dem Laden des Moduls prüfen
module load chem/orca/6.1.1
which orca

# Sollte zeigen:
# /opt/bwhpc/common/chem/orca/6.1.1/orca
```

### Job hängt in Queue

```bash
# Priorität und Queue-Status prüfen
squeue -u $USER -o "%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R"

# Partition-Info
sinfo -p cpu

# Job-Details
scontrol show job JOBID
```

### Zu wenig Speicher

```bash
# Im Submit-Skript erhöhen:
#SBATCH --mem=240G

# Oder High-Memory nutzen:
#SBATCH --partition=highmem
#SBATCH --mem=500G
```

### Timeout während SCF

```bash
# In CONTROL.txt:
enable_job_timeouts=no          # Deaktiviert Timeouts
enable_auto_recovery=yes        # Hilft bei SCF-Problemen
```

## Schritt 7: Beispiel-Workflows

### Schneller Test (dev_cpu)

```bash
# test.sh
#SBATCH --partition=dev_cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:20:00

module load chem/orca/6.1.1 devel/python/3.11
source ~/delfin_env/bin/activate
delfin
```

### Standard Redox-Berechnung

```bash
# Nutze: bwunicluster3_delfin.sh
sbatch bwunicluster3_delfin.sh
```

### Großes System mit High-Memory

```bash
# Nutze: bwunicluster3_delfin_highmem.sh
sbatch bwunicluster3_delfin_highmem.sh
```

## Kontakt und Support

- **BwUniCluster Wiki**: https://wiki.bwhpc.de/e/BwUniCluster3.0
- **DELFIN GitHub**: https://github.com/ComPlat/DELFIN
- **Cluster Support**: support@bwhpc.de

## Zusammenfassung

```bash
# 1. DELFIN installieren
module load devel/python/3.11
python -m venv ~/delfin_env
source ~/delfin_env/bin/activate
pip install delfin-complat

# 2. Projekt einrichten
mkdir mein_projekt && cd mein_projekt
delfin --define=struktur.xyz
# CONTROL.txt anpassen

# 3. Submit-Skript kopieren und anpassen
cp ~/ComPlat/DELFIN/examples/example_Job_Submission_Scripts/bwunicluster3_delfin.sh .
nano bwunicluster3_delfin.sh  # Module-Namen prüfen

# 4. Job starten
sbatch bwunicluster3_delfin.sh

# 5. Überwachen
squeue -u $USER
tail -f delfin_*.out
```

Viel Erfolg mit deinen Berechnungen auf dem BwUniCluster 3.0! 🚀
