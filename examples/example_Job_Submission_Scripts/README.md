# DELFIN auf BwUniCluster 3.0 - Komplette Anleitung

Diese Sammlung enthält alle benötigten Skripte für DELFIN auf dem **BwUniCluster 3.0**.

## 📁 Inhalt

- **`setup_delfin.sh`** - Automatische DELFIN-Installation
- **`submit_standard.sh`** - Standard Job (40 CPUs, 120GB, 48h)
- **`submit_highmem.sh`** - High-Memory Job (40 CPUs, 500GB, 72h)
- **`submit_test.sh`** - Test Job (8 CPUs, 16GB, 30min)
- **`CONTROL_template.txt`** - CONTROL.txt Template mit optimierten Einstellungen
- **`README.md`** - Diese Anleitung

## 🚀 Schnellstart

### 1. Setup durchführen (einmalig)

```bash
# Auf BwUniCluster einloggen
ssh username@uc3.scc.kit.edu

# Skripte holen (wenn noch nicht vorhanden)
cd ~
mkdir -p UniBW
cd UniBW
# (Skripte von lokalem System hochladen oder aus DELFIN-Repo kopieren)

# Setup ausführen
bash setup_delfin.sh
```

**WICHTIG:** Vor dem Setup die Module-Namen prüfen:
```bash
module avail python    # z.B. devel/python/3.11
module avail orca      # z.B. chem/orca/6.1.1
```

Dann in `setup_delfin.sh` anpassen (Zeile 18-19).

### 2. Projekt erstellen

```bash
# Neues Projekt anlegen
mkdir -p ~/delfin_projects/mein_molekuel
cd ~/delfin_projects/mein_molekuel

# Module laden und DELFIN aktivieren
delfin-modules
delfin-activate

# CONTROL.txt und input.txt erstellen
delfin --define=struktur.xyz

# Oder CONTROL.txt Template verwenden
cp ~/UniBW/CONTROL_template.txt CONTROL.txt
nano CONTROL.txt  # Anpassen (charge, mult, oxidation_steps, etc.)
```

### 3. Job einreichen

```bash
# Submit-Skript kopieren
cp ~/UniBW/submit_standard.sh .

# Module-Namen im Skript prüfen/anpassen
nano submit_standard.sh  # Zeile 25-26

# Job starten
sbatch submit_standard.sh
```

### 4. Job überwachen

```bash
# Job-Status
squeue -u $USER

# Output live anschauen
tail -f delfin_*.out

# Job abbrechen
scancel <job_id>
```

## 📝 Submit-Skripte im Detail

### submit_standard.sh
**Für die meisten Rechnungen**
- 40 CPUs, 120 GB RAM
- 48 Stunden max. Laufzeit
- BeeOND für schnelles I/O
- Partition: `cpu`

**Verwendung:**
```bash
cp ~/UniBW/submit_standard.sh .
sbatch submit_standard.sh
```

### submit_highmem.sh
**Für große Systeme**
- 40 CPUs, 500 GB RAM (bis 1 TB möglich)
- 72 Stunden max. Laufzeit
- Für große Moleküle (>200 Atome)
- Korrelierte Methoden (CCSD(T), NEVPT2)
- Partition: `highmem`

**Verwendung:**
```bash
cp ~/UniBW/submit_highmem.sh .
sbatch submit_highmem.sh
```

### submit_test.sh
**Für schnelle Tests**
- 8 CPUs, 16 GB RAM
- 30 Minuten max. Laufzeit
- Development Queue (schneller Durchsatz)
- Partition: `dev_cpu`

**Verwendung:**
```bash
cp ~/UniBW/submit_test.sh .
sbatch submit_test.sh
```

## ⚙️ CONTROL.txt Konfiguration

### Wichtigste Einstellungen

```ini
# Basis
charge=0                    # Molekül-Ladung
mult=1                      # Multiplizität
input_file=input.txt        # Geometrie-Datei

# Workflow
method=OCCUPIER
OCCUPIER_method=auto
oxidation_steps=1,2,3
reduction_steps=1,2,3

# Cluster-Optimierungen
pal_jobs=40                 # Wird von SLURM überschrieben
parallel_workflows=yes      # Bei ≥40 Cores empfohlen
enable_auto_recovery=yes    # Automatische Fehlerbehandlung
enable_job_timeouts=no      # Für lange Rechnungen

# XTB Pre-Optimierung
XTB_OPT=yes                 # Empfohlen für schnellere Konvergenz
```

### Template verwenden

```bash
cp ~/UniBW/CONTROL_template.txt CONTROL.txt
nano CONTROL.txt  # Anpassen
```

## 🔧 Häufige Anpassungen

### Module-Namen ändern

In allen `submit_*.sh` Skripten:

```bash
# Diese Zeilen anpassen:
module load chem/orca/6.1.1        # Deine ORCA-Version
module load devel/python/3.11      # Deine Python-Version
```

Module-Namen herausfinden:
```bash
module avail orca
module avail python
```

### Ressourcen anpassen

Im Submit-Skript die SBATCH-Zeilen ändern:

```bash
#SBATCH --cpus-per-task=40     # CPU-Anzahl
#SBATCH --mem=120G             # Speicher
#SBATCH --time=48:00:00        # Max. Laufzeit (HH:MM:SS)
#SBATCH --partition=cpu        # Partition
```

### Längere Laufzeit

```bash
#SBATCH --time=120:00:00       # 5 Tage

# In CONTROL.txt:
enable_job_timeouts=no         # Kein Software-Timeout
```

## 📊 Typische Workflows

### Standard Redox-Berechnung

```bash
# 1. Projekt erstellen
mkdir ~/delfin_projects/fe_complex && cd ~/delfin_projects/fe_complex

# 2. Setup
delfin-modules && delfin-activate
delfin --define=fe_complex.xyz

# 3. CONTROL.txt anpassen
nano CONTROL.txt
# charge=2, mult=1, oxidation_steps=1,2, reduction_steps=1,2

# 4. Job starten
cp ~/UniBW/submit_standard.sh .
sbatch submit_standard.sh

# 5. Überwachen
squeue -u $USER
tail -f delfin_*.out
```

### Großes System (>200 Atome)

```bash
# High-Memory Job nutzen
cp ~/UniBW/submit_highmem.sh .
nano submit_highmem.sh  # Ggf. --mem=800G für sehr große Systeme
sbatch submit_highmem.sh
```

### Quick Test

```bash
# Test-Job für schnelles Feedback
cp ~/UniBW/submit_test.sh .
nano CONTROL.txt  # oxidation_steps=1, reduction_steps=1
sbatch submit_test.sh
```

## 🐛 Troubleshooting

### Module nicht gefunden

```bash
# Verfügbare Module anzeigen
module avail

# Modul laden testen
module load chem/orca/6.1.1
which orca
```

### Job hängt in Queue

```bash
# Queue-Status prüfen
sinfo -p cpu

# Job-Details
scontrol show job <job_id>

# Eventuell andere Partition probieren
#SBATCH --partition=dev_cpu  # Für Tests
```

### ORCA Fehler

```bash
# In CONTROL.txt:
enable_auto_recovery=yes        # Automatische Fehlerbehandlung
max_recovery_attempts=2         # Mehr Versuche

# Oder weniger aggressive Methode:
method_name=PBE                 # Statt B3LYP
basis_set=def2-SVP              # Kleineres Basis-Set
```

### Zu wenig Speicher

```bash
# High-Memory Job nutzen
cp ~/UniBW/submit_highmem.sh .
sbatch submit_highmem.sh

# Oder in CONTROL.txt:
basis_set=def2-SVP              # Statt def2-TZVP
orca_parallel_strategy=threads  # Weniger Speicher als MPI
```

## 📚 Hilfreiche Kommandos

### SLURM

```bash
# Alle eigenen Jobs
squeue -u $USER

# Job-Details
scontrol show job <job_id>

# Job abbrechen
scancel <job_id>

# Alle eigenen Jobs abbrechen
scancel -u $USER

# Queue-Info
sinfo -p cpu
```

### DELFIN

```bash
# Module + Environment
delfin-modules && delfin-activate

# Projekt erstellen
delfin --define=struktur.xyz

# Hilfe
delfin --help

# Version
delfin --version

# Update
cd ~/DELFIN && git pull
```

### Module

```bash
# Verfügbare Module
module avail

# Geladene Module
module list

# Modul laden
module load chem/orca/6.1.1

# Alle entladen
module purge
```

## 🔄 Updates

### DELFIN updaten

```bash
cd ~/DELFIN
git pull
source .venv/bin/activate
pip install -e .
```

### Skripte updaten

```bash
cd ~/UniBW
# Neue Versionen von GitHub oder lokalem System holen
```

## 💡 Best Practices

1. **Immer zuerst testen**
   - Nutze `submit_test.sh` mit `oxidation_steps=1, reduction_steps=1`
   - Prüfe ob Setup funktioniert

2. **XTB Pre-Optimierung nutzen**
   - `XTB_OPT=yes` in CONTROL.txt
   - Spart Zeit und verbessert Konvergenz

3. **Auto-Recovery aktivieren**
   - `enable_auto_recovery=yes`
   - Vermeidet manuelle Interventionen bei SCF-Problemen

4. **BeeOND nutzen**
   - `#SBATCH --constraint=BEEOND`
   - Deutlich schnelleres I/O

5. **Ressourcen optimal nutzen**
   - Bei ≥40 Cores: `parallel_workflows=yes`
   - Bei <40 Cores: `parallel_workflows=no`

## 📞 Support

- **BwUniCluster Wiki**: https://wiki.bwhpc.de/e/BwUniCluster3.0
- **DELFIN GitHub**: https://github.com/ComPlat/DELFIN
- **Cluster Support**: support@bwhpc.de

## ✅ Checkliste

- [ ] `setup_delfin.sh` ausgeführt
- [ ] Module-Namen in Submit-Skripten angepasst
- [ ] Test-Job erfolgreich durchgelaufen
- [ ] CONTROL.txt für Projekt konfiguriert
- [ ] Submit-Skript ins Projektverzeichnis kopiert
- [ ] Job mit `sbatch` eingereicht
- [ ] Output mit `tail -f` überwacht

---

**Viel Erfolg mit deinen DELFIN-Berechnungen! 🚀**
