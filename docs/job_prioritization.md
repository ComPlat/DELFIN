# Job-Prioritisierung und Ressourcen-Boost in DELFIN

## Konzept

DELFIN erkennt automatisch **Bottleneck-Jobs** (Engstellen im Workflow) und optimiert sie zweifach:

1. **Priority-Boost**: Höhere Priorität → werden vorgelassen
2. **Resource-Boost**: Exklusive Bottlenecks bekommen alle verfügbaren Cores

## Wie es funktioniert

### 1. Bottleneck-Erkennung

Vor dem Start eines Workflows analysiert DELFIN alle Jobs:

```python
for job in all_jobs:
    downstream_count = count_downstream_jobs(job)

    if downstream_count >= 3:  # Bottleneck!
        job.priority = HIGH
```

**Bottleneck-Definition:** Ein Job ist ein Bottleneck wenn **>= 3 andere Jobs** (direkt oder indirekt) davon abhängen.

### 2. Prioritäts-Boost

Jobs mit hoher Priorität werden in der Queue **vorgelassen**:

```
Queue ohne Priorität:
[Job #5 (normal)] [Job #6 (normal)] [Bottleneck (normal)] [Job #7 (normal)]
→ Bottleneck wartet bis #5, #6 fertig sind

Queue mit Priorität:
[Bottleneck (HIGH)] [Job #5 (normal)] [Job #6 (normal)] [Job #7 (normal)]
→ Bottleneck startet sofort!
```

### 3. Normale Ausführung

Danach läuft alles normal:
- Jobs werden nach Priorität aus der Queue geholt
- Dependencies werden weiterhin beachtet
- Core-Allokation bleibt unverändert

## Resource-Boost für exklusive Bottlenecks

### Was ist ein exklusiver Bottleneck?

Ein Job ist ein **exklusiver Bottleneck** wenn ALLE anderen pending Jobs von ihm abhängen:

```python
pending_jobs = {job_A, job_B, job_C, bottleneck}

# Prüfe: Hängen A, B, C alle von bottleneck ab?
if {job_A, job_B, job_C}.issubset(downstream_of_bottleneck):
    # JA → exklusiver Bottleneck!
```

### Automatischer Resource-Boost mit dynamischer Erkennung

**Kontinuierliche Prüfung während der Laufzeit:**

Das System prüft in jedem Scheduling-Zyklus:
1. Ist ein exklusiver Bottleneck ready?
2. Laufen gerade andere Jobs?
3. **Wenn ja:** Warte bis andere fertig sind, dann starte Bottleneck mit max cores

**Warum warten?**
- ORCA kann cores zur Laufzeit nicht ändern
- Besser: Bottleneck wartet 2 Min → startet mit 64 cores → spart 10 Min
- Als: Bottleneck startet sofort mit 32 cores

**Beispiel-Timeline:**

```
T=0:  Job #1 (32c), Job #2 (32c) laufen
T=2:  occ_proc_initial wird ready
      Prüfung: Ist exklusiver Bottleneck? JA!
      Prüfung: Laufen andere Jobs? JA!
      → WARTE (starte nicht parallel)

[INFO] Waiting for running jobs to finish before starting exclusive bottleneck occ_proc_initial

T=5:  Job #1, #2 fertig
      Prüfung: Laufen andere Jobs? NEIN!
      → Starte occ_proc_initial mit 64 cores

[INFO] Job occ_proc_initial is exclusive bottleneck → allocating max cores (64)

T=5-15: occ_proc_initial läuft mit 64 cores (10 Min)

Gesamt: 15 Min (statt 20 Min wenn parallel mit 32 cores)
→ 25% schneller!
```

**Ohne dynamische Erkennung:**
```
T=0:  Job #1, #2 laufen (keine exklusiven Bottlenecks)
T=2:  occ_proc_initial wird ready
      → Startet mit 32 cores (parallel zu #1, #2)
T=2-22: occ_proc_initial läuft (20 Min)

Gesamt: 22 Min
```

## Beispiel

### Workflow:

```
occ_proc_initial
    ↓
    ├─→ occupier_initial (IMAG)
    ├─→ occ_proc_ox_1
    ├─→ occ_proc_red_1
    ├─→ occupier_absorption
    └─→ occupier_t1_state
```

**Analyse:**
- `occ_proc_initial` hat 5 downstream jobs → **HIGH priority**
- Andere Jobs: NORMAL priority

**Effekt:**
- Wenn mehrere Jobs ready sind, startet `occ_proc_initial` **zuerst**
- 5 nachfolgende Jobs werden früher unblocked
- **Gesamtdurchsatz steigt**

## Logging

```
[INFO] [occupier] Boosted priority for 1 bottleneck job(s)
[INFO] [priority] Job occ_proc_initial boosted to HIGH priority (5 downstream jobs depend on it)
```

## Konfiguration

**Bottleneck Threshold** (Standard: 3):
```python
# In job_priority.py
adjust_job_priorities(all_jobs, bottleneck_threshold=3)
```

Jobs mit >= 3 downstream Jobs bekommen HIGH priority.

## Vorteile

✅ **Einfach** - Nur 60 Zeilen Code
✅ **Automatisch** - Keine manuelle Konfiguration
✅ **Effektiv** - Reduziert Wartezeiten für kritische Jobs
✅ **Sicher** - Ändert nichts an Core-Allokation oder Dependencies
✅ **Transparent** - Alle Entscheidungen werden geloggt

## Technische Details

### Downstream-Zählung (BFS)

```python
def count_downstream_jobs(job_id, all_jobs):
    count = 0
    queue = [job_id]
    visited = set()

    while queue:
        current = queue.pop(0)
        visited.add(current)

        # Finde alle Jobs die von current abhängen
        for other_job in all_jobs:
            if current in other_job.dependencies:
                count += 1
                queue.append(other_job.job_id)

    return count
```

**Komplexität:** O(V + E) - V = Jobs, E = Dependencies

### Integration

Die Funktion wird automatisch in `_WorkflowManager.run()` aufgerufen:

```python
def run(self):
    # Adjust job priorities based on bottleneck detection
    self._adjust_priorities_for_bottlenecks()

    # Normal scheduling...
```

## Limitierungen

1. **Statisch**: Prioritäten werden nur **einmal** vor dem Start gesetzt, nicht dynamisch angepasst
2. **Heuristisch**: Threshold von 3 ist arbiträr (könnte adaptiv sein)
3. **Keine Kosten**: Berücksichtigt nicht die Laufzeit der Jobs (nur Anzahl Dependencies)

## Zukünftige Erweiterungen

- **Adaptive Thresholds**: Basierend auf Workflow-Größe
- **Gewichtete Prioritäten**: Jobs mit längerer erwarteter Laufzeit höher gewichten
- **Runtime-Anpassung**: Prioritäten während Laufzeit updaten
