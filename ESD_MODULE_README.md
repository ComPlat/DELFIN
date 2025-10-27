# ESD Module für DELFIN

## Übersicht

Das ESD (Excited State Dynamics) Modul ermöglicht die Berechnung von elektronischen Zuständen (S0, S1, T1, T2) und deren Übergängen (ISCs und ICs) in einem separaten `ESD/` Verzeichnis.

## Implementierte Dateien

1. **delfin/esd_input_generator.py** - Generiert ORCA Input-Dateien für:
   - Elektronische Zustände (S0, S1, T1, T2)
   - Intersystem Crossings (ISCs)
   - Internal Conversions (ICs)

2. **delfin/esd_module.py** - Hauptmodul mit Workflow-Logik:
   - Orchestriert alle ESD-Berechnungen
   - Verwaltet Dependencies zwischen Jobs
   - Parallelisiert Berechnungen über GlobalJobManager

3. **delfin/pipeline.py** - Integration in DELFIN Pipeline:
   - `run_esd_phase()` Funktion hinzugefügt
   - Kann nach classic/manually/OCCUPIER laufen

4. **delfin/cli.py** - CLI Integration:
   - ESD wird automatisch gestartet wenn `ESD_modul=yes`

## Konfiguration im CONTROL File

### Automatisches Template erstellen

Das ESD-Modul ist vollständig in `delfin --define` integriert:

```bash
delfin --define
```

Dies erstellt ein CONTROL.txt mit allen ESD-Feldern vorausgefüllt.

### Manuelle Konfiguration

Oder füge folgende Zeilen zu deinem CONTROL.txt hinzu:

```ini
------------------------------------
ESD Module:
ESD_modul=yes
states=S0,S1,T1,T2
ISCs=S1>T1,T1>S1,S1>T2,T2>S1,S0>T1,T1>S0
ICs=S1>S0,S0>S1,T1>T2,T2>T1
------------------------------------
deltaSCF:
deltaSCF_DOMOM=false
deltaSCF_PMOM=true
deltaSCF_keepinitialref=true
deltaSCF_SOSCFHESSUP=LBFGS
------------------------------------
```

### Parameter-Erklärung

- **ESD_modul**: `yes` oder `no` - Aktiviert/deaktiviert das ESD-Modul
- **states**: Komma-getrennte Liste der zu berechnenden Zustände (S0, S1, T1, T2)
- **ISCs**: Komma-getrennte Liste der Intersystem Crossings (Format: `InitialState>FinalState`)
- **ICs**: Komma-getrennte Liste der Internal Conversions (Format: `InitialState>FinalState`)

### deltaSCF Einstellungen (optional)

Diese werden für S1 und T2 Berechnungen verwendet:

- **deltaSCF_DOMOM**: `true`/`false` - DOMOM Methode aktivieren
- **deltaSCF_PMOM**: `true`/`false` - PMOM Methode aktivieren
- **deltaSCF_keepinitialref**: `true`/`false` - Initiale Referenz behalten
- **deltaSCF_SOSCFHESSUP**: SCF Hessian Update Methode (z.B. `LBFGS`)

## Verwendung

### Standalone (nur ESD, kein Redox)

```bash
# CONTROL.txt mit ESD_modul=yes und method nicht gesetzt oder leer
delfin
```

In diesem Fall muss `initial.inp` bereits existieren oder wird vom ESD-Modul erstellt.

### Nach Classic/Manually/OCCUPIER

```bash
# CONTROL.txt:
# method=classic  (oder manually oder OCCUPIER)
# ESD_modul=yes
delfin
```

Das ESD-Modul läuft automatisch nach der gewählten Methode.

## Workflow und Dependencies

### State Berechnungen

```
S0 (initial.inp → S0.out, S0.xyz, S0.gbw)
├─> S1 (benötigt S0.gbw für MOREAD)
├─> T1 (benötigt S0.gbw für MOREAD)
    └─> T2 (benötigt T1.gbw für MOREAD)
```

**Eigenschaften:**
- Alle States haben Multiplizität M=1
- Ladung wird aus CONTROL übernommen
- S0: UKS, OPT, FREQ
- S1: UKS, deltaSCF, OPT, FREQ (mit alphaconf 0,1 / betaconf 0)
- T1: UKS, OPT, FREQ
- T2: UKS, deltaSCF, OPT, FREQ (mit alphaconf 0,1 / betaconf 0)

### ISC/IC Berechnungen

ISCs und ICs benötigen beide beteiligten Zustände:

```
S1>T1 ISC: benötigt S1.out, S1.hess, T1.out, T1.hess
T1>S1 RISC: benötigt T1.out, T1.hess, S1.out, S1.hess
```

**DELE Berechnung:**
- Wird automatisch aus den final single point energies berechnet
- Konvertierung: Hartree → cm⁻¹ (Faktor: 219474.63)
- Wenn DELE nicht berechnet werden kann: Fallback auf 1692 cm⁻¹

**Temperatur:**
- Wird aus `temperature` Parameter im CONTROL übernommen (Standard: 298.15 K)

## Ordnerstruktur

```
Arbeitsverzeichnis/
├── CONTROL.txt
├── initial.inp (falls vorhanden)
├── initial.out
├── initial.xyz
├── initial.gbw
├── ESD/                    # Separater ESD-Ordner
    ├── S0.inp
    ├── S0.out
    ├── S0.xyz
    ├── S0.gbw
    ├── S0.hess
    ├── S1.inp
    ├── S1.out
    ├── S1.xyz
    ├── S1.gbw
    ├── S1.hess
    ├── T1.inp
    ├── T1.out
    ├── T1.xyz
    ├── T1.gbw
    ├── T1.hess
    ├── T2.inp
    ├── T2.out
    ├── T2.xyz
    ├── T2.gbw
    ├── T2.hess
    ├── S1_T1_ISC.inp
    ├── S1_T1_ISC.out
    ├── T1_S1_IC.inp
    ├── T1_S1_IC.out
    └── ... (weitere ISC/IC Jobs)
```

## Parallelisierung

Das ESD-Modul nutzt den globalen DynamicCorePool von DELFIN:

- **States** werden parallel berechnet wenn Dependencies erfüllt sind
- **ISCs/ICs** laufen parallel sobald beide States fertig sind
- Nutzt `parallel_workflows` Setting aus CONTROL
- Core-Zuteilung erfolgt automatisch via Workflow Manager

Beispiel Parallelität:
```
Zeit →
S0 ████████████████
   ├─> S1 ████████████████
   └─> T1 ████████████████
           └─> T2 ████████████████

Nach S1 und T1:
S1>T1 ISC ██████████
T1>S1 IC  ██████████
S1>S0 IC  ██████████
(alle parallel)
```

## Beispiel CONTROL.txt

```ini
input_file=input.txt
NAME=Benzophenone_ESD
charge=0
------------------------------------
Solvation:
implicit_solvation_model=CPCM
solvent=DMF
------------------------------------
Redox steps:
calc_initial=yes
method=classic
------------------------------------
Level of Theory:
functional=PBE0
disp_corr=D4
ri_jkx=RIJCOSX
aux_jk=def2/J
main_basisset=def2-SVP
temperature=298.15
maxiter=125
------------------------------------
deltaSCF:
deltaSCF_DOMOM=false
deltaSCF_PMOM=true
deltaSCF_keepinitialref=true
deltaSCF_SOSCFHESSUP=LBFGS
------------------------------------
ESD_modul=yes
states=S0,S1,T1,T2
ISCs=S1>T1,T1>S1
ICs=S1>S0
------------------------------------
Resource Settings:
PAL=32
maxcore=25000
parallel_workflows=yes
------------------------------------
```

## Fehlerbehebung

### S0.out existiert nicht
- Stelle sicher, dass `calc_initial=yes` gesetzt ist
- Oder kopiere `initial.out` manuell nach `ESD/S0.out`

### DELE kann nicht berechnet werden
- Prüfe ob beide State-Berechnungen erfolgreich waren
- Schaue in die `.out` Dateien nach "FINAL SINGLE POINT ENERGY"

### ISC/IC schlagen fehl
- Prüfe ob `.hess` Dateien existieren
- Frequenzberechnungen müssen für beide States erfolgreich sein

### Zu wenig Parallelität
- Erhöhe `PAL` im CONTROL
- Setze `parallel_workflows=yes`
- Prüfe `pal_jobs` Setting

## Technische Details

### ORCA Input-Generierung

**S0 Input:**
```
! PBE0 UKS def2-SVP D4 RIJCOSX def2/J CPCM(DMF) OPT FREQ
%base "S0"
%pal nprocs 32 end
%maxcore 25000
* xyzfile 0 1 initial.xyz
```

**S1 Input (mit deltaSCF):**
```
! PBE0 UKS def2-SVP D4 RIJCOSX def2/J CPCM(DMF) deltaSCF OPT FREQ NODIIS MOREAD
%base "S1"
%moinp "S0.gbw"
%scf
  DOMOM false
  pmom true
  keepinitialref true
  alphaconf 0,1
  betaconf 0
  SOSCFHESSUP LBFGS
end
%pal nprocs 32 end
%maxcore 25000
* xyzfile 0 1 S0.xyz
```

**ISC Input (S1>T1):**
```
! PBE0 RKS def2-SVP D4 RIJCOSX def2/J CPCM(DMF) ESD(ISC)
%base "S1_T1_ISC"
%TDDFT
  NROOTS  5
  SROOT   1
  TROOT   1
  TROOTSSL 0
  DOSOC   TRUE
END
%ESD
  ISCISHESS       "S1.hess"
  ISCFSHESS       "T1.hess"
  USEJ            TRUE
  DOHT            TRUE
  TEMP            298.15
  DELE            XXXX    # automatisch berechnet aus S1/T1 Energien
END
%pal nprocs 32 end
%maxcore 25000
* xyzfile 0 1 S1.xyz
```

## Erweiterungsmöglichkeiten

Das ESD-Modul kann leicht erweitert werden für:

- Weitere elektronische Zustände (S2, S3, T3, ...)
- Spin-Orbit Coupling (SOC) Berechnungen
- Vibronische Kopplungen
- Rate-Konstanten Berechnungen
- Spektren-Simulationen

## Lizenz und Autor

Teil von DELFIN - Implementiert 2025
