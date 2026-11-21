# Autoinstaller und Cluster-Laufzeit

## Gemeinsamer Mechanismus mit QM

Bio-Rezepte hängen an der zentralen Installer-Gruppe `bio` und verwenden
DELFINs Laufzeit- und Discovery-Infrastruktur. Eine passende vorhandene
Installation soll genutzt werden, bevor eine neue Umgebung entsteht.
Python-Pakete werden durch Import und exakte Distributionsversion geprüft;
Kommandozeilentools durch ihren Versionsaufruf. Eine reine Präsenzanzeige ist
schwächer als diese Laufzeitprüfung.

Cluster-Module werden im passenden Modulkontext geprüft und ausgeführt.
Dieser Kontext bleibt beim Worker erhalten. Modul-Python verwendet `-s -P`,
damit das Modul-PYTHONPATH nutzbar bleibt und User-Site beziehungsweise unsichere
Suchpfade ausgeschlossen werden; gewöhnliche separate Umgebungen verwenden
`-I`. Virtuelle Python-Umgebungen dürfen durch Symlink-Auflösung nicht versehentlich
auf den Basisinterpreter zurückfallen. Reale Modulkonfigurationen sind noch extern
abzunehmen.

## Vorhandene verwaltete Rezepte

| Tool | Standardversion | Prüfungen / zusätzliche Abhängigkeiten |
| --- | --- | --- |
| Vina | 1.2.5 | `vina --version` |
| Meeko | 0.7.1 | Import und Version; RDKit, gemmi, scipy; Vorbereitungs-/Exportkommandos |
| PDBFixer | 1.12 | Import und Version; OpenMM 8.4.0 |
| OpenMM | 8.4.0 | Import und Version |
| OpenMMForceFields | 0.15.0 | Generatoren, OpenFF Toolkit 0.17.0, OpenMM 8.4.0 und RDKit |
| MDAnalysis | 2.9.0 | MDAnalysis, Align-Modul und NumPy; Biopython |
| ProLIF | 2.2.3 | ProLIF, MDAnalysis und RDKit |
| PLIP | 3.0.1 | PLIP und Open Babel; PLIP-Kommando |
| PoseBusters | 0.5.1 | Paket und CLI-Modul, RDKit ≥2024.9, pandas, NumPy, PyYAML; `bust` |

Die maßgeblichen Rezepte stehen in `delfin/tools/bio/environments.py`.
Die angegebenen Pakete sind keine Behauptung, dass Solver und echte Engines hier
ausgeführt wurden. Weitere transitive Paketversionen werden beim Lösen bestimmt
und anschließend mit Build-/Hash-Provenienz aufgezeichnet.

## Provisionierung und Reparatur

Der Standardpfad ist `~/.local/share/delfin/bio_tools`; ein anderer Speicherort
kann mit `DELFIN_BIO_TOOLS_ROOT` gewählt werden. Rezept und Versionen bestimmen
einen eigenen Prefix unter `envs`. Die Installation verändert keine vorhandene
Cluster-Umgebung und nicht die DELFIN-Hauptumgebung.

Ein vorhandenes micromamba, mamba oder conda wird benötigt. Der Bio-Installer
bootstrapt keinen Paketmanager. Er nutzt conda-forge mit strikter Kanalpriorität
und Python 3.11. Ein Dateilock schützt die Prefix-Installation; ein Ownership-
Manifest dokumentiert die beabsichtigte Umgebung schon vor dem Installationslauf.
Beschädigte fremde oder nicht zuordenbare Prefixe werden nicht automatisch repariert.
Eine explizite Reparatur verlangt ein passendes Ownership-Manifest.

Der gemeinsame Schalter `DELFIN_AUTO_INSTALL_QM_TOOLS` gilt auch für diese
automatische Bio-Provisionierung. Read-only-Planung installiert niemals;
auch im Engine-Kindprozess gibt es keinen versteckten Installer. Ein gescheiterter
Versuch wird innerhalb der Sitzung nicht endlos wiederholt.

Die folgenden Kommandos sind für die spätere Zielumgebung dokumentiert; sie
wurden für diese Dokumentation nicht ausgeführt:

```sh
python -m delfin.installer --list
python -m delfin.installer --plan openmm,mdanalysis
python -m delfin.installer --install openmm mdanalysis
python -m delfin.installer --repair openmm
```

## Noch fehlende Installationsparität

Nicht jeder Adapter besitzt bereits ein verwaltetes Rezept. Insbesondere sind
GNINA, smina und ESMFold nicht durch die obige Rezeptliste abgedeckt.
ESMFold benötigt außerdem lokal bereitgestellte, hashgeprüfte Modellgewichte;
der Adapter lädt sie nicht selbst herunter. Für weitere Predictor-, Design- und
MD-Backends fehlen zusätzlich Ausführungsadapter oder vollständig geprüfte
Installationspfade.

Neue Rezepte benötigen Versions-/Importproben, vorhandene-Installation-zuerst,
separate Prefixe, nachvollziehbare Fehler, Reparaturtests und eine spätere echte
Installation auf dem Zielcluster. GPU-, Lizenz-, Modell- und Datenbankvoraussetzungen
müssen separat ausgewiesen werden. Ein Katalogeintrag darf diese Voraussetzungen
nicht verdecken.

PoseBusters 0.5.1 ist im [conda-forge-Paketbestand](https://anaconda.org/channels/conda-forge/packages/posebusters/files)
verfügbar. Abhängigkeiten und `bust`-Einstiegspunkt wurden anhand der
[versionierten Projektmetadaten](https://github.com/maabuu/posebusters/blob/v0.5.1/pyproject.toml)
abgeglichen. Das Rezept verwendet bewusst die zum bestehenden Adapterbeispiel
passende Version. Installation und wissenschaftlicher Lauf bleiben extern zu prüfen.

Auch ein explizit ausgewählter Interpreter und weitere Pakete eines Workers
müssen die Importprüfung bestehen. Vorhandene Versionsmetadaten allein reichen
nicht: eine fehlende native Bibliothek oder ein beschädigtes Python-Modul darf
nicht erst nach erfolgreichem Preflight als fehlgeschlagene Rechnung erscheinen.
