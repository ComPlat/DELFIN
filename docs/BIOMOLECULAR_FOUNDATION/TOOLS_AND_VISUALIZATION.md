# Werkzeuge, Processing und Visualisierung

Der aktuelle Katalog enthält 48 Verträge und 30 Ausführungsbindungen. Die Zahlen
beschreiben Software-Schnittstellen, nicht 30 wissenschaftlich validierte Engines.
Die Laufzeitabfrage des Katalogs ist maßgeblich; Grenzen einzelner Adapter sind
in der [Workflow-Referenz](../BIOMOLECULAR_WORKFLOWS.md) beschrieben.

## Ausführbare Bausteine und ihre Grenzen

| Bereich | Implementierter Weg | Wesentliche Grenze |
| --- | --- | --- |
| Liganden | RDKit-Vorbereitung, Tautomere, explizite Zustands-/Datensatzwahl, Deskriptoren, Scaffolds, 2D-SVG | Keine vollständige pH-Mikrozustandsverteilung oder ADMET-Validierung |
| Strukturen | PDB/mmCIF lesen, Modell/Kette/Altloc/Wasser auswählen, Mapping | Selektion ist keine chemische Reparatur |
| Rezeptor | PDBFixer/OpenMM-Vorbereitung und Meeko-Konversion | Eingeschränkte kanonische PDB-Teilbereiche; keine allgemeine Metall-/Cofaktorchemie |
| Docking | Vina, GNINA und smina; Meeko-Ligand-/Rezeptor-/Pose-Konversion | Versions-, Box-, Seed- und Ressourcenangaben erforderlich; keine automatische Affinitätsinterpretation |
| Qualitäts-/Kontaktanalyse | PoseBusters, PLIP, ProLIF, Pose-RMSD | Fehlerhafte Posen bleiben sichtbar; ProLIF benötigt passende Graphen/H-Atome und Mapping |
| MedChem-Daten | Exakter Assay-Join mit Zuständen, Einheiten, Bedingungen, Zensierung, Replikaten | Kein stilles Vermischen von IC50, Ki, Docking-Score oder inkompatiblen Einheiten |
| Vorhersage | ESMFold-v1-Adapter; Import, Konfidenz und Modellauswahl | Einzelprotein und lokale Gewichte; keine allgemeine Cofolding-Unterstützung |
| Systemaufbau | ff14SB-Protein mit optional einem SMIRNOFF-Liganden in Vakuum oder TIP3P/Na+/K+/Cl- mit PME | Explizite SDF-Ladungen und lokales OFFXML; orthorhombische Box |
| Dynamik | Begrenzte Minimierung und LangevinMiddle-Dynamik, Checkpoint-Fortsetzung | CPU/Reference, feste Box, kein NPT; kein wissenschaftlicher Äquilibrierungsnachweis |
| Trajektorien | MDAnalysis-basierte Ausrichtung, RMSD/RMSF und Trägheitsradius | Explizite Atomauswahl/Gewichtung; periodisch nur bestätigtes `as_stored`, kein Unwrapping |

## Visualisierung als Teil des Fundaments

Der Offline-3Dmol-Export enthält Struktur-/Poseansicht, Ketten-/Wassersteuerung
und Bindungsbox. JavaScript und Eingabeartefakte liegen beim Export; ein Webserver
ist für das dokumentierte lokale Beispiel nicht erforderlich. Zusätzliche PyMOL-
und ChimeraX-Einstiegsskripte ermöglichen Desktop-Weiterarbeit. Diese Skripte sind
kein Installer für die Desktop-Programme und kein Nachweis, dass deren Rendering
auf dem Zielrechner funktioniert.

Die 2D-Ligandendarstellung verwendet ein SVG-Raster. Die Protein-Chatansicht
unterstützt unter anderem Cartoon-Darstellung, Heterokomponenten, Kettenwahl,
Residuenmarkierung und Vollbild. Konfidenzinterpretation bleibt explizit: B-Faktoren
sind ohne entsprechende Angabe nicht automatisch pLDDT.

PoseBusters-Redocking verlangt genau einen Referenzliganden. Der Bericht muss
jede Poseposition genau einmal aus der tatsächlich übergebenen Pose-Datei
enthalten. Eine passende Zeilenzahl allein genügt nicht; doppelte Positionen,
ungültige Indizes und Berichte über andere Dateien werden abgelehnt. Negative
Prüfergebnisse bleiben als solche erhalten und sind kein Ausführungsfehler.

Bereits vorhanden sind Export- und Steuerlogiktests sowie Syntaxprüfungen für
Desktop-Skripte. Echte Browser-/WebGL- und Desktop-Abnahme bleibt erforderlich.
Trajektorienanimation, periodische Rekonstruktion, konsistente Kontakt-Overlays,
Ensemblevergleiche und vollständige vergleichende Berichte sind noch nicht als
fertige allgemeine Integrationen ausgewiesen.

Für die Abnahme muss jede Ansicht die angezeigten Quelldateien und deren
Identität nachvollziehbar machen. Fehlende Scores oder Konfidenz dürfen nicht als
Nullwerte erscheinen. Pose-, Modell-, Ketten- und Zustandswechsel dürfen Daten
nicht verwechseln. Ein Render-Erfolg belegt keine wissenschaftliche Güte.

## Katalogumfang ohne fertigen Ausführungspfad

AlphaFold 2/Multimer/3, ColabFold, Boltz-2, Chai-1 und RoseTTAFold All-Atom haben
Planungs-/Eingabeabdeckung, aber noch keinen vollständigen ausführbaren Pfad in
dieser Schicht. Native Eingabeexporte starten keine Modelle. Dasselbe gilt je nach
Backend für Molscrub, DiffDock, HADDOCK3, Glide, GOLD, GROMACS, Amber, NAMD,
RFdiffusion, ProteinMPNN und die Bio-QM/MM-Regionskonstruktion.

Diese Bereiche brauchen jeweils eigene, versionsgebundene Adapter und
Abnahmekriterien. Eine generische Shell-Kommando-Vorlage ersetzt keine Integration.
Die bestehende ORCA-QM/MM-Funktionalität ersetzt nicht die noch fehlende Auswahl
von Bio-Region, Grenzatomen, Caps, Ladung, Multiplizität und Embedding.

Details zum [Ligandenparameter-Vertrag](LIGAND_PARAMETERIZATION.md).
