# Ligandenparameter und Protein-Ligand-MD

## Implementierter Weg

`protein_ligand_md` verbindet vier Schritte: `bio_ligand_parameterize` →
`bio_openmm_complex_parameterize` → `bio_openmm_run` → `bio_mdanalysis`.
Die wissenschaftlichen Engines wurden hier nicht installiert oder ausgeführt.
Die Tests ersetzen ihre Prozessgrenzen; echte Installation und Referenzläufe
gehören zur Cluster-Abnahme.

Der Parameteradapter verwendet OpenMMForceFields 0.15.0, OpenFF Toolkit 0.17.0
und OpenMM 8.4.0. Das neue verwaltete Rezept `openmmforcefields` umfasst diese
Abhängigkeiten und RDKit. Discovery, bestehende Cluster-Module, isolierte
Installation und Reparatur verwenden dieselben Mechanismen wie die anderen
[Bio-Rezepte](INSTALLATION.md). Ein vorhandener Paketmanager ist erforderlich.

## Verbindliche Eingaben

- Ein einzelnes zusammenhängendes Molekül als 3D-SDF, einschließlich aller
  Wasserstoffe als Atome; festgelegte Stereochemie und formale Ladungen.
- Eine lokale SMIRNOFF-OFFXML-Datei mit expliziten vdW-/Elektrostatik-Handlern.
  Keine automatischen Downloads, Includes, virtuellen Zentren oder Plugins.
- Partialladungen in exakt derselben Atomreihenfolge wie im SDF. Ihre Summe
  muss innerhalb von 0.0001 Elementarladungen zur formalen Ligandenladung passen.
- Ein vorbereitetes kanonisches ff14SB-Protein mit Wasserstoffen und genau
  einem `LIG`-Residuum; dessen Atome heißen `L001`, `L002`, … in SDF-Reihenfolge.
  Dafür den vorhandenen Komplexassemblierungsadapter verwenden. Koordinaten
  müssen mit dem SDF übereinstimmen (PDB-Rundungstoleranz 0.0006 Å).

Die Ladungsdatei ist ein JSON-Objekt mit genau diesen Feldern:

| Feld | Inhalt |
| --- | --- |
| `schema_version` | `1` |
| `ligand_sha256` | SHA-256 der exakten SDF-Datei |
| `unit` | `elementary_charge` |
| `method` | Dokumentierte tatsächliche Ladungsmethode, einschließlich relevanter Einstellungen |
| `source` | Referenz auf den Quelllauf und seine Provenienz |
| `values` | Liste endlicher Zahlen, ein Wert pro SDF-Atom |

Es gibt keine Platzhalterladungen und keine automatische Wahl zwischen RESP,
AM1-BCC oder einer anderen Methode. Ladungen können aus einem externen oder
QM-Lauf übernommen werden, sofern dessen Atomzuordnung geprüft wurde. Ein
universeller QM-Ladungsexporter ist damit noch nicht implementiert.
Vollständig verschwindende Ladungen werden abgelehnt: der ausgewählte Generator
würde sie als fehlend behandeln. Implizite Ladungsberechnung, neue Konformere
und fraktionelle Bindungsordnungsberechnung werden im Worker ausdrücklich blockiert.

## Artefakte und Erhaltungsprüfungen

Das `ligand_parameters`-Paket enthält SDF, Ladungs-JSON, ursprüngliches OFFXML,
generiertes OpenMM-FFXML sowie ein Manifest mit Hashes, chemischer Identität
und Laufzeitprovenienz. Der Validator prüft Atomnamen/-reihenfolge, Elemente,
Bindungsgraph, vollständige Bindungs-/Winkelparameter und endliche Parameter.
Die gelieferten Ladungen müssen im generierten FFXML erhalten bleiben.

Der Systemadapter übernimmt das Paket in sein Systembundle. Er prüft die
Liganden-Atomzuordnung, einzelne Ladungen, Massen und Lennard-Jones-Parameter
bis zum serialisierten OpenMM-System. Er lehnt Protein-Ligand-Kovalenz ab.
Verschachtelte Dateien werden bei Weitergabe, Kopie und Resume hashgeprüft.
Dies belegt Datenkonsistenz, keine wissenschaftliche Eignung eines Kraftfelds.

Die Systemschemas 3 (Vakuum) und 4 (TIP3P/PME) erweitern die bestehenden
Protein-Schemas 1 und 2. Constraint-Politik bleibt ausdrücklich `none` oder
`HBonds`; der Generator entfernt eigene OFFXML-Constraints vor der Übergabe.
MD verwendet CPU/Reference, eine feste Box und die bestehenden geprüften
Checkpoint-/Trajektorienverträge. Periodische Analyse verlangt `as_stored`;
sie rekonstruiert keine über Boxgrenzen verteilten Moleküle.

## Plan für den Cluster erstellen

Der vorhandene Planbuilder erzeugt nur einen Plan, ohne Programme zu starten:

```sh
python -m examples.biomolecular_foundation.build_md_plan --help
```

Zu den expliziten Protein-/MD-/Analyseargumenten kommen gemeinsam hinzu:

```text
--structure /shared/prepared/complex.pdb
--ligand /shared/prepared/ligand.sdf
--charges-file /shared/prepared/charges.json
--ligand-forcefield /shared/forcefields/chosen.offxml
--openmmforcefields-version 0.15.0
--openff-version 0.17.0
```

Diese Zeilen sind zusätzliche Argumente, kein vollständiger ausführbarer Befehl.
`--charge` bezeichnet die erwartete Gesamtladung des fertigen Systems, nicht
nur die Ligandenladung. Für TIP3P/PME zusätzlich eine explizite Solvationsdatei
und `--periodic-policy as_stored` angeben. Vor neutralisierender Solvation muss
`solute_charge` in dieser Datei die Protein-plus-Ligand-Ladung enthalten.
Die CLI prüft Argumente und Planstruktur; Datei-/Chemieprüfung erfolgt beim Preflight
beziehungsweise beim zuständigen Ausführungsschritt.

## Grenzen und externe Abnahme

Unterstützt wird ein organischer, nichtkovalenter Ligand mit den Elementen
H, B, C, N, O, F, Si, P, S, Cl, Br und I, ohne Isotope oder Radikale.
Ob das gewählte Kraftfeld den konkreten Liganden abdeckt, entscheidet dessen
Parameterzuweisung. Mehrere Liganden, Metall-/Cofaktorchemie, kovalente Liganden,
GAFF, automatische Ladungserzeugung, NPT und andere MD-Engines benötigen weitere
Implementierungen.

Auf dem Cluster zuerst ein dokumentiertes Protein-Ligand-Referenzsystem wählen.
Installation/Module und echte Generierung prüfen; Atomzuordnung, Ladungen,
Massen, LJ- und gebundene Parameter sowie Energien gegen einen unabhängigen
Referenzaufbau vergleichen. Danach Vakuum und feste TIP3P-Box, kurze Dynamik,
Checkpoint-Fortsetzung und Analyse abnehmen. Auch veränderte SDF-/Ladungsdateien,
fehlende Parameter und unzulässige Chemie müssen kontrolliert scheitern.

API-Grundlage: [OpenMMForceFields 0.15.0](https://github.com/openmm/openmmforcefields/blob/0.15.0/openmmforcefields/generators/template_generators.py)
und [OpenFF Toolkit 0.17.0 Molecule](https://docs.openforcefield.org/projects/toolkit/en/0.17.0/api/generated/openff.toolkit.topology.Molecule.html).
