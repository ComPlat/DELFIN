# Review-Stand von PR #746

Dieser PR ist eine Implementierung des biomolekularen Grundgerüsts mit ausgewählten
ausführbaren Pfaden. Er beansprucht keine vollständige wissenschaftliche oder
funktionale QM-Parität. Die [weiteren Ausbauschritte](ACCEPTANCE.md) bleiben
sichtbar und werden nicht durch einen grünen Schnittstellentest erledigt.

## Abgegrenzter Inhalt

Reviewt werden die typisierten Pläne und Artefakte, 28 Ausführungsbindungen,
acht isolierte Installationsrezepte, vorhandene Cluster-Umgebungen, serialisierte
Ausführung und Wiederaufnahme sowie Processing, Viewer und die dokumentierten
eingeschränkten Engine-Adapter. Der OpenMM-Pfad ist Protein-only, ff14SB,
Vakuum oder feste orthorhombische TIP3P/PME-Box, CPU/Reference. Die periodische
Analyse rekonstruiert keine über Boxgrenzen geteilten Moleküle.

Liganden-Kraftfeldzuweisung, zusätzliche Predictor-/MD-/Design-Engineausführung,
automatische Modellbeschaffung, NPT, allgemeine Metall-/Cofaktorchemie,
Pro-Knoten-Scheduling und Bio-QM/MM-Regionsaufbau sind Folgearbeit. Diese Grenzen
müssen bei Planung, Ausführung und Ergebnisinterpretation erhalten bleiben.

## Vor dem Review geschlossene Blocker

- Vollständiges `main` bis `88cf9ed31a5b6042f4f8549eefa02bec59ce371e` integriert;
  keine Konflikte, auch die Binärartefakte sind vorhanden und hashgeprüft.
- PoseBusters-Provisionierung wird jetzt vom DAG-Ausführungspreflight erreicht;
  read-only-Preflight und fehlerhafte Pläne installieren weiterhin nichts.
- Die beiden reinen Bio-MCP-Planungswerkzeuge sind im Gate-Snapshot und in den
  Modusprüfungen berücksichtigt; Ausführungswerkzeuge erhalten dadurch keine
  allgemeine read-only-Freigabe.
- Der echte Vina-Smoke ist ein explizites Cluster-Beispiel statt eines vom
  Host abhängigen Tests. Installer-Testprozesse bekommen ein isoliertes HOME
  über den bestehenden `child_env`-Mechanismus.
- Fixture-Reexports passieren die CI-Lintprüfung. Die bisher bekannten
  Vertragsfehler bei ORCA-Thermochemie, QM/MM-Kategorie und Producer-Hinweisen
  sind korrigiert, ohne wissenschaftliche ORCA-Ausführung zu ändern.

## Prüfung und Umgebungsgrenzen

Die Prüfung umfasst alle durch den PR hinzugefügten oder geänderten Python-
Testdateien sowie angrenzende Installer-, Laufzeit-, Plattform- und CI-Gates.
Echte Engines, Paketinstallationen, Modell-Downloads und Slurm-Jobs werden
hier nicht ausgeführt. Die externe Abnahme bleibt in [ACCEPTANCE.md](ACCEPTANCE.md)
beschrieben.

Die vorhandene Python-Umgebung enthält alte `delfin-complat`-Metadaten mit
`numpy>=1.26`; der aktuelle Checkout deklariert `numpy>=1.26,<2.3`. Deshalb wurde
die PEP-517-Metadatenfunktion von setuptools für den aktuellen Checkout in einem
temporären Verzeichnis ausgeführt. Nur diese erzeugten Distributionsmetadaten
werden für die abschließende lokale Testsitzung zusätzlich auf den Python-Pfad
gelegt. Kein Paket wurde installiert oder aktualisiert, kein Test abgeschwächt.
Auf CI und Cluster muss DELFIN regulär aus dem zu prüfenden Checkout installiert
sein. Ohne diese passende Metadatenbasis scheitert lokal der vorhandene Test,
der das Weiterreichen der NumPy-Obergrenze an Familieninstaller prüft.

Der CI-Lintumfang mit ruff 0.16.3 (`F821,F811,F601,F402,F823`) ist fehlerfrei.
Die JavaScript-Viewer-Steuerprüfung und das Offline-Processing-Beispiel laufen
erfolgreich. Das belegt weder Browser/WebGL-Rendering noch wissenschaftliche
Qualität eines Docking- oder MD-Laufs.

Die abschließende Auswahl besteht mit **708 bestandenen und acht übersprungenen
Tests**, ohne Fehler. Sie umfasst 46 Testdateien; die Dateiliste steht in
[REVIEW_TESTS.txt](REVIEW_TESTS.txt). Zusätzlich liefen die 118 gezielten
CI-Regressionsprüfungen erfolgreich. Die vollständige GitHub-CI muss auf dem
veröffentlichten Merge-Commit erneut laufen; frühere grüne Einzeljobs werden
nicht als Prüfung dieses neuen Commits ausgegeben.

## Reihenfolge für Reviewer

1. Verträge, Bindings und Input-/Output-Prüfungen in `delfin/tools/bio/` lesen.
2. Installer/Discovery und gemeinsame Änderungen an `qm_runtime`, Installer,
   Step-Runner und Slurm-CLI gegen bestehende QM-Pfade prüfen.
3. Je Engine den dokumentierten Eingabeteilbereich und die zugehörigen
   Mock-/Negativtests vergleichen. Besonders Chemieerhalt, Einheiten,
   verschachtelte Hashes, Box und Checkpoint-Protokoll prüfen.
4. Processing-/Viewer-Beispiel öffnen und die Abnahmematrix auf einem
   geeigneten Cluster mit echten Versionen und Referenzsystemen durchführen.

Ein Merge setzt erfolgreiche verpflichtende GitHub-Checks und Review voraus.
Externe Engine-Abnahme muss vor produktiver wissenschaftlicher Verwendung des
jeweiligen Pfads erfolgen. Der Draft wird nicht allein wegen lokaler Tests
automatisch freigegeben oder gemergt.
