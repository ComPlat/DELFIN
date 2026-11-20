# Fundament-Abnahme und offene Arbeit

## Vier getrennte Statusstufen

| Stufe | Bedeutung | Kein Nachweis für |
| --- | --- | --- |
| Katalog/Export | Vertrag, Plan oder native Eingabe vorhanden | Engine-Ausführung |
| Implementiert | Adapter, Laufzeitweg und Ergebnisprüfung vorhanden | Installierbarkeit auf jedem Cluster |
| Schnittstelle getestet | Parser, Mocks, Fehler- und Identitätsfälle geprüft | Wissenschaftliche Richtigkeit eines echten Laufs |
| Extern abgenommen | Reale Installation und dokumentierter Referenzlauf auf benanntem System | Universelle Eignung für andere Chemie/Hardware |

Die bisherigen externen Engine-Grenztests nutzen überwiegend Mocks. Ein früherer
synthetischer Vina-I/O-Smoke ist kein Redocking-Benchmark. Hier werden keine
wissenschaftlichen Programme, Modell-Downloads oder Cluster-Jobs ausgeführt.

## Noch notwendige Arbeit bis zur umfassenden QM-Parität

| Priorität | Arbeit | Konkretes Abnahmekriterium |
| --- | --- | --- |
| 1 | Liganden-Parameterzuweisung | Maßgeblicher SDF-Graph, explizite Partialladungsmethode/-werte, lokale versionierte Kraftfelder, Atomzuordnung und Ladungsprüfung bis zum fertigen System |
| 1 | Installationsabdeckung | Ausführbare Adapter besitzen ein geprüftes Rezept oder dokumentierten externen Bereitstellungsvertrag; keine stillen Downloads |
| 1 | Echte Kernpfade | Vorbereitung → Docking → Validierung → Viewer und Systemaufbau → Dynamik → Analyse auf dem Zielcluster, inklusive kontrollierter Fehler |
| 1 | Visualisierungsabnahme | Browser/WebGL und Desktop-Exports mit richtigen Posen, Ketten, Einheiten, Konfidenzen und Dateien |
| 2 | Breitere Rezeptorchemie | Protonierungszustände, Cofaktoren, Metalle, Modifikationen und kovalente Bindungen mit explizitem Erhaltungs-/Änderungsprotokoll |
| 2 | Predictor-Ausführung | Je Engine überprüfte native Eingabe, Version, Assets, MSA-/Template-/Netzwerkpolitik, Ressourcen und geprüfter Output |
| 2 | Periodische Analyse/MD | Molekülrekonstruktion, Unwrapping, weitere Kraftfelder/Ensembles und NPT mit passenden Box-/Checkpoint-Verträgen |
| 2 | Ressourcen und Recovery | Pro-Knoten-Anforderungen, bounded Retry, Prozessabbruch und Wiederaufnahme ohne stillen Methoden- oder Chemiewechsel |
| 2 | QM/MM-Brücke | Auswahl, Grenzatome/Caps, Ladung/Multiplizität und Embedding mit überprüfbarer Bio↔QM-Atomzuordnung |
| 3 | Weitere Engines und Design | GROMACS/Amber/NAMD sowie Design-/weitere Docking-Backends mit eigenständigen Adaptern und Referenzen |
| 3 | Vergleichende Ergebnisse | Ensembles, Trajektorienansichten, Kontakt-Overlays und vergleichende Berichte mit konsistenter Identität |

Priorität 3 bedeutet Ausbau nach dem überprüfbaren Kern und ist keine Aussage,
dass diese Bereiche für das gewünschte umfassende Fundament entfallen.

## Externe Abnahmematrix

Jeder echte Test protokolliert Git-Revision, Betriebssystem, Scheduler,
Interpreter/Modulkontext, genaue Paket-/Engineversionen, Eingabe- und Asset-Hashes,
Ressourcen, Seeds, Ausgaben, Prüfungen und verbleibende Grenzen.

- **Installation:** passende Installation vorhanden; Tool fehlt; falsche Version;
  defekter verwalteter Prefix; Autoinstall deaktiviert; Netzwerk-/Solverfehler;
  zwei konkurrierende Installationsversuche; explizite Reparatur.
- **Cluster:** Modul stellt eigenes Python oder ein abhängiges Python bereit;
  Preflight und Worker importieren dieselben Pakete; Slurm-Exitcode, Timeout und
  Abbruch werden korrekt weitergegeben.
- **Docking:** kleiner dokumentierter Redocking-Fall mit Referenzpose,
  explizitem Rezeptor-/Ligandzustand, Box und Seed; Scoreeinheit, Mapping,
  RMSD und fehlgeschlagene Posen prüfen.
- **MD:** gültiges Referenzprotein im unterstützten Teilbereich; Vakuum und
  solvatisierte feste Box; Masse/Ladung/Constraints/PME/Einheiten prüfen;
  kurze Dynamik, Fortsetzung und absichtlich inkompatibler Checkpoint.
- **Analyse:** bekannte Translation/Rotation, bekannte RMSD-/RMSF-Werte,
  Auswahl-/Gewichtungswechsel und periodischer Grenzübertritt. Bis zur
  Unwrapping-Unterstützung muss dessen Einschränkung sichtbar bleiben.
- **Prediction:** lokale verifizierte Assets, passende Sequenz und
  Konfidenzzuordnung; fehlendes Asset und Versionsabweichung verweigern;
  Laufzeit/GPU-Speicher dokumentieren.
- **Viewer:** lokale Dateien ohne Netzwerk, Pose-/Ketten-/Wasserumschaltung,
  Bindungsbox, Konfidenz und Desktop-Skripte prüfen; Browser und Version nennen.
- **Resume:** erfolgreicher unveränderter Wiederanlauf; Änderungen an Eingabe,
  Einstellungen, Code, Laufzeit oder verschachtelten Artefakten werden erkannt.

## Integrationsstand

Am 20. November wurden die 46 gezielten System-/Periodik-/Dynamik-/Analysetests
erneut erfolgreich ausgeführt. Das sind Schnittstellentests ohne echte OpenMM-
oder MDAnalysis-Engineausführung.

Nach Ergänzung des PoseBusters-Installers und der Importprüfung für explizite
Interpreter ergibt die breitere Bio-/QM-Runtime-/CLI-/Tool-Vertragssuite 322 bestandene Tests
und drei bereits vorher reproduzierte Fehler: fehlende Pflichtparameterprüfung
für `orca_th_opi.output_file`, fehlende Vertragskategorie bei `orca_qmmm_opi`
und ein abgeschnittener Producer-Hinweis im Plattformtest. Diese Fehler werden
nicht als neue Periodik-Regressionen oder als bestandene Prüfungen ausgewiesen.

Der zuletzt integrierte `main`-Stand ist
`03a02d048f674a972d6c2e2551908eb5ece24b6c`. Der am 20. November abgefragte aktuelle
Stand `4e936b999de8cb6b60dcfc1a2b85c8b74171f3d8` ist noch nicht integriert:
direkter Git-Zugriff scheitert in der Entwicklungsumgebung an der DNS-Auflösung,
und die verfügbare GitHub-Dateianbindung kann ein benötigtes Binärartefakt nicht
übertragen. Vor externer Abnahme und vor einem Merge muss dieser Abgleich samt
Konfliktprüfung und passenden Tests nachgeholt werden. PR #746 bleibt Draft.
