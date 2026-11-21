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
| 1 | Liganden-Parameterzuweisung extern abnehmen | Der implementierte SMIRNOFF-Pfad mit gelieferten Ladungen benötigt echte Referenzläufe; automatische Ladungserzeugung und weitere Kraftfelder bleiben Ausbau |
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

Die zuvor bekannten Fehler bei `orca_th_opi.output_file`, der Kategorie von
`orca_qmmm_opi` und abgeschnittenen Producer-Hinweisen sind im Review-Stand
behoben. Die zugehörige Vertrags-/Plattform-/ORCA-Suite besteht mit 98 bestandenen
und sieben übersprungenen Tests.

`main` ist bis `2e3e2e1d083879fb079aa1ed3591d7db03be23dd` integriert, ohne
Merge-Konflikte. Der vollständige Snapshot einschließlich Binärdateien wurde
über GitHub übertragen und anhand der Git-Objekthashes geprüft. Der frühere
Übertragungsblocker ist damit aufgehoben. Die lokale Kopie hält diesen
Main-Commit als flache Historiengrenze; der veröffentlichte Merge referenziert
den echten Main-Commit als zweiten Elterncommit.

Die aktuelle Review-Abnahme und ihre Umgebungsgrenzen stehen in [REVIEW.md](REVIEW.md).
PR #746 bleibt bis zur erfolgreichen CI-/Review-/Cluster-Abnahme Draft.

Der echte Vina-I/O-Smoke gehört jetzt ausdrücklich zur manuellen Cluster-Abnahme:

```sh
python -m examples.biomolecular_foundation.vina_smoke --output /shared/new-vina-smoke
```

Er verlangt ein bereits installiertes Vina und einen neuen Ausgabeordner. Er
läuft nicht durch gewöhnliche pytest-Sammlung und wurde hier nicht gestartet.
Die regulären Tests verwenden kontrollierte Engine-Grenzen; ein zufällig auf
dem Testrechner vorhandenes Vina führt zu keiner echten Rechnung.
