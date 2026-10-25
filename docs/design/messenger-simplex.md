# Messenger for the DELFIN agents

Status: bridge, gateway and provisioning implemented; real-phone acceptance pending.

## Goal
Talk to the DELFIN agents from a phone: get told when an agent is done,
failed, has a question or needs an approval; answer questions and give
new tasks by replying. Open source, end-to-end encrypted, free.

## Choice: SimpleX Chat (user decision)
- Open source (AGPL): app, protocol, relays, CLI. End-to-end encrypted. Free.
- No phone number for anyone. The bot has no public address; it is paired
  once by a one-time link / QR code shown on the Pi.
- Ready aarch64 CLI: `simplex-chat-ubuntu-24_04-aarch64` (v7.0.3, 2026-09-26);
  pinned version + SHA-256 in the provisioning (repo rule, no curl | sh).
- Outbound only (connects to SMP relays); the Pi's firewall stays closed.
- Rejected for now: Signal -- the bot needs its own number for good, and
  signal-cli's libsignal has no linux-aarch64 build (0.14.9 bundles amd64
  .so only), so libsignal 0.103.0 would be built from source with Rust >=
  1.93 (Debian trixie has 1.85) and rebuilt regularly. WhatsApp -- not E2E
  on the Business API, needs a public webhook; unofficial libraries break
  the terms of service.

## Architecture
```
phone (SimpleX app) <-> SimpleX relays (E2E, no content) <-> Pi (outbound only)
  [simplex-chat daemon]  local API reachable by delfin-gw only (nftables meta skuid)
  [gateway]              user delfin-gw: no sudo, no access to delfin's home
        | unix socket; the bridge accepts a peer only if SO_PEERCRED uid == delfin-gw
  [bridge in the agent service]  user delfin: sessions, attention inbox
```

## What the phone can do (user decision: tasks + answers)
- free text -> task for the "phone session"; `@name text` -> a named session
- reply to a question (`1`, `2`, or text) -> answers the agent's question
- `status`, `sitzungen` / `sessions`, `stop` (stops all agents), `hilfe` / `help`
- the agent reports: done, failed, question, approval needed
- approvals (writes, commands, more rights) stay at dashboard/terminal with the
  one-time approval code; the message says where to approve

## Threats and measures
| threat | measure |
|---|---|
| strangers message the bot | no public address; only the paired contact; other contact requests refused and logged; security code compared at pairing |
| the agent forges "user" messages | bridge checks the peer uid (kernel); daemon API port only for delfin-gw (nftables); bot data 0700 under delfin-gw |
| the agent sends secrets to the phone | redaction (service Redactor), length limit, no file contents; the agent cannot choose a recipient |
| lost phone | `stop`; `delfin-gw revoke` on the Pi unpairs; approvals never went via phone |
| bug in the messenger program | runs as delfin-gw without rights; systemd hardening (NoNewPrivileges, ProtectSystem=strict, PrivateTmp, no capabilities); reaches only the narrow bridge |
| flooding, replay | rate limit, message ids processed once, size limit |
| tampered download | pinned version, SHA-256 checked |

## Phases
0. Feasibility (with the user's phone): simplex-chat on the Pi, bot API, one pairing.
1. Bridge in the agent service (product, PR from main): socket, peer check, commands,
   questions/results to the gateway. Tests.
2. Gateway (product): messenger <-> bridge, command parsing, redaction, limits;
   tested against a fake messenger API.
3. Provisioning `scripts/pi/pi_messenger.sh` (on the scripts/pi stack): user delfin-gw,
   pinned download, systemd units with hardening, nftables rule, pairing command,
   README checklist.
4. Red team (separate agent, not a builder): forge messages as delfin, reach the daemon
   port, replay, flood, get an approval through the phone, extract secrets. All must fail.
5. Acceptance on the Pi with the user's phone: task, answer a question, notification, stop.
   Then PRs; pushes only on the user's order.

Phases 1-3 run in parallel (swarm, separate worktrees); every result is checked on the Pi
before it goes anywhere.

## One-time for the user
Install the SimpleX app (Android / iOS / F-Droid), scan the QR code shown on the Pi,
compare the security code.

## Auftrag für DELFIN (diese Sitzung)

Setze die Phasen 1–3 um, in diesem Worktree auf dem Branch
`agent/messenger-simplex`. Phase 0 und 5 (echtes Handy) und Phase 4
(Red Team) machen andere; baue aber so, dass sie prüfbar sind.

Umfang:
1. Bridge im Agent-Service (`delfin/agent/service/`): ein Unix-Socket mit
   Rechten 0660 in einem Verzeichnis, das der Gruppe `delfin-gw` gehört;
   eine Verbindung wird nur angenommen, wenn `SO_PEERCRED` die UID des
   Nutzers `delfin-gw` liefert (konfigurierbar, im Test überschreibbar).
   Feste Befehle, JSON-Zeilen: `status`, `sessions`, `pending` (wartende
   Fragen und Freigaben aus der attention-Inbox), `answer` (beantwortet eine
   FRAGE, niemals eine Freigabe), `prompt` (Auftrag an eine Sitzung), `stop`
   (wie `delfin-agent stop-all`), und ein Strom ausgehender Ereignisse
   (fertig, gescheitert, Frage, Freigabe nötig). Jede Ausgabe läuft durch
   den vorhandenen Redactor (`service/redact.py`) und wird gekürzt.
2. Gateway (neues Paket `delfin/messenger/`): spricht mit dem simplex-chat-
   Daemon über dessen lokale WebSocket-API (nur 127.0.0.1) und mit der Bridge;
   nimmt nur Nachrichten des einen gekoppelten Kontakts an (Kontakt-ID
   gespeichert beim Koppeln); lehnt Kontaktanfragen ab; Befehle `status`,
   `sitzungen`/`sessions`, `stop`, `hilfe`/`help`, `@name text`, sonst
   Auftrag an die Handy-Sitzung; Antwort auf eine Frage per Nummer oder
   Text; Ratenlimit, Größenlimit, jede Nachrichten-ID nur einmal.
   Einstiegspunkte: `python -m delfin.messenger pair|run|revoke|status`.
3. Provisionierung `scripts/pi/pi_messenger.sh` (läuft als root, idempotent):
   Nutzer `delfin-gw` ohne sudo; simplex-chat `v7.0.3`
   `simplex-chat-ubuntu-24_04-aarch64` GEPINNT mit SHA-256 (die Summe beim
   Schreiben einmal ermitteln und eintragen; kein `curl | sh`, siehe
   tests/test_micromamba_is_fetched_as_one_pinned_release.py); systemd-
   Dienste für Daemon und Gateway mit Härtung (NoNewPrivileges,
   ProtectSystem=strict, ProtectHome, PrivateTmp, CapabilityBoundingSet=,
   RestrictAddressFamilies); eine nftables-Regel, die den API-Port nur für
   `meta skuid delfin-gw` öffnet; README-Abschnitt und Eintrag in der
   Checkliste "After the first boot". Keine festen Namen in /tmp
   (tests/test_a_pi_script_runs_nothing_from_a_fixed_tmp_name.py).

Regeln:
- Tests für alles, gegen einen nachgebauten Messenger und eine nachgebaute
  Bridge; kein Netz in den Tests. Mindestens: fremder Absender wird verworfen,
  falsche UID an der Bridge wird abgewiesen, `answer` auf eine Freigabe wird
  abgelehnt, Wiederholung derselben Nachrichten-ID wirkt nur einmal, Secrets
  in Ausgaben werden geschwärzt.
- Keine privaten Daten (Nummern, IPs, E-Mails, Namen) in Code oder Skripten.
- Kleine, thematisch saubere Commits, Nachrichten rein technisch, am Ende
  jeder Nachricht die Zeile: Co-Authored-By: Claude <noreply@anthropic.com>
- Nichts außerhalb dieses Worktrees ändern. Kein Push ohne Freigabe; wenn
  alles grün ist: Branch pushen und den Compare-Link nennen.

## Deployment and regression checks

The provisioner uses port 5227 for both the daemon and gateway. Install DELFIN
with the `messenger` extra (`pip install ".[agent,messenger]"` from the checkout) in `/opt/delfin-messenger/bin/venv`
before provisioning. The agent service must already run as `delfin` (override
with `MESSENGER_SERVICE_USER`); `MESSENGER_BRIDGE_PYTHON` can select that user's
Python installation. The script provisions a bridge unit as that user, with
primary group `delfin-gw`, and a separate 0750 runtime directory at
`/run/delfin-messenger-bridge`. Its socket is 0660. The agent service's own
runtime directory remains private. The bridge uses the host's `/tmp` to reach
the agent service and restarts after service connection failures.

Pair as `delfin-gw`, using the same state directory as the gateway:

```sh
sudo -u delfin-gw env DELFIN_MESSENGER_DIR=/var/lib/delfin-messenger/pairing /opt/delfin-messenger/bin/venv/bin/python -m delfin.messenger pair
```

Use the same environment and executable with `revoke` to remove the pairing.
The running gateway checks the pairing before accepting commands and sending
notifications; a command already handed to the agent service is not undone.
Restart the gateway after pairing a different phone. Approvals remain at the
terminal/dashboard.

Regression tests cover unsolicited bridge events, gateway notifications,
revocation during a running gateway, real service summary fields, named-session
routing, and agreement between generated units on port, socket and state path.
A real SimpleX daemon and phone are still required for deployment acceptance.
