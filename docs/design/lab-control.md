# Lab control: a safety kernel between the model and the hardware

Status: design + phase 2 core (simulator only). Nothing here is a
certification claim. IEC 61508, ISO 13849 and IEC 60204-1 are used for their
ideas (safe state, watchdog, two-channel stop, reset does not restart), not as
standards this software complies with. A safety function that protects people
is a hardware matter: an emergency-stop button that cuts power, a thermal fuse
and a fume-hood airflow switch. This software adds a second, independent
layer. It never replaces the first.

**The rule the whole feature rests on.** The model proposes actions, and a
deterministic kernel in the executor decides them. Nothing the model writes
can move a device outside its declared envelope: no argument, no permission
mode (including bypass/all_free), no prompt and no file content. The kernel
does not read the permission mode at all.

## 1. Scope

| In scope (phase 2, this branch) | Later phases |
|---|---|
| `delfin/lab/`: manifest, units, interlocks, the safety kernel, dead-man, e-stop, audit | real drivers (Shelly, Tasmota/MQTT, Home Assistant, libgpiod): phase 3 |
| simulator driver with physics, faults and a device-side timer | SiLA 2 and OPC UA LADS drivers, discovery: phase 4 |
| the coordinator–edge interface (signed envelopes, replay window, manifest pin) | an mTLS transport and node enrolment: phase 3 |
| modules, a sequential/graph workflow, sample hand-over, resume, run record | parallel scheduler across many modules: phase 5 |
| five agent tools; `delfin-agent lab ...`; systemd units for the Pi | Chemotion upload (export file only now); operator view in the dashboard: phase 5 |

The reference case is a modular Materials Acceleration Platform such as KIT's
ChemASAP: storage, synthesis, purification and analysis modules for automated,
parallel organic synthesis, with data planned to flow into Chemotion
[chemasap]. A Raspberry Pi 5 is the first edge node and the demo. It is not
the permanent home of the system.

## 2. Architecture

```
 clients: dashboard / TUI / phone / CLI   (operator view, approvals, e-stop)
        |                     \
 coordinator (any host)        \  e-stop: CLI + daemon on EVERY edge, no network needed
   plans, schedules, run record \
        | signed envelope (HMAC per node; mTLS in phase 3)
 edge node (Pi, lab PC, module controller)   -- one per module or several modules
   EdgeNode: MAC, replay window, manifest pin check
   SafetyKernel (the ONLY path to a driver write)
     manifest -> units -> range/step/rate -> lock (person+session) -> interlocks on
     FRESH readings -> approval -> intent -> write -> readback -> result
   dead-man daemon + e-stop daemon (separate processes / systemd units)
   drivers: simulator | mqtt | home_assistant | shelly_rpc | tasmota | gpiod | sila2 | opcua_lads
        |
 devices (separate VLAN, no cloud)
```

| Layer | File | What it decides |
|---|---|---|
| Driver | `delfin/lab/drivers/base.py`, `simulator.py` | how to reach a device; it decides nothing |
| Manifest | `delfin/lab/manifest.py` | what exists, its units, limits, safe state, interlocks, levels |
| Units | `delfin/lab/units.py` | a closed table; a missing, unknown or ambiguous unit is refused |
| Interlocks | `delfin/lab/expr.py` | a closed grammar with three-valued logic; UNKNOWN fails closed |
| Kernel | `delfin/lab/safety.py` | every write; nothing else imports a driver's `write` |
| Approvals | `delfin/lab/approvals.py` | single use, bound to the action, signed by the person |
| Dead-man / e-stop | `delfin/lab/watchdog.py`, `estop.py` | safe state on a missed heartbeat, a dead owner process or an e-stop |
| Edge interface | `delfin/lab/edge.py` | the coordinator never reaches a kernel without a valid envelope |
| Workflow | `delfin/lab/workflow.py`, `record.py` | modules, sample hand-over, resume, run record, export |
| Tools | `delfin/lab/tools.py` (wired in `api_client`) | five tools, advertised only where a manifest exists |
| CLI | `delfin/lab/cli.py` (`delfin-agent lab ...`) | the person's side: approve, clear, pin, e-stop, daemons, smoke |

**Modules.** A module is a group of devices with its own manifest section,
interlocks, positions and approval floor (`min_level`). A workflow is a
graph of module operations (`set`, `wait`, `measure`, `transfer`) with
sample hand-over (sample id, container id, position). The coordinator takes
the **module lock** for a run. A step's device write still goes through the
edge kernel's device lock and every check.

**Mapping to SiLA 2 and OPC UA LADS.** These are first-class targets, and
the interface was cut to fit them:

| DELFIN concept | SiLA 2 [sila-proto, sila-lock] | OPC UA LADS [lads] |
|---|---|---|
| node / edge | SiLA server (discovered via mDNS `_sila._tcp`) | LADS Device (OPC UA server) |
| module | a server's set of features | a FunctionalUnit (own state machine) |
| `read` with source timestamp | observable property | variable with SourceTimestamp |
| `write` -> command id | observable command -> CommandExecutionUUID + ExecutionInfo | method on a Function / Program -> result |
| device lock (person+session) | LockController `LockServer(id, timeout)`, id as metadata on every call | LADS uses OPC UA sessions/roles; no lock feature (UNVERIFIED whether a profile adds one) |
| simulator-first | SimulationController `StartSimulationMode` / `StartRealMode` | none in the base spec (UNVERIFIED) |
| authenticated channel | TLS mandatory (Part B, quoted second-hand: UNVERIFIED), AuthenticationService | SecureChannel, SignAndEncrypt, X.509 application certificates |

## 3. Threat model

Every row has a rule in the code and a test through the real tool/executor
path (`_DocToolExecutor.execute`) or the real process boundary. Each test has
a control that shows the allowed case passing.

| # | Threat | Rule that catches it | Test (tests/…) |
|---|---|---|---|
| T1 | wrong command: out of range, too big a step, too fast | manifest `min/max/max_step/min_interval_s`, checked before any write | `test_lab_safety_kernel.py::test_out_of_range_setpoint_is_refused` |
| T2 | hallucinated unit, or none | closed unit table; no unit means refused; other dimension refused | `::test_missing_or_wrong_unit_is_refused` |
| T3 | hallucinated device or channel | manifest lookup, with the known ids in the refusal | `::test_unknown_device_is_refused` |
| T4 | missing interlock (heater without stirrer or hood) | every interlock that names the channel is evaluated with the proposed value, both directions | `::test_missing_interlock_is_refused`, `::test_switching_off_a_required_device_is_refused` |
| T5 | stale or missing sensor | a reading older than `stale_after_s` is UNKNOWN, and UNKNOWN fails closed | `::test_stale_sensor_fails_closed` |
| T6 | a hung agent | dead-man: no heartbeat in `deadman_s` puts the device in its safe state and latches it | `test_lab_deadman_and_estop.py::test_missed_heartbeat_goes_safe` |
| T7 | a crashed agent or host | dead owner pid means safe state; a device-side timer (Shelly `toggle_after`, Tasmota `PulseTime`) holds without any process | `::test_killed_agent_mid_action_goes_safe`, `::test_device_side_watchdog_holds_without_any_process` |
| T8 | network loss / coordinator–edge partition | the edge's own dead-man; the edge never waits for the coordinator to be safe | `test_lab_edge_and_workflow.py::test_partition_goes_safe_on_the_edge` |
| T9 | malicious prompt or file content ("ignore the limits") | the kernel reads only the manifest; nothing in a tool argument widens a limit; bypass does not reach the kernel | `test_lab_safety_kernel.py::test_prompt_injection_in_a_file_does_not_widen_limits`, `::test_hazardous_needs_approval_even_in_bypass` |
| T10 | forbidden action | level `forbidden` is never executed, approval or not | `::test_forbidden_is_never_executed` |
| T11 | another person on a shared box | locks and approvals name the person; another person's lock is never taken and their action never approved | `::test_another_person_cannot_take_a_lock_or_approve` |
| T12 | the agent answers its own approval | `lab approve` / `lab clear` / `estop --clear` refused with the agent-session marker; the approval's MAC uses a key an agent child cannot obtain; `approval_guard` refuses the command text | `test_lab_tools_registered.py::test_the_agent_cannot_approve_or_clear` |
| T13 | replay (approval or envelope) | an approval is single use and bound to device, channel, value, person and session; an envelope has a nonce and a ±30 s window | `test_lab_safety_kernel.py::test_replayed_approval_is_refused`, `test_lab_edge_and_workflow.py::test_edge_refuses_bad_mac_replay_and_wrong_manifest` |
| T14 | race between two sessions | per-device lock under `flock`; the second session is refused, never queued silently | `test_lab_safety_kernel.py::test_two_sessions_race_for_one_device` |
| T15 | compromised device on the LAN (lies about its state) | readback mismatch puts it in safe state and latches it; interlocks should name an INDEPENDENT sensor (airflow switch, thermocouple), not the device's own report; VLAN isolation | `::test_readback_mismatch_goes_safe_and_latches` (the lie itself is caught only by a second channel) |
| T16 | compromised edge | an edge's key signs nothing another edge accepts; edges never accept edge-to-edge commands; a coordinator treats an edge's results as data (run record), not as permission | `test_lab_edge_and_workflow.py::test_an_edge_key_is_not_accepted_by_another_edge` |
| T17 | manifest tampered with on disk or in transit | the edge pins a checksum (`lab manifest pin`, the person's step); a different checksum, locally or in an envelope, refuses every write | `test_lab_safety_kernel.py::test_a_tampered_manifest_is_refused`, `test_lab_edge_and_workflow.py::test_edge_refuses_bad_mac_replay_and_wrong_manifest` |
| T18 | sample mix-up between modules | the registry knows where every sample is; a transfer from the wrong module, to a taken position, or with a container whose barcode does not match stops the run | `test_lab_edge_and_workflow.py::test_sample_mixup_stops_the_run` |
| T19 | two coordinators | the module lock names the coordinator actor; the second is refused | `::test_two_coordinators_cannot_hold_one_module` |
| T20 | e-stop while the agent is busy or hung | the e-stop is a separate process (CLI + daemon) that writes the stop flag first, then drives every device safe; writes are refused until a person clears it | `test_lab_deadman_and_estop.py::test_estop_while_the_agent_is_busy`, `::test_writes_blocked_after_estop_until_cleared` |

**Residual risks (not covered by software).** A device that reports a false
state and has no independent sensor beside it (T15). A person with the
account's password (same Unix user in people mode: the person is an
attribution and a kernel rule, not an OS boundary; see §9). A Pi that loses
power while a GPIO relay is energised: this is safe only if the relay is
de-energise-to-trip, so the wiring must make "off" the unpowered state.

## 4. Approval levels

| Level | Who decides | Bypass / all_free | Audit |
|---|---|---|---|
| read (free) | nobody | n/a | no (reads are frequent; the run record keeps them) |
| benign | the kernel's checks | same | intent, result, readback |
| hazardous | a person, for exactly this action (device, channel, value, unit, session), single use, expires in 300 s | **still asks**: the kernel does not read the mode | as benign + approval id + approver |
| forbidden | nobody, never | same | refusal |

A module's `min_level` raises every write in the module (an analysis module
with a laser might be `hazardous` throughout). How it works: `lab_set`
without an approval returns `approval_required` with a request id, and the
person runs `delfin-agent lab approve <id>`. That command is refused inside
an agent session and refused for anyone but the requesting person. The grant
is an HMAC under the person's trust-store key (`trust_store.derived_key`),
which returns nothing inside an agent session. The model then repeats the
call with `approval=<id>`. The kernel checks the MAC, the binding, the
expiry and that the grant is unused, consumes it, and only then writes.

A write of a channel's own **safe value** (a hotplate setpoint back to
20 °C, a pump to 0 mL/min) is treated as benign. It moves the device toward
its safe state, so a cool-down never waits for a person. Every other check
still applies, and the write is audited.

A workflow pauses at a hazardous step (`waiting_approval`, with the request
id in the run's checkpoint). After the person approves, `resume` repeats
that step with the grant. Approving a whole plan in one step is phase 5.

## 5. Dead-man

| Item | Design |
|---|---|
| obligation | every device with a hazardous write that is not in its safe state has a controlling actor (person, session, pid, host) |
| heartbeat | every lab tool call of that session, and `Kernel.heartbeat` (an operator view or a coordinator may call it); default `deadman_s` = 30 s per device |
| trip | heartbeat older than `deadman_s`, or the owner pid dead (same host), or the edge losing its coordinator (the edge sees no heartbeat) |
| action | driver `safe_state`, readback verified, device **latched** (writes refused until a person runs `lab clear <device>`), lock revoked, security event + audit |
| who checks | `delfin-agent lab daemon --role deadman` (systemd `delfin-lab-watchdog.service`); the agent loop never checks itself |
| device side | armed at every hazardous write and every heartbeat for `deadman_s`: the simulator has it; Shelly `toggle_after` [shelly-switch] and Tasmota `PulseTime` [tasmota] have it; Home Assistant, plain Zigbee plugs and libgpiod lines do not (use a hardware watchdog relay) |

Open issue: model latency. A slow model can take minutes per turn, and a turn
does not heartbeat. Under a 30 s dead-man a long heating step trips unless
something else heartbeats. The options are a per-device timeout, the
operator's view as the heartbeat (a person watching), or a coordinator that
heartbeats only while its run is healthy. This is a decision for the user,
listed in §11.

## 6. Emergency stop

| Item | Design |
|---|---|
| scopes | global, or one module |
| triggers | `delfin-agent lab estop [--module M]`, the `lab_estop` tool, `delfin.lab.estop.trigger()` for dashboard/TUI/phone, later a GPIO button on the edge |
| order | 1. stop flag written (under `flock`; every write is refused from here on) 2. every device in scope driven to its safe state, continuing past failures 3. locks revoked, pending approvals cancelled 4. audit + security event |
| a write in flight | the kernel checks the flag again after every driver write. If the stop arrived while a write was in flight (the stop's own pass may have run first), the write is undone at once by putting the device in its safe state, and the call reports `estop` |
| ownership | its own process: the CLI does the stop itself, and `delfin-lab-estop.service` re-asserts the safe state every 2 s while the stop is active, reporting a readback that does not hold. The agent loop is not on this path |
| reset | `lab estop --clear` by a person (refused inside an agent session). Clearing restarts nothing: devices stay in their safe state (IEC 60204-1 9.2.3.4: a reset "shall not restart the machinery but only permit restarting" [abb-estop]) |
| two channels | the software stop is channel 2. Channel 1 is a hardware e-stop that cuts power to the devices (category 0 stop [kollmorgen]); a two-channel safety relay is the established design [rs-relay] |

## 7. Audit

The installation keeps `<lab_dir>/audit.log`, written through DELFIN's
hash-chained, keyed `audit_log.append` (tool `lab`). Each write leaves three
records in order (`intent`, `result`, `readback`), plus `refused`,
`safe_state`, `estop`, `estop_cleared`, `deadman`, `approval_requested`,
`approval_granted`, `approval_consumed`, `latched` and `cleared`. The
security-event panel receives the alarming ones under the new kinds
`lab_refused`, `lab_safe_state`, `lab_estop`, `lab_estop_cleared`,
`lab_deadman`, `lab_readback_mismatch`, `lab_approval_refused` and
`lab_manifest_refused`. `tests/test_lab_tools_registered.py` pins those kinds.

## 8. Data: the run record and Chemotion

`delfin/lab/record.py` writes, per run, `runs/<run_id>/record.json`
(schema `delfin.lab.run/1`). It holds the plan and its digest, the manifest
version and checksum, the actor, every step's request, outcome and readback,
measurements, deviations, the samples and their positions, and FAIR metadata
(identifier, creator, licence, provenance and timestamps [fair]).
`ml_train` (other branch) reads these records.

Chemotion: `to_chemotion()` writes `chemotion_export.json`, a research-plan /
LabIMotion-shaped payload, and nothing is uploaded yet.

| Fact | Status |
|---|---|
| Grape API under `/api/v1`; Bearer API token (`chemtoken_…`) or JWT | verified in source [chemotion-api] |
| `POST authentication/token` (username, password) returns a token | verified in source; mount path and expiry UNVERIFIED |
| `/samples` CRUD, `/samples/import`; separate APIs for reactions, research plans, containers, attachments | verified in source |
| collection export/import as zip jobs | verified; zip layout UNVERIFIED |
| LabIMotion generic elements / segments / datasets | verified [labimotion] |
| Repository mints DataCite DOIs and exposes OAI-PMH | secondary source [nfdi4chem] |
| how analyses/containers/datasets nest | UNVERIFIED: inferred from file names only |

## 9. Deployment: edge on a Raspberry Pi 5, and a plain workstation

| Item | Pi 5 (Debian arm64) | workstation |
|---|---|---|
| e-stop owner | `delfin-lab-estop.service` | `delfin-agent lab daemon --role estop` (foreground, or `--detach`: double fork, pid and log file in the lab directory) |
| dead-man | `delfin-lab-watchdog.service` | `lab daemon --role deadman` |
| installed by | `scripts/pi/pi_install.sh` as systemd **user** units (agent user, linger enabled by the admin), idempotent, restarted on deploy | by hand |
| hardening | `NoNewPrivileges`, `ProtectSystem=strict`, `ProtectHome=read-only` + `ReadWritePaths=<lab_dir>`, `PrivateTmp`, `RestrictAddressFamilies=AF_UNIX`, `MemoryDenyWriteExecute`, `SystemCallArchitectures=native` | n/a |
| smoke | `python -m delfin.lab.smoke` in every deploy: set, interlock, readback, e-stop, clear on the simulator; prints `smoke lab ok` | the same |

GPIO on the Pi 5. The header pins sit on the RP1 chip. RPi.GPIO uses
`/dev/mem` register access and does not work there; this is reported by
secondary sources [rpi-gpio], and gpiozero's documentation says "only lgpio
works on the Pi 5" [gpiozero]. The phase 3 driver uses libgpiod v2 (`gpiod`
on PyPI, `request_lines` / `LineSettings`) [gpiod]. Since kernel PR #6144
(Aug 2024) the RP1 lines moved from gpiochip4 to gpiochip0 [rpi-6144], so
the driver finds the chip by its label `pinctrl-rp1` and never by number.

Admin step: `sudo usermod -aG gpio <agent user>`. Raspberry Pi OS's udev
rule gives `/dev/gpiochip*` to group `gpio`, mode 0660 [rpi-udev]. The group
grants the **whole chip**, not single lines. "Only the declared lines" is
enforced by the driver, which requests only the manifest's lines and holds
them for the process, so a second requester gets EBUSY. Against the agent's
own shell this is not a boundary. A deployment that needs one runs the edge
kernel as a separate `delfin-lab` user in the `gpio` group, with the agent
talking to it over a Unix socket (phase 3). This is a deviation from the
requested "agent user gets only the gpio group".

People mode. The manifest is per installation: the lab directory
resolves to `<account>/.delfin/lab` for a HOME of the form
`<account>/people/<name>`, never to the person's own HOME. Locks and
approvals name the person. One person's session never takes another
person's lock (`locked`, the holder named), and a grant for another
person's request is refused. On top of that, each person's grant is signed
with that person's own trust-store key, so it does not verify for someone
else's request. Test: `test_lab_safety_kernel.py::test_another_person_cannot_take_a_lock_or_approve`.

Failure and resume. A failed step stops the run (`failed`). The run keeps
its module locks so it can resume, and the record is written anyway,
because a failure is data. `resume` continues from the failed step and
never repeats a completed one. A device the edge latched (a readback
mismatch, a dead-man trip) stays latched until a person runs `lab clear
<device>`. `Coordinator.abandon(run_id)` releases a run that will not
resume.

## 10. Network isolation

| Advice | Why |
|---|---|
| devices on their own VLAN or segment; the edge is the only host that reaches it | a compromised device cannot reach the coordinator or the internet |
| no cloud accounts on devices (Shelly cloud off, Tasmota without cloud) | commands only come from the edge |
| MQTT: Mosquitto with `certfile/keyfile/cafile`, `require_certificate true`, `use_identity_as_username true`, an ACL file with one `user` per device and only its own `topic` lines, `allow_anonymous false` [mosquitto] | a device can publish only its own state |
| no plain MQTT for commands across hosts; coordinator–edge over mTLS (phase 3) plus the envelope MAC (now) | the envelope MAC holds even if TLS is misconfigured |
| Shelly Gen2 digest auth (SHA-256) on [shelly-auth]; Tasmota HTTP sends the password in clear, so Tasmota is used over TLS MQTT only [tasmota] | |
| Home Assistant: a long-lived token for one dedicated non-admin user (REST `Authorization: Bearer`, WebSocket `auth` message) [ha-rest, ha-ws] | |

## 11. Adding a device type

1. Write `delfin/lab/drivers/<name>.py`, a subclass of `Driver`: `read` (with
   the device's own timestamp), `write` (returns a command id), and
   `arm_watchdog` if the device has a timer.
2. Register it in `drivers.REGISTRY` and move its name from `DRIVERS_PLANNED`
   to `DRIVERS_IMPLEMENTED`.
3. Add a row to `DEVICE_WATCHDOG_SUPPORT` and to §5.
4. Give it a simulator twin, or a recorded-response fake, for tests. Tests
   never touch hardware.
5. The first run on real hardware is a separate step, and a person approves
   the plan.

## 12. Phases

| Phase | Deliverable |
|---|---|
| 1 | this document |
| 2 | edge kernel, simulator, dead-man, e-stop, audit, five tools, coordinator–edge interface (in-process transport), simulated two-module workflow (synthesis -> analysis with sample hand-over), Pi units + smoke |
| 3 | mTLS transport and node enrolment; separate `delfin-lab` user + socket; Shelly, Tasmota (TLS MQTT), Home Assistant and gpiod drivers; first real-device demo after a person approves the plan |
| 4 | SiLA 2 client driver (sila2 / UniteLabs), discovery, LockController and SimulationController mapping; OPC UA LADS via asyncua (LGPL-3.0) |
| 5 | scheduler for parallel runs across modules; operator view (all modules, per-module and global e-stop); Chemotion upload; ml_train loop |

## 13. Decisions open for the user

* Approval level per device class. Proposal: heaters above 40 °C, pumps
  moving solvent, anything with a laser or high voltage hazardous; plugs for
  lights and fans benign; mains switching of unknown loads forbidden.
* The dead-man timeout, and what heartbeats during a slow model turn (§5).
* What to buy or wire first for a demo. Proposal: one Shelly Plus Plug S
  (has `toggle_after`), one relay HAT on the Pi with de-energise-to-trip
  wiring and a lamp as the load, one DS18B20 sensor, and a mushroom e-stop
  in series with the relay supply.

## References

* [sila-proto] SiLA 2 framework protobuf, https://gitlab.com/SiLA2/sila_base/-/raw/master/protobuf/SiLAFramework.proto
* [sila-lock] LockController v2.0, https://gitlab.com/SiLA2/sila_base/-/raw/master/feature_definitions/org/silastandard/core/LockController-v2_0.sila.xml ; core features (AuthenticationService, AuthorizationService, SimulationController), https://gitlab.com/SiLA2/sila_base/-/tree/master/feature_definitions/org/silastandard/core ; Python `sila2` 0.14.0 (maintenance-only), https://pypi.org/project/sila2/
* [lads] OPC 30500 LADS 1.0.0 (2023-11-30), https://reference.opcfoundation.org/LADS/v100/docs/ ; model, https://reference.opcfoundation.org/specs/OPC-30500-1/4 ; security modes, https://reference.opcfoundation.org/Core/Part4/7.20
* [abb-estop] reset does not restart (IEC 60204-1 9.2.3.4, secondary), https://help.plc.abb.com/sf_emergencystop.html
* [kollmorgen] stop categories 0/1, https://www.kollmorgen.com/en-us/developer-network/stop-and-emergency-stop-function
* [rs-relay] two-channel e-stop with safety relay (secondary), https://docs.rs-online.com/15c9/A700000013733811.pdf
* [ha-rest] https://contributors.home-assistant.io/core/internals/rest-api/ ; [ha-ws] https://contributors.home-assistant.io/core/internals/websocket-api/
* [shelly-switch] https://shelly-api-docs.shelly.cloud/gen2/ComponentsAndServices/Switch ; [shelly-auth] https://shelly-api-docs.shelly.cloud/gen2/General/Authentication
* [tasmota] https://tasmota.github.io/docs/Commands/ (PulseTime: 1..111 in 0.1 s, 112..64900 = seconds + 100)
* [gpiozero] https://gpiozero.readthedocs.io/en/stable/api_pins.html ; [gpiod] https://pypi.org/project/gpiod/ ; [rpi-6144] https://github.com/raspberrypi/linux/pull/6144 ; [rpi-udev] https://github.com/RPi-Distro/raspberrypi-sys-mods ; [rpi-gpio] secondary, https://raspberry.tips/en/raspberrypi-tutorials/raspberry-pi-5-gpio-rpigpio-not-working-alternatives
* [mosquitto] https://mosquitto.org/man/mosquitto-conf-5.html ; MQTT 5.0, https://docs.oasis-open.org/mqtt/mqtt/v5.0/os/mqtt-v5.0-os.html
* [chemotion-api] https://github.com/ComPlat/chemotion_ELN/blob/main/app/api/api.rb ; https://github.com/ComPlat/chemotion_ELN/tree/main/app/api/chemotion ; [labimotion] https://docs.chemotion.net/labimotion/user ; [nfdi4chem] https://knowledgebase.nfdi4chem.de/knowledge_base/docs/chemotion_repository/
* [chemasap] https://www.chemasap.kit.edu/
* Self-driving labs: Nat. Synth. 2, 483 (2023), doi:10.1038/s44160-022-00231-0 ; Chem. Rev. 124, 9633 (2024), doi:10.1021/acs.chemrev.4c00055 ; safety monitoring preprints (weight accordingly): arXiv:2508.05148, arXiv:2602.15061
* [fair] Sci. Data 3, 160018 (2016), doi:10.1038/sdata.2016.18 ; AnIML https://animl.org/overview ; Allotrope https://allotropefoundation.org/
