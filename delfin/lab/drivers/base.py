"""The driver interface: the only code that talks to hardware.

A driver is called by the safety kernel and by nothing else. It knows how
to reach a device; it does not decide whether an action is allowed.

Methods (all per device of the manifest):

=====================  ======================================================
``read(device)``        ``{channel: Reading(value, unit, ts)}``; ``ts`` is when
                        the DEVICE produced the value (not when it was asked),
                        so a cached or frozen value shows as stale
``write(device, ch, v)``  one setpoint/state; returns a command id; raises
                        ``DriverError`` when the device refused or is gone
``safe_state(device)``  write the manifest's safe values, every channel, even
                        after an earlier one failed; returns the failures
``arm_watchdog(d, s)``  device-side timer: the DEVICE returns to its safe state
                        after ``s`` seconds unless re-armed (Shelly
                        ``toggle_after``, Tasmota ``PulseTime``, an OPC UA /
                        SiLA server-side lock timeout); False when unsupported
``disarm_watchdog(d)``  cancel the timer after a deliberate safe state
``health(device)``      ``Health(ok, detail)``; never raises
=====================  ======================================================

How the targets map onto it (design doc, section "Driver interface"):

* MQTT (Zigbee2MQTT, Tasmota): ``write`` publishes on the command topic,
  ``read`` takes the retained/last state message and its timestamp; the
  readback is the device's own state message, not our publish.
* Home Assistant: ``write`` = ``POST /api/services/<domain>/<service>``,
  ``read`` = ``GET /api/states/<entity_id>`` (``last_updated`` is ``ts``).
* Shelly Gen2: ``Switch.Set`` (with ``toggle_after`` as the device-side
  watchdog), ``Switch.GetStatus`` for the readback.
* libgpiod: a line request held for the process; the readback is the line
  value read back; no device-side timer (the kernel's dead-man and a
  hardware watchdog relay cover it).
* SiLA 2: ``write`` = an observable command (its execution UUID is the
  command id), ``read`` = observable properties; device-side lock with
  timeout = the Lock Controller feature.
* OPC UA LADS: ``write`` = a method call on a FunctionalUnit function or
  program, ``read`` = variable values with their source timestamps.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, NamedTuple, Optional


class DriverError(RuntimeError):
    """The device refused, timed out or could not be reached."""


class Reading(NamedTuple):
    value: Any
    unit: Optional[str]
    ts: float


@dataclass(frozen=True)
class Health:
    ok: bool
    detail: str = ""


@dataclass(frozen=True)
class DriverCapabilities:
    device_watchdog: bool = False     # a timer ON the device holds the safe state
    readback: bool = True             # the device reports what it is doing
    observable_commands: bool = False # long commands report progress (SiLA)
    device_lock: bool = False         # the device itself can be locked


class Driver(ABC):
    name: str = ""
    capabilities = DriverCapabilities()

    @abstractmethod
    def read(self, device) -> dict:
        ...

    @abstractmethod
    def write(self, device, channel: str, value: Any) -> str:
        ...

    def safe_state(self, device) -> list[str]:
        failures = []
        for ch, value in device.safe_state.items():
            try:
                self.write(device, ch, value)
            except Exception as exc:          # keep going: every channel
                failures.append(f"{ch}: {exc}")
        return failures

    def arm_watchdog(self, device, seconds: float) -> bool:
        return False

    def disarm_watchdog(self, device) -> None:
        return None

    def health(self, device) -> Health:
        return Health(True, "")
