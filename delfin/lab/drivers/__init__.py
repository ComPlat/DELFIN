"""Device drivers: simulator, Shelly Gen2 RPC and Tasmota are implemented."""

from __future__ import annotations

from .base import (Driver, DriverCapabilities, DriverError, Health,  # noqa: F401
                   Reading)

#: driver name -> class. Real drivers register here when they exist.
REGISTRY: dict = {}


def _register_builtin() -> None:
    from .simulator import SimulatorDriver
    from .shelly_tasmota import ShellyDriver, TasmotaDriver
    REGISTRY.setdefault("simulator", SimulatorDriver)
    REGISTRY.setdefault("shelly_rpc", ShellyDriver)
    REGISTRY.setdefault("tasmota", TasmotaDriver)


def make(name: str) -> Driver:
    _register_builtin()
    cls = REGISTRY.get(name)
    if cls is None:
        raise DriverError(f"no driver {name!r} on this node")
    return cls()


#: What each driver class can do on the device itself (design doc table).
DEVICE_WATCHDOG_SUPPORT = {
    "simulator": "yes (deadline in the simulated device)",
    "shelly_rpc": "yes: Switch.Set toggle_after",
    "tasmota": "yes: PulseTime",
    "mqtt": "only when the device firmware has a timer (Tasmota yes, "
            "Zigbee plugs mostly no) (planned)",
    "home_assistant": "no: HA automations run on the HA host, not on the "
                      "device (planned)",
    "gpiod": "no: a line keeps its value while the chip is powered; use a "
             "hardware watchdog relay (planned)",
    "sila2": "lock timeout yes; safe state on timeout depends on the "
             "server (planned)",
    "opcua_lads": "server-dependent (planned)",
}
