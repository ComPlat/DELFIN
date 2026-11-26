"""Simulated devices with physics, faults and a device-side timer.

The simulated world lives in ``<lab_dir>/sim_world.json`` under ``flock``,
so every process sees one world: the agent, the e-stop CLI, the watchdog
daemon, and a process killed mid-action. Time is evaluated lazily: the
state at time *t* is computed from the last write and *t*, so a device
"keeps running" (a hotplate heats, a timer expires) while no process is
alive, which is what a real device does.

Kinds (``model`` in the manifest):

==========  ================================  ===============================
kind        writes                            reads
==========  ================================  ===============================
plug        power on/off                      power, power_w
relay       power on/off                      power
hotplate    setpoint degC, power on/off       temperature (first-order lag,
                                              tau 60 s), setpoint, power
reactor     setpoint degC, duration min,      yield %, mean_temperature,
            run idle/start                    setpoint, duration, run
stirrer     speed rpm, power on/off           speed, power
sensor      --                                temperature (follows the device
                                              named in ``address.follows``)
pump        flow mL/min, power on/off         flow, power
analyser    run idle/measure                  run, signal, barcode
reader      --                                barcode (container at position)
==========  ================================  ===============================

``reactor`` is the one kind whose reading DEPENDS ON THE CONDITIONS: writing
``run=start`` computes ``yield`` from the current ``setpoint`` and
``duration`` with :mod:`delfin.lab.surface`, so a loop can measure a response
surface on it. Every other kind answers what it was last told.

Every write channel reads back under the same name. Device-side watchdog:
``arm_watchdog`` sets a deadline; past it the device applies its safe
state by itself (like Shelly ``toggle_after``) even if nothing is running.

Faults for tests (``inject``): ``stuck`` (writes are acknowledged and do
nothing: the readback shows it), ``stale`` (timestamps freeze), ``offline``
(every call raises ``DriverError``).
"""

from __future__ import annotations

import math
import time
from typing import Any

from .. import store as _store
from .. import surface as _surface
from .base import Driver, DriverCapabilities, DriverError, Health, Reading

AMBIENT_C = 20.0
TAU_S = 60.0

_KIND_CHANNELS = {
    "plug": {"power": "off"},
    "relay": {"power": "off"},
    "hotplate": {"setpoint": AMBIENT_C, "power": "off"},
    "reactor": {"setpoint": AMBIENT_C, "duration": 0.0, "run": "idle"},
    "stirrer": {"speed": 0.0, "power": "off"},
    "sensor": {},
    "pump": {"flow": 0.0, "power": "off"},
    "analyser": {"run": "idle"},
    "reader": {"container": ""},
}

#: Kinds whose reading is COMPUTED from the conditions rather than remembered
#: from the last write. Only these can be the target of a closed loop: for
#: every other kind a measurement carries no information about the setpoint.
CONDITIONAL_KINDS = frozenset({"reactor"})


def _world():
    return _store.locked_json("sim_world.json", {"devices": {}})


def _entry(world: dict, device) -> dict:
    devs = world.setdefault("devices", {})
    e = devs.get(device.id)
    if e is None:
        kind = device.model or "relay"
        if kind not in _KIND_CHANNELS:
            raise DriverError(f"simulator has no kind {kind!r}")
        t = time.time()
        values = dict(_KIND_CHANNELS[kind])
        if kind in CONDITIONAL_KINDS:
            # A reactor starts at the conditions the manifest declares as its
            # resting state, not at the kind's artificial default. Starting
            # at 20 degC would put the device outside the envelope its own
            # manifest declares (30..110 degC), and every first step would
            # then be judged against a state the envelope does not admit.
            for ch, sv in (device.safe_state or {}).items():
                if ch in values:
                    values[ch] = sv
        e = {"kind": kind, "values": values,
             "temp": AMBIENT_C, "temp_ts": t, "ts": t, "fault": None,
             "watchdog": None, "watchdog_armed": None,
             "safe": dict(device.safe_state),
             "address": dict(device.address)}
        devs[device.id] = e
    return e


def _advance(world: dict, e: dict, t: float) -> None:
    """Bring one device to time ``t``: thermal lag and the device timer."""
    wd = e.get("watchdog")
    armed = e.get("watchdog_armed")
    # t >= wd: deadline reached. t < armed: the clock jumped backward below
    # the moment the watchdog was armed, so the deadline can no longer be
    # reached by this clock; fail safe instead of waiting for real time to
    # catch up with a timestamp that has already been passed.
    if wd is not None and (t >= wd or (armed is not None and t < armed)):
        _settle_temp(e, wd)
        e["values"].update(e.get("safe") or {})
        e["watchdog"] = None
        e["watchdog_fired"] = wd
        e["ts"] = wd
    _settle_temp(e, t)


def _settle_temp(e: dict, t: float) -> None:
    if e["kind"] != "hotplate":
        return
    v = e["values"]
    target = float(v.get("setpoint", AMBIENT_C)) if v.get("power") == "on" \
        else AMBIENT_C
    dt = max(0.0, t - float(e.get("temp_ts", t)))
    e["temp"] = target + (float(e["temp"]) - target) * math.exp(-dt / TAU_S)
    e["temp_ts"] = t


def _reactor_reading(e: dict, out: dict, ts: float) -> None:
    """Add the reactor's computed readings to ``out``.

    ``yield`` is what the LAST STARTED run produced, and it stays readable
    after ``run`` returns to idle so a loop can read the result of the run it
    started. A reactor that has never been started reports ``yield`` 0.0 --
    nothing has been made yet, which is not the same as a failed run (a
    failure is a refusal at the kernel, not a reading).

    ``mean_temperature`` is reported alongside so a loop's record can say
    which mean the surface was asked about, not just the setpoint.
    """
    at = e.get("yield_at") or {}
    if at:
        got = _surface.yield_percent(float(at["setpoint"]), float(at["duration"]))
        mean = _surface.mean_temperature(float(at["setpoint"]),
                                         float(at["duration"]))
    else:
        got, mean = 0.0, AMBIENT_C
    out["yield"] = Reading(round(float(got), 6), "%", ts)
    out["mean_temperature"] = Reading(round(float(mean), 6), "degC", ts)


class SimulatorDriver(Driver):
    name = "simulator"
    capabilities = DriverCapabilities(device_watchdog=True, readback=True,
                                      observable_commands=False,
                                      device_lock=False)

    def __init__(self, clock=time.time):
        self.clock = clock

    def _check(self, e: dict) -> None:
        if e.get("fault") == "offline":
            raise DriverError("device does not answer (simulated offline)")

    def read(self, device) -> dict:
        t = self.clock()
        with _world() as world:
            e = _entry(world, device)
            self._check(e)
            _advance(world, e, t)
            ts = float(e.get("stale_ts", e["ts"])) \
                if e.get("fault") == "stale" else t
            out = {}
            for ch, val in e["values"].items():
                spec = device.channel(ch)
                out[ch] = Reading(val, spec.unit if spec else None, ts)
            kind = e["kind"]
            if kind == "hotplate":
                out["temperature"] = Reading(round(e["temp"], 3), "degC", ts)
            if kind == "plug":
                out["power_w"] = Reading(
                    35.0 if e["values"].get("power") == "on" else 0.0, "W", ts)
            if kind == "sensor":
                temp = AMBIENT_C
                follows = (e.get("address") or {}).get("follows")
                other = world.get("devices", {}).get(follows) if follows else None
                if other is not None:
                    _advance(world, other, t)
                    temp = other.get("temp", AMBIENT_C)
                out["temperature"] = Reading(round(temp, 3), "degC", ts)
            if kind == "analyser":
                out["signal"] = Reading(float(e.get("signal", 0.0)), "au", ts)
            if kind == "reactor":
                _reactor_reading(e, out, ts)
            if kind in ("analyser", "reader"):
                # The barcode a position reports IS the container it holds: for
                # a reader that is the slot a hand-over writes, for an analyser
                # what a `place` put there. Keeping them one value means a
                # transfer cannot succeed in the registry and fail in the world.
                bc = e.get("barcode")
                if bc is None or bc == "":
                    bc = e["values"].get("container", "")
                out["barcode"] = Reading(str(bc), None, ts)
            return out

    def write(self, device, channel: str, value: Any) -> str:
        t = self.clock()
        with _world() as world:
            e = _entry(world, device)
            self._check(e)
            _advance(world, e, t)
            if channel not in e["values"]:
                raise DriverError(f"{device.id} has no channel {channel}")
            if e.get("fault") != "stuck":
                e["values"][channel] = value
                if e["kind"] == "analyser" and channel == "run" \
                        and value == "measure":
                    e["signal"] = float(e.get("signal_next", 1.0))
                if e["kind"] == "reactor" and channel == "run" \
                        and value == "start":
                    # The one place a measurement is computed from the
                    # conditions. `yield_at` keeps the conditions the run was
                    # started at, so a later setpoint write cannot retroactively
                    # change a result that was already reported.
                    e["yield_at"] = {
                        "setpoint": float(e["values"].get("setpoint", 0.0)),
                        "duration": float(e["values"].get("duration", 0.0))}
            if e.get("fault") != "stale":
                e["ts"] = t
            e["cmd"] = int(e.get("cmd", 0)) + 1
            return f"sim-{device.id}-{e['cmd']}"

    def arm_watchdog(self, device, seconds: float) -> bool:
        with _world() as world:
            e = _entry(world, device)
            self._check(e)
            e["watchdog"] = self.clock() + float(seconds)
            e["watchdog_armed"] = self.clock()
            e["safe"] = dict(device.safe_state)
        return True

    def disarm_watchdog(self, device) -> None:
        with _world() as world:
            e = _entry(world, device)
            e["watchdog"] = None
            e["watchdog_armed"] = None

    def health(self, device) -> Health:
        try:
            with _world() as world:
                e = _entry(world, device)
                if e.get("fault"):
                    return Health(False, f"fault: {e['fault']}")
            return Health(True, "simulated")
        except Exception as exc:
            return Health(False, str(exc))


# -- test and demo helpers (never reached from an agent tool) ---------------

def inject(device_id: str, fault, **fields) -> None:
    """Set (or clear with None) a fault on a simulated device; ``stale_ts``
    sets the frozen timestamp of a ``stale`` fault."""
    with _world() as world:
        e = world.setdefault("devices", {}).get(device_id)
        if e is None:
            raise KeyError(device_id)
        e["fault"] = fault
        e.update(fields)


def place(device_id: str, **fields) -> None:
    """Put simulator facts in place: ``barcode``, ``signal_next``."""
    with _world() as world:
        e = world.setdefault("devices", {}).get(device_id)
        if e is None:
            raise KeyError(device_id)
        e.update(fields)


def raw(device_id: str) -> dict:
    """The simulated device's own state (what is physically true)."""
    doc = _store.read_json(_store.lab_dir() / "sim_world.json", {"devices": {}})
    return dict(doc.get("devices", {}).get(device_id) or {})
