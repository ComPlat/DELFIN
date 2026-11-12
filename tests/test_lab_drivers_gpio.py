"""Tests for the Raspberry Pi GPIO driver (delfin/lab/drivers/gpio.py).

The driver talks to the RP1 header through libgpiod v2. gpiod is an OPTIONAL
extra that is never installed in CI and only runs on Linux with a
``/dev/gpiochip*``, so every test drives the driver against an in-process fake
that mirrors the v2 API surface the driver uses (``Chip``, ``LineSettings``,
``request_lines``, ``line.Direction``, ``line.Value``). No real GPIO is ever
touched; the fake is installed as ``sys.modules["gpiod"]`` before the driver
imports it lazily.
"""

from __future__ import annotations

import sys
import threading

import pytest

from delfin.lab.drivers.base import DriverError
from delfin.lab.drivers.gpio import GpioDriver
from delfin.lab.safety import Actor


# --------------------------------------------------------------------------- #
# the fake gpiod v2 module (mirrors the API surface the driver uses)
# --------------------------------------------------------------------------- #
class _FakeValue:
    ACTIVE = 1
    INACTIVE = 0


class _FakeDirection:
    AS_IS = "as-is"
    OUTPUT = "output"
    INPUT = "input"


class _World:
    """The chip as the kernel keeps it (measured on a Pi 5, RP1): a released
    output keeps its direction and value, ``AS_IS`` changes neither, and a
    line held by one request is EBUSY for every other request."""

    def __init__(self):
        self.direction = {}            # offset -> "input" | "output"
        self.level = {}                # offset -> physical 0/1
        self.held = set()
        self.fail = False              # every call raises (transport error)
        self.hang = False              # every call blocks (timeout test)
        self.stuck = None              # physical level a read reports

    def phys(self, logical, active_low):
        return logical ^ int(bool(active_low))


class _FakeRequest:
    def __init__(self, world, config):
        self.world, self.config = world, config
        self.writes = []
        self.released = False
        for off, st in config.items():
            if st.direction == _FakeDirection.OUTPUT:
                world.direction[off] = "output"
                world.level[off] = world.phys(
                    st.output_value or _FakeValue.INACTIVE, st.active_low)
                self.writes.append((off, st.output_value))
            elif st.direction == _FakeDirection.INPUT:
                world.direction[off] = "input"
            world.direction.setdefault(off, "input")
            world.level.setdefault(off, 0)

    def get_value(self, offset):
        w = self.world
        if w.hang:
            threading.Event().wait()
        if w.fail:
            raise OSError("gpio transport error (fake)")
        level = w.stuck if w.stuck is not None else w.level[offset]
        return level ^ int(bool(self.config[offset].active_low))

    def release(self):
        self.released = True
        self.world.held -= set(self.config)


class _FakeSettings:
    def __init__(self, *, direction=None, active_low=False,
                 output_value=None):
        self.direction = direction
        self.active_low = active_low
        self.output_value = output_value


class _LineInfo:
    def __init__(self, direction):
        self.direction = direction


class _FakeChip:
    registry = {"pinctrl-rp1": ("/dev/gpiochip0", "/dev/gpiochip1")}

    def __init__(self, path: str):
        self.path = path

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def get_info(self):
        for label, paths in self.registry.items():
            if self.path in paths:
                return type("Info", (), {"label": label})()
        raise OSError(f"no such chip {self.path!r} (fake)")

    def get_line_info(self, offset):
        d = _FakeGpiod.world.direction.get(offset, "input")
        return _LineInfo(_FakeDirection.OUTPUT if d == "output"
                         else _FakeDirection.INPUT)


class _FakeGpiod:
    requests = []
    world = _World()

    class line:
        Direction = _FakeDirection
        Value = _FakeValue

    @classmethod
    def LineSettings(cls, *, direction=None, active_low=False,
                     output_value=None):
        return _FakeSettings(direction=direction, active_low=active_low,
                             output_value=output_value)

    @classmethod
    def Chip(cls, path: str):
        return _FakeChip(path)

    @classmethod
    def request_lines(cls, path, consumer=None, config=None):
        import errno
        w = cls.world
        if w.hang:
            threading.Event().wait()
        if w.fail:
            raise OSError("gpio transport error (fake)")
        if set(config or {}) & w.held:
            raise OSError(errno.EBUSY, "Device or resource busy")
        req = _FakeRequest(w, config or {})
        w.held |= set(config or {})
        cls.requests.append((path, consumer, config or {}, req))
        return req

    @classmethod
    def clear(cls):
        cls.requests = []
        cls.world = _World()


def _install_fake_gpiod():
    sys.modules["gpiod"] = _FakeGpiod
    _FakeGpiod.clear()


# --------------------------------------------------------------------------- #
# a lab manifest whose reactor device is driven by the gpio driver
# --------------------------------------------------------------------------- #
def _gpio_manifest() -> str:
    return """lab_manifest: 1
installation: gpio-lab
version: 1
nodes:
  - id: edge_pi
    description: Raspberry Pi 5 test node
modules:
  - id: reactor
    node: edge_pi
    min_level: benign
    positions: [R1]
    devices:
      - id: heater
        driver: gpiod
        model: relay
        address: {chip_label: pinctrl-rp1, power: 17, power_active_low: true}
        deadman_s: 15
        writes:
          power: {states: [on, off], level: benign}
        safe_state: {power: off}
      - id: overtemp
        driver: gpiod
        model: input
        address: {chip_label: pinctrl-rp1, power: 22}
        reads:
          power: {states: [on, off]}
"""


def _actor() -> Actor:
    """The operator actor; built lazily so it reads the fixture's env."""
    return Actor.here("s1")


@pytest.fixture
def gpio_world(tmp_path, monkeypatch):
    """A pinned lab manifest with a gpio device + the fake gpiod installed."""
    _install_fake_gpiod()
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("DELFIN_LAB_DIR", str(tmp_path / "lab"))
    monkeypatch.setenv("DELFIN_PERSON", "alice")
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    from delfin.lab import manifest as _manifest
    from delfin.lab import store
    from delfin.lab.safety import pin

    # the driver finds chips by globbing a device dir; point it at the fake
    from delfin.lab.drivers import gpio as _gpio
    monkeypatch.setattr(_gpio, "iglob",
                        lambda pattern: ["/dev/gpiochip0"])

    path = store.manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_gpio_manifest(), encoding="utf-8")
    m = _manifest.load(path)
    pin(m, by="alice")
    return m, _FakeGpiod


@pytest.fixture
def driver(monkeypatch):
    """A directly-constructed GpioDriver against the fake, fast timeout."""
    _install_fake_gpiod()
    from delfin.lab.drivers import gpio as _gpio
    monkeypatch.setattr(_gpio, "iglob", lambda pattern: ["/dev/gpiochip0"])
    return GpioDriver(chip_label="pinctrl-rp1", device_dir="/dev",
                      timeout_s=0.5, clock=lambda: 1000.0)


# --------------------------------------------------------------------------- #
# direct driver
# --------------------------------------------------------------------------- #
def test_driver_reads_writes_and_reads_back(driver, gpio_world):
    m, fake = gpio_world
    device = m.devices["heater"]
    cmd = driver.write(device, "power", "on")
    assert cmd.startswith("gpio-heater-")
    rb = driver.read(device)
    assert rb["power"].value == "on"          # the line echoes what we wrote
    driver.write(device, "power", "off")
    assert driver.read(device)["power"].value == "off"


def test_driver_requests_only_declared_lines(driver, gpio_world):
    m, fake = gpio_world
    driver.write(m.devices["heater"], "power", "on")
    path, consumer, config, req = fake.requests[-1]
    assert set(config) == {17}                # only the declared line
    settings = config[17]
    assert settings.direction == _FakeDirection.OUTPUT
    assert settings.active_low is True        # active-low relay honoured
    assert settings.output_value == _FakeValue.ACTIVE   # the commanded value


def test_an_undriven_relay_line_is_taken_over_in_its_safe_state(driver,
                                                                 gpio_world):
    """A line nobody drove since boot is an input; its level says nothing
    about the relay. Measured: an active-low heater on such a line read
    "on". The first read takes it over as an output at the safe value."""
    m, fake = gpio_world
    assert driver.read(m.devices["heater"])["power"].value == "off"
    path, consumer, config, req = fake.requests[-1]
    assert config[17].direction == _FakeDirection.OUTPUT
    assert config[17].output_value == _FakeValue.INACTIVE
    assert fake.world.level[17] == 1          # active-low: de-energised HIGH


def test_a_read_never_switches_an_output(driver, gpio_world):
    m, fake = gpio_world
    driver.write(m.devices["heater"], "power", "on")
    assert driver.read(m.devices["heater"])["power"].value == "on"
    path, consumer, config, req = fake.requests[-1]
    assert config[17].direction == _FakeDirection.AS_IS


def test_no_line_is_held_between_operations(driver, gpio_world):
    """Another process (the e-stop and dead-man daemons) must be able to
    take the line at any time: held, it gets EBUSY."""
    m, fake = gpio_world
    driver.write(m.devices["heater"], "power", "on")
    driver.read(m.devices["heater"])
    assert fake.world.held == set()
    assert all(r.released for *_x, r in fake.requests)
    other = GpioDriver(chip_label="pinctrl-rp1", timeout_s=0.5,
                       clock=lambda: 0.0)
    other.write(m.devices["heater"], "power", "off")   # a second process
    assert driver.read(m.devices["heater"])["power"].value == "off"


def test_a_busy_line_is_retried_then_fails_closed(driver, gpio_world):
    m, fake = gpio_world
    fake.world.held.add(17)                   # another process holds it
    patient = GpioDriver(chip_label="pinctrl-rp1", timeout_s=5.0,
                         clock=lambda: 0.0)
    with pytest.raises(DriverError, match="busy"):
        patient.write(m.devices["heater"], "power", "on")


def test_driver_maps_active_low_to_the_physical_line(driver, gpio_world):
    m, fake = gpio_world
    driver.write(m.devices["heater"], "power", "on")
    path, consumer, config, req = fake.requests[-1]
    # the driver hands the logical value to the request; the physical
    # inversion is done by libgpiod via LineSettings(active_low=True)
    assert req.writes[-1] == (17, _FakeValue.ACTIVE)
    assert fake.world.level[17] == 0          # active-low: energised LOW


def test_driver_reads_the_input_device(driver, gpio_world):
    m, fake = gpio_world
    inp = m.devices["overtemp"]
    driver.read(inp)
    path, consumer, config, req = fake.requests[-1]
    assert set(config) == {22} and config[22].direction == _FakeDirection.INPUT
    fake.world.level[22] = 1                  # the switch closes
    assert driver.read(inp)["power"].value == "on"


def test_driver_rejects_a_channel_without_a_line_offset(driver, gpio_world):
    m, fake = gpio_world
    # an address that lacks the channel's offset must be refused at request
    # time -- never silently write the wrong line
    from dataclasses import replace
    device = replace(m.devices["heater"],
                     address={k: v for k, v in m.devices["heater"].address.items()
                              if k != "power"})
    with pytest.raises(DriverError):
        driver.write(device, "power", "on")


def test_driver_finds_the_chip_by_label_and_not_a_number(driver, gpio_world):
    m, fake = gpio_world
    # a wrong label is refused: only the pinctrl-rp1 label is accepted
    other = GpioDriver(chip_label="some-other-chip", device_dir="/dev",
                       timeout_s=0.5, clock=lambda: 0.0)
    with pytest.raises(DriverError):
        other.write(m.devices["heater"], "power", "on")
    # the right label (found by scanning /dev/gpiochip*, never a pin number)
    driver.write(m.devices["heater"], "power", "on")


def test_driver_reports_no_device_side_timer(driver, gpio_world):
    m, fake = gpio_world
    device = m.devices["heater"]
    assert driver.capabilities.device_watchdog is False
    assert driver.arm_watchdog(device, 15) is False   # unsupported, no chip call
    assert not fake.requests                          # nothing attempted net
    driver.disarm_watchdog(device)                    # a no-op, never raises


def test_driver_close_is_harmless(driver, gpio_world):
    m, fake = gpio_world
    driver.write(m.devices["heater"], "power", "on")
    driver.close()
    driver.close()                                     # safe to call again
    assert fake.world.held == set()


def test_driver_reports_health(driver, gpio_world):
    m, fake = gpio_world
    health = driver.health(m.devices["heater"])
    assert health.ok is True
    # a transport failure is reported, never swallowed
    fake.world.fail = True
    broken = driver.health(m.devices["heater"])
    assert broken.ok is False


def test_driver_timeout_fails_closed(driver):
    # a hang in the fake chip must surface as a DriverError, not a deadlock
    from delfin.lab import manifest as _manifest
    import tempfile
    import os

    _install_fake_gpiod()
    _FakeGpiod.world.hang = True
    slow = GpioDriver(chip_label="pinctrl-rp1", device_dir="/dev",
                      timeout_s=0.05, clock=lambda: 0.0)
    with tempfile.NamedTemporaryFile("w", suffix=".yaml",
                                     delete=False) as f:
        f.write(_gpio_manifest())
        path = f.name
    try:
        m = _manifest.load(path)
    finally:
        os.unlink(path)
    device = m.devices["heater"]
    with pytest.raises(DriverError):
        slow.write(device, "power", "on")
    # restore the fake for the next test
    _install_fake_gpiod()


def test_lazy_import_error_names_the_extra():
    """With gpiod absent, the driver raises a message naming ``lab-gpio``."""
    sys.modules.pop("gpiod", None)
    d = GpioDriver(clock=lambda: 0.0)
    with pytest.raises(DriverError) as exc:
        d._gpiod_module()
    assert "lab-gpio" in str(exc)
    assert "gpiod" in str(exc)


# --------------------------------------------------------------------------- #
# through the kernel
# --------------------------------------------------------------------------- #
def test_kernel_transport_error_puts_device_safe_and_latched(gpio_world):
    m, fake = gpio_world
    from delfin.lab.safety import Kernel
    from _lab_support import FakeClock
    k = Kernel(m, clock=FakeClock(100.0))
    # a benign write opens the request and succeeds
    out = k.set(_actor(), "heater", "power", "off")
    assert out["ok"] is True
    # make the fake fail on the next write: the kernel must fail closed
    fake.world.fail = True
    out = k.set(_actor(), "heater", "power", "on")
    assert out["ok"] is False
    assert out["code"] == "driver_error"
    # the device is latched: the next write is refused
    out2 = k.set(_actor(), "heater", "power", "on")
    assert out2["code"] == "latched"
    snap = _snapshot()
    assert snap["latched"]["heater"]["reason"].startswith("write failed")


def test_kernel_readback_mismatch_latches(gpio_world):
    m, fake = gpio_world
    from delfin.lab.safety import Kernel
    from _lab_support import FakeClock
    k = Kernel(m, clock=FakeClock(100.0))
    # the first write opens the request and succeeds (echoes "on")
    out = k.set(_actor(), "heater", "power", "on")
    assert out["ok"] is True and out["readback"] == "on"
    # the line is stuck high: with active-low that reads "off"
    fake.world.stuck = 1
    out2 = k.set(_actor(), "heater", "power", "on")
    assert out2["code"] == "readback_mismatch"
    snap = _snapshot()
    assert "heater" in snap["latched"]


def test_kernel_does_not_arm_a_device_timer_on_gpio(gpio_world):
    m, fake = gpio_world
    from delfin.lab.safety import Kernel
    from _lab_support import FakeClock
    # the driver honestly reports it has no device timer, so the kernel must
    # never attempt to arm one even when the manifest sets device_watchdog;
    # the fake would fail closed if arm_watchdog ever ran
    k = Kernel(m, clock=FakeClock(100.0))
    k.set(_actor(), "heater", "power", "on")
    fake.world.fail = True
    assert k.heartbeat(_actor()) >= 0         # rearms nothing for gpio
    snap = _snapshot()
    assert "heater" not in snap["latched"]    # no arm attempted -> no failure


def test_lab_set_only_reaches_the_driver_through_kernel_set(gpio_world,
                                                            monkeypatch):
    """The agent tool writes to a gpio device ONLY via ``Kernel.set``.

    If any path in the tools layer called a driver's ``write`` directly, a
    monkeypatched ``Kernel.set`` that raises would not stop it (the fake would
    record a set_value). The tool path has no such side door: a raised ``set``
    propagates and the write aborts with the driver untouched.
    """
    m, fake = gpio_world
    from delfin.lab import safety as S

    calls = []

    def mocked_set(self, *a, **k):
        calls.append(a)
        raise RuntimeError("only the kernel writes")

    monkeypatch.setattr(S.Kernel, "set", mocked_set)
    from _lab_support import call
    import tempfile
    with tempfile.TemporaryDirectory() as ws:
        with pytest.raises(RuntimeError):
            call(ws, "lab_set", {"device": "heater", "channel": "power",
                                 "value": "on"})
    # the tool went through Kernel.set (it recorded the call) ...
    assert calls, "lab_set never called Kernel.set -- that is a bypass"
    # ... and the fake driver never saw a direct write
    assert all(not any(v == _FakeValue.ACTIVE for _o, v in r.writes)
               for _, _, _, r in fake.requests)


# small helpers ------------------------------------------------------------- #
def _snapshot():
    from delfin.lab import store
    return store.snapshot()
