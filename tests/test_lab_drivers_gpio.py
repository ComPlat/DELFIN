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
    OUTPUT = "output"
    INPUT = "input"


class _FakeRequest:
    def __init__(self, config):
        self.config = config
        self.values = {off: (s.output_value if s.output_value is not None
                             else _FakeValue.INACTIVE)
                       for off, s in config.items()}
        self.fail = False          # raise on the next call (transport error)
        self.hang = False          # block forever (timeout test)
        self.stuck = None          # get_value over-ride to simulate a stuck line
        self.closed = False
        self.writes = []           # (offset, value) recorded on set_value

    def set_value(self, offset, value):
        if self.hang:
            threading.Event().wait()
        if self.fail:
            raise OSError("gpio transport error (fake)")
        self.writes.append((offset, value))
        self.values[offset] = value

    def get_value(self, offset):
        if self.hang:
            threading.Event().wait()
        if self.fail:
            raise OSError("gpio transport error (fake)")
        if self.stuck is not None:
            return self.stuck
        return self.values[offset]

    def close(self):
        self.closed = True


class _FakeSettings:
    def __init__(self, *, direction=None, active_low=False,
                 output_value=None):
        self.direction = direction
        self.active_low = active_low
        self.output_value = output_value


class _FakeChip:
    registry = {"pinctrl-rp1": ("/dev/gpiochip0", "/dev/gpiochip1")}

    def __init__(self, path: str):
        self.path = path
        self.closed = False

    def get_label(self) -> str:
        for label, paths in self.registry.items():
            if self.path in paths:
                return label
        raise OSError(f"no such chip {self.path!r} (fake)")

    def close(self):
        self.closed = True


class _FakeGpiod:
    requests = []

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
        req = _FakeRequest(config or {})
        cls.requests.append((path, consumer, config or {}, req))
        return req

    @classmethod
    def clear(cls):
        cls.requests = []


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
    assert settings.output_value == _FakeValue.INACTIVE  # de-energise-to-trip


def test_driver_maps_active_low_to_the_physical_line(driver, gpio_world):
    m, fake = gpio_world
    driver.write(m.devices["heater"], "power", "on")
    path, consumer, config, req = fake.requests[-1]
    # the driver hands the logical value to the request; the physical
    # inversion is done by libgpiod via LineSettings(active_low=True)
    assert req.writes[-1] == (17, _FakeValue.ACTIVE)


def test_driver_reads_the_input_device(driver, gpio_world):
    m, fake = gpio_world
    inp = m.devices["overtemp"]
    # the first read opens the request for the input line
    driver.read(inp)
    path, consumer, config, req = fake.requests[-1]
    assert set(config) == {22} and config[22].direction == _FakeDirection.INPUT
    # pre-set the line value so the read reflects the device state
    req.values[22] = _FakeValue.ACTIVE
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


def test_driver_releases_lines_on_close(driver, gpio_world):
    m, fake = gpio_world
    driver.write(m.devices["heater"], "power", "on")
    path, consumer, config, req = fake.requests[-1]
    assert req.closed is False
    driver.close()
    assert req.closed is True                          # released on exit
    driver.close()                                     # safe to call again


def test_driver_reports_health(driver, gpio_world):
    m, fake = gpio_world
    health = driver.health(m.devices["heater"])
    assert health.ok is True
    # a transport failure is reported, never swallowed
    for p, c, cfg, r in fake.requests:
        r.fail = True
    broken = driver.health(m.devices["heater"])
    assert broken.ok is False


def test_driver_timeout_fails_closed(driver):
    # a hang in the fake chip must surface as a DriverError, not a deadlock
    from delfin.lab import manifest as _manifest
    import tempfile
    import os

    class Hanger(_FakeGpiod):
        @classmethod
        def request_lines(cls, path, consumer=None, config=None):
            req = _FakeRequest(config or {})
            req.hang = True
            cls.requests.append((path, consumer, config or {}, req))
            return req

    sys.modules["gpiod"] = Hanger
    Hanger.clear()
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
    for p, c, cfg, r in fake.requests:
        r.fail = True
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
    # make the line stuck "off": ready to verify, the device says it is off
    for p, c, cfg, r in fake.requests:
        r.stuck = _FakeValue.INACTIVE
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
    for p, c, cfg, r in fake.requests:
        r.fail = True
    n = k.heartbeat(_actor())                 # rearms nothing for gpio
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
    assert all(len(r.writes) == 0 for _, _, _, r in fake.requests)


# small helpers ------------------------------------------------------------- #
def _snapshot():
    from delfin.lab import store
    return store.snapshot()
