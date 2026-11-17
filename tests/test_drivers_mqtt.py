"""The MQTT device driver: a real driver behind the safety kernel.

Every test runs against an in-process fake broker (no network, no paho-mqtt
installed). The driver is exercised two ways:

* **directly** -- a ``MqttDriver`` with an injected fake transport, against a
  parsed manifest :class:`~delfin.lab.manifest.Device`: ``read``, ``write``
  with readback, QoS 1, retained state on connect, the device-side timer,
  the availability (LWT) rule that an offline device is unknown, health and
  credentials resolved only by reference.
* **through the kernel** -- the real ``_DocToolExecutor`` path a model's
  tool call takes: every transport error makes the kernel go safe (fail
  closed), and the agent's only write tool routes through ``Kernel.set``.

The manifest address for an MQTT device is flat (the schema forbids nested
addresses and inline secrets):

    address:
      host: 192.0.2.50            # broker (not a secret)
      port: 8883
      topic_base: lab/relay1      # topics derived below
      secret_ref: cred:mqtt_1     # username/password by reference
      tls: true
      tls_ca: /etc/delfin/ca.pem  # path, not a secret
      arm_topic: lab/relay1/arm   # optional device-side timer command topic

Topics per channel follow ``{topic_base}/cmd/<channel>`` (command, QoS 1,
not retained), ``{topic_base}/state/<channel>`` (the device's retained
state -- the readback), and ``{topic_base}/availability`` (the LWT; a
device that is offline is unknown).
"""

from __future__ import annotations

import inspect
import time

import pytest
import yaml

from delfin.lab import manifest
from delfin.lab.drivers.base import DriverError, Reading
from delfin.lab.drivers import mqtt

from _lab_support import lab_fixture  # noqa: F401  (the fixture is named `lab`)

#: topic fragments that separate command from state.
_CMD = "/cmd/"
_STATE = "/state/"


class FakeBroker:
    """A working in-process relay broker.

    Publishing to ``.../cmd/<ch>`` makes the device report that value on
    ``.../state/<ch>`` (a relay does what it is told and reports it back),
    which is what gives the kernel a readback. ``offline`` puts the device's
    LWT availability topic to False; ``publish_error`` makes every command
    fail with a transport error (the device stopped answering).
    """

    def __init__(self, clock):
        self.clock = clock
        self._retained: dict[str, str] = {}
        self._stamp: dict[str, float] = {}
        #: avail topic -> online (LWT). Default online.
        self.online: dict[str, bool] = {}
        self.publish_error: Exception | None = None
        self.published: list[tuple] = []        # (topic, payload, retain, qos)
        self.commands: list[str] = []
        self.subscribed: list[str] = []
        self._mid = 0

    # -- driver transport interface -----------------------------------------
    def open(self):
        pass

    def close(self):
        pass

    def subscribe(self, topic: str, qos: int = 1):
        self.subscribed.append(topic)

    def publish(self, topic: str, payload: str, qos: int = 1,
                retain: bool = False) -> str:
        if self.publish_error is not None:
            raise self.publish_error
        self._mid += 1
        self.published.append((topic, payload, retain, qos))
        mid = f"msg-{self._mid}"
        if _CMD in topic:
            # the relay obeys the command and reports its state back
            state_topic = topic.replace(_CMD, _STATE)
            self._retained[state_topic] = payload
            self._stamp[state_topic] = self.clock()
            self.commands.append(payload)
        return mid

    def retained(self, topic: str):
        """(payload, ts) of the device's last state on ``topic`` or None."""
        if topic in self._retained:
            return self._retained[topic], self._stamp[topic]
        return None

    def available(self, topic: str) -> bool:
        return self.online.get(topic, True)

    # -- test helpers ---------------------------------------------------------
    def device_offline(self, avail_topic: str):
        self.online[avail_topic] = False


class _Clock:
    def __init__(self, t0: float = 1000.0):
        self.t = t0

    def __call__(self) -> float:
        return self.t

    def advance(self, s: float) -> None:
        self.t += s


def _device_doc() -> dict:
    return {
        "id": "relay1", "driver": "mqtt", "model": "generic-relay",
        "address": {"host": "192.0.2.50", "port": 8883,
                    "topic_base": "lab/relay1", "secret_ref": "cred:mqtt_1"},
        "stale_after_s": 5, "deadman_s": 30,
        "reads": {"power": {"states": ["on", "off"]}},
        "writes": {"power": {"states": ["on", "off"], "level": "benign"}},
        "safe_state": {"power": "off"},
    }


def _manifest_with(*devices: dict):
    doc = {
        "lab_manifest": 1, "installation": "mqtt-lab", "version": 1,
        "nodes": [{"id": "edge"}],
        "modules": [{"id": "synth", "node": "edge", "min_level": "benign",
                     "devices": list(devices)}],
    }
    return manifest.parse_document(doc)


def _relay1(m=None):
    m = m or _manifest_with(_device_doc())
    return m.devices["relay1"]


# -- direct driver tests --------------------------------------------------------

def test_read_returns_the_devices_retained_state_with_its_timestamp():
    clock = _Clock()
    fake = FakeBroker(clock)
    # the device reported it is on at t0 and has not changed since
    fake.publish("lab/relay1/cmd/power", "on")
    clock.advance(2)                       # we ask two seconds later
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    got = drv.read(_relay1())
    # the Reading carries WHEN THE DEVICE produced the value (t0), not when
    # it was asked -- so a cached/frozen device shows as stale in the kernel
    assert got["power"] == Reading("on", None, 1000.0)


def test_subscribes_to_every_state_topic_and_availability_on_connect():
    clock = _Clock()
    fake = FakeBroker(clock)
    fake.publish("lab/relay1/cmd/power", "off")
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    drv.read(_relay1())
    for want in ("lab/relay1/state/power", "lab/relay1/availability"):
        assert want in fake.subscribed


def test_write_publishes_qos1_unretained_to_the_command_topic():
    clock = _Clock()
    fake = FakeBroker(clock)
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    dev = _relay1()
    cmd = drv.write(dev, "power", "on")
    assert cmd.startswith("msg-")
    topic, payload, retain, qos = fake.published[-1]
    assert topic == "lab/relay1/cmd/power"
    assert payload == "on"
    assert qos == 1
    assert retain is False
    # and the device reported its new state back -- the readback the kernel
    # will use to confirm the write
    rb = drv.read(dev)["power"]
    assert rb.value == "on"


def test_device_that_is_offline_is_unknown_not_valuable():
    clock = _Clock()
    fake = FakeBroker(clock)
    fake.publish("lab/relay1/cmd/power", "off")
    fake.device_offline("lab/relay1/availability")
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    with pytest.raises(DriverError, match="offline"):
        drv.read(_relay1())


def test_every_transport_error_raises_drivererror():
    clock = _Clock()
    fake = FakeBroker(clock)
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    fake.publish_error = DriverError("broker refused the command")
    with pytest.raises(DriverError, match="broker refused"):
        drv.write(_relay1(), "power", "on")
    fake.publish_error = DriverError("no answer")
    with pytest.raises(DriverError):
        drv.read(_relay1())


def test_arm_watchdog_sets_the_device_side_timer():
    clock = _Clock()
    fake = FakeBroker(clock)
    doc = _device_doc()
    doc["address"]["arm_topic"] = "lab/relay1/arm"
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    assert drv.arm_watchdog(_relay1(_manifest_with(doc)), 30) is True
    assert ("lab/relay1/arm", "30", False, 1) in fake.published


def test_arm_watchdog_is_false_when_the_firmware_has_no_timer():
    clock = _Clock()
    fake = FakeBroker(clock)
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    assert drv.arm_watchdog(_relay1(), 30) is False


def test_safe_state_writes_every_safe_channel_past_failures():
    clock = _Clock()
    fake = FakeBroker(clock)
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    fake.publish_error = DriverError("unit is off")
    # must not raise: every channel is attempted, failures are reported
    failures = drv.safe_state(_relay1())
    assert failures


def test_health_reports_the_device_available_and_never_raises():
    clock = _Clock()
    fake = FakeBroker(clock)
    fake.publish("lab/relay1/cmd/power", "off")
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    assert drv.health(_relay1()).ok
    fake.device_offline("lab/relay1/availability")
    h = drv.health(_relay1())
    assert not h.ok and "offline" in h.detail


def test_credentials_resolved_by_reference_only(monkeypatch):
    clock = _Clock()
    fake = FakeBroker(clock)
    seen = []
    monkeypatch.setattr(
        mqtt, "credentials",
        lambda name: seen.append(name) or
        '{"username": "edge", "password": "s3cr3t"}')
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    resolved = drv.credentials_for(_relay1())
    assert seen == ["mqtt_1"]                     # the store was asked once
    assert resolved["username"] == "edge"
    assert resolved["password"] == "s3cr3t"
    assert resolved["tls_ca"] is None


# -- through the kernel ---------------------------------------------------------

def _kernel_manifest() -> dict:
    return {
        "lab_manifest": 1, "installation": "mqtt-lab", "version": 1,
        "nodes": [{"id": "edge"}],
        "modules": [{"id": "synth", "node": "edge", "min_level": "benign",
                     "devices": [_device_doc()]}],
    }


def test_transport_error_on_a_write_makes_the_kernel_go_safe(lab, monkeypatch):
    from delfin.lab import safety, store
    from _lab_support import call

    fake = FakeBroker(lambda: time.time())
    monkeypatch.setattr(mqtt, "broker_factory",
                        lambda device, **kw: fake)
    fake.publish("lab/relay1/cmd/power", "off")

    path = store.manifest_path()
    path.write_text(yaml.safe_dump(_kernel_manifest()), encoding="utf-8")
    safety.pin(manifest.load(path), by="alice")

    out = call(lab, "lab_set", {"device": "relay1", "channel": "power",
                                "value": "on"})
    assert out.get("ok"), out

    # the device stops answering: the write fails and the kernel goes safe
    fake.publish_error = DriverError("relay gone")
    out2 = call(lab, "lab_set", {"device": "relay1", "channel": "power",
                                 "value": "off"})
    assert out2.get("code") == "driver_error", out2
    assert "safe state" in out2.get("reason", "")
    assert store.snapshot()["latched"].get("relay1")


def test_a_write_impossible_without_kernel_set_is_not_in_the_tool_surface():
    from delfin.lab import tools
    assert "lab_set" in tools.TOOL_NAMES
    # no tool name offers a direct write / publish / drive through the agent
    assert not (tools.TOOL_NAMES & {"lab_write", "lab_publish", "lab_drive",
                                    "lab_write_raw", "lab_publish_raw"})
    # the only write tool routes through Kernel.set and nothing else writes
    src = inspect.getsource(tools)
    assert "kernel.set(" in src
    assert ".publish(" not in src and "write_to" not in src


# -- connection per operation, settled before use -------------------------------

class _LateBroker(FakeBroker):
    """A broker whose retained messages arrive only after the subscribe, as
    on a real connection: ``wait_for`` is what delivers them."""

    def __init__(self, clock):
        super().__init__(clock)
        self.pending: dict[str, str] = {"lab/relay1/availability": "online",
                                        "lab/relay1/state/power": "off"}
        self.delivered: dict[str, str] = {}
        self.closed = 0
        self.waited: list = []

    def wait_for(self, topics, timeout_s):
        self.waited.append(("for", tuple(topics), timeout_s))
        self.delivered.update(self.pending)
        for t, v in self.pending.items():
            if "/state/" in t and t not in self._retained:
                self._retained[t], self._stamp[t] = v, self.clock()
        return True

    def wait_newer(self, topic, since, timeout_s):
        self.waited.append(("newer", topic, since))
        return self._stamp.get(topic, 0.0) > since

    def available(self, topic):
        return self.delivered.get(topic) == "online" and super().available(topic)

    def close(self):
        self.closed += 1


def test_a_read_waits_for_the_retained_availability_and_state():
    """Before the wait, the device looked offline and every write was
    refused, so an e-stop could not switch a reachable relay off."""
    clock = _Clock()
    fake = _LateBroker(clock)
    drv = mqtt.MqttDriver(clock=clock, transport_factory=lambda device: fake)
    assert drv.read(_relay1())["power"].value == "off"
    kind, topics, timeout = fake.waited[0]
    assert kind == "for" and "lab/relay1/availability" in topics
    assert "lab/relay1/state/power" in topics and timeout == mqtt._SETTLE_S


def test_a_write_waits_for_the_devices_answer():
    clock = _Clock()
    fake = _LateBroker(clock)
    drv = mqtt.MqttDriver(clock=clock, transport_factory=lambda device: fake)
    drv.read(_relay1())
    clock.advance(1)
    drv.write(_relay1(), "power", "on")
    assert ("newer", "lab/relay1/state/power", 1000.0) in fake.waited
    assert fake.retained("lab/relay1/state/power")[0] == "on"


def test_every_operation_closes_the_connection_it_built():
    clock = _Clock()
    fake = _LateBroker(clock)
    drv = mqtt.MqttDriver(clock=clock, transport_factory=lambda device: fake)
    d = _relay1()
    drv.read(d)
    drv.write(d, "power", "on")
    drv.health(d)
    assert fake.closed == 3


def test_an_injected_transport_is_not_closed():
    clock = _Clock()
    fake = _LateBroker(clock)
    drv = mqtt.MqttDriver(clock=clock, transport=fake)
    drv.read(_relay1())
    assert fake.closed == 0


def test_a_refused_login_is_reported_at_once_without_the_password(monkeypatch):
    """paho reports a refused login in on_connect; the wait ends there."""
    broker = mqtt._PahoBroker.__new__(mqtt._PahoBroker)
    broker._last_err = None
    broker._connected = False

    class _Code:
        is_failure = True

        def __str__(self):
            return "Not authorized"

    broker._on_connect(None, None, None, _Code())
    assert broker._connected is False and broker._last_err == "Not authorized"

    class _Ok:
        is_failure = False

    broker._on_connect(None, None, None, _Ok())
    assert broker._connected is True
