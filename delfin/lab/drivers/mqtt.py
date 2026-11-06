"""An MQTT device driver: generic topic-based relays and sensors.

The driver is the ONLY code that talks to the broker. It decides nothing:
a read returns what the device reported, a write publishes a command and
reports what the device says it did (the readback). Every transport error
raises :class:`DriverError`, so the safety kernel fails closed.

Transport
---------

``paho-mqtt`` is an optional dependency (extra ``lab-mqtt``). It is
imported lazily the first time a real broker is needed; tests inject a
fake broker and never import it. When paho is missing the driver raises a
:class:`DriverError` that names the ``lab-mqtt`` extra.

The broker seam an injection replaces (``broker_factory``) is a small
object protocol:

    open() -> None                      connect; raises DriverError
    subscribe(topic, qos=1) -> None     QoS 1 subscribe
    publish(topic, payload, qos=1, retain=False) -> str   QoS 1, returns an id
    retained(topic) -> (payload, ts) | None   the device's last state (readback)
    available(topic) -> bool            the (LWT) availability state
    close() -> None

Topics from the manifest (the address is flat; the schema forbids nested
addresses and inline secrets)::

    topic_base: lab/relay1         # from address.topic_base
      {topic_base}/cmd/<channel>   # command: publish, QoS 1, not retained
      {topic_base}/state/<channel> # retained device state -- the readback
      {topic_base}/availability    # LWT: a device that is offline is unknown
      address.arm_topic            # optional device-side timer command topic

Credentials
-----------

Resolved by reference only, from DELFIN's store
(:func:`delfin.agent.credentials.load_credential`), at driver start, never
by the kernel and never inline in a manifest. ``address.secret_ref`` is
``cred:<name>``; the stored value is either a password or a JSON object
with ``username``, ``password`` and optional ``ca``. A credential lookup
that returns nothing is a :class:`DriverError` whose text never names the
credential.

TLS
---

``address.tls: true`` enables it; ``address.tls_ca`` names the CA file (a
public path, not a secret). Per-device identity comes from the resolved
credentials, never from the manifest.

Safety
------

The device-side watchdog (``arm_watchdog``) is only armed when the
firmware has a timer the driver can reach, i.e. when ``address.arm_topic``
is set (Tasmota ``PulseTime`` and friends); else ``arm_watchdog`` returns
False and the manifest should set ``device_watchdog: false``. An MQTT
relay whose firmware has no timer stays safe through the kernel's own
dead-man and a hardware watchdog relay (design doc, section "Dead-man").
"""

from __future__ import annotations

import json
import time
from typing import Any, Callable, Optional

from .base import Driver, DriverCapabilities, DriverError, Health, Reading


def _credential(name: str) -> str:
    """Resolve one credential; an empty value reads as an error.

    Imported at call time so the module works without ``delfin.agent``.
    """
    from delfin.agent.credentials import load_credential
    return load_credential(name)


def _build_real_broker(device: Any) -> "_PahoBroker":
    """Create and configure the real paho-backed broker for ``device``."""
    return _PahoBroker(device, _resolved_credentials(device))


#: The broker seam. Tests replace this with an in-process fake.
broker_factory: Callable[..., Any] = _build_real_broker

#: The credentials seam. Tests replace this to observe resolution.
credentials: Callable[[str], str] = _credential

_QOS = 1
_CMD = "/cmd/"
_STATE = "/state/"
_AVAIL = "/availability"


def _addr(device: Any, key: str, default=None):
    return (device.address or {}).get(key, default)


def _cmd_topic(device: Any, channel: str) -> str:
    return f"{_addr(device, 'topic_base', '')}/cmd/{channel}"


def _state_topic(device: Any, channel: str) -> str:
    return f"{_addr(device, 'topic_base', '')}/state/{channel}"


def _avail_topic(device: Any) -> str:
    return f"{_addr(device, 'topic_base', '')}{_AVAIL}"


def _readable_channels(device: Any):
    names = set(device.reads)
    names.update(device.writes)     # a written channel reads back
    return names


def _resolved_credentials(device: Any) -> dict:
    """Resolve ``address.secret_ref`` to ``{username, password, tls_ca}``.

    The value is either a JSON object or a bare password; an empty or
    missing value raises a :class:`DriverError` that NEVER names the
    credential.
    """
    ref = _addr(device, "secret_ref", "") or ""
    if not ref.startswith("cred:"):
        raise DriverError("MQTT device has no secret_ref: 'cred:<name>'")
    value = credentials(ref[len("cred:"):])
    if not value:
        raise DriverError(
            "the MQTT broker credential could not be read (the store "
            "returned nothing); refusing to connect")
    try:
        parsed = json.loads(value)
        if isinstance(parsed, dict):
            return {"username": parsed.get("username", ""),
                    "password": parsed.get("password", ""),
                    "tls_ca": parsed.get("ca")}
    except ValueError:
        pass
    return {"username": "", "password": value, "tls_ca": None}


class _PahoBroker:
    """The real transport, wrapping ``paho.mqtt.client`` if importable.

    Built per device (per-device identity), connected on ``open()``.
    ``retained`` and ``available`` read what the DEVICE published; the
    command publishes are QoS 1 and never retained.
    """

    _extra = "lab-mqtt"

    def __init__(self, device: Any, resolved: dict):
        try:
            import paho.mqtt.client as mqtt_client  # noqa: PLC0415
        except ImportError as exc:                  # pragma: no cover
            raise DriverError(
                "the MQTT driver needs paho-mqtt; install DELFIN with the "
                f"{self._extra!r} extra (pip install 'delfin-complat[{self._extra}]')"
            ) from exc
        self._client = mqtt_client.Client(
            client_id=f"delfin-{device.id}", protocol=mqtt_client.MQTTv311)
        # in-flight QoS 1 allows one unacknowledged message: retries are the
        # client's buffered queue, never re-published by this driver
        self._client.max_inflight_messages_set(1)
        addr = device.address or {}
        if addr.get("tls") or addr.get("tls_ca"):
            self._client.tls_set(
                ca_certs=_addr(device, "tls_ca") or None)
        user = resolved.get("username") or ""
        self._client.username_pw_set(user, resolved.get("password") or "")
        self._client.on_connect = self._on_connect
        self._client.on_message = self._on_message
        self._client.on_disconnect = self._on_disconnect
        self._retained: dict[str, str] = {}
        self._stamp: dict[str, float] = {}
        self._connected = False
        self._connect_started = False
        self._broker_host = addr.get("host", "127.0.0.1")
        self._broker_port = int(addr.get("port", 1883))
        self._last_err: Optional[str] = None

    # -- paho callbacks -------------------------------------------------------
    def _on_connect(self, _c, _u, _f, rc, _p=None):
        self._connected = rc == 0

    def _on_disconnect(self, _c, _u, _rc):
        self._connected = False

    def _on_message(self, _c, _u, msg):
        self._retained[msg.topic] = msg.payload.decode("utf-8", "replace")
        self._stamp[msg.topic] = time.time()

    # -- the seam --------------------------------------------------------------
    def open(self):
        if self._connect_started:
            return
        self._connect_started = True
        try:
            self._client.connect(self._broker_host, self._broker_port,
                                 keepalive=60)
            self._client.loop_start()
            # QoS 1 CONNACK: give the broker a moment (timeout on the call)
            _wait(lambda: self._connected or self._last_err, timeout_s=10.0)
            if not self._connected:
                raise DriverError(
                    f"could not connect to the MQTT broker at "
                    f"{self._broker_host}:{self._broker_port}")
        except DriverError:
            raise
        except Exception as exc:
            raise DriverError(f"MQTT connect failed: {exc}") from exc

    def subscribe(self, topic: str, qos: int = _QOS):
        if self._connected:
            self._client.subscribe(topic, qos)

    def publish(self, topic: str, payload: str, qos: int = _QOS,
                retain: bool = False) -> str:
        self.open()
        try:
            info = self._client.publish(topic, payload, qos=qos,
                                        retain=retain)
            info.wait_for_publish(timeout=10.0)
        except Exception as exc:
            raise DriverError(f"MQTT publish to {topic} failed: {exc}") from exc
        return f"mqtt-{info.mid}"

    def retained(self, topic: str):
        if topic in self._retained:
            return self._retained[topic], self._stamp[topic]
        return None

    def available(self, topic: str) -> bool:
        # LWT: the retained availability payload is "online" while the device
        # is up; anything else -- especially the "offline" LWT -- means it is
        # not answering. An empty availability is also unknown (False).
        return self._retained.get(topic) == "online"

    def close(self):
        try:
            self._client.loop_stop()
            self._client.disconnect()
        except Exception:
            pass


def _wait(pred, *, timeout_s: float):
    """Spin until ``pred()`` is truthy or ``timeout_s`` elapses."""
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        if pred():
            return True
        time.sleep(0.05)
    return False


class MqttDriver(Driver):
    name = "mqtt"
    capabilities = DriverCapabilities(
        device_watchdog=True,    # only when address.arm_topic is configured
        readback=True,
        observable_commands=False,
        device_lock=False)

    def __init__(self, clock=time.time, transport: Optional[Any] = None,
                 transport_factory: Optional[Callable[..., Any]] = None):
        self.clock = clock
        self._transport_arg = transport  # an injected (fake) broker
        self._transport_factory = transport_factory or broker_factory
        self._transports: dict[str, Any] = {}

    def _transport(self, device: Any) -> Any:
        """One connected (open + subscribed) transport per device."""
        t = self._transports.get(device.id)
        if t is None:
            t = self._transport_arg or self._transport_factory(device)
            t.open()
            for ch in _readable_channels(device):
                t.subscribe(_state_topic(device, ch), _QOS)
            t.subscribe(_avail_topic(device), _QOS)
            self._transports[device.id] = t
        return t

    # -- credentials, resolved by reference only ------------------------------
    def credentials_for(self, device: Any) -> dict:
        """Resolve ``address.secret_ref`` to ``{username, password, tls_ca}``."""
        return _resolved_credentials(device)

    # -- the driver interface ---------------------------------------------------
    def read(self, device: Any) -> dict:
        t = self._transport(device)
        if not t.available(_avail_topic(device)):
            raise DriverError(
                f"{device.id} is offline (its availability/LWT topic does "
                "not report it online); its state is unknown and it is "
                "treated as failed closed")
        out: dict[str, Reading] = {}
        for ch in _readable_channels(device):
            got = t.retained(_state_topic(device, ch))
            if got is None:
                raise DriverError(
                    f"{device.id}.{ch} has not reported a state ({_state_topic(device, ch)})")
            payload, ts = got
            spec = device.channel(ch)
            out[ch] = Reading(payload, spec.unit if spec else None,
                              float(ts))
        return out

    def write(self, device: Any, channel: str, value: Any) -> str:
        t = self._transport(device)
        if not t.available(_avail_topic(device)):
            raise DriverError(
                f"{device.id} is offline; the command was not sent")
        return t.publish(_cmd_topic(device, channel), str(value), _QOS)

    def arm_watchdog(self, device: Any, seconds: float) -> bool:
        arm_topic = _addr(device, "arm_topic")
        if not arm_topic:
            return False
        # the firmware's timer command: ``arm_topic`` receives the deadline
        # in whole seconds (Tasmota ``PulseTime`` and friends)
        self._transport(device).publish(arm_topic, str(int(seconds)), _QOS)
        return True

    def disarm_watchdog(self, device: Any) -> None:
        arm_topic = _addr(device, "arm_topic")
        if arm_topic:
            self._transport(device).publish(arm_topic, "0", _QOS)

    def health(self, device: Any) -> Health:
        try:
            t = self._transport(device)
            if not t.available(_avail_topic(device)):
                return Health(False, f"{device.id} offline")
            return Health(True, "online")
        except Exception as exc:
            return Health(False, str(exc))



