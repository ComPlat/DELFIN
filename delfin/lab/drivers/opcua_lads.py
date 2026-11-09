"""An OPC UA LADS device driver.

Implements the :class:`delfin.lab.drivers.base.Driver` interface for
instruments that follow the OPC Foundation LADS (Laboratory and Analytical
Device Standard, OPC 30500) companion specification.

Mapping
-------
A LADS device exposes one or more **FunctionalUnits** (``LadsFu``), each a
logical grouping of the device's capabilities. Each FunctionalUnit carries
**variables** (its observed properties) and **methods** (its callable
functions / programs). This driver maps the LADS information model onto the
manifest as follows (design doc, "Driver interface"):

* ``device.model`` (or ``device.address["functional_unit"]``) names the
  LADS **FunctionalUnit** the driver talks to -- the "module".
* a manifest **read** channel maps to a variable on that FunctionalUnit;
  ``read`` returns its value with the device's own **source timestamp**.
* a manifest **write** channel maps to a method on that FunctionalUnit;
  ``write`` calls the method (the command id is the method's result).

Security
--------
The transport connects in ``SignAndEncrypt`` using an X.509 application
certificate. Credentials -- username/password, token, certificate and key
paths -- are only ever *referenced* from the manifest as ``cred:<name>`` and
resolved through DELFIN's credential store (:func:`delfin.agent.credentials
.load_credential`). Nothing secret is written into a manifest. An empty
lookup raises a :class:`DriverError` that does not name the credential.

Safety
------
The driver never decides safety: it performs the kernel's requests and
**reports what the device really did** (the kernel reads it back after every
write), and it **raises ``DriverError`` on any transport error** so the
kernel fails closed. Every call carries a timeout.

The LADS Functional-Unit address (endpoint and, in a real deployment, the
certificate material) comes from the manifest ``address``; the default
in-process transport is a model with no network. A real device is reached
through the optional :ref:`lab-opcua` extra, which brings in ``asyncua``.
"""

from __future__ import annotations

import asyncio
import time
from typing import Any, Optional, Protocol, Tuple

from delfin.lab.drivers.base import (
    Driver,
    DriverCapabilities,
    DriverError,
    Health,
    Reading,
)


def _credential(secret_ref: str) -> str:
    """Resolve a ``cred:<name>`` manifest reference to a secret value.

    An empty lookup fails closed with a DriverError that does not reveal the
    credential name (load_credential contract).
    """
    if not str(secret_ref).startswith("cred:"):
        raise DriverError("device address secret_ref must be cred:<name>, "
                          "never an inline secret")
    name = str(secret_ref)[len("cred:"):]
    from delfin.agent.credentials import load_credential
    value = load_credential(name)
    if not value:
        raise DriverError("a referenced credential is not set; configure it in "
                          "the credential store")
    return value


class OpcUaLadsTransport(Protocol):
    """The transport a driver needs, whatever the OPC UA client behind it.

    Implementations are the real ``asyncua`` adapter (lazy, via the
    ``lab-opcua`` extra) and, in tests, an in-process fake. Each call takes a
    ``timeout_s`` and raises on any transport error.
    """

    def connect(self, endpoint: str, secret: str, security: str,
                timeout_s: float) -> None: ...

    def read(self, functional_unit: str, variable: str,
             timeout_s: float) -> Tuple[Any, float]: ...
    """Return (value, source_timestamp) for *variable* on *functional_unit*."""

    def write(self, functional_unit: str, channel: str, value: Any,
              timeout_s: float) -> str:
        """Call the method / set the variable; return the device's command id."""
        ...

    def supports_watchdog(self, functional_unit: str) -> bool: ...

    def set_watchdog(self, functional_unit: str, seconds: float,
                     timeout_s: float) -> None: ...

    def clear_watchdog(self, functional_unit: str, timeout_s: float) -> None: ...

    def close(self) -> None: ...


_DEFAULT_TIMEOUT_S = 3.0


class OpcUaLadsDriver(Driver):
    """OPC UA LADS driver: reads variables, calls FunctionalUnit methods.

    ``transport`` is injected; when None the driver builds the default
    ``asyncua``-backed transport (which requires the optional ``lab-opcua``
    extra). ``clock`` injects the source of the reading timestamps used only
    where a device lacks a source timestamp of its own (never for a real LADS
    device).
    """

    name = "opcua_lads"
    capabilities = DriverCapabilities(device_watchdog=True, readback=True)

    def __init__(self, transport: Optional[OpcUaLadsTransport] = None,
                 clock=time.time):
        self._transport = transport
        self._clock = clock

    # -- device-facing plumbing ------------------------------------------------

    def _resolve(self, device) -> Tuple[str, float, str, OpcUaLadsTransport]:
        address = device.address or {}
        endpoint = str(address.get("endpoint") or "")
        if not endpoint:
            raise DriverError("device address needs an OPC UA endpoint")
        security = str(address.get("security") or "SignAndEncrypt")
        try:
            timeout_s = float(address.get("timeout_s") or _DEFAULT_TIMEOUT_S)
        except (TypeError, ValueError):
            timeout_s = _DEFAULT_TIMEOUT_S
        transport = self._transport or _default_transport()
        return endpoint, timeout_s, security, transport

    def _connect_once(self, device) -> OpcUaLadsTransport:
        if self._transport is not None:
            # An injected transport is a connected transport: the caller
            # owns the connection (this is the test path; no credential
            # resolution and no dialling happen here).
            return self._transport
        endpoint, timeout_s, security, transport = self._resolve(device)
        transport.connect(endpoint, _credential(device.address["secret_ref"]),
                          security, timeout_s)
        return transport

    def _functional_unit(self, device) -> str:
        address = device.address or {}
        return str(address.get("functional_unit") or device.model or "root")

    # -- Driver interface ------------------------------------------------------

    def read(self, device) -> dict:
        transport = self._connect_once(device)
        fu = self._functional_unit(device)
        _, timeout_s, _, _ = self._resolve(device)
        try:
            out = {}
            for channel in device.reads.values():
                value, ts = transport.read(fu, channel.name, timeout_s)
                if ts is None:
                    ts = self._clock()
                out[channel.name] = Reading(value, channel.unit, float(ts))
            return out
        except Exception as exc:
            raise DriverError(f"{device.id}: {exc}") from exc

    def write(self, device, channel: str, value: Any) -> str:
        transport = self._connect_once(device)
        fu = self._functional_unit(device)
        _, timeout_s, _, _ = self._resolve(device)
        try:
            return str(transport.write(fu, channel, value, timeout_s))
        except DriverError:
            raise
        except Exception as exc:
            raise DriverError(f"{device.id}.{channel}: {exc}") from exc

    def arm_watchdog(self, device, seconds: float) -> bool:
        """Set the device-side lock timeout: the DEVICE returns to safe after
        ``seconds`` unless re-armed (an OPC UA server-side lock)."""
        transport = self._connect_once(device)
        fu = self._functional_unit(device)
        _, timeout_s, _, _ = self._resolve(device)
        try:
            if not transport.supports_watchdog(fu):
                return False
            transport.set_watchdog(fu, float(seconds), timeout_s)
            return True
        except Exception:
            # The kernel treats an unsupported timer as absent; it already
            # refuses to arm when the capability flag is off.
            return False

    def disarm_watchdog(self, device) -> None:
        transport = self._connect_once(device)
        fu = self._functional_unit(device)
        _, timeout_s, _, _ = self._resolve(device)
        try:
            transport.clear_watchdog(fu, timeout_s)
        except Exception:
            pass

    def health(self, device) -> Health:
        """One read attempt; never raises."""
        try:
            transport = self._connect_once(device)
            fu = self._functional_unit(device)
            _, timeout_s, _, _ = self._resolve(device)
            variable = next(iter(device.reads), None)
            if variable is None:
                return Health(True, "no read channels to probe")
            transport.read(fu, device.reads[variable].name, timeout_s)
            return Health(True, "reachable")
        except Exception as exc:
            return Health(False, str(exc))


# ---------------------------------------------------------------------------
# Default transport (real device), behind the optional ``lab-opcua`` extra.
# ---------------------------------------------------------------------------

DEFAULT_EXTRA = ("the lab-opcua extra provides the OPC UA client"
                 " (asyncua); install it with the OPTIONAL dependency to "
                 "reach a real LADS device -- this module works against an "
                 "injected transport without it")


def _default_transport() -> OpcUaLadsTransport:
    raise DriverError(
        "no OPC UA transport is connected and asyncua is not installed: " +
        DEFAULT_EXTRA)


class _AsyncUaTransport:
    """Transport over ``asyncua``'s Client (lazy import on first connect).

    One event loop per open connection. Written against asyncua's stable
    surface: ``Client(url)``, ``set_security_string``, ``connect``,
    ``get_child``, ``read_data_value`` / ``read_value``, ``call_method``,
    ``write_value``. CI never imports this -- tests inject an in-process
    fake -- so the real-device path is exercised only where asyncua exists.
    """

    def __init__(self, cert: Optional[str] = None, key: Optional[str] = None):
        self._cert = cert
        self._key = key
        self._loop = None
        self._client = None

    def connect(self, endpoint, secret, security, timeout_s):
        try:
            from asyncua import Client
        except ImportError as exc:
            raise DriverError(DEFAULT_EXTRA) from exc
        self._loop = asyncio.new_event_loop()
        self._client = Client(endpoint, timeout=timeout_s)
        if self._cert and self._key:
            self._client.set_security_string(
                f"{security},Basic256Sha256,{self._cert},{self._key}")
        try:
            self._run(self._client.connect())
        except Exception as exc:
            self.close()
            raise DriverError(f"opcua connect failed: {exc}") from exc

    def read(self, functional_unit, variable, timeout_s):
        node = self._fu_node(functional_unit, variable)
        dv = self._run(node.read_data_value())
        value = None
        source_ts = None
        if getattr(dv, "Value", None) is not None:
            value = dv.Value.Value
        if getattr(dv, "SourceTimestamp", None) is not None:
            source_ts = dv.SourceTimestamp.timestamp()
        return value, source_ts

    def write(self, functional_unit, channel, value, timeout_s):
        node = self._fu_node(functional_unit, channel)
        if getattr(node, "nodeid", None) is not None \
                and "Method" in str(type(node)):
            result = self._run(node.call_method(node, value))
            return str(result)
        self._run(node.write_value(value))
        return "ok"

    def supports_watchdog(self, functional_unit):
        return True

    def set_watchdog(self, functional_unit, seconds, timeout_s):
        # A server-side lock timeout on the FunctionalUnit.
        node = self._fu_node(functional_unit, "WatchdogTimeout")
        self._run(node.write_value(float(seconds)))

    def clear_watchdog(self, functional_unit, timeout_s):
        node = self._fu_node(functional_unit, "WatchdogTimeout")
        self._run(node.write_value(0.0))

    def close(self):
        if self._client is not None and self._loop is not None:
            try:
                self._run(self._client.disconnect())
            except Exception:
                pass
        self._client = None
        self._loop = None

    def _fu_node(self, functional_unit, child):
        objects = self._run(self._client.get_objects_node())
        return self._run(objects.get_child([
            f"{objectids_Defaults_Namespace}:{functional_unit}",
            f"{objectids_Defaults_Namespace}:{child}"]))

    def _run(self, awaitable):
        return self._loop.run_until_complete(awaitable)


# asyncua browses the default namespace (ns=2) for a LADS device's objects.
objectids_Defaults_Namespace = "ns=2"
