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
                timeout_s: float, namespace: Optional[str] = None) -> None: ...

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
        ref = (device.address or {}).get("secret_ref")
        if not ref:
            raise DriverError(f"{device.id}: address needs secret_ref "
                              "(the client certificate, key and the "
                              "device's certificate)")
        transport.connect(endpoint, _credential(ref), security, timeout_s,
                          namespace=(device.address or {}).get("namespace"))
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
                if ts is None or channel.name in device.writes:
                    # A setpoint or state the device just reported is
                    # current although its source timestamp is the time it
                    # was last changed; only a pure measurement keeps the
                    # source timestamp, so a frozen sensor still reads stale.
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

#: Security modes a lab device is reached with. ``None`` (no signature, no
#: encryption) is refused: a setpoint must not travel unauthenticated.
SECURITY_MODES = ("SignAndEncrypt", "Sign")
_POLICY = "Basic256Sha256"


def _default_transport() -> OpcUaLadsTransport:
    return _AsyncUaTransport()


def _secret_fields(secret: str) -> dict:
    """The credential's JSON: ``cert``, ``key`` (the client's application
    certificate and private key), ``server_cert`` (the device's certificate,
    pinned) and optionally ``key_password``, ``username``, ``password``."""
    import json
    try:
        doc = json.loads(secret)
    except ValueError:
        doc = None
    if not isinstance(doc, dict):
        raise DriverError("the OPC UA credential must be a JSON object with "
                          "cert, key and server_cert")
    missing = [k for k in ("cert", "key", "server_cert") if not doc.get(k)]
    if missing:
        raise DriverError("the OPC UA credential lacks " + ", ".join(missing)
                          + "; a device is reached signed and encrypted, "
                          "with its certificate pinned")
    return doc


def _application_uri(cert_path: str) -> Optional[str]:
    """The URI in the client certificate's subjectAltName: the server
    refuses a session whose application URI differs from it."""
    try:
        from cryptography import x509
        data = open(cert_path, "rb").read()
        cert = (x509.load_pem_x509_certificate(data) if data.lstrip()
                .startswith(b"-----") else x509.load_der_x509_certificate(data))
        san = cert.extensions.get_extension_for_class(
            x509.SubjectAlternativeName).value
        uris = san.get_values_for_type(x509.UniformResourceIdentifier)
        return uris[0] if uris else None
    except Exception:
        return None


class _AsyncUaTransport:
    """Transport over ``asyncua``'s Client (lazy import on first use).

    ``connect`` checks and keeps the connection settings; every operation
    then opens one session, does its work and closes it (``asyncio.run`` in
    a worker thread, so it works whether or not the calling thread runs an
    event loop). Nothing stays open between kernel calls: the kernel loads
    a new driver for every call and never closes one.

    Security: ``Basic256Sha256`` with ``Sign`` or ``SignAndEncrypt``, the
    client's own certificate and key, and the device's certificate pinned
    (a server presenting another certificate fails the handshake). The
    application URI sent is the one in the client certificate. Browse paths
    are ``<ns>:<FunctionalUnit>/<ns>:<child>`` with ``ns`` the index of
    ``address.namespace`` (a URI), or 2 without one.
    """

    def __init__(self):
        self._cfg: Optional[dict] = None

    def connect(self, endpoint, secret, security, timeout_s, namespace=None):
        if security not in SECURITY_MODES:
            raise DriverError(
                f"OPC UA security {security!r} is refused; use one of "
                f"{', '.join(SECURITY_MODES)}")
        fields = _secret_fields(secret)
        try:
            import asyncua  # noqa: F401
        except ImportError as exc:
            raise DriverError(DEFAULT_EXTRA) from exc
        self._cfg = {"endpoint": endpoint, "security": security,
                     "timeout_s": float(timeout_s), "namespace": namespace,
                     **fields}

    def _session(self, op):
        cfg = self._cfg
        if cfg is None:
            raise DriverError("OPC UA transport used before connect")
        result: dict = {}

        async def run():
            from asyncua import Client
            client = Client(cfg["endpoint"], timeout=cfg["timeout_s"])
            uri = _application_uri(cfg["cert"])
            if uri:
                client.application_uri = uri
            key = cfg["key"]
            if cfg.get("key_password"):
                key = f"{key}::{cfg['key_password']}"
            await client.set_security_string(
                f"{_POLICY},{cfg['security']},{cfg['cert']},{key},"
                f"{cfg['server_cert']}")
            if cfg.get("username"):
                client.set_user(cfg["username"])
                client.set_password(cfg.get("password") or "")
            async with client:
                ns_uri = cfg.get("namespace")
                ns = (await client.get_namespace_index(ns_uri)
                      if ns_uri else 2)
                return await op(client, ns)

        def worker():
            try:
                result["value"] = asyncio.run(asyncio.wait_for(
                    run(), timeout=cfg["timeout_s"] * 4))
            except BaseException as exc:          # noqa: BLE001
                result["error"] = exc

        import threading
        t = threading.Thread(target=worker, daemon=True,
                             name="delfin-opcua-op")
        t.start()
        t.join(cfg["timeout_s"] * 4 + 5)
        if t.is_alive():
            raise DriverError("OPC UA operation timed out")
        if "error" in result:
            exc = result["error"]
            raise DriverError(f"OPC UA: {type(exc).__name__}: {exc}") from None
        return result.get("value")

    @staticmethod
    async def _child(client, ns, functional_unit, child):
        return await client.nodes.objects.get_child(
            [f"{ns}:{functional_unit}", f"{ns}:{child}"])

    def read(self, functional_unit, variable, timeout_s):
        async def op(client, ns):
            node = await self._child(client, ns, functional_unit, variable)
            dv = await node.read_data_value()
            ts = dv.SourceTimestamp or dv.ServerTimestamp
            return dv.Value.Value, (ts.timestamp() if ts else None)
        return self._session(op)

    def write(self, functional_unit, channel, value, timeout_s):
        async def op(client, ns):
            from asyncua import ua
            node = await self._child(client, ns, functional_unit, channel)
            if await node.read_node_class() == ua.NodeClass.Method:
                parent = await client.nodes.objects.get_child(
                    [f"{ns}:{functional_unit}"])
                return str(await parent.call_method(node, value))
            vtype = await node.read_data_type_as_variant_type()
            if vtype in (ua.VariantType.Double, ua.VariantType.Float):
                value_ = float(value)
            elif vtype == ua.VariantType.String:
                value_ = str(value)
            else:
                value_ = value
            await node.write_value(ua.DataValue(ua.Variant(value_, vtype)))
            return f"opcua-write-{functional_unit}.{channel}"
        return self._session(op)

    def supports_watchdog(self, functional_unit):
        async def op(client, ns):
            try:
                await self._child(client, ns, functional_unit,
                                  "WatchdogTimeout")
                return True
            except Exception:
                return False
        return bool(self._session(op))

    def set_watchdog(self, functional_unit, seconds, timeout_s):
        self.write(functional_unit, "WatchdogTimeout", float(seconds),
                   timeout_s)

    def clear_watchdog(self, functional_unit, timeout_s):
        self.write(functional_unit, "WatchdogTimeout", 0.0, timeout_s)

    def close(self):
        self._cfg = None
