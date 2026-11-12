"""Shelly Gen2 (RPC over HTTP) and Tasmota (HTTP) smart plugs and relays.

Both families are LAN-only devices: the driver speaks plain HTTP on the
local network (``address.host``/``address.port``) and never contacts a
cloud service. Transport is the standard library (``urllib.request``);
no third-party HTTP client is imported (see scratch/LICENSE_NOTES.md).

Shelly Gen2 (driver name ``shelly_rpc``)
    Every call is one JSON-RPC frame ``{"id", "src", "method", "params"}``
    POSTed to ``/rpc``; the answer carries ``result`` or ``error``
    (``{"code", "message"}``), both with HTTP 200. ``address.slot`` is the
    switch's ``id`` on the device (0 on a single plug).
    ``write``      ``Switch.Set {id, on}``; ``on`` is required by the device
    ``read``       ``Switch.GetStatus {id}`` -> ``output`` (on/off) and
                   ``apower`` (W); the readback channel is ``output``
    auth           HTTP digest as the device sends it: SHA-256 (Gen2
                   firmware) or MD5, user ``admin``, a fresh cnonce per
                   request
Tasmota (driver name ``tasmota``)
    ``write``      GET ``/cm?cmnd=Power On|Off``
    ``read``       GET ``/cm?cmnd=Power`` -> ``{"POWER": "ON"}``; the
                   readback channel is ``power``
    auth           Tasmota's web password travels as the ``user`` and
                   ``password`` query parameters of ``/cm`` (the device
                   answers HTTP 200 with a ``WARNING`` without them, not 401);
                   they never appear in an error text or a log line
Credentials (both)
    By reference only: ``address.secret_ref: cred:<name>`` names a value in
    DELFIN's secret store, either ``user:password`` or a bare password (user
    ``admin``, the fixed web user of both firmwares). An error never names
    the reference.
Device-side timer
    Shelly: ``toggle_after`` FLIPS the output after the delay, so it is only
    ever sent together with ``on: true``: on a switch that is off it would
    switch the load on. ``arm_watchdog`` on an output that is on re-sends
    ``Switch.Set {id, on: true, toggle_after}`` (the device restarts its
    timer); on an output that is off it only records the delay, and the next
    ``write(power, on)`` carries it. ``disarm_watchdog`` re-sends ``on: true``
    without the timer when the output is on.
    Tasmota: ``PulseTime`` -- 1..111 in 0.1 s units, 112..64900 as
    ``seconds + 100``; setting it restarts the running timer, at expiry the
    relay goes off, 0 clears it.

Safety stays in the kernel: this module only moves a device inside the
manifest envelope and reports what the device really did (readback). Every
transport failure -- connection refused, DNS, timeout, non-2xx answer, an
RPC error, a refused password -- raises ``DriverError`` so the kernel fails
closed; ``health`` is the one method that never raises. Every call carries
a timeout.
"""

from __future__ import annotations

import base64
import hashlib
import itertools
import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Optional

from .base import Driver, DriverCapabilities, DriverError, Health, Reading

DEFAULT_TIMEOUT_S = 5.0

#: The manifest address key that names the Shelly switch id ("slot").
SLOT_KEY = "slot"

#: The fixed web user of Shelly Gen2 and Tasmota.
DEFAULT_USER = "admin"

_SRC = "delfin"
_RPC_IDS = itertools.count(1)


def _base_url(device) -> str:
    addr = dict(device.address or {})
    host = addr.get("host")
    if not host:
        raise DriverError(f"{device.id}: address has no host")
    port = addr.get("port")
    netloc = f"{host}:{port}" if port else host
    return f"http://{netloc}"


def _switch_id(device) -> int:
    return int((device.address or {}).get(SLOT_KEY, 0))


def _auth_ref(device) -> Optional[tuple[str, str]]:
    """Resolve ``secret_ref: cred:<name>`` from the DELFIN secret store.

    ``user:password`` or a bare password (user ``admin``). No message
    names the reference.
    """
    ref = (device.address or {}).get("secret_ref")
    if not ref:
        return None
    scheme, _, name = str(ref).partition(":")
    if scheme != "cred" or not name:
        raise DriverError(
            f"{device.id}: secret_ref must be 'cred:<name>'")
    from delfin.agent import credentials as _cred
    try:
        value = _cred._stored_value(name)
    except Exception as exc:
        raise DriverError(f"{device.id}: the device credential could not be "
                          f"read ({type(exc).__name__})") from None
    if not value:
        raise DriverError(f"{device.id}: the device credential is empty; "
                          "refusing to connect")
    user, sep, pw = str(value).partition(":")
    if not sep:
        return (DEFAULT_USER, str(value))
    return (user or DEFAULT_USER, pw)


def _basic_auth(device, user: str, pw: str) -> str:
    token = base64.b64encode(f"{user}:{pw}".encode()).decode()
    return f"Basic {token}"


def _challenge_fields(challenge: str) -> dict[str, str]:
    """The ``key=value`` fields of a WWW-Authenticate challenge (quoted
    values may contain commas)."""
    parts: dict[str, str] = {}
    rest = challenge.split(None, 1)[1] if " " in challenge.strip() else ""
    for _ in range(64):
        k, sep, v = rest.partition("=")
        k = k.strip().strip(",").strip()
        if not sep or not k:
            break
        v = v.lstrip()
        if v.startswith('"'):
            v = v[1:]
            i = v.find('"')
            part, rest = (v[:i], v[i + 1:]) if i >= 0 else (v, "")
        else:
            part, _, rest = v.partition(",")
        parts[k.lower()] = part.strip()
        rest = rest.lstrip(", ")
    return parts


_DIGEST_HASHES = {"md5": hashlib.md5, "sha-256": hashlib.sha256}


def _digest_auth(user: str, pw: str, method: str, uri: str,
                 challenge: str) -> str:
    """RFC 7616 digest response to a WWW-Authenticate challenge, with the
    hash the challenge names (``algorithm=SHA-256`` on Shelly Gen2, MD5
    when absent) and a fresh cnonce, so one authenticated request is sent
    (no reliance on an HTTP client's 401-replay internals)."""
    f = _challenge_fields(challenge)
    algorithm = f.get("algorithm", "MD5")
    hfun = _DIGEST_HASHES.get(algorithm.lower())
    if hfun is None:
        raise DriverError(f"unsupported digest algorithm {algorithm!r}")

    def H(text: str) -> str:
        return hfun(text.encode()).hexdigest()

    realm, nonce = f.get("realm", ""), f.get("nonce", "")
    qop = f.get("qop", "")
    qop = "auth" if "auth" in [q.strip() for q in qop.split(",")] else ""
    ha1 = H(f"{user}:{realm}:{pw}")
    ha2 = H(f"{method}:{uri}")
    cnonce = os.urandom(8).hex()
    nc = "00000001"
    if qop:
        response = H(f"{ha1}:{nonce}:{nc}:{cnonce}:{qop}:{ha2}")
    else:
        response = H(f"{ha1}:{nonce}:{ha2}")
    hdr = (f'Digest username="{user}", realm="{realm}", nonce="{nonce}", '
           f'uri="{uri}", algorithm={algorithm}, response="{response}"')
    if qop:
        hdr += f', qop={qop}, nc={nc}, cnonce="{cnonce}"'
    if "opaque" in f:
        hdr += f', opaque="{f["opaque"]}"'
    return hdr


def _http_transport(device, method: str, path: str, body: Optional[bytes],
                    headers: dict) -> tuple[int, bytes, dict]:
    """Real HTTP transport. LAN-only: routes around the environment's egress
    proxy (the device is on the LAN and must be reached directly) and around
    any HTTP client's 401-replay (the caller computes the Authorization
    header itself). Raises ``DriverError`` on any transport error."""
    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({}),
        urllib.request.HTTPHandler(),
        urllib.request.HTTPSHandler(),
    )
    url = f"{_base_url(device)}{path}"
    req = urllib.request.Request(url, data=body, method=method, headers=headers)
    try:
        with opener.open(req, timeout=DEFAULT_TIMEOUT_S) as resp:
            return resp.status, resp.read(), dict(resp.headers.items())
    except urllib.error.HTTPError as exc:
        hdrs = dict(exc.headers.items()) if exc.headers else {}
        return exc.code, b"", hdrs
    except urllib.error.URLError as exc:
        raise DriverError(f"{device.id}: cannot reach device: {exc.reason}")
    except TimeoutError:
        raise DriverError(f"{device.id}: device timed out after "
                          f"{DEFAULT_TIMEOUT_S}s")
    except OSError as exc:
        raise DriverError(f"{device.id}: transport error: {exc}")


#: Injectable transport. Tests replace this with an in-process fake so the
#: driver is exercised end to end without any real device or socket; the
#: production value performs real LAN HTTP.
TRANSPORT = _http_transport


def _request(device, method: str, path: str,
             body: Optional[bytes] = None, *,
             shown: Optional[str] = None) -> tuple[int, dict]:
    """One HTTP request with explicit auth; raises ``DriverError`` on any
    transport error (kernel fails closed). Returns ``(status, parsed body)``.
    A 401 challenge is answered once with the computed Authorization header.
    *shown* is the path as an error may print it (without a secret query).
    """
    shown = shown or path
    user = pw = None
    if (device.address or {}).get("secret_ref"):
        user, pw = _auth_ref(device) or (None, None)
    headers = {"Content-Type": "application/json"} if body is not None else {}

    def _send() -> tuple[int, bytes, dict]:
        return TRANSPORT(device, method, path, body, headers)

    status, raw, resp_headers = _send()
    www_auth = {k.lower(): v for k, v in resp_headers.items()}.get(
        "www-authenticate")
    if status == 401 and user is not None and www_auth:
        challenge = www_auth.lstrip()
        if challenge.lower().startswith("basic"):
            headers["Authorization"] = _basic_auth(device, user, pw)
        elif challenge.lower().startswith("digest"):
            headers["Authorization"] = _digest_auth(
                user, pw, method, path, challenge)
        else:
            raise DriverError(f"{device.id}: unsupported auth scheme")
        # retry once with the computed Authorization header; the header dict is
        # the same object, so TRANSPORT sees it
        status, raw, resp_headers = _send()
    if status == 401:
        raise DriverError(f"{device.id}: the device refused the credentials "
                          f"(HTTP 401 from {_base_url(device)}{shown})")
    if status >= 400 and not raw:
        raise DriverError(f"{device.id}: HTTP {status} from "
                          f"{_base_url(device)}{shown}")
    try:
        return status, json.loads(raw.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise DriverError(f"{device.id}: device sent no JSON: {exc}")


def _rpc(device, method: str, params: dict) -> Any:
    """One Shelly RPC call as a JSON-RPC frame on ``/rpc``; returns its
    ``result``. An ``error`` answer raises ``DriverError`` with the
    device's code and message."""
    body = json.dumps({"id": next(_RPC_IDS), "src": _SRC, "method": method,
                       "params": params}).encode()
    _status, doc = _request(device, "POST", "/rpc", body)
    if not isinstance(doc, dict):
        raise DriverError(f"{device.id}: {method}: no RPC frame in the answer")
    err = doc.get("error")
    if err is not None:
        if isinstance(err, dict):
            raise DriverError(f"{device.id}: {method} refused: "
                              f"{err.get('message')} ({err.get('code')})")
        raise DriverError(f"{device.id}: {method} refused: {err}")
    if "result" not in doc:
        raise DriverError(f"{device.id}: {method}: answer has no result")
    return doc["result"]


def _get_cmnd(device, cmnd: str) -> dict:
    """One Tasmota HTTP command. Raises ``DriverError`` likewise; a
    password the device refuses (its ``WARNING`` answer) is an error."""
    query = {"cmnd": cmnd}
    shown = f"/cm?{urllib.parse.urlencode(query)}"
    cred = _auth_ref(device) if (device.address or {}).get("secret_ref") else None
    if cred is not None:
        query = {"user": cred[0], "password": cred[1], "cmnd": cmnd}
    path = f"/cm?{urllib.parse.urlencode(query)}"
    _status, doc = _request(device, "GET", path, shown=shown)
    if not isinstance(doc, dict):
        raise DriverError(f"{device.id}: {cmnd}: unexpected answer")
    if "WARNING" in doc:
        raise DriverError(f"{device.id}: the device refused the command "
                          f"({doc['WARNING']}); check its credential")
    return doc


def _pulsetime_value(seconds: float) -> int:
    """Tasmota PulseTime encoding: 0.1 s units up to 11.1 s, then s + 100."""
    s = float(seconds)
    if s * 10 <= 111:
        return max(1, int(round(s * 10)))
    return min(64900, max(112, int(round(s)) + 100))


def _out_state(flag: Any, spec) -> str:
    on = "on" if flag else "off"
    if spec and spec.states and on not in spec.states:
        raise DriverError(f"device reported {flag!r}, not a legal state")
    return on


class ShellyDriver(Driver):
    name = "shelly_rpc"
    capabilities = DriverCapabilities(device_watchdog=True, readback=True)

    def __init__(self, clock=time.time):
        self.clock = clock
        #: device id -> toggle_after seconds for the next switch-on
        self._timer: dict[str, float] = {}

    # -- Shelly Gen2 RPC ---------------------------------------------------
    def _status(self, device) -> dict:
        result = _rpc(device, "Switch.GetStatus", {"id": _switch_id(device)})
        if not isinstance(result, dict):
            raise DriverError(f"{device.id}: Switch.GetStatus returned no result")
        return result

    def _set(self, device, on: bool, toggle_after: Optional[float] = None):
        params: dict = {"id": _switch_id(device), "on": bool(on)}
        if on and toggle_after:
            params["toggle_after"] = float(toggle_after)
        return _rpc(device, "Switch.Set", params)

    def read(self, device) -> dict:
        result = self._status(device)
        t = self.clock()
        out = {}
        if "output" in result:
            out["output"] = Reading(
                _out_state(bool(result["output"]),
                           device.channel("output")),
                None, t)
        if "apower" in result:
            out["apower"] = Reading(float(result["apower"]), "W", t)
        if not out:
            raise DriverError(f"{device.id}: status had no usable fields")
        return out

    def write(self, device, channel: str, value: Any) -> str:
        spec = device.channel(channel)
        if spec is None:
            raise DriverError(f"{device.id} has no channel {channel}")
        if spec.states and value not in spec.states:
            raise DriverError(
                f"{device.id}.{channel}: {value!r} is not one of {spec.states}")
        if channel != "power":
            raise DriverError(
                f"{device.id}.{channel}: Shelly switch writes on/off only")
        on = value == "on"
        self._set(device, on, self._timer.get(device.id) if on else None)
        return f"shelly-{device.id}-{self.clock():.0f}"

    def arm_watchdog(self, device, seconds: float) -> bool:
        if float(seconds) <= 0:
            return False
        self._timer[device.id] = float(seconds)
        if self._status(device).get("output"):
            # restart the device's flip-back timer on the output that is on
            self._set(device, True, float(seconds))
        return True

    def disarm_watchdog(self, device) -> None:
        self._timer.pop(device.id, None)
        if self._status(device).get("output"):
            self._set(device, True)

    def health(self, device) -> Health:
        try:
            self.read(device)
            return Health(True, "Switch.GetStatus answered")
        except Exception as exc:
            return Health(False, str(exc))


class TasmotaDriver(Driver):
    name = "tasmota"
    capabilities = DriverCapabilities(device_watchdog=True, readback=True)

    def __init__(self, clock=time.time):
        self.clock = clock

    def read(self, device) -> dict:
        doc = _get_cmnd(device, "Power")
        state = doc.get("POWER")
        if state is None:
            raise DriverError(f"{device.id}: Power query returned no POWER")
        t = self.clock()
        return {"power": Reading(
            _out_state(str(state).lower() == "on", device.channel("power")),
            None, t)}

    def write(self, device, channel: str, value: Any) -> str:
        spec = device.channel(channel)
        if spec is None:
            raise DriverError(f"{device.id} has no channel {channel}")
        if channel != "power":
            raise DriverError(f"{device.id}.{channel}: Tasmota writes "
                              "Power on/off only")
        if value not in (spec.states or ["on", "off"]):
            raise DriverError(
                f"{device.id}.{channel}: {value!r} is not one of {spec.states}")
        doc = _get_cmnd(device, f"Power {'On' if value == 'on' else 'Off'}")
        if doc.get("POWER") is None:
            raise DriverError(f"{device.id}: Power command refused")
        return f"tasmota-{device.id}-{self.clock():.0f}"

    def arm_watchdog(self, device, seconds: float) -> bool:
        if float(seconds) <= 0:
            return False
        _get_cmnd(device, f"PulseTime {_pulsetime_value(seconds)}")
        return True

    def disarm_watchdog(self, device) -> None:
        _get_cmnd(device, "PulseTime 0")

    def health(self, device) -> Health:
        try:
            self.read(device)
            return Health(True, "Power query answered")
        except Exception as exc:
            return Health(False, str(exc))
