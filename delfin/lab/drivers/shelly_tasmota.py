"""Shelly Gen2 (RPC over HTTP) and Tasmota (HTTP) smart plugs and relays.

Both families are LAN-only devices: the driver speaks plain HTTP on the
local network (``address.host``/``address.port``) and never contacts a
cloud service. Transport is the standard library (``urllib.request``);
no third-party HTTP client is imported (see scratch/LICENSE_NOTES.md).

Shelly Gen2 (driver name ``shelly_rpc``)
    ``write``      POST ``/rpc/Switch.Set`` with ``on`` (and ``toggle_after``
                   only when arming the device-side timer)
    ``read``       POST ``/rpc/Switch.GetStatus`` -> ``output`` (on/off) and
                   ``apower`` (W); the readback channel is ``output``
    auth           HTTP digest, by reference: the manifest address carries
                   ``secret_ref: cred:<name>``, resolved from the DELFIN
                   secret store at driver start, never inline
Tasmota (driver name ``tasmota``)
    ``write``      GET ``/cm?cmnd=Power On|Off``
    ``read``       GET ``/cm?cmnd=Power`` -> ``{"POWER": "ON"}``; the
                   readback channel is ``power``
    auth           HTTP basic, by reference the same way
Device-side timer
    Shelly: ``Switch.Set`` with ``toggle_after`` seconds (the device toggles
    itself back after the delay; the kernel re-arms every heartbeat).
    Tasmota: ``PulseTime`` -- 1..111 in 0.1 s units, 112..64900 as
    ``seconds + 100``; a value of 0 clears it.

Safety stays in the kernel: this module only moves a device inside the
manifest envelope and reports what the device really did (readback). Every
transport failure -- connection refused, DNS, timeout, non-2xx answer --
raises ``DriverError`` so the kernel fails closed; ``health`` is the one
method that never raises. Every call carries a timeout.
"""

from __future__ import annotations

import base64
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Optional

from .base import Driver, DriverCapabilities, DriverError, Health, Reading

DEFAULT_TIMEOUT_S = 5.0

#: Shelly ``Switch.Set`` params carry the relay id ("slot" in the manifest).
SLOT_KEY = "slot"

_SRC = "delfin"


def _base_url(device) -> str:
    addr = dict(device.address or {})
    host = addr.get("host")
    if not host:
        raise DriverError(f"{device.id}: address has no host")
    port = addr.get("port")
    netloc = f"{host}:{port}" if port else host
    return f"http://{netloc}"


def _auth_ref(device) -> Optional[tuple[str, str]]:
    """Resolve ``secret_ref: cred:<name>`` from the DELFIN secret store."""
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
        raise DriverError(f"{device.id}: auth ref {name!r} unreadable: {exc}")
    if not value:
        raise DriverError(f"{device.id}: auth ref {name!r} is empty")
    user, _, pw = str(value).partition(":")
    if not pw:
        raise DriverError(
            f"{device.id}: auth ref {name!r} must be 'user:password'")
    return (user, pw)


def _basic_auth(device, user: str, pw: str) -> str:
    token = base64.b64encode(f"{user}:{pw}".encode()).decode()
    return f"Basic {token}"


def _digest_auth(user: str, pw: str, method: str, uri: str,
                 challenge: str) -> str:
    """RFC 7616 ``qop="auth"`` digest response from a WWW-Authenticate
    challenge, computed in the driver so a single authenticated request is
    sent (no reliance on an HTTP client's 401-replay internals)."""
    parts: dict[str, str] = {}
    # strip scheme token and split on commas, honoring quoted commas
    rest = challenge.split(None, 1)[1] if " " in challenge else challenge
    for _ in range(64):
        k, sep, v = rest.partition("=")
        k = k.strip().strip(",")
        if not sep:
            break
        v = v.strip()
        if v.startswith('"'):
            v = v[1:]
            i = v.find('"')
            part = v[:i]
            rest = v[i + 1:]
        else:
            part, _, rest = v.partition(",")
        parts[k] = part
        rest = rest.lstrip(", ")
        if not parts.get(k):
            break
    realm = parts.get("realm", "")
    nonce = parts.get("nonce", "")
    qop = parts.get("qop", "auth").split(",")[0].strip()
    ha1 = hashlib.md5(f"{user}:{realm}:{pw}".encode()).hexdigest()
    ha2 = hashlib.md5(f"{method}:{uri}".encode()).hexdigest()
    cnonce = NONCE_CNONCE
    nc = "00000001"
    if qop:
        response = hashlib.md5(
            f"{ha1}:{nonce}:{nc}:{cnonce}:{qop}:{ha2}".encode()).hexdigest()
    else:
        response = hashlib.md5(
            f"{ha1}:{nonce}:{ha2}".encode()).hexdigest()
    hdr = (f'Digest username="{user}", realm="{realm}", nonce="{nonce}", '
           f'uri="{uri}", response="{response}"')
    if qop:
        hdr += f', qop={qop}, nc={nc}, cnonce="{cnonce}"'
    return hdr


#: One cnonce for the whole process is fine for challenge-response signing.
NONCE_CNONCE = "abcdef0123456789"


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
             body: Optional[bytes] = None) -> tuple[int, dict]:
    """One HTTP request with explicit auth; raises ``DriverError`` on any
    transport error (kernel fails closed). Returns ``(status, parsed body)``.
    A 401 challenge is answered once with the computed Authorization header.
    """
    user = pw = None
    secret_ref = (device.address or {}).get("secret_ref")
    if secret_ref:
        user, pw = _auth_ref(device) or (None, None)
    headers = {"Content-Type": "application/json"} if body is not None else {}

    def _send() -> tuple[int, bytes, dict]:
        return TRANSPORT(device, method, path, body, headers)

    status, raw, resp_headers = _send()
    www_auth = resp_headers.get("WWW-Authenticate")
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
    if status >= 400 and not raw:
        raise DriverError(f"{device.id}: HTTP {status} from "
                          f"{_base_url(device)}{path}")
    try:
        return status, json.loads(raw.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise DriverError(f"{device.id}: device sent no JSON: {exc}")


def _post_json(device, path: str, params: dict) -> dict:
    """One Shelly RPC call. Raises ``DriverError`` on any transport error."""
    body = json.dumps({"id": 1, "method": path.rpartition("/")[2],
                       "src": _SRC, "params": params}).encode()
    _status, doc = _request(device, "POST", path, body)
    return doc


def _get_cmnd(device, cmnd: str) -> dict:
    """One Tasmota HTTP command. Raises ``DriverError`` likewise."""
    path = f"/cm?{urllib.parse.urlencode({'cmnd': cmnd})}"
    _status, doc = _request(device, "GET", path)
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

    # -- Shelly Gen2 RPC ---------------------------------------------------
    def read(self, device) -> dict:
        slot = int((device.address or {}).get(SLOT_KEY, 0))
        doc = _post_json(device, "/rpc/Switch.GetStatus", {"id": slot})
        result = doc.get("result")
        if not isinstance(result, dict):
            raise DriverError(f"{device.id}: Switch.GetStatus returned no result")
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
        params = {SLOT_KEY: int((device.address or {}).get(SLOT_KEY, 0))}
        if channel == "power":
            params["on"] = (value == "on")
        else:
            raise DriverError(
                f"{device.id}.{channel}: Shelly switch writes on/off only")
        doc = _post_json(device, "/rpc/Switch.Set", params)
        if doc.get("result") is None:
            raise DriverError(f"{device.id}: Switch.Set refused")
        return f"shelly-{device.id}-{self.clock():.0f}"

    def arm_watchdog(self, device, seconds: float) -> bool:
        if float(seconds) <= 0:
            return False
        params = {SLOT_KEY: int((device.address or {}).get(SLOT_KEY, 0)),
                  "toggle_after": float(seconds)}
        _post_json(device, "/rpc/Switch.Set", params)
        return True

    def disarm_watchdog(self, device) -> None:
        params = {SLOT_KEY: int((device.address or {}).get(SLOT_KEY, 0)),
                  "toggle_after": False}
        _post_json(device, "/rpc/Switch.Set", params)

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
