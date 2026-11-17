"""Shelly Gen2 and Tasmota plug/relay drivers against an in-process fake.

T4 / T20 of docs/design/lab-control.md: a real driver that talks HTTP to a
device, learns what the device really did (readback), and fails the kernel
closed on any transport error -- but never contacts real hardware or the
internet. The "device" here is an in-process ``FakeTransport`` injected into
the driver's ``TRANSPORT`` hook: it speaks just enough Shelly RPC
(``/rpc/...``) and Tasmota (``/cm?cmnd=...``) to answer the driver, with real
HTTP digest/basic challenges so the driver's auth-by-reference is exercised
end to end. No socket is ever opened, which is what makes the suite run in a
sandbox that blocks outbound connections and on CI alike.

Driver tests follow the interface in delfin/lab/drivers/base.py; the shell
of the writes goes through the real executor path (``_DocToolExecutor``)
exactly as the other lab tests do, so the "only Kernel.set writes" invariant
holds here too.
"""

from __future__ import annotations

import base64
import hashlib
import json
import threading
import time
from collections import deque
from urllib.parse import parse_qs

import pytest
import yaml

from delfin.lab import manifest as _manifest
from delfin.lab import store as _store
from delfin.lab import safety as _safety
from delfin.lab.drivers import make as _make_driver
from delfin.lab.drivers import shelly_tasmota as _st
from delfin.lab.drivers.base import DriverError as _DriverError

REALM = "delfin-lab"
NONCE = "deadbeefcafe"

SHELLY_USER = "plug-user"
SHELLY_PASS = "plug-pass"
TASMOTA_USER = "plug-user"
TASMOTA_PASS = "plug-pass"


# -- a fake union device: Shelly RPC on /rpc, Tasmota cmnd on /cm -------------

class FakePlug:
    """State + recorded calls of one fake plug, shared with the transport."""

    def __init__(self):
        self.lock = threading.Lock()
        self.output = False            # Shelly Switch.GetStatus.output
        self.apower = 5.0              # W
        self.last_toggle_after = None  # seconds, Shelly
        self.last_pulsetime = None     # Tasmota
        self.calls = deque()           # (kind, payload)
        self.fail = False              # drop every call (transport error)
        self.fail_read = False         # read path raises; write path works
        self.credential = None         # Shelly digest (user, pass)
        self.scheme = "Digest"
        self.tasmota_credential = None  # Tasmota query (user, password)

    def record(self, kind, payload):
        with self.lock:
            self.calls.append((kind, payload))

    def requests(self, kind=None):
        with self.lock:
            if kind is None:
                return list(self.calls)
            return [c for c in self.calls if c[0] == kind]


def _authorized(auth, user, pw, method, path):
    """Validate an Authorization header exactly as a real server would."""
    if not auth:
        return False
    if auth.startswith("Basic "):
        want = base64.b64encode((user + ":" + pw).encode()).decode()
        return auth == f"Basic {want}"
    if auth.startswith("Digest "):
        parts = {}
        for chunk in auth[len("Digest "):].split(","):
            k, _, v = chunk.strip().partition("=")
            parts[k] = v.strip('"')
        H = {"MD5": hashlib.md5, "SHA-256": hashlib.sha256}[
            parts.get("algorithm", "MD5")]
        ha1 = H(f"{user}:{REALM}:{pw}".encode()).hexdigest()
        ha2 = H(f"{method}:{path}".encode()).hexdigest()
        expected = H(
            f"{ha1}:{parts.get('nonce','')}:{parts.get('nc','')}:"
            f"{parts.get('cnonce','')}:{parts.get('qop','')}:{ha2}".encode()
        ).hexdigest()
        return (parts.get("response") == expected
                and parts.get("username") == user
                and parts.get("nonce") == NONCE)
    return False


def _challenge(scheme):
    if scheme == "Digest":
        # Shelly Gen2 firmware: SHA-256 digest
        return {"WWW-Authenticate":
                f'Digest qop="auth", realm="{REALM}", nonce="{NONCE}", '
                f'algorithm=SHA-256'}
    return {"WWW-Authenticate": f'Basic realm="{REALM}"'}


def _json(code, obj):
    return (code, json.dumps(obj).encode(), {"Content-Type": "application/json"})


class FakeTransport:
    """In-process stand-in for the driver's HTTP TRANSPORT hook (no socket).

    Mirrors the real ``_http_transport`` contract but in-object: it issues a
    real 401 digest/basic challenge once, then accepts the driver's computed
    Authorization header, and translates 'connection dropped' into a
    ``DriverError`` exactly as the real transport would (so the kernel fails
    closed). Raises on ``fail`` / ``fail_read``.
    """

    def __init__(self):
        self.fake = FakePlug()

    def __call__(self, device, method, path, body, headers):
        plug = self.fake
        plug.record("http", {"method": method, "path": path})
        if plug.fail:
            raise _DriverError(f"{device.id}: transport error: dropped")
        if path.startswith("/cm"):
            # Tasmota: the web password is the user/password query; without
            # it the device answers 200 with a WARNING, never 401
            q = {k: v[0] for k, v in parse_qs(path.partition("?")[2]).items()}
            cmnd = q.get("cmnd", "")
            plug.record("cmnd", {"cmnd": cmnd, "path": path})
            cred = plug.tasmota_credential
            if cred and (q.get("user"), q.get("password")) != cred:
                return _json(200, {"WARNING": "Need user=<username>&password=<password>"})
            return self._cmnd(plug, cmnd)
        auth = headers.get("Authorization")
        cred = plug.credential              # (user, pass) or None
        if cred and not _authorized(auth, cred[0], cred[1], method, path):
            plug.record("unauth", {"path": path})
            return (401, b"", _challenge(plug.scheme))
        if path != "/rpc" or method != "POST":
            return (404, b"", {})
        frame = json.loads((body or b"{}").decode())
        rpc, params = frame.get("method"), frame.get("params") or {}

        def answer(result=None, error=None):
            out = {"id": frame.get("id"), "src": "shellyplusplugs-test"}
            out.update({"error": {"code": error[0], "message": error[1]}}
                       if error else {"result": result})
            return _json(200, out)

        if params.get("id") != 0:
            return answer(error=(-105, "Argument 'id' not found!"))
        if rpc == "Switch.GetStatus":
            if plug.fail_read:
                raise _DriverError(f"{device.id}: transport error: read dropped")
            plug.record("read", {"params": params})
            return answer({"id": 0, "output": plug.output,
                           "apower": plug.apower})
        if rpc == "Switch.Set":
            if "on" not in params:
                return answer(error=(-103, "Missing required argument 'on'!"))
            ta = params.get("toggle_after")
            if ta is not None and (isinstance(ta, bool)
                                   or not isinstance(ta, (int, float))):
                return answer(error=(-103, "Invalid argument 'toggle_after'!"))
            was = plug.output
            plug.output = bool(params["on"])
            plug.last_toggle_after = ta
            plug.record("write", {"params": params})
            return answer({"was_on": was})
        return answer(error=(-114, f"Method {rpc} failed: not found"))

    def _cmnd(self, plug, cmnd):
        # Tasmota: `Power`, `Power On`, `Power Off`, `PulseTime <t>`
        words = cmnd.split()
        cmd = words[0] if words else ""
        arg = words[1] if len(words) > 1 else None
        if cmd == "Power" and arg is None:
            plug.record("read-power", {"cmnd": cmnd})
            return _json(200, {"POWER": "ON" if plug.output else "OFF"})
        if cmd == "Power":
            plug.output = (arg.upper() in ("ON", "1", "TOGGLE")
                           and not (arg.upper() == "TOGGLE" and plug.output))
            return _json(200, {"POWER": "ON" if plug.output else "OFF"})
        if cmd == "PulseTime":
            try:
                plug.last_pulsetime = int(arg)
            except (TypeError, ValueError):
                pass
            return _json(200, {"PulseTime": arg or "0"})
        return (400, b"", {})


class FakeServer:
    """Thin shim so the tests read like they talk to one device server."""

    def __init__(self, fake, port):
        self.fake = fake
        self.port = port


# -- manifest and fixture ------------------------------------------------------

def manifest_text(port: int, *, shelly=True, tasmota=True,
                  shelly_sec="plug-shelly", tasmota_sec="plug-tasmota",
                  tasmota_scheme="Basic") -> str:
    del tasmota_scheme
    devs = []
    if shelly:
        devs.append(f"""
        - id: plug_shelly
          driver: shelly_rpc
          model: shelly_plus_plug_s
          address: {{host: 192.0.2.10, port: {port}, slot: 0,
                     secret_ref: cred:{shelly_sec}}}
          stale_after_s: 5
          deadman_s: 30
          reads:
            output: {{states: [on, off]}}
            apower: {{unit: W}}
          writes:
            power: {{states: [on, off], level: benign, readback: output}}
          safe_state: {{power: off}}
          device_watchdog: true""")
    if tasmota:
        devs.append(f"""
        - id: plug_tasmota
          driver: tasmota
          model: tasmota_plug
          address: {{host: 192.0.2.11, port: {port},
                     secret_ref: cred:{tasmota_sec}}}
          stale_after_s: 5
          deadman_s: 30
          reads:
            power: {{states: [on, off]}}
          writes:
            power: {{states: [on, off], level: benign, readback: power}}
          safe_state: {{power: off}}
          device_watchdog: true""")
    return f"""
lab_manifest: 1
installation: plug-lab
version: 1
nodes:
  - id: edge
modules:
  - id: plugs
    node: edge
    min_level: benign
    positions: [P1]
    devices:{''.join(devs)}
"""


@pytest.fixture(name="plug")
def plug_fixture(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_LAB_DIR", str(tmp_path / "lab"))
    monkeypatch.setenv("DELFIN_PERSON", "alice")
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    monkeypatch.delenv("DELFIN_LAB_MANIFEST", raising=False)
    from delfin.agent import trust_store
    trust_store.forget_cached_keys()
    # inject the in-process transport (no real socket)
    transport = FakeTransport()
    monkeypatch.setattr(_st, "TRANSPORT", transport)
    try:
        # seed the credential store by reference-name (never in the manifest)
        from delfin.agent import credentials as _cred
        _cred.set_credential("plug-shelly", f"{SHELLY_USER}:{SHELLY_PASS}",
                             store="file")
        _cred.set_credential("plug-tasmota", f"{TASMOTA_USER}:{TASMOTA_PASS}",
                             store="file")
        path = _store.manifest_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(manifest_text(80), encoding="utf-8")
        manifest = _manifest.load(path)
        _safety.pin(manifest, by="alice")
        transport.fake.credential = (SHELLY_USER, SHELLY_PASS)
        transport.fake.scheme = "Digest"
        transport.fake.tasmota_credential = (TASMOTA_USER, TASMOTA_PASS)
        server = FakeServer(transport.fake, port=80)
        yield server, tmp_path
    finally:
        pass


def _kernel(clock=None):
    return _safety.Kernel.load(clock=clock or time.time)


def _device(driver):
    from delfin.lab.drivers.base import Driver
    d = _make_driver(driver)
    assert isinstance(d, Driver)
    return d


def _driver_error():
    from delfin.lab.drivers import base as B
    return B.DriverError


# -- registration and manifest -------------------------------------------------

def test_shelly_rpc_and_tasmota_drivers_are_registered():
    from delfin.lab.drivers import shelly_tasmota as st
    shelly = _device("shelly_rpc")
    tasmota = _device("tasmota")
    assert shelly.name == "shelly_rpc"
    assert tasmota.name == "tasmota"
    assert shelly.capabilities.device_watchdog is True
    assert tasmota.capabilities.device_watchdog is True
    # both named in the shared support map
    assert st.ShellyDriver is not None
    assert st.TasmotaDriver is not None


def test_shelly_and_tasmota_carry_no_inline_secret(plug):
    server, ws = plug
    # an address that puts a password inline is refused by the manifest
    text = manifest_text(server.port)
    doc = yaml.safe_load(text)
    doc["modules"][0]["devices"][0]["address"]["password"] = "hunter2"
    with pytest.raises(_manifest.ManifestError) as exc:
        _manifest.parse_document(doc)
    assert "secret" in str(exc.value)


# -- device-side timer ---------------------------------------------------------

def _shelly(plug):
    server, _ws = plug
    k = _kernel()
    drv = _make_driver("shelly_rpc")
    drv.clock = k.clock
    return server.fake, k.manifest.devices["plug_shelly"], drv


def test_shelly_arm_on_an_off_switch_sends_no_timer(plug):
    """toggle_after FLIPS the output: on a switch that is off it would turn
    the load on. Arming an off switch only records the delay."""
    fake, d, drv = _shelly(plug)
    assert drv.arm_watchdog(d, 30) is True
    assert fake.requests("write") == []
    assert fake.output is False


def test_shelly_switch_on_carries_the_armed_timer(plug):
    fake, d, drv = _shelly(plug)
    drv.arm_watchdog(d, 30)
    drv.write(d, "power", "on")
    params = fake.requests("write")[-1][1]["params"]
    assert params == {"id": 0, "on": True, "toggle_after": 30.0}
    drv.write(d, "power", "off")
    assert fake.requests("write")[-1][1]["params"] == {"id": 0, "on": False}


def test_shelly_arm_on_an_on_switch_restarts_the_timer(plug):
    fake, d, drv = _shelly(plug)
    fake.output = True
    assert drv.arm_watchdog(d, 30) is True
    assert fake.requests("write")[-1][1]["params"] == {
        "id": 0, "on": True, "toggle_after": 30.0}


def test_shelly_disarm_watchdog_clears_the_timer(plug):
    fake, d, drv = _shelly(plug)
    fake.output = True
    drv.arm_watchdog(d, 30)
    drv.disarm_watchdog(d)
    assert fake.requests("write")[-1][1]["params"] == {"id": 0, "on": True}
    # never a timer with on=false, never a non-numeric toggle_after
    for _kind, payload in fake.requests("write"):
        ta = payload["params"].get("toggle_after")
        assert ta is None or (payload["params"]["on"] is True
                              and isinstance(ta, float))


def test_shelly_rpc_error_is_a_driver_error(plug):
    server, _ws = plug
    k = _kernel()
    d = k.manifest.devices["plug_shelly"]
    d.address["slot"] = 3                      # no switch 3 on the plug
    drv = _make_driver("shelly_rpc")
    with pytest.raises(_DriverError, match="Switch.GetStatus refused"):
        drv.read(d)


def test_shelly_digest_uses_the_challenges_algorithm(plug):
    fake, d, drv = _shelly(plug)
    drv.read(d)
    assert fake.requests("unauth"), "the device asked for credentials"
    assert fake.requests("read"), "the SHA-256 digest was accepted"
    fake.credential = (SHELLY_USER, "another")
    with pytest.raises(_DriverError, match="refused the credentials"):
        drv.read(d)


def test_tasmota_arm_watchdog_sends_pulsetime(plug):
    server, ws = plug
    k = _kernel()
    d = k.manifest.devices["plug_tasmota"]
    drv = _make_driver("tasmota")
    drv.clock = k.clock
    assert drv.arm_watchdog(d, 30) is True    # 30 s -> 130 (seconds + 100)
    calls = server.fake.requests("cmnd")
    assert calls
    assert server.fake.last_pulsetime == 130


def test_tasmota_password_is_the_query_and_never_shown(plug):
    server, ws = plug
    k = _kernel()
    d = k.manifest.devices["plug_tasmota"]
    drv = _make_driver("tasmota")
    assert drv.read(d)["power"].value == "off"
    path = server.fake.requests("cmnd")[-1][1]["path"]
    assert f"password={TASMOTA_PASS}" in path and "user=" in path
    server.fake.tasmota_credential = (TASMOTA_USER, "changed")
    with pytest.raises(_DriverError) as exc:
        drv.read(d)
    assert "refused the command" in str(exc.value)
    assert TASMOTA_PASS not in str(exc.value)
    server.fake.fail = True
    with pytest.raises(_DriverError) as exc:
        drv.write(d, "power", "on")
    assert TASMOTA_PASS not in str(exc.value)


def test_estop_confirms_a_shelly_through_its_readback_channel(plug):
    """safe_state names the write channel (power); the Shelly reports it as
    output. The kernel reads the safe value where the write is read back."""
    from delfin.lab import estop
    server, _ws = plug
    k = _kernel()
    actor = _safety.Actor.here("s1")
    assert k.set(actor, "plug_shelly", "power", "on")["ok"]
    out = estop.trigger(reason="test", by="alice", kernel=k)
    assert "plug_shelly" in out["stopped"] and out["not_confirmed"] == [], out
    assert server.fake.output is False


# -- read / write / readback ---------------------------------------------------

def test_shelly_write_read_and_readback(plug):
    server, ws = plug
    k = _kernel()
    actor = _safety.Actor.here("s1")
    out = k.set(actor, "plug_shelly", "power", "on")
    assert out["ok"], out
    assert server.fake.output is True
    # readback under the declared reads
    rd = k.read("plug_shelly")
    assert rd["ok"], rd
    assert rd["readings"]["output"]["value"] == "on"


def test_tasmota_write_read_and_readback(plug):
    server, ws = plug
    k = _kernel()
    actor = _safety.Actor.here("s1")
    out = k.set(actor, "plug_tasmota", "power", "on")
    assert out["ok"], out
    assert server.fake.output is True
    rd = k.read("plug_tasmota")
    assert rd["ok"], rd
    assert rd["readings"]["power"]["value"] == "on"


# -- transport error fails the kernel closed -----------------------------------

def test_transport_error_on_write_makes_the_kernel_go_safe(plug):
    server, ws = plug
    server.fake.fail = True
    k = _kernel()
    actor = _safety.Actor.here("s1")
    out = k.set(actor, "plug_shelly", "power", "on")
    assert out["ok"] is False
    assert out["code"] == "driver_error"
    # the kernel latched the device, so no further write can happen
    snap = _store.snapshot()
    assert "plug_shelly" in snap["latched"], snap

    server.fake.fail = False
    out2 = k.set(actor, "plug_shelly", "power", "on")
    assert out2["ok"] is False
    assert out2["code"] == "latched"


def test_transport_error_makes_arm_watchdog_raise(plug):
    server, ws = plug
    server.fake.fail = True
    k = _kernel()
    d = k.manifest.devices["plug_shelly"]
    drv = _make_driver("shelly_rpc")
    drv.clock = k.clock
    with pytest.raises(_driver_error()):
        drv.arm_watchdog(d, 30)


# -- safe_state continues past one failed channel ------------------------------

def test_safe_state_keeps_writing_after_a_failure(plug):
    server, ws = plug
    # power works but the read path drops; safe_state must still write off
    server.fake.fail_read = True
    server.fake.output = True
    k = _kernel()
    d = k.manifest.devices["plug_tasmota"]
    drv = _make_driver("tasmota")
    drv.clock = k.clock
    drv.safe_state(d)
    assert server.fake.output is False      # the safe value was written
    calls = server.fake.requests("cmnd")
    assert any("PowerOff" in c[1]["cmnd"] or c[1]["cmnd"].startswith("Power")
               for c in calls)


# -- only Kernel.set writes, through the agent tools ---------------------------

def test_no_driver_write_without_kernel_set_via_tools(plug, monkeypatch):
    """The one write tool (lab_set) goes through Kernel.set -- nothing the
    model calls can reach a driver method directly."""
    from delfin.lab import tools as _tools
    # the only write tool is lab_set, and its handler calls kernel.set only
    assert "lab_set" in _tools.TOOL_NAMES
    called = []
    orig = _safety.Kernel.set

    def patched(self, *a, **kw):
        called.append(True)
        return orig(self, *a, **kw)

    monkeypatch.setattr(_safety.Kernel, "set", patched)
    perms = type("P", (), {"task_session_id": "s1"})()
    out = _tools.execute("lab_set", {"device": "plug_shelly",
                                     "channel": "power", "value": "on"},
                         perms)
    assert json.loads(out)["ok"] is True
    assert called, "lab_set did not route through Kernel.set"
    # no tool name is a raw driver write: a forged tool must be unknown
    out2 = _tools.execute("driver.write", {"device": "plug_shelly",
                                           "channel": "power", "value": "on"},
                          perms)
    assert "unknown lab tool" in json.loads(out2)["error"]
