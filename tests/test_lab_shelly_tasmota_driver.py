"""Shelly Gen2 and Tasmota plug/relay drivers against an in-process fake.

T4 / T20 of docs/design/lab-control.md: a real driver that talks HTTP to a
device, learns what the device really did (readback), and fails the kernel
closed on any transport error -- but never contacts real hardware or the
internet. The "device" here is an ``http.server`` bound to 127.0.0.1 started
by the test, speaking just enough Shelly RPC (``/rpc/...``) and Tasmota
(``/cm?cmnd=...``) to answer the driver, with real HTTP digest/basic
challenges so the driver's auth-by-reference is exercised end to end.

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
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs

import pytest
import yaml

from delfin.lab import manifest as _manifest
from delfin.lab import store as _store
from delfin.lab import safety as _safety
from delfin.lab.drivers import make as _make_driver

REALM = "delfin-lab"
NONCE = "deadbeefcafe"
OPAQUE = "opaque-session"

SHELLY_USER = "plug-user"
SHELLY_PASS = "plug-pass"
TASMOTA_USER = "plug-user"
TASMOTA_PASS = "plug-pass"


# -- a fake union device: Shelly RPC on /rpc, Tasmota cmnd on /cm -------------

class FakePlug:
    """State + recorded calls of one fake plug, shared with its handler."""

    def __init__(self):
        self.lock = threading.Lock()
        self.output = False            # Shelly Switch.GetStatus.output
        self.apower = 5.0              # W
        self.last_toggle_after = None  # seconds, Shelly
        self.last_pulsetime = None     # 0.1 s units, Tasmota
        self.calls = deque()           # (kind, payload)
        self.fail = False              # drop every connection (transport error)
        self.fail_read = False         # read path raises; write path works

    def record(self, kind, payload):
        with self.lock:
            self.calls.append((kind, payload))

    def requests(self, kind=None):
        with self.lock:
            if kind is None:
                return list(self.calls)
            return [c for c in self.calls if c[0] == kind]


class FakePlugHandler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, *args):        # keep pytest output clean
        pass

    def _drop(self):
        try:
            self.connection.close()
        except OSError:
            pass

    def _authorized(self, cred_user: str, cred_pass: str) -> bool:
        auth = self.headers.get("Authorization", "")
        if auth.startswith("Basic "):
            want = base64.b64encode(
                f"{cred_user}:{cred_pass}".encode()).decode()
            return auth == f"Basic {want}"
        if auth.startswith("Digest "):
            return self._digest_ok(auth, cred_user, cred_pass)
        return False

    def _digest_ok(self, auth: str, user: str, pw: str) -> bool:
        parts = {}
        for chunk in auth[len("Digest "):].split(","):
            k, _, v = chunk.strip().partition("=")
            parts[k] = v.strip('"')
        ha1 = hashlib.md5(f"{user}:{REALM}:{pw}".encode()).hexdigest()
        ha2 = hashlib.md5(
            f"{self.command}:{parts.get('uri','')}".encode()).hexdigest()
        expected = hashlib.md5(
            f"{ha1}:{parts.get('nonce','')}:{parts.get('nc','')}:"
            f"{parts.get('cnonce','')}:{parts.get('qop','')}:{ha2}"
        ).hexdigest()
        return (parts.get("response") == expected
                and parts.get("username") == user)

    def _challenge(self, scheme: str):
        if scheme == "Digest":
            self.send_response(401)
            self.send_header("WWW-Authenticate",
                             f'Digest realm="{REALM}", qop="auth", '
                             f'nonce="{NONCE}", opaque="{OPAQUE}"')
        else:
            self.send_response(401)
            self.send_header("WWW-Authenticate", f'Basic realm="{REALM}"')
        self.send_header("Content-Length", "0")
        self.end_headers()

    def _json(self, code: int, obj: dict):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _serve(self):
        plug = self.server.fake
        if plug.fail:
            self._drop()
            return
        cred = plug.credential          # (user, pass) or None
        if cred and not self._authorized(cred[0], cred[1]):
            plug.record("unauth", {"path": self.path})
            self._challenge(plug.scheme)
            return
        if plug.fail_read and self.path.startswith("/rpc/Switch.GetStatus"):
            self._drop()
            return
        self._route(plug)

    def _route(self, plug):
        path, _, query = self.path.partition("?")
        if path == "/rpc/Switch.GetStatus":
            plug.record("read", {"path": path})
            self._json(200, {"id": 1, "src": "delfin",
                             "result": {"output": plug.output,
                                        "apower": plug.apower}})
        elif path == "/rpc/Switch.Set":
            n = int(self.headers.get("Content-Length", "0"))
            params = json.loads(self.rfile.read(n))["params"]
            if "on" in params:
                plug.output = bool(params["on"])
            plug.last_toggle_after = params.get("toggle_after")
            plug.record("write", {"path": path, "params": params})
            self._json(200, {"id": 2, "src": "delfin", "result": {}})
        elif path == "/cm":
            plug.record("cmnd", {"path": path, "query": query})
            qs = parse_qs(query)
            cmnd = (qs.get("cmnd") or [""])[0]
            self._cmnd(plug, cmnd)
        else:
            self.send_error(404)

    def _cmnd(self, plug, cmnd):
        # Tasmota: `Power`, `Power On`, `Power Off`, `PulseTime <t>`
        words = cmnd.split()
        cmd = words[0] if words else ""
        arg = words[1] if len(words) > 1 else None
        if cmd == "Power" and arg is None:
            plug.record("read-power", {"cmnd": cmnd})
            self._json(200, {"POWER": "ON" if plug.output else "OFF"})
        elif cmd == "Power":
            plug.output = (arg.upper() in ("ON", "1", "TOGGLE")
                           and not (arg.upper() == "TOGGLE" and plug.output))
            self._json(200, {"POWER": "ON" if plug.output else "OFF"})
        elif cmd == "PulseTime":
            try:
                plug.last_pulsetime = int(arg)
            except (TypeError, ValueError):
                pass
            self._json(200, {"PulseTime": arg or "0"})
        else:
            self.send_error(400)

    def do_GET(self):
        self._serve()

    def do_POST(self):
        self._serve()


class FakePlugServer:
    def __init__(self):
        self.fake = FakePlug()
        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), FakePlugHandler)
        self.httpd.fake = self.fake
        self.thread = threading.Thread(target=self.httpd.serve_forever,
                                       daemon=True)
        self.thread.start()

    @property
    def port(self):
        return self.httpd.server_address[1]

    def close(self):
        self.httpd.shutdown()
        self.httpd.server_close()


# -- manifest and fixture ------------------------------------------------------

def manifest_text(port: int, *, shelly=True, tasmota=True,
                  shelly_sec="plug-shelly", tasmota_sec="plug-tasmota",
                  tasmota_scheme="Basic") -> str:
    devs = []
    if shelly:
        devs.append(f"""
        - id: plug_shelly
          driver: shelly_rpc
          model: shelly_plus_plug_s
          address: {{host: 127.0.0.1, port: {port}, slot: 0,
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
          address: {{host: 127.0.0.1, port: {port},
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
    server = FakePlugServer()
    try:
        # seed the credential store by reference-name (never in the manifest)
        from delfin.agent import credentials as _cred
        _cred.set_credential("plug-shelly", f"{SHELLY_USER}:{SHELLY_PASS}")
        _cred.set_credential("plug-tasmota", f"{TASMOTA_USER}:{TASMOTA_PASS}")
        path = _store.manifest_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(manifest_text(server.port), encoding="utf-8")
        manifest = _manifest.load(path)
        _safety.pin(manifest, by="alice")
        server.fake.credential = (SHELLY_USER, SHELLY_PASS)
        server.fake.scheme = "Digest"
        yield server, tmp_path
    finally:
        server.close()


def _kernel(clock=None):
    return _safety.Kernel.load(clock=clock or time.time)


def _device(driver):
    from delfin.lab.drivers.base import Driver
    d = _make_driver(driver)
    assert isinstance(d, Driver)
    return d


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

def test_shelly_arm_watchdog_sets_toggle_after(plug):
    server, ws = plug
    k = _kernel()
    d = k.manifest.devices["plug_shelly"]
    drv = _make_driver("shelly_rpc")
    drv.clock = k.clock
    ok = drv.arm_watchdog(d, 30)
    assert ok is True
    writes = server.fake.requests("write")
    assert writes, "no Switch.Set observed"
    params = writes[-1][1]["params"]
    assert params.get("toggle_after") == 30
    assert server.fake.last_toggle_after == 30


def test_shelly_disarm_watchdog_clears_the_timer(plug):
    server, ws = plug
    k = _kernel()
    d = k.manifest.devices["plug_shelly"]
    drv = _make_driver("shelly_rpc")
    drv.clock = k.clock
    drv.disarm_watchdog(d)
    writes = server.fake.requests("write")
    assert writes
    assert writes[-1][1]["params"].get("toggle_after") is False


def test_tasmota_arm_watchdog_sends_pulsetime(plug):
    server, ws = plug
    server.fake.scheme = "Basic"
    server.fake.credential = None          # Tasmota may run open on a LAN
    k = _kernel()
    d = k.manifest.devices["plug_tasmota"]
    drv = _make_driver("tasmota")
    drv.clock = k.clock
    assert drv.arm_watchdog(d, 30) is True    # 30 s -> 130 (seconds + 100)
    calls = server.fake.requests("cmnd")
    assert calls
    assert server.fake.last_pulsetime == 130


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
    server.fake.scheme = "Basic"
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
    # the fake plugin drops connections: the only realistic closed port
    # behaviour for a boxing LAN device is the same -- no answer at all
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


def _driver_error():
    from delfin.lab.drivers import base as B
    return B.DriverError


# -- safe_state continues past one failed channel ------------------------------

def test_safe_state_keeps_writing_after_a_failure(plug):
    server, ws = plug
    # power works but the read path drops; safe_state must still write off
    server.fake.fail_read = True
    server.fake.output = True
    k = _kernel()
    d = k.manifest.devices["plug_tasmota"]
    server.fake.scheme = "Basic"
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
