"""OPC UA LADS driver: behaviours required of a real-device driver.

Every behaviour is exercised against an in-process fake of an OPC UA LADS
server; no real device, no asyncua, no network. The driver's transport is
injected through :data:`delfin.lab.drivers.REGISTRY`, so these tests never
import asyncua and run on the CI interpreter without OPI.

The rules these protect: the driver never decides safety (the kernel does),
it reports what the device really did (readback with the device's own source
timestamp), it raises on any transport error so the kernel fails closed, and
the only way a value reaches a device is ``lab_set`` -> ``Kernel.set``.

Device model: one LADS FunctionalUnit ``pump``. ``flow`` is a measured
variable (read with a source timestamp) and also the setpoint written by a
LADS method; a server-side lock with a timeout stands in for the device
watchdog.
"""

from __future__ import annotations

import time
from pathlib import Path

import pytest

from tests._lab_support import approve_and_set, call

MANIFEST = """\
lab_manifest: 1
installation: opcua-test
version: 1
nodes:
  - id: edge_synth
    description: test node
modules:
  - id: synthesis
    node: edge_synth
    min_level: benign
    positions: [S1]
    devices:
      - id: pump
        driver: opcua_lads
        model: pump
        address:
          endpoint: opc.tcp://127.0.0.1:4840
          security: SignAndEncrypt
          timeout_s: 1.0
          secret_ref: "cred:opcua_test"
        deadman_s: 30
        reads:
          flow: {unit: mL/min}
        writes:
          flow: {unit: mL/min, min: 0, max: 10, max_step: 5,
                 level: hazardous, tolerance: 0.01}
        safe_state: {flow: 0}
"""


class FakeOpcUaServer:
    """An in-process model of a LADS device: one FunctionalUnit ``pump``.

    Variables: ``flow`` (the measured value, read with a source timestamp).
    Methods: ``set_flow(value)``. A server-side lock with a timeout stands
    in for the device watchdog. Adding a ``"<kind>:<channel>"`` key to
    :attr:`fail` makes the next such call raise like a transport error.
    """

    def __init__(self):
        self.flow = 0.0
        # The source timestamp must be near the kernel's wall clock, or the
        # kernel reads every value as stale and refuses a setpoint for that.
        self._ts = time.time()
        self.watchdog_armed = None   # seconds, or None when disarmed
        self.connect_count = 0
        self.fail = set()            # "<kind>:<channel>" that raises next call

    def _tick(self):
        self._ts += 1.0
        return self._ts

    def connect(self, endpoint, secret, security, timeout_s):
        self.connect_count += 1
        self._endpoint = endpoint
        self._security = security
        if "fail_on_connect" in self.fail:
            raise ConnectionError("server unreachable")

    def read(self, functional_unit, variable, timeout_s):
        self._maybe_fail("read", variable)
        if functional_unit != "pump":
            raise KeyError(f"no FunctionalUnit {functional_unit!r}")
        if variable == "flow":
            return self.flow, self._tick()
        raise KeyError(f"no variable {variable!r} on pump")

    def write(self, functional_unit, channel, value, timeout_s):
        self._maybe_fail("write", channel)
        if functional_unit != "pump":
            raise KeyError(f"no FunctionalUnit {functional_unit!r}")
        if channel == "flow":
            self.flow = float(value)
            return f"cmd-{int(self._tick())}"
        raise KeyError(f"no method {channel!r} on pump")

    def supports_watchdog(self, functional_unit):
        return functional_unit == "pump"

    def set_watchdog(self, functional_unit, seconds, timeout_s):
        self._maybe_fail("watchdog", "")
        self.watchdog_armed = seconds

    def clear_watchdog(self, functional_unit, timeout_s):
        self.watchdog_armed = None

    def close(self):
        pass

    def _maybe_fail(self, kind, channel):
        key = f"{kind}:{channel}" if channel else kind
        if key in self.fail:
            self.fail.remove(key)
            raise OSError(f"transport failure on {key}")


@pytest.fixture(name="lab_opcua")
def lab_opcua_fixture(tmp_path, monkeypatch):
    """A pinned custom manifest with an OPC UA pump + a fake server wired in."""
    from delfin.lab.drivers import REGISTRY
    from delfin.lab.drivers.opcua_lads import OpcUaLadsDriver

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_LAB_DIR", str(tmp_path / "lab"))
    monkeypatch.setenv("DELFIN_PERSON", "alice")
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    monkeypatch.delenv("DELFIN_LAB_MANIFEST", raising=False)
    from delfin.agent import trust_store
    trust_store.forget_cached_keys()

    server = FakeOpcUaServer()
    # The transport is injected so the driver never constructs a real client.
    monkeypatch.setitem(REGISTRY, "opcua_lads",
                        lambda: OpcUaLadsDriver(transport=server))

    from delfin.lab import manifest, store
    from delfin.lab.safety import pin
    path = store.manifest_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(MANIFEST, encoding="utf-8")
    pin(manifest.load(path), by="alice")
    ws = tmp_path / "ws"
    ws.mkdir()
    return server, ws


def test_read_returns_device_values_with_source_timestamp(lab_opcua):
    server, ws = lab_opcua
    server.flow = 3.5

    from delfin.lab.safety import Kernel
    kernel = Kernel.load()
    device = kernel.manifest.devices["pump"]
    drv = kernel.driver(device)

    got = drv.read(device)
    assert set(got) == {"flow"}
    r = got["flow"]
    assert r.value == 3.5
    assert r.unit == "mL/min"
    # The timestamp is when the DEVICE produced the value, not when asked.
    assert r.ts >= 1_700_000_000.0


def test_lab_read_tool_reports_device_value(lab_opcua):
    server, ws = lab_opcua
    server.flow = 2.0
    out = call(ws, "lab_read", {"device": "pump"})
    assert out["ok"] is True
    assert out["readings"]["flow"]["value"] == 2.0


def test_write_sets_device_and_readback_reports_what_it_did(lab_opcua):
    server, ws = lab_opcua
    out = approve_and_set(ws, {"device": "pump", "channel": "flow",
                               "value": 5, "unit": "mL/min"}, session="s1")
    assert out["ok"], out
    # The device really changed: the fake's own value reflects the write.
    assert server.flow == 5.0
    # The kernel read back the device's own value, not our command.
    assert out["readback"] == 5
    assert out["value"] == 5
    # The write went out as a LADS method call that returns a command id.
    assert str(out["command"]).startswith("cmd-")


def test_write_without_kernel_set_is_impossible_through_agent_tools(
        lab_opcua):
    server, ws = lab_opcua
    from delfin.lab import tools
    assert "lab_set" in tools.TOOL_NAMES
    # No tool reaches a device except lab_set; a direct write tool does not
    # exist, so the executor rejects it before any driver call.
    out = call(ws, "lab_nonexistent_write", {"device": "pump",
                                             "channel": "flow"})
    assert "Unknown tool" in out["error"]
    assert server.flow == 0.0
    # The only way flow changes is lab_set -> Kernel.set.
    approve_and_set(ws, {"device": "pump", "channel": "flow", "value": 4,
                         "unit": "mL/min"},
                    session="s1")
    assert server.flow == 4.0


def test_hazardous_write_arms_the_device_side_watchdog(lab_opcua):
    server, ws = lab_opcua
    out = approve_and_set(ws, {"device": "pump", "channel": "flow",
                               "value": 3, "unit": "mL/min"}, session="s1")
    assert out["ok"], out
    # The dead-man arms a server-side lock: the DEVICE returns to safe after
    # deadman_s seconds unless re-armed.
    assert server.watchdog_armed == 30


def test_transport_error_makes_the_kernel_go_safe(lab_opcua):
    server, ws = lab_opcua
    # A transport error during a write must race the kernel to safe (fail
    # closed), never silently continue.
    server.fail.add("write:flow")
    # value 5 is exactly max_step 5 from the safe value 0, so the refusal
    # comes from the transport, not the step gate.
    out = approve_and_set(ws, {"device": "pump", "channel": "flow",
                               "value": 5, "unit": "mL/min"}, session="s1")
    assert out["ok"] is False
    assert out["code"] == "driver_error"
    # The device was driven to its safe state (flow walked to 0).
    assert server.flow == 0.0


def test_health_reports_unreachable_device(lab_opcua):
    server, ws = lab_opcua
    from delfin.lab.safety import Kernel
    kernel = Kernel.load()
    device = kernel.manifest.devices["pump"]
    drv = kernel.driver(device)
    assert drv.health(device).ok is True
    # A transport error makes health report not-ok; it never raises.
    server.fail.add("read:flow")
    h = drv.health(device)
    assert h.ok is False
    assert h.detail


# -- the real transport: what it refuses before any network -------------------

def test_the_default_transport_is_the_asyncua_one():
    """It raised for every device: no real instrument could be reached."""
    from delfin.lab.drivers import opcua_lads
    assert isinstance(opcua_lads._default_transport(),
                      opcua_lads._AsyncUaTransport)


@pytest.mark.parametrize("security", ["None", "", "none"])
def test_an_unencrypted_endpoint_is_refused(security):
    from delfin.lab.drivers import opcua_lads
    from delfin.lab.drivers.base import DriverError
    with pytest.raises(DriverError, match="refused"):
        opcua_lads._AsyncUaTransport().connect(
            "opc.tcp://192.0.2.10:4840", '{"cert": "c", "key": "k", '
            '"server_cert": "s"}', security, 1.0)


@pytest.mark.parametrize("secret, missing", [
    ("not json", "JSON object"),
    ('{"cert": "c", "key": "k"}', "server_cert"),
    ('{"server_cert": "s"}', "cert, key"),
])
def test_the_credential_needs_both_certificates_and_the_key(secret, missing):
    """Without its own certificate the client cannot sign; without the
    device's it trusts whatever answers."""
    from delfin.lab.drivers import opcua_lads
    from delfin.lab.drivers.base import DriverError
    with pytest.raises(DriverError, match=missing):
        opcua_lads._AsyncUaTransport().connect(
            "opc.tcp://192.0.2.10:4840", secret, "SignAndEncrypt", 1.0)


def test_a_device_without_secret_ref_is_a_driver_error(lab_opcua):
    from delfin.lab.drivers.base import DriverError
    from delfin.lab.drivers.opcua_lads import OpcUaLadsDriver
    from delfin.lab.safety import Kernel
    device = Kernel.load().manifest.devices["pump"]
    device.address.pop("secret_ref")
    with pytest.raises(DriverError, match="secret_ref"):
        OpcUaLadsDriver()._connect_once(device)


def test_a_setpoint_reads_fresh_and_a_measurement_keeps_its_source_time(
        lab_opcua):
    """A setpoint the device just reported is current although its source
    timestamp is the time it last changed; a pure measurement keeps the
    source timestamp so a frozen sensor still reads stale."""
    server, _ws = lab_opcua
    from delfin.lab.drivers.opcua_lads import OpcUaLadsDriver
    from delfin.lab.safety import Kernel
    device = Kernel.load().manifest.devices["pump"]
    server._ts = 1_000.0                       # last changed long ago
    drv = OpcUaLadsDriver(transport=server, clock=lambda: 5_000.0)
    assert drv.read(device)["flow"].ts == 5_000.0      # flow is written too
    device.writes.pop("flow")                          # now a measurement
    assert drv.read(device)["flow"].ts == 1_002.0


def test_the_shipped_example_is_a_valid_manifest():
    """It used keys no manifest has (module, functions) and could not be
    pinned; power was written but not read, so its readback was unknown."""
    from delfin.lab import manifest
    path = (Path(__file__).resolve().parents[1] / "delfin" / "lab"
            / "examples" / "opcua_lads_example.yaml")
    m = manifest.load(path)
    dev = m.devices["opcua_pump"]
    assert dev.driver == "opcua_lads"
    assert set(dev.writes) <= set(dev.reads)
    assert dev.address["security"] == "SignAndEncrypt"
