"""The device manifest: what exists, what it may do, and what it needs.

Input: a YAML or JSON document (schema ``lab_manifest: 1``). Output: a
validated, immutable :class:`Manifest` with a checksum, or
``ManifestError`` naming every problem found.

Shape (all keys closed; an unknown key is an error, not ignored)::

    lab_manifest: 1
    installation: demo-lab          # one manifest per installation
    version: 3                      # raised by the person who edits it
    nodes:   [{id: edge_synth}]     # edge nodes that host drivers
    modules:
      - id: synthesis
        node: edge_synth
        min_level: benign           # floor for every write in the module
        positions: [A1, A2]         # sample positions (hand-over)
        devices:
          - id: hotplate
            driver: simulator
            model: hotplate          # simulator kind / vendor model
            address: {host: 192.0.2.10, secret_ref: "cred:shelly_hp"}
            stale_after_s: 5
            deadman_s: 30
            reads:  {temperature: {unit: degC}}
            writes:
              setpoint: {unit: degC, min: 20, max: 150, max_step: 30,
                         min_interval_s: 1, level: hazardous,
                         tolerance: 0.5}
              power:    {states: [on, off], level: benign}
            safe_state: {setpoint: 20, power: off}
            interlocks:
              - when: "hotplate.setpoint > 40"
                requires: "stirrer.power == on and fume_hood.power == on"

Rules checked here: identifiers are ``[a-z][a-z0-9_]{0,31}`` and unique per
installation; every unit is in :mod:`units`; ``min < max``; every writable
channel has a safe value inside its own range or states; levels are one of
benign / hazardous / forbidden (reads are free); interlock expressions
parse and name channels that exist; an address carries no secret inline
(a key that looks like one must be ``secret_ref: "cred:<name>"``, resolved
by the credential store at driver start, never by the kernel).

The checksum is SHA-256 over the canonical JSON of the parsed document
(sorted keys, no whitespace), so a YAML reformat does not change it and any
change of meaning does. The edge pins it (``pin``); a manifest whose
checksum differs from the pin is refused for every write.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from . import expr as _expr
from . import units as _units

LEVELS = ("benign", "hazardous", "forbidden")
_LEVEL_RANK = {"free": 0, "benign": 1, "hazardous": 2, "forbidden": 3}
_ID = re.compile(r"^[a-z][a-z0-9_]{0,31}$")
_SECRET_KEY = re.compile(r"(pass|secret|token|key|auth|cred)", re.I)
_SECRET_REF = re.compile(r"^cred:[A-Za-z0-9_.-]{1,64}$")
DRIVERS_IMPLEMENTED = ("simulator", "gpiod")
#: Named in the interface and the design doc; refused until written.
DRIVERS_PLANNED = ("mqtt", "home_assistant", "shelly_rpc", "tasmota",
                   "sila2", "opcua_lads", "serial", "modbus")


def state_str(v: Any) -> str:
    """A state value as text. YAML 1.1 reads bare ``on``/``off`` as
    booleans; they are the words a manifest means."""
    if v is True:
        return "on"
    if v is False:
        return "off"
    return str(v)


class ManifestError(ValueError):
    def __init__(self, problems: list[str]):
        self.problems = list(problems)
        super().__init__("; ".join(problems))


@dataclass(frozen=True)
class Channel:
    name: str
    unit: Optional[str] = None
    states: tuple = ()
    min: Optional[float] = None
    max: Optional[float] = None
    max_step: Optional[float] = None
    min_interval_s: float = 0.0
    level: str = "free"
    tolerance: float = 0.0
    readback: Optional[str] = None

    @property
    def is_state(self) -> bool:
        return bool(self.states)


@dataclass(frozen=True)
class Interlock:
    when: str
    requires: str
    device: str

    def text(self) -> str:
        return f"{self.when} requires {self.requires}"


@dataclass(frozen=True)
class Device:
    id: str
    module: str
    node: str
    driver: str
    model: str
    address: dict
    reads: dict
    writes: dict
    safe_state: dict
    interlocks: tuple
    stale_after_s: float
    deadman_s: float
    device_watchdog: bool

    def channel(self, name: str) -> Optional[Channel]:
        return self.writes.get(name) or self.reads.get(name)

    @property
    def hazardous(self) -> bool:
        return any(c.level == "hazardous" for c in self.writes.values())


@dataclass(frozen=True)
class Module:
    id: str
    node: str
    min_level: str
    positions: tuple
    devices: tuple


@dataclass(frozen=True)
class Manifest:
    installation: str
    version: int
    nodes: tuple
    modules: dict
    devices: dict
    checksum: str
    raw: dict = field(repr=False, default_factory=dict)

    def level(self, device: Device, channel: Channel) -> str:
        """The effective approval level: the channel's or the module floor."""
        floor = self.modules[device.module].min_level
        return max(channel.level, floor, key=lambda lv: _LEVEL_RANK[lv])

    def interlocks(self) -> list[Interlock]:
        out: list[Interlock] = []
        for d in self.devices.values():
            out.extend(d.interlocks)
        return out


def checksum_of(doc: dict) -> str:
    blob = json.dumps(doc, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def _closed(obj: Any, allowed: set, where: str, problems: list) -> dict:
    if not isinstance(obj, dict):
        problems.append(f"{where}: expected a mapping")
        return {}
    extra = sorted(set(obj) - allowed)
    if extra:
        problems.append(f"{where}: unknown key(s) {', '.join(extra)}")
    return obj


def _num(v, where, problems) -> Optional[float]:
    if v is None:
        return None
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        problems.append(f"{where}: expected a number, got {v!r}")
        return None
    return float(v)


def _channel(name, spec, where, writable, problems) -> Channel:
    allowed = {"unit", "states", "min", "max", "max_step", "min_interval_s",
               "level", "tolerance", "readback"} if writable else {
                   "unit", "states", "text"}
    spec = _closed(spec, allowed, where, problems)
    states = tuple(state_str(s) for s in (spec.get("states") or ()))
    unit = spec.get("unit")
    text = bool(spec.get("text")) and not writable
    if (states or text) and unit:
        problems.append(f"{where}: a state or text channel has no unit")
    if not states and not text:
        try:
            unit = _units.canonical(unit)
        except _units.UnitError as exc:
            problems.append(f"{where}: {exc}")
            unit = None
    level = str(spec.get("level", "benign" if writable else "free"))
    if writable and level not in LEVELS:
        problems.append(f"{where}: level must be one of {', '.join(LEVELS)}")
    lo = _num(spec.get("min"), where + ".min", problems)
    hi = _num(spec.get("max"), where + ".max", problems)
    if writable and not states and level != "forbidden" and (
            lo is None or hi is None):
        problems.append(f"{where}: a numeric setpoint needs min and max")
    if lo is not None and hi is not None and not lo < hi:
        problems.append(f"{where}: min must be below max")
    return Channel(
        name=name, unit=unit, states=states, min=lo, max=hi,
        max_step=_num(spec.get("max_step"), where + ".max_step", problems),
        min_interval_s=_num(spec.get("min_interval_s"), where, problems) or 0.0,
        level=level if writable else "free",
        tolerance=_num(spec.get("tolerance"), where, problems) or 0.0,
        readback=spec.get("readback"))


def _address(addr, where, problems) -> dict:
    if addr is None:
        return {}
    if not isinstance(addr, dict):
        problems.append(f"{where}: expected a mapping")
        return {}
    for k, v in addr.items():
        if isinstance(v, (dict, list)):
            problems.append(f"{where}.{k}: nested values are not allowed")
        if k == "secret_ref":
            if not (isinstance(v, str) and _SECRET_REF.match(v)):
                problems.append(f"{where}.secret_ref: must be 'cred:<name>' "
                                "(a credential-store reference)")
        elif _SECRET_KEY.search(str(k)):
            problems.append(f"{where}.{k}: looks like a secret; secrets are "
                            "never inline, use secret_ref: 'cred:<name>'")
    return dict(addr)


def parse_document(doc: Any) -> Manifest:
    problems: list[str] = []
    top = _closed(doc, {"lab_manifest", "installation", "version", "nodes",
                        "modules"}, "manifest", problems)
    if top.get("lab_manifest") != 1:
        problems.append("manifest: lab_manifest: 1 is required")
    inst = str(top.get("installation", "") or "")
    if not inst:
        problems.append("manifest: installation is required")
    version = top.get("version")
    if not isinstance(version, int) or isinstance(version, bool) or version < 1:
        problems.append("manifest: version must be a positive integer")
        version = 0
    nodes = []
    for i, n in enumerate(top.get("nodes") or []):
        n = _closed(n, {"id", "description"}, f"nodes[{i}]", problems)
        nid = str(n.get("id", ""))
        if not _ID.match(nid):
            problems.append(f"nodes[{i}]: bad id {nid!r}")
        nodes.append(nid)
    if not nodes:
        problems.append("manifest: at least one node is required")
    modules: dict[str, Module] = {}
    devices: dict[str, Device] = {}
    for mi, m in enumerate(top.get("modules") or []):
        m = _closed(m, {"id", "node", "min_level", "positions", "devices",
                        "description"}, f"modules[{mi}]", problems)
        mid = str(m.get("id", ""))
        if not _ID.match(mid) or mid in modules:
            problems.append(f"modules[{mi}]: bad or duplicate id {mid!r}")
        node = str(m.get("node", ""))
        if node not in nodes:
            problems.append(f"modules[{mi}]: node {node!r} is not declared")
        floor = str(m.get("min_level", "benign"))
        if floor not in LEVELS:
            problems.append(f"modules[{mi}]: min_level must be one of "
                            f"{', '.join(LEVELS)}")
            floor = "benign"
        dev_ids = []
        for di, d in enumerate(m.get("devices") or []):
            where = f"{mid}.devices[{di}]"
            d = _closed(d, {"id", "driver", "model", "address", "reads",
                            "writes", "safe_state", "interlocks",
                            "stale_after_s", "deadman_s", "device_watchdog",
                            "description"}, where, problems)
            did = str(d.get("id", ""))
            if not _ID.match(did) or did in devices:
                problems.append(f"{where}: bad or duplicate id {did!r}")
            where = f"{mid}.{did}"
            driver = str(d.get("driver", ""))
            if driver in DRIVERS_PLANNED:
                problems.append(f"{where}: driver {driver!r} is planned, not "
                                "implemented in this release")
            elif driver not in DRIVERS_IMPLEMENTED:
                problems.append(f"{where}: unknown driver {driver!r}")
            reads = {k: _channel(k, v, f"{where}.reads.{k}", False, problems)
                     for k, v in (d.get("reads") or {}).items()}
            writes = {k: _channel(k, v, f"{where}.writes.{k}", True, problems)
                      for k, v in (d.get("writes") or {}).items()}
            safe = d.get("safe_state") or {}
            if writes and not isinstance(safe, dict):
                problems.append(f"{where}: safe_state must be a mapping")
                safe = {}
            for ch_name, ch in writes.items():
                if ch.level == "forbidden":
                    continue
                if ch_name not in safe:
                    problems.append(f"{where}: safe_state has no value for "
                                    f"{ch_name}")
                    continue
                v = safe[ch_name]
                if ch.is_state and state_str(v) not in ch.states:
                    problems.append(f"{where}: safe value {v!r} for {ch_name} "
                                    "is not one of its states")
                if not ch.is_state and (
                        _num(v, f"{where}.safe_state.{ch_name}", problems)
                        is None or not ch.min <= float(v) <= ch.max):
                    problems.append(f"{where}: safe value {v!r} for {ch_name} "
                                    "is outside its range")
            extra_safe = set(safe) - set(writes)
            if extra_safe:
                problems.append(f"{where}: safe_state names unknown channel(s)"
                                f" {', '.join(sorted(extra_safe))}")
            locks = []
            for li, il in enumerate(d.get("interlocks") or []):
                il = _closed(il, {"when", "requires"},
                             f"{where}.interlocks[{li}]", problems)
                locks.append(Interlock(str(il.get("when", "") or ""),
                                       str(il.get("requires", "") or ""), did))
            devices[did] = Device(
                id=did, module=mid, node=node, driver=driver,
                model=str(d.get("model", "") or ""),
                address=_address(d.get("address"), f"{where}.address",
                                 problems),
                reads=reads, writes=writes,
                safe_state={k: (state_str(v) if k in writes and writes[k].is_state
                                else v) for k, v in safe.items()},
                interlocks=tuple(locks),
                stale_after_s=_num(d.get("stale_after_s", 5), where,
                                   problems) or 5.0,
                deadman_s=_num(d.get("deadman_s", 30), where, problems) or 30.0,
                device_watchdog=bool(d.get("device_watchdog", True)))
            dev_ids.append(did)
        modules[mid] = Module(mid, node, floor,
                              tuple(str(p) for p in (m.get("positions") or ())),
                              tuple(dev_ids))
    # Interlocks: parse and resolve every reference.
    for d in devices.values():
        for il in d.interlocks:
            for label, text in (("when", il.when), ("requires", il.requires)):
                try:
                    tree = _expr.parse(text)
                except (_expr.ExprError, _units.UnitError) as exc:
                    problems.append(f"{d.id}: interlock {label}: {exc}")
                    continue
                for ref in _expr.refs(tree):
                    other = devices.get(ref.device)
                    if other is None or other.channel(ref.channel) is None:
                        problems.append(f"{d.id}: interlock names unknown "
                                        f"channel {ref.name}")
    if problems:
        raise ManifestError(problems)
    return Manifest(inst, int(version), tuple(nodes), modules, devices,
                    checksum_of(top), raw=top)


def demo_text() -> str:
    """The simulated two-module demo lab (``examples/demo_lab.yaml``)."""
    from importlib import resources
    return (resources.files("delfin.lab") / "examples" / "demo_lab.yaml"
            ).read_text(encoding="utf-8")


def load(path: Path) -> Manifest:
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise ManifestError([f"cannot read {path.name}: {exc}"])
    if path.suffix == ".json":
        doc = json.loads(text)
    else:
        import yaml
        doc = yaml.safe_load(text)
    return parse_document(doc)
