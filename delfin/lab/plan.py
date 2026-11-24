"""A workflow document for the simulator, judged before anything runs.

Input: a validated :class:`~delfin.lab.manifest.Manifest` and a workflow
document of the shape :mod:`delfin.lab.workflow` executes. Output: a
:class:`Proposal` -- every structural problem found, the kernel's own verdict
on every write the plan contains, and the steps that ask a person.

WHY THIS EXISTS. ``Kernel.plan_check`` is a dry run of *writes*: it reads
``step["device"]`` and ``step["channel"]`` and answers a ``set``-shaped
verdict for each entry. Give it a real workflow's own steps and every
operation that is not a ``set`` comes back wrong. Measured on the demo
manifest, on a seven-step plan of which four are writes:

    hood  set       ok=True
    stir  set       ok=True
    hp    set       ok=True
    heat  set       ok=True
    hot   wait      ok=False not_writable  hotplate_probe has no writable
                                            channel 'temperature'
    move  transfer  ok=False unknown_device no device None in the manifest
    uv    measure   ok=False not_writable  uvvis has no writable channel None

Three of three, all fabricated: the coordinator drives ``uvvis.run`` for a
measurement, reads a probe for a wait and moves a sample for a transfer, and
none of those is a write of the shape this check models. And an EMPTY step
list reports ``ok: true``. A planner told "the kernel refused this" for a
plan the kernel would run perfectly is being taught to distrust the one
component that must be trusted -- and told a plan with no steps is fine.

So the structural half is judged here, against the manifest the kernel pins,
and the write half is still judged by the kernel itself -- ``Proposal``
carries ``Kernel.plan_check``'s own answer for exactly the ``set`` steps and
copying that answer is never done. Two judges, each on what it can see.

Nothing in this module writes a device, takes a lock or touches the kernel's
state: :func:`validate` is pure and :func:`run` hands a document it has
already judged to :class:`~delfin.lab.workflow.Coordinator`, which is the one
runner. There is no second scheduler here.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Optional

OPS = ("set", "wait", "transfer", "measure")
CMP_OPS = (">=", "<=", ">", "<", "==", "!=")

#: Keys a step may carry, per operation. A closed set for the same reason the
#: manifest's keys are closed: a key nobody reads is a plan that says
#: something the run will not do, and silently dropping it is how a
#: ``value`` meant for one channel ends up applied to another.
STEP_KEYS = {
    "set": {"id", "op", "device", "channel", "value", "unit", "after",
            "module", "note"},
    "wait": {"id", "op", "device", "channel", "cmp", "value", "unit",
             "timeout_s", "poll_s", "after", "module", "note"},
    "transfer": {"id", "op", "sample", "from", "to", "position", "after",
                 "module", "note"},
    "measure": {"id", "op", "device", "sample", "reader", "after", "module",
                "note"},
}


@dataclass
class Proposal:
    """A judged workflow document. Nothing here has been executed."""

    ok: bool
    workflow: dict
    problems: list[str] = field(default_factory=list)
    #: Kernel.plan_check's own answer for the ``set`` steps, in plan order.
    write_check: dict = field(default_factory=dict)
    #: ``{(step id or index): reason}`` for steps the kernel says need a person.
    needs_approval: dict = field(default_factory=dict)
    #: Devices the plan names, and the ones it must never reach.
    devices: list[str] = field(default_factory=list)
    forbidden_touched: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {"ok": self.ok, "problems": list(self.problems),
                "needs_approval": dict(self.needs_approval),
                "devices": list(self.devices),
                "forbidden_touched": list(self.forbidden_touched),
                "write_check": self.write_check,
                "touched_hardware": False}


def _step_ids(steps: list) -> tuple[dict, list[str]]:
    """Each step's id (falling back to its index) and the problems found."""
    seen: dict = {}
    problems: list[str] = []
    for i, st in enumerate(steps):
        if not isinstance(st, dict):
            problems.append(f"step {i}: is {type(st).__name__}, not an object")
            continue
        sid = str(st.get("id") or "")
        if not sid:
            problems.append(f"step {i}: no id; a step is named so a resume "
                            "can say which one it stopped at")
            sid = f"step-{i}"
        if sid in seen:
            problems.append(f"step {i}: id {sid!r} is used twice; a resume "
                            "addresses a step by its id")
        seen[sid] = st
    return seen, problems


def _check_set(st, manifest, where: str, problems: list) -> Optional[dict]:
    """The ``set`` step as the kernel's dry run wants it, or None."""
    dev = manifest.devices.get(str(st.get("device") or ""))
    if dev is None:
        problems.append(f"{where}: device {st.get('device')!r} is not in the "
                        f"manifest (known: {', '.join(sorted(manifest.devices))})")
        return None
    ch = dev.writes.get(str(st.get("channel") or ""))
    if ch is None:
        writable = ", ".join(sorted(dev.writes)) or "nothing"
        problems.append(f"{where}: {dev.id} has no writable channel "
                        f"{st.get('channel')!r} (writable: {writable})")
        return None
    if "value" not in st:
        problems.append(f"{where}: no value for {dev.id}.{ch.name}")
        return None
    return {"device": dev.id, "channel": ch.name, "value": st.get("value"),
            "unit": st.get("unit")}


def _check_wait(st, manifest, where: str, problems: list) -> None:
    dev = manifest.devices.get(str(st.get("device") or ""))
    if dev is None:
        problems.append(f"{where}: device {st.get('device')!r} is not in the "
                        "manifest")
        return
    name = str(st.get("channel") or "")
    if dev.channel(name) is None:
        readable = ", ".join(sorted(set(dev.reads) | set(dev.writes))) or "nothing"
        problems.append(f"{where}: {dev.id} has no channel {name!r} "
                        f"(has: {readable})")
    if str(st.get("cmp") or "") not in CMP_OPS:
        problems.append(f"{where}: cmp {st.get('cmp')!r} is not one of "
                        f"{', '.join(CMP_OPS)}")
    if "value" not in st:
        problems.append(f"{where}: no value to compare against")
    try:
        timeout = float(st.get("timeout_s", 600))
    except (TypeError, ValueError):
        problems.append(f"{where}: timeout_s {st.get('timeout_s')!r} is not a "
                        "number of seconds")
    else:
        if timeout <= 0:
            problems.append(f"{where}: timeout_s must be positive")


def _check_transfer(st, manifest, workflow, where: str, problems: list) -> None:
    sid = str(st.get("sample") or "")
    samples = workflow.get("samples") or {}
    if sid not in samples:
        problems.append(f"{where}: sample {sid!r} is not declared in the "
                        f"document's samples ({', '.join(sorted(samples)) or 'none'})")
    frm, to = str(st.get("from") or ""), str(st.get("to") or "")
    for side, mod in (("from", frm), ("to", to)):
        if mod not in manifest.modules:
            problems.append(f"{where}: {side} {mod!r} is not a module of the "
                            f"manifest (modules: {', '.join(sorted(manifest.modules))})")
    if sid in samples:
        declared = str((samples[sid] or {}).get("module") or "")
        if frm and declared and declared != frm:
            problems.append(f"{where}: sample {sid} is declared in {declared!r}"
                            f", not {frm!r}")
    if to in manifest.modules:
        positions = manifest.modules[to].positions
        if str(st.get("position") or "") not in positions:
            problems.append(f"{where}: position {st.get('position')!r} is not "
                            f"declared in module {to} "
                            f"(positions: {', '.join(positions) or 'none'})")


def _check_measure(st, manifest, where: str, problems: list) -> None:
    sid = str(st.get("sample") or "")
    if not sid:
        problems.append(f"{where}: no sample named")
    dev = manifest.devices.get(str(st.get("device") or ""))
    if dev is None:
        problems.append(f"{where}: device {st.get('device')!r} is not in the "
                        "manifest")
        return
    if "run" not in dev.writes:
        problems.append(f"{where}: {dev.id} has no 'run' channel, so the "
                        "coordinator cannot trigger a measurement on it "
                        f"(writable: {', '.join(sorted(dev.writes)) or 'nothing'})")
    reader = str(st.get("reader") or st.get("device") or "")
    rdev = manifest.devices.get(reader)
    if not reader:
        problems.append(f"{where}: no reader: a measurement checks the "
                        "container's barcode at the module first")
    elif rdev is None:
        problems.append(f"{where}: reader {reader!r} is not in the manifest")
    elif rdev.module != dev.module:
        problems.append(f"{where}: reader {reader} is in {rdev.module}, the "
                        f"analyser {dev.id} in {dev.module}; a barcode is read "
                        "at the module that holds the sample")


def validate(manifest, workflow, kernel=None) -> Proposal:
    """Judge *workflow* against *manifest*. Pure: nothing is written.

    ``kernel`` is optional and only supplies the write half: when given, its
    ``plan_check`` decides every ``set`` step (ranges, steps, rate limits,
    interlocks, e-stop) and this function copies that answer verbatim. Without
    one the proposal says so in ``write_check`` rather than guessing.
    """
    workflow = workflow if isinstance(workflow, dict) else {}
    problems: list[str] = []
    steps = workflow.get("steps")
    if not isinstance(steps, list) or not steps:
        problems.append("the document has no steps; a run with nothing to do "
                        "is not a run")
        steps = []
    if not str(workflow.get("id") or "").strip():
        problems.append("the document has no id; the run record is filed "
                        "under it")
    samples = workflow.get("samples") or {}
    if not isinstance(samples, dict):
        problems.append("samples must be an object {sample id: {container, "
                        "module, position}}")
        samples = {}

    named, id_problems = _step_ids(steps)
    problems.extend(id_problems)
    seen_ids = {st.get("id") for st in steps if isinstance(st, dict) and st.get("id")}

    ordered: list[tuple[str, dict]] = []
    write_steps: list[tuple[str, dict]] = []
    for i, st in enumerate(steps):
        if not isinstance(st, dict):
            continue
        op = str(st.get("op") or "")
        sid = str(st.get("id") or f"step-{i}")
        where = f"step {sid!r}"
        if op not in OPS:
            problems.append(f"{where}: op {op!r} is not one of "
                            f"{', '.join(OPS)}")
            continue
        ordered.append((sid, st))
        unknown = sorted(set(st) - STEP_KEYS[op])
        if unknown:
            problems.append(f"{where}: {op} does not take "
                            f"{', '.join(unknown)}; a key nothing reads is a "
                            "plan that says what the run will not do")
        for prev in st.get("after") or []:
            if prev not in seen_ids:
                problems.append(f"{where}: after names {prev!r}, which is not "
                                "a step of this document")
        if op == "set":
            shaped = _check_set(st, manifest, where, problems)
            if shaped is not None:
                write_steps.append((sid, shaped))
        elif op == "wait":
            _check_wait(st, manifest, where, problems)
        elif op == "transfer":
            _check_transfer(st, manifest, workflow, where, problems)
        elif op == "measure":
            _check_measure(st, manifest, where, problems)

    # The write half: the kernel judges it, and its answer is copied, not
    # summarised. A plan whose writes the kernel would refuse is not a plan.
    write_check: dict = {}
    needs_approval: dict = {}
    if write_steps:
        if kernel is None or getattr(kernel, "manifest", None) is None:
            problems.append("the write steps were not judged: no pinned "
                            "manifest is loaded, and a plan is not checked "
                            "against a guess")
        else:
            write_check = kernel.plan_check(
                None, [shaped for _sid, shaped in write_steps])
            for (sid, _shaped), verdict in zip(write_steps,
                                              write_check.get("steps") or []):
                if not verdict.get("ok"):
                    problems.append(
                        f"step {sid!r}: the kernel refuses this write "
                        f"({verdict.get('code')}): {verdict.get('reason')}")
                elif verdict.get("needs_approval"):
                    needs_approval[sid] = (
                        f"{verdict.get('level')}: a person approves exactly "
                        "this action")

    devices: list[str] = []
    for _sid, st in ordered:
        for key in ("device", "reader"):
            d = str(st.get(key) or "")
            if d and d not in devices:
                devices.append(d)
        if str(st.get("op")) == "measure" and st.get("reader") is None:
            d = str(st.get("device") or "")
            if d and d not in devices:
                devices.append(d)
    forbidden = [d for d in devices
                 if any(c.level == "forbidden"
                        for c in (manifest.devices[d].writes.values()
                                  if d in manifest.devices else ()))]
    if forbidden:
        problems.append("the plan names a device with a forbidden channel: "
                        + ", ".join(sorted(forbidden)))

    return Proposal(ok=not problems, workflow=workflow, problems=problems,
                    write_check=write_check, needs_approval=needs_approval,
                    devices=devices, forbidden_touched=forbidden)


def load(path) -> dict:
    """Read a workflow document from *path*, or raise ``ValueError``."""
    from pathlib import Path
    p = Path(path)
    try:
        text = p.read_text(encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"cannot read {p.name}: {exc}") from None
    try:
        doc = json.loads(text)
    except ValueError as exc:
        raise ValueError(f"{p.name} is not JSON: {exc}") from None
    if not isinstance(doc, dict):
        raise ValueError(f"{p.name} must hold one JSON object, not "
                         f"{type(doc).__name__}")
    return doc


def clients_for(manifest, *, coordinator: str, clock=None,
                nodes: Optional[list] = None) -> dict:
    """One in-process :class:`EdgeClient` per node the manifest declares.

    The same in-process transport :mod:`delfin.lab.edge` documents for phase
    2 (one host, the simulator): the edge's own kernel still verifies the
    envelope, its manifest pin and its replay window, so this is a transport
    choice and not a way round a check.
    """
    import time as _time
    from . import edge as _edge
    clock = clock or _time.time
    from .safety import Kernel
    out: dict = {}
    for node in (nodes or sorted({m.node for m in manifest.modules.values()})):
        kernel = Kernel.load(node=node, clock=clock)
        if kernel.manifest is None:
            raise ValueError(f"node {node} has no pinned manifest: "
                             f"{kernel.load_error}")
        key = _edge.provision_key(node)
        edge_node = _edge.EdgeNode(node, kernel=kernel, key=key, clock=clock)
        out[node] = _edge.EdgeClient(
            node, key, transport=edge_node.handle,
            manifest_checksum=manifest.checksum, coordinator=coordinator,
            clock=clock)
    return out


def run(manifest, proposal: Proposal, run_id: str, *, person: str,
        coordinator: str = "lab-specialist", nodes=None, clock=None,
        sleep=None) -> dict:
    """Run a proposal that was judged. Refuses one that was not.

    The refusal is the point of the separate call: a runner that accepts an
    unjudged document makes :func:`validate` advisory, and an advisory check
    on the path to a hotplate is a comment.
    """
    if not isinstance(proposal, Proposal):
        raise ValueError("run() needs a Proposal; call validate() first")
    if not proposal.ok:
        raise ValueError("this plan was not judged valid, so it is not run: "
                         + "; ".join(proposal.problems))
    from .workflow import Coordinator
    clients = clients_for(manifest, coordinator=coordinator, clock=clock,
                          nodes=nodes)
    coord = Coordinator(manifest, clients, person=person, name=coordinator,
                        **({"clock": clock} if clock else {}),
                        **({"sleep": sleep} if sleep else {}))
    return coord.run(proposal.workflow, run_id)


__all__ = ["CMP_OPS", "OPS", "Proposal", "STEP_KEYS", "clients_for", "load",
           "run", "validate"]
