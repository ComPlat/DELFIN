"""A chemist's experiment becomes a workflow, and only the kernel runs it.

MEASURED, and the reason this file exists. Before this change the
self-driving-lab loop had no agent-facing entrance at all: `Coordinator` --
the only thing that runs a workflow across modules -- was constructed in
`tests/` and nowhere else, and the five lab tools stopped at a single write.
The model could turn one device to `on` and nothing more.

And the check that was supposed to cover a plan could not. `lab_plan_check`
is a dry run of WRITES: given a real workflow's own steps on the demo
manifest it answered

    hot   wait      ok=False  not_writable   hotplate_probe has no writable
                                             channel 'temperature'
    move  transfer  ok=False  unknown_device no device None in the manifest
    uv    measure   ok=False  not_writable   uvvis has no writable channel None

Three of three fabricated -- the coordinator drives `uvvis.run` for a
measurement, reads a probe for a wait, and moves a sample for a transfer,
none of which is a write of the shape that check models -- while an EMPTY
step list answered `ok: True`. A planner taught to trust that would either
distrust the kernel or ship a plan with no steps.

What is asserted here, and each is a case that can FAIL:

  * a workflow that would run is judged ok, and one that would not is
    refused with the specific defect named (a module, a position, a channel,
    a forbidden device, a range the kernel itself refuses);
  * running a plan that was NOT judged ok is refused -- validation is not
    advisory;
  * a hazardous step stops the run and says so, and the device it names is
    NOT in its target state afterwards (the honest-report case: "waiting for
    approval" must never read as "done");
  * the whole thing is reached through the model's own dispatch
    (``_DocToolExecutor.execute``) and through ``delfin-agent lab run``, not
    by calling the planner directly.

Simulator only. Nothing here needs a model, a network or a device.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from _lab_support import call, lab_fixture  # noqa: F401

REPO = Path(__file__).resolve().parents[1]


def _wf(**over) -> dict:
    """The two-module experiment: heat, hand over, measure."""
    doc = {
        "id": "wf-suzuki",
        "samples": {"S-001": {"container": "vial-7", "module": "synthesis",
                              "position": "S1"}},
        "steps": [
            {"id": "hood", "op": "set", "device": "fume_hood",
             "channel": "power", "value": "on"},
            {"id": "stir", "op": "set", "device": "stirrer",
             "channel": "power", "value": "on"},
            {"id": "hp_on", "op": "set", "device": "hotplate",
             "channel": "power", "value": "on"},
            {"id": "move", "op": "transfer", "sample": "S-001",
             "from": "synthesis", "to": "analysis", "position": "A1"},
            {"id": "uv", "op": "measure", "device": "uvvis", "sample": "S-001",
             "reader": "reader"},
        ],
    }
    doc.update(over)
    return doc


def _kernel():
    from delfin.lab.safety import Kernel
    return Kernel.load()


def _place_the_vial(ws):
    """The person's part: the container is in the analyser's position."""
    from delfin.agent import api_client as A
    from delfin.lab.drivers import simulator
    perms = A.KitToolPermissions(workspace=ws)
    perms.mode = "bypassPermissions"
    perms.task_session_id = "s-place"
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    for dev in ("reader", "uvvis"):
        ex.execute("lab_read", {"device": dev}, perms)
    simulator.place("reader", barcode="vial-7")
    simulator.place("uvvis", signal_next=0.42)


# ---------------------------------------------------------------------------
# The structural half: what the write-only check cannot see
# ---------------------------------------------------------------------------

def test_a_runnable_workflow_is_judged_runnable(lab):
    """The control case: the same document the run below completes.

    RED on the base for the reason above -- on the baseline `lab_run` is not
    a tool at all, so this file cannot even ask. Once the tool exists, this
    is the case that says the judgement is not "refuse everything".
    """
    from delfin.lab import plan
    p = plan.validate(_kernel().manifest, _wf(), _kernel())
    assert p.ok, p.problems
    assert p.needs_approval == {}          # nothing hazardous in this one
    assert set(p.devices) == {"fume_hood", "stirrer", "hotplate", "uvvis",
                              "reader"}


def test_what_the_write_only_check_gets_wrong_is_not_reported_as_a_refusal(lab):
    """`lab_plan_check` on this workflow's own steps vs. this judgement.

    The finding, asserted directly: three steps the kernel's write check
    calls refusals are not refusals, and a plan is not refused for them.
    """
    k = _kernel()
    steps = [{"device": s.get("device"), "channel": s.get("channel"),
              "value": s.get("value"), "unit": s.get("unit")}
             for s in _wf()["steps"]]
    write_ok = k.plan_check(None, steps)
    assert write_ok["ok"] is False
    refused = [v for v in write_ok["steps"] if not v.get("ok")]
    assert len(refused) == 2           # the transfer and the measurement
    from delfin.lab import plan
    assert plan.validate(k.manifest, _wf(), k).ok is True


def test_an_empty_plan_is_not_ok(lab):
    """`lab_plan_check` answers ok: True for no steps. A run of nothing is
    not a run, and a planner must not be told it is."""
    from delfin.lab import plan
    k = _kernel()
    assert k.plan_check(None, [])["ok"] is True
    p = plan.validate(k.manifest, {"id": "wf-empty", "steps": []}, k)
    assert p.ok is False
    assert any("no steps" in prob for prob in p.problems)


@pytest.mark.parametrize("mutation,expected", [
    # a module that does not exist
    (lambda d: d["steps"][3].update(to="purification"), "purification"),
    # a position the destination does not declare
    (lambda d: d["steps"][3].update(position="Z9"), "Z9"),
    # a reader in another module: a barcode must be read where the sample is
    (lambda d: d["steps"][4].update(reader="hotplate_probe"), "hotplate_probe"),
    # a comparison operator nobody implements
    (lambda d: d["steps"].insert(4, {"id": "w", "op": "wait",
                                     "device": "hotplate_probe",
                                     "channel": "temperature", "cmp": "=>",
                                     "value": 55}), "=>"),
    # an operation that does not exist
    (lambda d: d["steps"].insert(1, {"id": "z", "op": "centrifuge"}), "centrifuge"),
    # a key the operation does not read: the plan says what will not happen
    (lambda d: d["steps"][0].update(colour="blue"), "colour"),
    # a device with a forbidden channel
    (lambda d: d["steps"].insert(1, {"id": "by", "op": "set",
                                     "device": "mains_bypass",
                                     "channel": "power", "value": "on"}),
     "mains_bypass"),
    # two steps with one id: a resume cannot address them
    (lambda d: d["steps"][1].update(id="hood"), "hood"),
    # a step that depends on one that is not there
    (lambda d: d["steps"][1].update(after=["nowhere"]), "nowhere"),
])
def test_a_plan_with_that_defect_is_refused_and_names_it(lab, mutation, expected):
    from delfin.lab import plan
    doc = _wf()
    mutation(doc)
    p = plan.validate(_kernel().manifest, doc, _kernel())
    assert p.ok is False
    assert any(expected in prob for prob in p.problems), p.problems


def test_a_range_the_kernel_refuses_is_refused_with_the_kernels_own_words(lab):
    """The write half is the KERNEL's answer, not a second opinion here.

    A setpoint of 200 degC is outside the manifest's 20..150; the refusal
    must quote the kernel's code, so a planner reads the same sentence
    `lab_set` would have given it.
    """
    from delfin.lab import plan
    doc = _wf()
    doc["steps"].insert(4, {"id": "hot", "op": "set", "device": "hotplate",
                            "channel": "setpoint", "value": 200,
                            "unit": "degC"})
    p = plan.validate(_kernel().manifest, doc, _kernel())
    assert p.ok is False
    assert any("out_of_range" in prob for prob in p.problems), p.problems


def test_the_write_half_is_the_kernels_own_verdict(lab):
    """Copied, not re-derived: the proposal carries plan_check's answer."""
    from delfin.lab import plan
    p = plan.validate(_kernel().manifest, _wf(), _kernel())
    assert p.write_check["touched_hardware"] is False
    assert len(p.write_check["steps"]) == 3          # hood, stir, hp_on
    assert all(v["ok"] for v in p.write_check["steps"])


def test_an_unjudged_plan_is_not_run(lab):
    """Validation is not advisory: the runner refuses what was not passed."""
    from delfin.lab import plan
    k = _kernel()
    doc = _wf()
    doc["steps"].insert(1, {"id": "z", "op": "centrifuge"})
    p = plan.validate(k.manifest, doc, k)
    assert p.ok is False
    with pytest.raises(ValueError, match="not run"):
        plan.run(k.manifest, p, "run-unjudged", person="alice")
    # and nothing was written: no run directory was created
    from delfin.lab.workflow import run_dir
    assert not (run_dir("run-unjudged") / "checkpoint.json").exists()


# ---------------------------------------------------------------------------
# The model's own entrance
# ---------------------------------------------------------------------------

def test_the_tool_reaches_the_planner_and_the_runner(lab):
    """Through ``_DocToolExecutor.execute`` -- the path a model's call takes.

    A test that imported ``plan`` and called ``run`` would measure that the
    module exists. This one measures that the model can reach it, which is
    the half that was missing.
    """
    _place_the_vial(lab)
    dry = call(lab, "lab_run", {"workflow": _wf()})
    assert "error" not in dry, (
        "the tool is advertised but not executable -- the three api_client.py "
        "lists (_LAB_TOOL_NAMES, the dispatch tuple, _PLAN_READONLY_TOOLS) are "
        f"in scratch/lab-run-dispatch.patch and must be applied: {dry}")
    assert dry["ran"] is False
    assert dry["ok"] is True
    assert "nothing was written" in dry["note"]
    # the dry run really wrote nothing
    from delfin.lab.drivers import simulator
    assert simulator.raw("stirrer").get("values", {}).get("power") in (None, "off")

    out = call(lab, "lab_run", {"workflow": _wf(), "dry_run": False,
                                "run_id": "run-tool-1"})
    assert "error" not in out, out
    assert out["ran"] is True
    assert out["status"] == "completed", out
    assert out["measurements"][0]["value"] == 0.42
    assert out["measurements"][0]["container"] == "vial-7"
    assert out["record_written"] is True
    from delfin.lab.workflow import run_dir
    rec = json.loads((run_dir("run-tool-1") / "record.json").read_text())
    assert rec["status"] == "completed"
    assert rec["plan"]["id"] == "wf-suzuki"


def test_the_tool_refuses_a_plan_it_judged_bad_and_writes_nothing(lab):
    doc = _wf()
    doc["steps"].insert(1, {"id": "by", "op": "set", "device": "mains_bypass",
                            "channel": "power", "value": "on"})
    out = call(lab, "lab_run", {"workflow": doc, "dry_run": False,
                                "run_id": "run-refused"})
    assert out["ran"] is False
    assert out["ok"] is False
    assert any("mains_bypass" in p for p in out["problems"])
    from delfin.lab.workflow import run_dir
    assert not (run_dir("run-refused") / "checkpoint.json").exists()
    # the forbidden device was not touched
    from delfin.lab.drivers import simulator
    assert simulator.raw("mains_bypass").get("values", {}).get("power") != "on"


# ---------------------------------------------------------------------------
# The hazardous step: an honest "not yet"
# ---------------------------------------------------------------------------

def test_a_hazardous_step_stops_the_run_and_does_not_claim_it_ran(lab):
    """The case this file exists for.

    A person is not present, so the setpoint must NOT be reached. What the
    tool answers has to say that: status waiting_approval, the step that
    waits, the request id, and a note that no record exists because the run
    has not finished. And the device must still be cold.
    """
    _place_the_vial(lab)
    doc = _wf()
    doc["steps"].insert(4, {"id": "heat", "op": "set", "device": "hotplate",
                            "channel": "setpoint", "value": 60,
                            "unit": "degC"})
    dry = call(lab, "lab_run", {"workflow": doc})
    assert "error" not in dry, (
        "lab_run is advertised but not executable: apply "
        f"scratch/lab-run-dispatch.patch. {dry}")
    assert dry["ok"] is True
    assert dry["needs_approval"] == {
        "heat": "hazardous: a person approves exactly this action"}

    out = call(lab, "lab_run", {"workflow": doc, "dry_run": False,
                                "run_id": "run-needs-person"})
    assert out["ran"] is True
    assert out["status"] == "waiting_approval"
    assert out["waiting_step"] == "heat"
    assert out["approval_request"].startswith("lab-")
    # the steps before the hazardous one, in order: the plan puts `heat`
    # after the hand-over, so `move` has happened and `uv` has not
    assert out["done"] == ["hood", "stir", "hp_on", "move"]
    assert "uv" not in out["done"]
    assert out["record_written"] is False
    assert "NOT finished" in out["note"]

    # the device is not at the setpoint the plan asked for
    from delfin.lab.drivers import simulator
    values = simulator.raw("hotplate")["values"]
    assert values["setpoint"] != 60


def test_a_person_approving_lets_the_same_run_finish(lab):
    """The other half of the stop: after the grant, the run resumes.

    The grant is made the way a person makes it -- outside an agent session,
    through `approvals.grant` -- and the run picks up at the step that
    waited rather than repeating the completed ones.
    """
    _place_the_vial(lab)
    from delfin.lab import approvals, plan
    doc = _wf()
    doc["steps"].insert(4, {"id": "heat", "op": "set", "device": "hotplate",
                            "channel": "setpoint", "value": 60,
                            "unit": "degC"})
    k = _kernel()
    first = plan.run(k.manifest, plan.validate(k.manifest, doc, k),
                     "run-approved", person="alice")
    assert first["status"] == "waiting_approval"
    approvals.grant(first["pending"]["approval"], approver="alice")
    again = plan.run(k.manifest, plan.validate(k.manifest, doc, k),
                     "run-approved", person="alice")
    assert again["status"] == "completed", again.get("deviations")
    assert again["done"].count("hood") == 1        # not repeated
    assert "heat" in again["done"]
    assert again["measurements"][0]["value"] == 0.42


# ---------------------------------------------------------------------------
# The person's own entrance: delfin-agent lab run
# ---------------------------------------------------------------------------

def _json_head(text: str) -> dict:
    """The JSON object at the start of the CLI's stdout.

    `lab run` prints the answer as JSON and then one human sentence (the
    house style of `delfin-agent lab`), so the document ends where the
    sentence begins.
    """
    dec = json.JSONDecoder()
    obj, _end = dec.raw_decode(text.lstrip())
    return obj


def _lab_cli(args, ws, tmp_path):
    """`delfin-agent lab ...` in a child process, as a person runs it."""
    from conftest import child_env
    env = child_env(tmp_path)
    # The code under test, not whatever `delfin` the host's interpreter
    # finds first. Measured: without this the child imported a different
    # checkout and answered with a CLI that has no `lab run` at all, so the
    # test failed for a reason it does not name. A test binds to DELFIN,
    # never to which tree the host's PYTHONPATH happens to point at.
    env["PYTHONPATH"] = str(REPO)
    env.pop("DELFIN_AGENT_SESSION", None)          # this IS the person
    env["DELFIN_LAB_DIR"] = os.environ["DELFIN_LAB_DIR"]
    env["DELFIN_PERSON"] = "alice"
    return subprocess.run([sys.executable, "-m", "delfin.agent.cli", "lab",
                           *args], cwd=str(ws), env=env, text=True,
                          capture_output=True, timeout=180)


def test_the_cli_judges_a_plan_and_writes_nothing(lab, tmp_path):
    from delfin.lab.drivers import simulator
    doc = tmp_path / "plan.json"
    doc.write_text(json.dumps(_wf()), encoding="utf-8")
    proc = _lab_cli(["run", str(doc), "--dry-run"], lab, tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = _json_head(proc.stdout)
    assert out["ok"] is True and out["ran"] is False
    assert "nothing was written" in proc.stdout
    assert simulator.raw("stirrer").get("values", {}).get("power") != "on"


def test_the_cli_refuses_a_bad_plan_with_a_nonzero_exit(lab, tmp_path):
    doc = _wf()
    doc["steps"].insert(1, {"id": "z", "op": "centrifuge"})
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    proc = _lab_cli(["run", str(path)], lab, tmp_path)
    assert proc.returncode == 1
    assert "refused" in proc.stderr
    assert "centrifuge" in proc.stdout


def test_the_cli_stops_at_the_hazardous_step_and_says_it_has_not_run(lab,
                                                                   tmp_path):
    _place_the_vial(lab)
    doc = _wf()
    doc["steps"].insert(4, {"id": "heat", "op": "set", "device": "hotplate",
                            "channel": "setpoint", "value": 60,
                            "unit": "degC"})
    path = tmp_path / "heat.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    proc = _lab_cli(["run", str(path), "--run-id", "cli-run"], lab, tmp_path)
    assert proc.returncode == 0, proc.stderr
    out = _json_head(proc.stdout)
    assert out["status"] == "waiting_approval"
    assert "needs a person" in proc.stderr
    assert "It has NOT run." in proc.stderr
    # the id the person needs is printed, and the resume command with it
    assert out["pending"]["approval"] in proc.stderr
    assert "--resume" in proc.stderr


# ---------------------------------------------------------------------------
# The specialist definition
# ---------------------------------------------------------------------------

def test_the_specialist_is_a_type_and_the_mode_routes_to_it():
    """"One definition, two entrances", the pattern `workflow-operator` set.

    The route names the PRESET, so the prompt and the tool surface come from
    the one file; a mode with its own copy is the drift that pattern ended.
    """
    from delfin.agent import subagents as sa
    from delfin.agent.prompt_loader import PromptLoader
    sa.reload_subagent_presets()
    try:
        assert "lab-specialist" in sa.subagent_type_names()
        assert sa.PRESET_LOAD_PROBLEMS.get(
            str(Path(sa.__file__).parent / "pack" / "agents"
                / "lab-specialist_subagent.md")) in (None, ())
        loader = PromptLoader(repo_dir=REPO / "delfin" / "agent")
        assert loader.load_mode("lab")["route"] == ["lab-specialist"]
        body = loader.load_role_prompt("lab-specialist")
        assert "lab_run" in body
    finally:
        sa.reload_subagent_presets()


def test_the_specialist_names_the_tools_it_needs_and_nothing_that_widens():
    """A preset that cannot call `lab_run` is a file, not a specialist.

    `lab_run` and `lab_devices` are the job; `bash` and `write_file` are
    not, and their absence is what keeps this specialist from reaching the
    hardware by another route.
    """
    from delfin.agent import subagents as sa
    sa.reload_subagent_presets()
    try:
        preset = sa.SUBAGENT_PRESETS["lab-specialist"]
        assert "lab_run" in preset.tools
        assert "lab_devices" in preset.tools
        assert "lab_set" not in preset.tools      # the plan is the unit
        assert "bash" not in preset.tools
        assert "write_file" not in preset.tools
        assert "edit_file" not in preset.tools
        # and the executor agrees: a bash call is refused for this role
        from delfin.agent import api_client as A
        assert A._role_allowlist("lab-specialist") == frozenset(
            str(t).rsplit("__", 1)[-1] for t in preset.tools)
    finally:
        sa.reload_subagent_presets()


def test_the_specialist_is_told_the_kernel_decides_and_how_to_report_a_stop():
    """The rule, in the file the model reads.

    A specialist that may run a workflow and is not told that a waiting run
    has not happened will report a temperature it never reached. That is the
    failure this assertion is for; the prose is where it is prevented.
    """
    from delfin.agent import subagents as sa
    sa.reload_subagent_presets()
    try:
        body = sa.SUBAGENT_PRESETS["lab-specialist"].system_prompt
    finally:
        sa.reload_subagent_presets()
    assert "waiting_approval" in body
    assert "has not happened" in body or "has NOT happened" in body
    assert "lab approve" in body
    assert "dry_run" in body
    assert "simulator" in body

# ---------------------------------------------------------------------------
# Plan mode is a permission profile, and lab_run must not talk it out of one
# ---------------------------------------------------------------------------

def test_plan_mode_judges_a_plan_and_does_not_run_it(lab):
    """The read-only profile reaches the dry run and nothing further.

    ``_PLAN_READONLY_TOOLS`` is a set of NAMES, so a tool whose acting form
    differs from its dry form by an ARGUMENT cannot be admitted by it. When
    `lab_run` was put in that set (wave 44, a1's scratch patch), this call
    -- plan mode, `dry_run: false` -- was admitted by the name-only gate,
    the workflow ran, and the simulator's stirrer came back `power: 'on'`.
    Plan mode then means "unless the tool has a second form", which is not
    a permission profile.

    The input that makes this fail: add "lab_run" to
    ``_PLAN_READONLY_TOOLS`` and drop the ``dry_run`` branch from
    ``_is_honoured_dry_run``. The acting call then returns a run, and the
    simulator reading below catches it.
    """
    from delfin.lab.drivers import simulator

    def stirrer_power():
        return simulator.raw("stirrer").get("values", {}).get("power")

    before = stirrer_power()

    # Judging is the read-only thing plan mode is FOR, and it is admitted
    # with the tool's own default -- no argument needed.
    judged = call(lab, "lab_run", {"workflow": _wf()}, mode="plan")
    assert "error" not in judged, judged
    assert judged["ran"] is False
    assert judged["ok"] is True, judged["problems"]

    # Acting in plan mode is refused, and the refusal names the tool.
    acted = call(lab, "lab_run", {"workflow": _wf(), "dry_run": False},
                 mode="plan")
    assert "error" in acted, acted
    assert "plan" in acted["error"].lower(), acted
    assert "ran" not in acted or acted["ran"] is False, acted

    # The measurement that would have caught the original hole: the device
    # did not move. A refusal that still moved the stirrer is not a refusal.
    assert stirrer_power() == before, (
        "plan mode refused the call and the device moved anyway")

    # Control case: the same acting call is admitted in a profile that can
    # act, so the refusal above is the mode's, not a broken tool.
    from _lab_support import prepare_heating
    prepare_heating(lab)
    ok = call(lab, "lab_run", {"workflow": _wf(), "dry_run": False},
              mode="bypassPermissions")
    assert "error" not in ok, ok
    assert ok["ran"] is True, ok

# ---------------------------------------------------------------------------
# Advertised is not executable -- the defect this wave's arms both hit
# ---------------------------------------------------------------------------

def test_every_tool_the_specialist_is_offered_actually_dispatches(lab):
    """A name in the preset's list is not a capability.

    Measured in wave 44: `lab_run` was in the preset's `tools:`, in
    ``tool_catalog._DOC_TOOLS_OPENAI`` and in ``TOOL_NAMES``, and every call
    to it answered ``{"error": "Unknown tool: lab_run"}`` -- because the
    executor's dispatch tuple and ``_LAB_TOOL_NAMES`` are separate lists and
    neither named it. Two changes in this wave added tool names to an
    allow-list, passed their own tests, and reached nothing.

    So this calls each tool the preset offers, through the executor, with
    arguments that are benign, and fails on the one answer that means
    "unreachable". Everything else -- a refusal, a validation complaint, an
    empty result -- proves the call ARRIVED at the tool.

    The input that makes it fail: drop "lab_run" from the dispatch tuple in
    ``_DocToolExecutor.execute`` (or from ``_LAB_TOOL_NAMES``). The name
    stays in the preset; the call answers "Unknown tool".
    """
    from delfin.agent import api_client as A
    from delfin.agent import subagents as sa

    sa.reload_subagent_presets()
    try:
        preset = sa.SUBAGENT_PRESETS["lab-specialist"]
    finally:
        sa.reload_subagent_presets()

    benign = {
        "lab_devices": {},
        "lab_read": {"device": "stirrer"},
        "lab_plan_check": {"steps": []},
        "lab_run": {"workflow": _wf()},
        "lab_estop": {},
        "ask_user_question": {"question": "q",
                              "options": [{"label": "a"}, {"label": "b"}]},
        "read_file": {"path": "missing-file.txt"},
        "grep_file": {"pattern": "x"},
        "list_files": {},
        "report_verdict": {"status": "approve"},
        "subagent_result": {"sa_id": "none"},
        "subagent_message": {"sa_id": "none", "message": "m"},
        "task_create": {"subject": "s"},
        "task_update": {"task_id": 1},
        "task_list": {},
    }

    unreachable = []
    for name in preset.tools:
        out = call(lab, name, benign.get(name, {}), mode="bypassPermissions")
        err = str(out.get("error") or "")
        if "Unknown tool" in err:
            unreachable.append((name, err))
    assert not unreachable, (
        "advertised to the model but not executable: " + repr(unreachable))
