"""Construction batch in the Submit Job tab: ``delfin cluster`` behind four buttons.

A list of ``ID;SMILES`` lines is built with MANTA, Architector, molSimplify or epic-MACE as a
sharded, resumable run (docs/CONSTRUCTION_BATCH.md).  The panel has no code path of its own:
every button turns the fields into the argument list of ``delfin cluster`` and runs it through
:func:`delfin.cluster_bench.cli.cbatch_dispatch`, so a run prepared here and one prepared on the
command line with the same settings are the same run directory.  The command line is shown under
the buttons, for the methods section and for reproducing a run without the dashboard.

On a Slurm backend Submit writes and submits the array script (partitions chosen the way every
DELFIN job gets them, ``DELFIN_MODULES`` loaded on the node); on the local backend the shards are
built one after the other in the background.
"""
from __future__ import annotations

import html
import os
import shlex
import subprocess
import sys
import threading
from pathlib import Path

import ipywidgets as widgets

from delfin.cluster_bench.provenance import MODES, TOOLS

#: Dropdown label -> ``delfin cluster`` tool name, in the order of the editor's buttons.
CB_TOOL_CHOICES = (("MANTA", "manta"), ("ARCHITECTOR", "architector"),
                   ("MOLSIMPLIFY", "molsimplify"), ("MACE", "mace"))
CB_RUNS_SUBDIR = "construction_batch"


def construction_batch_argv(action: str, f: dict) -> list:
    """The ``delfin cluster`` arguments of one button.  Options at their default are left out,
    as a user typing the command would leave them out."""
    run = str(f["run_dir"])
    if action == "prepare":
        argv = ["prepare", "--tool", f["tool"], "--input", str(f["input"]), "--run-dir", run]
        if f.get("mode") and f["mode"] != TOOLS[f["tool"]]["mode"]:
            argv += ["--mode", f["mode"]]
        for key, opt in (("select", "--select"), ("specs", "--specs"),
                         ("tool_python", "--tool-python")):
            if str(f.get(key) or "").strip():
                argv += [opt, str(f[key]).strip()]
        if int(f.get("shard_size") or 0) > 0:
            argv += ["--shard-size", str(int(f["shard_size"]))]
        if int(f.get("timeout") or 0) > 0 and int(f["timeout"]) != 21600:
            argv += ["--timeout", str(int(f["timeout"]))]
        if str(f.get("speed_factor") or "1.0").strip() not in ("1", "1.0"):
            argv += ["--speed-factor", str(f["speed_factor"]).strip()]
        if int(f.get("repeat") or 0) > 0:
            argv += ["--repeat", str(int(f["repeat"]))]
        return argv
    if action == "slurm":
        argv = ["slurm", run, "--submit"]
        if f.get("set", "main") != "main":
            argv += ["--set", f["set"]]
        if int(f.get("throttle", 40)) != 40:
            argv += ["--throttle", str(int(f["throttle"]))]
        if str(f.get("time_limit") or "72:00:00").strip() != "72:00:00":
            argv += ["--time", str(f["time_limit"]).strip()]
        for line in f.get("setup") or ():
            argv += ["--setup", line]
        return argv
    if action in ("status", "collect"):
        argv = [action, run]
        if f.get("set", "main") != "main":
            argv += ["--set", f["set"]]
        if action == "status":
            argv.append("-v")
        return argv
    if action == "repeat-stats":
        return ["repeat-stats", run]
    raise ValueError(f"unknown action {action!r}")


def construction_batch_command_line(argv) -> str:
    return "delfin cluster " + " ".join(shlex.quote(str(a)) for a in argv)


def construction_batch_call(argv) -> tuple:
    """Run ``delfin cluster ARGV`` in this process -> (exit code, its output as text)."""
    from delfin.cluster_bench.cli import cbatch_dispatch, cbatch_parser

    lines = []
    try:
        rc = cbatch_dispatch(cbatch_parser().parse_args([str(a) for a in argv]), emit=lines.append)
    except SystemExit as exc:          # the CLI's refusals: a message, nothing half-written
        rc = exc.code if isinstance(exc.code, int) else 1
        if not isinstance(exc.code, int) and exc.code is not None:
            lines.append(str(exc.code))
    except Exception as exc:  # noqa: BLE001 -- shown in the panel, the dashboard keeps running
        rc = 1
        lines.append(f"{type(exc).__name__}: {exc}")
    return rc, "\n".join(str(x) for x in lines)


def construction_batch_slurm_setup(backend) -> list:
    """Shell lines before the run on a node: the modules every DELFIN job of this site loads."""
    mods = os.environ.get("DELFIN_MODULES", "").strip()
    if not mods and backend is not None:
        site = getattr(backend, "_PROFILE_ENV", {}).get(getattr(backend, "slurm_profile", ""), {})
        mods = str(site.get("DELFIN_MODULES", "")).strip()
    return [f"module load {mods}"] if mods else []


def construction_batch_input_from_text(text: str) -> str:
    """The Batch SMILES field as an ``ID;SMILES`` list: ``name;SMILES;key=value`` keeps its
    first two fields; anything else is passed on unchanged, so the CLI's own check names it."""
    rows = []
    for ln in text.splitlines():
        if ln.strip():
            parts = ln.strip().split(";")
            rows.append(";".join(parts[:2]) if len(parts) >= 2 else ln.strip())
    return "".join(f"{ln}\n" for ln in rows)


def construction_batch_manifest(run_dir):
    from delfin.cluster_bench.prepare import cbatch_load_manifest

    try:
        return cbatch_load_manifest(run_dir)
    except SystemExit:
        return None


def construction_batch_run_locally(run_dir, sets, emit) -> threading.Thread:
    """Build every shard of ``sets`` one after the other, each in its own ``delfin cluster
    run-shard`` process (log in RUN/logs/local_<set>_<k>.log).  Returns the started thread."""
    man = construction_batch_manifest(run_dir)

    def _cb_local_work():
        for set_name in sets:
            for ent in man["sets"][set_name]["shards"]:
                k = ent["shard"]
                log = Path(run_dir) / "logs" / f"local_{set_name}_{k:04d}.log"
                with open(log, "a") as fh:
                    rc = subprocess.run([sys.executable, "-m", "delfin.cluster_bench", "run-shard",
                                         str(run_dir), "--shard", str(k), "--set", set_name],
                                        stdout=fh, stderr=subprocess.STDOUT).returncode
                emit(f"{set_name} shard {k}: exit {rc}")
        emit("local run finished -- press Collect")

    th = threading.Thread(target=_cb_local_work, daemon=True, name="construction-batch-local")
    th.start()
    return th


def create_construction_batch_panel(ctx, batch_text_widget) -> widgets.Accordion:
    """The panel; ``batch_text_widget`` is the tab's Batch SMILES field (used when no list
    file is given)."""
    style = {"description_width": "120px"}
    wide = widgets.Layout(width="100%")
    half = widgets.Layout(width="260px")

    tool = widgets.Dropdown(options=list(CB_TOOL_CHOICES), value="manta", description="Builder:",
                            style=style, layout=half)
    mode = widgets.Dropdown(options=list(MODES["manta"]), value=TOOLS["manta"]["mode"],
                            description="Mode:", style=style, layout=half)
    run_name = widgets.Text(value="", placeholder="name of the run", description="Run name:",
                            style=style, layout=half)
    list_path = widgets.Text(value="", placeholder="ID;SMILES file (empty: the Batch SMILES field)",
                             description="List file:", style=style, layout=wide)
    select_path = widgets.Text(value="", placeholder="optional: IDs to build, one per line",
                               description="Selection:", style=style, layout=wide)
    specs_path = widgets.Text(value="", placeholder="optional: specs JSONL of an earlier run",
                              description="Specs:", style=style, layout=wide)
    tool_python = widgets.Text(value="", description="Tool python:", style=style, layout=wide)
    shard_size = widgets.BoundedIntText(value=0, min=0, max=100000, description="Shard size:",
                                        style=style, layout=half)
    timeout = widgets.BoundedIntText(value=21600, min=1, max=10 ** 7, description="Limit/system s:",
                                     style=style, layout=half)
    speed = widgets.Text(value="1.0", description="Speed factor:", style=style, layout=half)
    throttle = widgets.BoundedIntText(value=40, min=0, max=100000, description="Throttle %N:",
                                      style=style, layout=half)
    repeat = widgets.BoundedIntText(value=0, min=0, max=100000, description="Repeat N:",
                                    style=style, layout=half)
    wall = widgets.Text(value="72:00:00", description="Job wall time:", style=style, layout=half)
    buttons = {name: widgets.Button(description=name, button_style=bs,
                                    layout=widgets.Layout(width="110px"))
               for name, bs in (("Prepare", "info"), ("Submit", "success"),
                                ("Status", ""), ("Collect", "warning"))}
    command = widgets.HTML("")
    out = widgets.Output(layout=widgets.Layout(max_height="320px", overflow_y="auto"))

    def _cb_tool_changed(*_):
        t = tool.value
        mode.options = list(MODES[t])
        mode.value = TOOLS[t]["mode"]
        if not TOOLS[t]["needs_specs"]:
            tool_python.placeholder = "MANTA runs in DELFIN's own interpreter"
            tool_python.disabled = True
            return
        tool_python.disabled = False
        try:
            from delfin.common import external_builders as eb

            tool_python.placeholder = "auto: " + eb.tool_python(t)
        except Exception:  # noqa: BLE001
            tool_python.placeholder = "auto"

    tool.observe(_cb_tool_changed, names="value")
    _cb_tool_changed()

    def _cb_fields() -> dict:
        name = run_name.value.strip()
        if not name or "/" in name or name.startswith("."):
            raise ValueError("Run name: a plain directory name is needed.")
        return {"tool": tool.value, "mode": mode.value,
                "run_dir": Path(ctx.calc_dir) / CB_RUNS_SUBDIR / name,
                "select": select_path.value, "specs": specs_path.value,
                "tool_python": tool_python.value, "shard_size": shard_size.value,
                "timeout": timeout.value, "speed_factor": speed.value, "repeat": repeat.value,
                "throttle": throttle.value, "time_limit": wall.value}

    def _cb_show(argv, rc, text):
        line = construction_batch_command_line(argv)
        command.value = f"<code style='font-size:12px'>{html.escape(line)}</code>"
        with out:
            print(f"$ {line}")
            if text:
                print(text)
            if rc:
                print(f"(exit {rc})")

    def _cb_prepared(f):
        man = construction_batch_manifest(f["run_dir"])
        if man is None:
            raise ValueError("Prepare first.")
        return man

    def _cb_prepare():
        f = _cb_fields()
        run = f["run_dir"]
        src = list_path.value.strip()
        if not src:
            text = construction_batch_input_from_text(batch_text_widget.value)
            if not text:
                raise ValueError("No list file and the Batch SMILES field is empty.")
            src = run.parent / f"{run.name}.input.txt"
            if src.exists() and src.read_text() != text:
                raise ValueError(f"{src} exists with another list; choose a new run name.")
            src.parent.mkdir(parents=True, exist_ok=True)
            src.write_text(text)
        f["input"] = src
        argv = construction_batch_argv("prepare", f)
        rc, text = construction_batch_call(argv)
        _cb_show(argv, rc, text)
        man = construction_batch_manifest(run) if rc == 0 else None
        if man:
            s = man["settings"]
            biggest = max(e["n"] for e in man["sets"]["main"]["shards"])
            bound_h = -(-biggest // s["workers"]) * s["timeout_s"] / 3600
            with out:
                print(f"{man['n_systems']} systems; "
                      + ", ".join(f"{k}: {v['n_shards']} shard(s)" for k, v in man["sets"].items())
                      + f"; limit {s['timeout_s']} s per system, {s['workers']} builds per node: "
                      f"a shard of {biggest} takes at most {bound_h:.1f} h even if every system "
                      "ran to the limit.")

    def _cb_submit():
        f = _cb_fields()
        sets = list(_cb_prepared(f)["sets"])
        if getattr(ctx.backend, "backend_name", "") == "SlurmJobBackend":
            setup = construction_batch_slurm_setup(ctx.backend)
            for set_name in sets:
                argv = construction_batch_argv("slurm", dict(f, set=set_name, setup=setup))
                _cb_show(argv, *construction_batch_call(argv))
            return
        with out:
            print("No Slurm here: the shards are built one after the other in the background "
                  f"({', '.join(sets)}); Status shows the progress.")
        construction_batch_run_locally(f["run_dir"], sets, lambda m: out.append_stdout(m + "\n"))

    def _cb_status():
        f = _cb_fields()
        for set_name in _cb_prepared(f)["sets"]:
            argv = construction_batch_argv("status", dict(f, set=set_name))
            _cb_show(argv, *construction_batch_call(argv))

    def _cb_collect():
        f = _cb_fields()
        man = _cb_prepared(f)
        for set_name in man["sets"]:
            argv = construction_batch_argv("collect", dict(f, set=set_name))
            _cb_show(argv, *construction_batch_call(argv))
        if "repeat" in man["sets"]:
            argv = construction_batch_argv("repeat-stats", f)
            _cb_show(argv, *construction_batch_call(argv))

    def _cb_on_click(fn):
        def _cb_click(_button):
            out.clear_output()
            try:
                fn()
            except Exception as exc:  # noqa: BLE001
                with out:
                    print(f"Error: {exc}")
        return _cb_click

    for name, fn in (("Prepare", _cb_prepare), ("Submit", _cb_submit), ("Status", _cb_status),
                     ("Collect", _cb_collect)):
        buttons[name].on_click(_cb_on_click(fn))

    row = dict(layout=widgets.Layout(gap="8px", flex_wrap="wrap"))
    body = widgets.VBox([
        widgets.HTML("Many SMILES with one builder, sharded and resumable "
                     "(<code>delfin cluster</code>, see docs/CONSTRUCTION_BATCH.md). "
                     f"Runs go to <code>{CB_RUNS_SUBDIR}/&lt;run name&gt;</code> in the "
                     "calculation folder."),
        widgets.HBox([tool, mode, run_name], **row),
        list_path, select_path, specs_path, tool_python,
        widgets.HBox([shard_size, timeout, speed], **row),
        widgets.HBox([throttle, repeat, wall], **row),
        widgets.HBox(list(buttons.values()), **row),
        command, out,
    ])
    acc = widgets.Accordion(children=[body])
    acc.set_title(0, "Construction batch (MANTA / ARCHITECTOR / MOLSIMPLIFY / MACE)")
    acc.selected_index = None
    acc.construction_batch_widgets = {
        "tool": tool, "mode": mode, "run_name": run_name, "list_path": list_path,
        "select_path": select_path, "specs_path": specs_path, "tool_python": tool_python,
        "shard_size": shard_size, "timeout": timeout, "speed": speed, "throttle": throttle,
        "repeat": repeat, "wall": wall, "buttons": buttons, "output": out, "command": command}
    return acc
