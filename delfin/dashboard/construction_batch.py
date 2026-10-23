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

Only Builder, Mode, Run name and List file are shown by default; the other options sit in a
collapsed Advanced section with their real defaults (shard size of the builder, the time limit in
hours, the computed wall time).  Defaults are left out of the command, as a user would leave them
out, so the panel never changes what the command line does.
"""
from __future__ import annotations

import datetime
import html
import json
import os
import shlex
import subprocess
import sys
import threading
from pathlib import Path

import ipywidgets as widgets

from delfin.cluster_bench.provenance import MODES, TOOLS, cbatch_effective_timeout

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
        if int(f.get("shard_size") or 0) not in (0, TOOLS[f["tool"]]["shard_size"]):
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


#: What each mode builds, in one line, shown under the Mode dropdown (docs/CONSTRUCTION_BATCH.md).
CB_MODE_HELP = {
    ("manta", "champion"): "the shipped MANTA construction: every isomer and conformer frame, "
                           "deterministic (recommended)",
    ("manta", "builder"): "MANTA's force-field-free builder with its base switches only, "
                          "without the shipped construction switches",
    ("architector", "full"): "10 symmetries x 10 conformers: every distinct isomer Architector finds",
    ("architector", "default"): "Architector's own defaults: one isomer per core geometry",
    ("molsimplify", "full"): "one structure per geometry of the coordination number",
    ("mace", "paper"): "octahedral and square-planar geometries, every stereomer, ten conformers each",
    ("mace", "extended"): "as paper, plus SPY, TBP, TET and SAN geometries",
}
#: ``delfin cluster slurm`` defaults the panel shows instead of leaving them implicit.
CB_DEFAULT_WALL = "72:00:00"
CB_DEFAULT_WALL_H = 72
CB_CPUS_PER_JOB = 48
_CB_HELP = "<span style='color:#757575;font-size:11.5px'>{}</span>"
_CB_WARN = "<span style='color:#e65100;font-size:11.5px'>{}</span>"


def construction_batch_field_rows(text: str) -> tuple:
    """The Batch SMILES field -> ([(id, smiles)], [reasons with the FIELD's line numbers]),
    by the CLI's own line rule (``name;SMILES;key=value`` keeps its first two fields)."""
    from delfin.cluster_bench.prepare import cbatch_check_line

    rows, errors, seen = [], [], set()
    for i, ln in enumerate(text.splitlines(), 1):
        if not ln.strip():
            continue
        parts = ln.strip().split(";")
        row = cbatch_check_line(";".join(parts[:2]) if len(parts) >= 2 else ln.strip())
        if isinstance(row, str):
            errors.append(f"line {i}: {row}")
        elif row[0] in seen:
            errors.append(f"line {i}: ID {row[0]} twice")
        else:
            seen.add(row[0])
            rows.append(row)
    return rows, errors


def construction_batch_default_run_name(tool: str, n, root, today=None) -> str:
    """``<builder>_<YYYYMMDD>_<n>mol``, with ``_2``, ``_3`` ... when that folder is taken."""
    day = (today or datetime.date.today()).strftime("%Y%m%d")
    base = f"{tool}_{day}_{n}mol" if n else f"{tool}_{day}"
    name, k = base, 2
    while (Path(root) / name).exists() or (Path(root) / f"{name}.input.txt").exists():
        name, k = f"{base}_{k}", k + 1
    return name


def construction_batch_plan(shard_ns, timeout_s, workers, throttle=40,
                            cpus=CB_CPUS_PER_JOB) -> dict:
    """Upper bounds of a run from its shard sizes: every molecule running to its limit."""
    shard_ns = [int(n) for n in shard_ns if int(n) > 0] or [0]
    bound_s = -(-max(shard_ns) // max(1, int(workers))) * int(timeout_s)
    n_jobs = len(shard_ns)
    parallel = min(n_jobs, int(throttle)) if int(throttle or 0) > 0 else n_jobs
    return {"n": sum(shard_ns), "n_shards": n_jobs, "biggest": max(shard_ns),
            "bound_h": bound_s / 3600, "parallel": parallel,
            "wall_h": -(-n_jobs // max(1, parallel)) * bound_s / 3600,
            "core_h": n_jobs * bound_s * cpus / 3600}


def construction_batch_summary(man: dict, run_dir, throttle=40) -> str:
    """The first line Prepare shows: what was prepared and what it can cost at most."""
    s = man["settings"]
    ns = [e["n"] for v in man["sets"].values() for e in v["shards"]]
    p = construction_batch_plan(ns, s["timeout_s"], s["workers"], throttle)
    rep = man["sets"].get("repeat")
    shards = (f"{man['sets']['main']['n_shards']} shard(s)"
              + (f" + {rep['n_shards']} repeat shard(s)" if rep else ""))
    return (f"{man['n_systems']} molecules, {shards}; at most {p['core_h']:,.0f} core-hours and "
            f"~{p['wall_h']:.1f} h wall with {p['parallel']} parallel job(s) (worst case: every "
            f"molecule runs to its limit of {s['timeout_s'] / 3600:.1f} h); output folder {run_dir}")


def construction_batch_status_line(st: dict) -> str:
    """One line per shard set: done/total, ok/timeout/failed, shards in progress."""
    by = st["by_class"]
    other = ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in sorted(by.items())
                      if k not in ("ok", "timeout", "fail") and v)
    running = sum(1 for r in st["shards"] if r["state"] in ("started", "running/partial"))
    return (f"{st['set']}: {st['n_done']}/{st['n_systems']} done -- ok {by.get('ok', 0)}, "
            f"timeout {by.get('timeout', 0)}, failed {by.get('fail', 0)}"
            + (f", {other}" if other else "")
            + f" -- shards {st['shards_done']}/{st['n_shards']} complete, {running} in progress")


def _cb_count_lines(path):
    try:
        return sum(1 for ln in Path(path).read_text(encoding="utf-8").splitlines() if ln.strip())
    except (OSError, UnicodeDecodeError):
        return None


def create_construction_batch_panel(ctx, batch_text_widget) -> widgets.Accordion:
    """The panel; ``batch_text_widget`` is the tab's Batch SMILES field (used when no list
    file is given).  Builder, Mode, Run name, List file and the buttons are always shown; the
    rest sits in a collapsed Advanced section, every field with its real default shown."""
    style = {"description_width": "90px"}
    adv_style = {"description_width": "250px"}
    wide = widgets.Layout(width="100%")
    half = widgets.Layout(width="300px")
    adv = widgets.Layout(width="400px")
    runs_root = Path(ctx.calc_dir) / CB_RUNS_SUBDIR

    def _help(text=""):
        return widgets.HTML(_CB_HELP.format(text))

    def _row(w, h):
        return widgets.HBox([w, h], layout=widgets.Layout(gap="10px", align_items="center",
                                                          flex_wrap="wrap"))

    tool = widgets.Dropdown(options=list(CB_TOOL_CHOICES), value="manta", description="Builder:",
                            style=style, layout=half,
                            tooltip="Which program builds the structures of every molecule.")
    mode = widgets.Dropdown(options=[], description="Mode:", style=style, layout=half,
                            tooltip="What the builder builds per molecule.")
    run_name = widgets.Text(value="", description="Run name:", style=style, layout=half,
                            tooltip="Folder of this run; empty = an automatic name.")
    list_path = widgets.Text(value="", placeholder="empty = the Batch SMILES field above",
                             description="List file:", style=style, layout=wide,
                             tooltip="A file with one ID;SMILES per line; empty = the Batch "
                                     "SMILES field above.")
    mode_help, run_help, list_help, plan = _help(), _help(), _help(), widgets.HTML("")

    select_path = widgets.Text(value="", placeholder="empty = all molecules",
                               description="Only these IDs (file):", style=adv_style, layout=adv,
                               tooltip="A file with the IDs to build, one per line.")
    specs_path = widgets.Text(value="", placeholder="empty = made from the SMILES",
                              description="Reuse specs (JSONL file):", style=adv_style, layout=adv,
                              tooltip="Ligand/geometry specs of an earlier run of an external "
                                      "builder, so both runs get the same input.")
    tool_python = widgets.Text(value="", description="Builder environment (python):",
                               style=adv_style, layout=adv,
                               tooltip="Python interpreter that has the external builder "
                                       "installed; empty = found automatically.")
    shard_size = widgets.BoundedIntText(value=TOOLS["manta"]["shard_size"], min=0, max=100000,
                                        description="Molecules per job (shard size):",
                                        style=adv_style, layout=adv,
                                        tooltip="How many molecules one Slurm array task builds.")
    timeout = widgets.BoundedIntText(value=21600, min=1, max=10 ** 7,
                                     description="Time limit per molecule (s):",
                                     style=adv_style, layout=adv,
                                     tooltip="A molecule still building after this is stopped and "
                                             "recorded as timeout. Not the job wall time.")
    speed = widgets.Text(value="1.0", description="Cluster slowness factor (x):",
                         style=adv_style, layout=adv,
                         tooltip="Multiplies the time limit, e.g. 1.3 on nodes 30 % slower than "
                                 "the reference machine.")
    throttle = widgets.BoundedIntText(value=40, min=0, max=100000,
                                      description="Max. parallel jobs:", style=adv_style,
                                      layout=adv, tooltip="Slurm array tasks running at once "
                                                          "(0 = no limit).")
    repeat = widgets.BoundedIntText(value=0, min=0, max=100000,
                                    description="Determinism check: rebuild N molecules:",
                                    style=adv_style, layout=adv,
                                    tooltip="N molecules are built a second time and compared "
                                            "(0 = off).")
    wall = widgets.Text(value="", placeholder="auto", description="Job wall time (hh:mm:ss):",
                        style=adv_style, layout=adv,
                        tooltip="Slurm time limit of each array task; empty = auto.")
    shard_help, timeout_help, wall_help = _help(), _help(), _help()

    buttons = {name: widgets.Button(description=name, button_style=bs, tooltip=tip,
                                    layout=widgets.Layout(width="110px"))
               for name, bs, tip in (
                   ("Prepare", "info", "Check the list and write the run folder (shards, "
                                       "settings); nothing is built yet."),
                   ("Submit", "success", "Start the prepared run (Slurm array, or one shard "
                                         "after the other here)."),
                   ("Status", "", "Progress of the run, with the details of every shard."),
                   ("Collect", "warning", "Merge the finished shards into one archive."))}
    refresh = widgets.Button(description="Refresh", icon="refresh",
                             tooltip="Update the progress line",
                             layout=widgets.Layout(width="100px"))
    status_html = widgets.HTML("")
    command = widgets.HTML("")
    out = widgets.Output(layout=widgets.Layout(max_height="320px", overflow_y="auto"))
    last_tool = {"tool": "manta"}

    def _cb_count():
        if list_path.value.strip():
            return _cb_count_lines(list_path.value.strip())
        return sum(1 for ln in str(batch_text_widget.value).splitlines() if ln.strip())

    def _cb_refresh(*_):
        t, n = tool.value, _cb_count()
        mode_help.value = _CB_HELP.format(html.escape(CB_MODE_HELP.get((t, mode.value), "")))
        name = run_name.value.strip()
        run_name.placeholder = "auto: " + construction_batch_default_run_name(t, n, runs_root)
        run = runs_root / (name or "_")
        if name and run.exists():
            run_help.value = _CB_WARN.format(
                f"{html.escape(str(run))} exists: Prepare will refuse (a run is never "
                "overwritten); Submit, Status and Collect work on it.")
        else:
            run_help.value = _CB_HELP.format(f"folder {html.escape(str(runs_root))}/&lt;run "
                                             "name&gt;; empty = the automatic name shown")
        if list_path.value.strip():
            list_help.value = (_CB_HELP.format(f"{n} molecule line(s), one ID;SMILES per line")
                               if n is not None else _CB_WARN.format("file not found"))
        else:
            list_help.value = _CB_HELP.format(f"{n} line(s) in the Batch SMILES field above "
                                              "(one ID;SMILES per line)")
        default_size = TOOLS[t]["shard_size"]
        shard_help.value = _CB_HELP.format(f"default for {t}: {default_size}; one shard = one "
                                           "job (Slurm array task)")
        try:
            tmo = cbatch_effective_timeout(int(timeout.value), speed.value.strip() or "1.0")
        except Exception:  # noqa: BLE001 -- shown, not raised
            timeout_help.value = _CB_WARN.format("the slowness factor must be a number, e.g. 1.3")
            plan.value = ""
            return
        timeout_help.value = _CB_HELP.format(
            f"= {int(timeout.value) / 3600:.1f} h; with the slowness factor {tmo / 3600:.1f} h "
            f"({tmo} s); not the job wall time")
        man = construction_batch_manifest(run) if name else None
        p = None
        if man:
            s = man["settings"]
            ns = [e["n"] for v in man["sets"].values() for e in v["shards"]]
            p = construction_batch_plan(ns, s["timeout_s"], s["workers"], throttle.value)
            source, tmo = "prepared run", s["timeout_s"]
        elif n:
            size = int(shard_size.value) or default_size
            ns = [min(size, n - k) for k in range(0, n, size)]
            if int(repeat.value) > 0:
                ns.append(min(int(repeat.value), n, size))
            p = construction_batch_plan(ns, tmo, TOOLS[t]["workers"], throttle.value)
            source = "estimate"
        if p is None:
            wall_help.value = _CB_HELP.format(f"auto = {CB_DEFAULT_WALL_H} h")
            plan.value = ""
            return
        need = p["bound_h"]
        too_long = need > CB_DEFAULT_WALL_H
        wall_help.value = (_CB_WARN if too_long else _CB_HELP).format(
            f"auto = {CB_DEFAULT_WALL_H} h (computed: {need:.1f} h for the biggest shard of "
            f"{p['biggest']})" + (" -- longer than the job: lower the shard size or set a longer "
                                  "wall time (an unfinished shard resumes when submitted again)"
                                  if too_long else ""))
        plan.value = _CB_HELP.format(
            f"<b>{source}:</b> {p['n']} molecules, {p['n_shards']} shard(s), at most "
            f"{p['core_h']:,.0f} core-hours, ~{p['wall_h']:.1f} h wall with {p['parallel']} "
            f"parallel job(s) (if every molecule runs to its {tmo / 3600:.1f} h limit)")

    def _cb_tool_changed(*_):
        t = tool.value
        old_default = TOOLS[last_tool["tool"]]["shard_size"]
        last_tool["tool"] = t
        mode.options = [(f"{m} (default)" if m == TOOLS[t]["mode"] else m, m) for m in MODES[t]]
        mode.value = TOOLS[t]["mode"]
        if int(shard_size.value) in (0, old_default):
            shard_size.value = TOOLS[t]["shard_size"]
        if not TOOLS[t]["needs_specs"]:
            tool_python.placeholder = "not needed: MANTA runs in DELFIN's own interpreter"
            tool_python.disabled = specs_path.disabled = True
        else:
            tool_python.disabled = specs_path.disabled = False
            try:
                from delfin.common import external_builders as eb

                tool_python.placeholder = "auto: " + eb.tool_python(t)
            except Exception:  # noqa: BLE001
                tool_python.placeholder = "auto"
        _cb_refresh()

    tool.observe(_cb_tool_changed, names="value")
    for w in (mode, run_name, list_path, shard_size, timeout, speed, throttle, repeat):
        w.observe(_cb_refresh, names="value")
    if hasattr(batch_text_widget, "observe"):
        batch_text_widget.observe(_cb_refresh, names="value")
    _cb_tool_changed()

    def _cb_fields() -> dict:
        name = run_name.value.strip()
        if not name:
            raise ValueError(f"Enter the run name (a folder under {runs_root}), or press "
                             "Prepare to start a new run with the automatic name.")
        if "/" in name or name.startswith("."):
            raise ValueError("Run name: a plain folder name is needed (no '/', not starting "
                             "with '.').")
        return {"tool": tool.value, "mode": mode.value, "run_dir": runs_root / name,
                "select": select_path.value, "specs": specs_path.value,
                "tool_python": tool_python.value, "shard_size": shard_size.value,
                "timeout": timeout.value, "speed_factor": speed.value, "repeat": repeat.value,
                "throttle": throttle.value, "time_limit": wall.value.strip() or CB_DEFAULT_WALL}

    def _cb_show(argv, rc, text):
        line = construction_batch_command_line(argv)
        command.value = (_CB_HELP.format("same as:")
                         + f" <code style='font-size:12px'>{html.escape(line)}</code>")
        with out:
            print(f"$ {line}")
            if text:
                print(text)
            if rc:
                print(f"(exit {rc})")

    def _cb_prepared(f):
        man = construction_batch_manifest(f["run_dir"])
        if man is None:
            raise ValueError(f"No prepared run {f['run_dir'].name!r} in {runs_root} -- press "
                             "Prepare first, or enter the name of an existing run.")
        return man

    def _cb_status_lines(run_dir, man) -> list:
        from delfin.cluster_bench.report import cbatch_status

        lines = []
        for set_name in man["sets"]:
            try:
                lines.append(construction_batch_status_line(cbatch_status(run_dir, set_name)))
            except SystemExit as exc:
                lines.append(f"{set_name}: {exc}")
        return lines

    def _cb_prepare():
        if not run_name.value.strip():
            run_name.value = construction_batch_default_run_name(tool.value, _cb_count(),
                                                                 runs_root)
            with out:
                print(f"Run name: {run_name.value} (automatic)")
        f = _cb_fields()
        run = f["run_dir"]
        if run.exists():
            raise ValueError(f"The folder {run} already exists. Prepare never overwrites a run: "
                             "choose another run name (Status and Collect work on the "
                             "existing one).")
        src = list_path.value.strip()
        if not src:
            rows, errors = construction_batch_field_rows(str(batch_text_widget.value))
            if errors:
                raise ValueError(f"{len(errors)} invalid line(s) in the Batch SMILES field -- "
                                 "nothing written:\n  " + "\n  ".join(errors[:20])
                                 + "\nEach line: ID;SMILES (further ;key=value fields are "
                                   "ignored).")
            if not rows:
                raise ValueError("The list is empty: enter ID;SMILES lines in the Batch SMILES "
                                 "field above, or give a list file.")
            text = construction_batch_input_from_text(str(batch_text_widget.value))
            src = run.parent / f"{run.name}.input.txt"
            if src.exists() and src.read_text() != text:
                raise ValueError(f"{src} exists with another list; choose a new run name.")
            src.parent.mkdir(parents=True, exist_ok=True)
            src.write_text(text)
        elif not Path(src).is_file():
            raise ValueError(f"List file {src} not found.")
        f["input"] = src
        argv = construction_batch_argv("prepare", f)
        rc, text = construction_batch_call(argv)
        man = construction_batch_manifest(run) if rc == 0 else None
        with out:
            if man:
                print(construction_batch_summary(man, run, throttle.value))
                print("Next: Submit.\n\nDetails:")
            else:
                head = (text.strip().splitlines() or ["see below"])[0]
                print(f"Prepare refused, nothing was written: {head}\n\nDetails:")
        _cb_show(argv, rc, text)
        _cb_refresh()

    def _cb_submit():
        f = _cb_fields()
        sets = list(_cb_prepared(f)["sets"])
        if getattr(ctx.backend, "backend_name", "") == "SlurmJobBackend":
            setup = construction_batch_slurm_setup(ctx.backend)
            for set_name in sets:
                argv = construction_batch_argv("slurm", dict(f, set=set_name, setup=setup))
                _cb_show(argv, *construction_batch_call(argv))
            with out:
                print("Submitted; Status or Refresh shows the progress.")
            return
        with out:
            print("No Slurm here: the shards are built one after the other in the background "
                  f"({', '.join(sets)}); Status or Refresh shows the progress.")
        construction_batch_run_locally(f["run_dir"], sets, lambda m: out.append_stdout(m + "\n"))

    def _cb_set_status(lines):
        status_html.value = "<br>".join(
            f"<code style='font-size:12px'>{html.escape(x)}</code>" for x in lines)

    def _cb_status():
        f = _cb_fields()
        man = _cb_prepared(f)
        lines = _cb_status_lines(f["run_dir"], man)
        _cb_set_status(lines)
        with out:
            print("\n".join(lines) + "\n\nDetails:")
        for set_name in man["sets"]:
            argv = construction_batch_argv("status", dict(f, set=set_name))
            _cb_show(argv, *construction_batch_call(argv))

    def _cb_refresh_status(_button=None):
        try:
            f = _cb_fields()
            _cb_set_status(_cb_status_lines(f["run_dir"], _cb_prepared(f)))
        except Exception as exc:  # noqa: BLE001
            status_html.value = _CB_WARN.format(html.escape(str(exc)))

    def _cb_collect():
        f = _cb_fields()
        man = _cb_prepared(f)
        for set_name in man["sets"]:
            label = man["label"] + ("" if set_name == "main" else f"_{set_name}_main")
            dest = Path(f["run_dir"]).resolve() / "collected" / label
            if dest.exists():
                with out:
                    print(f"{set_name}: already collected into {dest} (Collect never "
                          "overwrites; rename that folder to collect again).")
                continue
            argv = construction_batch_argv("collect", dict(f, set=set_name))
            rc, text = construction_batch_call(argv)
            try:
                s = json.loads(text)
            except ValueError:
                s = None
            with out:
                if isinstance(s, dict) and "n_xyz" in s:
                    left = (f"; shards not finished: {s['resubmit_array']} (Submit again, they "
                            "resume)") if s.get("resubmit_array") else ""
                    print(f"{set_name}: results in {dest / ('archive_' + s['label'])} -- "
                          f"{s['n_xyz']} of {s['n_systems_in_set']} molecules have frames"
                          f"{left}.")
                else:
                    head = (text.strip().splitlines() or ["see below"])[0]
                    print(f"{set_name}: Collect refused: {head}")
            _cb_show(argv, rc, text)
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
    refresh.on_click(_cb_refresh_status)

    advanced_body = widgets.VBox([
        _row(select_path, _help("a file with one ID per line")),
        _row(specs_path, _help("external builders only")),
        _row(tool_python, _help("external builders only")),
        _row(shard_size, shard_help),
        _row(timeout, timeout_help),
        _row(speed, _help("1.0 = the reference machine")),
        _row(throttle, _help("Slurm array tasks at once; 0 = no limit")),
        _row(repeat, _help("0 = off; rebuilt as set 'repeat', compared by Collect")),
        _row(wall, wall_help),
    ])
    advanced = widgets.Accordion(children=[advanced_body])
    advanced.set_title(0, "Advanced")
    advanced.selected_index = None
    basic = widgets.VBox([
        widgets.HTML(_CB_HELP.format(
            "Build many molecules with one builder, sharded and resumable: Prepare, Submit, "
            "then Status and Collect (docs/CONSTRUCTION_BATCH.md).")),
        _row(tool, _help("the program that builds the structures")),
        _row(mode, mode_help),
        _row(run_name, run_help),
        list_path, list_help, plan,
    ])
    body = widgets.VBox([
        basic, advanced,
        widgets.HBox(list(buttons.values()) + [refresh],
                     layout=widgets.Layout(gap="8px", flex_wrap="wrap")),
        status_html, command, out,
    ])
    acc = widgets.Accordion(children=[body])
    acc.set_title(0, "Construction batch (MANTA / ARCHITECTOR / MOLSIMPLIFY / MACE)")
    acc.selected_index = None
    acc.construction_batch_widgets = {
        "tool": tool, "mode": mode, "run_name": run_name, "list_path": list_path,
        "select_path": select_path, "specs_path": specs_path, "tool_python": tool_python,
        "shard_size": shard_size, "timeout": timeout, "speed": speed, "throttle": throttle,
        "repeat": repeat, "wall": wall, "buttons": buttons, "output": out, "command": command,
        "basic": basic, "advanced": advanced, "plan": plan, "status": status_html,
        "refresh": refresh, "mode_help": mode_help, "run_help": run_help,
        "timeout_help": timeout_help, "shard_help": shard_help, "wall_help": wall_help}
    return acc
