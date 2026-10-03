"""Adversarial reviewer tests for package F (cluster-bench node-local I/O).

Reviewer s4, wave 11.  These try to break the builder's fix (agent/s3-f11,
commits 7eca02d2 + 783db151) from outside; they do not repeat its own tests.

The three staging tests EXECUTE the bash the generated sbatch would run --
string assertions cannot catch a path the script never creates.  Everything
else is faked: a tiny stdlib worker instead of torch/MACE, tmp_path venvs, no
SLURM.  No test depends on wall-clock speed: the one wait loop a worker has is
bounded and only orders events, never times them.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.cluster_bench import prepare as cp
from delfin.cluster_bench import provenance as cprov
from delfin.cluster_bench import report as crep
from delfin.cluster_bench import runner as crun
from delfin.cluster_bench.slurm_script import cbatch_render_sbatch

CISPLATIN = "[Cl][Pt-2]([Cl])([NH3+])[NH3+]"


def _man(tool_python=None) -> dict:
    return {"tool": "architector", "label": "adv",
            "settings": {"tool_python": tool_python, "timeout_base_s": 21600,
                         "speed_factor": "1.0", "workers": 48, "threads": 1},
            "sets": {"main": {"n_shards": 2}}}


def _staging_block(text: str) -> str:
    """The staging part of a rendered sbatch, from its marker to the final joint export."""
    lines = text.splitlines()
    start = next(i for i, ln in enumerate(lines) if "node-local I/O" in ln)
    end = next(i for i, ln in enumerate(lines)
               if ln.startswith("export DELFIN_CLUSTER_LOCAL_ROOT DELFIN_CLUSTER_TOOL_PYTHON"))
    return "\n".join(lines[start:end + 1])


def _run_staging(tmp_path: Path, tool_python, tmpdir):
    """Execute the staging block under the sbatch's own `set -euo pipefail`,
    with a real tiny venv to pack; return (result, exported vars)."""
    venv = tmp_path / "ws_venv"
    (venv / "bin").mkdir(parents=True, exist_ok=True)
    (venv / "bin" / "python").write_text("#!/bin/sh\nexit 0\n")
    (venv / "bin" / "python").chmod(0o755)
    (venv / "pyvenv.cfg").write_text("home = /nowhere\n")
    site = venv / "lib" / "python3.9" / "site-packages"
    site.mkdir(parents=True, exist_ok=True)
    (site / "pkg_module.py").write_text("x = 1\n")
    text = cbatch_render_sbatch(tmp_path / "run", _man(tool_python=str(venv / "bin" / "python")))
    script = tmp_path / "staging.sh"
    script.write_text("#!/bin/bash\nset -euo pipefail\n" + _staging_block(text)
                      + "\nprintf 'LOCAL_ROOT=%s\\nTOOL_PYTHON=%s\\n' "
                        "\"$DELFIN_CLUSTER_LOCAL_ROOT\" \"$DELFIN_CLUSTER_TOOL_PYTHON\"\n")
    env = dict(os.environ)
    env["DELFIN_VENV_CACHE_DIR"] = str(tmp_path / "cache")
    if tmpdir is None:
        env.pop("TMPDIR", None)
    else:
        env["TMPDIR"] = tmpdir
    done = subprocess.run(["bash", str(script)], capture_output=True, text=True, env=env,
                          timeout=120)
    out = dict(ln.split("=", 1) for ln in done.stdout.splitlines() if "=" in ln)
    return done, out


def test_the_staging_block_exports_a_python_that_exists_and_runs(tmp_path):
    """The exported staged interpreter must be a file that actually exists: children are
    spawned with it 48 x thousands of times, so a wrong path kills every array task."""
    nodetmp = tmp_path / "nodetmp"
    nodetmp.mkdir()
    done, out = _run_staging(tmp_path, None, str(nodetmp))
    assert done.returncode == 0, done.stderr
    tp = out.get("TOOL_PYTHON", "")
    assert tp, "no staged interpreter exported"
    assert tp != str(tmp_path / "ws_venv" / "bin" / "python"), "exported the workspace venv"
    assert os.path.isfile(tp) and os.access(tp, os.X_OK), f"staged python not usable: {tp}"
    lr = Path(out["LOCAL_ROOT"])
    assert lr.is_dir()
    assert list(lr.rglob("pkg_module.py")), "venv content not unpacked under the local root"


def test_the_staging_block_falls_back_when_tmpdir_is_unset_not_crash(tmp_path):
    """The script runs under `set -euo pipefail`; with TMPDIR unset (no SLURM default, a
    bare node) `[ -n \"$TMPDIR\" ]` is an unbound variable: the block must keep today's
    behaviour with one warning line, not kill the array task."""
    done, out = _run_staging(tmp_path, None, None)
    assert done.returncode == 0, done.stderr
    assert out.get("TOOL_PYTHON", "") == ""
    assert "WARNING" in done.stdout


def _fake_workers_dir(tmp_path: Path) -> Path:
    """A tiny stdlib worker recording the archive/work/log paths it was given; writes an
    xyz + final meta.  `CB_FAIL_IDS` optionally makes it exit 7 instead (the runner sees
    rc=7 and records a crash -- no `final` meta for those)."""
    wd = tmp_path / "workers"
    wd.mkdir(exist_ok=True)
    (wd / "env_probe.py").write_text("import json\nprint(json.dumps({'packages': {}}))\n")
    (wd / "bench_worker.py").write_text(
        "import json, os, sys\n"
        "tool, ref, specs, archive, work = sys.argv[1:6]\n"
        "os.makedirs(os.path.join(archive, '_meta'), exist_ok=True)\n"
        "os.makedirs(os.path.join(work), exist_ok=True)\n"
        "open(os.environ['CB_PATHS'], 'a').write('%s\\t%s\\t%s\\n' % (archive, work, 'LOG'))\n"
        "open(os.path.join(work, 'scratch.bin'), 'wb').write(b'x' * 8)\n"
        "if ref in os.environ.get('CB_FAIL_IDS', '').split(','):\n"
        "    sys.exit(7)\n"
        "open(os.path.join(archive, ref + '.xyz'), 'w').write('1\\n%s frame0\\nFe 0 0 0\\n' % ref)\n"
        "json.dump({'refcode': ref, 'tool': tool, 'status': 'ok', 'n_frames': 1},\n"
        "          open(os.path.join(archive, '_meta', ref + '.json'), 'w'))\n")
    return wd


def _prepared_run(tmp_path, monkeypatch, tool_python, ids, fail_ids=(), size=50):
    wd = _fake_workers_dir(tmp_path)
    monkeypatch.setattr(cprov, "WORKERS_DIR", wd)
    monkeypatch.setattr(crun, "WORKERS_DIR", wd)
    inp = tmp_path / "in.txt"
    inp.write_text("".join(f"{rid};{CISPLATIN}\n" for rid in ids))
    specs = tmp_path / "specs.jsonl"
    specs.write_text("".join(json.dumps({"refcode": rid, "status": "ok", "cn": 4}) + "\n"
                             for rid in ids))
    cp.cbatch_prepare(tool="molsimplify", input_list=inp, run_dir=tmp_path / "run",
                      specs_file=specs, tool_python=tool_python, timeout_base=60,
                      size=size, check_tool=False)
    return tmp_path / "run"


def _interpreter_wrapper(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nprintf '%s\\t%s\\n' \"$0\" \"$1\" >> \"$CB_INTERP_MARKER\"\n"
                    "exec \"$CB_REALPY\" \"$@\"\n")
    path.chmod(0o755)
    return path


def _env(monkeypatch, tmp_path, tool_python):
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    _interpreter_wrapper(Path(tool_python))
    return tool_python


def test_exported_staged_python_and_runner_choice_agree(tmp_path, monkeypatch):
    """The path the staging block exports must be the path the runner would use: one wrong
    constant here sends 48 workers x thousands of systems at a nonexistent interpreter."""
    nodetmp = tmp_path / "nodetmp"
    nodetmp.mkdir()
    done, out = _run_staging(tmp_path, None, str(nodetmp))
    assert done.returncode == 0, done.stderr
    exported = out["TOOL_PYTHON"]
    lr = Path(out["LOCAL_ROOT"])
    assert lr.is_dir() and exported.startswith(str(lr)), \
        f"exported interpreter {exported} not under the exported local root"
    # ... and the runner picks up exactly this var
    man = {"settings": {"tool_python": "/ws/manifest/bin/python"}}
    monkeypatch.setenv("DELFIN_CLUSTER_TOOL_PYTHON", exported)
    assert crun.cbatch_tool_python(man) == exported


def test_children_run_from_the_staged_interpreter_and_everything_lies_under_the_local_root(
        tmp_path, monkeypatch):
    """With both vars set (as the staging block exports them), children must run from the
    staged interpreter and write work/log/archive under the local root -- nothing per-system
    on the workspace.  An ok system's log may stay local (it dies with the node)."""
    venv_python = str(tmp_path / "manifest_venv" / "bin" / "python")
    _env(monkeypatch, tmp_path, venv_python)
    staged = tmp_path / "nodetmp" / "cb_x" / "bin" / "python"
    _interpreter_wrapper(staged)
    monkeypatch.setenv("DELFIN_CLUSTER_TOOL_PYTHON", str(staged))
    local = tmp_path / "local"
    monkeypatch.setenv("DELFIN_CLUSTER_LOCAL_ROOT", str(local))
    run_dir = _prepared_run(tmp_path, monkeypatch, venv_python, ["M1", "M2"])

    assert crun.cbatch_run_shard(run_dir, 0, log=lambda m: None) == 0

    chunk = Path(run_dir) / "out" / "main_main" / "chunk_0000"
    recs = [ln.split("\t") for ln in (tmp_path / "paths").read_text().splitlines()]
    for archive, work, _log in recs:
        assert Path(archive) == local / "archive"
        assert Path(work).is_relative_to(local / "work")
    assert (chunk / "work").exists() is False
    assert (chunk / "logs").exists() is False
    # results copied back: every ok system's xyz + final meta on the workspace archive
    for rid in ("M1", "M2"):
        assert (chunk / "archive" / "_meta" / f"{rid}.json").exists()
        assert json.loads((chunk / "archive" / "_meta" / f"{rid}.json").read_text())["final"]
    s = crep.cbatch_collect(run_dir, set_name="main", run_name="main")
    assert s["by_class"] == {"ok": 2}


def test_a_failed_systems_log_survives_on_the_workspace_an_ok_ones_may_not(tmp_path, monkeypatch):
    """Operator review note: the per-system log is the only diagnosis of a fail/timeout, and
    it now dies with the node -- so the batched copy-back must carry NON-ok systems' logs to
    the workspace chunk, while ok logs stay local (they would only multiply files there)."""
    venv_python = str(tmp_path / "manifest_venv" / "bin" / "python")
    _env(monkeypatch, tmp_path, venv_python)
    local = tmp_path / "local"
    monkeypatch.setenv("DELFIN_CLUSTER_LOCAL_ROOT", str(local))
    run_dir = _prepared_run(tmp_path, monkeypatch, venv_python, ["OK1", "BAD1"])
    monkeypatch.setenv("CB_FAIL_IDS", "BAD1")

    assert crun.cbatch_run_shard(run_dir, 0, log=lambda m: None) == 0

    chunk = Path(run_dir) / "out" / "main_main" / "chunk_0000"
    assert json.loads((chunk / "archive" / "_meta" / "BAD1.json").read_text())[
        "status"].startswith("crash"), "the child's rc=7 must be visible in the final record"
    assert (chunk / "logs" / "BAD1.log").exists(), \
        "the failed system's log vanished with the node-local root"
    assert (chunk / "logs" / "OK1.log").exists() is False, \
        "the ok system's log was multiplied onto the workspace"
