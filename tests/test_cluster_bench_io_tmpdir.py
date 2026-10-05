"""Package F: the cluster benchmark must not hammer the shared file system.

New, pinned-behaviour tests for the node-local I/O a cluster run now does:
* the generated sbatch script stages the tool's venv once per job onto ``$TMPDIR`` and
  exports ``DELFIN_CLUSTER_TOOL_PYTHON`` pointing at the staged ``bin/python``; when
  ``$TMPDIR`` is missing or too small it falls back to today's behaviour with one warning;
* the runner runs every tool child from the staged interpreter when the job staged one,
  from the manifest's interpreter otherwise;
* per-system ``work/``, ``logs/`` and the *build* archive live under a node-local root
  (``DELFIN_CLUSTER_LOCAL_ROOT``, a ``$TMPDIR`` subdir) when the job staged one; finished
  systems are copied back to the workspace chunk archive in batches so the on-disk archive
  contract (``archive/<rid>.xyz`` + ``_meta/<rid>.json``) that resume and ``cbatch_collect``
  read is preserved;
* without a staged local root the runner keeps today's per-chunk layout exactly;
* resume after a wall-time kill still skips finished systems (final records on the workspace).

Everything is faked: no real torch/MACE venv is packed or unpacked, no SLURM is touched, and
the tool child is a tiny stdlib worker.  No test depends on wall-clock time.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from delfin.cluster_bench import prepare as cp
from delfin.cluster_bench import provenance as cprov
from delfin.cluster_bench import report as crep
from delfin.cluster_bench import runner as crun
from delfin.cluster_bench.slurm_script import cbatch_render_sbatch

CISPLATIN = "O[Pt](N)(N)(Cl)Cl"
FAKE_ENV_PROBE = (
    "import json, sys\n"
    "print(json.dumps({'python': sys.version_info[:3], 'executable': sys.executable,"
    " 'packages': {}}))\n"
)


def _sbatch_man(**settings):
    s = {"timeout_base_s": 21600, "speed_factor": "1.0", "workers": 48, "threads": 1}
    s.update(settings)
    return {"tool": "architector", "label": "arch_x",
            "settings": s, "sets": {"main": {"n_shards": 2}}}


# -------------------------------------------------------------------- the sbatch stages an env
def test_the_sbatch_stages_a_cached_venv_onto_tmpdir_once_per_job(tmp_path):
    text = cbatch_render_sbatch(tmp_path, _sbatch_man(tool_python="/ws/venv/bin/python"))
    assert "venv_cache_key" in text or "ensure_venv_tar" in text   # the shared cache logic
    assert "${TMPDIR}" in text                              # node-local root under $TMPDIR
    assert 'DELFIN_CLUSTER_LOCAL_ROOT="${TMPDIR}' in text   # but never triggers set -u (guarded)
    assert "export DELFIN_CLUSTER_TOOL_PYTHON=" in text            # children use the staged copy
    assert "SLURM_JOB_ID" in text                                  # unique per job


def test_the_staged_python_point_to_the_unpacked_copy_not_the_workspace_venv(tmp_path):
    text = cbatch_render_sbatch(tmp_path, _sbatch_man(tool_python="/ws/venv/bin/python"))
    assert 'DELFIN_CLUSTER_TOOL_PYTHON="$DELFIN_CLUSTER_LOCAL_ROOT/venv/bin/python"' in text
    assert "TOOL_PYTHON=/ws/venv/bin/python" in text   # the venv to pack is named, from the manifest


def test_the_rendered_staging_block_stages_and_runs_a_fake_venv(tmp_path):
    """EXECUTES the emitted bash (not a string assert): a tiny fake venv packed onto a temp
    TMPDIR must produce an interpreter that exists, is executable, is still a symlink (no
    sed over bin/*), and keeps the venv tree separate from work/logs/archive."""
    from delfin.cluster_bench.slurm_script import _cbatch_staging_lines

    venv = tmp_path / "toolenv"
    (venv / "bin").mkdir(parents=True)
    (venv / "bin" / "python").symlink_to(sys.executable)          # real venvs symlink the interpreter
    (venv / "bin" / "activate").write_text('export VIRTUAL_ENV="x"\n')
    (venv / "pyvenv.cfg").write_text("home = %s\n" % sys.prefix)
    (venv / "lib" / "python3.12" / "site-packages").mkdir(parents=True)
    (venv / "lib" / "python3.12" / "site-packages" / "marker.txt").write_text("x")

    node = tmp_path / "nodetmp"; node.mkdir()
    cache = tmp_path / "venvcache"; cache.mkdir()
    expr = tmp_path / "expr"; root = tmp_path / "root"

    block = "\n".join(_cbatch_staging_lines(str(venv / "bin" / "python")))
    script = (
        "set -euo pipefail\n"
        "export TMPDIR=%s\n"
        "export DELFIN_VENV_CACHE_DIR=%s\n"
        "%s\n"
        "printf '%%s\\n' \"$DELFIN_CLUSTER_TOOL_PYTHON\" > %s\n"
        "printf '%%s\\n' \"$DELFIN_CLUSTER_LOCAL_ROOT\" > %s\n"
    ) % (node, cache, block, expr, root)
    sh_file = tmp_path / "stage.sh"
    sh_file.write_text(script)

    res = subprocess.run(["bash", str(sh_file)], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    tp = Path(expr.read_text().strip())
    assert tp.exists() and os.access(tp, os.X_OK), tp          # staged interpreter present + executable
    assert tp.is_symlink() and tp.resolve() == Path(sys.executable).resolve(), \
        "bin/python must stay a symlink (no sed over bin/*)"
    lr = Path(root.read_text().strip())
    assert (lr / "venv" / "bin" / "python").exists()
    # the venv tree is its own subdir: no runner dirs inside it
    for sub in ("work", "logs", "archive"):
        assert not (lr / "venv" / sub).exists(), sub


def test_the_sbatch_falls_back_with_a_warning_when_tmpdir_is_unusable(tmp_path):
    text = cbatch_render_sbatch(tmp_path, _sbatch_man(tool_python="/ws/venv/bin/python"))
    assert "WARNING:" in text
    # the fallback leaves the exports unset so the runner uses the workspace env
    assert 'DELFIN_CLUSTER_TOOL_PYTHON=""' in text
    assert "shared file system" in text.lower() or "workspace" in text.lower()


def test_the_sbatch_renders_no_staging_when_the_run_has_no_tool_python(tmp_path):
    text = cbatch_render_sbatch(tmp_path, _sbatch_man())
    assert "venv_cache_key" not in text and "DELFIN_CLUSTER_TOOL_PYTHON" not in text


# ----------------------------------------------------------- fake tool environment + a run dir
def _interpreter_wrapper(path):
    """A fake ``bin/python``: appends ``<invoked-path>\t<argv1>`` (the interpreter path the
    runner chose and the script it ran) to a marker file, then runs the real interpreter so the
    stdlib worker executes.  The provenance probe runs the manifest interpreter against
    ``env_probe.py``; the tool child runs the chosen interpreter against ``bench_worker.py`` --
    the tests filter the child lines by that script name."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nprintf '%s\\t%s\\n' \"$0\" \"$1\" >> \"$CB_INTERP_MARKER\"\n"
                    "exec \"$CB_REALPY\" \"$@\"\n")
    path.chmod(0o755)
    return path


def _child_interpreters(interp_marker):
    """The interpreter path each tool child ran from (probe lines filtered out)."""
    out = set()
    for ln in interp_marker.read_text().splitlines():
        interp, script = ln.split("\t")
        if script.rsplit("/", 1)[-1] == "bench_worker.py":
            out.add(interp)
    return out


def _fake_workers_dir(tmp_path):
    wd = tmp_path / "workers"
    wd.mkdir()
    (wd / "env_probe.py").write_text(FAKE_ENV_PROBE)
    (wd / "bench_worker.py").write_text(
        "import json, os, sys\n"
        "tool, ref, specs, archive, work = sys.argv[1:6]\n"
        "os.makedirs(os.path.join(archive, '_meta'), exist_ok=True)\n"
        "open(os.environ['CB_PATHS'], 'a').write('%s\\t%s\\n' % (archive, work))\n"
        "open(os.path.join(archive, ref + '.xyz'), 'w').write('1\\n%s frame0\\nFe 0 0 0\\n' % ref)\n"
        "json.dump({'refcode': ref, 'tool': tool, 'status': 'ok', 'n_frames': 1},\n"
        "          open(os.path.join(archive, '_meta', ref + '.json'), 'w'))\n")
    return wd


def _prepared_run(tmp_path, monkeypatch, tool_python, ids, repeat=0, size=50):
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
                      size=size, repeat=repeat, check_tool=False)
    return tmp_path / "run"


# --------------------------------- the runner runs children from the staged interpreter
def test_a_staged_interpreter_runs_the_children_when_set(tmp_path, monkeypatch):
    interp_marker = tmp_path / "interp"
    monkeypatch.setenv("CB_INTERP_MARKER", str(interp_marker))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    venv_python = tmp_path / "manifest_venv" / "bin" / "python"   # the manifest's interpreter
    _interpreter_wrapper(venv_python)
    staged_python = tmp_path / "staged_env" / "bin" / "python"     # the staged copy
    _interpreter_wrapper(staged_python)
    run_dir = _prepared_run(tmp_path, monkeypatch, str(venv_python), ["M1", "M2"])
    monkeypatch.setenv("DELFIN_CLUSTER_TOOL_PYTHON", str(staged_python))

    assert crun.cbatch_run_shard(run_dir, 0, log=lambda m: None) == 0

    assert _child_interpreters(interp_marker) == {str(staged_python)}  # children ran staged


def test_the_manifest_interpreter_runs_when_no_staged_one_is_set(tmp_path, monkeypatch):
    interp_marker = tmp_path / "interp"
    monkeypatch.setenv("CB_INTERP_MARKER", str(interp_marker))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    venv_python = tmp_path / "manifest_venv" / "bin" / "python"
    _interpreter_wrapper(venv_python)
    staged_python = tmp_path / "staged_env" / "bin" / "python"
    _interpreter_wrapper(staged_python)
    run_dir = _prepared_run(tmp_path, monkeypatch, str(venv_python), ["M1", "M2"])
    monkeypatch.delenv("DELFIN_CLUSTER_TOOL_PYTHON", raising=False)

    assert crun.cbatch_run_shard(run_dir, 0, log=lambda m: None) == 0

    assert _child_interpreters(interp_marker) == {str(venv_python)}  # fallback = manifest interp


def test_tool_python_helper_prefers_the_staged_copy_when_set(monkeypatch):
    man = {"settings": {"tool_python": "/ws/venv/bin/python"}}
    monkeypatch.delenv("DELFIN_CLUSTER_TOOL_PYTHON", raising=False)
    assert crun.cbatch_tool_python(man) == "/ws/venv/bin/python"
    monkeypatch.setenv("DELFIN_CLUSTER_TOOL_PYTHON", "/tmp/staged/bin/python")
    assert crun.cbatch_tool_python(man) == "/tmp/staged/bin/python"


# ------------------------------------------------- work/log/build-archive under the local root
def test_work_logs_and_build_archive_lie_under_the_local_root_when_staged(tmp_path, monkeypatch):
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    venv_python = tmp_path / "manifest_venv" / "bin" / "python"
    _interpreter_wrapper(venv_python)
    run_dir = _prepared_run(tmp_path, monkeypatch, str(venv_python), ["M1", "M2"])
    local = tmp_path / "local"
    monkeypatch.setenv("DELFIN_CLUSTER_LOCAL_ROOT", str(local))

    assert crun.cbatch_run_shard(run_dir, 0, log=lambda m: None) == 0

    paths = [ln.split("\t") for ln in (tmp_path / "paths").read_text().splitlines()]
    ws_archive = Path(run_dir) / "out" / "main_main" / "chunk_0000" / "archive"
    for archive, work in paths:
        assert Path(archive).is_relative_to(local / "archive")     # child wrote into the local build
        assert Path(work).is_relative_to(local / "work")           # scratch under the local root
    assert (local / "logs" / "M1.log").exists()                     # per-system log is local
    # the workspace chunk archive still ends up complete via the batched copy-back
    for rid in ("M1", "M2"):
        assert (ws_archive / f"{rid}.xyz").exists()
        assert (ws_archive / "_meta" / f"{rid}.json").exists()


def test_work_logs_and_archive_stay_on_the_chunk_without_a_local_root(tmp_path, monkeypatch):
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    venv_python = tmp_path / "manifest_venv" / "bin" / "python"
    _interpreter_wrapper(venv_python)
    run_dir = _prepared_run(tmp_path, monkeypatch, str(venv_python), ["M1", "M2"])
    monkeypatch.delenv("DELFIN_CLUSTER_LOCAL_ROOT", raising=False)

    assert crun.cbatch_run_shard(run_dir, 0, log=lambda m: None) == 0

    chunk = Path(run_dir) / "out" / "main_main" / "chunk_0000"
    paths = [ln.split("\t") for ln in (tmp_path / "paths").read_text().splitlines()]
    for archive, work in paths:
        assert Path(archive) == chunk / "archive"    # today's per-chunk archive
        assert Path(work).is_relative_to(chunk / "work")
    assert (chunk / "logs" / "M1.log").exists()
    assert (chunk / "archive" / "_meta" / "M1.json").exists()


# ----------------------------------------------- the chunk copies finished results back (batched)
def test_finished_systems_are_collected_from_the_workspace_archive_after_a_local_run(
        tmp_path, monkeypatch):
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    venv_python = tmp_path / "manifest_venv" / "bin" / "python"
    _interpreter_wrapper(venv_python)
    run_dir = _prepared_run(tmp_path, monkeypatch, str(venv_python), ["M1", "M2"])
    local = tmp_path / "local"
    monkeypatch.setenv("DELFIN_CLUSTER_LOCAL_ROOT", str(local))

    assert crun.cbatch_run_shard(run_dir, 0, log=lambda m: None) == 0

    # cbatch_collect (report.py, unchanged) reads the workspace archive and sees every system
    s = crep.cbatch_collect(run_dir, set_name="main", run_name="main")
    assert s["n_problems"] == 0, s["problems"]
    assert s["by_class"] == {"ok": 2}


# ------------------------------------------------------------- resume after a kill skips finals
def test_resume_after_a_kill_skips_finished_systems_with_a_local_root(tmp_path, monkeypatch):
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    venv_python = tmp_path / "manifest_venv" / "bin" / "python"
    _interpreter_wrapper(venv_python)
    run_dir = _prepared_run(tmp_path, monkeypatch, str(venv_python), ["M1", "M2"])
    local = tmp_path / "local"
    monkeypatch.setenv("DELFIN_CLUSTER_LOCAL_ROOT", str(local))
    chunk = Path(run_dir) / "out" / "main_main" / "chunk_0000"
    chunk.mkdir(parents=True, exist_ok=True)
    (chunk / "archive" / "_meta").mkdir(parents=True, exist_ok=True)
    # M1 already has a final record on the workspace archive (built before the kill)
    (chunk / "archive" / "_meta" / "M1.json").write_text(
        json.dumps({"refcode": "M1", "final": True, "status": "ok"}))

    assert crun.cbatch_run_shard(run_dir, 0, log=lambda m: None) == 0

    # only M2 was rebuilt: the child never ran for M1 (its path is not in the marker)
    build = (chunk / "archive").read_text() if (chunk / "archive").is_file() else ""
    ws_meta = chunk / "archive" / "_meta" / "M2.json"
    assert ws_meta.exists() and json.loads(ws_meta.read_text()).get("final") is True
    assert (local / "work" / "M1") .exists() is False    # M1 was not rebuilt on the local root
    assert (local / "logs" / "M1.log").exists() is False


def test_empty_shard_and_none_workers_still_work(tmp_path, monkeypatch):
    monkeypatch.setenv("CB_INTERP_MARKER", str(tmp_path / "interp"))
    monkeypatch.setenv("CB_REALPY", sys.executable)
    monkeypatch.setenv("CB_PATHS", str(tmp_path / "paths"))
    venv_python = tmp_path / "manifest_venv" / "bin" / "python"
    _interpreter_wrapper(venv_python)
    # 100 systems, workers=None -> the manifest default (48) is used; must all finish and collect
    run_dir = _prepared_run(tmp_path, monkeypatch, str(venv_python),
                            [f"S{i:03d}" for i in range(100)], size=1000)
    monkeypatch.setenv("DELFIN_CLUSTER_TOOL_PYTHON", str(venv_python))
    assert crun.cbatch_run_shard(run_dir, 0, workers=None, log=lambda m: None) == 0
    s = crep.cbatch_collect(run_dir, set_name="main", run_name="main")
    assert s["by_class"] == {"ok": 100}
