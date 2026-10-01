"""The worker of one array task: build every system of one shard, one child process per system.

Common to all tools:
* the shard (and its specs) must match the sha256 of the manifest, and the code state and the
  environments must match the provenance recorded by ``prepare`` -- otherwise nothing is built
  (exit 2).  A run whose chunks were built with different code is not one result;
* per-system wall limit = ceil(timeout_base x speed_factor); at the limit the whole process
  group of that system is killed and the system is recorded as ``timeout``.  The limit decides
  only WHICH systems finish, never the content of a finished one;
* resumable: a system with a final record is skipped, so resubmitting the same array index
  continues a task that hit the job's wall time or lost its node;
* exit 0 when every system of the shard has a final record, 4 otherwise (resubmit the index).

MANTA (``manta_child.py``): IDs with identical SMILES are built once, sequentially in one worker;
the others get a copy of the frames with their own ID in the header (``cached``).  A build that
did not end ``ok`` is not served to the others -- each is then built itself.
External builders (``tool_workers/``): run with the tool's own interpreter; one worker process
per system writes ``<ID>.xyz`` and ``_meta/<ID>.json``; the runner adds wall time, limit and peak
RSS and marks the record final.  Start order: highest CN first (shorter tail; no effect on any
output).
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from delfin.cluster_bench.prepare import cbatch_chunk_dir, cbatch_load_manifest, cbatch_shard_rows
from delfin.cluster_bench.provenance import (HERE, REPO_ROOT, WORKERS_DIR, cbatch_effective_timeout,
                                             cbatch_provenance, cbatch_provenance_mismatch,
                                             cbatch_sha256_file, cbatch_tool_child_env)

_SANITISE_PREFIXES = ("DELFIN_", "WEDDELL_", "LOOP_", "EYE_")
_SANITISE_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                  "NUMEXPR_NUM_THREADS", "PYTHONHASHSEED", "PYTHONPATH", "PYTHONSTARTUP")
MANTA_STATUSES = ("ok", "empty", "timeout", "fail")


def cbatch_manta_env() -> dict:
    """The build environment: inherited variables that could change a build are removed
    (every DELFIN_/WEDDELL_/LOOP_/EYE_ switch, thread counts, hash seed, PYTHONPATH), then the
    fixed settings.  The construction switches come from the child itself."""
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(_SANITISE_PREFIXES) and k not in _SANITISE_KEYS}
    env["DELFIN_DETERMINISTIC"] = "1"
    env["PYTHONHASHSEED"] = "0"
    env["PYTHONNOUSERSITE"] = "1"
    for t in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[t] = "1"
    return env


def cbatch_run_child(cmd, env, timeout):
    """-> (returncode, stdout, stderr, timed_out); at the limit the process group is killed."""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
                         start_new_session=True)
    try:
        out, err = p.communicate(timeout=timeout)
        return p.returncode, out, err, False
    except subprocess.TimeoutExpired:
        try:
            os.killpg(p.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        p.communicate()
        return None, "", "", True


def cbatch_relabel_xyz(src: Path, dst: Path, rid: str) -> None:
    """Copy a multi-frame xyz with ``rid`` as the first word of every comment line."""
    lines = src.read_text().splitlines()
    parts, i = [], 0
    while i < len(lines):
        try:
            n = int(lines[i].split()[0])
        except (ValueError, IndexError):
            i += 1
            continue
        hdr = lines[i + 1].split(None, 1) if i + 1 < len(lines) else []
        rest = hdr[1] if len(hdr) > 1 else ""
        parts.append(f"{n}\n{rid} {rest}\n" + "\n".join(lines[i + 2:i + 2 + n]) + "\n")
        i += 2 + n
    tmp = dst.with_name(f".{dst.name}.part")
    tmp.write_text("".join(parts))
    os.replace(tmp, dst)


def cbatch_read_status(path: Path) -> dict:
    done = {}
    if path.exists():
        for ln in path.read_text().splitlines():
            try:
                rec = json.loads(ln)
                done[rec["rid"]] = rec
            except (ValueError, KeyError):
                pass
    return done


def _manta_shard(man, rows, chunk: Path, workers, timeout, log) -> dict:
    s = man["settings"]
    archive = chunk / "archive"
    archive.mkdir(parents=True, exist_ok=True)
    (chunk / "logs").mkdir(exist_ok=True)
    jl = chunk / "status.jsonl"
    done = cbatch_read_status(jl)
    env = cbatch_manta_env()
    child = str(HERE / "manta_child.py")
    lock = threading.Lock()

    def _record_status(rec):
        with lock:
            with open(jl, "a") as fh:
                fh.write(json.dumps(rec, sort_keys=True) + "\n")
            done[rec["rid"]] = rec
        log(f"{rec['rid']:<12} {rec['status']:<8} {rec.get('time_s', 0):9.1f}s"
            + (" cached" if rec.get("cached") else ""))

    def _build_one_system(rid, smi):
        t0 = time.time()
        rec = {"rid": rid, "status": "fail"}
        cmd = [s["tool_python"], child, str(REPO_ROOT), rid, smi, str(archive), s["mode"],
               str(s["threads"])]
        try:
            rc, so, se, to = cbatch_run_child(cmd, env, timeout)
            if to:
                rec["status"] = "timeout"
            else:
                for ln in so.splitlines():
                    if ln.startswith("{"):
                        try:
                            j = json.loads(ln)
                        except ValueError:
                            continue
                        rec.update({k: j[k] for k in ("niso", "peak_gb") if k in j})
                        rec["status"] = j.get("status", "fail")
                        break
                if rec["status"] not in MANTA_STATUSES:
                    rec["status"] = "fail"
                if rec["status"] == "fail":
                    rec["returncode"] = rc
                    (chunk / "logs" / f"{rid}.log").write_text((se or "") + "\n--- STDOUT ---\n" + (so or ""))
        except Exception as e:  # noqa: BLE001
            rec["error"] = f"{type(e).__name__}: {e}"
        rec["time_s"] = round(time.time() - t0, 2)
        _record_status(rec)
        return rec

    def _build_smiles_group(smi, rids):
        first = None
        for rid in rids:
            if rid in done:
                if first is None and done[rid]["status"] == "ok" and (archive / f"{rid}.xyz").exists():
                    first = rid
                continue
            if first is not None:
                t0 = time.time()
                cbatch_relabel_xyz(archive / f"{first}.xyz", archive / f"{rid}.xyz", rid)
                _record_status({"rid": rid, "status": "ok", "cached": 1, "served_from": first,
                                "time_s": round(time.time() - t0, 2)})
                continue
            if _build_one_system(rid, smi)["status"] == "ok":
                first = rid

    groups: dict = {}
    for rid, smi in rows:
        groups.setdefault(smi, []).append(rid)
    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        for f in as_completed([ex.submit(_build_smiles_group, smi, rids) for smi, rids in groups.items()]):
            f.result()
    ids = [r for r, _ in rows]
    st = {r: done[r]["status"] for r in ids if r in done}
    (chunk / "build.json").write_text(json.dumps(st, indent=1, sort_keys=True))
    (chunk / "buildtime.json").write_text(json.dumps(
        {r: done[r].get("time_s") for r in ids if r in done}, indent=0, sort_keys=True))
    (chunk / "buildmem.json").write_text(json.dumps(
        {r: done[r]["peak_gb"] for r in ids if r in done and "peak_gb" in done[r]}, indent=1,
        sort_keys=True))
    return st


def cbatch_tool_one(man, rid, specs, archive: Path, workroot: Path, logroot: Path, timeout) -> dict:
    s = man["settings"]
    tool = man["tool"]
    work = workroot / rid
    if tool == "mace":
        cmd = [s["tool_python"], str(WORKERS_DIR / "mace_worker.py"), rid, str(specs), str(archive),
               str(work), "--geoms", s["mode"]]
    else:
        cmd = [s["tool_python"], str(WORKERS_DIR / "bench_worker.py"), tool, rid, str(specs),
               str(archive), str(work), "--mode", s["mode"]]
    env = cbatch_tool_child_env({})
    t0 = time.time()
    killed = {"v": False}
    with open(logroot / f"{rid}.log", "w") as log:
        p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                             env=env)

        def _cbatch_kill_on_limit():
            killed["v"] = True
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

        tm = threading.Timer(timeout, _cbatch_kill_on_limit)
        tm.start()
        _, status, ru = os.wait4(p.pid, 0)
        tm.cancel()
        p.returncode = os.waitstatus_to_exitcode(status)
        if killed["v"]:
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
    wall = round(time.time() - t0, 2)
    mp = archive / "_meta" / f"{rid}.json"
    if killed["v"] or not mp.exists():
        meta = {"refcode": rid, "tool": tool, "mode": s["mode"],
                "status": "timeout" if killed["v"] else f"crash(rc={p.returncode})", "n_frames": 0}
        xyz = archive / f"{rid}.xyz"
        if xyz.exists():
            xyz.unlink()
    else:
        meta = json.loads(mp.read_text())
    meta["wall_total_s"] = wall
    meta["timeout_s"] = timeout
    meta["rss_gb"] = round(ru.ru_maxrss / 1048576.0, 3)
    meta["final"] = True
    tmp = mp.with_name(f".{mp.name}.part")
    tmp.write_text(json.dumps(meta, indent=1))
    os.replace(tmp, mp)
    return meta


def cbatch_read_final_metas(archive: Path, ids) -> dict:
    out = {}
    for rid in ids:
        mp = archive / "_meta" / f"{rid}.json"
        if mp.exists():
            try:
                m = json.loads(mp.read_text())
            except ValueError:
                continue
            if m.get("final"):
                out[rid] = m
    return out


def _tool_shard(man, rows, specs: Path, chunk: Path, workers, timeout, log) -> dict:
    archive = chunk / "archive"
    (archive / "_meta").mkdir(parents=True, exist_ok=True)
    workroot = Path(os.environ.get("DELFIN_CLUSTER_WORK_ROOT") or (chunk / "work"))
    workroot.mkdir(parents=True, exist_ok=True)
    logroot = chunk / "logs"
    logroot.mkdir(exist_ok=True)
    ids = [r for r, _ in rows]
    finished = cbatch_read_final_metas(archive, ids)
    cn = {}
    for ln in specs.read_text().splitlines():
        if ln.strip():
            rec = json.loads(ln)
            cn[rec["refcode"]] = int(rec.get("cn") or 0)
    todo = [r for r in ids if r not in finished]
    todo.sort(key=lambda r: -cn.get(r, 0))
    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex, \
            open(chunk / "run_summary.jsonl", "a") as summ:
        futs = [ex.submit(cbatch_tool_one, man, r, specs, archive, workroot, logroot, timeout)
                for r in todo]
        for f in as_completed(futs):
            m = f.result()
            summ.write(json.dumps(m) + "\n")
            summ.flush()
            finished[m["refcode"]] = m
            log(f"{m['refcode']:<12} {str(m['status'])[:40]:<40} frames={m['n_frames']:<4} "
                f"wall={m['wall_total_s']:8.1f}s rss={m['rss_gb']:.2f}G")
    return {r: finished[r].get("status") for r in ids if r in finished}


def cbatch_run_shard(run_dir, shard, *, set_name="main", run_name="main", workers=None,
                     speed_factor=None, allow_env_change=False, log=print) -> int:
    run_dir = Path(run_dir).resolve()
    man = cbatch_load_manifest(run_dir)
    sets = man["sets"]
    if set_name not in sets:
        log(f"FATAL: shard set {set_name!r} not in this run ({sorted(sets)})")
        return 2
    table = sets[set_name]["shards"]
    shard = int(shard)
    if not 0 <= shard < len(table):
        log(f"FATAL: shard {shard} out of range 0..{len(table) - 1}")
        return 2
    ent = table[shard]
    sdir = run_dir / "shards" / set_name
    if cbatch_sha256_file(sdir / ent["file"]) != ent["sha256"]:
        log(f"FATAL: {ent['file']} sha256 differs from the manifest")
        return 2
    specs = sdir / f"specs_{shard:04d}.jsonl"
    if "specs_sha256" in ent and cbatch_sha256_file(specs) != ent["specs_sha256"]:
        log(f"FATAL: {specs.name} sha256 differs from the manifest")
        return 2
    s = man["settings"]
    chunk = cbatch_chunk_dir(run_dir, set_name, run_name, shard)
    chunk.mkdir(parents=True, exist_ok=True)
    prov = cbatch_provenance(man["tool"], s["tool_python"], s["mode"])
    bad = cbatch_provenance_mismatch(man["provenance"], prov)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    if bad and not allow_env_change:
        (chunk / f"refused_{stamp}.json").write_text(json.dumps(
            {"differs": bad, "manifest": {k: man["provenance"].get(k) for k in bad},
             "here": {k: prov.get(k) for k in bad}}, indent=1))
        log(f"FATAL: code or environment differs from the manifest in {bad}; "
            f"see {chunk}/refused_{stamp}.json (override: --allow-env-change)")
        return 2
    sf = str(speed_factor) if speed_factor is not None else s["speed_factor"]
    timeout = cbatch_effective_timeout(s["timeout_base_s"], sf)
    tj = chunk / "timeout.json"
    tmo = {"timeout_base_s": s["timeout_base_s"], "speed_factor": sf, "timeout_s": timeout}
    if tj.exists() and json.loads(tj.read_text()) != tmo and not allow_env_change:
        log(f"FATAL: this chunk was started with {tj.read_text().strip()}, now {tmo}; "
            "one chunk, one limit (override: --allow-env-change)")
        return 2
    tj.write_text(json.dumps(tmo))
    workers = int(workers or s["workers"])
    rows = cbatch_shard_rows(run_dir, set_name, shard)
    (chunk / f"run_meta_{stamp}.json").write_text(json.dumps(
        {"tool": man["tool"], "set": set_name, "run": run_name, "shard": shard, "workers": workers,
         "threads": s["threads"], "mode": s["mode"], **tmo, "host": os.uname().nodename,
         "slurm_job": os.environ.get("SLURM_JOB_ID"),
         "slurm_array_task": os.environ.get("SLURM_ARRAY_TASK_ID"),
         "python": sys.executable, "provenance": prov, "provenance_differs": bad,
         "n_shard": len(rows), "started": stamp}, indent=1))
    log(f"[{time.strftime('%Y-%m-%dT%H:%M:%S')}] {man['tool']} {set_name}/{run_name} shard {shard}: "
        f"{len(rows)} systems, {workers} workers, limit {timeout} s")
    if man["tool"] == "manta":
        st = _manta_shard(man, rows, chunk, workers, timeout, log)
    else:
        st = _tool_shard(man, rows, specs, chunk, workers, timeout, log)
    missing = [r for r, _ in rows if r not in st]
    if missing:
        log(f"INCOMPLETE: {len(missing)} of {len(rows)} systems without a final record -- "
            "resubmit this index")
        return 4
    tally: dict = {}
    for v in st.values():
        tally[v] = tally.get(v, 0) + 1
    (chunk / "DONE.json").write_text(json.dumps({"n": len(rows), "by_status": tally,
                                                 "finished": time.strftime("%Y-%m-%dT%H:%M:%S")}))
    log(f"DONE shard {shard}: {tally}")
    return 0
