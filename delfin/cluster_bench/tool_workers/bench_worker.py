#!/usr/bin/env python
"""Build ONE system with ONE external builder and write it in DELFIN's archive format.

    bench_worker.py {architector|molsimplify} REFCODE SPECS.jsonl ARCHIVE_DIR WORK_DIR [--mode M]

Writes ARCHIVE_DIR/<REFCODE>.xyz (multi-frame: count line, comment line, 'Sym x y z' lines;
comment = '<REF> frame<i> <label> tool=<tool> E=<energy or na>') and
ARCHIVE_DIR/_meta/<REFCODE>.json (status, n_frames, wall time, error).  Atom order is whatever
the tool emits; compare structures by graph isomorphism, not by atom index.  Runs in the Architector/molSimplify environment, single-threaded.

Architector modes:
  full    : n_conformers = n_symmetries = 10  -> every distinct symmetry (isomer) Architector
            builds per core geometry is relaxed and returned (its own RMSD/energy dedup applies)
  default : Architector defaults (n_symmetries 10, n_conformers 1 -> one isomer per core geometry)
molSimplify: one build per molSimplify geometry of the CN (bench_convert.MS_GEOMS).
"""
import argparse
import contextlib
import glob
import json
import os
import shutil
import sys
import time
import traceback

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
os.environ.setdefault("OMP_STACKSIZE", "1G")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bench_convert as bc  # noqa: E402


def bench_write_archive(path, ref, frames, tool):
    with open(path, "w") as fh:
        for i, (label, syms, xyz, e) in enumerate(frames):
            fh.write("%d\n" % len(syms))
            fh.write("%s frame%d %s tool=%s E=%s\n"
                     % (ref, i, label, tool, "na" if e is None else "%.6f" % e))
            for s, (x, y, z) in zip(syms, xyz):
                fh.write("%-2s %14.6f %14.6f %14.6f\n" % (s, x, y, z))


def bench_run_architector(rec, work, mode):
    from architector import build_complex
    params = {"temp_prefix": work.rstrip("/") + "/"}
    if mode == "full":
        params.update({"n_symmetries": 10, "n_conformers": 10})
    inp, err = bc.bench_to_architector(rec, params)
    if err:
        return [], "convert:" + err, {}
    out = build_complex(inp)
    frames = []
    for key, v in out.items():
        at = v["ase_atoms"]
        frames.append((key, at.get_chemical_symbols(), at.get_positions().tolist(),
                       float(v.get("energy")) if v.get("energy") is not None else None))
    return frames, ("ok" if frames else "no_structure"), {"keys": list(out.keys())}


def bench_run_molsimplify(rec, work):
    from molSimplify.Scripts.generator import startgen_pythonic
    runs, err = bc.bench_to_molsimplify(rec)
    if err:
        return [], "convert:" + err, {}
    frames, per_geo = [], {}
    for geo, d in runs:
        rd = os.path.join(work, geo)
        name = "%s_%s" % (rec["refcode"], geo)
        d = dict(d)
        d["-rundir"] = rd + "/"
        d["-name"] = name
        try:
            startgen_pythonic(d, write=True)
        except (Exception, SystemExit) as e:  # one geometry failing (molSimplify also quit()s) must not hide the others
            per_geo[geo] = "exception:%s:%s" % (type(e).__name__, str(e)[:120])
            continue
        xyzs = sorted(glob.glob(os.path.join(rd, "**", name + ".xyz"), recursive=True))
        if not xyzs:
            per_geo[geo] = "no_xyz"
            continue
        per_geo[geo] = "ok"
        for f in xyzs:
            ln = open(f).read().splitlines()
            n = int(ln[0].split()[0])
            syms, xyz = [], []
            for l in ln[2:2 + n]:
                p = l.split()
                syms.append(p[0])
                xyz.append([float(p[1]), float(p[2]), float(p[3])])
            frames.append((geo, syms, xyz, None))
    st = "ok" if frames else ("no_structure:" + ";".join("%s=%s" % kv for kv in per_geo.items()))
    return frames, st, {"per_geometry": per_geo}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("tool", choices=["architector", "molsimplify"])
    ap.add_argument("ref")
    ap.add_argument("specs")
    ap.add_argument("archive")
    ap.add_argument("work")
    ap.add_argument("--mode", default="full")
    a = ap.parse_args()
    os.makedirs(os.path.join(a.archive, "_meta"), exist_ok=True)
    if os.path.isdir(a.work):
        shutil.rmtree(a.work)
    os.makedirs(a.work)
    rec = None
    for l in open(a.specs):
        if l.startswith('{"refcode": "%s"' % a.ref):
            rec = json.loads(l)
            break
    t0 = time.time()
    meta = {"refcode": a.ref, "tool": a.tool, "mode": a.mode}
    log = open(os.path.join(a.work, "tool.log"), "w")
    try:
        with contextlib.redirect_stdout(log):
            if rec is None:
                frames, st, extra = [], "no_spec", {}
            elif a.tool == "architector":
                frames, st, extra = bench_run_architector(rec, a.work, a.mode)
            else:
                frames, st, extra = bench_run_molsimplify(rec, a.work)
    except Exception as e:
        frames, st, extra = [], "exception:%s:%s" % (type(e).__name__, str(e)[:200]), {}
        meta["traceback"] = traceback.format_exc()[-2000:]
    log.close()
    meta.update({"status": st, "n_frames": len(frames), "wall_s": round(time.time() - t0, 2)})
    meta.update(extra)
    if frames:
        bench_write_archive(os.path.join(a.archive, a.ref + ".xyz"), a.ref, frames, a.tool)
    json.dump(meta, open(os.path.join(a.archive, "_meta", a.ref + ".json"), "w"), indent=1)
    if st == "ok" and os.environ.get("BENCH_KEEP_WORK", "0") != "1":
        shutil.rmtree(a.work, ignore_errors=True)
