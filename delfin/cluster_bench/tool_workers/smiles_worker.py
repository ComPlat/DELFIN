#!/usr/bin/env python
"""Build ONE system with rxembed, Molassembler or MetalloGen and write it in DELFIN's archive format.

    smiles_worker.py {rxembed|molassembler|metallogen} REFCODE SPECS.jsonl ARCHIVE_DIR WORK_DIR
                     [--mode M]

Same archive format as bench_worker.py: ARCHIVE_DIR/<REF>.xyz (multi-frame; comment
'<REF> frame<i> <label> tool=<tool> E=<energy or na>') and ARCHIVE_DIR/_meta/<REF>.json.  The build
is the one the dashboard and CONTROL run (smiles_tools.py, imported as a sibling file; nothing
else from DELFIN), with the same defaults.  rxembed and Molassembler take the SMILES of the spec
line as written; MetalloGen takes the spec (DELFIN's split) and needs ORCA on PATH.  Atom order is
whatever the tool emits; compare structures by graph isomorphism.  Runs in the tool's own
environment, single-threaded.

Modes: rxembed ``default`` (rx.embed(rx.isomers(smiles))); Molassembler ``default`` (every
stereoisomer, 10 conformers each); MetalloGen ``common`` (the usual polyhedra of the site count)
or ``all`` (every geometry MetalloGen knows for it).
"""
import argparse
import contextlib
import json
import os
import shutil
import sys
import time
import traceback

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import smiles_tools as st  # noqa: E402

MODES = {"rxembed": ("default",), "molassembler": ("default",), "metallogen": ("common", "all")}


def sw_write_archive(path, ref, frames, tool):
    with open(path, "w") as fh:
        for i, f in enumerate(frames):
            e = f.get("energy")
            fh.write("%d\n" % len(f["symbols"]))
            fh.write("%s frame%d %s tool=%s E=%s\n"
                     % (ref, i, f["label"], tool, "na" if e is None else "%.6f" % e))
            for s, (x, y, z) in zip(f["symbols"], f["coords"]):
                fh.write("%-2s %14.6f %14.6f %14.6f\n" % (s, x, y, z))


def sw_run(tool, rec, work, mode):
    """``(frames, status, extra)`` of one spec line."""
    if tool in st.TAKES_SPEC:
        if rec.get("status") != "ok":
            return [], "convert:not_expressible:" + str(rec.get("bail_msg") or rec.get("bail")), {}
        # Each finished arrangement is appended there; a run killed at its limit keeps it.
        source, options = rec, {"geometries": mode,
                                "partial_xyz": os.path.join(work, "partial.xyz")}
    else:
        source, options = rec["smiles"], {}
    try:
        frames, extra = st.RUNNERS[tool](source, work, options)
    except st.NotExpressible as exc:
        return [], "convert:not_expressible:" + str(exc)[:200], {}
    return frames, ("ok" if frames else "no_structure"), extra


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("tool", choices=sorted(MODES))
    ap.add_argument("ref")
    ap.add_argument("specs")
    ap.add_argument("archive")
    ap.add_argument("work")
    ap.add_argument("--mode", default=None)
    a = ap.parse_args()
    mode = a.mode or MODES[a.tool][0]
    if mode not in MODES[a.tool]:
        ap.error("mode %r not valid for %s: %s" % (mode, a.tool, MODES[a.tool]))
    os.makedirs(os.path.join(a.archive, "_meta"), exist_ok=True)
    if os.path.isdir(a.work):
        shutil.rmtree(a.work)
    os.makedirs(a.work)
    rec = None
    for line in open(a.specs):
        if line.startswith('{"refcode": "%s"' % a.ref):
            rec = json.loads(line)
            break
    t0 = time.time()
    meta = {"refcode": a.ref, "tool": a.tool, "mode": mode}
    log = open(os.path.join(a.work, "tool.log"), "w")
    try:
        with contextlib.redirect_stdout(log):
            if rec is None:
                frames, status, extra = [], "no_spec", {}
            else:
                frames, status, extra = sw_run(a.tool, rec, a.work, mode)
    except (Exception, SystemExit) as e:  # MetalloGen also exit()s
        frames, status, extra = [], "exception:%s:%s" % (type(e).__name__, str(e)[:200]), {}
        meta["traceback"] = traceback.format_exc()[-2000:]
    log.close()
    if status.startswith("convert:not_expressible:"):
        meta["not_expressible_reason"] = status.split(":", 2)[2]
    meta.update({"status": status, "n_frames": len(frames), "wall_s": round(time.time() - t0, 2)})
    meta.update(extra)
    if frames:
        sw_write_archive(os.path.join(a.archive, a.ref + ".xyz"), a.ref, frames, a.tool)
    json.dump(meta, open(os.path.join(a.archive, "_meta", a.ref + ".json"), "w"), indent=1)
    if status == "ok" and os.environ.get("BENCH_KEEP_WORK", "0") != "1":
        shutil.rmtree(a.work, ignore_errors=True)
