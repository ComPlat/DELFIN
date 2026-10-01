#!/usr/bin/env python
"""Build ONE system with epic-MACE and write it in DELFIN's archive format.

    mace_worker.py REFCODE SPECS.jsonl ARCHIVE_DIR WORK_DIR [--geoms paper|extended]
                   [--num-confs 10] [--max-attempts 10]

Same archive format as bench_worker.py: ARCHIVE_DIR/<REF>.xyz (multi-frame; comment
'<REF> frame<i> <label> tool=mace E=<UFF energy of MACE's force field>') and
ARCHIVE_DIR/_meta/<REF>.json.  Label = '<GEOM>-iso<k>-conf<j>': k = MACE stereomer index (the
order of Complex.GetStereomers), j = conformer rank by MACE's MM energy.  MACE's centroid dummies
of hapto ligands (element X in its xyz) are dropped; atom order otherwise as MACE emits it
(the eye maps by graph isomorphism).  Runs in the epic-MACE environment, single-threaded.

Settings = the defaults of epic-MACE's CLI (mace/__main__.py: prepare_complexes +
run_mace_for_system), with one deliberate exception (enantiomers are kept):
  * input: mace.ComplexFromLigands(ligand SMILES with mapped donors, central atom, geom) from
    mace_convert.py; all donor map numbers are reset to 1 and the complex is rebuilt
    (ComplexFromMol), exactly as the CLI does before a stereomer search;
  * stereomers: GetStereomers(regime='all', dropEnantiomers=False, minTransCycle=None,
    merRule=False): every arrangement at the metal AND every unassigned ligand stereocentre
    (the held-out SMILES carry no stereo), enantiomers kept (the crystal is one of them;
    MANTA and the eye distinguish them; CLI default would drop one of each pair).
    minTransCycle=None = CLI default --trans-cycle unset (no chelate ring spans trans
    positions); merRule=False = CLI default (--mer-rule not given).  The library API default
    merRule=True was tried first and rejected: its empirical "rigid X-Y-Z only mer" rule
    returns ZERO stereomers for fac-only tripods (scorpionate/tripodal tetradentates in the
    pilot) and dropped 4 of 10 stereomers of a hexacoordinate pilot system;
  * 3D: AddConformers(numConfs=10, maxAttempts=10, rmsThresh=-1) per stereomer (the CLI default
    num-confs 10; its default rms-thresh 0.0 never removes a conformer, -1 is the same result
    without the RMS computations), then OrderConfsByEnergy.  Each conformer: RDKit distance
    geometry with MACE's bounds/coordMap (enforceChirality, random coordinates), UFF with MACE's
    metal-centre constraints (maxIts 1000), rejected if the metal centre's chirality is wrong;
  * no representative-conformer selection (the CLI default).  All conformers are written.
  * random seed: MACE does not expose one (RDKit EmbedParameters.randomSeed stays -1).
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
import mace_convert as mc  # noqa: E402

REGIME, DROP_ENANTIOMERS, MIN_TRANS_CYCLE, MER_RULE = "all", False, None, False


def mw_write_archive(path, ref, frames):
    with open(path, "w") as fh:
        for i, (label, syms, xyz, e) in enumerate(frames):
            fh.write("%d\n" % len(syms))
            fh.write("%s frame%d %s tool=mace E=%s\n"
                     % (ref, i, label, "na" if e is None else "%.6f" % e))
            for s, (x, y, z) in zip(syms, xyz):
                fh.write("%-2s %14.6f %14.6f %14.6f\n" % (s, x, y, z))


def mw_frames_of(X, geom, k):
    out = []
    for j in range(X.GetNumConformers()):
        block = X.ToXYZBlock(j).splitlines()
        n = int(block[0])
        e = json.loads(block[1]).get("E")
        syms, xyz = [], []
        for l in block[2:2 + n]:
            p = l.split()
            if p[0] == "X":  # MACE's centroid dummy of a hapto ligand: not an atom
                continue
            syms.append(p[0])
            xyz.append([float(p[1]), float(p[2]), float(p[3])])
        out.append(("%s-iso%d-conf%d" % (geom, k, j), syms, xyz, e))
    return out


def mw_run(rec, geoms, num_confs, max_attempts, meta):
    import mace
    job, err = mc.mace_convert(rec, geoms)
    if job is not None:
        meta["convert_info"] = job["info"]
    if err:
        return [], "convert:" + err
    meta["mace_input"] = {"ligands": job["ligands"], "CA": job["CA"]}
    frames, per = [], {}
    for geom in job["geoms"]:
        g = {}
        per[geom] = g
        t0 = time.time()
        try:
            X = mace.ComplexFromLigands(job["ligands"], job["CA"], geom)
            g["eta_detected"] = sorted(len(v) for v in getattr(X, "_haptic_DAs", {}).values())
            for idx in X._DAs:  # as the CLI (prepare_complexes) before a stereomer search
                X.mol.GetAtomWithIdx(idx).SetAtomMapNum(1)
                X.mol.GetAtomWithIdx(idx).SetIsotope(1)
            X = mace.ComplexFromMol(X.mol, X.geom)
            Xs = X.GetStereomers(REGIME, DROP_ENANTIOMERS, MIN_TRANS_CYCLE, MER_RULE)
        except Exception as e:
            g["status"] = "stereo_exception:%s:%s" % (type(e).__name__,
                                                      str(e).strip().splitlines()[0][:120]
                                                      if str(e).strip() else "")
            g["t_stereo_s"] = round(time.time() - t0, 2)
            continue
        g["t_stereo_s"] = round(time.time() - t0, 2)
        g["n_stereomers"] = len(Xs)
        # classes up to mirror image (MACE's own IsEnantiomer), for comparison with counts
        # reported without enantiomers
        rep = []
        for x in Xs:
            try:
                if not any(x.IsEqual(r) or x.IsEnantiomer(r) for r in rep):
                    rep.append(x)
            except Exception:
                rep.append(x)
        g["n_stereomers_no_enantiomers"] = len(rep)
        t1 = time.time()
        noconf, errs = [], {}
        for k, x in enumerate(Xs):
            try:
                x.AddConformers(numConfs=num_confs, maxAttempts=max_attempts, rmsThresh=-1)
                if x.GetNumConformers():
                    x.OrderConfsByEnergy()
                    frames += mw_frames_of(x, geom, k)
                else:
                    noconf.append(k)
            except Exception as e:
                noconf.append(k)
                errs[k] = "%s:%s" % (type(e).__name__, str(e).strip()[:120])
        g["t_embed_s"] = round(time.time() - t1, 2)
        g["stereomers_without_conformer"] = noconf
        if errs:
            g["embed_exceptions"] = errs
        g["n_frames"] = sum(1 for f in frames if f[0].startswith(geom + "-"))
        g["status"] = "ok" if g["n_frames"] else "no_conformer"
    meta["per_geometry"] = per
    if frames:
        return frames, "ok"
    return [], "no_structure:" + ";".join("%s=%s" % (k, v.get("status")) for k, v in per.items())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("specs")
    ap.add_argument("archive")
    ap.add_argument("work")
    ap.add_argument("--geoms", default="paper", choices=sorted(mc.GEOMS))
    ap.add_argument("--num-confs", type=int, default=10)
    ap.add_argument("--max-attempts", type=int, default=10)
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
    meta = {"refcode": a.ref, "tool": "mace", "mode": a.geoms,
            "settings": {"regime": REGIME, "dropEnantiomers": DROP_ENANTIOMERS,
                         "minTransCycle": MIN_TRANS_CYCLE, "merRule": MER_RULE,
                         "numConfs": a.num_confs, "maxAttempts": a.max_attempts,
                         "rmsThresh": -1, "maxResonanceStructures": 1}}
    log = open(os.path.join(a.work, "tool.log"), "w")
    try:
        with contextlib.redirect_stdout(log):
            if rec is None:
                frames, st = [], "no_spec"
            else:
                frames, st = mw_run(rec, a.geoms, a.num_confs, a.max_attempts, meta)
    except Exception as e:
        frames, st = [], "exception:%s:%s" % (type(e).__name__, str(e)[:200])
        meta["traceback"] = traceback.format_exc()[-2000:]
    log.close()
    if st.startswith("convert:not_expressible:"):
        meta["not_expressible_reason"] = st.split(":", 2)[2]
    meta.update({"status": st, "n_frames": len(frames), "wall_s": round(time.time() - t0, 2)})
    per = meta.get("per_geometry", {})
    meta["n_stereoisomers"] = sum(v.get("n_stereomers", 0) for v in per.values())
    meta["n_stereoisomers_no_enantiomers"] = sum(v.get("n_stereomers_no_enantiomers", 0)
                                                 for v in per.values())
    if frames:
        mw_write_archive(os.path.join(a.archive, a.ref + ".xyz"), a.ref, frames)
    json.dump(meta, open(os.path.join(a.archive, "_meta", a.ref + ".json"), "w"), indent=1)
    if st == "ok" and os.environ.get("BENCH_KEEP_WORK", "0") != "1":
        shutil.rmtree(a.work, ignore_errors=True)
