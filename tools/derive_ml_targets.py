#!/usr/bin/env python
"""Derive the bond-valence metal-donor target tables of MANTA (DELFIN_FFFREE_ML_TARGETS_BV).

The shipped table (delfin/manta/data/ml_targets_cod.tsv) comes from the Crystallography Open
Database (CC0) and is rebuilt with the ``records-cod`` and ``table`` steps below.  A table from the
Cambridge Structural Database is NOT shipped: anyone with a CSD licence derives it locally with the
``records-pairs`` and ``table`` steps from crystal coordinates exported from their own CSD
installation (tools/export_csd_crystals.py, run under the CSD Python API), and points MANTA at it
with DELFIN_FFFREE_ML_TARGETS_BV_TABLE.  Its rows are statistics of licensed data -- keep them out of
any public repository.

Held-out systems never enter a derivation: ``--exclude`` takes id;SMILES (or id|SMILES, or id) files
and drops every structure of the same refcode family (ABCDEF, ABCDEF01, ...) and every structure with
the same molecular formula.

    derive_ml_targets.py records-pairs PAIRS --xyz-dir DIR [--xyz-dir DIR ...] --out REC.jsonl
                         [--exclude F1,F2] [--workers N] [--max N] [--all-members]
    derive_ml_targets.py records-cod --xyz-dir COD_XYZ --cif-dir COD_CIF --out REC.jsonl
                         [--exclude F1,F2] [--workers N] [--max N]
    derive_ml_targets.py perception-check PAIRS --xyz-dir DIR [...] --out CHECK.jsonl [--max N]
    derive_ml_targets.py table REC.jsonl --out TABLE.tsv --source TEXT [--consensus 0.5]
                         [--fit] [--b B --rho RHO] [--min-n 10] [--boot 200]

PAIRS: one ``id;SMILES`` or ``id|SMILES`` per line (the SMILES as MANTA reads it).  Crystal files are
``<DIR>/<id>.xyz`` with the single complex (largest metal-bearing component, primary disorder site).
records-pairs keeps one crystal per refcode family (the first id in sorted order that has a file)
unless --all-members.  See delfin/manta/ml_targets_derive.py for every rule.
"""
from __future__ import annotations

import argparse
import json
import signal
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")
from delfin.manta import ml_targets_derive as MD  # noqa: E402

TIMEOUT_S = 60


class _StructureTimeout(Exception):
    pass


def _structure_alarm(*_):
    raise _StructureTimeout()


def _run_one_structure(fn, *args):
    signal.signal(signal.SIGALRM, _structure_alarm)
    signal.alarm(TIMEOUT_S)
    try:
        return fn(*args)
    except _StructureTimeout:
        return {"id": args[0], "st": "timeout"}
    except Exception as e:                                   # one structure never stops a run
        return {"id": args[0], "st": "err", "err": str(e)[:160]}
    finally:
        signal.alarm(0)


def _find_xyz(dirs, rid):
    for d in dirs:
        p = Path(d) / f"{rid}.xyz"
        if p.exists():
            return p
    return None


def _pairs_record(rid, smi, path):
    fr = MD.read_first_xyz_frame(path)
    if fr is None:
        return {"id": rid, "st": "crystal_parse"}
    return MD.record_from_smiles_crystal(rid, smi, *fr)


def _check_record(rid, smi, path):
    from delfin.manta import _ml_targets_bv as BV
    from delfin.manta._oxidation_state import _OX_GROUP
    fr = MD.read_first_xyz_frame(path)
    if fr is None:
        return {"id": rid, "st": "crystal_parse"}
    el, X = fr
    mol = BV._mol_from_smiles(smi)
    if mol is None:
        return {"id": rid, "st": "smiles"}
    mis = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in _OX_GROUP]
    if len(mis) != 1:
        return {"id": rid, "st": "not_one_tm"}
    cd = BV.complex_donors(mol, mis[0])
    if cd is None:
        return {"id": rid, "st": "cut_none"}
    xm = [i for i, e in enumerate(el) if MD.EL.is_metal(e)]
    if len(xm) != 1:
        return {"id": rid, "st": "crystal_metal"}
    qsum, pdon, _nh, st = MD.perceive_complex(el, X, xm[0])
    if pdon is None:
        return {"id": rid, "st": "perc_" + st}
    qtot = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    return {"id": rid, "st": "ok", "metal": cd[0], "ox": cd[1], "ox_p": qtot - qsum,
            "bvs": MD.bvs_literature(cd[0], [(e, d) for _g, e, _c, d in pdon])}


def _cod_record(cid, path, cif_dir, excl_forms):
    fr = MD.read_first_xyz_frame(path)
    if fr is None:
        return {"id": cid, "st": "parse"}
    fs = MD.cif_formula_sum(Path(cif_dir) / f"{cid}.cif")
    return MD.record_from_cod(cid, fr[0], fr[1], fs, excl_forms)


def _job_pairs(args):
    return _run_one_structure(_pairs_record, *args)


def _job_check(args):
    return _run_one_structure(_check_record, *args)


def _job_cod(args):
    path, cif_dir, excl_forms = args
    return _run_one_structure(_cod_record, Path(path).stem, path, cif_dir, excl_forms)


def _read_pairs(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            sep = ";" if ";" in ln else ("|" if "|" in ln else None)
            if not sep:
                continue
            rid, smi = ln.split(sep, 1)
            out.append((rid.strip(), smi.strip()))
    return out


def _pair_jobs(a, keep_one_per_family=True):
    ex_roots, ex_forms = MD.load_exclusions(a.exclude.split(",") if a.exclude else [])
    fam, n_ex = {}, Counter()
    for rid, smi in _read_pairs(a.pairs):
        if MD.root_of(rid) in ex_roots:
            n_ex["family"] += 1
            continue
        if MD.smiles_formula(smi) in ex_forms:
            n_ex["formula"] += 1
            continue
        key = MD.root_of(rid) if keep_one_per_family else rid
        fam.setdefault(key, []).append((rid, smi))
    jobs = []
    for k in sorted(fam):
        for rid, smi in sorted(fam[k]):
            p = _find_xyz(a.xyz_dir, rid)
            if p is not None:
                jobs.append((rid, smi, str(p)))
                if keep_one_per_family:
                    break
    print(f"# excluded {dict(n_ex)}; {len(fam)} families/ids; {len(jobs)} with a crystal", flush=True)
    return jobs[:a.max] if a.max else jobs


def _run_jobs(jobs, fn, out, workers):
    st = Counter()
    with Pool(workers) as pl, open(out, "w") as fh:
        for r in pl.imap_unordered(fn, jobs, chunksize=8):
            st[r.get("st")] += 1
            fh.write(json.dumps(r) + "\n")
    print("# status", dict(st.most_common()), flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("records-pairs", "perception-check"):
        p = sub.add_parser(name)
        p.add_argument("pairs")
        p.add_argument("--xyz-dir", action="append", required=True)
        p.add_argument("--out", required=True)
        p.add_argument("--exclude", default="")
        p.add_argument("--workers", type=int, default=8)
        p.add_argument("--max", type=int, default=0)
        p.add_argument("--all-members", action="store_true")
    p = sub.add_parser("records-cod")
    p.add_argument("--xyz-dir", required=True)
    p.add_argument("--cif-dir", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--exclude", default="")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--max", type=int, default=0)
    p = sub.add_parser("table")
    p.add_argument("records")
    p.add_argument("--out", required=True)
    p.add_argument("--source", required=True)
    p.add_argument("--consensus", type=float, default=None,
                   help="keep a record only if |BVS(literature R0) - ox| <= this (COD rows)")
    p.add_argument("--fit", action="store_true", help="choose b and rho on half the records")
    p.add_argument("--b", type=float, default=0.37)
    p.add_argument("--rho", type=float, default=1.0)
    p.add_argument("--min-n", type=int, default=10)
    p.add_argument("--boot", type=int, default=200)
    p.add_argument("--exclude", default="", help="drop records of these families/formulas (id;SMILES files)")
    p.add_argument("--jt", action="store_true", help="add Cu(II) Jahn-Teller axial/equatorial rows (J)")
    p.add_argument("--spin", action="store_true",
                   help="label spin states (3d d4-d7, CN >= 5) and add per-spin rows and low-spin shares")
    a = ap.parse_args(argv)

    if a.cmd == "records-pairs":
        _run_jobs(_pair_jobs(a, not a.all_members), _job_pairs, a.out, a.workers)
    elif a.cmd == "perception-check":
        _run_jobs(_pair_jobs(a, not a.all_members), _job_check, a.out, a.workers)
    elif a.cmd == "records-cod":
        _ex_roots, ex_forms = MD.load_exclusions(a.exclude.split(",") if a.exclude else [])
        files = sorted(Path(a.xyz_dir).glob("*.xyz"))
        files = files[:a.max] if a.max else files
        _run_jobs([(str(f), a.cif_dir, ex_forms) for f in files], _job_cod, a.out, a.workers)
    else:
        recs = [json.loads(ln) for ln in open(a.records)]
        recs = [r for r in recs if r.get("st") == "ok"]
        if a.exclude:
            ex_roots, ex_forms = MD.load_exclusions(a.exclude.split(","))
            recs = [r for r in recs if MD.root_of(r["id"]) not in ex_roots and r.get("formula") not in ex_forms]
        n0 = len(recs)
        if a.consensus is not None:
            recs = [r for r in recs if MD.consensus_ok(r, a.consensus)]
        b, rho = a.b, a.rho
        head = [f"Bond-valence M-D target table -- {a.source}",
                f"{len(recs)} structures ({n0} before the oxidation-state consensus)"
                + (f", consensus |BVS - ox| <= {a.consensus}" if a.consensus is not None else ""),
                f"built by tools/derive_ml_targets.py table; rows with >= {a.min_n} structures"]
        if a.fit:
            f = MD.fit_bv(recs)
            b, rho = f["b"], f["rho"]
            head.append(f"b, rho fitted on half A, scored on half B of these structures: b = {b}, rho = {rho}, "
                        f"median |R-level error| {f['score']:.4f} A (n = {f['n_scored']}); literature "
                        f"b = 0.37, rho = 1: {f['literature_score']:.4f} A (n = {f['literature_n']})")
            print("# fit", {k: v for k, v in f.items() if k != "grid"}, flush=True)
        rows = MD.build_table(recs, b, rho, a.min_n, a.boot)
        if a.jt:
            rows += MD.jt_rows(recs, "Cu", 2, a.min_n, a.boot)
            rows += MD.jt_rows_sp5(recs, "Cu", 2, a.min_n, a.boot)
            head.append("J rows: Cu(II) CN6 axial/equatorial M-D per element and Jahn-Teller form (elongated / compressed)")
        if a.spin:
            summ = MD.label_spins(recs)
            nl = Counter(r.get("spin") for r in recs if r.get("spin"))
            head.append(f"spin labels (two-component split, posterior >= {MD.SPIN_POSTERIOR}): {dict(nl)}")
            for k in sorted(summ):
                head.append(f"  spin split {k[0]}({k[1]}): {summ[k]}")
            rows += MD.spin_rows(recs, a.min_n, a.boot)
        MD.write_target_table(a.out, rows, b, rho, head)
        print(f"# {len(recs)} structures -> {len(rows)} rows -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
