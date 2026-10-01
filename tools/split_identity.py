#!/usr/bin/env python
"""Identity harness for the MANTA constructor split.

Builds a fixed set of SMILES (``tools/split_identity_smiles.tsv``) exactly the
way the ``delfin manta`` CLI and the dashboard build them -- the shipped
champion configuration from ``delfin.cli_manta.construction_env`` applied to a
fresh subprocess with ``PYTHONHASHSEED=0`` (``delfin.common.manta_build``) --
and compares the emitted multi-frame XYZ, the frame labels, the error string
and the FF-free ISO trace byte for byte against a recorded reference.

Three ways to use it::

    # 1. record the reference from the working tree (writes the manifest)
    python tools/split_identity.py --record

    # 2. compare the working tree against the recorded manifest
    python tools/split_identity.py

    # 3. compare the working tree against a reference commit built on the fly
    python tools/split_identity.py --ref-commit <sha>

The manifest (``tools/split_identity_manifest.json``) stores the sha256 of
every artefact per SMILES plus the commit it was recorded at.  A compare run
rebuilds every SMILES on the tree under test and reports N/N identical; any
difference is listed with the first differing line.  Exit status 0 only when
everything is identical.

The interpreter that runs this script is the one the builds use (it should be
the environment the loop builds with).  The tree under test is put in front of
``PYTHONPATH`` so an installed DELFIN elsewhere is never picked up.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_SMILES = HERE / "split_identity_smiles.tsv"
DEFAULT_MANIFEST = HERE / "split_identity_manifest.json"

# The build-side script: the same four lines delfin.common.manta_build runs in
# its subprocess, plus the frame serialisation the CLI writes (count line,
# label line, atom lines).
_WORKER = r"""
import json, sys
from delfin.smiles_converter import smiles_to_xyz_isomers
req = json.loads(sys.stdin.read())
res, err = smiles_to_xyz_isomers(req['smiles'], **req['kwargs'])
frames = []
labels = []
for xyz, label in (res or []):
    lines = [ln for ln in xyz.splitlines() if ln.strip()]
    if len(lines) >= 2 and lines[0].strip().isdigit():
        lines = lines[2:]
    frames.append("%d\n%s\n%s\n" % (len(lines), label or "", "\n".join(lines)))
    labels.append(label)
out = {"xyz": "".join(frames), "labels": labels, "error": err}
sys.stdout.write("__SPLIT_IDENTITY__" + json.dumps(out))
"""

# How the CLI calls the builder when no option is given (delfin.cli_manta.main):
# every isomer, no label collapse, binding-mode isomers on, UFF on,
# deterministic, quality extreme.  ``max_isomers`` is read from the tree.
_ENV_PROBE = r"""
import json
from delfin import cli_manta
print(json.dumps({
    "env": cli_manta.construction_env("champion"),
    "max_isomers": cli_manta._ALL_ISOMERS,
}))
"""


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_smiles(path: Path) -> list[tuple[str, str]]:
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("id\t"):
            continue
        parts = line.split("\t")
        rows.append((parts[0], parts[1]))
    return rows


def tree_env(tree: Path) -> dict:
    """Environment for a subprocess that must import DELFIN from ``tree``."""
    env = dict(os.environ)
    env["PYTHONPATH"] = str(tree) + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    return env


def probe_construction(tree: Path) -> dict:
    out = subprocess.run([sys.executable, "-c", _ENV_PROBE], env=tree_env(tree),
                         capture_output=True, text=True, check=True, cwd=str(tree))
    return json.loads(out.stdout.strip().splitlines()[-1])


def build_one(tree: Path, construction: dict, sid: str, smiles: str,
              out_dir: Path, timeout: float | None) -> dict:
    """Build one SMILES on ``tree`` and write xyz / labels / trace files."""
    env = tree_env(tree)
    env.update(construction["env"])
    env["PYTHONHASHSEED"] = "0"
    trace_path = out_dir / f"{sid}.trace"
    if trace_path.exists():
        trace_path.unlink()
    env["DELFIN_FFFREE_ISO_TRACE"] = str(trace_path)
    kwargs = {
        "max_isomers": construction["max_isomers"],
        "collapse_label_variants": False,
        "include_binding_mode_isomers": True,
        "apply_uff": True,
        "deterministic": True,
        "quality_mode": "extreme",
    }
    payload = json.dumps({"smiles": smiles, "kwargs": kwargs})
    t0 = time.time()
    proc = subprocess.Popen([sys.executable, "-c", _WORKER], stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True, env=env, cwd=str(tree), start_new_session=True)
    try:
        out, err = proc.communicate(input=payload, timeout=timeout)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, 9)
        except Exception:
            pass
        return {"id": sid, "status": "timeout", "seconds": time.time() - t0}
    seconds = time.time() - t0
    marker = "__SPLIT_IDENTITY__"
    result = None
    for line in out.splitlines():
        idx = line.find(marker)
        if idx >= 0:
            result = json.loads(line[idx + len(marker):])
    if proc.returncode != 0 or result is None:
        tail = (err or "").splitlines()[-8:]
        return {"id": sid, "status": "failed", "seconds": seconds,
                "detail": f"exit {proc.returncode}: " + " | ".join(tail)}
    xyz = result["xyz"].encode()
    labels = json.dumps(result["labels"], ensure_ascii=False).encode()
    error = json.dumps(result["error"]).encode()
    (out_dir / f"{sid}.xyz").write_bytes(xyz)
    (out_dir / f"{sid}.labels.json").write_bytes(labels)
    (out_dir / f"{sid}.error.json").write_bytes(error)
    if not trace_path.exists():
        trace_path.write_bytes(b"")
    trace = trace_path.read_bytes()
    n_atoms = int(xyz.split(b"\n", 1)[0]) if xyz else 0
    return {
        "id": sid, "status": "ok", "seconds": round(seconds, 1),
        "n_frames": len(result["labels"]), "n_atoms": n_atoms,
        "sha_xyz": _sha(xyz), "sha_labels": _sha(labels),
        "sha_error": _sha(error), "sha_trace": _sha(trace),
    }


def build_set(tree: Path, rows, out_dir: Path, workers: int, timeout) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    construction = probe_construction(tree)
    results = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(build_one, tree, construction, sid, smi, out_dir, timeout): sid
                for sid, smi in rows}
        for fut in concurrent.futures.as_completed(futs):
            r = fut.result()
            results[r["id"]] = r
            print(f"  built {r['id']:18s} {r['status']:8s} "
                  f"{r.get('n_frames', '-')!s:>4} frames  {r['seconds']:6.1f} s",
                  flush=True)
    return {"construction": construction, "results": results}


def first_diff(a: bytes, b: bytes) -> str:
    la, lb = a.splitlines(), b.splitlines()
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            return f"line {i + 1}: {x[:70]!r} != {y[:70]!r}"
    if len(la) != len(lb):
        return f"length {len(la)} vs {len(lb)} lines"
    return "identical"


def compare(ref: dict, new: dict, ref_dir: Path | None, new_dir: Path) -> int:
    keys = ("sha_xyz", "sha_labels", "sha_error", "sha_trace")
    n_ok = 0
    n_total = 0
    for sid in sorted(new):
        n_total += 1
        r, n = ref.get(sid), new[sid]
        if r is None:
            print(f"  {sid:18s} NOT IN REFERENCE")
            continue
        if r.get("status") != "ok" or n.get("status") != "ok":
            same = r.get("status") == n.get("status")
            print(f"  {sid:18s} {'same' if same else 'DIFF'} status "
                  f"{r.get('status')} / {n.get('status')} {n.get('detail', '')}")
            n_ok += int(same)
            continue
        bad = [k for k in keys if r[k] != n[k]]
        if not bad:
            n_ok += 1
            continue
        print(f"  {sid:18s} DIFF in {', '.join(bad)}")
        if ref_dir is not None:
            for k, suffix in (("sha_xyz", ".xyz"), ("sha_labels", ".labels.json"),
                              ("sha_error", ".error.json"), ("sha_trace", ".trace")):
                if k in bad:
                    a = (ref_dir / f"{sid}{suffix}").read_bytes()
                    b = (new_dir / f"{sid}{suffix}").read_bytes()
                    print(f"      {suffix}: {first_diff(a, b)}")
    print(f"\nidentity: {n_ok}/{n_total} identical")
    return 0 if n_ok == n_total and n_total > 0 else 1


def export_commit(commit: str, dest: Path) -> None:
    """``git archive`` of ``commit`` into ``dest`` (never a checkout)."""
    dest.mkdir(parents=True, exist_ok=True)
    archive = subprocess.run(["git", "archive", "--format=tar", commit], cwd=str(ROOT),
                             capture_output=True, check=True)
    subprocess.run(["tar", "-x", "-C", str(dest)], input=archive.stdout, check=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--smiles", type=Path, default=DEFAULT_SMILES)
    ap.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    ap.add_argument("--tree", type=Path, default=ROOT,
                    help="tree under test (default: the repo this script lives in)")
    ap.add_argument("--ref-commit", default=None,
                    help="build this commit too and compare against it")
    ap.add_argument("--record", action="store_true",
                    help="write the manifest from the tree under test")
    ap.add_argument("--out", type=Path, default=None,
                    help="where outputs go (default: a temporary directory)")
    ap.add_argument("--keep-ref-outputs", type=Path, default=None,
                    help="with --record: keep the reference artefacts here; "
                         "in a compare: read them from here for the first-diff lines")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=3600.0)
    ap.add_argument("--only", default=None, help="comma-separated ids to build")
    ap.add_argument("--quick", type=float, default=None, metavar="SECONDS",
                    help="only the ids whose recorded build took at most this long")
    args = ap.parse_args(argv)

    rows = read_smiles(args.smiles)
    if args.only:
        wanted = set(args.only.split(","))
        rows = [r for r in rows if r[0] in wanted]
    if args.quick is not None:
        recorded = json.loads(args.manifest.read_text())["results"]
        rows = [r for r in rows
                if recorded.get(r[0], {}).get("seconds", 1e9) <= args.quick]
    print(f"{len(rows)} SMILES, interpreter {sys.executable}")

    tmp = None
    out = args.out
    if out is None:
        tmp = tempfile.mkdtemp(prefix="split_identity_")
        out = Path(tmp)
    out.mkdir(parents=True, exist_ok=True)
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(args.tree),
                              capture_output=True, text=True).stdout.strip()
        new_dir = out / "tree"
        print(f"building on tree {args.tree} (HEAD {head[:12]})")
        new = build_set(args.tree, rows, new_dir, args.workers, args.timeout)

        if args.record:
            manifest = {
                "recorded_at_commit": head,
                "construction": new["construction"],
                "results": new["results"],
            }
            args.manifest.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n")
            n_ok = sum(1 for r in new["results"].values() if r["status"] == "ok")
            print(f"recorded {n_ok}/{len(rows)} ok builds into {args.manifest}")
            if args.keep_ref_outputs is not None:
                if args.keep_ref_outputs.exists():
                    shutil.rmtree(args.keep_ref_outputs)
                shutil.copytree(new_dir, args.keep_ref_outputs)
            return 0 if n_ok == len(rows) else 1

        if args.ref_commit:
            ref_tree = out / "ref_tree"
            export_commit(args.ref_commit, ref_tree)
            ref_dir = out / "ref"
            print(f"building reference commit {args.ref_commit}")
            ref = build_set(ref_tree, rows, ref_dir, args.workers, args.timeout)
            print(f"\ncompare tree vs commit {args.ref_commit}:")
            return compare(ref["results"], new["results"], ref_dir, new_dir)

        manifest = json.loads(args.manifest.read_text())
        ref_dir = args.keep_ref_outputs
        print(f"\ncompare tree vs manifest (recorded at {manifest['recorded_at_commit'][:12]}):")
        if manifest["construction"] != new["construction"]:
            print("  WARNING: construction env differs from the recorded one")
        return compare(manifest["results"], new["results"], ref_dir, new_dir)
    finally:
        if tmp is not None and os.environ.get("SPLIT_IDENTITY_KEEP") != "1":
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())
