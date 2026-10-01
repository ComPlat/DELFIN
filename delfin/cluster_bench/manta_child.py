"""Build ONE MANTA system in its own process -- the child of ``delfin cluster run-shard``.

    python manta_child.py REPO_ROOT ID SMILES OUT_DIR CONFIG THREADS

Run by path, never imported: the construction switches must be in the environment BEFORE
``delfin`` is imported, and a fresh process per system keeps every build independent of the one
before it (module caches, RDKit state, peak memory).

What the build is -- the same call the campaign harness of the MANTA paper makes for every
system, so a cluster run and a local reference run give the same bytes:

* the construction switches of ``delfin.cli_manta.construction_env(CONFIG)`` (``champion`` =
  the shipped MANTA construction).  If they cannot be read the build STOPS (exit 3): a guessed
  configuration would give numbers that look like an answer and measure something else;
* every ``DELFIN_*`` switch whose name marks a post-construction pass (REPAIR, RESTORE, FIXUP,
  POSTHOC) is set to 0, so only the raw construction is measured (none is set by the champion;
  this only guards against an inherited switch);
* ``smiles_to_xyz_isomers(SMILES, apply_uff=True, collapse_label_variants=False,
  include_binding_mode_isomers=True, deterministic=True, max_isomers=100000,
  quality_mode="extreme")`` -- the complete manifold, as the ``delfin-manta`` CLI ships it;
* threads: ``DELFIN_MAX_THREAD_WORKERS`` = THREADS, BLAS/OpenMP pinned to 1 (multi-threaded
  reductions re-associate float sums, and the conformer dedup decides on an RMSD threshold);
  ``DELFIN_MAX_PROCESS_WORKERS`` = 1.

Output: ``OUT_DIR/<ID>.xyz`` (one multi-frame xyz, header ``<ID> frame<k> <label>``) and one
JSON line on stdout: ``{"rid", "status": "ok"|"empty", "niso", "peak_gb"}``.  Nothing is
written before the build is complete, so a killed build leaves no file.
"""
import json
import os
import sys

REPO, RID, SMI, OUT, CONFIG, THREADS = sys.argv[1:7]

# Started by path, Python puts this directory first on sys.path; its module names (cli, runner,
# ...) must never shadow anything the build imports.
_here = os.path.dirname(os.path.abspath(__file__))
sys.path[:] = [p for p in sys.path if os.path.abspath(p or ".") != _here]
sys.path.insert(0, REPO)

os.environ["DELFIN_UI_INLINE"] = "1"
os.environ["DELFIN_MAX_PROCESS_WORKERS"] = "1"
os.environ["DELFIN_MAX_THREAD_WORKERS"] = str(int(THREADS))
for _t in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_t, "1")

from rdkit import RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")
try:
    from delfin.cli_manta import _apply_construction_env  # noqa: E402

    _apply_construction_env(CONFIG)
except Exception as _e:  # noqa: BLE001
    sys.stderr.write("FATAL: cannot read the construction configuration %r "
                     "(delfin.cli_manta.construction_env): %s: %s\n" % (CONFIG, type(_e).__name__, _e))
    raise SystemExit(3)

for _k in list(os.environ):
    if _k.startswith("DELFIN_") and any(_t in _k for _t in ("REPAIR", "RESTORE", "FIXUP", "POSTHOC")):
        os.environ[_k] = "0"

import delfin.smiles_converter as sc  # noqa: E402


def _cb_child_peak_gb():
    try:
        import resource
        return round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0, 3)
    except Exception:  # noqa: BLE001
        return None


iso, _err = sc.smiles_to_xyz_isomers(
    SMI, apply_uff=True, collapse_label_variants=False, include_binding_mode_isomers=True,
    deterministic=True, max_isomers=100000, quality_mode="extreme")
if not iso:
    print(json.dumps({"rid": RID, "status": "empty", "peak_gb": _cb_child_peak_gb()}))
    raise SystemExit(0)
os.makedirs(OUT, exist_ok=True)
_drop_identical = os.environ.get("DELFIN_FFFREE_DROP_IDENTICAL_FRAMES", "0") == "1"
_tmp = os.path.join(OUT, ".%s.xyz.part" % RID)
with open(_tmp, "w") as fh:
    _seen = set()
    _wi = 0
    for xyz, lbl in iso:
        atoms = [ln for ln in str(xyz).splitlines() if len(ln.split()) == 4]
        if _drop_identical:
            _blk = "\n".join(atoms)
            if _blk in _seen:
                continue
            _seen.add(_blk)
        fh.write(f"{len(atoms)}\n{RID} frame{_wi} {lbl}\n" + "\n".join(atoms) + "\n")
        _wi += 1
os.replace(_tmp, os.path.join(OUT, RID + ".xyz"))
print(json.dumps({"rid": RID, "status": "ok", "niso": len(iso), "peak_gb": _cb_child_peak_gb()}))
