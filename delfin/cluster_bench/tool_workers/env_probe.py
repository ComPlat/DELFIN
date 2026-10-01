#!/usr/bin/env python
"""Print the provenance of the interpreter that runs a tool, as one JSON object.

    python env_probe.py

Runs in every tool environment (python >= 3.7, no DELFIN import).  Reports the interpreter,
glibc, and the version of every package a worker may use.  For epic-MACE also the sha256 over
its installed ``mace/*.py``: its GitHub state carries the same version string as the release,
so only the files tell them apart.  ``delfin cluster`` records this in the manifest when the run
is prepared and again in every chunk, and refuses a chunk whose environment differs.
"""
import hashlib
import json
import os
import platform
import sys

try:
    import importlib.metadata as _md
except ImportError:  # python 3.7
    _md = None

PACKAGES = ("rdkit", "rdkit-pypi", "numpy", "scipy", "architector", "molSimplify", "mendeleev",
            "openbabel", "openbabel-wheel", "xtb", "tblite", "ase", "numba", "epic-mace", "pyyaml",
            "networkx", "stk")


def probe_version(name):
    if _md is not None:
        try:
            return _md.version(name)
        except Exception:
            return None
    try:
        import pkg_resources
        return pkg_resources.get_distribution(name).version
    except Exception:
        return None


def probe_mace_files():
    try:
        import mace
    except Exception:
        return None
    root = os.path.dirname(os.path.abspath(mace.__file__))
    h = hashlib.sha256()
    for name in sorted(os.listdir(root)):
        if name.endswith(".py"):
            h.update(name.encode())
            with open(os.path.join(root, name), "rb") as fh:
                h.update(hashlib.sha256(fh.read()).hexdigest().encode())
    return h.hexdigest()


def main():
    rep = {"python": platform.python_version(), "executable": sys.executable,
           "glibc": "-".join(platform.libc_ver()), "machine": platform.machine()}
    rep["packages"] = {p: v for p, v in ((p, probe_version(p)) for p in PACKAGES) if v}
    try:
        from rdkit import rdBase
        rep["rdkit_runtime"] = rdBase.rdkitVersion
    except Exception:
        rep["rdkit_runtime"] = None
    mf = probe_mace_files()
    if mf:
        rep["mace_files_sha256"] = mf
    print(json.dumps(rep, sort_keys=True))


if __name__ == "__main__":
    main()
