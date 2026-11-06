"""Agent-side MANTA checks: parity gates, load reproducibility, raw geometry audit.

Phase 1 of the mt1 task: tools only, no constructor or eye change.

This module ports the *semantics* of the private audit scripts into clean,
configurable DELFIN code.  All paths come from arguments; no user names,
archive names or home paths are hard-coded, and no script data is shipped.

Everything calls MANTA only through its public entry points: the CLI
(``delfin-manta`` or ``python -m delfin.cli_manta``) or
:func:`delfin.cli_manta.construction_env`.  No private copy of the
construction switches is made here and no ``DELFIN_*`` construction variable
is ever set globally in this process -- the environment dicts are handed to
child processes only.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

__all__ = [
    "blocks_of",
    "compare_build_to_archive",
    "parity_gate",
    "parity_env",
    "load_repro",
    "geom_audit_frames",
    "geom_compare_table",
    "main",
]


# --------------------------------------------------------------------------
# Frame/block parsing (parity_compare semantics)
# --------------------------------------------------------------------------

def blocks_of(text: str) -> Tuple[str, ...]:
    """XYZ blocks of a multi-frame file as joined coordinate-line tuples.

    A block starts at a line that is a lone integer (the atom count), then one
    comment line, then that many coordinate lines.  Each coordinate line is
    normalised by collapsing internal whitespace so that byte-noise does not
    produce phantom differences.
    """
    out: List[str] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        tok = lines[i].split()
        if len(tok) == 1 and tok[0].isdigit():
            n = int(tok[0])
            out.extend(
                " ".join(x.split()) for x in lines[i + 2:i + 2 + n]
            )
            i += 2 + n
        else:
            i += 1
    return tuple(out)


def _manifest_frames(out_dir: Path) -> List[str]:
    """Per-isomer coordinate lines of one build directory, manifest order."""
    man = out_dir / "manifest.json"
    if not man.exists():
        return []
    lines: List[str] = []
    for e in json.loads(man.read_text())["isomers"]:
        f = out_dir / e["file"]
        if f.exists():
            lines.extend(blocks_of(f.read_text()))
    return lines


def compare_build_to_archive(
    out_dir: Path, archive_xyz: Path
) -> Dict[str, Any]:
    """Classify one built system against one reference archive frame set.

    Returns a dict with keys ``system``, ``verdict`` (one of ``identical``,
    ``same_set_other_order``, ``different``, ``missing``), ``n_cli``,
    ``n_ref`` and ``n_common``.
    """
    rid = archive_xyz.stem
    man = out_dir / "manifest.json"
    if not man.exists() or not archive_xyz.exists():
        return {"system": rid, "verdict": "missing", "n_cli": 0, "n_ref": 0,
                "n_common": 0}
    cli = _manifest_frames(out_dir)
    ref = blocks_of(archive_xyz.read_text())
    if cli == list(ref):
        k = "identical"
    elif sorted(cli) == sorted(ref):
        k = "same_set_other_order"
    else:
        k = "different"
    return {
        "system": rid,
        "verdict": k,
        "n_cli": len(cli),
        "n_ref": len(ref),
        "n_common": len(set(cli) & set(ref)),
    }


def _pool_lines(ids_file: Path, pool_file: Path) -> List[str]:
    """Pool records (``rid|smiles``) whose rid appears in the id list."""
    ids = set(ids_file.read_text().split())
    out = []
    for ln in pool_file.read_text().splitlines():
        rid = ln.replace(";", "|").split("|")[0].strip()
        if rid and rid in ids:
            smi = ln.split("|", 1)[1].strip() if "|" in ln else ""
            out.append(rid + "|" + smi)
    return out


def _inputs_from_ids(ids_file: Path, pool_file: Optional[Path]) -> List[Tuple[str, str]]:
    """``(rid, smiles)`` pairs: from the pool when given, else SMILES per line.

    With a pool file, the record whose rid matches the id list wins.  Without
    one, each non-empty line of the id file is ``rid`` and the SMILES is
    expected after the first ``|`` separator (``rid|smiles``); a line without a
    separator is treated as a bare SMILES with a sequential ``sysN`` id.
    """
    if pool_file is not None:
        out = []
        for ln in _pool_lines(ids_file, pool_file):
            rid, _, smi = ln.partition("|")
            out.append((rid.strip(), smi.strip()))
        return out
    out = []
    for k, ln in enumerate(ids_file.read_text().splitlines()):
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        rid, sep, smi = ln.partition("|")
        out.append((rid.strip() if sep else f"sys{k}", smi.strip() if sep else ln))
    return out


def _run_cli_build(
    command: Sequence[str],
    smiles: str,
    out_dir: Path,
    env_overrides: Dict[str, str],
    log_path: Path,
    nice: Optional[int],
) -> int:
    """One public-CLI build as a child process; returns the exit code.

    The construction environment (a dict from ``construction_env`` or user
    VAR=VALUE pairs) is applied to the child environment only -- never to this
    process.
    """
    out_dir.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update(env_overrides)
    cmd = list(command)
    if nice:
        cmd = ["nice", "-n", str(nice)] + cmd
    proc = subprocess.run(
        cmd + [smiles, "-o", str(out_dir)],
        env=env, capture_output=True, text=True, timeout=None,
    )
    log_path.write_text(proc.stdout + proc.stderr)
    return proc.returncode


def _run_parity_batch(
    inputs: Sequence[Tuple[str, str]],
    work: Path,
    command_builder: Callable[[Path], Sequence[str]],
    env_overrides: Dict[str, str],
    parallel: int,
    nice: Optional[int],
    ref_archive: Path,
) -> List[Dict[str, Any]]:
    """Run builds over the inputs, then classify each against the archive."""
    rows: List[Dict[str, Any]] = []

    def build(item: Tuple[str, str]) -> None:
        rid, smi = item
        out_dir = work / "out" / rid
        if (out_dir / "manifest.json").exists():
            return
        cmd = command_builder(out_dir)
        _run_cli_build(cmd, smi, out_dir, env_overrides,
                       work / "out" / f"{rid}.log", nice)

    parallel = max(1, int(parallel))
    if parallel == 1:
        for item in inputs:
            build(item)
    else:
        with ThreadPoolExecutor(max_workers=parallel) as ex:
            list(ex.map(build, inputs))
    for rid, _smi in inputs:
        rows.append(compare_build_to_archive(work / "out" / rid,
                                             ref_archive / f"{rid}.xyz"))
    return rows


def _summary(rows: Iterable[Dict[str, Any]]) -> Dict[str, int]:
    cnt = {"identical": 0, "same_set_other_order": 0, "different": 0,
           "missing": 0}
    for r in rows:
        cnt[r["verdict"]] = cnt.get(r["verdict"], 0) + 1
    return cnt


# --------------------------------------------------------------------------
# Tool 1: parity_gate -- fresh clone + clean venv + public CLI vs an archive
# --------------------------------------------------------------------------

def parity_gate(
    ref: str,
    ids_file: Path,
    ref_archive: Path,
    work: Path,
    *,
    parallel: int = 1,
    base_repo: Optional[str] = None,
    python_for_venv: Optional[str] = None,
    pool_file: Optional[Path] = None,
    extra_env: Optional[Dict[str, str]] = None,
    nice: Optional[int] = 10,
    clone_dir: Optional[Path] = None,
    install_cmd: Optional[Sequence[str]] = None,
    allow_install: bool = False,
    command: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Plan and (on ``allow_install=True``) run the clone/install/parity gate.

    Dry-run by default: with ``allow_install=False`` nothing is cloned or
    installed and nothing is built; the plan plus the comparison of what is
    already in ``work`` against the archive is reported.

    * ``ref``: a commit hash or branch name checked out in the clone.
    * ``ids_file``: systems to build (rid list, or ``rid|smiles`` lines, or a
      pool file via ``pool_file``).
    * ``ref_archive``: directory of ``<rid>.xyz`` multi-frame reference files.
    * ``work``: scratch directory for the clone, venv and build outputs.
    * ``base_repo``: the DELFIN checkout to clone from (default: the parent of
      this file's tree, ``git rev-parse`` from cwd when missing).
    * ``python_for_venv``: interpreter used to create the clean venv.
    * ``install_cmd``: full command used to install the clone into the venv
      (default ``[pip, install, ./clone]``).
    * ``command``: build command override (tests pass a fake CLI here; the
      real gate runs ``<venv>/bin/delfin-manta`` -- the public CLI).
    * ``allow_install``: permission flag; clone+install+builds run only when
      it is true.
    """
    work = Path(work)
    ids_file = Path(ids_file)
    ref_archive = Path(ref_archive)
    plan = {
        "ref": ref,
        "ids_file": str(ids_file),
        "ref_archive": str(ref_archive),
        "work": str(work),
        "parallel": parallel,
        "allow_install": allow_install,
        "base_repo": base_repo,
        "python_for_venv": python_for_venv,
    }
    inputs = _inputs_from_ids(ids_file, pool_file)
    plan["n_systems"] = len(inputs)
    if not allow_install:
        plan["mode"] = "dry_run"
        rows = [
            compare_build_to_archive(work / "out" / rid,
                                     ref_archive / f"{rid}.xyz")
            for rid, _smi in inputs
        ]
        plan["rows"] = rows
        plan["summary"] = _summary(rows)
        return plan

    clone = Path(clone_dir) if clone_dir else work / "clone"
    plan["clone_dir"] = str(clone)
    if not clone.exists():
        base_repo = base_repo or _infer_base_repo()
        subprocess.run(
            ["git", "clone", "-q", "--no-checkout", base_repo, str(clone)],
            check=True, capture_output=True, text=True)
        subprocess.run(
            ["git", "-C", str(clone), "checkout", "-q", ref],
            check=True, capture_output=True, text=True)
    venv_py = _ensure_venv(work, python_for_venv, install_cmd, clone)
    build_cmd = list(command) if command is not None else \
        [str(venv_py.parent / "delfin-manta")]
    plan["mode"] = "executed"
    rows = _run_parity_batch(
        inputs, work,
        command_builder=lambda out: build_cmd,
        env_overrides=dict(extra_env or {}),
        parallel=parallel, nice=nice, ref_archive=ref_archive)
    plan["rows"] = rows
    plan["summary"] = _summary(rows)
    return plan


def _infer_base_repo() -> str:
    """Repo root of the checkout this module runs from (no hard-coded path)."""
    import delfin

    return str(Path(delfin.__file__).resolve().parent.parent)


def _ensure_venv(
    work: Path,
    python_for_venv: Optional[str],
    install_cmd: Optional[Sequence[str]],
    clone: Path,
) -> Path:
    """Create ``work/venv`` and install the clone into it; return python path."""
    venv = work / "venv"
    py = venv / "bin" / "python"
    if not py.exists():
        base = python_for_venv or sys.executable
        subprocess.run([base, "-m", "venv", str(venv)],
                       check=True, capture_output=True, text=True)
    cmd = list(install_cmd or [str(venv / "bin" / "pip"), "install",
                               str(clone)])
    if cmd and cmd[0].endswith("pip"):
        pass  # explicit pip command as given
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return py


# --------------------------------------------------------------------------
# Tool 2: parity_env -- same comparison on a chosen interpreter, no clone
# --------------------------------------------------------------------------

def parity_env(
    ids_file: Path,
    ref_archive: Path,
    work: Path,
    *,
    pybin: Optional[str] = None,
    source_dir: Optional[str] = None,
    parallel: int = 1,
    pool_file: Optional[Path] = None,
    extra_env: Optional[Dict[str, str]] = None,
    nice: Optional[int] = 10,
    allow_builds: bool = True,
    command: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Run the public CLI (``python -m delfin.cli_manta``) from a source tree
    under a chosen interpreter and compare against the reference archive.

    No clone, no install: the source tree is put on the child's PYTHONPATH.
    * ``command``: override of the build command (tests pass a fake CLI); the
      default is the public entry ``python -m delfin.cli_manta``.
    With ``allow_builds=False`` the call only classifies what already exists
    in ``work/out`` (dry-run flavour).
    """
    ids_file = Path(ids_file)
    ref_archive = Path(ref_archive)
    work = Path(work)
    plan = {
        "mode": "executed" if allow_builds else "dry_run",
        "ids_file": str(ids_file),
        "ref_archive": str(ref_archive),
        "work": str(work),
        "pybin": pybin or sys.executable,
        "source_dir": source_dir,
        "parallel": parallel,
    }
    inputs = _inputs_from_ids(ids_file, pool_file)
    plan["n_systems"] = len(inputs)
    if not allow_builds:
        rows = [
            compare_build_to_archive(work / "out" / rid,
                                     ref_archive / f"{rid}.xyz")
            for rid, _smi in inputs
        ]
    else:
        env = {"PYTHONPATH": source_dir} if source_dir else {}
        env.update(extra_env or {})
        base_cmd = list(command) if command is not None else \
            [pybin or sys.executable, "-m", "delfin.cli_manta"]
        rows = _run_parity_batch(
            inputs, work,
            command_builder=lambda out: base_cmd,
            env_overrides=env,
            parallel=parallel, nice=nice, ref_archive=ref_archive)
    plan["rows"] = rows
    plan["summary"] = _summary(rows)
    return plan


# --------------------------------------------------------------------------
# Tool 3: load_repro -- N parallel identical builds; empty frames + distinct
# --------------------------------------------------------------------------

_LOAD_REPRO_TEMPLATE = "tests/test_cli_dashboard_parity.py"  # side-script shape


def load_repro(
    smiles: str,
    work: Path,
    *,
    n: int = 2,
    config: str = "champion",
    command: Optional[Sequence[str]] = None,
    parallel: Optional[int] = None,
    timeout_s: float = 900.0,
) -> Dict[str, Any]:
    """Run ``n`` (<= 16) identical public-CLI builds in parallel.

    Returns the number of empty builds (no manifest / no isomers) and the
    number of distinct coordinate-line sets across the builds.  Deterministic
    construction must give one distinct result; a run with more is a
    reproducibility failure.

    ``command`` defaults to ``[sys.executable, -m, delfin.cli_manta]``; tests
    pass a fake CLI here.  The construction environment for each child comes
    from :func:`delfin.cli_manta.construction_env` applied to the child
    environment only.
    """
    work = Path(work)
    n = int(n)
    if n < 1 or n > 16:
        raise ValueError("load_repro needs 1 <= n <= 16")
    workers = parallel or n
    from delfin import cli_manta  # public entry point, switches come from here

    env_base = {k: v for k, v in os.environ.items()
                if not k.startswith("DELFIN_")}

    def one(i: int) -> Dict[str, Any]:
        out = work / f"build_{i}"
        cmd = list(command or [sys.executable, "-m", "delfin.cli_manta"])
        child_env = dict(env_base)
        child_env.update(
            cli_manta.construction_env(config, environ=env_base))
        out.mkdir(parents=True, exist_ok=True)
        proc = subprocess.run(
            cmd + [smiles, "-o", str(out)], env=child_env,
            capture_output=True, text=True, timeout=timeout_s)
        (work / f"build_{i}.log").write_text(proc.stdout + proc.stderr)
        empty = not (out / "manifest.json").exists()
        frames: List[str] = []
        if not empty:
            frames = _manifest_frames(out)
            empty = not frames
        return {"index": i, "exit": proc.returncode, "empty": empty,
                "frames": frames}

    with ThreadPoolExecutor(max_workers=workers) as ex:
        runs = list(ex.map(one, range(n)))
    distinct = len({tuple(r["frames"]) for r in runs})
    return {
        "n": n,
        "n_empty": sum(1 for r in runs if r["empty"]),
        "n_distinct": distinct,
        "distinct_results": distinct,
        "n_failed": sum(1 for r in runs if r["exit"] != 0),
        "runs": [{"index": r["index"], "exit": r["exit"], "empty": r["empty"]}
                 for r in runs],
    }


# --------------------------------------------------------------------------
# Tool 4: geom_audit -- RAW geometric measures only (no CSD table needed)
# --------------------------------------------------------------------------

def _vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _vec_norm(a):
    return math.sqrt(_vec_dot(a, a))


def _vec_angle(a, b) -> float:
    """Angle between two vectors in degrees."""
    na, nb = _vec_norm(a), _vec_norm(b)
    if na == 0 or nb == 0:
        return float("nan")
    c = _vec_dot(a, b) / (na * nb)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def _dist(a, b) -> float:
    return _vec_norm(_vec_sub(a, b))


# covalent radii used for the raw contact measure; extendable by the caller
_RCOV = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
    "Sc": 1.70, "Ti": 1.60, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32,
    "Co": 1.26, "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Y": 1.90, "Zr": 1.75,
    "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42, "Pd": 1.39,
    "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "La": 2.07, "Ce": 2.04,
    "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36,
    "Hg": 1.32, "Pb": 1.46,
}
_VDW = {"H": 1.20, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
        "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98}


def parse_xyz(text: str) -> Tuple[List[str], List[Tuple[float, float, float]]]:
    """First XYZ frame of ``text``: element list and coordinates."""
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        tok = ln.split()
        if len(tok) == 1 and tok[0].isdigit():
            n = int(tok[0])
            syms, xyz = [], []
            for cl in lines[i + 2:i + 2 + n]:
                t = cl.split()
                syms.append(t[0])
                xyz.append((float(t[1]), float(t[2]), float(t[3])))
            return syms, xyz
    raise ValueError("no XYZ frame found")


_GEOM_UNITS = ("bond_len", "m_contact", "nb_clash", "angle")


def geom_audit_frames(
    frames: Iterable[Tuple[str, Sequence[str], Sequence[Tuple[float, float, float]]]],
    *,
    bonds: Optional[Dict[int, List[int]]] = None,
    rcov: Optional[Dict[str, float]] = None,
    vdw: Optional[Dict[str, float]] = None,
    use_rdkit: bool = True,
) -> List[Dict[str, Any]]:
    """RAW geometric measures of every frame -- no threshold table involved.

    Each frame is a ``(frame_id, symbols, coords)`` triple.  Returns unit rows
    (dicts) with ``frame_id``, ``quantity`` and ``value``:

    * ``bond_len`` -- bonded heavy-pair distance normalised by the covalent
      radii sum (class key ``bond_len``, two-sided quantity).
    * ``m_contact`` -- heavy atom not bonded to a metal / metal distances,
      normalised by the covalent radii sum, low side.
    * ``nb_clash`` -- heavy pairs >= 4 bonds apart normalised by the van der
      Waals sum, low side.
    * ``angle`` -- bond angle at each middle atom over two of its bonded
      neighbours (hydrogens included; class key ``angle``).

    When ``bonds`` is not given, bonds come from RDKit (if installed and
    ``use_rdkit``) reading a SMILES-free distance graph; without RDKit a
    covalent-distance graph is built from the coordinates themselves.  The
    measures are pure geometry: no chemistry class assignment beyond element
    identity, and the optional COD comparison is a separate call.
    """
    rows: List[Dict[str, Any]] = []
    rcov = dict(_RCOV if rcov is None else rcov)
    vdw = dict(_VDW if vdw is None else vdw)
    metal_hint = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                  "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                  "In", "Sn", "La", "Ce", "W", "Re", "Os", "Ir", "Pt", "Au",
                  "Hg", "Pb"}
    for frame_id, syms, xyz in frames:
        n = len(syms)
        if n < 2:
            continue
        bmap = _bonds_for(syms, xyz, bonds, use_rdkit)
        heavy = [i for i in range(n) if syms[i] != "H"]
        # bond_len: bonded heavy pairs
        for i in heavy:
            for j in bmap.get(i, []):
                if j <= i or syms[j] == "H":
                    continue
                d = _dist(xyz[i], xyz[j])
                s = rcov.get(syms[i], 0.8) + rcov.get(syms[j], 0.8)
                if s > 0:
                    rows.append({"frame_id": frame_id, "quantity": "bond_len",
                                 "value": d / s, "class": "bond_len",
                                 "i": i, "j": j})
        # m_contact: heavy non-metal atom to metal distances (not bonded)
        metals = [i for i in range(n) if syms[i] in metal_hint]
        for m in metals:
            for x in heavy:
                if x == m or syms[x] in metal_hint:
                    continue
                if x in bmap.get(m, []) or m in bmap.get(x, []):
                    continue
                d = _dist(xyz[m], xyz[x])
                s = rcov.get(syms[m], 1.3) + rcov.get(syms[x], 0.8)
                if s > 0:
                    rows.append({"frame_id": frame_id, "quantity": "m_contact",
                                 "value": d / s, "class": "m_contact",
                                 "i": m, "j": x})
        # nb_clash: heavy pairs >= 4 bonds apart (graph hops), vdW normalised
        hops = _hops(n, bmap)
        for a_i, x in enumerate(heavy):
            for y in heavy[a_i + 1:]:
                if hops.get((x, y), hops.get((y, x))) is not None and \
                        hops.get((x, y), hops.get((y, x))) < 4:
                    continue
                d = _dist(xyz[x], xyz[y])
                s = vdw.get(syms[x], 1.7) + vdw.get(syms[y], 1.7)
                if s > 0:
                    rows.append({"frame_id": frame_id, "quantity": "nb_clash",
                                 "value": d / s, "class": "nb_clash",
                                 "i": x, "j": y})
        # angle: bond angle at each atom over bonded-neighbour pairs (H
        # included -- a raw geometric measure needs no chemistry classification).
        for i in range(n):
            nb = sorted(bmap.get(i, []))
            if len(nb) < 2:
                continue
            for a1 in range(len(nb)):
                for a2 in range(a1 + 1, len(nb)):
                    x, y = nb[a1], nb[a2]
                    v1 = _vec_sub(xyz[x], xyz[i])
                    v2 = _vec_sub(xyz[y], xyz[i])
                    ang = _vec_angle(v1, v2)
                    if not math.isnan(ang):
                        rows.append({"frame_id": frame_id, "quantity": "angle",
                                     "value": ang, "class": "angle",
                                     "i": i, "j": x, "k": y})
    return rows


def _bonds_for(syms, xyz, bonds, use_rdkit):
    if bonds is not None:
        return bonds
    if use_rdkit:
        try:
            bmap = _rdkit_bonds(syms, xyz)
            if bmap:
                return bmap
        except Exception:
            pass
    return _distance_bonds(syms, xyz)


def _rdkit_bonds(syms, xyz):
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    body_lines = [f"{s} {c[0]:.4f} {c[1]:.4f} {c[2]:.4f}"
                  for s, c in zip(syms, xyz)]
    text = f"{len(syms)}\ngeom_audit\n" + "\n".join(body_lines) + "\n"
    mol = Chem.MolFromXYZBlock(text)
    if mol is None:
        raise ValueError("RDKit could not parse the XYZ frame")
    bmap: Dict[int, List[int]] = {}
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bmap.setdefault(i, []).append(j)
        bmap.setdefault(j, []).append(i)
    return bmap


def _distance_bonds(syms, xyz):
    bmap: Dict[int, List[int]] = {}
    n = len(syms)
    for i in range(n):
        for j in range(i + 1, n):
            s = _RCOV.get(syms[i], 0.8) + _RCOV.get(syms[j], 0.8)
            if s > 0 and _dist(xyz[i], xyz[j]) <= 1.25 * s:
                bmap.setdefault(i, []).append(j)
                bmap.setdefault(j, []).append(i)
    return bmap


def _hops(n, bmap) -> Dict[Tuple[int, int], int]:
    """BFS hop distances between all heavy-atom-relevant pairs (small graphs)."""
    hops: Dict[Tuple[int, int], int] = {}
    for src in range(n):
        seen = {src: 0}
        frontier = [src]
        while frontier:
            nxt = []
            for u in frontier:
                for v in bmap.get(u, []):
                    if v not in seen:
                        seen[v] = seen[u] + 1
                        nxt.append(v)
            frontier = nxt
        for tgt, d in seen.items():
            hops[(src, tgt)] = d
    return hops


def geom_compare_table(
    rows: Sequence[Dict[str, Any]],
    table_path: Path,
    *,
    min_n: int = 300,
) -> Dict[str, Any]:
    """OPTIONAL comparison of raw unit rows against a user-supplied table.

    The table (JSON, derived from COD or the user's own data -- never CSD)
    maps ``quantity -> class -> {"p99": x, "p999": y}`` or directly a list of
    unit values from which those percentiles are computed.  A frame is
    flagged ``sev`` when any of its unit values crosses the severe side.
    Direction follows the quantity: ``m_contact``/``nb_clash`` flag low
    values, everything else flags high values.
    """
    table = json.loads(Path(table_path).read_text())
    q_rows: Dict[str, List[Dict[str, Any]]] = {}
    for r in rows:
        q_rows.setdefault(r["quantity"], []).append(r)

    def th(q: str, cls: str) -> Optional[Tuple[float, float]]:
        e = table.get(q, {}).get(cls)
        if e is None:
            e = table.get(q, {}).get("*")
        if e is None:
            return None
        if isinstance(e, dict):
            return (e.get("p99"), e.get("p999"))
        vals = sorted(float(v) for v in e)
        if len(vals) < min_n:
            return None
        return (vals[int(0.99 * (len(vals) - 1))],
                vals[int(0.999 * (len(vals) - 1))])

    low_side = {"m_contact", "nb_clash"}
    frames: Dict[str, Dict[str, Any]] = {}
    for q, rs in q_rows.items():
        for r in rs:
            cls = r.get("class", q)
            t = th(q, cls)
            if t is None:
                continue
            p99, p999 = t
            hit = None
            if q in low_side and r["value"] < p999:
                hit = "sev" if r["value"] < p99 else "warn"
            elif q not in low_side and r["value"] > p999:
                hit = "sev" if r["value"] > p99 * 1.0 else "warn"
            if hit:
                frames.setdefault(r["frame_id"], {})[q] = hit
    return {
        "n_units": len(rows),
        "n_frames_flagged": len(frames),
        "frames": frames,
        "table": str(table_path),
    }


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        prog="delfin-manta-checks",
        description="MANTA parity / load-repro / raw-geometry agent checks")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("parity_gate", help="clone+venv+CLI vs archive "
                                            "(dry-run by default)")
    p.add_argument("--ref", required=True)
    p.add_argument("--ids", required=True, type=Path)
    p.add_argument("--archive", required=True, type=Path)
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--pool", type=Path, default=None)
    p.add_argument("--parallel", type=int, default=1)
    p.add_argument("--base-repo", default=None)
    p.add_argument("--python-for-venv", default=None)
    p.add_argument("--run", action="store_true",
                   help="execute clone/install/builds (default: dry run)")

    p = sub.add_parser("parity_env", help="chosen interpreter, no clone")
    p.add_argument("--ids", required=True, type=Path)
    p.add_argument("--archive", required=True, type=Path)
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--pool", type=Path, default=None)
    p.add_argument("--pybin", default=None)
    p.add_argument("--source-dir", default=None)
    p.add_argument("--parallel", type=int, default=1)
    p.add_argument("--dry-run", action="store_true")

    p = sub.add_parser("load_repro", help="N <= 16 identical builds")
    p.add_argument("--smiles", required=True)
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--n", type=int, default=2)
    p.add_argument("--config", default="champion")

    p = sub.add_parser("geom_audit", help="raw geometric measures")
    p.add_argument("frames", nargs="+", type=Path, help=".xyz frame files")
    p.add_argument("--bonds-from-rdkit", action="store_true", default=True)
    p.add_argument("--no-rdkit", dest="bonds_from_rdkit", action="store_false")
    p.add_argument("--compare-table", type=Path, default=None)
    p.add_argument("--json", action="store_true")

    args = ap.parse_args(argv)
    if args.cmd == "parity_gate":
        res = parity_gate(args.ref, args.ids, args.archive, args.work,
                          parallel=args.parallel, pool_file=args.pool,
                          base_repo=args.base_repo,
                          python_for_venv=args.python_for_venv,
                          allow_install=args.run)
    elif args.cmd == "parity_env":
        res = parity_env(args.ids, args.archive, args.work,
                         pybin=args.pybin, source_dir=args.source_dir,
                         parallel=args.parallel, pool_file=args.pool,
                         allow_builds=not args.dry_run)
    elif args.cmd == "load_repro":
        res = load_repro(args.smiles, args.work, n=args.n,
                         config=args.config)
    else:
        frames = []
        for p in args.frames:
            syms, xyz = parse_xyz(Path(p).read_text())
            frames.append((Path(p).stem, syms, xyz))
        rows = geom_audit_frames(frames, use_rdkit=args.bonds_from_rdkit)
        res = {"n_units": len(rows), "units": rows}
        if args.compare_table:
            res["comparison"] = geom_compare_table(rows, args.compare_table)
    json.dump(res, sys.stdout, indent=2, sort_keys=True, default=str)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
