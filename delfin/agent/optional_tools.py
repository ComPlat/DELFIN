"""Optional chemistry programs: which are present, and does each one compute.

Input: the interpreter this module runs in and the binaries DELFIN's tool
resolver finds (``qm_runtime.find_tool_executable``, PATH as fallback).

Output, two forms over one table (:data:`GROUPS`):

* :func:`availability` -- present or not, without starting a computation
  (``importlib.util.find_spec`` and the resolver; no import of torch or
  pyscf). ``delfin-agent doctor`` shows one row per group from it.
* :func:`smoke` -- one tiny computation per group, about a second each on a
  Raspberry Pi 5: water RHF/STO-3G in PySCF, the energy of one harmonic bond
  in an OpenMM context, a 2x2 eigenproblem in torch with MACE's model class
  imported, ``xtb --version`` and a GFN2 water single point, an ORCA input
  written by OPI (ORCA itself is not needed for that). ``pi_install.sh``
  runs it after every deploy; ``python -m delfin.agent.optional_tools smoke``
  runs it by hand.

The groups are the optional install groups of ``scripts/pi/pi_install.sh``
(``DELFIN_PI_*`` switches), named here so the deploy, the doctor and the
README speak of the same things. Rejected: an import probe as the smoke
test. An import succeeds for a PySCF whose BLAS is missing a symbol, or an
xtb whose parameter files are not where the binary looks; a computation
against a known value does not.

Each Python computation runs in a child interpreter with one thread and a
timeout, so a crash, a hang or a 4-core BLAS start-up in one group cannot
take the others down, and a group's failure is reported, never raised.
"""

from __future__ import annotations

import importlib.util
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

__all__ = [
    "Group", "GROUPS", "OK", "PARTIAL", "MISSING", "FAILED", "OFF",
    "availability", "smoke", "off_text", "pi_hint", "format_smoke_line",
    "parse_smoke", "main",
]

OK = "ok"
PARTIAL = "partial"
MISSING = "missing"
FAILED = "failed"
OFF = "off"
STATUSES = (OK, PARTIAL, MISSING, FAILED, OFF)

#: Water at the experimental geometry (Angstrom), shared by every check.
WATER = (("O", 0.0, 0.0, 0.1173), ("H", 0.0, 0.7572, -0.4692),
         ("H", 0.0, -0.7572, -0.4692))

# Reference values, measured 2026-11-03 on x86-64 (PySCF 2.14.0, OpenMM
# 8.6.1, xtb 6.7.1 conda-forge). The tolerances admit another BLAS or CPU,
# not another method.
_RHF_STO3G_EH = -74.963023
_GFN2_EH = -5.070371
_BOND_KJ = 0.2          # 0.5 * 1000 kJ/mol/nm^2 * (0.10 - 0.12 nm)^2

_CHILD_TIMEOUT_S = 300   # torch + e3nn import on a Pi 5 SD card: tens of s
_BINARY_TIMEOUT_S = 60


@dataclass(frozen=True)
class Group:
    """One optional install group.

    ``modules`` are imports, ``binaries`` executables, ``dist`` the
    distribution whose version is reported. ``switch`` is the Pi's
    variable, ``install`` the line for any other host.
    """

    name: str
    purpose: str
    switch: str
    install: str
    modules: tuple[str, ...] = ()
    binaries: tuple[str, ...] = ()
    dist: str = ""
    _code: str = field(default="", repr=False)


_PYSCF = f"""
from pyscf import gto, scf
atoms = "; ".join(f"{{a}} {{x}} {{y}} {{z}}" for a, x, y, z in {WATER!r})
mol = gto.M(atom=atoms, basis="sto-3g", verbose=0)
e = scf.RHF(mol).kernel()
assert abs(e - ({_RHF_STO3G_EH})) < 1e-3, e
print(json.dumps({{"detail": f"water RHF/STO-3G {{e:.6f}} Eh"}}))
"""

_OPENMM = f"""
import openmm, openmm.unit as unit
system = openmm.System()
system.addParticle(1.0); system.addParticle(1.0)
bond = openmm.HarmonicBondForce(); bond.addBond(0, 1, 0.12, 1000.0)
system.addForce(bond)
platform = openmm.Platform.getPlatformByName("Reference")
ctx = openmm.Context(system, openmm.VerletIntegrator(0.001), platform)
ctx.setPositions([openmm.Vec3(0, 0, 0), openmm.Vec3(0.1, 0, 0)])
e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
assert abs(e - {_BOND_KJ}) < 1e-6, e
print(json.dumps({{"detail": f"2-atom bond {{e:.4f}} kJ/mol (OpenMM {{openmm.__version__}})"}}))
"""

_TORCH = """
import torch
torch.set_num_threads(1)
w = torch.linalg.eigvalsh(torch.tensor([[2.0, 1.0], [1.0, 2.0]], dtype=torch.float64))
assert [round(float(v), 9) for v in w] == [1.0, 3.0], w
from mace.modules import MACE  # noqa: F401  the model class, not a download
import mace
ver = getattr(mace, "__version__", "")
print(json.dumps({"detail": f"torch {torch.__version__} eigvalsh ok, MACE {ver} model class importable"}))
"""

_OPI = f"""
import tempfile, pathlib
from opi.core import Calculator
from opi.input import Input
from opi.input.simple_keywords import BasisSet, Method
from opi.input.structures import Atom, Structure
work = pathlib.Path(tempfile.mkdtemp(prefix="delfin-opi-smoke-"))
calc = Calculator("water", working_dir=work, version_check=False)
calc.input = Input()
calc.input.add_simple_keywords(Method.HF, BasisSet.DEF2_SVP)
calc.structure = Structure([Atom(a, coordinates=(x, y, z)) for a, x, y, z in {WATER!r}],
                           charge=0, multiplicity=1)
calc.write_input(force=True)
text = (work / "water.inp").read_text()
for f in work.iterdir():
    f.unlink()
work.rmdir()
assert "!hf" in text.lower() and "* xyz 0 1" in text, text
from importlib import metadata
print(json.dumps({{"detail": f"OPI {{metadata.version('orca-pi')}} wrote an HF/def2-SVP water input"}}))
"""

GROUPS: tuple[Group, ...] = (
    Group("opi", "ORCA's Python interface: builds ORCA inputs, reads ORCA outputs "
          "(no ORCA binary exists for ARM)",
          switch="DELFIN_PI_SKIP_EXTRAS",
          install='pip install "delfin-complat[orca]"',
          modules=("opi",), dist="orca-pi", _code=_OPI),
    Group("pyscf", "PySCF: HF/DFT/MP2/CC on small molecules",
          switch="DELFIN_PI_PYSCF", install="pip install pyscf",
          modules=("pyscf",), dist="pyscf", _code=_PYSCF),
    Group("openmm", "OpenMM: molecular dynamics, CPU platform",
          switch="DELFIN_PI_OPENMM", install="pip install openmm",
          modules=("openmm",), dist="openmm", _code=_OPENMM),
    Group("torch", "torch (CPU) + MACE: ML potentials",
          switch="DELFIN_PI_TORCH", install="pip install torch mace-torch",
          modules=("torch", "mace"), dist="torch", _code=_TORCH),
    Group("xtb", "xtb, CREST, tblite, dftd4 (conda-forge builds)",
          switch="DELFIN_PI_XTB", install="python -m delfin.installer --install xtb crest",
          binaries=("xtb", "crest", "tblite", "dftd4")),
)

_BY_NAME = {g.name: g for g in GROUPS}


def _find_binary(name: str) -> str | None:
    """The binary the way DELFIN runs it (qm_tools, *_BINARY, PATH)."""
    try:
        from delfin import qm_runtime
        found = qm_runtime.find_tool_executable(name)
    except Exception:
        found = None
    return found or shutil.which(name)


def _dist_version(dist: str) -> str:
    from importlib import metadata
    try:
        return metadata.version(dist)
    except metadata.PackageNotFoundError:
        return ""


def _module_present(module: str) -> bool:
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):
        return False


def availability(group: Group) -> tuple[str, str]:
    """(status, detail) without running anything: ok, partial or missing."""
    have = [m for m in group.modules if _module_present(m)]
    have += [b for b in group.binaries if _find_binary(b)]
    wanted = list(group.modules) + list(group.binaries)
    lacking = [w for w in wanted if w not in have]
    if not have:
        return MISSING, "not installed"
    version = _dist_version(group.dist) if group.dist else ""
    found = " ".join(have) + (f" {version}" if version else "")
    if lacking:
        return PARTIAL, f"{found}; missing: {' '.join(lacking)}"
    return OK, found


def _child_env() -> dict:
    env = dict(os.environ)
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[var] = "1"
    return env


def _run_python(code: str) -> tuple[bool, str]:
    """Run ``code`` in a child interpreter; (passed, detail)."""
    program = "import json\n" + code
    try:
        proc = subprocess.run(
            [sys.executable, "-c", program], capture_output=True, text=True,
            timeout=_CHILD_TIMEOUT_S, env=_child_env(), check=False,
            cwd=tempfile.gettempdir(),
        )
    except subprocess.TimeoutExpired:
        return False, f"no result within {_CHILD_TIMEOUT_S} s"
    except OSError as exc:
        return False, f"could not start {sys.executable}: {exc}"
    if proc.returncode == 0:
        for line in reversed(proc.stdout.strip().splitlines()):
            try:
                return True, str(json.loads(line)["detail"])
            except (ValueError, KeyError, TypeError):
                continue
        return False, "finished without a result line"
    tail = [ln for ln in (proc.stderr or proc.stdout).strip().splitlines() if ln.strip()]
    return False, (tail[-1] if tail else f"exit code {proc.returncode}")[:300]


def _version_line(path: str) -> str:
    try:
        proc = subprocess.run([path, "--version"], capture_output=True, text=True,
                              timeout=_BINARY_TIMEOUT_S, env=_child_env(), check=False)
    except (subprocess.TimeoutExpired, OSError) as exc:
        raise RuntimeError(f"--version: {exc}") from exc
    text = (proc.stdout or "") + (proc.stderr or "")
    m = re.search(r"[Vv]ersion\s+([0-9][0-9.]*)", text)
    if not m:
        raise RuntimeError(f"--version answered no version (exit {proc.returncode})")
    return m.group(1).rstrip(".")


def _xtb_water(path: str) -> float:
    with tempfile.TemporaryDirectory(prefix="delfin-xtb-smoke-") as work:
        xyz = Path(work) / "water.xyz"
        xyz.write_text("3\nwater\n" + "".join(
            f"{a} {x} {y} {z}\n" for a, x, y, z in WATER))
        proc = subprocess.run([path, str(xyz), "--sp", "--gfn", "2"], cwd=work,
                              capture_output=True, text=True, check=False,
                              timeout=_BINARY_TIMEOUT_S, env=_child_env())
    m = re.search(r"TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh", proc.stdout or "")
    if proc.returncode != 0 or not m:
        tail = [ln for ln in (proc.stderr or proc.stdout or "").splitlines() if ln.strip()]
        raise RuntimeError("water single point: " + (tail[-1].strip() if tail else
                                                     f"exit {proc.returncode}"))
    energy = float(m.group(1))
    if abs(energy - _GFN2_EH) > 1e-3:
        raise RuntimeError(f"water GFN2 {energy:.6f} Eh, expected {_GFN2_EH:.6f}")
    return energy


def _smoke_binaries(group: Group) -> tuple[str, str]:
    parts, failed, lacking = [], [], []
    for name in group.binaries:
        path = _find_binary(name)
        if not path:
            lacking.append(name)
            continue
        try:
            version = _version_line(path)
            if name == "xtb":
                energy = _xtb_water(path)
                parts.append(f"xtb {version} water GFN2 {energy:.6f} Eh")
            else:
                parts.append(f"{name} {version}")
        except (RuntimeError, subprocess.TimeoutExpired, OSError) as exc:
            failed.append(f"{name}: {exc}")
    if failed:
        return FAILED, "; ".join(failed + parts)
    if not parts:
        return MISSING, "not installed"
    detail = ", ".join(parts)
    if lacking:
        return PARTIAL, f"{detail}; missing: {' '.join(lacking)}"
    return OK, detail


def smoke(group: Group) -> tuple[str, str]:
    """(status, detail) of one tiny computation; never raises."""
    if group.binaries:
        return _smoke_binaries(group)
    status, detail = availability(group)
    if status == MISSING:
        return MISSING, detail
    if status == PARTIAL:
        return PARTIAL, detail
    passed, detail = _run_python(group._code)
    return (OK if passed else FAILED), detail


def off_text(group: Group) -> str:
    """How the Pi switches ``group`` off."""
    if group.switch == "DELFIN_PI_SKIP_EXTRAS":
        return "orca in DELFIN_PI_SKIP_EXTRAS"
    return f"{group.switch}=0"


def pi_hint(group: Group) -> str:
    """How the Pi's deploy installs ``group``."""
    if group.switch == "DELFIN_PI_SKIP_EXTRAS":
        return "installed unless orca is in DELFIN_PI_SKIP_EXTRAS"
    return f"{group.switch}=1 (the default)"


def format_smoke_line(group: str, status: str, detail: str, seconds: float | None = None) -> str:
    """``smoke <group> <status> <detail>``: the line pi_install.sh prints."""
    took = f" ({seconds:.1f} s)" if seconds is not None else ""
    return f"smoke {group:<7} {status:<8} {detail}{took}"


_LINE = re.compile(r"^smoke (\S+)\s+(" + "|".join(STATUSES) + r")\s+(.*)$")


def parse_smoke(text: str) -> dict[str, tuple[str, str]]:
    """{group: (status, detail)} from the lines :func:`format_smoke_line`
    writes; any other line is ignored."""
    out: dict[str, tuple[str, str]] = {}
    for line in text.splitlines():
        m = _LINE.match(line.rstrip())
        if m:
            out[m.group(1)] = (m.group(2), m.group(3))
    return out


def _run_smoke(names: list[str], off: set[str]) -> int:
    counts = dict.fromkeys(STATUSES, 0)
    for name in names:
        group = _BY_NAME[name]
        if name in off:
            status, detail, took = OFF, f"switched off ({off_text(group)})", None
        else:
            start = time.monotonic()
            try:
                status, detail = smoke(group)
            except Exception as exc:  # the report must complete
                status, detail = FAILED, f"check crashed: {type(exc).__name__}: {exc}"
            took = time.monotonic() - start
        counts[status] += 1
        print(format_smoke_line(name, status, detail, took), flush=True)
    print("smoke summary: " + ", ".join(f"{counts[s]} {s}" for s in STATUSES if counts[s]),
          flush=True)
    return 0


def main(argv: list[str] | None = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(prog="python -m delfin.agent.optional_tools",
                                     description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    sm = sub.add_parser("smoke", help="one tiny computation per group")
    sm.add_argument("groups", nargs="*", help=f"default: all ({' '.join(_BY_NAME)})")
    sm.add_argument("--off", default="", help="comma/space separated groups reported as off")
    sub.add_parser("list", help="present or not, nothing is run")
    args = parser.parse_args(argv)
    if args.cmd == "list":
        for group in GROUPS:
            status, detail = availability(group)
            print(format_smoke_line(group.name, status, detail))
        return 0
    names = args.groups or list(_BY_NAME)
    unknown = [n for n in names if n not in _BY_NAME]
    if unknown:
        parser.error(f"unknown group(s): {' '.join(unknown)}; known: {' '.join(_BY_NAME)}")
    off = {w for w in re.split(r"[,\s]+", args.off) if w}
    return _run_smoke(names, off)


if __name__ == "__main__":
    sys.exit(main())
