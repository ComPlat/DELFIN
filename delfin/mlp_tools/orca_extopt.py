"""ORCA <-> MLP bridge for the ORCA external-optimizer interface (ExtOpt).

This module couples DELFIN's MLP backends (delfin.mlp_tools.calculators)
to ORCA's optimizer infrastructure so that ``! ExtOpt`` jobs -- geometry
optimizations, NEB-TS, IRC, GOAT and relaxed surface scans -- can run on
machine-learned potential gradients.

Interface contract, verified against the ORCA 6.1.1 manual, chapter
4.1.10 "Optimizing with External Methods" (ORCA manual text
``structurereactivity/optimizations.txt``, sections 4.1.10, lines
1537-1701):

* The ORCA input activates the interface with ``! ExtOpt`` and names the
  wrapper in a ``%method`` block::

      ! ExtOpt Opt
      %method
        ProgExt "/full/path/to/wrapperscript"
        Ext_Params "optional command line arguments"   # optional
      end

  ``ProgExt`` may also be given as a file named ``otool_external`` next to
  the ORCA executables, or via the ``EXTOPTEXE`` environment variable.
  All electronic-structure information in the input is discarded.

* At every step ORCA writes ``basename_EXT.extinp.tmp``::

      basename_EXT.xyz # xyz filename: string, ending in '.xyz'
      0                # charge: integer
      1                # multiplicity: positive integer
      1                # NCores: positive integer
      0                # do gradient: 0 or 1
      pointcharges.pc  # point charge filename: string (optional)

  and ``basename_EXT.xyz`` in the working directory.  ``#`` comments run
  to end of line.  The energy is always requested; the gradient only when
  needed (``do gradient`` = 1).

* ORCA calls ``scriptname basename_EXT.extinp.tmp [args]`` where
  ``[args]`` is the verbatim ``Ext_Params`` string.

* The wrapper writes ``basename_EXT.engrad`` (same basename as the xyz
  file)::

      #
      # Number of atoms: must match the XYZ
      #
      <n>
      #
      # The current total energy in Eh
      #
      <energy>
      #
      # The current gradient in Eh/bohr: Atom1X, Atom1Y, Atom1Z, Atom2X, ...
      #
      <g1x> <g1y> <g1z> ...

  ``#`` comments and comment-only lines are ignored.

Unit conventions: DELFIN's MLP layer (ASE) reports energies in eV and
forces in eV/Angstrom; ORCA expects Hartree and Hartree/Bohr.  The
conversion factors below are CODATA 2018 (1 Eh = 27.211386245988 eV,
1 Bohr = 0.529177210903 Angstrom), the same values ase.units uses.

The OPI classes ``opi.external_methods.interface.ExtoptInterface``
(read_extopt_input / write_orca_input) implement the same file formats;
the format constants here were checked against that source.

Thread safety on shared hosts: the wrapper honours the machine rule that
MLP runs stay on CPU with at most 2 threads (OMP/MKL/torch), set before
torch is imported (see :func:`main`).
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np

from delfin.mlp_tools.calculators import create_calculator  # noqa: E402

#: Hartree per eV (CODATA 2018: 1 Eh = 27.211386245988 eV).
EV_TO_EH: float = 1.0 / 27.211386245988

#: Hartree/Bohr per eV/Angstrom: 1 Eh/Bohr = 1 Eh / (Bohr) = 27.211386245988 eV
#: / 0.529177210903 A = 51.422067480645 eV/A, so eV/A -> Eh/Bohr DIVIDES by it.
EV_PER_ANG_TO_EH_PER_BOHR: float = 27.211386245988 / 0.529177210903

#: CPU thread cap for MLP evaluation on shared hosts (machine rule:
#: real MLP runs here only on CPU, <= 2 threads, no GPU).
MAX_CPU_THREADS = 2

_EXTINP_SUFFIX = ".extinp.tmp"
_ENGRAD_SUFFIX = ".engrad"


@dataclass(frozen=True)
class ExtinpInfo:
    """Parsed ``basename_EXT.extinp.tmp`` (manual 4.1.10)."""

    xyz_filename: str
    charge: int
    multiplicity: int
    ncores: int
    do_gradient: bool
    pointcharge_filename: Optional[str] = None


def parse_extinp(path: Path | str) -> ExtinpInfo:
    """Parse the input file ORCA writes for the external program.

    Lines: xyz filename, charge, multiplicity, NCores, do_gradient
    (0/1), optional pointcharge filename.  ``#`` starts a comment.
    """
    lines = [
        line.split("#", 1)[0].strip()
        for line in Path(path).read_text().splitlines()
        if line.split("#", 1)[0].strip()
    ]
    if len(lines) < 5:
        raise ValueError(
            f"{path}: external-tool input must have at least 5 fields, "
            f"got {len(lines)} (xyz, charge, multiplicity, NCores, do_gradient)"
        )
    xyz_filename = lines[0]
    charge = int(lines[1])
    multiplicity = int(lines[2])
    ncores = int(lines[3])
    do_gradient = int(lines[4])
    if do_gradient not in (0, 1):
        raise ValueError(f"{path}: do_gradient must be 0 or 1, got {do_gradient}")
    if multiplicity < 1:
        raise ValueError(f"{path}: multiplicity must be >= 1, got {multiplicity}")
    if ncores < 1:
        raise ValueError(f"{path}: NCores must be >= 1, got {ncores}")
    pointcharge_filename = lines[5] if len(lines) >= 6 else None
    return ExtinpInfo(
        xyz_filename=xyz_filename,
        charge=charge,
        multiplicity=multiplicity,
        ncores=ncores,
        do_gradient=bool(do_gradient),
        pointcharge_filename=pointcharge_filename,
    )


def write_engrad(
    path: Path | str,
    *,
    n_atoms: int,
    energy_eh: float,
    gradient_eh_bohr: Optional[np.ndarray] = None,
) -> None:
    """Write an ``.engrad`` file in the format ORCA expects (manual 4.1.10).

    ``energy_eh`` is the total energy in Hartree; ``gradient_eh_bohr`` the
    gradient flattened atom-major (Atom1X, Atom1Y, Atom1Z, Atom2X, ...) in
    Hartree/Bohr, omitted entirely when None (energy-only request).
    """
    out = []
    out.append("#\n")
    out.append("# Number of atoms: must match the XYZ\n")
    out.append("#\n")
    out.append(f"{n_atoms}\n")
    out.append("#\n")
    out.append("# The current total energy in Eh\n")
    out.append("#\n")
    out.append(f"{energy_eh:.12f}\n")
    if gradient_eh_bohr is not None:
        grad = np.asarray(gradient_eh_bohr, dtype=float).reshape(-1)
        if grad.size != 3 * n_atoms:
            raise ValueError(
                f"gradient has {grad.size} components, expected {3 * n_atoms} "
                f"for {n_atoms} atoms"
            )
        out.append("#\n")
        out.append("# The current gradient in Eh/bohr: Atom1X, Atom1Y, "
                   "Atom1Z, Atom2X, etc.\n")
        out.append("#\n")
        out.extend(f"{g:.12f}\n" for g in grad)
    Path(path).write_text("".join(out))


def _limit_cpu_threads() -> None:
    """Cap OMP/MKL/torch threads at MAX_CPU_THREADS before torch loads."""
    for var in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "TORCH_NUM_THREADS",
    ):
        os.environ.setdefault(var, str(MAX_CPU_THREADS))
    try:  # pragma: no cover - torch may be absent in CI
        import torch

        torch.set_num_threads(MAX_CPU_THREADS)
    except ImportError:
        pass


def read_xyz(path: Path | str):
    """Minimal xyz reader: returns (symbols, positions ndarray in Angstrom)."""
    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path}: not an xyz file")
    try:
        n_atoms = int(lines[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: first line is not an atom count") from exc
    symbols = []
    positions = []
    for line in lines[2 : 2 + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"{path}: malformed atom line {line!r}")
        symbols.append(parts[0])
        positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
    if len(symbols) != n_atoms:
        raise ValueError(
            f"{path}: declares {n_atoms} atoms but has {len(symbols)} coordinate lines"
        )
    return symbols, np.array(positions, dtype=float)


def run_wrapper(
    extinp_path: Path | str,
    *,
    backend: str = "ani2x",
    calculator_factory: Optional[Callable] = None,
    device: str = "cpu",
    max_threads: int = MAX_CPU_THREADS,
) -> Path:
    """One external-program call: read extinp + xyz, evaluate the MLP,
    write the ``.engrad`` file ORCA reads back.

    ``calculator_factory`` defaults to
    ``delfin.mlp_tools.calculators.create_calculator``; tests inject a
    fake.  A missing backend package or model must surface as the
    backend's own clear error (available=False + reason + install hint in
    the calculator layer) -- it is never silently swallowed here.
    """
    extinp_path = Path(extinp_path)
    info = parse_extinp(extinp_path)

    xyz_path = extinp_path.parent / info.xyz_filename
    if not xyz_path.exists():
        raise FileNotFoundError(
            f"ORCA wrote {extinp_path.name} but the xyz file "
            f"{xyz_path} is missing"
        )
    symbols, positions = read_xyz(xyz_path)

    if calculator_factory is None:
        calculator_factory = create_calculator

    atoms = _atoms_from_xyz(symbols, positions, info.charge, info.multiplicity)
    calc = calculator_factory(
        backend, device=device, charge=info.charge, mult=info.multiplicity
    )
    atoms.calc = calc
    energy_ev = float(atoms.get_potential_energy())
    energy_eh = energy_ev * EV_TO_EH
    if info.do_gradient:
        forces_ev_ang = np.asarray(atoms.get_forces(), dtype=float)
        # ORCA wants the gradient (= -force) in Eh/bohr, atom-major.
        gradient_eh_bohr = -forces_ev_ang.reshape(-1) * EV_PER_ANG_TO_EH_PER_BOHR
    else:
        gradient_eh_bohr = None

    engrad_path = extinp_path.with_name(
        extinp_path.name.replace(_EXTINP_SUFFIX, _ENGRAD_SUFFIX)
    )
    if _EXTINP_SUFFIX not in extinp_path.name:
        # Fallback: same stem + .engrad (the manual only guarantees the
        # basename association, not the exact tmp suffix).
        engrad_path = extinp_path.with_suffix(_ENGRAD_SUFFIX)
    write_engrad(
        engrad_path,
        n_atoms=len(symbols),
        energy_eh=energy_eh,
        gradient_eh_bohr=gradient_eh_bohr,
    )
    return engrad_path


def _atoms_from_xyz(symbols, positions, charge: int, mult: int):
    """ASE Atoms from the xyz (lazy ASE import so unit tests stay light)."""
    from ase import Atoms
    from ase.units import Bohr  # noqa: F401  (documented unit anchor)

    atoms = Atoms(symbols=symbols, positions=positions)
    atoms.info["charge"] = charge
    atoms.info["multiplicity"] = mult
    return atoms


#: Supported ORCA run types driven through the external-optimizer
#: infrastructure (manual 4.1.10: "Other calculations which exclusively
#: rely on energies and gradients are also supported"; the examples name
#: Opt, NEB-TS and GOAT; IRC and relaxed scans are energy+gradient runs).
EXTOPT_RUN_TYPES: dict[str, tuple[str, ...]] = {
    # run type -> simple-keyword line (first token) written for it
    "opt": ("Opt",),
    "optts": ("OptTS",),
    "neb-ts": ("NEB-TS",),
    "neb": ("NEB",),
    "irc": ("NEB-TS", "IRC"),
    "goat": ("GOAT",),
    "scan": ("Opt",),  # relaxed scan via %geom Scan block
}


def build_extopt_input(
    run_type: str,
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int = 0,
    mult: int = 1,
    wrapper_path: str = "delfin-mlp-extopt",
    backend: str = "ani2x",
    extra_keywords: Sequence[str] = (),
    scan: Optional[dict] = None,
    base_name: str = "extopt",
    work_dir: Path | str = ".",
) -> Path:
    """Write a valid ``! ExtOpt`` ORCA input for a chosen run type.

    Parameters
    ----------
    run_type : one of :data:`EXTOPT_RUN_TYPES` ("opt", "optts", "neb-ts",
        "neb", "irc", "goat", "scan").
    elements/coordinates : the geometry (xyz order, Angstrom).
    charge/mult : molecular charge and spin multiplicity.
    wrapper_path : path ORCA passes as ``ProgExt`` (the wrapper script).
    backend : MLP backend name, passed to the wrapper through
        ``Ext_Params`` as ``--backend <name>`` (manual 4.1.10: Ext_Params
        arguments are forwarded verbatim to the script).
    extra_keywords : additional simple keywords (e.g. constraints).
    scan : for ``run_type="scan"``, a dict for the ``%geom Scan`` block
        with keys ``mode`` ("A"/"B"/"D"), ``atom1``, ``atom2``, optional
        ``atom3``/``atom4``, ``startvalue``, ``endvalue``, ``steps``
        (ORCA manual, surface scans; OPI class
        ``opi.input.blocks.block_geom.Scan``).  Angstrom / degrees, as
        ORCA expects scan values in the input geometry units.

    The OPI classes ``ExternalTools.EXTOPT`` (simple keyword ``extopt``),
    ``BlockMethod`` (fields ``ProgExt``/``Ext_Params``) and ``Scan``
    express this input natively; when OPI is importable the file is
    written through ``opi.core.Calculator``, otherwise a plain-text
    writer produces the same content (identical keywords, verified in
    tests/test_opi_mlp_extopt_builder.py for parity).

    Returns the written ``<base_name>.inp`` path.
    """
    run_type = run_type.lower().strip()
    if run_type not in EXTOPT_RUN_TYPES:
        raise ValueError(
            f"run_type {run_type!r} not one of {sorted(EXTOPT_RUN_TYPES)}"
        )
    if run_type == "scan" and scan is None:
        raise ValueError("run_type='scan' requires the scan= dict")
    keywords = list(EXTOPT_RUN_TYPES[run_type])
    keywords.extend(extra_keywords)

    ext_params = f"--backend {backend} --device cpu"
    xyz_name = f"{base_name}.xyz"

    try:
        return _build_with_opi(
            keywords,
            elements,
            coordinates,
            charge,
            mult,
            wrapper_path,
            ext_params,
            scan,
            base_name,
            work_dir,
        )
    except ImportError:
        return _write_plain_input(
            keywords,
            elements,
            coordinates,
            charge,
            mult,
            wrapper_path,
            ext_params,
            scan,
            xyz_name,
            Path(work_dir) / f"{base_name}.inp",
        )


def _build_with_opi(keywords, elements, coordinates, charge, mult,
                    wrapper_path, ext_params, scan, base_name, work_dir) -> Path:
    """Build the input through OPI (py3.12 venv with the ``opi`` package)."""
    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks.block_method import BlockMethod
    from opi.input.structures import Atom, Structure
    from opi.input.simple_keywords.base import SimpleKeyword

    def _kw(name: str) -> str:
        # OPI stores simple keywords lower-case; pass the canonical form.
        return str(SimpleKeyword(name.lower()))

    inp = Input()
    for kw in keywords:
        inp.add_simple_keywords(_kw(kw))
    inp.add_simple_keywords("extopt")  # ExternalTools.EXTOPT (OPI 6.1.1)
    inp.ncores = 1
    inp.memory = 1000

    method = BlockMethod()
    method.ProgExt = wrapper_path  # type: ignore[assignment]
    method.Ext_Params = ext_params
    inp.add_blocks(method)

    atoms = [
        Atom(str(el), coordinates=(float(x), float(y), float(z)))
        for el, (x, y, z) in zip(elements, coordinates)
    ]
    structure = Structure(atoms, charge=charge, multiplicity=mult)

    calc = Calculator(base_name, working_dir=work_dir, version_check=False)
    calc.input = inp
    calc.structure = structure
    if scan is not None:
        from opi.input.blocks.block_geom import BlockGeom, Scan

        geom = BlockGeom()
        geom.Scan = Scan(
            mode=str(scan["mode"]).upper(),
            atom1=int(scan["atom1"]),
            atom2=int(scan["atom2"]),
            atom3=scan.get("atom3"),
            atom4=scan.get("atom4"),
            startvalue=scan.get("startvalue"),
            endvalue=scan.get("endvalue"),
            steps=scan.get("steps"),
            valuelist=scan.get("valuelist"),
        )
        calc.input.add_blocks(geom)
    calc.write_input(force=True)
    return Path(work_dir) / f"{base_name}.inp"


def _write_plain_input(keywords, elements, coordinates, charge, mult,
                       wrapper_path, ext_params, scan, xyz_name, out_path) -> Path:
    """Plain-text fallback identical to the OPI output (no OPI needed)."""
    lines = [f"! {' '.join(keywords)} ExtOpt", ""]
    lines.append("%method")
    lines.append(f'  ProgExt "{wrapper_path}"')
    lines.append(f'  Ext_Params "{ext_params}"')
    lines.append("end")
    lines.append("")
    if scan is not None:
        lines.append("%geom")
        parts = [f"Scan {scan['mode']}"]
        for key in ("atom1", "atom2", "atom3", "atom4"):
            if scan.get(key) is not None:
                parts.append(str(scan[key]))
        if scan.get("valuelist") is not None:
            parts.append(" ".join(str(v) for v in scan["valuelist"]))
        else:
            for key in ("startvalue", "endvalue", "steps"):
                parts.append(str(scan[key]))
        lines.append(f"  {parts[0]} {' '.join(parts[1:])}")
        lines.append("end")
        lines.append("")
    lines.append(f"*xyzfile {charge} {mult} {xyz_name}")
    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry point; ORCA calls ``script <extinp> [args]``.

    ``--backend`` selects the MLP backend (default ``ani2x``); the rest of
    ``args`` (ORCA's ``Ext_Params``) is accepted and ignored, so arbitrary
    backend options can be appended later without breaking ORCA.
    """
    _limit_cpu_threads()
    parser = argparse.ArgumentParser(
        prog="delfin-mlp-extopt",
        description="ORCA ExtOpt wrapper: MLP energy/gradient for ORCA's "
        "external optimizer (manual 4.1.10).",
    )
    parser.add_argument("extinp", help="basename_EXT.extinp.tmp written by ORCA")
    parser.add_argument(
        "--backend", default="ani2x",
        help="MLP backend name for delfin.mlp_tools.calculators.create_calculator"
        " (default: ani2x)",
    )
    parser.add_argument(
        "--device", default="cpu",
        help="MLP device; this host is CPU-only (default: cpu)",
    )
    args, _unknown = parser.parse_known_args(argv)
    engrad = run_wrapper(args.extinp, backend=args.backend, device=args.device)
    print(f"delfin-mlp-extopt: wrote {engrad}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
