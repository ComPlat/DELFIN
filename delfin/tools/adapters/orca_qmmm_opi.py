"""DELFIN's ORCA QM/MM building block on OPI: an additive QM/MM single point
for a solute in explicit solvent or a cluster.

This module builds an additive QM/MM ORCA calculation entirely with the ORCA
Python Interface (OPI, orca-pi 2.0.0): the ``!QMMM`` simple keyword and the
``%qmmm`` block come from OPI's input classes (``opi.input``), the geometry
and the QM-region charge/multiplicity from ``opi.input.structures``, and the
``.inp`` file is written by OPI's ``Calculator.write_input()``.  The binary is
launched by DELFIN's existing runner (``delfin.orca.run_orca``) and every
keyword / block option is checked against DELFIN's ORCA knowledge store
(``delfin.orca_knowledge``) before the input is written.

Method family
-------------
The block serves ORCA's additive QM/MM method for a single QM region plus a
force-field-described MM region (``!QMMM``).  ORCA manual 6.1,
multiscalesimulations/qmmm-molecules.txt section 6.2.1 documents the basic
usage::

    ! QMMM

    %qmmm
      QMAtoms { 0 1 2 27 28 } end
      ORCAFFFilename "UBQ.ORCAFF.prms"
    end

    * xyzfile 0 1 UBQ.xyz

The QM region is selected by atom indices (``QMAtoms``), the force field is
given **by reference** to an ORCA force-field parameter file
(``ORCAFFFilename``), and the charge / multiplicity of the QM part is the
coordinate-section charge / multiplicity (``* xyz <charge> <mult>``), which
for a single QM region defines the high-level region (manual
multiscalesimulations/qmmm-general.txt, section 6.2.4 "System Charges for
Stacking or Supermolecular QM:QM2 =  MM Calculations").  ORCA 6.1 places link
atoms automatically at QM-MM boundary bonds (qmmm-general.txt, section 6.2.2
"Boundary Region"); the force-field file supplies the boundary / capping
descriptors without which ORCA cannot build the QM-MM boundary.

Refusals (each returns an error string from :func:`validate_qmmm`)
---------------------------------------------------------------
- **No force-field file**: ``orcaff_filename`` is empty (/chemistry would
  silently drop the MM region), or the referenced file does not exist on
  disk.  ORCA requires ``ORCAFFFilename`` for QM/MM (manual section 6.2.1
  note: "ORCAFF prms files need to be provided for QM/MM ... but not for
  QM/QM2"); a missing file is a prerequisite error, not something we repair.
- **QM region cuts through the boundary without link-atom handling / leaves
  no MM region**: a non-empty ``qmatoms`` selection that covers the whole
  system leaves no MM region for the force field -- the boundary (and the
  link-atom / capping descriptors it needs) cannot be built, which silently
  turns the run into a pure-QM calculation.  We also refuse an empty
  ``qmatoms`` (no QM region) and an index out of the given atom range (a
  boundary cut at an atom that does not exist).  These are the operational,
  testable forms of "a QM region that cuts through a bond without link-atom
  handling"; we deliberately do not force method / charge / multiplicity.
- **Physically wrong charge / multiplicity**: non-integer charge, or a
  multiplicity that is not a positive integer.

Reading
-------
Results come back as typed values with units on :class:`QmmmResult`.  The
**total** energy is read through ``delfin.orca_read.read(...,
Quantity.FINAL_SINGLE_POINT_ENERGY)`` (an additive QM/MM run prints one
``FINAL SINGLE POINT ENERGY`` for the whole system).  The **QM-region**
energy is the converged SCF electronic energy of the high-level subsystem,
printed as ``Total Energy : <float> Eh`` in the converged SCF block -- the
same line the existing reader matches for ``final_energy`` (manual
singlepoint section: the SCF cycle reports the high-level electronic energy
in the electrostatic field of the MM point charges).  The converged geometry
of an additive QM/MM single point is the job geometry (the input
elements + coordinates), returned on the result.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from delfin.orca_read import Quantity, _scf_converged_text, _terminated_normally_text, read

# ---------------------------------------------------------------------------
# defaults (each chosen default is documented against the ORCA manual section
# cited in the module docstring and in the validator's inline comments)
# ---------------------------------------------------------------------------
DEFAULT_EMBEDDING = "electrostatic"
VALID_EMBEDDINGS = ("electrostatic", "mechanical")
VALID_COUPLINGS = ("additive", "subtractive")

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.orca import run_orca

_TOTAL_ENERGY_PATTERN = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
_QM_SCF_PATTERN = re.compile(r"Total Energy\s*:\s*(-?\d+\.\d+)")


@dataclass
class QmmmResult:
    """Typed outcome of an additive QM/MM single point.

    ``meaningful`` is True only when the run terminated normally and the SCF
    converged.  ``total_energy_Eh`` / ``qm_energy_Eh`` are floats with unit
    Hartree when read, ``None`` otherwise.  ``geometry`` is a list of
    ``(element, (x, y, z))`` tuples -- the (converged) geometry of a QM/MM
    single point is the job geometry.
    """

    ok: bool = False
    work_dir: Optional[str] = None
    output_file: Optional[str] = None
    terminated_normally: bool = False
    scf_converged: bool = False
    meaningful: bool = False
    total_energy_Eh: Optional[float] = None
    qm_energy_Eh: Optional[float] = None
    geometry: List[Tuple[str, Tuple[float, float, float]]] = field(
        default_factory=list
    )
    findings: List[str] = field(default_factory=list)
    failure_lines: List[str] = field(default_factory=list)
    error: Optional[str] = None


def validate_qmmm(
    *,
    qmatoms: Sequence[int],
    orcaff_filename: str,
    n_atoms: int,
    charge: int,
    mult: int,
    embedding: str = DEFAULT_EMBEDDING,
    coupling: str = "additive",
    extra_keywords: Sequence[str] = (),
    qmmm_options: Optional[Dict[str, object]] = None,
) -> Tuple[List[str], List[str]]:
    """Validate an additive QM/MM input; returns ``(errors, warnings)``.

    Never "fixes" method, charge, multiplicity or state -- an error always
    means the user must change the input.  Individual refuse conditions are
    documented in the module docstring.
    """
    errors: List[str] = []
    warnings: List[str] = []

    # --- QM region by atom indices -----------------------------------------
    if not qmatoms:
        errors.append(
            "no QM region: QMAtoms 'qmatoms' is empty -- an additive QM/MM "
            "run needs at least one QM atom (manual 6.2.1)."
        )
    covered = set(qmatoms)
    if qmatoms and any(i < 0 or i >= n_atoms for i in qmatoms):
        errors.append(
            f"QM region cuts through the system: a 'qmatoms' index is outside "
            f"[0, {n_atoms - 1}] -- no atom exists there for the boundary "
            f"link-atom handling (manual 6.2.2)."
        )
    elif qmatoms and len(covered) >= n_atoms:
        errors.append(
            "QM region covers the whole system: no MM region remains for the "
            "force field, so the QM-MM boundary (and its link-atom handling) "
            "cannot be built -- this silently turns the run into a pure-QM "
            "calculation (manual 6.2.1)."
        )

    # --- force-field file by reference --------------------------------------
    if not orcaff_filename:
        errors.append(
            "missing force-field file: additive QM/MM requires "
            "ORCAFFFilename by reference (manual 6.2.1)."
        )
    else:
        ff = Path(orcaff_filename)
        if not ff.is_file():
            errors.append(
                f"force-field file '{orcaff_filename}' not found on disk -- "
                "an additive QM/MM run cannot build its MM / boundary "
                "descriptors without it (manual 6.2.1)."
            )

    # --- physically wrong charge / multiplicity ------------------------------
    if not isinstance(charge, int) or isinstance(charge, bool):
        errors.append(f"charge must be an integer, got {charge!r}.")
    if not isinstance(mult, int) or isinstance(mult, bool) or mult < 1:
        errors.append(f"multiplicity must be a positive integer, got {mult!r}.")

    # --- block options / keywords must live in the knowledge store ----------
    from delfin.orca_knowledge import known_block_option, known_keyword

    if embedding not in VALID_EMBEDDINGS:
        errors.append(
            f"unknown embedding {embedding!r}; expected one of "
            f"{VALID_EMBEDDINGS} (manual 6.2.1 EmbeddingType)."
        )
    if coupling not in VALID_COUPLINGS:
        errors.append(
            f"unknown coupling {coupling!r}; expected one of {VALID_COUPLINGS}."
        )
    if not known_block_option("qmmm", "qmatoms"):
        errors.append("the 'qmmm' block does not know the 'qmatoms' option.")
    if not known_block_option("qmmm", "orcafffilename"):
        errors.append("the 'qmmm' block does not know 'orcafffilename'.")
    if not known_block_option("qmmm", "embedding"):
        errors.append("the 'qmmm' block does not know 'embedding'.")
    if not known_block_option("qmmm", "coupling"):
        errors.append("the 'qmmm' block does not know 'coupling'.")

    for kw in extra_keywords:
        if kw and not known_keyword(kw):
            errors.append(f"unknown ORCA simple keyword {kw!r}.")
    qmmm_options = qmmm_options or {}
    for opt in qmmm_options:
        if not known_block_option("qmmm", opt):
            errors.append(f"the 'qmmm' block does not know option {opt!r}.")

    return errors, warnings


def _total_energy_text(out_path: Path) -> Optional[float]:
    """Last ``FINAL SINGLE POINT ENERGY`` (the additive QM/MM total)."""
    try:
        text = out_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    match = [m for m in _TOTAL_ENERGY_PATTERN.finditer(text)]
    if not match:
        return None
    return float(match[-1].group(1))


def _qm_energy_text(out_path: Path) -> Optional[float]:
    """Last converged SCF ``Total Energy :`` (high-level QM electronic)."""
    try:
        text = out_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    match = list(_QM_SCF_PATTERN.finditer(text))
    if not match:
        return None
    return float(match[-1].group(1))


def read_qmmm_output(work_dir: str, base_name: str) -> QmmmResult:
    """Read total + QM energies from a finished additive QM/MM run.

    The total energy is read through ``delfin.orca_read.read`` with
    ``Quantity.FINAL_SINGLE_POINT_ENERGY`` so it obeys the reader's OPI /
    convergence gates; the QM-region energy is parsed here from the same
    converged SCF line the reader matches.  Units: Hartree.
    """
    result = QmmmResult(work_dir=work_dir, output_file=f"{base_name}.out")
    out_path = Path(work_dir) / f"{base_name}.out"
    if not out_path.is_file():
        result.error = f"output file {out_path.name} not found"
        result.findings.append(result.error)
        return result

    total = _total_energy_text(out_path)
    text = out_path.read_text(encoding="utf-8", errors="replace")
    scf_ok = _scf_converged_text(text)
    normal = _terminated_normally_text(out_path)
    result.terminated_normally = normal
    result.scf_converged = scf_ok
    result.total_energy_Eh = total
    result.qm_energy_Eh = _qm_energy_text(out_path)
    result.meaningful = bool(normal and scf_ok and total is not None)
    result.ok = bool(normal and total is not None)
    if not result.ok:
        result.error = (
            "run did not terminate normally" if not normal
            else "no FINAL SINGLE POINT ENERGY found"
        )
        result.failure_lines.extend(
            line.strip() for line in text.splitlines()[-8:]
        )

    # Gate consistency with the reader's view of the same file.
    try:
        reading = read(Path(work_dir), base_name, Quantity.FINAL_SINGLE_POINT_ENERGY)
        if reading is not None and reading.value is not None:
            total_rd = float(reading.value)
            if result.total_energy_Eh is not None and abs(
                result.total_energy_Eh - total_rd
            ) > 1e-6:
                result.findings.append(
                    "total energy differs from the reader's "
                    "FINAL_SINGLE_POINT_ENERGY view; reported the reader value"
                )
                result.total_energy_Eh = total_rd
    except Exception as exc:  # reader is best-effort here
        result.findings.append(f"reader check skipped: {exc}")

    return result


def build_qmmm_input(
    work_dir: str,
    *,
    base_name: str,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    qmatoms: Sequence[int],
    orcaff_filename: str,
    embedding: str = DEFAULT_EMBEDDING,
    coupling: str = "additive",
    qmmm_options: Optional[Dict[str, object]] = None,
    ncores: int = 1,
) -> str:
    """Build an additive QM/MM input with OPI's typed objects.

    Returns the absolute path of the written ``.inp``.  OPI is required (the
    caller gates on it).  All emitted options are validated by
    :func:`validate_qmmm` first (the caller); this function performs the
    mechanical write.
    """
    errors, _ = validate_qmmm(
        qmatoms=qmatoms,
        orcaff_filename=orcaff_filename,
        n_atoms=len(elements),
        charge=charge,
        mult=mult,
        embedding=embedding,
        coupling=coupling,
        qmmm_options=qmmm_options,
    )
    if errors:
        raise ValueError("; ".join(errors))

    from opi.core import Calculator
    from opi.input import Input
    from opi.input.blocks import BlockQmmm
    from opi.input.blocks.util import InputFilePath
    from opi.input.structures import Atom, Structure

    inp = Input()
    inp.add_simple_keywords("QMMM")
    inp.ncores = ncores

    block_args: Dict[str, object] = {
        "qmatoms": list(qmatoms),
        "orcafffilename": InputFilePath.from_string(orcaff_filename),
        "embedding": embedding,
        "coupling": coupling,
    }
    for opt, val in (qmmm_options or {}).items():
        block_args[opt] = val
    inp.add_blocks(BlockQmmm(**block_args))

    atoms = [
        Atom(str(el), coordinates=tuple(map(float, c)))
        for el, c in zip(elements, coordinates)
    ]
    structure = Structure(atoms, charge=charge, multiplicity=mult)

    work_dir_path = Path(work_dir)
    work_dir_path.mkdir(parents=True, exist_ok=True)
    calc = Calculator(base_name, working_dir=str(work_dir_path), version_check=False)
    calc.input = inp
    calc.structure = structure
    calc.write_input()
    return str(work_dir_path / f"{base_name}.inp")


class OrcaQmmmOpiAdapter(StepAdapter):
    """OPI-backed additive QM/MM single-point step (``orca_qmmm_opi``)."""

    name = "orca_qmmm_opi"

    def __init__(self) -> None:
        pass

    def method(self, properties: dict) -> str:
        from delfin.orca_knowledge import known_keyword

        errors, _ = validate_qmmm(
            qmatoms=properties.get("qmatoms", []),
            orcaff_filename=properties.get("orcaff_filename", ""),
            n_atoms=_n_atoms(properties),
            charge=_as_int(properties.get("charge", 0)),
            mult=_as_int(properties.get("mult", 1)),
            embedding=properties.get("embedding", DEFAULT_EMBEDDING),
            coupling=properties.get("coupling", "additive"),
            extra_keywords=properties.get("extra_keywords", ()),
            qmmm_options=properties.get("qmmm_options"),
        )
        if errors:
            raise ValueError("; ".join(errors))
        for kw in properties.get("extra_keywords", ()):
            if kw and not known_keyword(kw):
                raise ValueError(f"unknown ORCA simple keyword {kw!r}.")
        b = "QMMM"
        if properties.get("coupling", "additive") == "subtractive":
            b = f"QM/{properties.get('method', 'XTB')}/MM"
        return b

    def execute(self, run_spec: dict) -> QmmmResult:
        return run_qmmm(**run_spec)


register(OrcaQmmmOpiAdapter())


def _n_atoms(properties: dict) -> int:
    coords = properties.get("coordinates")
    if coords:
        return len(coords)
    elements = properties.get("elements")
    return len(elements) if elements else 0


def _as_int(value: object) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return -1


def run_qmmm(
    *,
    working_dir: str,
    base_name: str,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    mult: int,
    qmatoms: Sequence[int],
    orcaff_filename: str,
    embedding: str = DEFAULT_EMBEDDING,
    coupling: str = "additive",
    qmmm_options: Optional[Dict[str, object]] = None,
    ncores: int = 1,
    timeout: Optional[int] = None,
    properties: Optional[dict] = None,
) -> QmmmResult:
    """End-to-end additive QM/MM single point (validate / build / run / read).

    ``properties``, when given, is the control-extracted job properties dict
    the adapter passes in; it takes precedence over the keyword arguments for
    ``elements``/``coordinates``/``charge``/``mult``.
    """
    props = dict(properties or {})
    if "elements" in props:
        elements = props["elements"]
    if "coordinates" in props:
        coordinates = props["coordinates"]
    if "charge" in props:
        charge = _as_int(props["charge"])
    if "mult" in props:
        mult = _as_int(props["mult"])
    if "qmatoms" in props:
        qmatoms = props["qmatoms"]
    if "orcaff_filename" in props:
        orcaff_filename = props["orcaff_filename"]

    errors, _ = validate_qmmm(
        qmatoms=qmatoms,
        orcaff_filename=orcaff_filename,
        n_atoms=len(elements),
        charge=charge,
        mult=mult,
        embedding=embedding,
        coupling=coupling,
        extra_keywords=props.get("extra_keywords", ()),
        qmmm_options=props.get("qmmm_options") or qmmm_options,
    )
    if errors:
        return QmmmResult(ok=False, error="; ".join(errors), output_file=f"{base_name}.out")

    inp_path = build_qmmm_input(
        working_dir,
        base_name=base_name,
        elements=elements,
        coordinates=coordinates,
        charge=charge,
        mult=mult,
        qmatoms=qmatoms,
        orcaff_filename=orcaff_filename,
        embedding=embedding,
        coupling=coupling,
        qmmm_options=props.get("qmmm_options") or qmmm_options,
        ncores=ncores,
    )
    out_path = inp_path.with_suffix(".out")
    if not run_orca(str(inp_path), str(out_path), timeout=timeout,
                    working_dir=Path(working_dir)):
        result = QmmmResult(ok=False, work_dir=working_dir,
                            output_file=out_path.name,
                            error="ORCA run failed")
        result.findings.append("ORCA terminated abnormally (run_orca returned False)")
        return result
    result = read_qmmm_output(working_dir, base_name)
    result.geometry = [
        (str(el), (float(a), float(b_), float(c)))
        for el, (a, b_, c) in zip(elements, coordinates)
    ]
    return result
