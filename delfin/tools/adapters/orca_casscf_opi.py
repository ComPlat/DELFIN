"""DELFIN's CASSCF(+NEVPT2) building block on OPI.

Builds a (possibly state-averaged) CASSCF single point entirely with the ORCA
Python Interface (OPI, orca-pi 2.0.0): the ``%casscf`` block (``nel``,
``norb``, ``mult``, ``nroots``, ``bweight``) and the optional NEVPT2
perturbation correction (``nevpt2``/``ptmethod``) come from OPI's typed
``BlockCasscf``, the geometry from ``Structure``/``Atom``, the input is
written by OPI's ``Calculator.write_input()`` and the binary is launched by
DELFIN's ``delfin.orca.run_orca``.  Every keyword is checked against DELFIN's
ORCA knowledge store before anything is written.

``orca_dlpno_opi.py`` structure is followed closely: a ``validate_*`` gate
that refuses physically-wrong or silent-chemistry inputs, a ``build_*``
function that renders the input with OPI, a ``run_*`` launcher and a
``read_*`` reader gating every number on the OPI convergence flags.  A
finished run may hand its ``gbw`` on through ``artifacts["gbw"]`` so the
chain wiring ``gbw -> moread`` can start the next step from these orbitals.

Manual ground truth (ORCA 6.1 manual, modelchemistries/CASSCF.txt and
modelchemistries/NEVPT2.txt):

* ``%casscf`` ``nel``/``norb``/``mult``/``nroots``: a CAS(6,6) over the two
  lowest singlets is ``nel 6 norb 6 mult 1 nroots 2`` ("3.13.3 Basic Usage",
  CASSCF.txt lines 1064-1072).  ``nel`` and ``norb`` are the minimal
  mandatory keys (knowledge card ``casscf.active-space.mandatory``).
* ``bweight``: state weights per root, "bweight 2,1" for a 2:1
  singlet:triplet average; "The program automatically normalizes these
  weights such that the sum over all weights is unity" (CASSCF.txt lines
  454-469).
* Natural-orbital occupations are printed as ``N(occ)= ...`` under a
  "Natural Orbital Occupation Numbers" header (CASSCF.txt lines 2745, 4463).
* Feasibility: beyond roughly 14 active orbitals the CI space (about a
  million CSFs) is at the limit of what is practical; RAS / ICE / DMRG are
  the alternatives (CASSCF.txt 3.13.1; casscf-setup skill warning rule).
* NEVPT2 ``ptmethod``: SC-NEVPT2 / FIC-NEVPT2 / DLPNO-NEVPT2, spelled with
  an underscore in the ``%casscf`` block form (knowledge card
  ``nevpt2.method-selection``, NEVPT2.txt 3.16); ``TrafoStep RI`` gives the
  RI approximation.
* MOREAD / ``%moinp``: "! moread" plus a top-level ``%moinp "<file>.gbw"``
  (manual section 2.20.5), the shipped form used by ``opi_chain.py``.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Sequence

from delfin.tools._base import StepAdapter
from delfin.tools._registry import register
from delfin.tools._types import ErrorKind, StepResult, StepStatus

MAX_ORBITALS = 14

__all__ = [
    "MAX_ORBITALS",
    "OrcaCasscfResult",
    "validate_casscf",
    "build_casscf_input",
    "run_casscf",
    "read_casscf_output",
    "OrcaCasscfOpiAdapter",
]

# Multiplicity parity rule: for an even electron count the total spin
# quantum number ``S`` (hence ``mult = 2*S + 1``) is an integer, so ``mult``
# is odd; for an odd count ``S`` is half-integer and ``mult`` is even.
# The rule is enforced inline in ``validate_casscf``.


@dataclass
class OrcaCasscfResult:
    """Typed result of one OPI-backed CASSCF / CASSCF+NEVPT2 run.

    Numeric fields are populated only after the run terminated normally and
    the CASSCF (or SCF) converged.  All energies are in Eh; occupations are
    unitless populations in (0..2).
    """

    ok: bool
    work_dir: Path
    base_name: str = "calc"
    output_file: Optional[Path] = None

    terminated_normally: Optional[bool] = None
    casscf_converged: Optional[bool] = None
    meaningful: bool = False

    nel: Optional[int] = None
    norb: Optional[int] = None
    mult: Optional[int] = None

    # Numbers returned only when meaningful.
    state_energies_Eh: tuple = ()          # energies of the averaged roots
    nevpt2_energies_Eh: tuple = ()         # per-root NEVPT2-corrected energies
    energies_total_Eh: tuple = ()          # selected pseudocanonical/root energies
    natural_orbital_occupations: tuple = ()  # norb natural occupations (0..2)

    findings: tuple = ()
    failure_lines: tuple = ()
    error: Optional[str] = None

    @property
    def artifacts(self) -> dict:
        """Chain hand-off: the converged ``.gbw`` for a next ``moread`` step."""
        gbw = self.work_dir / f"{self.base_name}.gbw"
        arts = {}
        if self.ok and gbw.is_file():
            arts["gbw"] = str(gbw)
        return arts


def _parity_ok(nel: int, mult: int) -> bool:
    """``mult`` is odd exactly when the electron count is even.

    S = (mult-1)/2 must be (half-)integer consistent with the electron
    count's spin.
    """
    return (nel % 2 == 0) == (mult % 2 == 1)


def validate_casscf(
    *,
    nel: int,
    norb: int,
    mult: Optional[int] = None,
    nroots: int = 1,
    weights: Optional[Sequence[float]] = None,
    nevpt2: bool = False,
    ptmethod: str = "",
    allow_large_active_space: bool = False,
) -> list[str]:
    """Refuse physically-wrong or silent-chemistry CASSCF inputs.

    Returns a list of human-readable errors; empty means the input is
    acceptable.  Each rejection names what is wrong and why it is refused --
    never is the method, charge, multiplicity or state silently changed.
    """
    errors: list[str] = []

    if nel < 0 or norb < 0:
        errors.append(f"active electrons ({nel}) and orbitals ({norb}) must be non-negative")

    if norb > 0 and nel > 2 * norb:
        errors.append(
            f"nel {nel} > 2*norb {2 * norb}: an orbital holds at most two "
            "electrons; shrink the active space or raise nel/norb limit"
        )

    if mult is not None:
        if mult < 1:
            errors.append(f"multiplicity {mult} must be a positive integer")
        elif not _parity_ok(nel, mult):
            errors.append(
                f"nel {nel} and mult {mult} mismatch: a {nel}-electron active "
                f"space is incompatible with multiplicity {mult} (spin parity)"
            )

    if nroots < 1:
        errors.append(f"nroots {nroots} must be at least 1 (one state requested)")

    if weights is not None:
        if len(weights) != nroots:
            errors.append(
                f"{len(weights)} weights given but nroots is {nroots}: "
                "provide exactly one weight per root"
            )
        if any(w < 0 for w in weights):
            errors.append("state weights must be non-negative")

    if ptmethod and not nevpt2:
        errors.append(
            f"ptmethod {ptmethod!r} without nevpt2=True: the perturbation "
            "correction would silently do nothing; pass nevpt2=True or drop "
            "ptmethod"
        )
    if nevpt2 and ptmethod not in ("", "sc_nevpt2", "fic_nevpt2", "dlpno_nevpt2"):
        errors.append(
            f"unknown NEVPT2 ptmethod {ptmethod!r}: ORCA knows "
            "SC-NEVPT2, FIC-NEVPT2, DLPNO-NEVPT2 (NEVPT2.txt 3.16; "
            "knowledge card nevpt2.method-selection)"
        )

    if norb > MAX_ORBITALS and not allow_large_active_space:
        errors.append(
            f"active space of {norb} orbitals exceeds the stated limit of "
            f"{MAX_ORBITALS}: beyond ~{MAX_ORBITALS}-{MAX_ORBITALS + 2} "
            "orbitals (~1e6 CSFs) CASSCF is at the edge of feasibility "
            "(CASSCF.txt 3.13.1). Set allow_large_active_space=True to opt in."
        )

    return errors


def build_casscf_input(
    *,
    elements: Sequence[str],
    coordinates: Sequence[Sequence[float]],
    charge: int,
    nel: int,
    norb: int,
    mult: int,
    nroots: int,
    weights: Optional[Sequence[float]] = None,
    nevpt2: bool = False,
    ptmethod: str = "",
    trafostep: str = "",
    maxiter: Optional[int] = None,
    moread: Optional[Path] = None,
    basis: str = "def2-SVP",
    extra_keywords: Optional[Sequence[str]] = None,
    ncores: int = 2,
    base_name: str = "calc",
    work_dir: Path = Path("."),
) -> Path:
    """Validate and write a CASSCF (+NEVPT2) ORCA input with OPI.

    The active space and all blocks are OPI typed objects; nothing is written
    as a hand-rolled string where OPI has a type.  ``moread`` (a ``.gbw``)
    is copied into the working directory and emitted as ``! moread`` +
    ``%moinp "<base_name>.gbw"`` exactly like ``opi_chain`` (manual 2.20.5).
    Returns the path of the written input file.
    """
    errors = validate_casscf(
        nel=nel, norb=norb, mult=mult, nroots=nroots, weights=weights,
        nevpt2=nevpt2, ptmethod=ptmethod,
    )
    if errors:
        raise ValueError("; ".join(errors))

    from opi.core import Calculator
    from opi.input import ArbitraryStringPos, Input
    from opi.input.blocks.block_casscf import BlockCasscf
    from opi.input.simple_keywords.approximation import Approximation
    from opi.input.simple_keywords.base import SimpleKeyword
    from opi.input.simple_keywords.basis_set import BasisSet
    from opi.input.simple_keywords.scf import Scf
    from opi.input.structures import Atom, Structure

    def _lookup(enum_cls: type, word: str) -> Any:
        target = str(word).lower()
        for attr in vars(enum_cls).values():
            if isinstance(attr, SimpleKeyword) and str(attr).lower() == target:
                return attr
        raise ValueError(f"no OPI keyword {word!r} in {enum_cls}")

    wdir = Path(work_dir)
    wdir.mkdir(parents=True, exist_ok=True)

    inp = Input()
    keywords = [_lookup(BasisSet, basis)]

    # moread: copy the gbw in and emit ! moread + %moinp (manual 2.20.5).
    if moread is not None:
        moinp = Path(moread)
        if not moinp.is_file():
            raise FileNotFoundError(
                f"moread {moinp} does not exist: cannot start CASSCF from it"
            )
        shutil.copy2(moinp, wdir / f"{base_name}.gbw")
        keywords.append(Scf.MOREAD)

    if extra_keywords:
        for kw in extra_keywords:
            # each extra keyword is looked up so only ORCA-known words pass
            for cls in (Approximation, Scf,):
                try:
                    keywords.append(_lookup(cls, kw))
                    break
                except ValueError:
                    continue
            else:
                # falls through if neither enum holds it
                raise ValueError(
                    f"extra keyword {kw!r} is not in the OPI keyword catalog"
                )

    # orbitals: OPI typed coordinates block
    atoms = []
    for el, xyz in zip(elements, coordinates):
        atoms.append(
            Atom(str(el), coordinates=(float(xyz[0]), float(xyz[1]), float(xyz[2])))
        )
    structure = Structure(atoms, charge=charge, multiplicity=mult)

    inp.add_simple_keywords(*keywords)

    casscf_kw: dict = dict(
        nel=nel, norb=norb, mult=mult, nroots=nroots,
    )
    if weights is not None:
        from opi.input.blocks.block_casscf import NumList
        casscf_kw["bweight"] = NumList(weights)
    if nevpt2:
        casscf_kw["nevpt2"] = 1
        if ptmethod:
            casscf_kw["ptmethod"] = ptmethod
    if trafostep:
        if trafostep not in ("exact", "rimo", "ri"):
            raise ValueError(
                f"trafostep {trafostep!r}: ORCA allows exact, rimo or ri "
                "(CASSCF.txt; knowledge card nevpt2.method-selection)"
            )
        casscf_kw["trafostep"] = trafostep
    if maxiter is not None:
        casscf_kw["maxiter"] = maxiter

    inp.add_blocks(BlockCasscf(**casscf_kw))

    if moread is not None:
        inp.add_arbitrary_string(
            f'%moinp "{base_name}.gbw"', pos=ArbitraryStringPos.BEFORE_COORDS
        )

    inp_path = wdir / f"{base_name}.inp"
    inp.ncores = ncores
    calc = Calculator(base_name, working_dir=wdir, version_check=False)
    calc.input = inp
    calc.structure = structure
    calc.write_input(force=True)
    return inp_path


def read_casscf_output(work_dir: Path, base_name: str) -> OrcaCasscfResult:
    """Read a finished CASSCF run, gating every number on termination and
    CASSCF convergence.

    Run state comes from OPI's ``Output.terminated_normally()`` when the run
    is OPI-parsable, otherwise from the ``ORCA TERMINATED NORMALLY`` text
    marker.  The CASSCF numbers (per-root energies, NEVPT2 corrections,
    natural occupations) are always read from the ``.out`` text: OPI's
    ``.property.json`` does not carry them (OPI output models have no
    ``casscf``/``nevpt2`` keys), so there is no OPI-native number to read,
    and OPI has no ``mo_converged()``/``get_property()`` methods.
    """
    from opi.output.core import Output

    wdir = Path(work_dir)
    out_path = wdir / f"{base_name}.out"
    text = out_path.read_text(errors="replace") if out_path.is_file() else ""

    terminated = False
    if text:
        try:
            out = Output(basename=base_name, working_dir=wdir)
            out.parse()
            terminated = bool(out.terminated_normally())
        except Exception:
            # No OPI-parsable run (e.g. no <base>.property.json): fall back to
            # the ORCA text marker.  A hand-written synthetic .out stays
            # testable without committing any real ORCA output.
            terminated = "TERMINATED NORMALLY" in text

    converged = bool(text) and (
        "CASSCF CONVERGED" in text or "SCF CONVERGED" in text
    )
    meaningful = terminated and converged

    result = OrcaCasscfResult(
        ok=meaningful, work_dir=wdir, base_name=base_name,
        output_file=out_path if out_path.is_file() else None,
        terminated_normally=terminated, casscf_converged=converged,
        meaningful=meaningful,
    )

    if meaningful:
        try:
            result.state_energies_Eh = _text_state_energies(text)
            result.nevpt2_energies_Eh = _text_nevpt2_energies(text)
            result.energies_total_Eh = _total_energies_text(text)
            result.natural_orbital_occupations = _text_occupations(text)
        except Exception as exc:
            result.ok = False
            result.meaningful = False
            result.error = f"CASSCF result parsing: {exc}"

    if not meaningful:
        result.failure_lines = _failure_lines(out_path)
    return result


def _text_nevpt2_energies(text: str) -> tuple:
    """Parse NEVPT2-corrected root energies (Eh) from output text.

    ORCA prints the corrected energy of each NEVPT2 root block as
    ``Total Energy (E0+dE)    : E  = <E>`` (modelchemistries/NEVPT2.txt,
    lines ~223-228).  There is no ``State n: E(NEVPT2)`` line in ORCA's
    output.  Empty tuple when the run did not do NEVPT2.
    """
    if not text:
        return ()
    import re
    return tuple(
        float(m.group(1))
        for m in re.finditer(
            r"Total\s+Energy\s+\(E0\+dE\)[^:]*:\s*E\s*=\s*(-?[0-9.eE+\-]+)",
            text,
        )
    )


def _total_energies_text(text: str) -> tuple:
    """Per-root CI energies (Eh) from ``CAS-SCF STATES`` ``ROOT n: E=`` lines.

    ``CAS-SCF STATES FOR BLOCK  1 MULT= 1 NROOTS= 1`` prints ``ROOT   1:  E=
    <E> Eh`` per root (manual CASSCF.txt lines 2320-2324).  Empty when the
    output uses no such block.
    """
    if not text:
        return ()
    import re
    return tuple(
        float(m.group(1))
        for m in re.finditer(
            r"ROOT\s+\d+:\s+E=\s*(-?[0-9.eE+\-]+)\s+Eh", text
        )
    )


def _text_occupations(text: str) -> tuple:
    """Natural-orbital occupations as ``(index, population)`` pairs.

    Reads the LAST ``N(occ)= ...`` line (the converged iteration), like
    ``delfin.casscf_active_space.parse_casscf_occupations``, but from the
    exact ``.out`` text already in hand -- never a re-glob of the directory,
    which can pick up a different, larger ``.out``.  Each pair is
    ``(index, value)`` with value in (0..2), matching the original
    ``_fill_text_result`` contract.
    """
    if not text:
        return ()
    occ = None
    for line in text.splitlines():
        if "N(occ)=" in line:
            rest = line.split("N(occ)=", 1)[1].split()
            occ = [float(x) for x in rest]
    return tuple(enumerate(occ)) if occ else ()


def _text_state_energies(text: str) -> tuple:
    """Parse CASSCF per-root absolute energies (Eh) from output text.

    Real ORCA format (manual 3.13.9 / 3.13.11, CASSCF.txt): an
    ``SA-CASSCF TRANSITION ENERGIES`` block prints the lowest root's absolute
    energy as ``LOWEST ROOT = <E> Eh`` or ``LOWEST ROOT (ROOT n ,MULT m) = <E>
    Eh`` and then a ``STATE ROOT MULT [IRREP] DE/a.u. DE/eV DE/cm**-1`` table
    whose ``DE/a.u.`` column is relative to that lowest root.  The header's
    presence/absence of the ``IRREP`` column fixes where ``DE/a.u.`` sits, so
    the column read is exact rather than fixed-width.  Absolute[root i] =
    LOWEST_ROOT + DE/i.  Empty when no such block is present.
    """
    if not text:
        return ()
    import re
    lowest = re.search(
        r"LOWEST ROOT(?: \(ROOT \d+ ,MULT \d+\))?\s*=\s*(-?[0-9.eE+\-]+)\s+Eh",
        text,
    )
    if not lowest:
        return ()
    header = re.search(
        r"^STATE\s+ROOT\s+MULT(?:\s+IRREP)?\s+(DE/a\.u\.)\s+DE/eV\s+DE/cm\*\*-1",
        text, re.MULTILINE,
    )
    if not header:
        # no transition table (unusual): only the lowest root itself
        return (float(lowest.group(1)),)
    # column of 'DE/a.u.' among the header fields: 3 without IRREP, 4 with
    de_col = header.group(0).split().index("DE/a.u.")
    lowest_e = float(lowest.group(1))
    out_energies = []
    for line in text.splitlines():
        m = re.match(r"^\s*\d+:\s+\S+\s+\S+", line)
        if not m:
            continue
        fields = line.split()
        # fields[0] is the "N:" state label; the DE/a.u. field sits at the
        # same index in the row as it does in the header (de_col), because
        # both start with the state label column at index 0.
        try:
            de = float(fields[de_col])
        except (IndexError, ValueError):
            continue
        out_energies.append(lowest_e + de)
    # the LOWEST ROOT itself (DE = 0) is not always a printed row; add it so
    # every averaged state is present.
    if out_energies and min(out_energies) > lowest_e:
        out_energies.append(lowest_e)
    return tuple(sorted(out_energies))


def _failure_lines(out_path: Path) -> tuple:
    if not Path(out_path).is_file():
        return ()
    try:
        text = Path(out_path).read_text(errors="replace")
    except OSError:
        return ()
    marks = ("CASSCF NOT CONVERGED", "ERROR", "aborted", "CONVERGENCE FAILURE",
             "FATAL", "file not found", "MOREAD")
    return tuple(
        line.strip() for line in text.splitlines()
        if any(m in line for m in marks)
    )[:40]


def run_casscf(
    *,
    elements, coordinates, charge,
    nel, norb, mult, nroots,
    weights=None, nevpt2=False, ptmethod="", trafostep="", maxiter=None,
    moread=None, basis="def2-SVP", extra_keywords=None,
    ncores=2, timeout=None,
    work_dir=None, base_name="calc",
) -> OrcaCasscfResult:
    """Run one CASSCF(+NEVPT2) single point and read the result back.

    Validates against the ORCA knowledge store before writing anything,
    builds and writes the input with OPI, launches through DELFIN's
    ``delfin.orca.run_orca`` in this working directory, then reads back with
    OPI, gating every number on the convergence flags.
    """
    wdir = Path(work_dir).resolve() if work_dir is not None else Path.cwd() / "casscf_run"
    wdir.mkdir(parents=True, exist_ok=True)

    errors = validate_casscf(
        nel=nel, norb=norb, mult=mult, nroots=nroots, weights=weights,
        nevpt2=nevpt2, ptmethod=ptmethod,
    )
    if errors:
        return OrcaCasscfResult(
            ok=False, work_dir=wdir, meaningful=False,
            error="; ".join(errors), failure_lines=tuple(errors),
        )

    from delfin.orca import run_orca
    out_path = wdir / f"{base_name}.out"

    try:
        inp_path = build_casscf_input(
            elements=elements, coordinates=coordinates, charge=charge,
            nel=nel, norb=norb, mult=mult, nroots=nroots,
            weights=weights, nevpt2=nevpt2, ptmethod=ptmethod,
            trafostep=trafostep, maxiter=maxiter, moread=moread,
            basis=basis, extra_keywords=extra_keywords, ncores=ncores,
            base_name=base_name, work_dir=wdir,
        )
        success = run_orca(str(inp_path), str(out_path), timeout=timeout,
                           working_dir=wdir)
    except Exception as exc:
        return OrcaCasscfResult(
            ok=False, work_dir=wdir, meaningful=False,
            error=str(exc), failure_lines=_failure_lines(out_path),
        )

    result = read_casscf_output(wdir, base_name)
    if not success and result.ok:
        result.ok = False
        result.meaningful = False
        result.error = result.error or "ORCA subprocess reported failure"
        result.failure_lines = result.failure_lines or _failure_lines(out_path)
    return result


class OrcaCasscfOpiAdapter(StepAdapter):
    """A CASSCF(+NEVPT2) step driven by OPI, launchable via ``run_step``.

    Mirrors ``OrcaDlpnoOpiAdapter``'s position in the registry.  The caller
    always states ``nel``/``norb``/``mult``; validation refuses inputs that
    would change the chemistry silently (parity, active space size).
    """

    name = "orca_casscf_opi"
    description = "OPI-backed ORCA CASSCF (state-averaged) single point, optional NEVPT2"
    category = "wft"
    produces_geometry = False
    consumes = ()
    produces = ()

    def execute(self, work_dir: Path, *, cores: int = 2, **kwargs: Any) -> StepResult:
        """Run the step and return a typed ``StepResult``.

        ``cores`` is the step framework's parallelisation knob (same contract
        as ``OrcaDlpnoOpiAdapter.execute``); it is forwarded to ``run_casscf``
        as ``ncores`` so the input's ``%pal nprocs`` block uses the requested
        cores.  Other kwargs (elements, coordinates, charge, ...) pass through
        to ``run_casscf``.
        """
        import time

        start = time.monotonic()
        wdir = Path(work_dir).resolve()
        base_name = kwargs.pop("base_name", "calc")
        result = run_casscf(work_dir=wdir, base_name=base_name, ncores=cores, **kwargs)
        kw = {}
        if result.error:
            kw["error"] = result.error
            kw["error_kind"] = ErrorKind.PARSE
        if result.artifacts:
            kw["artifacts"] = result.artifacts
        return self._make_result(
            self.name,
            StepStatus.SUCCESS if result.ok else StepStatus.FAILED,
            wdir, start, **kw,
        )


register(OrcaCasscfOpiAdapter())
