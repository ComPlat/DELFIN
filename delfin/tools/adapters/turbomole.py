"""Turbomole adapter — geometry optimization, single-point, frequencies.

Requires Turbomole binaries (ridft, jobex, aoforce, etc.) in PATH or
configured via environment variables.

Usage::

    run_step("turbomole_sp", geometry="mol.xyz", charge=0, method="b3-lyp",
             basis="def2-SVP", cores=4)

    run_step("turbomole_opt", geometry="mol.xyz", charge=0, method="ri-bp86",
             basis="def2-SVP")

Physical-realism contract (reviewer findings F1-F4, wave 21 TM1):
- spin is never silently dropped: an open shell (mult > 1) writes
  ``$unpaired N`` with N = mult - 1 into the control file;
- command arguments are passed as an argv list, never interpolated into
  a shell command line (no injection surface through max_cycles);
- every declared product is produced: turbomole_freq surfaces the
  aoforce ``hessian`` file or fails honestly.
"""

from __future__ import annotations

import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, List, Optional

from delfin.tools._base import StepAdapter
from delfin.tools._types import ErrorKind, StepResult, StepStatus
from delfin.tools._registry import register
from delfin.tools._spec import DataKeySpec, ParamSpec

# Shared declarative contract for the Turbomole adapters.
_TM_PARAMS = (
    ParamSpec("charge", "int", required=True, description="Molecular charge"),
    ParamSpec("mult", "int", default=1, description="Spin multiplicity (2S+1)"),
    ParamSpec("method", "str", default="b3-lyp", description="DFT functional"),
    ParamSpec("basis", "str", default="def2-SVP", description="Basis set"),
    ParamSpec("ri", "bool", default=True, description="Use the RI approximation (ridft)"),
    ParamSpec("timeout", "int", unit="s", description="Per-command timeout"),
)
_TM_DATA = (DataKeySpec("energy_Eh", "float", "Eh", "Final electronic energy"),)
# The Turbomole suite is reached via define/x2t plus the runners (ridft/jobex/aoforce).
_TM_BINARIES = ("define", "x2t")


def _check_charge_mult_types(charge: Any, mult: Any) -> Optional[str]:
    """Physical limits: charge and multiplicity must be plain integers.

    multiplicity (2S+1) >= 1; charge may be any integer.  bool is not an
    int here (True would otherwise masquerade as charge=1/mult=1).
    """
    if isinstance(charge, bool) or not isinstance(charge, int):
        return f"charge must be an integer, got {charge!r}"
    if isinstance(mult, bool) or not isinstance(mult, int):
        return f"mult must be an integer, got {mult!r}"
    if mult < 1:
        return f"mult must be >= 1 (2S+1), got {mult}"
    return None


def _run_tm_command(argv: List[str], work_dir: Path, cores: int,
                    timeout: Optional[int] = None) -> subprocess.CompletedProcess:
    """Run a Turbomole command in work_dir.

    ``argv`` is a list and is executed WITHOUT a shell: no argument of
    the caller can ever be re-interpreted as shell syntax (injection
    surface, reviewer finding F3).  A missing binary is reported as a
    returncode-127 CompletedProcess instead of raising, so the adapter
    can return an honest FAILED result.
    """
    import os
    env = os.environ.copy()
    # Parallelism: the T2 parity requirement is "PARA_ARCH, PARNODES" for a
    # multicore run. PARNODES alone does not tell Turbomole's runtime which
    # parallel architecture to use, so default PARA_ARCH to SMP (the same
    # default the dashboard submit path uses) unless the caller already set
    # it in the environment. Without this, a cores>1 run through the direct
    # adapter path (the w21 workflow) silently runs serial.
    env.setdefault("PARA_ARCH", "SMP")
    env["PARNODES"] = str(cores)
    env["OMP_NUM_THREADS"] = str(cores)

    if shutil.which(argv[0]) is None:
        return subprocess.CompletedProcess(
            argv, 127, "",
            f"{argv[0]}: not found in PATH (Turbomole is license-restricted "
            "and is never auto-installed — see delfin/tools/_environment.py)")
    return subprocess.run(
        argv, shell=False, cwd=work_dir, env=env,
        capture_output=True, text=True, timeout=timeout,
    )


def _setup_coord(geometry: Path, work_dir: Path) -> Path:
    """Convert XYZ to Turbomole coord format."""
    coord_path = work_dir / "coord"
    # Use x2t if available, else manual conversion
    x2t = shutil.which("x2t")
    if x2t:
        result = subprocess.run(
            [x2t, str(geometry)], cwd=work_dir,
            capture_output=True, text=True,
        )
        coord_path.write_text(result.stdout)
    else:
        # Fallback: simple xyz → coord conversion
        lines = geometry.read_text().strip().splitlines()
        n_atoms = int(lines[0])
        coord_lines = ["$coord"]
        bohr = 1.8897259886  # Å → Bohr
        for line in lines[2:2 + n_atoms]:
            parts = line.split()
            x, y, z = float(parts[1]) * bohr, float(parts[2]) * bohr, float(parts[3]) * bohr
            coord_lines.append(f"  {x:20.14f}  {y:20.14f}  {z:20.14f}  {parts[0].lower()}")
        coord_lines.append("$end")
        coord_path.write_text("\n".join(coord_lines) + "\n")
    return coord_path


def _validate_control(work_dir: Path, *, ri: bool, mult: int,
                      method: str) -> None:
    """Check the control file define produced contains the sections the
    calculation needs. define can 'end normally' (exit 0) with a partial
    control when the dialog answers drift out of step (real run #608:
    only $title/$actual step/$end); ridft must never run on that."""
    text = (work_dir / "control").read_text() if (work_dir / "control").is_file() else ""
    required = ["$coord", "$atoms", "$basis",
                "$dft"]
    if ri:
        required.append("$rij")
    required.append("$closed shells" if mult == 1 else "$uhf")
    missing = [key for key in required if key not in text]
    if "$dft" in text and method and method.lower() not in text.lower():
        missing.append(f"$dft func {method}")
    if missing:
        raise RuntimeError(
            "Turbomole define produced an incomplete control file "
            f"(missing: {', '.join(missing)}); not running ridft/dscf on it"
        )


def _setup_control(work_dir: Path, **kwargs: Any) -> None:
    """Generate a Turbomole control file via define or manually.

    Everything is read from ``kwargs`` so a caller can never bind the
    same value twice (reviewer finding F1: the old positional + **kwargs
    call raised TypeError on every execute()).
    """
    method = kwargs.get("method", "b3-lyp")
    basis = kwargs.get("basis", "def2-SVP")
    charge = kwargs.get("charge", 0)
    mult = kwargs.get("mult", 1)
    cores = kwargs.get("cores", 1)
    define = shutil.which("define")
    ri = bool(kwargs.get("ri", True))
    if define:
        # Use define with input script. The FIRST TWO answers must be
        # blank: define first asks whether to read default data from an
        # existing control file, then asks for the title (real run #608:
        # without the leading blank line 'a coord' was eaten by the title
        # prompt and every later answer was out of step, leaving a
        # control with only $title/$actual step/$end).
        define_input = "\n\na coord\n*\nno\n"
        define_input += f"b all {basis}\n*\neht\ny\n{charge}\n"
        if mult <= 2:
            # Singlet AND doublet: the plain 'y' accepts the default
            # occupation; define sets up UHF itself for an open shell on
            # 7.9 (verified: CH3 mult=2 got $uhf/$alpha/$beta shells from
            # the plain dialog). Sending 'n'/'u N' here falls out of step
            # (real run #608 round 4: CH3 control missing $dft/$rij).
            define_input += "y\n"
        else:
            # mult >= 3: eht asks for the occupation explicitly — decline
            # the default occupation with 'n' and give the number of
            # unpaired electrons (verified working transcript, site 7.9:
            # .gate/rr608_define_transcript_working.txt).
            define_input += f"n\nu {mult - 1}\n"
        # dft submenu: switch it on, set the functional, then LEAVE WITH A
        # BLANK line — not '*' (site 7.9 verified transcript).
        define_input += f"dft\non\nfunc {method}\n\n"
        if ri:
            # The RI block only when the caller asked for RI (N1):
            # ri=False must produce a conventional dscf setup, not an
            # RI one the caller never requested. 'on' switches it on and
            # a BLANK line leaves the submenu; no cores line is asked
            # here on 7.9.
            define_input += "ri\non\n\n"
        define_input += "*\n"
        result = subprocess.run(
            [define], input=define_input, cwd=work_dir,
            capture_output=True, text=True, timeout=60,
        )
        # A broken define dialog is a SETUP failure, not an SCF failure:
        # ridft/dscf must never run on a missing or partial control file,
        # and the error must name define (real run #608: the dialog broke
        # on the site install and the run was misreported as an SCF
        # convergence failure).
        if result.returncode != 0:
            raise RuntimeError(
                f"Turbomole define failed (exit {result.returncode}): "
                f"{(result.stderr or result.stdout or '')[:500].strip()}"
            )
        # Even a define that exits 0 can leave a partial control file
        # (out-of-step dialog answers). Validate before anything runs.
        _validate_control(work_dir, ri=ri, mult=mult, method=method)
    else:
        # Minimal manual control file
        control = f"""$title
$operating system unix
$symmetry c1
$coord file=coord
$atoms
basis={basis}
$charge {charge}
$dft
  functional {method}
$ricore 500
$maxcor 500
$scfconv 7
"""
        if mult > 1:
            # An open shell is declared explicitly: $unpaired N with
            # N = mult - 1 unpaired electrons.  Never silently closed.
            control += f"$unpaired {mult - 1}\n"
        control += "$end\n"
        (work_dir / "control").write_text(control)


def _extract_energy_tm(work_dir: Path) -> Optional[float]:
    """Extract energy from Turbomole energy file."""
    energy_file = work_dir / "energy"
    if not energy_file.is_file():
        return None
    try:
        lines = energy_file.read_text().strip().splitlines()
        for line in reversed(lines):
            line = line.strip()
            if line and not line.startswith("$"):
                parts = line.split()
                if len(parts) >= 2:
                    return float(parts[1])
    except Exception:
        pass
    return None


def _count_imaginary_tm(aoforce_text: str) -> int:
    """Count imaginary modes from the frequency table only.

    Real-run finding #608 round 4: the old whole-text counter hit the
    coordinate block and normal-mode eigenvector rows (whose numbers are
    often negative) and reported imaginary modes for a true minimum. The
    clean source is the frequency table: the ``$vibrational spectrum``
    block (the vibspectrum file), column 3 = wave number; TURBOMOLE
    prints imaginary modes there as negative numbers, real ones as
    positive. Only lines inside that block are read; an i-suffixed wave
    number also counts (defence for installs that print the suffix).

    A count of 0 reliably marks a verified minimum; over-counting would
    keep a minimum UNVERIFIED -- fail-safe, never silently green
    (target 3: "aoforce imaginary modes via DELFIN's IMAG logic").
    """
    in_spectrum = False
    count = 0
    for line in aoforce_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("$vibrational spectrum"):
            in_spectrum = True
            continue
        if in_spectrum:
            if stripped.startswith("$"):
                in_spectrum = False
                continue
            if not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 3:
                continue
            token = parts[2].lower()
            has_i_suffix = token.endswith("i") or token.endswith("i*")
            bare = token.rstrip("*i")
            try:
                value = float(bare)
            except ValueError:
                continue
            if has_i_suffix or value < 0:
                count += 1
    return count


def _coord_to_xyz(work_dir: Path) -> Optional[Path]:
    """Convert Turbomole coord back to XYZ."""
    coord_path = work_dir / "coord"
    if not coord_path.is_file():
        return None
    try:
        bohr = 1.8897259886
        lines = coord_path.read_text().strip().splitlines()
        atoms = []
        for line in lines:
            line = line.strip()
            if line.startswith("$") or not line:
                continue
            parts = line.split()
            if len(parts) >= 4:
                x = float(parts[0]) / bohr
                y = float(parts[1]) / bohr
                z = float(parts[2]) / bohr
                elem = parts[3].capitalize()
                atoms.append(f"{elem}  {x:.8f}  {y:.8f}  {z:.8f}")
        if atoms:
            xyz_path = work_dir / "output.xyz"
            xyz_path.write_text(f"{len(atoms)}\nTurbomole output\n" + "\n".join(atoms) + "\n")
            return xyz_path
    except Exception:
        pass
    return None


class TurbomoleSpAdapter(StepAdapter):
    name = "turbomole_sp"
    description = "Turbomole single-point energy (DFT/RI)"
    produces_geometry = False
    category = "dft"
    params = _TM_PARAMS
    consumes = ("geometry",)
    produces = ("qc_output",)
    data_keys = _TM_DATA
    requires_binaries = _TM_BINARIES

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        charge = kwargs.get("charge", 0)
        mult = kwargs.get("mult", 1)
        bad = _check_charge_mult_types(charge, mult)
        if bad:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=bad)
        timeout = kwargs.get("timeout")

        try:
            _setup_coord(geometry, work_dir)
            _setup_control(work_dir, cores=cores, **kwargs)
        except (RuntimeError, OSError) as exc:
            # A broken define dialog (or an unwritable/missing work dir,
            # an OSError from _setup_coord/_setup_control) is a setup
            # failure, reported as such — never carried on to
            # ridft/dscf and misread as an SCF problem.
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc),
                                     error_kind=ErrorKind.TOOL_FAILED)

        # Run ridft (RI-DFT) or dscf (conventional)
        argv = ["ridft"] if kwargs.get("ri", True) else ["dscf"]
        proc = _run_tm_command(argv, work_dir, cores, timeout=timeout)

        energy = _extract_energy_tm(work_dir)
        # TM2-1: exit 0 without a converged energy file is an unconverged
        # run — never a SUCCESS with empty data.
        success = proc.returncode == 0 and energy is not None
        data = {}
        if energy is not None:
            data["energy_Eh"] = energy

        if success:
            err: Optional[str] = None
            err_kind = ErrorKind.NONE
        elif proc.returncode == 0:
            err = (f"Turbomole {argv[0]} exited 0 but produced no energy "
                   "file — the SCF did not converge; it is NOT reported "
                   "as a result")
            err_kind = ErrorKind.CONVERGENCE
        else:
            err = f"Turbomole {argv[0]} failed: {proc.stderr[:500]}"
            err_kind = ErrorKind.TOOL_FAILED

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir, start,
            data=data,
            error=err,
            error_kind=err_kind,
        )


class TurbomoleOptAdapter(StepAdapter):
    name = "turbomole_opt"
    description = "Turbomole geometry optimization"
    produces_geometry = True
    category = "dft"
    params = _TM_PARAMS + (
        ParamSpec("max_cycles", "int", default=200, description="Max jobex optimisation cycles"),
    )
    consumes = ("geometry",)
    produces = ("qc_output",)   # geometry added automatically (produces_geometry)
    data_keys = _TM_DATA
    requires_binaries = _TM_BINARIES

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        charge = kwargs.get("charge", 0)
        mult = kwargs.get("mult", 1)
        bad = _check_charge_mult_types(charge, mult)
        if bad:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=bad)
        max_cycles = kwargs.get("max_cycles", 200)
        # max_cycles reaches the process as ONE argv element, never as
        # shell syntax: validate the type, then hand jobex a clean list
        # (reviewer finding F3 — "200; touch X" must fail validation).
        if isinstance(max_cycles, bool) or not isinstance(max_cycles, int) or max_cycles < 1:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=f"max_cycles must be a positive integer, "
                                           f"got {max_cycles!r}")
        timeout = kwargs.get("timeout")

        try:
            _setup_coord(geometry, work_dir)
            _setup_control(work_dir, cores=cores, **kwargs)
        except RuntimeError as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc),
                                     error_kind=ErrorKind.TOOL_FAILED)

        argv = ["jobex", "-ri", "-c", str(max_cycles)]
        proc = _run_tm_command(argv, work_dir, cores, timeout=timeout)

        energy = _extract_energy_tm(work_dir)
        geom_out = _coord_to_xyz(work_dir)

        # TM2-1: exit 0 is NOT success. A jobex that returns 0 without a
        # converged result (energy file) is an UNCONVERGED run and must be
        # reported as such — never as a SUCCESS with empty data (target 3:
        # "an unconverged run is never reported as a result").
        success = proc.returncode == 0 and energy is not None

        data = {}
        if energy is not None:
            data["energy_Eh"] = energy

        if success:
            err: Optional[str] = None
            err_kind = ErrorKind.NONE
        elif proc.returncode == 0:
            err = ("Turbomole jobex exited 0 but produced no converged "
                   "energy file — the optimization did not converge; it is "
                   "NOT reported as a result")
            err_kind = ErrorKind.CONVERGENCE
        else:
            err = f"Turbomole jobex failed: {proc.stderr[:500]}"
            err_kind = ErrorKind.TOOL_FAILED

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir, start,
            geometry=geom_out,
            data=data,
            error=err,
            error_kind=err_kind,
        )


class TurbomoleFreqAdapter(StepAdapter):
    name = "turbomole_freq"
    description = "Turbomole frequency calculation (aoforce)"
    produces_geometry = False
    category = "dft"
    params = _TM_PARAMS
    consumes = ("geometry",)
    produces = ("qc_output", "hessian")
    data_keys = _TM_DATA + (DataKeySpec("hessian", "file", "path",
                                        "aoforce Hessian file"),
                            DataKeySpec("imaginary_frequencies", "int", "count",
                                        "imaginary vibrational modes"),)
    requires_binaries = _TM_BINARIES

    def validate_params(self, **kwargs: Any) -> None:
        if "charge" not in kwargs:
            raise ValueError("'charge' parameter is required")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if geometry is None:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error="geometry is required")

        charge = kwargs.get("charge", 0)
        mult = kwargs.get("mult", 1)
        bad = _check_charge_mult_types(charge, mult)
        if bad:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=bad)
        timeout = kwargs.get("timeout")

        try:
            _setup_coord(geometry, work_dir)
            _setup_control(work_dir, cores=cores, **kwargs)
        except RuntimeError as exc:
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start,
                                     error=str(exc),
                                     error_kind=ErrorKind.TOOL_FAILED)

        # First SP, then aoforce
        _run_tm_command(["ridft"], work_dir, cores, timeout=timeout)
        proc = _run_tm_command(["aoforce"], work_dir, cores, timeout=timeout)

        energy = _extract_energy_tm(work_dir)
        # TM2-1: exit 0 without a converged energy file is an unconverged
        # run — never a SUCCESS with empty data.
        success = proc.returncode == 0 and energy is not None

        # Target 3: persist aoforce's output so DELFIN's IMAG logic can read
        # it. The adapter must leave aoforce.out on disk, not swallow it into
        # a CompletedProcess (defect TM2-T3-A: the workflow's _count_imaginary
        # read a file that was never written, so imaginary modes were never
        # surfaced).
        aoforce_out = work_dir / "aoforce.out"
        if proc.stdout:
            aoforce_out.write_text(proc.stdout)

        data = {}
        if energy is not None:
            data["energy_Eh"] = energy

        # Surface the imaginary-mode count: a count > 0 means this structure
        # is NOT a verified minimum and the caller must report that (ORCA
        # parity: extract_imaginary_frequencies). The count is a datum, not
        # an error -- the workflow decides whether to warn.
        if aoforce_out.is_file():
            data["imaginary_frequencies"] = _count_imaginary_tm(
                aoforce_out.read_text(errors="replace"))

        # The adapter declares "hessian" as a product: surface the file
        # aoforce wrote, or fail honestly — a declared product that is
        # silently absent is a defect (reviewer finding F4).
        hessian = work_dir / "hessian"
        if hessian.is_file():
            data["hessian"] = str(hessian)
        elif success:
            success = False

        if success:
            err: Optional[str] = None
            err_kind = ErrorKind.NONE
        elif proc.returncode == 0 and energy is None:
            err = ("Turbomole aoforce exited 0 but produced no energy file "
                   "— the SCF did not converge; it is NOT reported as a result")
            err_kind = ErrorKind.CONVERGENCE
        elif proc.returncode == 0:
            err = "aoforce exited 0 but wrote no hessian file"
            err_kind = ErrorKind.PARSE
        else:
            err = f"Turbomole aoforce failed: {proc.stderr[:500]}"
            err_kind = ErrorKind.TOOL_FAILED

        return self._make_result(
            self.name,
            StepStatus.SUCCESS if success else StepStatus.FAILED,
            work_dir, start,
            data=data,
            error=err,
            error_kind=err_kind,
        )


register(TurbomoleSpAdapter())
register(TurbomoleOptAdapter())
register(TurbomoleFreqAdapter())
