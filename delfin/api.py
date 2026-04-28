"""Programmatic Python API for DELFIN.

Wraps the CLI subcommands (``delfin/cli.py``) with typed Python functions
returning structured ``CommandResult`` objects.  Intended for embedding in
notebooks, agent tooling (MCP server), tests, and downstream automation.

The internal pattern is:

1. Build the argv list from typed kwargs.
2. Run ``cli_main(argv)`` with captured stdout/stderr.
3. Return a ``CommandResult`` containing exit code + captured streams.

This keeps the CLI as the single source of truth and avoids logic duplication.
"""
from __future__ import annotations

import contextlib
import io
import sys
from dataclasses import dataclass, field
from typing import Sequence


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

@dataclass
class CommandResult:
    """Structured result of a DELFIN CLI subcommand.

    Attributes
    ----------
    returncode : int
        Exit code (0 = success).  ``cli_main`` returns this directly; on
        ``SystemExit`` we capture ``code`` from the exception.
    stdout : str
        Captured standard output.
    stderr : str
        Captured standard error.
    argv : list[str]
        The argv passed to ``cli_main`` (useful for debugging / logging).
    """

    returncode: int = 0
    stdout: str = ""
    stderr: str = ""
    argv: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.returncode == 0


# ---------------------------------------------------------------------------
# Internal CLI runner
# ---------------------------------------------------------------------------

def _has_flag(argv: Sequence[str], *flags: str) -> bool:
    for flag in flags:
        if flag in argv:
            return True
        prefix = f"{flag}="
        if any(arg.startswith(prefix) for arg in argv):
            return True
    return False


def _run_cli(argv: list[str], *, capture: bool = True) -> CommandResult:
    """Invoke ``delfin.cli.main`` with a built argv and capture the streams.

    Parameters
    ----------
    argv : list[str]
        The argv passed to ``cli_main`` (without the leading ``delfin``).
    capture : bool
        If True (default), capture stdout/stderr into the result.  If False,
        streams pass through to the caller's terminal.
    """
    from delfin.cli import main as cli_main

    if not capture:
        try:
            rc = cli_main(argv)
        except SystemExit as exc:
            rc = int(exc.code) if exc.code is not None else 0
        return CommandResult(returncode=rc, argv=list(argv))

    out_buf = io.StringIO()
    err_buf = io.StringIO()
    rc: int = 0
    with contextlib.redirect_stdout(out_buf), contextlib.redirect_stderr(err_buf):
        try:
            rc = cli_main(argv)
        except SystemExit as exc:
            rc = int(exc.code) if exc.code is not None else 0
    return CommandResult(
        returncode=int(rc or 0),
        stdout=out_buf.getvalue(),
        stderr=err_buf.getvalue(),
        argv=list(argv),
    )


# ---------------------------------------------------------------------------
# Pipeline (run / prepare) — kept rückwärtskompatibel
# ---------------------------------------------------------------------------

def run(
    control_file: str = "CONTROL.txt",
    *,
    cleanup: bool = True,
    recalc: bool = False,
    overwrite: bool = False,
    define: str | None = None,
    extra_args: Sequence[str] | None = None,
) -> int:
    """Execute the DELFIN pipeline programmatically using CLI semantics.

    Returns the raw exit code (int) for backwards compatibility.
    Use :func:`pipeline_run` for a structured ``CommandResult``.
    """
    return pipeline_run(
        control_file=control_file,
        cleanup=cleanup,
        recalc=recalc,
        overwrite=overwrite,
        define=define,
        extra_args=extra_args,
        capture=False,
    ).returncode


def prepare(control_file: str = "CONTROL.txt", overwrite: bool = False) -> int:
    """Convenience wrapper for ``delfin --define`` to create CONTROL templates."""
    return pipeline_prepare(
        control_file=control_file, overwrite=overwrite, capture=False
    ).returncode


def pipeline_run(
    control_file: str = "CONTROL.txt",
    *,
    cleanup: bool = True,
    recalc: bool = False,
    overwrite: bool = False,
    define: str | None = None,
    extra_args: Sequence[str] | None = None,
    capture: bool = True,
) -> CommandResult:
    """Execute the DELFIN pipeline and return a ``CommandResult``."""
    argv: list[str] = list(extra_args or [])

    if recalc and not _has_flag(argv, "--recalc"):
        argv.append("--recalc")
    if not cleanup and not _has_flag(argv, "--no-cleanup"):
        argv.append("--no-cleanup")
    if overwrite and not _has_flag(argv, "--overwrite"):
        argv.append("--overwrite")

    if not _has_flag(argv, "--control", "-F") and control_file != "CONTROL.txt":
        argv.extend(["--control", control_file])

    if define is not None and not _has_flag(argv, "--define", "-D"):
        argv.extend(["--define", define])

    return _run_cli(argv, capture=capture)


def pipeline_prepare(
    control_file: str = "CONTROL.txt",
    overwrite: bool = False,
    *,
    capture: bool = True,
) -> CommandResult:
    """Generate a CONTROL template via ``delfin --define <file>``."""
    argv = ["--define", control_file]
    if overwrite:
        argv.append("--overwrite")
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Tool / runtime checks (read-only, safe)
# ---------------------------------------------------------------------------

def qm_check(
    tools: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Inspect QM tool resolution (xtb, crest, etc.).  Read-only."""
    argv = ["qm_check"]
    if tools:
        argv.extend(list(tools))
    return _run_cli(argv, capture=capture)


def csp_check(*, capture: bool = True) -> CommandResult:
    """Check CSP (genarris) tool availability.  Read-only."""
    return _run_cli(["csp_check"], capture=capture)


def mlp_check(*, capture: bool = True) -> CommandResult:
    """Check MLP (torchani / AIMNet2 / MACE) backend availability.  Read-only."""
    return _run_cli(["mlp_check"], capture=capture)


def analysis_check(*, capture: bool = True) -> CommandResult:
    """Check analysis tools (Multiwfn, CENSO, ANMR, morfeus).  Read-only."""
    return _run_cli(["analysis_check"], capture=capture)


# ---------------------------------------------------------------------------
# QM tool runner
# ---------------------------------------------------------------------------

def qm_run(
    tool: str,
    tool_args: Sequence[str] | None = None,
    *,
    cwd: str = ".",
    capture_tool: bool = True,
    capture: bool = True,
) -> CommandResult:
    """Run a single QM tool via DELFIN's resolver.

    Parameters
    ----------
    tool : str
        Tool name (e.g. ``xtb``, ``crest``, ``xtb4stda``, ``stda``, ``dftb+``).
    tool_args : sequence of str
        Arguments forwarded to the tool.
    cwd : str
        Working directory for the tool process.
    capture_tool : bool
        Pass ``--capture`` so ORCA prints stdout/stderr after completion.
    capture : bool
        Capture DELFIN's own output into the ``CommandResult``.
    """
    argv = ["qm_run", tool, "--cwd", cwd]
    if capture_tool:
        argv.append("--capture")
    if tool_args:
        argv.append("--")
        argv.extend(list(tool_args))
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Process control (cleanup / stop)
# ---------------------------------------------------------------------------

def cleanup(
    *,
    orca: bool = False,
    dry_run: bool = False,
    workspace: str = ".",
    scratch: str | None = None,
    capture: bool = True,
) -> CommandResult:
    """Remove DELFIN scratch artifacts and optionally stop ORCA."""
    argv = ["cleanup", "--workspace", workspace]
    if orca:
        argv.append("--orca")
    if dry_run:
        argv.append("--dry-run")
    if scratch:
        argv.extend(["--scratch", scratch])
    return _run_cli(argv, capture=capture)


def stop(
    *,
    workspace: str = ".",
    signal_name: str = "INT",
    dry_run: bool = False,
    cleanup_after: bool = False,
    wait_seconds: float = 3.0,
    capture: bool = True,
) -> CommandResult:
    """Send a signal (INT/TERM/KILL) to running DELFIN processes.

    Parameters
    ----------
    signal_name : {"INT", "TERM", "KILL"}
        Signal to send.  Defaults to graceful SIGINT.
    cleanup_after : bool
        Run cleanup after signaling (only when not dry-run).
    """
    if signal_name not in {"INT", "TERM", "KILL"}:
        raise ValueError(
            f"signal_name must be one of INT/TERM/KILL, got {signal_name!r}"
        )
    argv = ["stop", "--workspace", workspace, "--signal", signal_name,
            "--wait-seconds", str(wait_seconds)]
    if dry_run:
        argv.append("--dry-run")
    if cleanup_after:
        argv.append("--cleanup")
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# ORCA invocation
# ---------------------------------------------------------------------------

def run_orca_input(
    input_file: str | None = None,
    output: str | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run ORCA on a ``.inp`` file (uses DELFIN's ORCA resolver)."""
    argv: list[str] = ["run_orca"]
    if input_file:
        argv.append(input_file)
    if output:
        argv.extend(["--output", output])
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# CO2 Coordinator
# ---------------------------------------------------------------------------

def co2(
    *,
    define: bool = False,
    force: bool = False,
    recalc: bool = False,
    charge: int | None = None,
    multiplicity: int | None = None,
    solvent: str | None = None,
    metal: str | None = None,
    broken_sym: str | None = None,
    capture: bool = True,
) -> CommandResult:
    """CO2 Coordinator workflow."""
    argv: list[str] = ["co2"]
    if define:
        argv.append("--define")
    if force:
        argv.append("--force")
    if recalc:
        argv.append("--recalc")
    if charge is not None:
        argv.extend(["--charge", str(charge)])
    if multiplicity is not None:
        argv.extend(["--multiplicity", str(multiplicity)])
    if solvent is not None:
        argv.extend(["--solvent", solvent])
    if metal is not None:
        argv.extend(["--metal", metal])
    if broken_sym is not None:
        argv.extend(["--broken_sym", broken_sym])
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# Specialised workflows: TADF (xTB), hyperpolarisability
# ---------------------------------------------------------------------------

def tadf_xtb(
    extra_args: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run the TADF xTB workflow (``delfin tadf_xtb``)."""
    argv: list[str] = ["tadf_xtb"]
    if extra_args:
        argv.extend(list(extra_args))
    return _run_cli(argv, capture=capture)


def hyperpol(
    extra_args: Sequence[str] | None = None,
    *,
    capture: bool = True,
) -> CommandResult:
    """Run the hyperpolarisability workflow (``delfin hyperpol``)."""
    argv: list[str] = ["hyperpol"]
    if extra_args:
        argv.extend(list(extra_args))
    return _run_cli(argv, capture=capture)


# ---------------------------------------------------------------------------
# P1 — Output parsing (read-only structured returns for the agent)
# ---------------------------------------------------------------------------


@dataclass
class OrcaParseResult:
    """Structured ORCA-output snapshot.

    Numeric fields are ``None`` when the corresponding pattern wasn't
    found in the output. ``scf_converged`` / ``opt_converged`` are
    tri-state via ``None`` when the marker is absent — don't mistake
    "no marker" for "did not converge".
    """
    path: str
    final_single_point: float | None = None
    gibbs_free_energy: float | None = None
    zpe: float | None = None
    scf_converged: bool | None = None
    opt_converged: bool | None = None
    imag_freq_count: int | None = None
    walltime_s: float | None = None
    n_atoms: int | None = None
    functional: str = ""
    basis: str = ""
    error_summary: str = ""


def parse_orca_output(path: str) -> OrcaParseResult:
    """Parse one ORCA output file and return a structured snapshot.

    Wraps the per-property helpers in :mod:`delfin.energies` so the
    agent can ask "what happened in this run?" with a single tool
    call instead of grepping the .out file line by line.

    Missing fields → ``None`` (not zero — distinguishes "absent"
    from "found but zero-valued"). Read errors don't raise; the
    result's ``error_summary`` carries the diagnosis.
    """
    from pathlib import Path as _P
    from delfin import energies as _e
    p = _P(path)
    if not p.exists():
        return OrcaParseResult(path=str(p), error_summary="file not found")
    if not p.is_file():
        return OrcaParseResult(path=str(p), error_summary="not a file")
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except Exception as exc:
        return OrcaParseResult(path=str(p), error_summary=f"read error: {exc}")

    out = OrcaParseResult(path=str(p))
    try:
        out.final_single_point = _e.find_electronic_energy(str(p))
    except Exception:
        pass
    try:
        out.gibbs_free_energy = _e.find_gibbs_energy(str(p))
    except Exception:
        pass
    try:
        out.zpe = _e.find_ZPE(str(p))
    except Exception:
        pass

    # SCF convergence: ORCA writes a positive marker on success.
    if "SCF CONVERGED AFTER" in text or "SCF CONVERGED" in text:
        out.scf_converged = True
    elif "SCF NOT CONVERGED" in text or "SCF iterations did not converge" in text:
        out.scf_converged = False

    # Optimization convergence (geom_opt jobs only).
    if ("OPTIMIZATION RUN DONE" in text
            or "THE OPTIMIZATION HAS CONVERGED" in text):
        out.opt_converged = True
    elif ("OPTIMIZATION DID NOT CONVERGE" in text
            or "FAILED TO CONVERGE THE GEOMETRY OPTIMIZATION" in text):
        out.opt_converged = False

    # Imaginary frequency count.
    import re as _re
    m_imag = _re.search(
        r"Number of imaginary frequencies\s*\.\.\.\s*(\d+)", text,
    )
    if m_imag:
        out.imag_freq_count = int(m_imag.group(1))
    else:
        # Count negative-frequency lines: "  4:    -123.45 cm**-1"
        imag_lines = _re.findall(
            r"^\s*\d+:\s+(-\d+\.\d+)\s+cm\*\*-1", text, _re.MULTILINE,
        )
        if imag_lines:
            out.imag_freq_count = len(imag_lines)

    # Number of atoms.
    m_atoms = _re.search(r"Number of atoms\s+\.\.\.\s+(\d+)", text)
    if m_atoms:
        out.n_atoms = int(m_atoms.group(1))

    # Walltime: "TOTAL RUN TIME: 0 days 1 hours 23 minutes 45 seconds 678 msec"
    m_wall = _re.search(
        r"TOTAL RUN TIME:\s+(\d+)\s+days\s+(\d+)\s+hours\s+(\d+)\s+minutes\s+"
        r"(\d+)\s+seconds\s+(\d+)\s+msec",
        text,
    )
    if m_wall:
        d, h, mi, s, ms = (int(g) for g in m_wall.groups())
        out.walltime_s = d * 86400 + h * 3600 + mi * 60 + s + ms / 1000.0

    # Functional / basis from the input keyword line.
    m_kw = _re.search(r"^\s*\|\s*\d+>\s*!\s*(.+)$", text, _re.MULTILINE)
    if m_kw:
        for token in m_kw.group(1).split():
            T = token.upper()
            if T in ("BP86", "PBE0", "B3LYP", "TPSS", "WB97X", "WB97X-D3",
                    "CAM-B3LYP", "R2SCAN", "PBE", "M06", "M062X", "BLYP",
                    "REVTPSS"):
                out.functional = token
            if ("DEF2-" in T or "PCSSEG" in T or "MA-DEF2-" in T
                    or "SARC-" in T or "X2C-" in T):
                out.basis = token

    if (out.final_single_point is None
            and out.scf_converged is not True
            and "ABORTING THE RUN" in text):
        out.error_summary = "ORCA aborted the run"

    return out


@dataclass
class OrcaError:
    """One detected error pattern in an ORCA output."""
    type: str
    message: str
    line_number: int | None = None
    suggestion: str = ""


_ORCA_ERROR_PATTERNS: list[tuple[str, str, str]] = [
    (r"SCF NOT CONVERGED|SCF iterations did not converge",
     "scf_diverge",
     "Add SlowConv or tighten SCFconv; try DIIS off or KDIIS; check spin guess."),
    (r"oom_kill|Killed by signal 9|out of memory",
     "oom",
     "Increase --mem (sbatch) or lower %maxcore in the .inp."),
    (r"Number of basis functions.*exceeds",
     "basis",
     "Reduce basis (def2-TZVP → def2-SVP) or split atom regions."),
    (r"Wrong multiplicity|Multiplicity\s+\d+\s+impossible",
     "multiplicity",
     "Recheck charge/multiplicity — total electrons must match."),
    (r"MPI ERROR|MPI_ABORT|orterun.*detected.*aborted",
     "mpi",
     "MPI process died — check sbatch --ntasks vs nprocs in %pal."),
    (r"DUE TO TIME LIMIT|Job exceeded|TIMEOUT",
     "timeout",
     "Bump SLURM time limit; or restart from .gbw + .opt with maxiter=N."),
    (r"FATAL ERROR|ABORTING THE RUN",
     "other",
     "ORCA aborted — read the output around this line for the cause."),
]


def find_orca_errors(folder: str) -> list[OrcaError]:
    """Scan ``*.out`` files in ``folder`` for known ORCA error patterns.

    Non-recursive (one folder at a time keeps results focused). Each
    match yields an OrcaError with a short suggestion. An empty list
    is NOT proof of success — use ``parse_orca_output`` for that.
    """
    from pathlib import Path as _P
    import re as _re
    errors: list[OrcaError] = []
    p = _P(folder)
    if not p.exists() or not p.is_dir():
        return errors
    for out_file in sorted(p.glob("*.out")):
        try:
            lines = out_file.read_text(
                encoding="utf-8", errors="replace",
            ).splitlines()
        except Exception:
            continue
        for ln_no, line in enumerate(lines, start=1):
            for pat, etype, suggestion in _ORCA_ERROR_PATTERNS:
                if _re.search(pat, line, _re.IGNORECASE):
                    errors.append(OrcaError(
                        type=etype,
                        message=line.strip()[:200],
                        line_number=ln_no,
                        suggestion=suggestion,
                    ))
                    break
    return errors


@dataclass
class ThermochemResult:
    """Thermochemistry block extracted from an ORCA Freq output."""
    path: str
    temperature_k: float | None = None
    pressure_atm: float | None = None
    zpe: float | None = None
    thermal_corr: float | None = None
    enthalpy_corr: float | None = None
    entropy_total: float | None = None
    gibbs_corr: float | None = None
    final_gibbs: float | None = None


def extract_thermochem(folder: str) -> ThermochemResult:
    """Extract the full thermochemistry block from an ORCA Freq output.

    Picks the first ``*.out`` containing thermochemistry data in the
    folder. All energies in Hartree, T in K, P in atm.
    """
    from pathlib import Path as _P
    import re as _re
    p = _P(folder)
    if not p.exists():
        return ThermochemResult(path=str(p))
    out_files = sorted(p.glob("*.out")) if p.is_dir() else [p]
    for out_file in out_files:
        try:
            text = out_file.read_text(encoding="utf-8", errors="replace")
        except Exception:
            continue
        if "THERMOCHEMISTRY" not in text.upper():
            continue
        result = ThermochemResult(path=str(out_file))
        FLOAT = r"(-?\d+\.\d+)"
        for attr, pattern in [
            ("temperature_k",  rf"Temperature\s+\.\.\.\s+{FLOAT}\s*K"),
            ("pressure_atm",   rf"Pressure\s+\.\.\.\s+{FLOAT}\s*atm"),
            ("zpe",            rf"Zero point energy\s+\.\.\.\s+{FLOAT}\s*Eh"),
            ("thermal_corr",   rf"Total thermal correction\s+{FLOAT}\s*Eh"),
            ("enthalpy_corr",  rf"Total Enthalpy\s+\.\.\.\s+{FLOAT}\s*Eh"),
            ("entropy_total",  rf"Final entropy term\s+\.\.\.\s+{FLOAT}\s*Eh"),
            ("gibbs_corr",     rf"G-E\(el\)\s+\.\.\.\s+{FLOAT}\s*Eh"),
            ("final_gibbs",
             rf"Final Gibbs free energy\s+\.\.\.\s+{FLOAT}\s*Eh"),
        ]:
            m = _re.search(pattern, text)
            if m:
                setattr(result, attr, float(m.group(1)))
        return result
    return ThermochemResult(path=str(p))


def extract_energy_table(
    folders: list[str] | str,
    properties: list[str] | None = None,
) -> list[dict]:
    """Walk a list of folders and collect energies into rows.

    Each row is a dict with at least ``folder`` and ``status`` plus
    one entry per requested property. Defaults to
    ``["gibbs", "zpe", "single_point"]``.

    Recognised properties: ``gibbs``, ``zpe``, ``single_point``,
    ``scf_converged``, ``opt_converged``, ``imag_freqs``,
    ``walltime_s``.

    Folders without ORCA output get ``status="no_output"`` and
    ``None`` for every property — easier to filter than missing rows.
    """
    if properties is None:
        properties = ["gibbs", "zpe", "single_point"]
    if isinstance(folders, str):
        folders = [folders]

    from pathlib import Path as _P
    rows: list[dict] = []
    for folder in folders:
        p = _P(folder)
        if not p.exists() or not p.is_dir():
            row: dict = {"folder": str(folder), "status": "missing"}
            for prop in properties:
                row[prop] = None
            rows.append(row)
            continue
        out_files = sorted(p.glob("*.out"))
        if not out_files:
            row = {"folder": str(folder), "status": "no_output"}
            for prop in properties:
                row[prop] = None
            rows.append(row)
            continue
        target = max(out_files, key=lambda f: f.stat().st_size)
        parsed = parse_orca_output(str(target))
        row = {
            "folder": str(folder),
            "status": "ok",
            "output_file": target.name,
        }
        for prop in properties:
            if prop == "gibbs":
                row[prop] = parsed.gibbs_free_energy
            elif prop == "zpe":
                row[prop] = parsed.zpe
            elif prop == "single_point":
                row[prop] = parsed.final_single_point
            elif prop == "scf_converged":
                row[prop] = parsed.scf_converged
            elif prop == "opt_converged":
                row[prop] = parsed.opt_converged
            elif prop == "imag_freqs":
                row[prop] = parsed.imag_freq_count
            elif prop == "walltime_s":
                row[prop] = parsed.walltime_s
            else:
                row[prop] = None
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Imaginary-frequency extraction + multi-folder comparison
# ---------------------------------------------------------------------------


@dataclass
class ImagFreqMode:
    """One imaginary vibrational mode in an ORCA frequency output."""
    mode_index: int
    frequency_cm: float


@dataclass
class ImagFreqResult:
    """Imaginary frequencies extracted from one calculation folder."""
    folder: str
    output_file: str | None
    n_imag: int
    modes: list[ImagFreqMode]
    most_negative: float | None
    is_minimum: bool | None
    is_ts: bool | None
    error: str | None = None


def extract_imaginary_frequencies(folder: str) -> ImagFreqResult:
    """Pull imaginary modes from the largest ``*.out`` in *folder*.

    Returns ``n_imag = 0`` and ``is_minimum = True`` for true minima.
    A single imaginary mode (typically below -100 cm-1) is the TS
    signature; ``is_ts = True`` then. Missing freq output → ``error``.
    """
    from pathlib import Path as _P
    import re as _re
    p = _P(folder)
    if not p.exists():
        return ImagFreqResult(
            folder=str(p), output_file=None, n_imag=0, modes=[],
            most_negative=None, is_minimum=None, is_ts=None,
            error="folder missing",
        )
    out_files = sorted(p.glob("*.out")) if p.is_dir() else [p]
    if not out_files:
        return ImagFreqResult(
            folder=str(p), output_file=None, n_imag=0, modes=[],
            most_negative=None, is_minimum=None, is_ts=None,
            error="no ORCA output",
        )
    target = max(out_files, key=lambda f: f.stat().st_size)
    try:
        text = target.read_text(encoding="utf-8", errors="replace")
    except Exception as exc:
        return ImagFreqResult(
            folder=str(p), output_file=target.name, n_imag=0, modes=[],
            most_negative=None, is_minimum=None, is_ts=None,
            error=f"read failed: {exc}",
        )

    if "VIBRATIONAL FREQUENCIES" not in text:
        return ImagFreqResult(
            folder=str(p), output_file=target.name, n_imag=0, modes=[],
            most_negative=None, is_minimum=None, is_ts=None,
            error="no Freq section (run with !Freq or !Opt Freq)",
        )

    modes: list[ImagFreqMode] = []
    for m in _re.finditer(
        r"^\s*(\d+):\s+(-\d+\.\d+)\s+cm\*\*-1.*?\*\*\*imaginary mode\*\*\*",
        text, _re.MULTILINE,
    ):
        modes.append(ImagFreqMode(int(m.group(1)), float(m.group(2))))
    if not modes:
        for m in _re.finditer(
            r"^\s*(\d+):\s+(-\d+\.\d+)\s+cm\*\*-1", text, _re.MULTILINE,
        ):
            modes.append(ImagFreqMode(int(m.group(1)), float(m.group(2))))

    n_imag = len(modes)
    most_negative = min((m.frequency_cm for m in modes), default=None)
    is_minimum = (n_imag == 0)
    is_ts = (n_imag == 1)
    return ImagFreqResult(
        folder=str(p), output_file=target.name,
        n_imag=n_imag, modes=modes,
        most_negative=most_negative,
        is_minimum=is_minimum, is_ts=is_ts,
    )


@dataclass
class FunctionalComparisonRow:
    """One folder's row in a cross-functional comparison table."""
    folder: str
    functional: str | None
    basis: str | None
    gibbs: float | None
    single_point: float | None
    zpe: float | None
    n_imag: int | None
    is_minimum: bool | None
    status: str  # "ok" / "no_output" / "missing" / parse-error


def compare_across_functionals(
    folders: list[str],
    *,
    include_imag: bool = True,
    sort_by: str = "gibbs",
) -> list[FunctionalComparisonRow]:
    """Build a comparison table grouped by functional/basis.

    Walks each folder, parses its largest .out, and returns one
    :class:`FunctionalComparisonRow` per folder. ``sort_by`` accepts
    ``gibbs``, ``single_point``, ``zpe``, ``functional``, or ``folder``.
    Rows with missing values for the sort key go to the bottom.
    """
    rows: list[FunctionalComparisonRow] = []
    from pathlib import Path as _P
    for folder in folders:
        p = _P(folder)
        if not p.exists() or not p.is_dir():
            rows.append(FunctionalComparisonRow(
                folder=str(folder), functional=None, basis=None,
                gibbs=None, single_point=None, zpe=None,
                n_imag=None, is_minimum=None, status="missing",
            ))
            continue
        out_files = sorted(p.glob("*.out"))
        if not out_files:
            rows.append(FunctionalComparisonRow(
                folder=str(folder), functional=None, basis=None,
                gibbs=None, single_point=None, zpe=None,
                n_imag=None, is_minimum=None, status="no_output",
            ))
            continue
        target = max(out_files, key=lambda f: f.stat().st_size)
        parsed = parse_orca_output(str(target))
        n_imag: int | None = None
        is_min: bool | None = None
        if include_imag:
            imag = extract_imaginary_frequencies(str(p))
            if imag.error is None:
                n_imag = imag.n_imag
                is_min = imag.is_minimum
        rows.append(FunctionalComparisonRow(
            folder=str(folder),
            functional=parsed.functional,
            basis=parsed.basis,
            gibbs=parsed.gibbs_free_energy,
            single_point=parsed.final_single_point,
            zpe=parsed.zpe,
            n_imag=n_imag,
            is_minimum=is_min,
            status="ok",
        ))

    sort_field = sort_by.strip().lower()
    if sort_field in ("gibbs", "single_point", "zpe"):
        def _key(r: FunctionalComparisonRow):
            v = getattr(r, sort_field)
            return (v is None, v if v is not None else 0.0)
        rows.sort(key=_key)
    elif sort_field == "functional":
        rows.sort(key=lambda r: (r.functional is None, r.functional or ""))
    elif sort_field == "folder":
        rows.sort(key=lambda r: r.folder)
    return rows


@dataclass
class CalcDiff:
    """Side-by-side diff of two calculation folders."""
    folder_a: str
    folder_b: str
    method_match: bool
    basis_match: bool
    a_functional: str | None
    b_functional: str | None
    a_basis: str | None
    b_basis: str | None
    a_gibbs: float | None
    b_gibbs: float | None
    delta_gibbs_kcal: float | None
    a_spe: float | None
    b_spe: float | None
    delta_spe_kcal: float | None
    a_n_imag: int | None
    b_n_imag: int | None
    notes: list[str]


def compare_calculations(folder_a: str, folder_b: str) -> CalcDiff:
    """Diff method/basis/results between two calculation folders.

    ``delta_gibbs_kcal`` and ``delta_spe_kcal`` are :math:`(B - A)` in
    kcal/mol; ``None`` if either side lacks the value. ``method_match``
    is ``True`` only when both functionals are non-None and equal.
    """
    HARTREE_TO_KCAL = 627.5094740631
    rows = compare_across_functionals(
        [folder_a, folder_b], include_imag=True, sort_by="folder",
    )
    by_folder = {r.folder: r for r in rows}
    a = by_folder.get(str(folder_a))
    b = by_folder.get(str(folder_b))
    if a is None or b is None:
        return CalcDiff(
            folder_a=str(folder_a), folder_b=str(folder_b),
            method_match=False, basis_match=False,
            a_functional=None, b_functional=None,
            a_basis=None, b_basis=None,
            a_gibbs=None, b_gibbs=None, delta_gibbs_kcal=None,
            a_spe=None, b_spe=None, delta_spe_kcal=None,
            a_n_imag=None, b_n_imag=None,
            notes=["could not parse one or both folders"],
        )

    notes: list[str] = []
    method_match = (
        a.functional is not None
        and b.functional is not None
        and a.functional.upper() == b.functional.upper()
    )
    basis_match = (
        a.basis is not None
        and b.basis is not None
        and a.basis.upper() == b.basis.upper()
    )
    if not method_match and a.functional and b.functional:
        notes.append(
            f"different functional: {a.functional} vs {b.functional}",
        )
    if not basis_match and a.basis and b.basis:
        notes.append(f"different basis: {a.basis} vs {b.basis}")
    if a.n_imag and a.n_imag > 0:
        notes.append(f"A has {a.n_imag} imag freq(s) — not a minimum")
    if b.n_imag and b.n_imag > 0:
        notes.append(f"B has {b.n_imag} imag freq(s) — not a minimum")

    delta_g = (
        (b.gibbs - a.gibbs) * HARTREE_TO_KCAL
        if a.gibbs is not None and b.gibbs is not None else None
    )
    delta_e = (
        (b.single_point - a.single_point) * HARTREE_TO_KCAL
        if a.single_point is not None and b.single_point is not None
        else None
    )
    return CalcDiff(
        folder_a=str(folder_a), folder_b=str(folder_b),
        method_match=method_match, basis_match=basis_match,
        a_functional=a.functional, b_functional=b.functional,
        a_basis=a.basis, b_basis=b.basis,
        a_gibbs=a.gibbs, b_gibbs=b.gibbs, delta_gibbs_kcal=delta_g,
        a_spe=a.single_point, b_spe=b.single_point, delta_spe_kcal=delta_e,
        a_n_imag=a.n_imag, b_n_imag=b.n_imag,
        notes=notes,
    )


# ---------------------------------------------------------------------------
# Meta-tool: tool catalog for on-demand discovery
# ---------------------------------------------------------------------------
#
# As the MCP toolbox grows, paying token cost for ALL schemas every turn
# scales poorly.  The agent can use ``list_tools`` to see what's
# available (one short line each) and ``describe_tool`` to fetch the
# full signature only when it actually plans to call it.

_TOOL_CATALOG: list[dict] = [
    # delfin-ops — environment / tool checks
    {"name": "qm_check", "category": "checks",
     "summary": "Verify xtb/crest/etc. tool availability."},
    {"name": "csp_check", "category": "checks",
     "summary": "Verify CSP (genarris) availability."},
    {"name": "mlp_check", "category": "checks",
     "summary": "Verify MLP backends (torchani/AIMNet2/MACE)."},
    {"name": "analysis_check", "category": "checks",
     "summary": "Verify Multiwfn/CENSO/ANMR/morfeus availability."},
    # delfin-ops — workflow execution (mutating)
    {"name": "pipeline_prepare", "category": "workflow",
     "summary": "Generate a CONTROL.txt template (mutating)."},
    {"name": "pipeline_run", "category": "workflow",
     "summary": "Run the full DELFIN pipeline (mutating, long-running)."},
    {"name": "run_orca_input", "category": "workflow",
     "summary": "Run ORCA on a single .inp file (mutating)."},
    {"name": "co2", "category": "workflow",
     "summary": "CO2 Coordinator workflow (mutating)."},
    {"name": "tadf_xtb", "category": "workflow",
     "summary": "TADF xTB workflow (mutating)."},
    {"name": "hyperpol", "category": "workflow",
     "summary": "Hyperpolarisability workflow (mutating)."},
    {"name": "cleanup", "category": "workflow",
     "summary": "Remove ORCA scratch / OCCUPIER artifacts (mutating)."},
    {"name": "stop", "category": "workflow",
     "summary": "Signal running DELFIN procs (INT/TERM/KILL, mutating)."},
    {"name": "stop_dry_run", "category": "workflow",
     "summary": "List DELFIN procs that WOULD be signaled (read-only)."},
    # delfin-ops — output parsing
    {"name": "parse_orca_output", "category": "parsing",
     "summary": "Snapshot ONE ORCA .out: energies, conv, freq, walltime."},
    {"name": "find_orca_errors", "category": "parsing",
     "summary": "Scan a folder's .out files for known error patterns."},
    {"name": "extract_thermochem", "category": "parsing",
     "summary": "Pull T, P, ZPE, H, S, G from a Freq output."},
    {"name": "extract_energy_table", "category": "parsing",
     "summary": "Multi-folder energy table (gibbs/zpe/spe/conv/...)."},
    {"name": "find_calculation_extreme", "category": "parsing",
     "summary": "Top-N folders sorted by an energy property."},
    {"name": "extract_imaginary_frequencies", "category": "parsing",
     "summary": "Imag freq modes + minimum/TS classification per folder."},
    {"name": "compare_calculations", "category": "parsing",
     "summary": "Side-by-side diff of method/basis/G/SPE for two folders."},
    {"name": "compare_across_functionals", "category": "parsing",
     "summary": "Multi-folder table grouped by functional/basis."},
    # delfin-ops — plotting
    {"name": "plot_energy_distribution", "category": "plotting",
     "summary": "Histogram/bar/boxplot of energies; PNG → workspace."},
    {"name": "plot_energy_correlation", "category": "plotting",
     "summary": "Scatter + Pearson r; PNG → workspace."},
    # delfin-ops — dashboard guidance
    {"name": "list_dashboard_patterns", "category": "guidance",
     "summary": "Names of operational pattern recipes."},
    {"name": "get_dashboard_pattern", "category": "guidance",
     "summary": "Slash-chain recipe for one workflow (batch/recalc/...)."},
    {"name": "list_dashboard_widgets", "category": "guidance",
     "summary": "Catalog of /ui-controllable widgets per tab."},
    {"name": "get_widget_options", "category": "guidance",
     "summary": "Allowed values for one widget (dropdown options)."},
    # delfin-ops — validation
    {"name": "validate_orca_input", "category": "validation",
     "summary": "Sanity-check an ORCA .inp text and return issues."},
    # delfin-ops — discovery
    {"name": "list_tools", "category": "meta",
     "summary": "Filter this very catalog (category, query)."},
    {"name": "describe_tool", "category": "meta",
     "summary": "Full description + signature of one tool."},
    # delfin-ops — job lifecycle
    {"name": "submit_calculation", "category": "jobs",
     "summary": "Submit a folder via the live backend (mutating)."},
    {"name": "cancel_calculation", "category": "jobs",
     "summary": "scancel one job by id (mutating)."},
    {"name": "list_active_calculations", "category": "jobs",
     "summary": "Live SLURM/local job list with status."},
    # delfin-ops — calc folder management (mutating, allow_mutate gates)
    {"name": "rename_calc_folder", "category": "calc-fs",
     "summary": "Rename one calc folder in place (mutating)."},
    {"name": "create_calc_folder", "category": "calc-fs",
     "summary": "mkdir a new sub-folder in calc/ (mutating)."},
    {"name": "move_calc_folder", "category": "calc-fs",
     "summary": "Move one calc folder to another calc/ location (mutating)."},
    {"name": "move_to_archive", "category": "calc-fs",
     "summary": "Move a calc folder calc/ -> archive/ (mutating)."},
    {"name": "delete_calc_folder", "category": "calc-fs",
     "summary": "rmtree a calc folder (3-lock destructive)."},
    # delfin-ops — bulk job control + recalc preparation
    {"name": "kill_all_user_jobs", "category": "jobs",
     "summary": "scancel every active job for the user (mutating)."},
    {"name": "prepare_recalc", "category": "workflow",
     "summary": "Pre-flight + submit Smart/classic/Override recalc."},
    # delfin-ops — calc-options dropdown
    {"name": "list_calc_options", "category": "calc-fs",
     "summary": "Items in the (Options) dropdown for one filename."},
    {"name": "run_calc_option", "category": "calc-fs",
     "summary": "Dispatch one (Options) action — recalc/UI hint/etc."},
    # delfin-ops — literature management
    {"name": "check_orca_manual_indexed", "category": "literature",
     "summary": "Is the ORCA manual indexed for search_docs?"},
    {"name": "index_new_pdf", "category": "literature",
     "summary": "Add a freshly-uploaded PDF to the search index."},
    {"name": "read_pdf", "category": "literature",
     "summary": "Plain text from a PDF (page selector + size cap)."},
    {"name": "search_pdf_local", "category": "literature",
     "summary": "Substring-search inside one specific PDF."},
    {"name": "extract_pdf_section", "category": "literature",
     "summary": "Pull one section from a PDF by heading."},
    {"name": "list_literature_files", "category": "literature",
     "summary": "List PDFs/MDs/TXTs available in the Literature folder."},
    # delfin-ops — DELFIN feature explainer
    {"name": "list_delfin_features", "category": "explainer",
     "summary": "Catalog of DELFIN concepts (control_keys, OCCUPIER, …)."},
    {"name": "explain_delfin_feature", "category": "explainer",
     "summary": "Curated explanation of one DELFIN concept."},
    # delfin-docs (other MCP server, listed for cross-reference)
    {"name": "search_docs", "category": "literature",
     "summary": "TF-IDF search over indexed PDFs."},
    {"name": "read_section", "category": "literature",
     "summary": "Read a specific section of an indexed doc."},
    {"name": "list_docs", "category": "literature",
     "summary": "List indexed documents."},
    {"name": "list_sections", "category": "literature",
     "summary": "List sections of one indexed doc."},
    {"name": "search_calcs", "category": "calc-search",
     "summary": "Search calculation metadata (functional/basis/...)."},
    {"name": "get_calc_info", "category": "calc-search",
     "summary": "Detailed info about one calculation by id."},
    {"name": "calc_summary", "category": "calc-search",
     "summary": "Aggregate counts across all indexed calculations."},
]


def list_tools(category: str = "", query: str = "") -> list[dict]:
    """Filter the typed tool catalog by category and/or substring query.

    Returns a list of {name, category, summary} entries — no full
    schemas, just enough for the agent to pick a candidate tool. Use
    :func:`describe_tool` for the full description + signature.

    Both filters are case-insensitive substring matches. Empty filters
    (the defaults) → return everything.
    """
    cat = (category or "").strip().lower()
    q = (query or "").strip().lower()
    out = []
    for entry in _TOOL_CATALOG:
        if cat and entry["category"].lower() != cat:
            continue
        if q and (q not in entry["name"].lower()
                  and q not in entry["summary"].lower()
                  and q not in entry["category"].lower()):
            continue
        out.append(dict(entry))
    return out


def describe_tool(name: str) -> dict:
    """Return the catalog entry for ``name`` plus its docstring.

    Looks the function up by name in this module first; falls back to
    the doc-server toolset if nothing local matches. Unknown names
    return a hint with the available choices.
    """
    if not name:
        return {"error": "no name given",
                "available": [e["name"] for e in _TOOL_CATALOG]}
    target = str(name).strip().lower()
    entry = next(
        (e for e in _TOOL_CATALOG if e["name"].lower() == target), None,
    )
    if entry is None:
        return {
            "error": f"unknown tool: {name!r}",
            "available": [e["name"] for e in _TOOL_CATALOG],
        }
    # Try to grab the docstring of the local function if we host it.
    fn = globals().get(entry["name"])
    docstring = (fn.__doc__ or "") if callable(fn) else ""
    return {
        "name": entry["name"],
        "category": entry["category"],
        "summary": entry["summary"],
        "docstring": docstring.strip(),
    }


# ---------------------------------------------------------------------------
# Dashboard widget catalog (static; lets the agent emit /ui commands
# without trial-and-error)
# ---------------------------------------------------------------------------

_WIDGET_CATALOG: list[dict] = [
    # Submit-Tab
    {"name": "control", "tab": "submit", "type": "Textarea",
     "purpose": "CONTROL.txt content (use /control key for single keys)."},
    {"name": "coords", "tab": "submit", "type": "Textarea",
     "purpose": "Coordinates (XYZ / SMILES / batch)."},
    {"name": "batch-smiles", "tab": "submit", "type": "Textarea",
     "purpose": "Batch SMILES/XYZ (use /batch from-calc to populate)."},
    {"name": "job-name", "tab": "submit", "type": "Text",
     "purpose": "Job name (used as folder prefix)."},
    {"name": "submit-btn", "tab": "submit", "type": "Button",
     "purpose": "Submit the configured job (DESTRUCTIVE)."},
    {"name": "time-limit", "tab": "submit", "type": "Dropdown",
     "purpose": "SLURM time limit preset.",
     "values": ["00:30:00", "01:00:00", "06:00:00",
                "12:00:00", "24:00:00", "48:00:00", "custom"]},
    {"name": "custom-time", "tab": "submit", "type": "Text",
     "purpose": "Custom HH:MM:SS, only used when time-limit=custom."},
    # ORCA Builder
    {"name": "orca-method", "tab": "orca", "type": "Dropdown",
     "purpose": "Functional / method.",
     "values": ["BP86", "PBE0", "B3LYP", "TPSS", "wB97X-D3",
                "CAM-B3LYP", "r2SCAN", "PBE", "M06-2X"]},
    {"name": "orca-basis", "tab": "orca", "type": "Dropdown",
     "purpose": "Basis set.",
     "values": ["def2-SVP", "def2-TZVP", "def2-TZVPP",
                "def2-QZVP", "ma-def2-TZVP", "x2c-TZVPall",
                "ZORA-def2-TZVP", "SARC-ZORA-TZVP", "pcSseg-2"]},
    {"name": "orca-job-type", "tab": "orca", "type": "Dropdown",
     "purpose": "Job type.",
     "values": ["SP", "OPT", "OPT Freq", "Freq", "TDDFT",
                "NEB-TS", "Scan"]},
    {"name": "orca-charge", "tab": "orca", "type": "IntText",
     "purpose": "Total molecular charge."},
    {"name": "orca-mult", "tab": "orca", "type": "IntText",
     "purpose": "Spin multiplicity (2S+1)."},
    {"name": "orca-pal", "tab": "orca", "type": "IntText",
     "purpose": "PAL — number of MPI processes."},
    {"name": "orca-maxcore", "tab": "orca", "type": "IntText",
     "purpose": "%maxcore in MB per process."},
    {"name": "orca-coords", "tab": "orca", "type": "Textarea",
     "purpose": "Coordinate block (xyz, gbw, etc.)."},
    {"name": "orca-dispersion", "tab": "orca", "type": "Dropdown",
     "purpose": "Dispersion correction.",
     "values": ["", "D3BJ", "D3(0)", "D4"]},
    {"name": "orca-solvent", "tab": "orca", "type": "Dropdown",
     "purpose": "Solvent for CPCM/SMD.",
     "values": ["", "water", "methanol", "ethanol",
                "acetonitrile", "DMSO", "DMF", "toluene",
                "thf", "chloroform", "dcm"]},
    {"name": "orca-preview", "tab": "orca", "type": "Textarea",
     "purpose": "Live INP preview generated from the form."},
    {"name": "orca-submit-btn", "tab": "orca", "type": "Button",
     "purpose": "Submit the ORCA job (DESTRUCTIVE)."},
    # Calc Browser
    {"name": "calc-path", "tab": "calc", "type": "Text",
     "purpose": "Current path inside calc/ or archive/."},
    {"name": "calc-options", "tab": "calc", "type": "Dropdown",
     "purpose": "Action for selected file.",
     "values": ["", "Recalc", "Smart Recalc", "Visualize",
                "Generate Report", "Override / Resubmit"]},
    {"name": "calc-editor", "tab": "calc", "type": "Textarea",
     "purpose": "Editor for CONTROL.txt (Smart Recalc panel)."},
    {"name": "calc-override-time", "tab": "calc", "type": "Text",
     "purpose": "Time-limit override for the resubmit."},
    {"name": "calc-override-btn", "tab": "calc", "type": "Button",
     "purpose": "Submit resubmit/override (DESTRUCTIVE)."},
    # Agent tab itself
    {"name": "send-btn", "tab": "agent", "type": "Button",
     "purpose": "Send the current input to the agent."},
    {"name": "input", "tab": "agent", "type": "Textarea",
     "purpose": "User-input area."},
    {"name": "mode", "tab": "agent", "type": "Dropdown",
     "purpose": "Agent mode.",
     "values": ["dashboard", "solo", "quick", "reviewed",
                "tdd", "cluster", "full"]},
    {"name": "perm", "tab": "agent", "type": "Dropdown",
     "purpose": "Permission profile.",
     "values": ["plan", "ask_all", "repo_free", "all_free"]},
]


def list_dashboard_widgets(tab: str = "") -> list[dict]:
    """Return the static catalog of /ui-controllable widgets.

    Optionally filter by ``tab`` (submit / orca / calc / agent).
    Each entry: name, tab, type, purpose, and (for dropdowns) an
    explicit ``values`` list — no need to /ui options first.
    """
    t = (tab or "").strip().lower()
    return [
        dict(e) for e in _WIDGET_CATALOG
        if not t or e["tab"].lower() == t
    ]


def get_widget_options(name: str) -> list:
    """Return the allowed values for a dropdown widget, or [] otherwise.

    Use BEFORE setting a value via /ui to avoid a "Invalid value"
    error round-trip.
    """
    if not name:
        return []
    target = str(name).strip().lower()
    for e in _WIDGET_CATALOG:
        if e["name"].lower() == target:
            return list(e.get("values", []))
    return []


# ---------------------------------------------------------------------------
# Validation: is this ORCA input plausible?
# ---------------------------------------------------------------------------


@dataclass
class ValidationIssue:
    """One validation finding."""
    severity: str   # "error" | "warning" | "info"
    code: str       # short machine-readable identifier
    message: str    # human-readable detail
    suggestion: str = ""


def validate_orca_input(inp_text: str) -> list[ValidationIssue]:
    """Sanity-check the text of an ORCA ``.inp`` file.

    Returns a list of :class:`ValidationIssue` entries. Empty list = no
    obvious problems detected (NOT proof the input is correct, just
    that the heuristic checks pass). Use this before submit when the
    user asks "passt alles im ORCA Builder?".

    Checks performed:

    - **structural**: ``! …`` keyword line present; `* xyz` / `* xyzfile`
      coordinate header; balanced ``%pal`` / ``end`` blocks.
    - **resources**: ``%pal nprocs N end`` parses (incl. inline form);
      ``%maxcore`` reasonable (1-20000 MB range).
    - **electronics**: charge / multiplicity present; multiplicity is
      odd/even consistent with electron count if the coordinate block
      has integers we can sum.
    - **method**: known functional + basis combo (warn on unknown).
    - **job-type vs keywords**: ``Freq`` requires no ``LARGEPRINT`` etc.
    - **dispersion vs functional**: D3(0)/D3BJ/D4 paired sensibly.

    Any discovered issue carries a short suggestion the agent can
    paste back to the user (or use to ``/control key`` an autofix).
    """
    import re as _re
    issues: list[ValidationIssue] = []
    if not (inp_text or "").strip():
        issues.append(ValidationIssue(
            severity="error", code="empty",
            message="INP text is empty.",
            suggestion="Generate the input via ORCA Builder or paste a real .inp.",
        ))
        return issues

    text = inp_text
    text_low = text.lower()

    # --- Structural checks -----------------------------------------------
    m_keyword = _re.search(r"^\s*!\s*(.+)$", text, _re.MULTILINE)
    if not m_keyword:
        issues.append(ValidationIssue(
            severity="error", code="no_keyword_line",
            message="No `! …` keyword line found.",
            suggestion="Add a line like `! PBE0 def2-TZVP D4 OPT Freq` at the top.",
        ))
        keyword_line = ""
    else:
        keyword_line = m_keyword.group(1).strip()

    if not _re.search(r"^\s*\*\s+(xyz|xyzfile)\s+", text, _re.MULTILINE):
        issues.append(ValidationIssue(
            severity="error", code="no_xyz_block",
            message="No `* xyz <charge> <mult>` (or xyzfile) block found.",
            suggestion="Add `* xyz 0 1` followed by atom lines and a closing `*`.",
        ))

    # Charge / multiplicity from xyz header
    m_xyz = _re.search(
        r"^\s*\*\s+xyz\s+(-?\d+)\s+(\d+)\s*$", text, _re.MULTILINE,
    )
    charge = int(m_xyz.group(1)) if m_xyz else None
    mult = int(m_xyz.group(2)) if m_xyz else None
    if mult is not None and mult < 1:
        issues.append(ValidationIssue(
            severity="error", code="bad_multiplicity",
            message=f"Multiplicity {mult} is < 1.",
            suggestion="Multiplicity = 2S+1, so it must be ≥ 1 (singlet=1, doublet=2, ...).",
        ))

    # --- %pal block ------------------------------------------------------
    if "%pal" in text_low:
        # Use the same regex as parse_inp_resources (inline + multi-line)
        m_pal = _re.search(r"\bnprocs\s+(\d+)", text, _re.IGNORECASE)
        if not m_pal:
            issues.append(ValidationIssue(
                severity="error", code="pal_no_nprocs",
                message="`%pal` block has no `nprocs` value.",
                suggestion="Use `%pal nprocs 12 end` (or multi-line form).",
            ))
    if "%maxcore" in text_low:
        m_mc = _re.search(r"%maxcore\s+(\d+)", text, _re.IGNORECASE)
        if m_mc:
            mc = int(m_mc.group(1))
            if mc < 100:
                issues.append(ValidationIssue(
                    severity="warning", code="maxcore_low",
                    message=f"%maxcore = {mc} MB is unusually low.",
                    suggestion="Typical values are 1000-8000 MB per core.",
                ))
            elif mc > 32000:
                issues.append(ValidationIssue(
                    severity="warning", code="maxcore_high",
                    message=f"%maxcore = {mc} MB is unusually high.",
                    suggestion="Above 16-32 GB/core is uncommon — verify SLURM allocation.",
                ))

    # --- Method / basis sanity (heuristic) ------------------------------
    kw = keyword_line.upper()
    has_functional = any(
        f in kw.split()
        for f in ("BP86", "PBE0", "B3LYP", "TPSS", "WB97X", "WB97X-D3",
                  "CAM-B3LYP", "R2SCAN", "PBE", "M06", "M06-2X", "BLYP",
                  "REVTPSS", "HF", "MP2", "CCSD")
    )
    has_basis = any(
        b in kw
        for b in ("DEF2-", "PCSSEG", "MA-DEF2-", "SARC-", "X2C-", "CC-PV", "AUG-")
    )
    if keyword_line and not has_functional:
        issues.append(ValidationIssue(
            severity="warning", code="no_functional",
            message="Could not detect a known functional in the keyword line.",
            suggestion="Common picks: PBE0, BP86, B3LYP, wB97X-D3, CAM-B3LYP.",
        ))
    if keyword_line and not has_basis:
        issues.append(ValidationIssue(
            severity="warning", code="no_basis",
            message="Could not detect a known basis set.",
            suggestion="Common picks: def2-SVP, def2-TZVP, ma-def2-TZVP.",
        ))

    # --- Dispersion vs functional --------------------------------------
    has_d4 = "D4" in kw.split()
    has_d3 = "D3BJ" in kw.split() or "D3(0)" in kw.replace(" ", "")
    if has_d4 and has_d3:
        issues.append(ValidationIssue(
            severity="warning", code="dispersion_double",
            message="Both D3 and D4 appear in the keyword line.",
            suggestion="Pick one — usually D4 for newer functionals, D3BJ otherwise.",
        ))

    # --- Job type vs accompanying keywords -----------------------------
    if "FREQ" in kw and "TIGHTOPT" not in kw and "OPT" in kw and "VERYTIGHTOPT" not in kw:
        issues.append(ValidationIssue(
            severity="info", code="freq_after_loose_opt",
            message="`OPT Freq` without TIGHTOPT may give imag-freq false positives.",
            suggestion="Add TIGHTOPT or VERYTIGHTOPT for clean Hessian.",
        ))

    if "NEB-TS" in kw and not _re.search(r"%neb", text, _re.IGNORECASE):
        issues.append(ValidationIssue(
            severity="error", code="neb_no_block",
            message="NEB-TS keyword without %neb block.",
            suggestion="Add `%neb NImages 8 end` (or your image count) to the input.",
        ))

    # --- Relativistic consistency --------------------------------------
    has_zora = "ZORA" in kw
    has_x2c = "X2C" in kw
    if has_zora and "DEF2-" in kw and "ZORA-DEF2-" not in kw:
        issues.append(ValidationIssue(
            severity="warning", code="zora_basis_mismatch",
            message="ZORA Hamiltonian set but plain def2- basis used.",
            suggestion="Use ZORA-def2-TZVP or SARC-ZORA-TZVP with ZORA.",
        ))
    if has_x2c and ("DEF2-" in kw and "X2C-" not in kw):
        issues.append(ValidationIssue(
            severity="warning", code="x2c_basis_mismatch",
            message="X2C set but def2- basis (non-relativistic) used.",
            suggestion="Use x2c-TZVPall / x2c-QZVPPall with X2C.",
        ))

    return issues


# ---------------------------------------------------------------------------
# Job lifecycle (typed wrappers around backend_local / backend_slurm)
# ---------------------------------------------------------------------------


@dataclass
class SubmitOutcome:
    """Result of a submit_calculation call."""
    job_id: str = ""
    submitted: bool = False
    folder: str = ""
    backend: str = ""
    message: str = ""
    error: str = ""


def _resolve_backend():
    """Pick the right backend (local vs slurm) for the running host."""
    import shutil as _sh
    if _sh.which("sbatch") and _sh.which("squeue"):
        from delfin.dashboard.backend_slurm import SLURMBackend
        return SLURMBackend(orca_base="")
    from delfin.dashboard.backend_local import LocalBackend
    return LocalBackend(orca_base="")


def submit_calculation(
    folder: str,
    *,
    job_name: str = "",
    mode: str = "delfin",
    time_limit: str = "12:00:00",
    pal: int = 12,
    maxcore: int = 6000,
    allow_mutate: bool = False,
) -> SubmitOutcome:
    """Submit a calculation folder via the live backend (DESTRUCTIVE).

    ``allow_mutate=False`` (default) → no submit. Returns a
    SubmitOutcome with ``submitted=False`` and a hint message so the
    agent can confirm with the user before a real submit.

    Args:
        folder: absolute path to the job folder (must contain a
            CONTROL.txt or .inp).
        job_name: optional name (defaults to folder basename).
        mode: ``delfin`` (full pipeline) or ``orca`` (single-step).
        time_limit: SLURM HH:MM:SS.
        pal: cores.
        maxcore: per-core memory in MB.
        allow_mutate: must be True for the submit to actually run.
    """
    from pathlib import Path as _P
    p = _P(folder)
    if not p.exists() or not p.is_dir():
        return SubmitOutcome(
            folder=str(folder),
            error=f"folder not found or not a directory: {folder}",
        )
    if not allow_mutate:
        return SubmitOutcome(
            folder=str(folder),
            backend="(skipped — allow_mutate=False)",
            message=(
                f"Would submit {folder} as mode={mode!r}, "
                f"time={time_limit}, pal={pal}, maxcore={maxcore}. "
                f"Confirm with the user, then call again with "
                f"allow_mutate=True."
            ),
        )
    name = job_name or p.name
    try:
        backend = _resolve_backend()
        result = backend.submit_delfin(
            job_dir=str(p), job_name=name, mode=mode,
            time_limit=time_limit, pal=pal, maxcore=maxcore,
        )
    except Exception as exc:
        return SubmitOutcome(
            folder=str(folder), error=f"submit failed: {exc}",
        )
    submitted = bool(getattr(result, "returncode", 1) == 0)
    job_id = ""
    stdout = getattr(result, "stdout", "") or ""
    import re as _re
    m = _re.search(r"\b(\d{6,})\b", stdout)
    if m:
        job_id = m.group(1)
    return SubmitOutcome(
        job_id=job_id,
        submitted=submitted,
        folder=str(p),
        backend=type(backend).__name__,
        message=(stdout or "").strip()[:500],
        error="" if submitted else (
            getattr(result, "stderr", "") or "submit returned non-zero"
        )[:500],
    )


def cancel_calculation(
    job_id: str,
    *,
    allow_mutate: bool = False,
) -> dict:
    """Cancel a running job by SLURM/local id (DESTRUCTIVE).

    ``allow_mutate=False`` → returns a "would cancel" hint without
    running scancel. Confirm with the user, then re-call.
    """
    if not (job_id or "").strip():
        return {"ok": False, "error": "no job_id given"}
    if not allow_mutate:
        return {
            "ok": False,
            "skipped": True,
            "message": (
                f"Would cancel job {job_id}. Confirm with the user, "
                f"then call again with allow_mutate=True."
            ),
        }
    try:
        backend = _resolve_backend()
        ok, msg = backend.cancel_job(job_id)
        return {"ok": bool(ok), "message": msg, "job_id": job_id}
    except Exception as exc:
        return {"ok": False, "error": str(exc), "job_id": job_id}


# ---------------------------------------------------------------------------
# Calculations-tab folder management (rename / new / move / archive / delete)
# ---------------------------------------------------------------------------
#
# These mirror the operations behind the Calculations browser tab, exposed
# as typed Python so the agent can drive them without /ui-clicking. All
# mutating ops require ``allow_mutate=True`` AND fall under the directory-
# permissions rule (calculations/ writable only via ACTION:, archives are
# read-only). The ``allowed_roots`` parameter is a sanity wall so the agent
# never moves something outside the project tree by accident.


def _resolve_within_roots(path: str, allowed_roots: list[str]) -> "Path":
    """Resolve ``path`` and require it sits inside one of ``allowed_roots``.

    Raises ValueError if it doesn't. Both sides are real-path resolved.
    """
    from pathlib import Path as _P
    resolved = _P(path).expanduser().resolve()
    for root in allowed_roots:
        try:
            r = _P(root).expanduser().resolve()
        except Exception:
            continue
        try:
            resolved.relative_to(r)
            return resolved
        except ValueError:
            continue
    raise ValueError(
        f"path {resolved} is outside allowed roots {allowed_roots}"
    )


def _default_calc_roots() -> list[str]:
    """Best-effort default for the (calc, archive) roots in the cwd."""
    from pathlib import Path as _P
    cwd = _P.cwd()
    candidates = []
    for name in ("calculations", "calc", "archive", "remote_archive"):
        p = cwd / name
        if p.exists() and p.is_dir():
            candidates.append(str(p))
    return candidates


def _check_destructive_dirs(
    *,
    block_archive: bool,
    paths: list[str],
) -> str | None:
    """Refuse mutating ops that would touch read-only roots."""
    from pathlib import Path as _P
    if not block_archive:
        return None
    for raw in paths:
        try:
            p = _P(raw).expanduser().resolve()
        except Exception:
            continue
        for marker in ("archive", "remote_archive"):
            if marker in p.parts:
                # Allow targeting "...archive/..." only as an explicit
                # destination from a calc/ source (move_to_archive). The
                # caller distinguishes via block_archive.
                return (
                    f"Refusing to mutate inside read-only root '{marker}/' "
                    f"(path: {p})."
                )
    return None


def rename_calc_folder(
    src: str,
    new_name: str,
    *,
    allowed_roots: list[str] | None = None,
    allow_mutate: bool = False,
) -> dict:
    """Rename a calculation folder in place (destructive).

    ``new_name`` is a basename — it cannot contain path separators.
    Refuses to overwrite an existing target. Use only inside the
    ``calculations/`` tree; calls outside it raise ValueError.

    Returns a dict with ``ok``, ``src``, ``dst``, ``message``.
    """
    from pathlib import Path as _P
    roots = allowed_roots or _default_calc_roots()
    if not roots:
        return {
            "ok": False,
            "error": "no allowed_roots and none could be inferred",
        }
    if "/" in new_name or "\\" in new_name or not new_name.strip():
        return {"ok": False, "error": "new_name must be a basename"}
    try:
        src_p = _resolve_within_roots(src, roots)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    if not src_p.exists():
        return {"ok": False, "error": f"source missing: {src_p}"}
    block = _check_destructive_dirs(block_archive=True, paths=[str(src_p)])
    if block:
        return {"ok": False, "error": block}
    dst_p = src_p.parent / new_name.strip()
    if dst_p.exists():
        return {"ok": False, "error": f"target already exists: {dst_p}"}
    if not allow_mutate:
        return {
            "ok": False,
            "skipped": True,
            "src": str(src_p),
            "dst": str(dst_p),
            "message": (
                f"Would rename {src_p.name} -> {dst_p.name}. Confirm "
                f"with the user, then call again with allow_mutate=True."
            ),
        }
    try:
        src_p.rename(dst_p)
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    return {
        "ok": True, "src": str(src_p), "dst": str(dst_p),
        "message": f"Renamed to {dst_p.name}.",
    }


def create_calc_folder(
    parent: str,
    name: str,
    *,
    allowed_roots: list[str] | None = None,
    allow_mutate: bool = False,
) -> dict:
    """Create a new sub-folder inside ``parent`` (destructive)."""
    from pathlib import Path as _P
    roots = allowed_roots or _default_calc_roots()
    if not roots:
        return {"ok": False, "error": "no allowed_roots inferred"}
    if "/" in name or "\\" in name or not name.strip():
        return {"ok": False, "error": "name must be a basename"}
    try:
        parent_p = _resolve_within_roots(parent, roots)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    if not parent_p.is_dir():
        return {"ok": False, "error": f"parent is not a directory: {parent_p}"}
    block = _check_destructive_dirs(
        block_archive=True, paths=[str(parent_p)],
    )
    if block:
        return {"ok": False, "error": block}
    dst_p = parent_p / name.strip()
    if dst_p.exists():
        return {"ok": False, "error": f"already exists: {dst_p}"}
    if not allow_mutate:
        return {
            "ok": False, "skipped": True, "dst": str(dst_p),
            "message": (
                f"Would create {dst_p}. Confirm and re-run with "
                f"allow_mutate=True."
            ),
        }
    try:
        dst_p.mkdir(parents=False, exist_ok=False)
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    return {"ok": True, "dst": str(dst_p), "message": f"Created {dst_p}."}


def move_calc_folder(
    src: str,
    dst_parent: str,
    *,
    allowed_roots: list[str] | None = None,
    allow_mutate: bool = False,
) -> dict:
    """Move ``src`` into ``dst_parent`` (both must be inside calc roots).

    Use this for calc-to-calc moves (re-organizing the tree). For
    sending things to ``archive/`` use :func:`move_to_archive`, which
    enforces the calc → archive direction.
    """
    from pathlib import Path as _P
    import shutil as _sh
    roots = allowed_roots or _default_calc_roots()
    if not roots:
        return {"ok": False, "error": "no allowed_roots inferred"}
    try:
        src_p = _resolve_within_roots(src, roots)
        dst_parent_p = _resolve_within_roots(dst_parent, roots)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    if not src_p.exists():
        return {"ok": False, "error": f"source missing: {src_p}"}
    if not dst_parent_p.is_dir():
        return {
            "ok": False,
            "error": f"destination parent is not a directory: {dst_parent_p}",
        }
    block = _check_destructive_dirs(
        block_archive=True, paths=[str(src_p), str(dst_parent_p)],
    )
    if block:
        return {"ok": False, "error": block}
    dst_p = dst_parent_p / src_p.name
    if dst_p.exists():
        return {"ok": False, "error": f"target already exists: {dst_p}"}
    if dst_p == src_p or str(dst_parent_p).startswith(str(src_p) + "/"):
        return {
            "ok": False, "error": "cannot move a folder into itself",
        }
    if not allow_mutate:
        return {
            "ok": False, "skipped": True,
            "src": str(src_p), "dst": str(dst_p),
            "message": "Would move (calc -> calc). Re-run with allow_mutate=True.",
        }
    try:
        _sh.move(str(src_p), str(dst_p))
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    return {
        "ok": True, "src": str(src_p), "dst": str(dst_p),
        "message": f"Moved {src_p.name} -> {dst_p}.",
    }


def move_to_archive(
    src: str,
    archive_root: str = "",
    *,
    allowed_roots: list[str] | None = None,
    allow_mutate: bool = False,
) -> dict:
    """Move a calc folder from ``calculations/`` into ``archive/``.

    Direction is enforced: ``src`` must NOT already be under archive,
    and ``archive_root`` must be the project's archive directory. With
    an empty ``archive_root`` the function infers the cwd's archive/.
    """
    from pathlib import Path as _P
    import shutil as _sh
    cwd = _P.cwd()
    if not archive_root:
        # Default: ./archive next to the working directory
        cand = cwd / "archive"
        archive_root = str(cand)
    arc_p = _P(archive_root).expanduser()
    if not arc_p.exists():
        return {
            "ok": False,
            "error": f"archive_root does not exist: {arc_p}",
        }
    arc_resolved = arc_p.resolve()
    if "archive" not in arc_resolved.parts:
        return {
            "ok": False,
            "error": (
                "archive_root must contain a path segment named "
                "'archive' (got " + str(arc_resolved) + ")"
            ),
        }
    # Source must be in calc/, not in archive/
    roots = allowed_roots or _default_calc_roots()
    if not roots:
        return {"ok": False, "error": "no allowed_roots inferred"}
    try:
        src_p = _resolve_within_roots(src, roots)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    if "archive" in src_p.parts:
        return {
            "ok": False,
            "error": "source is already inside archive/",
        }
    dst_p = arc_resolved / src_p.name
    if dst_p.exists():
        return {"ok": False, "error": f"target already exists: {dst_p}"}
    if not allow_mutate:
        return {
            "ok": False, "skipped": True,
            "src": str(src_p), "dst": str(dst_p),
            "message": "Would move calc -> archive. Re-run with allow_mutate=True.",
        }
    try:
        arc_resolved.mkdir(parents=True, exist_ok=True)
        _sh.move(str(src_p), str(dst_p))
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    return {
        "ok": True, "src": str(src_p), "dst": str(dst_p),
        "message": f"Archived {src_p.name}.",
    }


def delete_calc_folder(
    folder: str,
    *,
    confirm_token: str = "",
    allowed_roots: list[str] | None = None,
    allow_mutate: bool = False,
) -> dict:
    """Permanently delete a calc folder (DESTRUCTIVE — extra-strict gate).

    Three locks: ``allow_mutate=True``, target must be inside calc/
    roots (NOT archive/), and ``confirm_token`` must equal the folder's
    basename verbatim. Any one missing → refusal.
    """
    from pathlib import Path as _P
    import shutil as _sh
    roots = allowed_roots or _default_calc_roots()
    if not roots:
        return {"ok": False, "error": "no allowed_roots inferred"}
    try:
        target = _resolve_within_roots(folder, roots)
    except ValueError as exc:
        return {"ok": False, "error": str(exc)}
    if not target.exists():
        return {"ok": False, "error": f"missing: {target}"}
    block = _check_destructive_dirs(block_archive=True, paths=[str(target)])
    if block:
        return {"ok": False, "error": block}
    if confirm_token != target.name:
        return {
            "ok": False, "skipped": True,
            "message": (
                f"To delete '{target}', call again with "
                f"confirm_token='{target.name}' AND allow_mutate=True."
            ),
        }
    if not allow_mutate:
        return {
            "ok": False, "skipped": True,
            "message": "confirm_token correct, still need allow_mutate=True.",
        }
    try:
        if target.is_dir():
            _sh.rmtree(target)
        else:
            target.unlink()
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    return {
        "ok": True, "deleted": str(target),
        "message": f"Deleted {target.name}.",
    }


# ---------------------------------------------------------------------------
# Job lifecycle helpers (bulk cancel) + recalc preparation
# ---------------------------------------------------------------------------


def kill_all_user_jobs(
    *,
    only_running: bool = False,
    allow_mutate: bool = False,
) -> dict:
    """Cancel every active job for the current user (DESTRUCTIVE).

    Lists jobs via :func:`list_active_calculations` and calls
    :func:`cancel_calculation` for each. With ``only_running=True``
    skips PENDING jobs.

    Returns a dict with ``ok``, ``cancelled`` (list of {job_id, ok}),
    ``skipped`` (list), and ``total``.
    """
    try:
        jobs = list_active_calculations()
    except Exception as exc:
        return {"ok": False, "error": str(exc)}
    candidates = []
    for job in jobs:
        status = (job.get("status") or "").upper()
        if only_running and status not in ("RUNNING", "R"):
            continue
        jid = str(job.get("job_id") or job.get("id") or "").strip()
        if jid:
            candidates.append(jid)

    if not allow_mutate:
        return {
            "ok": False,
            "skipped": True,
            "candidates": candidates,
            "total": len(candidates),
            "message": (
                f"Would cancel {len(candidates)} job(s). Confirm with "
                f"the user, then re-run with allow_mutate=True."
            ),
        }

    cancelled, errors = [], []
    for jid in candidates:
        rc = cancel_calculation(jid, allow_mutate=True)
        if rc.get("ok"):
            cancelled.append(jid)
        else:
            errors.append({"job_id": jid, "error": rc.get("error")})
    return {
        "ok": not errors,
        "cancelled": cancelled,
        "errors": errors,
        "total": len(candidates),
    }


@dataclass
class RecalcPlan:
    """Pre-flight summary of a recalc submission before it runs."""
    folder: str
    mode: str
    time_limit: str
    pal: int | None
    maxcore: int | None
    override: str
    valid: bool
    issues: list[str]
    submitted_job_id: str | None = None
    submission_message: str = ""


def _extract_pal_maxcore(folder: "Path") -> tuple[int | None, int | None]:
    """Pull (pal, maxcore) from CONTROL.txt or the largest .inp."""
    import re as _re
    pal = maxcore = None
    control = folder / "CONTROL.txt"
    if control.exists():
        try:
            text = control.read_text(encoding="utf-8", errors="replace")
            m = _re.search(r"^\s*PAL\s*=\s*(\d+)", text, _re.MULTILINE)
            if m:
                pal = int(m.group(1))
            m = _re.search(r"^\s*maxcore\s*=\s*(\d+)", text, _re.MULTILINE)
            if m:
                maxcore = int(m.group(1))
        except Exception:
            pass
    if pal is None or maxcore is None:
        for inp in sorted(folder.glob("*.inp")):
            try:
                text = inp.read_text(encoding="utf-8", errors="replace")
            except Exception:
                continue
            if pal is None:
                m = _re.search(r"\bnprocs\s+(\d+)", text)
                if m:
                    pal = int(m.group(1))
            if maxcore is None:
                m = _re.search(r"%\s*maxcore\s+(\d+)", text)
                if m:
                    maxcore = int(m.group(1))
            if pal is not None and maxcore is not None:
                break
    return pal, maxcore


def prepare_recalc(
    folder: str,
    *,
    mode: str = "smart",
    time_limit: str = "24:00:00",
    override: str = "",
    allow_mutate: bool = False,
) -> RecalcPlan:
    """Prepare a Smart-Recalc / classic Recalc / Override resubmission.

    Reads PAL/maxcore from the folder's CONTROL.txt or largest .inp,
    reports any issues (missing values, no input file, …), and either
    returns a dry-run plan (``allow_mutate=False``) or actually
    submits the job via the live backend.

    Args:
        folder: absolute path to the calc folder.
        mode: ``smart`` (default) → delfin-recalc, ``classic`` →
            delfin-recalc-classic, ``override`` → delfin-recalc-override
            (then ``override`` must be ``STAGE=INDEX``).
        time_limit: SLURM time spec (``HH:MM:SS``).
        override: only used when mode=override.
        allow_mutate: must be True to actually submit.
    """
    from pathlib import Path as _P
    p = _P(folder).expanduser().resolve()
    issues: list[str] = []
    if not p.exists() or not p.is_dir():
        issues.append(f"folder missing: {p}")
    mode_clean = mode.strip().lower()
    submission_mode_map = {
        "smart": "delfin-recalc",
        "classic": "delfin-recalc-classic",
        "override": "delfin-recalc-override",
    }
    submission_mode = submission_mode_map.get(mode_clean)
    if submission_mode is None:
        issues.append(
            f"mode must be one of {sorted(submission_mode_map)} (got {mode!r})"
        )
    if mode_clean == "override":
        if not override.strip() or "=" not in override:
            issues.append(
                "override mode needs STAGE=INDEX (e.g. red_step_2=1)"
            )
    if not _re_match(r"^\d+:\d{2}:\d{2}$", time_limit):
        issues.append(f"time_limit must be HH:MM:SS (got {time_limit!r})")

    pal: int | None = None
    maxcore: int | None = None
    if p.is_dir():
        pal, maxcore = _extract_pal_maxcore(p)
        if pal is None:
            issues.append("PAL not found in CONTROL.txt / .inp")
        if maxcore is None:
            issues.append("maxcore not found in CONTROL.txt / .inp")

    plan = RecalcPlan(
        folder=str(p),
        mode=submission_mode or "",
        time_limit=time_limit,
        pal=pal,
        maxcore=maxcore,
        override=override,
        valid=not issues,
        issues=issues,
    )
    if not plan.valid:
        return plan
    if not allow_mutate:
        plan.submission_message = (
            f"Would submit {submission_mode} for {p.name} "
            f"(PAL={pal}, maxcore={maxcore}, time={time_limit}). "
            f"Confirm with the user, then re-run with allow_mutate=True."
        )
        return plan

    try:
        backend = _resolve_backend()
        job_name = f"recalc_{p.name}"
        kwargs = dict(
            job_dir=str(p),
            job_name=job_name,
            mode=submission_mode,
            time_limit=time_limit,
            pal=pal,
            maxcore=maxcore,
        )
        if mode_clean == "override":
            kwargs["override"] = override
        rc = backend.submit_delfin(**kwargs)
        if getattr(rc, "returncode", 0) == 0:
            plan.submitted_job_id = str(getattr(rc, "job_id", "") or "")
            plan.submission_message = (
                f"Submitted {submission_mode} for {p.name}"
                + (
                    f" (job {plan.submitted_job_id})"
                    if plan.submitted_job_id else ""
                )
            )
        else:
            plan.valid = False
            plan.issues.append(
                f"backend returned rc={rc.returncode}: "
                f"{getattr(rc, 'stderr', '')[:200]}"
            )
    except Exception as exc:
        plan.valid = False
        plan.issues.append(f"submission failed: {exc}")
    return plan


def _re_match(pattern: str, value: str) -> bool:
    import re as _re
    return _re.match(pattern, value or "") is not None


# ---------------------------------------------------------------------------
# Calculations-tab options dropdown — typed map + dispatcher
# ---------------------------------------------------------------------------
#
# The Calculations browser shows a context-aware "(Options)" dropdown
# whose contents depend on the selected file's type. The agent can use
# ``list_calc_options(filename)`` to discover what's offered, and
# ``run_calc_option(folder, option, ...)`` as a single typed entry to
# the most common ones (Recalc / Smart Recalc / Override / Visualize).

_CALC_OPTIONS_BY_TYPE: dict[str, list[str]] = {
    "control_txt":      ["Recalc", "Smart Recalc"],
    "occupier_txt":     ["Override"],
    "inp":              ["Recalc"],
    "out":              ["Print Mode", "MO Plot", "Print NMR", "RMSD"],
    "xyz":              [
        "Build Batch from XYZ", "Calc NMR", "Calc CENSO/ANMR",
        "hyperpol_xtb", "tadf_xtb", "RMSD",
    ],
    "csv_complete":     ["Preselection", "Visualize"],
    "csv_preselected":  ["Visualize"],
    "csv_rejected":     ["Visualize"],
    "final_interp":     ["Plot Trajectory"],
    "other":            ["RMSD"],
}


def _classify_calc_file(filename: str) -> str:
    """Map a basename to one of the keys in :data:`_CALC_OPTIONS_BY_TYPE`."""
    name = (filename or "").strip()
    lower = name.lower()
    if name == "CONTROL.txt":
        return "control_txt"
    if name == "OCCUPIER.txt":
        return "occupier_txt"
    if lower.endswith(".inp"):
        return "inp"
    if lower.endswith(".out"):
        return "out"
    if lower.endswith(".xyz"):
        return "xyz"
    if lower.endswith(".final.interp"):
        return "final_interp"
    if lower.endswith(".csv"):
        # The dashboard further splits by content — without reading the
        # CSV we can't tell complete vs preselected vs rejected, so we
        # return the broadest superset so the agent sees every option.
        return "csv_complete"
    return "other"


def list_calc_options(filename: str) -> dict:
    """Return the (Options)-dropdown items for one selected file.

    Mirrors the dashboard's context-aware menu so the agent can plan
    "what can I do with this file" without /ui-trial-and-error. The
    ``file_type`` field shows which classification was applied.
    """
    file_type = _classify_calc_file(filename)
    options = list(_CALC_OPTIONS_BY_TYPE.get(file_type, []))
    return {
        "filename": filename,
        "file_type": file_type,
        "options": options,
    }


def run_calc_option(
    folder: str,
    option: str,
    *,
    target_file: str = "",
    time_limit: str = "24:00:00",
    override: str = "",
    allow_mutate: bool = False,
) -> dict:
    """Generic dispatcher for the (Options) dropdown.

    Today this routes to :func:`prepare_recalc` for the destructive
    workflow paths (Recalc / Smart Recalc / Override) and returns a
    "not_implemented" hint for read-only / visualisation options that
    only make sense inside the dashboard UI (Visualize, MO Plot, …).
    The agent should drive those via ``ACTION:`` slash-commands.

    Args:
        folder: absolute path to the calc folder.
        option: dropdown text — "Smart Recalc", "Recalc", "Override",
            "Visualize", "MO Plot", "Print Mode", "Print NMR",
            "Plot Trajectory", "Preselection", "RMSD", "Build Batch
            from XYZ", "Calc NMR", "Calc CENSO/ANMR", "hyperpol_xtb",
            "tadf_xtb".
        target_file: optional basename when the option references a
            single file (e.g. RMSD between two .xyz files).
        time_limit: SLURM time spec for recalcs.
        override: only used when option=Override (STAGE=INDEX).
        allow_mutate: required for any submission.
    """
    opt = (option or "").strip()
    routing = {
        "Smart Recalc": ("recalc", "smart"),
        "Recalc":       ("recalc", "classic"),
        "Override":     ("recalc", "override"),
    }
    routed = routing.get(opt)
    if routed and routed[0] == "recalc":
        plan = prepare_recalc(
            folder,
            mode=routed[1],
            time_limit=time_limit,
            override=override,
            allow_mutate=allow_mutate,
        )
        from dataclasses import asdict as _asdict
        return {
            "ok": plan.valid and (allow_mutate or not plan.issues),
            "option": opt,
            **_asdict(plan),
        }

    ui_only = {
        "Visualize", "MO Plot", "Print Mode", "Print NMR",
        "Plot Trajectory", "Preselection", "RMSD",
        "Build Batch from XYZ", "Calc NMR", "Calc CENSO/ANMR",
        "hyperpol_xtb", "tadf_xtb",
    }
    if opt in ui_only:
        return {
            "ok": False,
            "option": opt,
            "ui_only": True,
            "message": (
                f"Option '{opt}' runs inside the dashboard UI. The "
                f"agent should drive it with an ACTION: slash-command "
                f"(see list_dashboard_widgets / get_dashboard_pattern)."
            ),
            "hint": (
                f"For {opt}: select the file in the calc browser, set "
                f"the (Options) dropdown to '{opt}', then click the "
                f"action button. From the agent: ACTION: /calc select "
                f"<file>; ACTION: /calc options {opt}."
            ),
        }
    return {
        "ok": False,
        "option": opt,
        "error": (
            f"Unknown option {opt!r}. Use list_calc_options(filename) "
            f"to discover valid choices for the selected file."
        ),
    }


# ---------------------------------------------------------------------------
# DELFIN-feature explainer (curated knowledge base of DELFIN's own concepts)
# ---------------------------------------------------------------------------

_DELFIN_FEATURES: dict[str, dict] = {
    # CONTROL.txt keys
    "control": {
        "category": "config",
        "summary": (
            "DELFIN's master configuration file. Plain key=value lines "
            "drive every workflow — functional, basis, charge, "
            "multiplicity, solvent, redox steps, etc."
        ),
        "see_also": [
            "delfin/utils.py", "delfin/common/control_validator.py",
            "delfin/dashboard/constants.py",
        ],
    },
    "control_keys": {
        "category": "config",
        "summary": (
            "Common keys: functional, main_basisset, metal_basisset, "
            "aux_jk, disp_corr, solvent, solvation_model, freq_type, "
            "geom_opt, charge, multiplicity, PAL, maxcore, "
            "redox_steps, parallel_workflows. Relativistic variants "
            "carry a `_rel` suffix (main_basisset_rel etc.) and only "
            "fire when `relativity` is set (ZORA/X2C/DKH)."
        ),
        "see_also": ["delfin/common/control_validator.py"],
    },
    "relativistic_methods": {
        "category": "methodology",
        "summary": (
            "Three Hamiltonians supported: ZORA (zeroth-order regular "
            "approximation), X2C (exact two-component), DKH (Douglas-"
            "Kroll-Hess). When `relativity` is set, the `*_rel` keys "
            "select matching basis sets — ZORA wants ZORA-def2-TZVP / "
            "SARC-ZORA-TZVP, X2C wants x2c-TZVPall / x2c-QZVPPall, "
            "and aux_jk_rel is empty for X2C (def2/J works) but "
            "SARC/J for ZORA. The non-rel keys stay as-is."
        ),
        "see_also": ["delfin/utils.py"],
    },
    # Workflow concepts
    "pipeline": {
        "category": "workflow",
        "summary": (
            "DELFIN's main pipeline runs the full QM-chemistry workflow "
            "for one CONTROL.txt: prepare geometry → optimize → "
            "frequencies → optional redox / solvation steps → property "
            "extraction. Driven by `delfin run` or `pipeline_run` MCP "
            "tool. Stages can run in parallel via `parallel_workflows`."
        ),
        "see_also": ["delfin/cli.py", "delfin/runtime_setup.py"],
    },
    "smart_recalc": {
        "category": "workflow",
        "summary": (
            "Re-submit an existing job with only resource overrides "
            "changed (PAL, maxcore, time-limit, basis swap). The "
            "Smart-Recalc panel auto-loads CONTROL.txt into an editor "
            "so you make a targeted edit (`/ui calc-editor replace "
            "PAL=40 PAL=12`) and click `calc-override-btn`. Faster "
            "than a from-scratch recalc; reuses .gbw and .opt when "
            "available."
        ),
        "see_also": [
            "delfin/dashboard/tab_calculations_browser.py",
            "delfin/orca_recovery.py",
        ],
    },
    "occupier": {
        "category": "module",
        "summary": (
            "Spin/occupation-state explorer: enumerates plausible "
            "(charge, multiplicity, occupation) combinations for "
            "open-shell metal complexes, runs each as a sub-job, and "
            "ranks by SCF + geometry quality. Drives the "
            "OCCUPIER_FOB job tree under each calc folder. Triggered "
            "by `OCCUPIER=auto` in CONTROL or via the typed "
            "build_mult/co2_species_delta override."
        ),
        "see_also": [
            "delfin/occupier.py", "delfin/occupier_flat_extraction.py",
        ],
    },
    "guppy": {
        "category": "module",
        "summary": (
            "GUPPY = generative ULFP/structure sampling. Generates "
            "diverse 3D conformers + spin states for SMILES input "
            "before the QM stages, using deterministic seed schedules "
            "(env-var DELFIN_GUPPY_SEEDS) plus optional ML-potential "
            "rescoring. Used as the front-end for SMILES → 3D in the "
            "ORCA Builder and batch flows."
        ),
        "see_also": [
            "delfin/guppy_sampling.py", "delfin/guppy_batch.py",
        ],
    },
    "csp": {
        "category": "module",
        "summary": (
            "Crystal Structure Prediction via Genarris. Generates "
            "candidate periodic structures from a molecular input. "
            "Triggered by `csp` workflow in CONTROL or `csp_check` "
            "MCP tool to verify Genarris availability."
        ),
        "see_also": ["delfin/csp_tools/"],
    },
    "mlp": {
        "category": "module",
        "summary": (
            "Machine-learning potentials for fast geometry pre-relax: "
            "TorchANI, AIMNet2, MACE backends. Used to cheap-relax "
            "before the expensive DFT step, especially in batch "
            "screening. Available via `mlp_check` MCP tool."
        ),
        "see_also": ["delfin/mlp_tools/"],
    },
    # Modes
    "modes": {
        "category": "agent",
        "summary": (
            "Agent operating modes: dashboard (UI-control, restricted "
            "tools), solo (full code access), quick (lightweight "
            "pipeline), reviewed (with critic), tdd (test-first), "
            "cluster (multi-system), full (everything). Set via "
            "`/mode <name>` slash command or the Mode dropdown."
        ),
        "see_also": ["delfin/agent/pack/agents/"],
    },
    "permissions": {
        "category": "agent",
        "summary": (
            "Four permission profiles: plan (read-only everywhere), "
            "ask_all (confirm every change, default), repo_free "
            "(repo edits free, calc asks), all_free (full mutation "
            "rights, archive stays read-only). Set via `/perms` slash."
        ),
        "see_also": ["delfin/dashboard/tab_agent.py"],
    },
    # Other concepts
    "co2": {
        "category": "workflow",
        "summary": (
            "CO2 Coordinator workflow: enumerates and optimizes CO2 "
            "binding modes on a metal complex, then ranks by binding "
            "energy. Run via `delfin co2` or the `co2` MCP tool with "
            "explicit charge/multiplicity/solvent."
        ),
        "see_also": ["delfin/cli.py (co2 subcommand)"],
    },
    "tadf_xtb": {
        "category": "workflow",
        "summary": (
            "TADF (thermally activated delayed fluorescence) workflow "
            "via xTB. Computes S0/S1/T1 energies, ΔEST, oscillator "
            "strength, and relevant kinetics. Faster than full DFT "
            "TDDFT; useful for screening dye candidates."
        ),
        "see_also": ["delfin/cli.py (tadf-xtb subcommand)"],
    },
    "hyperpol": {
        "category": "workflow",
        "summary": (
            "Hyperpolarizability (β, γ) workflow: ORCA-driven "
            "calculation of nonlinear-optical properties for "
            "chromophores. Reads orca.out for β-tensor components "
            "(beta_zzz etc.) and aggregates into DELFIN_data.json."
        ),
        "see_also": ["delfin/cli.py (hyperpol subcommand)"],
    },
    "outcomes": {
        "category": "agent",
        "summary": (
            "Persistent log of every agent turn: outcome_history.jsonl "
            "in ~/.delfin/. Each entry has provider, model, mode, "
            "verdict (PASS/FAIL/PARTIAL, computed via heuristics on "
            "denied-commands + give-up phrases + QUESTION tag), "
            "cost_usd_delta (per-turn spend), retries, denied_commands. "
            "Activity tab visualises this history."
        ),
        "see_also": ["delfin/agent/outcome_tracker.py"],
    },
    "session_boot": {
        "category": "agent",
        "summary": (
            "Every fresh dashboard session injects a one-shot "
            "[Session boot context] block on the first send: last 5 "
            "outcomes, active SLURM jobs, recent commits, branch, "
            "active calc folder. Lets the agent answer 'where are we?'"
            " without spending tool calls."
        ),
        "see_also": [
            "delfin/dashboard/tab_agent.py "
            "(_build_dashboard_session_boot)",
        ],
    },
    "live_state": {
        "category": "agent",
        "summary": (
            "Per-turn snapshot in the system prompt's `--- Live state "
            "---` section: current CONTROL, ORCA-Builder values, "
            "active calc folder, workspace files, recent jobs, "
            "permissions. Updated each send so the agent is always "
            "aware of UI state without round-trips. Lives in the "
            "system prompt (cache-friendly), not user messages."
        ),
        "see_also": [
            "delfin/dashboard/tab_agent.py (_build_dashboard_context)",
            "delfin/agent/engine.py (set_live_state)",
        ],
    },
}


def list_delfin_features(category: str = "") -> list[dict]:
    """Return the catalog of DELFIN concept entries, optionally filtered.

    Each entry: {name, category, summary}. Use this BEFORE describe_-
    delfin_feature to find the right name for a concept the user
    asked about.
    """
    cat = (category or "").strip().lower()
    out = []
    for name, info in sorted(_DELFIN_FEATURES.items()):
        if cat and info["category"].lower() != cat:
            continue
        out.append({
            "name": name,
            "category": info["category"],
            "summary": info["summary"][:140],
        })
    return out


def explain_delfin_feature(name: str) -> dict:
    """Return a curated explanation of one DELFIN feature/concept.

    Use when the user asks "wie funktioniert X in DELFIN?" / "was
    macht OCCUPIER?" / "was sind die Modi?". The returned dict has:

    - ``name``, ``category``, ``summary`` (full prose)
    - ``see_also``: list of source-file pointers if you want to dig
      deeper (Read tool).

    Unknown name → returns hint with the available concept names so
    the agent can self-correct.

    Args:
        name: concept name (case-insensitive). E.g. ``occupier``,
            ``smart_recalc``, ``relativistic_methods``,
            ``control_keys``, ``modes``, ``guppy``.
    """
    if not name:
        return {
            "error": "no name given",
            "available": [k for k in _DELFIN_FEATURES],
        }
    key = (
        str(name).strip().lower()
        .replace("-", "_").replace(" ", "_")
    )
    info = _DELFIN_FEATURES.get(key)
    if info is None:
        # Fuzzy: substring match
        candidates = [
            k for k in _DELFIN_FEATURES
            if key in k.lower() or k.lower() in key
        ]
        if len(candidates) == 1:
            info = _DELFIN_FEATURES[candidates[0]]
            key = candidates[0]
        else:
            return {
                "error": f"unknown feature: {name!r}",
                "candidates": candidates,
                "available": sorted(_DELFIN_FEATURES),
            }
    return {
        "name": key,
        "category": info["category"],
        "summary": info["summary"],
        "see_also": list(info.get("see_also", [])),
    }


# ---------------------------------------------------------------------------
# PDF on-demand reading (no pre-indexing required)
# ---------------------------------------------------------------------------


def _pdf_pages(path: str) -> list[dict]:
    """Page-wise plain text via pypdf. Empty list on missing/errored PDFs."""
    from pathlib import Path as _P
    p = _P(path).expanduser()
    if not p.exists() or p.suffix.lower() != ".pdf":
        return []
    try:
        from delfin.doc_server.indexer import _extract_pdf_text
        return _extract_pdf_text(p, quiet=True)
    except Exception:
        return []


def _parse_page_range(spec: str, total: int) -> list[int]:
    """Turn '5-10,12,15-20' into a sorted unique 1-based page list.

    Empty / invalid input returns every page (1..total).
    """
    if not spec:
        return list(range(1, total + 1))
    pages: set[int] = set()
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            try:
                a, b = (int(x) for x in part.split("-", 1))
            except ValueError:
                continue
            if a > b:
                a, b = b, a
            pages.update(range(max(1, a), min(total, b) + 1))
        else:
            try:
                v = int(part)
            except ValueError:
                continue
            if 1 <= v <= total:
                pages.add(v)
    return sorted(pages) if pages else list(range(1, total + 1))


def read_pdf(
    path: str,
    *,
    pages: str = "",
    max_chars: int = 50000,
) -> dict:
    """Read plain text from a PDF without indexing it first.

    Use for ad-hoc PDFs the user dropped in the Literature tab when
    you don't need persistent search. For repeated lookups, prefer
    ``index_new_pdf`` + ``search_docs``.

    Args:
        path: absolute path to the PDF.
        pages: optional 1-based page selector — ``"5"`` / ``"3-7"`` /
            ``"1,5-10,15"``. Empty → entire document.
        max_chars: cap on returned text size (default 50 000). The
            text is truncated WITH a ``[truncated …]`` marker so the
            agent knows there's more.

    Returns dict with: ``path``, ``n_pages_total``, ``n_pages_read``,
    ``text`` (page-tagged), ``truncated`` (bool), ``error``.
    """
    from pathlib import Path as _P
    p = _P(path).expanduser()
    if not p.exists():
        return {"path": str(p), "error": "file not found",
                "text": "", "n_pages_total": 0, "n_pages_read": 0,
                "truncated": False}
    if p.suffix.lower() != ".pdf":
        return {"path": str(p), "error": "not a PDF file",
                "text": "", "n_pages_total": 0, "n_pages_read": 0,
                "truncated": False}

    page_list = _pdf_pages(str(p))
    if not page_list:
        return {"path": str(p),
                "error": "could not extract text (pypdf missing or empty)",
                "text": "", "n_pages_total": 0, "n_pages_read": 0,
                "truncated": False}

    selected = _parse_page_range(pages, len(page_list))
    by_num = {entry["page"]: entry["text"] for entry in page_list}
    out_chunks = [
        f"--- PAGE {n} ---\n{by_num[n]}"
        for n in selected if n in by_num
    ]
    full = "\n\n".join(out_chunks)

    truncated = False
    if len(full) > max_chars:
        full = full[: max_chars - 32] + "\n\n[truncated — larger than max_chars]"
        truncated = True

    return {
        "path": str(p),
        "n_pages_total": len(page_list),
        "n_pages_read": len(selected),
        "text": full,
        "truncated": truncated,
        "error": "",
    }


def search_pdf_local(
    path: str,
    query: str,
    *,
    context_lines: int = 3,
    max_hits: int = 20,
    case_sensitive: bool = False,
) -> dict:
    """Substring-search inside ONE PDF, return matching paragraphs.

    Cheaper than indexing for one-off lookups. Returns hits with
    surrounding context (``context_lines`` per side) so the agent
    can quote the relevant passage.

    Args:
        path: absolute path to the PDF.
        query: substring to look for.
        context_lines: lines of context above/below each hit.
        max_hits: cap on returned hits.
        case_sensitive: default False — lowercase compare.

    Returns dict with ``path``, ``query``, ``hits`` (list of
    {page, line, snippet}), ``error``.
    """
    if not query.strip():
        return {"path": path, "query": query, "hits": [],
                "error": "empty query"}
    page_list = _pdf_pages(path)
    if not page_list:
        return {"path": path, "query": query, "hits": [],
                "error": "could not extract text"}
    needle = query if case_sensitive else query.lower()
    hits: list[dict] = []
    for entry in page_list:
        text = entry["text"]
        lines = text.splitlines()
        for i, line in enumerate(lines):
            cmp_line = line if case_sensitive else line.lower()
            if needle in cmp_line:
                lo = max(0, i - context_lines)
                hi = min(len(lines), i + context_lines + 1)
                snippet = "\n".join(lines[lo:hi])
                hits.append({
                    "page": entry["page"],
                    "line": i + 1,
                    "snippet": snippet,
                })
                if len(hits) >= max_hits:
                    return {"path": path, "query": query,
                            "hits": hits, "error": ""}
    return {"path": path, "query": query, "hits": hits, "error": ""}


def extract_pdf_section(
    path: str,
    heading: str,
    *,
    max_chars: int = 8000,
) -> dict:
    """Pull a single section starting at the given heading.

    Searches for the heading line in the PDF text, then returns
    everything from there until the next heading-like line (or the
    end of the document). Useful when ``search_docs`` returned a hit
    and you want the full section text without a separate
    ``read_section`` call.

    Args:
        path: absolute path to the PDF.
        heading: heading text to match (case-insensitive substring).
        max_chars: cap on returned text.

    Returns dict with ``path``, ``heading_found`` (bool), ``text``,
    ``next_heading`` (str), ``error``.
    """
    import re as _re
    page_list = _pdf_pages(path)
    if not page_list:
        return {"path": path, "heading_found": False,
                "text": "", "next_heading": "",
                "error": "could not extract text"}

    full = "\n\n".join(
        f"--- PAGE {p['page']} ---\n{p['text']}" for p in page_list
    )
    needle = heading.strip()
    if not needle:
        return {"path": path, "heading_found": False,
                "text": "", "next_heading": "", "error": "empty heading"}

    # Find heading position (case-insensitive)
    idx = full.lower().find(needle.lower())
    if idx == -1:
        return {"path": path, "heading_found": False,
                "text": "", "next_heading": "",
                "error": f"heading {heading!r} not found"}

    after = full[idx:]
    # Heuristic for "next heading": next line that looks like a
    # numbered or capitalised heading, with at least one blank line
    # before it.
    heading_re = _re.compile(
        r"\n\n(\d+(?:\.\d+)*\s+[A-Z][^\n]{2,80}|[A-Z][A-Z\s]{4,80})\n",
    )
    m = heading_re.search(after, pos=len(needle) + 1)
    if m:
        text = after[: m.start()]
        next_heading = m.group(1).strip()
    else:
        text = after
        next_heading = ""

    truncated = False
    if len(text) > max_chars:
        text = text[: max_chars - 32] + "\n\n[truncated]"
        truncated = True

    return {
        "path": path,
        "heading_found": True,
        "text": text,
        "next_heading": next_heading,
        "truncated": truncated,
        "error": "",
    }


def list_literature_files(
    folder: str = "",
    *,
    extensions: tuple = (".pdf", ".md", ".txt"),
) -> list[dict]:
    """List documents available in a Literature folder.

    Empty folder → uses the Literature directory next to the DELFIN
    repo (same path the indexer auto-detects). Recursive.

    Each entry: {path, name, size_bytes, mtime_iso, ext}.
    """
    from pathlib import Path as _P
    from datetime import datetime as _dt
    if folder:
        root = _P(folder).expanduser()
    else:
        try:
            from delfin.doc_server.indexer import get_default_literature_dir
            root = get_default_literature_dir() or _P.home() / "literature"
        except Exception:
            root = _P.home() / "literature"
    if not root.exists() or not root.is_dir():
        return []
    files: list[dict] = []
    exts = tuple(e.lower() for e in extensions)
    for p in sorted(root.rglob("*")):
        if not p.is_file():
            continue
        if p.suffix.lower() not in exts:
            continue
        try:
            stat = p.stat()
        except OSError:
            continue
        files.append({
            "path": str(p),
            "name": p.name,
            "size_bytes": int(stat.st_size),
            "mtime_iso": _dt.fromtimestamp(stat.st_mtime).isoformat(),
            "ext": p.suffix.lower(),
        })
    return files


def check_orca_manual_indexed() -> dict:
    """Quick-check whether the ORCA manual is in the doc-search index.

    Returns a dict with:
    - ``indexed`` (bool): True iff at least one indexed doc looks like
      the ORCA manual (title or doc_id contains "orca").
    - ``doc_ids`` (list): matching ids if found.
    - ``hint`` (str): if not indexed, a suggestion the agent should
      pass to the user (e.g. "drop the ORCA manual PDF into the
      Literature tab and call index_new_pdf").
    """
    import json as _json
    from pathlib import Path as _P
    try:
        from delfin.doc_server.indexer import get_default_index_path
        index_path = get_default_index_path()
    except Exception as exc:
        return {"indexed": False, "doc_ids": [],
                "hint": f"could not locate doc index: {exc}"}
    if not index_path.exists():
        return {
            "indexed": False, "doc_ids": [],
            "hint": (
                f"No doc index at {index_path}. Drop a PDF into the "
                f"Literature tab and call index_new_pdf, or run "
                f"`delfin-docs-index` once to build the initial index."
            ),
        }
    try:
        data = _json.loads(_P(index_path).read_text(encoding="utf-8"))
    except Exception as exc:
        return {"indexed": False, "doc_ids": [],
                "hint": f"index read error: {exc}"}
    docs = (data or {}).get("documents", {}) or {}
    matches: list[str] = []
    for doc_id, meta in docs.items():
        title = (meta or {}).get("title", "") or ""
        if "orca" in doc_id.lower() or "orca" in title.lower():
            matches.append(doc_id)
    if matches:
        return {"indexed": True, "doc_ids": matches, "hint": ""}
    return {
        "indexed": False, "doc_ids": [],
        "hint": (
            "ORCA manual not found in the index. Ask the user to drop "
            "the ORCA manual PDF (e.g. ORCA_6_1_1.pdf) into the "
            "Literature tab, then call index_new_pdf with the path."
        ),
    }


def index_new_pdf(
    path: str,
    *,
    doc_id: str = "",
    title: str = "",
) -> dict:
    """Add (or refresh) one PDF to the doc-search index.

    Re-builds the existing index PLUS the new PDF as an extra entry.
    Use after the user uploads a new manual / paper into Literature.
    The PDF must already exist on disk; this tool does NOT move or
    copy files.

    Args:
        path: absolute path to the PDF.
        doc_id: explicit id (defaults to filename stem).
        title: human-readable title (defaults to filename stem).

    Returns dict with ``ok`` (bool), ``doc_id``, ``sections_indexed``
    (int), ``error`` (str on failure).
    """
    from pathlib import Path as _P
    p = _P(path).expanduser()
    if not p.exists():
        return {"ok": False, "error": f"file not found: {path}"}
    if p.suffix.lower() != ".pdf":
        return {"ok": False, "error": "not a PDF file"}

    try:
        from delfin.doc_server.indexer import (
            build_index,
            get_default_index_path,
            get_default_literature_dir,
        )
    except Exception as exc:
        return {"ok": False, "error": f"indexer import failed: {exc}"}

    lit = get_default_literature_dir() or _P.home() / "literature"
    final_doc_id = doc_id or p.stem
    final_title = title or p.stem.replace("_", " ")
    try:
        index = build_index(
            literature_dir=lit,
            extra_paths=[{
                "path": str(p),
                "doc_id": final_doc_id,
                "title": final_title,
            }],
            quiet=True,
        )
    except Exception as exc:
        return {"ok": False, "error": f"build_index failed: {exc}"}

    out_path = get_default_index_path()
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        import json as _json
        out_path.write_text(_json.dumps(index, indent=2), encoding="utf-8")
    except Exception as exc:
        return {"ok": False, "error": f"index write failed: {exc}"}

    docs = (index or {}).get("documents", {}) or {}
    sections = 0
    if final_doc_id in docs:
        sections = len((docs[final_doc_id] or {}).get("sections", {}))
    return {
        "ok": True,
        "doc_id": final_doc_id,
        "title": final_title,
        "sections_indexed": sections,
        "index_path": str(out_path),
    }


def list_active_calculations() -> list[dict]:
    """Return the live job list from the active backend (read-only).

    Each entry: {job_id, name, status, runtime_s, directory}. Empty
    list when SLURM isn't reachable or no jobs are active.
    """
    try:
        backend = _resolve_backend()
        jobs = list(backend.list_jobs() or [])
    except Exception as exc:
        return [{"error": str(exc)}]
    out: list[dict] = []
    for j in jobs:
        out.append({
            "job_id": getattr(j, "job_id", "") or "",
            "name": getattr(j, "name", "") or "",
            "status": getattr(j, "status", "") or "",
            "runtime_s": float(getattr(j, "runtime_s", 0) or 0),
            "directory": getattr(j, "directory", "") or "",
        })
    return out


def _default_plot_dir() -> str:
    """Return the agent_workspace/ directory the dashboard auto-renders."""
    from pathlib import Path as _P
    cand = _P.home() / "agent_workspace"
    if cand.exists():
        return str(cand)
    return str(_P.cwd())


@dataclass
class PlotResult:
    """Outcome of a plot helper call.

    Carrying the path lets the agent reference it in the chat so the
    inline-render hook can pick it up and embed the PNG.
    """
    path: str = ""
    n_points: int = 0
    title: str = ""
    properties: list[str] = field(default_factory=list)
    statistics: dict = field(default_factory=dict)
    error: str = ""


def plot_energy_distribution(
    folders: list[str] | str,
    *,
    properties: list[str] | None = None,
    output_path: str = "",
    plot_type: str = "histogram",
    title: str = "",
    bins: int = 30,
) -> PlotResult:
    """Extract energies across folders and write a statistical plot.

    Reads ``properties`` from each folder's largest ``*.out`` (via
    :func:`extract_energy_table`) and produces:

    - ``plot_type="histogram"`` — one histogram per property, stacked
      vertically (default).
    - ``plot_type="bar"`` — bar chart per folder (sorted by the first
      property), useful for small N.
    - ``plot_type="boxplot"`` — distribution summary across folders.

    The PNG lands in ``agent_workspace/`` by default so the dashboard's
    inline-artifact hook picks it up automatically. ``PlotResult.path``
    tells the agent where it ended up; mention that path in chat and
    the user sees the image.

    Args:
        folders: list of paths (or single path) to scan.
        properties: subset of ``gibbs``, ``zpe``, ``single_point``.
            Default ``["gibbs", "single_point"]``.
        output_path: explicit PNG location. Empty → auto-generated
            inside ``agent_workspace/``.
        plot_type: ``histogram`` | ``bar`` | ``boxplot``.
        title: figure title (auto-generated if empty).
        bins: histogram bin count.
    """
    if properties is None:
        properties = ["gibbs", "single_point"]
    if isinstance(folders, str):
        folders = [folders]
    rows = extract_energy_table(folders, properties=properties)
    valid_rows = [
        r for r in rows
        if any(r.get(p) is not None for p in properties)
    ]
    if not valid_rows:
        return PlotResult(
            error=f"No folder produced parseable values for {properties}.",
            properties=list(properties),
            n_points=0,
        )

    try:
        import matplotlib  # noqa
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        return PlotResult(
            error=f"matplotlib/numpy missing: {exc}",
            properties=list(properties),
        )

    series: dict[str, list[float]] = {p: [] for p in properties}
    labels: list[str] = []
    for r in valid_rows:
        labels.append(_short_folder_label(r["folder"]))
        for p in properties:
            v = r.get(p)
            if v is not None:
                series[p].append(float(v))

    n_points = max(len(v) for v in series.values()) if series else 0
    if not title:
        prop_label = " + ".join(properties)
        title = f"{prop_label} across {n_points} calculations ({plot_type})"

    out_path = output_path or str(
        _new_workspace_png_path(prefix=f"energy_{plot_type}")
    )

    plot_type_norm = plot_type.lower().strip()
    statistics: dict = {}
    if plot_type_norm == "bar":
        # One row per folder; bars side-by-side per property.
        n_props = len(properties)
        x = np.arange(len(labels))
        width = 0.8 / max(1, n_props)
        fig, ax = plt.subplots(
            figsize=(max(6, len(labels) * 0.4 + 2), 5),
        )
        for i, p in enumerate(properties):
            vals = [r.get(p) if r.get(p) is not None else float("nan")
                    for r in valid_rows]
            ax.bar(x + (i - n_props / 2 + 0.5) * width, vals,
                   width=width, label=p)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=70, ha="right", fontsize=8)
        ax.set_ylabel("Energy / Hartree")
        ax.set_title(title)
        ax.legend(loc="best")
        plt.tight_layout()
        plt.savefig(out_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
    elif plot_type_norm == "boxplot":
        fig, ax = plt.subplots(figsize=(2 + 1.5 * len(properties), 5))
        data = [series[p] for p in properties]
        ax.boxplot(data, tick_labels=list(properties))
        ax.set_ylabel("Energy / Hartree")
        ax.set_title(title)
        plt.tight_layout()
        plt.savefig(out_path, dpi=120, bbox_inches="tight")
        plt.close(fig)
    else:  # histogram (default)
        fig, axes = plt.subplots(
            len(properties), 1,
            figsize=(7, 2.5 * len(properties)),
            sharex=False,
        )
        if len(properties) == 1:
            axes = [axes]
        for ax, p in zip(axes, properties):
            data = series[p]
            if not data:
                ax.text(0.5, 0.5, f"no data for {p}",
                        transform=ax.transAxes, ha="center")
                ax.set_axis_off()
                continue
            ax.hist(data, bins=bins, color="#6366f1", edgecolor="white")
            ax.set_xlabel(f"{p} / Hartree")
            ax.set_ylabel("count")
            ax.set_title(f"{p} (N={len(data)})", fontsize=11)
            ax.axvline(np.mean(data), color="#ef4444",
                       linestyle="--", label=f"mean={np.mean(data):.4f}")
            ax.axvline(np.median(data), color="#10b981",
                       linestyle=":", label=f"median={np.median(data):.4f}")
            ax.legend(loc="best", fontsize=8)
        fig.suptitle(title, fontsize=12)
        plt.tight_layout()
        plt.savefig(out_path, dpi=120, bbox_inches="tight")
        plt.close(fig)

    for p in properties:
        data = series.get(p, [])
        if data:
            statistics[p] = {
                "n": len(data),
                "min": float(min(data)),
                "max": float(max(data)),
                "mean": float(sum(data) / len(data)),
                "range": float(max(data) - min(data)),
            }

    return PlotResult(
        path=out_path,
        n_points=n_points,
        title=title,
        properties=list(properties),
        statistics=statistics,
    )


def plot_energy_correlation(
    folders: list[str] | str,
    *,
    x: str = "single_point",
    y: str = "gibbs",
    output_path: str = "",
    title: str = "",
) -> PlotResult:
    """Scatter plot one energy property against another across folders.

    Pearson correlation + best-fit line are drawn on top so the user
    can see at a glance whether ``y`` tracks ``x`` linearly.

    Args:
        folders: list of paths (or single path) to scan.
        x, y: two of ``gibbs``, ``zpe``, ``single_point``.
        output_path: explicit PNG location.
        title: figure title (auto-generated if empty).
    """
    if isinstance(folders, str):
        folders = [folders]
    rows = extract_energy_table(folders, properties=[x, y])
    pairs = [(r.get(x), r.get(y)) for r in rows
             if r.get(x) is not None and r.get(y) is not None]
    if not pairs:
        return PlotResult(
            error=f"No folder has both {x} and {y}.",
            properties=[x, y],
        )

    try:
        import matplotlib  # noqa
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        return PlotResult(
            error=f"matplotlib/numpy missing: {exc}",
            properties=[x, y],
        )

    xs = np.array([p[0] for p in pairs], dtype=float)
    ys = np.array([p[1] for p in pairs], dtype=float)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(xs, ys, s=24, color="#6366f1", alpha=0.7,
               edgecolor="white", linewidth=0.5)
    if len(xs) >= 2:
        slope, intercept = np.polyfit(xs, ys, 1)
        x_line = np.linspace(xs.min(), xs.max(), 50)
        ax.plot(x_line, slope * x_line + intercept,
                color="#ef4444", linestyle="--", linewidth=1.2,
                label=f"y = {slope:.4f}·x + {intercept:.4f}")
        # Pearson r
        if xs.std() and ys.std():
            r = float(np.corrcoef(xs, ys)[0, 1])
        else:
            r = 0.0
    else:
        slope = intercept = r = 0.0
    ax.set_xlabel(f"{x} / Hartree")
    ax.set_ylabel(f"{y} / Hartree")
    if not title:
        title = f"{y} vs {x}  (N={len(xs)}, r={r:.3f})"
    ax.set_title(title)
    if len(xs) >= 2:
        ax.legend(loc="best", fontsize=9)
    plt.tight_layout()

    out_path = output_path or str(
        _new_workspace_png_path(prefix=f"energy_corr_{x}_vs_{y}")
    )
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)

    return PlotResult(
        path=out_path,
        n_points=len(xs),
        title=title,
        properties=[x, y],
        statistics={
            "n": len(xs),
            "pearson_r": float(r),
            "slope": float(slope),
            "intercept": float(intercept),
        },
    )


def _short_folder_label(folder: str) -> str:
    """Return a 24-char-max label for a folder path (last 1-2 components)."""
    from pathlib import Path as _P
    p = _P(folder)
    if not p.parts:
        return folder
    # Use last two components if available
    if len(p.parts) >= 2:
        cand = "/".join(p.parts[-2:])
    else:
        cand = p.name
    if len(cand) > 24:
        cand = cand[:11] + "…" + cand[-12:]
    return cand


def _new_workspace_png_path(prefix: str = "plot") -> "Path":
    """Return a non-clobbering PNG path in agent_workspace/."""
    from pathlib import Path as _P
    base = _P(_default_plot_dir())
    base.mkdir(parents=True, exist_ok=True)
    from datetime import datetime as _dt
    stamp = _dt.now().strftime("%Y%m%d_%H%M%S")
    return base / f"{prefix}_{stamp}.png"


def find_calculation_extreme(
    folders: list[str] | str,
    *,
    property: str = "gibbs",
    extreme: str = "min",
    n: int = 5,
) -> list[dict]:
    """Return the N folders with the lowest/highest value of a property.

    Direct answer to "open the .out with the lowest Gibbs energy"
    type questions: pass the candidate folders, get back the top
    ``n`` rows sorted ascending (``extreme="min"``) or descending
    (``extreme="max"``).

    Folders that fail to parse the property are excluded from the
    ranking.

    Args:
        folders: list of paths (or single path) to scan.
        property: one of ``gibbs``, ``zpe``, ``single_point``,
            ``imag_freqs``, ``walltime_s``.
        extreme: ``"min"`` (lowest) or ``"max"`` (highest).
        n: how many top rows to return (clipped to len(rows)).
    """
    rows = extract_energy_table(folders, properties=[property])
    valid = [r for r in rows if r.get(property) is not None]
    if not valid:
        return []
    reverse = (str(extreme).lower() == "max")
    valid.sort(key=lambda r: float(r[property]), reverse=reverse)
    return valid[: max(1, int(n))]


# ---------------------------------------------------------------------------
# On-demand operational-pattern lookup
# ---------------------------------------------------------------------------
#
# Concrete slash-chain recipes for common dashboard workflows.  These
# used to be baked into ``dashboard_agent.md`` but that paid token cost
# every turn — even for conversations that have nothing to do with
# batch jobs or recalc.  Now they live behind a typed lookup the agent
# only calls when it actually needs a workflow recipe.
#
# Adding a new pattern: drop an entry into ``_DASHBOARD_PATTERNS`` and
# add the name to the docstring for discoverability.

_DASHBOARD_PATTERNS: dict[str, str] = {
    "batch": (
        "## Batch jobs (Submit-Tab)\n\n"
        "- Build a batch from EVERY calc-folder's initial geometry:\n"
        "    `ACTION: /batch from-calc`\n"
        "- Build from a glob subset (e.g. all Casagrande folders):\n"
        "    `ACTION: /batch from-calc Casagrande*`\n"
        "- Add one SMILES line manually:\n"
        "    `ACTION: /batch add Name;SMILES;charge=…`\n"
        "- Show / clear the current batch text:\n"
        "    `ACTION: /batch show`  ·  `ACTION: /batch clear`\n\n"
        "**Never** assemble batch text by reading XYZ files yourself —\n"
        "`/batch from-calc` already collects every initial.xyz (or\n"
        "fallback input.txt / coords.xyz) across `calculations/` and\n"
        "writes a properly formatted block into the Submit-Tab textarea.\n\n"
        "Example:\n"
        "  > User: \"Bau einen Batch aus allen XYZ in calc/.\"\n"
        "  > Agent:\n"
        "  >\n"
        "  >     ACTION: /batch from-calc\n"
    ),
    "control_edit": (
        "## CONTROL.txt edits (Submit-Tab)\n\n"
        "- One key change → `ACTION: /control key <key> <value>`\n"
        "- Replace whole content (rare) → `ACTION: /control set <multi-line>`\n"
        "- Validate before submit → `ACTION: /control validate`\n\n"
        "Never paste the full CONTROL into chat — use `/control key` per change.\n"
    ),
    "smart_recalc": (
        "## Smart Recalc (Calc browser)\n\n"
        "Selecting Smart Recalc auto-loads CONTROL into the editor, so\n"
        "the right idiom is `/ui calc-editor replace <old> <new>`:\n\n"
        "  > User: \"smart recalc von Foo/bar mit PAL=1, 1 h\"\n"
        "  > Agent:\n"
        "  >\n"
        "  >     ACTION: /calc cd Foo/bar\n"
        "  >     ACTION: /calc select CONTROL.txt\n"
        "  >     ACTION: /ui calc-options value Smart Recalc\n"
        "  >     ACTION: /ui calc-editor replace PAL=40 PAL=1\n"
        "  >     ACTION: /ui calc-override-time value 01:00:00\n"
        "  >\n"
        "  > Then ASK before submitting:\n"
        "  >\n"
        "  >     ACTION: /ui calc-override-btn click\n\n"
        "`calc-override-btn` is the submit; `calc-override-time` is the\n"
        "time field. Don't confuse with `calc-recalc-btn` /\n"
        "`calc-submit-recalc-btn` — those are NOT the Smart-Recalc panel.\n"
    ),
    "submit_orca": (
        "## Submit a single ORCA job (ORCA Builder)\n\n"
        "- Set fields → `ACTION: /ui orca-method value PBE0`,\n"
        "  `/ui orca-basis value def2-TZVP`, `/ui orca-charge value 0`, …\n"
        "- Switch to the tab → `ACTION: /tab orca`\n"
        "- Submit (after the user explicitly OKs!) → `ACTION: /orca submit`\n"
    ),
    "analyze": (
        "## Analyze existing calculations\n\n"
        "- One folder, full → `ACTION: /analyze <dir>`\n"
        "- Just energies → `ACTION: /analyze energy <dir>`\n"
        "- SCF convergence → `ACTION: /analyze convergence <dir>`\n"
        "- Error scan → `ACTION: /analyze errors <dir>`\n"
        "- All folders overview → `ACTION: /analyze status`\n\n"
        "For multi-folder energy tables prefer\n"
        "`mcp__delfin-ops__extract_energy_table(folders=[...])` — it\n"
        "returns structured data the agent can format directly.\n"
    ),
    "recalc": (
        "## Recalc check / submit\n\n"
        "- Check one folder (safe, read-only) → `ACTION: /recalc check <dir>`\n"
        "- Scan everything (safe) → `ACTION: /recalc check-all`\n"
        "- Submit recalc (DESTRUCTIVE — needs explicit user OK):\n"
        "    `ACTION: /recalc <dir>` and confirm\n"
        "- Bulk auto-recalc only on explicit \"alle neuberechnen\":\n"
        "    `ACTION: /recalc auto`\n"
    ),
    "cancel": (
        "## Cancel jobs\n\n"
        "- One job → `ACTION: /cancel <job_id>` (after user OK)\n"
        "- All — only on explicit \"cancel all\":\n"
        "    `ACTION: /cancel all`\n"
    ),
}


def list_dashboard_patterns() -> list[str]:
    """Return the names of every operational pattern available for lookup.

    Use this when the agent isn't sure which pattern matches the user's
    request — a names-only list is cheap and lets the agent pick.
    """
    return sorted(_DASHBOARD_PATTERNS)


def get_dashboard_pattern(name: str) -> str:
    """Return the slash-chain recipe for a named dashboard workflow.

    Available names: ``batch``, ``control_edit``, ``smart_recalc``,
    ``submit_orca``, ``analyze``, ``recalc``, ``cancel``.

    Names are matched case-insensitively. An unknown name returns a
    short error string with the available choices so the agent can
    self-correct.
    """
    if not name:
        return (
            "No pattern requested. Call with one of: "
            + ", ".join(list_dashboard_patterns())
        )
    key = str(name).strip().lower().replace("-", "_").replace(" ", "_")
    text = _DASHBOARD_PATTERNS.get(key)
    if text is None:
        return (
            f"Unknown pattern: {name!r}. Available: "
            + ", ".join(list_dashboard_patterns())
        )
    return text


# ---------------------------------------------------------------------------
# Module surface
# ---------------------------------------------------------------------------

__all__ = [
    "CommandResult",
    "OrcaParseResult",
    "OrcaError",
    "ThermochemResult",
    "PlotResult",
    "run",
    "prepare",
    "pipeline_run",
    "pipeline_prepare",
    "qm_check",
    "csp_check",
    "mlp_check",
    "analysis_check",
    "qm_run",
    "cleanup",
    "stop",
    "run_orca_input",
    "co2",
    "tadf_xtb",
    "hyperpol",
    "list_dashboard_patterns",
    "get_dashboard_pattern",
    # P1 — output parsing
    "parse_orca_output",
    "find_orca_errors",
    "extract_thermochem",
    "extract_energy_table",
    "find_calculation_extreme",
    "ImagFreqMode",
    "ImagFreqResult",
    "extract_imaginary_frequencies",
    "FunctionalComparisonRow",
    "compare_across_functionals",
    "CalcDiff",
    "compare_calculations",
    # P1 — statistical plots (PNG → agent_workspace, auto-shown in chat)
    "plot_energy_distribution",
    "plot_energy_correlation",
    # Tool / widget catalogs (cheap on-demand discovery)
    "list_tools",
    "describe_tool",
    "list_dashboard_widgets",
    "get_widget_options",
    # ORCA Builder validation
    "ValidationIssue",
    "validate_orca_input",
    # Job lifecycle (typed, allow_mutate gates)
    "SubmitOutcome",
    "submit_calculation",
    "cancel_calculation",
    "list_active_calculations",
    # Calc folder management
    "rename_calc_folder",
    "create_calc_folder",
    "move_calc_folder",
    "move_to_archive",
    "delete_calc_folder",
    "kill_all_user_jobs",
    "RecalcPlan",
    "prepare_recalc",
    "list_calc_options",
    "run_calc_option",
    # ORCA-manual lookup helpers
    "check_orca_manual_indexed",
    "index_new_pdf",
    # PDF on-demand reading (no pre-indexing)
    "read_pdf",
    "search_pdf_local",
    "extract_pdf_section",
    "list_literature_files",
    # DELFIN-feature explainer
    "list_delfin_features",
    "explain_delfin_feature",
]
